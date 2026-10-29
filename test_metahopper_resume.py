"""Completion/retry regressions using tiny files and simulated external tools."""
import argparse
from collections import Counter
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import metahopper as pipeline
from metahopper_state import Checkpoints, output_lock


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


class StateTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.root = Path(self.temp.name)
        self.source = write(self.root / 'input', 'input')
        self.output = write(self.root / 'output', 'output')

    def tearDown(self): self.temp.cleanup()

    def state(self, **kwargs): return Checkpoints(self.root, log=lambda _: None, **kwargs)

    def complete(self, state):
        self.assertIsNone(state.get('test', inputs=[self.source], required=[self.output]))
        state.done('test', [self.output])

    def test_only_matching_success_with_intact_outputs_is_reused(self):
        state = self.state(); self.complete(state)
        self.assertIsNotNone(self.state().get('test', inputs=[self.source], required=[self.output]))
        write(self.output, 'truncated')
        self.assertIsNone(self.state().get('test', inputs=[self.source], required=[self.output]))

    def test_changed_inputs_options_and_force_invalidate_completion(self):
        self.complete(self.state())
        self.assertIsNone(self.state().peek('test', inputs=[self.source], options={'new': True}))
        self.assertIsNone(self.state(force=True).peek('test', inputs=[self.source]))
        write(self.source, 'different input')
        self.assertIsNone(self.state().peek('test', inputs=[self.source]))

    def test_running_and_failed_steps_are_retryable_and_visible_at_root(self):
        state = self.state(); state.get('test', inputs=[self.source])
        self.assertEqual(json.loads((self.root / 'pipeline_status.json').read_text())['current_step'], 'test')
        state = self.state()
        self.assertEqual(state.status['steps']['test']['status'], 'interrupted')
        self.assertIsNone(state.get('test', inputs=[self.source]))
        state.fail(error='simulated failure')
        self.assertEqual(state.status['steps']['test']['status'], 'failed')
        self.assertIsNone(self.state().peek('test', inputs=[self.source]))
        self.assertIn('simulated failure', (self.root / 'pipeline.log').read_text())

    def test_corrupt_marker_and_directory_alone_never_prove_completion(self):
        self.assertIsNone(self.state().peek('test', inputs=[self.source], required=[self.output]))
        write(self.root / 'work/checkpoints/test.json', 'broken json')
        self.assertIsNone(self.state().peek('test', inputs=[self.source], required=[self.output]))

    def test_an_active_output_directory_cannot_be_written_by_a_second_run(self):
        with output_lock(self.root):
            with self.assertRaises(RuntimeError), output_lock(self.root): pass
        with output_lock(self.root): pass

    def test_old_success_marker_requires_saved_compatible_context(self):
        marker = write(self.root / 'work/checkpoints/test.json', json.dumps({'files': [str(self.output)]}))
        self.assertIsNone(self.state().peek('test', inputs=[self.source]))
        state = self.state(contexts={'test': {'mode': 1}}, previous_contexts={'test': {'mode': 1}})
        self.assertIsNotNone(state.get('test', inputs=[self.source]))
        self.assertIn('signature', json.loads(marker.read_text()))

    def test_inspection_reuses_success_and_force_bypasses_the_tree_cache(self):
        write(self.root / 'final/assembly/consolidated_contigs.fasta', '>c1\n' + 'A'*1000 + '\n')
        write(self.root / 'final/classification/contig_classification.tsv', 'contig\nc1\n')
        calls = []
        def inspect(root, **kwargs):
            calls.append(kwargs['force'])
            result = dict(errors=[], trees=[])
            write(root / 'inspection/inspection.json', json.dumps(result))
            return result
        module = argparse.Namespace(inspect_run=inspect)
        spec = argparse.Namespace(loader=argparse.Namespace(exec_module=lambda _: None))
        args = argparse.Namespace(inspect='all', references=None, threads=1, force_overwrite=False)
        with patch('importlib.util.spec_from_file_location', return_value=spec), \
                patch('importlib.util.module_from_spec', return_value=module), contextlib.redirect_stdout(io.StringIO()):
            state = self.state()
            pipeline.stage_inspect(args, self.root, 'genus', [self.source], state)
            pipeline.stage_inspect(args, self.root, 'genus', [self.source], state)
            args.force_overwrite = True
            pipeline.stage_inspect(args, self.root, 'genus', [self.source], self.state(force=True))
        self.assertEqual(calls, [False, True])


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.root = Path(self.temp.name)
        self.out = self.root / 'run'
        self.raw = write(self.root / 'reads.fastq', '@read\nACGT\n+\nIIII\n')
        self.assembly = write(self.root / 'assembly.fasta', '>c1\n' + 'A' * 1000 + '\n>c2\n' + 'G' * 1000 + '\n')
        database = write(self.root / 'nr.dmnd', 'fake database')
        taxdump = self.root / 'taxdump'
        write(taxdump / 'nodes.dmp', 'fixture'); write(taxdump / 'names.dmp', 'fixture')
        self.config = write(self.root / 'config.json', json.dumps({'diamond_db': str(database), 'taxdump': str(taxdump)}))
        self.tokens = ['-1', str(self.raw), '--assembly', str(self.assembly), '-o', str(self.out), '-c', str(self.config),
                       '--prefilter', 'none', '--skip-phyloflash', '--inspect', 'off', '--ranks', 'phylum,genus',
                       '--min-bin-length', '1000', '--checkm']
        self.counts = Counter(); self.fail_tool = None; self.fail_label = None
        self.patches = [patch.object(pipeline, 'validate'), patch.object(pipeline, 'Taxonomy', return_value=object()),
                        patch.object(pipeline.Runner, 'run', self.run_tool),
                        patch.object(pipeline, 'stream_mapping', side_effect=self.mapping),
                        patch.object(pipeline.tb, 'classify_hits', side_effect=self.votes),
                        patch.object(pipeline.tb, 'ordinate', side_effect=self.ordinate)]
        for item in self.patches: item.start()

    def tearDown(self):
        for item in reversed(self.patches): item.stop()
        self.temp.cleanup()

    def execute(self, extra=(), resume=False):
        tokens = ['--resume', str(self.out)] if resume else self.tokens
        with contextlib.redirect_stdout(io.StringIO()): pipeline.cli(tokens + list(extra))

    def run_tool(self, tool, arguments, label, **kwargs):
        self.counts[tool] += 1
        if tool == self.fail_tool or label == self.fail_label:
            self.fail_tool = None; self.fail_label = None
            raise RuntimeError('Simulated ' + tool + ' failure')
        args = list(map(str, arguments))
        def value(flag): return Path(args[args.index(flag) + 1])
        if tool == 'fastp':
            for flag in ('-o', '--json', '--html'): write(value(flag), 'fixture\n')
        elif tool == 'diamond' and args[0] == 'blastp':
            write(value('--out'), 'fixture hits\n')
        elif tool == 'prodigal':
            seqs = pipeline.load_sequences(value('-i'))
            write(value('-a'), ''.join(f'>{cid}_1 # 1 # 900 # 1 # ID=x\nMMMM\n' for cid in seqs))
            write(value('-o'), '##gff-version 3\n')
        elif tool == 'unicycler':
            sequence = 'A' * 1000 + 'G' * 1000
            write(value('-o') / 'assembly.fasta', '>1\n' + sequence + '\n')
            write(value('-o') / 'assembly.gfa', 'S\t1\t' + sequence + '\n')
        elif tool == 'minimap2':
            source = pipeline.load_sequences(Path(args[-1]))
            write(kwargs['stdout'], ''.join(f'{cid}\t{len(seq)}\t0\t{len(seq)}\t+\ttarget\t2000\t0\t{len(seq)}\t{len(seq)}\t{len(seq)}\t60\n'
                                            for cid, seq in source.items()))
        elif tool == 'checkm':
            ids = sorted(path.stem for path in Path(args[-2]).glob('*.fasta'))
            write(value('-f'), 'Bin Id\tCompleteness\tContamination\tStrain heterogeneity\n' +
                  ''.join(f'{cid}\t95\t1\t0\n' for cid in ids))

    def mapping(self, runner, assembly, reads, threads, label, **kwargs):
        self.counts[label] += 1
        seqs = pipeline.load_sequences(assembly)
        pools = kwargs.get('pools')
        if pools:
            for name in set(kwargs['contig_bins'].values()):
                pools.counts[name] = 4; pools.bases[name] = 4000
                write(pools.directory / (name + '.SE.fq.gz'), 'fixture reads')
        return ({cid: 2.0 for cid in seqs}, {}, {cid: 1.0 for cid in seqs})

    def votes(self, *args):
        self.counts['votes'] += 1
        names = dict(domain='Bacteria', kingdom='Unclassified', phylum='Pseudomonadota',
                     **{'class': 'Gammaproteobacteria'}, order='Enterobacterales', family='Enterobacteriaceae',
                     genus='Sodalis', species='Sodalis_test')
        return ({cid: {'calls': {rank: (names[rank], 1.0) for rank in pipeline.tb.ORDER}, 'orfs_with_hits': 1}
                 for cid in ('c1', 'c2')}, {'hits': 2})

    def ordinate(self, seqs, *args, **kwargs):
        self.counts['ordination'] += 1
        return {'4PC1': {cid: float(i) for i, cid in enumerate(seqs)},
                '4PC2': {cid: float(i) for i, cid in enumerate(seqs)}}

    def test_existing_output_auto_reuses_heavy_steps_including_final_and_checkm(self):
        self.execute(); first = self.counts.copy()
        self.assertEqual(first['unicycler'], 1); self.assertEqual(first['final_mapping'], 1)
        self.assertFalse(list((self.out / 'work/read_pools').glob('*.fq.gz')))
        self.execute()
        self.assertEqual(self.counts, first)
        state = json.loads((self.out / 'pipeline_status.json').read_text())
        self.assertEqual(state['status'], 'completed')
        self.assertEqual(state['steps']['final']['status'], 'reused')
        self.assertEqual(state['steps']['checkm']['status'], 'reused')
        self.assertIn('95.00', (self.out / 'metahopper_report.html').read_text())

    def test_force_overwrite_reruns_and_is_not_sticky_on_resume(self):
        self.execute(); first = self.counts.copy()
        self.execute(['--force-overwrite'], resume=True)
        for key in ('fastp', 'diamond', 'prodigal', 'unicycler', 'mapping', 'final_mapping', 'ordination', 'checkm', 'votes'):
            self.assertEqual(self.counts[key], first[key] * 2, key)
        forced = self.counts.copy(); self.execute(resume=True)
        self.assertEqual(self.counts, forced)

    def test_missing_final_output_repeats_final_work_without_repeating_diamond_or_reassembly(self):
        self.execute(); first = self.counts.copy()
        (self.out / 'final/binarena/binarena_input.tsv').unlink()
        self.execute()
        self.assertEqual(self.counts['diamond'], first['diamond'])
        self.assertEqual(self.counts['unicycler'], first['unicycler'])
        self.assertEqual(self.counts['final_mapping'], first['final_mapping'] + 1)

    def test_failed_checkm_is_retried_without_other_heavy_steps(self):
        self.fail_tool = 'checkm'; self.execute(); first = self.counts.copy()
        self.assertFalse((self.out / 'work/checkpoints/checkm.json').exists())
        self.execute()
        self.assertEqual(self.counts['checkm'], first['checkm'] + 1)
        self.assertEqual(self.counts['diamond'], first['diamond'])
        self.assertEqual(self.counts['final_mapping'], first['final_mapping'])

    def test_failed_diamond_records_failure_then_resumes_from_that_step(self):
        self.fail_label = 'preliminary_diamond'
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit): self.execute()
        state = json.loads((self.out / 'pipeline_status.json').read_text())
        self.assertEqual(state['status'], 'failed')
        self.assertEqual(state['steps']['preliminary_diamond']['status'], 'failed')
        first = self.counts.copy(); self.execute()
        self.assertEqual(self.counts['fastp'], first['fastp'])
        self.assertEqual(self.counts['prodigal'], first['prodigal'] + 1)  # final assembly ORFs only
        self.assertEqual(self.counts['diamond'], first['diamond'] + 1)

    def test_failed_candidate_keeps_pools_and_retries_instead_of_becoming_permanent_fallback(self):
        self.fail_tool = 'unicycler'; self.execute(); first = self.counts.copy()
        self.assertTrue(list((self.out / 'work/read_pools').glob('*.fq.gz')))
        self.execute()
        self.assertEqual(self.counts['mapping'], first['mapping'])
        self.assertEqual(self.counts['unicycler'], first['unicycler'] + 1)
        result = json.loads((self.out / 'reassembly/genus/Sodalis/result.json').read_text())
        self.assertEqual(result['outcome'], 'accepted')

    def test_unfinished_bin_with_missing_pool_regenerates_mapping(self):
        self.fail_tool = 'unicycler'; self.execute(); first = self.counts.copy()
        for path in (self.out / 'work/read_pools').glob('*.fq.gz'): path.unlink()
        self.execute()
        self.assertEqual(self.counts['mapping'], first['mapping'] + 1)
        self.assertEqual(self.counts['unicycler'], first['unicycler'] + 1)

    def test_force_is_passed_to_per_bin_annotation(self):
        def annotate(root, config, rank, *args, **kwargs):
            calls.append(kwargs['force'])
            summary = write(root / 'final/annotation' / rank / 'summary.tsv', 'bin\tstatus\nSodalis\tcompleted\n')
            return dict(bins=[dict(status='completed')], summary=str(summary), warnings=[])
        calls = []
        with patch('metahopper_annotate.annotate_run', side_effect=annotate):
            self.execute(['--annotate', 'bakta'])
            self.execute(['--force-overwrite'], resume=True)
        self.assertEqual(calls, [False, True])

    def test_changed_input_invalidates_downstream_results(self):
        self.execute(); first = self.counts.copy()
        self.assembly.write_text(self.assembly.read_text().replace('A' * 1000, 'T' * 1000))
        self.execute()
        self.assertEqual(self.counts['fastp'], first['fastp'])
        self.assertEqual(self.counts['diamond'], first['diamond'] + 1)
        self.assertEqual(self.counts['mapping'], first['mapping'] + 1)

    def test_old_complete_run_and_markers_can_be_adopted_without_repeating_alignments(self):
        self.execute(); first = self.counts.copy()
        for path in (self.out / 'work/checkpoints').glob('*.json'):
            if path.stem in ('dbinfo', 'qc', 'assembly', 'preliminary_prodigal', 'preliminary_diamond', 'mapping'):
                old = json.loads(path.read_text())
                path.write_text(json.dumps({'files': old['files'], 'time': old['time']}))
            else: path.unlink()
        (self.out / 'pipeline_status.json').unlink()
        (self.out / 'classification/votes.json.gz').unlink()
        self.execute()
        for key in ('fastp', 'diamond', 'prodigal', 'unicycler', 'mapping', 'final_mapping', 'ordination', 'checkm'):
            self.assertEqual(self.counts[key], first[key], key)

    def test_resume_accepts_explicit_new_options_and_new_defaults(self):
        self.execute()
        manifest = json.loads((self.out / 'run.json').read_text())
        manifest['arguments'].pop('force_overwrite', None)
        manifest['arguments'].pop('annotation_mode', None)
        (self.out / 'run.json').write_text(json.dumps(manifest))
        with patch.object(pipeline, 'main') as main:
            pipeline.cli(['--resume', str(self.out), '-t32', '--annotate', 'bakta'])
            args = main.call_args.args[0]
        self.assertEqual(args.threads, 32)
        self.assertEqual(args.annotate, 'bakta')
        self.assertEqual(args.annotation_mode, 'meta')
        self.assertFalse(args.force_overwrite)


if __name__ == '__main__': unittest.main()
