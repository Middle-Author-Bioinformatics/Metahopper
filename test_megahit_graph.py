"""MEGAHIT graph orientation, source scope, atomic writes and assembly reuse."""
import argparse
import contextlib
import io
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import MetaHopper as pipeline
import metahopper_graph as graph
import metahopper_report as report
from metahopper_state import Checkpoints


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def assembly(path, records):
    return write(path, ''.join(f'>{name}\n{seq}\n' for name, seq in records.items()))


class MegahitGraphTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.source = self.root / 'megahit/final.contigs.fa'
        self.dest = self.root / 'megahit/assembly.gfa'
        self.manifest = self.root / 'megahit/assembly_graph.json'

    def tearDown(self):
        self.tmp.cleanup()

    def export(self, records, **kwargs):
        assembly(self.source, records)
        result = graph.export_megahit_graph(self.root, **kwargs)
        data, nodes, _, _ = graph.load_gfa(self.dest, 'test')
        return result, data, nodes

    def test_end_links_equal_brute_force_oriented_overlaps_and_keep_ids(self):
        records = {'k3_1': 'AACCGT', 'k3_2': 'CGTTAA', 'k3_3': 'TACGGA',
                   'k3_4': 'CCGTAA', 'k5_1': 'AACCGTTAA', 'k3_5': 'NNNAAA', 'k3_6': 'GGGACG'}
        result, data, nodes = self.export(records)
        flip = {'+': '-', '-': '+'}
        expected = set()
        for left, seq in records.items():
            for right, target in records.items():
                k = int(left.split('_')[0][1:])
                if right.split('_')[0] != left.split('_')[0]:
                    continue
                for a, aseq in (('+', seq), ('-', seq.translate(graph.COMPLEMENT)[::-1])):
                    for b, bseq in (('+', target), ('-', target.translate(graph.COMPLEMENT)[::-1])):
                        end = aseq[-k:]
                        if len(aseq) > k and len(bseq) > k and set(end) <= set('ACGT') and end == bseq[:k]:
                            edge = (left, a, right, b)
                            opposite = (right, flip[b], left, flip[a])
                            expected.add(min(edge, opposite))
        actual = {(x['source'], x['source_orient'], x['target'], x['target_orient']) for x in data['links']}
        self.assertEqual(actual, expected)
        self.assertTrue(any(a[1] != a[3] for a in actual))
        self.assertEqual(set(nodes), set(records))
        self.assertEqual(result['links'], len(expected))
        self.assertEqual(data['construction'], 'megahit_kmer_overlaps')
        self.assertEqual(data['source_scope'], 'retained_contigs')

    def test_intermediates_keep_short_neighbours_and_cache_survives_cleanup(self):
        records = {'k3_1': 'AACCGT', 'k1_1': 'AT'}
        assembly(self.source, records)
        directory = self.root / 'work/megahit/intermediate_contigs'
        assembly(directory / 'k3.contigs.fa', {'k3_1': 'AACCGT', 'k3_2': 'CGTA'})
        assembly(directory / 'k1.final.contigs.fa', {'k1_1': 'AT'})
        result = graph.export_megahit_graph(self.root)
        self.assertEqual(result['source_scope'], 'intermediates')
        self.assertEqual(result['segments'], 3)
        for path in directory.iterdir():
            path.unlink()
        before = self.dest.stat().st_mtime_ns
        self.assertTrue(graph.export_megahit_graph(self.root)['cached'])
        self.assertEqual(self.dest.stat().st_mtime_ns, before)

    def test_original_final_contigs_map_by_recorded_id_and_graph_scope_is_visible(self):
        self.export({'k3_1': 'AACCGT', 'k3_2': 'CGTTAA'})
        write(self.root / 'final/assembly/contig_provenance.tsv',
              'final_contig\tsource_rank\tsource_bin\tsource_type\tsource_record\n'
              'MH_original_1\tgenus\tTest\toriginal\tk3_1\n')
        catalog = graph.build_catalog(self.root)
        self.assertEqual(catalog['contigs']['MH_original_1'][0]['steps'], [['k3_1', '+']])
        self.assertIn('retained contigs only', ' '.join(catalog['warnings']))

    def test_modified_source_or_graph_invalidates_export_cache(self):
        self.export({'k3_1': 'AACCGT'})
        write(self.dest, 'broken')
        self.assertFalse(graph.export_megahit_graph(self.root)['cached'])
        assembly(self.source, {'k3_2': 'CGTTAA'})
        result = graph.export_megahit_graph(self.root)
        self.assertFalse(result['cached'])
        self.assertEqual(set(graph.load_gfa(self.dest, 'test')[1]), {'k3_2'})

    def test_link_limit_and_bad_ids_preserve_prior_complete_graph(self):
        self.export({'k3_1': 'AACCGT'})
        before = (self.dest.read_bytes(), self.manifest.read_bytes())
        assembly(self.source, {'k3_1': 'AAAA', 'k3_2': 'AAAA'})
        with self.assertRaisesRegex(ValueError, 'exceeds'):
            graph.export_megahit_graph(self.root, max_links=1)
        self.assertEqual((self.dest.read_bytes(), self.manifest.read_bytes()), before)
        assembly(self.source, {'renamed_contig': 'AACCGT'})
        with self.assertRaisesRegex(ValueError, 'Cannot infer'):
            graph.export_megahit_graph(self.root)
        self.assertEqual((self.dest.read_bytes(), self.manifest.read_bytes()), before)
        self.assertFalse(list(self.dest.parent.glob('.megahit_graph_*')))

    def test_manual_graph_is_not_overwritten_without_force(self):
        assembly(self.source, {'k3_1': 'AACCGT'})
        write(self.dest, 'manual graph')
        with self.assertRaisesRegex(ValueError, 'already exists'):
            graph.export_megahit_graph(self.root)
        self.assertEqual(self.dest.read_text(), 'manual graph')
        self.assertFalse(graph.export_megahit_graph(self.root, force=True)['cached'])

    def test_duplicate_final_ids_are_rejected_and_intermediate_conflicts_are_namespaced(self):
        write(self.source, '>k3_1\nAACCGT\n>k3_1\nCGTTAA\n')
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            graph.export_megahit_graph(self.root)
        assembly(self.source, {'k3_1': 'AACCGT'})
        assembly(self.root / 'work/megahit/intermediate_contigs/k3.contigs.fa', {'k3_1': 'CGTTAA'})
        result = graph.export_megahit_graph(self.root)
        self.assertEqual(result['renamed_intermediates'], 1)
        nodes = graph.load_gfa(self.dest, 'test')[1]
        self.assertIn('k3_1', nodes)
        self.assertEqual(len(nodes), 2)

    def test_existing_run_cli_exports_a_loadable_catalogue_without_external_tools(self):
        assembly(self.source, {'k3_1': 'AACCGT', 'k3_2': 'CGTTAA'})
        completed = subprocess.run([sys.executable, str(Path(graph.__file__)), '-i', str(self.root),
                                    '--bin-set', 'preliminary', '--export-megahit-graph'],
                                   text=True, capture_output=True, check=True)
        self.assertIn('2 segments', completed.stdout)
        catalog = json.loads((self.root / 'assembly_graphs.preliminary.json').read_text())
        self.assertEqual(catalog['counts']['matched_contigs'], 2)
        self.assertEqual(catalog['graphs'][0]['source_scope'], 'retained_contigs')

    def test_report_embeds_exported_original_graph_and_its_mapping(self):
        self.export({'k3_1': 'AACCGT', 'k3_2': 'CGTTAA'})
        assembly(self.root / 'final/assembly/consolidated_contigs.fasta', {'MH_original_1': 'AACCGT'})
        assembly(self.root / 'final/bins/genus/Test.fasta', {'MH_original_1': 'AACCGT'})
        write(self.root / 'final/assembly/contig_provenance.tsv',
              'final_contig\tsource_rank\tsource_bin\tsource_type\tsource_record\n'
              'MH_original_1\tgenus\tTest\toriginal\tk3_1\n')
        write(self.root / 'final/binarena/binarena_input.tsv',
              'contig\tlength\tGC\tcoverage\tbin\tkmer_4_PC1\tkmer_4_PC2\n'
              'MH_original_1\t6\t50\t2\tTest\t0.1\t0.2\n')
        args = argparse.Namespace(no_fasta_stats=False, min_length=1, max_points=100,
                                  include_uncoordinated=False, title=None)
        html = report.build_report(report.RunLayout(self.root, 'final'), args)
        payload = json.loads(re.search(r'var DATA = (.*?);\s*/\*', html, re.S)[1])
        data = payload['assemblyGraphs']
        self.assertEqual(data['counts']['matched_contigs'], 1)
        self.assertEqual(data['contigs']['MH_original_1'][0]['steps'], [['k3_1', '+']])
        self.assertEqual(data['graphs'][0]['construction'], 'megahit_kmer_overlaps')


class AssemblyGraphIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.reads = [write(self.root / 'reads.fq', '@r\nACGT\n+\nIIII\n')]
        self.out = self.root / 'run'
        self.ck = Checkpoints(self.out, log=lambda _: None)
        self.args = argparse.Namespace(assembly=None, assembly_graph=None, min_contig_len=20,
                                       threads=1, memory=1, force_overwrite=False)
        self.calls = 0
        outer = self
        class FakeRunner:
            def run(self, tool, arguments, label):
                outer.calls += 1
                args = list(map(str, arguments))
                directory = Path(args[args.index('-o') + 1])
                records = {'k15_1': 'A'*15 + 'T' + 'C'*15}
                assembly(directory / 'final.contigs.fa', records)
                assembly(directory / 'intermediate_contigs/k15.contigs.fa',
                         dict(records, k15_2='C'*15 + 'G'))
        self.runner = FakeRunner()

    def tearDown(self):
        self.tmp.cleanup()

    def run_stage(self):
        with contextlib.redirect_stdout(io.StringIO()):
            return pipeline.stage_assembly(self.runner, self.args, self.out, self.reads, self.ck)

    def test_fresh_run_exports_before_cleanup_and_resume_does_not_assemble_again(self):
        self.run_stage()
        self.assertEqual(self.calls, 1)
        result = json.loads((self.out / 'megahit/assembly_graph.json').read_text())
        self.assertEqual(result['segments'], 2)
        self.assertEqual(result['source_scope'], 'intermediates')
        self.assertFalse((self.out / 'work/megahit').exists())
        self.run_stage()
        self.assertEqual(self.calls, 1)

    def test_export_failure_preserves_intermediates_and_retries_without_reassembly(self):
        with patch.object(graph, 'export_megahit_graph', side_effect=ValueError('simulated')):
            self.run_stage()
        self.assertTrue((self.out / 'work/megahit/intermediate_contigs/k15.contigs.fa').is_file())
        self.assertEqual(self.ck.status['steps']['megahit_graph']['status'], 'failed')
        self.run_stage()
        self.assertEqual(self.calls, 1)
        self.assertFalse((self.out / 'work/megahit').exists())

    def test_old_completed_assembly_backfills_graph_without_rerunning_megahit(self):
        self.run_stage()
        for filename in ('assembly.gfa', 'assembly_graph.json'):
            (self.out / 'megahit' / filename).unlink()
        self.run_stage()
        self.assertEqual(self.calls, 1)
        result = json.loads((self.out / 'megahit/assembly_graph.json').read_text())
        self.assertEqual(result['source_scope'], 'retained_contigs')
        self.assertEqual(result['segments'], 1)


if __name__ == '__main__':
    unittest.main()
