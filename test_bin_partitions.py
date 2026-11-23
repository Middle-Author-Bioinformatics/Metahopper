import argparse
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import metahopper_report as report
from metahopper_taxonomy import assembly_partition, is_unbinned, taxonomy_summary
from metahopper_report_figures import assembly_context
import metahopper_postprocess as postprocess
import MetaHopper as pipeline


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


class BinPartitionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name) / 'Sample'
        self.sequences = {'a': 'A' * 1000, 'b': 'C' * 300, 'c': 'G' * 700, 'extra': 'T' * 200}
        self.assembly = write(self.root / 'final/assembly/consolidated_contigs.fasta',
                              ''.join('>' + k + '\n' + v + '\n' for k, v in self.sequences.items()))
        self.bin('genus', 'Sodalis', ['a'])
        self.bin('genus', 'unclassified_Bacteria', ['b'])
        self.bin('genus', 'Unclassified', ['c'])
        self.bin('phylum', 'Pseudomonadota', ['a', 'b'])
        self.bin('phylum', 'Unclassified', ['c'])
        self.layout = report.RunLayout(self.root, 'final')

    def tearDown(self):
        self.tmp.cleanup()

    def bin(self, rank, name, ids):
        return write(self.root / 'final/bins' / rank / (name + '.fasta'),
                     ''.join('>' + k + '\n' + self.sequences[k] + '\n' for k in ids))

    def rows(self):
        return report.collect_bins(self.layout, True)

    def partition(self):
        return assembly_partition(self.layout, self.rows(), report.iter_fasta, report.read_tsv)

    def test_each_rank_conserves_full_assembly_and_missing_unbinned_contigs_are_included(self):
        data = self.partition()
        self.assertEqual(data['total_bp'], 2200)
        for rank in ('phylum', 'genus'):
            parts = data['ranks'][rank]['parts']
            self.assertEqual(sum(p['bp'] for p in parts), 2200)
            self.assertEqual(parts[-1], dict(bin='Unbinned', bp=900, contigs=2, unbinned=True))
        self.assertIn('unclassified_Bacteria', [p['bin'] for p in data['ranks']['genus']['parts']])

    def test_assembly_scope_describes_the_actual_fasta_without_changing_totals(self):
        data = self.partition()
        html = report.render_bins(self.rows(), report.rank_totals(self.rows()), partition=data)
        self.assertIn('id="assembly-scope-note"', html)
        self.assertIn('Removed host and other excluded contigs are not shown.', html)
        self.assertIn('"assembly_label":"Retained assembly"', html)
        self.assertEqual(data['total_bp'], 2200)
        self.assertNotIn('assembly_label', data)
        for path in ('megahit/final.contigs.fa', 'assembly/preliminary_contigs.fasta'):
            label, note = assembly_context(dict(source=path))
            self.assertEqual(label, 'Initial assembly')
            self.assertIn('Before taxonomic filtering', note)
        self.assertEqual(assembly_context(dict(source='unknown.fasta'))[0], 'Displayed assembly')
        self.assertEqual(assembly_context(None)[0], 'Displayed assembly')

    def test_generic_catchalls_hidden_but_domain_assignments_remain_and_quality_is_preserved(self):
        for name in ('Unclassified', 'unbinned', 'Unassigned'):
            self.assertTrue(is_unbinned(name))
        for name in ('unclassified_Bacteria', 'unclassified_Viruses', 'unclassified_Archaea'):
            self.assertFalse(is_unbinned(name))
        write(self.root / 'final/checkm.tsv', 'Bin Id\tCompleteness\tContamination\n'
                                            'genus__Sodalis\t95.2\t1.3\n')
        rows = self.rows()
        html = report.render_bins(rows, report.rank_totals(rows), partition=self.partition())
        self.assertNotIn('data-bin-name="Unclassified"', html)
        self.assertIn('data-bin-name="unclassified_Bacteria"', html)
        self.assertIn('>95.20<', html)
        self.assertIn('>1.30<', html)
        for removed in ('bin-quality-figure', 'bin-size-figure', 'class="stats"', 'Refinement compares preliminary'):
            self.assertNotIn(removed, html)
        self.assertIn('xlog:true,ylog:true', html)
        self.assertIn("p.unbinned?'#000000'", html)
        index = taxonomy_summary('Sample', 'final', rows)
        self.assertFalse(any(b['bin'] == 'Unclassified' for b in index['bins']))
        self.assertEqual(next(b['completeness_percent'] for b in index['bins'] if b['bin'] == 'Sodalis'), 95.2)

    def test_overlap_and_sequence_mismatch_fail_instead_of_double_counting(self):
        self.bin('genus', 'Overlap', ['a'])
        self.assertIn('Overlapping', self.partition()['ranks']['genus']['error'])
        (self.root / 'final/bins/genus/Overlap.fasta').unlink()
        write(self.root / 'final/bins/genus/Sodalis.fasta', '>a\nAAAA\n')
        self.assertIn('do not match', self.partition()['ranks']['genus']['error'])

    def test_membership_is_verified_against_summary_and_legacy_assembly_export_is_accepted(self):
        write(self.root / 'final/bins/genus/bin_membership.tsv', 'contig\tbin\na\tSodalis\n')
        self.assertEqual(self.partition()['ranks']['genus']['parts'][0]['bp'], 1000)
        self.assembly.rename(self.assembly.with_name('contigs.fasta'))
        self.layout = report.RunLayout(self.root, 'final')
        self.assertEqual(self.partition()['total_bp'], 2200)


class CheckmOnlyTests(unittest.TestCase):
    def test_quality_only_preflight_does_not_require_bakta_or_pathwaygenie(self):
        result = argparse.Namespace(stdout='--reduced_tree --pplacer_threads --tab_table')
        with patch.object(postprocess.annotation, 'validate_config') as bakta, \
             patch.object(postprocess.pathways, 'preflight_runtime') as pathways, \
             patch.object(postprocess.subprocess, 'run', return_value=result) as command:
            config = {'tools': {'checkm': ['conda', 'run', '-n', 'checkm_env', 'checkm']}}
            postprocess.validate(config, checkm_only=True)
            bakta.assert_not_called(); pathways.assert_not_called()
            self.assertEqual(command.call_args.args[0][-2:], ['lineage_wf', '-h'])

    def test_existing_bins_are_unchanged_other_rank_scores_survive_and_success_is_reused(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder) / 'Sample'
            fasta = write(root / 'final/bins/genus/Sodalis.fasta', '>a\n' + 'A' * 1000 + '\n')
            write(root / 'final/assembly/consolidated_contigs.fasta', fasta.read_text())
            for relative in ('classification/contig_classification.tsv', 'binarena/binarena_input.tsv'):
                write(root / 'final' / relative, 'contig\na\n')
            summary = root / 'final/bins/genus/summary.tsv'
            pipeline.tsv(summary, [dict(bin='Sodalis', **report.fasta_stats(fasta))], pipeline.tb.SUMMARY_FIELDS)
            config = write(root / 'config.json', '{}')
            options = pipeline.parse_args(['-1', 'reads.fq', '-o', str(root), '-c', str(config)])
            write(root / 'run.json', json.dumps(dict(status='completed', arguments=vars(options), config={}), default=str))
            write(root / 'final/checkm.tsv', 'Bin Id\tCompleteness\tContamination\tStrain heterogeneity\n'
                                            'phylum__Pseudomonadota\t99\t20\t0\n')
            before = postprocess.sequence_snapshot(root)
            calls = []
            def checkm(tool, arguments, label):
                self.assertEqual(tool, 'checkm')
                calls.append(tool)
                arguments = list(map(str, arguments))
                output = Path(arguments[arguments.index('-f') + 1])
                write(output, 'Bin Id\tCompleteness\tContamination\tStrain heterogeneity\n'
                              'genus__Sodalis\t95\t1\t0\n')
            with patch.object(postprocess, 'validate'), \
                 patch.object(pipeline.Runner, 'run', side_effect=checkm), \
                 patch.object(postprocess, 'prepare_report', side_effect=lambda *a: write(root / 'metahopper_report.html', 'updated')) as refresh, \
                 patch.object(postprocess.annotation, 'annotate_run') as bakta, \
                 patch.object(postprocess.pathways, 'run_pathways') as pathways:
                postprocess.run(root, config, checkm_only=True)
                postprocess.run(root, config, checkm_only=True)
                postprocess.run(root, config, checkm_only=True, refresh_report=False)
            self.assertEqual(calls, ['checkm'])
            bakta.assert_not_called(); pathways.assert_not_called()
            self.assertEqual(refresh.call_count, 2)
            self.assertEqual(before, postprocess.sequence_snapshot(root))
            table = report.read_tsv(root / 'final/checkm.tsv')
            self.assertEqual({r['Bin Id'] for r in table}, {'phylum__Pseudomonadota', 'genus__Sodalis'})
            self.assertEqual(report.read_tsv(summary)[0]['completeness_percent'], '95')
            self.assertEqual(json.loads((root / 'postprocess_status.json').read_text())['status'], 'completed')


if __name__ == '__main__':
    unittest.main()
