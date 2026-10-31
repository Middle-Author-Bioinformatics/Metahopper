"""Data-free regressions for graph provenance and per-bin annotation orchestration."""
import argparse
import json
import sys
import tempfile
import unittest
from pathlib import Path

import metahopper_graph as graph
import metahopper_annotate as ann
import metahopper_report as report


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def assembly(path, records):
    return write(path, ''.join('>'+cid+'\n'+seq+'\n' for cid, seq in records.items()))


class GraphTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def fixture(self, sequences, aliases='', kind='reassembled'):
        directory = self.root / 'reassembly/genus/Test'
        assembly(directory / 'assembly.fasta', sequences)
        if aliases:
            write(directory / 'contig_graph_aliases.tsv', 'contig\tgfa_record\n'+aliases)
        prov = 'final_contig\tsource_rank\tsource_bin\tsource_type\tsource_fasta\tsource_record\n'
        for i, cid in enumerate(sequences):
            prov += f'MH_test_{i}\tgenus\tTest\t{kind}\tunused.fasta\t{cid}\n'
        write(self.root / 'final/assembly/contig_provenance.tsv', prov)
        return directory

    def test_exact_sequence_and_reverse_complement_recover_old_names(self):
        directory = self.fixture({'renamed_a':'AACCGGTA', 'renamed_b':'TTAACCGG'})
        write(directory/'assembly.gfa', 'H\tVN:Z:1.0\nS\t1\tAACCGGTA\nS\t2\tCCGGTTAA\nL\t1\t+\t2\t-\t0M\n')
        data = graph.build_catalog(self.root)
        self.assertEqual(data['contigs']['MH_test_0'][0]['steps'], [['1','+']])
        self.assertEqual(data['contigs']['MH_test_1'][0]['steps'], [['2','-']])
        self.assertEqual(data['graphs'][0]['links'][0]['target_orient'], '-')
        self.assertNotIn('AACCGGTA', json.dumps(data))

    def test_duplicate_sequences_are_ambiguous_not_assigned_by_order(self):
        directory=self.fixture({'Test_reassembled_00001':'AACCGGTA'})
        write(directory/'assembly.gfa','S\t1\tAACCGGTA\nS\t2\tAACCGGTA\n')
        data=graph.build_catalog(self.root)
        self.assertNotIn('MH_test_0',data['contigs'])
        self.assertEqual(data['unmapped']['MH_test_0'],'ambiguous_exact_sequence')

    def test_preserved_alias_resolves_duplicate_sequences(self):
        directory=self.fixture({'renamed':'AACCGGTA'},'renamed\t2\n')
        write(directory/'assembly.gfa','S\t1\tAACCGGTA\nS\t2\tAACCGGTA\n')
        self.assertEqual(graph.build_catalog(self.root)['contigs']['MH_test_0'][0]['steps'],[['2','+']])

    def test_rejected_candidate_graph_is_not_used_for_original_fallback(self):
        directory=self.fixture({'original':'AACCGGTA'},kind='original')
        write(directory/'assembly.gfa','S\t1\tAACCGGTA\n')
        data=graph.build_catalog(self.root)
        self.assertFalse(data['contigs'])
        self.assertIn('preliminary_graph',data['unmapped']['MH_test_0'])

    def test_paths_walks_and_bidirected_link_deduplication(self):
        p=write(self.root/'test.gfa','H\tVN:Z:1.1\nS\ta\tAA\nS\tb\tGG\nL\ta\t+\tb\t-\t0M\nL\tb\t+\ta\t-\t0M\nP\tp\ta+,b-\t0M\nW\tsample\t0\tchr1\t0\t4\t>a<b\n')
        data,nodes,paths,_=graph.load_gfa(p,'test')
        self.assertEqual(len(data['links']),1)
        self.assertEqual(paths['p'],[['a','+'],['b','-']])
        self.assertEqual(paths['sample#0#chr1#0#4'],[['a','+'],['b','-']])

    def test_missing_nodes_and_gfa2_fail_explicitly(self):
        for text in ['H\tVN:Z:2.0\n','S\ta\tAA\nL\ta\t+\tmissing\t+\t0M\n']:
            p=write(self.root/'bad.gfa',text)
            with self.assertRaises(ValueError): graph.load_gfa(p,'bad')

    def test_graph_limit_does_not_silently_drop_segments(self):
        p=write(self.root/'large.gfa','S\ta\tAA\nS\tb\tGG\n')
        with self.assertRaisesRegex(ValueError,'exceeds'):
            graph.load_gfa(p,'large',max_nodes=1)

    def test_read_counts_are_not_misreported_as_depth(self):
        p=write(self.root/'counts.gfa','S\ta\tACGT\tRC:i:40\tKC:i:60\n')
        data,_,_,_=graph.load_gfa(p,'counts')
        self.assertIsNone(data['nodes'][0]['depth'])
        self.assertEqual(data['nodes'][0]['read_count'],40)
        self.assertEqual(data['nodes'][0]['kmer_count'],60)

    def test_preliminary_report_uses_preliminary_ids_even_with_final(self):
        assembly(self.root/'assembly/preliminary_contigs.fasta',{'old':'AAA'})
        assembly(self.root/'final/assembly/consolidated_contigs.fasta',{'new':'AAA'})
        (self.root/'final/bins').mkdir(parents=True)
        self.assertEqual(report.RunLayout(self.root,'preliminary').default_contigs(),self.root/'assembly/preliminary_contigs.fasta')


FAKE_BAKTA = '''import sys,json,pathlib
if '--version' in sys.argv:
 print('bakta-test 1.0');raise SystemExit(0)
args=sys.argv
output=pathlib.Path(args[args.index('--output')+1]);output.mkdir(parents=True)
with (output.parent/'executions.txt').open('a') as handle:handle.write('run\\n')
source=pathlib.Path(args[-1]);cid=source.read_text().splitlines()[0][1:]
if 'Fail' in source.stem: raise SystemExit(7)
payload={'version':{'bakta':'test'},'features':[{'type':'cds','contig':cid,'start':1,'stop':9,'strand':'+','locus':'MH_TEST_1','product':'test enzyme','db_xrefs':['EC:1.1.1.1']},{'type':'cds','contig':cid,'start':10,'stop':18,'strand':'+','locus':'MH_TEST_2','product':'hypothetical protein','hypothetical':True},{'type':'trna','contig':cid,'start':20,'stop':30}]}
(output/'annotation.json').write_text(json.dumps(payload))
for suffix in ['.gff3','.gbff','.faa','.ffn','.tsv']:(output/('annotation'+suffix)).write_text('test output\\n')
'''


class AnnotationTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.db=self.root/'db';self.db.mkdir();write(self.db/'version.json','{"major":5,"minor":0,"date":"2026-01-01","type":"full"}')
        tool=write(self.root/'fake_bakta.py',FAKE_BAKTA)
        self.config=dict(bakta_db=str(self.db),tools={'bakta':[sys.executable,str(tool)]})
        self.bin=self.root/'final/bins/genus/Sodalis.fasta';assembly(self.bin,{'contig_a':'A'*1200})
        self.classification=self.root/'final/classification/contig_classification.tsv'
        write(self.classification,'contig\tdomain\tgenus\tlength_bp\ncontig_a\tBacteria\tSodalis\t1200\n')

    def tearDown(self): self.tmp.cleanup()

    def run_annotation(self,**kw): return ann.annotate_run(self.root,self.config,log=lambda _:None,**kw)

    def test_force_runs_bakta_again_after_successful_cache_hit(self):
        self.run_annotation()
        self.assertTrue(self.run_annotation()['bins'][0]['cached'])
        self.assertFalse(self.run_annotation(force=True)['bins'][0]['cached'])
        executions=self.root/'final/annotation/genus/Sodalis/executions.txt'
        self.assertEqual(executions.read_text().splitlines(),['run','run'])

    def test_report_does_not_present_annotation_of_changed_bin_as_current(self):
        self.run_annotation()
        assembly(self.bin,{'new_contig':'G'*1300})
        html=report.render_annotation(report.RunLayout(self.root,'final'))
        self.assertIn('stale',html)
        self.assertIn('does not match the current bin',html)
        self.assertNotIn('>GenBank</a>',html)

    def test_outputs_identity_features_and_cache(self):
        result=self.run_annotation();row=result['bins'][0]
        self.assertEqual(row['status'],'completed');self.assertEqual(row['cds'],2)
        self.assertEqual(row['named_proteins'],1);self.assertEqual(row['hypothetical_proteins'],1)
        directory=self.root/'final/annotation/genus/Sodalis'
        features=graph.table(directory/'features.tsv');self.assertEqual(features[0]['contig'],'contig_a')
        self.assertEqual(features[0]['db_xrefs'],'EC:1.1.1.1')
        modules=json.loads((directory/'modules.json').read_text());self.assertEqual(modules['biosynthetic_clusters']['status'],'not_run')
        self.assertTrue(self.run_annotation()['bins'][0]['cached'])
        assembly(self.bin,{'contig_a':'A'*1201})
        self.assertFalse(self.run_annotation()['bins'][0]['cached'])

    def test_database_version_change_invalidates_cache(self):
        self.run_annotation();write(self.db/'version.json','{"major":5,"minor":1,"date":"2026-01-01","type":"full"}')
        self.assertFalse(self.run_annotation()['bins'][0]['cached'])

    def test_non_bacterial_and_domain_catchall_skipped(self):
        assembly(self.bin.parent/'unclassified_Bacteria.fasta',{'u':'A'*100})
        write(self.classification,'contig\tdomain\tgenus\ncontig_a\tArchaea\tSodalis\nu\tBacteria\tUnclassified\n')
        rows=self.run_annotation()['bins'];self.assertEqual({r['status'] for r in rows},{'skipped'})
        self.assertIn('bacterial',rows[0]['reason'])

    def test_code4_hint_and_explicit_metadata(self):
        self.bin.rename(self.bin.with_name('Spiroplasma.fasta'))
        write(self.classification,'contig\tdomain\tgenus\ncontig_a\tBacteria\tSpiroplasma\n')
        row=self.run_annotation()['bins'][0];self.assertEqual(row['translation_table'],4)
        meta=write(self.root/'annotation_metadata.tsv','bin\ttranslation_table\nSpiroplasma\t25\n')
        self.config['annotation_metadata']=str(meta)
        row=self.run_annotation()['bins'][0];self.assertEqual(row['translation_table'],25)
        self.assertEqual(row['translation_source'],'metadata')

    def test_per_bin_failure_does_not_erase_success(self):
        assembly(self.bin.parent/'Fail.fasta',{'f':'A'*1100})
        with open(self.classification,'a') as f:f.write('f\tBacteria\tFail\t1100\n')
        data=self.run_annotation();status={r['bin']:r['status'] for r in data['bins']}
        self.assertEqual(status,{'Fail':'failed','Sodalis':'completed'})
        self.assertEqual(len(data['warnings']),1)

    def test_managed_options_cannot_be_overridden(self):
        self.config['bakta_args']=['--output=/elsewhere']
        with self.assertRaisesRegex(ValueError,'managed'):ann.validate_config(self.config)

    def test_subset_preserves_other_successful_bins_and_flags_changed_inputs(self):
        second=self.bin.parent/'Other.fasta';assembly(second,{'other':'G'*1200})
        with open(self.classification,'a') as f:f.write('other\tBacteria\tOther\t1200\n')
        self.run_annotation()
        rows=self.run_annotation(selected_bins=['Sodalis'])['bins']
        self.assertEqual({r['bin'] for r in rows},{'Sodalis','Other'})
        self.assertEqual(next(r for r in rows if r['bin']=='Other')['status'],'completed')
        assembly(second,{'other':'G'*1300})
        rows=self.run_annotation(selected_bins=['Sodalis'])['bins']
        self.assertEqual(next(r for r in rows if r['bin']=='Other')['status'],'stale')

    def test_mixed_translation_code_hints_require_metadata(self):
        assembly(self.bin,{'contig_a':'A'*1200,'contig_b':'G'*1200})
        write(self.classification,'contig\tdomain\tgenus\ncontig_a\tBacteria\tSpiroplasma\ncontig_b\tBacteria\tSodalis\n')
        row=self.run_annotation()['bins'][0]
        self.assertEqual(row['status'],'skipped')
        self.assertIn('mixed translation',row['reason'])

    def test_changed_taxonomic_eligibility_does_not_leave_completed_module_status(self):
        self.run_annotation()
        write(self.classification,'contig\tdomain\tgenus\ncontig_a\tEukaryota\tSodalis\n')
        self.run_annotation()
        path=self.root/'final/annotation/genus/Sodalis/modules.json'
        self.assertEqual(json.loads(path.read_text())['bakta']['status'],'skipped')

    def test_annotation_report_has_real_results_and_no_fake_specialist_scores(self):
        self.run_annotation()
        html=report.render_annotation(report.RunLayout(self.root,'final'))
        self.assertIn('GenBank',html);self.assertIn('Specialized pathway, BGC',html)
        self.assertNotIn('virulence factors: 0',html)


class ReportTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        assembly(self.root/'final/bins/genus/Test.fasta',{'f1':'A'*1000,'f2':'G'*1000,'f3':'T'*1000})
        write(self.root/'final/bins/genus/summary.tsv','bin\tnum_contigs\ttotal_length_bp\tN50\tGC_percent\tcompleteness_percent\tcontamination_percent\nTest\t3\t3000\t1000\t33.3\tNA\tNA\n')
        write(self.root/'bins/genus/summary.tsv','bin\tnum_contigs\ttotal_length_bp\tN50\nTest\t8\t2800\t500\n')
        write(self.root/'bins/genus/reassembly_summary.tsv','bin\toutcome\tcontigs_before\ttotal_length_before_bp\tcontigs_after\ttotal_length_after_bp\nTest\taccepted\t8\t2800\t5\t3100\n')
        self.layout=report.RunLayout(self.root,'final')

    def tearDown(self): self.tmp.cleanup()

    def test_compact_refinement_uses_actual_final_metrics(self):
        rows=report.collect_bins(self.layout,True)
        html=report.render_bins(rows,report.rank_totals(rows))
        self.assertIn('8 &rarr; 3 contigs (5 fewer)',html)
        self.assertIn('+200 bp',html)
        self.assertNotIn('8 &rarr; 5 contigs',html)
        self.assertNotIn('<th>recovered</th>',html)

    def test_raw_checkm_is_used_when_summary_scores_are_missing(self):
        write(self.root/'final/checkm.tsv','Bin Id\tCompleteness\tContamination\tStrain heterogeneity\n'+'genus__Test\t95.2\t3.1\t0\n')
        rows=report.collect_bins(self.layout,True)
        self.assertEqual(rows[0]['completeness_percent'],95.2)
        self.assertEqual(rows[0]['contamination_percent'],3.1)
        html=report.render_bins(rows,report.rank_totals(rows))
        self.assertIn('completeness %',html);self.assertIn('contamination %',html)
        self.assertIn('95.20',html);self.assertIn('3.10',html)

    def test_searchable_bins_retain_rank_quality_size_and_escape_names(self):
        row=dict(rank='species',bin='Sodalis_test_"<name>"',num_contigs=3,total_length_bp=3000,
                 completeness_percent=95.2,contamination_percent=3.1,unclassified=False)
        html=report.render_bins([row],report.rank_totals([row]))
        self.assertIn('id="bin-search"',html)
        self.assertIn('data-rank="species" data-bin-name="Sodalis_test_&quot;&lt;name&gt;&quot;"',html)
        self.assertIn('95.20',html);self.assertIn('3.10',html)
        self.assertIn('3.00 kb',html)
        self.assertNotIn('data-bin-name="Sodalis_test_"<name>"',html)

    def test_missing_contamination_does_not_pass_overview_quality_filter(self):
        rows=report.collect_bins(self.layout,True);rows[0]['completeness_percent']=95
        html=report.render_overview(self.layout,rows,report.rank_totals(rows),None,{})
        self.assertIn('0 / 1',html)

    def test_no_expansion_tab_or_pane(self):
        opts=argparse.Namespace(no_fasta_stats=False,min_length=1000,max_points=25000,include_uncoordinated=False,title=None)
        html=report.build_report(self.layout,opts)
        self.assertNotIn('pane-expansion',html)
        self.assertNotIn('"label": "Expansion"',html)
        self.assertIn('refinement-cell',html)

    def test_only_requested_tabs_and_first_available_rank_when_phylum_absent(self):
        opts=argparse.Namespace(no_fasta_stats=False,min_length=1000,max_points=25000,include_uncoordinated=False,title=None)
        html=report.build_report(self.layout,opts)
        self.assertIn('id="pane-bins" class="pane on"',html)
        self.assertIn('<option value="genus" selected>Genus</option>',html)
        for name in ('overview','inspection','files'):
            self.assertNotIn('id="pane-'+name+'"',html)
        self.assertNotIn('>all ranks</option>',html)
        self.assertIn('Classified bins',html)
        self.assertIn('Assembly size',html)
        self.assertIn('Contigs plotted',html)


if __name__=='__main__':unittest.main()

