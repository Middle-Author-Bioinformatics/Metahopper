import csv
import hashlib
import json
from pathlib import Path
import re
import tempfile
import unittest

import metahopper_annotation_view as view
import metahopper_report as report


class AnnotationViewTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.directory = self.root / 'final/annotation/genus/Test'
        self.directory.mkdir(parents=True)
        self.bin = self.root / 'final/bins/genus/Test.fasta'
        self.bin.parent.mkdir(parents=True)
        self.bin.write_text('>contig_1\nATGAAATAA\n')
        self.row = dict(rank='genus', bin='Test', status='completed', cds='2', trna='1',
                        output_directory='final/annotation/genus/Test/bakta')
        (self.directory / 'result.json').write_text(json.dumps(dict(status='completed',
            signature={'input_sha256':hashlib.sha256(self.bin.read_bytes()).hexdigest()})))
        self.fields = ['rank', 'bin'] + view.FIELDS
        self.features = [dict(rank='genus',bin='Test',locus_tag='MH_1',gene='dnaA',
            product='Chromosomal replication initiator',type='cds',contig='contig_1',start='1',end='9',strand='+',db_xrefs='COG:0593'),
            dict(rank='genus',bin='Test',locus_tag='MH_2',gene='',product='hypothetical protein',type='cds',contig='contig_1',start='2',end='8',strand='-',db_xrefs=''),
            dict(rank='genus',bin='Test',locus_tag='MH_3',gene='',product='tRNA-Leu',type='trna',contig='contig_1',start='3',end='7',strand='+',db_xrefs='')]
        self.write_features()
        with (self.directory.parent / 'summary.tsv').open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(self.row),delimiter='\t');writer.writeheader();writer.writerow(self.row)
        self.layout=report.RunLayout(self.root,'final')

    def tearDown(self):
        self.tmp.cleanup()

    def write_features(self):
        with (self.directory/'features.tsv').open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=self.fields,delimiter='\t');writer.writeheader();writer.writerows(self.features)

    def test_features_retain_genes_products_coordinates_and_empty_gene_names(self):
        data=view.collect_features(self.layout,[self.row])
        item=data['bins'][0]
        self.assertEqual(item['key'],'genus/Test')
        self.assertEqual(item['features'][0],['MH_1','dnaA','Chromosomal replication initiator','cds','contig_1','1','9','+','COG:0593'])
        self.assertEqual(item['features'][1][1],'')
        self.assertEqual(item['features'][2][3],'trna')
        self.assertEqual(item['protein_path'],'final/annotation/genus/Test/bakta/annotation.faa')

    def test_rendered_summary_and_bins_link_to_the_same_gene_browser(self):
        content=report.render_annotation(self.layout)
        self.assertIn('data-annotation-bin="genus/Test"',content)
        self.assertIn('id="annotation-feature-table"',content)
        data=json.loads(re.search(r'id="annotation-feature-data">(.*?)</script>',content,re.S)[1])
        self.assertEqual(data['bins'][0]['features'][0][1],'dnaA')
        bins=report.render_bins([dict(rank='genus',bin='Test')],[dict(rank='genus',num_classified_bins=1,total_bp=9)],annotation_keys={('genus','Test')})
        self.assertIn('data-annotation-bin="genus/Test"',bins)

    def test_changed_bin_is_stale_and_its_old_genes_are_not_embedded(self):
        self.bin.write_text('>contig_1\nATGCCCTAA\n')
        rows=report.collect_annotation_rows(self.layout)
        self.assertEqual(rows[0]['status'],'stale')
        self.assertEqual(view.collect_features(self.layout,rows)['bins'],[])

    def test_feature_rows_from_another_bin_are_rejected(self):
        self.features[1]['bin']='Another'
        self.write_features()
        item=view.collect_features(self.layout,[self.row])['bins'][0]
        self.assertEqual(item['features'],[])
        self.assertIn('different bin',item['error'])

    def test_missing_feature_file_is_explained_without_fabricated_gene_rows(self):
        (self.directory/'features.tsv').unlink()
        item=view.collect_features(self.layout,[self.row])['bins'][0]
        self.assertTrue(item['error'])
        self.assertEqual(item['features'],[])

    def test_annotation_strings_cannot_end_the_embedded_json_script(self):
        self.features[0]['product']='</script><img src=x onerror=alert(1)>'
        self.write_features()
        data=view.collect_features(self.layout,[self.row])
        content=view.render_browser(data)
        self.assertNotIn('</script><img',content)
        parsed=json.loads(re.search(r'id="annotation-feature-data">(.*?)</script>',content,re.S)[1])
        self.assertEqual(parsed['bins'][0]['features'][0][2],self.features[0]['product'])

    def test_symlinked_features_outside_the_sample_are_rejected(self):
        feature=self.directory/'features.tsv';outside=self.root.parent/(self.root.name+'.outside.tsv')
        outside.write_bytes(feature.read_bytes());feature.unlink();feature.symlink_to(outside)
        try:
            item=view.collect_features(self.layout,[self.row])['bins'][0]
            self.assertIn('outside',item['error']);self.assertEqual(item['features'],[])
        finally:outside.unlink()

    def write_json(self, features):
        path=self.directory/'bakta/annotation.json'
        path.parent.mkdir(exist_ok=True)
        path.write_text(json.dumps({'features':features}))

    def test_bakta_json_preserves_sequences_and_overlapping_categories(self):
        self.write_json([
            dict(type='cds',sequence='contig_1',locus='B_1',gene='dnaA',product='Replication initiator',start=1,stop=9,strand='+',aa='MK*',nt='ATGAAATAA',db_xrefs=['COG:0593']),
            dict(type='cds',sequence='contig_1',locus='B_2',product='hypothetical protein',start=2,stop=8,strand='-',pseudogene={'description':'Frameshift'},aa='M',nt='ATG'),
            dict(type='cds',sequence='contig_1',locus='B_3',product='ATPase',pseudogene={}),
            dict(type='trna',sequence='contig_1',locus='B_4',product='tRNA-Leu')])
        item=view.collect_features(self.layout,[self.row])['bins'][0]
        self.assertEqual(item['source'],'Bakta JSON')
        self.assertEqual(item['features'][0],['B_1','dnaA','Replication initiator','cds','contig_1','1','9','+','COG:0593'])
        self.assertEqual(item['details'][0]['aa'],'MK*')
        self.assertEqual(item['details'][0]['nt'],'ATGAAATAA')
        self.assertEqual(item['details'][0]['aa_length'],2)
        self.assertIs(item['details'][0]['pseudogene'],False)
        self.assertIs(item['details'][1]['hypothetical'],True)
        self.assertIs(item['details'][1]['pseudogene'],True)
        self.assertEqual(item['details'][1]['pseudogene_note'],'Frameshift')
        self.assertIs(item['details'][2]['pseudogene'],True)
        self.assertIsNone(item['details'][3]['hypothetical'])

    def test_old_table_flags_are_unknown_and_sequences_match_by_locus(self):
        bakta=self.directory/'bakta';bakta.mkdir()
        (bakta/'annotation.faa').write_text('>MH_2 hypothetical\nMQQ*\n>OTHER unrelated\nABC\n')
        (bakta/'annotation.ffn').write_text('>MH_1 dnaA\nATGAAA\n')
        item=view.collect_features(self.layout,[self.row])['bins'][0]
        self.assertIsNone(item['details'][0]['pseudogene'])
        self.assertEqual(item['details'][0]['aa'],'')
        self.assertEqual(item['details'][0]['nt'],'ATGAAA')
        self.assertEqual(item['details'][1]['aa'],'MQQ*')
        self.assertEqual(item['details'][1]['aa_length'],3)
        self.assertIs(item['details'][1]['hypothetical'],True)

    def test_json_foreign_contig_is_rejected_without_tsv_substitution(self):
        self.write_json([dict(type='cds',sequence='other_sample',locus='B_1')])
        item=view.collect_features(self.layout,[self.row])['bins'][0]
        self.assertEqual(item['features'],[])
        self.assertIn('outside this bin',item['error'])

    def test_malformed_json_does_not_silently_become_tsv_data(self):
        self.write_json([])
        (self.directory/'bakta/annotation.json').write_text('{broken')
        item=view.collect_features(self.layout,[self.row])['bins'][0]
        self.assertEqual(item['features'],[])
        self.assertTrue(item['error'])

    def test_absent_products_and_unavailable_flags_do_not_become_hypothetical_or_pseudo(self):
        detail=view.feature_detail(dict(pseudogene='unavailable'),'', 'cds')
        self.assertIsNone(detail['hypothetical'])
        self.assertIsNone(detail['pseudogene'])
        rows,details=view.records_from_json({'features':[dict(type='cds')]})
        self.assertIsNone(details[0]['hypothetical'])

    def test_table_precedes_figures_and_download_labels_are_unbroken(self):
        bins=report.render_bins([dict(rank='phylum',bin='Test')],[dict(rank='phylum',num_classified_bins=1,total_bp=9)])
        self.assertLess(bins.index('id="bin-table"'),bins.index('id="bin-figures"'))
        annotation=report.render_annotation(self.layout)
        self.assertLess(annotation.index('id="annotation-summary-table"'),annotation.index('id="annotation-overview-figure"'))
        self.assertIn('white-space:nowrap!important',annotation)
        self.assertIn('.asset-link[hidden]{display:none}',annotation)


if __name__=='__main__':
    unittest.main()
