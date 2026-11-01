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


if __name__=='__main__':
    unittest.main()
