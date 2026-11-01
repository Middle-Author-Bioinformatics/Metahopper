"""Synthetic feature records for the browser regression; never used in reports."""
from metahopper_annotation_view import FIELDS, render_browser

rows = [['MH_1','dnaA','Chromosomal replication initiator','cds','contig_1','1','100','+','COG:0593'],
        ['MH_2','','<img src=x onerror=alert(1)>','cds','contig_2','1','100','-','']]
rows += [[f'MH_{i}','','hypothetical protein','cds','contig_1','1','100','+',''] for i in range(3,206)]
rows += [['MH_RNA','','tRNA-Leu','trna','contig_1','10','80','+','']]
first = dict(key='genus/Test',rank='genus',bin='Test',features=rows,error='',
             features_path='final/annotation/genus/Test/features.tsv',
             protein_path='final/annotation/genus/Test/bakta/annotation.faa',fasta_path='final/bins/genus/Test.fasta')
second = dict(first,key='genus/Second',bin='Second',features=[rows[0]],
              protein_path='final/annotation/genus/Second/bakta/annotation.faa')
if __name__ == '__main__':
    print(render_browser(dict(schema='metahopper-annotation-features-1',run_name='Test',stage='final',fields=FIELDS,bins=[first,second])))
