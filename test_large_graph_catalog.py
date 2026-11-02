"""Large source graphs retain complete display neighbourhoods and sequence proof."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import metahopper_graph as graph


def node(name):
    return dict(id=str(name),length=9,depth=None,sequence_present=False)


def edge(a,b,source='+',target='+',kind='L'):
    return dict(source=str(a),source_orient=source,target=str(b),target_orient=target,overlap='0M',kind=kind)


def catalogue(seed_ids=(0,),size=9):
    return dict(graphs=[dict(id='source.gfa',nodes=[node(i) for i in range(size)],
                            links=[edge(i,i+1) for i in range(size-1)])],
                contigs={f'contig_{i}':[dict(graph='source.gfa',steps=[[str(i),'+']])] for i in seed_ids},warnings=[])


class CompactGraphTests(unittest.TestCase):
    def test_complete_two_link_neighbourhood_excludes_only_more_distant_nodes(self):
        data=catalogue()
        graph.compact_catalog(data,min_nodes=2)
        source=data['graphs'][0]
        self.assertEqual([n['id'] for n in source['nodes']],['0','1','2'])
        self.assertEqual(source['links'],[edge(0,1),edge(1,2)])
        self.assertEqual(source['catalogue_scope'],dict(kind='mapped_contig_neighbourhoods',hops=2,
                         source_nodes=9,source_links=8,retained_nodes=3,retained_links=2))

    def test_union_of_all_mapped_contigs_preserves_both_neighbourhoods(self):
        data=catalogue(seed_ids=(0,8))
        before=deepcopy(data['contigs'])
        graph.compact_catalog(data,min_nodes=2)
        self.assertEqual({n['id'] for n in data['graphs'][0]['nodes']},{'0','1','2','6','7','8'})
        self.assertEqual(data['contigs'],before)
        self.assertTrue(all(e['source'] in {'0','1','2','6','7','8'} and e['target'] in {'0','1','2','6','7','8'} for e in data['graphs'][0]['links']))

    def test_orientation_jump_records_and_edges_between_retained_nodes_are_preserved(self):
        data=catalogue()
        data['graphs'][0]['links'].append(edge(2,0,'-','+','J'))
        graph.compact_catalog(data,min_nodes=2)
        self.assertIn(edge(2,0,'-','+','J'),data['graphs'][0]['links'])

    def test_equal_node_ids_in_different_sources_stay_separate(self):
        data=catalogue()
        other=deepcopy(data['graphs'][0]);other['id']='other.gfa';data['graphs'].append(other)
        data['contigs']['other']=[dict(graph='other.gfa',steps=[['8','-']])]
        graph.compact_catalog(data,min_nodes=2)
        self.assertEqual([n['id'] for n in data['graphs'][0]['nodes']],['0','1','2'])
        self.assertEqual([n['id'] for n in data['graphs'][1]['nodes']],['6','7','8'])

    def test_small_graphs_are_unchanged_and_unmapped_large_sources_are_explained(self):
        data=catalogue();before=deepcopy(data)
        graph.compact_catalog(data,min_nodes=20)
        self.assertEqual(data,before)
        data['contigs']={};graph.compact_catalog(data,min_nodes=2)
        self.assertEqual(data['graphs'],[])
        self.assertIn('full GFA remains available',data['warnings'][0])
        with self.assertRaises(ValueError):graph.compact_catalog(catalogue(),hops=1)

    def test_filtered_sequence_index_keeps_all_ambiguous_and_reverse_complement_matches(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'source.gfa';path.write_text('S\ta\tAAAC\nS\tb\tGTTT\nS\tc\tAGCG\n')
            _,nodes,_,index=graph.load_gfa(path,'source',sequence_queries={graph.sequence_key('AAAC')})
            self.assertEqual(len(nodes),3)
            self.assertEqual(list(index),[graph.sequence_key('AAAC')])
            self.assertEqual(index[graph.sequence_key('AAAC')],[['a','+'],['b','-']])
            self.assertEqual(graph.match_sequence('AAAC',index),(None,'ambiguous_exact_sequence'))

    def test_name_match_with_wrong_sequence_is_rejected_even_when_filtered_index_is_empty(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            for sub,text in {
                'final/assembly/consolidated_contigs.fasta':'>final_1\nAAAA\n',
                'final/assembly/contig_provenance.tsv':'final_contig\tsource_type\tsource_record\nfinal_1\toriginal\told\n',
                'megahit/final.contigs.fa':'>old\nAAAA\n',
                'megahit/assembly.gfa':'S\told\tCCCC\n',
            }.items():
                path=root/sub;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text)
            data=graph.build_catalog(root)
            self.assertNotIn('final_1',data['contigs'])
            self.assertIn('final_1',data['unmapped'])

    def test_source_above_old_report_limit_is_imported_then_compacted_without_changing_gfa(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);(root/'megahit').mkdir();(root/'final/assembly').mkdir(parents=True)
            (root/'final/assembly/consolidated_contigs.fasta').write_text('>final_1\nAAACCCGGG\n')
            (root/'final/assembly/contig_provenance.tsv').write_text('final_contig\tsource_type\tsource_record\nfinal_1\toriginal\t0\n')
            (root/'megahit/final.contigs.fa').write_text('>0\nAAACCCGGG\n')
            path=root/'megahit/assembly.gfa'
            with path.open('w') as handle:
                for i in range(1000001):handle.write(f'S\t{i}\t*\tLN:i:9\n')
                for i in range(3):handle.write(f'L\t{i}\t+\t{i+1}\t+\t0M\n')
            before=path.stat()
            data=graph.build_catalog(root)
            self.assertEqual(data['counts']['matched_contigs'],1)
            self.assertEqual(len(data['graphs'][0]['nodes']),3)
            self.assertEqual(data['graphs'][0]['catalogue_scope']['source_nodes'],1000001)
            self.assertEqual((before.st_size,before.st_mtime_ns),(path.stat().st_size,path.stat().st_mtime_ns))
            self.assertLess(len(json.dumps(data)),5000)


if __name__=='__main__':unittest.main()
