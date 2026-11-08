"""Small, explicitly synthetic evidence for pathway integration regressions."""
import csv
import hashlib
import json
from pathlib import Path


def write_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer=csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)


def fixture(root, output=None):
    root=Path(root);output=Path(output) if output else root/'final/pathways/genus'
    output.mkdir(parents=True,exist_ok=True)
    inputs=[];qc=[];summary=[];components=[];hits=[]
    config={'pathways':[dict(id='p1',label='Synthetic panel',scored_components=['c1','c2'],route_components=['c1','c2','c3'],extra_components=[],target_genes=['g1','g2','AAT'],notes='Synthetic fixture only; unmeasured c3.')],
            'components':{'c1':{'assay_scope':'targeted'},'c2':{'assay_scope':'targeted'},'c3':{'assay_scope':'not_assayed'}}}
    config_bytes=(json.dumps(config,indent=2)+'\n').encode();(output/'config_used.json').write_bytes(config_bytes)
    definitions=[dict(pathway='p1',component=c,label=c,assay_scope='not_assayed' if c=='c3' else 'targeted',in_heatmap_denominator=int(c!='c3'),in_route_definition=1,notes='Synthetic') for c in ['c1','c2','c3']]
    for name in ['A','Zero','Failed','Stale']:
        source=root/'final/bins/genus'/(name+'.fasta');source.parent.mkdir(parents=True,exist_ok=True);source.write_text('>'+name+'_contig\nATGAAATAA\n')
        digest=hashlib.sha256(source.read_bytes()).hexdigest()
        inputs.append(dict(bin_id=name,path=str(source),sha256='0'*64 if name=='Stale' else digest))
        qc.append(dict(bin_id=name,status='failed' if name=='Failed' else 'ok',source=str(source),contigs=1,bases=9,orfs=1,error='Synthetic failure' if name=='Failed' else '',folder=''))
        summary.append(dict(bin_id=name,pathway='p1',panel_percent='NA' if name=='Failed' else '0' if name=='Zero' else '100',panel_status='bin_failed' if name=='Failed' else 'not_detected_in_panel' if name=='Zero' else 'complete_target_panel',route_status='bin_failed' if name=='Failed' else 'not_detected_in_panel' if name=='Zero' else 'unresolved_unmeasured_steps',detected_scored_components='NA' if name=='Failed' else 0 if name=='Zero' else 2,total_scored_components='NA' if name=='Failed' else 2,detected_components='c1;c2' if name in ['A','Stale'] else '',missing_scored_components='c1;c2' if name=='Zero' else '',candidate_components='c1' if name=='Zero' else '',unmeasured_route_steps='c3',operon_support='NA',pathway_notes='Synthetic panel; not a demonstrated route',error=''))
        if name=='Failed':continue
        for c in ['c1','c2','c3']:
            components.append(dict(bin_id=name,component=c,state='not_assayed' if c=='c3' else 'not_detected' if name=='Zero' else 'detected',value='NA' if c=='c3' else 0 if name=='Zero' else 1,copies=0 if c=='c3' or name=='Zero' else 1,proteins='c0000001_1' if name!='Zero' and c!='c3' else '',candidate_proteins='c0000001_1' if name=='Zero' and c=='c1' else '',assay_scope='not_assayed' if c=='c3' else 'targeted',notes='Synthetic fixture'))
        folder=output/'bins'/(name+'_test')
        folder.mkdir(parents=True,exist_ok=True)
        (folder/'bin_status.json').write_text(json.dumps(dict(bin_id=name,status='ok')))
        gene=dict(protein='c0000001_1',contig='c0000001',original_contig=name+'_contig',start=1,end=9,strand=1,partial='00',length=2,rank=0)
        write_csv(folder/'orfs.csv',[gene],list(gene))
        (folder/'proteins.faa').write_text('>c0000001_1\nMK*\n');(folder/'genes.fna').write_text('>c0000001_1\nATGAAATAA\n')
        for component in (['c1','c2'] if name!='Zero' else ['']):
            hits.append(dict(bin_id=name,protein='c0000001_1',contig=name+'_contig',start=1,end=9,strand=1,orf_rank=0,partial='00',protein_length=2,model_accession='SYNTHETIC_'+(component or 'AAT'),source_database='Synthetic',source_product='Synthetic product',target_gene='AAT' if name=='Zero' else 'g1' if component=='c1' else 'g2',target_call='shared_family' if component=='c2' else 'specific',source_group='context' if name=='Zero' else 'function_models',status='context_candidate' if name=='Zero' else 'accepted',reason='Synthetic fixture only',components=component,sequence_bit_score=100,domain_bit_score=95,sequence_evalue='1e-20',independent_evalue='2e-20',hmm_coverage=.8,hmm_intervals='[[1,2]]',protein_intervals='[[1,2]]',normalized_cutoff_margin=.2))
    manifest=dict(program='PathwayGenie',version='synthetic-test',finished_utc='2026-01-01T00:00:00+00:00',bin_inputs=inputs,config_sha256=hashlib.sha256(config_bytes).hexdigest(),arguments={},versions={})
    (output/'run_manifest.json').write_text(json.dumps(manifest))
    write_csv(output/'bin_qc.csv',qc,list(qc[0]))
    write_csv(output/'pathway_summary.long.csv',summary,list(summary[0]))
    write_csv(output/'component_definitions.csv',definitions,list(definitions[0]))
    write_csv(output/'component_evidence.long.csv',components,list(components[0]))
    write_csv(output/'all_model_hits.csv',hits,list(hits[0]))
    write_csv(output/'pathway_heatmap.csv',[dict(bin_id='A',p1=100)],['bin_id','p1'])
    write_csv(output/'neighborhoods.csv',[],['bin_id','pathway','neighborhood_id','contig','start','end'])
    write_csv(output/'operon_candidates.csv',[],['bin_id','pathway','neighborhood_id','template_id','score'])
    return output
