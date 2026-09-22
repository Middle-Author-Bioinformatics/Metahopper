#!/usr/bin/env python3
"""MetaHopper coverage, terminal-junction and endosymbiont/reference inspection.

Standalone: python metahopper_inspect.py -i RUN --references REFERENCES -t 24
One genome per FASTA (multi-contig allowed); group references in genus subfolders.
See MetaHopper_inspection_README.md for evidence limits and reference metadata.
"""
import argparse
import csv
import gzip
import hashlib
import html
import json
import math
import re
import shutil
import statistics
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path

# Screening labels, not a claim that every host should contain every organism.
PANEL = {
    'Sulcia': (('sulcia', 'karelsulcia'), 11, 'nutritional symbiont'),
    'Vidania': (('vidania',), 4, 'nutritional symbiont'),
    'Nasuia': (('nasuia',), 4, 'nutritional symbiont'),
    'Zinderia': (('zinderia',), 4, 'nutritional symbiont; also relevant to spittlebugs'),
    'Baumannia': (('baumannia',), 11, 'nutritional symbiont'),
    'Arsenophonus': (('arsenophonus',), 11, 'association varies by lineage'),
    'Sodalis': (('sodalis',), 11, 'association varies by lineage'),
    'Wolbachia': (('wolbachia',), 11, 'facultative/intracellular associate'),
    'Rickettsia': (('rickettsia',), 11, 'facultative/intracellular associate'),
    'Cardinium': (('cardinium',), 11, 'facultative/intracellular associate'),
    'Spiroplasma': (('spiroplasma',), 4, 'facultative associate'),
}
FASTA_SUFFIXES = ('.fa', '.fna', '.fasta', '.fas', '.fa.gz', '.fna.gz', '.fasta.gz')


def fasta(path):
    records = {}
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt') as fh:
        name, parts = None, []
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    records[name] = ''.join(parts).upper()
                name = line[1:].split()[0]
                if name in records:
                    raise ValueError(f'Duplicate FASTA identifier {name} in {path}')
                parts = []
            else:
                parts.append(line.strip())
        if name is not None:
            if name in records:
                raise ValueError(f'Duplicate FASTA identifier {name} in {path}')
            records[name] = ''.join(parts).upper()
    return records


def write_fasta(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w') as fh:
        for name, seq in records.items():
            fh.write(f'>{name}\n')
            for start in range(0, len(seq), 80):
                fh.write(seq[start:start + 80] + '\n')


def tsv(path):
    if not path.is_file():
        return []
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write_tsv(path, rows, fields=None):
    fields = fields or list(dict.fromkeys(k for row in rows for k in row))
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def run(cmd, log, output=None):
    with open(log, 'a') as lf:
        lf.write('COMMAND: ' + ' '.join(map(str, cmd)) + '\n'); lf.flush()
        if output:
            with open(output, 'w') as out:
                p = subprocess.run(list(map(str, cmd)), stdout=out, stderr=lf)
        else:
            p = subprocess.run(list(map(str, cmd)), stdout=lf, stderr=lf)
    if p.returncode:
        raise RuntimeError(f'{cmd[0]} exited {p.returncode}; see {log}')


def stream(cmd, log):
    with open(log, 'a') as lf:
        p = subprocess.Popen(list(map(str, cmd)), stdout=subprocess.PIPE, stderr=lf, text=True)
        try:
            yield from p.stdout
        finally:
            p.stdout.close()
            if p.poll() is None:
                p.wait()
        if p.returncode:
            raise RuntimeError(f'{cmd[0]} exited {p.returncode}; see {log}')


def group_name(text):
    words = set(re.findall(r'[a-z]+', text.lower()))
    matched = [name for name, (aliases, _, _) in PANEL.items() if words.intersection(aliases)]
    return matched[0] if len(matched) == 1 else None


def code_for(group):
    return PANEL.get(group, ((), 11, 'user-defined group'))[1]


def safe(text):
    return re.sub(r'[^A-Za-z0-9_.-]', '_', text)


def map_reads(reference, r1, r2, work, threads, local=False):
    work.mkdir(parents=True, exist_ok=True)
    index = work / 'index'
    log = work / 'mapping.log'
    run(['bowtie2-build', '--threads', threads, reference, index], log)
    bam = work / 'reads.sorted.bam'
    with open(log, 'a') as lf:
        bt = subprocess.Popen(['bowtie2', '--very-sensitive-local' if local else '--very-sensitive',
             '--no-unal', '-X', '1000', '-p', str(threads), '-x', str(index), '-1', str(r1), '-2', str(r2)],
             stdout=subprocess.PIPE, stderr=lf)
        sort = subprocess.Popen(['samtools', 'sort', '-@', str(threads), '-o', str(bam), '-'],
                                stdin=bt.stdout, stderr=lf)
        bt.stdout.close()
        sr, br = sort.wait(), bt.wait()
    if sr or br:
        raise RuntimeError(f'Read mapping failed; see {log}')
    run(['samtools', 'index', bam], log)
    return bam


def coverage_from_lines(lines, lengths, window=1000):
    """Zero-fill positions omitted by samtools depth; bounded window arrays, no per-base arrays."""
    data = {cid: {'length': n, 'sum': 0, 'sum2': 0, 'covered': 0, 'ge10': 0,
                  'hist': Counter(), 'positions': 0, 'windows': [0] * ((n + window - 1)//window),
                  'window_covered': [0] * ((n + window - 1)//window)} for cid, n in lengths.items()}
    for line in lines:
        f = line.rstrip().split('\t')
        if len(f) < 3 or f[0] not in data:
            continue
        cid, pos, depth = f[0], int(f[1]) - 1, int(f[2]); d = data[cid]
        if not 0 <= pos < d['length']:
            raise ValueError('Depth coordinate outside reference')
        d['sum'] += depth; d['sum2'] += depth*depth; d['positions'] += 1
        d['covered'] += depth > 0; d['ge10'] += depth >= 10; d['hist'][depth] += 1
        d['windows'][pos//window] += depth
        d['window_covered'][pos//window] += depth > 0
    rows, windows = {}, []
    for cid, d in data.items():
        n = d['length']
        d['hist'][0] += n - d['positions']
        middle = ((n - 1)//2, n//2); c = 0; med = []
        for depth, count in sorted(d['hist'].items()):
            for m in middle:
                if c <= m < c + count:
                    med.append(depth)
            c += count
        mean = d['sum']/max(1,n)
        cv = math.sqrt(max(0, d['sum2']/max(1,n)-mean*mean))/mean if mean else None
        rows[cid] = dict(contig=cid, length_bp=n, mean_depth=mean,
                         median_depth=sum(med)/max(1,len(med)), breadth_1x=d['covered']/max(1,n),
                         breadth_10x=d['ge10']/max(1,n), depth_cv=cv)
        for i, total in enumerate(d['windows']):
            start, end = i*window, min(n,(i+1)*window)
            avg = total/max(1,end-start)
            windows.append(dict(contig=cid, start=start, end=end, mean_depth=avg,
                covered_fraction=d['window_covered'][i]/max(1,end-start),
                flag='zero_coverage' if total == 0 else ('low_relative_depth' if mean and avg < .2*mean else ('high_relative_depth' if mean and avg > 3*mean else ''))))
    return rows, windows


def terminal_overlap(seq, minimum=30, maximum=10000):
    """Longest exact suffix-prefix overlap within 20% of sequence length."""
    size = min(maximum, len(seq)//5)
    if size < minimum:
        return 0
    text = seq[:size] + '$' + seq[-size:]
    pi = [0]*len(text)
    for i in range(1,len(text)):
        j = pi[i-1]
        while j and text[i] != text[j]:
            j = pi[j-1]
        if text[i] == text[j]:
            j += 1
        pi[i] = j
    length = pi[-1]
    return length if length >= minimum and set(seq[:length]) <= set('ACGT') else 0


def junction_records(seqs):
    records, meta = {}, {}
    for cid, seq in seqs.items():
        if len(seq) < 1000:
            continue
        overlap = terminal_overlap(seq)
        core = seq[:-overlap] if overlap else seq
        flank = min(500, len(core)//3)
        jid = f'MHJ_{len(records):07d}'
        records[jid] = core[-flank:] + core[:flank]
        meta[jid] = dict(contig=cid, terminal_overlap_bp=overlap, junction_coordinate=flank,
                         tested_core_length_bp=len(core))
    return records, meta


def spans_junction(position, cigar, boundary, anchor=25):
    """Require a continuous aligned block across the join: gaps cannot fake support."""
    cursor = position
    for count, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        n = int(count)
        if op in 'M=X':
            if cursor <= boundary-anchor and cursor+n >= boundary+anchor:
                return True
            cursor += n
        elif op in 'DN':
            cursor += n
    return False


def junction_anchor_quality(position, cigar, boundary, qualities, anchor=25):
    if qualities == '*':
        return False
    query = 0
    for count, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        n = int(count)
        if op in 'M=X':
            if position <= boundary-anchor and position+n >= boundary+anchor:
                start = query + boundary - position - anchor
                segment = qualities[start:start+2*anchor]
                return len(segment) == 2*anchor and all(ord(c)-33 >= 20 for c in segment)
            query += n; position += n
        elif op in 'IS': query += n
        elif op in 'DN': position += n
    return False


def junction_support(lines, meta):
    reads, starts, pairs = defaultdict(set), defaultdict(set), defaultdict(set)
    for line in lines:
        f = line.rstrip().split('\t')
        if len(f) < 11 or f[2] not in meta:
            continue
        flag, mq, pos = int(f[1]), int(f[4]), int(f[3])-1
        if flag & (4|256|512|1024|2048) or mq < 20:
            continue
        aligned = sum(int(n) for n, op in re.findall(r'(\d+)([MI=X])', f[5]))
        nm = next((int(t[5:]) for t in f[11:] if t.startswith('NM:i:')), None)
        if not aligned or nm is None or nm/aligned > .02:
            continue
        jid = f[2]; boundary = meta[jid]['junction_coordinate']
        if spans_junction(pos, f[5], boundary) and junction_anchor_quality(pos, f[5], boundary, f[10]):
            reads[jid].add(f[0]); starts[jid].add((pos, bool(flag&16)))
        # Proper pairs bracket boundary without either mate necessarily crossing it.
        if flag&2 and not flag&8 and f[6] == '=' and 0 < int(f[8]) <= 1000:
            if pos < boundary-25 and int(f[7])-1 >= boundary+25:
                pairs[jid].add(f[0])
    return [dict(**info, spanning_templates=len(reads[jid]), distinct_alignment_starts=len(starts[jid]),
        bracketing_pairs=len(pairs[jid]), status=('junction_read_supported' if len(reads[jid]) >= 5 and len(starts[jid]) >= 3
        else 'limited_junction_support' if reads[jid] or pairs[jid] else 'no_unique_junction_support')) for jid, info in meta.items()]


def reference_inputs(folder):
    rows, issues = [], []
    manifest = {r['file']: r for r in tsv(folder/'references.tsv')}
    for path in sorted(folder.rglob('*')):
        if not path.is_file() or not path.name.lower().endswith(FASTA_SUFFIXES):
            continue
        rel = str(path.relative_to(folder)); entry = manifest.get(rel, {})
        op = gzip.open if path.name.endswith('.gz') else open
        with op(path,'rt') as fh:
            header = next((x.strip() for x in fh if x.startswith('>')), '')
        group = entry.get('group') or group_name(rel + ' ' + header)
        if not group and path.parent != folder:
            group = path.relative_to(folder).parts[0]
        if not group:
            issues.append(dict(file=rel, reason='Unassigned group: use a genus subfolder or references.tsv'))
            continue
        group = group_name(group) or group
        code = int(entry.get('translation_table') or code_for(group))
        if code not in (4,11):
            raise ValueError(f'{rel}: this bacterial inspection supports translation tables 4 or 11')
        seqs = fasta(path)
        if not seqs or any(set(seq)-set('ACGTRYSWKMBDHVN') for seq in seqs.values()):
            raise ValueError(f'{path}: expected nucleotide FASTA')
        rows.append(dict(path=path, group=group, label=entry.get('label') or path.name,
                         translation_table=code, kind='reference', length_bp=sum(map(len,seqs.values()))))
    return rows, issues


def reciprocal_markers(lines, owners, anchor, taxa):
    """Strict reciprocal best hits; reject near-tied alternatives (within 10%)."""
    hits = defaultdict(dict)
    for line in lines:
        q,s,identity,qcov,scov,score = line.rstrip().split('\t')
        if q not in owners or s not in owners or owners[q] == owners[s]:
            continue
        if float(identity)<30 or min(float(qcov),float(scov))<60:
            continue
        key=(q,owners[s]);hits[key][s]=max(float(score),hits[key].get(s,0))
    unique={}
    for key,values in hits.items():
        ordered=sorted(values.items(),key=lambda x:x[1],reverse=True)
        if len(ordered)==1 or ordered[1][1] < .9*ordered[0][1]:
            unique[key]=ordered[0][0]
    markers={}
    for gene,owner in owners.items():
        if owner != anchor:
            continue
        members={anchor:gene}
        for taxon in taxa:
            other=unique.get((gene,taxon))
            if other and unique.get((other,anchor))==gene:
                members[taxon]=other
        if len(members)>=max(3,math.ceil(.8*len(taxa))):
            markers[gene]=members
    return markers


def trim_alignment(alignment, taxa):
    lengths={len(seq) for seq in alignment.values()}
    if len(lengths)!=1:
        raise ValueError('MAFFT returned unequal sequence lengths')
    n=next(iter(lengths),0)
    allowed=set('ACDEFGHIKLMNPQRSTVWY')
    keep=[i for i in range(n) if sum(i<len(alignment.get(t,'')) and alignment[t][i] in allowed for t in taxa)>=math.ceil(.8*len(taxa))]
    return {t:''.join(alignment.get(t,'-'*n)[i] for i in keep) for t in taxa}


def newick_svg(newick, labels):
    """Render IQ-TREE's safe-ID Newick; orientation is arbitrary, not an inferred root."""
    tokens=re.findall(r'[(),:;]|[^(),:;\s]+',newick);cursor=0
    def node():
        nonlocal cursor
        children=[];label='';length=0.0
        if tokens[cursor]=='(':
            cursor+=1
            while True:
                children.append(node())
                if tokens[cursor]==',':cursor+=1;continue
                if tokens[cursor]!=')':raise ValueError('Malformed Newick')
                cursor+=1;break
        if cursor<len(tokens) and tokens[cursor] not in ':,);':label=tokens[cursor];cursor+=1
        if cursor<len(tokens) and tokens[cursor]==':':
            cursor+=1;length=max(0,float(tokens[cursor]));cursor+=1
        return dict(children=children,label=label,length=length)
    tree=node();leaves=[];nodes=[]
    def positions(n,x):
        n['x']=x+n['length'];nodes.append(n)
        for c in n['children']:positions(c,n['x'])
        if n['children']:n['y']=sum(c['y'] for c in n['children'])/len(n['children'])
        else:n['y']=30+24*len(leaves);leaves.append(n)
    positions(tree,0);maximum=max((n['x'] for n in nodes),default=0) or 1
    width=max(1000,600+7*max((len(labels.get(n['label'],n['label'])) for n in leaves),default=0))
    height=max(100,len(leaves)*24+60)
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" style="width:100%;height:auto;min-width:700px" role="img" aria-label="Unrooted phylogenetic tree, arbitrary display orientation">']
    def draw(n):
        x=20+550*n['x']/maximum;y=n['y']
        if n['children']:
            ys=[c['y'] for c in n['children']]
            out.append(f'<path d="M{x},{min(ys)} V{max(ys)}" stroke="#334155" fill="none"/>')
            for c in n['children']:
                cx=20+550*c['x']/maximum
                out.append(f'<path d="M{x},{c["y"]} H{cx}" stroke="#334155" fill="none"/>');draw(c)
            if n['label']:out.append(f'<text x="{x+2}" y="{y-3}" font-size="9">{html.escape(n["label"])}</text>')
        else:
            out.append(f'<text x="{x+5}" y="{y+4}" font-size="12">{html.escape(labels.get(n["label"],n["label"]))}</text>')
    draw(tree)
    out.append(f'<text x="20" y="{height-12}" font-size="11">Branch scale: full horizontal span = {maximum:.4g} substitutions/site; internal labels = bootstrap support</text></svg>')
    return ''.join(out)


def build_tree(group, genomes, work, threads):
    result=dict(group=group, genomes=len(genomes), status='not_run')
    if not any(g['kind']=='bin' for g in genomes) or not any(g['kind']=='reference' for g in genomes):
        return dict(result,status='requires_bins_and_references')
    if len(genomes)<4:
        return dict(result,status='requires_at_least_four_genomes')
    work.mkdir(parents=True,exist_ok=True);log=work/'phylogeny.log'
    proteins,owners,labels,counts={}, {}, {}, {}
    genome_rows=[]
    for i,g in enumerate(genomes):
        gid=f'G{i:05d}';labels[gid]=f'{g["kind"]}: {g["label"]}'
        dna=work/f'{gid}.fna';write_fasta(dna,fasta(g['path']))
        faa=work/f'{gid}.faa'
        # Single-genome training plus explicit code; this avoids meta-mode code auto-selection.
        run(['prodigal','-i',dna,'-a',faa,'-o',work/f'{gid}.gff','-f','gff',
             '-p','single','-g',g['translation_table'],'-q'],log)
        predicted=fasta(faa);counts[gid]=len(predicted)
        for j,seq in enumerate(predicted.values()):
            pid=f'{gid}_P{j:06d}';proteins[pid]=seq.rstrip('*');owners[pid]=gid
        genome_rows.append(dict(id=gid,label=g['label'],kind=g['kind'],source=str(g['path']),
                                translation_table=g['translation_table'],proteins=len(predicted)))
    write_tsv(work/'genomes.tsv',genome_rows)
    refs=[row['id'] for row in genome_rows if row['kind']=='reference']
    anchor=max(refs,key=lambda gid:counts[gid]);allfaa=work/'proteins.faa';write_fasta(allfaa,proteins)
    db=work/'proteins';hits=work/'homology.tsv'
    run(['diamond','makedb','--in',allfaa,'-d',db,'--threads',threads],log)
    run(['diamond','blastp','-q',allfaa,'-d',db,'-o',hits,'--very-sensitive','--evalue','1e-5',
         '--max-target-seqs','0','--max-hsps','1','--threads',threads,'--outfmt','6',
         'qseqid','sseqid','pident','qcovhsp','scovhsp','bitscore'],log)
    with open(hits) as fh:markers=reciprocal_markers(fh,owners,anchor,list(labels))
    result['markers']=len(markers)
    if len(markers)<10:return dict(result,status='insufficient_single_copy_markers')
    occupancy={gid:sum(gid in members for members in markers.values())/len(markers) for gid in labels}
    write_tsv(work/'marker_occupancy.tsv',[dict(id=gid,label=labels[gid],marker_fraction=fraction) for gid,fraction in occupancy.items()])
    if min(occupancy.values())<.7:return dict(result,status='insufficient_genome_marker_occupancy')
    concat={gid:'' for gid in labels};marker_rows=[]
    for i,(anchor_gene,members) in enumerate(markers.items()):
        src=work/f'marker_{i:05d}.faa';aln=work/f'marker_{i:05d}.aln.faa'
        write_fasta(src,{gid:proteins[pid] for gid,pid in members.items()})
        run(['mafft','--auto','--thread',threads,src],log,aln)
        trimmed=trim_alignment(fasta(aln),list(labels));length=len(next(iter(trimmed.values())))
        if length<30:continue
        start=len(next(iter(concat.values())))+1
        for gid in concat:concat[gid]+=trimmed[gid]
        marker_rows.append(dict(marker=anchor_gene,start=start,end=start+length-1,genomes=len(members)))
    sites=len(next(iter(concat.values())));result.update(aligned_sites=sites,markers=len(marker_rows))
    write_tsv(work/'markers.tsv',marker_rows,['marker','start','end','genomes'])
    if sites<1000 or len(marker_rows)<10:return dict(result,status='insufficient_trimmed_alignment')
    # Require >=70% observed amino acids per genome after trimming too.
    if any(sum(c in 'ACDEFGHIKLMNPQRSTVWY' for c in seq)/sites<.7 for seq in concat.values()):
        return dict(result,status='insufficient_aligned_genome_occupancy')
    alignment=work/'concatenated.faa';write_fasta(alignment,concat)
    iq=shutil.which('iqtree2') or shutil.which('iqtree')
    prefix=work/'tree'
    run([iq,'-s',alignment,'-st','AA','-m','LG+F+G4','-bb','1000','-bnni','-nt',threads,
         '-seed','42','-pre',prefix,'-redo'],log)
    tree=Path(str(prefix)+'.treefile');svg=newick_svg(tree.read_text(),labels)
    (work/'tree.svg').write_text(svg)
    return dict(result,status='tree_built',newick=tree.read_text(),svg=svg,
                method='Within-group reference-anchored reciprocal-best-hit protein concatenation; LG+F+G4; 1000 ultrafast bootstraps; unrooted; exploratory')


def inspect_run(root, mode='endosymbionts', references=None, threads=8, rank=None,
                assembly=None, r1=None, r2=None, bam=None, assembly_sequences=None):
    started=time.monotonic();root=Path(root).resolve();out=root/'inspection';out.mkdir(parents=True,exist_ok=True)
    base=root/'final' if (root/'final/bins').is_dir() else root
    ranks=[p.name for p in (base/'bins').iterdir() if p.is_dir()] if (base/'bins').is_dir() else []
    rank=rank or next((r for r in ('genus','species','family') if r in ranks),next(iter(ranks),None))
    if not rank or not (base/'bins'/rank).is_dir():raise ValueError('No usable bin rank found')
    result=dict(version=1,rank=rank,mode=mode,bins=[],coverage=[],windows=[],junctions=[],trees=[],errors=[],
        notes=['Taxon labels are screening evidence, not proof of host association or absence.',
               'Depth uses MAPQ >=20, base quality >=20 and samtools overlap suppression; ambiguous repeats can look undercovered.',
               'Terminal junctions only: internal assembly joins are not independently validated.',
               'Junction support does not by itself establish a complete circular genome.',
               'Trees are unrooted exploratory within-group protein trees. Genus mixtures, long branches and reduced genomes require review.'])
    refrows,refissues=reference_inputs(Path(references).resolve()) if references else ([],[])
    result['reference_issues']=refissues
    custom={row['group'] for row in refrows}
    binmeta={row['bin']:row for row in tsv(Path(references)/'bins.tsv')} if references else {}
    bins=[];used=set()
    for path in sorted((base/'bins'/rank).glob('*.fasta')):
        entry=binmeta.get(path.stem,{})
        group=entry.get('group') or group_name(path.stem) or next((g for g in custom if safe(g).lower()==path.stem.lower()),None)
        if path.stem.lower()=='unclassified' or (mode=='endosymbionts' and not group):continue
        group=group_name(group) or group if group else None
        code=int(entry.get('translation_table') or code_for(group))
        if code not in (4,11):raise ValueError('Bin translation table must be 4 or 11')
        seqs=fasta(path);used.update(seqs)
        bins.append(dict(path=path,label=path.stem,group=group or path.stem,
                         translation_table=code,kind='bin',source='final_bin' if base!=root else 'preliminary_bin',
                         length_bp=sum(map(len,seqs.values())),seqs=seqs))
    if assembly is None:
        assembly=next((p for p in [base/'assembly/consolidated_contigs.fasta',root/'megahit/final.contigs.fa'] if p.is_file()),None)
    manifest={}
    if (root/'run_manifest.json').is_file():manifest=json.loads((root/'run_manifest.json').read_text())
    if assembly is None and manifest.get('input') and Path(manifest['input']).is_file():assembly=Path(manifest['input'])
    # Use a whole-assembly competitor reference, not just selected symbiont sequences.
    if assembly is None:
        allseqs={}
        for path in sorted((base/'bins'/rank).glob('*.fasta')):
            for cid,seq in fasta(path).items():
                if cid in allseqs:raise ValueError(f'Duplicate original contig ID {cid} across bins')
                allseqs[cid]=seq
        assembly=out/'all_bins_reference.fasta';write_fasta(assembly,allseqs)
    else:
        assembly=Path(assembly);allseqs=assembly_sequences if assembly_sequences is not None else fasta(assembly)
    # Recover inspection candidates below the bin-size filter without promoting them
    # to accepted genomes or changing the pipeline's bins.
    extras=defaultdict(dict)
    for row in tsv(base/'classification/contig_classification.tsv'):
        cid=row.get('contig');group=group_name(row.get(rank,'')+' '+row.get('genus','')+' '+row.get('species',''))
        excluded=any(str(row.get(k,'')).lower() in ('1','true','yes') for k in ('excluded_animal_plant','excluded_gene_density'))
        if group and cid in allseqs and cid not in used and not excluded:
            extras[group][cid]=allseqs[cid]
    for group,seqs in extras.items():
        path=out/'taxonomic_candidates'/f'{safe(group)}.fasta';write_fasta(path,seqs)
        bins.append(dict(path=path,label=group+'_unbinned_candidate',group=group,
                         translation_table=code_for(group),kind='bin',source='unbinned_taxonomic_candidate',
                         length_bp=sum(map(len,seqs.values())),seqs=seqs))
    selected={cid:seq for b in bins for cid,seq in b['seqs'].items()}
    result['panel']=[dict(taxon=g,status='candidate_detected' if any(b['group']==g for b in bins) else 'not_detected_in_taxonomic_bins',
                          translation_table_hint=code,role=role) for g,(_,code,role) in PANEL.items()]
    r1=Path(r1 or manifest['r1']) if r1 or manifest.get('r1') else None
    r2=Path(r2 or manifest['r2']) if r2 or manifest.get('r2') else None
    reads_available=bool(r1 and r2 and r1.is_file() and r2.is_file())
    if bam is None:
        candidate=base/'classification/coverage/reads_to_contigs.sorted.bam'
        bam=candidate if candidate.is_file() else None
    coverage={}
    if selected and (bam or reads_available):
        try:
            needed=['samtools']+([] if bam else ['bowtie2','bowtie2-build'])
            missing=[x for x in needed if not shutil.which(x)]
            if missing:raise RuntimeError('Missing coverage tools: '+', '.join(missing))
            bam=Path(bam) if bam else map_reads(assembly,r1,r2,out/'coverage_mapping',threads,True)
            header={}
            for line in stream(['samtools','view','-H',bam],out/'coverage.log'):
                if line.startswith('@SQ\t'):
                    tags=dict(x.split(':',1) for x in line.rstrip().split('\t')[1:]);header[tags['SN']]=int(tags['LN'])
            if any(header.get(cid)!=len(seq) for cid,seq in selected.items()):
                raise ValueError('BAM reference IDs/lengths do not match selected contigs; remap the current assembly')
            # BED limits depth output while retaining the whole-assembly mapping context.
            bed=out/'selected_contigs.bed'
            bed.write_text(''.join(f'{cid}\t0\t{len(seq)}\n' for cid,seq in selected.items()))
            lines=stream(['samtools','depth','-s','-q','20','-Q','20','-G','3844','-b',bed,bam],out/'coverage.log')
            coverage,windows=coverage_from_lines(lines,{cid:len(seq) for cid,seq in selected.items()})
            result['coverage']=list(coverage.values());result['windows']=windows
            write_tsv(out/'coverage.tsv',result['coverage'])
            write_tsv(out/'coverage_windows.tsv',windows)
        except (RuntimeError,ValueError,OSError) as exc:result['errors'].append(str(exc))
    elif selected:result['errors'].append('Coverage unavailable: no reads or matching BAM supplied')
    junctions,meta=junction_records(selected)
    if junctions and reads_available:
        try:
            missing=[x for x in ('bowtie2','bowtie2-build','samtools') if not shutil.which(x)]
            if missing:raise RuntimeError('Missing junction tools: '+', '.join(missing))
            # Distinct IDs for original contigs; junction mappings compete against all
            # original sequence, avoiding false uniqueness within a small target subset.
            competitor={f'MHO_{i:09d}':seq for i,seq in enumerate(allseqs.values())}
            competitor.update(junctions);ref=out/'junction_reference.fasta';write_fasta(ref,competitor)
            jb=map_reads(ref,r1,r2,out/'junction_mapping',threads,False)
            result['junctions']=junction_support(stream(['samtools','view',jb],out/'junction.log'),meta)
        except (RuntimeError,ValueError,OSError) as exc:
            result['errors'].append(str(exc))
            result['junctions']=[dict(**x,status='not_tested',reason=str(exc)) for x in meta.values()]
    else:
        result['junctions']=[dict(**x,status='not_tested',reason='Paired reads unavailable') for x in meta.values()]
    for cid,seq in selected.items():
        if len(seq)<1000:result['junctions'].append(dict(contig=cid,status='not_tested',reason='Contig shorter than 1000 bp'))
    write_tsv(out/'junctions.tsv',result['junctions'],['contig','terminal_overlap_bp','tested_core_length_bp',
              'junction_coordinate','spanning_templates','distinct_alignment_starts','bracketing_pairs','status','reason'])
    for b in bins:
        covs=[coverage[cid] for cid in b['seqs'] if cid in coverage];total=b['length_bp']
        refs=[g['length_bp'] for g in refrows if g['group']==b['group']]
        median=statistics.median(refs) if refs else None
        row={k:v for k,v in b.items() if k not in ('seqs','path','kind')}
        row.update(contigs=len(b['seqs']),gc_percent=100*sum(s.count('G')+s.count('C') for s in b['seqs'].values())/max(1,total),
            mean_depth=sum(c['mean_depth']*c['length_bp'] for c in covs)/max(1,total) if len(covs)==len(b['seqs']) else None,
            breadth_1x=sum(c['breadth_1x']*c['length_bp'] for c in covs)/max(1,total) if len(covs)==len(b['seqs']) else None,
            reference_median_bp=median,reference_size_ratio=total/median if median else None,
            junctions_supported=sum(j['status']=='junction_read_supported' and j['contig'] in b['seqs'] for j in result['junctions']))
        flags=[]
        if row['breadth_1x'] is not None and row['breadth_1x']<.95:flags.append('coverage_breadth_below_95pct')
        if median and not .5<=total/median<=1.5:flags.append('size_differs_from_supplied_references')
        if b['source']=='unbinned_taxonomic_candidate':flags.append('unbinned_candidate_not_validated_genome')
        row['review_flags']=';'.join(flags);result['bins'].append(row)
    write_tsv(out/'endosymbionts.tsv',result['bins'])
    if refrows:
        groups=defaultdict(list)
        for g in refrows+bins:groups[g['group']].append(g)
        dependencies=['prodigal','diamond','mafft']
        missing=[x for x in dependencies if not shutil.which(x)]
        if not(shutil.which('iqtree2') or shutil.which('iqtree')):missing.append('iqtree2')
        for group,genomes in sorted(groups.items()):
            try:
                if missing:raise RuntimeError('Missing tree tools: '+', '.join(missing))
                work=out/'phylogeny'/safe(group)
                # Exact-content cache: changing a genome/code or tool binary invalidates it.
                signature={'algorithm':1,'inputs':[
                    [str(g['path']),hashlib.sha256(g['path'].read_bytes()).hexdigest(),g['translation_table'],g['kind'],g['label']]
                    for g in genomes],'tools':[[x,shutil.which(x),Path(shutil.which(x)).stat().st_mtime_ns]
                        for x in dependencies+(['iqtree2'] if shutil.which('iqtree2') else ['iqtree'])]}
                cache=work/'result.json';sigpath=work/'inputs.json'
                if cache.is_file() and sigpath.is_file() and json.loads(sigpath.read_text())==signature:
                    tree=json.loads(cache.read_text());tree['cached']=True
                else:
                    if work.exists():shutil.rmtree(work)
                    tree=build_tree(group,genomes,work,threads)
                    if tree['status']=='tree_built':
                        cache.write_text(json.dumps(tree));sigpath.write_text(json.dumps(signature))
                result['trees'].append(tree)
            except (RuntimeError,ValueError,OSError) as exc:
                result['trees'].append(dict(group=group,status='failed',reason=str(exc)))
    result['elapsed_seconds']=time.monotonic()-started
    (out/'inspection.json').write_text(json.dumps(result,indent=2,allow_nan=False))
    write_tsv(out/'tree_summary.tsv',[{k:v for k,v in t.items() if k not in ('svg','newick')} for t in result['trees']],
              ['group','status','genomes','markers','aligned_sites','cached','reason','method'])
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i','--input',type=Path,required=True,help='Existing MetaHopper run directory')
    p.add_argument('--references',type=Path,help='Reference-genome folder; one genome per FASTA')
    p.add_argument('--inspect',choices=['endosymbionts','all'],default='endosymbionts')
    p.add_argument('-t','--threads',type=int,default=8)
    p.add_argument('--rank',help='Bin rank; default genus, then species')
    p.add_argument('--assembly',type=Path)
    p.add_argument('-1','--r1',type=Path);p.add_argument('-2','--r2',type=Path)
    p.add_argument('--bam',type=Path,help='Whole-assembly BAM for depth; junctions still require FASTQs')
    a=p.parse_args()
    if a.threads<1:p.error('Threads must be positive')
    if bool(a.r1)!=bool(a.r2):p.error('Supply both -1 and -2')
    if a.references and not a.references.is_dir():p.error('Reference folder does not exist')
    data=inspect_run(a.input,a.inspect,a.references,a.threads,a.rank,a.assembly,a.r1,a.r2,a.bam)
    # Regenerate the existing self-contained report without rerunning assembly.
    try:
        import metahopper_report as report
        opts=argparse.Namespace(no_fasta_stats=False,min_length=1000,max_points=25000,
                               include_uncoordinated=False,title=None)
        (a.input/'metahopper_report.html').write_text(report.build_report(report.RunLayout(a.input,'auto'),opts))
    except (ImportError,OSError,ValueError) as exc:print(f'Report not regenerated: {exc}')
    print(f'Inspection saved to {a.input / "inspection"}; {len(data["bins"])} candidate bins inspected')
    for error in data['errors']:print('Inspection note:',error)
    for tree in data['trees']:print('Tree:',tree['group'],tree['status'],tree.get('reason',''))


if __name__=='__main__':main()
