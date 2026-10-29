#!/usr/bin/env python3
"""Build a sequence-free GFA1 graph catalogue for MetaHopper reports.

Contig-to-segment links use recorded assembler identifiers, matching names in a
known source graph, or unique exact sequence matches (either strand). Numeric
contig suffixes and taxonomic co-membership are never used to invent links.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

SCHEMA = 'metahopper-graphs-1'
COMPLEMENT = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')
MEGAHIT_GRAPH_SCHEMA = 1


def opn(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def fasta(path):
    with opn(path) as f:
        name, parts = None, []
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(parts).upper()
                name, parts = line[1:].split()[0], []
            elif name is not None:
                parts.append(line.strip())
        if name is not None:
            yield name, ''.join(parts).upper()


def table(path):
    if not Path(path).is_file():
        return []
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def export_megahit_graph(root, intermediate_dir=None, force=False, max_links=2000000):
    """Save a GFA of same-k exact contig-end overlaps, keeping MEGAHIT IDs.

    This follows MEGAHIT's contig2fastg endpoint-matching algorithm, rather
    than claiming to recover its internal succinct de Bruijn graph. Fresh
    runs include final-round intermediate unitigs and earlier finalized
    contigs. Old runs can use their retained final contigs; absent unitigs
    and connections through them cannot be restored. Different k values
    are separate graphs within the GFA, with no invented cross-k links.
    """
    root = Path(root).resolve()
    source = root / 'megahit/final.contigs.fa'
    if not source.is_file():
        raise ValueError(f'No retained MEGAHIT contigs: {source}')
    if max_links < 1:
        raise ValueError('MEGAHIT graph link limit must be positive')
    dest = source.parent / 'assembly.gfa'
    manifest = source.parent / 'assembly_graph.json'

    def stamp(path):
        stat = Path(path).stat()
        return dict(size=stat.st_size, mtime_ns=stat.st_mtime_ns)

    source_stamp = stamp(source)
    try:
        previous = json.loads(manifest.read_text())
    except (OSError, ValueError):
        previous = {}
    if dest.exists() and previous.get('export_schema') != MEGAHIT_GRAPH_SCHEMA and not force:
        raise ValueError(f'{dest} already exists without an exporter record; preserve it or explicitly force export')
    if (not force and previous.get('export_schema') == MEGAHIT_GRAPH_SCHEMA
            and previous.get('status') == 'completed' and previous.get('input') == source_stamp
            and dest.is_file() and previous.get('output') == stamp(dest)):
        return dict(previous, cached=True)

    directory = Path(intermediate_dir) if intermediate_dir else root / 'work/megahit/intermediate_contigs'
    candidates = []
    for path in directory.glob('k*.contigs.fa'):
        match = re.fullmatch(r'k(\d+)\.contigs\.fa', path.name)
        if match and path.stat().st_size:
            candidates.append((int(match[1]), path))
    parts = [max(candidates, key=lambda x: x[0])[1]] if candidates else []
    parts += sorted(p for p in directory.glob('k*.final.contigs.fa') if p.stat().st_size)
    parts = list(dict.fromkeys(parts))
    scope = 'intermediates' if parts else 'retained_contigs'
    nodes = {}
    for path in parts + [source]:
        seen = set()
        for name, seq in fasta(path):
            match = re.fullmatch(r'k(\d+)_\d+', name)
            if not match or int(match[1]) < 1 or int(match[1]) % 2 != 1:
                raise ValueError(f'Cannot infer a valid MEGAHIT k value from contig ID {name!r} in {path}')
            if name in seen:
                raise ValueError(f'Duplicate contig ID {name!r} in {path}')
            seen.add(name)
            if not seq or not re.fullmatch(r'[ACGTRYSWKMBDHVN]+', seq):
                raise ValueError(f'Invalid DNA sequence for {name!r} in {path}')
            record = (seq, int(match[1]))
            if name in nodes and nodes[name] != record:
                raise ValueError(f'Conflicting MEGAHIT sequences for {name!r}')
            nodes[name] = record
    if not nodes:
        raise ValueError('No MEGAHIT contigs available for graph export')

    def reverse(seq):
        return seq.translate(COMPLEMENT)[::-1]

    starts = defaultdict(list)
    for name, (seq, k) in nodes.items():
        if len(seq) <= k:
            continue
        for end, orient in ((seq[:k], '+'), (reverse(seq[-k:]), '-')):
            if re.fullmatch(r'[ACGT]+', end):
                starts[(k, end)].append((name, orient))
    flip = {'+': '-', '-': '+'}
    links = 0
    # Never expose a partial graph or overwrite a valid graph on export failure.
    with tempfile.TemporaryDirectory(prefix='.megahit_graph_', dir=dest.parent) as temporary:
        partial = Path(temporary) / dest.name
        with partial.open('w') as handle:
            handle.write(f'H\tVN:Z:1.0\tTS:Z:megahit_kmer_overlaps\tSC:Z:{scope}\n')
            for name, (seq, k) in nodes.items():
                handle.write(f'S\t{name}\t{seq}\tLN:i:{len(seq)}\tKK:i:{k}\n')
            for name, (seq, k) in nodes.items():
                if len(seq) <= k:
                    continue
                for end, orient in ((seq[-k:], '+'), (reverse(seq[:k]), '-')):
                    for target, target_orient in starts.get((k, end), ()):
                        edge = (name, orient, target, target_orient)
                        opposite = (target, flip[target_orient], name, flip[orient])
                        if edge > opposite:
                            continue
                        links += 1
                        if links > max_links:
                            raise ValueError(f'MEGAHIT overlap graph exceeds {max_links:,} links; no partial graph saved')
                        handle.write(f'L\t{name}\t{orient}\t{target}\t{target_orient}\t{k}M\n')
        result = dict(export_schema=MEGAHIT_GRAPH_SCHEMA, status='completed', input=source_stamp,
                      output=stamp(partial), construction='megahit_kmer_overlaps', source_scope=scope,
                      source_files=[str(p.relative_to(root)) if p.is_relative_to(root) else str(p) for p in parts + [source]],
                      segments=len(nodes), links=links, kmer_sizes=sorted({k for _, k in nodes.values()}))
        record = Path(temporary) / manifest.name
        record.write_text(json.dumps(result, indent=2) + '\n')
        partial.replace(dest)
        record.replace(manifest)
    return dict(result, cached=False)


def sequence_key(seq):
    return (len(seq), hashlib.sha256(seq.encode('ascii')).hexdigest())


def path_steps(text):
    # GFA1.2 semicolons mark jumps; a path still identifies these segments.
    parts = re.split(r'(?<=[+-])[,;]', text)
    if not parts or any(len(x) < 2 or x[-1] not in '+-' for x in parts):
        raise ValueError('Malformed GFA path segment list')
    return [[x[:-1], x[-1]] for x in parts]


def load_gfa(path, graph_id, max_nodes=200000):
    nodes, links, paths, seq_index = {}, [], {}, defaultdict(list)
    metadata = {}
    unsupported = Counter()
    with opn(path) as f:
        for number, line in enumerate(f, 1):
            p = line.rstrip('\r\n').split('\t')
            kind = p[0]
            if not kind or kind.startswith('#'):
                continue
            try:
                if kind == 'H':
                    if any(x.startswith('VN:Z:2') for x in p[1:]):
                        raise ValueError('GFA2 is unsupported; provide GFA1')
                    tags = {x[:2]: x[5:] for x in p[1:] if len(x) >= 5 and x[2:5] == ':Z:'}
                    if tags.get('TS') == 'megahit_kmer_overlaps':
                        metadata.update(construction=tags['TS'], source_scope=tags.get('SC', 'unknown'))
                elif kind == 'S':
                    if len(p) < 3 or p[1] in nodes:
                        raise ValueError('Missing or duplicate segment identifier')
                    if p[2].isdigit() and len(p) >= 4:
                        raise ValueError('GFA2 segment record; provide GFA1')
                    if len(nodes) >= max_nodes:
                        raise ValueError(f'Graph exceeds {max_nodes:,} segments')
                    seq = p[2].upper()
                    tags = {x[:2]: x[5:] for x in p[3:] if len(x) >= 5 and x[2] == ':' and x[4] == ':'}
                    length = len(seq) if seq != '*' else int(tags.get('LN', 0))
                    if length < 0:
                        raise ValueError('Negative segment length')
                    depth = None
                    for tag in ('dp', 'DP'):
                        if tag in tags:
                            try:
                                value = float(tags[tag])
                                depth = value
                                if not (0 <= depth < float('inf')):
                                    depth = None
                                break
                            except ValueError:
                                pass
                    nodes[p[1]] = dict(id=p[1], length=length, depth=depth, sequence_present=seq != '*')
                    for tag, field in [('KC', 'kmer_count'), ('RC', 'read_count')]:
                        try:
                            if tag in tags:
                                count = int(tags[tag])
                                if count >= 0:
                                    nodes[p[1]][field] = count
                        except ValueError:
                            pass
                    if seq != '*':
                        seq_index[sequence_key(seq)].append([p[1], '+'])
                        seq_index[sequence_key(seq.translate(COMPLEMENT)[::-1])].append([p[1], '-'])
                elif kind in ('L', 'J'):
                    if len(p) < 6 or p[2] not in ('+', '-') or p[4] not in ('+', '-'):
                        raise ValueError('Malformed link/jump record')
                    links.append(dict(source=p[1], source_orient=p[2], target=p[3], target_orient=p[4],
                                      overlap=p[5], kind=kind))
                elif kind == 'P':
                    if len(p) < 4 or p[1] in paths:
                        raise ValueError('Malformed or duplicate path')
                    paths[p[1]] = path_steps(p[2])
                elif kind == 'W':
                    if len(p) < 7:
                        raise ValueError('Malformed walk')
                    steps = [[m[1], '+' if m[0] == '>' else '-'] for m in re.findall(r'([><])([^><]+)', p[6])]
                    if not steps or ''.join(('>' if s[1] == '+' else '<') + s[0] for s in steps) != p[6]:
                        raise ValueError('Malformed walk steps')
                    # SeqId alone is ambiguous across samples, haplotypes or intervals.
                    key = '#'.join(p[1:6])
                    if key in paths:
                        raise ValueError('Duplicate walk identifier')
                    paths[key] = steps
                else:
                    unsupported[kind] += 1
            except (ValueError, IndexError) as e:
                raise ValueError(f'{path.name}:{number}: {e}') from e
    if not nodes:
        raise ValueError('No GFA1 segments found')
    missing = sum(x['source'] not in nodes or x['target'] not in nodes for x in links)
    if missing:
        raise ValueError(f'{missing} links reference absent segments')
    for name, steps in paths.items():
        if any(s[0] not in nodes for s in steps):
            raise ValueError(f'Path {name} references absent segments')
    # Opposite-strand representations of the same physical link are duplicates.
    seen, unique = set(), []
    flip = {'+': '-', '-': '+'}
    for edge in links:
        a = (edge['source'], edge['source_orient'], edge['target'], edge['target_orient'])
        b = (edge['target'], flip[edge['target_orient']], edge['source'], flip[edge['source_orient']])
        key = (min(a, b), edge['kind'], edge['overlap'])
        if key not in seen:
            seen.add(key)
            unique.append(edge)
    graph = dict(id=graph_id, source=graph_id, nodes=list(nodes.values()), links=unique,
                 ignored_record_types=dict(unsupported), **metadata)
    return graph, nodes, paths, seq_index


def match_sequence(seq, seq_index):
    hits = seq_index.get(sequence_key(seq), [])
    distinct = {h[0] for h in hits}
    if len(distinct) != 1:
        return None, ('ambiguous_exact_sequence' if distinct else 'no_exact_segment_match')
    return [hits[0]], 'unique_exact_sequence'


def build_catalog(root, stage='final', contig_ids=None, graph_paths=(), max_nodes=200000):
    root = Path(root).resolve()
    if not root.is_dir() or stage not in ('final', 'preliminary'):
        raise ValueError('Graph input must be an existing run directory with a final or preliminary bin set')
    base = root / 'final' if stage == 'final' else root
    result = dict(schema=SCHEMA, stage=stage, run_name=root.name, graphs=[], contigs={}, unmapped={}, warnings=[])
    wanted = set(contig_ids) if contig_ids is not None else None
    prov = table(root / 'final/assembly/contig_provenance.tsv') if stage == 'final' else []
    by_source = defaultdict(list)
    if prov:
        for row in prov:
            cid = row.get('final_contig')
            if not cid or (wanted is not None and cid not in wanted):
                continue
            by_source[(row.get('source_type'), row.get('source_rank'), row.get('source_bin'))].append(row)
    else:
        source = next((p for p in [base / 'assembly/consolidated_contigs.fasta', root / 'assembly/preliminary_contigs.fasta',
                                   root / 'megahit/final.contigs.fa'] if p.is_file()), None)
        if source:
            for cid, _ in fasta(source):
                if wanted is None or cid in wanted:
                    by_source[('original', '', '')].append(dict(final_contig=cid, source_record=cid))
    if wanted is None:
        wanted = {row['final_contig'] for rows in by_source.values() for row in rows}
    result['lengths'] = {r['contig']: int(r['length_bp']) for r in table(base / 'classification/contig_classification.tsv')
                         if r.get('contig') in wanted and str(r.get('length_bp', '')).isdigit()}
    loaded = {}

    def get_graph(path):
        path = Path(path).resolve()
        if path in loaded:
            return loaded[path]
        try:
            graph_id = path.relative_to(root).as_posix()
        except ValueError:
            graph_id = 'external/' + hashlib.sha256(str(path).encode()).hexdigest()[:12] + '/' + path.name
        try:
            data = load_gfa(path, graph_id, max_nodes)
        except (ValueError, OSError) as e:
            result['warnings'].append(f'{graph_id}: {e}')
            loaded[path] = None
            return None
        loaded[path] = data
        result['graphs'].append(data[0])
        if data[0].get('construction') == 'megahit_kmer_overlaps':
            message = graph_id + ': MEGAHIT links represent exact contig-end overlaps within each k value.'
            if data[0].get('source_scope') == 'retained_contigs':
                message += ' Rebuilt from retained contigs only; missing intermediate segments and their links are not recoverable.'
            result['warnings'].append(message)
        if data[0]['ignored_record_types']:
            result['warnings'].append(graph_id + ': unsupported records omitted: ' +
                                      ', '.join(f'{k}={v}' for k, v in data[0]['ignored_record_types'].items()))
        return data

    preliminary_graphs = list(map(Path, graph_paths))
    for folder in (root / 'assembly', root / 'megahit'):
        if folder.is_dir():
            preliminary_graphs.extend(sorted(folder.glob('*.gfa')))
            preliminary_graphs.extend(sorted(folder.glob('*.gfa.gz')))
    preliminary_graphs = list(dict.fromkeys(p.resolve() for p in preliminary_graphs))

    original_ids = {r.get('source_record') for (kind, _, _), rows in by_source.items() if kind == 'original' for r in rows}
    original_sequences = {}
    if preliminary_graphs and original_ids:
        source = next((p for p in (root / 'assembly/preliminary_contigs.fasta', root / 'megahit/final.contigs.fa') if p.is_file()), None)
        if source:
            original_sequences = {c: s for c, s in fasta(source) if c in original_ids}
        elif (root / 'final/assembly/consolidated_contigs.fasta').is_file():
            final_to_original = {r['final_contig']: r.get('source_record') for (kind, _, _), rows in by_source.items()
                                 if kind == 'original' for r in rows}
            original_sequences = {final_to_original[c]: s for c, s in fasta(root / 'final/assembly/consolidated_contigs.fasta')
                                  if c in final_to_original}

    def matching(data, name, seq=None, alias=None):
        graph, nodes, paths, seq_index = data
        key = alias or name
        if key in paths:
            return paths[key], 'recorded_path' if alias else 'source_path_name'
        if key in nodes:
            # Known-source names may be trusted, but reject a detectable content mismatch.
            if seq is not None and seq_index:
                hits = seq_index.get(sequence_key(seq), [])
                if any(h[0] == key for h in hits):
                    return [next(h for h in hits if h[0] == key)], 'recorded_segment' if alias else 'source_segment_name'
                if nodes[key]['length'] and nodes[key]['length'] != len(seq):
                    return None, 'source_name_length_mismatch'
                # Equal length is not sufficient proof when sequence is available.
                if nodes[key]['sequence_present']:
                    return None, 'source_name_sequence_mismatch'
            return [[key, '+']], 'recorded_segment' if alias else 'source_segment_name'
        if seq is not None:
            return match_sequence(seq, seq_index)
        return None, 'no_recorded_graph_identifier'

    for (kind, rank, bin_name), rows in by_source.items():
        if kind == 'reassembled':
            # Never associate original fallbacks with graphs of rejected candidates.
            directory = root / 'reassembly' / str(rank) / str(bin_name)
            gpath = next((p for p in [directory / 'assembly.gfa', directory / 'assembly.gfa.gz'] if p.is_file()), None)
            data = get_graph(gpath) if gpath else None
            seqs = dict(fasta(directory / 'assembly.fasta')) if data and (directory / 'assembly.fasta').is_file() else {}
            aliases = {r['contig']: r['gfa_record'] for r in table(directory / 'contig_graph_aliases.tsv')
                       if r.get('contig') and r.get('gfa_record')}
            for row in rows:
                cid, original = row['final_contig'], row.get('source_record', '')
                if not data:
                    result['unmapped'][cid] = 'source_graph_unavailable'
                    continue
                # Old pipeline records were renamed; only a preserved alias or exact
                # sequence can recover numeric Unicycler identifiers safely.
                steps, method = matching(data, original, seqs.get(original), aliases.get(original))
                if steps:
                    result['contigs'][cid] = [dict(graph=data[0]['id'], steps=steps, method=method)]
                else:
                    result['unmapped'][cid] = method
        else:
            datas = [d for p in preliminary_graphs if (d := get_graph(p)) is not None]
            for row in rows:
                cid, original = row['final_contig'], row.get('source_record', '')
                matches = []
                for data in datas:
                    steps, method = matching(data, original, original_sequences.get(original))
                    if steps:
                        matches.append(dict(graph=data[0]['id'], steps=steps, method=method))
                if len(matches) == 1:
                    result['contigs'][cid] = matches
                else:
                    result['unmapped'][cid] = 'multiple_source_graph_matches' if matches else 'preliminary_graph_unavailable_or_unmatched'
    for cid in wanted:
        if cid not in result['contigs'] and cid not in result['unmapped']:
            result['unmapped'][cid] = 'provenance_unavailable'
    result['counts'] = dict(graphs=len(result['graphs']), matched_contigs=len(result['contigs']),
                            unmatched_contigs=len(result['unmapped']))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-i', '--input', type=Path, required=True)
    p.add_argument('-o', '--output', type=Path)
    p.add_argument('--bin-set', choices=['auto', 'preliminary', 'final'], default='auto')
    p.add_argument('--graph', type=Path, action='append', default=[], help='Additional GFA1 for preliminary contigs')
    p.add_argument('--max-nodes', type=int, default=200000)
    p.add_argument('--export-megahit-graph', action='store_true',
                   help='Export megahit/assembly.gfa before building the catalogue; old runs use retained contigs')
    p.add_argument('--force-megahit-graph', action='store_true', help='Rebuild the MEGAHIT GFA even when cached')
    p.add_argument('--megahit-max-links', type=int, default=2000000, help='Maximum exported MEGAHIT overlap links')
    a = p.parse_args()
    if a.max_nodes < 1:
        p.error('--max-nodes must be positive')
    if a.megahit_max_links < 1:
        p.error('--megahit-max-links must be positive')
    if a.force_megahit_graph and not a.export_megahit_graph:
        p.error('--force-megahit-graph requires --export-megahit-graph')
    if a.export_megahit_graph:
        try:
            exported = export_megahit_graph(a.input, force=a.force_megahit_graph, max_links=a.megahit_max_links)
        except (OSError, ValueError) as exc:
            p.exit(1, f'ERROR: {exc}\n')
        print(f'MEGAHIT GFA: {a.input / "megahit/assembly.gfa"}; '
              f'{exported["segments"]:,} segments, {exported["links"]:,} links; '
              f'source scope: {exported["source_scope"]}; cached: {exported["cached"]}')
    stage = 'final' if a.bin_set == 'auto' and (a.input / 'final/bins').is_dir() else a.bin_set
    if stage == 'auto':
        stage = 'preliminary'
    data = build_catalog(a.input, stage, graph_paths=a.graph, max_nodes=a.max_nodes)
    dest = a.output or a.input / f'assembly_graphs.{stage}.json'
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(data, separators=(',', ':'), allow_nan=False))
    print(f'Graph catalogue: {dest}; {data["counts"]}')
    for warning in data['warnings']:
        print('WARNING:', warning)


VIEWER_HTML = r'''
<dialog id="graph-dialog" class="graph-dialog" aria-labelledby="graph-heading">
  <div class="row" style="justify-content:space-between;align-items:center">
    <h2 id="graph-heading" style="margin:0">Assembly graph for selected contigs</h2>
    <button id="graph-close" type="button">Close</button>
  </div>
  <p class="note">Segments and links come from source GFA files. Contigs can map to one or more segments. Connections are never inferred from taxonomy or proximity in the scatter plot.</p>
  <div class="row">
    <div><label for="graph-hops">Neighbourhood</label><select id="graph-hops"><option value="0">Selected segments only</option><option value="1" selected>One link away</option><option value="2">Two links away</option></select></div>
    <div><label for="graph-limit">Node limit</label><select id="graph-limit"><option value="200">200</option><option value="500" selected>500</option><option value="1000">1000</option></select></div>
    <div><label for="graph-labels">Segment labels</label><select id="graph-labels"><option value="yes">Show labels</option><option value="no">Hide labels</option></select></div>
    <button id="graph-fit" type="button">Fit graph</button>
    <button id="graph-download" type="button" disabled>Download SVG</button>
  </div>
  <p id="graph-status" class="note" role="status"></p>
  <div id="graph-stage" class="graph-stage"></div>
  <p class="note">Drag the background to pan, scroll to zoom, or drag a node. Selected segments use the current scatter colours; grey nodes are neighbours. Hatched nodes belong to selected contigs with different colours. Dashed edges are GFA jump records. Node spacing is a layout choice, not genomic distance.</p>
  <p id="graph-detail" class="note graph-detail"></p>
  <details id="graph-unmapped"><summary>Contigs without a usable graph mapping</summary><p id="graph-missing" class="note graph-detail"></p></details>
  <details><summary>Load graph data for an existing report</summary>
    <p class="note">On the server, run <code>python metahopper_graph.py -i RUN --bin-set final</code>, then select the resulting assembly_graphs.final.json file below. Regenerating the report with the updated metahopper_report.py embeds these data automatically. A final report needs a final catalogue from the same run.</p>
    <input id="graph-file" type="file" accept=".json,application/json">
  </details>
</dialog>
'''

VIEWER_JS = r'''
  var graphData=DATA.assemblyGraphs||{schema:'metahopper-graphs-1',graphs:[],contigs:{},unmapped:{},warnings:[]};
  var graphSelection=[],graphSvg=null,graphGroup=null,graphNodes=[],graphEdges=[],graphTransform={x:0,y:0,k:1};
  var graphNodeMap=new Map(),graphAdj=new Map(),graphIndex=new Map(),graphDrag=null;
  var SVG_NS='http://www.w3.org/2000/svg';
  function svgElement(tag,attrs,text){
    var element=document.createElementNS(SVG_NS,tag);
    Object.keys(attrs||{}).forEach(function(key){element.setAttribute(key,attrs[key]);});
    if(text!==undefined)element.textContent=text;
    return element;
  }
  function graphKey(graph,node){return JSON.stringify([graph,node]);}
  function graphMessage(message){el('graph-status').textContent=message;}
  function graphIndexData(data){
    if(!data||data.schema!=='metahopper-graphs-1'||!Array.isArray(data.graphs)||!data.contigs||typeof data.contigs!=='object')throw Error('Not a MetaHopper graph catalogue');
    if(data.stage&&DATA.assemblyGraphs&&DATA.assemblyGraphs.stage&&data.stage!==DATA.assemblyGraphs.stage)throw Error('Catalogue bin set does not match this report');
    if(data.run_name&&DATA.runName&&data.run_name!==DATA.runName)throw Error('Catalogue run name does not match this report');
    var index=new Map();
    data.graphs.forEach(function(graph){
      if(typeof graph.id!=='string'||index.has(graph.id)||!Array.isArray(graph.nodes)||!Array.isArray(graph.links))throw Error('Invalid or duplicate graph source');
      var nodes=new Map();
      graph.nodes.forEach(function(node){
        if(typeof node.id!=='string'||nodes.has(node.id)||!Number.isFinite(node.length)||node.length<0)throw Error('Invalid graph segment');
        nodes.set(node.id,node);
      });
      graph.links.forEach(function(edge){
        if(!nodes.has(edge.source)||!nodes.has(edge.target)||['+','-'].indexOf(edge.source_orient)<0||['+','-'].indexOf(edge.target_orient)<0||['L','J'].indexOf(edge.kind)<0)throw Error('Graph contains an invalid link');
      });
      index.set(graph.id,{data:graph,nodes:nodes});
    });
    Object.keys(data.contigs).forEach(function(cid){
      var mappings=data.contigs[cid];
      if(!Array.isArray(mappings)||!mappings.length)throw Error('Invalid contig mapping');
      mappings.forEach(function(mapping){
        var graph=index.get(mapping.graph);
        if(!graph||!Array.isArray(mapping.steps)||!mapping.steps.length||mapping.steps.some(function(step){return !Array.isArray(step)||!graph.nodes.has(step[0])||['+','-'].indexOf(step[1])<0;}))throw Error('Contig mapping references an absent segment');
      });
      var i=IDX[cid],length=data.lengths&&data.lengths[cid];
      if(i!==undefined&&length!==undefined&&num('length')[i]!==length)throw Error('Catalogue contig lengths do not match this report: '+cid);
    });
    graphIndex=index;
  }
  function graphApplyTransform(){
    if(graphGroup)graphGroup.setAttribute('transform','translate('+graphTransform.x+' '+graphTransform.y+') scale('+graphTransform.k+')');
  }
  function graphFit(){
    if(!graphNodes.length||!graphSvg)return;
    var w=el('graph-stage').clientWidth||1000,h=540,minX=Infinity,minY=Infinity,maxX=-Infinity,maxY=-Infinity;
    graphNodes.forEach(function(n){minX=Math.min(minX,n.x-n.width/2-15);maxX=Math.max(maxX,n.x+n.width/2+15);minY=Math.min(minY,n.y-35);maxY=Math.max(maxY,n.y+45);});
    graphTransform.k=Math.min(2,(w-40)/Math.max(1,maxX-minX),(h-40)/Math.max(1,maxY-minY));
    graphTransform.x=(w-(maxX+minX)*graphTransform.k)/2;graphTransform.y=(h-(maxY+minY)*graphTransform.k)/2;
    graphApplyTransform();
  }
  function graphEdgePath(edge){
    var a=graphNodeMap.get(edge.a),b=graphNodeMap.get(edge.b);
    var ax=a.x+(edge.source_orient==='+'?a.width/2:-a.width/2),bx=b.x+(edge.target_orient==='+'?-b.width/2:b.width/2);
    if(a===b){var bend=45;return 'M '+ax+' '+a.y+' C '+(ax+40)+' '+(a.y-bend*2)+' '+(bx-40)+' '+(b.y-bend*2)+' '+bx+' '+b.y;}
    var offset=Math.max(40,Math.abs(bx-ax)/2),ac=ax+(edge.source_orient==='+'?offset:-offset),bc=bx+(edge.target_orient==='+'?-offset:offset);
    return 'M '+ax+' '+a.y+' C '+ac+' '+a.y+' '+bc+' '+b.y+' '+bx+' '+b.y;
  }
  function graphUpdatePositions(){
    graphNodes.forEach(function(n){n.element.setAttribute('transform','translate('+n.x+' '+n.y+')');});
    graphEdges.forEach(function(e){e.element.setAttribute('d',graphEdgePath(e));});
  }
  function graphLayout(){
    var visited=new Set(),components=[];
    graphNodes.slice().sort(function(a,b){return (b.owners.length-a.owners.length)||((graphAdj.get(b.key)||[]).length-(graphAdj.get(a.key)||[]).length)||a.key.localeCompare(b.key);}).forEach(function(seed){
      if(visited.has(seed.key))return;
      var queue=[[seed.key,0]],layers=[];visited.add(seed.key);
      for(var q=0;q<queue.length;q++){
        var key=queue[q][0],level=queue[q][1];(layers[level]=layers[level]||[]).push(graphNodeMap.get(key));
        (graphAdj.get(key)||[]).forEach(function(other){if(!visited.has(other)){visited.add(other);queue.push([other,level+1]);}});
      }
      components.push(layers);
    });
    var cx=0,cy=0,rowHeight=0;
    components.forEach(function(layers){
      var width=Math.max(190,layers.length*190),height=Math.max(100,Math.max.apply(null,layers.map(function(x){return x.length;}))*75+30);
      if(cx&&cx+width>1400){cx=0;cy+=rowHeight+40;rowHeight=0;}
      layers.forEach(function(layer,level){layer.forEach(function(n,i){n.x=cx+90+level*190;n.y=cy+40+(i+.5)*height/layer.length;});});
      cx+=width+40;rowHeight=Math.max(rowHeight,height);
    });
  }
  function graphDetails(node){
    var methods=[];
    node.owners.forEach(function(cid){(graphData.contigs[cid]||[]).forEach(function(m){if(m.graph===node.graph)methods.push(m.method);});});
    el('graph-detail').textContent='Source: '+node.graph+'\nSegment: '+node.id+'; '+node.length.toLocaleString()+' bp'+(node.depth!==null&&node.depth!==undefined?'; graph depth '+node.depth.toFixed(2)+'x':'')+
      (node.read_count!==undefined?'; read count '+node.read_count:'')+(node.kmer_count!==undefined?'; k-mer count '+node.kmer_count:'')+
      '\nSelected contigs: '+(node.owners.join(', ')||'Neighbour only')+'\nMapping: '+Array.from(new Set(methods)).join(', ');
  }
  function renderGraph(){
    el('graph-stage').innerHTML='';el('graph-detail').textContent='';el('graph-download').disabled=true;graphSvg=null;graphNodes=[];graphEdges=[];
    try{graphIndexData(graphData);}catch(error){graphMessage(error.message);return;}
    var seeds=new Set(),owners=new Map(),missing=[];
    graphSelection.forEach(function(cid){
      var mappings=graphData.contigs[cid];
      if(!mappings){missing.push(cid+': '+((graphData.unmapped||{})[cid]||'no graph mapping'));return;}
      mappings.forEach(function(mapping){mapping.steps.forEach(function(step){var key=graphKey(mapping.graph,step[0]);seeds.add(key);if(!owners.has(key))owners.set(key,[]);if(owners.get(key).indexOf(cid)<0)owners.get(key).push(cid);});});
    });
    el('graph-missing').textContent=missing.join('\n')||'All selected contigs have graph mappings.';
    el('graph-unmapped').open=missing.length>0;
    var warn=(graphData.warnings||[]).length?' '+graphData.warnings.join(' '):'';
    if(!seeds.size){graphMessage('No source graph mapping is available for the '+graphSelection.length+' selected contig(s). Rebuild this report on the run directory or load its graph catalogue below.'+warn);return;}
    var limit=Number(el('graph-limit').value)||500;
    if(seeds.size>limit){graphMessage(seeds.size+' selected segments exceed the '+limit+' node limit. Select fewer contigs or increase the node limit. No selected segments have been silently dropped.');return;}
    var chosen=new Set(seeds),adj=new Map(),edgePool=[],activeSources=new Set();
    seeds.forEach(function(key){activeSources.add(JSON.parse(key)[0]);});
    graphIndex.forEach(function(graph,gid){if(!activeSources.has(gid))return;graph.data.links.forEach(function(edge){
      var a=graphKey(gid,edge.source),b=graphKey(gid,edge.target);
      if(!adj.has(a))adj.set(a,[]);if(!adj.has(b))adj.set(b,[]);adj.get(a).push(b);adj.get(b).push(a);
      edgePool.push(Object.assign({a:a,b:b,graph:gid},edge));
    });});
    var frontier=Array.from(seeds),omitted=new Set();
    for(var hop=0;hop<Number(el('graph-hops').value);hop++){
      var next=[];
      frontier.forEach(function(key){(adj.get(key)||[]).forEach(function(other){
        if(chosen.has(other))return;if(chosen.size>=limit){omitted.add(other);return;}chosen.add(other);next.push(other);
      });});frontier=next;
    }
    graphNodeMap=new Map();graphAdj=new Map();
    chosen.forEach(function(key){
      var parts=JSON.parse(key),original=graphIndex.get(parts[0]).nodes.get(parts[1]);
      var node=Object.assign({},original,{key:key,graph:parts[0],owners:owners.get(key)||[],width:Math.max(42,Math.min(135,25+19*Math.log10(Math.max(10,original.length))))});
      graphNodes.push(node);graphNodeMap.set(key,node);graphAdj.set(key,[]);
    });
    graphEdges=edgePool.filter(function(edge){return chosen.has(edge.a)&&chosen.has(edge.b);});
    graphEdges.forEach(function(edge){graphAdj.get(edge.a).push(edge.b);graphAdj.get(edge.b).push(edge.a);});
    graphLayout();
    var width=el('graph-stage').clientWidth||1000;
    graphSvg=svgElement('svg',{xmlns:SVG_NS,viewBox:'0 0 '+width+' 540',role:'img','aria-label':'Assembly graph from GFA segments and links'});
    var defs=svgElement('defs'),pattern=svgElement('pattern',{id:'graph-shared',width:8,height:8,patternUnits:'userSpaceOnUse'});
    pattern.appendChild(svgElement('rect',{width:8,height:8,fill:'#e2e8f0'}));pattern.appendChild(svgElement('path',{d:'M0 8L8 0',stroke:'#334155','stroke-width':2}));defs.appendChild(pattern);graphSvg.appendChild(defs);
    graphGroup=svgElement('g');graphSvg.appendChild(graphGroup);
    graphEdges.forEach(function(edge){
      var path=svgElement('path',{fill:'none',stroke:edge.kind==='J'?'#b45309':'#94a3b8','stroke-width':1.7,'stroke-dasharray':edge.kind==='J'?'6 4':''});
      path.appendChild(svgElement('title',{},edge.source+edge.source_orient+' -> '+edge.target+edge.target_orient+'; '+edge.kind+' '+edge.overlap+'; '+edge.graph));
      edge.element=path;graphGroup.appendChild(path);
    });
    graphNodes.forEach(function(node){
      var group=svgElement('g',{'class':'graph-node','data-segment':node.id,'data-source':node.graph});
      var colours=Array.from(new Set(node.owners.map(function(cid){return IDX[cid]===undefined?'#64748b':pointColour(IDX[cid]);})));
      var fill=colours.length>1?'url(#graph-shared)':colours[0]||'#d1d9e0';
      group.appendChild(svgElement('rect',{x:-node.width/2,y:-12,width:node.width,height:24,rx:5,fill:fill,stroke:node.owners.length?'#16222b':'#94a3b8','stroke-width':node.owners.length?2:1}));
      group.appendChild(svgElement('title',{},node.graph+' / '+node.id+'; '+node.length+' bp; '+(node.owners.join(', ')||'neighbour')));
      if(el('graph-labels').value==='yes')group.appendChild(svgElement('text',{x:0,y:30,'text-anchor':'middle','font-size':11,fill:'#334155'},node.id));
      group.onpointerdown=function(event){event.stopPropagation();graphDrag={node:node,x:event.clientX,y:event.clientY};graphSvg.setPointerCapture(event.pointerId);graphDetails(node);};
      group.onclick=function(){graphDetails(node);};node.element=group;graphGroup.appendChild(group);
    });
    graphSvg.onpointerdown=function(event){graphDrag={x:event.clientX,y:event.clientY};graphSvg.setPointerCapture(event.pointerId);};
    graphSvg.onpointermove=function(event){
      if(!graphDrag)return;var dx=event.clientX-graphDrag.x,dy=event.clientY-graphDrag.y;graphDrag.x=event.clientX;graphDrag.y=event.clientY;
      if(graphDrag.node){graphDrag.node.x+=dx/graphTransform.k;graphDrag.node.y+=dy/graphTransform.k;graphUpdatePositions();}
      else{graphTransform.x+=dx;graphTransform.y+=dy;graphApplyTransform();}
    };
    graphSvg.onpointerup=graphSvg.onpointercancel=function(){graphDrag=null;};
    graphSvg.addEventListener('wheel',function(event){
      event.preventDefault();var rect=graphSvg.getBoundingClientRect(),x=event.clientX-rect.left,y=event.clientY-rect.top,old=graphTransform.k;
      graphTransform.k=Math.max(.02,Math.min(20,old*Math.exp(-event.deltaY*.0015)));
      graphTransform.x=x-(x-graphTransform.x)*graphTransform.k/old;graphTransform.y=y-(y-graphTransform.y)*graphTransform.k/old;graphApplyTransform();
    },{passive:false});
    el('graph-stage').appendChild(graphSvg);graphUpdatePositions();graphFit();el('graph-download').disabled=false;
    var mapped=graphSelection.length-missing.length,sources=new Set(graphNodes.map(function(n){return n.graph;}));
    graphMessage(mapped+' / '+graphSelection.length+' selected contigs mapped; '+sources.size+' source graph(s); '+graphNodes.length+' segments; '+graphEdges.length+' links.'+
      (omitted.size?' Neighbourhood truncated at '+limit+' nodes.':'')+(!graphEdges.length?' No links connect the displayed segments in the source GFA.':'')+warn);
  }
  el('graph-open').onclick=function(){graphSelection=Array.from(selected);el('graph-dialog').showModal();renderGraph();};
  el('graph-close').onclick=function(){el('graph-dialog').close();};
  ['graph-hops','graph-limit','graph-labels'].forEach(function(id){el(id).onchange=renderGraph;});
  el('graph-fit').onclick=graphFit;
  el('graph-download').onclick=function(){
    if(!graphSvg)return;var svg=new XMLSerializer().serializeToString(graphSvg);download(new Blob([svg],{type:'image/svg+xml'}),'selected_assembly_graph.svg');
  };
  el('graph-file').onchange=async function(){
    var file=this.files&&this.files[0];if(!file)return;
    if(file.size>100*1024*1024){graphMessage('Catalogue exceeds the 100 MB local-load limit. Build a smaller source catalogue.');return;}
    try{var data=JSON.parse(await file.text());graphIndexData(data);graphData=data;renderGraph();}catch(error){graphMessage('Graph catalogue could not be loaded: '+error.message);}
  };
'''


if __name__ == '__main__':
    main()

