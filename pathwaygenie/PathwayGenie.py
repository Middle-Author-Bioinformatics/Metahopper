#!/usr/bin/env python3
"""PathwayGenie: Prodigal -> HMMER -> activity matrices and flexible neighborhoods.

Python 3.9+ standard library. External dependencies: Prodigal and hmmsearch.
This is a standalone workflow inspired by the supplied GasGenie example.
"""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from dataclasses import dataclass, field
from datetime import datetime, timezone
from functools import lru_cache
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

VERSION = '1.0.0'


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def write_json(path, value):
    path = Path(path)
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    temp.replace(path)


def write_csv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp')
    with temp.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)


def fasta_records(path):
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    header = None
    sequence = []
    with opener(path, 'rt') as handle:
        for number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(sequence)
                header = line[1:]
                if not header.strip():
                    raise ValueError(f'Empty FASTA identifier at {path}:{number}')
                sequence = []
            else:
                if header is None:
                    raise ValueError(f'Sequence before FASTA header at {path}:{number}')
                sequence.append(''.join(line.split()))
    if header is not None:
        yield header, ''.join(sequence)


@dataclass
class ORF:
    protein: str
    contig: str
    original_contig: str
    start: int
    end: int
    strand: int
    partial: str
    length: int
    rank: int = 0


@dataclass
class Hit:
    protein: str
    accession: str
    target: str
    source_group: str
    call_type: str
    seq_score: float
    domain_score: float
    seq_evalue: float
    independent_evalue: float
    hmm_coverage: float
    protein_intervals: list
    hmm_intervals: list
    margin: float
    components: list
    status: str = 'accepted'
    reason: str = ''


def union_length(intervals):
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(end, merged[-1][1])
    return sum(end - start + 1 for start, end in merged)


def overlap_fraction(left, right):
    intersections = [(max(a, c), min(b, d)) for a, b in left for c, d in right
                     if max(a, c) <= min(b, d)]
    denominator = min(union_length(left), union_length(right))
    return union_length(intersections) / denominator if denominator else 0.0


def model_headers(path):
    result = {}
    record = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith('HMMER3/'):
                record = {}
            elif line.strip() == '//':
                accession = record.get('ACC')
                if not accession or accession in result:
                    raise ValueError(f'Missing/duplicate model accession in {path}: {accession}')
                if record.get('ALPH') != 'amino':
                    raise ValueError(f'Non-protein HMM in {path}: {accession}')
                result[accession] = record
                record = {}
            elif line.startswith(('NAME ', 'ACC ', 'LENG ', 'ALPH ', 'TC ', 'GA ')):
                key, value = line.split(None, 1)
                record[key] = value.strip()
    return result


def load_panel(root, config, include_context=True):
    root = Path(root)
    mapping_path = root / 'gene_model_map.tsv'
    manifest_path = root / 'download_manifest.json'
    with mapping_path.open() as handle:
        mapping = list(csv.DictReader(handle, delimiter='\t'))
    manifest = json.loads(manifest_path.read_text())
    original = {m['accession']: m for m in manifest['models']}
    specs = [('main', 'libraries/pathway_gene_models.hmm', 'trusted')]
    if include_context:
        specs += [('context_ncbi', 'libraries/context_ncbi.hmm', 'trusted'),
                  ('context_pfam', 'libraries/context_pfam.hmm', 'gathering')]
    libraries = []
    models = {}
    for group, relative, cutoff in specs:
        path = root / relative
        headers = model_headers(path)
        if not headers:
            raise ValueError(f'Empty HMM library: {path}')
        for accession, header in headers.items():
            if accession in models or accession not in original:
                raise ValueError(f'Unknown/duplicate library model: {accession}')
            if accession not in config['model_rules']:
                raise ValueError(f'Missing config rule for {accession}')
            native = original[accession]
            if (native['group'] == 'context_only') != bool(config['model_rules'][accession]['context_only']):
                raise ValueError('Configuration/source context classification mismatch: ' + accession)
            source_path = root / native['relative_file']
            if sha256(source_path) != native['sha256']:
                raise ValueError(f'Source HMM checksum mismatch: {source_path}')
            # Bundles may normalize NAME only; parameters must match their source files.
            normalized = lambda s: re.sub(r'^NAME\s+.*$', 'NAME ACCESSION', s, flags=re.M).strip()
            text = source_path.read_text()
            models[accession] = dict(native, header=header, search_group=group,
                                     cutoff_mode=cutoff)
            if header['ACC'] != native['accession'] or int(header['LENG']) != native['length']:
                raise ValueError(f'Library/source metadata mismatch: {accession}')
        # Validate all model bodies too, not just bundle headers.
        bundle_records = re.split(r'(?m)^//\s*$', path.read_text())
        for text in bundle_records:
            if not text.strip():
                continue
            match = re.search(r'^ACC\s+(\S+)', text, re.M)
            source_text = (root / original[match.group(1)]['relative_file']).read_text()
            if normalized(text) != normalized(source_text.replace('\n//', '\n')):
                raise ValueError(f'Bundle parameters differ from source: {match.group(1)}')
        libraries.append({'group': group, 'path': str(path.resolve()),
                          'cutoff_mode': cutoff, 'sha256': sha256(path)})
    mapped = {r['model_accession'] for r in mapping}
    if set(models) - mapped:
        raise ValueError('Library models missing from gene_model_map.tsv')
    return models, mapping, libraries, manifest


def validate_config(config):
    if config.get('schema_version') != 1:
        raise ValueError('Unsupported configuration schema')
    components = set(config['components'])
    for accession, rule in config['model_rules'].items():
        if rule['context_only'] and rule['components']:
            raise ValueError('Context-only models cannot confirm components: ' + accession)
    ids = [p['id'] for p in config['pathways']]
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate pathway IDs')
    for pathway in config['pathways']:
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', pathway['id']):
            raise ValueError('Unsafe pathway ID: ' + pathway['id'])
        for key in ['scored_components', 'route_components', 'extra_components']:
            values = pathway[key]
            if len(values) != len(set(values)) or set(values) - components:
                raise ValueError(f'Unknown/duplicate components in {pathway["id"]}/{key}')
        if not pathway['scored_components']:
            raise ValueError('Every pathway needs scored components')
        for component in pathway['scored_components']:
            if config['components'][component]['assay_scope'] != 'targeted':
                raise ValueError('Unassayed component cannot enter heatmap denominator: ' + component)
            if not any(component in r['components'] and not r['context_only']
                       for r in config['model_rules'].values()):
                raise ValueError('Scored component has no confirming model: ' + component)
    for accession, rule in config['model_rules'].items():
        if set(rule['components']) - components:
            raise ValueError('Unknown model components: ' + accession)
        if not 0 <= rule.get('min_hmm_coverage', 0) <= 1:
            raise ValueError('Invalid per-model coverage: ' + accession)
    for template in config['operon_templates']:
        if template['pathway'] not in ids:
            raise ValueError('Unknown operon pathway')
        if any(set(slot) - components or not slot for slot in template['slots']):
            raise ValueError('Invalid operon slot: ' + template['id'])
        if not 2 <= template['min_slots'] <= len(template['slots']):
            raise ValueError('Invalid minimum slot count: ' + template['id'])
        if not 0 <= template['min_fraction'] <= 1 or template['max_missing'] < 0:
            raise ValueError('Invalid operon tolerance: ' + template['id'])


def normalize_bin(input_path, output_path, contig_table):
    contigs = {}
    originals = set()
    rows = []
    temp = output_path.with_suffix('.tmp')
    with temp.open('w') as handle:
        for index, (header, sequence) in enumerate(fasta_records(input_path), 1):
            original = header.split()[0]
            if original in originals:
                raise ValueError(f'Duplicate contig ID: {original}')
            originals.add(original)
            sequence = sequence.upper()
            if not sequence or set(sequence) - set('ACGTRYSWKMBDHVN'):
                raise ValueError(f'Empty/non-DNA sequence in contig {original}')
            cid = f'c{index:07d}'
            contigs[cid] = {'original': original, 'length': len(sequence)}
            rows.append({'contig': cid, 'original_contig': original,
                         'length': len(sequence), 'original_header': header})
            handle.write('>' + cid + '\n')
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + '\n')
    if not contigs:
        temp.unlink(missing_ok=True)
        raise ValueError(f'No FASTA records in {input_path}')
    temp.replace(output_path)
    write_csv(contig_table, rows, ['contig', 'original_contig', 'length', 'original_header'])
    return contigs


def load_contigs(path):
    with path.open() as handle:
        return {r['contig']: {'original': r['original_contig'], 'length': int(r['length'])}
                for r in csv.DictReader(handle)}


def parse_proteins(path, contigs, cleaned=None):
    genes = {}
    sequences = []
    for header, sequence in fasta_records(path):
        parts = header.split(' # ')
        if len(parts) < 5:
            raise ValueError('Prodigal protein header lacks coordinates: ' + header)
        protein = parts[0].split()[0]
        contig = protein.rsplit('_', 1)[0]
        if protein in genes or contig not in contigs:
            raise ValueError('Duplicate protein/unknown contig: ' + protein)
        start, end, strand = map(int, parts[1:4])
        if not 1 <= start <= end <= contigs[contig]['length'] or strand not in [-1, 1]:
            raise ValueError('Invalid ORF coordinates: ' + protein)
        partial = re.search(r'(?:^|;)partial=([01]{2})(?:;|$)', parts[4])
        sequence = sequence.upper().rstrip('*')
        if not sequence or '*' in sequence:
            raise ValueError('Empty protein/internal stop: ' + protein)
        genes[protein] = ORF(protein, contig, contigs[contig]['original'], start, end,
                             strand, partial.group(1) if partial else 'unknown', len(sequence))
        sequences.append((header, sequence))
    by_contig = defaultdict(list)
    for gene in genes.values():
        by_contig[gene.contig].append(gene)
    for group in by_contig.values():
        for rank, gene in enumerate(sorted(group, key=lambda g: (g.start, g.end, g.protein))):
            gene.rank = rank
    if cleaned:
        temp = cleaned.with_suffix('.tmp')
        with temp.open('w') as handle:
            for header, sequence in sequences:
                handle.write('>' + header + '\n' + sequence + '\n')
        temp.replace(cleaned)
    return genes


def parse_domtbl(path, models, genes, config, min_coverage):
    grouped = defaultdict(list)
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.split(maxsplit=22)
        if len(fields) < 22:
            raise ValueError(f'Malformed domtblout at {path}:{number}')
        protein, accession = fields[0], fields[4] if fields[4] != '-' else fields[3]
        if protein not in genes or accession not in models:
            raise ValueError(f'Unknown protein/model at {path}:{number}: {protein}/{accession}')
        gene, model = genes[protein], models[accession]
        if int(fields[2]) != gene.length or int(fields[5]) != model['length']:
            raise ValueError(f'Protein/model length mismatch at {path}:{number}')
        seq_e, seq_score, dom_e, dom_score = map(float, [fields[6], fields[7], fields[12], fields[13]])
        if not all(math.isfinite(x) for x in [seq_e, seq_score, dom_e, dom_score]):
            raise ValueError('Non-finite HMMER score')
        hmm_start, hmm_end, ali_start, ali_end = map(int, fields[15:19])
        if not (1 <= hmm_start <= hmm_end <= model['length'] and
                1 <= ali_start <= ali_end <= gene.length):
            raise ValueError('Invalid HMMER domain coordinates')
        cutoff = 'GA' if model['cutoff_mode'] == 'gathering' else 'TC'
        seq_cut, domain_cut = [float(v.rstrip(';')) for v in model['header'][cutoff].split()[:2]]
        if seq_score < seq_cut or dom_score < domain_cut:
            continue
        grouped[(protein, accession)].append((seq_score, dom_score, seq_e, dom_e,
                                               (hmm_start, hmm_end), (ali_start, ali_end),
                                               seq_cut, domain_cut))
    hits = []
    for (protein, accession), domains in sorted(grouped.items()):
        model = models[accession]
        rule = config['model_rules'][accession]
        hmm_intervals = [r[4] for r in domains]
        protein_intervals = [r[5] for r in domains]
        coverage = union_length(hmm_intervals) / model['length']
        best = max(domains, key=lambda r: r[1])
        margin = min((best[0] - best[6]) / max(abs(best[6]), 1),
                     (best[1] - best[7]) / max(abs(best[7]), 1))
        hit = Hit(protein, accession, rule['target'], model['group'], rule['target_call'],
                  best[0], best[1], best[2], min(d[3] for d in domains), coverage,
                  protein_intervals, hmm_intervals, margin, list(rule['components']))
        required = max(min_coverage, rule.get('min_hmm_coverage', 0))
        if coverage < required:
            hit.status = 'fragment_candidate'
            hit.reason = f'HMM coverage {coverage:.3f} below required {required:.3f}'
        elif rule['context_only']:
            hit.status = 'context_candidate'
            hit.reason = 'Supporting/alternative model; never confirms a component by default'
        hits.append(hit)
    return hits


def resolve_hits(hits, ambiguity_margin=0.10, overlap=0.5):
    by_protein = defaultdict(list)
    for hit in hits:
        if hit.status == 'accepted':
            by_protein[hit.protein].append(hit)
    for group in by_protein.values():
        retained = []
        for hit in sorted(group, key=lambda h: (-h.margin, -h.hmm_coverage, h.accession)):
            conflicting = []
            for earlier in retained:
                if earlier.status == 'competing_hit':
                    continue
                if earlier.target == hit.target:
                    continue
                if overlap_fraction(earlier.protein_intervals, hit.protein_intervals) < overlap:
                    continue
                if set(earlier.components) & set(hit.components):
                    # Preserve reaction evidence; an overlapping shared family is not
                    # independent evidence of two differently named isozymes.
                    if earlier.call_type == hit.call_type == 'specific':
                        if abs(earlier.margin - hit.margin) <= ambiguity_margin:
                            earlier.call_type = hit.call_type = 'shared_family'
                    continue
                conflicting.append(earlier)
            if conflicting:
                close = [h for h in conflicting if h.margin - hit.margin <= ambiguity_margin]
                if close:
                    hit.status = 'ambiguous'
                    hit.reason = 'Overlapping incompatible models have similar normalized cutoff margins'
                    for earlier in close:
                        earlier.status = 'ambiguous'
                        earlier.reason = hit.reason
                else:
                    hit.status = 'competing_hit'
                    hit.reason = 'Stronger incompatible overlapping model: ' + conflicting[0].accession
            retained.append(hit)
    return hits


def component_evidence(hits, config):
    accepted = defaultdict(set)
    candidates = defaultdict(set)
    for hit in hits:
        for component in hit.components:
            destination = accepted if hit.status == 'accepted' else candidates
            destination[component].add(hit.protein)
    result = {}
    for component, definition in config['components'].items():
        proteins = accepted[component]
        if proteins:
            state, value, count = 'detected', 1, len(proteins)
        elif definition['assay_scope'] != 'targeted':
            state = 'not_assayed' if definition['assay_scope'] == 'not_assayed' else 'fusion_only_not_detected'
            value, count = 'NA', 'NA'
        else:
            state = 'candidate_only' if candidates[component] else 'not_detected'
            value, count = 0, 0
        result[component] = {'state': state, 'value': value, 'count': count,
                             'proteins': proteins, 'candidates': candidates[component]}
    return result


def target_evidence(hits, targets):
    result = {}
    for target in targets:
        selected = [h for h in hits if h.target == target]
        good = [h for h in selected if h.status == 'accepted']
        calls = {h.call_type for h in good}
        if 'specific' in calls:
            state = 'detected'
        elif 'shared_family' in calls:
            state = 'shared_family'
        elif 'domain_only' in calls:
            state = 'domain_only'
        elif selected:
            state = 'candidate_only'
        else:
            state = 'not_detected'
        result[target] = state
    return result


def neighborhood_windows(genes, max_intervening, max_bp):
    """Bounded windows: never bridge contigs or ignore intervening unmarked ORFs."""
    by_contig = defaultdict(list)
    for gene in genes:
        by_contig[gene.contig].append(gene)
    windows = []
    for contig_genes in by_contig.values():
        ordered = sorted(contig_genes, key=lambda g: g.rank)
        for start in range(len(ordered)):
            window = [ordered[start]]
            for gene in ordered[start + 1:]:
                if gene.rank - window[-1].rank - 1 > max_intervening:
                    break
                if max(g.end for g in window + [gene]) - window[0].start + 1 > max_bp:
                    break
                window.append(gene)
            if len(window) < 2:
                continue
            ids = {g.protein for g in window}
            if any(ids <= {g.protein for g in prior} for prior in windows):
                continue
            windows.append(window)
    return windows


def ordered_slot_match(ordered_components, slots):
    """LCS with missing slots/extra genes; a verified fusion may fill adjacent slots."""
    @lru_cache(None)
    def solve(i, j):
        if i == len(ordered_components) or j == len(slots):
            return ()
        alternatives = [solve(i + 1, j), solve(i, j + 1)]
        if set(ordered_components[i]) & set(slots[j]):
            alternatives += [((i, j),) + solve(i + 1, j + 1)]
            # Same ORF may encode distinct activities. Do not reuse the same
            # component to satisfy repeated biological slots.
            for k in range(j + 1, len(slots)):
                if (set(ordered_components[i]) & set(slots[k]) and
                        not set(slots[j]) & set(slots[k])):
                    alternatives.append(((i, j), (i, k)) + solve(i + 1, k + 1))
        return max(alternatives, key=len)
    return solve(0, 0)


def operon_metrics(window, component_by_protein, template, args):
    slots = template['slots']
    observed = {j for j, slot in enumerate(slots)
                if any(set(slot) & component_by_protein[g.protein] for g in window)}
    if len(observed) < template.get('min_slots', 2):
        return None
    min_fraction = args.operon_min_fraction if args.operon_min_fraction is not None else template['min_fraction']
    max_missing = args.operon_max_missing if args.operon_max_missing is not None else template['max_missing']
    coverage = len(observed) / len(slots)
    if coverage < min_fraction or len(slots) - len(observed) > max_missing:
        return None
    options = []
    for orientation in [1, -1]:
        ordered = sorted(window, key=lambda g: g.rank, reverse=orientation == -1)
        if args.operon_order == 'ignore':
            matches = tuple((next(i for i, g in enumerate(ordered)
                                  if set(slots[j]) & component_by_protein[g.protein]), j)
                            for j in sorted(observed))
        else:
            matches = ordered_slot_match([tuple(sorted(component_by_protein[g.protein])) for g in ordered], slots)
        if not matches:
            continue
        matched_genes = {ordered[i].protein for i, j in matches}
        if len(matched_genes) < 2:
            continue
        strand_score = sum(g.strand == orientation for g in ordered if g.protein in matched_genes) / len(matched_genes)
        order_score = len(matches) / len(observed)
        if args.operon_order == 'strict' and order_score < 1:
            continue
        if args.operon_strand == 'strict' and strand_score < 1:
            continue
        score = coverage
        if args.operon_order != 'ignore':
            score *= order_score
        if args.operon_strand != 'ignore':
            score *= strand_score
        options.append({'score': score, 'coverage': coverage, 'order_score': order_score,
                        'strand_score': strand_score, 'orientation': orientation,
                        'matched_slots': ';'.join(str(j + 1) for i, j in matches),
                        'matched_proteins': ';'.join(ordered[i].protein for i, j in matches),
                        'observed_slots': len(observed), 'total_slots': len(slots)})
    return max(options, key=lambda r: (r['score'], r['order_score'], r['strand_score'])) if options else None


def external_version(executable, flag):
    result = subprocess.run([executable, flag], stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, timeout=20)
    text = result.stdout.strip()
    if not text:
        raise RuntimeError('Cannot obtain tool version: ' + executable)
    return '\n'.join(text.splitlines()[:5])


def run_command(command, logfile):
    with open(logfile, 'w') as handle:
        handle.write('COMMAND_JSON: ' + json.dumps(command) + '\n')
        handle.flush()
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode:
        tail = '\n'.join(Path(logfile).read_text(errors='replace').splitlines()[-8:])
        raise RuntimeError(f'{Path(command[0]).name} failed ({result.returncode}); {logfile}\n{tail}')


def cache_valid(checkpoint, signature, files):
    try:
        old = json.loads(Path(checkpoint).read_text())
        return old['signature'] == signature and all(
            p.is_file() and old['files'].get(p.name) == sha256(p) for p in files)
    except (OSError, ValueError, KeyError):
        return False


def save_checkpoint(path, signature, files):
    write_json(path, {'signature': signature,
                     'files': {p.name: sha256(p) for p in files}})


def analyze_neighborhoods(bin_id, genes, evidence, config, args):
    by_protein = defaultdict(set)
    for component, record in evidence.items():
        for protein in record['proteins']:
            by_protein[protein].add(component)
    neighborhoods, operons = [], []
    for pathway in config['pathways']:
        pid = pathway['id']
        relevant = set(pathway['scored_components'] + pathway['route_components'] + pathway['extra_components'])
        selected = [genes[p] for p, cs in by_protein.items() if cs & relevant]
        windows = neighborhood_windows(selected, args.max_intervening_genes, args.max_neighborhood_bp)
        templates = [t for t in config['operon_templates'] if t['pathway'] == pid]
        for index, window in enumerate(windows, 1):
            components = set().union(*(by_protein[g.protein] & relevant for g in window))
            if len(components) < 2:
                continue
            nid = f'{pid}_n{index:04d}'
            base = {'bin_id': bin_id, 'pathway': pid, 'neighborhood_id': nid,
                    'contig': window[0].original_contig,
                    'start': min(g.start for g in window), 'end': max(g.end for g in window),
                    'protein_count': len(window),
                    'proteins': ';'.join(g.protein for g in window),
                    'components': ';'.join(sorted(components)),
                    'max_intervening_orfs': max(b.rank - a.rank - 1 for a, b in zip(window, window[1:])),
                    'observed_order': ' | '.join(','.join(sorted(by_protein[g.protein] & relevant)) +
                                                ('(+)' if g.strand == 1 else '(-)') for g in window)}
            neighborhoods.append(base)
            for template in templates:
                metric = operon_metrics(window, by_protein, template, args)
                if metric:
                    operons.append(dict(base, **metric, template_id=template['id'],
                                        reference_taxon=template['reference_taxon'],
                                        source_url=template['source_url'], template_notes=template['notes']))
    return neighborhoods, operons


def run_bin(item, out, models, libraries, config, args, versions):
    bin_id, source = item
    slug = re.sub(r'[^A-Za-z0-9_.-]+', '_', bin_id).strip('.')[:70] or 'bin'
    folder = out / 'bins' / (slug + '_' + fingerprint(bin_id)[:8])
    folder.mkdir(parents=True, exist_ok=True)
    normalized, contig_file = folder / 'contigs.fna', folder / 'contigs.csv'
    raw, proteins = folder / 'proteins.raw.faa', folder / 'proteins.faa'
    gff, coding = folder / 'genes.gff', folder / 'genes.fna'
    orf_files = [normalized, contig_file, raw, proteins, gff, coding]
    orf_signature = fingerprint({'input_sha256': sha256(source), 'mode': args.prodigal_mode,
                                 'prodigal': versions['prodigal'], 'pipeline_orf_schema': 1})
    reused_orfs = not args.force and cache_valid(folder / 'orfs.done.json', orf_signature, orf_files)
    if reused_orfs:
        contigs = load_contigs(contig_file)
        genes = parse_proteins(proteins, contigs)
    else:
        contigs = normalize_bin(source, normalized, contig_file)
        temporary = [p.with_name(p.name + '.tmp') for p in [raw, gff, coding]]
        command = [args.prodigal, '-i', str(normalized), '-a', str(temporary[0]),
                   '-o', str(temporary[1]), '-d', str(temporary[2]), '-f', 'gff',
                   '-p', args.prodigal_mode, '-q']
        run_command(command, folder / 'prodigal.log')
        for temp, final in zip(temporary, [raw, gff, coding]):
            if not temp.is_file():
                raise RuntimeError('Prodigal did not create expected file: ' + str(temp))
            temp.replace(final)
        genes = parse_proteins(raw, contigs, cleaned=proteins)
        save_checkpoint(folder / 'orfs.done.json', orf_signature, orf_files)
    write_csv(folder / 'orfs.csv', [vars(g) for g in genes.values()],
              ['protein', 'contig', 'original_contig', 'start', 'end', 'strand', 'partial', 'length', 'rank'])
    hits = []
    reused_searches = 0
    for library in libraries:
        stem = folder / ('hmm_' + library['group'])
        dom, table, report = Path(str(stem) + '.domtbl'), Path(str(stem) + '.tbl'), Path(str(stem) + '.txt')
        checkpoint = Path(str(stem) + '.done.json')
        signature = fingerprint({'proteins': sha256(proteins), 'hmm': library['sha256'],
                                 'mode': library['cutoff_mode'], 'hmmer': versions['hmmsearch']})
        cached = not args.force and cache_valid(checkpoint, signature, [dom, table, report])
        if cached:
            reused_searches += 1
        elif genes:
            temps = [p.with_name(p.name + '.tmp') for p in [dom, table, report]]
            flag = '--cut_ga' if library['cutoff_mode'] == 'gathering' else '--cut_tc'
            command = [args.hmmsearch, '--cpu', str(args.cpus_per_job), flag, '--noali', '--notextw',
                       '--domtblout', str(temps[0]), '--tblout', str(temps[1]),
                       '-o', str(temps[2]), library['path'], str(proteins)]
            run_command(command, Path(str(stem) + '.command.log'))
            for temp, final in zip(temps, [dom, table, report]):
                if not temp.is_file():
                    raise RuntimeError('HMMER did not create expected output: ' + str(temp))
                temp.replace(final)
            save_checkpoint(checkpoint, signature, [dom, table, report])
        else:
            for path in [dom, table, report]:
                path.write_text('# No predicted proteins; search not run.\n')
            save_checkpoint(checkpoint, signature, [dom, table, report])
        hits.extend(parse_domtbl(dom, models, genes, config, args.min_hmm_coverage))
    resolve_hits(hits, args.ambiguity_margin, args.overlap_fraction)
    evidence = component_evidence(hits, config)
    target_ids = sorted({rule['target'] for rule in config['model_rules'].values()})
    targets = target_evidence(hits, target_ids)
    for target in target_ids:
        if not any(rule['target'] == target and not rule['context_only']
                   for accession, rule in config['model_rules'].items() if accession in models):
            targets[target] = 'candidate_only' if any(h.target == target for h in hits) else 'not_assayed'
    if any(h.status == 'accepted' and h.accession == 'TIGR00118.1' for h in hits):
        for target in ['ilvB', 'ilvG']:
            if targets.get(target) == 'not_detected':
                targets[target] = 'shared_family'
    neighborhoods, operons = analyze_neighborhoods(bin_id, genes, evidence, config, args)
    result = {'bin_id': bin_id, 'source': str(source), 'status': 'ok' if genes else 'no_predicted_genes',
              'error': '', 'contigs': len(contigs), 'bases': sum(c['length'] for c in contigs.values()),
              'orfs': len(genes), 'genes': genes, 'hits': hits, 'evidence': evidence,
              'targets': targets, 'neighborhoods': neighborhoods, 'operons': operons,
              'folder': str(folder), 'orfs_reused': reused_orfs, 'searches_reused': reused_searches}
    write_json(folder / 'bin_status.json', {k: result[k] for k in ['bin_id', 'source', 'status', 'error',
               'contigs', 'bases', 'orfs', 'folder', 'orfs_reused', 'searches_reused']})
    return result


HIT_FIELDS = ['bin_id', 'protein', 'contig', 'start', 'end', 'strand', 'orf_rank', 'partial',
              'protein_length', 'model_accession', 'source_database', 'source_product',
              'target_gene', 'target_call', 'source_group', 'status', 'reason', 'components',
              'sequence_bit_score', 'domain_bit_score', 'sequence_evalue', 'independent_evalue',
              'hmm_coverage', 'hmm_intervals', 'protein_intervals', 'normalized_cutoff_margin']
NEIGHBOR_FIELDS = ['bin_id', 'pathway', 'neighborhood_id', 'contig', 'start', 'end',
                   'protein_count', 'proteins', 'components', 'max_intervening_orfs', 'observed_order']
OPERON_FIELDS = NEIGHBOR_FIELDS + ['template_id', 'reference_taxon', 'score', 'coverage', 'order_score',
                                  'strand_score', 'orientation', 'matched_slots', 'matched_proteins',
                                  'observed_slots', 'total_slots', 'source_url', 'template_notes']


def hit_row(result, hit, models, mapping):
    gene = result['genes'][hit.protein]
    source = mapping[hit.accession]
    return {'bin_id': result['bin_id'], 'protein': gene.protein, 'contig': gene.original_contig,
            'start': gene.start, 'end': gene.end, 'strand': gene.strand, 'orf_rank': gene.rank,
            'partial': gene.partial, 'protein_length': gene.length, 'model_accession': hit.accession,
            'source_database': models[hit.accession]['source_database'],
            'source_product': source['source_product'], 'target_gene': hit.target,
            'target_call': hit.call_type, 'source_group': hit.source_group,
            'status': hit.status, 'reason': hit.reason, 'components': ';'.join(hit.components),
            'sequence_bit_score': hit.seq_score, 'domain_bit_score': hit.domain_score,
            'sequence_evalue': hit.seq_evalue, 'independent_evalue': hit.independent_evalue,
            'hmm_coverage': round(hit.hmm_coverage, 4), 'hmm_intervals': json.dumps(hit.hmm_intervals),
            'protein_intervals': json.dumps(hit.protein_intervals), 'normalized_cutoff_margin': round(hit.margin, 4)}


def summarize_pathway(result, pathway, config):
    pid, bid = pathway['id'], result['bin_id']
    if result['status'] == 'failed':
        return {'bin_id': bid, 'pathway': pid, 'panel_percent': 'NA', 'panel_status': 'bin_failed',
                'route_status': 'bin_failed', 'operon_support': 'NA', 'error': result['error']}
    evidence = result['evidence']
    scored, route = pathway['scored_components'], pathway['route_components']
    detected = [c for c in scored if evidence[c]['value'] == 1]
    missing = [c for c in scored if evidence[c]['value'] != 1]
    unresolved = [c for c in route if evidence[c]['value'] == 'NA']
    candidates = [c for c in scored if evidence[c]['state'] == 'candidate_only']
    percent = round(100 * len(detected) / len(scored), 3)
    panel_status = 'complete_target_panel' if not missing else 'partial_target_panel' if detected else 'not_detected_in_panel'
    if all(evidence[c]['value'] == 1 for c in route):
        route_status = 'supported_selected_route'
    elif not detected:
        route_status = 'not_detected_in_panel'
    elif unresolved:
        route_status = 'unresolved_unmeasured_steps'
    else:
        route_status = 'partial_evidence'
    templates = [t for t in config['operon_templates'] if t['pathway'] == pid]
    candidates_op = [op['score'] for op in result['operons'] if op['pathway'] == pid]
    support = round(max(candidates_op, default=0), 4) if templates else 'NA'
    return {'bin_id': bid, 'pathway': pid, 'panel_percent': percent, 'panel_status': panel_status,
            'route_status': route_status, 'detected_scored_components': len(detected),
            'total_scored_components': len(scored), 'detected_components': ';'.join(detected),
            'missing_scored_components': ';'.join(missing), 'candidate_components': ';'.join(candidates),
            'unmeasured_route_steps': ';'.join(unresolved), 'operon_support': support,
            'pathway_notes': pathway['notes'], 'error': ''}


def write_outputs(out, results, config, models, mapping):
    results = sorted(results, key=lambda r: r['bin_id'])
    pids = [p['id'] for p in config['pathways']]
    summaries = [summarize_pathway(r, p, config) for r in results for p in config['pathways']]
    by_pair = {(r['bin_id'], r['pathway']): r for r in summaries}
    for filename, fieldname in [('pathway_heatmap.csv', 'panel_percent'),
                                ('pathway_status.csv', 'panel_status'),
                                ('route_status.csv', 'route_status'),
                                ('operon_support_heatmap.csv', 'operon_support')]:
        rows = [dict(bin_id=r['bin_id'], **{p: by_pair[(r['bin_id'], p)][fieldname] for p in pids})
                for r in results]
        write_csv(out / filename, rows, ['bin_id'] + pids)
    summary_fields = ['bin_id', 'pathway', 'panel_percent', 'panel_status', 'route_status',
                      'detected_scored_components', 'total_scored_components', 'detected_components',
                      'missing_scored_components', 'candidate_components', 'unmeasured_route_steps',
                      'operon_support', 'pathway_notes', 'error']
    write_csv(out / 'pathway_summary.long.csv', summaries, summary_fields)
    target_ids = sorted({r['target'] for r in config['model_rules'].values()})
    write_csv(out / 'gene_target_evidence.csv', [dict(bin_id=r['bin_id'], **{
        g: r['targets'].get(g, 'not_detected') if r['status'] != 'failed' else 'bin_failed'
        for g in target_ids}) for r in results], ['bin_id'] + target_ids)
    for pathway in config['pathways']:
        pid = pathway['id']
        comps = list(dict.fromkeys(pathway['scored_components'] + pathway['route_components'] + pathway['extra_components']))
        for suffix, fieldname in [('components', 'value'), ('copy_counts', 'count'), ('component_status', 'state')]:
            rows = [dict(bin_id=r['bin_id'], **{c: r['evidence'][c][fieldname] if r['status'] != 'failed'
                        else 'bin_failed' if fieldname == 'state' else 'NA' for c in comps}) for r in results]
            write_csv(out / 'pathways' / f'{pid}.{suffix}.csv', rows, ['bin_id'] + comps)
        genes = list(dict.fromkeys(pathway['target_genes']))
        write_csv(out / 'pathways' / f'{pid}.gene_evidence.csv', [dict(bin_id=r['bin_id'], **{
            g: r['targets'].get(g, 'not_detected') if r['status'] != 'failed' else 'bin_failed'
            for g in genes}) for r in results], ['bin_id'] + genes)
    source = {r['model_accession']: r for r in mapping}
    rows = [hit_row(r, h, models, source) for r in results if r['status'] != 'failed' for h in r['hits']]
    write_csv(out / 'all_model_hits.csv', rows, HIT_FIELDS)
    write_csv(out / 'accepted_hits.csv', [r for r in rows if r['status'] == 'accepted'], HIT_FIELDS)
    write_csv(out / 'review_hits.csv', [r for r in rows if r['status'] != 'accepted' or r['target_call'] != 'specific'], HIT_FIELDS)
    components = []
    for result in results:
        if result['status'] == 'failed':
            continue
        for cid, e in result['evidence'].items():
            components.append({'bin_id': result['bin_id'], 'component': cid, 'state': e['state'],
                               'value': e['value'], 'copies': e['count'],
                               'proteins': ';'.join(sorted(e['proteins'])),
                               'candidate_proteins': ';'.join(sorted(e['candidates'])),
                               'assay_scope': config['components'][cid]['assay_scope'],
                               'notes': config['components'][cid]['notes']})
    write_csv(out / 'component_evidence.long.csv', components,
              ['bin_id', 'component', 'state', 'value', 'copies', 'proteins', 'candidate_proteins', 'assay_scope', 'notes'])
    write_csv(out / 'neighborhoods.csv', [n for r in results for n in r.get('neighborhoods', [])], NEIGHBOR_FIELDS)
    write_csv(out / 'operon_candidates.csv', [n for r in results for n in r.get('operons', [])], OPERON_FIELDS)
    write_csv(out / 'bin_qc.csv', results, ['bin_id', 'source', 'status', 'contigs', 'bases', 'orfs',
                                           'orfs_reused', 'searches_reused', 'folder', 'error'])
    definitions = []
    for pathway in config['pathways']:
        for c in dict.fromkeys(pathway['scored_components'] + pathway['route_components'] + pathway['extra_components']):
            definitions.append({'pathway': pathway['id'], 'component': c,
                                'in_heatmap_denominator': int(c in pathway['scored_components']),
                                'in_route_definition': int(c in pathway['route_components']),
                                **config['components'][c]})
    write_csv(out / 'component_definitions.csv', definitions, ['pathway', 'component', 'label',
              'assay_scope', 'in_heatmap_denominator', 'in_route_definition', 'notes'])


def discover_bins(folder, extension=None, recursive=False, exclude=None):
    folder = Path(folder).resolve()
    if not folder.is_dir():
        raise ValueError('Bin directory does not exist: ' + str(folder))
    extensions = [extension.lstrip('.').lower()] if extension else ['fa', 'fna', 'fasta', 'fas', 'fsa']
    result = []
    for path in sorted(folder.rglob('*') if recursive else folder.iterdir()):
        if exclude is not None and (path.resolve() == exclude or exclude in path.resolve().parents):
            continue
        if not path.is_file():
            continue
        filename = path.name[:-3] if path.name.lower().endswith('.gz') else path.name
        match = next((ext for ext in extensions if filename.lower().endswith('.' + ext)), None)
        if not match:
            continue
        relative = path.relative_to(folder)
        bid = str(relative.parent / filename[:-(len(match) + 1)])
        result.append((bid, path))
    if not result:
        raise ValueError('No matching nucleotide FASTA bins found')
    ids = [item[0] for item in result]
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate bin IDs after removing FASTA suffixes; rename these bins')
    return result


def make_parser():
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-bin_dir', '--bin-dir', required=True, help='Directory of nucleotide bins')
    parser.add_argument('-bin_ext', '--bin-ext', help='Optional single FASTA suffix without period; also accepts its .gz form')
    parser.add_argument('-out', '--out', default='pathwaygenie_out', help='Output directory')
    parser.add_argument('--hmm-dir', type=Path, default=base / 'Pathway_Gene_HMMs', help='Extracted HMM panel directory')
    parser.add_argument('--config', type=Path, default=base / 'pathwaygenie_config.json', help='Editable pathway/activity/operon configuration')
    parser.add_argument('-t', '--threads', type=int, default=8, help='Total HMMER worker CPU budget')
    parser.add_argument('-j', '--jobs', type=int, default=1, help='Number of bins processed concurrently; divides --threads')
    parser.add_argument('--prodigal-mode', choices=['meta', 'single'], default='meta')
    parser.add_argument('--meta', action='store_const', const='meta', dest='prodigal_mode', help='GasGenie-style alias for --prodigal-mode meta')
    parser.add_argument('--prodigal', default='prodigal', help='Prodigal executable')
    parser.add_argument('--hmmsearch', default='hmmsearch', help='HMMER3 hmmsearch executable')
    parser.add_argument('--min-hmm-coverage', type=float, default=0.60, help='Minimum union HMM coverage; protected model-specific floors also apply')
    parser.add_argument('--ambiguity-margin', type=float, default=0.10, help='Normalized cutoff-margin tolerance for incompatible overlapping hits')
    parser.add_argument('--overlap-fraction', type=float, default=0.5, help='Minimum overlap relative to shorter alignment to compare competitors')
    parser.add_argument('-d', '--max-intervening-genes', type=int, default=5, help='Maximum unmarked ORFs BETWEEN successive marked genes')
    parser.add_argument('--max-neighborhood-bp', type=int, default=20000, help='Maximum total neighborhood span in bp')
    parser.add_argument('--operon-order', choices=['prefer', 'strict', 'ignore'], default='prefer')
    parser.add_argument('--operon-strand', choices=['prefer', 'strict', 'ignore'], default='prefer')
    parser.add_argument('--operon-max-missing', type=int, help='Override per-template number of missing slots allowed')
    parser.add_argument('--operon-min-fraction', type=float, help='Override per-template detected-slot fraction required')
    parser.add_argument('--skip-context', action='store_true', help='Skip supporting/ambiguous HMM searches')
    parser.add_argument('--recursive', action='store_true', help='Search nested bin directories')
    parser.add_argument('--force', action='store_true', help='Rerun Prodigal/HMMER; otherwise use verified content-based caches')
    parser.add_argument('--version', action='version', version=VERSION)
    return parser


def main(argv=None):
    parser = make_parser()
    args = parser.parse_args(argv)
    if args.threads < 1 or not 1 <= args.jobs <= args.threads:
        parser.error('Require 1 <= jobs <= threads')
    if not 0 < args.min_hmm_coverage <= 1 or not 0 < args.overlap_fraction <= 1:
        parser.error('Coverage and overlap must lie in (0,1]')
    if args.ambiguity_margin < 0 or args.max_intervening_genes < 0 or args.max_neighborhood_bp < 1:
        parser.error('Invalid distance or ambiguity threshold')
    if args.operon_max_missing is not None and args.operon_max_missing < 0:
        parser.error('Operon missing-slot allowance must be nonnegative')
    if args.operon_min_fraction is not None and not 0 <= args.operon_min_fraction <= 1:
        parser.error('Operon fraction must lie in [0,1]')
    args.cpus_per_job = max(1, args.threads // args.jobs)
    try:
        out = Path(args.out).resolve()
        if out == Path(args.bin_dir).resolve():
            raise ValueError('Output directory must differ from input bin directory')
        bins = discover_bins(args.bin_dir, args.bin_ext, args.recursive, exclude=out)
        config = json.loads(args.config.read_text())
        validate_config(config)
        models, mapping, libraries, manifest = load_panel(args.hmm_dir.resolve(), config, not args.skip_context)
        for key in ['prodigal', 'hmmsearch']:
            executable = shutil.which(getattr(args, key))
            if not executable:
                raise ValueError(f'{key} executable not found: {getattr(args, key)}')
            setattr(args, key, executable)
        versions = {'prodigal': external_version(args.prodigal, '-v'),
                    'hmmsearch': external_version(args.hmmsearch, '-h')}
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.SubprocessError) as error:
        parser.exit(2, 'ERROR: ' + str(error) + '\n')
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    lock = out / '.pathwaygenie.lock'
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError:
        parser.exit(2, 'ERROR: Output is locked: ' + str(lock) + '. If a prior process crashed, verify it has stopped before removing this lock.\n')
    os.write(descriptor, str(os.getpid()).encode())
    os.close(descriptor)
    try:
        if args.config.resolve() != (out / 'config_used.json').resolve():
            shutil.copyfile(args.config, out / 'config_used.json')
        run_info = {'program': 'PathwayGenie', 'version': VERSION,
                    'started_utc': datetime.now(timezone.utc).isoformat(), 'versions': versions,
                    'arguments': {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                    'config_sha256': sha256(args.config),
                    'metadata_sha256': sha256(args.hmm_dir / 'gene_model_map.tsv'),
                    'ncbi_release': manifest.get('ncbi_release'), 'libraries': libraries,
                    'bin_inputs': [{'bin_id': bid, 'path': str(path), 'sha256': sha256(path)} for bid, path in bins]}
        write_json(out / 'run_manifest.json', run_info)
        print(f'PathwayGenie {VERSION}: {len(bins)} bins; {args.jobs} jobs; {args.cpus_per_job} HMMER CPUs/job.', flush=True)
        results = []
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(run_bin, item, out, models, libraries, config, args, versions): item for item in bins}
            for future in as_completed(futures):
                bid, source = futures[future]
                try:
                    result = future.result()
                    print(f'[{len(results)+1}/{len(bins)}] {bid}: {result["orfs"]} ORFs; '
                          f'{sum(h.status == "accepted" for h in result["hits"])} accepted model matches; '
                          f'ORFs reused={result["orfs_reused"]}.', flush=True)
                except Exception as error:
                    result = {'bin_id': bid, 'source': str(source), 'status': 'failed', 'error': str(error)}
                    print(f'FAILED {bid}: {error}', file=sys.stderr, flush=True)
                results.append(result)
        write_outputs(out, results, config, models, mapping)
        run_info['finished_utc'] = datetime.now(timezone.utc).isoformat()
        run_info['failed_bins'] = [r['bin_id'] for r in results if r['status'] == 'failed']
        write_json(out / 'run_manifest.json', run_info)
        print('Results: ' + str(out / 'pathway_heatmap.csv'), flush=True)
        return 1 if run_info['failed_bins'] else 0
    finally:
        lock.unlink(missing_ok=True)


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
