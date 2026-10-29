#!/usr/bin/env python3
"""Optional, resumable Bakta annotation of one MetaHopper bin partition.

One rank is annotated, preventing duplicated calls on nested rank partitions.
Original contig headers and raw Bakta outputs are retained for later modules.
"""
from __future__ import annotations
import argparse
import contextlib
import csv
import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import time
from collections import Counter
from pathlib import Path

from metahopper_graph import fasta, table

SUMMARY_FIELDS = ['rank', 'bin', 'status', 'reason', 'length_bp', 'num_contigs', 'translation_table',
                  'translation_source', 'bacterial_fraction', 'cds', 'sorf', 'hypothetical_proteins',
                  'named_proteins', 'trna', 'rrna', 'elapsed_seconds', 'cached', 'output_directory']
FEATURE_FIELDS = ['rank', 'bin', 'contig', 'locus_tag', 'type', 'start', 'end', 'strand', 'gene', 'product', 'db_xrefs']
CODE4 = {'vidania', 'nasuia', 'zinderia', 'spiroplasma'}
MANAGED = {'--db', '-d', '--output', '-o', '--prefix', '-p', '--threads', '-t', '--locus-tag',
           '--translation-table', '--meta', '--complete', '--force', '-f', '--keep-contig-headers', '--compliant'}


def write_tsv(path, rows, fields):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, data):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False, default=str) + '\n')
    tmp.replace(path)


def resolve_config(path):
    path = Path(path).resolve()
    config = json.loads(path.read_text())
    for key in ('bakta_db', 'annotation_metadata'):
        if config.get(key):
            value = Path(os.path.expandvars(config[key])).expanduser()
            config[key] = str((path.parent / value).resolve())
    return config


def validate_config(config):
    db = Path(config.get('bakta_db', ''))
    if not config.get('bakta_db') or not db.is_dir():
        raise ValueError('Bakta annotation requires config.bakta_db pointing to a prepared Bakta database directory')
    command = config.get('tools', {}).get('bakta', ['bakta'])
    if not isinstance(command, list) or not command or not all(isinstance(x, str) for x in command):
        raise ValueError('tools.bakta must be a nonempty executable-prefix array')
    if not shutil.which(command[0]):
        raise ValueError(f'Missing Bakta executable/prefix: {command[0]}')
    extra = config.get('bakta_args', ['--skip-plot'])
    if not isinstance(extra, list) or not all(isinstance(x, str) for x in extra):
        raise ValueError('bakta_args must be a JSON array of strings')
    if any(x.split('=', 1)[0] in MANAGED or re.match(r'^-[doptf][^-].+', x) for x in extra):
        raise ValueError('bakta_args cannot override managed input, identity, completeness, database, output or translation settings')
    if config.get('annotation_metadata') and not Path(config['annotation_metadata']).is_file():
        raise ValueError('annotation_metadata file does not exist')
    return db, command, extra


def run_command(command, log, timeout_h):
    with open(log, 'a') as f:
        f.write('COMMAND: ' + json.dumps(list(map(str, command))) + '\n')
        f.flush()
        proc = subprocess.Popen(list(map(str, command)), stdout=f, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = proc.wait(timeout=timeout_h * 3600 if timeout_h else None)
        except BaseException:
            with contextlib.suppress(ProcessLookupError):
                os.killpg(proc.pid, signal.SIGTERM)
            try:
                proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                with contextlib.suppress(ProcessLookupError):
                    os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
            raise
    if code:
        raise RuntimeError(f'Bakta exited {code}; see {log}')


def summarize_json(path, output, rank, bin_name, valid_contigs):
    data = json.loads(path.read_text())
    features = data.get('features')
    if not isinstance(features, list):
        raise ValueError('Bakta JSON has no features array')
    counts = Counter()
    normalized = []
    for feature in features:
        kind = str(feature.get('type', '')).lower()
        counts[kind] += 1
        product = str(feature.get('product') or '')
        if kind in ('cds', 'sorf'):
            hypothetical = bool(feature.get('hypothetical')) or not product or product.lower().startswith('hypothetical protein')
            counts['hypothetical_proteins' if hypothetical else 'named_proteins'] += 1
        cid = feature.get('contig') or feature.get('sequence') or feature.get('sequence_id') or ''
        if isinstance(cid, dict):
            cid = cid.get('id', '')
        if not cid or cid not in valid_contigs:
            raise ValueError(f'Bakta changed a contig ID despite --keep-contig-headers: {cid}')
        refs = feature.get('db_xrefs') or feature.get('dbxref') or []
        if isinstance(refs, str):
            refs = [refs]
        normalized.append(dict(rank=rank, bin=bin_name, contig=cid,
                               locus_tag=feature.get('locus') or feature.get('locus_tag') or feature.get('id') or '',
                               type=kind, start=feature.get('start', ''), end=feature.get('stop', feature.get('end', '')),
                               strand=feature.get('strand', ''), gene=feature.get('gene', ''), product=product,
                               db_xrefs=';'.join(map(str, refs))))
    write_tsv(output, normalized, FEATURE_FIELDS)
    return {k: counts[k] for k in ('cds', 'sorf', 'hypothetical_proteins', 'named_proteins', 'trna', 'rrna')}, data.get('version', {})


def annotate_run(root, config, rank='genus', stage='final', threads=16, mode='meta', timeout_h=6,
                 force=False, translation_table=11, selected_bins=None, runner=None, log=print):
    root = Path(root).resolve()
    if rank in ('.', '..') or Path(rank).name != rank:
        raise ValueError('Annotation rank must be a single bin-folder name')
    base = root / 'final' if stage == 'final' else root
    bindir = base / 'bins' / rank
    if not bindir.is_dir():
        raise ValueError(f'No {stage} bins at rank {rank}: {bindir}')
    if threads < 1 or timeout_h < 0 or translation_table not in (4, 11, 25) or mode not in ('meta', 'single'):
        raise ValueError('Invalid annotation threads, timeout, translation table or mode')
    db, prefix, extra = validate_config(config)
    out = base / 'annotation' / rank
    out.mkdir(parents=True, exist_ok=True)
    # Never infer successful annotation from the presence of an output directory.
    try:
        version = subprocess.run(prefix + ['--version'], text=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, timeout=60, check=True).stdout.strip()
    except (subprocess.SubprocessError, OSError) as e:
        raise ValueError(f'Could not determine Bakta version: {e}') from e
    db_versions = {p.name: p.read_text() for p in (db / 'version.json', db / 'version.txt') if p.is_file()}
    db_identity = dict(path=str(db), version_files=db_versions, directory_mtime_ns=db.stat().st_mtime_ns)
    metadata = {}
    for row in table(Path(config['annotation_metadata'])) if config.get('annotation_metadata') else []:
        name = row.get('bin')
        if name and (not row.get('rank') or row['rank'] == rank):
            if name in metadata:
                raise ValueError(f'Duplicate annotation metadata for {rank}/{name}')
            metadata[name] = dict(row)
            if row.get('translation_table'):
                code = int(row['translation_table'])
                if code not in (4, 11, 25):
                    raise ValueError(f'Unsupported translation table for {name}: {code}')
                metadata[name]['translation_table'] = code
            else:
                metadata[name].pop('translation_table', None)
    taxonomy = {r.get('contig'): r for r in table(base / 'classification/contig_classification.tsv')}
    selected = set(selected_bins) if selected_bins else None
    all_fastas = sorted(bindir.glob('*.fasta'))
    if selected and selected - {p.stem for p in all_fastas}:
        raise ValueError('Requested annotation bins not found: ' + ', '.join(sorted(selected - {p.stem for p in all_fastas})))
    rows, warnings = [], []
    for path in all_fastas:
        name = path.stem
        if name in ('.', '..'):
            raise ValueError('Invalid bin filename stem')
        if selected and name not in selected:
            continue
        t0 = time.monotonic()
        seqs = dict(fasta(path))
        length = sum(map(len, seqs.values()))
        domains = {taxonomy[c].get('domain') for c in seqs if c in taxonomy and taxonomy[c].get('domain') not in ('', None, 'Unclassified')}
        bacterial_bp = sum(len(s) for c, s in seqs.items() if taxonomy.get(c, {}).get('domain') == 'Bacteria')
        hint_codes = {4 if CODE4.intersection(re.findall(r'[a-z]+', taxonomy.get(c, {}).get('genus', '').lower())) else 11
                      for c in seqs if taxonomy.get(c, {}).get('domain') == 'Bacteria'
                      and taxonomy.get(c, {}).get('genus') not in (None, '', 'Unclassified')}
        # Explicit metadata can establish a bacterial bin and override taxon-based code hints.
        meta = metadata.get(name, {})
        effective_domains = {meta['domain']} if meta.get('domain') else domains
        hint = bool(CODE4.intersection(re.findall(r'[a-z]+', name.lower())))
        code = meta.get('translation_table', 4 if hint or hint_codes == {4} else translation_table)
        code_source = 'metadata' if meta.get('translation_table') else 'taxon_hint' if hint or hint_codes == {4} else 'default'
        row = dict(rank=rank, bin=name, status='', reason='', length_bp=length, num_contigs=len(seqs),
                   translation_table=code, translation_source=code_source,
                   bacterial_fraction=round(bacterial_bp / max(1, length), 4), cached=False,
                   output_directory=(out / name / 'bakta').relative_to(root).as_posix())
        reason = None
        if name == 'Unclassified' or name.startswith('unclassified_'):
            reason = 'unclassified or domain-only catch-all bin'
        elif not seqs:
            reason = 'empty bin'
        elif effective_domains != {'Bacteria'}:
            reason = 'Bakta requires bacterial bins; observed domains: ' + ','.join(sorted(effective_domains or {'unknown'}))
        elif len(hint_codes) > 1 and not meta.get('translation_table'):
            reason = 'mixed translation-code hints; supply per-bin annotation_metadata'
        if reason:
            row.update(status='skipped', reason=reason, elapsed_seconds=round(time.monotonic() - t0, 2))
            skipped_dir = out / name
            skipped_dir.mkdir(exist_ok=True)
            atomic_json(skipped_dir / 'result.json', dict(status='skipped', summary=row))
            atomic_json(skipped_dir / 'modules.json', dict(schema=1, rank=rank, bin=name,
                        input_sha256=sha256(path), bakta=dict(status='skipped', reason=reason),
                        metabolic_pathways=dict(status='not_run'), biosynthetic_clusters=dict(status='not_run'),
                        virulence_screen=dict(status='not_run')))
            rows.append(row)
            continue
        digest = sha256(path)
        locus = 'MH' + hashlib.sha256((root.name + rank + name + digest).encode()).hexdigest()[:12].upper()
        signature = dict(schema=1, input_sha256=digest, rank=rank, bin=name, tool_version=version,
                         command_prefix=prefix, db=db_identity, args=extra, translation_table=code,
                         mode=mode, locus_tag=locus)
        bdir = out / name
        bdir.mkdir(exist_ok=True)
        record = bdir / 'result.json'
        final = bdir / 'bakta'
        required = [final / ('annotation' + suffix) for suffix in ('.json', '.gff3', '.gbff', '.faa', '.ffn', '.tsv')]
        cached = json.loads(record.read_text()) if record.is_file() else {}
        if not force and cached.get('status') == 'completed' and cached.get('signature') == signature and all(p.is_file() for p in required) and (bdir / 'features.tsv').is_file():
            row.update({k: cached['summary'].get(k, 0) for k in ('cds', 'sorf', 'hypothetical_proteins', 'named_proteins', 'trna', 'rrna')})
            row.update(status='completed', cached=True, elapsed_seconds=0.0)
            rows.append(row)
            log(f'Bakta {rank}/{name}: cached')
            continue
        work = bdir / 'bakta.partial'
        if work.exists():
            shutil.rmtree(work)
        # Keep a prior successful annotation until a new one is complete.
        cmd = ['--db', str(db), '--output', str(work), '--prefix', 'annotation', '--threads', str(threads),
               '--locus-tag', locus, '--translation-table', str(code), '--keep-contig-headers']
        if mode == 'meta':
            cmd.append('--meta')
        cmd += extra + [str(path)]
        logfile = bdir / 'bakta_runner.log'
        log(f'Bakta {rank}/{name}: {length:,} bp; translation table {code} ({code_source})')
        try:
            if runner is None:
                run_command(prefix + cmd, logfile, timeout_h)
            else:
                runner.run('bakta', cmd, 'bakta_' + rank + '_' + name, timeout=timeout_h * 3600 or None)
            missing = [p.name for p in required if not (work / p.name).is_file()]
            if missing:
                raise ValueError('Bakta exited without required outputs: ' + ', '.join(missing))
            stats, tool_details = summarize_json(work / 'annotation.json', bdir / 'features.partial.tsv', rank, name, set(seqs))
            if final.exists():
                shutil.rmtree(final)
            work.replace(final)
            (bdir / 'features.partial.tsv').replace(bdir / 'features.tsv')
            row.update(stats, status='completed', reason='')
            modules = dict(schema=1, rank=rank, bin=name, input_fasta=path.relative_to(root).as_posix(),
                           input_sha256=digest, coordinate_system='1-based inclusive',
                           bakta=dict(status='completed', version=version, details=tool_details,
                                      translation_table=code, translation_source=code_source, locus_tag=locus, mode=mode,
                                      directory=final.relative_to(root).as_posix(), features=(bdir / 'features.tsv').relative_to(root).as_posix()),
                           metabolic_pathways=dict(status='not_run'), biosynthetic_clusters=dict(status='not_run'),
                           virulence_screen=dict(status='not_run'))
            atomic_json(bdir / 'modules.json', modules)
        except (RuntimeError, ValueError, OSError, subprocess.TimeoutExpired) as e:
            row.update(status='failed', reason=str(e))
            warnings.append(f'Bakta {rank}/{name}: {e}')
            atomic_json(bdir / 'modules.json', dict(schema=1, rank=rank, bin=name, input_sha256=digest,
                        bakta=dict(status='failed', reason=str(e), previous_output_preserved=final.exists()),
                        metabolic_pathways=dict(status='not_run'), biosynthetic_clusters=dict(status='not_run'),
                        virulence_screen=dict(status='not_run')))
        row['elapsed_seconds'] = round(time.monotonic() - t0, 2)
        atomic_json(record, dict(status=row['status'], signature=signature, summary=row))
        rows.append(row)
        write_tsv(out / 'summary.tsv', rows, SUMMARY_FIELDS)
    if selected:
        # Preserve valid annotations of other bins when only a subset is requested.
        # Changed inputs are labelled stale instead of reusing old feature counts.
        for path in all_fastas:
            if path.stem in selected:
                continue
            record = out / path.stem / 'result.json'
            if not record.is_file():
                continue
            previous = json.loads(record.read_text())
            old = dict(previous.get('summary', {}))
            if previous.get('status') == 'completed':
                if previous.get('signature', {}).get('input_sha256') == sha256(path):
                    old.update(reason='Existing annotation retained; not selected in this invocation', cached=True)
                else:
                    old = dict(rank=rank, bin=path.stem, status='stale', reason='Bin changed since its annotation; rerun this bin')
                rows.append(old)
        rows.sort(key=lambda r: r['bin'])
    write_tsv(out / 'summary.tsv', rows, SUMMARY_FIELDS)
    atomic_json(out / 'annotation_manifest.json', dict(schema=1, module='bakta', stage=stage, rank=rank,
                tool_version=version, database=db_identity, bins=rows, warnings=warnings,
                specialized_modules=dict(metabolic_pathways='not_run', biosynthetic_clusters='not_run', virulence_screen='not_run')))
    return dict(bins=rows, warnings=warnings, summary=str(out / 'summary.tsv'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-i', '--input', type=Path, required=True)
    p.add_argument('-c', '--config', type=Path, required=True)
    p.add_argument('--rank', default='genus')
    p.add_argument('--bin-set', choices=['auto', 'final', 'preliminary'], default='auto')
    p.add_argument('-t', '--threads', type=int, default=16)
    p.add_argument('--mode', choices=['meta', 'single'], default='meta')
    p.add_argument('--translation-table', type=int, choices=[4, 11, 25], default=11)
    p.add_argument('--timeout-h', type=float, default=6)
    p.add_argument('--bins', help='Comma-separated bin filename stems to annotate')
    p.add_argument('--force', action='store_true', help='Reannotate, ignoring successful cache entries')
    a = p.parse_args()
    stage = 'final' if a.bin_set == 'auto' and (a.input / 'final/bins').is_dir() else a.bin_set
    if stage == 'auto':
        stage = 'preliminary'
    try:
        data = annotate_run(a.input, resolve_config(a.config), a.rank, stage, a.threads, a.mode, a.timeout_h,
                            a.force, a.translation_table, a.bins.split(',') if a.bins else None)
    except (ValueError, OSError) as e:
        p.exit(1, f'ERROR: {e}\n')
    print('Annotation summary:', data['summary'])
    # Refresh existing reports without importing main-pipeline dependencies.
    import metahopper_report as report
    opts = argparse.Namespace(no_fasta_stats=False, min_length=1000, max_points=25000,
                             include_uncoordinated=False, title=None)
    try:
        (a.input / 'metahopper_report.html').write_text(report.build_report(report.RunLayout(a.input, stage), opts))
    except (OSError, ValueError, RuntimeError) as e:
        print('WARNING: annotation saved but report refresh failed:', e)
    if any(row['status'] == 'failed' for row in data['bins']):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
