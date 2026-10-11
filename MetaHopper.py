#!/usr/bin/env python3
"""MetaHopper Lean 2.0.0: taxonomy-first rank partitions, BinaRena tables and the 7.4 report.

Workflow
  1  phyloFlash SSU survey                     (context only)
  2  fastp read QC
  3  MEGAHIT (or a supplied assembly) + contig length filter
  4  Prodigal ORFs + Tiara domain prefilter (eukarya/organelle contigs are not searched)
  5  DIAMOND vs NR on the first/spaced N ORFs per contig
  6  7.4 bitscore-weighted per-rank vote; Metazoa/Viridiplantae exclusion; fallback names
  7  one stringent BWA pass: coverage/breadth for every contig + competitive read pools
  8  rank partitions (phylum..species) -> bins/<rank>/*.fasta, summaries, BinaRena table
  9  guarded per-bin Unicycler reassembly at one rank, with the 7.4 retention gate
 10  consolidated final assembly, inherited taxonomy, final bins/BinaRena, optional CheckM
 11  metahopper_report.py (7.4 report, unchanged)

No coverage/GC/k-mer splitting or recruitment is done: composition is exported to BinaRena
for manual curation.
"""
from __future__ import annotations
import argparse, contextlib, csv, json, os, shutil, signal, subprocess, sys, time, zlib
from collections import Counter, defaultdict
from pathlib import Path

for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(key, '1')
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from core import Taxonomy, fasta, tsv, phyloflash_abundance, Pools, stream_mapping, downsample_pool, opn  # noqa: E402
import taxbin as tb  # noqa: E402

VERSION = '2.0.0'
PROKARYOTIC_TIARA = {'bacteria', 'archaea', 'prokarya'}
DEFAULT_DIAMOND_ARGS = ['--evalue', '1e-5', '--max-target-seqs', '25']


def note(message):
    print(time.strftime('[%H:%M:%S] ') + message, flush=True)


# ---------------------------------------------------------------------------------------
# Command runner (process-group kill on timeout / Ctrl-C)
# ---------------------------------------------------------------------------------------
class Runner:
    def __init__(self, config, out, memory):
        self.config, self.logdir = config, out / 'logs'
        self.logdir.mkdir(parents=True, exist_ok=True)
        self.env = dict(os.environ)
        self.env['JAVA_TOOL_OPTIONS'] = f'-Xms256m -Xmx{max(1, memory // 4)}g'
        self.commands = self.logdir / 'commands.jsonl'

    def command(self, tool, args):
        prefix = self.config.get('tools', {}).get(tool, [tool])
        if not isinstance(prefix, list) or not prefix: raise ValueError(f'tools.{tool} must be a nonempty JSON string array')
        return prefix + list(map(str, args))

    def record(self, cmd):
        with open(self.commands, 'a') as f:
            f.write(json.dumps(dict(time=time.strftime('%FT%T%z'), command=cmd)) + '\n')

    def run(self, tool, args, label, cwd=None, timeout=None, stdout=None):
        cmd = self.command(tool, args); self.record(cmd)
        log = self.logdir / (label + '.log')
        with open(log, 'wb') as f, contextlib.ExitStack() as stack:
            out = stack.enter_context(open(stdout, 'wb')) if stdout else f
            proc = subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT if not stdout else f,
                                    cwd=cwd, env=self.env, start_new_session=True)
            try:
                code = proc.wait(timeout=timeout)
            except BaseException as e:
                with contextlib.suppress(ProcessLookupError): os.killpg(proc.pid, signal.SIGTERM)
                try: proc.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    with contextlib.suppress(ProcessLookupError): os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
                if isinstance(e, subprocess.TimeoutExpired):
                    raise RuntimeError(f'{tool} exceeded {timeout / 3600:.1f} h; killed. See {log}') from None
                raise
        if code: raise RuntimeError(f'{tool} exited {code}; see {log}')


class Checkpoints:
    """work/checkpoints/<stage>.json markers; a stage is reused only on --resume."""

    def __init__(self, out, enabled):
        self.dir = out / 'work' / 'checkpoints'; self.dir.mkdir(parents=True, exist_ok=True)
        self.enabled = enabled

    def get(self, stage):
        p = self.dir / f'{stage}.json'
        if self.enabled and p.is_file():
            data = json.loads(p.read_text())
            if all(Path(x).exists() for x in data.get('files', [])): return data
        return None

    def done(self, stage, files=(), **info):
        (self.dir / f'{stage}.json').write_text(json.dumps(dict(info, files=[str(f) for f in files],
                                                                 time=time.strftime('%FT%T%z')), indent=2))


# ---------------------------------------------------------------------------------------
# FASTA helpers
# ---------------------------------------------------------------------------------------
def write_fasta(path, records, wrap=70):
    """7.4 write_fasta: plain, 70-column."""
    with open(path, 'w') as fh:
        for sid, seq in records:
            fh.write(f'>{sid}\n')
            for i in range(0, len(seq), wrap): fh.write(seq[i:i + wrap] + '\n')


def length_filter(source, dest, min_len):
    """Keep contigs >= min_len; ID = first header token (MEGAHIT names kept)."""
    stats = Counter(); seen = set()
    with open(dest, 'w') as f:
        for name, seq in fasta(source):
            sid = name.split()[0]
            if sid in seen: raise ValueError(f'Duplicate contig ID in assembly: {sid}')
            seen.add(sid)
            if len(seq) < min_len:
                stats['removed'] += 1; stats['removed_bp'] += len(seq); continue
            stats['kept'] += 1; stats['kept_bp'] += len(seq)
            f.write(f'>{sid}\n{seq}\n')
    return dict(stats, min_len=min_len)


def load_sequences(path):
    return {name.split()[0]: seq.upper() for name, seq in fasta(path)}


# ---------------------------------------------------------------------------------------
# ORFs, Tiara prefilter, DIAMOND query
# ---------------------------------------------------------------------------------------
def orf_table(faa):
    coding, orfs = Counter(), Counter()
    with open(faa) as f:
        for line in f:
            if line.startswith('>'):
                p = line[1:].split(' # ')
                c = tb.contig_of(p[0])
                coding[c] += abs(int(p[2]) - int(p[1])) + 1; orfs[c] += 1
    return coding, orfs


def tiara_calls(path):
    calls = {}
    with open(path) as f:
        next(f, None)
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 3: calls[p[0].split()[0]] = (p[1], p[2])
    return calls


def prefilter_decisions(seqs, tiara, coding, orfs, rescue, rescue_eukarya, audit_fraction):
    """Tiara bacteria/archaea/prokarya -> search. Tiara eukarya -> removed (unless
    --rescue-eukarya and coding density > --rescue-coding-density). Unknown/organelle ->
    searched when coding density > threshold. A deterministic audit sample of removed
    contigs is searched in the same DIAMOND run and classified normally."""
    cut = int(audit_fraction * 2 ** 32); rows = {}
    for c, seq in seqs.items():
        n = len(seq); c1, c2 = tiara.get(c, ('unknown', 'n/a'))
        dens = min(1.0, coding[c] / n) if n else 0.0
        if c1 in PROKARYOTIC_TIARA: d = 'search'
        elif c1 == 'eukarya' and not rescue_eukarya: d = 'skip_eukarya'
        elif dens > rescue: d = 'search_rescued_coding_density'
        elif c1 == 'organelle': d = 'skip_organelle'
        elif c1 == 'eukarya': d = 'skip_eukarya'
        else: d = 'skip_unknown'
        audit = d.startswith('skip') and orfs[c] > 0 and zlib.crc32(c.encode(), 17) < cut
        rows[c] = dict(contig=c, length=n, tiara=c1, tiara_stage2=c2, orfs=orfs[c],
                       coding_density=round(dens, 4), decision=d, audited=audit)
    return rows


def select_orfs(n, cap, mode):
    if cap is None or n <= cap: return None   # all
    if mode == 'first': return set(range(1, cap + 1))
    if cap == 1: return {1 + (n - 1) // 2}
    return {1 + round(i * (n - 1) / (cap - 1)) for i in range(cap)}


def write_query(faa, dest, contigs, orf_counts, cap, mode):
    total = sent = 0; chosen = {}
    with open(faa) as f, open(dest, 'w') as o:
        write = False
        for line in f:
            if line.startswith('>'):
                total += 1
                orf = line[1:].split(' # ', 1)[0]; c, idx = orf.rsplit('_', 1)
                if contigs is not None and c not in contigs: write = False
                else:
                    if c not in chosen: chosen[c] = select_orfs(orf_counts[c], cap, mode)
                    sel = chosen[c]; write = sel is None or int(idx) in sel
                sent += write
            if write: o.write(line)
    return total, sent


def orf_limit(value):
    if str(value).lower() == 'all': return 'all'
    try: n = int(value)
    except ValueError: raise argparse.ArgumentTypeError("must be a positive integer or 'all'")
    if n < 1: raise argparse.ArgumentTypeError("must be a positive integer or 'all'")
    return n


# ---------------------------------------------------------------------------------------
# Tables (7.4 formats)
# ---------------------------------------------------------------------------------------
CLASS_HEADER = (['contig', 'length_bp', 'GC_percent', 'coding_density', 'n_cds', 'mean_depth', 'covered_fraction',
                 'triage_call', 'excluded_gene_density', 'refinement_status', 'n_orfs', 'n_orfs_with_hits']
                + [x for r in tb.ORDER for x in (r, f'{r}_support')]
                + ['excluded_animal_plant', 'tiara_class', 'prefilter_decision', 'diamond_searched'])


def fmt(v, spec, na='NA'):
    return na if v is None else format(v, spec)


def write_classification(path, ids, seqs, calls, cov, breadth, coding, orfs, hits, pre, excluded_kingdom, prefilter_excluded):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n'); w.writerow(CLASS_HEADER)
        for c in ids:
            s = seqs[c]; n = len(s); p = pre.get(c, {})
            res = calls.get(c, {}).get('calls') or tb.empty_calls()
            gc = 100.0 * (s.count('G') + s.count('C')) / n if n else 0.0
            row = [c, n, f'{gc:.3f}', fmt(min(1.0, coding[c] / n) if n and c in coding else None, '.5f'),
                   orfs.get(c, 'NA'), fmt(cov.get(c), '.5f'), fmt(breadth.get(c), '.5f'),
                   p.get('decision', 'not_evaluated'), c in prefilter_excluded, 'not_evaluated',
                   orfs.get(c, 0), hits.get(c, 0)]
            for r in tb.ORDER: row += list(res[r])
            row += [c in excluded_kingdom, p.get('tiara', 'NA'), p.get('decision', 'NA'),
                    bool(p) and (p.get('decision', '').startswith('search') or p.get('audited', False))]
            w.writerow(row)


def write_binarena(path, ids, seqs, calls, labels_by_rank, bin_rank, part, cov, breadth, coding, orfs,
                   pre, excluded, columns, bin_ranks):
    names = list(columns)
    header = (['ID', 'length', 'GC', 'coverage', 'covered_fraction', 'coding_density', 'n_cds']
              + [f'taxon_{r}' for r in bin_ranks]
              + ['bin', 'taxon_support', 'triage', 'refinement_status', 'gc_outlier', 'coverage_outlier',
                 'joint_outlier', 'retained'] + names)
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n'); w.writerow(header)
        for c in ids:
            s = seqs[c]; n = len(s)
            res = calls.get(c, {}).get('calls') or tb.empty_calls()
            d = cov.get(c); b = breadth.get(c)
            row = [c, n, f'{100.0 * tb.gc_degenerate(s) / n:.2f}' if n else '0.00',
                   f'{d:.4f}' if d is not None else '', f'{b:.4f}' if b is not None else '',
                   f'{min(1.0, coding[c] / n):.4f}' if n and c in coding else 'NA', orfs.get(c, 'NA')]
            row += [labels_by_rank[r].get(c, tb.UNCLASSIFIED) for r in bin_ranks]
            row += [tb.sanitize(part.get(c, 'excluded') if c not in excluded else 'excluded'),
                    f'{res[bin_rank][1]:.4f}', pre.get(c, {}).get('decision', 'NA'), 'not_evaluated',
                    'false', 'false', 'false', 'false' if c in excluded else 'true']
            for nm in names:
                v = columns[nm].get(c); row.append(f'{v:.6f}' if v is not None else '')
            w.writerow(row)


def write_rank_bins(base, bin_ranks, partitions, seqs, unclassified_shared):
    """bins/<rank>/<bin>.fasta + bin_membership.tsv + summary.tsv; Unclassified shared via symlink."""
    gc = {c: s.count('G') + s.count('C') for c, s in seqs.items()}
    lengths = {c: len(s) for c, s in seqs.items()}
    unclassified_sets = {r: frozenset(c for c, l in partitions[r].items() if l == tb.UNCLASSIFIED) for r in bin_ranks}
    common = None
    if len(set(unclassified_sets.values())) == 1 and next(iter(unclassified_sets.values())):
        common = next(iter(unclassified_sets.values()))
        write_fasta(unclassified_shared, ((c, seqs[c]) for c in seqs if c in common))
    summaries = {}
    for r in bin_ranks:
        d = base / r
        if d.exists(): shutil.rmtree(d)
        d.mkdir(parents=True)
        groups = defaultdict(list)
        for c, l in partitions[r].items(): groups[tb.sanitize(l)].append(c)
        rows = []
        for name in sorted(groups):
            members = groups[name]; path = d / f'{name}.fasta'
            if name == tb.UNCLASSIFIED and common is not None:
                path.symlink_to(os.path.relpath(unclassified_shared, d))
            else:
                write_fasta(path, ((c, seqs[c]) for c in members))
            st = tb.bin_stats(members, lengths, gc)
            rows.append(dict(bin=name, **st, completeness_percent='NA', contamination_percent='NA',
                             strain_heterogeneity_percent='NA', marker_lineage='NA'))
        with open(d / 'bin_membership.tsv', 'w') as fh:
            fh.write('contig\tbin\n')
            for c, l in partitions[r].items(): fh.write(f'{c}\t{tb.sanitize(l)}\n')
        tsv(d / 'summary.tsv', rows, tb.SUMMARY_FIELDS)
        summaries[r] = {row['bin']: row for row in rows}
        note(f'     {r:<8} {len(groups):>4} bins ({sum(1 for g in groups if g != tb.UNCLASSIFIED)} named), '
             f'Unclassified {sum(lengths[c] for c in groups.get(tb.UNCLASSIFIED, [])) / 1e6:,.2f} Mbp')
    return summaries


def paf_retention(paf, original_lengths, min_identity=0.95):
    """7.4 paf_retention: union of aligned original spans at >= min_identity."""
    intervals = defaultdict(list)
    with open(paf) as fh:
        for line in fh:
            f = line.rstrip().split('\t')
            if len(f) < 12 or f[0] not in original_lengths: continue
            if int(f[10]) <= 0 or int(f[9]) / int(f[10]) < min_identity: continue
            lo, hi = max(0, int(f[2])), min(original_lengths[f[0]], int(f[3]))
            if hi > lo: intervals[f[0]].append((lo, hi))
    aligned = 0
    for spans in intervals.values():
        end = 0
        for lo, hi in sorted(spans):
            aligned += max(0, hi - max(lo, end)); end = max(end, hi)
    return aligned / max(1, sum(original_lengths.values()))


def candidate_decision(before, after, retained, min_retained=0.95, max_inflation=1.25):
    """7.4 candidate_decision without marker comparison."""
    ratio = after['total_length_bp'] / max(1, before['total_length_bp'])
    if retained < min_retained or ratio < min_retained: return 'rejected_sequence_loss'
    if ratio > max_inflation: return 'rejected_length_inflation'
    if after['num_contigs'] <= before['num_contigs'] * 0.9 or after['N50'] >= max(1, before['N50']) * 1.2:
        return 'accepted'
    return 'rejected_no_gain'


# ---------------------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------------------
def validate(args, config):
    for p in [args.read1, args.read2, args.assembly]:
        if p and not Path(p).is_file(): raise ValueError(f'Input file does not exist: {p}')
    for key in ('diamond_db', 'taxdump'):
        if not config.get(key): raise ValueError(f'Missing config key: {key}')
        if not Path(config[key]).exists(): raise ValueError(f'{key} does not exist: {config[key]}')
    for name in ('nodes.dmp', 'names.dmp'):
        if not (Path(config['taxdump']) / name).is_file(): raise ValueError(f'Missing taxdump/{name}')
    bad = [r for r in args.ranks if r not in tb.ORDER[2:]]
    if bad: raise ValueError(f'--ranks must be drawn from {",".join(tb.ORDER[2:])}; got {bad}')
    if not args.skip_reassembly and args.reassemble_rank not in args.ranks:
        raise ValueError('--reassemble-rank must be one of --ranks')
    required = ['fastp', 'prodigal', 'diamond', 'bwa']
    if not args.skip_phyloflash:
        required.append('phyloFlash.pl')
        if not config.get('phyloflash_db') or not Path(config['phyloflash_db']).exists():
            raise ValueError('phyloflash_db missing from config (or use --skip-phyloflash)')
    if args.prefilter == 'tiara': required.append('tiara')
    if not args.assembly: required.append('megahit')
    if not args.skip_reassembly: required += ['unicycler', 'minimap2']
    if args.checkm: required.append('checkm')
    for name in required:
        command = config.get('tools', {}).get(name, [name])
        if not isinstance(command, list) or not command or not all(isinstance(x, str) for x in command):
            raise ValueError(f'tools.{name} must be an executable/prefix array')
        if not shutil.which(command[0]): raise ValueError(f'Missing executable for {name}: {command[0]}')
    if args.inspect != 'off' and not (HERE / 'metahopper_inspect.py').is_file():
        raise ValueError('metahopper_inspect.py must sit next to metahopper.py (or use --inspect off)')
    if not (HERE / 'metahopper_report.py').is_file() and not args.skip_report:
        raise ValueError('metahopper_report.py must sit next to metahopper.py (or use --skip-report)')
    if args.threads < 1 or args.memory < 4: raise ValueError('Use at least 1 thread and 4 GB memory.')
    import pysam, sklearn, numpy  # noqa: F401  (fail early if missing)


# ---------------------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------------------
def stage_phyloflash(runner, args, config, out, ck):
    if ck.get('phyloflash'): return
    if args.skip_phyloflash: note('1    phyloFlash skipped'); return
    note('1    phyloFlash SSU survey')
    directory = out / 'work' / 'phyloflash'
    if directory.exists(): shutil.rmtree(directory)
    directory.mkdir(parents=True)
    with opn(args.read1) as f:
        f.readline(); length = len(f.readline().strip())
    cmd = ['-lib', 'sample', '-read1', args.read1, '-dbhome', config['phyloflash_db'], '-CPUs', args.threads,
           '-readlength', min(max(length, 50), 500), '-skip_spades', '-noemirge']
    if args.read2: cmd += ['-read2', args.read2]
    runner.run('phyloFlash.pl', cmd, 'phyloflash', cwd=directory)
    source = directory / 'sample.phyloFlash.NTUfull_abundance.csv'
    if not source.exists(): raise RuntimeError('phyloFlash produced no NTU abundance CSV; see logs/phyloflash.log')
    abundance = phyloflash_abundance(source)
    (out / 'phyloflash').mkdir(exist_ok=True)
    shutil.copy(source, out / 'phyloflash' / 'NTUfull_abundance.csv')
    tsv(out / 'phyloflash' / 'phyloflash_all.tsv', abundance,
        ['lineage', 'ssu_count', 'retained', 'reason', 'rank_retained', 'fraction_all_ssu', 'fraction_retained_ssu'])
    shutil.rmtree(directory)
    ck.done('phyloflash', [out / 'phyloflash' / 'phyloflash_all.tsv'])


def stage_qc(runner, args, out, ck):
    qdir = out / 'qc'; qc1, qc2 = qdir / 'qc.R1.fq.gz', qdir / 'qc.R2.fq.gz'
    reads = [qc1] + ([qc2] if args.read2 else [])
    if ck.get('qc'): return reads
    note('2    Read QC (fastp)')
    qdir.mkdir(exist_ok=True)
    cmd = ['-i', args.read1, '-o', qc1, '--thread', min(args.threads, 16), '--qualified_quality_phred', '20',
           '--unqualified_percent_limit', '40', '--length_required', '75', '--cut_right', '--cut_right_mean_quality', '20',
           '--json', qdir / 'fastp.json', '--html', qdir / 'fastp.html']
    if args.read2: cmd += ['-I', args.read2, '-O', qc2, '--detect_adapter_for_pe']
    runner.run('fastp', cmd, 'fastp')
    ck.done('qc', reads)
    return reads


def stage_assembly(runner, args, out, reads, ck):
    if args.assembly: path = out / 'assembly' / 'preliminary_contigs.fasta'
    else: path = out / 'megahit' / 'final.contigs.fa'
    if ck.get('assembly'): return path
    path.parent.mkdir(parents=True, exist_ok=True)
    if args.assembly:
        note('3    Preliminary assembly: supplied FASTA'); source = Path(args.assembly)
    else:
        note('3    Preliminary assembly (MEGAHIT)')
        mdir = out / 'work' / 'megahit'
        if mdir.exists(): shutil.rmtree(mdir)
        cmd = ['-1', reads[0], '-2', reads[1]] if len(reads) == 2 else ['-r', reads[0]]
        runner.run('megahit', cmd + ['-o', mdir, '-t', args.threads, '-m', args.memory * 1024 ** 3,
                                     '--min-contig-len', args.min_contig_len], 'megahit')
        source = mdir / 'final.contigs.fa'
    tmp = path.with_suffix('.tmp')
    stats = length_filter(source, tmp, args.min_contig_len)
    if not stats.get('kept'): raise RuntimeError(f'No preliminary contigs >= {args.min_contig_len} bp.')
    tmp.replace(path)
    if not args.assembly: shutil.rmtree(out / 'work' / 'megahit')
    (path.parent / 'length_filter.json').write_text(json.dumps(stats, indent=2) + '\n')
    note(f'     length filter (>= {args.min_contig_len:,} bp): kept {stats.get("kept", 0):,} contigs '
         f'({stats.get("kept_bp", 0) / 1e6:,.1f} Mbp), removed {stats.get("removed", 0):,} ({stats.get("removed_bp", 0) / 1e6:,.1f} Mbp)')
    ck.done('assembly', [path])
    return path


def run_diamond(runner, args, config, query, hits, temp, label):
    temp.mkdir(parents=True, exist_ok=True)
    tmp_hits = hits.with_suffix('.partial')
    extra = config.get('diamond_args', DEFAULT_DIAMOND_ARGS)
    runner.run('diamond', ['blastp', '--db', config['diamond_db'], '--query', query, '--out', tmp_hits,
                           '--threads', args.threads, *extra, '--tmpdir', temp,
                           '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore', 'staxids'], label)
    tmp_hits.replace(hits)
    shutil.rmtree(temp, ignore_errors=True)


def classify_assembly(runner, args, config, out, assembly, seqs, tax, ck, warnings, label='preliminary'):
    """Prodigal -> Tiara prefilter -> DIAMOND -> per-rank vote. Returns everything downstream needs."""
    cdir = out / ('classification' if label == 'preliminary' else 'final/classification')
    cdir.mkdir(parents=True, exist_ok=True)
    faa = cdir / 'orfs.faa'
    if not ck.get(f'{label}_prodigal'):
        note(f'4    {label}: Prodigal ORFs')
        runner.run('prodigal', ['-i', assembly, '-a', faa, '-o', cdir / 'orfs.gff', '-f', 'gff', '-p', 'meta', '-q'],
                   f'{label}_prodigal')
        ck.done(f'{label}_prodigal', [faa])
    coding, orfs = orf_table(faa)
    pre, selected = {}, None
    if args.prefilter == 'tiara' and label == 'preliminary':
        tpath = cdir / 'tiara.tsv'
        if not ck.get('tiara'):
            note('4    Tiara domain prefilter')
            runner.run('tiara', ['-i', assembly, '-o', tpath, '-m', args.min_contig_len, '-t', args.threads, '--probabilities'],
                       'tiara')
            ck.done('tiara', [tpath])
        pre = prefilter_decisions(seqs, tiara_calls(tpath), coding, orfs, args.rescue_coding_density,
                                  args.rescue_eukarya, args.prefilter_audit)
        selected = {c for c, r in pre.items() if r['decision'].startswith('search') or r['audited']}
        by = defaultdict(lambda: [0, 0])
        for r in pre.values(): by[r['decision']][0] += 1; by[r['decision']][1] += r['length']
        for k, (n, bp) in sorted(by.items()): note(f'       {k}: {n:,} contigs, {bp / 1e6:,.2f} Mbp')
        note(f'       audit sample: {sum(r["audited"] for r in pre.values()):,} removed contigs searched anyway')
    hits = (out / 'diamond' / 'hits.tsv') if label == 'preliminary' else (out / 'final' / 'diamond' / 'hits.tsv')
    hits.parent.mkdir(parents=True, exist_ok=True)
    if not ck.get(f'{label}_diamond'):
        query = cdir / 'diamond_query.faa'
        cap = None if args.max_orfs_per_contig == 'all' else int(args.max_orfs_per_contig)
        total, sent = write_query(faa, query, selected, orfs, cap, args.orf_selection)
        note(f'5    {label}: DIAMOND on {sent:,}/{total:,} proteins ({100 * sent / max(1, total):.1f}%), '
             f'ORF cap {args.max_orfs_per_contig} ({args.orf_selection})')
        if sent: run_diamond(runner, args, config, query, hits, out / 'work' / f'{label}_diamond_tmp', f'{label}_diamond')
        else: hits.write_text('')
        query.unlink(missing_ok=True)
        ck.done(f'{label}_diamond', [hits], proteins_total=total, proteins_sent=sent)
    note(f'6    {label}: per-rank bitscore vote')
    calls, stats = tb.classify_hits(hits, tax, tb.ORDER, args.bitscore_range, args.max_hits_per_orf,
                                    args.min_support, args.support_denominator)
    if stats.get('noncontiguous_orfs'):
        warnings.append(f'{label} DIAMOND output had {stats["noncontiguous_orfs"]} non-contiguous ORF groups.')
    if stats.get('hits') and not calls:
        raise RuntimeError('DIAMOND hits carried no usable taxids: the database needs taxonomy (makedb --taxonmap/--taxonnodes).')
    (cdir / 'vote_stats.json').write_text(json.dumps(stats, indent=2) + '\n')
    return calls, pre, coding, orfs


def main(args):
    config_path = Path(args.config).resolve(); config = json.loads(config_path.read_text())
    for key in ('diamond_db', 'taxdump', 'phyloflash_db', 'bakta_db'):
        if config.get(key):
            p = Path(os.path.expandvars(config[key])).expanduser()
            config[key] = str((config_path.parent / p).resolve())
    for key in ('read1', 'read2', 'assembly'):
        if getattr(args, key): setattr(args, key, str(Path(getattr(args, key)).resolve()))
    validate(args, config)
    out = Path(args.out).resolve()
    if not args.resume and out.exists() and any(out.iterdir()):
        raise ValueError('Output directory must be new or empty (use --resume RUN_DIR to continue a run).')
    out.mkdir(parents=True, exist_ok=True)
    runner = Runner(config, out, args.memory)
    ck = Checkpoints(out, bool(args.resume))
    manifest = dict(version=VERSION, arguments={k: v for k, v in vars(args).items()}, config=config,
                    start=time.strftime('%FT%T%z'), status='running')
    (out / 'run.json').write_text(json.dumps(manifest, indent=2, default=str))
    warnings = []
    t0 = time.time()
    try:
        if not ck.get('dbinfo'):
            runner.run('diamond', ['dbinfo', '--db', config['diamond_db']], 'database_preflight'); ck.done('dbinfo')
        stage_phyloflash(runner, args, config, out, ck)
        reads = stage_qc(runner, args, out, ck)
        assembly = stage_assembly(runner, args, out, reads, ck)
        seqs = load_sequences(assembly)
        ids = list(seqs); lengths = {c: len(s) for c, s in seqs.items()}
        tax = Taxonomy(config['taxdump'])
        calls, pre, coding, orfs = classify_assembly(runner, args, config, out, assembly, seqs, tax, ck, warnings)

        # Exclusions: Tiara-removed eukarya/organelles (unless DIAMOND on the audit sample
        # showed otherwise) and Metazoa/Viridiplantae by kingdom vote (7.4 rule).
        exclude_kingdoms = {k.strip() for k in args.exclude_kingdoms.split(',') if k.strip()}
        excluded_kingdom = {c for c, r in calls.items() if r['calls']['kingdom'][0] in exclude_kingdoms}
        prefilter_excluded = set()
        audit_rows = []
        for c, r in pre.items():
            if r['decision'] in ('skip_eukarya', 'skip_organelle'):
                if r['audited']:
                    dom = calls.get(c, {}).get('calls', tb.empty_calls())['domain'][0]
                    audit_rows.append(dict(r, audit_domain=dom,
                                           audit_genus=calls.get(c, {}).get('calls', tb.empty_calls())['genus'][0]))
                    if dom in ('Bacteria', 'Archaea'): continue   # DIAMOND evidence overrides Tiara
                prefilter_excluded.add(c)
            elif r['decision'] == 'skip_unknown' and r['audited']:
                dom = calls.get(c, {}).get('calls', tb.empty_calls())['domain'][0]
                audit_rows.append(dict(r, audit_domain=dom,
                                       audit_genus=calls.get(c, {}).get('calls', tb.empty_calls())['genus'][0]))
        excluded = excluded_kingdom | prefilter_excluded
        if pre:
            tsv(out / 'classification' / 'prefilter.tsv', list(pre.values()),
                ['contig', 'length', 'tiara', 'tiara_stage2', 'orfs', 'coding_density', 'decision', 'audited'])
            if audit_rows:
                tsv(out / 'classification' / 'prefilter_audit.tsv', audit_rows,
                    ['contig', 'length', 'tiara', 'coding_density', 'decision', 'audit_domain', 'audit_genus'])
            removed = [r for r in pre.values() if r['decision'] in ('skip_eukarya', 'skip_organelle')]
            aud = [r for r in audit_rows if r['decision'] in ('skip_eukarya', 'skip_organelle')]
            prok = [r for r in aud if r['audit_domain'] in ('Bacteria', 'Archaea')]
            if aud:
                frac = sum(r['length'] for r in prok) / max(1, sum(r['length'] for r in aud))
                est = frac * sum(r['length'] for r in removed)
                msg = (f'Tiara audit: {len(prok)}/{len(aud)} sampled removed contigs were prokaryotic by DIAMOND '
                       f'(~{est / 1e6:.2f} Mbp of removed sequence estimated prokaryotic; top genera: '
                       f'{", ".join(g for g, _ in Counter(r["audit_genus"] for r in prok).most_common(4))}).')
                note('     ' + msg)
                if prok: warnings.append(msg + ' Consider --rescue-eukarya or --prefilter none for this sample.')
        if excluded:
            write_fasta(out / 'classification' / 'excluded_animal_plant_contamination.fasta',
                        ((c, seqs[c]) for c in ids if c in excluded_kingdom))
            write_fasta(out / 'classification' / 'excluded_tiara_eukarya_organelle.fasta',
                        ((c, seqs[c]) for c in ids if c in prefilter_excluded))

        binnable = {c: (calls[c]['calls'] if c in calls else tb.empty_calls()) for c in ids if c not in excluded}
        partitions = {r: tb.partition(binnable, r, lengths, args.min_bin_length, args.min_bin_contigs) for r in args.ranks}
        labels = {r: {c: tb.fallback_label(v, r) for c, v in binnable.items()} for r in args.ranks}
        hits_n = {c: calls[c]['orfs_with_hits'] for c in calls}
        R = args.reassemble_rank
        part = partitions.get(R) or partitions[args.ranks[-1]]

        # Eligible bins for reassembly: named, prokaryotic majority, under the size cap.
        groups = defaultdict(list)
        for c, l in part.items(): groups[tb.sanitize(l)].append(c)
        eligible = {}
        if not args.skip_reassembly:
            cats = {x.strip() for x in args.reassemble_categories.split(',')}
            for name, members in groups.items():
                if name == tb.UNCLASSIFIED: continue
                if name.startswith('unclassified_'):   # domain-only catch-all: a mixture by construction
                    eligible[name] = dict(members=members, bp=sum(lengths[c] for c in members), domain='mixed',
                                          skip='domain-only catch-all bin (mixed lineages)'); continue
                bp = sum(lengths[c] for c in members); dom = tb.majority_domain(members, binnable, lengths)
                eligible[name] = dict(members=members, bp=bp, domain=dom,
                                      skip=(f'domain {dom} not in --reassemble-categories' if dom not in cats else
                                            f'{bp / 1e6:.1f} Mbp exceeds --max-bin-mbp {args.max_bin_mbp:g}'
                                            if bp > args.max_bin_mbp * 1e6 else None))
        contig_bins = {c: n for n, e in eligible.items() if not e['skip'] for c in e['members']}

        # 7: one stringent mapping pass for coverage, breadth and competitive pools.
        covfile = out / 'work' / 'coverage.tsv'; poolfile = out / 'work' / 'pool_counts.json'
        pool = Pools(out / 'work' / 'read_pools') if contig_bins else None
        cpt = ck.get('mapping')
        if cpt:
            cov, breadth = {}, {}
            with open(covfile) as f:
                next(f)
                for line in f:
                    c, d, b = line.rstrip('\n').split('\t'); cov[c] = float(d); breadth[c] = float(b)
            pc = json.loads(poolfile.read_text()) if poolfile.exists() else {}
            if pool: pool.counts.update(pc.get('counts', {})); pool.bases.update(pc.get('bases', {}))
            routing = pc.get('routing', {})
        else:
            note(f'7    Stringent BWA mapping of all QC reads: coverage + read pools for {len(set(contig_bins.values()))} bins')
            cov, routing, breadth = stream_mapping(runner, assembly, reads, args.threads, 'mapping', strict=True,
                                                   contig_bins=contig_bins, pools=pool, breadth=True)
            with open(covfile, 'w') as f:
                f.write('contig\tmean_depth\tcovered_fraction\n')
                for c in ids: f.write(f'{c}\t{cov.get(c, 0.0)}\t{breadth.get(c, 0.0)}\n')
            poolfile.write_text(json.dumps(dict(counts=dict(pool.counts) if pool else {},
                                                bases=dict(pool.bases) if pool else {}, routing=routing)))
            ck.done('mapping', [covfile])
        write_classification(out / 'classification' / 'contig_classification.tsv', ids, seqs, calls, cov, breadth,
                             coding, orfs, hits_n, pre, excluded_kingdom, prefilter_excluded)

        note('8    Taxonomic rank partitions (fallback names; small bins cascade upward)')
        write_rank_bins(out / 'bins', args.ranks, partitions, {c: seqs[c] for c in binnable},
                        out / 'classification' / 'unclassified_contigs.fasta')
        bdir = out / 'binarena'; bdir.mkdir(exist_ok=True)
        orient = {c: seqs[c] for c in ids if (args.binarena_include_excluded or c not in excluded)}
        columns = tb.ordinate(orient, args.binarena_kmers, args.binarena_methods, args.binarena_min_length,
                              args.binarena_max_contigs, args.binarena_force, args.binarena_seed, outdir=bdir, log=note)
        write_binarena(bdir / 'binarena_input.tsv', ids, seqs, calls, labels, R if R in args.ranks else args.ranks[-1],
                       part, cov, breadth, coding, orfs, pre, excluded, columns, args.ranks)

        # 9: guarded Unicycler reassembly at one rank.
        accepted = {}
        if not args.skip_reassembly:
            accepted = stage_reassembly(runner, args, out, R, eligible, groups, seqs, pool, cov, warnings)
        else:
            note('9    Reassembly skipped')

        # 10: consolidated final assembly.
        stage_final(runner, args, config, out, R, groups, seqs, binnable, calls, accepted, cov, breadth,
                    coding, orfs, pre, excluded, reads, warnings, ck)
        if args.inspect != 'off': stage_inspect(args, out, R, reads)
        manifest.update(status='completed')
    except BaseException as e:
        manifest.update(status='failed', error=str(e) or type(e).__name__)
        raise
    finally:
        manifest.update(end=time.strftime('%FT%T%z'), hours=round((time.time() - t0) / 3600, 3), warnings=warnings)
        (out / 'run.json').write_text(json.dumps(manifest, indent=2, default=str))
        if warnings: (out / 'logs' / 'warnings.txt').write_text('\n'.join(warnings) + '\n')
    if not args.skip_report: write_report(out, args)
    for w in warnings: note('WARNING: ' + w)
    note(f'Done in {(time.time() - t0) / 3600:.2f} h: {out}')


def stage_reassembly(runner, args, out, R, eligible, groups, seqs, pool, cov, warnings):
    note(f'9    Unicycler reassembly at {R} rank (one bin at a time)')
    rroot = out / 'reassembly' / R; rroot.mkdir(parents=True, exist_ok=True)
    started = time.time(); rows = {}; accepted = {}
    order = sorted(eligible, key=lambda n: eligible[n]['bp'])
    lengths = {c: len(s) for c, s in seqs.items()}
    gcn = {c: s.count('G') + s.count('C') for c, s in seqs.items()}
    for i, name in enumerate(order, 1):
        e = eligible[name]; wdir = rroot / name; done = wdir / 'result.json'
        if done.is_file():
            row = json.loads(done.read_text()); rows[name] = row
            if row['outcome'] == 'accepted': accepted[name] = wdir / 'assembly.fasta'
            continue
        wdir.mkdir(exist_ok=True)
        before = tb.bin_stats(e['members'], lengths, gcn)
        row = dict(bin=name, outcome='', elapsed_seconds=0.0, contigs_before=before['num_contigs'],
                   total_length_before_bp=before['total_length_bp'], N50_before=before['N50'],
                   contigs_after='', total_length_after_bp='', N50_after='', original_sequence_retained='',
                   recovered_fraction='', marker_check='not evaluated', extension_rounds=0,
                   seed_templates=(pool.counts.get(name, 0) // (2 if args.read2 else 1)) if pool else 0, accepted_templates='',
                   recruitment_seconds='', assembly_seconds='', polish_seconds='', polish_outcome='not_run', reason='')
        note(f'  {i}/{len(order)} {name}: {before["total_length_bp"] / 1e6:.2f} Mbp, {before["num_contigs"]:,} contigs')
        t = time.time(); udir = out / 'work' / f'unicycler_{name}'
        try:
            if e['skip']: raise ValueError('skipped: ' + e['skip'])
            if not pool or not pool.counts.get(name): raise ValueError('skipped: no unambiguous recruited read pairs')
            suffixes = ['R1', 'R2'] if args.read2 else ['SE']
            depth = pool.bases.get(name, 0) / max(1, before['total_length_bp'])
            templates = pool.counts.get(name, 0) // len(suffixes)
            marker = wdir / 'downsampled.json'
            if marker.is_file(): templates = json.loads(marker.read_text())['templates']
            elif args.target_coverage and depth > args.target_coverage:
                templates = downsample_pool(pool.directory, name, suffixes, args.target_coverage / depth)
                marker.write_text(json.dumps(dict(templates=templates, depth=depth)))
                note(f'     recruited depth ~{depth:,.0f}x -> downsampled to ~{args.target_coverage:g}x ({templates:,} templates)')
            row['accepted_templates'] = templates
            ureads = (['-1', pool.directory / f'{name}.R1.fq.gz', '-2', pool.directory / f'{name}.R2.fq.gz']
                      if args.read2 else ['-s', pool.directory / f'{name}.SE.fq.gz'])
            if udir.exists(): shutil.rmtree(udir)
            runner.run('unicycler', ureads + ['-o', udir, '-t', args.threads, '--mode', args.mode, '--keep', '0',
                                              '--spades_options', f'-m {args.memory}'], f'{name}_unicycler',
                       timeout=args.reassembly_timeout_h * 3600 or None)
            row['assembly_seconds'] = round(time.time() - t, 1)
            cand = [(f'{name}_reassembled_{k:05d}', s.upper()) for k, (n_, s) in
                    enumerate((x for x in fasta(udir / 'assembly.fasta') if len(x[1]) >= args.min_contig_len), 1)]
            if not cand: raise ValueError(f'no reassembled contigs >= {args.min_contig_len} bp')
            write_fasta(wdir / 'assembly.fasta', cand)
            if (udir / 'assembly.gfa').exists(): shutil.copy(udir / 'assembly.gfa', wdir / 'assembly.gfa')
            write_fasta(wdir / 'original.fasta', ((c, seqs[c]) for c in e['members']))
            runner.run('minimap2', ['-x', 'asm5', '-c', '--secondary=no', '-t', args.threads,
                                    wdir / 'assembly.fasta', wdir / 'original.fasta'], f'{name}_retention',
                       stdout=wdir / 'original_to_candidate.paf')
            retained = paf_retention(wdir / 'original_to_candidate.paf', {c: lengths[c] for c in e['members']})
            clen = {c: len(s) for c, s in cand}; cgc = {c: s.count('G') + s.count('C') for c, s in cand}
            after = tb.bin_stats(list(clen), clen, cgc)
            row.update(contigs_after=after['num_contigs'], total_length_after_bp=after['total_length_bp'],
                       N50_after=after['N50'], original_sequence_retained=round(retained, 4),
                       recovered_fraction=round(after['total_length_bp'] / max(1, before['total_length_bp']), 4))
            row['outcome'] = candidate_decision(before, after, retained, args.min_retained, args.max_length_inflation)
            row['reason'] = f'retained {retained:.1%}; length ratio {row["recovered_fraction"]:.2f}'
            (wdir / 'original.fasta').unlink(missing_ok=True)
            if row['outcome'] == 'accepted': accepted[name] = wdir / 'assembly.fasta'
            note(f'     {row["outcome"]}: {row["reason"]}; contigs {before["num_contigs"]}->{after["num_contigs"]}, '
                 f'N50 {before["N50"]:,}->{after["N50"]:,}')
        except (RuntimeError, ValueError, FileNotFoundError) as ex:
            msg = str(ex)
            row['outcome'] = 'skipped' if msg.startswith('skipped') else 'failed_fallback_original'
            row['reason'] = msg
            if not msg.startswith('skipped'): warnings.append(f'{name}: {msg}')
            note(f'     {msg}' if row['outcome'] == 'skipped' else f'     {row["outcome"]}: {msg}')
        finally:
            if udir.exists(): shutil.rmtree(udir)
        row['elapsed_seconds'] = round(time.time() - t, 1)
        for side in ('before', 'after'):   # 7.4 historical aliases
            row[f'total_length_{side}'] = row[f'total_length_{side}_bp']
        done.write_text(json.dumps(row, indent=2)); rows[name] = row
        if pool:   # pools are removed only once the bin has a recorded outcome (Ctrl-C keeps them for --resume)
            for p in pool.directory.glob(f'{name}.*.fq.gz'): p.unlink()
    fields = list(dict.fromkeys(k for r in rows.values() for k in r))
    with open(out / 'bins' / R / 'reassembly_summary.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fields, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows.values())
    (rroot / 'rank_timing.json').write_text(json.dumps(dict(
        elapsed_seconds=time.time() - started, candidate_quality_seconds=0,
        note='Rank time covers Unicycler and retention checks; the shared mapping pass is excluded.'), indent=2))
    return accepted


def stage_final(runner, args, config, out, R, groups, seqs, binnable, calls, accepted, cov, breadth,
                coding, orfs, pre, excluded, reads, warnings, ck):
    note('10   Consolidated final assembly and final bins')
    fdir = out / 'final'; adir = fdir / 'assembly'; adir.mkdir(parents=True, exist_ok=True)
    records, prov, final_calls = [], [], {}
    consensus = {n: tb.bin_consensus(m, binnable) for n, m in groups.items() if n in accepted}
    for name in sorted(groups):
        safe = tb.sanitize(name)
        if name in accepted:
            src = list((n_.split()[0], s.upper()) for n_, s in fasta(accepted[name])); typ = 'reassembled'
            src_path = accepted[name]
        else:
            src = [(c, seqs[c]) for c in groups[name]]; typ = 'original'; src_path = out / 'bins' / R / f'{safe}.fasta'
        for k, (sid, s) in enumerate(src, 1):
            fid = f'MH_{typ}_{safe}_{k:07d}'
            records.append((fid, s)); prov.append([fid, R, name, typ, str(src_path), sid])
            if typ == 'original': final_calls[fid] = binnable[sid]
            else: final_calls[fid] = consensus[name]
    write_fasta(adir / 'consolidated_contigs.fasta', records)
    with open(adir / 'contig_provenance.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['final_contig', 'source_rank', 'source_bin', 'source_type', 'source_fasta', 'source_record'])
        w.writerows(prov)
    fseqs = dict(records); fids = list(fseqs); flen = {c: len(s) for c, s in fseqs.items()}
    src_of = {p[0]: p[5] for p in prov if p[3] == 'original'}
    # Coverage: reuse preliminary values unless reassembled sequence exists (then remap).
    if accepted:
        note('     remapping QC reads to the consolidated assembly (reassembled bins present)')
        fcov, _, fbreadth = stream_mapping(runner, adir / 'consolidated_contigs.fasta', reads, args.threads,
                                           'final_mapping', strict=True, breadth=True)
        reas = [(c, fseqs[c]) for c in fids if c not in src_of]
        rdir = out / 'work' / 'final_prodigal'; rdir.mkdir(parents=True, exist_ok=True)
        write_fasta(rdir / 're.fasta', reas)
        runner.run('prodigal', ['-i', rdir / 're.fasta', '-a', rdir / 're.faa', '-o', rdir / 're.gff', '-f', 'gff',
                                '-p', 'meta', '-q'], 'final_prodigal')
        rcod, rorf = orf_table(rdir / 're.faa'); shutil.rmtree(rdir)
    else:
        fcov = {c: cov.get(src_of[c]) for c in fids}; fbreadth = {c: breadth.get(src_of[c]) for c in fids}
        rcod, rorf = Counter(), Counter()
    fcoding = Counter({c: (coding[src_of[c]] if c in src_of else rcod[c]) for c in fids})
    forfs = {c: (orfs.get(src_of[c], 0) if c in src_of else rorf.get(c, 0)) for c in fids}
    fpre = {c: pre[src_of[c]] for c in fids if c in src_of and src_of[c] in pre}
    fcalls = {c: dict(calls=v, orfs_with_hits=(calls.get(src_of[c], {}).get('orfs_with_hits', 0) if c in src_of else 0))
              for c, v in final_calls.items()}
    write_classification(fdir / 'classification' / 'contig_classification.tsv', fids, fseqs, fcalls, fcov, fbreadth,
                         fcoding, forfs, {c: v['orfs_with_hits'] for c, v in fcalls.items()}, fpre, set(), set())
    fparts = {r: tb.partition(final_calls, r, flen, args.min_bin_length, args.min_bin_contigs) for r in args.ranks}
    flabels = {r: {c: tb.fallback_label(v, r) for c, v in final_calls.items()} for r in args.ranks}
    summaries = write_rank_bins(fdir / 'bins', args.ranks, fparts, fseqs, fdir / 'classification' / 'unclassified_contigs.fasta')
    fb = fdir / 'binarena'; fb.mkdir(exist_ok=True)
    columns = tb.ordinate(fseqs, args.binarena_kmers, args.binarena_methods, args.binarena_min_length,
                          args.binarena_max_contigs, args.binarena_force, args.binarena_seed, outdir=fb, log=note)
    write_binarena(fb / 'binarena_input.tsv', fids, fseqs, fcalls, flabels, R if R in args.ranks else args.ranks[-1],
                   fparts.get(R) or fparts[args.ranks[-1]], fcov, fbreadth, fcoding, forfs, fpre, set(), columns, args.ranks)
    if args.checkm: run_checkm(runner, args, out, fdir / 'bins', args.ranks, summaries, warnings)


def run_checkm(runner, args, out, bins_dir, ranks, summaries, warnings):
    """One CheckM run over every distinct named final bin (Unclassified excluded)."""
    note('     CheckM (lineage_wf, reduced tree) on named final bins')
    inp = out / 'work' / 'checkm_input'; dest = out / 'work' / 'checkm_work'
    for d in (inp, dest):
        if d.exists(): shutil.rmtree(d)
    inp.mkdir(parents=True)
    keymap = {}
    for r in ranks:
        for name in summaries[r]:
            if name == tb.UNCLASSIFIED: continue
            key = f'{r}__{name}'; keymap[key] = (r, name)
            (inp / f'{key}.fasta').symlink_to(bins_dir / r / f'{name}.fasta')
    if not keymap: return
    table = out / 'final' / 'checkm.tsv'
    try:
        runner.run('checkm', ['lineage_wf', '-x', 'fasta', '-t', args.threads, '--pplacer_threads', min(args.threads, 4),
                              '--reduced_tree', '--tab_table', '-f', table, inp, dest], 'checkm')
    except RuntimeError as e:
        warnings.append(f'CheckM failed: {e}'); return
    finally:
        shutil.rmtree(inp, ignore_errors=True)
    with open(table) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            r, name = keymap.get(row.get('Bin Id'), (None, None))
            if not r: continue
            summaries[r][name].update(completeness_percent=row.get('Completeness', 'NA'),
                                      contamination_percent=row.get('Contamination', 'NA'),
                                      strain_heterogeneity_percent=row.get('Strain heterogeneity', 'NA'),
                                      marker_lineage=row.get('Marker lineage', 'NA'))
    for r in ranks: tsv(bins_dir / r / 'summary.tsv', list(summaries[r].values()), tb.SUMMARY_FIELDS)
    shutil.rmtree(dest, ignore_errors=True)


def stage_inspect(args, out, R, reads):
    """7.4 inspection (metahopper_inspect.py, unchanged): coverage profiles, terminal-junction
    support and optional reference trees for the endosymbiont panel + bins >=10% complete."""
    import importlib.util
    note(f'11   Inspection ({args.inspect}) with metahopper_inspect.py')
    try:
        spec = importlib.util.spec_from_file_location('metahopper_inspect', HERE / 'metahopper_inspect.py')
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        fasm = out / 'final' / 'assembly' / 'consolidated_contigs.fasta'
        result = module.inspect_run(out, mode=args.inspect, references=args.references, threads=args.threads,
                                    rank=R, assembly=fasm, r1=reads[0], r2=reads[1] if len(reads) > 1 else None,
                                    assembly_sequences={n: s for n, s in fasta(fasm)})
        for issue in result.get('errors', []): note(f'WARNING: inspection: {issue}')
    except (RuntimeError, ValueError, OSError) as e:
        note(f'WARNING: inspection could not complete: {e}')
        (out / 'inspection').mkdir(exist_ok=True)
        (out / 'inspection' / 'inspection.json').write_text(json.dumps({'errors': [str(e)]}))


def write_report(out, args):
    """Render metahopper_report.html with the unchanged 7.4 metahopper_report.py."""
    import importlib.util
    script = HERE / 'metahopper_report.py'
    try:
        spec = importlib.util.spec_from_file_location('metahopper_report', script)
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        layout = module.RunLayout(out, 'auto')
        opts = argparse.Namespace(input=out, output=out / 'metahopper_report.html', bin_set='auto', title=None,
                                  min_length=args.binarena_min_length, max_points=args.report_max_points,
                                  include_uncoordinated=False, no_fasta_stats=False, apply_split=None,
                                  split_outdir=None, contigs=None)
        (out / 'metahopper_report.html').write_text(module.build_report(layout, opts), encoding='utf-8')
        note(f'12   HTML report: {out / "metahopper_report.html"}')
    except Exception as e:  # a report failure must not fail a finished run
        note(f'WARNING: report not written ({type(e).__name__}: {e}); run metahopper_report.py -i {out}')


def cli(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0], formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-1', '--read1'); p.add_argument('-2', '--read2')
    p.add_argument('-o', '--out'); p.add_argument('-c', '--config')
    p.add_argument('--assembly', help='Skip MEGAHIT and use this FASTA')
    p.add_argument('--resume', metavar='RUN_DIR', help='Continue a run; finished stages are reused')
    p.add_argument('-t', '--threads', type=int, default=16); p.add_argument('--memory', type=int, default=64, help='GB')
    g = p.add_argument_group('contigs and prefilter')
    g.add_argument('--min-contig-len', type=int, default=1000, help='Contig length filter (bp), preliminary and reassembled')
    g.add_argument('--prefilter', choices=['tiara', 'none'], default='tiara')
    g.add_argument('--rescue-coding-density', type=float, default=0.40,
                   help='Tiara unknown/organelle contigs are still searched if coding density is greater than this')
    g.add_argument('--rescue-eukarya', action='store_true', help='Also apply --rescue-coding-density to Tiara eukarya')
    g.add_argument('--prefilter-audit', type=float, default=0.02,
                   help='Fraction of removed contigs searched anyway (same DIAMOND run) to measure losses')
    g = p.add_argument_group('DIAMOND and taxonomy vote (7.4 algorithm)')
    g.add_argument('--max-orfs-per-contig', type=orf_limit, default=10, help="ORFs searched per contig, or 'all'")
    g.add_argument('--orf-selection', choices=['spaced', 'first'], default='spaced',
                   help='Which ORFs when a contig has more than the cap: evenly spaced along the contig, or the first N')
    g.add_argument('--bitscore-range', type=float, default=0.9)
    g.add_argument('--max-hits-per-orf', type=int, default=5)
    g.add_argument('--min-support', type=float, default=0.5)
    g.add_argument('--support-denominator', choices=['all', 'named'], default='all')
    g.add_argument('--exclude-kingdoms', default='Metazoa,Viridiplantae')
    g = p.add_argument_group('bins')
    g.add_argument('--ranks', type=lambda s: [x.strip() for x in s.split(',') if x.strip()],
                   default=list(tb.DEFAULT_BIN_RANKS), help='Comma list of partition ranks')
    g.add_argument('--min-bin-length', type=int, default=50000, help='Smaller bins cascade to the next called higher rank')
    g.add_argument('--min-bin-contigs', type=int, default=1)
    g = p.add_argument_group('reassembly')
    g.add_argument('--skip-reassembly', action='store_true')
    g.add_argument('--reassemble-rank', default='genus')
    g.add_argument('--reassemble-categories', default='Bacteria,Archaea')
    g.add_argument('--max-bin-mbp', type=float, default=20)
    g.add_argument('--target-coverage', type=float, default=200)
    g.add_argument('--reassembly-timeout-h', type=float, default=4)
    g.add_argument('--mode', choices=['bold', 'normal', 'conservative'], default='bold')
    g.add_argument('--min-retained', type=float, default=0.95)
    g.add_argument('--max-length-inflation', type=float, default=1.25)
    g = p.add_argument_group('BinaRena and report')
    g.add_argument('--binarena-kmers', type=lambda s: [int(x) for x in s.split(',')], default=[4, 5, 6])
    g.add_argument('--binarena-methods', type=lambda s: [x.strip() for x in s.split(',')], default=['pca'])
    g.add_argument('--binarena-min-length', type=int, default=1000)
    g.add_argument('--binarena-max-contigs', type=int, default=20000)
    g.add_argument('--binarena-force', action='store_true')
    g.add_argument('--binarena-seed', type=int, default=42)
    g.add_argument('--binarena-include-excluded', action='store_true',
                   help='Also ordinate Tiara/kingdom-excluded contigs (they are always listed in the table)')
    g.add_argument('--inspect', choices=['endosymbionts', 'all', 'off'], default='endosymbionts',
                   help='7.4 coverage/junction inspection (needs bowtie2 + samtools; trees need --references)')
    g.add_argument('--references', type=Path, help='Reference-genome folder for inspection trees')
    g.add_argument('--report-max-points', type=int, default=25000)
    g.add_argument('--skip-report', action='store_true')
    g.add_argument('--skip-phyloflash', action='store_true')
    g.add_argument('--checkm', action='store_true')
    args = p.parse_args(argv)
    try:
        if args.resume:
            prev = json.loads((Path(args.resume) / 'run.json').read_text())['arguments']
            prev.update(resume=str(Path(args.resume).resolve()), out=str(Path(args.resume).resolve()))
            args = argparse.Namespace(**prev)
        elif not (args.read1 and args.out and args.config):
            p.error('-1, -o and -c are required unless using --resume RUN_DIR')
        main(args)
    except (ValueError, RuntimeError, OSError, ImportError) as e:
        p.exit(1, 'ERROR: ' + str(e) + '\n')


if __name__ == '__main__':
    cli()
