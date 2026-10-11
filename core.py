"""New taxonomy-first binning and streaming read recruitment implementation."""
from __future__ import annotations
import csv, gzip, itertools, math, re, subprocess, shutil
from collections import Counter, defaultdict, OrderedDict
from pathlib import Path
import numpy as np


def opn(path, mode='rt'):
    return gzip.open(path, mode) if str(path).endswith('.gz') else open(path, mode)


def fasta(path):
    with opn(path) as f:
        name, chunks = None, []
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, ''.join(chunks).upper()


def tsv(path, rows, fields=None):
    rows = list(rows)
    fields = fields or (list(rows[0]) if rows else [])
    with opn(path, 'wt') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore', lineterminator='\n')
        w.writeheader(); w.writerows(rows)


def gzcopy(src, dest):
    with open(src, 'rb') as fi, gzip.open(dest, 'wb', compresslevel=3) as fo:
        shutil.copyfileobj(fi, fo)


class Taxonomy:
    def __init__(self, directory):
        directory = Path(directory)
        self.parent, self.rank, self.names, self.merged = {}, {}, {}, {}
        # PERF: lineages were recomputed by walking the tree on every call (millions of
        # times in parse_hits/bin_contigs). Cache them; results are identical.
        self._lineage, self._lineage_set = {}, {}
        with open(directory / 'nodes.dmp') as f:
            for line in f:
                p = [x.strip() for x in line.split('|')]
                self.parent[int(p[0])] = int(p[1]); self.rank[int(p[0])] = p[2]
        with open(directory / 'names.dmp') as f:
            for line in f:
                p = [x.strip() for x in line.split('|')]
                if p[3] == 'scientific name': self.names[int(p[0])] = p[1]
        if (directory / 'merged.dmp').exists():
            with open(directory / 'merged.dmp') as f:
                for line in f:
                    p = line.split('|'); self.merged[int(p[0])] = int(p[1])

    def resolve(self, tid):
        seen = set()
        while tid in self.merged and tid not in seen:
            seen.add(tid); tid = self.merged[tid]
        return tid if tid in self.parent else 0

    def lineage(self, tid):
        cached = self._lineage.get(tid)
        if cached is not None: return cached
        start = tid
        tid, result, seen = self.resolve(tid), [], set()
        while tid and tid not in seen:
            result.append(tid); seen.add(tid)
            tid = self.parent.get(tid, 0)
        result = tuple(result)
        self._lineage[start] = result
        return result

    def lineage_set(self, tid):
        cached = self._lineage_set.get(tid)
        if cached is None:
            cached = self._lineage_set[tid] = frozenset(self.lineage(tid))
        return cached

    def lca(self, tids):
        # PERF: NR hits carry many duplicate taxids; the LCA of a set equals the LCA of
        # the list, so deduplicate before intersecting cached lineage sets.
        unique = [t for t in dict.fromkeys(tids) if self.lineage(t)]
        if not unique: return 0
        common = self.lineage_set(unique[0]).intersection(*(self.lineage_set(t) for t in unique[1:]))
        return next((t for t in self.lineage(unique[0]) if t in common), 0)

    def category(self, tid):
        lin = self.lineage_set(tid)
        if 33208 in lin or 3193 in lin: return 'excluded_animal_or_land_plant'
        if 2 in lin: return 'Bacteria'
        if 2157 in lin: return 'Archaea'
        if 4751 in lin: return 'Fungi'
        if 2759 in lin: return 'Other_eukaryote_putative_protist'
        return 'Unknown_or_other'

    def anchor(self, tid):
        lin = self.lineage(tid)
        for rank in ('genus', 'family', 'order'):
            x = next((t for t in lin if self.rank.get(t) == rank), 0)
            if x: return x
        return 0

    def consensus(self, votes, threshold=.7):
        """One unit vote per informative ORF; no best-hit score domination."""
        votes = [self.resolve(t) for t in votes if self.resolve(t)]
        if not votes: return 0, 0.
        support = Counter(t for v in votes for t in self.lineage(v))
        candidates = [t for t, n in support.items() if n / len(votes) >= threshold]
        tid = max(candidates, key=lambda t: len(self.lineage(t)), default=0)
        return tid, support[tid] / len(votes) if tid else 0.


def diamond_rows(path, stats):
    """Preserve empty TSV fields; fail on damaged/wrong-format rows with a line number."""
    with open(path) as f:
        for number, line in enumerate(f, 1):
            if not line.strip(): continue
            # rstrip() without arguments removes the TAB marking an empty final taxid.
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != 3:
                raise RuntimeError(f'DIAMOND row {number} in {path}: expected 3 tab-separated '
                                   f'columns (qseqid, bitscore, staxids), found {len(fields)}. '
                                   'File retained; inspect this row before rerunning DIAMOND.')
            gene, score_text, taxids = fields
            try:
                score = float(score_text)
            except ValueError:
                raise RuntimeError(f'DIAMOND row {number} in {path}: invalid bitscore {score_text!r}') from None
            if not gene or not math.isfinite(score) or score < 0:
                raise RuntimeError(f'DIAMOND row {number} in {path}: invalid query ID or bitscore')
            tids = []
            for token in taxids.split(';'):
                token = token.strip()
                if token in ('', '0', 'N/A', 'NA', '*'): continue
                if not token.isdigit():
                    raise RuntimeError(f'DIAMOND row {number} in {path}: invalid taxonomy ID {token!r}')
                tids.append(int(token))
            stats['hit_rows'] += 1
            if not tids: stats['rows_without_taxid'] += 1
            yield gene, score, tids


def parse_hits(path, tax, diagnostics=None):
    votes = defaultdict(list)
    stats = dict(hit_rows=0, rows_without_taxid=0, unresolved_taxids=0,
                 queries_with_hits=0, informative_queries=0)
    try:
        for gene, group in itertools.groupby(diamond_rows(path, stats), lambda p: p[0]):
            hits = list(group); stats['queries_with_hits'] += 1
            best = max(p[1] for p in hits)
            tids = []
            for _, score, raw_tids in hits:
                if score >= best * .9:
                    for raw_tid in raw_tids:
                        tid = tax.resolve(raw_tid)
                        if tid: tids.append(tid)
                        else: stats['unresolved_taxids'] += 1
            tid = tax.lca(tids)
            if tid:
                stats['informative_queries'] += 1
                votes[gene.rsplit('_', 1)[0]].append(tid)
        if stats['queries_with_hits'] and not stats['informative_queries']:
            raise RuntimeError('DIAMOND hits contain no usable taxids. Check the taxonomy-enabled '
                               'database and matching taxdump before rebuilding or rerunning.')
    finally:
        if diagnostics is not None:
            import json
            Path(diagnostics).write_text(json.dumps(stats, indent=2)+'\n')
    result = {}
    for contig, tids in votes.items():
        tid, confidence = tax.consensus(tids)
        result[contig] = dict(taxid=tid, taxonomy=tax.names.get(tid, 'unclassified'),
                             tax_support=confidence, informative_orfs=len(tids), category=tax.category(tid))
    return result


# Reverse-complement-collapsed tetranucleotide frequencies (136 dimensions).
_KMERS = [''.join(k) for k in itertools.product('ACGT', repeat=4)]
_COMP = str.maketrans('ACGT', 'TGCA')
_CANON = sorted({min(k, k.translate(_COMP)[::-1]) for k in _KMERS})
_KINDEX = {k: _CANON.index(min(k, k.translate(_COMP)[::-1])) for k in _KMERS}


def sequence_features(seq):
    valid = sum(seq.count(x) for x in 'ACGT')
    gc = 100 * (seq.count('G') + seq.count('C')) / valid if valid else 0.
    counts = np.zeros(136, dtype=np.float32)
    # Count non-overlapping occurrences at all four offsets in C regex engine.
    for offset in range(4):
        for k, n in Counter(re.findall('.{4}', seq[offset:])).items():
            if k in _KINDEX: counts[_KINDEX[k]] += n
    return gc, counts / counts.sum() if counts.sum() else counts


def feature_table(path, calls, depths):
    rows, vectors = [], []
    for name, seq in fasta(path):
        gc, vector = sequence_features(seq)
        row = dict(contig=name, length=len(seq), gc=gc, coverage=depths.get(name, 0.),
                   taxid=0, taxonomy='unclassified', tax_support=0., informative_orfs=0,
                   category='Unknown_or_other', bin='unassigned', assignment='unassigned')
        row.update(calls.get(name, {})); rows.append(row); vectors.append(vector)
    return rows, np.asarray(vectors, dtype=np.float32)


def bin_contigs(rows, frequencies, tax, max_bin_bp=20_000_000, warnings=None):
    """Taxonomic hard partitions, conservative within-partition splits, bounded recruitment."""
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from sklearn.cluster import MiniBatchKMeans
    from sklearn.metrics import silhouette_score
    n = len(rows)
    if not n: return {}
    dimensions = min(6, n, frequencies.shape[1])
    if n > 1 and np.any(np.ptp(frequencies, axis=0) > 0):
        pc = PCA(n_components=dimensions, svd_solver='randomized', random_state=17).fit_transform(np.sqrt(frequencies))
    else:
        pc = np.zeros((n, dimensions))
    raw = np.column_stack(([math.log1p(r['coverage']) for r in rows], [r['gc'] for r in rows], pc))
    x = StandardScaler().fit_transform(raw)
    # Coverage and GC each get one unit; all TNF axes together get one unit.
    x[:, 2:] /= math.sqrt(dimensions)
    for i, r in enumerate(rows):
        r.update(pc1=float(pc[i, 0]), pc2=float(pc[i, 1]) if dimensions > 1 else 0.)
    groups = defaultdict(list)
    for i, r in enumerate(rows):
        anchor = tax.anchor(r['taxid'])
        r['anchor_taxid'] = anchor
        if r['category'].startswith('excluded'):
            r['bin'] = 'excluded'; r['assignment'] = 'excluded'; continue
        if anchor and r['informative_orfs'] >= 2 and r['tax_support'] >= .7:
            groups[anchor].append(i)
    seeds = []
    for tid, ids in sorted(groups.items()):
        labels = np.zeros(len(ids), dtype=int)
        # Cap fitting and silhouette sample size; never form a whole-metagenome distance matrix.
        if len(ids) >= 20:
            rng = np.random.default_rng(17)
            sample = rng.choice(len(ids), min(5000, len(ids)), replace=False)
            baseline = .5
            for k in range(2, min(6, len(ids) // 10) + 1):
                model = MiniBatchKMeans(n_clusters=k, random_state=17, n_init=3, batch_size=1024)
                model.fit(x[np.array(ids)[sample]])
                lab = model.predict(x[ids])
                sizes = np.bincount(lab, minlength=k)
                if min(sizes) < 5 or len(set(lab[sample])) < 2: continue
                score = silhouette_score(x[np.array(ids)[sample]], lab[sample], sample_size=min(1000, len(sample)), random_state=17)
                if score > baseline:
                    baseline, labels = score, lab
        for label in sorted(set(labels)):
            members = [ids[j] for j in range(len(ids)) if labels[j] == label]
            center = np.mean(x[members], axis=0)
            d = np.linalg.norm(x[members] - center, axis=1)
            radius = max(.75, float(np.quantile(d, .95)) * 1.5)
            seeds.append((tid, members, center, radius))
    bins = {}
    for number, (tid, ids, center, radius) in enumerate(seeds, 1):
        bid = f'bin_{number:04d}'
        bins[bid] = dict(taxid=tid, name=tax.names.get(tid, str(tid)), members=list(ids), center=center, radius=radius)
        for i in ids: rows[i].update(bin=bid, assignment='taxonomy_seed')
    # Only assign compatible ancestors/descendants, never override conflicting taxonomy.
    # PERF: the original looped over every (unassigned contig x bin) pair in Python,
    # recomputing two lineages per pair: O(contigs * bins) interpreter work, which is
    # days for ~1e6 contigs and thousands of taxonomy seeds. Distances are now screened
    # as matrix blocks; exact per-pair norms/tie rules are kept for the few bins in range.
    if not bins: return bins
    # GUARD (v1.0.2): contigs whose consensus is only root / 'cellular organisms'
    # (typical of host contigs with mixed eukaryotic + HGT-like ORF votes) used to count
    # as 'compatible' with every bin. Ancestors must now lie inside the bin's domain.
    domains = (2, 2157, 2759, 10239)
    def compatible_tax(tid, btid):
        if not tid: return True
        if tid in tax.lineage_set(btid):      # contig taxon is an ancestor of the bin taxon
            return any(d in tax.lineage_set(tid) for d in domains)
        return btid in tax.lineage_set(tid)  # contig taxon is a descendant
    bids = list(bins)
    centers = np.stack([bins[b]['center'] for b in bids]).astype(np.float64)
    radii = np.array([bins[b]['radius'] for b in bids], dtype=np.float64)
    btax = [bins[b]['taxid'] for b in bids]
    csq = (centers ** 2).sum(1)
    todo = [i for i, r in enumerate(rows) if r['bin'] == 'unassigned']
    proposals = defaultdict(list)
    step = max(1, min(4096, 2_000_000 // max(1, len(bids))))
    for start in range(0, len(todo), step):
        chunk = todo[start:start+step]
        xc = x[chunk]
        d2 = (xc ** 2).sum(1)[:, None] + csq[None, :] - 2 * xc @ centers.T
        near = d2 <= (radii[None, :] + 1e-6) ** 2 + 1e-6   # generous screen; exact test below
        for row_index, i in enumerate(chunk):
            js = np.flatnonzero(near[row_index])
            if not len(js): continue
            tid = rows[i]['taxid']
            compatible = []
            for j in js:
                if not compatible_tax(tid, btax[j]): continue
                dist = float(np.linalg.norm(x[i] - centers[j]))
                if dist <= radii[j]: compatible.append((dist, bids[j]))
            compatible.sort()
            if not compatible: continue
            if len(compatible) > 1 and compatible[0][0] > .75 * compatible[1][0]: continue
            proposals[compatible[0][1]].append((compatible[0][0], i))
    # GUARD (v1.0.2): a prokaryotic bin cannot be hundreds of Mbp. Accept recruits
    # nearest-first only while the bin stays under max_bin_bp; report what was refused.
    for bid in bids:
        b = bins[bid]; size = sum(rows[i]['length'] for i in b['members'])
        refused = refused_bp = 0
        for dist, i in sorted(proposals.get(bid, ())):
            n = rows[i]['length']
            if size + n > max_bin_bp:
                refused += 1; refused_bp += n; rows[i]['assignment'] = 'recruitment_refused_size_cap'; continue
            size += n
            rows[i].update(bin=bid, assignment='compatible_feature_recruitment')
            b['members'].append(i)
        b['length_bp'] = size
        if refused and warnings is not None:
            warnings.append(f'{bid} ({b["name"]}): refused {refused:,} recruited contigs '
                            f'({refused_bp/1e6:.1f} Mbp) at the {max_bin_bp/1e6:.0f} Mbp bin-size cap.')
    return bins


def phyloflash_abundance(path):
    """phyloFlash's headerless full NTU CSV: lineage,count."""
    rows = []
    with open(path, newline='') as f:
        for p in csv.reader(f):
            if len(p) < 2: continue
            try: count = float(p[-1])
            except ValueError: continue
            lineage = ','.join(p[:-1]).strip()
            parts = {v.strip().casefold() for v in lineage.split(';')}
            excluded = bool(parts & {'metazoa', 'animalia', 'embryophyta'})
            allowed = bool(parts & {'bacteria','archaea','eukaryota','eukaryotes','eukarya'}) and not excluded
            rows.append(dict(lineage=lineage, ssu_count=count, retained=allowed,
                             reason='retained' if allowed else 'animal/land-plant-or-unresolved-lineage'))
    total = sum(r['ssu_count'] for r in rows)
    retained = sum(r['ssu_count'] for r in rows if r['retained'])
    rows.sort(key=lambda r: -r['ssu_count'])
    rank = 0
    for r in rows:
        if r['retained']: rank += 1
        r.update(rank_retained=rank if r['retained'] else '',
                 fraction_all_ssu=r['ssu_count']/total if total else 0.,
                 fraction_retained_ssu=r['ssu_count']/retained if r['retained'] and retained else 0.)
    return rows


def corroboration(tid, tax, abundance):
    # Exact full lineage-token matches, not fuzzy substring claims.
    names = {tax.names[t].casefold() for t in tax.lineage(tid)
             if tax.rank.get(t) in ('species','genus','family','order') and t in tax.names}
    matches = [r for r in abundance if r['retained'] and names.intersection(x.strip().casefold() for x in r['lineage'].split(';'))]
    return sum(r['ssu_count'] for r in matches), min((r['rank_retained'] for r in matches), default='')


class BufferedWriters:
    """Many output files without many open handles or per-record open/close.

    PERF: the original kept an LRU of 64 open gzip streams. With more than ~32 bins
    (R1+R2 handles each) and reads arriving in random bin order, nearly every read
    closed one gzip stream and opened another: two filesystem metadata operations, a
    zlib init/flush and a new gzip member per read -- ~1e8 times per pass, and far
    worse on NFS/Lustre. Records are now buffered per file and appended in large
    blocks. File contents (decompressed) and record order are unchanged.
    """
    def __init__(self, compress=True, per_file=4 << 20, total=512 << 20):
        self.compress, self.per_file, self.total_limit = compress, per_file, total
        self.buffers, self.sizes, self.total = defaultdict(list), Counter(), 0
        self.touched = set()

    def write(self, path, text):
        self.buffers[path].append(text); n = len(text)
        self.sizes[path] += n; self.total += n
        if self.sizes[path] >= self.per_file: self.flush(path)
        elif self.total >= self.total_limit: self.flush_all()

    def flush(self, path):
        data = ''.join(self.buffers.pop(path, ()))
        self.total -= self.sizes.pop(path, 0)
        if not data: return
        mode = 'a' if path in self.touched else 'w'
        self.touched.add(path)
        if self.compress:
            with gzip.open(path, mode+'t', compresslevel=1) as h: h.write(data)
        else:
            with open(path, mode) as h: h.write(data)

    def flush_all(self):
        for path in list(self.buffers): self.flush(path)


class Pools:
    """Per-bin FASTQ pools written through buffered appends."""
    def __init__(self, directory, max_open=64):
        self.directory = Path(directory); self.directory.mkdir(parents=True, exist_ok=True)
        self.writers = BufferedWriters(compress=True); self.counts = Counter(); self.bases = Counter()

    def write(self, bid, suffix, record):
        seq, qual = record.get_forward_sequence(), record.get_forward_qualities()
        if seq is None or qual is None: raise RuntimeError('Alignment lacks FASTQ sequence/quality.')
        # PERF: C-level quality encoding instead of a Python generator per base.
        self.writers.write(self.directory / f'{bid}.{suffix}.fq.gz',
                           f'@{record.query_name}\n{seq}\n+\n{_qualstring(qual)}\n')
        self.counts[bid] += 1; self.bases[bid] += len(seq)

    def close(self):
        self.writers.flush_all()


def downsample_pool(directory, bid, suffixes, fraction, seed=17):
    """Keep a deterministic read-name-hashed fraction; mates stay paired across files."""
    import zlib
    cutoff = int(fraction * 2**32); kept = 0
    for suffix in suffixes:
        src = Path(directory) / f'{bid}.{suffix}.fq.gz'; tmp = src.with_suffix('.tmp.gz')
        with gzip.open(src, 'rt') as fi, gzip.open(tmp, 'wt', compresslevel=1) as fo:
            for header in fi:
                block = header + fi.readline() + fi.readline() + fi.readline()
                if zlib.crc32(header.split()[0].encode(), seed) < cutoff:
                    fo.write(block); kept += 1
        tmp.replace(src)
    return kept // len(suffixes)


def _qualstring(qual):
    import pysam
    return pysam.qualities_to_qualitystring(qual)


def route_group(records, paired, contig_bins, pools):
    primary = [r for r in records if not r.is_secondary and not r.is_supplementary]
    good = [r for r in primary if not r.is_unmapped and r.mapping_quality >= 10]
    if paired:
        if len(primary) != 2: raise RuntimeError('Expected two primary records per QNAME; check paired FASTQ synchronization and duplicate IDs.')
        if len(good) != 2: return 'unmapped_or_ambiguous'
        left = next((r for r in good if r.is_read1), None)
        right = next((r for r in good if r.is_read2), None)
        if left is None or right is None: raise RuntimeError('Invalid read-pair flags.')
        a, b = contig_bins.get(left.reference_name), contig_bins.get(right.reference_name)
        if not a or a != b: return 'cross_bin'
        pools.write(a, 'R1', left); pools.write(a, 'R2', right)
    else:
        if len(primary) != 1: raise RuntimeError('Duplicate QNAME in single-end reads.')
        if not good: return 'unmapped_or_ambiguous'
        bid = contig_bins.get(good[0].reference_name)
        if not bid: return 'unassigned'
        pools.write(bid, 'SE', good[0])
    return 'accepted'


def stream_mapping(runner, reference, reads, threads, logname, strict=False, contig_bins=None, pools=None,
                   breadth=False, keep_index=False):
    """BWA emits each template's records together; stream without a disk SAM/BAM.

    Returns ({contig: mean depth}, routing counts) or, with breadth=True,
    ({contig: mean depth}, routing counts, {contig: covered fraction}). Depth counts
    M/=/X bases of primary alignments with MAPQ >= 10; breadth marks their reference spans.
    """
    import pysam
    if not keep_index or not Path(str(reference)+'.bwt').exists():
        runner.run('bwa', ['index', reference], logname+'_index')
    args = ['mem','-t',str(threads)]
    if strict: args += ['-k','50','-B','10','-O','10','-T','90']
    args += [str(reference)] + list(map(str, reads))
    cmd = runner.command('bwa', args)
    runner.record(cmd)
    bases, counts = Counter(), Counter()
    log = runner.logdir / (logname+'.log')
    with open(log, 'wb') as err:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err, env=runner.env)
        try:
            with pysam.AlignmentFile(proc.stdout, 'r', check_sq=True) as sam:
                lengths = dict(zip(sam.references, sam.lengths))
                if breadth:
                    offset, total = {}, 0
                    for name in sam.references: offset[name] = total; total += lengths[name]
                    covered = bytearray(total); ones = memoryview(b'\x01' * max(lengths.values(), default=1))
                for _, group in itertools.groupby(sam, lambda r: r.query_name):
                    records = list(group)
                    for r in records:
                        if not (r.is_unmapped or r.is_secondary or r.is_supplementary) and r.mapping_quality >= 10:
                            bases[r.reference_name] += sum(n for op,n in (r.cigartuples or []) if op in (0,7,8))
                            if breadth:
                                a = offset[r.reference_name] + r.reference_start; b = offset[r.reference_name] + r.reference_end
                                covered[a:b] = ones[:b-a]
                    if pools: counts[route_group(records, len(reads)==2, contig_bins, pools)] += 1
            code = proc.wait()
            if code: raise RuntimeError(f'BWA failed ({code}); see {log}')
        except BaseException:
            proc.kill(); proc.wait(); raise
        finally:
            if pools: pools.close()
            if not keep_index:
                for suffix in ('.amb','.ann','.bwt','.pac','.sa'):
                    Path(str(reference)+suffix).unlink(missing_ok=True)
    depth = {c: bases[c]/l for c,l in lengths.items()}
    if not breadth: return depth, dict(counts)
    import numpy as np
    names = list(lengths)
    if not names or not total: return depth, dict(counts), {}
    arr = np.frombuffer(covered, dtype=np.uint8)   # no copy of the (assembly-sized) mask
    starts = np.array([offset[n] for n in names], dtype=np.int64)
    sums = np.add.reduceat(arr, starts, dtype=np.int64)
    return depth, dict(counts), {n: float(sums[i]) / lengths[n] for i, n in enumerate(names)}


def metrics(path, depths, calls, tax):
    lengths, gc_bases, valid_bases, depth_bases = [], 0, 0, 0.
    weighted = Counter(); informative = 0
    for name, seq in fasta(path):
        n = len(seq); lengths.append(n)
        gc_bases += seq.count('G')+seq.count('C'); valid_bases += sum(seq.count(b) for b in 'ACGT')
        depth_bases += n * depths.get(name,0.)
        tid = calls.get(name,{}).get('taxid',0)
        if tid:
            informative += n
            for anc in tax.lineage(tid): weighted[anc] += n
    total = sum(lengths); running = 0; n50 = l50 = 0
    for i,n in enumerate(sorted(lengths,reverse=True),1):
        running += n
        if running >= total/2: n50, l50 = n,i; break
    tids = [t for t,v in weighted.items() if informative and v/informative >= .7]
    tid = max(tids,key=lambda t:len(tax.lineage(t)),default=0)
    return dict(length_bp=total,contigs=len(lengths),n50=n50,l50=l50,
                contigs_per_mbp=len(lengths)*1e6/total if total else 0.,
                gc=100*gc_bases/valid_bases if valid_bases else 0.,
                coverage=depth_bases/total if total else 0.,taxid=tid,
                taxonomy=tax.names.get(tid,'unclassified'), category=tax.category(tid),
                classified_bp_fraction=informative/total if total else 0.,
                tax_support=weighted[tid]/informative if informative and tid else 0.,lengths=lengths)
