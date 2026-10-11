"""Taxonomy-first contig classification, rank partitions and BinaRena features.

The per-rank vote is the MetaHopper 7.4 algorithm (classify_contig): for every ORF, keep
DIAMOND hits within --bitscore-range of the best hit (at most --max-hits-per-orf), and add
each kept hit's bitscore to the taxon it names at every rank (split evenly when a hit
carries several taxa). A rank is called when the winning taxon holds at least
--min-support of the rank's weight. Ranks are voted independently.

What is new relative to 7.4:
  * hit lineages come from DIAMOND staxids + the NCBI taxdump (any rank, incl. class/order);
  * partitions use *fallback names*: a contig not called at the bin rank is named after its
    nearest called higher rank ("Enterobacteriaceae bacterium", "Sodalis sp.") instead of
    being dumped into Unclassified;
  * bins smaller than --min-bin-length cascade upward into those fallback names instead of
    into Unclassified;
  * no composition/coverage splitting or recruitment: that is left to BinaRena.
"""
from __future__ import annotations
from collections import Counter, defaultdict
from pathlib import Path

ORDER = ['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
DEFAULT_BIN_RANKS = ['phylum', 'class', 'order', 'family', 'genus', 'species']
SUFFIX = {'Bacteria': 'bacterium', 'Archaea': 'archaeon', 'Eukaryota': 'eukaryote', 'Viruses': 'virus'}
UNCLASSIFIED = 'Unclassified'


def sanitize(name):
    """Identical to MetaHopper 7.4 sanitize(): bin file names."""
    if not name: name = UNCLASSIFIED
    out = ''.join(ch if (ch.isalnum() or ch in '-._') else '_' for ch in name).strip('_')
    return out or UNCLASSIFIED


# ---------------------------------------------------------------------------------------
# Hit lineages
# ---------------------------------------------------------------------------------------
class RankNames:
    """staxids field -> {rank: tuple(unique names)}; cached."""

    def __init__(self, tax):
        self.tax, self._tid, self._field = tax, {}, {}

    def of_taxid(self, tid):
        hit = self._tid.get(tid)
        if hit is not None: return hit
        out, tax = {}, self.tax
        lin = tax.lineage(tid)
        for t in lin:
            r = tax.rank.get(t)
            if r == 'superkingdom': r = 'domain'
            if r in ORDER and r not in out: out[r] = tax.names.get(t, str(t))
        if 'domain' not in out and 10239 in tax.lineage_set(tid): out['domain'] = 'Viruses'
        self._tid[tid] = out
        return out

    def of_field(self, field):
        hit = self._field.get(field)
        if hit is not None: return hit
        per = defaultdict(dict)
        for x in field.replace(',', ';').split(';'):
            x = x.strip()
            if not x.isdigit() or x == '0': continue
            for r, n in self.of_taxid(int(x)).items(): per[r][n] = None
        out = {r: tuple(v) for r, v in per.items()}
        self._field[field] = out
        return out


def contig_of(orf_id):
    return orf_id.rsplit('_', 1)[0]


def classify_hits(hits_path, tax, ranks=ORDER, bitscore_range=0.9, max_hits=5,
                  min_support=0.5, denominator='all'):
    """Stream DIAMOND outfmt-6 (qseqid sseqid pident length evalue bitscore staxids).

    Returns {contig: {'calls': {rank: (taxon, support)}, 'orfs_with_hits': n}} plus stats.
    Hits for one ORF must be contiguous (DIAMOND writes them that way).
    """
    names = RankNames(tax)
    weight = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    total = defaultdict(lambda: defaultdict(float))
    with_hits = Counter()
    stats = Counter()

    def flush(orf, hits):
        if not hits: return
        hits.sort(key=lambda h: h[0], reverse=True)  # stable, as 7.4
        best = hits[0][0]
        if best <= 0: return
        kept = [h for h in hits if h[0] >= bitscore_range * best][:max_hits]
        if not kept: return
        c = contig_of(orf); with_hits[c] += 1
        W, T = weight[c], total[c]
        for bs, lin in kept:
            for r in ranks:
                taxa = lin.get(r, ())
                if taxa or denominator == 'all': T[r] += bs
                if taxa:
                    share = bs / len(taxa)
                    for t in taxa: W[r][t] += share

    current, buf, done = None, [], set()
    with open(hits_path) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 7: stats['malformed'] += 1; continue
            q = p[0]
            if q != current:
                flush(current, buf)
                if current is not None: done.add(current)
                if q in done: stats['noncontiguous_orfs'] += 1
                current, buf = q, []
            try: bs = float(p[5])
            except ValueError: stats['malformed'] += 1; continue
            stats['hits'] += 1
            buf.append((bs, names.of_field(p[6])))
    flush(current, buf)
    result = {}
    for c in set(weight) | set(total):
        calls = {}
        for r in ranks:
            W, T = weight[c].get(r), total[c].get(r, 0.)
            if not W or T <= 0:
                calls[r] = (UNCLASSIFIED, 0.0); continue
            taxon, w = max(W.items(), key=lambda kv: kv[1])
            s = w / T
            calls[r] = (taxon, round(s, 4)) if s >= min_support else (UNCLASSIFIED, round(s, 4))
        result[c] = dict(calls=calls, orfs_with_hits=with_hits[c])
    stats['contigs_with_hits'] = len(result)
    return result, dict(stats)


def empty_calls():
    return {r: (UNCLASSIFIED, 0.0) for r in ORDER}


# ---------------------------------------------------------------------------------------
# Fallback names and rank partitions
# ---------------------------------------------------------------------------------------
def suffix_for(calls):
    return SUFFIX.get(calls.get('domain', (UNCLASSIFIED,))[0], 'organism')


def named_at(calls, start, bin_index=None):
    """Nearest called rank at or above ORDER[start]; returns (label, rank index).
    Labels above the bin rank get fallback names ("Enterobacteriaceae bacterium",
    "Sodalis sp.", "unclassified Bacteria")."""
    bin_index = start if bin_index is None else bin_index
    for i in range(start, -1, -1):
        t = calls.get(ORDER[i], (UNCLASSIFIED,))[0]
        if t == UNCLASSIFIED: continue
        if i == bin_index: return t, i
        if ORDER[i] == 'domain': return f'unclassified {t}', i
        if ORDER[i] == 'genus': return f'{t} sp.', i
        return f'{t} {suffix_for(calls)}', i
    return UNCLASSIFIED, -1


def fallback_label(calls, rank):
    return named_at(calls, ORDER.index(rank))[0]


def partition(contigs, rank, lengths, min_bin_length=50000, min_bin_contigs=1):
    """contigs: {id: calls}. Returns {id: label}.

    Bins smaller than min_bin_length (or with fewer than min_bin_contigs) cascade: each
    member is renamed after its next called higher rank, and groups are re-evaluated until
    stable. Domain-level groups ("unclassified Bacteria") are kept whatever their size, so
    only contigs with no call at any rank end up in Unclassified."""
    bi = ORDER.index(rank)
    level, label = {}, {}
    for c, calls in contigs.items():
        label[c], level[c] = named_at(calls, bi, bi)
    for _ in range(len(ORDER) + 1):
        size, count, top = Counter(), Counter(), {}
        for c, l in label.items(): size[l] += lengths[c]; count[l] += 1; top[l] = level[c]
        small = {l for l in size if l != UNCLASSIFIED and top[l] > 0
                 and (size[l] < min_bin_length or count[l] < min_bin_contigs)}
        if not small: break
        for c, l in label.items():
            if l in small: label[c], level[c] = named_at(contigs[c], level[c] - 1, bi)
    return label


def bin_consensus(members, contigs):
    """Per-rank unanimity among called members (7.4 bin_consensus_lineages)."""
    out = {}
    for r in ORDER:
        called = {contigs[c][r][0] for c in members if contigs[c][r][0] != UNCLASSIFIED}
        out[r] = (next(iter(called)), 1.0) if len(called) == 1 else (UNCLASSIFIED, 0.0)
    return out


def majority_domain(members, contigs, lengths):
    bp = Counter()
    for c in members: bp[contigs[c]['domain'][0]] += lengths[c]
    return bp.most_common(1)[0][0] if bp else UNCLASSIFIED


# ---------------------------------------------------------------------------------------
# Stats (7.4 summary.tsv fields)
# ---------------------------------------------------------------------------------------
SUMMARY_FIELDS = ['bin', 'num_contigs', 'total_length_bp', 'largest_contig_bp', 'N50', 'L50',
                  'GC_percent', 'completeness_percent', 'contamination_percent',
                  'strain_heterogeneity_percent', 'marker_lineage']


def n50(lengths):
    ls = sorted(lengths, reverse=True); total = sum(ls); run = 0
    for i, l in enumerate(ls, 1):
        run += l
        if run * 2 >= total: return l, i
    return 0, 0


def bin_stats(members, lengths, gc_counts):
    ls = [lengths[c] for c in members]; total = sum(ls)
    v, l50 = n50(ls)
    gc = sum(gc_counts[c] for c in members)
    return dict(num_contigs=len(ls), total_length_bp=total, largest_contig_bp=max(ls, default=0),
                N50=v, L50=l50, GC_percent=round(100 * gc / total, 2) if total else 0.0)


# ---------------------------------------------------------------------------------------
# BinaRena features (7.4 semantics, numpy k-mer counting)
# ---------------------------------------------------------------------------------------
def _lut():  # noqa
    import numpy as np
    lut = np.full(256, -1, dtype=np.int8)
    for i, b in enumerate(b'ACGT'): lut[b] = i
    for i, b in enumerate(b'acgt'): lut[b] = i
    return lut


_LUT = None


def kmer_counts(seq, k):
    """Canonical-strand k-mer counts identical to 7.4 count_kmers (both strands counted,
    windows containing non-ACGT skipped)."""
    import numpy as np
    global _LUT
    if _LUT is None: _LUT = _lut()
    size = 4 ** k
    a = _LUT[np.frombuffer(seq.encode(), dtype=np.uint8)]
    n = a.size
    if n < k: return np.zeros(size, dtype=np.int64)
    valid = a >= 0
    bad = np.concatenate(([0], np.cumsum(~valid)))
    ok = (bad[k:] - bad[:-k]) == 0
    c = np.where(valid, a, 0).astype(np.int64)
    m = n - k + 1
    fwd = np.zeros(m, dtype=np.int64); rev = np.zeros(m, dtype=np.int64)
    for j in range(k):
        s = c[j:j + m]
        fwd = (fwd << 2) + s
        rev += (3 - s) << (2 * j)
    return np.bincount(fwd[ok], minlength=size) + np.bincount(rev[ok], minlength=size)


def gc_degenerate(seq):
    """7.4 count_gc_degenerate: G+C equivalents (IUPAC-weighted count, not a fraction)."""
    res = 0
    for c, w in (('G', 6), ('C', 6), ('S', 6), ('R', 3), ('Y', 3), ('K', 3), ('M', 3), ('N', 3),
                 ('D', 2), ('H', 2), ('B', 4), ('V', 4)):
        res += seq.count(c) * w
    return res / 6


METHOD_SUFFIXES = {'pca': ('PC1', 'PC2'), 'tsne': ('tsne1', 'tsne2'), 'umap': ('UM1', 'UM2')}


def ordinate(seqs, ks=(4, 5, 6), methods=('pca',), min_len=1000, max_contigs=20000, force=False,
             seed=42, perplexity=30.0, neighbors=15, outdir=None, log=print):
    """7.4 extract_binarena_features + reduce_kmer_dimensions.

    seqs: {id: upper-case sequence} of contigs eligible for ordination. Returns
    {column_name: {id: value}} and writes kmer_<k>.<method>.tsv like 7.4."""
    import numpy as np
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    columns = {}
    eligible = [c for c, s in seqs.items() if len(s) >= min_len]
    if not eligible:
        log(f'     BinaRena: no contigs >= {min_len} bp; skipping k-mer composition'); return columns
    requested = list(methods)
    if len(eligible) > max_contigs and not force:
        dropped = [m for m in requested if m in ('tsne', 'umap')]
        if dropped: log(f'     BinaRena: {len(eligible):,} contigs > {max_contigs:,}; PCA only (skipping {"/".join(dropped)})')
        requested = [m for m in requested if m == 'pca']
    for k in ks:
        ids, rows = [], []
        for c in eligible:
            v = kmer_counts(seqs[c], k); t = v.sum()
            if t <= 0: continue
            ids.append(c); rows.append((v / t).astype(np.float64))
        if not ids: continue
        data = np.vstack(rows); del rows
        n = data.shape[0]
        scaled = StandardScaler().fit_transform(data); del data
        results = {}
        n_pre = int(min(50, n - 1, scaled.shape[1]))
        pre = PCA(n_components=n_pre, random_state=seed).fit_transform(scaled) if n_pre >= 2 and \
            any(m in requested for m in ('tsne', 'umap')) else scaled
        if 'pca' in requested and n >= 2:
            results['pca'] = PCA(n_components=2, random_state=seed).fit_transform(scaled)
        if 'tsne' in requested and n >= 10:
            from sklearn.manifold import TSNE
            p = float(min(perplexity, max(5.0, (n - 1) / 3.0)))
            results['tsne'] = TSNE(n_components=2, perplexity=p, init='pca', random_state=seed).fit_transform(pre)
        if 'umap' in requested and n >= 10:
            try:
                import umap
                results['umap'] = umap.UMAP(n_components=2, n_neighbors=int(min(neighbors, n - 1)),
                                            random_state=seed).fit_transform(pre)
            except ImportError: log('     BinaRena: umap-learn unavailable; skipping UMAP')
        for method, xy in results.items():
            s1, s2 = METHOD_SUFFIXES[method]
            if outdir is not None:
                with open(Path(outdir) / f'kmer_{k}.{method}.tsv', 'w') as fh:
                    fh.write(f'ID\t{k}{s1}\t{k}{s2}\n')
                    for c, (x, y) in zip(ids, xy): fh.write(f'{c}\t{x:.6f}\t{y:.6f}\n')
            columns[f'{k}{s1}'] = {c: float(x) for c, (x, _) in zip(ids, xy)}
            columns[f'{k}{s2}'] = {c: float(y) for c, (_, y) in zip(ids, xy)}
        log(f'     BinaRena: {k}-mer {"/".join(results) or "no"} ordination on {n:,} contigs')
    return columns
