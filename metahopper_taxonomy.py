"""Compact bin summaries and exact, single-rank assembly partitions for reports."""
from collections import defaultdict


def is_unbinned(name):
    """Domain-only assignments remain bins; only the generic catch-all is unbinned."""
    return str(name or '').strip().casefold() in {'', 'unclassified', 'unbinned', 'unassigned'}


def taxonomy_summary(sample, stage, rows):
    fields = ('rank', 'bin', 'total_length_bp', 'num_contigs', 'N50', 'GC_percent',
              'completeness_percent', 'contamination_percent', 'strain_heterogeneity_percent')
    return dict(schema='metahopper-taxonomy-1', sample=sample, stage=stage,
                bins=[{key: row.get(key) for key in fields} for row in rows
                      if not is_unbinned(row.get('bin'))
                      and row.get('rank') in {'phylum', 'class', 'order', 'family', 'genus', 'species'}])


def assembly_partition(layout, rows, iter_fasta, read_tsv):
    """Account for each assembly base exactly once at each rank, including unbinned.

    Summary totals alone cannot prove a partition. Use sequence IDs and lengths from
    the full displayed assembly and membership tables (or the bin FASTAs). Reject
    overlapping or incompatible inputs instead of drawing a misleading flow.
    """
    source = layout.default_contigs()
    result = dict(schema='metahopper-partition-1', total_bp=None, source=None, ranks={})
    if source is None:
        result['error'] = 'The full assembly FASTA is unavailable. Bin tables remain available.'
        return result
    lengths = {}
    for identifier, sequence in iter_fasta(source):
        if identifier in lengths:
            result['error'] = 'The full assembly contains duplicate contig IDs; a partition cannot be verified.'
            return result
        lengths[identifier] = len(sequence)
    result.update(total_bp=sum(lengths.values()), source=source.relative_to(layout.root).as_posix())
    grouped = defaultdict(list)
    for row in rows:
        grouped[row['rank']].append(row)
    for rank, bins in grouped.items():
        members = defaultdict(list)
        for item in read_tsv(layout.bins_dir / rank / 'bin_membership.tsv'):
            members[item.get('bin')].append(item.get('contig'))
        assigned, binned, parts, error = set(), set(), [], None
        for row in bins:
            name = row['bin']
            if name in members:
                records = [(cid, lengths.get(cid)) for cid in members[name]]
            else:
                records = [(cid, len(seq)) for cid, seq in iter_fasta(layout.bins_dir / rank / (name + '.fasta'))]
            ids = [cid for cid, _ in records]
            if len(set(ids)) != len(ids) or assigned.intersection(ids):
                error = 'Overlapping contig assignments prevent a verified partition at this rank.'
                break
            if any(cid not in lengths or size != lengths[cid] for cid, size in records):
                error = 'Bin contigs do not match the full assembly FASTA at this rank. Regenerate the final assembly export.'
                break
            bp = sum(lengths[cid] for cid in ids)
            expected = row.get('total_length_bp')
            if expected is not None and bp != expected:
                error = 'Bin summary totals do not match assembly membership at this rank. Refresh the bin summaries.'
                break
            assigned.update(ids)
            if not is_unbinned(name):
                binned.update(ids)
                if bp:
                    parts.append(dict(bin=name, bp=bp, contigs=len(ids), unbinned=False))
        if error:
            result['ranks'][rank] = dict(error=error)
            continue
        parts.sort(key=lambda p: (-p['bp'], p['bin']))
        remainder = set(lengths) - binned
        # Keep the black destination even for a completely binned assembly.
        parts.append(dict(bin='Unbinned', bp=sum(lengths[cid] for cid in remainder),
                          contigs=len(remainder), unbinned=True))
        result['ranks'][rank] = dict(parts=parts)
    return result
