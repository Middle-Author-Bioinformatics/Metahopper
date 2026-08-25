# MetaHopper

End-to-end pipeline for taxonomic classification and binning of metagenomic contigs, starting from either raw paired-end FASTQs or an already-assembled contigs FASTA.

```
R1.fastq, R2.fastq  (raw paired-end reads; optional -- skip straight to contigs.fasta)
    -> [optional] fastp poly-G trimming (--trim-polyg)
    -> QC: Trimmomatic (adapter clip + quality trim) + FLASH (merge overlapping pairs)
    -> MEGAHIT assembly
contigs.fasta
    -> Prodigal (ORF/gene prediction, meta mode)
    -> DIAMOND blastp vs. a DIAMOND protein database (nr.dmnd or similar)
    -> per-ORF top hits -> per-contig taxonomic classification
    -> bin FASTA files at genus / species level (or whatever --ranks you ask for)
    -> per-bin-set summary.tsv (QUAST contiguity + CheckM completeness/contamination)
```

## Why

Classifying a *contig* from many *ORF* hits is the same problem CAT/BAT, MEGAN's LCA, and Kraken-style consensus callers solve. MetaHopper's approach:

1. For each ORF, keep the DIAMOND hits within `--bitscore-range` (default 90%) of that ORF's best bitscore -- a "bit-score competitive set," not just the single top hit.
2. Each kept hit's organism name is parsed straight out of its `stitle` (the trailing `[Genus species]`, standard NCBI nr header format), resolved to a full NCBI lineage, and contributes its bitscore as a vote for every taxon in that lineage.
3. Votes are tallied independently at each requested rank. A taxon wins a rank only if its share of the total weighted vote clears `--min-support` (default 0.5, i.e. a majority; set to `0` for a pure plurality/"most votes wins" call). Otherwise the contig is `Unclassified` at that rank.

Because organism names come from `stitle`, **no `--taxonmap`/`--taxonnodes` embedding is required in the DIAMOND database** -- any `nr.dmnd` built with a plain `diamond makedb --in nr.faa` works.

## Usage

From raw paired-end reads:

```bash
MetaHopper.py --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz \
  --trimmomatic-folder /path/to/Trimmomatic-0.39 \
  -d /path/to/nr.dmnd -o metahop_out -t 16 \
  --taxdump-dir /path/to/taxdump \
  --ranks domain,phylum,family,genus,species
```

From an existing assembly:

```bash
MetaHopper.py -i contigs.fasta -d /path/to/nr.dmnd -o metahop_out -t 16 \
  --taxdump-dir /path/to/taxdump
```

Run `MetaHopper.py --help` for the full option list (QC/MEGAHIT tuning, `--min-support`, `--min-bin-contigs`/`--min-bin-length` to collapse small/noisy bins, `--reuse-prodigal`/`--reuse-diamond` to re-classify at a different rank without re-running the expensive steps, etc.).

## Requirements

**Only if starting from raw reads (`--r1`/`--r2`):**
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (skip with `--skip-qc`)
- [FLASH](https://ccb.jhu.edu/software/FLASH/) (skip with `--skip-qc`)
- [pigz](https://zlib.net/pigz/) (falls back to `gzip` if missing)
- [fastp](https://github.com/OpenGene/fastp) (only with `--trim-polyg`)
- [MEGAHIT](https://github.com/voutcn/megahit)

**Always:**
- [Prodigal](https://github.com/hyattpd/Prodigal)
- [DIAMOND](https://github.com/bbuchfink/diamond)
- [QUAST](https://github.com/ablab/quast) (optional; falls back to a built-in N50/L50/GC calculator if missing or `--skip-quast`)
- [CheckM](https://github.com/Ecogenomics/CheckM) `lineage_wf` (optional; skipped with a warning if missing or `--skip-checkm`)

**Taxonomy lookups** -- pick one:
- `--taxdump-dir <dir>`: a local NCBI taxdump (`nodes.dmp` + `names.dmp`, e.g. an extracted [`taxdump.tar.gz`](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)). Fast, offline, no extra Python dependency. **Recommended.**
- Otherwise falls back to `ete3.NCBITaxa`, which downloads and builds its own local sqlite taxonomy db on first use (`pip install ete3 six`) -- requires outbound internet access.

## Files

- `MetaHopper.py` -- the pipeline
- `full_lineage.py` -- standalone helper: expand an NCBI taxid to its complete ranked lineage (domain through species/strain) from a local taxdump, e.g. for inspecting `staxids` pulled from a taxonomy-embedded DIAMOND database

## License

Add a license of your choice (e.g. MIT) before publishing.
