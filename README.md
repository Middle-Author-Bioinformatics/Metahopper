# MetaHopper

End-to-end pipeline for taxonomic classification and binning of metagenomic contigs, starting from either raw paired-end FASTQs or an already-assembled contigs FASTA -- with an optional targeted reassembly module for pulling fragmented, low-abundance genomes (e.g. host-restricted symbionts) back together.

```
R1.fastq, R2.fastq  (raw paired-end reads; optional -- skip straight to contigs.fasta)
    -> [optional] fastp poly-G trimming (--trim-polyg)
    -> QC: Trimmomatic (adapter clip + quality trim) + FLASH (merge overlapping pairs)
    -> MEGAHIT assembly
contigs.fasta
    -> Prodigal (ORF/gene prediction, meta mode)
    -> DIAMOND blastp vs. a DIAMOND protein database (nr.dmnd or similar)
    -> per-ORF top hits -> per-contig taxonomic classification
    -> drop host-animal/plant contamination (Metazoa/Viridiplantae, configurable)
    -> bin FASTA files at genus / species level (or whatever --ranks you ask for)
    -> per-bin-set summary.tsv (QUAST contiguity + CheckM completeness/contamination)
    -> [optional] targeted per-bin reassembly (--reassemble-bins): recruit raw reads
       back onto each bin's own contigs, extend recruitment with exact-kmer frontier
       scans, reassemble the recruited pool alone with SPAdes/metaSPAdes, report a
       before/after contiguity comparison
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
MetaHopper.v2.py --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz \
  --trimmomatic-folder /path/to/Trimmomatic-0.39 \
  -d /path/to/nr.dmnd -o metahop_out -t 16 \
  --taxdump-dir /path/to/taxdump \
  --ranks domain,phylum,family,genus,species
```

From an existing assembly:

```bash
MetaHopper.v2.py -i contigs.fasta -d /path/to/nr.dmnd -o metahop_out -t 16 \
  --taxdump-dir /path/to/taxdump
```

With targeted bin reassembly, to pull fragmented low-abundance/symbiont genomes back together after the initial binning pass:

```bash
MetaHopper.v2.py --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz \
  --trimmomatic-folder /path/to/Trimmomatic-0.39 \
  -d /path/to/nr.dmnd -o metahop_out -t 16 \
  --taxdump-dir /path/to/taxdump \
  --ranks genus,species \
  --reassemble-bins --reassemble-ranks genus,species
```

Run `MetaHopper.v2.py --help` for the full option list (QC/MEGAHIT tuning, `--min-support`, `--min-bin-contigs`/`--min-bin-length` to collapse small/noisy bins, `--reuse-prodigal`/`--reuse-diamond` to re-classify at a different rank without re-running the expensive steps, the full `--reassemble-*` group, etc.).

## Targeted bin reassembly (`--reassemble-bins`)

Whole-metagenome co-assembly (MEGAHIT) tends to fragment low-abundance or fast-diverging genomes -- especially host-restricted symbionts (*Wolbachia*, *Rickettsia*, etc.) that share conserved genes/k-mers with the rest of the community and rarely have a close enough reference genome to assemble against directly.

Instead of requiring a reference, this module uses each bin's **own already-classified contigs** as the seed:

1. Raw reads are mapped back onto that bin's contigs with `bowtie2 --very-sensitive-local`.
2. Confidently-mapping read pairs (+ their mates, however the mate aligned) are recruited.
3. Recruitment is optionally extended with cheap exact-kmer "frontier" scans (BBDuk) to catch divergent regions the initial contigs missed entirely -- guarded by growth-rate and accepted-fraction stop conditions so it can't snowball into an unrelated, similar-coverage genome.
4. The resulting small, mostly-single-organism read pool is QC'd and reassembled alone with SPAdes/metaSPAdes.

A subset-only assembly graph is far simpler than the whole-community graph, so it can often resolve tangles (shared genes, similar-coverage strains) that fragmented the original bin. An optional `--anchor-db` (repeatable) lets you add external reference sequences -- a related genome, conserved marker genes -- as supplementary seeds, useful for recruiting reads from genome regions the current fragmented bin has zero contigs for at all.

This **does not overwrite** the original bin FASTA. It writes a separate reassembled contig set plus a before/after comparison so you can judge whether it actually helped:

```
<outdir>/reassembly/<rank>/<bin>/reassembled.fasta
<outdir>/reassembly/<rank>/<bin>/recruitment.tsv       (per-round recruitment metrics)
<outdir>/bins/<rank>/reassembly_summary.tsv            (contigs/N50/length, before vs. after)
```

**Caveat:** this maps the *entire* raw read set against each bin's contigs, once per bin, so with many bins across multiple ranks it can get slow and disk-heavy. Use `--reassemble-ranks` and `--reassemble-min-bin-contigs` to narrow scope to the bins that actually look fragmented.

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

**Only with `--reassemble-bins`:**
- [bowtie2](https://github.com/BenLangmead/bowtie2) (`bowtie2` + `bowtie2-build`)
- [samtools](https://github.com/samtools/samtools) (needs `samtools view -N` support)
- [BBMap](https://sourceforge.net/projects/bbmap/) (`bbduk.sh`)
- [SPAdes](https://github.com/ablab/spades) (`spades.py`)
- Trimmomatic + FLASH (same as above -- required even if `--skip-qc` was used for the main assembly, since recruited reads always get QC'd)

**Taxonomy lookups** -- pick one:
- `--taxdump-dir <dir>`: a local NCBI taxdump (`nodes.dmp` + `names.dmp`, e.g. an extracted [`taxdump.tar.gz`](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)). Fast, offline, no extra Python dependency. **Recommended.**
- Otherwise falls back to `ete3.NCBITaxa`, which downloads and builds its own local sqlite taxonomy db on first use (`pip install ete3 six`) -- requires outbound internet access.

## Files

- `MetaHopper.v2.py` -- the current pipeline (raw reads or contigs in, classification + binning out, optional targeted bin reassembly)
- `MetaHopper.py` -- earlier version, kept for reference
- `MetaHopper_contigs_only.py` -- stripped-down variant that only accepts an already-assembled `contigs.fasta` (no raw-read QC/MEGAHIT support)
- `full_lineage.py` -- standalone helper: expand an NCBI taxid to its complete ranked lineage (domain through species/strain) from a local taxdump, e.g. for inspecting `staxids` pulled from a taxonomy-embedded DIAMOND database

## License

Add a license of your choice (e.g. MIT) before publishing.
