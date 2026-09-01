# MetaHopper

**Metagenome assembly, microbial-contig retention, taxonomic binning, and targeted extension**

MetaHopper accepts paired reads, assembled contigs, or both. It builds or uses a complete metagenome assembly, predicts proteins with Prodigal, and classifies contigs against the taxonomy-enabled `nr-tax.dmnd` DIAMOND database. Animal and plant contigs are separated by default, while microbial contigs are retained and organized into domain-, phylum-, family-, genus-, and species-level bins.

When reads are available, MetaHopper automatically uses the classified bins as seeds for targeted read recruitment. Reads map competitively across all bins, ambiguous ties are excluded, and guarded exact-k-mer extension attempts to walk beyond the original contigs. Recruited reads are assembled with SPAdes/metaSPAdes, passed through Unicycler for graph resolution and possible circularization, and polished with Pilon. The refined assembly is then completely reclassified and assessed with QUAST and CheckM.

```mermaid
flowchart LR
    A[Reads, contigs, or both] --> B[Full assembly]
    B --> C[Prodigal and nr-tax.dmnd]
    C --> D[Microbial taxonomic bins]
    D --> E[Competitive seed extension]
    E --> F[SPAdes, Unicycler, and Pilon]
    F --> G[Final bins, QUAST, and CheckM]
```

Seed-and-extension is enabled by default whenever paired reads are supplied. Contigs-only runs stop after classification, binning, and quality assessment because no reads are available for extension.

## Input behavior

| Input | Behavior |
|---|---|
| Paired reads only | QC → MEGAHIT assembly → classification → seed extension → focused reassembly → final classification and assessment |
| Contigs plus paired reads | Uses the supplied contigs, skips MEGAHIT, and uses the reads for seed extension and refinement |
| Contigs only | Classifies, retains, bins, and assesses the supplied assembly; no extension or circularization |

## Requirements

MetaHopper uses Python 3, Prodigal, DIAMOND ≥2.1.17, MEGAHIT, Bowtie2, Samtools, BBMap/BBDuk, SPAdes, Unicycler, Pilon, Trimmomatic, FLASH2, QUAST, and CheckM. fastp is optional for poly-G trimming.

Create and activate an environment:

```bash
mamba create -n metahopper \
    -c conda-forge \
    -c bioconda \
    python prodigal diamond megahit bowtie2 samtools bbmap spades \
    unicycler pilon trimmomatic flash2 fastp pigz quast checkm-genome \
    --yes

mamba activate metahopper
chmod +x MetaHopper.v2.py
```

## Build `nr-tax.dmnd`

MetaHopper requires a taxonomy-enabled NCBI nr database. Download a compatible set of:

- `nr.gz`
- `prot.accession2taxid.FULL.gz`
- `nodes.dmp`
- `names.dmp`

Build the database:

```bash
diamond makedb \
    --in nr.gz \
    --db nr-tax \
    --taxonmap prot.accession2taxid.FULL.gz \
    --taxonnodes nodes.dmp \
    --taxonnames names.dmp \
    --threads 24
```

This creates `nr-tax.dmnd`. The accession-to-taxid mapping, hierarchical NCBI taxonomy, and scientific names are baked into the database, so MetaHopper does not require a separate taxdump directory or `ete3` at runtime.

## Quick start

### Paired reads only

```bash
./MetaHopper.v2.py \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    --trimmomatic-folder /path/to/Trimmomatic-0.39 \
    -d /path/to/nr-tax.dmnd \
    -o metahopper_sample \
    --ranks domain,phylum,family,genus,species \
    -t 24
```

### Existing contigs plus paired reads

```bash
./MetaHopper.v2.py \
    -i assembly.fasta \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    --trimmomatic-folder /path/to/Trimmomatic-0.39 \
    -d /path/to/nr-tax.dmnd \
    -o metahopper_refined \
    --ranks domain,phylum,family,genus,species \
    -t 24
```

### Existing contigs only

```bash
./MetaHopper.v2.py \
    -i assembly.fasta \
    -d /path/to/nr-tax.dmnd \
    -o metahopper_contigs \
    --ranks domain,phylum,family,genus,species \
    -t 24
```

## Important behavior

- Metazoa and Viridiplantae contigs are excluded from microbial bins by default and written to a separate FASTA.
- DIAMOND protein hits are combined across all ORFs using a rank-specific, bitscore-weighted consensus.
- Competitive seed mapping assigns each template to one uniquely best-scoring bin; equal-score ties are excluded.
- BBDuk extension uses newly recruited reads to walk outward, with growth safeguards to prevent runaway recruitment.
- Unicycler attempts graph resolution and circularization; Pilon polishes the resulting sequence.
- Failed or skipped refinements fall back to the original preliminary-bin contigs.
- Refined contigs receive a new Prodigal/DIAMOND classification before final bins are created.
- `--skip-reassembly` disables read recruitment and the second classification pass.

## Primary outputs

```text
<output>/
├── classification/contig_classification.tsv
├── bins/<rank>/<taxon>.fasta
├── reassembly/<rank>/competitive_seed/competitive_mapping.tsv
├── reassembly/<rank>/<bin>/reassembled.fasta
├── final/assembly/consolidated_contigs.fasta
├── final/assembly/contig_provenance.tsv
├── final/classification/contig_classification.tsv
└── final/bins/<rank>/summary.tsv
```

The `final/` directory is produced when seed-and-extension runs. In contigs-only mode, `bins/<rank>/` contains the assessed final bins.

Run `./MetaHopper.v2.py --help` for all classification, recruitment, assembly, and filtering options.
