# MetaHopper: coverage, junctions and endosymbiont reference trees

Keep `MetaHopper.py`, `metahopper_report.py` and `metahopper_inspect.py` together.
This update preserves the earlier conservative bin-expansion defaults and adds an
Inspection tab to the HTML report. It does not automatically trim, circularize,
merge or promote candidate sequences into accepted bins.

## Usage

For a new run, add `--references reference_genomes` to your usual MetaHopper
command if you want phylogenies. Reference trees are optional. The only new
workflow controls are:

- `--inspect endosymbionts` (default): inspect recognized endosymbiont bins and
  groups assigned through reference metadata.
- `--inspect all`: inspect all classified bins at the selected rank.
- `--inspect off`: disable inspection in the main pipeline.
- `--references PATH`: add within-group reference phylogenies.

To inspect an existing run without repeating binning or assembly:

```bash
python metahopper_inspect.py -i metahopper_output --references reference_genomes -t 24
```

Omit `--references` for coverage/junction inspection alone. The standalone script
recovers assembly and read paths from the run when available. For moved files or
older runs, supply `--assembly assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz`.
A matching whole-assembly BAM can be supplied with `--bam reads.sorted.bam`;
junction testing still needs paired FASTQ reads. `--rank genus` selects a bin rank.
The standalone script regenerates `metahopper_report.html`.

## Reference folder

Provide one nucleotide genome per FASTA file; multiple contigs per genome are
allowed. Supported suffixes are `.fa`, `.fna`, `.fasta`, `.fas`, optionally `.gz`.
For example, put Sulcia references in `reference_genomes/Sulcia/` and Vidania
references in `reference_genomes/Vidania/`. Use at least three suitable references
if comparing a single reconstructed bin: trees require at least four genomes,
including a bin and a reference. Separate taxonomic groups produce separate trees.
Do not place phylogenetically distant symbionts into a single catch-all group.

Known taxon names in paths or the first FASTA header are recognized. Otherwise,
the first subfolder supplies the group; flat accession-only files need metadata.
Optional `references.tsv` in the reference folder overrides automatic assignments:

```tsv
file	group	translation_table	label
Sulcia/GCF_example.fna	Sulcia	11	Sulcia reference A
Vidania/GCF_example.fna	Vidania	4	Vidania reference A
```

The example accessions are placeholders. Paths in `file` are relative to the
reference folder. An optional `bins.tsv` explicitly assigns reconstructed bins:

```tsv
bin	group	translation_table
my_bin_name	Vidania	4
```

Use the bin filename without `.fasta`. Verify translation tables against reference
metadata. Default hints are code 4 for Vidania, Nasuia, Zinderia and Spiroplasma,
and code 11 for the other panel groups. Only codes 4 and 11 are supported here.
These settings affect tree-specific gene calling; they do not retroactively fix
translation choices used by the upstream classification pipeline.

## What is inspected

The bacterial screening panel includes Sulcia/Karelsulcia, Vidania, Nasuia,
Zinderia, Baumannia, Arsenophonus, Sodalis, Wolbachia, Rickettsia, Cardinium and
Spiroplasma. They are candidate taxa across different host lineages, not an
expected checklist for every insect. This is a taxonomy-based screen, not a new
sensitive homology search: unknown or misclassified symbionts can still be missed.
Fungal replacements are outside this panel.

Small taxonomically labelled fragments excluded from bin-size retention are
included as inspection-only candidates. Host/gene-density excluded records are
not rescued. Candidate presence is not proof of symbiosis; lack of a taxonomic
bin is not proof of absence. Supplied reference lengths provide a size comparison,
not an estimate of genome completeness.

**Coverage:** mean and median depth, breadth at 1x/10x, coefficient of variation,
and 1 kb coverage profiles. Zero, unusually low (<20% of contig mean), and high
(>3 times mean) windows are flagged for review. Depth uses base and mapping
qualities >=20, excludes unmapped/secondary/supplementary/QC-failed/duplicate
records, and suppresses overlapping-mate double counting. Duplicate filtering
only works where input BAM records are already marked as duplicates. Whole-assembly
mapping retains competing sequences; ambiguous repeats may appear undercovered.
BAM reference names and lengths are checked, but sequence identity cannot be
proven from those fields: use a BAM mapped to the exact current assembly.

**Terminal junctions:** exact end overlaps (30 bp to the smaller of 10 kb or 20%
of contig length) are recorded. For contigs >=1 kb, synthetic end-to-start junctions
are tested against paired reads, competitively with the full assembly. Overlap
removal is confined to the test sequence. A supported junction requires >=5 read
templates with >=3 distinct alignment start/strand combinations, each continuously
spanning the join by >=25 bp on each side, MAPQ >=20, anchor base qualities >=20,
and <=2% alignment edits. Paired fragments bracketing the join are reported but
are not sufficient on their own. Repeats can prevent unique support. These checks
do not validate internal joins or prove a complete circular chromosome.

**Reference trees:** Prodigal uses an explicit translation table for each genome.
DIAMOND reciprocal-best-hit proteins are anchored to the reference with the most
predicted proteins. Matches require >=30% identity and >=60% bilateral coverage;
near-tied alternative hits are excluded. Markers need >=80% taxon occupancy,
with >=70% marker occupancy per genome. MAFFT alignments retain columns with >=80%
canonical amino acids. At least 10 markers, 1,000 retained sites, and >=70%
observed amino acids per genome are required. IQ-TREE uses LG+F+G4 and 1,000
ultrafast bootstraps with `-bnni`. Outputs include Newick, SVG, alignments, genome
metadata and marker tables. Missing tools or insufficient evidence produce an
explicit status rather than a fabricated tree.

Trees are exploratory, unrooted, within-group protein trees—not a curated universal
marker analysis or species assignment. Genome reduction, composition bias,
paralogy and strain mixtures can distort placement. Short candidates may fail
Prodigal single-genome training. The display orientation does not establish a root.

## Runtime and dependencies

Coverage reuses the pipeline BAM when available; otherwise it maps the reads once.
Junction inspection adds one competitive mapping of the paired read library, which
can be substantial for large datasets. It does not reassemble reads. Reference
trees run only when a reference folder is supplied; successful trees are cached
against input content, metadata and tool binary path/mtime. Coverage and junction
mapping are currently rerun on each inspection invocation. No speedup or biological
quality gain has been measured on your data.

External inspection dependencies: Bowtie2, samtools, Prodigal, DIAMOND, MAFFT and
IQ-TREE (`iqtree2` or `iqtree`). MAFFT/IQ-TREE are only needed for reference trees.
The Python inspection module itself uses the standard library. Existing pipeline
dependencies remain necessary for full pipeline runs.

Outputs live under `inspection/`: `inspection.json`, `endosymbionts.tsv`,
`coverage.tsv`, `coverage_windows.tsv`, `junctions.tsv`, `tree_summary.tsv`,
mapping logs/BAMs, and per-group `phylogeny/` files. The report embeds profiles and
trees for review. Missing results and failures are reported in the JSON and HTML.

## Relation to your supplied workflows

This adds the read-depth and end-overlap review emphasized by the Anna and
Verdanus workflows, while making evidence explicit and keeping original sequences
intact. Translation-table choices are explicit in the new tree stage. Manual
self-overlap trimming remains a reviewed decision. Core assessment still uses the
existing CheckM workflow; this update does not substitute CheckM2 or claim that
marker completeness is reliable for every reduced endosymbiont genome.

Validation: all 27 unit and mocked integration checks passed, including coverage
zero filling, junction CIGAR/quality filtering, reference grouping, paralog rejection,
report rendering and a mocked complete tree workflow. Actual aligners/assemblers
were not installed in the validation environment; no real-data biological benchmark
or end-to-end external-tool execution was performed.

Tool documentation: [samtools depth](https://www.htslib.org/doc/samtools-depth.html),
[DIAMOND options](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options),
[MAFFT](https://mafft.cbrc.jp/alignment/software/manual/manual.html), and
[IQ-TREE](https://iqtree.github.io/doc/Command-Reference).
