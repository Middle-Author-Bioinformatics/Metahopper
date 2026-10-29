# MetaHopper Lean 2.2

MetaHopper assembles or accepts a short-read metagenome assembly alongside reads,
classifies contigs, writes taxonomic partitions, and attempts guarded Unicycler
reassembly at one rank. The HTML report supports manual bin curation, source
assembly-graph inspection and optional per-bin Bakta annotation summaries.

Taxonomic bins can contain several genomes or strains. Neither a bin label nor an
annotation establishes that a bin is a pure, complete genome.

## Install and keep the modules together

The active files are `metahopper.py`, `core.py`, `taxbin.py`,
`metahopper_report.py`, `metahopper_inspect.py`, `metahopper_graph.py`,
`metahopper_annotate.py`, and `metahopper_state.py`. Copy all eight into your MetaHopper installation. Do not
replace your working `config.json` with the example configuration.

The main pipeline requires Python with numpy, scikit-learn and pysam, plus fastp,
Prodigal, DIAMOND and BWA. MEGAHIT is required unless an assembly is supplied.
Default stages also use Tiara, phyloFlash and its database, Unicycler and minimap2.
Inspection uses Bowtie2 and samtools. Reference trees additionally use MAFFT and
IQ-TREE. CheckM and Bakta are optional stages with their own prepared databases.
The report, graph catalogue and annotation orchestration modules use the Python
standard library; external Bakta still requires its own installed dependencies.

## Main workflow

```bash
python metahopper.py \
  -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  -c config.json -o sample_metahopper -t 24 --memory 128
```

Supply `--assembly assembly.fasta` to skip MEGAHIT. Reads are still required by the
main workflow. `--skip-reassembly` skips Unicycler but still produces `final/`.
`--skip-phyloflash`, `--prefilter none`, and `--inspect off` disable the corresponding
optional stages. `--checkm` enables final-bin quality assessment.

The default partitions are phylum, class, order, family, genus and species. Bins
without sufficient rank resolution receive higher-rank fallback names; small
bins merge upward. Reassembly uses genus by default. The candidate must retain
at least 95% of original sequence at 95% identity, have 95-125% of original total
length, and reduce contig count by at least 10% or increase N50 by at least 20%.
Otherwise original contigs are kept. Final taxonomy is inherited rather than
searched again with DIAMOND.

## Existing output directories and completion logs

Repeating a command with an existing `-o RUN` now automatically checks completion
records and reuses successful stages with matching inputs, settings and outputs.
You can also restore a run's saved arguments with:

```bash
python metahopper.py --resume RUN
python metahopper.py --resume RUN --annotate bakta --annotation-rank genus -t 24
```

Explicit options supplied with `--resume` now override saved arguments. A flag
such as `--checkm` can therefore add an optional stage without rerunning completed
assembly or classification. Use the standalone report command below when only
the HTML needs refreshing.

The parent output directory contains two small status files:

| File | Contents |
| --- | --- |
| `pipeline.log` | Append-only text timeline: step, time, running/completed/reused/failed/skipped state |
| `pipeline_status.json` | Current step, overall run status and the latest state of each step |

For example, `tail -n 20 RUN/pipeline.log` shows recent progress while a run is
active. Detailed tool logs remain in `logs/`; successful output metadata and
settings are retained under `work/checkpoints/`.

Reuse covers QC, assembly, Prodigal, Tiara, DIAMOND, taxonomy voting, coverage/read
pool mapping, preliminary bins/ordination, each reassembly candidate, consolidated
final-bin work, CheckM, inspection and Bakta's per-bin cache. Failed or interrupted
steps retry. A failed reassembly preserves its recruited reads for retry; pools
already removed after completed candidates do not trigger remapping. Missing
pools needed by unfinished candidates are regenerated. The HTML is always rebuilt.

Older runs can reuse compatible success checkpoints. A saved completed
`run.json` plus the required outputs can also verify older binning/final-bin and
CheckM stages that did not have individual checkpoints. Plain output directories
and free-text tool logs alone are not treated as proof of success. Stages without
a verifiable completion record rerun once and receive the new records.

To recompute all **enabled** stages, bypassing completion and per-bin caches:

```bash
python metahopper.py --resume RUN --force-overwrite
```

The flag applies only to that invocation. Input checks use paths, sizes and
nanosecond modification times, plus relevant settings and database/version
metadata; large FASTQs/databases are not hashed in full. Bakta additionally hashes
bin FASTAs. If files or software were edited in place while preserving those
identifiers, use `--force-overwrite`.

## Report navigation

The report has exactly three tabs: **Bins**, **BinaRena**, and **Annotation**.
Bins opens first with classified-bin, assembly-size and plotted-contig cards above
the table. The rank selector starts at Phylum and offers one rank at a time;
there is no combined all-ranks view. If Phylum was not produced, the highest
available rank opens instead. Backend inspection and detailed audit files remain
available in the run directory.
Use **Search taxonomic name** above the Bins table to find matching bin names
across every available rank. Matching ignores capitalization and treats spaces
and underscores alike; partial names work. Each matching row retains its rank,
size, contig statistics, CheckM scores when available, and refinement summary.
The rank selector pauses while a search is active. **Clear**, Escape, or deleting
the query restores the previous rank selection. No matches are shown explicitly.
The classified-bin card retains the total across available ranks; these groups
are nested, so it is not a count of distinct genomes. Annotation rows are marked
stale if their source bin changed or disappeared since annotation.

## Fix colors and add graph review to an existing run

Replace the scripts, then regenerate the report on the server containing the
run's actual files:

```bash
python /path/to/Metahopper/metahopper_report.py \
  -i /path/to/APHALN_SEG35A-5
```

No assembly or DIAMOND rerun is needed. The updated report preserves the plotted
data, uses the same opaque colors for dots and legend entries, and offers every
category in an expandable legend. Colors remain stable when filtering. Manual
assignments no longer silently replace taxonomy colors: select **Manual sub-bin**
in **Colour by** to view those assignments.

Lasso contigs in BinaRena and click **View assembly graph**. The dialog displays
selected segments and zero, one or two links of neighbouring context. It supports
pan, zoom, node dragging, segment details and SVG export. Grey nodes are context;
hatched nodes belong to selected contigs with conflicting colors. Links retain
GFA orientation, and jump records are dashed. Layout positions do not represent
genomic distances.

### Graph sources and mappings

Graph generation uses GFA1 `S`, `L`, `J`, `P` and `W` records. Unsupported record
types are counted and reported. GFA2 is rejected. Source files are kept separate;
segment `1` in Sodalis is never confused with segment `1` in another source graph.

For final reassembled contigs, the report reads
`final/assembly/contig_provenance.tsv`, the source bin's `assembly.fasta`, and
`reassembly/<rank>/<bin>/assembly.gfa`. New runs also preserve
`contig_graph_aliases.tsv`, recording original Unicycler identifiers before
renaming. Older runs can be mapped by a unique exact segment sequence match on
either strand. Duplicate sequence matches remain ambiguous; contig-number order
is never used as evidence of a mapping.

Original fallback contigs need their preliminary graph. A GFA from a rejected
candidate is not used as the graph of the original final sequence. For supplied
assemblies, new runs can retain a matching graph with:

```bash
python metahopper.py --assembly assembly.fasta --assembly-graph assembly.gfa \
  -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -c config.json -o sample_run
```

For an existing run, pass `--graph preliminary.gfa` to the report generator.
Missing graph files, unsupported mappings and ambiguous matches are reported
explicitly. A FASTA alone cannot supply assembly links. Exact sequence recovery
does not handle rotated circular sequences or reconstruction from unnamed graph
paths; those cases need preserved identifiers or an explicit named path.

The graph catalogue contains topology and metadata, not nucleotide sequences.
Graphs are embedded automatically up to 200,000 segments per source and 32 MB of
combined JSON. The interactive view defaults to 500 nodes and refuses to drop
selected segments silently. A context neighbourhood may be explicitly truncated.
Larger or separately transferred catalogues can be loaded locally in the dialog:

```bash
python metahopper_graph.py -i /path/to/RUN --bin-set final
# Writes RUN/assembly_graphs.final.json
```

Use a catalogue from the same run and bin set. The viewer checks run names and
available contig lengths. `--skip-graphs` omits automatic graph embedding, and
`--graph-max-nodes` adjusts the per-source import limit.

## Compact refinement and CheckM results

The Bins table includes a short refinement summary with the recorded outcome,
contigs before/after, and absolute and percent size change. It uses displayed
final-bin metrics rather than substituting candidate-assembly statistics.
Provenance identifies bins containing accepted reassembled sequence, including
bins renamed or regrouped at the final stage. The detailed candidate audit files
remain on disk; there is no separate Expansion tab.

When available, CheckM completeness, contamination and strain heterogeneity are
shown in the Bins table. Scores are read from each rank's summary and can fall
back to the selected bin set's checkm.tsv using rank-qualified bin IDs. Missing
scores remain NA. A missing contamination score is not treated as a zero or a
quality pass. New runs still require --checkm to perform CheckM assessment.

## Optional Bakta annotation

Bakta targets bacterial genomes and plasmids. This module processes bacterial bins
at one rank, skips literal Unclassified and domain-only catch-all bins, and skips
bins with conflicting or non-bacterial domain calls. Unresolved contigs may remain
within a bacterial taxonomic group; inspect the recorded bacterial fraction and
bin composition before interpreting functions as one organism's capabilities.

Add the standard prepared Bakta database to your existing configuration:

```json
{
  "bakta_db": "/data/databases/bakta/db",
  "bakta_args": ["--skip-plot"]
}
```

These are additions to the complete existing JSON object. They do not replace
`diamond_db`, `taxdump` or your other settings. To invoke Bakta in a separate
conda environment, merge this entry into `tools`:

```json
"bakta": ["conda", "run", "--no-capture-output", "-n", "bakta_env", "bakta"]
```

### Annotate an existing run

```bash
python /path/to/Metahopper/metahopper_annotate.py \
  -i /path/to/APHALN_SEG35A-5 -c /path/to/config.json \
  --rank genus --bin-set final -t 24
```

This requires only Bakta and its prepared database, not a repeat of the upstream
analysis. It also refreshes the HTML report. `--bins Sodalis,Symbiopectobacterium`
limits the requested bins; valid existing annotations of other bins are retained.
Changed unselected inputs are labelled stale. `--force` reannotates requested
bins regardless of the successful cache. A six-hour per-bin timeout is the
default; `--timeout-h 0` disables it.

### Annotate during a new analysis

```bash
python metahopper.py \
  -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  -c config.json -o sample_metahopper -t 24 --memory 128 \
  --annotate bakta --annotation-rank genus
```

Only one rank is annotated, avoiding duplicate annotation across nested taxonomic
partitions. Bins run sequentially using the requested threads. A failed bin is
recorded and later bins continue. The standalone command exits nonzero if a
requested annotation fails; the main pipeline records annotation warnings.

Bakta uses metagenome CDS prediction by default (`--meta`), preserves contig
headers, and does not mark all contigs complete or circular. Use standalone
`--mode single` or main `--annotation-mode single` for single-genome training.
The default translation table is 11. Vidania, Nasuia, Zinderia and Spiroplasma
labels provide code-4 hints. Conflicting code hints require explicit metadata.
Verify code choices for your genomes; taxon hints do not prove the genetic code.

For per-bin overrides, add `"annotation_metadata": "annotation_metadata.tsv"`
to the config. Relative database and metadata paths resolve from that config.

```tsv
rank	bin	translation_table
genus	Spiroplasma	4
genus	Sodalis	11
```

Create this as a tab-delimited file.
The bin value is the filename stem, without `.fasta`. An optional `domain`
column can explicitly establish a bin as `Bacteria` when its taxonomy is unresolved.
Codes 4, 11 and 25 are supported. The standalone `--translation-table` or main
`--annotation-translation-table` sets the general default. Managed parameters
cannot be replaced through `bakta_args`; other Bakta options can be supplied as
an array. The standard Bakta database is required now; additional specialist
databases are not required by this integration.

### Outputs and later modules

For final bins, outputs are under `final/annotation/<rank>/`:

| Output | Contents |
| --- | --- |
| `summary.tsv` | Per-bin status, reason, genetic code, sequence size and feature counts |
| `annotation_manifest.json` | Tool/database identity, bin results and module status |
| `<bin>/bakta/annotation.*` | Raw Bakta GenBank, GFF3, proteins, feature DNA, JSON, TSV and other outputs |
| `<bin>/features.tsv` | Normalized contig/locus IDs, feature type, 1-based inclusive coordinates, strand, gene, product and db cross-references |
| `<bin>/result.json` | Input/configuration fingerprint and completed/failed/skipped state |
| `<bin>/modules.json` | Per-bin entry points and status for downstream modules |

The report's **Annotation** tab summarizes actual completed results and links to
files when the full run folder is available. Download links will not work if
only the HTML has been uploaded to a report hosting service.

The cache considers input bytes, Bakta version, command prefix, database version
metadata, selected code, mode and extra options. Database content is not hashed
in full; use standalone `--force` (or main-pipeline `--force-overwrite`) if a database was changed in place without updating its
version metadata. Previous successful outputs remain available until a new
annotation completes, while a failed rerun is marked failed rather than current.

Pathway reconstruction, BGC detection and a dedicated virulence screen remain
`not_run`, with no fabricated zero counts. Future modules can consume the same
`.faa`, `.gbff`, coordinates and locus IDs through `modules.json`. Bakta's own
expert annotations are retained in its raw output; they are separate from a
future dedicated screening workflow.

## Validation

```bash
python -m unittest discover -p 'test_metahopper*.py' -v
node test_report_frontend.cjs /path/to/metahopper_report_updated.html
```

The Python suite exercises graph parsing, provenance, exact/reverse-complement
matching, ambiguous and missing graphs, code handling, caching and per-bin failure
isolation with a fake Bakta executable. JavaScript tests exercise report logic
with a DOM/canvas harness and the supplied 13,375-contig report. A simulated full
pipeline exercises automatic reuse, missing outputs, failed-step retries, legacy
checkpoints and force overwrite without executing real aligners/assemblers. These checks do
not constitute execution of real Bakta, biological validation, or browser visual QA.

Primary specifications: [GFA1](https://gfa-spec.github.io/GFA-spec/GFA1.html) and
[Bakta](https://github.com/oschwengers/bakta).
