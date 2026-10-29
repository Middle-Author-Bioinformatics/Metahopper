# MetaHopper Lean 2.2.1

- Added a taxonomic-name search above the Bins table. Nonempty searches find matching bin names across ranks and retain each row's rank, size, CheckM scores when available, and refinement information. Matching ignores case and normalizes spaces/underscores. Clear or Escape restores the previously selected rank.

# MetaHopper Lean 2.2.0

- Report navigation now contains only Bins, BinaRena and Annotation. The three requested summary cards appear above Bins, which opens at Phylum with no combined all-ranks option.
- Existing output directories automatically reuse verified successful stages. Added parent-directory pipeline.log and pipeline_status.json, input/settings/output checks, compatibility with older completion proofs, and --force-overwrite to bypass enabled-stage caches.
- Cached final-bin work, taxonomy voting, ordination, CheckM and inspection. Failed reassemblies retain their read pools and retry; already completed candidates do not remap just because their pools were deliberately cleaned up.
- --resume now merges explicit CLI options with stored arguments and current defaults. Force overwrite is never inherited by the next invocation.

# MetaHopper Lean 2.1.0

- Expansion outcomes now appear as compact Bins-table summaries with actual before/final contig and size changes. The separate Expansion tab is removed. CheckM scores remain visible when present, including fallback to rank-qualified raw CheckM records. Missing contamination no longer passes the Overview quality badge.

- One shared, opaque colour mapping for scatter points, graph nodes and the legend. The legend can show every category. Manual assignments use an explicit Manual sub-bin colour mode.
- Lasso-selected contigs can open a GFA1 graph neighbourhood, with source isolation, orientation-aware links, pan/zoom, node details and SVG export.
- Report generation discovers source graphs and resolves final contigs through provenance, preserved assembler aliases or unique exact sequence matches. Unavailable or ambiguous mappings remain explicit.
- New runs preserve Unicycler contig-to-graph identifiers. --assembly-graph retains an optional preliminary GFA.
- Optional per-bin Bakta annotation at one rank, with a standalone existing-run command, genetic-code metadata, content/configuration caching, per-bin failure handling and raw output preservation.
- Annotation summary tab and normalized features.tsv / modules.json outputs for later specialist modules. Dedicated pathway, BGC and virulence-screening modules remain not_run.
- Preliminary reports now select the preliminary FASTA for manual extraction, even when final/ is present.

# MetaHopper Lean 2.0.0

## What 2.0.0 restores from 7.4
- **Classification.** It uses 7.4's bitscore-weighted vote for each rank: hits within 0.9 of the best bitscore, at most 5 hits per ORF, minimum support 0.5, and an `all` or `named` denominator. On 400 random test contigs it gave the same result as `classify_contig` in 7.4. Hit lineages now come from DIAMOND `staxids` plus the NCBI taxdump, so class and order are available.
- **DIAMOND settings.** These match 7.4: `--evalue 1e-5 --max-target-seqs 25`, default sensitivity. You can override them with `diamond_args` in the config.
- **Kingdom exclusion.** Contigs voted Metazoa or Viridiplantae at kingdom level are excluded, as in 7.4. They go to `classification/excluded_animal_plant_contamination.fasta`.
- **Output layout.** The layout is 7.4's: `classification/`, `bins/<rank>/` (plain 70-column FASTA, `bin_membership.tsv`, `summary.tsv`, `reassembly_summary.tsv`), `binarena/`, `reassembly/`, `final/assembly/`, `final/classification/`, `final/bins/`, `final/binarena/` and `inspection/`.
- **Report and inspection.** `metahopper_report.py` and `metahopper_inspect.py` are byte-identical to the 7.4 copies. Inspection runs by default with `--inspect endosymbionts`, as in 7.4.
- **BinaRena.** The BinaRena table has the same columns as 7.4. It includes 4-, 5- and 6-mer PCA (with t-SNE and UMAP as options), GC, coverage, covered fraction and coding density. K-mer counting now uses numpy and gives the same results as 7.4's `count_kmers`.

## New in 2.0.0
- **Fallback names.** A contig with no call at the bin rank is named after its nearest called higher rank: "Enterobacteriaceae bacterium", "Sodalis sp.", "unclassified Bacteria". Only contigs with no call at any rank go to `Unclassified`.
- **Small bins merge upward.** A bin under `--min-bin-length` (50 kb) is renamed after the next called higher rank, so a small Escherichia bin joins "Enterobacteriaceae bacterium". This repeats until every bin is big enough. Anything still too small ends in one "unclassified \<domain\>" bin per domain.
- **Ranks.** The default ranks are `phylum,class,order,family,genus,species`; set them with `--ranks`.
- **No composition splitting.** Bins are not split or recruited by coverage, GC or k-mers. Composition is exported to BinaRena for manual splitting.
- **Plain FASTA bins.** `Unclassified.fasta` is one shared file, linked from every rank folder.
- **One mapping pass.** A single stringent BWA pass gives per-contig depth, covered fraction and read pools, with every contig present as a competitor.
- **Reassembly.**
  - It runs at one rank (`--reassemble-rank genus`) and skips domain-only catch-all bins.
  - The 1.0.2 limits are kept: 20 Mbp bin cap, downsampling to 200×, 4 h timeout, Bacteria/Archaea only.
  - Unicycler output must pass the 7.4 retention check: at least 95% of the original sequence aligned at 95% identity, length between 95% and 125% of the original, and fewer contigs or a higher N50. Otherwise the original contigs are kept.
- **Final taxonomy.** It is inherited from the preliminary calls (7.4 bin consensus), with no second DIAMOND run. Reads are remapped only when a bin was accepted.
- **ORF cap.** At most `--max-orfs-per-contig` ORFs per contig (default 10, or `all`) go to DIAMOND. They're evenly spaced along the contig by default; `--orf-selection first` takes the first N instead.
- **Tiara prefilter.**
  - Eukarya contigs are removed from the search by default.
  - Unknown and organelle contigs above 40% coding density are still searched (`--rescue-coding-density`).
  - `--rescue-eukarya` applies the same threshold to eukarya.
  - A 2% audit sample of removed contigs is searched anyway. If DIAMOND calls an audited contig Bacteria or Archaea, it is binned.
- **Length filter.** `--min-contig-len 1000` applies to preliminary and reassembled contigs.
- **Resume.** `--resume RUN_DIR` reuses QC, assembly, Prodigal, Tiara, DIAMOND, mapping and each finished Unicycler bin.

## Tools
fastp, MEGAHIT, Prodigal, DIAMOND (database with taxonomy), bwa, Tiara and phyloFlash.
- Reassembly: Unicycler and minimap2.
- Inspection: bowtie2 and samtools.
- `--checkm`: CheckM.
