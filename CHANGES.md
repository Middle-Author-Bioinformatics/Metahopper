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
