# PathwayGenie

A standalone, Genie-style workflow based on the supplied GasGenie example:

**Nucleotide bins → Prodigal → HMMER → enzyme/activity evidence → pathway CSVs → flexible neighborhoods and reference-operon support.**

The original 114-profile HMM collection and its metadata are included. No further HMM download is needed. The supplied GasGenie.py is an example only; this implementation does not require MagicLamp environment variables, DIAMOND, R, pandas, or Biopython.

## Install and run

When this folder is nested under **MetaHopper**, the script, rules and complete
HMM panel stay together here. Use MetaHopper's existing environment; no separate
checkout, environment or HMM download is needed:

```bash
conda run --no-capture-output -n metahopper-lean \
  python ~/bin/MetaHopper/metahopper_pathways.py check

conda run --no-capture-output -n metahopper-lean \
  python ~/bin/MetaHopper/metahopper_pathways.py run \
  -i /path/to/MetaHopper_sample --rank genus --bin-set final -t 24 -j 3
```

The wrapper resolves the script, `Pathway_Gene_HMMs/` and configuration relative
to its own installation and uses the active Python. Its results are written to
`final/pathways/genus/` and discovered automatically by the report generator.
Add `--panel-only` to `check` to validate the resources without external tools.

Requires Python 3.9+, Prodigal 2.6.x, and HMMER 3.x. For example:

```bash
mamba install -n metahopper-lean --override-channels \
  -c conda-forge -c bioconda \
  --file ~/bin/MetaHopper/pathwaygenie/runtime-conda-spec.txt
```

Install these executables only if missing. This adds them to the existing
environment; no Python package is required by this module. The specification
records the Prodigal 2.6.3/HMMER 3.4 versions used by the supplied successful run.

To use the nested program directly on an arbitrary bin folder, run from this
PathwayGenie directory:

```bash
python PathwayGenie.py \
  -bin_dir /path/to/bins \
  -out /path/to/pathwaygenie_results \
  -t 24 -j 3
```

This processes three bins concurrently with eight HMMER worker CPUs per bin. Prodigal is single-threaded. Minor coordinator overhead is additional to the HMMER worker budget. `-t` is the total worker budget, not a per-bin budget.

FASTA suffixes are detected automatically: `.fa`, `.fna`, `.fasta`, `.fas`, `.fsa`, and compressed `.gz` versions. To select one suffix:

```bash
python PathwayGenie.py -bin_dir /path/to/bins -bin_ext fna \
  -out pathwaygenie_results -t 24 -j 3
```

Prodigal runs in **meta mode by default**, appropriate for fragmented bins. `--prodigal-mode single` selects its training-based genome mode; small or fragmented bins may not provide enough sequence to train it. Use `--recursive` for nested bin directories. The output directory is excluded from recursive input discovery. Input contig IDs and bin names may contain underscores; whitespace descriptions are preserved in the mapping. Duplicate contig IDs or duplicate bin IDs are rejected.

Executable overrides are available as `--prodigal /path/to/prodigal` and `--hmmsearch /path/to/hmmsearch`.

## Main outputs

Every matrix has **one row per bin**. Pathway/component names are columns. CSV quoting is handled correctly.

| File | Meaning |
|---|---|
| `pathway_heatmap.csv` | Percentage of **scored biochemical components detected**, from 0 to 100, for each pathway. Suitable for a heatmap. This is neither relative abundance nor percent of the genome. |
| `pathway_status.csv` | `complete_target_panel`, `partial_target_panel`, `not_detected_in_panel`, or `bin_failed`. |
| `route_status.csv` | Whether all components of the configured conventional route are supported, or whether there are missing/unmeasured steps. |
| `pathway_summary.long.csv` | Numerators, denominators, detected/missing components, unresolved route steps, operon support, and interpretation notes. |
| `pathways/PATHWAY.components.csv` | Granular component evidence: `1` detected, `0` no accepted evidence, `NA` not assayed or only incompletely assayed. |
| `pathways/PATHWAY.copy_counts.csv` | Unique predicted proteins supporting each component; alternative HMMs and multiple domains do not multiply the copy count. |
| `pathways/PATHWAY.component_status.csv` | Detailed labels: detected, candidate-only, not-detected, not-assayed, fusion-only-not-detected, or bin-failed. |
| `pathways/PATHWAY.gene_evidence.csv` | Evidence for the original gene labels, including shared-family and domain-only results. This is categorical because not every profile identifies an exact gene symbol. |
| `gene_target_evidence.csv` | All original target labels in one categorical matrix. |
| `operon_support_heatmap.csv` | Best reference-template support per bin/pathway, from 0 to 1; `NA` means no reference template is supplied or the bin failed. This is a heuristic support score, not a probability. |
| `neighborhoods.csv` | Same-contig neighborhoods, including coordinates, strand/order signatures, intervening-ORF counts, and components. |
| `operon_candidates.csv` | Neighborhoods matching reference-operon projections, with coverage/order/strand scores and source links. |
| `all_model_hits.csv` | HMMER hits above native model score cutoffs, including fragments, ambiguities, and context candidates. |
| `accepted_hits.csv` | Accepted component-level evidence, with original contig coordinates and source model metadata. |
| `review_hits.csv` | Fragments, conflicting hits, broad/supporting profiles, shared-family calls, and domain-only evidence. |
| `component_evidence.long.csv` | Direct trace from every component to its accepted and candidate protein IDs. |
| `component_definitions.csv` | Which components enter each heatmap denominator and each route definition. |
| `bin_qc.csv` | Input size, contig/ORF counts, success/failure, and cache reuse. |
| `config_used.json`, `run_manifest.json` | Exact configuration, model/input hashes, software versions, and run parameters. |

Per-bin `bins/` folders retain normalized nucleotide FASTA, contig-ID mappings, Prodigal GFF, coding sequences and proteins, ORF coordinate tables, raw HMMER output, logs, and verified cache records. ORF coordinates are **1-based, inclusive**; `orf_rank` is **0-based within each contig**.

### Reading the heatmap correctly

For example, 80 means four out of five configured, scored components have accepted evidence. It does **not** prove an 80%-functional pathway. A value of 100 means the **scored target panel** is complete; unmeasured reactions may still prevent a full route call. Read `route_status.csv` and `pathway_summary.long.csv` alongside the heatmap.

`supported_selected_route` means the configured conventional reaction set has evidence. It does not establish expression, metabolite production, flux, physiology, or coverage of every alternative route. `not_detected_in_panel` means no accepted marker evidence, not proven biological absence. Failed bins are `NA`, never silently converted to all-zero rows.

## Biological rules

- Prodigal predicts genes afresh; this workflow does not reuse the Bakta annotations. Results can differ because gene calling and annotation thresholds differ.
- NCBI/TIGRFAM models use native **sequence AND domain trusted cutoffs**. Pfam supporting models use native gathering cutoffs. A universal E-value is not substituted for curated scores.
- Default minimum union coverage of a model is **0.60**. Coverage is calculated from the union of above-threshold HMM-domain coordinates, not summed overlapping domains.
- Whole-fusion models used to establish multiple activities have a protected **0.85** coverage floor. A truncated fusion does not automatically establish all of its reactions. Lowering the global coverage flag cannot bypass that protection; edit a per-model floor only after reviewing the evidence.
- Multiple models supporting the same activity on the same protein count once. A protein can legitimately contribute different activities if a verified fusion or nonoverlapping domains support them.
- HisB phosphatase and dehydratase are scored separately. ThrA contributes aspartokinase and homoserine dehydrogenase activities. PheA/TyrA chorismate-mutase domains do not establish their dehydratase/dehydrogenase activities.
- General acetolactate synthase, DAHP synthase, shikimate kinase, and aspartokinase profiles support their activities without automatically identifying every isozyme name.
- Overlapping incompatible functions are compared by normalized margins above their model-specific sequence/domain cutoffs. A close conflict is marked ambiguous; a clearly weaker conflicting hit is not counted. The default tolerance is 0.10. This is a conservative assignment heuristic, not a statistically calibrated probability. Nonoverlapping fusion domains and compatible/shared activities are retained.
- `ilvY`, `tyrR`, and `metR` are reported as regulators but excluded from biochemical coverage denominators.
- AAT is still unresolved. Its candidate models are reported for review and never complete phenylalanine synthesis by default. OAT-like hits establish a homologous ornithine aminotransferase activity, not an insect origin for a bacterial-bin sequence.
- Missing standalone TrpD/TrpF, HisC/HisI, BioH, DapC, LysA, IlvA, or a specific aromatic aminotransferase cannot be called absent using this panel. Relevant unassayed components are explicitly `NA`; the selected TrpG/TrpD and TrpC/TrpF fusion models can supply positive evidence for the otherwise unassayed fusion-associated activities.
- Pathway labels cover conventional selected reactions, not all possible biosynthetic routes. The selected lysine route is the succinylated DAP route; MetH, ArgJ, type-I AroD, and other unrepresented alternatives need additional models/configuration. Central substrate supply is not assumed to have been independently demonstrated.

## Flexible neighborhoods and reference operons

Component evidence is aggregated across the entire bin. Neighborhoods are evaluated separately and **do not veto a pathway-component detection**. Dispersed biosynthetic genes can remain valid evidence.

Defaults allow:

- Up to **5 intervening predicted ORFs** between successive marked genes.
- A total neighborhood span of at most **20,000 bp**.
- Missing template slots: up to **2**, with at least **50%** of slots observed and at least two matched slots on **two distinct ORFs**. Some templates are therefore effectively stricter because their total size differs.
- Relative order and same-strand organization are **preferred**, rather than required.

Neighborhoods never bridge separate contigs. Distance counts all intervening predicted ORFs, not just HMM-positive genes. Reverse-strand/reverse-coordinate reference operons are recognized. A single bifunctional protein is not called an operon. Circular-origin wraparound, transcription, and read-supported contig linking are not evaluated.

Six reference templates are included:

| Reference pattern | Interpretation and source |
|---|---|
| `trpEDCBA_projected` | E. coli trpE–(TrpG/TrpD)–(TrpC/TrpF)–trpB–trpA projected onto assayed activities. [Yanofsky, 2007](https://rnajournal.cshlp.org/content/13/8/1141.full.pdf). |
| `hisGDCBHAFI_projected` | Histidine operon projected onto the assayed activities; not all original genes are measured. [Biosynthesis of Histidine](https://doi.org/10.1128/ecosalplus.3.6.1.9). |
| `ilvGMEDA_projected` | Large-subunit/E/D projection; this does not identify every gene or guarantee an IlvG-specific match. [Lawther et al., 1987](https://pubmed.ncbi.nlm.nih.gov/3550695/). |
| `argCBH` | The argCBH transcription unit, without requiring divergent argE to share its strand. [Piette et al., 1982](https://pubmed.ncbi.nlm.nih.gov/6761650/). |
| `carAB` | Shared carbamoyl-phosphate synthase subunit pair. [Piette et al., 1984](https://pubmed.ncbi.nlm.nih.gov/6330744/). |
| `thrABC` | ThrA represented by its homoserine-dehydrogenase activity, followed by ThrB/ThrC. [Comparative organization study](https://pmc.ncbi.nlm.nih.gov/articles/PMC1885841/). |

These are **reference patterns**, mostly from E. coli, not universal operon structures. Other pathways still receive generic neighborhood analysis even without a preset operon. Each preset retains its taxon, projection limits, and source URL in the JSON and output.

Support score = observed-slot fraction × order agreement × strand agreement, with ignored terms omitted. Order agreement is the fraction of observed slots that fit the reference relative order, allowing extra ORFs, missing slots, and supported fusions. Strand agreement is the fraction of matched ORFs oriented in the candidate transcriptional direction. The best forward/reverse score is retained. No score is added to biochemical completeness.

For stricter synteny:

```bash
python PathwayGenie.py -bin_dir /path/to/bins -out strict_results -t 24 -j 3 \
  -d 2 --max-neighborhood-bp 12000 \
  --operon-order strict --operon-strand strict --operon-max-missing 1
```

For flexible reduced/fragmented genomes:

```bash
python PathwayGenie.py -bin_dir /path/to/bins -out flexible_results -t 24 -j 3 \
  -d 8 --max-neighborhood-bp 30000 \
  --operon-order ignore --operon-strand prefer \
  --operon-max-missing 4 --operon-min-fraction 0.4
```

Here `strict` order means strict **relative order of observed slots**; the configured missing-slot allowance still applies. Neither setting joins different contigs.

## Editing the configuration

`pathwaygenie_config.json` contains:

1. `model_rules`: versioned HMM accession → original target label, supported activities, gene-call specificity, context-only status, and coverage floor.
2. `components`: activity IDs, assay scope (`targeted`, `fusion_only`, `not_assayed`), and notes.
3. `pathways`: scored components, full selected-route components, optional/context components, and original gene labels.
4. `operon_templates`: ordered slots, minimum coverage/slot count, missing-slot allowance, reference taxon, and source.

All scored components are required for a 100% target-panel score. Alternative gene families can support the **same activity ID** through different model rules. Full route definitions are explicit AND lists; multiple alternative routes should be represented as separate definitions rather than silently merged into one permissive route. A slot lists alternative component IDs, e.g. the HisB locus can be represented by either of its separately scored activities. It is positional evidence, not a claim that both activities are measured.

Use `--config my_configuration.json` to run an edited copy. Adding new HMMs also requires updating the HMM package's libraries, gene-to-model table, and download manifest; integrity checks reject inconsistent model bodies, lengths, or source checksums. Context-only source profiles cannot be silently promoted by changing just a configuration flag.

## Resume and failure handling

Rerun the same command to resume. ORF calls and HMM searches are reused only when their input hashes, parameters, tool versions, and output hashes match. Changes to DNA, Prodigal mode, HMM library, cutoff mode, or tools invalidate the appropriate cache. Changes to pathway rules, neighborhood tolerances, or post-search coverage filters regenerate summaries from retained raw searches.

Use `--force` to rerun both external stages. Logs and per-bin intermediates are retained. A failed bin gets `NA` matrix entries; successful bins are still summarized, and the process returns a nonzero exit code. Check `bin_qc.csv` and the relevant per-bin command log.

An output lock prevents concurrent runs writing the same directory. If a process is killed and leaves `.pathwaygenie.lock`, verify that the old process has stopped before deleting that one lock file and resuming.

## Validation and example data

Run the included unit checks with:

```bash
python -m unittest discover -s tests -v
```

`VALIDATION.md` records the checks and real-tool integration runs. `example_results/` contains example CSVs generated from the supplied Rickettsia assembly and the public E. coli K12 MG1655 assembly. These examples are pipeline demonstrations, not experimental validation of metabolic capacity. Inspect the detailed hit tables and route notes before drawing biological conclusions.

## Model provenance

The nested `Pathway_Gene_HMMs/` directory preserves the original collection, source URLs, per-gene metadata, license/attribution notices, and its separate downloader. Its original README describes the model-collection stage; this top-level README describes the new analysis workflow. TIGRFAM, NCBIfam/PRK, and Pfam source terms remain in force.

Primary software/method references:

- [Prodigal](https://github.com/hyattpd/Prodigal) and [its output documentation](https://github.com/hyattpd/Prodigal/wiki/Understanding-the-Prodigal-Output).
- [HMMER](http://hmmer.org/) and [model-specific cutoff options](https://i5k.nal.usda.gov/training/webapp/static/hmmer/pdf/Userguide.pdf).
- [NCBI protein annotation evidence](https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/).
