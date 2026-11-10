# Pathway gene HMM panel

Prepared 6 November 2026. Source targets: the supplied pathway screenshot, not every gene in the attached genome annotations.

## Contents

- **70 unique target labels**: 69 from panel (a), plus host OAT from panel (b).
- **114 downloaded HMMER3 models**: 79 TIGRFAM, 33 NCBIfam/PRK, 2 Pfam.
- **106 function-directed models** and **8 context/alternative models**, held in separate libraries.
- Original individual HMMs, native thresholds, a gene-to-model mapping, source metadata, UniProt references, direct-download URLs, SHA-256 checksums, and a standard-library Python downloader.

Start with **index.html** for the searchable gene catalog, **gene_summary.tsv** for one row per target, and **gene_model_map.tsv** for one row per target/model relationship.

## How the models were selected

Curated UniProt entries were used to check enzyme functions and NCBIfam/Pfam cross-references. The majority of reference entries are reviewed E. coli K12 proteins; additional reviewed bacterial and insect entries resolve missing functions and nomenclature. The Drosophila AAT example Q8IPY3 is unreviewed and is labeled as such. Reference examples are not claims about the identity of the insect in the figure.

TIGRFAM profiles were preferred where they describe the intended function. NCBIfam/PRK profiles were used for gaps, variants, fusions, or more specific isozyme evidence. Pfam domain profiles are included only as supporting AAT/OAT candidates, not promoted to gene-specific calls. Source functions, taxonomic range, and family types were checked; genes were not mapped solely by their symbols. Not every selected profile is cross-referenced by the limited reference-protein sample: the mapping explicitly distinguishes direct UniProt cross-reference support from selection using NCBI functional metadata.

NCBI source: PGAP HMM release **20.0**, pinned download URLs. This is the maintained NCBI collection, rather than the archived TIGRFAM 15.0 collection. Pfam profiles are versioned in their ACC fields and checksummed; the API's current-version endpoint may change in the future.

## Quick start: HMMER3

Run from the extracted folder, with proteins in FASTA format:

```bash
hmmsearch --cpu 16 --cut_tc \
  --tblout pathway_hits.tbl \
  --domtblout pathway_domains.tbl \
  libraries/pathway_gene_models.hmm proteins.faa > pathway_search.txt
```

For supplementary NCBI context and alternative-symbol models:

```bash
hmmsearch --cpu 16 --cut_tc \
  --domtblout context_ncbi.domtbl \
  libraries/context_ncbi.hmm proteins.faa > context_ncbi.txt
```

For broad Pfam candidate domains, use their native gathering thresholds:

```bash
hmmsearch --cpu 16 --cut_ga \
  --domtblout context_pfam.domtbl \
  libraries/context_pfam.hmm proteins.faa > context_pfam.txt
```

`hmmsearch` accepts multiple HMMs without hmmpress. To use `hmmscan`, first run `hmmpress` on the desired library. Bundled NAME fields are set to the versioned model accession to prevent duplicate-name problems. All emission/transition parameters and thresholds are preserved. Original individual files remain unmodified.

NCBI trusted cutoffs include sequence and domain scores. Do not replace them with one universal E-value. Domain boundaries, HMM coverage, sequence length, taxonomic applicability, and competing family hits still need review. No universal coverage cutoff was invented for this panel.

## Interpretation limits that matter

1. **AAT is unresolved.** The figure supplies only an acronym. NF006719 (aromatic transaminase), NF006228 (aspartate aminotransferase family), and PF00155 (broad aminotransferase domain) are provisional context-only candidates. NF006719 is cross-referenced by human and Drosophila aspartate aminotransferases, but its curated NCBI range is Bacteria. None proves which host AAT or substrate specificity the figure means. Obtain the host protein accession or intended reaction before making a definitive call.
2. **aroH and aroQ have multiple meanings.** Main-panel aroH is DAHP synthase; main-panel aroQ is type II dehydroquinate dehydratase, based on their placement in the chorismate group. Chorismate-mutase meanings are supplied only in the context library. Confirm against the source paper before treating these interpretations as definitive.
3. **ribC is interpreted by reaction, not nomenclature alone.** Bacterial riboflavin synthase (often ribE) and archaeal-type riboflavin synthase are different families; lumazine synthase (ribH, sometimes ribE) is a different reaction again.
4. **ilvY, tyrR, and metR are regulators**, not catalytic pathway requirements. Shared acetolactate synthase models may not distinguish IlvB/IlvG/IlvI or all small-subunit isozymes. General aspartokinase and shikimate-kinase models do not uniquely establish a specific regulatory isozyme name.
5. **Bifunctional proteins require activity-level interpretation.** HisB, ThrA, TrpC/TrpF, PheA, TyrA, and TrpG/TrpD may be split or fused. A single domain hit does not establish all reactions performed by a bifunctional protein. The table lists source model types and fusion descriptions.
6. **The screenshot is not a complete universal pathway definition.** Its repeated trpC, cysK, and hisA rows were preserved in screenshot_gene_rows.tsv and collapsed in the unique target list. They were not silently changed to other gene names. Missing labels such as trpD/trpF, hisC/hisI, lysA, metH, or alternative pathway enzymes were not added as independent targets. Some listed fusion models incidentally cover extra activities. The figure's annotations and low resolution leave unresolved naming questions.
7. **A negative hit is not proof of absence.** Several profiles are taxon-specific; this panel does not exhaust every clade, nonorthologous replacement, or alternative route. Review assembly completeness, pseudogenes, contig ends, and uncalled genes. Protein searches alone cannot recover an uncalled ORF.

No presence/absence matrix or biological performance benchmark is claimed here. This deliverable collects and checks the models; it does not claim to have validated pathway function in the uploaded Rickettsia bin.

## Download again, or download selected genes

The HMMs are already included; no download is needed to search them. To verify/recreate the local panel:

```bash
python3 download_hmms.py
```

To download only selected targets into a separate directory:

```bash
python3 download_hmms.py --genes trpA trpB carA carB OAT --output selected_models
```

The script needs only Python 3.9+ and internet access. It uses official source URLs, decompresses Pfam responses, checks SHA-256, and rebuilds the three libraries. A source version change triggers an error rather than a silent substitution. NCBI may eventually retire individual-file directories; the release tarball is the archival fallback at https://ftp.ncbi.nlm.nih.gov/hmm/20.0/hmm_PGAP.HMM.tgz .

## Sources and attribution

- UniProtKB: https://www.uniprot.org/ ; individual reference links in index.html and accession lists in provenance/uniprot_references.tsv. UniProt data: CC BY 4.0 (https://www.uniprot.org/help/license).
- NCBI PGAP HMM collection: https://ftp.ncbi.nlm.nih.gov/hmm/20.0/ ; source metadata subset and release notes included.
- NCBI protein annotation evidence documentation: https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/
- TIGRFAMs at NCBI: https://www.ncbi.nlm.nih.gov/refseq/annotation_prok/tigrfams/ . TIGRFAM data: CC BY-SA 4.0. Original attribution/license notice included in provenance/TIGRFAM_license.txt. Individual source profiles are unmodified; bundled NAME identifiers alone were normalized to accessions. This compiled collection and its mapping are supplied under CC BY-SA 4.0; third-party source licenses remain in force.
- Li et al. (2021), RefSeq: expanding the Prokaryotic Genome Annotation Pipeline reach with protein family model curation. Nucleic Acids Research 49:D1020-D1028. https://doi.org/10.1093/nar/gkaa1105
- Pfam via InterPro: https://www.ebi.ac.uk/interpro/ . Pfam downloadable data: CC0 1.0, https://interpro-documentation.readthedocs.io/en/latest/license.html .

See VALIDATION.txt for file and HMM parser checks. Integrity checks do not establish biological sensitivity or substrate specificity.
