#!/usr/bin/env python3
"""
metatax_binner.py

End-to-end pipeline for assembly, microbial-contig retention, taxonomic genome binning,
seed-and-extension, assembly refinement, and quality assessment. Inputs may be reads alone,
contigs alone, or contigs plus their paired reads.

    R1.fastq, R2.fastq  (reads-only mode)
        -> [optional] fastp poly-G trimming (--trim-polyg -- for NextSeq/NovaSeq
           two-channel-chemistry poly-G tail artifacts)
        -> QC: Trimmomatic (adapter clip + quality trim) + FLASH (merge overlapping pairs)
           [skippable with --skip-qc]
        -> MEGAHIT assembly
    initial contigs.fasta  (supplied directly, or produced above)
        -> preliminary Prodigal (ORF/gene prediction, meta mode)
        -> gene-density triage: quarantine long, low-coding-density eukaryotic-like
           contigs before DIAMOND (default; disable with --skip-contig-triage)
        -> preliminary DIAMOND blastp vs. taxonomy-enabled nr-tax.dmnd; rank names are
           emitted directly from taxonomy embedded in the database
        -> preliminary per-ORF hits -> per-contig taxonomic classification
        -> drop host-animal/plant contamination (Metazoa/Viridiplantae kingdoms, by
           default -- see --exclude-kingdoms); Bacteria, Archaea, Fungi, and protists
           are kept
        -> [when reads exist] map all reads to contigs; calculate GC and depth; prevent
           strong joint GC+coverage outliers from acting as extension seeds
        -> preliminary seed-bin FASTAs at the requested ranks
        -> [default whenever reads exist] one competitive Bowtie2 seed mapping against all
           bins at a source rank; uniquely winning templates enter that bin and ties are
           excluded, followed by guarded BBDuk frontier extension
        -> focused SPAdes/metaSPAdes, Unicycler circularization attempt, and Pilon polishing
        -> consolidate successful reassemblies and original fallbacks from one source rank
        -> final Prodigal + DIAMOND classification of the consolidated contigs
        -> remap all reads and repeat GC+coverage coherence refinement; joint outliers
           are demoted to Unclassified at the selected rank and more-specific ranks
        -> reapply animal/plant filtering and rebuild final multirank FASTA bins
        -> final QUAST contiguity + CheckM completeness/contamination summaries

The default seed-and-extension stage is aimed at
           fragmented, low-abundance, or fast-diverging genomes (e.g. host symbionts) that
the whole-metagenome MEGAHIT co-assembly split into many pieces. Pass
--skip-reassembly to stop after the preliminary bins and assess those directly.

CLASSIFICATION METHOD (what this script does and why)
-------------------------------------------------------
Classifying a *contig* from many *ORF* hits is the same problem CAT/BAT, MEGAN's LCA,
and Kraken-style consensus callers solve, and the approach used here follows that lineage:

  1. For each ORF, keep the diamond hits within `--bitscore-range` (default 90%) of that
     ORF's best bitscore (i.e. a "bit-score competitive set", not just the single top hit).
     This avoids over-trusting one alignment when several equally-good references exist.
  2. Each kept hit's organism name is parsed from its stitle (the "[Genus species]" at the
     hit's rank names are read directly from nr-tax.dmnd and its bitscore contributes a
     vote at each available rank.
  3. Votes are tallied *independently at each rank* (genus, species by default -- ask for
     domain/phylum/family too via --ranks) across every ORF on the contig. At each rank,
     the taxon with the largest share of total bitscore wins IF its share clears
     `--min-support` (default 0.5, i.e. a majority of the weighted evidence). Set
     --min-support 0 for a pure plurality ("most votes wins") call instead. Otherwise the
     contig is "Unclassified" at that rank.
  4. This gives a bitscore-weighted majority (or plurality) vote per rank -- simpler to
     implement/debug than a full weighted-LCA tree walk, but captures the same core idea.

This is a reasonable default for whole-contig binning. If you later want the more rigorous
tree-based LCA that CAT/BAT itself implements, swap `classify_contig()` for a call to CAT,
and feed its per-contig output into `bin_contigs()` / `summarize_bins()` below unchanged.

EXTERNAL TOOLS REQUIRED ON $PATH
-------------------------------------------------------
  Only needed if starting from raw reads (--r1/--r2):
    trimmomatic         http://www.usadellab.org/cms/?page=trimmomatic  (skip with --skip-qc)
    flash               https://ccb.jhu.edu/software/FLASH/             (skip with --skip-qc)
    pigz                https://zlib.net/pigz/                          (falls back to gzip)
    fastp               https://github.com/OpenGene/fastp               (only with --trim-polyg)
    megahit             https://github.com/voutcn/megahit

  Always needed:
    prodigal            https://github.com/hyattpd/Prodigal
    diamond >=2.1.17    https://github.com/bbuchfink/diamond
                         (requires nr-tax.dmnd built with --taxonmap, --taxonnodes,
                         and --taxonnames)
    quast.py            https://github.com/ablab/quast          (optional; falls back to a
                         built-in N50/L50/GC calculator if missing or --skip-quast is given)
    checkm (lineage_wf) https://github.com/Ecogenomics/CheckM    (optional; skipped with a
                         warning if missing or --skip-checkm is given)

  Needed for read-derived coverage refinement (disable with --skip-bin-refinement):
    bowtie2, bowtie2-build   https://github.com/BenLangmead/bowtie2
    samtools                 https://github.com/samtools/samtools

  Additionally needed for seed-and-extension (disable with --skip-reassembly):
    samtools                 (needs `samtools view -N`)
    bbduk.sh (BBMap)         https://sourceforge.net/projects/bbmap/
    spades.py                https://github.com/ablab/spades
    unicycler                https://github.com/rrwick/Unicycler
    pilon                    https://github.com/broadinstitute/pilon
    trimmomatic, flash       (same as above -- required even if --skip-qc was used for the
                              main assembly, since recruited reads always get QC'd)

EXAMPLES
-------------------------------------------------------
  # From an existing assembly:
  python metatax_binner.py \\
      -i contigs.fasta \\
      -d /dbs/nr-tax.dmnd \\
      -o metatax_out \\
      -t 16

  # From raw paired-end reads (full default workflow, including seed-and-extension):
  python metatax_binner.py \\
      --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz \\
      --trimmomatic-folder /path/to/Trimmomatic-0.39 \\
      -d /dbs/nr-tax.dmnd \\
      -o metatax_out \\
      -t 16

OUTPUT LAYOUT
-------------------------------------------------------
  <outdir>/
    qc/proteins... (Trimmomatic + FLASH intermediates; only with --r1/--r2 and no --skip-qc)
    megahit/final.contigs.fa, megahit.log            (only with --r1/--r2)
    prodigal/proteins.faa, genes.gff
    diamond/hits.tsv
    classification/contig_classification.tsv
    classification/bin_refinement.tsv, coverage/reads_to_contigs.sorted.bam  (with reads)
    classification/excluded_eukaryotic_like_gene_density.fasta  (when triage quarantines any)
    classification/seed_refinement_outliers.fasta               (when refinement finds any)
    classification/excluded_animal_plant_contamination.fasta  (only if --exclude-kingdoms
                                                                 dropped anything)
    bins/genus/<taxon>.fasta ... bin_membership.tsv        (preliminary seed bins)
    reassembly/<rank>/<bin>/reassembled.fasta, recruitment.tsv
    final/assembly/consolidated_contigs.fasta, contig_provenance.tsv
    final/prodigal/proteins.faa, genes.gff
    final/diamond/hits.tsv
    final/classification/contig_classification.tsv
    final/classification/bin_refinement.tsv, refinement_outliers.fasta
    final/bins/<rank>/<taxon>.fasta, bin_membership.tsv, summary.tsv
"""

import argparse
import csv
import gzip
import logging
import math
import re
import shlex
import shutil
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

# --------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------

# Both "domain" and "superkingdom" are requested because taxonomy-enabled DIAMOND
# databases can expose either top-rank label. "kingdom" lets the default retention filter
# distinguish Fungi from Metazoa and Viridiplantae within Eukaryota.
WANTED_RANKS = ["domain", "superkingdom", "kingdom", "phylum", "family", "genus", "species"]
BIN_RANKS_DEFAULT = ["genus", "species"]

# Eukaryotic kingdoms dropped from binning by default (--exclude-kingdoms). Metazoa and
# Viridiplantae are the only two of NCBI's formal "kingdom"-rank taxa besides Fungi --
# everything else under Eukaryota (the various protist lineages: SAR, Excavata,
# Amoebozoa, etc.) has no kingdom-rank ancestor at all in NCBI's taxonomy, so it's never
# matched by this filter and passes through untouched, same as Fungi.
DEFAULT_EXCLUDED_KINGDOMS = ["Metazoa", "Viridiplantae"]

# Deliberately conservative heuristic/refinement defaults. Gene-density triage only
# quarantines long contigs with very little Prodigal-predicted coding sequence. GC and
# coverage refinement requires a joint outlier: either signal by itself is reported but
# never sufficient to reject a seed or demote a final taxonomic assignment.
TRIAGE_MIN_LENGTH_BP = 3000
TRIAGE_MAX_EUKARYOTIC_CODING_DENSITY = 0.35
REFINEMENT_MIN_CONTIG_LENGTH_BP = 2000
REFINEMENT_MIN_REFERENCE_CONTIGS = 4
REFINEMENT_MIN_TAXON_SUPPORT = 0.70
REFINEMENT_GC_ABSOLUTE_FLOOR = 0.08
REFINEMENT_COVERAGE_LOG2_FLOOR = 2.0
REFINEMENT_ROBUST_Z = 4.0
REFINEMENT_RANK_ORDER = [
    "domain", "superkingdom", "kingdom", "phylum", "family", "genus", "species",
]

DIAMOND_RANK_FIELDS = {
    "domain": "sdomain",
    "superkingdom": "ssuperkingdom",
    "kingdom": "skingdom",
    "phylum": "sphylum",
    "family": "sfamily",
    "genus": "sgenus",
    "species": "sspecies",
}

# DIAMOND >=2.1.17 exposes generic sRANK output fields. The accession-to-taxid map,
# NCBI tree, and scientific names used to resolve these fields are baked into nr-tax.dmnd
# at `diamond makedb` time, so no external taxdump is needed at runtime.
DIAMOND_FIELDS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "staxids", *DIAMOND_RANK_FIELDS.values(), "stitle",
]


log = logging.getLogger("metatax_binner")


# --------------------------------------------------------------------------------------
# Small utilities
# --------------------------------------------------------------------------------------

def setup_logging(verbose: bool = False) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )


class StepCounter:
    """Tiny helper so log lines read 'Step i/N: ...' even though N varies depending on
    whether QC/MEGAHIT ran (raw-reads mode) or we started straight from a contigs FASTA."""

    def __init__(self, total: int):
        self.total = total
        self.i = 0

    def next(self, label: str) -> None:
        self.i += 1
        log.info("Step %d/%d: %s", self.i, self.total, label)


def which_or_die(tool: str) -> str:
    path = shutil.which(tool)
    if path is None:
        log.error("Required tool '%s' not found on PATH.", tool)
        sys.exit(1)
    return path


def require_diamond_taxonomy_fields() -> None:
    """Require the DIAMOND release that introduced generic ``sRANK`` fields."""
    proc = subprocess.run(["diamond", "version"], capture_output=True, text=True)
    version_text = f"{proc.stdout} {proc.stderr}".strip()
    match = re.search(r"(\d+)\.(\d+)\.(\d+)", version_text)
    if proc.returncode != 0 or not match:
        log.error("Could not determine DIAMOND version: %s", version_text or "no output")
        sys.exit(1)
    version = tuple(int(part) for part in match.groups())
    if version < (2, 1, 17):
        log.error(
            "DIAMOND %s is too old. MetaHopper requires >=2.1.17 for taxonomy fields "
            "such as sphylum, sfamily, and sgenus.",
            ".".join(str(part) for part in version),
        )
        sys.exit(1)


def run_cmd(cmd, log_file: Path = None, cwd: Path = None) -> None:
    log.info("Running: %s", " ".join(str(c) for c in cmd))
    with open(log_file, "a") if log_file else open("/dev/null", "w") as lf:
        proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, cwd=cwd, text=True)
    if proc.returncode != 0:
        tail = ""
        if log_file and Path(log_file).exists():
            tail = "\n".join(Path(log_file).read_text().splitlines()[-30:])
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(str(c) for c in cmd)}\n{tail}")


def run_pipeline(cmds, stdout_path: Path = None, log_file: Path = None) -> None:
    """Runs cmds[0] | cmds[1] | ... | cmds[-1], like a bash pipe. The last stage's stdout is
    written to stdout_path if given, otherwise discarded; every stage's stderr is appended
    to log_file. Used for the bowtie2 | samtools and samtools collate | samtools fastq
    pipelines in the targeted bin-reassembly module (Step 7)."""
    opened = []
    lf = open(log_file, "a") if log_file else subprocess.DEVNULL
    if lf is not subprocess.DEVNULL:
        opened.append(lf)
    procs = []
    prev_stdout = None
    try:
        for i, cmd in enumerate(cmds):
            is_last = i == len(cmds) - 1
            if is_last:
                out = open(stdout_path, "wb") if stdout_path else subprocess.DEVNULL
                if stdout_path:
                    opened.append(out)
            else:
                out = subprocess.PIPE
            proc = subprocess.Popen(cmd, stdin=prev_stdout, stdout=out, stderr=lf)
            if prev_stdout is not None:
                prev_stdout.close()
            prev_stdout = proc.stdout
            procs.append(proc)
        for proc in procs:
            proc.wait()
    finally:
        for fh in opened:
            fh.close()
    for cmd, proc in zip(cmds, procs):
        if proc.returncode != 0:
            raise RuntimeError(f"Pipeline stage failed ({proc.returncode}): {' '.join(str(c) for c in cmd)}")


def count_fastq_reads(path: Path) -> int:
    """Fast read count via `zcat -f | wc -l` (handles both gzipped and plain FASTQ)."""
    if not path or not Path(path).exists() or Path(path).stat().st_size == 0:
        return 0
    proc = subprocess.run(f"zcat -f -- {shlex.quote(str(path))} | wc -l",
                           shell=True, capture_output=True, text=True, check=True)
    lines = int((proc.stdout or "0").strip() or 0)
    return lines // 4


def sanitize(name: str) -> str:
    if not name:
        name = "Unclassified"
    keep = []
    for ch in name:
        keep.append(ch if (ch.isalnum() or ch in "-._") else "_")
    out = "".join(keep).strip("_")
    return out or "Unclassified"


def read_fasta(path: Path) -> dict:
    """Minimal FASTA reader: {seq_id: sequence}. seq_id = header up to first whitespace."""
    seqs = {}
    header = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            seqs[header] = "".join(chunks)
    return seqs


def write_fasta(path: Path, records: dict, wrap: int = 70) -> None:
    with open(path, "w") as fh:
        for seq_id, seq in records.items():
            fh.write(f">{seq_id}\n")
            for i in range(0, len(seq), wrap):
                fh.write(seq[i:i + wrap] + "\n")


# --------------------------------------------------------------------------------------
# Step 0a: poly-G trimming (fastp) -- optional, only for raw-reads (--r1/--r2) input
# --------------------------------------------------------------------------------------

def trim_poly_g(r1: Path, r2: Path, outdir: Path, threads: int,
                 fastp_cmd: str = "fastp", poly_g_min_len: int = 10) -> tuple:
    """Trims poly-G tails with fastp's dedicated detector, and nothing else.

    Poly-G runs are a well-known artifact of two-channel Illumina chemistry (NextSeq,
    NovaSeq): those instruments call a dark/no-signal cycle as 'G', so reads that run
    past the end of a short insert (or into a quality dropout) pick up a run of spurious
    G's that isn't part of the biological sequence. Trimmomatic's adapter/quality-trim
    steps don't reliably catch this -- poly-G runs are often called at deceptively high
    quality and don't match the ILLUMINACLIP adapter sequences -- so this runs before
    everything else in the QC chain. Every other fastp filter is disabled here (adapter
    trimming, quality filtering, length filtering) so it doesn't duplicate/interfere with
    what Trimmomatic + FLASH already do downstream.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    out1 = outdir / "polyg_trimmed_1.fastq"
    out2 = outdir / "polyg_trimmed_2.fastq"
    log_file = outdir / "fastp_polyg.log"
    cmd = [
        fastp_cmd, "-i", str(r1), "-I", str(r2), "-o", str(out1), "-O", str(out2),
        "--trim_poly_g", "--poly_g_min_len", str(poly_g_min_len),
        "--disable_adapter_trimming", "--disable_quality_filtering", "--disable_length_filtering",
        "--thread", str(threads),
        "--json", str(outdir / "fastp_polyg.json"), "--html", str(outdir / "fastp_polyg.html"),
    ]
    run_cmd(cmd, log_file=log_file)
    return out1, out2


# --------------------------------------------------------------------------------------
# Step 0b: QC (Trimmomatic + FLASH) -- optional, only for raw-reads (--r1/--r2) input
# --------------------------------------------------------------------------------------

def run_qc(r1: Path, r2: Path, outdir: Path, threads: int,
           trimmomatic_cmd: str, trimmomatic_folder: Path,
           flash_cmd: str, flash_max_overlap: int, pigz_cmd: str,
           qc_quality: int, qc_minlen: int, keep_tmp: bool = False) -> tuple:
    """Adapter-trim (Trimmomatic ILLUMINACLIP), quality-trim (SLIDINGWINDOW), and merge
    overlapping pairs (FLASH) -- the same recipe as QC_reads.sh (Trimmomatic PE adapter
    clip -> quality-trim the orphaned singles -> FLASH-merge the still-paired reads ->
    quality-trim the merged fragments and the unmerged pairs -> pool every singleton/
    merged read into one file -> gzip). Returns (paired_r1.gz, paired_r2.gz, unpaired.gz).
    """
    outdir.mkdir(parents=True, exist_ok=True)
    log_file = outdir / "qc.log"
    base = outdir / "reads"

    adapters = Path(trimmomatic_folder) / "adapters" / "TruSeq3-PE-2.fa"
    if not adapters.exists():
        log.error("Adapter file not found: %s (check --trimmomatic-folder)", adapters)
        sys.exit(1)

    # 1. Adapter clip (paired)
    run_cmd([
        trimmomatic_cmd, "PE", "-threads", str(threads), "-baseout", str(base),
        str(r1), str(r2), f"ILLUMINACLIP:{adapters}:2:30:10",
    ], log_file=log_file)
    p1, u1 = Path(f"{base}_1P"), Path(f"{base}_1U")
    p2, u2 = Path(f"{base}_2P"), Path(f"{base}_2U")

    # 2. Quality-trim the reads that lost their mate during adapter clipping
    u1_qt, u2_qt = Path(f"{u1}.qual_trimmed"), Path(f"{u2}.qual_trimmed")
    for u_in, u_out in ((u1, u1_qt), (u2, u2_qt)):
        run_cmd([trimmomatic_cmd, "SE", "-threads", str(threads), str(u_in), str(u_out),
                 f"SLIDINGWINDOW:4:{qc_quality}", f"MINLEN:{qc_minlen}"], log_file=log_file)

    # 3. Merge overlapping pairs
    run_cmd([
        flash_cmd, "--threads", str(threads), "--output-prefix", "flash",
        "--max-overlap", str(flash_max_overlap), str(p1), str(p2),
        "--output-directory", str(outdir),
    ], log_file=log_file)
    merged = outdir / "flash.extendedFrags.fastq"
    nc1 = outdir / "flash.notCombined_1.fastq"
    nc2 = outdir / "flash.notCombined_2.fastq"

    # 4. Quality-trim the merged fragments
    merged_final = outdir / "merged.final.fastq"
    run_cmd([trimmomatic_cmd, "SE", "-threads", str(threads), str(merged), str(merged_final),
             f"SLIDINGWINDOW:4:{qc_quality}", f"MINLEN:{qc_minlen}"], log_file=log_file)

    # 5. Quality-trim the pairs FLASH couldn't merge
    nc_base = outdir / "notcombined.final"
    run_cmd([trimmomatic_cmd, "PE", "-threads", str(threads), str(nc1), str(nc2),
             "-baseout", str(nc_base), f"SLIDINGWINDOW:4:{qc_quality}", f"MINLEN:{qc_minlen}"],
            log_file=log_file)
    nc_p1, nc_u1 = Path(f"{nc_base}_1P"), Path(f"{nc_base}_1U")
    nc_p2, nc_u2 = Path(f"{nc_base}_2P"), Path(f"{nc_base}_2U")

    # 6. Pool every orphan/singleton/merged read into one unpaired file
    unpaired = outdir / "unpaired.fq"
    with open(unpaired, "w") as out_fh:
        for part in (merged_final, nc_u1, nc_u2, u1_qt, u2_qt):
            if part.exists() and part.stat().st_size > 0:
                with open(part) as in_fh:
                    shutil.copyfileobj(in_fh, out_fh)

    # 7. Compress the final paired + unpaired reads (pigz if available, else gzip)
    gzip_cmd = pigz_cmd if shutil.which(pigz_cmd) else "gzip"
    gzip_opts = ["--best", "--processes", str(threads)] if gzip_cmd == pigz_cmd else ["-9"]
    run_cmd([gzip_cmd, "-f", *gzip_opts, str(nc_p1), str(nc_p2), str(unpaired)], log_file=log_file)
    final_r1, final_r2, final_u = Path(f"{nc_p1}.gz"), Path(f"{nc_p2}.gz"), Path(f"{unpaired}.gz")

    if not keep_tmp:
        for f in (p1, u1, p2, u2, u1_qt, u2_qt, merged, nc1, nc2, merged_final, nc_u1, nc_u2):
            Path(f).unlink(missing_ok=True)

    log.info("QC done: %s, %s, %s", final_r1, final_r2, final_u)
    return final_r1, final_r2, final_u


# --------------------------------------------------------------------------------------
# Step 0c: MEGAHIT assembly -- optional, only for raw-reads (--r1/--r2) input
# --------------------------------------------------------------------------------------

def run_megahit(r1: Path, r2: Path, unpaired: Path, outdir: Path, threads: int,
                 megahit_cmd: str = "megahit", min_contig_len: int = None,
                 extra_args: str = None) -> Path:
    """Runs MEGAHIT on paired (+ optional unpaired) reads. MEGAHIT refuses to write into
    an existing output directory, so it's removed first if present (e.g. on a rerun)."""
    if outdir.exists():
        log.warning("Removing existing MEGAHIT output dir: %s", outdir)
        shutil.rmtree(outdir)
    outdir.parent.mkdir(parents=True, exist_ok=True)

    cmd = [megahit_cmd, "-1", str(r1), "-2", str(r2), "-t", str(threads), "-o", str(outdir)]
    if unpaired:
        cmd += ["-r", str(unpaired)]
    if min_contig_len:
        cmd += ["--min-contig-len", str(min_contig_len)]
    if extra_args:
        cmd += shlex.split(extra_args)

    log_file = outdir.parent / "megahit.log"
    run_cmd(cmd, log_file=log_file)

    contigs = outdir / "final.contigs.fa"
    if not contigs.exists():
        raise RuntimeError(f"MEGAHIT did not produce {contigs}")
    return contigs


# --------------------------------------------------------------------------------------
# Step 1: Prodigal
# --------------------------------------------------------------------------------------

def run_prodigal(fasta: Path, outdir: Path, mode: str = "meta") -> tuple:
    outdir.mkdir(parents=True, exist_ok=True)
    faa = outdir / "proteins.faa"
    gff = outdir / "genes.gff"
    log_file = outdir / "prodigal.log"
    cmd = [
        "prodigal", "-i", str(fasta), "-a", str(faa), "-f", "gff", "-o", str(gff),
        "-p", mode, "-q",
    ]
    run_cmd(cmd, log_file=log_file)
    return faa, gff


def parse_orf_to_contig(faa: Path) -> dict:
    """Prodigal names ORFs '<contig_id>_<gene_number>'. Map orf_id -> contig_id."""
    mapping = {}
    with open(faa) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            orf_id = line[1:].split()[0]
            idx = orf_id.rfind("_")
            if idx == -1 or not orf_id[idx + 1:].isdigit():
                log.warning("Could not parse contig id from ORF header: %s", orf_id)
                continue
            contig_id = orf_id[:idx]
            mapping[orf_id] = contig_id
    return mapping


def compute_contig_metrics(contig_seqs: dict, gff: Path) -> dict:
    """Calculate sequence and Prodigal gene-density metrics for every contig.

    Coding density is the fraction of contig bases covered by the union of predicted CDS
    intervals. Using the interval union prevents overlapping calls from producing values
    above one. Coverage fields are populated later when paired reads are available.
    """
    cds_intervals = defaultdict(list)
    if gff and Path(gff).exists():
        with open(gff) as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5 or (len(fields) >= 3 and fields[2] != "CDS"):
                    continue
                try:
                    start = max(0, int(fields[3]) - 1)
                    end = int(fields[4])
                except ValueError:
                    continue
                if end > start:
                    cds_intervals[fields[0]].append((start, end))

    metrics = {}
    for contig_id, sequence in contig_seqs.items():
        seq = sequence.upper()
        length = len(seq)
        acgt = sum(seq.count(base) for base in "ACGT")
        gc_fraction = ((seq.count("G") + seq.count("C")) / acgt) if acgt else 0.0
        intervals = sorted(cds_intervals.get(contig_id, []))
        merged = []
        for start, end in intervals:
            start, end = min(start, length), min(end, length)
            if end <= start:
                continue
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        coding_bp = sum(end - start for start, end in merged)
        cds_lengths = [min(end, length) - min(start, length) for start, end in intervals]
        cds_lengths = [value for value in cds_lengths if value > 0]
        metrics[contig_id] = {
            "length_bp": length,
            "gc_fraction": gc_fraction,
            "coding_bp": coding_bp,
            "coding_density": (coding_bp / length) if length else 0.0,
            "n_cds": len(intervals),
            "mean_cds_length_bp": (
                sum(cds_lengths) / len(cds_lengths) if cds_lengths else 0.0
            ),
            "mean_depth": None,
            "covered_fraction": None,
        }
    return metrics


def triage_contigs(contig_metrics: dict, disabled: bool = False) -> tuple:
    """Return ``(call_by_contig, quarantined_ids)`` from coding density.

    Short contigs are retained as uncertain because coding density is unstable at short
    lengths. The threshold is intentionally stringent: this is a preliminary screen, not
    a substitute for the later DIAMOND taxonomy classification.
    """
    calls = {}
    quarantined = set()
    for contig_id, metric in contig_metrics.items():
        if disabled:
            call = "disabled"
        elif metric["length_bp"] < TRIAGE_MIN_LENGTH_BP:
            call = "uncertain_short_retained"
        elif metric["coding_density"] < TRIAGE_MAX_EUKARYOTIC_CODING_DENSITY:
            call = "eukaryotic_like_quarantined"
            quarantined.add(contig_id)
        else:
            call = "prokaryotic_like_retained"
        calls[contig_id] = call
    return calls, quarantined


def write_candidate_proteins(faa: Path, orf_to_contig: dict, excluded_contigs: set,
                             out_path: Path) -> Path:
    """Write the Prodigal proteins whose parent contigs passed preliminary triage."""
    if not excluded_contigs:
        out_path.unlink(missing_ok=True)
        return faa
    proteins = read_fasta(faa)
    retained = {
        orf_id: sequence for orf_id, sequence in proteins.items()
        if orf_to_contig.get(orf_id) not in excluded_contigs
    }
    write_fasta(out_path, retained)
    return out_path


# --------------------------------------------------------------------------------------
# Step 2: DIAMOND
# --------------------------------------------------------------------------------------

def run_diamond(query_faa: Path, db: Path, outdir: Path, threads: int, evalue: float,
                 max_target_seqs: int) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    hits_tsv = outdir / "hits.tsv"
    log_file = outdir / "diamond.log"
    if not query_faa.exists() or query_faa.stat().st_size == 0:
        log.warning("No proteins passed contig triage; writing an empty DIAMOND result.")
        hits_tsv.write_text("")
        return hits_tsv
    cmd = [
        "diamond", "blastp",
        "-q", str(query_faa),
        "-d", str(db),
        "-o", str(hits_tsv),
        "-e", str(evalue),
        "-k", str(max_target_seqs),
        "--threads", str(threads),
        "--outfmt", "6", *DIAMOND_FIELDS,
    ]
    run_cmd(cmd, log_file=log_file)
    return hits_tsv


class Hit:
    __slots__ = ("sseqid", "pident", "length", "evalue", "bitscore", "lineage")

    def __init__(self, sseqid, pident, length, evalue, bitscore, lineage):
        self.sseqid = sseqid
        self.pident = pident
        self.length = length
        self.evalue = evalue
        self.bitscore = bitscore
        self.lineage = lineage


def parse_diamond_taxa(value: str) -> tuple:
    """Return the unique taxa represented by one DIAMOND taxonomy field.

    A protein accession can be associated with more than one taxid. DIAMOND may render
    the resulting names with semicolon or ``<>`` separators. Keeping every unique value
    and splitting the hit's vote avoids making an arbitrary first-name assignment.
    """
    if not value or value in {"N/A", "*", "0"}:
        return ()
    return tuple(dict.fromkeys(
        item.strip() for item in re.split(r"\s*(?:<>|;)\s*", value) if item.strip()
    ))


def parse_diamond_hits(hits_tsv: Path) -> dict:
    """Return ``{orf_id: [Hit, ...]}`` with hierarchy read from nr-tax.dmnd."""
    hits_by_orf = defaultdict(list)
    idx = {f: i for i, f in enumerate(DIAMOND_FIELDS)}
    any_taxonomy = False
    with open(hits_tsv) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(DIAMOND_FIELDS):
                continue
            lineage = {
                rank: parse_diamond_taxa(f[idx[field]])
                for rank, field in DIAMOND_RANK_FIELDS.items()
            }
            if any(lineage.values()):
                any_taxonomy = True
            hits_by_orf[f[idx["qseqid"]]].append(Hit(
                sseqid=f[idx["sseqid"]],
                pident=float(f[idx["pident"]]),
                length=int(f[idx["length"]]),
                evalue=float(f[idx["evalue"]]),
                bitscore=float(f[idx["bitscore"]]),
                lineage=lineage,
            ))
    if hits_by_orf and not any_taxonomy:
        log.warning(
            "DIAMOND returned protein hits but no rank-resolved taxonomy. Confirm that "
            "-d points to nr-tax.dmnd built with --taxonmap, --taxonnodes, and "
            "--taxonnames, and that DIAMOND is version 2.1.17 or newer."
        )
    return hits_by_orf


# --------------------------------------------------------------------------------------
# Step 3: Taxonomy is read directly from nr-tax.dmnd DIAMOND output fields.
# --------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------
# Step 4: Per-contig classification (bitscore-weighted majority/plurality vote per rank)
# --------------------------------------------------------------------------------------

def classify_contig(orf_ids, hits_by_orf: dict, ranks,
                     bitscore_range: float = 0.9, max_hits_per_orf: int = 5,
                     min_support: float = 0.5) -> dict:
    rank_weights = {r: defaultdict(float) for r in ranks}
    total_weight = {r: 0.0 for r in ranks}
    n_orfs_with_hits = 0

    for orf_id in orf_ids:
        hits = hits_by_orf.get(orf_id)
        if not hits:
            continue
        hits_sorted = sorted(hits, key=lambda h: h.bitscore, reverse=True)
        best_bs = hits_sorted[0].bitscore
        if best_bs <= 0:
            continue
        kept = [h for h in hits_sorted if h.bitscore >= bitscore_range * best_bs][:max_hits_per_orf]
        if kept:
            n_orfs_with_hits += 1
        for h in kept:
            for r in ranks:
                total_weight[r] += h.bitscore
                names = h.lineage.get(r, ())
                if names:
                    per_name_weight = h.bitscore / len(names)
                    for name in names:
                        rank_weights[r][name] += per_name_weight

    result = {}
    for r in ranks:
        weights = rank_weights[r]
        tw = total_weight[r]
        if tw <= 0 or not weights:
            result[r] = ("Unclassified", 0.0)
            continue
        top_taxon, top_w = max(weights.items(), key=lambda kv: kv[1])
        support = top_w / tw
        result[r] = (top_taxon, round(support, 4)) if support >= min_support else ("Unclassified", round(support, 4))
    result["_n_orfs_total"] = len(orf_ids)
    result["_n_orfs_with_hits"] = n_orfs_with_hits
    return result


def classify_all_contigs(contig_ids, orf_to_contig, hits_by_orf, ranks,
                          bitscore_range, max_hits_per_orf, min_support) -> dict:
    contig_to_orfs = defaultdict(list)
    for orf_id, contig_id in orf_to_contig.items():
        contig_to_orfs[contig_id].append(orf_id)

    classifications = {}
    for i, contig_id in enumerate(contig_ids, 1):
        orfs = contig_to_orfs.get(contig_id, [])
        classifications[contig_id] = classify_contig(
            orfs, hits_by_orf, ranks, bitscore_range, max_hits_per_orf, min_support
        )
        if i % 500 == 0:
            log.info("Classified %d/%d contigs...", i, len(contig_ids))
    return classifications


def write_classification_table(classifications: dict, ranks, out_path: Path,
                                excluded_ids: set = None, contig_metrics: dict = None,
                                triage_calls: dict = None, triage_excluded_ids: set = None,
                                refinement_decisions: dict = None) -> None:
    excluded_ids = excluded_ids or set()
    contig_metrics = contig_metrics or {}
    triage_calls = triage_calls or {}
    triage_excluded_ids = triage_excluded_ids or set()
    refinement_decisions = refinement_decisions or {}
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        header = [
            "contig", "length_bp", "GC_percent", "coding_density", "n_cds",
            "mean_depth", "covered_fraction", "triage_call", "excluded_gene_density",
            "refinement_status", "n_orfs", "n_orfs_with_hits",
        ]
        for r in ranks:
            header += [r, f"{r}_support"]
        header.append("excluded_animal_plant")
        w.writerow(header)
        for contig_id, res in classifications.items():
            metric = contig_metrics.get(contig_id, {})
            decision = refinement_decisions.get(contig_id, {})
            mean_depth = metric.get("mean_depth")
            covered_fraction = metric.get("covered_fraction")
            row = [
                contig_id,
                metric.get("length_bp", "NA"),
                (f"{100.0 * metric.get('gc_fraction', 0.0):.3f}"
                 if metric else "NA"),
                (f"{metric.get('coding_density', 0.0):.5f}" if metric else "NA"),
                metric.get("n_cds", res["_n_orfs_total"]),
                (f"{mean_depth:.5f}" if mean_depth is not None else "NA"),
                (f"{covered_fraction:.5f}" if covered_fraction is not None else "NA"),
                triage_calls.get(contig_id, "not_evaluated"),
                contig_id in triage_excluded_ids,
                decision.get("status", "not_evaluated"),
                res["_n_orfs_total"],
                res["_n_orfs_with_hits"],
            ]
            for r in ranks:
                taxon, support = res[r]
                row += [taxon, support]
            row.append(contig_id in excluded_ids)
            w.writerow(row)


def split_excluded_eukaryotes(classifications: dict, exclude_kingdoms) -> tuple:
    """Splits off contigs whose (always-computed) 'kingdom' classification confidently
    matches one of exclude_kingdoms -- by default Metazoa/Viridiplantae, i.e. host-animal
    or plant contamination. Bacteria, Archaea, Fungi, and protist lineages (which mostly
    have no formal kingdom-rank ancestor in NCBI taxonomy, so 'kingdom' comes back
    Unclassified for them) are never matched here and always pass through.

    Returns (kept_classifications, excluded_contig_ids).
    """
    if not exclude_kingdoms:
        return classifications, set()
    exclude_set = set(exclude_kingdoms)
    kept = {}
    excluded = set()
    for contig_id, res in classifications.items():
        kingdom, _support = res.get("kingdom", ("Unclassified", 0.0))
        if kingdom in exclude_set:
            excluded.add(contig_id)
        else:
            kept[contig_id] = res
    return kept, excluded


def _median_absolute_deviation(values, center=None) -> float:
    if not values:
        return 0.0
    center = statistics.median(values) if center is None else center
    return statistics.median(abs(value - center) for value in values)


def choose_refinement_rank(ranks, requested: str = "auto") -> str:
    """Choose a genome-oriented rank for GC/coverage coherence checks."""
    if requested and requested != "auto":
        if requested not in ranks:
            raise ValueError(
                f"--refinement-rank '{requested}' must occur in --ranks "
                f"(available: {','.join(ranks)})."
            )
        return requested
    # Prefer species so distinct species within one genus are not mistaken for GC/depth
    # outliers. Fall back toward broader ranks only when finer calls were not requested.
    for rank in ("species", "genus", "family", "phylum", "kingdom", "superkingdom", "domain"):
        if rank in ranks:
            return rank
    return ranks[-1]


def refine_taxonomic_bins(classifications: dict, contig_metrics: dict, rank: str,
                          enabled: bool = True) -> tuple:
    """Detect strong within-taxon GC *and* coverage outliers.

    Robust bin centers are learned only from long, strongly classified contigs. A contig
    is called discordant only when both GC and log2-depth exceed conservative, MAD-based
    thresholds (with absolute floors). This deliberately avoids treating strain-level GC
    variation, repeats, plasmids, or coverage noise as enough evidence on their own.
    """
    decisions = {}
    groups = defaultdict(list)
    for contig_id, res in classifications.items():
        taxon, support = res.get(rank, ("Unclassified", 0.0))
        decision = {
            "rank": rank,
            "taxon": taxon,
            "support": support,
            "status": "not_evaluated",
            "gc_outlier": False,
            "coverage_outlier": False,
            "joint_outlier": False,
        }
        decisions[contig_id] = decision
        if taxon == "Unclassified":
            decision["status"] = "unclassified_not_refined"
        else:
            groups[taxon].append(contig_id)

    outliers = set()
    for taxon, contig_ids in groups.items():
        references = [
            contig_id for contig_id in contig_ids
            if contig_metrics.get(contig_id, {}).get("length_bp", 0)
            >= REFINEMENT_MIN_CONTIG_LENGTH_BP
            and classifications[contig_id][rank][1] >= REFINEMENT_MIN_TAXON_SUPPORT
            and contig_metrics.get(contig_id, {}).get("mean_depth") is not None
        ]
        if not enabled:
            for contig_id in contig_ids:
                decisions[contig_id]["status"] = "disabled"
            continue
        if len(references) < REFINEMENT_MIN_REFERENCE_CONTIGS:
            for contig_id in contig_ids:
                decisions[contig_id]["status"] = "insufficient_reference_contigs"
            continue

        gc_values = [contig_metrics[contig_id]["gc_fraction"] for contig_id in references]
        depth_values = [contig_metrics[contig_id]["mean_depth"] for contig_id in references]
        positive_depths = [value for value in depth_values if value > 0]
        if not positive_depths:
            for contig_id in contig_ids:
                decisions[contig_id]["status"] = "no_coverage_signal"
            continue
        depth_floor = max(0.01, statistics.median(positive_depths) * 0.01)
        log_depth_values = [math.log2(value + depth_floor) for value in depth_values]
        gc_center = statistics.median(gc_values)
        depth_center = statistics.median(log_depth_values)
        gc_threshold = max(
            REFINEMENT_GC_ABSOLUTE_FLOOR,
            REFINEMENT_ROBUST_Z * 1.4826 * _median_absolute_deviation(gc_values, gc_center),
        )
        depth_threshold = max(
            REFINEMENT_COVERAGE_LOG2_FLOOR,
            REFINEMENT_ROBUST_Z * 1.4826
            * _median_absolute_deviation(log_depth_values, depth_center),
        )

        for contig_id in contig_ids:
            metric = contig_metrics.get(contig_id, {})
            decision = decisions[contig_id]
            decision.update({
                "gc_center": gc_center,
                "log2_depth_center": depth_center,
                "gc_threshold": gc_threshold,
                "log2_depth_threshold": depth_threshold,
            })
            if metric.get("length_bp", 0) < REFINEMENT_MIN_CONTIG_LENGTH_BP:
                decision["status"] = "short_contig_not_refined"
                continue
            depth = metric.get("mean_depth")
            if depth is None:
                decision["status"] = "coverage_unavailable"
                continue
            gc_delta = abs(metric.get("gc_fraction", 0.0) - gc_center)
            depth_delta = abs(math.log2(depth + depth_floor) - depth_center)
            decision["gc_delta"] = gc_delta
            decision["log2_depth_delta"] = depth_delta
            decision["gc_outlier"] = gc_delta > gc_threshold
            decision["coverage_outlier"] = depth_delta > depth_threshold
            decision["joint_outlier"] = (
                decision["gc_outlier"] and decision["coverage_outlier"]
            )
            if decision["joint_outlier"]:
                decision["status"] = "joint_gc_coverage_outlier"
                outliers.add(contig_id)
            elif decision["gc_outlier"]:
                decision["status"] = "gc_only_outlier_retained"
            elif decision["coverage_outlier"]:
                decision["status"] = "coverage_only_outlier_retained"
            else:
                decision["status"] = "coherent"
    return decisions, outliers


def write_refinement_table(decisions: dict, contig_metrics: dict, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "contig", "rank", "taxon", "taxon_support", "length_bp", "GC_percent",
        "mean_depth", "covered_fraction", "GC_center_percent", "GC_delta_percent",
        "GC_threshold_percent", "log2_depth_center", "log2_depth_delta",
        "log2_depth_threshold", "GC_outlier", "coverage_outlier", "joint_outlier",
        "status",
    ]
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(fields)
        for contig_id, decision in decisions.items():
            metric = contig_metrics.get(contig_id, {})
            writer.writerow([
                contig_id, decision.get("rank", "NA"), decision.get("taxon", "NA"),
                decision.get("support", "NA"), metric.get("length_bp", "NA"),
                (f"{100.0 * metric.get('gc_fraction', 0.0):.3f}" if metric else "NA"),
                (f"{metric['mean_depth']:.5f}"
                 if metric.get("mean_depth") is not None else "NA"),
                (f"{metric['covered_fraction']:.5f}"
                 if metric.get("covered_fraction") is not None else "NA"),
                (f"{100.0 * decision['gc_center']:.3f}"
                 if "gc_center" in decision else "NA"),
                (f"{100.0 * decision['gc_delta']:.3f}"
                 if "gc_delta" in decision else "NA"),
                (f"{100.0 * decision['gc_threshold']:.3f}"
                 if "gc_threshold" in decision else "NA"),
                (f"{decision['log2_depth_center']:.5f}"
                 if "log2_depth_center" in decision else "NA"),
                (f"{decision['log2_depth_delta']:.5f}"
                 if "log2_depth_delta" in decision else "NA"),
                (f"{decision['log2_depth_threshold']:.5f}"
                 if "log2_depth_threshold" in decision else "NA"),
                decision.get("gc_outlier", False),
                decision.get("coverage_outlier", False),
                decision.get("joint_outlier", False), decision.get("status", "NA"),
            ])


def demote_refinement_outliers(classifications: dict, outlier_ids: set, rank: str) -> dict:
    """Demote a final joint outlier at ``rank`` and its more-specific descendants."""
    refined = {contig_id: result.copy() for contig_id, result in classifications.items()}
    try:
        start = REFINEMENT_RANK_ORDER.index(rank)
    except ValueError:
        return refined
    demoted_ranks = set(REFINEMENT_RANK_ORDER[start:])
    for contig_id in outlier_ids:
        if contig_id not in refined:
            continue
        for candidate_rank in demoted_ranks:
            if candidate_rank in refined[contig_id]:
                _taxon, support = refined[contig_id][candidate_rank]
                refined[contig_id][candidate_rank] = ("Unclassified", support)
    return refined


# --------------------------------------------------------------------------------------
# Step 5: Binning
# --------------------------------------------------------------------------------------

def bin_contigs(contig_seqs: dict, classifications: dict, ranks, outdir: Path,
                 include_unclassified: bool = True) -> dict:
    """Writes <outdir>/<rank>/<taxon>.fasta for every rank in `ranks`.
    Returns {rank: {bin_name: [contig_ids]}}."""
    membership = {}
    for r in ranks:
        rank_dir = outdir / r
        rank_dir.mkdir(parents=True, exist_ok=True)
        bins = defaultdict(list)
        for contig_id, res in classifications.items():
            taxon, _support = res[r]
            if taxon == "Unclassified" and not include_unclassified:
                continue
            bins[sanitize(taxon)].append(contig_id)
        for bin_name, contig_ids in bins.items():
            records = {cid: contig_seqs[cid] for cid in contig_ids if cid in contig_seqs}
            write_fasta(rank_dir / f"{bin_name}.fasta", records)
        membership[r] = bins
        with open(rank_dir / "bin_membership.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["bin", "contig"])
            for bin_name, contig_ids in bins.items():
                for cid in contig_ids:
                    w.writerow([bin_name, cid])
        log.info("Rank '%s': wrote %d bin FASTA files to %s", r, len(bins), rank_dir)
    return membership


def filter_small_bins(contig_seqs: dict, membership: dict, ranks, bins_outdir: Path,
                       min_bin_contigs: int, min_bin_length: int,
                       include_unclassified: bool) -> dict:
    """Merges bins smaller than the given thresholds into Unclassified (or drops them
    entirely if include_unclassified is False), rewriting the bin FASTAs/membership.tsv
    in place. Returns the updated {rank: {bin_name: [contig_ids]}}."""
    if min_bin_contigs <= 1 and min_bin_length <= 0:
        return membership

    for r in ranks:
        rank_dir = bins_outdir / r
        bins = membership[r]
        kept = {}
        merged = []
        n_dropped_bins = 0
        for bin_name, contig_ids in bins.items():
            if bin_name == "Unclassified":
                kept[bin_name] = list(contig_ids)
                continue
            total_len = sum(len(contig_seqs[c]) for c in contig_ids if c in contig_seqs)
            if len(contig_ids) < min_bin_contigs or total_len < min_bin_length:
                (rank_dir / f"{bin_name}.fasta").unlink(missing_ok=True)
                merged.extend(contig_ids)
                n_dropped_bins += 1
            else:
                kept[bin_name] = list(contig_ids)

        if merged and include_unclassified:
            kept.setdefault("Unclassified", [])
            kept["Unclassified"].extend(merged)
            records = {cid: contig_seqs[cid] for cid in kept["Unclassified"] if cid in contig_seqs}
            write_fasta(rank_dir / "Unclassified.fasta", records)

        with open(rank_dir / "bin_membership.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["bin", "contig"])
            for bin_name, contig_ids in kept.items():
                for cid in contig_ids:
                    w.writerow([bin_name, cid])

        log.info(
            "Rank '%s': merged %d small/short bins (%d contigs) -> %d bins remain.",
            r, n_dropped_bins, len(merged), len(kept),
        )
        membership[r] = kept
    return membership


# --------------------------------------------------------------------------------------
# Step 6: Per-bin-set summaries (QUAST + CheckM)
# --------------------------------------------------------------------------------------

def basic_assembly_stats(fasta_path: Path) -> dict:
    seqs = read_fasta(fasta_path)
    lengths = sorted((len(s) for s in seqs.values()), reverse=True)
    total = sum(lengths)
    gc_count = sum(s.upper().count("G") + s.upper().count("C") for s in seqs.values())
    gc_pct = round(100.0 * gc_count / total, 2) if total else 0.0
    n50 = l50 = 0
    cum = 0
    for i, l in enumerate(lengths, 1):
        cum += l
        if cum >= total / 2 and n50 == 0:
            n50, l50 = l, i
    return {
        "num_contigs": len(lengths),
        "total_length_bp": total,
        "largest_contig_bp": lengths[0] if lengths else 0,
        "N50": n50,
        "L50": l50,
        "GC_percent": gc_pct,
    }


def run_quast_multi(bin_fastas: dict, outdir: Path, threads: int, quast_cmd: str = "quast.py") -> dict:
    """bin_fastas: {bin_name: Path}. Returns {bin_name: {stat: value}}."""
    if not bin_fastas:
        return {}
    quast_argv = shlex.split(quast_cmd)
    if shutil.which(quast_argv[0]) is None:
        log.warning("%s not found on PATH -- falling back to built-in assembly stats.", quast_cmd)
        return {name: basic_assembly_stats(p) for name, p in bin_fastas.items()}

    outdir.mkdir(parents=True, exist_ok=True)
    names = list(bin_fastas.keys())
    paths = [str(bin_fastas[n]) for n in names]
    cmd = [
        *quast_argv, "-o", str(outdir), "--threads", str(threads),
        "--min-contig", "0", "--silent", "--labels", ",".join(names), *paths,
    ]
    try:
        run_cmd(cmd, log_file=outdir / "quast.log")
    except RuntimeError as exc:
        log.warning("QUAST failed (%s); falling back to built-in assembly stats.", exc)
        return {name: basic_assembly_stats(p) for name, p in bin_fastas.items()}

    report = outdir / "transposed_report.tsv"
    if not report.exists():
        log.warning("QUAST did not produce transposed_report.tsv; falling back to built-in stats.")
        return {name: basic_assembly_stats(p) for name, p in bin_fastas.items()}

    stats = {}
    with open(report) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        col = {h: i for i, h in enumerate(header)}
        for row in reader:
            name = row[col["Assembly"]]

            def g(key, cast=float, default=0):
                for k in col:
                    if k.startswith(key):
                        try:
                            return cast(row[col[k]])
                        except ValueError:
                            return default
                return default

            stats[name] = {
                "num_contigs": g("# contigs (>= 0", int),
                "total_length_bp": g("Total length (>= 0", int),
                "largest_contig_bp": g("Largest contig", int),
                "N50": g("N50", int),
                "L50": g("L50", int),
                "GC_percent": g("GC (%)", float),
            }
    return stats


def run_checkm(bin_dir: Path, outdir: Path, threads: int, checkm_cmd: str = "checkm") -> dict:
    """Runs `checkm lineage_wf`, reporting only completeness/contamination per bin (no
    plots, no extended stats -- just `-x fasta --tab_table -f <results.tsv>`, same as
    the minimal `checkm lineage_wf -t T -x fasta <bins_dir> <out_dir>` invocation).
    Returns {bin_name: {stat: value}}."""
    checkm_argv = shlex.split(checkm_cmd)
    if shutil.which(checkm_argv[0]) is None:
        log.warning("%s not found on PATH -- skipping completeness/contamination (will be NA).", checkm_cmd)
        return {}
    fastas = list(bin_dir.glob("*.fasta"))
    if not fastas:
        return {}
    outdir.mkdir(parents=True, exist_ok=True)
    results_tsv = outdir / "checkm_results.tsv"
    cmd = [
        *checkm_argv, "lineage_wf", "-x", "fasta", "--tab_table", "-f", str(results_tsv),
        "-t", str(threads), str(bin_dir), str(outdir),
    ]
    try:
        run_cmd(cmd, log_file=outdir / "checkm.log")
    except RuntimeError as exc:
        log.warning("CheckM failed (%s); completeness/contamination will be NA for this bin set.", exc)
        return {}

    stats = {}
    if results_tsv.exists():
        with open(results_tsv) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                bin_id = row.get("Bin Id")
                if not bin_id:
                    continue
                stats[bin_id] = {
                    "completeness_percent": row.get("Completeness"),
                    "contamination_percent": row.get("Contamination"),
                }
    return stats


def summarize_bin_set(rank: str, rank_dir: Path, threads: int, skip_quast: bool,
                       skip_checkm: bool, quast_cmd: str = "quast.py",
                       checkm_cmd: str = "checkm") -> None:
    bin_fastas = {p.stem: p for p in sorted(rank_dir.glob("*.fasta"))}
    if not bin_fastas:
        log.warning("No bins found for rank '%s'; skipping summary.", rank)
        return

    if skip_quast:
        quast_stats = {name: basic_assembly_stats(p) for name, p in bin_fastas.items()}
    else:
        quast_stats = run_quast_multi(bin_fastas, rank_dir / "quast_out", threads, quast_cmd)

    checkm_stats = {} if skip_checkm else run_checkm(rank_dir, rank_dir / "checkm_out", threads, checkm_cmd)

    out_path = rank_dir / "summary.tsv"
    fields = ["bin", "num_contigs", "total_length_bp", "largest_contig_bp", "N50", "L50",
              "GC_percent", "completeness_percent", "contamination_percent"]
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(fields)
        for name in bin_fastas:
            q = quast_stats.get(name, {})
            c = checkm_stats.get(name, {})
            w.writerow([
                name,
                q.get("num_contigs", "NA"),
                q.get("total_length_bp", "NA"),
                q.get("largest_contig_bp", "NA"),
                q.get("N50", "NA"),
                q.get("L50", "NA"),
                q.get("GC_percent", "NA"),
                c.get("completeness_percent", "NA"),
                c.get("contamination_percent", "NA"),
            ])
    log.info("Wrote %s", out_path)


# --------------------------------------------------------------------------------------
# Step 7: Targeted bin reassembly -- competitive seed mapping + ITSME-inspired BBDuk
# frontier extension + focused SPAdes/metaSPAdes + Unicycler/Pilon refinement.
# --------------------------------------------------------------------------------------
#
# Whole-metagenome co-assembly (MEGAHIT, Step 0c) fragments low-abundance or
# fast-diverging genomes -- especially host-restricted symbionts (Wolbachia, Rickettsia,
# etc.) that share conserved genes/k-mers with other community members and rarely have a
# close enough reference genome to assemble against directly. Rather than requiring a
# reference, this step combines every mutually exclusive bin at a selected rank into one
# seed index. Reads map competitively: only a unique best-scoring bin wins the template,
# while equal best-score ties are excluded. Winning read pools are optionally extended
# outward with cheap exact-kmer "frontier" scans (BBDuk) to catch divergent regions the
# initial contigs missed entirely, and the resulting small, mostly-single-organism read
# pool is reassembled alone with SPAdes/metaSPAdes. A subset-only assembly graph is far
# simpler than the whole-community graph, so it can often resolve tangles (shared genes,
# similar-coverage strains) that fragmented the original bin. Frontier extension is
# guarded by growth-rate/accepted-fraction stop conditions so recruitment cannot snowball.
# Optional external anchors are named BIN=FASTA so they participate for only their intended
# competitor. Unicycler then attempts graph circularization and Pilon polishes the result.
#
# This does NOT replace the original bin FASTA. It writes a separate reassembled contig
# set plus a before/after comparison table so you can judge whether it actually helped.

def build_bowtie2_index(fasta: Path, index_prefix: Path, threads: int, log_file: Path) -> None:
    run_cmd(["bowtie2-build", "--threads", str(threads), str(fasta), str(index_prefix)], log_file=log_file)


def map_reads_to_index(r1: Path, r2: Path, index_prefix: Path, out_bam: Path, threads: int,
                        score_min: str, max_insert: int, log_file: Path,
                        report_multiple: int = 1) -> None:
    """Maps ALL raw read pairs (mapped and unmapped alike -- no --no-unal) against
    index_prefix. Keeping unmapped records in the BAM is what lets extract_templates()
    later pull out reads recruited purely by frontier k-mer matching, not just direct
    alignment, using this same BAM."""
    bt2 = ["bowtie2", "--very-sensitive-local", "--score-min", score_min,
           "-X", str(max_insert), "-p", str(threads), "-x", str(index_prefix),
           "-1", str(r1), "-2", str(r2)]
    if report_multiple > 1:
        bt2 += ["-k", str(report_multiple)]
    sam = ["samtools", "view", "-@", str(threads), "-b", "-o", str(out_bam), "-"]
    run_pipeline([bt2, sam], log_file=log_file)


def estimate_contig_coverage(fasta: Path, r1: Path, r2: Path, outdir: Path,
                             threads: int, max_insert: int = 1000) -> dict:
    """Competitively map the complete read set and return per-contig mean depth.

    Bowtie2 reports one best placement by default, so a multi-mapping read does not add
    depth to every similar contig. ``samtools depth -aa`` is streamed and aggregated to
    avoid materializing a potentially very large depth table.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    index_prefix = outdir / "contigs_index"
    build_bowtie2_index(fasta, index_prefix, threads, outdir / "bowtie2-build.log")
    sorted_bam = outdir / "reads_to_contigs.sorted.bam"
    bt2 = [
        "bowtie2", "--very-sensitive-local", "--no-unal", "-X", str(max_insert),
        "-p", str(threads), "-x", str(index_prefix), "-1", str(r1), "-2", str(r2),
    ]
    view = ["samtools", "view", "-@", str(threads), "-b", "-F", "4", "-"]
    sort = ["samtools", "sort", "-@", str(threads), "-o", str(sorted_bam), "-"]
    run_pipeline([bt2, view, sort], log_file=outdir / "bowtie2.log")

    lengths = {contig_id: len(sequence) for contig_id, sequence in read_fasta(fasta).items()}
    depth_sum = defaultdict(float)
    covered = defaultdict(int)
    proc = subprocess.Popen(
        ["samtools", "depth", "-aa", str(sorted_bam)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3:
            continue
        try:
            depth = float(fields[2])
        except ValueError:
            continue
        depth_sum[fields[0]] += depth
        if depth > 0:
            covered[fields[0]] += 1
    stderr = proc.stderr.read() if proc.stderr is not None else ""
    returncode = proc.wait()
    if returncode != 0:
        raise RuntimeError(f"samtools depth failed ({returncode}): {stderr.strip()}")

    coverage = {}
    for contig_id, length in lengths.items():
        coverage[contig_id] = {
            "mean_depth": (depth_sum[contig_id] / length) if length else 0.0,
            "covered_fraction": (covered[contig_id] / length) if length else 0.0,
        }
    return coverage


def add_coverage_to_metrics(contig_metrics: dict, coverage: dict) -> None:
    for contig_id, values in coverage.items():
        if contig_id in contig_metrics:
            contig_metrics[contig_id].update(values)


def aligned_bases_from_cigar(cigar: str) -> int:
    return sum(int(n) for n, _op in re.findall(r"(\d+)([MI=X])", cigar))


def query_bases_from_cigar(cigar: str) -> int:
    """Query length represented by a CIGAR, including soft-clipped sequence."""
    return sum(int(n) for n, _op in re.findall(r"(\d+)([MIS=X])", cigar))


def build_competitive_seed_reference(bin_fastas, out_fasta: Path,
                                     anchors_by_bin: dict = None,
                                     excluded_contig_ids: set = None) -> dict:
    """Combine all rank bins into one reference and return ``reference -> bin``.

    Prefixing every record with a generated ID makes the owning bin unambiguous even
    when input FASTAs reuse contig names. Named external anchors, when present, join only
    their specified bin instead of being duplicated across every competitor.
    """
    anchors_by_bin = anchors_by_bin or {}
    excluded_contig_ids = excluded_contig_ids or set()
    ref_to_bin = {}
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(out_fasta, "w") as out_fh:
        for bin_n, bin_fasta in enumerate(bin_fastas, 1):
            bin_name = bin_fasta.stem
            sources = [("contig", bin_fasta)]
            sources.extend(("anchor", p) for p in anchors_by_bin.get(bin_name, []))
            record_n = 0
            for source_kind, source_fasta in sources:
                for original_id, seq in read_fasta(source_fasta).items():
                    if source_kind == "contig" and original_id in excluded_contig_ids:
                        continue
                    record_n += 1
                    ref_id = f"MHSEED_{bin_n:06d}_{record_n:09d}_{source_kind}"
                    ref_to_bin[ref_id] = bin_name
                    out_fh.write(f">{ref_id}\n")
                    for i in range(0, len(seq), 80):
                        out_fh.write(seq[i:i + 80] + "\n")
    if not ref_to_bin:
        raise RuntimeError("Competitive seed reference contains no sequences.")
    return ref_to_bin


def competitively_assign_templates(bam: Path, ref_to_bin: dict, min_aligned_fraction: float,
                                   min_identity: float, threads: int) -> tuple:
    """Assign each template to its unique best-scoring bin across the combined seed index.

    Secondary alignments are retained so conserved reads can expose cross-bin competition.
    The best AS score for each mate is summed within each candidate bin. Equal best scores
    across bins are called ambiguous and excluded from every bin.
    """
    proc = subprocess.run(
        ["samtools", "view", "-@", str(threads), "-F", "2052", str(bam)],
        capture_output=True, text=True, check=True,
    )
    per_template = defaultdict(lambda: defaultdict(dict))
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 11 or f[2] not in ref_to_bin or f[5] == "*":
            continue
        flag = int(f[1])
        aligned = aligned_bases_from_cigar(f[5])
        query_len = query_bases_from_cigar(f[5]) or len(f[9])
        nm = 0
        alignment_score = None
        for tag in f[11:]:
            if tag.startswith("NM:i:"):
                nm = int(tag.split(":", 2)[2])
            elif tag.startswith("AS:i:"):
                alignment_score = int(tag.split(":", 2)[2])
        if not query_len or not aligned or alignment_score is None:
            continue
        if aligned / query_len < min_aligned_fraction:
            continue
        if (aligned - nm) / aligned < min_identity:
            continue
        mate = 1 if flag & 64 else (2 if flag & 128 else 0)
        bin_name = ref_to_bin[f[2]]
        prior = per_template[f[0]][bin_name].get(mate)
        if prior is None or alignment_score > prior:
            per_template[f[0]][bin_name][mate] = alignment_score

    assignments = defaultdict(set)
    ambiguous = 0
    for template, bin_mates in per_template.items():
        scores = {
            bin_name: sum(mate_scores.values())
            for bin_name, mate_scores in bin_mates.items()
        }
        if not scores:
            continue
        best_score = max(scores.values())
        winners = [bin_name for bin_name, score in scores.items() if score == best_score]
        if len(winners) != 1:
            ambiguous += 1
            continue
        assignments[winners[0]].add(template)
    return {bin_name: names for bin_name, names in assignments.items()}, len(per_template), ambiguous


def recruit_read_names(bam: Path, min_aligned_fraction: float, min_identity: float, threads: int) -> set:
    """Primary, mapped (not secondary/supplementary) alignments only (-F 2308), filtered by
    aligned-fraction-of-read and approximate identity -- the same seed-hit criteria ITSME
    uses for its strict seed mapping."""
    proc = subprocess.run(["samtools", "view", "-@", str(threads), "-F", "2308", str(bam)],
                           capture_output=True, text=True, check=True)
    names = set()
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 11:
            continue
        seq, cigar = f[9], f[5]
        query_len = len(seq)
        if query_len == 0 or cigar == "*":
            continue
        aligned = aligned_bases_from_cigar(cigar)
        nm = 0
        for tag in f[11:]:
            if tag.startswith("NM:i:"):
                nm = int(tag.split(":", 2)[2])
                break
        aligned_fraction = aligned / query_len if query_len else 0.0
        identity = (aligned - nm) / aligned if aligned else 0.0
        if aligned_fraction >= min_aligned_fraction and identity >= min_identity:
            names.add(f[0])
    return names


def extract_templates(bam: Path, names: set, out_r1: Path, out_r2: Path, out_single: Path,
                       threads: int, tmp_prefix: Path) -> None:
    """Pulls the given read names (+ their mates, however the mate mapped) out of `bam`
    into paired/single FASTQs -- a Python port of ITSME's extract_templates_from_bam."""
    out_r1.parent.mkdir(parents=True, exist_ok=True)
    if not names:
        for p in (out_r1, out_r2, out_single):
            p.write_bytes(b"")
        return
    names_file = Path(f"{tmp_prefix}.names.txt")
    names_file.write_text("\n".join(sorted(names)) + "\n")
    selected_bam = Path(f"{tmp_prefix}.selected.bam")
    run_cmd(["samtools", "view", "-@", str(threads), "-b", "-F", "2304", "-N", str(names_file),
             "-o", str(selected_bam), str(bam)])
    other_fq = Path(f"{tmp_prefix}.other.fastq.gz")
    single_fq = Path(f"{tmp_prefix}.singleton.fastq.gz")
    collate = ["samtools", "collate", "-@", str(threads), "-u", "-O", str(selected_bam)]
    fastq = ["samtools", "fastq", "-@", str(threads), "-c", "1", "-n",
             "-1", str(out_r1), "-2", str(out_r2), "-0", str(other_fq), "-s", str(single_fq), "-"]
    run_pipeline([collate, fastq], log_file=Path(f"{tmp_prefix}.fastq.log"))
    with open(out_single, "wb") as out_fh:
        for p in (other_fq, single_fq):
            if p.exists() and p.stat().st_size > 0:
                out_fh.write(p.read_bytes())
    if not out_single.exists():
        out_single.write_bytes(b"")
    for p in (selected_bam, other_fq, single_fq, names_file):
        p.unlink(missing_ok=True)


def make_frontier_baits(fastq_paths, out_fasta: Path, word_size: int, entropy: float,
                         bbtools_memory: str, log_file: Path) -> None:
    """Converts a pool of recruited reads into a bait FASTA (dropping low-complexity
    sequence via BBDuk's entropy filter) to seed the next round's exact-kmer frontier scan."""
    raw = out_fasta.with_suffix(".unfiltered.fasta")
    tag = 0
    with open(raw, "w") as out_fh:
        for fq in fastq_paths:
            if not fq or not Path(fq).exists() or Path(fq).stat().st_size == 0:
                continue
            tag += 1
            opener = gzip.open if str(fq).endswith((".gz", ".bgz")) else open
            with opener(fq, "rt") as fh:
                record = 0
                for i, line in enumerate(fh):
                    mod = i % 4
                    if mod == 0:
                        stripped = line[1:].strip()
                        record += 1
                        name = stripped.split()[0] if stripped else str(record)
                        out_fh.write(f">bait{tag:02d}|{record:09d}|{name}\n")
                    elif mod == 1:
                        out_fh.write(line.strip().upper() + "\n")
    run_cmd(["bbduk.sh", f"-Xmx{bbtools_memory}", f"in={raw}", f"out={out_fasta}",
             "overwrite=t", f"minlen={word_size}", f"entropy={entropy}",
             "entropywindow=50", "entropyk=5"], log_file=log_file)
    raw.unlink(missing_ok=True)
    if not out_fasta.exists() or out_fasta.stat().st_size == 0:
        raise RuntimeError(f"No frontier sequences survived bait preparation; inspect {log_file}")


def scan_raw_with_baits(bait_fasta: Path, r1: Path, r2: Path, word_size: int,
                         min_word_hits: int, bbtools_memory: str, threads: int,
                         log_file: Path) -> set:
    """Exact-kmer scan of the ORIGINAL raw reads against a bait FASTA (BBDuk), returning
    template names with at least `min_word_hits` exact word matches -- cheap compared to
    another bowtie2 pass, and how the frontier extends beyond direct alignment recruitment."""
    cmd = ["bbduk.sh", f"-Xmx{bbtools_memory}", f"in1={r1}", f"in2={r2}", "outm=stdout.fq",
           f"ref={bait_fasta}", f"k={word_size}", "hdist=0", f"minkmerhits={min_word_hits}",
           "mm=f", "rcomp=t", f"t={threads}", "overwrite=t", "ordered=f"]
    with open(log_file, "w") as lf:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=lf, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"BBDuk frontier scan failed ({proc.returncode}); inspect {log_file}")
    names = set()
    for i, line in enumerate(proc.stdout.splitlines()):
        if i % 4 == 0 and len(line) > 1:
            name = re.sub(r"/[12]$", "", line[1:].split()[0])
            if name:
                names.add(name)
    return names


def combine_fastqs(output: Path, inputs) -> None:
    """Concatenates gzip FASTQs by raw byte concatenation (a valid multi-member gzip stream,
    same trick ITSME/`cat *.gz` rely on -- avoids a decompress/recompress round trip)."""
    with open(output, "wb") as out_fh:
        for p in inputs:
            if p and Path(p).exists() and Path(p).stat().st_size > 0:
                out_fh.write(Path(p).read_bytes())
    if not output.exists():
        output.write_bytes(b"")


def run_spades_targeted(r1: Path, r2: Path, single: Path, outdir: Path, threads: int,
                         memory_gb: int, mode: str, kmers: str, log_file: Path) -> Path:
    """Runs spades.py (metaSPAdes if mode == 'meta') on a small recruited read pool.
    Returns the path to contigs.fasta."""
    pairs = count_fastq_reads(r1)
    singles = count_fastq_reads(single)
    if pairs == 0 and singles == 0:
        raise RuntimeError("No recruited reads remain for reassembly.")
    cmd = ["spades.py", "--only-assembler", "-o", str(outdir), "-t", str(threads), "-m", str(memory_gb)]
    if pairs > 0:
        cmd += ["-1", str(r1), "-2", str(r2)]
    if singles > 0:
        cmd += ["-s", str(single)]
    if mode == "meta":
        if pairs == 0:
            raise RuntimeError("metaSPAdes (--reassemble-mode meta) requires paired reads.")
        cmd += ["--meta"]
    if kmers != "auto":
        cmd += ["-k", kmers]
    run_cmd(cmd, log_file=log_file)
    contigs = outdir / "contigs.fasta"
    if not contigs.exists():
        raise RuntimeError(f"SPAdes did not produce {contigs}; inspect {log_file}")
    return contigs


def run_unicycler_and_pilon(r1: Path, r2: Path, single: Path, workdir: Path,
                            threads: int, unicycler_mode: str) -> tuple:
    """Attempt graph circularization with Unicycler, then polish with Pilon.

    Unicycler performs the circularization attempt from the recruited read pool. Pilon
    does not circularize; it remaps the same reads to the Unicycler assembly and corrects
    bases/indels. The caller can retain the SPAdes assembly if this optional refinement
    fails, and can retain the Unicycler assembly if only Pilon fails.
    """
    pairs = count_fastq_reads(r1)
    singles = count_fastq_reads(single)
    if pairs == 0 and singles == 0:
        return None, "spades", 0

    unicycler_dir = workdir / "unicycler"
    cmd = [
        "unicycler", "-o", str(unicycler_dir), "--threads", str(threads),
        "--mode", unicycler_mode,
    ]
    if pairs:
        cmd += ["-1", str(r1), "-2", str(r2)]
    if singles:
        cmd += ["-s", str(single)]
    try:
        run_cmd(cmd, log_file=workdir / "unicycler.log")
    except RuntimeError as exc:
        log.warning("Unicycler circularization failed (%s); retaining focused SPAdes output.", exc)
        return None, "spades", 0

    unicycler_assembly = unicycler_dir / "assembly.fasta"
    if not unicycler_assembly.exists() or unicycler_assembly.stat().st_size == 0:
        log.warning("Unicycler produced no assembly; retaining focused SPAdes output.")
        return None, "spades", 0
    circular_count = sum(
        1 for line in unicycler_assembly.read_text().splitlines()
        if line.startswith(">") and "circular=true" in line.lower()
    )

    polish_dir = workdir / "pilon"
    polish_dir.mkdir(parents=True, exist_ok=True)
    index_prefix = polish_dir / "unicycler_index"
    bam = polish_dir / "reads_to_unicycler.sorted.bam"
    try:
        build_bowtie2_index(
            unicycler_assembly, index_prefix, threads, polish_dir / "bowtie2-build.log",
        )
        bt2 = [
            "bowtie2", "--very-sensitive", "-p", str(threads), "-x", str(index_prefix),
        ]
        if pairs:
            bt2 += ["-1", str(r1), "-2", str(r2)]
        if singles:
            bt2 += ["-U", str(single)]
        run_pipeline(
            [bt2, ["samtools", "sort", "-@", str(threads), "-o", str(bam), "-"]],
            log_file=polish_dir / "bowtie2.log",
        )
        run_cmd(["samtools", "index", "-@", str(threads), str(bam)])
        run_cmd([
            "pilon", "--genome", str(unicycler_assembly), "--frags", str(bam),
            "--output", "pilon_polished", "--outdir", str(polish_dir),
            "--threads", str(threads), "--fix", "all",
        ], log_file=workdir / "pilon.log")
    except RuntimeError as exc:
        log.warning("Pilon polishing failed (%s); retaining the Unicycler assembly.", exc)
        return unicycler_assembly, "unicycler", circular_count

    polished = polish_dir / "pilon_polished.fasta"
    if not polished.exists() or polished.stat().st_size == 0:
        log.warning("Pilon produced no polished FASTA; retaining the Unicycler assembly.")
        return unicycler_assembly, "unicycler", circular_count
    return polished, "pilon", circular_count


def reassemble_one_bin(bin_fasta: Path, bin_name: str, rank: str, r1_raw: Path, r2_raw: Path,
                        competitive_bam: Path, seed_names: set, blocked_seed_names: set,
                        outdir: Path, threads: int,
                        trimmomatic_folder: Path, qc_quality: int, qc_minlen: int,
                        flash_max_overlap: int, max_rounds: int, word_size: int,
                        min_word_hits: int,
                        bait_min_entropy: float, max_round_growth: float,
                        max_accepted_fraction: float, min_new_templates: int, min_growth: float,
                        bbtools_memory: str, spades_mode: str, spades_memory_gb: int,
                        kmers: str, circularize: bool, unicycler_mode: str) -> Path:
    """Runs one bin through seed-and-extend recruitment + reassembly. Returns the path to
    the reassembled contigs FASTA, or None if recruitment/reassembly didn't produce one."""
    workdir = outdir / "reassembly" / rank / bin_name
    seed_dir, recruit_dir = workdir / "seed", workdir / "recruitment"
    for d in (seed_dir, recruit_dir):
        d.mkdir(parents=True, exist_ok=True)

    total_templates = count_fastq_reads(r1_raw)
    accepted_names = set(seed_names)
    if not accepted_names:
        log.warning("Bin '%s' (%s): no reads won competitive seed mapping; skipping reassembly.",
                    bin_name, rank)
        return None
    accepted_count = len(accepted_names)
    accepted_fraction = accepted_count / total_templates if total_templates else 0.0

    metrics_path = workdir / "recruitment.tsv"
    metrics_rows = [
        ["round", "frontier_templates", "candidate_templates", "new_templates",
         "accepted_total", "growth_fraction", "accepted_fraction", "decision"],
        [0, accepted_count, accepted_count, accepted_count, accepted_count,
         "NA", f"{accepted_fraction:.8f}", "seed"],
    ]

    seed_raw_r1 = seed_dir / "accepted_raw_R1.fastq.gz"
    seed_raw_r2 = seed_dir / "accepted_raw_R2.fastq.gz"
    seed_raw_single = seed_dir / "accepted_raw_single.fastq.gz"
    extract_templates(competitive_bam, accepted_names, seed_raw_r1, seed_raw_r2, seed_raw_single,
                       threads, seed_dir / "extract")
    log.info("Bin '%s' (%s): competitive seed mapping assigned %d/%d templates (%.4f); running QC.",
              bin_name, rank, accepted_count, total_templates, accepted_fraction)
    seed_r1, seed_r2, seed_single = run_qc(
        seed_raw_r1, seed_raw_r2, workdir / "qc" / "seed", threads,
        trimmomatic_cmd="trimmomatic", trimmomatic_folder=trimmomatic_folder,
        flash_cmd="flash", flash_max_overlap=flash_max_overlap, pigz_cmd="pigz",
        qc_quality=qc_quality, qc_minlen=qc_minlen, keep_tmp=False,
    )
    for p in (seed_raw_r1, seed_raw_r2, seed_raw_single):
        p.unlink(missing_ok=True)
    assembly_r1_files, assembly_r2_files, assembly_single_files = [seed_r1], [seed_r2], [seed_single]

    stop_reason = "extension disabled (--reassemble-max-rounds 0)"
    rounds_accepted = 0
    if max_rounds > 0:
        frontier_names = accepted_names
        frontier_baits = recruit_dir / "round_00_frontier_baits.fasta"
        make_frontier_baits([seed_r1, seed_r2, seed_single], frontier_baits, word_size,
                             bait_min_entropy, bbtools_memory, recruit_dir / "round_00_bait_filter.log")

        stop_reason = "maximum recruitment rounds reached"
        for round_n in range(1, max_rounds + 1):
            round_dir = recruit_dir / f"round_{round_n:02d}"
            round_dir.mkdir(parents=True, exist_ok=True)
            frontier_count_this_round = len(frontier_names)

            candidate_names = scan_raw_with_baits(
                frontier_baits, r1_raw, r2_raw, word_size, min_word_hits, bbtools_memory,
                threads, round_dir / "bbduk.log",
            )
            # Never let frontier extension steal a template that competitive seed mapping
            # assigned to another bin. Previously unassigned reads remain eligible to bridge
            # outward from this bin's newest frontier.
            new_names = candidate_names - accepted_names - blocked_seed_names
            growth = (len(new_names) / accepted_count) if accepted_count else 0.0
            proposed_count = accepted_count + len(new_names)
            proposed_fraction = proposed_count / total_templates if total_templates else 0.0

            if growth > max_round_growth:
                metrics_rows.append([round_n, frontier_count_this_round, len(candidate_names),
                                      len(new_names), accepted_count, f"{growth:.8f}",
                                      f"{accepted_fraction:.8f}", "rejected_growth"])
                stop_reason = f"round-{round_n} growth {growth:.4f} exceeded {max_round_growth}; recruitment rejected"
                log.warning("Bin '%s' (%s): %s.", bin_name, rank, stop_reason)
                break
            if proposed_fraction > max_accepted_fraction:
                metrics_rows.append([round_n, frontier_count_this_round, len(candidate_names),
                                      len(new_names), accepted_count, f"{growth:.8f}",
                                      f"{accepted_fraction:.8f}", "rejected_total_fraction"])
                stop_reason = (f"round-{round_n} accepted fraction {proposed_fraction:.4f} exceeded "
                                f"{max_accepted_fraction}; recruitment rejected")
                log.warning("Bin '%s' (%s): %s.", bin_name, rank, stop_reason)
                break
            if not new_names:
                metrics_rows.append([round_n, frontier_count_this_round, len(candidate_names), 0,
                                      accepted_count, "0.00000000", f"{accepted_fraction:.8f}",
                                      "converged"])
                stop_reason = "no new templates were recruited"
                break

            accepted_names = accepted_names | new_names
            accepted_count = proposed_count
            accepted_fraction = proposed_fraction
            rounds_accepted += 1

            raw_r1 = round_dir / "new_raw_R1.fastq.gz"
            raw_r2 = round_dir / "new_raw_R2.fastq.gz"
            raw_single = round_dir / "new_raw_single.fastq.gz"
            extract_templates(competitive_bam, new_names, raw_r1, raw_r2, raw_single, threads,
                               round_dir / "extract")
            log.info("Bin '%s' (%s): round %d QC of %d newly accepted templates.",
                      bin_name, rank, round_n, len(new_names))
            clean_r1, clean_r2, clean_single = run_qc(
                raw_r1, raw_r2, round_dir / "qc", threads,
                trimmomatic_cmd="trimmomatic", trimmomatic_folder=trimmomatic_folder,
                flash_cmd="flash", flash_max_overlap=flash_max_overlap, pigz_cmd="pigz",
                qc_quality=qc_quality, qc_minlen=qc_minlen, keep_tmp=False,
            )
            for p in (raw_r1, raw_r2, raw_single):
                p.unlink(missing_ok=True)
            assembly_r1_files.append(clean_r1)
            assembly_r2_files.append(clean_r2)
            assembly_single_files.append(clean_single)

            metrics_rows.append([round_n, frontier_count_this_round, len(candidate_names),
                                  len(new_names), accepted_count, f"{growth:.8f}",
                                  f"{accepted_fraction:.8f}", "accepted"])
            log.info("Bin '%s' (%s): round %d accepted %d new templates; total=%d, growth=%.4f.",
                      bin_name, rank, round_n, len(new_names), accepted_count, growth)

            frontier_names = new_names
            frontier_baits = round_dir / "frontier_baits.fasta"
            make_frontier_baits([clean_r1, clean_r2, clean_single], frontier_baits, word_size,
                                 bait_min_entropy, bbtools_memory, round_dir / "bait_filter.log")

            if len(new_names) < min_new_templates:
                stop_reason = f"new-template count {len(new_names)} fell below {min_new_templates}"
                break
            if growth < min_growth:
                stop_reason = f"growth {growth:.6f} fell below {min_growth}"
                break
    else:
        log.info("Bin '%s' (%s): frontier extension disabled; assembling seed recruitment directly.",
                  bin_name, rank)

    log.info("Bin '%s' (%s): recruitment stopped (%s) after %d extension round(s); assembling.",
              bin_name, rank, stop_reason, rounds_accepted)

    final_r1 = workdir / "accepted_R1.fastq.gz"
    final_r2 = workdir / "accepted_R2.fastq.gz"
    final_single = workdir / "accepted_single.fastq.gz"
    combine_fastqs(final_r1, assembly_r1_files)
    combine_fastqs(final_r2, assembly_r2_files)
    combine_fastqs(final_single, assembly_single_files)

    with open(metrics_path, "w", newline="") as mf:
        csv.writer(mf, delimiter="\t").writerows(metrics_rows)

    spades_dir = workdir / "spades"
    try:
        contigs = run_spades_targeted(final_r1, final_r2, final_single, spades_dir, threads,
                                       spades_memory_gb, spades_mode, kmers,
                                       workdir / "spades.log")
    except RuntimeError as exc:
        log.warning("Bin '%s' (%s): reassembly SPAdes run failed (%s); keeping original bin contigs.",
                    bin_name, rank, exc)
        return None

    chosen_contigs = contigs
    chosen_stage = "spades"
    circular_count = 0
    if circularize:
        refined, chosen_stage, circular_count = run_unicycler_and_pilon(
            final_r1, final_r2, final_single, workdir, threads, unicycler_mode,
        )
        if refined is not None:
            chosen_contigs = refined

    with open(workdir / "assembly_stage.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["selected_stage", "unicycler_circular_contigs", "selected_fasta"])
        w.writerow([chosen_stage, circular_count, chosen_contigs])

    reassembled = workdir / "reassembled.fasta"
    shutil.copyfile(chosen_contigs, reassembled)
    return reassembled


def run_bin_reassembly(outdir: Path, ranks, r1_raw: Path, r2_raw: Path, threads: int,
                        anchors_by_bin: dict, trimmomatic_folder: Path, qc_quality: int, qc_minlen: int,
                        flash_max_overlap: int, seed_score_min: str, min_read_aligned: float,
                        min_read_identity: float, max_insert: int, max_rounds: int,
                        word_size: int, min_word_hits: int, bait_min_entropy: float,
                        max_round_growth: float, max_accepted_fraction: float,
                        min_new_templates: int, min_growth: float, bbtools_memory: str,
                        spades_mode: str, spades_memory_gb: int, kmers: str,
                        min_bin_contigs_to_reassemble: int, include_unclassified: bool,
                        circularize: bool, unicycler_mode: str,
                        seed_excluded_ids: set = None) -> dict:
    """Drives reassemble_one_bin() over every bin FASTA at each rank in `ranks`, and writes
    a before/after comparison table (contig count, N50, total length) per rank.

    Returns {rank: {bin_name: reassembled_fasta}} for successful reassemblies. Bins that
    were skipped, failed, or recruited no reads are absent and can therefore fall back to
    their original preliminary-bin FASTA during final-assembly consolidation.
    """
    successful = defaultdict(dict)
    seed_excluded_ids = seed_excluded_ids or set()
    for r in ranks:
        rank_dir = outdir / "bins" / r
        all_bin_fastas = sorted(rank_dir.glob("*.fasta"))
        if not all_bin_fastas:
            continue
        bin_fastas = [
            p for p in all_bin_fastas
            if p.name != "Unclassified.fasta" or include_unclassified
        ]

        # Map once against every mutually exclusive bin at this rank. Even bins that are
        # not scheduled for reassembly remain in the index so their reads cannot be
        # spuriously awarded to a different bin.
        competitive_dir = outdir / "reassembly" / r / "competitive_seed"
        competitive_dir.mkdir(parents=True, exist_ok=True)
        combined_seed = competitive_dir / "all_bins.fasta"
        ref_to_bin = build_competitive_seed_reference(
            all_bin_fastas, combined_seed, anchors_by_bin=anchors_by_bin,
            excluded_contig_ids=seed_excluded_ids,
        )
        seed_index = competitive_dir / "all_bins_index"
        build_bowtie2_index(
            combined_seed, seed_index, threads, competitive_dir / "bowtie2-build.log",
        )
        competitive_bam = competitive_dir / "reads_to_all_bins.bam"
        map_reads_to_index(
            r1_raw, r2_raw, seed_index, competitive_bam, threads, seed_score_min,
            max_insert, competitive_dir / "bowtie2.log", report_multiple=20,
        )
        assignments, passing_templates, ambiguous_templates = competitively_assign_templates(
            competitive_bam, ref_to_bin, min_read_aligned, min_read_identity, threads,
        )
        all_seed_assigned = set().union(*assignments.values()) if assignments else set()
        metrics_path = competitive_dir / "competitive_mapping.tsv"
        with open(metrics_path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["rank", "bin", "uniquely_assigned_templates"])
            for bin_fasta in all_bin_fastas:
                w.writerow([r, bin_fasta.stem, len(assignments.get(bin_fasta.stem, set()))])
            w.writerow([r, "__PASSING_ALIGNMENTS__", passing_templates])
            w.writerow([r, "__AMBIGUOUS_TIES_EXCLUDED__", ambiguous_templates])
        log.info(
            "Rank '%s': competitive seed mapping assigned %d/%d passing templates; "
            "%d equal-score ties were excluded.",
            r, sum(len(v) for v in assignments.values()), passing_templates,
            ambiguous_templates,
        )

        summary_rows = []
        for bin_fasta in bin_fastas:
            bin_name = bin_fasta.stem
            before = basic_assembly_stats(bin_fasta)
            if before["num_contigs"] < min_bin_contigs_to_reassemble:
                log.info("Bin '%s' (%s): only %d contig(s); skipping reassembly.",
                          bin_name, r, before["num_contigs"])
                continue
            log.info("Reassembling bin '%s' (rank %s, %d contigs, %d bp)...",
                      bin_name, r, before["num_contigs"], before["total_length_bp"])
            try:
                reassembled = reassemble_one_bin(
                    bin_fasta, bin_name, r, r1_raw, r2_raw, competitive_bam,
                    assignments.get(bin_name, set()),
                    all_seed_assigned - assignments.get(bin_name, set()),
                    outdir, threads, trimmomatic_folder,
                    qc_quality, qc_minlen, flash_max_overlap, max_rounds, word_size,
                    min_word_hits, bait_min_entropy, max_round_growth,
                    max_accepted_fraction, min_new_templates, min_growth, bbtools_memory,
                    spades_mode, spades_memory_gb, kmers, circularize, unicycler_mode,
                )
            except RuntimeError as exc:
                log.warning("Bin '%s' (%s): reassembly failed (%s); keeping original bin untouched.",
                            bin_name, r, exc)
                reassembled = None
            if reassembled is None:
                continue
            successful[r][bin_name] = reassembled
            after = basic_assembly_stats(reassembled)
            summary_rows.append({
                "bin": bin_name, "rank": r,
                "contigs_before": before["num_contigs"], "contigs_after": after["num_contigs"],
                "total_length_before_bp": before["total_length_bp"],
                "total_length_after_bp": after["total_length_bp"],
                "N50_before": before["N50"], "N50_after": after["N50"],
                "largest_contig_before_bp": before["largest_contig_bp"],
                "largest_contig_after_bp": after["largest_contig_bp"],
            })
            log.info("Bin '%s' (%s): %d -> %d contigs, N50 %d -> %d bp.",
                      bin_name, r, before["num_contigs"], after["num_contigs"],
                      before["N50"], after["N50"])

        if not summary_rows:
            continue
        summary_path = outdir / "bins" / r / "reassembly_summary.tsv"
        fields = ["bin", "rank", "contigs_before", "contigs_after",
                  "total_length_before_bp", "total_length_after_bp",
                  "N50_before", "N50_after", "largest_contig_before_bp", "largest_contig_after_bp"]
        with open(summary_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(summary_rows)
        log.info("Rank '%s': wrote reassembly comparison summary to %s", r, summary_path)
    return {rank: dict(paths) for rank, paths in successful.items()}


def choose_final_source_rank(ranks, reassemble_ranks, requested: str = "auto") -> str:
    """Selects the single preliminary partition used to construct the consolidated final
    assembly. A single source rank is essential: combining nested family/genus/species
    reassemblies would represent the same original contigs and reads multiple times.

    Auto mode prefers genus as a practical genome-oriented compromise, then species,
    family, phylum, kingdom, superkingdom, and domain. The chosen rank must be both a
    requested binning rank and one of the ranks actually sent through reassembly.
    """
    available = [r for r in ranks if r in reassemble_ranks]
    if not available:
        raise ValueError("No shared rank exists between --ranks and --reassemble-ranks.")
    if requested and requested != "auto":
        if requested not in available:
            raise ValueError(
                f"--final-source-rank '{requested}' must occur in both --ranks and "
                f"--reassemble-ranks (available: {','.join(available)})."
            )
        return requested
    for rank in ("genus", "species", "family", "phylum", "kingdom", "superkingdom", "domain"):
        if rank in available:
            return rank
    return available[-1]


def build_consolidated_final_assembly(outdir: Path, source_rank: str,
                                      successful_reassemblies: dict) -> Path:
    """Builds one nonredundant-by-partition final contig set from a preliminary rank.

    For each preliminary bin at `source_rank`, a successful targeted reassembly replaces
    that bin's original contigs. Skipped or failed bins contribute their original contigs.
    The source rank partitions each retained initial contig exactly once, avoiding the
    duplication that would result from pooling nested reassemblies across several ranks.

    New globally unique FASTA identifiers are assigned and a provenance table records the
    source bin, source type, source FASTA, and original/reassembled record identifier.
    """
    rank_dir = outdir / "bins" / source_rank
    bin_fastas = sorted(rank_dir.glob("*.fasta"))
    if not bin_fastas:
        raise RuntimeError(
            f"No preliminary bin FASTAs found at source rank '{source_rank}' in {rank_dir}."
        )

    final_assembly_dir = outdir / "final" / "assembly"
    final_assembly_dir.mkdir(parents=True, exist_ok=True)
    final_fasta = final_assembly_dir / "consolidated_contigs.fasta"
    provenance_tsv = final_assembly_dir / "contig_provenance.tsv"

    success_at_rank = successful_reassemblies.get(source_rank, {})
    final_records = {}
    provenance_rows = []
    n_reassembled_bins = 0
    n_fallback_bins = 0

    for bin_fasta in bin_fastas:
        bin_name = bin_fasta.stem
        reassembled = success_at_rank.get(bin_name)
        if reassembled and Path(reassembled).exists() and Path(reassembled).stat().st_size > 0:
            source_fasta = Path(reassembled)
            source_type = "reassembled"
            n_reassembled_bins += 1
        else:
            source_fasta = bin_fasta
            source_type = "original"
            n_fallback_bins += 1

        source_records = read_fasta(source_fasta)
        safe_bin = sanitize(bin_name)
        for record_n, (source_id, seq) in enumerate(source_records.items(), 1):
            final_id = f"MH_{source_type}_{safe_bin}_{record_n:07d}"
            # A repeated/stale filename should never collide, but guard explicitly so a
            # final FASTA record can never be silently overwritten.
            collision_n = 1
            candidate = final_id
            while candidate in final_records:
                collision_n += 1
                candidate = f"{final_id}_{collision_n}"
            final_id = candidate
            final_records[final_id] = seq
            provenance_rows.append([
                final_id, source_rank, bin_name, source_type, str(source_fasta), source_id,
            ])

    if not final_records:
        raise RuntimeError("Final-assembly consolidation produced no contigs.")

    write_fasta(final_fasta, final_records)
    with open(provenance_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "final_contig", "source_rank", "source_bin", "source_type",
            "source_fasta", "source_record",
        ])
        w.writerows(provenance_rows)

    log.info(
        "Consolidated final assembly from rank '%s': %d contigs; %d reassembled bin(s), "
        "%d original-fallback bin(s).",
        source_rank, len(final_records), n_reassembled_bins, n_fallback_bins,
    )
    return final_fasta


# --------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Assemble, classify, and bin metagenomic contigs. Whenever reads are supplied, "
                    "the default is preliminary microbial retention -> competitive ITSME-style "
                    "seed-and-extend "
                    "reassembly -> consolidated final Prodigal/DIAMOND classification -> "
                    "multirank bins with QUAST/CheckM summaries.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-i", "--input", type=Path, default=None,
                    help="Optional input contigs FASTA. Use alone for classification/binning, or "
                         "together with -1/-2 to use these contigs as the initial assembly and "
                         "enable read recruitment without running MEGAHIT.")
    p.add_argument("-1", "--r1", type=Path, default=None, help="Raw/forward paired-end FASTQ (R1)")
    p.add_argument("-2", "--r2", type=Path, default=None, help="Raw/reverse paired-end FASTQ (R2)")
    p.add_argument("-d", "--diamond-db", required=True, type=Path,
                    help="Path to taxonomy-enabled nr-tax.dmnd (DIAMOND >=2.1.17)")
    p.add_argument("-o", "--outdir", required=True, type=Path, help="Output directory")
    p.add_argument("-t", "--threads", type=int, default=8)
    p.add_argument("--skip-qc", action="store_true",
                    help="Feed -1/-2 straight to MEGAHIT, skipping the Trimmomatic/FLASH QC step "
                         "(--trim-polyg, if given, still runs).")

    pg = p.add_argument_group("Poly-G trimming (fastp; optional, runs before everything else)")
    pg.add_argument("--trim-polyg", action="store_true",
                     help="Trim poly-G tails with fastp before adapter clipping. Poly-G runs are a "
                          "known artifact of two-channel Illumina chemistry (NextSeq/NovaSeq) where "
                          "dark/no-signal cycles get miscalled as 'G' -- Trimmomatic's adapter/quality "
                          "trimming doesn't reliably catch these. Turn this on if your reads are known "
                          "to have long poly-G tails.")
    pg.add_argument("--poly-g-min-len", type=int, default=10,
                     help="Minimum length of a 3' G-run to trim (fastp --poly_g_min_len, default 10)")

    qc = p.add_argument_group("QC options (Trimmomatic + FLASH; used with -1/-2 unless --skip-qc)")
    qc.add_argument("--trimmomatic-folder", type=Path, default=None,
                     help="Path to the Trimmomatic install folder containing "
                          "adapters/TruSeq3-PE-2.fa (required unless --skip-qc)")
    qc.add_argument("--flash-max-overlap", type=int, default=150)
    qc.add_argument("--qc-quality", type=int, default=20,
                     help="Trimmomatic SLIDINGWINDOW:4:<qc_quality> for all quality-trim passes")
    qc.add_argument("--qc-minlen", type=int, default=50, help="Trimmomatic MINLEN for quality-trim passes")
    qc.add_argument("--keep-qc-tmp", action="store_true", help="Keep intermediate QC files")

    mh = p.add_argument_group("MEGAHIT options (used with -1/-2)")
    mh.add_argument("--megahit-min-contig-len", type=int, default=None,
                     help="MEGAHIT --min-contig-len (default: MEGAHIT's own default, 200bp)")
    mh.add_argument("--megahit-extra", default=None,
                     help="Extra raw arguments passed through to MEGAHIT verbatim, "
                          "e.g. --megahit-extra '--k-list 21,41,61'")

    p.add_argument("--prodigal-mode", choices=["single", "meta"], default="meta",
                    help="Prodigal procedure (default: meta, appropriate for mixed organisms).")

    tr = p.add_argument_group("Contig triage and bin refinement")
    tr.add_argument("--skip-contig-triage", action="store_true",
                    help="Disable the default pre-DIAMOND gene-density screen. By default, "
                         f"contigs >= {TRIAGE_MIN_LENGTH_BP} bp with < "
                         f"{TRIAGE_MAX_EUKARYOTIC_CODING_DENSITY * 100:.0f}%% coding density are "
                         "quarantined as eukaryotic-like; short/uncertain contigs are retained.")
    tr.add_argument("--skip-bin-refinement", action="store_true",
                    help="Disable conservative GC+coverage coherence checks. With reads, joint "
                         "GC-and-depth outliers are excluded from preliminary seed references "
                         "and demoted to Unclassified at the selected rank after reassembly. "
                         "Single-signal outliers are only reported.")
    tr.add_argument("--refinement-rank", default="auto",
                    help="Rank used for within-bin GC+coverage coherence (default auto: species, "
                         "then genus/family/...); must also occur in --ranks.")

    p.add_argument("-e", "--evalue", type=float, default=1e-5, help="DIAMOND e-value cutoff")
    p.add_argument("--max-target-seqs", type=int, default=25, help="DIAMOND -k (hits kept per ORF)")

    p.add_argument("--ranks", default=",".join(BIN_RANKS_DEFAULT),
                    help="Comma-separated ranks to classify/bin at (default: genus,species). "
                         "domain/superkingdom/phylum/family are also supported.")
    p.add_argument("--bitscore-range", type=float, default=0.9,
                    help="Keep hits within this fraction of an ORF's best bitscore (default 0.9)")
    p.add_argument("--max-hits-per-orf", type=int, default=5)
    p.add_argument("--min-support", type=float, default=0.5,
                    help="Minimum bitscore-weighted vote share required to assign a taxon at a "
                         "rank (0.5 = majority). Use 0 for a pure plurality/'most votes wins' call.")

    p.add_argument("--min-bin-contigs", type=int, default=1,
                    help="Drop bins with fewer than this many contigs (default 1, i.e. keep all). "
                         "Contigs from dropped bins fall back to 'Unclassified'.")
    p.add_argument("--min-bin-length", type=int, default=0,
                    help="Drop bins whose total length (bp) is below this (default 0, i.e. keep all).")
    p.add_argument("--exclude-unclassified-bins", action="store_true",
                    help="Do not write an Unclassified.fasta bin (default: include it)")

    p.add_argument("--exclude-kingdoms", default=",".join(DEFAULT_EXCLUDED_KINGDOMS),
                    help="Comma-separated Eukaryota kingdoms to drop from binning entirely "
                         "(default: Metazoa,Viridiplantae -- removes host-animal/plant "
                         "contamination while keeping Bacteria, Archaea, Fungi, and protists). "
                         "Dropped contigs are written to "
                         "<outdir>/classification/excluded_animal_plant_contamination.fasta "
                         "and flagged in contig_classification.tsv, not silently discarded. "
                         "Pass '' to disable this filter and keep everything.")

    p.add_argument("--reuse-prodigal", action="store_true",
                    help="Skip Prodigal if <outdir>/prodigal/proteins.faa already exists (reuse it).")
    p.add_argument("--reuse-diamond", action="store_true",
                    help="Skip DIAMOND if <outdir>/diamond/hits.tsv already exists (reuse it). "
                         "Combine with --reuse-prodigal and a new --ranks/--min-support to "
                         "re-classify/re-bin at a different granularity without re-running "
                         "Prodigal+DIAMOND.")

    p.add_argument("--skip-quast", action="store_true")
    p.add_argument("--skip-checkm", action="store_true")

    rb = p.add_argument_group(
        "Targeted bin reassembly (default whenever reads are supplied; competitive mapping, "
        "frontier extension, focused SPAdes/metaSPAdes, Unicycler, and Pilon)"
    )
    rb_toggle = rb.add_mutually_exclusive_group()
    rb_toggle.add_argument("--reassemble-bins", dest="reassemble_bins", action="store_true",
                     help="Explicitly enable the default workflow used whenever -1/-2 are present: "
                          "competitively map reads against all seed bins, exclude tied assignments, "
                          "extend winning pools, run focused SPAdes/metaSPAdes, then attempt "
                          "Unicycler circularization and Pilon polishing.")
    rb_toggle.add_argument("--skip-reassembly", dest="reassemble_bins", action="store_false",
                     help="Disable the default seed-and-extension/final-reclassification workflow "
                          "and stop after preliminary classification, binning, and QUAST/CheckM. "
                          "This is the implicit behavior only when contigs are supplied without reads.")
    p.set_defaults(reassemble_bins=None)
    rb.add_argument("--reassemble-ranks", default=None,
                     help="Comma-separated subset of --ranks to reassemble. By default, only the "
                          "automatically selected final source rank is reassembled (genus is "
                          "preferred). Additional ranks are diagnostic and are not pooled into "
                          "the consolidated final assembly. Each rank must also appear in --ranks.")
    rb.add_argument("--final-source-rank", default="auto",
                     help="Single preliminary rank whose mutually exclusive bins are used to "
                          "construct the consolidated final assembly. Successful reassemblies "
                          "replace their original bin contigs; skipped/failed bins fall back to "
                          "the originals. Default 'auto' prefers genus, then species, family, "
                          "phylum, kingdom, superkingdom, or domain. The rank must also occur in "
                          "--reassemble-ranks when that option is given.")
    rb.add_argument("--reassemble-min-bin-contigs", type=int, default=2,
                     help="Skip reassembly for bins with fewer contigs than this -- nothing to "
                          "gain from reassembling an already-single-contig bin (default 2).")
    rb.add_argument("--reassemble-include-unclassified", action="store_true",
                     help="Also attempt reassembly of the catch-all 'Unclassified' bin (default: "
                          "skipped, since it's a mixed leftover pool, not one coherent genome).")
    rb.add_argument("--anchor-db", action="append", default=[], metavar="BIN=FASTA",
                     help="Optional external seed assigned to exactly one competing bin; repeatable. "
                          "BIN is a bin FASTA stem, e.g. Escherichia_coli=reference.fasta.")
    rb.add_argument("--reassemble-seed-score-min", default="G,20,8",
                     help="Bowtie2 --score-min for the seed-recruitment mapping (default G,20,8).")
    rb.add_argument("--reassemble-min-read-aligned", type=float, default=0.75,
                     help="Minimum aligned fraction of a read for it to be recruited (default 0.75).")
    rb.add_argument("--reassemble-min-read-identity", type=float, default=0.85,
                     help="Minimum approximate identity for a recruited read (default 0.85).")
    rb.add_argument("--reassemble-max-insert", type=int, default=1000,
                     help="Maximum bowtie2 fragment length during recruitment (default 1000).")
    rb.add_argument("--reassemble-max-rounds", type=int, default=5,
                     help="Exact-kmer frontier-extension rounds after competitive seed "
                          "recruitment (default 5; pass 0 to disable extension).")
    rb.add_argument("--reassemble-word-size", type=int, default=31,
                     help="Exact k-mer/word length for frontier extension (default 31; 15-31).")
    rb.add_argument("--reassemble-min-word-hits", type=int, default=3,
                     help="Minimum exact-word hits required to recruit a read during extension "
                          "(default 3).")
    rb.add_argument("--reassemble-bait-min-entropy", type=float, default=0.45,
                     help="Excludes low-complexity bait sequence from frontier extension "
                          "(default 0.45).")
    rb.add_argument("--reassemble-max-round-growth", type=float, default=0.25,
                     help="Reject an extension round if new/accepted exceeds this fraction -- "
                          "guards against snowballing into an unrelated, similar-coverage genome "
                          "(default 0.25).")
    rb.add_argument("--reassemble-max-accepted-fraction", type=float, default=0.05,
                     help="Reject an extension round if accepted/total-raw-reads would exceed "
                          "this fraction (default 0.05).")
    rb.add_argument("--reassemble-min-new-templates", type=int, default=10,
                     help="Stop extension once a round recruits fewer new templates than this "
                          "(default 10).")
    rb.add_argument("--reassemble-min-growth", type=float, default=0.0001,
                     help="Stop extension once round-over-round growth falls below this fraction "
                          "(default 0.0001).")
    rb.add_argument("--reassemble-bbtools-memory", default="16g",
                     help="Java heap for BBDuk during frontier extension, e.g. 16g (default 16g).")
    rb.add_argument("--reassemble-mode", choices=["standard", "meta"], default="meta",
                     help="spades.py mode for the per-bin reassembly: 'meta' (metaSPAdes, "
                          "default -- more forgiving of residual strain heterogeneity/uneven "
                          "coverage in the recruited pool) or 'standard' (plain SPAdes, for a "
                          "very clean single-strain recruitment).")
    rb.add_argument("--reassemble-memory-gb", type=int, default=32,
                     help="SPAdes memory limit in GB for each bin's reassembly (default 32).")
    rb.add_argument("--reassemble-kmers", default="auto",
                     help="SPAdes k-mer list or 'auto' (default auto).")
    rb.add_argument("--skip-circularization", dest="circularize", action="store_false",
                     help="Keep focused SPAdes/metaSPAdes output; skip the default Unicycler "
                          "circularization attempt and Pilon polishing.")
    rb.add_argument("--unicycler-mode", choices=["conservative", "normal", "bold"],
                     default="normal", help="Unicycler bridging mode (default normal).")
    p.set_defaults(circularize=True)

    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def parse_anchor_specs(specs) -> dict:
    """Parse repeatable ``BIN=FASTA`` anchors for competitive seed mapping."""
    anchors = defaultdict(list)
    for spec in specs or []:
        if "=" not in spec:
            raise ValueError(f"--anchor-db must be BIN=FASTA, got: {spec}")
        raw_bin, raw_path = spec.split("=", 1)
        bin_name = sanitize(raw_bin.strip())
        anchor_path = Path(raw_path.strip())
        if not raw_bin.strip() or not raw_path.strip():
            raise ValueError(f"--anchor-db must be BIN=FASTA, got: {spec}")
        if not anchor_path.is_file():
            raise ValueError(f"Anchor FASTA does not exist: {anchor_path}")
        anchors[bin_name].append(anchor_path)
    return dict(anchors)


def main(argv=None):
    args = parse_args(argv)
    setup_logging(args.verbose)

    using_reads = args.r1 is not None or args.r2 is not None
    using_contigs = args.input is not None
    if using_reads and (args.r1 is None or args.r2 is None):
        log.error("Both -1 and -2 are required together.")
        sys.exit(1)
    if not using_reads and not using_contigs:
        log.error("Provide contigs (-i), paired reads (-1/-2), or both.")
        sys.exit(1)

    # Seed-and-extension is the default whenever reads are available. Contigs-only
    # mode implicitly disables it unless the user explicitly requested it, which is an
    # error because no reads exist to recruit.
    reassembly_was_explicit = args.reassemble_bins is True
    if args.reassemble_bins is None:
        args.reassemble_bins = using_reads
    if reassembly_was_explicit and not using_reads:
        log.error("--reassemble-bins requires paired reads (-1/-2); contigs alone cannot extend.")
        sys.exit(1)

    reads_only = using_reads and not using_contigs
    if reads_only and not args.skip_qc and args.trimmomatic_folder is None:
        log.error("--trimmomatic-folder is required for QC (or pass --skip-qc to bypass QC).")
        sys.exit(1)
    if args.reassemble_bins and args.trimmomatic_folder is None:
        log.error("Seed-and-extension is enabled by default whenever reads are supplied and requires "
                  "--trimmomatic-folder because every recruited batch is QC processed. "
                  "Provide the folder or pass --skip-reassembly.")
        sys.exit(1)
    if args.reassemble_bins and not (15 <= args.reassemble_word_size <= 31):
        log.error("--reassemble-word-size must be between 15 and 31.")
        sys.exit(1)
    try:
        anchors_by_bin = parse_anchor_specs(args.anchor_db)
    except ValueError as exc:
        log.error("%s", exc)
        sys.exit(1)

    ranks = [r.strip() for r in args.ranks.split(",") if r.strip()]
    if not ranks:
        log.error("--ranks must contain at least one supported rank.")
        sys.exit(1)
    for r in ranks:
        if r not in WANTED_RANKS:
            log.error("Unsupported rank '%s'. Supported: %s", r, WANTED_RANKS)
            sys.exit(1)
    try:
        refinement_rank = choose_refinement_rank(ranks, args.refinement_rank)
    except ValueError as exc:
        log.error("%s", exc)
        sys.exit(1)
    reassemble_ranks = []
    final_source_rank = None
    if args.reassemble_bins:
        if args.reassemble_ranks:
            reassemble_ranks = [r.strip() for r in args.reassemble_ranks.split(",") if r.strip()]
        else:
            try:
                auto_source = choose_final_source_rank(ranks, ranks, args.final_source_rank)
            except ValueError as exc:
                log.error("%s", exc)
                sys.exit(1)
            reassemble_ranks = [auto_source]
        for r in reassemble_ranks:
            if r not in ranks:
                log.error("--reassemble-ranks '%s' was not classified/binned (--ranks was '%s').",
                          r, args.ranks)
                sys.exit(1)
        try:
            final_source_rank = choose_final_source_rank(
                ranks, reassemble_ranks, args.final_source_rank,
            )
        except ValueError as exc:
            log.error("%s", exc)
            sys.exit(1)
        log.info(
            "Seed-and-extension enabled at rank(s) %s; consolidated final assembly source rank: %s.",
            reassemble_ranks, final_source_rank,
        )

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    faa = outdir / "prodigal" / "proteins.faa"
    gff = outdir / "prodigal" / "genes.gff"
    reuse_prodigal = (
        args.reuse_prodigal and faa.exists()
        and (args.skip_contig_triage or gff.exists())
    )
    hits_tsv = outdir / "diamond" / "hits.tsv"
    reuse_diamond = args.reuse_diamond and hits_tsv.exists()

    if using_reads:
        if args.trim_polyg:
            which_or_die("fastp")
        if reads_only and not args.skip_qc:
            which_or_die("trimmomatic")
            which_or_die("flash")
        if reads_only:
            which_or_die("megahit")
    if args.reassemble_bins or not reuse_prodigal:
        which_or_die("prodigal")
    if args.reassemble_bins or not reuse_diamond:
        which_or_die("diamond")
        require_diamond_taxonomy_fields()
    if not args.skip_quast and shutil.which("quast.py") is None:
        log.warning("quast.py not on PATH; will use built-in assembly stats instead.")
    if not args.skip_checkm and shutil.which("checkm") is None:
        log.warning("checkm not on PATH; completeness/contamination will be NA.")
    needs_coverage_mapping = using_reads and not args.skip_bin_refinement
    if args.reassemble_bins or needs_coverage_mapping:
        which_or_die("bowtie2")
        which_or_die("bowtie2-build")
        which_or_die("samtools")
    if args.reassemble_bins:
        which_or_die("trimmomatic")
        which_or_die("flash")
        which_or_die("bbduk.sh")
        which_or_die("spades.py")
        if args.circularize:
            which_or_die("unicycler")
            which_or_die("pilon")

    n_qc_steps = (
        (1 if using_reads and args.trim_polyg else 0)
        + (1 if reads_only and not args.skip_qc else 0)
        + (1 if reads_only else 0)
    )
    refinement_steps = (2 if args.reassemble_bins else 1) if needs_coverage_mapping else 0
    steps = StepCounter(n_qc_steps + (11 if args.reassemble_bins else 5) + refinement_steps)

    # 0a/0b/0c. Poly-G may prepare reads in either read mode; QC+MEGAHIT are reads-only.
    assembly_fasta = args.input
    if using_reads:
        r1_in, r2_in = args.r1, args.r2
        if args.trim_polyg:
            steps.next("Trimming poly-G tails (fastp)...")
            r1_in, r2_in = trim_poly_g(
                r1_in, r2_in, outdir / "polyg", args.threads,
                poly_g_min_len=args.poly_g_min_len,
            )

        if reads_only:
            if not args.skip_qc:
                steps.next("Running QC (Trimmomatic adapter/quality trim + FLASH merge)...")
                r1_final, r2_final, u_final = run_qc(
                    r1_in, r2_in, outdir / "qc", args.threads,
                    trimmomatic_cmd="trimmomatic", trimmomatic_folder=args.trimmomatic_folder,
                    flash_cmd="flash", flash_max_overlap=args.flash_max_overlap, pigz_cmd="pigz",
                    qc_quality=args.qc_quality, qc_minlen=args.qc_minlen,
                    keep_tmp=args.keep_qc_tmp,
                )
            else:
                r1_final, r2_final, u_final = r1_in, r2_in, None

            steps.next("Running MEGAHIT assembly...")
            assembly_fasta = run_megahit(
                r1_final, r2_final, u_final, outdir / "megahit", args.threads,
                min_contig_len=args.megahit_min_contig_len, extra_args=args.megahit_extra,
            )
            log.info("MEGAHIT assembly: %s", assembly_fasta)
        else:
            log.info(
                "Using supplied contigs as the initial assembly; reads are reserved for "
                "competitive seed-and-extension."
            )

    contig_seqs = read_fasta(assembly_fasta)

    # 1. Preliminary Prodigal and default pre-DIAMOND gene-density triage.
    if reuse_prodigal:
        steps.next(f"Reusing preliminary Prodigal output: {faa}")
    else:
        steps.next("Running preliminary Prodigal...")
        faa, gff = run_prodigal(assembly_fasta, outdir / "prodigal", mode=args.prodigal_mode)
    orf_to_contig = parse_orf_to_contig(faa)
    contig_metrics = compute_contig_metrics(contig_seqs, gff)
    triage_calls, triage_excluded_ids = triage_contigs(
        contig_metrics, disabled=args.skip_contig_triage,
    )
    classification_dir = outdir / "classification"
    classification_dir.mkdir(parents=True, exist_ok=True)
    triage_excluded_path = classification_dir / "excluded_eukaryotic_like_gene_density.fasta"
    triage_excluded_path.unlink(missing_ok=True)
    if triage_excluded_ids:
        write_fasta(
            triage_excluded_path,
            {contig_id: contig_seqs[contig_id] for contig_id in triage_excluded_ids},
        )
        log.info(
            "Gene-density triage quarantined %d long, low-coding-density contig(s) before DIAMOND.",
            len(triage_excluded_ids),
        )
    diamond_query_faa = write_candidate_proteins(
        faa, orf_to_contig, triage_excluded_ids,
        outdir / "prodigal" / "proteins.prokaryotic_candidates.faa",
    )
    log.info("Predicted %d ORFs.", len(orf_to_contig))

    # 2. Preliminary DIAMOND
    if reuse_diamond:
        steps.next(f"Reusing preliminary DIAMOND output: {hits_tsv}")
    else:
        steps.next(f"Running preliminary DIAMOND blastp vs {args.diamond_db}...")
        hits_tsv = run_diamond(
            diamond_query_faa, args.diamond_db, outdir / "diamond", args.threads, args.evalue,
            args.max_target_seqs,
        )
    hits_by_orf = parse_diamond_hits(hits_tsv)
    log.info("Got hits for %d/%d ORFs.", len(hits_by_orf), len(orf_to_contig))

    # 3. Preliminary classification and microbial-contig retention
    # "domain" and "kingdom" are always classified internally (even if not in --ranks) so
    # the animal/plant-contamination filter below can always run; they're only written as
    # bins if you actually asked for them in --ranks.
    internal_ranks = ranks + [r for r in ("domain", "kingdom") if r not in ranks]
    steps.next(f"Preliminary contig classification/retention at ranks: {internal_ranks}...")
    classifications = classify_all_contigs(
        list(contig_seqs.keys()), orf_to_contig, hits_by_orf, internal_ranks,
        args.bitscore_range, args.max_hits_per_orf, args.min_support,
    )

    exclude_kingdoms = [k.strip() for k in args.exclude_kingdoms.split(",") if k.strip()]
    triage_retained = {
        contig_id: result for contig_id, result in classifications.items()
        if contig_id not in triage_excluded_ids
    }
    classifications_for_binning, excluded_ids = split_excluded_eukaryotes(
        triage_retained, exclude_kingdoms,
    )
    preliminary_excluded_path = classification_dir / "excluded_animal_plant_contamination.fasta"
    preliminary_excluded_path.unlink(missing_ok=True)
    if excluded_ids:
        log.info(
            "Excluding %d contig(s) classified as %s (host-animal/plant contamination).",
            len(excluded_ids), exclude_kingdoms,
        )
        excluded_records = {cid: contig_seqs[cid] for cid in excluded_ids if cid in contig_seqs}
        write_fasta(
            preliminary_excluded_path,
            excluded_records,
        )

    preliminary_refinement_decisions = {}
    seed_refinement_outliers = set()
    if needs_coverage_mapping:
        steps.next(
            f"Estimating preliminary coverage and checking {refinement_rank}-bin coherence..."
        )
        coverage = estimate_contig_coverage(
            assembly_fasta, r1_in, r2_in, classification_dir / "coverage",
            args.threads, args.reassemble_max_insert,
        )
        add_coverage_to_metrics(contig_metrics, coverage)
        preliminary_refinement_decisions, seed_refinement_outliers = refine_taxonomic_bins(
            classifications_for_binning, contig_metrics, refinement_rank,
            enabled=not args.skip_bin_refinement,
        )
        write_refinement_table(
            preliminary_refinement_decisions, contig_metrics,
            classification_dir / "bin_refinement.tsv",
        )
        preliminary_outlier_path = classification_dir / "seed_refinement_outliers.fasta"
        preliminary_outlier_path.unlink(missing_ok=True)
        if seed_refinement_outliers:
            write_fasta(
                preliminary_outlier_path,
                {contig_id: contig_seqs[contig_id] for contig_id in seed_refinement_outliers},
            )
            if args.reassemble_bins:
                log.info(
                    "Excluding %d joint GC/coverage outlier(s) from seed references; "
                    "the contigs remain in preliminary bins and the final assembly.",
                    len(seed_refinement_outliers),
                )
            else:
                log.info(
                    "Demoting %d joint GC/coverage outlier(s) at %s and lower ranks.",
                    len(seed_refinement_outliers), refinement_rank,
                )
                classifications = demote_refinement_outliers(
                    classifications, seed_refinement_outliers, refinement_rank,
                )
                classifications_for_binning = {
                    contig_id: classifications[contig_id]
                    for contig_id in classifications_for_binning
                }

    write_classification_table(
        classifications, internal_ranks, classification_dir / "contig_classification.tsv",
        excluded_ids=excluded_ids, contig_metrics=contig_metrics,
        triage_calls=triage_calls, triage_excluded_ids=triage_excluded_ids,
        refinement_decisions=preliminary_refinement_decisions,
    )

    # 4. Preliminary binning (excluded contigs never become seed bins)
    steps.next("Writing preliminary seed-bin FASTA files...")
    membership = bin_contigs(
        contig_seqs, classifications_for_binning, ranks, outdir / "bins",
        include_unclassified=not args.exclude_unclassified_bins,
    )
    if args.min_bin_contigs > 1 or args.min_bin_length > 0:
        log.info(
            "Merging bins smaller than %d contigs / %d bp into Unclassified...",
            args.min_bin_contigs, args.min_bin_length,
        )
        filter_small_bins(
            contig_seqs, membership, ranks, outdir / "bins",
            args.min_bin_contigs, args.min_bin_length,
            include_unclassified=not args.exclude_unclassified_bins,
        )

    # 5. Default targeted seed-and-extension reassembly, then a complete final pass.
    if args.reassemble_bins:
        assembler_label = "metaSPAdes" if args.reassemble_mode == "meta" else "SPAdes"
        steps.next(f"Targeted bin reassembly at ranks {reassemble_ranks} "
                   f"(competitive seed-and-extend + focused {assembler_label}"
                   f"{' + Unicycler/Pilon' if args.circularize else ''})...")
        successful_reassemblies = run_bin_reassembly(
            outdir, reassemble_ranks, r1_in, r2_in, args.threads, anchors_by_bin,
            args.trimmomatic_folder, args.qc_quality, args.qc_minlen, args.flash_max_overlap,
            args.reassemble_seed_score_min, args.reassemble_min_read_aligned,
            args.reassemble_min_read_identity, args.reassemble_max_insert,
            args.reassemble_max_rounds, args.reassemble_word_size, args.reassemble_min_word_hits,
            args.reassemble_bait_min_entropy, args.reassemble_max_round_growth,
            args.reassemble_max_accepted_fraction, args.reassemble_min_new_templates,
            args.reassemble_min_growth, args.reassemble_bbtools_memory, args.reassemble_mode,
            args.reassemble_memory_gb, args.reassemble_kmers,
            args.reassemble_min_bin_contigs, args.reassemble_include_unclassified,
            args.circularize, args.unicycler_mode,
            seed_excluded_ids=seed_refinement_outliers,
        )

        steps.next(f"Consolidating final assembly from preliminary rank '{final_source_rank}'...")
        final_assembly_fasta = build_consolidated_final_assembly(
            outdir, final_source_rank, successful_reassemblies,
        )

        # 6. Final Prodigal, gene-density triage, and DIAMOND pass.
        steps.next("Running final Prodigal on consolidated contigs...")
        final_faa, final_gff = run_prodigal(
            final_assembly_fasta, outdir / "final" / "prodigal", mode=args.prodigal_mode,
        )
        final_contig_seqs = read_fasta(final_assembly_fasta)
        final_orf_to_contig = parse_orf_to_contig(final_faa)
        final_contig_metrics = compute_contig_metrics(final_contig_seqs, final_gff)
        final_triage_calls, final_triage_excluded_ids = triage_contigs(
            final_contig_metrics, disabled=args.skip_contig_triage,
        )
        final_classification_dir = outdir / "final" / "classification"
        final_classification_dir.mkdir(parents=True, exist_ok=True)
        final_triage_excluded_path = (
            final_classification_dir / "excluded_eukaryotic_like_gene_density.fasta"
        )
        final_triage_excluded_path.unlink(missing_ok=True)
        if final_triage_excluded_ids:
            write_fasta(
                final_triage_excluded_path,
                {
                    contig_id: final_contig_seqs[contig_id]
                    for contig_id in final_triage_excluded_ids
                },
            )
            log.info(
                "Final gene-density triage quarantined %d contig(s) before DIAMOND.",
                len(final_triage_excluded_ids),
            )
        final_diamond_query_faa = write_candidate_proteins(
            final_faa, final_orf_to_contig, final_triage_excluded_ids,
            outdir / "final" / "prodigal" / "proteins.prokaryotic_candidates.faa",
        )
        log.info("Final assembly: predicted %d ORFs.", len(final_orf_to_contig))

        steps.next(f"Running final DIAMOND blastp vs {args.diamond_db}...")
        final_hits_tsv = run_diamond(
            final_diamond_query_faa, args.diamond_db, outdir / "final" / "diamond", args.threads,
            args.evalue, args.max_target_seqs,
        )
        final_hits_by_orf = parse_diamond_hits(final_hits_tsv)
        log.info(
            "Final assembly: got hits for %d/%d ORFs.",
            len(final_hits_by_orf), len(final_orf_to_contig),
        )

        # 7. Reclassify and reapply the animal/plant filter because frontier extension
        # can introduce contigs whose taxonomy differs from the preliminary seed bin.
        steps.next(f"Final contig classification/retention at ranks: {internal_ranks}...")
        final_classifications = classify_all_contigs(
            list(final_contig_seqs.keys()), final_orf_to_contig, final_hits_by_orf,
            internal_ranks, args.bitscore_range, args.max_hits_per_orf,
            args.min_support,
        )
        final_triage_retained = {
            contig_id: result for contig_id, result in final_classifications.items()
            if contig_id not in final_triage_excluded_ids
        }
        final_classifications_for_binning, final_excluded_ids = split_excluded_eukaryotes(
            final_triage_retained, exclude_kingdoms,
        )
        final_excluded_path = final_classification_dir / "excluded_animal_plant_contamination.fasta"
        final_excluded_path.unlink(missing_ok=True)
        if final_excluded_ids:
            log.info(
                "Final pass: excluding %d contig(s) classified as %s.",
                len(final_excluded_ids), exclude_kingdoms,
            )
            final_excluded_records = {
                cid: final_contig_seqs[cid]
                for cid in final_excluded_ids if cid in final_contig_seqs
            }
            write_fasta(final_excluded_path, final_excluded_records)

        final_refinement_decisions = {}
        final_refinement_outliers = set()
        if needs_coverage_mapping:
            steps.next(
                f"Remapping all reads and refining final {refinement_rank}-level bins..."
            )
            final_coverage = estimate_contig_coverage(
                final_assembly_fasta, r1_in, r2_in,
                final_classification_dir / "coverage", args.threads,
                args.reassemble_max_insert,
            )
            add_coverage_to_metrics(final_contig_metrics, final_coverage)
            final_refinement_decisions, final_refinement_outliers = refine_taxonomic_bins(
                final_classifications_for_binning, final_contig_metrics, refinement_rank,
                enabled=not args.skip_bin_refinement,
            )
            write_refinement_table(
                final_refinement_decisions, final_contig_metrics,
                final_classification_dir / "bin_refinement.tsv",
            )
            final_refinement_path = final_classification_dir / "refinement_outliers.fasta"
            final_refinement_path.unlink(missing_ok=True)
            if final_refinement_outliers:
                write_fasta(
                    final_refinement_path,
                    {
                        contig_id: final_contig_seqs[contig_id]
                        for contig_id in final_refinement_outliers
                    },
                )
                log.info(
                    "Final refinement demoted %d joint GC/coverage outlier(s) at %s "
                    "and more-specific ranks.",
                    len(final_refinement_outliers), refinement_rank,
                )
                final_classifications = demote_refinement_outliers(
                    final_classifications, final_refinement_outliers, refinement_rank,
                )
                final_classifications_for_binning = {
                    contig_id: final_classifications[contig_id]
                    for contig_id in final_classifications_for_binning
                }
        write_classification_table(
            final_classifications, internal_ranks,
            final_classification_dir / "contig_classification.tsv",
            excluded_ids=final_excluded_ids, contig_metrics=final_contig_metrics,
            triage_calls=final_triage_calls,
            triage_excluded_ids=final_triage_excluded_ids,
            refinement_decisions=final_refinement_decisions,
        )

        # 8. Final multirank bins and quality assessment.
        steps.next("Writing final taxonomic bin FASTA files...")
        final_bins_dir = outdir / "final" / "bins"
        final_membership = bin_contigs(
            final_contig_seqs, final_classifications_for_binning, ranks, final_bins_dir,
            include_unclassified=not args.exclude_unclassified_bins,
        )
        if args.min_bin_contigs > 1 or args.min_bin_length > 0:
            filter_small_bins(
                final_contig_seqs, final_membership, ranks, final_bins_dir,
                args.min_bin_contigs, args.min_bin_length,
                include_unclassified=not args.exclude_unclassified_bins,
            )

        steps.next("Assessing final bins with QUAST and CheckM...")
        for r in ranks:
            summarize_bin_set(
                r, final_bins_dir / r, args.threads, args.skip_quast, args.skip_checkm,
            )
        log.info("Final classified and assessed bins: %s", final_bins_dir)
    else:
        # With contigs alone, or with --skip-reassembly, the preliminary assembly is the
        # final assembly and receives the one requested quality-assessment pass.
        steps.next("Assessing bins with QUAST and CheckM (reassembly disabled)...")
        for r in ranks:
            summarize_bin_set(
                r, outdir / "bins" / r, args.threads, args.skip_quast, args.skip_checkm,
            )

    log.info("Done. Results in %s", outdir)


if __name__ == "__main__":
    main()
