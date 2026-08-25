#!/usr/bin/env python3
"""
metatax_binner.py

End-to-end pipeline for taxonomic classification and binning of metagenomic contigs.
Input is either a pair of raw FASTQs (QC'd and assembled first) or an already-assembled
contigs FASTA.

    R1.fastq, R2.fastq  (raw paired-end reads; optional -- skip straight to contigs.fasta)
        -> [optional] fastp poly-G trimming (--trim-polyg -- for NextSeq/NovaSeq
           two-channel-chemistry poly-G tail artifacts)
        -> QC: Trimmomatic (adapter clip + quality trim) + FLASH (merge overlapping pairs)
           [skippable with --skip-qc]
        -> MEGAHIT assembly
    contigs.fasta
        -> Prodigal (ORF/gene prediction, meta mode)
        -> DIAMOND blastp vs. a pre-built nr.dmnd database (any nr.dmnd works -- no
           --taxonmap/--taxonnodes needed; organism names are parsed straight out of
           the trailing "[Organism name]" in each hit's stitle, standard NCBI nr format)
        -> per-ORF top hits -> per-contig taxonomic classification
        -> bin FASTA files at genus / species level (or whatever --ranks you ask for)
        -> per-bin-set summary.tsv (QUAST contiguity + CheckM completeness/contamination)

CLASSIFICATION METHOD (what this script does and why)
-------------------------------------------------------
Classifying a *contig* from many *ORF* hits is the same problem CAT/BAT, MEGAN's LCA,
and Kraken-style consensus callers solve, and the approach used here follows that lineage:

  1. For each ORF, keep the diamond hits within `--bitscore-range` (default 90%) of that
     ORF's best bitscore (i.e. a "bit-score competitive set", not just the single top hit).
     This avoids over-trusting one alignment when several equally-good references exist.
  2. Each kept hit's organism name is parsed from its stitle (the "[Genus species]" at the
     end), resolved to an NCBI lineage by name (local taxdump or ete3), and contributes its
     full bitscore as a "vote" for every taxon in that lineage.
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
    diamond             https://github.com/bbuchfink/diamond   (any nr.dmnd built with
                         `diamond makedb --in nr.faa[.gz] --db nr` works -- taxonomy is NOT
                         read from the database; organism names come from stitle instead, so
                         headers must look like standard NCBI nr: "... [Genus species]")
    quast.py            https://github.com/ablab/quast          (optional; falls back to a
                         built-in N50/L50/GC calculator if missing or --skip-quast is given)
    checkm (lineage_wf) https://github.com/Ecogenomics/CheckM    (optional; skipped with a
                         warning if missing or --skip-checkm is given)

PYTHON DEPENDENCIES
-------------------------------------------------------
  Only needed if NOT using --taxdump-dir (see below):
    pip install ete3 six --break-system-packages
    (ete3.NCBITaxa downloads/builds a local NCBI taxonomy sqlite db on first use -- this
    can take several minutes and ~1-2 GB of disk the very first time the script is run,
    and requires outbound internet access.)

TAXONOMY LOOKUPS
-------------------------------------------------------
  --taxdump-dir <dir>  Point this at a local NCBI taxdump (nodes.dmp + names.dmp, e.g. an
                        extracted taxdump.tar.gz from
                        https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz). Fast,
                        offline, no ete3 dependency. Recommended, and required if the
                        machine running this has no outbound internet access.
                        Without it, falls back to ete3.NCBITaxa (see above).

EXAMPLES
-------------------------------------------------------
  # From an existing assembly:
  python metatax_binner.py \\
      -i contigs.fasta \\
      -d /dbs/nr.dmnd \\
      -o metatax_out \\
      -t 16 \\
      --taxdump-dir /dbs/taxdump

  # From raw paired-end reads (QC -> MEGAHIT -> classify -> bin):
  python metatax_binner.py \\
      --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz \\
      --trimmomatic-folder /path/to/Trimmomatic-0.39 \\
      -d /dbs/nr.dmnd \\
      -o metatax_out \\
      -t 16 \\
      --taxdump-dir /dbs/taxdump

OUTPUT LAYOUT
-------------------------------------------------------
  <outdir>/
    qc/proteins... (Trimmomatic + FLASH intermediates; only with --r1/--r2 and no --skip-qc)
    megahit/final.contigs.fa, megahit.log            (only with --r1/--r2)
    prodigal/proteins.faa, genes.gff
    diamond/hits.tsv
    classification/contig_classification.tsv
    bins/genus/<taxon>.fasta ... summary.tsv, bin_membership.tsv
    bins/species/...
"""

import argparse
import csv
import logging
import re
import shlex
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

# --------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------

# Both "domain" and "superkingdom" are included because NCBI renamed this top rank from
# "superkingdom" to "domain" in their taxdump around 2023 (reflecting the three-domain
# system). Whichever label your nodes.dmp actually uses will match; the other is a
# harmless no-op.
WANTED_RANKS = ["domain", "superkingdom", "phylum", "family", "genus", "species"]
BIN_RANKS_DEFAULT = ["genus", "species"]

DIAMOND_FIELDS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "stitle",
]

# Matches the trailing "[Organism name]" NCBI nr convention, e.g.
# "chromosomal replication initiator protein DnaA [Escherichia coli]" -> "Escherichia coli"
ORGANISM_RE = re.compile(r"\[([^\[\]]+)\]\s*$")


def parse_organism(stitle: str) -> str:
    """Pull the organism name out of a diamond stitle's trailing [brackets]."""
    if not stitle:
        return None
    m = ORGANISM_RE.search(stitle.strip())
    return m.group(1).strip() if m else None


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


def run_cmd(cmd, log_file: Path = None, cwd: Path = None) -> None:
    log.info("Running: %s", " ".join(str(c) for c in cmd))
    with open(log_file, "a") if log_file else open("/dev/null", "w") as lf:
        proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, cwd=cwd, text=True)
    if proc.returncode != 0:
        tail = ""
        if log_file and Path(log_file).exists():
            tail = "\n".join(Path(log_file).read_text().splitlines()[-30:])
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(str(c) for c in cmd)}\n{tail}")


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


# --------------------------------------------------------------------------------------
# Step 2: DIAMOND
# --------------------------------------------------------------------------------------

def run_diamond(query_faa: Path, db: Path, outdir: Path, threads: int, evalue: float,
                 max_target_seqs: int) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    hits_tsv = outdir / "hits.tsv"
    log_file = outdir / "diamond.log"
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
    __slots__ = ("sseqid", "pident", "length", "evalue", "bitscore", "organism")

    def __init__(self, sseqid, pident, length, evalue, bitscore, organism):
        self.sseqid = sseqid
        self.pident = pident
        self.length = length
        self.evalue = evalue
        self.bitscore = bitscore
        self.organism = organism


def parse_diamond_hits(hits_tsv: Path) -> dict:
    """Returns {orf_id: [Hit, ...]}"""
    hits_by_orf = defaultdict(list)
    idx = {f: i for i, f in enumerate(DIAMOND_FIELDS)}
    any_organism = False
    with open(hits_tsv) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(DIAMOND_FIELDS):
                continue
            organism = parse_organism(f[idx["stitle"]])
            if organism:
                any_organism = True
            hits_by_orf[f[idx["qseqid"]]].append(Hit(
                sseqid=f[idx["sseqid"]],
                pident=float(f[idx["pident"]]),
                length=int(f[idx["length"]]),
                evalue=float(f[idx["evalue"]]),
                bitscore=float(f[idx["bitscore"]]),
                organism=organism,
            ))
    if hits_by_orf and not any_organism:
        log.warning(
            "Could not parse an organism name (trailing '[...]') out of any DIAMOND "
            "stitle. Your nr.dmnd may not have been built from headers in standard "
            "NCBI nr format ('... [Genus species]'), so taxonomic classification will "
            "be empty."
        )
    return hits_by_orf


# --------------------------------------------------------------------------------------
# Step 3: Taxonomy lookup (local NCBI taxdump, or ete3 as a fallback)
# --------------------------------------------------------------------------------------

class LocalTaxdump:
    """Pure-python reader for a local NCBI taxdump (nodes.dmp + names.dmp) -- no network
    access, no ete3/sqlite build required. Point --taxdump-dir at a directory containing
    both files (e.g. the extraction of taxdump.tar.gz from
    https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)."""

    def __init__(self, taxdump_dir: Path):
        nodes_path = Path(taxdump_dir) / "nodes.dmp"
        names_path = Path(taxdump_dir) / "names.dmp"
        if not nodes_path.exists() or not names_path.exists():
            log.error("--taxdump-dir %s must contain both nodes.dmp and names.dmp", taxdump_dir)
            sys.exit(1)

        log.info("Loading local NCBI taxdump from %s (nodes.dmp + names.dmp)...", taxdump_dir)
        self.parent = {}   # taxid -> parent taxid
        self.rank = {}     # taxid -> rank string
        self.name_of = {}  # taxid -> scientific name
        self.taxid_of = {}  # lowercased scientific name -> taxid (first one wins)

        with open(nodes_path) as fh:
            for line in fh:
                f = line.split("|")
                taxid = int(f[0].strip())
                parent = int(f[1].strip())
                rank = f[2].strip()
                self.parent[taxid] = parent
                self.rank[taxid] = rank

        with open(names_path) as fh:
            for line in fh:
                f = line.split("|")
                taxid = int(f[0].strip())
                name = f[1].strip()
                name_class = f[3].strip()
                if name_class == "scientific name":
                    self.name_of[taxid] = name
                    self.taxid_of.setdefault(name.lower(), taxid)

        log.info("Loaded taxdump: %d taxa.", len(self.parent))

    def lineage(self, taxid) -> dict:
        result = {}
        seen = set()
        cur = taxid
        while cur is not None and cur not in seen:
            seen.add(cur)
            r = self.rank.get(cur)
            if r in WANTED_RANKS and r not in result:
                result[r] = self.name_of.get(cur)
            parent = self.parent.get(cur)
            if parent is None or parent == cur:
                break
            cur = parent
        return result


class TaxonomyLookup:
    """Lazily-initialized organism name -> {rank: name} lookup.

    Uses a local NCBI taxdump (nodes.dmp/names.dmp, via --taxdump-dir) if given -- fast,
    offline, no extra dependency. Otherwise falls back to ete3.NCBITaxa, which downloads
    and builds its own local sqlite taxonomy db (~/.etetoolkit/taxa.sqlite) on first use;
    that requires outbound internet access from wherever the script runs.

    Organism names come from DIAMOND stitle (the trailing "[Genus species]"), not from
    staxids, so the DIAMOND database itself does not need --taxonmap/--taxonnodes at all.
    """

    def __init__(self, taxdump_dir: Path = None):
        self._taxdump = LocalTaxdump(taxdump_dir) if taxdump_dir else None
        self._ncbi = None
        self._cache = {}

    def _ensure_ncbi(self):
        if self._ncbi is None:
            try:
                from ete3 import NCBITaxa
            except ImportError:
                log.error("ete3 is required for taxonomy lookups: pip install ete3 six --break-system-packages")
                sys.exit(1)
            log.info("Initializing NCBI taxonomy database (ete3) -- first run may take a while...")
            self._ncbi = NCBITaxa()

    def lineage_by_name(self, organism: str) -> dict:
        """Resolve an organism name (e.g. 'Escherichia coli') to its NCBI lineage."""
        if not organism:
            return {}
        if organism in self._cache:
            return self._cache[organism]

        # Names with strain/isolate qualifiers ("Escherichia coli str. K-12") usually
        # aren't in the taxdump as scientific names -- fall back to just the binomial.
        candidates = [organism]
        parts = organism.split()
        if len(parts) > 2:
            candidates.append(" ".join(parts[:2]))

        result = {}
        if self._taxdump is not None:
            for name in candidates:
                taxid = self._taxdump.taxid_of.get(name.lower())
                if taxid is not None:
                    result = self._taxdump.lineage(taxid)
                    break
        else:
            self._ensure_ncbi()
            try:
                for name in candidates:
                    taxids = self._ncbi.get_name_translator([name]).get(name)
                    if taxids:
                        taxid = taxids[0]
                        lineage_ids = self._ncbi.get_lineage(taxid)
                        ranks = self._ncbi.get_rank(lineage_ids)
                        names = self._ncbi.get_taxid_translator(lineage_ids)
                        for tid in lineage_ids:
                            r = ranks.get(tid)
                            if r in WANTED_RANKS:
                                result[r] = names.get(tid)
                        break
            except Exception as exc:
                log.debug("Lineage lookup failed for organism '%s': %s", organism, exc)
                result = {}

        self._cache[organism] = result
        return result


# --------------------------------------------------------------------------------------
# Step 4: Per-contig classification (bitscore-weighted majority/plurality vote per rank)
# --------------------------------------------------------------------------------------

def classify_contig(orf_ids, hits_by_orf: dict, taxlookup: TaxonomyLookup, ranks,
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
            if not h.organism:
                continue
            lineage = taxlookup.lineage_by_name(h.organism)
            for r in ranks:
                total_weight[r] += h.bitscore
                name = lineage.get(r)
                if name:
                    rank_weights[r][name] += h.bitscore

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


def classify_all_contigs(contig_ids, orf_to_contig, hits_by_orf, taxlookup, ranks,
                          bitscore_range, max_hits_per_orf, min_support) -> dict:
    contig_to_orfs = defaultdict(list)
    for orf_id, contig_id in orf_to_contig.items():
        contig_to_orfs[contig_id].append(orf_id)

    classifications = {}
    for i, contig_id in enumerate(contig_ids, 1):
        orfs = contig_to_orfs.get(contig_id, [])
        classifications[contig_id] = classify_contig(
            orfs, hits_by_orf, taxlookup, ranks, bitscore_range, max_hits_per_orf, min_support
        )
        if i % 500 == 0:
            log.info("Classified %d/%d contigs...", i, len(contig_ids))
    return classifications


def write_classification_table(classifications: dict, ranks, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        header = ["contig", "n_orfs", "n_orfs_with_hits"]
        for r in ranks:
            header += [r, f"{r}_support"]
        w.writerow(header)
        for contig_id, res in classifications.items():
            row = [contig_id, res["_n_orfs_total"], res["_n_orfs_with_hits"]]
            for r in ranks:
                taxon, support = res[r]
                row += [taxon, support]
            w.writerow(row)


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


def run_quast_multi(bin_fastas: dict, outdir: Path, threads: int) -> dict:
    """bin_fastas: {bin_name: Path}. Returns {bin_name: {stat: value}}."""
    if not bin_fastas:
        return {}
    if shutil.which("quast.py") is None:
        log.warning("quast.py not found on PATH -- falling back to built-in assembly stats.")
        return {name: basic_assembly_stats(p) for name, p in bin_fastas.items()}

    outdir.mkdir(parents=True, exist_ok=True)
    names = list(bin_fastas.keys())
    paths = [str(bin_fastas[n]) for n in names]
    cmd = [
        "quast.py", "-o", str(outdir), "--threads", str(threads),
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


def run_checkm(bin_dir: Path, outdir: Path, threads: int) -> dict:
    """Runs `checkm lineage_wf` on a directory of bin FASTAs. Returns {bin_name: {stat: value}}."""
    if shutil.which("checkm") is None:
        log.warning("checkm not found on PATH -- skipping completeness/contamination (will be NA).")
        return {}
    fastas = list(bin_dir.glob("*.fasta"))
    if not fastas:
        return {}
    outdir.mkdir(parents=True, exist_ok=True)
    results_tsv = outdir / "checkm_results.tsv"
    cmd = [
        "checkm", "lineage_wf", "-x", "fasta", "--tab_table", "-f", str(results_tsv),
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
                    "strain_heterogeneity": row.get("Strain heterogeneity"),
                }
    return stats


def summarize_bin_set(rank: str, rank_dir: Path, threads: int, skip_quast: bool,
                       skip_checkm: bool) -> None:
    bin_fastas = {p.stem: p for p in sorted(rank_dir.glob("*.fasta"))}
    if not bin_fastas:
        log.warning("No bins found for rank '%s'; skipping summary.", rank)
        return

    if skip_quast:
        quast_stats = {name: basic_assembly_stats(p) for name, p in bin_fastas.items()}
    else:
        quast_stats = run_quast_multi(bin_fastas, rank_dir / "quast_out", threads)

    checkm_stats = {} if skip_checkm else run_checkm(rank_dir, rank_dir / "checkm_out", threads)

    out_path = rank_dir / "summary.tsv"
    fields = ["bin", "num_contigs", "total_length_bp", "largest_contig_bp", "N50", "L50",
              "GC_percent", "completeness_percent", "contamination_percent", "strain_heterogeneity"]
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
                c.get("strain_heterogeneity", "NA"),
            ])
    log.info("Wrote %s", out_path)


# --------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Taxonomically classify and bin metagenomic contigs (optionally starting "
                     "from raw paired-end reads: QC + MEGAHIT -> Prodigal + DIAMOND + "
                     "weighted-vote LCA-style classifier + QUAST/CheckM bin summaries).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-i", "--input", type=Path, default=None,
                    help="Input contigs FASTA file (skips QC/MEGAHIT entirely). "
                         "Mutually exclusive with --r1/--r2.")
    p.add_argument("-d", "--diamond-db", required=True, type=Path, help="Path to nr.dmnd (DIAMOND protein DB)")
    p.add_argument("-o", "--outdir", required=True, type=Path, help="Output directory")
    p.add_argument("-t", "--threads", type=int, default=8)

    reads = p.add_argument_group("Raw reads input (runs QC + MEGAHIT to build contigs first)")
    reads.add_argument("--r1", type=Path, default=None, help="Raw/forward paired-end FASTQ (R1)")
    reads.add_argument("--r2", type=Path, default=None, help="Raw/reverse paired-end FASTQ (R2)")
    reads.add_argument("--skip-qc", action="store_true",
                        help="Feed --r1/--r2 straight to MEGAHIT, skipping the Trimmomatic/FLASH QC step "
                             "(--trim-polyg, if given, still runs).")

    pg = p.add_argument_group("Poly-G trimming (fastp; optional, runs before everything else)")
    pg.add_argument("--trim-polyg", action="store_true",
                     help="Trim poly-G tails with fastp before adapter clipping. Poly-G runs are a "
                          "known artifact of two-channel Illumina chemistry (NextSeq/NovaSeq) where "
                          "dark/no-signal cycles get miscalled as 'G' -- Trimmomatic's adapter/quality "
                          "trimming doesn't reliably catch these. Turn this on if your reads are known "
                          "to have long poly-G tails.")
    pg.add_argument("--fastp-cmd", default="fastp", help="fastp executable (only used with --trim-polyg)")
    pg.add_argument("--poly-g-min-len", type=int, default=10,
                     help="Minimum length of a 3' G-run to trim (fastp --poly_g_min_len, default 10)")

    qc = p.add_argument_group("QC options (Trimmomatic + FLASH; used with --r1/--r2 unless --skip-qc)")
    qc.add_argument("--trimmomatic-cmd", default="trimmomatic", help="Trimmomatic executable")
    qc.add_argument("--trimmomatic-folder", type=Path, default=None,
                     help="Path to the Trimmomatic install folder containing "
                          "adapters/TruSeq3-PE-2.fa (required unless --skip-qc)")
    qc.add_argument("--flash-cmd", default="flash", help="FLASH executable")
    qc.add_argument("--flash-max-overlap", type=int, default=150)
    qc.add_argument("--pigz-cmd", default="pigz", help="pigz executable (falls back to gzip if missing)")
    qc.add_argument("--qc-quality", type=int, default=20,
                     help="Trimmomatic SLIDINGWINDOW:4:<qc_quality> for all quality-trim passes")
    qc.add_argument("--qc-minlen", type=int, default=50, help="Trimmomatic MINLEN for quality-trim passes")
    qc.add_argument("--keep-qc-tmp", action="store_true", help="Keep intermediate QC files")

    mh = p.add_argument_group("MEGAHIT options (used with --r1/--r2)")
    mh.add_argument("--megahit-cmd", default="megahit", help="MEGAHIT executable")
    mh.add_argument("--megahit-min-contig-len", type=int, default=None,
                     help="MEGAHIT --min-contig-len (default: MEGAHIT's own default, 200bp)")
    mh.add_argument("--megahit-extra", default=None,
                     help="Extra raw arguments passed through to MEGAHIT verbatim, "
                          "e.g. --megahit-extra '--k-list 21,41,61'")

    p.add_argument("--prodigal-mode", choices=["single", "meta"], default="meta")

    p.add_argument("-e", "--evalue", type=float, default=1e-5, help="DIAMOND e-value cutoff")
    p.add_argument("--max-target-seqs", type=int, default=25, help="DIAMOND -k (hits kept per ORF)")

    p.add_argument("--taxdump-dir", type=Path, default=None,
                    help="Directory containing a local NCBI taxdump (nodes.dmp + names.dmp, "
                         "e.g. an extracted taxdump.tar.gz). If given, taxonomy lookups use "
                         "this directly -- no network access or ete3 sqlite build needed. "
                         "Otherwise falls back to ete3.NCBITaxa (downloads its own copy).")

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

    p.add_argument("--reuse-prodigal", action="store_true",
                    help="Skip Prodigal if <outdir>/prodigal/proteins.faa already exists (reuse it).")
    p.add_argument("--reuse-diamond", action="store_true",
                    help="Skip DIAMOND if <outdir>/diamond/hits.tsv already exists (reuse it). "
                         "Combine with --reuse-prodigal and a new --ranks/--min-support to "
                         "re-classify/re-bin at a different granularity without re-running "
                         "Prodigal+DIAMOND.")

    p.add_argument("--skip-quast", action="store_true")
    p.add_argument("--skip-checkm", action="store_true")

    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    setup_logging(args.verbose)

    using_reads = args.r1 is not None or args.r2 is not None
    if using_reads and args.input is not None:
        log.error("Pass either -i/--input OR --r1/--r2, not both.")
        sys.exit(1)
    if using_reads and (args.r1 is None or args.r2 is None):
        log.error("Both --r1 and --r2 are required together.")
        sys.exit(1)
    if not using_reads and args.input is None:
        log.error("Provide either -i/--input (pre-assembled contigs) or --r1/--r2 (raw paired FASTQs).")
        sys.exit(1)
    if using_reads and not args.skip_qc and args.trimmomatic_folder is None:
        log.error("--trimmomatic-folder is required for QC (or pass --skip-qc to bypass QC).")
        sys.exit(1)

    ranks = [r.strip() for r in args.ranks.split(",") if r.strip()]
    for r in ranks:
        if r not in WANTED_RANKS:
            log.error("Unsupported rank '%s'. Supported: %s", r, WANTED_RANKS)
            sys.exit(1)

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    faa = outdir / "prodigal" / "proteins.faa"
    reuse_prodigal = args.reuse_prodigal and faa.exists()
    hits_tsv = outdir / "diamond" / "hits.tsv"
    reuse_diamond = args.reuse_diamond and hits_tsv.exists()

    if using_reads:
        if args.trim_polyg:
            which_or_die(args.fastp_cmd)
        if not args.skip_qc:
            which_or_die(args.trimmomatic_cmd)
            which_or_die(args.flash_cmd)
        which_or_die(args.megahit_cmd)
    if not reuse_prodigal:
        which_or_die("prodigal")
    if not reuse_diamond:
        which_or_die("diamond")
    if not args.skip_quast and shutil.which("quast.py") is None:
        log.warning("quast.py not on PATH; will use built-in assembly stats instead.")
    if not args.skip_checkm and shutil.which("checkm") is None:
        log.warning("checkm not on PATH; completeness/contamination will be NA.")

    n_qc_steps = (
        (1 if using_reads and args.trim_polyg else 0)
        + (1 if using_reads and not args.skip_qc else 0)
        + (1 if using_reads else 0)
    )
    steps = StepCounter(n_qc_steps + 5)

    # 0a/0b/0c. Poly-G trim + QC + MEGAHIT assembly (raw-reads mode only)
    assembly_fasta = args.input
    if using_reads:
        r1_in, r2_in = args.r1, args.r2
        if args.trim_polyg:
            steps.next("Trimming poly-G tails (fastp)...")
            r1_in, r2_in = trim_poly_g(
                r1_in, r2_in, outdir / "polyg", args.threads,
                args.fastp_cmd, args.poly_g_min_len,
            )

        if not args.skip_qc:
            steps.next("Running QC (Trimmomatic adapter/quality trim + FLASH merge)...")
            r1_final, r2_final, u_final = run_qc(
                r1_in, r2_in, outdir / "qc", args.threads,
                args.trimmomatic_cmd, args.trimmomatic_folder,
                args.flash_cmd, args.flash_max_overlap, args.pigz_cmd,
                args.qc_quality, args.qc_minlen, keep_tmp=args.keep_qc_tmp,
            )
        else:
            r1_final, r2_final, u_final = r1_in, r2_in, None

        steps.next("Running MEGAHIT assembly...")
        assembly_fasta = run_megahit(
            r1_final, r2_final, u_final, outdir / "megahit", args.threads,
            args.megahit_cmd, args.megahit_min_contig_len, args.megahit_extra,
        )
        log.info("MEGAHIT assembly: %s", assembly_fasta)

    # 1. Prodigal
    if reuse_prodigal:
        steps.next(f"Reusing existing Prodigal output: {faa}")
    else:
        steps.next("Running Prodigal...")
        faa, _gff = run_prodigal(assembly_fasta, outdir / "prodigal", mode=args.prodigal_mode)
    orf_to_contig = parse_orf_to_contig(faa)
    log.info("Predicted %d ORFs.", len(orf_to_contig))

    # 2. DIAMOND
    if reuse_diamond:
        steps.next(f"Reusing existing DIAMOND output: {hits_tsv}")
    else:
        steps.next(f"Running DIAMOND blastp vs {args.diamond_db}...")
        hits_tsv = run_diamond(
            faa, args.diamond_db, outdir / "diamond", args.threads, args.evalue,
            args.max_target_seqs,
        )
    hits_by_orf = parse_diamond_hits(hits_tsv)
    log.info("Got hits for %d/%d ORFs.", len(hits_by_orf), len(orf_to_contig))

    # 3. Classification
    steps.next(f"Classifying contigs (bitscore-weighted vote at ranks: {ranks})...")
    contig_seqs = read_fasta(assembly_fasta)
    taxlookup = TaxonomyLookup(taxdump_dir=args.taxdump_dir)
    classifications = classify_all_contigs(
        list(contig_seqs.keys()), orf_to_contig, hits_by_orf, taxlookup, ranks,
        args.bitscore_range, args.max_hits_per_orf, args.min_support,
    )
    write_classification_table(classifications, ranks, outdir / "classification" / "contig_classification.tsv")

    # 4. Binning
    steps.next("Writing bin FASTA files...")
    membership = bin_contigs(
        contig_seqs, classifications, ranks, outdir / "bins",
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

    # 5. Per-bin-set summaries
    steps.next("Summarizing bins (QUAST contiguity + CheckM completeness/contamination)...")
    for r in ranks:
        summarize_bin_set(r, outdir / "bins" / r, args.threads, args.skip_quast, args.skip_checkm)

    log.info("Done. Results in %s", outdir)


if __name__ == "__main__":
    main()
