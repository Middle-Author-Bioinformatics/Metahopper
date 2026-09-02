#!/usr/bin/env python3
"""Summarize MetaHopper bins across samples and taxonomic ranks.

The expected layout is:

    metahop_SAMPLE/bins/RANK/TAXON.fasta

If a completed two-pass run contains ``final/bins``, ``--bin-set auto`` uses that
directory instead. The script uses only the Python standard library.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import statistics
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path


RANK_ORDER = {
    "domain": 0,
    "superkingdom": 1,
    "kingdom": 2,
    "phylum": 3,
    "class": 4,
    "order": 5,
    "family": 6,
    "genus": 7,
    "species": 8,
}

FASTA_ENDINGS = (
    ".fasta", ".fa", ".fna", ".fas",
    ".fasta.gz", ".fa.gz", ".fna.gz", ".fas.gz",
)

BIN_FIELDS = [
    "sample",
    "bin_set",
    "rank",
    "bin_id",
    "taxon_label",
    "is_unclassified",
    "num_contigs",
    "total_bp",
    "smallest_contig_bp",
    "largest_contig_bp",
    "mean_contig_bp",
    "median_contig_bp",
    "N50_bp",
    "L50_contigs",
    "N90_bp",
    "L90_contigs",
    "GC_percent",
    "contigs_per_Mbp",
    "largest_contig_fraction",
    "N50_fraction",
    "fasta",
]

RANK_FIELDS = [
    "sample",
    "bin_set",
    "rank",
    "num_bins",
    "num_classified_bins",
    "total_contigs",
    "total_bp",
    "classified_contigs",
    "classified_bp",
    "unclassified_contigs",
    "unclassified_bp",
    "classified_fraction_contigs",
    "classified_fraction_bp",
    "total_contigs_per_Mbp",
    "median_classified_bin_bp",
    "median_classified_bin_contigs",
    "median_classified_bin_N50_bp",
    "largest_classified_bin",
    "largest_classified_bin_bp",
    "most_fragmented_classified_bin",
    "highest_contigs_per_Mbp",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate per-bin and per-rank size and fragmentation statistics from "
            "one or more MetaHopper result directories."
        )
    )
    parser.add_argument(
        "-i", "--input-dir", type=Path, default=Path("."),
        help="Directory containing metahop_* result folders, or one result folder (default: .)",
    )
    parser.add_argument(
        "-o", "--output-dir", type=Path, default=Path("metahopper_bin_summary"),
        help="Output directory (default: metahopper_bin_summary)",
    )
    parser.add_argument(
        "--sample-glob", default="metahop_*",
        help="Glob used to find sample result directories (default: metahop_*)",
    )
    parser.add_argument(
        "--strip-prefix", default="metahop_",
        help="Prefix removed from result-directory names to form sample IDs (default: metahop_)",
    )
    parser.add_argument(
        "--bin-set", choices=("auto", "preliminary", "final", "both"), default="auto",
        help=(
            "Bins to inspect: auto prefers final/bins when present; preliminary uses bins; "
            "final uses final/bins; both reports both (default: auto)"
        ),
    )
    parser.add_argument(
        "--ranks", default=None,
        help="Optional comma-separated ranks; default discovers all rank directories",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=4,
        help="FASTA files processed in parallel (default: 4; use 1 to minimize memory)",
    )
    return parser.parse_args()


def is_fasta(path: Path) -> bool:
    return path.is_file() and path.name.lower().endswith(FASTA_ENDINGS)


def has_fastas(bin_root: Path) -> bool:
    if not bin_root.is_dir():
        return False
    return any(is_fasta(path) for rank_dir in bin_root.iterdir() if rank_dir.is_dir()
               for path in rank_dir.iterdir())


def sample_id(sample_dir: Path, strip_prefix: str) -> str:
    name = sample_dir.name
    return name[len(strip_prefix):] if strip_prefix and name.startswith(strip_prefix) else name


def discover_samples(input_dir: Path, sample_glob: str) -> list[Path]:
    input_dir = input_dir.resolve()
    if has_fastas(input_dir / "bins") or has_fastas(input_dir / "final" / "bins"):
        return [input_dir]
    samples = sorted(path for path in input_dir.glob(sample_glob) if path.is_dir())
    return [path for path in samples
            if has_fastas(path / "bins") or has_fastas(path / "final" / "bins")]


def select_bin_roots(sample_dir: Path, mode: str) -> list[tuple[str, Path]]:
    preliminary = sample_dir / "bins"
    final = sample_dir / "final" / "bins"
    if mode == "auto":
        if has_fastas(final):
            return [("final", final)]
        return [("preliminary", preliminary)] if has_fastas(preliminary) else []
    if mode == "preliminary":
        return [("preliminary", preliminary)] if has_fastas(preliminary) else []
    if mode == "final":
        return [("final", final)] if has_fastas(final) else []
    roots = []
    if has_fastas(preliminary):
        roots.append(("preliminary", preliminary))
    if has_fastas(final):
        roots.append(("final", final))
    return roots


def bin_id_from_path(path: Path) -> str:
    name = path.name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    for ending in (".fasta", ".fna", ".fas", ".fa"):
        if name.lower().endswith(ending):
            return name[:-len(ending)]
    return Path(name).stem


def rank_sort_key(rank: str) -> tuple[int, str]:
    return RANK_ORDER.get(rank.lower(), 999), rank.lower()


def collect_tasks(samples: list[Path], args: argparse.Namespace) -> list[dict]:
    wanted_ranks = None
    if args.ranks:
        wanted_ranks = {rank.strip().lower() for rank in args.ranks.split(",") if rank.strip()}

    tasks = []
    for sample_dir in samples:
        sid = sample_id(sample_dir, args.strip_prefix)
        for bin_set, bin_root in select_bin_roots(sample_dir, args.bin_set):
            rank_dirs = sorted((path for path in bin_root.iterdir() if path.is_dir()),
                               key=lambda path: rank_sort_key(path.name))
            for rank_dir in rank_dirs:
                rank = rank_dir.name
                if wanted_ranks is not None and rank.lower() not in wanted_ranks:
                    continue
                for fasta in sorted(path for path in rank_dir.iterdir() if is_fasta(path)):
                    tasks.append({
                        "sample": sid,
                        "bin_set": bin_set,
                        "rank": rank,
                        "bin_id": bin_id_from_path(fasta),
                        "fasta": str(fasta.resolve()),
                    })
    return tasks


def open_fasta(path: Path):
    return gzip.open(path, "rt") if path.name.lower().endswith(".gz") else open(path)


def nx_lx(lengths_desc: list[int], total_bp: int, fraction: float) -> tuple[int, int]:
    if not lengths_desc or total_bp <= 0:
        return 0, 0
    target = total_bp * fraction
    cumulative = 0
    for index, length in enumerate(lengths_desc, 1):
        cumulative += length
        if cumulative >= target:
            return length, index
    return lengths_desc[-1], len(lengths_desc)


def analyze_fasta(task: dict) -> dict:
    path = Path(task["fasta"])
    lengths: list[int] = []
    current_length = 0
    total_gc = 0
    saw_header = False

    with open_fasta(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if saw_header:
                    lengths.append(current_length)
                saw_header = True
                current_length = 0
            else:
                if not saw_header:
                    raise ValueError(f"Sequence encountered before a FASTA header in {path}")
                sequence = line.upper()
                current_length += len(sequence)
                total_gc += sequence.count("G") + sequence.count("C")
    if saw_header:
        lengths.append(current_length)

    lengths.sort(reverse=True)
    num_contigs = len(lengths)
    total_bp = sum(lengths)
    n50, l50 = nx_lx(lengths, total_bp, 0.50)
    n90, l90 = nx_lx(lengths, total_bp, 0.90)
    largest = lengths[0] if lengths else 0
    smallest = lengths[-1] if lengths else 0
    mean_length = total_bp / num_contigs if num_contigs else 0.0
    median_length = statistics.median(lengths) if lengths else 0.0
    contigs_per_mbp = num_contigs / (total_bp / 1_000_000) if total_bp else 0.0
    bin_id = task["bin_id"]

    return {
        "sample": task["sample"],
        "bin_set": task["bin_set"],
        "rank": task["rank"],
        "bin_id": bin_id,
        "taxon_label": bin_id.replace("_", " "),
        "is_unclassified": str(bin_id.lower() == "unclassified").lower(),
        "num_contigs": num_contigs,
        "total_bp": total_bp,
        "smallest_contig_bp": smallest,
        "largest_contig_bp": largest,
        "mean_contig_bp": round(mean_length, 2),
        "median_contig_bp": round(float(median_length), 2),
        "N50_bp": n50,
        "L50_contigs": l50,
        "N90_bp": n90,
        "L90_contigs": l90,
        "GC_percent": round(100 * total_gc / total_bp, 4) if total_bp else 0.0,
        "contigs_per_Mbp": round(contigs_per_mbp, 4),
        "largest_contig_fraction": round(largest / total_bp, 6) if total_bp else 0.0,
        "N50_fraction": round(n50 / total_bp, 6) if total_bp else 0.0,
        "fasta": str(path),
    }


def safe_fraction(numerator: int, denominator: int) -> float:
    return round(numerator / denominator, 6) if denominator else 0.0


def median_value(rows: list[dict], field: str) -> float:
    return round(float(statistics.median(row[field] for row in rows)), 2) if rows else 0.0


def build_rank_summaries(bin_rows: list[dict]) -> list[dict]:
    grouped: dict[tuple[str, str, str], list[dict]] = {}
    for row in bin_rows:
        grouped.setdefault((row["sample"], row["bin_set"], row["rank"]), []).append(row)

    summaries = []
    for (sample, bin_set, rank), rows in sorted(
        grouped.items(), key=lambda item: (item[0][0], item[0][1], rank_sort_key(item[0][2]))
    ):
        classified = [row for row in rows if row["is_unclassified"] == "false"]
        unclassified = [row for row in rows if row["is_unclassified"] == "true"]
        total_contigs = sum(row["num_contigs"] for row in rows)
        total_bp = sum(row["total_bp"] for row in rows)
        classified_contigs = sum(row["num_contigs"] for row in classified)
        classified_bp = sum(row["total_bp"] for row in classified)
        unclassified_contigs = sum(row["num_contigs"] for row in unclassified)
        unclassified_bp = sum(row["total_bp"] for row in unclassified)
        largest = max(classified, key=lambda row: row["total_bp"], default=None)
        most_fragmented = max(classified, key=lambda row: row["contigs_per_Mbp"], default=None)

        summaries.append({
            "sample": sample,
            "bin_set": bin_set,
            "rank": rank,
            "num_bins": len(rows),
            "num_classified_bins": len(classified),
            "total_contigs": total_contigs,
            "total_bp": total_bp,
            "classified_contigs": classified_contigs,
            "classified_bp": classified_bp,
            "unclassified_contigs": unclassified_contigs,
            "unclassified_bp": unclassified_bp,
            "classified_fraction_contigs": safe_fraction(classified_contigs, total_contigs),
            "classified_fraction_bp": safe_fraction(classified_bp, total_bp),
            "total_contigs_per_Mbp": round(total_contigs / (total_bp / 1_000_000), 4)
            if total_bp else 0.0,
            "median_classified_bin_bp": median_value(classified, "total_bp"),
            "median_classified_bin_contigs": median_value(classified, "num_contigs"),
            "median_classified_bin_N50_bp": median_value(classified, "N50_bp"),
            "largest_classified_bin": largest["bin_id"] if largest else "NA",
            "largest_classified_bin_bp": largest["total_bp"] if largest else 0,
            "most_fragmented_classified_bin": most_fragmented["bin_id"]
            if most_fragmented else "NA",
            "highest_contigs_per_Mbp": most_fragmented["contigs_per_Mbp"]
            if most_fragmented else 0.0,
        })
    return summaries


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_matrix(path: Path, rows: list[dict], value_field: str) -> None:
    row_keys = sorted({(row["sample"], row["bin_set"]) for row in rows})
    columns = sorted(
        {(row["rank"], row["bin_id"]) for row in rows},
        key=lambda item: (rank_sort_key(item[0]), item[1].lower()),
    )
    values = {
        (row["sample"], row["bin_set"], row["rank"], row["bin_id"]): row[value_field]
        for row in rows
    }
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "bin_set"] + [f"{rank}::{bin_id}" for rank, bin_id in columns])
        for sample, bin_set in row_keys:
            writer.writerow(
                [sample, bin_set]
                + [values.get((sample, bin_set, rank, bin_id), 0) for rank, bin_id in columns]
            )


def main() -> int:
    args = parse_args()
    if args.threads < 1:
        print("ERROR: --threads must be at least 1", file=sys.stderr)
        return 2
    samples = discover_samples(args.input_dir, args.sample_glob)
    if not samples:
        print(
            f"ERROR: no MetaHopper result directories found under {args.input_dir}",
            file=sys.stderr,
        )
        return 1

    tasks = collect_tasks(samples, args)
    if not tasks:
        print("ERROR: result directories were found, but no matching bin FASTAs were found", file=sys.stderr)
        return 1

    print(
        f"Found {len(samples)} sample(s) and {len(tasks)} bin FASTA file(s); "
        f"processing with {args.threads} worker(s).",
        file=sys.stderr,
    )
    if args.threads == 1:
        bin_rows = [analyze_fasta(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            bin_rows = list(executor.map(analyze_fasta, tasks))

    bin_rows.sort(
        key=lambda row: (
            row["sample"], row["bin_set"], rank_sort_key(row["rank"]),
            row["is_unclassified"] == "true", row["bin_id"].lower(),
        )
    )
    rank_rows = build_rank_summaries(bin_rows)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output_dir / "bin_metrics.tsv", bin_rows, BIN_FIELDS)
    write_tsv(args.output_dir / "rank_summary.tsv", rank_rows, RANK_FIELDS)
    write_matrix(args.output_dir / "bin_size_bp_matrix.tsv", bin_rows, "total_bp")
    write_matrix(args.output_dir / "bin_contig_count_matrix.tsv", bin_rows, "num_contigs")
    write_matrix(args.output_dir / "bin_N50_bp_matrix.tsv", bin_rows, "N50_bp")
    write_matrix(
        args.output_dir / "bin_contigs_per_Mbp_matrix.tsv", bin_rows, "contigs_per_Mbp"
    )

    print(f"Wrote summaries to {args.output_dir.resolve()}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
