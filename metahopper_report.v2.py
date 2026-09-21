#!/usr/bin/env python3
"""Build a self-contained interactive HTML report from one MetaHopper output directory.

    metahopper_report.py -i metahop_out -o report.html
    metahopper_report.py -i metahop_out --apply-split manual_bins.tsv

Point ``-i`` at a MetaHopper run directory. Everything the report needs is discovered
inside it; nothing has to be pre-summarised. Files that a given run did not produce are
simply omitted from the report, so contigs-only runs, ``--skip-reassembly`` runs, and
full two-pass runs all work.

Files read (all optional)
-------------------------
    binarena/binarena_input.tsv        or final/binarena/binarena_input.tsv
    classification/contig_classification.tsv   or final/classification/...
    classification/bin_refinement.tsv          or final/classification/...
    bins/<rank>/summary.tsv                    or final/bins/<rank>/summary.tsv
    bins/<rank>/bin_membership.tsv             or final/bins/<rank>/...
    bins/<rank>/reassembly_summary.tsv
    final/assembly/contig_provenance.tsv

The BinaRena tab is an interactive canvas scatter. Choose any two numeric columns as
axes (GC, coverage, and every k-mer PCA/t-SNE/UMAP axis) and colour by any categorical
column, then lasso contigs and assign the selection to a named sub-bin.

    Axis scaling   Each axis has an independent scale: cube root, square root, linear,
                   square, cube, or log10. Powers are applied as sign(v)*|v|^p so that
                   ordination axes, which straddle zero, transform symmetrically rather
                   than folding. log10 keeps only positive values and reports how many
                   points it hid. Tick labels always read in the column's original
                   units; only the positions are transformed.

    Search         Substring match over the contig name and every categorical column.
                   Comma-separated terms are OR'd ("Wolbachia, Symbiopectobacterium")
                   and a leading '-' excludes ("-Unclassified"). It filters rather than
                   highlights, so a lasso only ever picks up what is currently visible.

Sub-bins can be exported two ways. "Download manual_bins.tsv" writes the assignment
table, which this same script turns back into FASTA on the server:

    metahopper_report.py -i metahop_out --apply-split manual_bins.tsv

Or select the run's contig FASTA in the report and "Write FASTA per sub-bin" produces
the sequences directly in the browser -- one .fasta for a single sub-bin, a .zip for
several. The file is read locally with the File API in 8 MB chunks and never uploaded,
so a multi-GB assembly streams without being held in memory. Sequences are deliberately
not embedded in the HTML, which is what keeps the report small. Record names are matched
on the first whitespace-delimited token, so MEGAHIT- and Unicycler-decorated headers work
unchanged; note that a final/ run renames contigs, so it needs
final/assembly/consolidated_contigs.fasta. Gzipped FASTA cannot be read in the browser --
use --apply-split for that.

Only the Python standard library is used, and the HTML has no external dependencies, so
the report works offline and on an air-gapped cluster.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import sys
from datetime import datetime
from pathlib import Path

RANK_ORDER = {
    "domain": 0, "superkingdom": 1, "kingdom": 2, "phylum": 3,
    "class": 4, "order": 5, "family": 6, "genus": 7, "species": 8,
}

# Columns treated as categorical (colour/filter) rather than numeric (axes), regardless
# of whether their values happen to parse as numbers.
FORCED_CATEGORICAL = {
    "bin", "triage", "refinement_status", "gc_outlier", "coverage_outlier",
    "joint_outlier", "retained", "source_type", "source_bin",
}

PREFERRED_X = ["4PC1", "5PC1", "6PC1", "GC", "4tsne1", "4UM1", "length"]
PREFERRED_Y = ["4PC2", "5PC2", "6PC2", "coverage", "4tsne2", "4UM2", "GC"]


# ------------------------------------------------------------------ small helpers

def rank_key(rank: str):
    return (RANK_ORDER.get(rank.lower(), 99), rank.lower())


def read_tsv(path: Path):
    """Read a TSV into a list of dicts, tolerating the trailing \\r our writers emit."""
    if not path or not path.is_file():
        return []
    with open(path, newline="") as fh:
        rows = []
        for row in csv.DictReader(fh, delimiter="\t"):
            rows.append({
                (k.strip() if k else k): (v.strip() if isinstance(v, str) else v)
                for k, v in row.items()
            })
        return rows


def to_float(value):
    if value is None:
        return None
    text = str(value).strip()
    if text == "" or text.upper() in ("NA", "NAN", "NONE", "-"):
        return None
    try:
        result = float(text)
    except ValueError:
        return None
    if result != result or result in (float("inf"), float("-inf")):
        return None
    return result


def open_fasta(path: Path):
    return gzip.open(path, "rt") if str(path).lower().endswith(".gz") else open(path)


def iter_fasta(path: Path):
    name, chunks = None, []
    with open_fasta(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def human_bp(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return "NA"
    for unit, size in (("Gb", 1e9), ("Mb", 1e6), ("kb", 1e3)):
        if abs(value) >= size:
            return f"{value / size:.2f} {unit}"
    return f"{value:.0f} bp"


# ------------------------------------------------------------------ discovery

class RunLayout:
    """Locate the files a MetaHopper run produced, preferring final/ when present."""

    def __init__(self, root: Path, prefer: str = "auto"):
        self.root = root.resolve()
        self.has_final = (self.root / "final" / "bins").is_dir()
        self.has_reassembly = (self.root / "reassembly").is_dir()
        if prefer == "final" and not self.has_final:
            raise SystemExit(
                f"--bin-set final was requested but {self.root/'final'/'bins'} does not "
                "exist. This run had no seed-and-extension stage (no reads, or "
                "--skip-reassembly)."
            )
        self.stage = "final" if (prefer == "auto" and self.has_final) or prefer == "final" \
            else "preliminary"
        base = self.root / "final" if self.stage == "final" else self.root
        self.base = base
        self.bins_dir = base / "bins"
        self.binarena = self._first(base / "binarena" / "binarena_input.tsv")
        self.classification = self._first(
            base / "classification" / "contig_classification.tsv")
        self.refinement = self._first(base / "classification" / "bin_refinement.tsv")
        self.provenance = self._first(self.root / "final" / "assembly" / "contig_provenance.tsv")
        self.consolidated = self._first(
            self.root / "final" / "assembly" / "consolidated_contigs.fasta")
        self.megahit = self._first(self.root / "megahit" / "final.contigs.fa")

    @staticmethod
    def _first(path: Path):
        return path if path.is_file() else None

    def ranks(self):
        if not self.bins_dir.is_dir():
            return []
        return sorted(
            (p.name for p in self.bins_dir.iterdir() if p.is_dir()), key=rank_key,
        )

    def default_contigs(self):
        return self.consolidated or self.megahit

    def bin_fasta_partition(self):
        """One rank's bin FASTAs, which together contain every binned contig once.

        Used as a contig source for --apply-split when neither a consolidated assembly
        nor a MEGAHIT assembly is present, i.e. when contigs were supplied with -i from
        somewhere outside the run directory.
        """
        for rank in self.ranks():
            fastas = sorted((self.bins_dir / rank).glob("*.fasta"))
            if fastas:
                return fastas
        return []


# ------------------------------------------------------------------ bin metrics

def fasta_stats(path: Path):
    lengths, gc = [], 0
    for _name, seq in iter_fasta(path):
        seq = seq.upper()
        lengths.append(len(seq))
        gc += seq.count("G") + seq.count("C")
    lengths.sort(reverse=True)
    total = sum(lengths)
    n50 = l50 = 0
    cumulative = 0
    for i, length in enumerate(lengths, 1):
        cumulative += length
        if cumulative >= total / 2:
            n50, l50 = length, i
            break
    return {
        "num_contigs": len(lengths),
        "total_length_bp": total,
        "largest_contig_bp": lengths[0] if lengths else 0,
        "N50": n50, "L50": l50,
        "GC_percent": round(100.0 * gc / total, 2) if total else 0.0,
    }


def collect_bins(layout: RunLayout, compute_missing: bool):
    """Per-bin rows, from summary.tsv where available and the FASTA otherwise."""
    rows = []
    for rank in layout.ranks():
        rank_dir = layout.bins_dir / rank
        summary = {r.get("bin"): r for r in read_tsv(rank_dir / "summary.tsv")}
        reassembly = {r.get("bin"): r for r in read_tsv(layout.root / "bins" / rank / "reassembly_summary.tsv")}
        for fasta in sorted(rank_dir.glob("*.fasta")):
            name = fasta.stem
            record = {
                "rank": rank, "bin": name,
                "unclassified": name.lower() == "unclassified",
            }
            src = summary.get(name)
            if src:
                for field in ("num_contigs", "total_length_bp", "largest_contig_bp",
                              "N50", "L50", "GC_percent", "completeness_percent",
                              "contamination_percent", "strain_heterogeneity_percent"):
                    record[field] = to_float(src.get(field))
                record["marker_lineage"] = src.get("marker_lineage") or "NA"
            if record.get("total_length_bp") in (None, 0) and compute_missing:
                record.update(fasta_stats(fasta))
            record.setdefault("completeness_percent", None)
            record.setdefault("contamination_percent", None)
            reas = reassembly.get(name)
            if reas:
                record["reassembly_outcome"] = reas.get("outcome") or "NA"
                record["recovered_fraction"] = to_float(reas.get("recovered_fraction"))
                record["contigs_before"] = to_float(reas.get("contigs_before"))
                record["total_length_before_bp"] = to_float(reas.get("total_length_before_bp"))
            total = record.get("total_length_bp") or 0
            contigs = record.get("num_contigs") or 0
            record["contigs_per_Mbp"] = round(contigs / (total / 1e6), 2) if total else 0.0
            rows.append(record)
    rows.sort(key=lambda r: (rank_key(r["rank"]), r["unclassified"],
                             -(r.get("total_length_bp") or 0)))
    return rows


def rank_totals(bin_rows):
    grouped = {}
    for row in bin_rows:
        grouped.setdefault(row["rank"], []).append(row)
    out = []
    for rank in sorted(grouped, key=rank_key):
        rows = grouped[rank]
        classified = [r for r in rows if not r["unclassified"]]
        total_bp = sum(r.get("total_length_bp") or 0 for r in rows)
        cls_bp = sum(r.get("total_length_bp") or 0 for r in classified)
        out.append({
            "rank": rank,
            "num_bins": len(rows),
            "num_classified_bins": len(classified),
            "total_bp": total_bp,
            "classified_bp": cls_bp,
            "unclassified_bp": total_bp - cls_bp,
            "classified_fraction": round(cls_bp / total_bp, 4) if total_bp else 0.0,
            "total_contigs": sum(r.get("num_contigs") or 0 for r in rows),
            "largest_bin": classified[0]["bin"] if classified else "NA",
            "largest_bin_bp": classified[0].get("total_length_bp") if classified else 0,
        })
    return out


# ------------------------------------------------------------------ scatter dataset

def build_scatter(layout: RunLayout, min_length: int, max_points: int,
                  require_coords: bool):
    """Columnar dataset for the BinaRena tab.

    Categorical columns are dictionary-encoded and numeric columns are emitted as plain
    arrays, which keeps the embedded JSON an order of magnitude smaller than a list of
    per-contig objects. Contigs are prioritised by length so that trimming to
    ``max_points`` keeps the sequence that matters.
    """
    if not layout.binarena:
        return None
    rows = read_tsv(layout.binarena)
    if not rows:
        return None

    header = list(rows[0].keys())
    id_col = header[0]
    coord_cols = [c for c in header if any(
        c.endswith(suffix) for suffix in ("PC1", "PC2", "tsne1", "tsne2", "UM1", "UM2"))]

    kept = []
    n_total = len(rows)
    for row in rows:
        length = to_float(row.get("length")) or 0
        if length < min_length:
            continue
        if require_coords and coord_cols and not any(
                to_float(row.get(c)) is not None for c in coord_cols):
            continue
        kept.append((length, row))
    n_eligible = len(kept)
    kept.sort(key=lambda item: -item[0])
    trimmed = n_eligible > max_points
    kept = kept[:max_points]
    rows = [row for _length, row in kept]

    numeric, categorical = {}, {}
    for col in header:
        if col == id_col:
            continue
        values = [row.get(col) for row in rows]
        is_taxon_label = col.startswith("taxon_") and col != "taxon_support"
        if col in FORCED_CATEGORICAL or is_taxon_label:
            categorical[col] = values
            continue
        floats = [to_float(v) for v in values]
        if sum(1 for f in floats if f is not None) >= max(1, 0.5 * len(floats)):
            numeric[col] = floats
        else:
            categorical[col] = values

    cat_encoded = {}
    for col, values in categorical.items():
        levels, codes = [], []
        index = {}
        for value in values:
            value = value if value not in (None, "") else "NA"
            if value not in index:
                index[value] = len(levels)
                levels.append(value)
            codes.append(index[value])
        cat_encoded[col] = {"levels": levels, "codes": codes}

    numeric_names = sorted(numeric)
    def pick(preferred, fallback_index):
        for name in preferred:
            if name in numeric:
                return name
        return numeric_names[min(fallback_index, len(numeric_names) - 1)] \
            if numeric_names else None

    return {
        "ids": [row.get(id_col, "") for row in rows],
        "numeric": numeric,
        "categorical": cat_encoded,
        "numericNames": numeric_names,
        "categoricalNames": sorted(cat_encoded),
        "defaultX": pick(PREFERRED_X, 0),
        "defaultY": pick(PREFERRED_Y, 1),
        "defaultColor": "bin" if "bin" in cat_encoded else (
            sorted(cat_encoded)[0] if cat_encoded else None),
        "nTotal": n_total,
        "nEligible": n_eligible,
        "nShown": len(rows),
        "trimmed": trimmed,
        "minLength": min_length,
        "source": str(layout.binarena),
        "fastaHint": fasta_hint(layout),
    }


def fasta_hint(layout: RunLayout) -> str:
    """The contig FASTA whose record names match the BinaRena table's IDs.

    A final/ run renames every contig to MH_<origin>_<bin>_<n>, so only the consolidated
    assembly matches; a preliminary run's IDs are the assembler's own names.
    """
    default = layout.default_contigs()
    if default:
        return str(default)
    partition = layout.bin_fasta_partition()
    if partition:
        return f"the bin FASTAs under {layout.bins_dir} (select them all)"
    return ""


# ------------------------------------------------------------------ split application

def apply_split(layout: RunLayout, split_tsv: Path, contigs: Path, outdir: Path):
    rows = read_tsv(split_tsv)
    if not rows:
        raise SystemExit(f"No rows found in {split_tsv}")
    key = "contig" if "contig" in rows[0] else list(rows[0].keys())[0]
    label_key = None
    for candidate in ("manual_bin", "sub_bin", "label", "bin"):
        if candidate in rows[0]:
            label_key = candidate
            break
    if label_key is None:
        raise SystemExit(
            f"{split_tsv} needs a 'manual_bin' column (found: {list(rows[0])})")

    wanted = {}
    for row in rows:
        contig = (row.get(key) or "").strip()
        label = (row.get(label_key) or "").strip()
        if contig and label:
            wanted[contig] = label
    if not wanted:
        raise SystemExit(f"{split_tsv} contained no usable contig/label pairs")

    sources = []
    if contigs is not None:
        if not Path(contigs).is_file():
            raise SystemExit(f"--contigs not found: {contigs}")
        sources = [Path(contigs)]
    else:
        default = layout.default_contigs()
        if default:
            sources = [default]
        else:
            # No assembly FASTA inside the run directory (contigs were supplied with -i
            # from elsewhere). One rank's bin FASTAs partition every binned contig.
            sources = layout.bin_fasta_partition()
            if sources:
                print(f"No assembly FASTA in the run directory; reading contigs from "
                      f"{len(sources)} bin FASTA(s) under {layout.bins_dir}.",
                      file=sys.stderr)
    if not sources:
        raise SystemExit(
            "Could not find any contig source. Pass --contigs explicitly; expected "
            "final/assembly/consolidated_contigs.fasta, megahit/final.contigs.fa, or "
            "bin FASTAs under bins/<rank>/."
        )

    outdir.mkdir(parents=True, exist_ok=True)
    handles, counts, lengths = {}, {}, {}
    found = 0
    seen = set()
    for source in sources:
        for name, seq in iter_fasta(Path(source)):
            label = wanted.get(name)
            if label is None or name in seen:
                continue
            seen.add(name)
            found += 1
            safe = "".join(c if (c.isalnum() or c in "._-") else "_" for c in label)
            if safe not in handles:
                handles[safe] = open(outdir / f"{safe}.fasta", "w")
                counts[safe] = 0
                lengths[safe] = 0
            handles[safe].write(f">{name}\n")
            for i in range(0, len(seq), 70):
                handles[safe].write(seq[i:i + 70] + "\n")
            counts[safe] += 1
            lengths[safe] += len(seq)
    for handle in handles.values():
        handle.close()

    print(f"Source contigs : {', '.join(str(x) for x in sources[:3])}"
          + (f" (+{len(sources)-3} more)" if len(sources) > 3 else ""),
          file=sys.stderr)
    print(f"Assignments    : {len(wanted)} contig(s) in {split_tsv}", file=sys.stderr)
    print(f"Matched        : {found}", file=sys.stderr)
    missing = len(wanted) - found
    if missing:
        print(
            f"WARNING        : {missing} assigned contig(s) were not found in the FASTA. "
            "If the report was built from final/, split against "
            "final/assembly/consolidated_contigs.fasta, whose records are renamed "
            "MH_<type>_<bin>_<n>.",
            file=sys.stderr,
        )
    for safe in sorted(handles):
        print(f"  {outdir/(safe + '.fasta')}: {counts[safe]} contigs, "
              f"{lengths[safe]} bp", file=sys.stderr)
    return 0


# ------------------------------------------------------------------ HTML template

HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title>
<style>
:root{--bg:#f5f7f8;--card:#fff;--ink:#16222b;--mut:#63757f;--line:#dae3e7;
--acc:#136b78;--accs:#e2f0f2;--warn:#b4453a;--ok:#2f7a44;}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);line-height:1.45;
font-family:Inter,ui-sans-serif,system-ui,-apple-system,"Segoe UI",sans-serif}
.wrap{max-width:1500px;margin:0 auto;padding:22px}
h1{margin:0;font-size:1.6rem;font-weight:660}
h2{margin:0 0 12px;font-size:1rem;font-weight:660}
.sub{color:var(--mut);margin:6px 0 0;font-size:.92rem}
.card{background:var(--card);border:1px solid var(--line);border-radius:11px;padding:16px;margin-bottom:16px}
.tabs{display:flex;gap:6px;margin:18px 0 14px;flex-wrap:wrap}
.tab{padding:8px 15px;border:1px solid var(--line);border-radius:8px;background:var(--card);
cursor:pointer;font-size:.92rem;font-weight:560;color:var(--mut)}
.tab.on{background:var(--acc);border-color:var(--acc);color:#fff}
.pane{display:none}.pane.on{display:block}
table{border-collapse:collapse;width:100%;font-size:.86rem}
th,td{padding:6px 9px;border-bottom:1px solid var(--line);text-align:right;white-space:nowrap}
th:first-child,td:first-child,th.l,td.l{text-align:left}
th{background:var(--accs);font-weight:620;cursor:pointer;position:sticky;top:0}
tbody tr:hover{background:#f0f6f7}
.scroll{max-height:62vh;overflow:auto;border:1px solid var(--line);border-radius:8px}
.stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin-bottom:16px}
.stat{background:var(--card);border:1px solid var(--line);border-radius:11px;padding:13px}
.stat .k{color:var(--mut);font-size:.78rem;text-transform:uppercase;letter-spacing:.05em}
.stat .v{font-size:1.3rem;font-weight:660;margin-top:3px}
label{display:block;font-size:.78rem;color:var(--mut);margin-bottom:3px;font-weight:560}
select,input,button{font:inherit;padding:6px 8px;border:1px solid var(--line);
border-radius:7px;background:var(--card);color:var(--ink);width:100%}
button{cursor:pointer;font-weight:580;width:auto}
button.p{background:var(--acc);border-color:var(--acc);color:#fff}
button.d{background:#fff;border-color:var(--warn);color:var(--warn)}
button:disabled{opacity:.45;cursor:not-allowed}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:11px}
.row{display:flex;gap:9px;align-items:flex-end;flex-wrap:wrap}
.plotwrap{display:grid;grid-template-columns:1fr 310px;gap:15px}
@media(max-width:1100px){.plotwrap{grid-template-columns:1fr}}
canvas{border:1px solid var(--line);border-radius:8px;background:#fff;cursor:crosshair;
display:block;width:100%}
.legend{display:flex;flex-wrap:wrap;gap:5px 12px;margin-top:9px;font-size:.78rem;color:var(--mut)}
.legend i{width:9px;height:9px;border-radius:50%;display:inline-block;margin-right:4px}
.note{color:var(--mut);font-size:.84rem}
.warn{color:var(--warn);font-weight:560}
.pill{display:inline-block;padding:1px 7px;border-radius:9px;font-size:.75rem;font-weight:600}
.pill.ok{background:#e4f2e8;color:var(--ok)}.pill.no{background:#fbeceb;color:var(--warn)}
.sb{display:flex;align-items:center;gap:7px;padding:6px 0;border-bottom:1px solid var(--line);font-size:.85rem}
.sb i{width:11px;height:11px;border-radius:3px;flex:0 0 auto}
.sb span{flex:1;overflow:hidden;text-overflow:ellipsis}
code{background:var(--accs);padding:1px 5px;border-radius:4px;font-size:.85em}
.tip{position:fixed;pointer-events:none;background:#16222b;color:#fff;padding:7px 10px;
border-radius:7px;font-size:.78rem;display:none;z-index:9;max-width:320px;line-height:1.4}
</style></head><body><div class="wrap">
<header><h1>__TITLE__</h1>
<p class="sub">__SUBTITLE__</p></header>
<div class="tabs" id="tabs"></div>
<div id="pane-overview" class="pane on">__OVERVIEW__</div>
<div id="pane-bins" class="pane">__BINS__</div>
<div id="pane-expansion" class="pane">__EXPANSION__</div>
<div id="pane-binarena" class="pane">__BINARENA__</div>
<div id="pane-files" class="pane">__FILES__</div>
</div>
<div class="tip" id="tip"></div>
<script>
var DATA = __DATA__;
/* ---------------- tabs ---------------- */
(function(){
  var names = __TABS__;
  var host = document.getElementById('tabs');
  names.forEach(function(n, i){
    var b = document.createElement('button');
    b.className = 'tab' + (i === 0 ? ' on' : '');
    b.textContent = n.label; b.dataset.id = n.id;
    b.onclick = function(){
      Array.prototype.forEach.call(document.querySelectorAll('.tab'), function(t){
        t.classList.toggle('on', t === b); });
      Array.prototype.forEach.call(document.querySelectorAll('.pane'), function(p){
        p.classList.toggle('on', p.id === 'pane-' + n.id); });
      if (n.id === 'binarena' && window.__draw) window.__draw();
    };
    host.appendChild(b);
  });
})();
/* ---------------- sortable tables ---------------- */
Array.prototype.forEach.call(document.querySelectorAll('table.sortable'), function(tb){
  Array.prototype.forEach.call(tb.querySelectorAll('th'), function(th, idx){
    var dir = 1;
    th.onclick = function(){
      var body = tb.tBodies[0];
      var rows = Array.prototype.slice.call(body.rows);
      rows.sort(function(a, b){
        var x = a.cells[idx].dataset.v, y = b.cells[idx].dataset.v;
        var nx = parseFloat(x), ny = parseFloat(y);
        if (!isNaN(nx) && !isNaN(ny)) return (nx - ny) * dir;
        return String(x || '').localeCompare(String(y || '')) * dir;
      });
      dir = -dir;
      rows.forEach(function(r){ body.appendChild(r); });
    };
  });
});
/* ---------------- bins tab rank filter ---------------- */
(function(){
  var sel = document.getElementById('rank-filter');
  if (!sel) return;
  sel.onchange = function(){
    var v = sel.value;
    Array.prototype.forEach.call(document.querySelectorAll('#bin-table tbody tr'), function(tr){
      tr.style.display = (v === '*' || tr.dataset.rank === v) ? '' : 'none';
    });
  };
})();
/* ---------------- BinaRena scatter ---------------- */
(function(){
  if (!DATA.scatter) return;
  var S = DATA.scatter, N = S.ids.length;
  var cv = document.getElementById('cv'), ctx = cv.getContext('2d');
  var tip = document.getElementById('tip');
  var PAL = ['#2E75B6','#4E9B5B','#8B5CB8','#E28833','#16888C','#C0504D','#B08A2E',
             '#6C7A89','#D06A9C','#3F8F8A','#9A6A3A','#5B7FBF'];
  var SPAL = ['#E8873A','#8B5CB8','#2F7A44','#B4453A','#1F6FEB','#B08A2E','#0F8C8C','#C2478F'];
  var sub = {}, order = [], selected = new Set(), lasso = null, view = [];
  var IDX = {};
  for (var _i = 0; _i < N; _i++) IDX[S.ids[_i]] = _i;   // O(1) id -> row lookup
  var el = function(id){ return document.getElementById(id); };

  function opts(sel, list, cur){
    sel.innerHTML = '';
    list.forEach(function(n){
      var o = document.createElement('option'); o.value = n; o.textContent = n;
      if (n === cur) o.selected = true; sel.appendChild(o);
    });
  }
  // BinaRena-style axis scaling: signed powers from a cube root up to a cube, plus log.
  // Powers are applied as sign(v)*|v|^p so that ordination axes, which are centred on
  // zero and half negative, transform symmetrically instead of collapsing.
  var SCALES = [['cbrt', 'cube root', 1 / 3, 'cube root'],
                ['sqrt', 'square root', 0.5, 'sqrt'],
                ['lin', 'linear', 1, ''],
                ['sq', 'square', 2, 'squared'],
                ['cube', 'cube', 3, 'cubed'],
                ['log', 'log10 (positive only)', 0, 'log10']];
  function scaleOf(id){
    var v = el(id).value;
    for (var k = 0; k < SCALES.length; k++) if (SCALES[k][0] === v) return SCALES[k];
    return SCALES[2];
  }
  function tf(v, sc){
    if (v === null) return null;
    if (sc[0] === 'log') return v > 0 ? Math.log(v) / Math.LN10 : null;
    if (sc[2] === 1) return v;
    return (v < 0 ? -1 : 1) * Math.pow(Math.abs(v), sc[2]);
  }
  function inv(t, sc){
    if (sc[0] === 'log') return Math.pow(10, t);
    if (sc[2] === 1) return t;
    return (t < 0 ? -1 : 1) * Math.pow(Math.abs(t), 1 / sc[2]);
  }
  function scaleOpts(sel){
    sel.innerHTML = '';
    SCALES.forEach(function(sc){
      var o = document.createElement('option');
      o.value = sc[0]; o.textContent = sc[1];
      if (sc[0] === 'lin') o.selected = true;
      sel.appendChild(o);
    });
  }
  scaleOpts(el('xs')); scaleOpts(el('ys'));

  opts(el('x'), S.numericNames, S.defaultX);
  opts(el('y'), S.numericNames, S.defaultY);
  opts(el('c'), S.categoricalNames, S.defaultColor);
  var binLevels = S.categorical.bin ? S.categorical.bin.levels : [];
  opts(el('bf'), ['*'].concat(binLevels), '*');
  el('bf').options[0].textContent = 'all bins';

  function num(name){ return S.numeric[name] || []; }
  function catOf(name){ return S.categorical[name]; }

  // Lazily built lowercase search key per contig: its id plus every categorical value.
  var HAY = null;
  function haystack(){
    if (HAY) return HAY;
    HAY = new Array(N);
    var cats = S.categoricalNames.map(catOf).filter(Boolean);
    for (var i = 0; i < N; i++){
      var parts = [S.ids[i]];
      for (var c = 0; c < cats.length; c++) parts.push(cats[c].levels[cats[c].codes[i]]);
      HAY[i] = parts.join(' ').toLowerCase();
    }
    return HAY;
  }
  // "a, b" keeps rows matching a OR b; a leading '-' on a term excludes instead.
  function parseQuery(raw){
    var inc = [], exc = [];
    (raw || '').toLowerCase().split(',').forEach(function(t){
      t = t.trim();
      if (!t) return;
      if (t.charAt(0) === '-'){ if (t.length > 1) exc.push(t.slice(1)); }
      else inc.push(t);
    });
    return (inc.length || exc.length) ? {inc: inc, exc: exc} : null;
  }
  var TX = null, TY = null, nHiddenByScale = 0;

  function rebuild(){
    var xs = num(el('x').value), ys = num(el('y').value);
    var xsc = scaleOf('xs'), ysc = scaleOf('ys');
    var bf = el('bf').value, bc = catOf('bin');
    var minL = parseFloat(el('ml').value) || 0;
    var lens = num('length');
    var q = parseQuery(el('q').value), hay = q ? haystack() : null;
    TX = new Array(N); TY = new Array(N);
    nHiddenByScale = 0;
    view = [];
    for (var i = 0; i < N; i++){
      if (xs[i] === null || ys[i] === null) continue;
      if (minL && lens.length && lens[i] !== null && lens[i] < minL) continue;
      if (bf !== '*' && bc && bc.levels[bc.codes[i]] !== bf) continue;
      if (q){
        var h = hay[i], ok = q.inc.length === 0, j;
        for (j = 0; j < q.inc.length && !ok; j++) if (h.indexOf(q.inc[j]) >= 0) ok = true;
        if (ok) for (j = 0; j < q.exc.length; j++) if (h.indexOf(q.exc[j]) >= 0){ ok = false; break; }
        if (!ok) continue;
      }
      var tx = tf(xs[i], xsc), ty = tf(ys[i], ysc);
      if (tx === null || ty === null){ nHiddenByScale++; continue; }   // log of <= 0
      TX[i] = tx; TY[i] = ty;
      view.push(i);
    }
    draw();
  }

  function extent(arr, idx){
    var lo = Infinity, hi = -Infinity;
    for (var k = 0; k < idx.length; k++){
      var v = arr[idx[k]];
      if (v === null) continue;
      if (v < lo) lo = v; if (v > hi) hi = v;
    }
    if (lo === Infinity){ lo = 0; hi = 1; }
    if (lo === hi){ lo -= 0.5; hi += 0.5; }
    var pad = (hi - lo) * 0.05;
    return [lo - pad, hi + pad];
  }

  var M = {l: 70, r: 16, t: 14, b: 46}, sx, sy, ex, ey;

  function draw(){
    var w = cv.clientWidth, h = 520, dpr = window.devicePixelRatio || 1;
    cv.width = w * dpr; cv.height = h * dpr;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    var xsc = scaleOf('xs'), ysc = scaleOf('ys');
    ex = extent(TX, view); ey = extent(TY, view);
    var pw = w - M.l - M.r, ph = h - M.t - M.b;
    sx = function(v){ return M.l + (v - ex[0]) / (ex[1] - ex[0]) * pw; };
    sy = function(v){ return M.t + ph - (v - ey[0]) / (ey[1] - ey[0]) * ph; };
    // axes
    ctx.strokeStyle = '#c8d4d9'; ctx.lineWidth = 1; ctx.fillStyle = '#63757f';
    ctx.font = '11px system-ui'; ctx.beginPath();
    ctx.moveTo(M.l, M.t); ctx.lineTo(M.l, M.t + ph); ctx.lineTo(M.l + pw, M.t + ph);
    ctx.stroke();
    for (var g = 0; g <= 4; g++){
      var vx = ex[0] + (ex[1] - ex[0]) * g / 4, vy = ey[0] + (ey[1] - ey[0]) * g / 4;
      ctx.strokeStyle = '#eef3f5'; ctx.beginPath();
      ctx.moveTo(sx(vx), M.t); ctx.lineTo(sx(vx), M.t + ph); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(M.l, sy(vy)); ctx.lineTo(M.l + pw, sy(vy)); ctx.stroke();
      ctx.fillStyle = '#63757f'; ctx.textAlign = 'center';
      ctx.fillText(fmt(inv(vx, xsc)), sx(vx), M.t + ph + 15);
      ctx.textAlign = 'right'; ctx.fillText(fmt(inv(vy, ysc)), M.l - 6, sy(vy) + 3);
    }
    ctx.textAlign = 'center'; ctx.fillStyle = '#16222b'; ctx.font = '600 12px system-ui';
    var xlab = el('x').value + (xsc[3] ? '  [' + xsc[3] + ']' : '');
    var ylab = el('y').value + (ysc[3] ? '  [' + ysc[3] + ']' : '');
    ctx.fillText(xlab, M.l + pw / 2, h - 8);
    ctx.save(); ctx.translate(17, M.t + ph / 2); ctx.rotate(-Math.PI / 2);
    ctx.fillText(ylab, 0, 0); ctx.restore();
    // points
    var cat = catOf(el('c').value), lens = num('length');
    var maxL = 1;
    for (var k = 0; k < view.length; k++){
      var L = lens.length ? lens[view[k]] : 1; if (L > maxL) maxL = L;
    }
    for (var k = 0; k < view.length; k++){
      var i = view[k], id = S.ids[i];
      var L = lens.length && lens[i] !== null ? lens[i] : 1;
      var r = 2.2 + 5.5 * Math.sqrt(L / maxL);
      var col = sub[id] !== undefined ? SPAL[order.indexOf(sub[id]) % SPAL.length]
              : (cat ? PAL[cat.codes[i] % PAL.length] : '#2E75B6');
      ctx.beginPath(); ctx.arc(sx(TX[i]), sy(TY[i]), r, 0, 6.2832);
      ctx.fillStyle = col; ctx.globalAlpha = selected.has(id) ? 1 : 0.72; ctx.fill();
      if (selected.has(id)){
        ctx.globalAlpha = 1; ctx.strokeStyle = '#16222b'; ctx.lineWidth = 1.8; ctx.stroke();
      }
    }
    ctx.globalAlpha = 1;
    if (lasso && lasso.length > 1){
      ctx.beginPath(); ctx.moveTo(lasso[0][0], lasso[0][1]);
      for (var k = 1; k < lasso.length; k++) ctx.lineTo(lasso[k][0], lasso[k][1]);
      ctx.closePath(); ctx.strokeStyle = '#E8873A'; ctx.lineWidth = 1.8;
      ctx.setLineDash([5, 4]); ctx.stroke(); ctx.setLineDash([]);
    }
    legend(cat); status();
  }
  window.__draw = draw;

  function fmt(v){
    var a = Math.abs(v);
    if (a >= 1e6) return (v / 1e6).toFixed(1) + 'M';
    if (a >= 1e3) return (v / 1e3).toFixed(1) + 'k';
    if (a >= 10) return v.toFixed(0);
    return v.toFixed(2);
  }

  function legend(cat){
    var host = el('lg'); host.innerHTML = '';
    if (!cat) return;
    var counts = {};
    for (var k = 0; k < view.length; k++){
      var c = cat.codes[view[k]]; counts[c] = (counts[c] || 0) + 1;
    }
    Object.keys(counts).sort(function(a, b){ return counts[b] - counts[a]; })
      .slice(0, 14).forEach(function(c){
        var s = document.createElement('span');
        s.innerHTML = '<i style="background:' + PAL[c % PAL.length] + '"></i>' +
          cat.levels[c] + ' (' + counts[c] + ')';
        host.appendChild(s);
      });
  }

  function status(){
    var lens = num('length'), gc = num('GC'), cov = num('coverage');
    var bp = 0, n = 0, gsum = 0, gn = 0, csum = 0, cn = 0;
    selected.forEach(function(id){
      var i = IDX[id]; if (i === undefined) return; n++;
      if (lens.length && lens[i] !== null) bp += lens[i];
      if (gc.length && gc[i] !== null){ gsum += gc[i]; gn++; }
      if (cov.length && cov[i] !== null){ csum += cov[i]; cn++; }
    });
    el('shown').textContent = view.length + ' of ' + S.nShown + ' plotted contigs' +
      (nHiddenByScale ? ' (' + nHiddenByScale + ' hidden: not positive on a log axis)' : '');
    el('selinfo').innerHTML = n === 0 ? '<span class="note">Drag on the plot to lasso contigs.</span>'
      : '<b>' + n + '</b> selected &middot; ' + (bp / 1e6).toFixed(3) + ' Mb' +
        (gn ? ' &middot; GC ' + (gsum / gn).toFixed(1) + '%' : '') +
        (cn ? ' &middot; cov ' + (csum / cn).toFixed(1) + '&times;' : '');
    el('assign').disabled = n === 0;
    var host = el('sblist'); host.innerHTML = '';
    if (!order.length){ host.innerHTML = '<p class="note">No sub-bins yet.</p>'; }
    var agg = {};
    Object.keys(sub).forEach(function(id){
      var label = sub[id];
      if (!agg[label]) agg[label] = [0, 0];
      agg[label][0] += 1;
      var i = IDX[id];
      if (i !== undefined && lens.length && lens[i] !== null) agg[label][1] += lens[i];
    });
    order.forEach(function(label, li){
      var cnt = (agg[label] || [0, 0])[0], tot = (agg[label] || [0, 0])[1];
      var d = document.createElement('div'); d.className = 'sb';
      d.innerHTML = '<i style="background:' + SPAL[li % SPAL.length] + '"></i>' +
        '<span>' + label + '</span><span style="flex:0 0 auto;color:#63757f">' + cnt +
        ' &middot; ' + (tot / 1e6).toFixed(2) + ' Mb</span>';
      var x = document.createElement('button'); x.className = 'd'; x.textContent = '\u00d7';
      x.style.padding = '1px 7px';
      x.onclick = function(){
        Object.keys(sub).forEach(function(id){ if (sub[id] === label) delete sub[id]; });
        order.splice(order.indexOf(label), 1); draw();
      };
      d.appendChild(x); host.appendChild(d);
    });
    el('exp').disabled = !order.length;
    var faSel = el('fa') && el('fa').files && el('fa').files.length;
    if (el('expfa')) el('expfa').disabled = !(order.length && faSel);
  }

  function inPoly(px, py, poly){
    var hit = false;
    for (var i = 0, j = poly.length - 1; i < poly.length; j = i++){
      var xi = poly[i][0], yi = poly[i][1], xj = poly[j][0], yj = poly[j][1];
      if (((yi > py) !== (yj > py)) && (px < (xj - xi) * (py - yi) / (yj - yi) + xi)) hit = !hit;
    }
    return hit;
  }

  var drag = false;
  function pos(e){
    var r = cv.getBoundingClientRect();
    return [e.clientX - r.left, e.clientY - r.top];
  }
  cv.addEventListener('mousedown', function(e){
    drag = true; lasso = [pos(e)]; });
  cv.addEventListener('mousemove', function(e){
    if (drag){ lasso.push(pos(e)); draw(); return; }
    var p = pos(e), best = -1, bd = 81;
    for (var k = 0; k < view.length; k++){
      var i = view[k], dx = sx(TX[i]) - p[0], dy = sy(TY[i]) - p[1], d = dx * dx + dy * dy;
      if (d < bd){ bd = d; best = i; }
    }
    if (best < 0){ tip.style.display = 'none'; return; }
    var lines = ['<b>' + S.ids[best] + '</b>'];
    ['length','GC','coverage','coding_density'].forEach(function(f){
      var a = num(f); if (a.length && a[best] !== null) lines.push(f + ': ' + a[best]);
    });
    S.categoricalNames.forEach(function(f){
      var c = catOf(f); if (c) lines.push(f + ': ' + c.levels[c.codes[best]]);
    });
    tip.innerHTML = lines.join('<br>');
    tip.style.display = 'block';
    tip.style.left = (e.clientX + 14) + 'px'; tip.style.top = (e.clientY + 12) + 'px';
  });
  cv.addEventListener('mouseleave', function(){ tip.style.display = 'none'; });
  window.addEventListener('mouseup', function(e){
    if (!drag) return;
    drag = false;
    if (lasso && lasso.length > 2){
      if (!e.shiftKey) selected.clear();
      for (var k = 0; k < view.length; k++){
        var i = view[k];
        if (inPoly(sx(TX[i]), sy(TY[i]), lasso)) selected.add(S.ids[i]);
      }
    }
    lasso = null; draw();
  });

  el('assign').onclick = function(){
    var label = (el('lbl').value || '').trim();
    if (!label){ alert('Give the sub-bin a name first.'); return; }
    if (order.indexOf(label) < 0) order.push(label);
    selected.forEach(function(id){ sub[id] = label; });
    selected.clear(); el('lbl').value = ''; draw();
  };
  el('clr').onclick = function(){ selected.clear(); draw(); };
  el('exp').onclick = function(){
    var bc = catOf('bin');
    var out = ['contig\tsource_bin\tmanual_bin'];
    Object.keys(sub).forEach(function(id){
      var i = IDX[id];
      var src = (bc && i !== undefined) ? bc.levels[bc.codes[i]] : 'NA';
      out.push(id + '\t' + src + '\t' + sub[id]);
    });
    var blob = new Blob([out.join('\n') + '\n'], {type: 'text/tab-separated-values'});
    var a = document.createElement('a');
    a.href = URL.createObjectURL(blob); a.download = 'manual_bins.tsv';
    document.body.appendChild(a); a.click(); a.remove();
  };
  // ---------------------------------------------------------------- FASTA export
  // Sequences are never embedded in this file -- a whole assembly would be far too
  // large -- so the user points us at the FASTA on disk and we read it locally with
  // the File API. Nothing leaves the machine.
  var CRC = (function(){
    var t = new Int32Array(256);
    for (var n = 0; n < 256; n++){
      var c = n;
      for (var k = 0; k < 8; k++) c = (c & 1) ? (0xEDB88320 ^ (c >>> 1)) : (c >>> 1);
      t[n] = c;
    }
    return function(bytes){
      var c = -1;
      for (var i = 0; i < bytes.length; i++) c = t[(c ^ bytes[i]) & 0xFF] ^ (c >>> 8);
      return (c ^ -1) >>> 0;
    };
  })();

  // Minimal STORE-only ZIP writer, so several FASTA files can come down as one archive
  // without pulling in a compression library.
  function makeZip(files){
    var enc = new TextEncoder(), parts = [], central = [], offset = 0;
    function u16(v){ return [v & 0xFF, (v >>> 8) & 0xFF]; }
    function u32(v){ return [v & 0xFF, (v >>> 8) & 0xFF, (v >>> 16) & 0xFF, (v >>> 24) & 0xFF]; }
    files.forEach(function(f){
      var name = enc.encode(f.name), body = enc.encode(f.text), crc = CRC(body);
      var local = [].concat(u32(0x04034B50), u16(20), u16(0), u16(0), u16(0), u16(0),
                            u32(crc), u32(body.length), u32(body.length),
                            u16(name.length), u16(0));
      parts.push(new Uint8Array(local), name, body);
      central.push([].concat(u32(0x02014B50), u16(20), u16(20), u16(0), u16(0), u16(0),
                             u16(0), u32(crc), u32(body.length), u32(body.length),
                             u16(name.length), u16(0), u16(0), u16(0), u16(0), u32(0),
                             u32(offset)));
      central[central.length - 1].nameBytes = name;
      offset += local.length + name.length + body.length;
    });
    var cdStart = offset, cdSize = 0;
    central.forEach(function(c){
      parts.push(new Uint8Array(c), c.nameBytes);
      cdSize += c.length + c.nameBytes.length;
    });
    parts.push(new Uint8Array([].concat(u32(0x06054B50), u16(0), u16(0),
      u16(files.length), u16(files.length), u32(cdSize), u32(cdStart), u16(0))));
    return new Blob(parts, {type: 'application/zip'});
  }

  function download(blob, name){
    var a = document.createElement('a');
    a.href = URL.createObjectURL(blob); a.download = name;
    document.body.appendChild(a); a.click();
    setTimeout(function(){ URL.revokeObjectURL(a.href); a.remove(); }, 1000);
  }

  // Streams a FASTA in chunks, keeping only the wanted records, so a multi-GB assembly
  // never has to sit in memory at once. Record names are the first whitespace-delimited
  // token, matching how MEGAHIT and Unicycler decorate their headers.
  function readFastaSubset(file, wanted, onSeq, onDone, onErr){
    var CHUNK = 8 << 20, pos = 0, tail = '', keep = false, name = null, buf = [];
    function flush(){
      if (keep && name !== null) onSeq(name, buf.join(''));
      buf = []; keep = false; name = null;
    }
    function handleLine(line){
      if (line.charAt(0) === '>'){
        flush();
        name = line.slice(1).split(/[\s]/)[0];
        keep = wanted.has(name);
      } else if (keep){
        buf.push(line.trim());
      }
    }
    function step(){
      if (pos >= file.size){
        tail.split('\n').forEach(function(l){ if (l !== '') handleLine(l); });
        flush(); onDone(); return;
      }
      var slice = file.slice(pos, Math.min(pos + CHUNK, file.size));
      pos += CHUNK;
      var fr = new FileReader();
      fr.onerror = function(){ onErr('Could not read ' + file.name); };
      fr.onload = function(){
        var text = tail + fr.result, lines = text.split('\n');
        tail = lines.pop();                       // may be a partial line
        for (var i = 0; i < lines.length; i++){
          var l = lines[i];
          if (l !== '' && l !== '\r') handleLine(l.charAt(l.length - 1) === '\r' ? l.slice(0, -1) : l);
        }
        el('fastat').textContent = 'Reading ' + file.name + '... ' +
          Math.min(100, Math.round(100 * pos / file.size)) + '%';
        setTimeout(step, 0);                      // yield so the UI stays responsive
      };
      fr.readAsText(slice);
    }
    step();
  }

  function wrap(seq){
    var out = [];
    for (var i = 0; i < seq.length; i += 70) out.push(seq.slice(i, i + 70));
    return out.join('\n');
  }

  el('fa').onchange = function(){
    el('expfa').disabled = !(this.files && this.files.length && order.length);
    el('fastat').textContent = this.files && this.files.length
      ? this.files.length + ' file(s) selected.' : '';
  };

  el('expfa').onclick = function(){
    var files = Array.prototype.slice.call(el('fa').files || []);
    if (!files.length){ alert('Choose the contig FASTA first.'); return; }
    if (files.some(function(f){ return /\.gz$/i.test(f.name); })){
      el('fastat').textContent = 'Gzipped FASTA cannot be read in the browser. ' +
        'Decompress it, or use --apply-split on the server.';
      return;
    }
    var wanted = new Set(Object.keys(sub));
    var seqs = {}, found = 0, bytes = 0;
    el('expfa').disabled = true;
    var fi = 0;
    function nextFile(){
      if (fi >= files.length){ finish(); return; }
      var f = files[fi++];
      readFastaSubset(f, wanted, function(name, seq){
        if (seqs[name] === undefined){ seqs[name] = seq; found++; bytes += seq.length; }
      }, nextFile, function(msg){
        el('fastat').textContent = msg; el('expfa').disabled = false;
      });
    }
    function finish(){
      if (!found){
        el('fastat').innerHTML = '<b>No matching records.</b> None of the ' + wanted.size +
          ' assigned contig names were found in the selected file(s). This usually means ' +
          'the wrong FASTA was chosen &mdash; a <code>final/</code> run needs ' +
          '<code>consolidated_contigs.fasta</code>, whose records are renamed ' +
          '<code>MH_...</code>.';
        el('expfa').disabled = false; return;
      }
      var byBin = {};
      Object.keys(sub).forEach(function(id){
        if (seqs[id] === undefined) return;
        (byBin[sub[id]] = byBin[sub[id]] || []).push(id);
      });
      var out = [];
      order.forEach(function(label){
        var ids = byBin[label];
        if (!ids || !ids.length) return;
        var safe = label.replace(/[^A-Za-z0-9._-]/g, '_');
        var text = ids.map(function(id){ return '>' + id + '\n' + wrap(seqs[id]); }).join('\n') + '\n';
        out.push({name: safe + '.fasta', text: text});
      });
      if (!out.length){ el('fastat').textContent = 'Nothing to write.';
                        el('expfa').disabled = false; return; }
      if (out.length === 1) download(new Blob([out[0].text], {type: 'text/plain'}), out[0].name);
      else download(makeZip(out), 'manual_bins_fasta.zip');
      var missing = wanted.size - found;
      el('fastat').innerHTML = 'Wrote ' + out.length + ' FASTA file(s), ' + found +
        ' contigs, ' + (bytes / 1e6).toFixed(2) + ' Mb.' +
        (missing > 0 ? ' <b>' + missing + ' assigned contig(s) were not found</b> in the ' +
         'selected file(s) and were skipped.' : '');
      el('expfa').disabled = false;
    }
    nextFile();
  };

  ['x','y','c','bf','ml','xs','ys'].forEach(function(id){ el(id).onchange = rebuild; });
  var qt = null;
  el('q').addEventListener('input', function(){
    clearTimeout(qt); qt = setTimeout(rebuild, 160);   // debounce so typing stays smooth
  });
  window.addEventListener('resize', function(){ if (view.length) draw(); });
  rebuild();
})();
</script></body></html>
"""


# ------------------------------------------------------------------ rendering

def esc(value):
    return html.escape("NA" if value is None else str(value))


def cell(value, fmt=None, klass=""):
    if value is None or value == "":
        return f'<td class="{klass}" data-v="">NA</td>'
    text = fmt(value) if fmt else esc(value)
    return f'<td class="{klass}" data-v="{esc(value)}">{text}</td>'


def num_fmt(value):
    return f"{value:,.0f}" if abs(value) >= 1000 else f"{value:g}"


def render_overview(layout: RunLayout, bin_rows, ranks, scatter, provenance_counts):
    total_bp = max((r["total_bp"] for r in ranks), default=0)
    n_bins = sum(r["num_classified_bins"] for r in ranks)
    checkm = [r for r in bin_rows if r.get("completeness_percent") is not None]
    good = [r for r in checkm
            if (r["completeness_percent"] or 0) >= 70
            and (r.get("contamination_percent") or 0) <= 10]
    stats = [
        ("Bin set shown", layout.stage),
        ("Ranks", str(len(ranks))),
        ("Classified bins", str(n_bins)),
        ("Assembly size", human_bp(total_bp)),
    ]
    if checkm:
        stats.append((">=70% compl, <=10% contam", f"{len(good)} / {len(checkm)}"))
    if scatter:
        stats.append(("Contigs plotted", f"{scatter['nShown']:,}"))

    out = ['<div class="stats">']
    for key, value in stats:
        out.append(f'<div class="stat"><div class="k">{esc(key)}</div>'
                   f'<div class="v">{esc(value)}</div></div>')
    out.append("</div>")

    out.append('<div class="card"><h2>Per-rank totals</h2><div class="scroll">'
               '<table class="sortable"><thead><tr>'
               '<th class="l">rank</th><th>bins</th><th>classified bins</th>'
               '<th>contigs</th><th>total bp</th><th>classified bp</th>'
               '<th>classified %</th><th class="l">largest bin</th>'
               '<th>largest bin bp</th></tr></thead><tbody>')
    for row in ranks:
        out.append(
            "<tr>"
            + cell(row["rank"], klass="l")
            + cell(row["num_bins"], num_fmt) + cell(row["num_classified_bins"], num_fmt)
            + cell(row["total_contigs"], num_fmt)
            + cell(row["total_bp"], lambda v: human_bp(v))
            + cell(row["classified_bp"], lambda v: human_bp(v))
            + cell(row["classified_fraction"] * 100, lambda v: f"{v:.1f}%")
            + cell(row["largest_bin"], klass="l")
            + cell(row["largest_bin_bp"], lambda v: human_bp(v))
            + "</tr>")
    out.append("</tbody></table></div>")
    out.append('<p class="note" style="margin-top:10px">Bins at different ranks are '
               'nested: the same contig appears in its family, genus and species bin. '
               'Totals are therefore per rank, not additive across ranks.</p></div>')

    if provenance_counts:
        out.append('<div class="card"><h2>Final-assembly provenance</h2>'
                   '<div class="scroll"><table class="sortable"><thead><tr>'
                   '<th class="l">source bin</th><th class="l">source</th>'
                   '<th>contigs</th></tr></thead><tbody>')
        for (bin_name, source_type), count in sorted(
                provenance_counts.items(), key=lambda kv: -kv[1]):
            pill = f'<span class="pill">{esc(source_type)}</span>' 
            out.append("<tr>" + cell(bin_name, klass="l")
                       + f'<td class="l" data-v="{esc(source_type)}">{pill}</td>'
                       + cell(count, num_fmt) + "</tr>")
        out.append("</tbody></table></div></div>")
    return "\n".join(out)


def render_bins(bin_rows, ranks):
    has_checkm = any(r.get("completeness_percent") is not None for r in bin_rows)
    has_reas = any(r.get("reassembly_outcome") for r in bin_rows)
    head = ['<th class="l">rank</th>', '<th class="l">bin</th>', "<th>contigs</th>",
            "<th>total bp</th>", "<th>N50</th>", "<th>largest</th>", "<th>GC %</th>",
            "<th>contigs/Mb</th>"]
    if has_checkm:
        head += ["<th>compl %</th>", "<th>contam %</th>", "<th>strain het %</th>"]
    if has_reas:
        head += ['<th class="l">reassembly</th>', "<th>recovered</th>"]
    out = ['<div class="card"><div class="row" style="margin-bottom:12px">'
           '<div style="max-width:220px"><label for="rank-filter">Rank</label>'
           '<select id="rank-filter"><option value="*">all ranks</option>']
    for row in ranks:
        out.append(f'<option value="{esc(row["rank"])}">{esc(row["rank"])}</option>')
    out.append("</select></div></div>")
    out.append('<div class="scroll"><table class="sortable" id="bin-table"><thead><tr>'
               + "".join(head) + "</tr></thead><tbody>")
    for row in bin_rows:
        cells = ["<tr" + f' data-rank="{esc(row["rank"])}">',
                 cell(row["rank"], klass="l"), cell(row["bin"], klass="l"),
                 cell(row.get("num_contigs"), num_fmt),
                 cell(row.get("total_length_bp"), lambda v: human_bp(v)),
                 cell(row.get("N50"), num_fmt),
                 cell(row.get("largest_contig_bp"), num_fmt),
                 cell(row.get("GC_percent"), lambda v: f"{v:.2f}"),
                 cell(row.get("contigs_per_Mbp"), lambda v: f"{v:.1f}")]
        if has_checkm:
            cells += [cell(row.get("completeness_percent"), lambda v: f"{v:.2f}"),
                      cell(row.get("contamination_percent"), lambda v: f"{v:.2f}"),
                      cell(row.get("strain_heterogeneity_percent"), lambda v: f"{v:.2f}")]
        if has_reas:
            outcome = row.get("reassembly_outcome")
            if outcome == "accepted":
                pill = '<span class="pill ok">accepted</span>'
            elif outcome == "linked":
                pill = '<span class="pill ok">linked</span>'
            elif outcome:
                pill = f'<span class="pill no">{esc(outcome)}</span>'
            else:
                pill = "&ndash;"
            cells.append(f'<td class="l" data-v="{esc(outcome or "")}">{pill}</td>')
            cells.append(cell(row.get("recovered_fraction"), lambda v: f"{100*v:.0f}%"))
        cells.append("</tr>")
        out.append("".join(cells))
    out.append("</tbody></table></div>")
    if has_reas:
        out.append('<p class="note" style="margin-top:10px">A <b>rejected</b> reassembly '
                   'retained less than <code>--reassemble-min-recovered-fraction</code> '
                   'of the bin\'s original length, so the original contigs were kept in '
                   'the consolidated assembly.</p>')
    out.append("</div>")
    return "\n".join(out)


def render_binarena(scatter):
    if not scatter:
        return ('<div class="card"><h2>No BinaRena table found</h2>'
                '<p class="note">This run has no <code>binarena/binarena_input.tsv</code>. '
                'Re-run without <code>--skip-binarena</code> to generate it.</p></div>')
    warn = ""
    if scatter["trimmed"]:
        warn = (f'<p class="note warn">Showing the {scatter["nShown"]:,} longest of '
                f'{scatter["nEligible"]:,} eligible contigs. Raise '
                f'<code>--max-points</code> or <code>--min-length</code> to change this.</p>')
    return f"""
<div class="card">
  <div class="grid" style="margin-bottom:11px">
    <div><label for="x">X axis</label><select id="x"></select></div>
    <div><label for="xs">X scale</label><select id="xs"></select></div>
    <div><label for="y">Y axis</label><select id="y"></select></div>
    <div><label for="ys">Y scale</label><select id="ys"></select></div>
  </div>
  <div class="grid" style="margin-bottom:13px">
    <div><label for="c">Colour by</label><select id="c"></select></div>
    <div><label for="bf">Restrict to bin</label><select id="bf"></select></div>
    <div><label for="ml">Min contig length</label>
      <input id="ml" type="number" value="{scatter['minLength']}" step="500" min="0"></div>
    <div style="grid-column:span 2"><label for="q">Search</label>
      <input id="q" type="search" placeholder="contig name or taxon; comma-separated = OR; -term excludes"></div>
  </div>
  <div class="plotwrap">
    <div>
      <canvas id="cv" height="520"></canvas>
      <div class="legend" id="lg"></div>
      <p class="note" style="margin-top:8px">
        Drag to lasso. Hold <b>Shift</b> while dragging to add to the selection.
        Hover a point for its values. <span id="shown"></span>
      </p>
      {warn}
    </div>
    <div>
      <h2>Split a bin</h2>
      <p id="selinfo" class="note" style="margin-bottom:11px"></p>
      <label for="lbl">Sub-bin name</label>
      <input id="lbl" placeholder="e.g. Symbiopectobacterium_cluster_01">
      <div class="row" style="margin:10px 0 14px">
        <button class="p" id="assign" disabled>Assign selection</button>
        <button id="clr">Clear</button>
      </div>
      <h2>Sub-bins</h2>
      <div id="sblist"></div>
      <div class="row" style="margin-top:12px">
        <button class="p" id="exp" disabled>Download manual_bins.tsv</button>
      </div>
      <h2 style="margin-top:18px">Download FASTA</h2>
      <p class="note">Pick the contig FASTA so the sequences can be written here in the
         browser. Nothing is uploaded.</p>
      <p class="note" style="margin-top:6px">Expected:<br><code>{esc(scatter['fastaHint']) or 'no assembly FASTA found in this run'}</code></p>
      <input id="fa" type="file" multiple accept=".fasta,.fa,.fna,.fasta.gz,.fa.gz,.txt"
             style="margin-top:8px">
      <div class="row" style="margin-top:10px">
        <button class="p" id="expfa" disabled>Write FASTA per sub-bin</button>
      </div>
      <p class="note" id="fastat" style="margin-top:9px"></p>
      <p class="note" style="margin-top:12px">Or do it from the TSV on the server:</p>
      <p><code>metahopper_report.py -i &lt;outdir&gt; --apply-split manual_bins.tsv</code></p>
    </div>
  </div>
  <p class="note" style="margin-top:12px">Source:
     <code>{esc(scatter['source'])}</code> &middot; {scatter['nTotal']:,} rows,
     {scatter['nEligible']:,} at or above {scatter['minLength']:,} bp with ordination
     coordinates.</p>
</div>"""


def render_files(layout: RunLayout):
    checks = [
        ("BinaRena table", layout.binarena),
        ("Contig classification", layout.classification),
        ("Bin refinement (step 9)", layout.refinement),
        ("Final-assembly provenance", layout.provenance),
        ("Consolidated contigs", layout.consolidated),
        ("MEGAHIT contigs", layout.megahit),
    ]
    out = ['<div class="card"><h2>Files discovered</h2>'
           '<div class="scroll"><table><thead><tr><th class="l">what</th>'
           '<th class="l">status</th><th class="l">path</th></tr></thead><tbody>']
    for label, path in checks:
        pill = ('<span class="pill ok">found</span>' if path
                else '<span class="pill no">absent</span>')
        out.append(f'<tr><td class="l">{esc(label)}</td><td class="l">{pill}</td>'
                   f'<td class="l">{esc(path) if path else "&ndash;"}</td></tr>')
    for rank in layout.ranks():
        rank_dir = layout.bins_dir / rank
        n = len(list(rank_dir.glob("*.fasta")))
        out.append(f'<tr><td class="l">bins/{esc(rank)}</td>'
                   f'<td class="l"><span class="pill ok">{n} FASTA</span></td>'
                   f'<td class="l">{esc(rank_dir)}</td></tr>')
    out.append("</tbody></table></div>")
    out.append(f'<p class="note" style="margin-top:11px">Run directory: '
               f'<code>{esc(layout.root)}</code>. Seed-and-extension '
               f'{"ran" if layout.has_reassembly else "did not run"}; '
               f'<code>final/</code> {"present" if layout.has_final else "absent"}.</p>')
    out.append("</div>")
    return "\n".join(out)


def render_expansion(layout):
    """Compare source-bin candidates separately from final taxonomically re-binned output."""
    root = layout.root
    tables = sorted((root / "bins").glob("*/reassembly_summary.tsv"))
    if not tables:
        return '<div class="card"><h2>Expansion</h2><p>No expansion comparison is available. A bin-only run does not perform expansion.</p></div>'
    out = ['<div class="card"><h2>Expansion: cost and outcome</h2>',
           '<p>Candidate metrics describe the proposed bin before final reclassification/refinement. '
           'Final metrics describe the final bin of the same name; renamed or split bins need manual review. '
           'A larger bin or higher N50 alone does not establish better genome recovery. '
           'Link mode adds existing contigs and does not join sequence.</p>']
    headers = ['rank', 'bin', 'outcome', 'bin minutes', 'contigs before', 'candidate contigs',
               'final contigs', 'bp before', 'candidate bp', 'final bp', 'N50 before',
               'candidate N50', 'final N50', 'original sequence retained %',
               'completeness before %', 'candidate completeness %', 'final completeness %',
               'contamination before %', 'candidate contamination %', 'final contamination %',
               'marker check', 'rounds', 'seed templates', 'accepted templates',
               'recruitment minutes', 'assembly minutes', 'polish minutes', 'polish outcome', 'reason / stop']
    out.append('<div class="scroll"><table class="sortable"><thead><tr>' +
               ''.join('<th>' + esc(x) + '</th>' for x in headers) + '</tr></thead><tbody>')
    timing_notes = []
    for path in tables:
        rank = path.parent.name
        before = {r.get('bin'): r for r in read_tsv(path.parent / 'summary.tsv')}
        after = {r.get('bin'): r for r in read_tsv(root / 'final' / 'bins' / rank / 'summary.tsv')}
        timing = root / 'reassembly' / rank / 'rank_timing.json'
        if timing.is_file():
            try:
                t = json.loads(timing.read_text())
                timing_notes.append(f"{rank}: {float(t.get('elapsed_seconds', 0)) / 60:.2f} minutes including shared mapping/extraction and candidate quality assessment")
            except (ValueError, TypeError):
                pass
        for row in read_tsv(path):
            b, a = before.get(row.get('bin'), {}), after.get(row.get('bin'), {})
            def n(key):
                return to_float(row.get(key))
            def minutes(key):
                value = n(key)
                return None if value is None else value / 60
            def baseline(key):
                return to_float(row.get(key + '_before')) if to_float(row.get(key + '_before')) is not None else to_float(b.get(key))
            retained = n('original_sequence_retained')
            values = [rank, row.get('bin'), row.get('outcome'), minutes('elapsed_seconds'),
                      n('contigs_before'), n('contigs_after'), to_float(a.get('num_contigs')),
                      n('total_length_before_bp'), n('total_length_after_bp'), to_float(a.get('total_length_bp')),
                      n('N50_before'), n('N50_after'), to_float(a.get('N50')),
                      None if retained is None else retained * 100,
                      baseline('completeness_percent'), n('completeness_percent_candidate'), to_float(a.get('completeness_percent')),
                      baseline('contamination_percent'), n('contamination_percent_candidate'), to_float(a.get('contamination_percent')),
                      row.get('marker_check', 'not evaluated'), n('extension_rounds'), n('seed_templates'), n('accepted_templates'),
                      minutes('recruitment_seconds'), minutes('assembly_seconds'), minutes('polish_seconds'), row.get('polish_outcome'),
                      row.get('reason') or row.get('recruitment_stop')]
            out.append('<tr>' + ''.join(cell(v, (lambda x: f'{x:,.2f}') if isinstance(v, (int, float)) else None) for v in values) + '</tr>')
    out.append('</tbody></table></div>')
    out.append('<p class="note">Per-bin time excludes shared mapping, extraction, and batch CheckM. Missing marker evidence is not a quality pass. Candidate marker scores precede optional polishing.</p>')
    for note in timing_notes:
        out.append('<p>' + esc(note) + '</p>')
    out.append('</div>')
    return '\n'.join(out)


def build_report(layout: RunLayout, args) -> str:
    bin_rows = collect_bins(layout, compute_missing=not args.no_fasta_stats)
    ranks = rank_totals(bin_rows)
    scatter = build_scatter(layout, args.min_length, args.max_points,
                            not args.include_uncoordinated)
    provenance_counts = {}
    if layout.provenance:
        for row in read_tsv(layout.provenance):
            key = (row.get("source_bin", "NA"), row.get("source_type", "NA"))
            provenance_counts[key] = provenance_counts.get(key, 0) + 1

    tabs = [{"id": "overview", "label": "Overview"}, {"id": "bins", "label": "Bins"},
            {"id": "expansion", "label": "Expansion"},
            {"id": "binarena", "label": "BinaRena"}, {"id": "files", "label": "Files"}]
    title = args.title or f"MetaHopper report \u2014 {layout.root.name}"
    subtitle = (f"{layout.stage} bin set &middot; "
                f"{len(ranks)} rank(s) &middot; built "
                f"{datetime.now().strftime('%Y-%m-%d %H:%M')}")

    payload = json.dumps({"scatter": scatter}, separators=(",", ":"),
                         allow_nan=False, default=str)
    payload = payload.replace("</", "<\\/")

    return (HTML
            .replace("__TITLE__", esc(title))
            .replace("__SUBTITLE__", subtitle)
            .replace("__TABS__", json.dumps(tabs))
            .replace("__OVERVIEW__", render_overview(layout, bin_rows, ranks, scatter,
                                                     provenance_counts))
            .replace("__BINS__", render_bins(bin_rows, ranks))
            .replace("__EXPANSION__", render_expansion(layout))
            .replace("__BINARENA__", render_binarena(scatter))
            .replace("__FILES__", render_files(layout))
            .replace("__DATA__", payload))


# ------------------------------------------------------------------ CLI

def parse_args():
    p = argparse.ArgumentParser(
        description="Interactive HTML report for one MetaHopper output directory.")
    p.add_argument("-i", "--input", type=Path, default=Path("."),
                   help="MetaHopper run directory (default: .)")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="Output HTML file (default: <input>/metahopper_report.html)")
    p.add_argument("--bin-set", choices=("auto", "preliminary", "final"), default="auto",
                   help="Which bin set to report: auto prefers final/ when present.")
    p.add_argument("--title", default=None, help="Report title")
    p.add_argument("--min-length", type=int, default=1000,
                   help="Minimum contig length shown in the BinaRena tab (default 1000)")
    p.add_argument("--max-points", type=int, default=25000,
                   help="Maximum contigs embedded in the scatter, longest first "
                        "(default 25000). Keeps the HTML small enough to open.")
    p.add_argument("--include-uncoordinated", action="store_true",
                   help="Also plot contigs with no k-mer ordination coordinates "
                        "(they can still be placed on GC/coverage axes).")
    p.add_argument("--no-fasta-stats", action="store_true",
                   help="Do not fall back to reading bin FASTAs when a rank has no "
                        "summary.tsv (faster on very large runs).")
    p.add_argument("--apply-split", type=Path, default=None, metavar="TSV",
                   help="Instead of building a report, read a manual_bins.tsv exported "
                        "from the BinaRena tab and write one FASTA per sub-bin.")
    p.add_argument("--split-outdir", type=Path, default=None,
                   help="Destination for --apply-split (default: <input>/manual_bins)")
    p.add_argument("--contigs", type=Path, default=None,
                   help="Contig FASTA used by --apply-split. Defaults to "
                        "final/assembly/consolidated_contigs.fasta, else "
                        "megahit/final.contigs.fa.")
    return p.parse_args()


def main():
    args = parse_args()
    root = args.input
    if not root.is_dir():
        raise SystemExit(f"Not a directory: {root}")
    layout = RunLayout(root, args.bin_set)

    if args.apply_split:
        outdir = args.split_outdir or (layout.root / "manual_bins")
        return apply_split(layout, args.apply_split, args.contigs, outdir)

    if not layout.bins_dir.is_dir():
        raise SystemExit(
            f"No bins directory under {layout.base}. Point -i at a MetaHopper run "
            "directory (the one containing bins/, classification/, prodigal/)."
        )
    out_path = args.output or (layout.root / "metahopper_report.html")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(build_report(layout, args))
    size_kb = out_path.stat().st_size / 1024
    print(f"Wrote {out_path} ({size_kb:,.0f} kB)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
