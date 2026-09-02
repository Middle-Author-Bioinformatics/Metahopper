#!/usr/bin/env python3
"""Build a self-contained interactive HTML report from MetaHopper bin metrics.

Expected directory layout
-------------------------

Parent directory containing multiple samples::

    metahop_SAMPLE_A/summary/bin_metrics.tsv
    metahop_SAMPLE_B/summary/bin_metrics.tsv

or one sample directory::

    metahop_SAMPLE_A/summary/bin_metrics.tsv

The report reads only ``summary/bin_metrics.tsv``. It deliberately ignores
``bins/*/checkm_out``, ``bins/*/quast_out``, and ``bins/*/summary.tsv``.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import sys
from datetime import datetime, timezone
from pathlib import Path


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__REPORT_TITLE__</title>
<style>
  :root {
    color-scheme: light dark;
    --bg: #f4f7f9;
    --surface: #ffffff;
    --surface-soft: #edf3f5;
    --text: #17222a;
    --muted: #60717c;
    --border: #d5e0e5;
    --accent: #176b78;
    --accent-soft: #dceef0;
    --track: #e2e9ec;
    --shadow: 0 8px 24px rgba(28, 54, 66, 0.08);
  }
  @media (prefers-color-scheme: dark) {
    :root {
      --bg: #10171b;
      --surface: #172126;
      --surface-soft: #1e2c32;
      --text: #edf5f6;
      --muted: #a7b8bf;
      --border: #31434b;
      --accent: #65c1ca;
      --accent-soft: #203d43;
      --track: #29383f;
      --shadow: none;
    }
  }
  * { box-sizing: border-box; }
  body {
    margin: 0;
    background: var(--bg);
    color: var(--text);
    font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    line-height: 1.45;
  }
  .page { max-width: 1420px; margin: 0 auto; padding: 28px; }
  h1 { margin: 0; font-size: clamp(1.55rem, 2.5vw, 2.2rem); font-weight: 650; }
  h2 { margin: 0 0 14px; font-size: 1.05rem; font-weight: 650; }
  p { margin: 0; }
  header { margin-bottom: 22px; }
  .subtitle { color: var(--muted); margin-top: 6px; max-width: 900px; }
  .source-note { color: var(--muted); margin-top: 8px; font-size: 0.88rem; }
  .panel, .stat {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 12px;
    box-shadow: var(--shadow);
  }
  .controls { padding: 16px; margin-bottom: 16px; }
  .control-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
    gap: 12px 16px;
    align-items: end;
  }
  label { display: block; font-size: 0.8rem; color: var(--muted); margin-bottom: 5px; }
  select, input[type="range"] {
    width: 100%;
    color: var(--text);
    background: var(--surface-soft);
    border: 1px solid var(--border);
    border-radius: 8px;
    min-height: 38px;
    padding: 7px 10px;
    font: inherit;
  }
  input[type="range"] { padding: 0; }
  .check-wrap {
    display: flex;
    align-items: center;
    min-height: 38px;
    gap: 8px;
    color: var(--text);
  }
  .check-wrap label { margin: 0; color: var(--text); font-size: 0.92rem; }
  .rank-label { margin-top: 15px; color: var(--muted); font-size: 0.8rem; }
  .rank-tabs { display: flex; flex-wrap: wrap; gap: 7px; margin-top: 7px; }
  .rank-tab {
    border: 1px solid var(--border);
    color: var(--text);
    background: var(--surface-soft);
    border-radius: 999px;
    padding: 7px 13px;
    font: inherit;
    cursor: pointer;
  }
  .rank-tab.active { background: var(--accent); color: #fff; border-color: var(--accent); }
  .stats {
    display: grid;
    grid-template-columns: repeat(4, minmax(150px, 1fr));
    gap: 12px;
    margin-bottom: 16px;
  }
  .stat { padding: 14px 16px; }
  .stat-label { color: var(--muted); font-size: 0.78rem; }
  .stat-value { margin-top: 3px; font-size: 1.35rem; font-weight: 650; font-variant-numeric: tabular-nums; }
  .stat-detail { margin-top: 2px; color: var(--muted); font-size: 0.76rem; }
  .main-grid { display: grid; grid-template-columns: minmax(0, 1.15fr) minmax(0, 0.85fr); gap: 16px; }
  .panel { padding: 18px; min-width: 0; }
  .chart-note { color: var(--muted); font-size: 0.82rem; margin: -8px 0 14px; }
  .bars { display: grid; gap: 9px; }
  .bar-row { display: grid; grid-template-columns: minmax(150px, 240px) minmax(120px, 1fr) 115px; gap: 10px; align-items: center; }
  .bar-label { overflow: hidden; text-overflow: ellipsis; white-space: nowrap; font-size: 0.84rem; }
  .bar-track { height: 14px; background: var(--track); border-radius: 7px; overflow: hidden; }
  .bar-fill { height: 100%; min-width: 2px; border-radius: 7px; }
  .bar-value { text-align: right; font-size: 0.82rem; font-variant-numeric: tabular-nums; }
  .legend { display: flex; align-items: center; gap: 9px; margin-top: 14px; color: var(--muted); font-size: 0.78rem; }
  .gradient { width: 120px; height: 9px; border-radius: 5px; background: linear-gradient(90deg, hsl(8 68% 51%), hsl(55 68% 47%), hsl(118 52% 42%)); }
  #scatter { display: block; width: 100%; min-height: 330px; }
  .axis { stroke: var(--border); stroke-width: 1; }
  .gridline { stroke: var(--border); stroke-width: 1; opacity: 0.55; }
  .axis-text { fill: var(--muted); font-size: 11px; }
  .axis-title { fill: var(--text); font-size: 12px; }
  .point { stroke: var(--surface); stroke-width: 1.4; cursor: pointer; }
  .table-panel { margin-top: 16px; }
  .table-wrap { overflow-x: auto; }
  table { width: 100%; border-collapse: collapse; font-size: 0.83rem; }
  th, td { border-bottom: 1px solid var(--border); padding: 8px 9px; text-align: left; white-space: nowrap; }
  th { color: var(--muted); font-weight: 600; background: var(--surface-soft); }
  td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
  .continuity-cell { display: flex; align-items: center; gap: 7px; }
  .continuity-dot { width: 10px; height: 10px; border-radius: 50%; flex: 0 0 auto; }
  .explain {
    margin-top: 16px;
    color: var(--muted);
    font-size: 0.86rem;
    display: grid;
    grid-template-columns: repeat(3, minmax(0, 1fr));
    gap: 14px;
  }
  .explain strong { color: var(--text); }
  .empty { padding: 36px 8px; text-align: center; color: var(--muted); }
  footer { margin-top: 18px; color: var(--muted); font-size: 0.76rem; }
  .tooltip {
    position: fixed;
    display: none;
    z-index: 10;
    pointer-events: none;
    max-width: 300px;
    background: var(--surface);
    color: var(--text);
    border: 1px solid var(--border);
    border-radius: 8px;
    box-shadow: var(--shadow);
    padding: 9px 11px;
    font-size: 0.78rem;
  }
  @media (max-width: 980px) {
    .main-grid { grid-template-columns: 1fr; }
    .stats { grid-template-columns: repeat(2, minmax(140px, 1fr)); }
  }
  @media (max-width: 650px) {
    .page { padding: 16px; }
    .stats { grid-template-columns: 1fr 1fr; }
    .bar-row { grid-template-columns: minmax(105px, 150px) 1fr; }
    .bar-value { grid-column: 2; }
    .explain { grid-template-columns: 1fr; }
  }
</style>
</head>
<body>
<main class="page">
  <header>
    <h1>__REPORT_TITLE__</h1>
    <p class="subtitle">Explore taxonomic bin size and assembly fragmentation across samples. Rank, sample, and bin-set filters update every view.</p>
    <p class="source-note">This report uses only <code>summary/bin_metrics.tsv</code>. CheckM output is not read or displayed.</p>
  </header>

  <section class="panel controls" aria-label="Report controls">
    <div class="control-grid">
      <div>
        <label for="sample-select">Sample</label>
        <select id="sample-select"></select>
      </div>
      <div id="bin-set-control">
        <label for="bin-set-select">Bin set</label>
        <select id="bin-set-select"></select>
      </div>
      <div>
        <label for="sort-select">Order bins by</label>
        <select id="sort-select">
          <option value="totalBp">Bin size</option>
          <option value="n50">N50</option>
          <option value="contigs">Contig count</option>
          <option value="cpm">Contigs per Mb</option>
        </select>
      </div>
      <div>
        <label for="top-range">Bins shown: <span id="top-value">30</span></label>
        <input id="top-range" type="range" min="5" max="100" step="5" value="30">
      </div>
      <div class="check-wrap">
        <input id="unclassified-check" type="checkbox">
        <label for="unclassified-check">Include Unclassified</label>
      </div>
    </div>
    <div class="rank-label">Taxonomic rank</div>
    <div id="rank-tabs" class="rank-tabs" role="group" aria-label="Taxonomic rank"></div>
  </section>

  <section class="stats" aria-live="polite">
    <div class="stat"><div class="stat-label">Displayed bins</div><div id="stat-bins" class="stat-value">—</div><div id="stat-samples" class="stat-detail"></div></div>
    <div class="stat"><div class="stat-label">Combined bin size</div><div id="stat-bp" class="stat-value">—</div><div class="stat-detail">Within the current filters</div></div>
    <div class="stat"><div class="stat-label">Median N50</div><div id="stat-n50" class="stat-value">—</div><div class="stat-detail">Higher generally means more contiguous</div></div>
    <div class="stat"><div class="stat-label">Median contigs per Mb</div><div id="stat-cpm" class="stat-value">—</div><div class="stat-detail">Higher means more fragmented</div></div>
  </section>

  <section class="main-grid">
    <div class="panel">
      <h2>Bin sizes</h2>
      <p class="chart-note">Bar length is total assembled sequence. Color is continuity: N50 divided by total bin size.</p>
      <div id="bars" class="bars"></div>
      <div class="legend"><span>Fragmented</span><span class="gradient"></span><span>More contiguous</span></div>
    </div>
    <div class="panel">
      <h2>Size versus fragmentation</h2>
      <p class="chart-note">Points farther right are larger; points higher up have more contigs per Mb.</p>
      <svg id="scatter" viewBox="0 0 620 360" role="img" aria-label="Scatter plot of bin size and contigs per megabase"></svg>
    </div>
  </section>

  <section class="panel table-panel">
    <h2>Bin details</h2>
    <div class="table-wrap">
      <table>
        <thead><tr>
          <th>Sample</th><th>Taxon</th><th class="num">Bin size</th><th class="num">Contigs</th>
          <th class="num">N50</th><th class="num">L50</th><th class="num">Contigs/Mb</th>
          <th class="num">Largest contig</th><th>Continuity</th>
        </tr></thead>
        <tbody id="table-body"></tbody>
      </table>
    </div>
  </section>

  <section class="explain" aria-label="Metric interpretation">
    <p><strong>Bin size</strong><br>Total base pairs across all contigs assigned to that taxon at the selected rank.</p>
    <p><strong>N50 and L50</strong><br>Higher N50 and lower L50 indicate that more of the bin is contained in fewer, longer contigs.</p>
    <p><strong>Contigs per Mb</strong><br>Contig count normalized for bin size. Higher values generally indicate greater fragmentation.</p>
  </section>

  <footer>Generated __GENERATED_AT__ from __SOURCE_COUNT__ bin-metrics file(s). No CheckM data included.</footer>
</main>
<div id="tooltip" class="tooltip" role="tooltip"></div>

<script id="bin-data" type="application/json">__DATA_JSON__</script>
<script>
(() => {
  "use strict";
  const data = JSON.parse(document.getElementById("bin-data").textContent);
  const rankOrder = ["domain", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"];
  const ranks = [...new Set(data.map(d => d.rank))].sort((a, b) => {
    const ai = rankOrder.indexOf(a.toLowerCase());
    const bi = rankOrder.indexOf(b.toLowerCase());
    return (ai < 0 ? 999 : ai) - (bi < 0 ? 999 : bi) || a.localeCompare(b);
  });
  const samples = [...new Set(data.map(d => d.sample))].sort();
  const binSets = [...new Set(data.map(d => d.binSet))].sort();
  const state = {
    rank: ranks.includes("genus") ? "genus" : ranks[0],
    sample: "__all__",
    binSet: "__all__",
    includeUnclassified: false,
    sort: "totalBp",
    topN: 30
  };

  const $ = id => document.getElementById(id);
  const sampleSelect = $("sample-select");
  const binSetSelect = $("bin-set-select");
  const sortSelect = $("sort-select");
  const topRange = $("top-range");
  const unclassifiedCheck = $("unclassified-check");
  const tooltip = $("tooltip");

  function option(select, value, label) {
    const node = document.createElement("option");
    node.value = value;
    node.textContent = label;
    select.appendChild(node);
  }
  option(sampleSelect, "__all__", "All samples");
  samples.forEach(s => option(sampleSelect, s, s));
  option(binSetSelect, "__all__", "All bin sets");
  binSets.forEach(s => option(binSetSelect, s, s[0].toUpperCase() + s.slice(1)));
  if (binSets.length <= 1) $("bin-set-control").style.display = "none";

  const rankTabs = $("rank-tabs");
  ranks.forEach(rank => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = "rank-tab";
    button.textContent = rank[0].toUpperCase() + rank.slice(1);
    button.dataset.rank = rank;
    button.addEventListener("click", () => { state.rank = rank; render(); });
    rankTabs.appendChild(button);
  });

  sampleSelect.addEventListener("change", () => { state.sample = sampleSelect.value; render(); });
  binSetSelect.addEventListener("change", () => { state.binSet = binSetSelect.value; render(); });
  sortSelect.addEventListener("change", () => { state.sort = sortSelect.value; render(); });
  topRange.addEventListener("input", () => {
    state.topN = Number(topRange.value);
    $("top-value").textContent = topRange.value;
    render();
  });
  unclassifiedCheck.addEventListener("change", () => {
    state.includeUnclassified = unclassifiedCheck.checked;
    render();
  });

  function formatNumber(value, digits = 0) {
    return Number(value).toLocaleString(undefined, { maximumFractionDigits: digits });
  }
  function formatBp(value) {
    const n = Number(value);
    if (n >= 1e9) return `${(n / 1e9).toFixed(2)} Gb`;
    if (n >= 1e6) return `${(n / 1e6).toFixed(2)} Mb`;
    if (n >= 1e3) return `${(n / 1e3).toFixed(1)} kb`;
    return `${formatNumber(n)} bp`;
  }
  function median(values) {
    if (!values.length) return 0;
    const sorted = [...values].sort((a, b) => a - b);
    const mid = Math.floor(sorted.length / 2);
    return sorted.length % 2 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
  }
  function clamp(value, min, max) { return Math.max(min, Math.min(max, value)); }
  function continuityColor(fraction) {
    const f = Math.sqrt(clamp(Number(fraction) || 0, 0, 1));
    const hue = 8 + f * 110;
    return `hsl(${hue} 62% 46%)`;
  }
  function filteredRows() {
    return data.filter(d =>
      d.rank === state.rank &&
      (state.sample === "__all__" || d.sample === state.sample) &&
      (state.binSet === "__all__" || d.binSet === state.binSet) &&
      (state.includeUnclassified || !d.unclassified)
    );
  }
  function sortedRows(rows) {
    const key = state.sort;
    return [...rows].sort((a, b) => Number(b[key]) - Number(a[key]) || a.taxon.localeCompare(b.taxon));
  }
  function labelFor(row) {
    return state.sample === "__all__" ? `${row.sample} · ${row.taxon}` : row.taxon;
  }
  function detailText(row) {
    return `${row.sample} — ${row.taxon}\nSize: ${formatBp(row.totalBp)}\nContigs: ${formatNumber(row.contigs)}\nN50: ${formatBp(row.n50)}\nL50: ${formatNumber(row.l50)}\nContigs/Mb: ${formatNumber(row.cpm, 2)}`;
  }

  function renderStats(rows) {
    $("stat-bins").textContent = formatNumber(rows.length);
    $("stat-samples").textContent = `${new Set(rows.map(d => d.sample)).size} sample(s)`;
    $("stat-bp").textContent = formatBp(rows.reduce((sum, d) => sum + d.totalBp, 0));
    $("stat-n50").textContent = formatBp(median(rows.map(d => d.n50)));
    $("stat-cpm").textContent = formatNumber(median(rows.map(d => d.cpm)), 2);
  }

  function renderBars(rows) {
    const container = $("bars");
    container.replaceChildren();
    if (!rows.length) {
      const empty = document.createElement("div");
      empty.className = "empty";
      empty.textContent = "No bins match the current filters.";
      container.appendChild(empty);
      return;
    }
    const maxBp = Math.max(...rows.map(d => d.totalBp), 1);
    rows.forEach(row => {
      const wrapper = document.createElement("div");
      wrapper.className = "bar-row";
      wrapper.title = detailText(row);
      const label = document.createElement("div");
      label.className = "bar-label";
      label.textContent = labelFor(row);
      const track = document.createElement("div");
      track.className = "bar-track";
      const fill = document.createElement("div");
      fill.className = "bar-fill";
      fill.style.width = `${Math.max(0.5, 100 * row.totalBp / maxBp)}%`;
      fill.style.background = continuityColor(row.n50Fraction);
      track.appendChild(fill);
      const value = document.createElement("div");
      value.className = "bar-value";
      value.textContent = formatBp(row.totalBp);
      wrapper.append(label, track, value);
      container.appendChild(wrapper);
    });
  }

  function svgNode(name, attrs = {}) {
    const node = document.createElementNS("http://www.w3.org/2000/svg", name);
    Object.entries(attrs).forEach(([key, value]) => node.setAttribute(key, value));
    return node;
  }
  function renderScatter(rows) {
    const svg = $("scatter");
    svg.replaceChildren();
    const width = 620, height = 360;
    const margin = { top: 16, right: 20, bottom: 55, left: 70 };
    const plotW = width - margin.left - margin.right;
    const plotH = height - margin.top - margin.bottom;
    if (!rows.length) {
      const text = svgNode("text", { x: width / 2, y: height / 2, "text-anchor": "middle", class: "axis-text" });
      text.textContent = "No bins match the current filters.";
      svg.appendChild(text);
      return;
    }
    const xValues = rows.map(d => Math.log10(Math.max(d.totalBp, 1)));
    const yValues = rows.map(d => Math.log10(Math.max(d.cpm, 0) + 1));
    let xMin = Math.min(...xValues), xMax = Math.max(...xValues);
    let yMin = Math.min(...yValues), yMax = Math.max(...yValues);
    if (xMin === xMax) { xMin -= 0.25; xMax += 0.25; }
    if (yMin === yMax) { yMin = Math.max(0, yMin - 0.25); yMax += 0.25; }
    const padX = (xMax - xMin) * 0.06;
    const padY = (yMax - yMin) * 0.08;
    xMin -= padX; xMax += padX; yMin = Math.max(0, yMin - padY); yMax += padY;
    const xScale = value => margin.left + (value - xMin) / (xMax - xMin) * plotW;
    const yScale = value => margin.top + plotH - (value - yMin) / (yMax - yMin) * plotH;

    for (let i = 0; i <= 4; i++) {
      const xLog = xMin + i * (xMax - xMin) / 4;
      const x = xScale(xLog);
      svg.appendChild(svgNode("line", { x1: x, y1: margin.top, x2: x, y2: margin.top + plotH, class: "gridline" }));
      const label = svgNode("text", { x, y: height - 29, "text-anchor": "middle", class: "axis-text" });
      label.textContent = formatBp(10 ** xLog);
      svg.appendChild(label);
      const yLog = yMin + i * (yMax - yMin) / 4;
      const y = yScale(yLog);
      svg.appendChild(svgNode("line", { x1: margin.left, y1: y, x2: width - margin.right, y2: y, class: "gridline" }));
      const yLabel = svgNode("text", { x: margin.left - 9, y: y + 4, "text-anchor": "end", class: "axis-text" });
      yLabel.textContent = formatNumber((10 ** yLog) - 1, 1);
      svg.appendChild(yLabel);
    }
    svg.appendChild(svgNode("line", { x1: margin.left, y1: margin.top + plotH, x2: width - margin.right, y2: margin.top + plotH, class: "axis" }));
    svg.appendChild(svgNode("line", { x1: margin.left, y1: margin.top, x2: margin.left, y2: margin.top + plotH, class: "axis" }));
    const xTitle = svgNode("text", { x: margin.left + plotW / 2, y: height - 5, "text-anchor": "middle", class: "axis-title" });
    xTitle.textContent = "Total bin size (log scale)";
    svg.appendChild(xTitle);
    const yTitle = svgNode("text", { x: 15, y: margin.top + plotH / 2, "text-anchor": "middle", class: "axis-title", transform: `rotate(-90 15 ${margin.top + plotH / 2})` });
    yTitle.textContent = "Contigs per Mb (log scale)";
    svg.appendChild(yTitle);

    rows.forEach(row => {
      const point = svgNode("circle", {
        cx: xScale(Math.log10(Math.max(row.totalBp, 1))),
        cy: yScale(Math.log10(Math.max(row.cpm, 0) + 1)),
        r: 5.5,
        fill: continuityColor(row.n50Fraction),
        class: "point"
      });
      const title = svgNode("title");
      title.textContent = detailText(row);
      point.appendChild(title);
      point.addEventListener("mousemove", event => {
        tooltip.style.display = "block";
        tooltip.style.left = `${event.clientX + 14}px`;
        tooltip.style.top = `${event.clientY + 14}px`;
        tooltip.textContent = detailText(row);
      });
      point.addEventListener("mouseleave", () => { tooltip.style.display = "none"; });
      svg.appendChild(point);
    });
  }

  function renderTable(rows) {
    const body = $("table-body");
    body.replaceChildren();
    rows.forEach(row => {
      const tr = document.createElement("tr");
      const continuity = Math.max(0, Math.min(1, row.n50Fraction || 0));
      const values = [
        [row.sample, ""], [row.taxon, ""], [formatBp(row.totalBp), "num"],
        [formatNumber(row.contigs), "num"], [formatBp(row.n50), "num"],
        [formatNumber(row.l50), "num"], [formatNumber(row.cpm, 2), "num"],
        [formatBp(row.largest), "num"]
      ];
      values.forEach(([value, cls]) => {
        const td = document.createElement("td");
        if (cls) td.className = cls;
        td.textContent = value;
        tr.appendChild(td);
      });
      const continuityTd = document.createElement("td");
      const continuityWrap = document.createElement("div");
      continuityWrap.className = "continuity-cell";
      const dot = document.createElement("span");
      dot.className = "continuity-dot";
      dot.style.background = continuityColor(continuity);
      const text = document.createElement("span");
      text.textContent = `${(100 * continuity).toFixed(1)}%`;
      continuityWrap.append(dot, text);
      continuityTd.appendChild(continuityWrap);
      tr.appendChild(continuityTd);
      body.appendChild(tr);
    });
  }

  function render() {
    [...rankTabs.children].forEach(button => button.classList.toggle("active", button.dataset.rank === state.rank));
    const filtered = filteredRows();
    const shown = sortedRows(filtered).slice(0, state.topN);
    renderStats(filtered);
    renderBars(shown);
    renderScatter(shown);
    renderTable(shown);
  }
  render();
})();
</script>
</body>
</html>
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a self-contained HTML explorer from MetaHopper "
            "summary/bin_metrics.tsv files. CheckM output is ignored."
        )
    )
    parser.add_argument(
        "-i", "--input", type=Path, default=Path("."),
        help=(
            "Parent directory containing metahop_* samples, one metahop_* sample directory, "
            "a summary directory, or a bin_metrics.tsv file (default: .)"
        ),
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("metahopper_bin_report.html"),
        help="Output HTML file (default: metahopper_bin_report.html)",
    )
    parser.add_argument(
        "--sample-glob", default="metahop_*",
        help="Sample-directory glob when --input is a parent directory (default: metahop_*)",
    )
    parser.add_argument(
        "--title", default="MetaHopper Bin Explorer",
        help="Report title (default: MetaHopper Bin Explorer)",
    )
    return parser.parse_args()


def discover_metric_files(input_path: Path, sample_glob: str) -> list[Path]:
    input_path = input_path.resolve()
    if input_path.is_file():
        return [input_path] if input_path.name == "bin_metrics.tsv" else []

    direct_candidates = [
        input_path / "bin_metrics.tsv",
        input_path / "summary" / "bin_metrics.tsv",
        input_path / "metahopper_bin_summary" / "bin_metrics.tsv",
    ]
    for candidate in direct_candidates:
        if candidate.is_file():
            return [candidate]

    candidates = sorted(
        path / "summary" / "bin_metrics.tsv"
        for path in input_path.glob(sample_glob)
        if path.is_dir()
    )
    return [path for path in candidates if path.is_file()]


def as_int(row: dict[str, str], field: str) -> int:
    try:
        return int(float(row.get(field, "0") or 0))
    except ValueError:
        return 0


def as_float(row: dict[str, str], field: str) -> float:
    try:
        return float(row.get(field, "0") or 0)
    except ValueError:
        return 0.0


def fallback_sample(metric_file: Path) -> str:
    if metric_file.parent.name == "summary":
        name = metric_file.parent.parent.name
        return name[len("metahop_"):] if name.startswith("metahop_") else name
    return metric_file.parent.name


def read_metrics(metric_files: list[Path]) -> list[dict]:
    records = []
    for metric_file in metric_files:
        with open(metric_file, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"rank", "bin_id", "num_contigs", "total_bp", "N50_bp"}
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise ValueError(
                    f"{metric_file} is missing required columns: {', '.join(sorted(missing))}"
                )
            for row in reader:
                total_bp = as_int(row, "total_bp")
                n50 = as_int(row, "N50_bp")
                bin_id = row.get("bin_id", "") or "Unclassified"
                records.append({
                    "sample": row.get("sample", "") or fallback_sample(metric_file),
                    "binSet": row.get("bin_set", "") or "preliminary",
                    "rank": row.get("rank", "") or "unknown",
                    "bin": bin_id,
                    "taxon": row.get("taxon_label", "") or bin_id.replace("_", " "),
                    "unclassified": (
                        (row.get("is_unclassified", "").lower() == "true")
                        or bin_id.lower() == "unclassified"
                    ),
                    "contigs": as_int(row, "num_contigs"),
                    "totalBp": total_bp,
                    "largest": as_int(row, "largest_contig_bp"),
                    "n50": n50,
                    "l50": as_int(row, "L50_contigs"),
                    "gc": round(as_float(row, "GC_percent"), 4),
                    "cpm": round(as_float(row, "contigs_per_Mbp"), 4),
                    "n50Fraction": round(
                        as_float(row, "N50_fraction")
                        if row.get("N50_fraction", "")
                        else (n50 / total_bp if total_bp else 0.0),
                        6,
                    ),
                })
    return records


def escape_json_for_html(data: list[dict]) -> str:
    return json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")


def build_html(records: list[dict], title: str, source_count: int) -> str:
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    return (
        HTML_TEMPLATE
        .replace("__REPORT_TITLE__", html.escape(title))
        .replace("__GENERATED_AT__", generated)
        .replace("__SOURCE_COUNT__", str(source_count))
        .replace("__DATA_JSON__", escape_json_for_html(records))
    )


def main() -> int:
    args = parse_args()
    metric_files = discover_metric_files(args.input, args.sample_glob)
    if not metric_files:
        print(
            "ERROR: no summary/bin_metrics.tsv files found. Point --input at the parent "
            "directory, one metahop_* directory, a summary directory, or bin_metrics.tsv.",
            file=sys.stderr,
        )
        return 1
    try:
        records = read_metrics(metric_files)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    if not records:
        print("ERROR: the discovered bin_metrics.tsv files contain no bin rows", file=sys.stderr)
        return 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(build_html(records, args.title, len(metric_files)), encoding="utf-8")
    sample_count = len({record["sample"] for record in records})
    rank_count = len({record["rank"] for record in records})
    print(
        f"Wrote {args.output.resolve()} with {len(records)} bins, "
        f"{sample_count} sample(s), and {rank_count} rank(s).",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
