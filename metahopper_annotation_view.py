"""Embed the per-bin gene/protein explorer, preserving Bakta's recorded evidence."""
from __future__ import annotations
import csv
import html
import io
import json
from pathlib import Path
from metahopper_report_figures import CHART_JS, FIGURE_CSS, json_script

FIELDS = ['locus_tag', 'gene', 'product', 'type', 'contig', 'start', 'end', 'strand', 'db_xrefs']


def bin_key(rank, name):
    return rank + '/' + name


def bin_button(rank, name):
    return ('<button type="button" class="annotation-bin-link" data-annotation-bin="'
            + html.escape(bin_key(rank, name), quote=True) + '">' + html.escape(name) + '</button>')


def safe_text(path, root):
    path = Path(path)
    if not path.resolve().is_relative_to(Path(root).resolve()):
        raise ValueError('Annotation file points outside this sample')
    before = path.stat()
    text = path.read_text(encoding='utf-8')
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError('Annotation changed while building the report; retry after annotation finishes')
    return text


def flag(value):
    if isinstance(value, str):
        return value.strip().lower() not in ('', 'false', 'no', '0', 'none', 'null')
    if isinstance(value, dict):
        return True
    return bool(value)


def optional_flag(value):
    if value is None or isinstance(value, str) and value.strip().lower() in ('', 'na', 'n/a', 'unknown', 'unavailable', 'null', 'none'):
        return None
    return flag(value)


def feature_detail(feature, product, kind, authoritative=False):
    protein = kind in ('cds', 'sorf')
    hypothetical = (optional_flag(feature['hypothetical']) if 'hypothetical' in feature else
                    product.lower().startswith('hypothetical protein') if product else
                    None)
    if not protein:
        hypothetical = None
    key = next((key for key in ('pseudogene', 'pseudo') if key in feature), None)
    pseudo = (flag(feature[key]) if authoritative else optional_flag(feature[key])) if key else False if authoritative else None
    value = feature.get(key) if key else None
    note = ''
    if isinstance(value, dict):
        note = str(value.get('description') or value.get('reason') or '')
        if not note and value:
            note = json.dumps(value, ensure_ascii=False, sort_keys=True)
    elif isinstance(value, str) and value.lower() not in ('true', 'false', 'yes', 'no', '0', '1'):
        note = value
    aa, nt = feature.get('aa') or '', feature.get('nt') or ''
    if not isinstance(aa, str) or not isinstance(nt, str):
        raise ValueError('Bakta feature sequences must be strings')
    return dict(hypothetical=hypothetical, pseudogene=pseudo, pseudogene_note=note,
                truncated=str(feature.get('truncated') or ''), aa=aa, nt=nt,
                aa_length=len(aa.rstrip('*')) if aa else None, nt_length=len(nt) if nt else None)


def records_from_json(data, valid_contigs=None):
    features = data.get('features') if isinstance(data, dict) else None
    if not isinstance(features, list):
        raise ValueError('Bakta JSON has no features array')
    rows, details = [], []
    for number, feature in enumerate(features, 1):
        if not isinstance(feature, dict):
            raise ValueError(f'Bakta feature {number} is not an object')
        kind = str(feature.get('type') or '').lower()
        contig = feature.get('sequence') or feature.get('contig') or feature.get('sequence_id') or ''
        if isinstance(contig, dict):
            contig = contig.get('id', '')
        contig = str(contig)
        if valid_contigs is not None and contig not in valid_contigs:
            raise ValueError(f'Bakta feature {number} belongs to a contig outside this bin: {contig}')
        product = str(feature.get('product') or '')
        refs = feature.get('db_xrefs') or feature.get('dbxref') or []
        refs = refs if isinstance(refs, str) else ';'.join(map(str, refs))
        rows.append([str(feature.get('locus') or feature.get('locus_tag') or feature.get('id') or ''),
                     str(feature.get('gene') or ''), product, kind, contig,
                     str(feature.get('start') or ''), str(feature.get('stop', feature.get('end', ''))),
                     str(feature.get('strand') or ''), refs])
        details.append(feature_detail(feature, product, kind, authoritative=True))
    return rows, details


def fasta_records(text):
    current, pieces = None, []
    for line in text.splitlines():
        if line.startswith('>'):
            if current is not None:
                yield current, ''.join(pieces)
            current, pieces = line[1:].split()[0], []
        elif current is not None:
            pieces.append(''.join(line.split()))
    if current is not None:
        yield current, ''.join(pieces)


def collect_features(layout, rows):
    """Only use completed rows that passed the report's source-bin freshness check."""
    result = []
    for summary in rows:
        rank, name = str(summary.get('rank', '')), str(summary.get('bin', ''))
        if summary.get('status') != 'completed':
            continue
        if any(not x or Path(x).name != x or x in ('.', '..') for x in (rank, name)):
            continue
        directory = layout.base / 'annotation' / rank / name
        relative = (directory / 'features.tsv').relative_to(layout.root).as_posix()
        item = dict(key=bin_key(rank, name), rank=rank, bin=name, features_path=relative,
                    fasta_path=(layout.base / 'bins' / rank / (name + '.fasta')).relative_to(layout.root).as_posix(),
                    protein_path=(directory / 'bakta/annotation.faa').relative_to(layout.root).as_posix(),
                    json_path=(directory / 'bakta/annotation.json').relative_to(layout.root).as_posix(),
                    features=[], details=[], expected_contigs=[], source='', error='')
        try:
            bin_path = layout.root / item['fasta_path']
            if bin_path.is_file():
                item['expected_contigs'] = [key for key, _ in fasta_records(safe_text(bin_path, layout.root))]
            path = layout.root / item['json_path']
            if path.is_file():
                data = json.loads(safe_text(path, layout.root))
                item['features'], item['details'] = records_from_json(data, set(item['expected_contigs']) if item['expected_contigs'] else None)
                item['source'] = 'Bakta JSON'
            else:
                reader = csv.DictReader(io.StringIO(safe_text(layout.root / relative, layout.root)), delimiter='\t')
                missing = set(FIELDS) - set(reader.fieldnames or [])
                if missing:
                    raise ValueError('Feature table is missing columns: ' + ', '.join(sorted(missing)))
                for number, feature in enumerate(reader, 2):
                    if feature.get('rank') not in (None, '', rank) or feature.get('bin') not in (None, '', name):
                        raise ValueError(f'Feature row {number} belongs to a different bin')
                    item['features'].append([str(feature.get(field) or '') for field in FIELDS])
                    item['details'].append(feature_detail(feature, str(feature.get('product') or ''), str(feature.get('type') or '').lower()))
                item['source'] = 'Feature table; pseudogene flags unavailable unless explicitly recorded'
            for suffix, field, length in [('.faa', 'aa', 'aa_length'), ('.ffn', 'nt', 'nt_length')]:
                path = directory / 'bakta' / ('annotation' + suffix)
                if path.is_file() and any(not detail[field] for detail in item['details']):
                    sequences = dict(fasta_records(safe_text(path, layout.root)))
                    for row, detail in zip(item['features'], item['details']):
                        if not detail[field] and row[0] in sequences:
                            detail[field] = sequences[row[0]]
                            detail[length] = len(detail[field].rstrip('*')) if field == 'aa' else len(detail[field])
        except (OSError, ValueError, TypeError, UnicodeError, csv.Error) as exc:
            item.update(features=[], details=[], error=str(exc))
        result.append(item)
    return dict(schema='metahopper-annotation-features-1', run_name=layout.root.name,
                stage=layout.stage, fields=FIELDS, bins=result)


def render_browser(data):
    script = Path(__file__).with_name('metahopper_annotation_browser.js').read_text(encoding='utf-8')
    return ('<style>' + FIGURE_CSS + BROWSER_CSS + '</style>' + BROWSER_HTML
            + json_script('annotation-feature-data', data) + '<script>' + CHART_JS + script + '</script>')


BROWSER_CSS = '''
.annotation-bin-link{border:0;background:transparent;color:var(--acc);padding:0;text-align:left;text-decoration:underline;text-underline-offset:3px;font:inherit;cursor:pointer}
.annotation-bin-link:focus-visible{outline:2px solid var(--acc);outline-offset:3px}
.annotation-controls{display:grid;grid-template-columns:minmax(180px,1.1fr) minmax(220px,2fr) minmax(160px,1fr) minmax(190px,1fr);gap:12px;margin:12px 0}
.annotation-chips{display:flex;gap:8px;flex-wrap:wrap;margin:12px 0}.annotation-chips button{font-size:.84rem;padding:7px 11px}.annotation-chips button.active{background:var(--accs);border-color:var(--acc)}
.annotation-actions{display:flex;gap:7px;flex-wrap:wrap;align-items:center;margin:12px 0}.annotation-actions button,.annotation-actions a{white-space:nowrap;font-size:.8rem}
#annotation-feature-table{min-width:1500px}#annotation-feature-table th{white-space:nowrap}#annotation-feature-table th button{font:inherit;padding:0;border:0;background:none;white-space:nowrap}
#annotation-feature-table td{vertical-align:top;white-space:nowrap;padding:8px 10px}
#annotation-feature-table td.product-cell{white-space:normal;min-width:260px;max-width:420px;overflow-wrap:anywhere}
#annotation-feature-table td.reference-cell{white-space:normal;max-width:240px;min-width:140px;overflow-wrap:anywhere}
.feature-tag{display:inline-block;font-size:.72rem;border-radius:5px;padding:2px 5px;background:#edf2f4;white-space:nowrap}
.feature-tag.hypothetical{background:#fff0d9;color:#8c5b15}.feature-tag.pseudogene{background:#eee7f7;color:#66428c}.feature-tag.named{background:#e3f2ef;color:#246559}
.annotation-feature-link{border:0;background:none;color:var(--acc);font:inherit;text-decoration:underline;padding:0;white-space:nowrap}
.annotation-help{font-size:.81rem;color:var(--mut);margin:6px 0}.annotation-import{font-size:.81rem}.annotation-import input{font-size:.8rem;max-width:340px;padding:4px}
.annotation-dialog{width:min(1000px,95vw);max-height:90vh;border:1px solid var(--line);border-radius:12px;padding:20px;color:var(--ink)}
.annotation-dialog::backdrop{background:rgba(15,23,42,.45)}.annotation-dialog .dialog-head{display:flex;justify-content:space-between;align-items:flex-start;gap:15px}
.annotation-dialog dl{display:grid;grid-template-columns:130px minmax(0,1fr);gap:7px 15px;font-size:.86rem}.annotation-dialog dt{font-weight:600;color:var(--mut)}.annotation-dialog dd{margin:0;overflow-wrap:anywhere}
.sequence-preview{font-family:ui-monospace,SFMono-Regular,Consolas,monospace;font-size:.8rem;white-space:pre-wrap;overflow-wrap:anywhere;max-height:220px;overflow:auto;background:#f4f7f8;border:1px solid var(--line);border-radius:7px;padding:12px;line-height:1.6}
@media(max-width:1000px){.annotation-controls{grid-template-columns:repeat(2,minmax(0,1fr))}}
@media(max-width:600px){.annotation-controls{grid-template-columns:1fr}.annotation-dialog dl{grid-template-columns:1fr}}
'''

BROWSER_HTML = '''
<section class="card" id="annotation-browser">
 <h2 id="annotation-browser-title">Gene and protein explorer</h2>
 <p class="annotation-help">Browse all annotated features or choose a bin. Click a locus tag for its predicted product, sequence and annotation details.</p>
 <div class="annotation-controls">
  <div><label for="annotation-bin">Bin</label><select id="annotation-bin"></select></div>
  <div><label for="annotation-search">Search genes, products, locus tags or contigs</label><input id="annotation-search" type="search" placeholder="e.g. dnaA, transporter, ribosomal" autocomplete="off"></div>
  <div><label for="annotation-type">Feature type</label><select id="annotation-type"></select></div>
  <div><label for="annotation-category">Annotation category</label><select id="annotation-category"><option value="all">All categories</option><option value="named">Non-hypothetical proteins</option><option value="hypothetical">Hypothetical proteins</option><option value="unassigned">Unassigned protein products</option><option value="pseudogene">Pseudogenes (flagged)</option><option value="pseudo-unknown">Pseudogene status unavailable</option></select></div>
 </div>
 <div class="annotation-chips" id="annotation-category-chips"></div>
 <p class="annotation-help">Pseudogene and hypothetical categories can overlap. Pseudogene status follows Bakta flags; older summary tables may not record those flags.</p>
 <div class="annotation-actions"><button type="button" id="annotation-clear">Reset filters</button><button type="button" id="annotation-export">Filtered features TSV</button><button type="button" id="annotation-export-proteins">Filtered proteins FASTA</button><button type="button" id="annotation-export-genes">Filtered genes FASTA</button><a id="annotation-feature-download" class="asset-link" download>All features TSV</a><a id="annotation-protein-download" class="asset-link" download>Proteins FASTA</a><a id="annotation-fasta-download" class="asset-link" download>Bin FASTA</a></div>
 <p class="note" id="annotation-feature-status" role="status" aria-live="polite"></p>
 <div class="scroll"><table id="annotation-feature-table"><thead><tr>
  <th class="l"><button data-feature-sort="bin">Bin</button></th><th class="l"><button data-feature-sort="locus_tag">Locus tag</button></th><th class="l"><button data-feature-sort="gene">Gene</button></th><th class="l"><button data-feature-sort="product">Product</button></th><th class="l"><button data-feature-sort="type">Type</button></th><th class="l">Product category</th><th class="l"><button data-feature-sort="pseudogene">Pseudogene</button></th><th class="l"><button data-feature-sort="contig">Contig</button></th><th><button data-feature-sort="start">Start</button></th><th><button data-feature-sort="end">End</button></th><th class="l">Strand</th><th><button data-feature-sort="aa_length">Length (aa)</button></th><th class="l">Database references</th>
 </tr></thead><tbody id="annotation-feature-body"></tbody></table></div>
 <div class="row" style="margin-top:12px"><button type="button" id="annotation-prev">Previous</button><span class="note" id="annotation-page"></span><button type="button" id="annotation-next">Next</button></div>
 <details class="annotation-import" style="margin-top:14px"><summary>Open an existing Bakta JSON file locally</summary><p class="annotation-help">Select its bin above, then choose that bin’s annotation.json. It is read in this browser. Rebuilding the report on the run directory embeds these records automatically.</p><label for="annotation-load-json">Bakta JSON for the selected bin</label><input id="annotation-load-json" type="file" accept=".json,application/json"><p id="annotation-import-status" class="annotation-help" role="status"></p></details>
</section>
<div class="report-figures" id="annotation-selection-figures">
 <section class="card"><h2>Feature types in the selected bin(s)</h2><div id="annotation-types-figure" class="figure-host"></div><p class="figure-note">Click a bar to filter the feature table. Counts use the loaded gene records.</p></section>
 <section class="card"><h2>Protein product categories</h2><div id="annotation-products-figure" class="figure-host"></div><p class="figure-note">Pseudogene flags may overlap either protein category and are counted separately above.</p></section>
</div>
<dialog class="annotation-dialog" id="annotation-feature-dialog">
 <div class="dialog-head"><div><h2 id="annotation-detail-title">Feature details</h2><p id="annotation-detail-product"></p></div><button type="button" id="annotation-detail-close">Close</button></div>
 <dl id="annotation-detail-fields"></dl>
 <div class="annotation-actions"><button type="button" id="annotation-detail-protein">Download protein FASTA</button><button type="button" id="annotation-detail-gene">Download gene FASTA</button></div>
 <h2>Predicted protein sequence</h2><pre class="sequence-preview" id="annotation-detail-aa"></pre>
 <h2>Nucleotide sequence</h2><pre class="sequence-preview" id="annotation-detail-nt"></pre>
</dialog>
'''
