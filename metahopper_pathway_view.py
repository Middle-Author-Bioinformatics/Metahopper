"""Read PathwayGenie evidence into the offline report without changing its calls."""
from __future__ import annotations
import csv
import fcntl
import hashlib
import html
import io
import json
import math
from pathlib import Path
import re

from metahopper_report_figures import json_script

SCHEMA = 'metahopper-pathways-1'
TABLES = {
    'summary': ('pathway_summary.long.csv', {'bin_id', 'pathway', 'panel_percent', 'panel_status', 'route_status'}),
    'bins': ('bin_qc.csv', {'bin_id', 'status'}),
    'definitions': ('component_definitions.csv', {'pathway', 'component', 'assay_scope'}),
    'components': ('component_evidence.long.csv', {'bin_id', 'component', 'state', 'value', 'copies'}),
    'hits': ('all_model_hits.csv', {'bin_id', 'protein', 'contig', 'start', 'end', 'status', 'target_call'}),
    'neighborhoods': ('neighborhoods.csv', {'bin_id', 'pathway', 'neighborhood_id', 'contig', 'start', 'end'}),
    'operons': ('operon_candidates.csv', {'bin_id', 'pathway', 'neighborhood_id', 'template_id', 'score'}),
}
FILE_SUFFIXES = {'.csv', '.json', '.faa', '.fna', '.gff', '.domtbl', '.tbl', '.txt', '.log'}


def safe_name(value):
    value = str(value)
    if not value or value in ('.', '..') or any(c in value for c in '/\\') or any(ord(c) < 32 for c in value):
        raise ValueError('Invalid PathwayGenie bin/rank identifier')
    return value


def safe_path(path, boundary):
    path = Path(path)
    if not path.resolve().is_relative_to(Path(boundary).resolve()):
        raise ValueError('Pathway file points outside its sample: ' + str(path))
    return path


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def finite(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def read_text(path, boundary, observed):
    path = safe_path(path, boundary)
    before = path.stat()
    content = path.read_bytes()
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError('Pathway outputs changed while reading; retry after the run finishes')
    observed[path.resolve()] = hashlib.sha256(content).hexdigest()
    return content.decode('utf-8-sig')


def csv_rows(path, boundary, observed, required=()):
    reader = csv.DictReader(io.StringIO(read_text(path, boundary, observed)))
    missing = set(required) - set(reader.fieldnames or [])
    if missing:
        raise ValueError(path.name + ' is missing columns: ' + ', '.join(sorted(missing)))
    rows = list(reader)
    if any(None in row for row in rows):
        raise ValueError('Malformed CSV row in ' + path.name)
    return rows


def fasta_records(path, boundary, observed):
    current, sequence = None, []
    for line in read_text(path, boundary, observed).splitlines():
        if line.startswith('>'):
            if current is not None:
                yield current, ''.join(sequence)
            current = line[1:].split()[0]
            sequence = []
        elif current is not None:
            sequence.append(''.join(line.split()))
    if current is not None:
        yield current, ''.join(sequence)


def result_directories(layout, rank='genus', explicit=None):
    rank = safe_name(rank)
    if explicit is not None:
        return [safe_path(Path(explicit), layout.root)]
    return [p for p in [layout.base / 'pathways' / rank,
                       layout.base / 'bins' / 'pathwaygenie',
                       layout.base / 'annotation' / rank / 'pathwaygenie',
                       layout.root / 'pathwaygenie']
            if (p / 'pathway_summary.long.csv').is_file()]


def read_results(directory, sample_root, stage='final', rank='genus', verify_inputs=True, check_wrapper=True):
    """Use long tables, not the optionally zero-filtered display heatmap CSV."""
    root = Path(sample_root).resolve()
    directory = safe_path(directory, root)
    wrapper_lock = directory.parent / ('.' + directory.name + '.metahopper-pathways.lock')
    if check_wrapper and wrapper_lock.is_file():
        with wrapper_lock.open() as handle:
            try:
                fcntl.flock(handle, fcntl.LOCK_SH | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise ValueError('Pathway run/import is active; finish it before generating the report') from exc
            finally:
                fcntl.flock(handle, fcntl.LOCK_UN)
    rank = safe_name(rank)
    base = root / 'final' if stage == 'final' else root
    if (directory / '.pathwaygenie.lock').exists():
        raise ValueError('PathwayGenie output is locked. Finish the run before building the report.')
    observed = {}
    manifest = json.loads(read_text(directory / 'run_manifest.json', root, observed))
    if manifest.get('program') != 'PathwayGenie' or not manifest.get('finished_utc'):
        raise ValueError('PathwayGenie run has no finished manifest; incomplete outputs are not displayed')
    config_path = directory / 'config_used.json'
    config = json.loads(read_text(config_path, root, observed)) if config_path.is_file() else {}
    if config_path.is_file() and manifest.get('config_sha256') and observed[config_path.resolve()] != manifest['config_sha256']:
        raise ValueError('PathwayGenie config_used.json does not match the recorded run checksum')
    data = {key: csv_rows(directory / filename, root, observed, required)
            for key, (filename, required) in TABLES.items()}
    warnings, inputs = [], []
    original_inputs = {}
    for item in manifest.get('bin_inputs', []):
        name = safe_name(item['bin_id'])
        if name in original_inputs:
            raise ValueError('Duplicate bin in PathwayGenie manifest: ' + name)
        original_inputs[name] = item
    bin_map = {}
    for row in data['bins']:
        name = safe_name(row['bin_id'])
        if name in bin_map:
            raise ValueError('Duplicate bin in bin_qc.csv: ' + name)
        expected = original_inputs.get(name, {}).get('sha256')
        candidates = [base / 'bins' / rank / (name + suffix) for suffix in ('.fasta', '.fna', '.fa', '.fas', '.fsa')]
        candidates += [Path(str(p) + '.gz') for p in candidates]
        current = next((safe_path(p, root) for p in candidates if p.is_file()), None)
        state = 'preview_unverified' if not verify_inputs else 'missing_bin' if current is None else 'unverified'
        if verify_inputs and current is not None:
            digest = sha256(current)
            inputs.append(dict(path=current.relative_to(root).as_posix(), sha256=digest))
            if expected:
                state = 'current' if digest == expected else 'stale_input'
        row.update(source_state=state, rank=rank)
        if state not in ('current', 'preview_unverified'):
            warnings.append(name + ': ' + state.replace('_', ' ') + '; prior pathway calls are hidden')
        bin_map[name] = row
    trusted = {name for name, row in bin_map.items()
               if row['source_state'] in ('current', 'preview_unverified') and row['status'] != 'failed'}
    seen = set()
    for row in data['summary']:
        name, pid = row['bin_id'], row['pathway']
        if name not in bin_map or (name, pid) in seen:
            raise ValueError('Unknown/duplicate bin-pathway pair in pathway summary: ' + name + '/' + pid)
        seen.add((name, pid))
        number = finite(row['panel_percent'])
        if number is not None and not 0 <= number <= 100:
            raise ValueError('Pathway panel percentage is outside 0–100')
        row['panel_percent'] = number
        row['source_state'] = bin_map[name]['source_state']
        if name not in trusted:
            status = 'bin_failed' if bin_map[name]['status'] == 'failed' else row['source_state']
            row.update(panel_percent=None, panel_status=status, route_status=status,
                       detected_scored_components=None, total_scored_components=None,
                       detected_components='', missing_scored_components='', candidate_components='',
                       unmeasured_route_steps='', operon_support=None)
    for key in ('hits', 'components', 'neighborhoods', 'operons'):
        if any(row['bin_id'] not in bin_map for row in data[key]):
            raise ValueError('Unknown bin in PathwayGenie ' + key)
        data[key] = [row for row in data[key] if row['bin_id'] in trusted]
    component_keys = [(r['bin_id'], r['component']) for r in data['components']]
    if len(component_keys) != len(set(component_keys)):
        raise ValueError('Duplicate bin/component evidence rows')
    for row in data['bins']:
        hits = [h for h in data['hits'] if h['bin_id'] == row['bin_id']]
        row['accepted_models'] = sum(h['status'] == 'accepted' for h in hits)
        row['accepted_proteins'] = len({h['protein'] for h in hits if h['status'] == 'accepted'})
        row['review_hits'] = sum(h['status'] != 'accepted' or h['target_call'] != 'specific' for h in hits)
    proteins = {}
    folders = {}
    for path in sorted((directory / 'bins').glob('*/bin_status.json')):
        record = json.loads(read_text(path, root, observed))
        name = record.get('bin_id')
        if name in folders:
            raise ValueError('Duplicate PathwayGenie per-bin output: ' + str(name))
        folders[name] = path.parent
    for name in trusted:
        folder = folders.get(name)
        if folder is None:
            warnings.append(name + ': no per-bin ORF/sequence folder found')
            continue
        wanted = {h['protein'] for h in data['hits'] if h['bin_id'] == name}
        neighborhoods = [n for n in data['neighborhoods'] if n['bin_id'] == name]
        orf_path = folder / 'orfs.csv'
        for hit in data['hits']:
            if hit['bin_id'] == name:
                proteins[name + '\x1f' + hit['protein']] = dict(protein=hit['protein'],bin_id=name,
                    original_contig=hit['contig'],start=hit['start'],end=hit['end'],strand=hit.get('strand',''),
                    partial=hit.get('partial',''),length=hit.get('protein_length',''),aa='',nt='')
        if not orf_path.is_file():
            warnings.append(name + ': orfs.csv is unavailable; using coordinates recorded in the hit table')
        for row in csv_rows(orf_path, root, observed, {'protein', 'original_contig', 'start', 'end', 'strand'}) if orf_path.is_file() else []:
            selected = [n for n in neighborhoods if n['contig'] == row['original_contig']
                        and int(row['end']) >= int(n['start']) and int(row['start']) <= int(n['end'])]
            if selected:
                wanted.add(row['protein'])
            if row['protein'] not in wanted:
                continue
            key = name + '\x1f' + row['protein']
            proteins[key] = dict(row, bin_id=name, aa='', nt='')
            for neighborhood in selected:
                neighborhood.setdefault('orfs', []).append(row['protein'])
        for suffix, field in [('proteins.faa', 'aa'), ('genes.fna', 'nt')]:
            path = folder / suffix
            if path.is_file():
                for protein, sequence in fasta_records(path, root, observed):
                    key = name + '\x1f' + protein
                    if key in proteins:
                        proteins[key][field] = sequence
    ordered = []
    for pathway in config.get('pathways', []):
        if any(r['pathway'] == pathway['id'] for r in data['summary']):
            ordered.append(pathway)
    for pid in dict.fromkeys(r['pathway'] for r in data['summary']):
        if not any(p['id'] == pid for p in ordered):
            ordered.append(dict(id=pid, label=pid.replace('_', ' '), notes=''))
    if not verify_inputs:
        warnings.append('Preview of the supplied output files. Input-bin checksums have not been compared with current MetaHopper FASTAs.')
    assets = []
    for path in sorted(directory.rglob('*')):
        if not path.is_file() or any(part.startswith('.') for part in path.relative_to(directory).parts) or path.suffix not in FILE_SUFFIXES:
            continue
        path = safe_path(path, root)
        digest = sha256(path)
        if path.resolve() in observed and digest != observed[path.resolve()]:
            raise ValueError('Pathway outputs changed during report generation; retry after the run finishes')
        assets.append(dict(path=path.relative_to(root).as_posix(), label=path.relative_to(directory).as_posix(),
                           size=path.stat().st_size, sha256=digest))
    if (directory / '.pathwaygenie.lock').exists() or sha256(directory / 'run_manifest.json') != observed[(directory / 'run_manifest.json').resolve()]:
        raise ValueError('PathwayGenie run changed during report generation; retry when it finishes')
    arguments = manifest.get('arguments', {})
    parameters = {key: arguments.get(key) for key in ('prodigal_mode', 'min_hmm_coverage', 'ambiguity_margin',
                  'overlap_fraction', 'max_intervening_genes', 'max_neighborhood_bp', 'operon_order',
                  'operon_strand', 'operon_max_missing', 'operon_min_fraction', 'skip_context')}
    return dict(data, schema=SCHEMA, run_name=root.name, stage=stage, rank=rank,
                directory=directory.relative_to(root).as_posix(), pathways=ordered, proteins=proteins,
                files=assets, inputs=inputs, warnings=warnings, parameters=parameters,
                version=manifest.get('version', ''), finished_utc=manifest.get('finished_utc', ''),
                tool_versions=manifest.get('versions', {}), preview=not verify_inputs)


def collect_pathways(layout, rank='genus', explicit=None):
    directories = result_directories(layout, rank, explicit)
    if not directories:
        return None
    data = read_results(directories[0], layout.root, layout.stage, rank)
    if len(directories) > 1:
        data['warnings'].append('Multiple result folders found; using ' + data['directory'] + '. Use --pathways-dir to select another.')
    return data


def render_pathways(data=None, error=''):
    if data is None:
        return ('<section class="card"><h2>HMM pathway annotation</h2><p class="note">'
                + html.escape(error or 'No PathwayGenie results found. Run metahopper_pathways.py, or import an existing PathwayGenie output folder, then regenerate the report.') + '</p></section>')
    script = Path(__file__).with_name('metahopper_pathway_browser.js').read_text()
    return ('<style>' + PATHWAY_CSS + '</style>' + PATHWAY_HTML
            + json_script('metahopper-pathways-data', data) + '<script>' + script + '</script>')


def render_annotation_sections(genes, pathways, default='genes'):
    return ('<style>.annotation-section-tabs{display:flex;gap:7px;margin:0 0 14px}.annotation-section-tabs button.active{background:var(--acc);color:white;border-color:var(--acc)}.annotation-section[hidden]{display:none!important}</style>'
            '<div class="annotation-section-tabs" role="tablist" aria-label="Annotation results">'
            '<button type="button" id="annotation-genes-tab" role="tab" aria-controls="annotation-genes-section">Genes &amp; proteins</button>'
            '<button type="button" id="annotation-pathways-tab" role="tab" aria-controls="annotation-pathways-section">Pathways</button></div>'
            '<div class="annotation-section" id="annotation-genes-section" role="tabpanel" aria-labelledby="annotation-genes-tab">' + genes + '</div>'
            '<div class="annotation-section" id="annotation-pathways-section" role="tabpanel" aria-labelledby="annotation-pathways-tab">' + pathways + '</div><script>'
            + r'''(function(){function show(name){['genes','pathways'].forEach(function(id){var on=id===name,button=document.getElementById('annotation-'+id+'-tab');document.getElementById('annotation-'+id+'-section').hidden=!on;button.setAttribute('aria-selected',String(on));button.classList.toggle('active',on);});}['genes','pathways'].forEach(function(id){document.getElementById('annotation-'+id+'-tab').addEventListener('click',function(){show(id);});});document.addEventListener('metahopper:annotation-bin',function(){show('genes');},true);document.addEventListener('click',function(e){if(e.target.closest&&e.target.closest('[data-annotation-bin]'))show('genes');},true);show(''' + json.dumps(default) + ''');})();</script>''')


PATHWAY_CSS = '''
.pathway-controls{display:grid;grid-template-columns:minmax(180px,1fr) minmax(180px,1fr) minmax(230px,1.5fr);gap:12px;margin:12px 0}
.pathway-actions{display:flex;gap:7px;flex-wrap:wrap;align-items:center;margin:10px 0}.pathway-actions a,.pathway-actions button{font-size:.8rem;white-space:nowrap}
.pathway-actions .pathway-file-link{display:inline-block;max-width:100%;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.pathway-check{display:flex;gap:7px;align-items:center;font-size:.83rem;color:var(--mut)}.pathway-check input{width:auto}
.pathway-scroll{overflow:auto;max-height:420px;border:1px solid var(--line);border-radius:8px}.pathway-scroll table{width:max-content;min-width:100%}
.pathway-scroll th,.pathway-scroll td{padding:7px 10px;vertical-align:top}.pathway-scroll td.wrap-text{white-space:normal;min-width:190px;max-width:320px}
.pathway-scroll .pathway-text-button{font:inherit;background:none;border:0;padding:0;color:var(--acc);text-decoration:underline;text-align:left;white-space:nowrap}
.pathway-tag{display:inline-block;white-space:nowrap;font-size:.75rem;border-radius:5px;padding:2px 5px;background:#eef2f4;color:#51616c}
.pathway-tag.ok{background:#e1f1e9;color:#246442}.pathway-tag.partial{background:#fff0d9;color:#8b5d16}.pathway-tag.bad{background:#f9e5e2;color:#9d3b30}
.pathway-map table{font-size:.77rem}.pathway-map th{min-width:95px;white-space:normal;max-width:115px;text-align:center}.pathway-map th:first-child{min-width:230px;text-align:left;white-space:nowrap;position:sticky;left:0;z-index:2}
.pathway-map td{padding:2px;min-width:95px}.pathway-map button{width:100%;border:0;border-radius:3px;padding:8px 4px;font-size:.78rem;font-variant-numeric:tabular-nums}.pathway-map button:focus-visible{outline:2px solid #16222b;outline-offset:-3px}
.pathway-map tbody th{background:var(--card);font-weight:500}.pathway-key{display:flex;gap:12px;flex-wrap:wrap;font-size:.78rem;color:var(--mut);margin:8px 0}.pathway-key i{display:inline-block;width:12px;height:12px;margin-right:4px;vertical-align:middle;border:1px solid var(--line)}
#pathway-neighborhood-figure{overflow:auto;margin-top:10px}#pathway-neighborhood-figure svg{display:block;width:100%;min-width:650px;max-height:180px}
.pathway-warning{background:#fff6e9;border-left:3px solid #dca54f;padding:9px 12px;font-size:.83rem;margin:10px 0}
.pathway-dialog{width:min(1000px,95vw);max-height:90vh;border:1px solid var(--line);border-radius:12px;padding:20px;color:var(--ink)}.pathway-dialog::backdrop{background:rgba(15,23,42,.45)}
.pathway-sequence{max-height:200px;overflow:auto;white-space:pre-wrap;overflow-wrap:anywhere;background:#f4f7f8;padding:12px;border:1px solid var(--line);border-radius:7px;font-size:.8rem}
.pathway-dialog dl{display:grid;grid-template-columns:140px minmax(0,1fr);gap:8px;font-size:.85rem}.pathway-dialog dd{margin:0;overflow-wrap:anywhere}.pathway-dialog dt{font-weight:600}
@media(max-width:800px){.pathway-controls{grid-template-columns:1fr}.pathway-dialog dl{grid-template-columns:1fr}}
'''

PATHWAY_HTML = '''
<section class="card" id="pathway-browser">
 <h2>PathwayGenie functional annotation</h2>
 <p class="note">HMM evidence from independently predicted Prodigal genes. Panel coverage measures the configured, scored components; it is not pathway activity or abundance. Unmeasured steps can remain at 100% panel coverage. Taxonomic bins can contain multiple genomes, so combined component evidence need not belong to one organism.</p>
 <div id="pathway-warnings"></div><div class="stats" id="pathway-stats"></div>
 <div class="pathway-controls"><div><label for="pathway-bin">Bin</label><select id="pathway-bin"></select></div><div><label for="pathway-pathway">Pathway</label><select id="pathway-pathway"></select></div><div><label for="pathway-search">Search bins, pathways, targets or model hits</label><input id="pathway-search" type="search" placeholder="e.g. Wolbachia, biotin, bioB, TIGR" autocomplete="off"></div></div>
 <div class="pathway-actions"><label class="pathway-check"><input type="checkbox" id="pathway-hide-zero" checked>Hide bins with all-zero panel results</label><button type="button" id="pathway-reset">Reset filters</button><button type="button" id="pathway-export-summary">Filtered pathway CSV</button><a id="pathway-summary-download" class="asset-link" download>Original pathway CSV</a><button type="button" id="pathway-bakta">Bakta annotations for this bin</button></div>
 <p class="note" id="pathway-selection-status" role="status" aria-live="polite"></p>
 <div class="pathway-scroll"><table id="pathway-summary-table"><thead><tr><th class="l">Bin</th><th class="l">Pathway</th><th>Panel coverage</th><th>Detected / scored</th><th class="l">Panel status</th><th class="l">Selected route</th><th>Operon support</th><th class="l">Input / run</th></tr></thead><tbody id="pathway-summary-body"></tbody></table></div>
 <div class="pathway-actions"><button type="button" id="pathway-summary-prev">Previous</button><span class="note" id="pathway-summary-page"></span><button type="button" id="pathway-summary-next">Next</button></div>
</section>
<section class="card"><h2>Scored component coverage by bin</h2><div class="pathway-key"><span><i style="background:#eff3f4"></i>0%: no accepted panel evidence</span><span><i style="background:#8bbcc1"></i>Partial panel</span><span><i style="background:#176d79"></i>100%: scored panel complete</span><span><i style="background:#ddd7e2"></i>Unavailable / failed / stale</span></div><div class="pathway-scroll pathway-map" id="pathway-heatmap"></div><p class="note">Click a cell to inspect that bin/pathway. Gray cells are unavailable, not zero. By default, all-zero bins are hidden; the full QC table remains available below.</p></section>
<section class="card" id="pathway-components-card"><h2>Component evidence</h2><p class="note" id="pathway-component-status"></p><div class="pathway-scroll"><table><thead><tr><th class="l">Bin</th><th class="l">Pathway / component</th><th class="l">Assay scope</th><th>Scored panel</th><th>Selected route</th><th class="l">Evidence</th><th>Protein copies</th><th class="l">Accepted proteins</th><th class="l">Candidate proteins</th><th class="l">Notes</th></tr></thead><tbody id="pathway-component-body"></tbody></table></div><div class="pathway-actions"><button type="button" id="pathway-component-prev">Previous</button><span class="note" id="pathway-component-page"></span><button type="button" id="pathway-component-next">Next</button><button type="button" id="pathway-export-components">Filtered component CSV</button></div></section>
<section class="card" id="pathway-hits-card"><h2>Supporting HMM hits</h2><p class="note">Accepted model evidence can still identify a shared family rather than an exact gene. The review filter includes these accepted non-specific calls as well as fragments and context candidates.</p><div class="pathway-actions"><label for="pathway-hit-filter">Evidence filter</label><select id="pathway-hit-filter" style="width:auto"><option value="all">All model hits</option><option value="accepted">Accepted component evidence</option><option value="review">Review / non-specific evidence</option></select><button type="button" id="pathway-export-hits">Filtered hits CSV</button><button type="button" id="pathway-export-proteins">Filtered protein FASTA</button></div>
<p class="note" id="pathway-hit-status" role="status"></p><div class="pathway-scroll"><table><thead><tr><th class="l">Bin</th><th class="l">Protein</th><th class="l">Target</th><th class="l">Product / model</th><th class="l">Call specificity</th><th class="l">Status / reason</th><th>HMM coverage</th><th>Sequence bits</th><th>Domain bits</th><th>Independent E-value</th><th class="l">Contig</th><th>Start</th><th>End</th><th class="l">Strand</th></tr></thead><tbody id="pathway-hit-body"></tbody></table></div><div class="pathway-actions"><button type="button" id="pathway-hit-prev">Previous</button><span class="note" id="pathway-hit-page"></span><button type="button" id="pathway-hit-next">Next</button></div></section>
<section class="card"><h2>Gene neighborhoods and reference-operon support</h2><p class="note">Neighborhoods stay within one contig. Operon support is a reference-pattern heuristic from 0 to 1, not a probability or a measure of biochemical completion. Click a neighborhood to inspect gene order and orientation.</p><div class="pathway-scroll"><table><thead><tr><th class="l">Bin</th><th class="l">Pathway</th><th class="l">Neighborhood</th><th class="l">Contig</th><th>Start</th><th>End</th><th>Proteins</th><th class="l">Observed order</th><th class="l">Reference template</th><th>Support</th></tr></thead><tbody id="pathway-neighborhood-body"></tbody></table></div><div id="pathway-neighborhood-figure"></div><p class="note" id="pathway-neighborhood-details"></p><div class="pathway-actions"><button type="button" id="pathway-export-neighborhoods">Filtered neighborhoods CSV</button><button type="button" id="pathway-export-operons">Filtered operon CSV</button></div></section>
<section class="card"><h2>Run quality and original files</h2><p class="note" id="pathway-run-details"></p><div class="pathway-scroll"><table><thead><tr><th class="l">Bin</th><th class="l">Run status</th><th class="l">Input check</th><th>Contigs</th><th>Size (bp)</th><th>Predicted ORFs</th><th>Accepted model hits</th><th>Proteins with accepted hits</th><th>Review hits</th><th class="l">Reason</th></tr></thead><tbody id="pathway-qc-body"></tbody></table></div><p class="note">Original downloads preserve PathwayGenie’s output, including records from any bins flagged as stale above. Bakta locus tags and PathwayGenie protein IDs belong to different gene-calling runs.</p><details><summary>Analysis settings</summary><dl id="pathway-parameters"></dl></details><details style="margin-top:12px"><summary>Download original PathwayGenie files</summary><label for="pathway-file-search">Filter file names</label><input type="search" id="pathway-file-search" placeholder="e.g. accepted_hits, proteins.faa, biotin" style="max-width:500px"><div class="pathway-actions" id="pathway-file-links"></div><p class="note" id="pathway-file-status"></p></details></section>
<dialog class="pathway-dialog" id="pathway-protein-dialog"><div class="row" style="justify-content:space-between"><h2 id="pathway-protein-title"></h2><button type="button" id="pathway-protein-close">Close</button></div><dl id="pathway-protein-fields"></dl><div class="pathway-actions"><button type="button" id="pathway-protein-aa-download">Protein FASTA</button><button type="button" id="pathway-protein-nt-download">Coding-sequence FASTA</button><button type="button" id="pathway-locate-contig">Locate contig in BinaRena</button></div><h2>Predicted protein</h2><pre class="pathway-sequence" id="pathway-protein-aa"></pre><h2>Coding sequence</h2><pre class="pathway-sequence" id="pathway-protein-nt"></pre><p class="note" id="pathway-protein-message"></p></dialog>
'''
