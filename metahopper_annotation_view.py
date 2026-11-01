"""Embed searchable per-bin Bakta features in the offline MetaHopper report."""
from __future__ import annotations

import csv
import html
import json
from pathlib import Path

FIELDS = ['locus_tag', 'gene', 'product', 'type', 'contig', 'start', 'end', 'strand', 'db_xrefs']


def bin_key(rank, name):
    return rank + '/' + name


def bin_button(rank, name):
    return ('<button type="button" class="annotation-bin-link" data-annotation-bin="'
            + html.escape(bin_key(rank, name), quote=True) + '">'
            + html.escape(name) + '</button>')


def collect_features(layout, rows):
    """Only expose features whose source bin passed the report's freshness check."""
    result = []
    for row in rows:
        rank, name = str(row.get('rank', '')), str(row.get('bin', ''))
        if row.get('status') != 'completed':
            continue
        if any(not x or Path(x).name != x or x in ('.', '..') for x in (rank, name)):
            continue
        relative = (layout.base / 'annotation' / rank / name / 'features.tsv').relative_to(layout.root).as_posix()
        item = dict(key=bin_key(rank, name), rank=rank, bin=name, features_path=relative,
                    fasta_path=(layout.base / 'bins' / rank / (name + '.fasta')).relative_to(layout.root).as_posix(),
                    protein_path=(layout.base / 'annotation' / rank / name / 'bakta/annotation.faa').relative_to(layout.root).as_posix(),
                    features=[], error='')
        try:
            path = layout.root / relative
            if not path.resolve().is_relative_to(layout.root.resolve()):
                raise ValueError('Feature file points outside this sample')
            before = path.stat()
            with path.open(newline='', encoding='utf-8') as handle:
                reader = csv.DictReader(handle, delimiter='\t')
                missing = set(FIELDS) - set(reader.fieldnames or [])
                if missing:
                    raise ValueError('Feature table is missing columns: ' + ', '.join(sorted(missing)))
                for number, feature in enumerate(reader, 2):
                    if feature.get('rank') not in (None, '', rank) or feature.get('bin') not in (None, '', name):
                        raise ValueError(f'Feature row {number} belongs to a different bin')
                    item['features'].append([str(feature.get(field) or '') for field in FIELDS])
            after = path.stat()
            if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                raise ValueError('Feature table changed while building the report; retry after annotation finishes')
        except (OSError, ValueError, UnicodeError, csv.Error) as exc:
            item.update(features=[], error=str(exc))
        result.append(item)
    return dict(schema='metahopper-annotation-features-1', run_name=layout.root.name,
                stage=layout.stage, fields=FIELDS, bins=result)


def render_browser(data):
    payload = json.dumps(data, separators=(',', ':'), ensure_ascii=False).replace('<', '\\u003c').replace('&', '\\u0026')
    return BROWSER_HTML + '<script type="application/json" id="annotation-feature-data">' + payload + '</script><script>' + BROWSER_JS + '</script>'


BROWSER_HTML = '''
<div class="card" id="annotation-browser">
  <h2 id="annotation-browser-title">Genes and products</h2>
  <p class="note">Click an annotated bin name in Bins or Annotation to inspect its features. Gene names are shown when assigned; products describe the predicted functions.</p>
  <div class="row" style="margin-bottom:12px">
    <div style="min-width:220px;flex:1"><label for="annotation-bin">Annotated bin</label><select id="annotation-bin"></select></div>
    <div style="min-width:230px;flex:2"><label for="annotation-search">Search genes, products, locus tags or contigs</label><input id="annotation-search" type="search" placeholder="e.g. dnaA, transporter, ribosomal" autocomplete="off"></div>
    <div style="min-width:190px"><label for="annotation-type">Feature type</label><select id="annotation-type"><option value="proteins">Proteins (CDS + sORF)</option><option value="all">All features</option><option value="cds">CDS</option><option value="sorf">sORF</option><option value="trna">tRNA</option><option value="rrna">rRNA</option><option value="ncrna">ncRNA</option></select></div>
    <button type="button" id="annotation-clear">Clear search</button>
  </div>
  <p class="note" id="annotation-feature-status" role="status" aria-live="polite"></p>
  <div class="row" style="margin-bottom:12px"><button type="button" id="annotation-export">Download filtered TSV</button><a id="annotation-feature-download" download>All features TSV</a><a id="annotation-protein-download" download>Proteins FASTA</a><a id="annotation-fasta-download" download>Bin FASTA</a></div>
  <div class="scroll"><table id="annotation-feature-table"><thead><tr><th class="l">Locus tag</th><th class="l">Gene</th><th class="l">Product</th><th class="l">Type</th><th class="l">Contig</th><th>Start</th><th>End</th><th class="l">Strand</th><th class="l">Database references</th></tr></thead><tbody id="annotation-feature-body"></tbody></table></div>
  <div class="row" style="margin-top:12px"><button type="button" id="annotation-prev">Previous</button><span class="note" id="annotation-page"></span><button type="button" id="annotation-next">Next</button></div>
</div>
<style>.annotation-bin-link{border:0;background:transparent;color:var(--acc);padding:0;text-align:left;text-decoration:underline;text-underline-offset:3px;font:inherit;cursor:pointer}.annotation-bin-link:focus-visible{outline:2px solid var(--acc);outline-offset:3px}#annotation-feature-table td:nth-child(3),#annotation-feature-table td:nth-child(9){white-space:normal;min-width:230px;overflow-wrap:anywhere}#annotation-browser .row a{font-size:.84rem}#annotation-feature-status{overflow-wrap:anywhere}</style>
'''


BROWSER_JS = r'''
(function(){
  'use strict';
  var dataNode=document.getElementById('annotation-feature-data');
  if(!dataNode)return;
  var data=JSON.parse(dataNode.textContent), bins=data.bins||[], fields=data.fields, index=Object.create(null);
  fields.forEach(function(field,i){index[field]=i;});
  var el=function(id){return document.getElementById(id);}, select=el('annotation-bin'), search=el('annotation-search'), type=el('annotation-type');
  var current=null, matches=[], page=0, pageSize=100;
  bins.forEach(function(bin){var option=document.createElement('option');option.value=bin.key;option.textContent=bin.rank+' / '+bin.bin;select.appendChild(option);});
  function link(id,path){var node=el(id);node.hidden=!path;if(path){node.setAttribute('data-sample-asset',path);node.setAttribute('href',path.split('/').map(encodeURIComponent).join('/'));}}
  function matchesQuery(row,query,featureType){
    var kind=row[index.type].toLowerCase();
    if(featureType==='proteins'&&kind!=='cds'&&kind!=='sorf')return false;
    if(featureType!=='all'&&featureType!=='proteins'&&kind!==featureType)return false;
    return !query||row.some(function(value){return value.toLocaleLowerCase().indexOf(query)>=0;});
  }
  function render(){
    var body=el('annotation-feature-body');body.textContent='';
    var query=search.value.trim().toLocaleLowerCase();
    matches=current?(current.features||[]).filter(function(row){return matchesQuery(row,query,type.value);}):[];
    var pages=Math.max(1,Math.ceil(matches.length/pageSize));page=Math.min(Math.max(0,page),pages-1);
    matches.slice(page*pageSize,(page+1)*pageSize).forEach(function(row){var tr=document.createElement('tr');fields.forEach(function(field,i){var td=document.createElement('td');td.textContent=row[i]||'—';if(['start','end'].indexOf(field)<0)td.className='l';tr.appendChild(td);});body.appendChild(tr);});
    el('annotation-browser-title').textContent=current?'Genes and products — '+current.bin:'Genes and products';
    var message=!current?'No current completed annotations are available.':current.error?'Gene list unavailable: '+current.error:matches.length.toLocaleString()+' matching features of '+current.features.length.toLocaleString()+' annotated features in '+current.rank+' / '+current.bin+'.';
    if(current&&!current.error&&!matches.length)message+=' No features match the current search and type filter.';
    el('annotation-feature-status').textContent=message;
    el('annotation-page').textContent='Page '+(page+1)+' of '+pages+' · '+pageSize+' rows per page';
    el('annotation-prev').disabled=page===0;el('annotation-next').disabled=page+1>=pages;el('annotation-export').disabled=!matches.length;
    select.disabled=!bins.length;search.disabled=!current;type.disabled=!current;
    link('annotation-feature-download',current&&current.features_path);link('annotation-protein-download',current&&current.protein_path);link('annotation-fasta-download',current&&current.fasta_path);
  }
  function openBin(key,focus){
    var bin=bins.find(function(item){return item.key===key;});if(!bin)return false;
    current=bin;select.value=key;page=0;search.value='';type.value='proteins';render();
    if(focus){var tab=document.querySelector('.tab[data-id="annotation"]');if(tab)tab.click();var panel=el('annotation-browser');if(panel.scrollIntoView)panel.scrollIntoView({block:'start'});}
    return true;
  }
  select.addEventListener('change',function(){openBin(select.value,false);});
  search.addEventListener('input',function(){page=0;render();});type.addEventListener('change',function(){page=0;render();});
  el('annotation-clear').addEventListener('click',function(){search.value='';page=0;render();search.focus();});
  el('annotation-prev').addEventListener('click',function(){page--;render();});el('annotation-next').addEventListener('click',function(){page++;render();});
  el('annotation-export').addEventListener('click',function(){
    if(!current||!matches.length)return;
    function tsvCell(value){value=String(value);return /[\t\r\n"]/.test(value)?'"'+value.replace(/"/g,'""')+'"':value;}
    var lines=[['rank','bin'].concat(fields).map(tsvCell).join('\t')];
    matches.forEach(function(row){lines.push([current.rank,current.bin].concat(row).map(tsvCell).join('\t'));});
    var url=URL.createObjectURL(new Blob([lines.join('\n')+'\n'],{type:'text/tab-separated-values;charset=utf-8'}));
    var anchor=document.createElement('a');anchor.href=url;anchor.download=current.bin.replace(/[^A-Za-z0-9_.-]/g,'_')+'.filtered_features.tsv';document.body.appendChild(anchor);anchor.click();anchor.remove();setTimeout(function(){URL.revokeObjectURL(url);},1000);
  });
  document.addEventListener('click',function(event){var button=event.target.closest&&event.target.closest('[data-annotation-bin]');if(button){event.preventDefault();openBin(button.getAttribute('data-annotation-bin'),true);}});
  if(bins.length)openBin(bins[0].key,false);else render();
})();
'''
