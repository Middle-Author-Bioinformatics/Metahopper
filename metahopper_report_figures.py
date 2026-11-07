"""Small, offline SVG figures for MetaHopper bin and annotation summaries."""
import json


def json_script(identifier, data):
    payload = json.dumps(data, separators=(',', ':'), ensure_ascii=False, allow_nan=False, default=str)
    payload = payload.replace('<', '\\u003c').replace('&', '\\u0026')
    return '<script type="application/json" id="' + identifier + '">' + payload + '</script>'


DOWNLOAD_LAYOUT_CSS = '''
/* Keep annotation download controls compact in current and older reports. */
#pane-annotation .scroll{max-width:100%;overflow-x:auto}
#pane-annotation table.sortable{width:max-content!important;min-width:100%;table-layout:auto!important}
#pane-annotation table.sortable th,#pane-annotation table.sortable td{padding:6px 10px;line-height:1.35;vertical-align:middle;white-space:nowrap}
#pane-annotation td.annotation-path{width:auto!important;min-width:370px!important;max-width:none!important;white-space:nowrap!important;word-break:normal!important;overflow-wrap:normal!important}
#pane-annotation td.annotation-path .asset-links{display:inline-flex!important;flex-wrap:nowrap!important;gap:6px;align-items:center;min-width:max-content;white-space:nowrap!important}
#pane-annotation td.annotation-path a{display:inline-block;flex:none;max-width:none!important;white-space:nowrap!important;word-break:normal!important;overflow-wrap:normal!important;border:1px solid var(--line,#dae3e7);border-radius:5px;padding:3px 7px;font-size:.79rem;line-height:1.25;color:var(--acc,#136b78);text-decoration:none;background:var(--card,#fff)}
#pane-annotation td.annotation-path a:hover{background:var(--accs,#e2f0f2)}
#pane-annotation td.annotation-path a:focus-visible{outline:2px solid var(--acc,#136b78);outline-offset:2px}
#pane-annotation table.sortable td.annotation-reason,#pane-annotation table.sortable td:nth-last-child(2){white-space:normal;min-width:90px;max-width:260px}
@media(min-width:900px){
 #pane-annotation table.sortable td.annotation-path{position:sticky;right:0;z-index:1;background:var(--card,#fff);box-shadow:-4px 0 7px rgba(22,34,43,.06)}
 #pane-annotation table.sortable th:last-child{position:sticky;right:0;z-index:3;background:var(--accs,#e2f0f2)}
}
'''


FIGURE_CSS = '''
.report-figures{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px;margin-bottom:16px}
.report-figures .card{margin:0;min-width:0}.report-figures .wide{grid-column:1/-1}
.figure-host{max-height:500px;overflow:auto;min-height:120px}.figure-host svg{width:100%;height:auto;display:block;font-family:inherit}
.figure-host text{fill:var(--ink);font-size:12px}.figure-host .axis-label{fill:var(--mut);font-size:11px}
.figure-host .chart-pick{cursor:pointer}.figure-host .chart-pick:focus{outline:2px solid var(--acc)}
.figure-key{display:flex;gap:12px;flex-wrap:wrap;color:var(--mut);font-size:.8rem;margin:6px 0 12px}
.figure-key i{display:inline-block;width:11px;height:11px;margin-right:5px;border-radius:2px}
.figure-note{margin:7px 0 0;color:var(--mut);font-size:.8rem}.figure-empty{padding:24px 5px;color:var(--mut);font-size:.85rem}
.annotation-summary-scroll{max-height:440px}.annotation-summary-table{width:max-content;min-width:100%;table-layout:auto}
.annotation-summary-table th,.annotation-summary-table td{padding:9px 12px;white-space:nowrap}
.annotation-summary-table .annotation-path{min-width:325px;white-space:nowrap;overflow-wrap:normal}
.annotation-summary-table .annotation-reason{white-space:normal;min-width:90px;max-width:280px}
.asset-links{display:inline-flex;flex-wrap:nowrap;gap:6px;align-items:center;min-width:max-content}
.asset-link{display:inline-flex;align-items:center;justify-content:center;white-space:nowrap!important;word-break:normal!important;overflow-wrap:normal!important;border:1px solid var(--line);border-radius:6px;padding:4px 8px;font-size:.79rem;color:var(--acc);text-decoration:none;background:#fff}
.asset-link:hover{background:var(--accs)}
.asset-link[hidden]{display:none}
@media(min-width:900px){
 .annotation-summary-table tr>:last-child,#bin-table tr>:last-child{position:sticky;right:0;z-index:1;background:var(--card);box-shadow:-4px 0 7px rgba(22,34,43,.06)}
 .annotation-summary-table thead th:last-child,#bin-table thead th:last-child{z-index:3;background:var(--accs)}
}
@media(max-width:800px){.report-figures{grid-template-columns:1fr}.report-figures .wide{grid-column:auto}.figure-host{max-height:430px}.figure-host svg{min-width:640px}}
''' + DOWNLOAD_LAYOUT_CSS

CHART_JS = r'''
(function(){
  'use strict';
  if(window.MHCharts)return;
  var NS='http://www.w3.org/2000/svg';
  function node(tag,attrs,text){var n=document.createElementNS(NS,tag);Object.keys(attrs||{}).forEach(function(k){n.setAttribute(k,String(attrs[k]));});if(text!==undefined)n.textContent=text;return n;}
  function clear(host){host.textContent='';}
  function empty(host,message){clear(host);var n=document.createElement('p');n.className='figure-empty';n.textContent=message;host.appendChild(n);}
  function number(value){return value===null||value===undefined||value===''?null:Number.isFinite(Number(value))?Number(value):null;}
  function fmt(v){return Number(v).toLocaleString(undefined,{maximumFractionDigits:2});}
  function action(n,callback){if(!callback)return;n.classList.add('chart-pick');n.setAttribute('tabindex','0');n.setAttribute('role','button');n.addEventListener('click',callback);n.addEventListener('keydown',function(e){if(e.key==='Enter'||e.key===' '){e.preventDefault();callback();}});}
  function bar(host,rows,options){
    options=options||{};if(!rows.length){empty(host,options.empty||'No values are available for this selection.');return;}
    clear(host);var left=235,right=82,top=35,rowH=30,width=900,plot=width-left-right;
    var max=Math.max.apply(null,rows.map(function(r){return r.values.reduce(function(a,b){return a+Math.max(0,Number(b)||0);},0);}))||1;
    var svg=node('svg',{viewBox:'0 0 '+width+' '+(top+rows.length*rowH+34),role:'img','aria-label':options.title||'Bar chart'});
    svg.appendChild(node('title',{},options.title||'Bar chart'));
    for(var tick=0;tick<=4;tick++){var x=left+plot*tick/4;svg.appendChild(node('line',{x1:x,x2:x,y1:top-5,y2:top+rows.length*rowH,stroke:'#e3eaed'}));svg.appendChild(node('text',{x:x,y:top-13,'text-anchor':'middle',class:'axis-label'},(options.format||fmt)(max*tick/4)));}
    rows.forEach(function(row,i){var y=top+i*rowH,total=0,group=node('g',{});var label=String(row.label);var title=row.title||label;
      var text=node('text',{x:left-12,y:y+18,'text-anchor':'end'},label.length>35?label.slice(0,33)+'…':label);text.appendChild(node('title',{},title));action(text,options.onPick?function(){options.onPick(row,null);}:null);svg.appendChild(text);
      row.values.forEach(function(value,j){value=Math.max(0,Number(value)||0);if(!value)return;var rect=node('rect',{x:left+total/max*plot,y:y+3,width:Math.max(1,value/max*plot),height:21,rx:2,fill:(row.colors||options.colors||['#207b86','#dc9747','#9baab2'])[j%3]});
        rect.appendChild(node('title',{},title+' — '+((options.names||[])[j]||options.unit||'value')+': '+(options.format||fmt)(value)));action(rect,options.onPick?function(){options.onPick(row,j);}:null);group.appendChild(rect);total+=value;});
      svg.appendChild(group);svg.appendChild(node('text',{x:left+plot+9,y:y+18},(options.format||fmt)(total)));
    });svg.appendChild(node('text',{x:left+plot/2,y:top+rows.length*rowH+24,'text-anchor':'middle',class:'axis-label'},options.unit||''));host.appendChild(svg);
  }
  function scatter(host,rows,options){
    var points=rows.filter(function(r){return number(r.x)!==null&&number(r.y)!==null;});
    if(!points.length){empty(host,options.empty||'No paired values are available for this selection.');return;}
    clear(host);var w=690,h=350,L=74,R=28,T=24,B=62,pw=w-L-R,ph=h-T-B;
    var xmax=Math.max(options.xminRange||1,Math.max.apply(null,points.map(function(p){return Number(p.x);}))*1.08),ymax=Math.max(options.yminRange||1,Math.max.apply(null,points.map(function(p){return Number(p.y);}))*1.08);
    var svg=node('svg',{viewBox:'0 0 '+w+' '+h,role:'img','aria-label':options.title});svg.appendChild(node('title',{},options.title));
    for(var i=0;i<=4;i++){var x=L+pw*i/4,y=T+ph-ph*i/4;
      svg.appendChild(node('line',{x1:x,x2:x,y1:T,y2:T+ph,stroke:'#e3eaed'}));svg.appendChild(node('line',{x1:L,x2:L+pw,y1:y,y2:y,stroke:'#e3eaed'}));
      svg.appendChild(node('text',{x:x,y:T+ph+22,'text-anchor':'middle',class:'axis-label'},(options.xformat||fmt)(xmax*i/4)));
      svg.appendChild(node('text',{x:L-12,y:y+4,'text-anchor':'end',class:'axis-label'},(options.yformat||fmt)(ymax*i/4)));}
    points.forEach(function(p){var dot=node('circle',{cx:L+Number(p.x)/xmax*pw,cy:T+ph-Number(p.y)/ymax*ph,r:6,fill:p.color||'#207b86',stroke:'#fff','stroke-width':1.5});dot.appendChild(node('title',{},p.label+'\n'+options.xlabel+': '+fmt(p.x)+'\n'+options.ylabel+': '+fmt(p.y)));action(dot,options.onPick?function(){options.onPick(p);}:null);svg.appendChild(dot);});
    svg.appendChild(node('text',{x:L+pw/2,y:h-12,'text-anchor':'middle',class:'axis-label'},options.xlabel));svg.appendChild(node('text',{x:17,y:T+ph/2,transform:'rotate(-90 17 '+(T+ph/2)+')','text-anchor':'middle',class:'axis-label'},options.ylabel));host.appendChild(svg);
  }
  window.MHCharts={bar:bar,scatter:scatter,empty:empty,number:number,format:fmt};
})();
'''


BIN_FIGURES_HTML = '''
<div id="bin-figures" class="report-figures">
 <section class="card wide"><h2>Assembly size by bin</h2><p class="note" id="bin-figures-status"></p><div id="bin-size-figure" class="figure-host"></div><p class="figure-note">Figures follow the rank and taxonomic search above. Bins at different ranks overlap; their sizes should not be added together.</p></section>
 <section class="card"><h2>Completeness and contamination</h2><div id="bin-quality-figure" class="figure-host"></div><p class="figure-note">Hover over a point for its bin and exact values. Missing CheckM results are omitted. Taxonomic bins can contain multiple genomes.</p></section>
 <section class="card"><h2>Contiguity and assembly size</h2><div id="bin-contiguity-figure" class="figure-host"></div><p class="figure-note">N50 is the sequence length at which half the assembly lies in contigs of that length or longer.</p></section>
</div>
'''

BIN_FIGURES_JS = r'''
(function(){
  var source=document.getElementById('bin-figure-data');if(!source||!window.MHCharts)return;
  var data=JSON.parse(source.textContent),C=window.MHCharts;
  function refresh(){var rank=document.getElementById('rank-filter'),search=document.getElementById('bin-search');
    function normal(v){return String(v||'').toLowerCase().replace(/[_\s]+/g,' ').trim();}
    var query=normal(search&&search.value),rows=data.filter(function(r){return query?normal(r.bin).indexOf(query)>=0:!rank||r.rank===rank.value;});
    function label(r){return r.rank+' / '+r.bin.replace(/_/g,' ');}
    C.bar(document.getElementById('bin-size-figure'),rows.filter(function(r){return C.number(r.total_length_bp)!==null;}).sort(function(a,b){return b.total_length_bp-a.total_length_bp;}).map(function(r){return {label:label(r),values:[r.total_length_bp/1e6]};}),{title:'Assembly size by bin',unit:'Assembly size (Mb)',format:function(v){return Number(v).toLocaleString(undefined,{maximumFractionDigits:3});}});
    C.scatter(document.getElementById('bin-quality-figure'),rows.map(function(r){return {label:label(r),x:r.contamination_percent,y:r.completeness_percent};}),{title:'Completeness and contamination',xlabel:'Contamination (%)',ylabel:'Completeness (%)',yminRange:100,empty:'No CheckM completeness and contamination values are available for these bins.'});
    C.scatter(document.getElementById('bin-contiguity-figure'),rows.map(function(r){return {label:label(r),x:r.total_length_bp===null?null:r.total_length_bp/1e6,y:r.N50===null?null:r.N50/1e3};}),{title:'Contiguity and assembly size',xlabel:'Assembly size (Mb)',ylabel:'N50 (kb)'});
    document.getElementById('bin-figures-status').textContent=rows.length+' bin'+(rows.length===1?'':'s')+' in the current table selection.';
  }
  ['rank-filter','bin-search','bin-search-clear'].forEach(function(id){var n=document.getElementById(id);if(n)['change','input','click'].forEach(function(event){n.addEventListener(event,function(){setTimeout(refresh,0);});});});
  refresh();
})();
'''


def render_bin_figures(rows):
    fields = ['rank', 'bin', 'total_length_bp', 'num_contigs', 'N50', 'GC_percent',
              'completeness_percent', 'contamination_percent']
    data = [{field: row.get(field) for field in fields} for row in rows]
    return ('<style>' + FIGURE_CSS + '</style>' + BIN_FIGURES_HTML
            + json_script('bin-figure-data', data) + '<script>' + CHART_JS + BIN_FIGURES_JS + '</script>')


def render_annotation_figures(rows):
    data = []
    for row in rows:
        if row.get('status') != 'completed':
            continue
        values = []
        for field in ('named_proteins', 'hypothetical_proteins'):
            try:
                values.append(max(0, int(row.get(field) or 0)))
            except (ValueError, TypeError):
                values.append(0)
        data.append(dict(label=str(row.get('bin', '')).replace('_', ' '),
                         key=str(row.get('rank', '')) + '/' + str(row.get('bin', '')), values=values))
    return ('<div class="card"><h2>Protein product assignments by bin</h2>'
            '<div class="figure-key"><span><i style="background:#207b86"></i>Non-hypothetical products</span>'
            '<span><i style="background:#dc9747"></i>Hypothetical proteins</span></div>'
            '<div id="annotation-overview-figure" class="figure-host"></div>'
            '<p class="figure-note">Counts come from completed Bakta summaries. Click a bin or a bar to open its feature explorer. A product assignment is a prediction, not experimental validation.</p></div>'
            + json_script('annotation-overview-data', data)
            + '<script>' + CHART_JS + r'''
(function(){var host=document.getElementById('annotation-overview-figure'),rows=JSON.parse(document.getElementById('annotation-overview-data').textContent);window.MHCharts.bar(host,rows,{title:'Protein product assignments by bin',unit:'Predicted protein features',names:['non-hypothetical products','hypothetical proteins'],onPick:function(row,segment){document.dispatchEvent(new CustomEvent('metahopper:annotation-bin',{detail:{key:row.key,category:segment===1?'hypothetical':'all'}}));}});})();
''' + '</script>')
