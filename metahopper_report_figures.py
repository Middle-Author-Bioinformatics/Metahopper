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
.figure-host .flow:hover path,.figure-host .flow:focus path{fill-opacity:.55}
.partition-total{font-size:1.2rem;font-weight:650;margin:3px 0 10px}.partition-total span{font-size:.8rem;font-weight:400;color:var(--mut)}
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
@media(max-width:1100px){.report-figures{grid-template-columns:1fr}.report-figures .wide{grid-column:auto}}
@media(max-width:800px){.figure-host{max-height:430px}.figure-host svg{min-width:640px}}
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
    var points=rows.filter(function(r){return number(r.x)!==null&&number(r.y)!==null&&(!options.xlog||Number(r.x)>0)&&(!options.ylog||Number(r.y)>0);});
    if(!points.length){empty(host,options.empty||'No paired values are available for this selection.');return;}
    clear(host);var w=690,h=350,L=74,R=28,T=24,B=62,pw=w-L-R,ph=h-T-B;
    function axis(key,log,minRange){
      var values=points.map(function(p){return Number(p[key]);}),lo=0,hi=Math.max(minRange||1,Math.max.apply(null,values)*1.08);
      if(log){lo=Math.floor(Math.log10(Math.min.apply(null,values)));hi=Math.ceil(Math.log10(Math.max.apply(null,values)));if(lo===hi){lo-=.5;hi+=.5;}}
      return {position:function(v){return ((log?Math.log10(v):v)-lo)/(hi-lo);},tick:function(i){var v=lo+(hi-lo)*i/4;return log?Math.pow(10,v):v;}};
    }
    var xa=axis('x',options.xlog,options.xminRange),ya=axis('y',options.ylog,options.yminRange);
    function tick(v){return Number(v).toLocaleString(undefined,{maximumSignificantDigits:3});}
    var svg=node('svg',{viewBox:'0 0 '+w+' '+h,role:'img','aria-label':options.title});svg.appendChild(node('title',{},options.title));
    for(var i=0;i<=4;i++){var x=L+pw*i/4,y=T+ph-ph*i/4;
      svg.appendChild(node('line',{x1:x,x2:x,y1:T,y2:T+ph,stroke:'#e3eaed'}));svg.appendChild(node('line',{x1:L,x2:L+pw,y1:y,y2:y,stroke:'#e3eaed'}));
      svg.appendChild(node('text',{x:x,y:T+ph+22,'text-anchor':'middle',class:'axis-label'},(options.xformat||tick)(xa.tick(i))));
      svg.appendChild(node('text',{x:L-12,y:y+4,'text-anchor':'end',class:'axis-label'},(options.yformat||tick)(ya.tick(i))));}
    points.forEach(function(p){var dot=node('circle',{cx:L+xa.position(Number(p.x))*pw,cy:T+ph-ya.position(Number(p.y))*ph,r:6,fill:p.color||'#207b86',stroke:'#fff','stroke-width':1.5});dot.appendChild(node('title',{},p.label+'\n'+options.xlabel+': '+tick(p.x)+'\n'+options.ylabel+': '+tick(p.y)));action(dot,options.onPick?function(){options.onPick(p);}:null);svg.appendChild(dot);});
    svg.appendChild(node('text',{x:L+pw/2,y:h-12,'text-anchor':'middle',class:'axis-label'},options.xlabel));svg.appendChild(node('text',{x:17,y:T+ph/2,transform:'rotate(-90 17 '+(T+ph/2)+')','text-anchor':'middle',class:'axis-label'},options.ylabel));host.appendChild(svg);
  }
  function color(name){var colors=['#207b86','#708c58','#9b72aa','#c88340','#4f81ad','#b46f78','#658f92','#a0914b'],hash=0;for(var i=0;i<name.length;i++)hash=(hash*31+name.charCodeAt(i))|0;return colors[(hash>>>0)%colors.length];}
  function flow(host,parts,total,options){
    if(!total||!parts.length){empty(host,'The assembly has no sequence to partition.');return;}
    clear(host);var scale=280/total,y=65,sourceY=65,sourceOffset=0,width=820;
    var positions=parts.map(function(p){var size=p.bp*scale,box=Math.max(32,size+12),pos={p:p,y:y+(box-size)/2,h:size};y+=box;return pos;});
    var height=Math.max(380,y+20),svg=node('svg',{viewBox:'0 0 '+width+' '+height,role:'img','aria-label':'Full assembly partition at the '+options.rank+' rank'});
    svg.appendChild(node('title',{},'Full assembly to '+options.rank+' bins, including unbinned sequence'));
    svg.appendChild(node('desc',{},'Ribbon thickness represents sequence length. Each assembly base is counted once at this rank. Black represents unbinned sequence.'));
    svg.appendChild(node('text',{x:25,y:25},'Full assembly'));
    svg.appendChild(node('text',{x:25,y:44,class:'axis-label'},options.bp(total)+' · '+options.rank+' rank'));
    svg.appendChild(node('rect',{x:25,y:sourceY,width:14,height:280,rx:3,fill:'#465e69'}));
    positions.forEach(function(pos){var p=pos.p,c=p.unbinned?'#000000':color(p.bin),s=sourceY+sourceOffset,t=pos.y,h=pos.h;
      var group=node('g',{class:'flow',tabindex:0});
      group.appendChild(node('title',{},p.bin.replace(/_/g,' ')+' — '+options.bp(p.bp)+' ('+fmt(p.bp/total*100)+'%); '+fmt(p.contigs)+' contigs'));
      if(h>0){group.appendChild(node('path',{d:'M39 '+s+' C205 '+s+' 205 '+t+' 390 '+t+' L390 '+(t+h)+' C205 '+(t+h)+' 205 '+(s+h)+' 39 '+(s+h)+' Z',fill:c,'fill-opacity':p.unbinned ? 1 : .25}));
        group.appendChild(node('rect',{x:390,y:t,width:11,height:Math.max(.7,h),rx:1,fill:c}));}
      var label=p.bin.replace(/_/g,' ');group.appendChild(node('text',{x:414,y:t+h/2-2,style:'font-size:14px'},label.length>47?label.slice(0,45)+'…':label));
      group.appendChild(node('text',{x:414,y:t+h/2+14,class:'axis-label',style:'font-size:12px'},options.bp(p.bp)+' · '+fmt(p.bp/total*100)+'%'));
      if(!p.unbinned)action(group,function(){options.onPick(p);});svg.appendChild(group);sourceOffset+=h;
    });host.appendChild(svg);
  }
  window.MHCharts={bar:bar,scatter:scatter,flow:flow,color:color,empty:empty,number:number,format:fmt};
})();
'''


BIN_FIGURES_HTML = '''
<div id="bin-figures" class="report-figures">
 <section class="card"><h2>Assembly partition</h2><p id="bin-partition-total" class="partition-total"></p><div id="bin-partition-figure" class="figure-host"></div><p class="figure-note" id="bin-partition-note"></p></section>
 <section class="card"><h2>Contiguity and assembly size</h2><p class="note" id="bin-figures-status"></p><div id="bin-contiguity-figure" class="figure-host"></div><p class="figure-note">Both axes use log10 scales. Hover for original values. Bins with missing or nonpositive size or N50 are omitted. N50 is the contig length at which half the bin sequence is in contigs of that length or longer.</p></section>
</div>
'''

BIN_FIGURES_JS = r'''
(function(){
  var source=document.getElementById('bin-figure-data');if(!source||!window.MHCharts)return;
  var data=JSON.parse(source.textContent),partition=JSON.parse(document.getElementById('bin-partition-data').textContent),C=window.MHCharts;
  function bp(v){return v>=1e6?C.format(v/1e6)+' Mb':v>=1e3?C.format(v/1e3)+' kb':C.format(v)+' bp';}
  function refresh(){var rank=document.getElementById('rank-filter'),search=document.getElementById('bin-search');
    function normal(v){return String(v||'').toLowerCase().replace(/[_\s]+/g,' ').trim();}
    var query=normal(search&&search.value),rows=data.filter(function(r){return query?normal(r.bin).indexOf(query)>=0:!rank||r.rank===rank.value;});
    function label(r){return r.rank+' / '+r.bin.replace(/_/g,' ');}
    var host=document.getElementById('bin-partition-figure'),entry=(partition.ranks||{})[rank&&rank.value];
    document.getElementById('bin-partition-total').textContent=partition.total_bp===null?'':bp(partition.total_bp);
    document.getElementById('bin-partition-note').textContent='One rank at a time; the full assembly is retained while searching. Black is unbinned sequence. Hover for sizes and percentages; click a named bin to find it in the table.'+(partition.source?' Source: '+partition.source+'.':'');
    if(partition.error||!entry||entry.error)C.empty(host,partition.error||(entry&&entry.error)||'No verified partition is available at this rank.');
    else C.flow(host,entry.parts,partition.total_bp,{rank:rank.value,bp:bp,onPick:function(p){search.value=p.bin;search.dispatchEvent(new Event('input'));document.getElementById('bin-table').scrollIntoView({block:'start',behavior:'smooth'});}});
    C.scatter(document.getElementById('bin-contiguity-figure'),rows.map(function(r){return {label:label(r),x:r.total_length_bp===null?null:r.total_length_bp/1e6,y:r.N50===null?null:r.N50/1e3,color:C.color(r.bin)};}),{title:'Contiguity and assembly size',xlabel:'Bin size (Mb, log10)',ylabel:'N50 (kb, log10)',xlog:true,ylog:true});
    document.getElementById('bin-figures-status').textContent=rows.length+' bin'+(rows.length===1?'':'s')+' in the current table selection.';
  }
  ['rank-filter','bin-search','bin-search-clear'].forEach(function(id){var n=document.getElementById(id);if(n)['change','input','click'].forEach(function(event){n.addEventListener(event,function(){setTimeout(refresh,0);});});});
  refresh();
})();
'''


def render_bin_figures(rows, partition=None):
    fields = ['rank', 'bin', 'total_length_bp', 'num_contigs', 'N50', 'GC_percent',
              'completeness_percent', 'contamination_percent']
    data = [{field: row.get(field) for field in fields} for row in rows]
    return ('<style>' + FIGURE_CSS + '</style>' + BIN_FIGURES_HTML
            + json_script('bin-figure-data', data)
            + json_script('bin-partition-data', partition or {'total_bp': None, 'ranks': {}})
            + '<script>' + CHART_JS + BIN_FIGURES_JS + '</script>')


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
