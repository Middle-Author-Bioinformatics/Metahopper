/* One report's browsing state only; manual bin assignments are never reset. */
(function(){
 'use strict';
 if(typeof document==='undefined'||window.MHReportControls)return;
 function init(){
  if(!document.querySelector||!document.createElement||!document.getElementById)return;
  if(window.MHReportControls||document.getElementById('mh-report-toolbar'))return;
  var host=document.querySelector('.wrap');if(!host||typeof host.insertBefore!=='function')return;
  var E=id=>document.getElementById(id),Q=selector=>document.querySelector(selector);
  function json(id,fallback){try{return JSON.parse(E(id).textContent);}catch(_){return fallback;}}
  function node(tag,value){var n=document.createElement(tag);if(value!==undefined)n.textContent=value;return n;}
  function event(id,type){var element=E(id);if(element)element.dispatchEvent(new Event(type||'change',{bubbles:true}));}
  var bar=node('nav');bar.id='mh-report-toolbar';bar.className='mh-report-toolbar';bar.setAttribute('aria-label','Report navigation');
  var back=node('button','← Back');back.id='mh-back-view';back.type='button';back.disabled=true;
  var reset=node('button','Reset filters');reset.id='mh-reset-view';reset.type='button';
  bar.appendChild(back);bar.appendChild(reset);
  var downloads=node('div');downloads.className='mh-run-downloads';bar.appendChild(downloads);
  host.insertBefore(bar,E('tabs')||host.firstChild);
  var manifest=json('metahopper-assets',null),files=manifest&&manifest.files||{},support=json('metahopper-run-downloads',null);
  [['Run info',support?support.info:'report_support/run_info.json'],['Full log',support?support.log:'report_support/full_run.log']].forEach(function(pair){
   if(!pair[1])return;
   if(!(manifest?files[pair[1]]:support&&Object.values(support).includes(pair[1])))return;
   var a=node('a',pair[0]);a.href=pair[1];a.setAttribute('download','');a.setAttribute('data-sample-asset',pair[1]);downloads.appendChild(a);
  });
  // Consolidate the repeated feature export row without removing its functionality.
  var exporters=['annotation-export','annotation-export-proteins','annotation-export-genes'].filter(id=>E(id));
  if(exporters.length){var menu=node('details');menu.id='mh-feature-exports';menu.hidden=true;menu.appendChild(node('summary','Download features ▾'));var list=node('div');
   exporters.forEach(function(id){var b=node('button',{'annotation-export':'Features TSV','annotation-export-proteins':'Proteins FASTA','annotation-export-genes':'Genes FASTA'}[id]);b.type='button';b.dataset.exportId=id;b.addEventListener('click',function(){E(id).click();menu.open=false;});list.appendChild(b);});menu.appendChild(list);downloads.appendChild(menu);}
  var tooltip=node('div');tooltip.id='mh-bin-tooltip';tooltip.hidden=true;tooltip.setAttribute('role','tooltip');document.body.appendChild(tooltip);
  function hideTip(){tooltip.hidden=true;}
  function format(v,d){return Number(v).toLocaleString(undefined,{maximumFractionDigits:d===undefined?2:d});}
  function bp(v){return v>=1e6?format(v/1e6)+' Mb':v>=1e3?format(v/1e3)+' kb':format(v,0)+' bp';}
  function known(v){return v!==null&&v!==undefined&&v!==''&&Number.isFinite(Number(v));}
  function showTip(row,point){
   tooltip.textContent='';tooltip.appendChild(node('strong',String(row.bin).replace(/_/g,' ')));
   var facts=[['Rank',row.rank],['Bin size',known(row.total_length_bp)?bp(row.total_length_bp):null],['N50',known(row.N50)?bp(row.N50):null],['Contigs',known(row.num_contigs)?format(row.num_contigs,0):null],['GC',known(row.GC_percent)?format(row.GC_percent)+'%':null],['Completeness',known(row.completeness_percent)?format(row.completeness_percent)+'%':null],['Contamination',known(row.contamination_percent)?format(row.contamination_percent)+'%':null]];
   var dl=node('dl');facts.forEach(function(pair){if(pair[1]!==null){dl.appendChild(node('dt',pair[0]));dl.appendChild(node('dd',pair[1]));}});tooltip.appendChild(dl);tooltip.hidden=false;
   var rect=point.getBoundingClientRect(),x=rect.right+12,y=rect.top;
   tooltip.style.left=Math.max(8,Math.min(x,window.innerWidth-tooltip.offsetWidth-8))+'px';tooltip.style.top=Math.max(8,Math.min(y,window.innerHeight-tooltip.offsetHeight-8))+'px';
  }
  function bindPoints(){
   var data=json('bin-figure-data',[]),figure=E('bin-contiguity-figure');if(!figure)return;
   figure.querySelectorAll('circle').forEach(function(point){
    if(point.dataset.mhTooltip)return;
    var title=point.querySelector('title'),name=title&&title.textContent.split('\n')[0];
    var row=data.find(r=>r.rank+' / '+String(r.bin).replace(/_/g,' ')===name);if(!row)return;
    point.dataset.mhTooltip='true';if(title)title.remove();point.setAttribute('tabindex','0');point.setAttribute('role','button');
    point.setAttribute('aria-label',row.bin+', '+row.rank+', bin size '+bp(row.total_length_bp)+', N50 '+bp(row.N50));point.setAttribute('aria-describedby','mh-bin-tooltip');
    point.addEventListener('pointerenter',()=>showTip(row,point));point.addEventListener('focus',()=>showTip(row,point));point.addEventListener('pointerleave',hideTip);point.addEventListener('blur',hideTip);
    function pick(){hideTip();if(E('bin-search')){E('bin-search').value=row.bin;event('bin-search','input');E('bin-table')?.scrollIntoView({block:'start',behavior:'smooth'});}}
    point.addEventListener('click',pick);point.addEventListener('keydown',function(e){if(e.key==='Escape')hideTip();if(e.key==='Enter'||e.key===' '){e.preventDefault();pick();schedule(e);}});
   });
  }
  var ids=['rank-filter','bin-search','x','y','xs','ys','c','bf','ml','q','annotation-bin','annotation-search','annotation-type','annotation-category','pathway-bin','pathway-pathway','pathway-search','pathway-hit-filter','pathway-hide-zero'];
  function state(){
   var values={};ids.forEach(function(id){var el=E(id);if(el)values[id]=el.type==='checkbox'?el.checked:el.value;});
   var tab=Q('.tab.on'),sub=Q('[role="tab"][aria-selected="true"][id^="annotation-"]');
   return {values:values,tab:tab&&tab.dataset.id||'bins',sub:sub&&sub.id||'annotation-genes-tab',selection:window.MHBinarenaView?window.MHBinarenaView.capture():[],scroll:window.scrollY||0};
  }
  function identity(s){return JSON.stringify([s.values,s.tab,s.sub,s.selection]);}
  var defaults=state(),current=defaults,history=[],restoring=false,queued=false,lastEdit='',lastEditTime=0;
  function update(){back.disabled=!history.length;var menu=E('mh-feature-exports');if(menu)menu.hidden=!(Q('#pane-annotation.on')&&E('annotation-genes-section')?.hidden!==true);}
  function commit(e){
   if(restoring)return;
   var next=state();if(identity(next)===identity(current)){compact();return;}
   if(next.values['rank-filter']!==current.values['rank-filter'])history=[];
   else {var edit=e&&e.type==='input'&&e.target&&e.target.id||'',now=Date.now();if(!edit||edit!==lastEdit||now-lastEditTime>900)history.push(current);lastEdit=edit;lastEditTime=now;if(history.length>25)history.shift();}
   current=next;update();compact();
  }
  function schedule(e){if(restoring||queued)return;queued=true;queueMicrotask(function(){queued=false;commit(e);});}
  function restore(s){
   restoring=true;hideTip();
   if(E('annotation-bin')&&s.values['annotation-bin']!==undefined){E('annotation-bin').value=s.values['annotation-bin'];event('annotation-bin');}
   Object.keys(s.values).forEach(function(id){if(id==='rank-filter')return;var el=E(id);if(el){if(el.type==='checkbox')el.checked=s.values[id];else el.value=s.values[id];}});
   event('bin-search','input');event('bf');event('q','input');event('annotation-search','input');event('pathway-search','input');
   Q('.tab[data-id="'+s.tab+'"]')?.click();E(s.sub)?.click();
   if(window.MHBinarenaView)window.MHBinarenaView.restore(s.selection);
   document.querySelectorAll('.mh-search-target').forEach(el=>el.classList.remove('mh-search-target'));
   current=state();restoring=false;lastEdit='';update();compact();window.scrollTo({top:s.scroll,behavior:'instant'});
  }
  back.addEventListener('click',function(){if(!history.length)return;restore(history.pop());document.dispatchEvent(new CustomEvent('metahopper:view-restored'));});
  reset.addEventListener('click',function(){
   var before=state();if(identity(before)!==identity(current))current=before;
   restoring=true;hideTip();
   if(E('bin-search')){E('bin-search').value='';event('bin-search','input');}
   ['x','y','xs','ys','c','bf','ml','q'].forEach(function(id){if(E(id)&&defaults.values[id]!==undefined)E(id).value=defaults.values[id];});
   event('bf');event('q','input');if(window.MHBinarenaView)window.MHBinarenaView.reset();else E('clr')?.click();
   if(E('annotation-bin')){E('annotation-bin').value='all';event('annotation-bin');}E('annotation-clear')?.click();E('pathway-reset')?.click();
   document.querySelectorAll('.mh-search-target').forEach(el=>el.classList.remove('mh-search-target'));
   var after=state();if(identity(before)!==identity(after))history.push(before);current=after;restoring=false;lastEdit='';update();compact();
   document.dispatchEvent(new CustomEvent('metahopper:reset-view'));
  });
  function compact(){
   var table=E('annotation-summary-table');if(table){var headers=Array.from(table.querySelectorAll('thead th'));headers.forEach(function(th,i){if(/^(hypothetical proteins|non-hypothetical proteins)$/i.test(th.textContent.trim()))table.querySelectorAll('tr').forEach(tr=>{if(tr.children[i])tr.children[i].classList.add('mh-redundant-column');});});}
   var status=E('annotation-feature-status');if(status)status.classList.toggle('mh-routine-status',/^[\d,]+ \/ [\d,]+ features$/.test(status.textContent));
   var parameters=E('pathway-parameters');if(parameters)parameters.closest('details')?.setAttribute('hidden','');
   ['annotation','pathway-summary','pathway-hit'].forEach(function(id){var prev=E(id+'-prev'),next=E(id+'-next');if(prev&&next&&prev.parentElement===next.parentElement)prev.parentElement.hidden=prev.disabled&&next.disabled;});
   document.querySelectorAll('[data-export-id]').forEach(b=>b.disabled=E(b.dataset.exportId).disabled);update();bindPoints();
  }
  ['click','input','change','pointerup','keyup'].forEach(type=>document.addEventListener(type,function(e){if(restoring||e.target?.closest?.('#mh-report-toolbar'))return;if(!queued)current.scroll=window.scrollY||0;schedule(e);},true));
  document.addEventListener('metahopper:bins-filter',function(){hideTip();bindPoints();schedule();});
  document.addEventListener('keydown',e=>{if(e.key==='Escape')hideTip();});
  if(typeof MutationObserver!=='undefined'){
   ['bin-contiguity-figure','annotation-feature-body','pathway-summary-body','pathway-hit-body'].forEach(function(id){var target=E(id);if(target)new MutationObserver(function(){hideTip();compact();}).observe(target,{childList:true});});
  }
  window.MHReportControls={reset:()=>reset.click(),back:()=>back.click(),capture:state};
  compact();
 }
 if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',init);else init();
})();
