/* Inlined report explorer; no remote scripts or external UI dependencies. */
(function(){
 'use strict';
 var source=document.getElementById('metahopper-pathways-data');if(!source)return;
 var D=JSON.parse(source.textContent),E=function(id){return document.getElementById(id);},bins=D.bins||[],paths=D.pathways||[],pages={summary:0,hit:0},sorts={},pageSize=100,current={},selectedProtein=null;
 var bin=E('pathway-bin'),path=E('pathway-pathway'),search=E('pathway-search'),hide=E('pathway-hide-zero');
 var labels={complete_target_panel:'All panel targets detected',partial_target_panel:'Some targets detected',not_detected_in_panel:'No panel genes detected',supported_selected_route:'Selected route supported',unresolved_unmeasured_steps:'Unmeasured route steps',partial_evidence:'Partial route evidence',bin_failed:'Bin failed',stale_input:'Stale input',missing_bin:'Current bin missing',unverified:'Input not verified',current:'Input verified',preview_unverified:'Preview: input not checked',detected:'Detected',candidate_only:'Candidate only',not_detected:'Not detected',not_assayed:'Not assayed',fusion_only_not_detected:'Fusion-only: no accepted evidence',accepted:'Accepted',context_candidate:'Context candidate',fragment_candidate:'Fragment candidate',shared_family:'Shared family',domain_only:'Domain only',specific:'Specific target',ambiguous:'Ambiguous',competing_hit:'Competing hit',targeted:'Targeted',fusion_only:'Selected fusions only',ok:'Completed',no_predicted_genes:'No predicted genes',failed:'Failed'};
 Object.assign(labels,{no_hits:'No hits',no_accepted_hit:'No accepted hit',fusion_only_not_detected:'No fusion hit'});
 function label(v){return labels[v]||String(v??'').replace(/_/g,' ');}
 function n(v){return v===null||v===undefined||v===''||!Number.isFinite(Number(v))?null:Number(v);}
 function number(v,digits){var x=n(v);return x===null?'NA':x.toLocaleString(undefined,{maximumFractionDigits:digits??3});}
 function text(tag,value,cls){var e=document.createElement(tag);if(value!==undefined)e.textContent=String(value);if(cls)e.className=cls;return e;}
 function key(b,p){return b+'\x1f'+p;}
 function split(v){return String(v||'').split(';').filter(Boolean);}
 function normalized(v){return String(v??'').toLowerCase().replace(/[_\s]+/g,' ').trim();}
 function includes(row,q){return !q||Object.values(row).some(function(v){return typeof v!=='object'&&normalized(v).includes(normalized(q));});}
 function option(select,value,title){var o=text('option',title);o.value=value;select.appendChild(o);}
 function button(title,callback){var b=text('button',title,'pathway-text-button');b.type='button';b.addEventListener('click',callback);return b;}
 function badge(value){var good=['accepted','detected','complete_target_panel','supported_selected_route','current','ok'].includes(value),bad=['bin_failed','stale_input','missing_bin','failed'].includes(value),missing=['no_hits','no_accepted_hit','fusion_only_not_detected'].includes(value),warn=['partial_target_panel','unresolved_unmeasured_steps','partial_evidence','candidate_only','context_candidate','fragment_candidate','shared_family','domain_only','ambiguous','competing_hit','unverified','preview_unverified'].includes(value);return text('span',label(value),'pathway-tag '+(good?'ok':bad?'bad':missing?'missing':warn?'partial':''));}
 function cell(tr,value,cls){var td=text('td',undefined,cls||'l');if(value&&typeof value==='object'&&value.nodeType)td.appendChild(value);else td.textContent=value===null||value===undefined||value===''?'—':String(value);tr.appendChild(td);return td;}
 function resetPages(){Object.keys(pages).forEach(function(k){pages[k]=0;});}
 function pick(b,p){bin.value=b;path.value=p||'all';search.value='';resetPages();refresh();E('pathway-hits-card').scrollIntoView({block:'start'});}
 function showProtein(b,p){selectedProtein={bin:b,protein:p};var record=D.proteins[key(b,p)]||{},hits=D.hits.filter(function(h){return h.bin_id===b&&h.protein===p;}),first=hits[0]||{};
  E('pathway-protein-title').textContent=b+' / '+p;var dl=E('pathway-protein-fields');dl.textContent='';
  [['Gene caller','Prodigal (PathwayGenie); IDs are independent of Bakta locus tags'],['Targets',Array.from(new Set(hits.map(function(h){return h.target_gene;}))).join('; ')],['Products',Array.from(new Set(hits.map(function(h){return h.source_product;}))).join('; ')],['Models',hits.map(function(h){return h.model_accession+' — '+label(h.status)+'; '+label(h.target_call);}).join('\n')],['Contig',record.original_contig||first.contig],['Coordinates',(record.start||first.start||'?')+'–'+(record.end||first.end||'?')+' ('+((record.strand||first.strand)==='-1'?'-':'+')+')'],['Partial-gene flag',record.partial??first.partial],['Protein length',record.length||first.protein_length],['HMM coverage',hits.map(function(h){return h.model_accession+': '+number(n(h.hmm_coverage)===null?null:n(h.hmm_coverage)*100,2)+'%';}).join('; ')],['Model / protein intervals',hits.map(function(h){return h.model_accession+': HMM '+h.hmm_intervals+'; protein '+h.protein_intervals;}).join('\n')]].forEach(function(pair){dl.appendChild(text('dt',pair[0]));dl.appendChild(text('dd',pair[1]||'—'));});
  E('pathway-protein-aa').textContent=record.aa?record.aa.match(/.{1,70}/g).join('\n'):'No protein sequence embedded for this record.';
  E('pathway-protein-nt').textContent=record.nt?record.nt.match(/.{1,70}/g).join('\n'):'No coding sequence embedded for this record.';
  E('pathway-protein-aa-download').disabled=!record.aa;E('pathway-protein-nt-download').disabled=!record.nt;
  var contig=record.original_contig||first.contig,scatter=window.DATA&&window.DATA.scatter;E('pathway-locate-contig').disabled=!(scatter&&scatter.ids&&scatter.ids.includes(contig));
  E('pathway-protein-message').textContent=E('pathway-locate-contig').disabled?'This contig is not in the report’s plotted contig set.':'';
  E('pathway-protein-dialog').showModal();
 }
 function paginated(name,rows,columns){
  var body=E('pathway-'+name+'-body'),table=body.closest('table');
  var sort=sorts[name];if(sort)rows=rows.slice().sort(function(a,b){var x=columns[sort.column].value(a),y=columns[sort.column].value(b);if(x===null||x===undefined)return y===null||y===undefined?0:1;if(y===null||y===undefined)return -1;return (typeof x==='number'&&typeof y==='number'?x-y:String(x).localeCompare(String(y),undefined,{numeric:true}))*sort.direction;});
  current[name]=rows;var totalPages=Math.max(1,Math.ceil(rows.length/pageSize));pages[name]=Math.max(0,Math.min(pages[name],totalPages-1));body.textContent='';
  rows.slice(pages[name]*pageSize,(pages[name]+1)*pageSize).forEach(function(row){var tr=text('tr',undefined,row.missing?'pathway-missing':'');columns.forEach(function(c){cell(tr,c.render?c.render(row):c.value(row),c.cls||'l');});body.appendChild(tr);});
  E('pathway-'+name+'-page').textContent=number(rows.length,0)+' rows · page '+(pages[name]+1)+' of '+totalPages;E('pathway-'+name+'-prev').disabled=pages[name]===0;E('pathway-'+name+'-next').disabled=pages[name]+1>=totalPages;
  Array.from(table.querySelectorAll('thead th')).forEach(function(th,i){if(th.dataset.pathwaySort)return;th.dataset.pathwaySort='1';th.tabIndex=0;th.setAttribute('role','button');th.title='Sort by '+th.textContent;function sortRows(){sorts[name]={column:i,direction:sorts[name]&&sorts[name].column===i?-sorts[name].direction:1};pages[name]=0;refresh();}th.addEventListener('click',sortRows);th.addEventListener('keydown',function(e){if(e.key==='Enter'||e.key===' '){e.preventDefault();sortRows();}});});
 }
 var definitionsByPath=Object.create(null);D.definitions.forEach(function(d){(definitionsByPath[d.pathway]||(definitionsByPath[d.pathway]=[])).push(d);});
 function hitInPath(hit,pid){if(hit.pathway)return hit.pathway===pid;var p=paths.find(p=>p.id===pid);return (definitionsByPath[pid]||[]).some(d=>split(hit.components).includes(d.component))||!!(p&&(p.target_genes||[]).includes(hit.target_gene));}
 option(bin,'all','All analyzed bins');bins.forEach(function(b){option(bin,b.bin_id,b.bin_id.replace(/_/g,' '));});
 option(path,'all','All pathways');paths.forEach(function(p){option(path,p.id,(p.label||p.id).replace(/_/g,' '));});
 (D.warnings||[]).forEach(function(w){E('pathway-warnings').appendChild(text('p',w,'pathway-warning'));});
 var summaryColumns=[
  {value:r=>r.bin_id,render:r=>button(r.bin_id,function(){pick(r.bin_id,path.value);} )},
  {value:r=>r.pathway,render:r=>button(r.pathway.replace(/_/g,' '),function(){pick(r.bin_id,r.pathway);})},
  {value:r=>n(r.panel_percent),render:r=>n(r.panel_percent)===null?'NA':number(r.panel_percent)+'%'},
  {value:r=>n(r.detected_scored_components),render:r=>n(r.detected_scored_components)===null?'NA':r.detected_scored_components+' / '+r.total_scored_components},
  {value:r=>r.panel_status,render:r=>badge(r.panel_status)},
  {value:r=>r.pathway,render:r=>button('Panel genes',function(){showPanel(r.pathway);})}];
 function hitPaths(hit){return paths.filter(function(p){return hitInPath(hit,p.id);});}
 function matchStatus(hit){return hit.status==='accepted'?(hit.target_call==='specific'?'accepted':hit.target_call):hit.status;}
 function panelRows(){return (D.panel_genes||D.definitions.map(function(d){return {pathway:d.pathway,component:d.component,gene:d.label||d.component,models:'',counted_in_coverage:d.in_heatmap_denominator==='1',measured:d.assay_scope!=='not_assayed',notes:d.notes||''};}));}
 // Missing rows come from verified component evidence, never from the current
 // hit filter. Unmeasured steps and failed/stale bins cannot become negatives.
 var missingRows=[];
 (D.component_states||[]).forEach(function(e){
  if(!['not_detected','candidate_only','fusion_only_not_detected'].includes(e.state))return;
  panelRows().filter(p=>p.component===e.component&&p.measured).forEach(function(p){
   missingRows.push({bin_id:e.bin_id,pathway:p.pathway,components:e.component,target_gene:p.gene,source_product:'',protein:'',model_accession:p.models,
    status:e.state==='candidate_only'?'candidate_only':e.has_hits?'no_accepted_hit':e.state==='fusion_only_not_detected'?'fusion_only_not_detected':'no_hits',
    independent_evalue:null,hmm_coverage:null,missing:true});
  });
 });
 function panelNote(row){return (row.measured?(row.counted_in_coverage?'Counts toward coverage.':'Supporting component.'):'Not measured by this panel.')+(row.notes?' '+row.notes:'');}
 function showPanel(pid){var p=paths.find(function(p){return p.id===pid;});E('pathway-panel-title').textContent=(p?p.label||p.id:pid)+' — panel genes';E('pathway-panel-note').textContent=p&&p.notes||'Genes and models configured for this pathway.';var body=E('pathway-panel-body');body.textContent='';panelRows().filter(function(r){return r.pathway===pid;}).forEach(function(r){var tr=text('tr');cell(tr,r.gene);cell(tr,r.models||(r.measured?'Models not recorded':'—'));cell(tr,panelNote(r),'l wrap-text');body.appendChild(tr);});E('pathway-panel-dialog').showModal();}
 var hitColumns=[{value:r=>r.bin_id},
  {value:r=>hitPaths(r).map(p=>p.id).join('; '),render:r=>{var span=text('span');hitPaths(r).forEach(function(p,i){if(i)span.appendChild(document.createTextNode('; '));span.appendChild(button(p.label||p.id,function(){showPanel(p.id);}));});return span;}},
  {value:r=>r.target_gene,render:r=>{var span=text('span',r.target_gene||'—');if(r.source_product)span.appendChild(text('div',r.source_product,'note'));return span;},cls:'l wrap-text'},
  {value:r=>r.protein,render:r=>r.protein?button(r.protein,function(){showProtein(r.bin_id,r.protein);}):'—'},
  {value:r=>r.model_accession,cls:'l wrap-text'},
  {value:r=>matchStatus(r),render:r=>badge(matchStatus(r))},
  {value:r=>n(r.independent_evalue),render:r=>r.independent_evalue},
  {value:r=>n(r.hmm_coverage),render:r=>n(r.hmm_coverage)===null?'—':number(n(r.hmm_coverage)*100,2)+'%'}];
 function heatmap(rows){var host=E('pathway-heatmap');host.textContent='';if(!rows.length){host.appendChild(text('p','No matching pathway records.','note'));return;}
  var shownBins=Array.from(new Set(rows.map(r=>r.bin_id))),shownPaths=paths.filter(p=>rows.some(r=>r.pathway===p.id)),lookup=Object.create(null);rows.forEach(r=>{lookup[key(r.bin_id,r.pathway)]=r;});
  var table=text('table'),head=text('thead'),tr=text('tr');tr.appendChild(text('th','Bin'));shownPaths.forEach(p=>tr.appendChild(text('th',(p.label||p.id).replace(/_/g,' '))));head.appendChild(tr);table.appendChild(head);var body=text('tbody');
  shownBins.forEach(function(b){var tr=text('tr');tr.appendChild(text('th',b));shownPaths.forEach(function(p){var row=lookup[key(b,p.id)],value=row?n(row.panel_percent):null,td=text('td'),control=text('button',value===null?'NA':number(value)+'%');control.type='button';
   var f=(value||0)/100;control.style.background=value===null?'#ddd7e2':value===0?'#eff3f4':'rgb('+Math.round(228-205*f)+','+Math.round(241-132*f)+','+Math.round(241-120*f)+')';control.style.color=value!==null&&value>=65?'#fff':'#163b43';
   control.title=b+' / '+p.id+' — '+(value===null?'Unavailable':number(value)+'% measured genes detected')+'; '+(row?label(row.route_status):'No record');control.setAttribute('aria-label',control.title);control.addEventListener('click',function(){pick(b,p.id);});td.appendChild(control);tr.appendChild(td);});body.appendChild(tr);});table.appendChild(body);host.appendChild(table);
 }
 function refresh(){
  var q=search.value.trim().toLowerCase(),chosenPaths=new Set(paths.filter(p=>path.value==='all'||p.id===path.value).map(p=>p.id));
  var selectedSummaries=D.summary.filter(r=>(bin.value==='all'||r.bin_id===bin.value)&&chosenPaths.has(r.pathway));
  var visibleBins=new Set(selectedSummaries.map(r=>r.bin_id)),hiddenZeros=0;
  if(hide.checked&&bin.value==='all')Array.from(visibleBins).forEach(function(b){var rows=selectedSummaries.filter(r=>r.bin_id===b);if(rows.length&&rows.every(r=>n(r.panel_percent)===0)){visibleBins.delete(b);hiddenZeros++;}});
  var queryHits=D.hits.concat(missingRows).filter(h=>includes(h,q)),matchingPairs=new Set();if(q)queryHits.forEach(function(h){chosenPaths.forEach(function(pid){if(hitInPath(h,pid))matchingPairs.add(key(h.bin_id,pid));});});
  var summaryRows=selectedSummaries.filter(r=>visibleBins.has(r.bin_id)&&(!q||includes(r,q)||matchingPairs.has(key(r.bin_id,r.pathway))));
  var activeBins=new Set(summaryRows.map(r=>r.bin_id));
  E('pathway-selection-status').textContent=number(activeBins.size,0)+' bins · '+number(new Set(summaryRows.map(r=>r.pathway)).size,0)+' pathways'+(hiddenZeros?' · '+hiddenZeros+' zero-coverage bins hidden':'');
  paginated('summary',summaryRows,summaryColumns);heatmap(summaryRows);
  var filter=E('pathway-hit-filter').value;
  function selected(h){return activeBins.has(h.bin_id)&&(path.value==='all'||hitInPath(h,path.value))&&(!q||includes(h,q)||hitPaths(h).some(p=>includes(p,q)));}
  var hits=filter==='missing'?[]:D.hits.filter(h=>selected(h)&&(!['accepted','overview'].includes(filter)||h.status==='accepted')&&(filter!=='review'||h.status!=='accepted'||h.target_call!=='specific'));
  var missing=['overview','missing'].includes(filter)?missingRows.filter(selected):[];
  var evidence=hits.concat(missing);if(filter==='overview')evidence.sort((a,b)=>a.bin_id.localeCompare(b.bin_id)||hitPaths(a).map(p=>p.id).join(';').localeCompare(hitPaths(b).map(p=>p.id).join(';'))||a.target_gene.localeCompare(b.target_gene));
  paginated('hit',evidence,hitColumns);var unique=new Set(hits.map(h=>key(h.bin_id,h.protein)));E('pathway-hit-status').textContent=number(hits.length,0)+' hits · '+number(unique.size,0)+' proteins'+(missing.length?' · '+number(missing.length,0)+' missing components':'');
  var qc=E('pathway-qc-body');qc.textContent='';bins.filter(b=>bin.value==='all'||b.bin_id===bin.value).forEach(function(b){var tr=text('tr');[b.bin_id,badge(b.status),badge(b.source_state),number(b.bases,0),number(b.orfs,0),b.error].forEach(v=>cell(tr,v));qc.appendChild(tr);});
  var annotationKey=D.rank+'/'+bin.value;E('pathway-bakta').disabled=bin.value==='all'||!Array.from(document.querySelectorAll('[data-annotation-bin]')).some(b=>b.getAttribute('data-annotation-bin')===annotationKey);
 }
 function tsvDownload(rows,columns,filename){function escape(v){v=String(v??'');return /["\t\r\n]/.test(v)?'"'+v.replace(/"/g,'""')+'"':v;}var lines=[columns.map(c=>escape(c[0])).join('\t')].concat(rows.map(r=>columns.map(c=>escape(c[1](r))).join('\t')));download(lines.join('\n')+'\n',filename,'text/tab-separated-values;charset=utf-8');}
 function download(content,name,mime){var url=URL.createObjectURL(new Blob([content],{type:mime||'text/plain;charset=utf-8'})),a=text('a');a.href=url;a.download=name;document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);}
 function fasta(b,p,field){var g=D.proteins[key(b,p)];return g&&g[field]?'>'+b+'|'+p+' '+(g.original_contig||'')+'\n'+g[field].match(/.{1,70}/g).join('\n')+'\n':'';}
 ['summary','hit'].forEach(function(name){E('pathway-'+name+'-prev').addEventListener('click',()=>{pages[name]--;refresh();});E('pathway-'+name+'-next').addEventListener('click',()=>{pages[name]++;refresh();});});
 [bin,path,hide,E('pathway-hit-filter')].forEach(el=>el.addEventListener('change',()=>{resetPages();refresh();}));search.addEventListener('input',()=>{resetPages();refresh();});
 E('pathway-reset').addEventListener('click',()=>{bin.value='all';path.value='all';search.value='';hide.checked=false;E('pathway-hit-filter').value='overview';resetPages();sorts={};refresh();});
 E('pathway-bakta').addEventListener('click',()=>document.dispatchEvent(new CustomEvent('metahopper:annotation-bin',{detail:{key:D.rank+'/'+bin.value}})));
 E('pathway-export-summary').addEventListener('click',()=>tsvDownload(current.summary,[['bin',r=>r.bin_id],['pathway',r=>r.pathway],['coverage_percent',r=>r.panel_percent],['detected',r=>r.detected_scored_components],['measured',r=>r.total_scored_components],['status',r=>label(r.panel_status)]],'pathways.tsv'));
 E('pathway-export-hits').addEventListener('click',()=>tsvDownload(current.hit,[['bin',r=>r.bin_id],['pathway',r=>hitPaths(r).map(p=>p.id).join('; ')],['gene',r=>r.target_gene],['product',r=>r.source_product],['protein',r=>r.protein],['hmm',r=>r.model_accession],['match',r=>label(matchStatus(r))],['evalue',r=>r.independent_evalue],['hmm_coverage',r=>r.hmm_coverage],['contig',r=>r.contig],['start',r=>r.start],['end',r=>r.end]],'hmm_hits.tsv'));
 E('pathway-export-panels').addEventListener('click',()=>tsvDownload(panelRows().filter(r=>path.value==='all'||r.pathway===path.value),[['pathway',r=>r.pathway],['gene',r=>r.gene],['component',r=>r.component],['hmm_models',r=>r.models],['counted_in_coverage',r=>r.counted_in_coverage],['measured',r=>r.measured],['notes',r=>r.notes]],'pathway_panels.tsv'));
 E('pathway-panel-close').addEventListener('click',()=>E('pathway-panel-dialog').close());
 E('pathway-protein-close').addEventListener('click',()=>E('pathway-protein-dialog').close());
 [['aa','faa'],['nt','fna']].forEach(pair=>E('pathway-protein-'+pair[0]+'-download').addEventListener('click',()=>{if(selectedProtein)download(fasta(selectedProtein.bin,selectedProtein.protein,pair[0]),selectedProtein.bin.replace(/[^A-Za-z0-9_.-]/g,'_')+'.'+selectedProtein.protein+'.'+pair[1]);}));
 E('pathway-locate-contig').addEventListener('click',()=>{if(!selectedProtein)return;var g=D.proteins[key(selectedProtein.bin,selectedProtein.protein)]||{},h=D.hits.find(h=>h.bin_id===selectedProtein.bin&&h.protein===selectedProtein.protein),contig=g.original_contig||(h&&h.contig);var q=document.getElementById('q'),tab=document.querySelector('.tab[data-id="binarena"]');if(q&&tab&&contig){E('pathway-protein-dialog').close();tab.click();q.value=contig;q.dispatchEvent(new Event('input',{bubbles:true}));q.scrollIntoView({block:'center'});}});
 E('pathway-run-details').textContent='PathwayGenie '+D.version+' · rank '+D.rank+' · finished '+D.finished_utc+' · '+D.directory+'.';
 Object.keys(D.parameters||{}).forEach(k=>{E('pathway-parameters').appendChild(text('dt',k.replace(/_/g,' ')));E('pathway-parameters').appendChild(text('dd',D.parameters[k]===null?'Per-template default':String(D.parameters[k])));});
 bin.value='all';path.value='all';hide.checked=false;E('pathway-hit-filter').value='overview';refresh();
})();
