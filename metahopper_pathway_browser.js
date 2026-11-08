/* Inlined report explorer; no remote scripts or external UI dependencies. */
(function(){
 'use strict';
 var source=document.getElementById('metahopper-pathways-data');if(!source)return;
 var D=JSON.parse(source.textContent),E=function(id){return document.getElementById(id);},bins=D.bins||[],paths=D.pathways||[],pages={summary:0,component:0,hit:0},sorts={},pageSize=100,current={},selectedProtein=null;
 var bin=E('pathway-bin'),path=E('pathway-pathway'),search=E('pathway-search'),hide=E('pathway-hide-zero');
 var labels={complete_target_panel:'Scored panel complete',partial_target_panel:'Partial panel',not_detected_in_panel:'No accepted panel evidence',supported_selected_route:'Selected route supported',unresolved_unmeasured_steps:'Unmeasured route steps',partial_evidence:'Partial route evidence',bin_failed:'Bin failed',stale_input:'Stale input',missing_bin:'Current bin missing',unverified:'Input not verified',current:'Input verified',preview_unverified:'Preview: input not checked',detected:'Detected',candidate_only:'Candidate only',not_detected:'Not detected',not_assayed:'Not assayed',fusion_only_not_detected:'Fusion-only: no accepted evidence',accepted:'Accepted',context_candidate:'Context candidate',fragment_candidate:'Fragment candidate',shared_family:'Shared family',domain_only:'Domain only',specific:'Specific target',ambiguous:'Ambiguous',competing_hit:'Competing hit',targeted:'Targeted',fusion_only:'Selected fusions only',ok:'Completed',no_predicted_genes:'No predicted genes',failed:'Failed'};
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
 function badge(value){var good=['accepted','detected','complete_target_panel','supported_selected_route','current','ok'].includes(value),bad=['bin_failed','stale_input','missing_bin','failed'].includes(value),warn=['partial_target_panel','unresolved_unmeasured_steps','partial_evidence','candidate_only','context_candidate','fragment_candidate','shared_family','domain_only','ambiguous','competing_hit','unverified','preview_unverified'].includes(value);return text('span',label(value),'pathway-tag '+(good?'ok':bad?'bad':warn?'partial':''));}
 function cell(tr,value,cls){var td=text('td',undefined,cls||'l');if(value&&typeof value==='object'&&value.nodeType)td.appendChild(value);else td.textContent=value===null||value===undefined||value===''?'—':String(value);tr.appendChild(td);return td;}
 function fileLink(asset,title){var a=text('a',title||asset.label,'asset-link pathway-file-link');a.title=asset.label;a.href=asset.path.split('/').map(encodeURIComponent).join('/');a.setAttribute('data-sample-asset',asset.path);a.download='';return a;}
 function resetPages(){Object.keys(pages).forEach(function(k){pages[k]=0;});}
 function pick(b,p){bin.value=b;path.value=p||'all';search.value='';resetPages();refresh();E('pathway-components-card').scrollIntoView({block:'start'});}
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
 function proteinButtons(b,value){var span=text('span');split(value).forEach(function(p,i){if(i)span.appendChild(document.createTextNode('; '));span.appendChild(button(p,function(){showProtein(b,p);}));});if(!span.childNodes.length)span.textContent='—';return span;}
 function paginated(name,rows,columns){
  var body=E('pathway-'+name+'-body'),table=body.closest('table');
  var sort=sorts[name];if(sort)rows=rows.slice().sort(function(a,b){var x=columns[sort.column].value(a),y=columns[sort.column].value(b);if(x===null||x===undefined)return y===null||y===undefined?0:1;if(y===null||y===undefined)return -1;return (typeof x==='number'&&typeof y==='number'?x-y:String(x).localeCompare(String(y),undefined,{numeric:true}))*sort.direction;});
  current[name]=rows;var totalPages=Math.max(1,Math.ceil(rows.length/pageSize));pages[name]=Math.max(0,Math.min(pages[name],totalPages-1));body.textContent='';
  rows.slice(pages[name]*pageSize,(pages[name]+1)*pageSize).forEach(function(row){var tr=text('tr');columns.forEach(function(c){cell(tr,c.render?c.render(row):c.value(row),c.cls||'l');});body.appendChild(tr);});
  E('pathway-'+name+'-page').textContent=number(rows.length,0)+' rows · page '+(pages[name]+1)+' of '+totalPages;E('pathway-'+name+'-prev').disabled=pages[name]===0;E('pathway-'+name+'-next').disabled=pages[name]+1>=totalPages;
  Array.from(table.querySelectorAll('thead th')).forEach(function(th,i){if(th.dataset.pathwaySort)return;th.dataset.pathwaySort='1';th.tabIndex=0;th.setAttribute('role','button');th.title='Sort by '+th.textContent;function sortRows(){sorts[name]={column:i,direction:sorts[name]&&sorts[name].column===i?-sorts[name].direction:1};pages[name]=0;refresh();}th.addEventListener('click',sortRows);th.addEventListener('keydown',function(e){if(e.key==='Enter'||e.key===' '){e.preventDefault();sortRows();}});});
 }
 var componentByKey=Object.create(null),definitionsByPath=Object.create(null);
 D.components.forEach(function(c){componentByKey[key(c.bin_id,c.component)]=c;});D.definitions.forEach(function(d){(definitionsByPath[d.pathway]||(definitionsByPath[d.pathway]=[])).push(d);});
 function hitInPath(hit,pid){var p=paths.find(p=>p.id===pid);return (definitionsByPath[pid]||[]).some(d=>split(hit.components).includes(d.component))||!!(p&&(p.target_genes||[]).includes(hit.target_gene));}
 option(bin,'all','All analyzed bins');bins.forEach(function(b){option(bin,b.bin_id,b.bin_id.replace(/_/g,' '));});
 option(path,'all','All pathways');paths.forEach(function(p){option(path,p.id,(p.label||p.id).replace(/_/g,' '));});
 (D.warnings||[]).forEach(function(w){E('pathway-warnings').appendChild(text('p',w,'pathway-warning'));});
 var uniqueProteins=new Set(D.hits.filter(function(h){return h.status==='accepted';}).map(function(h){return key(h.bin_id,h.protein);}));
 [['Analyzed bins',bins.length],['Pathway panels',paths.length],['Proteins with accepted HMM hits',uniqueProteins.size],['Bins with failed runs',bins.filter(function(b){return b.status==='failed';}).length]].forEach(function(item){var card=text('div',undefined,'stat');card.appendChild(text('div',item[0],'k'));card.appendChild(text('div',number(item[1],0),'v'));E('pathway-stats').appendChild(card);});
 var summaryFile=D.files.find(function(f){return f.label==='pathway_summary.long.csv';});if(summaryFile&&!D.preview){var a=E('pathway-summary-download');a.href=summaryFile.path.split('/').map(encodeURIComponent).join('/');a.setAttribute('data-sample-asset',summaryFile.path);}else E('pathway-summary-download').hidden=true;
 var summaryColumns=[
  {value:r=>r.bin_id,render:r=>button(r.bin_id,function(){pick(r.bin_id,path.value);} )},
  {value:r=>r.pathway,render:r=>button(r.pathway.replace(/_/g,' '),function(){pick(r.bin_id,r.pathway);})},
  {value:r=>n(r.panel_percent),render:r=>n(r.panel_percent)===null?'NA':number(r.panel_percent)+'%'},
  {value:r=>n(r.detected_scored_components),render:r=>n(r.detected_scored_components)===null?'NA':r.detected_scored_components+' / '+r.total_scored_components},
  {value:r=>r.panel_status,render:r=>badge(r.panel_status)},
  {value:r=>r.route_status,render:r=>badge(r.route_status)},
  {value:r=>n(r.operon_support),render:r=>number(r.operon_support)},
  {value:r=>r.source_state,render:r=>badge(r.source_state)}];
 var componentColumns=[{value:r=>r.bin_id},{value:r=>r.pathway+' / '+(r.label||r.component)},{value:r=>r.assay_scope,render:r=>badge(r.assay_scope)},
  {value:r=>r.in_heatmap_denominator,render:r=>r.in_heatmap_denominator==='1'?'Yes':'No'},
  {value:r=>r.in_route_definition,render:r=>r.in_route_definition==='1'?'Yes':'No'},
  {value:r=>r.state,render:r=>badge(r.state)},{value:r=>n(r.copies),render:r=>number(r.copies,0)},
  {value:r=>r.proteins,render:r=>proteinButtons(r.bin_id,r.proteins)},{value:r=>r.candidate_proteins,render:r=>proteinButtons(r.bin_id,r.candidate_proteins)},
  {value:r=>r.notes,cls:'l wrap-text'}];
 var hitColumns=[{value:r=>r.bin_id},{value:r=>r.protein,render:r=>button(r.protein,function(){showProtein(r.bin_id,r.protein);})},{value:r=>r.target_gene},
  {value:r=>r.source_product,render:r=>{var span=text('span',r.source_product+'\n'+r.model_accession+' ('+r.source_database+')');span.style.whiteSpace='pre-line';return span;},cls:'l wrap-text'},
  {value:r=>r.target_call,render:r=>badge(r.target_call)},
  {value:r=>r.status,render:r=>{var span=text('span');span.appendChild(badge(r.status));if(r.reason)span.appendChild(text('div',r.reason));return span;},cls:'l wrap-text'},
  {value:r=>n(r.hmm_coverage),render:r=>n(r.hmm_coverage)===null?'NA':number(n(r.hmm_coverage)*100,2)+'%'},
  {value:r=>n(r.sequence_bit_score),render:r=>number(r.sequence_bit_score)},
  {value:r=>n(r.domain_bit_score),render:r=>number(r.domain_bit_score)},
  {value:r=>n(r.independent_evalue),render:r=>r.independent_evalue},{value:r=>r.contig},
  {value:r=>n(r.start),render:r=>number(r.start,0)},{value:r=>n(r.end),render:r=>number(r.end,0)},{value:r=>r.strand,render:r=>r.strand==='-1'?'-':'+'}];
 function heatmap(rows){var host=E('pathway-heatmap');host.textContent='';if(!rows.length){host.appendChild(text('p','No matching pathway records.','note'));return;}
  var shownBins=Array.from(new Set(rows.map(r=>r.bin_id))),shownPaths=paths.filter(p=>rows.some(r=>r.pathway===p.id)),lookup=Object.create(null);rows.forEach(r=>{lookup[key(r.bin_id,r.pathway)]=r;});
  var table=text('table'),head=text('thead'),tr=text('tr');tr.appendChild(text('th','Bin'));shownPaths.forEach(p=>tr.appendChild(text('th',(p.label||p.id).replace(/_/g,' '))));head.appendChild(tr);table.appendChild(head);var body=text('tbody');
  shownBins.forEach(function(b){var tr=text('tr');tr.appendChild(text('th',b));shownPaths.forEach(function(p){var row=lookup[key(b,p.id)],value=row?n(row.panel_percent):null,td=text('td'),control=text('button',value===null?'NA':number(value)+'%');control.type='button';
   var f=(value||0)/100;control.style.background=value===null?'#ddd7e2':value===0?'#eff3f4':'rgb('+Math.round(228-205*f)+','+Math.round(241-132*f)+','+Math.round(241-120*f)+')';control.style.color=value!==null&&value>=65?'#fff':'#163b43';
   control.title=b+' / '+p.id+' — '+(value===null?'Unavailable':number(value)+'% scored components')+'; '+(row?label(row.route_status):'No record');control.setAttribute('aria-label',control.title);control.addEventListener('click',function(){pick(b,p.id);});td.appendChild(control);tr.appendChild(td);});body.appendChild(tr);});table.appendChild(body);host.appendChild(table);
 }
 function neighborhoodFigure(row){var host=E('pathway-neighborhood-figure');host.textContent='';var NS='http://www.w3.org/2000/svg';function svgEl(tag,attrs,content){var e=document.createElementNS(NS,tag);Object.keys(attrs).forEach(k=>e.setAttribute(k,String(attrs[k])));if(content!==undefined)e.textContent=content;return e;}
  var genes=(row.orfs||[]).map(p=>D.proteins[key(row.bin_id,p)]).filter(Boolean).sort((a,b)=>Number(a.start)-Number(b.start));
  if(!genes.length){E('pathway-neighborhood-details').textContent='No per-gene coordinates are available for this neighborhood.';return;}
  var left=Math.min(Number(row.start),...genes.map(g=>Number(g.start))),right=Math.max(Number(row.end),...genes.map(g=>Number(g.end))),width=1000,padding=35,scale=(width-2*padding)/Math.max(1,right-left);
  var svg=svgEl('svg',{viewBox:'0 0 1000 140',role:'img','aria-label':row.bin_id+' / '+row.neighborhood_id+' gene neighborhood'});svg.appendChild(svgEl('line',{x1:padding,x2:width-padding,y1:66,y2:66,stroke:'#aebdc4'}));
  genes.forEach(function(g,i){var a=padding+(Number(g.start)-left)*scale,b=padding+(Number(g.end)-left)*scale,tip=Math.min(9,Math.max(2,(b-a)/3)),minus=String(g.strand)==='-1',supported=D.hits.some(h=>h.bin_id===row.bin_id&&h.protein===g.protein&&h.status==='accepted');
   var points=minus?[[a,66],[a+tip,49],[b,49],[b,83],[a+tip,83]]:[[a,49],[b-tip,49],[b,66],[b-tip,83],[a,83]];
   var polygon=svgEl('polygon',{points:points.map(p=>p.join(',')).join(' '),fill:supported?'#207b86':'#cad4d9',stroke:'#fff','stroke-width':1,tabindex:0,role:'button'});
   polygon.appendChild(svgEl('title',{},g.protein+'; '+g.start+'–'+g.end+'; strand '+(minus?'-':'+')));polygon.addEventListener('click',()=>showProtein(row.bin_id,g.protein));polygon.addEventListener('keydown',e=>{if(e.key==='Enter'){showProtein(row.bin_id,g.protein);}});svg.appendChild(polygon);
   svg.appendChild(svgEl('text',{x:(a+b)/2,y:i%2?115:36,'text-anchor':'middle','font-size':10,fill:'#51616c'},g.protein));
  });svg.appendChild(svgEl('text',{x:padding,y:136,'font-size':10,fill:'#63757f'},number(left,0)+' bp'));svg.appendChild(svgEl('text',{x:width-padding,y:136,'font-size':10,'text-anchor':'end',fill:'#63757f'},number(right,0)+' bp'));host.appendChild(svg);
  var ops=current.operons.filter(o=>o.bin_id===row.bin_id&&o.pathway===row.pathway&&o.neighborhood_id===row.neighborhood_id);E('pathway-neighborhood-details').textContent=row.bin_id+' / '+row.pathway+' / '+row.contig+' · '+row.observed_order+' · teal: protein with accepted model evidence; gray: other predicted ORF. '+ops.map(o=>o.template_id+': score '+number(o.score)+'; order '+number(o.order_score)+'; strand '+number(o.strand_score)+'; '+o.reference_taxon+'; '+o.source_url).join(' | ');
 }
 function refresh(){
  var q=search.value.trim().toLowerCase(),chosenPaths=new Set(paths.filter(p=>path.value==='all'||p.id===path.value).map(p=>p.id));
  var selectedSummaries=D.summary.filter(r=>(bin.value==='all'||r.bin_id===bin.value)&&chosenPaths.has(r.pathway));
  var visibleBins=new Set(selectedSummaries.map(r=>r.bin_id)),hiddenZeros=0;
  if(hide.checked&&bin.value==='all')Array.from(visibleBins).forEach(function(b){var rows=selectedSummaries.filter(r=>r.bin_id===b);if(rows.length&&rows.every(r=>n(r.panel_percent)===0)){visibleBins.delete(b);hiddenZeros++;}});
  var queryHits=D.hits.filter(h=>includes(h,q)),matchingPairs=new Set();if(q)queryHits.forEach(function(h){chosenPaths.forEach(function(pid){if(hitInPath(h,pid))matchingPairs.add(key(h.bin_id,pid));});});
  var summaryRows=selectedSummaries.filter(r=>visibleBins.has(r.bin_id)&&(!q||includes(r,q)||matchingPairs.has(key(r.bin_id,r.pathway))));
  var activeBins=new Set(summaryRows.map(r=>r.bin_id));
  E('pathway-selection-status').textContent=number(activeBins.size,0)+' bins shown · '+number(summaryRows.length,0)+' bin/pathway records · '+hiddenZeros+' all-zero bins hidden. Coverage is scored component detection, not metabolic flux.';
  paginated('summary',summaryRows,summaryColumns);heatmap(summaryRows);
  var componentRows=[];summaryRows.forEach(function(row){(definitionsByPath[row.pathway]||[]).forEach(function(def){var c=componentByKey[key(row.bin_id,def.component)];var valid=row.source_state==='current'||row.source_state==='preview_unverified';componentRows.push(Object.assign({},def,c||{}, {bin_id:row.bin_id,pathway:row.pathway,component:def.component,label:def.label||def.component,state:!valid?row.source_state:c?c.state:row.panel_status==='bin_failed'?'bin_failed':'unavailable',copies:c?c.copies:null}));});});
  if(q)componentRows=componentRows.filter(c=>includes(c,q)||includes({bin_id:c.bin_id,pathway:c.pathway},q)||queryHits.some(h=>h.bin_id===c.bin_id&&split(h.components).includes(c.component)));
  paginated('component',componentRows,componentColumns);
  var chosen=paths.find(p=>p.id===path.value);E('pathway-component-status').textContent=chosen?(chosen.notes||'')+' '+(chosen.route_scope||''):'Select a heatmap cell or choose a pathway/bin to focus its scored, unmeasured and supporting components.';
  var filter=E('pathway-hit-filter').value;
  var hits=D.hits.filter(h=>activeBins.has(h.bin_id)&&(path.value==='all'||hitInPath(h,path.value))&&includes(h,q)&&(filter!=='accepted'||h.status==='accepted')&&(filter!=='review'||h.status!=='accepted'||h.target_call!=='specific'));
  paginated('hit',hits,hitColumns);var unique=new Set(hits.map(h=>key(h.bin_id,h.protein)));E('pathway-hit-status').textContent=number(hits.length,0)+' model hits on '+number(unique.size,0)+' predicted proteins. Multiple models on one protein are not extra gene copies.';
  E('pathway-export-proteins').disabled=!hits.some(h=>(D.proteins[key(h.bin_id,h.protein)]||{}).aa);
  current.neighborhoods=D.neighborhoods.filter(r=>activeBins.has(r.bin_id)&&chosenPaths.has(r.pathway)&&(!q||includes(r,q)||summaryRows.some(s=>s.bin_id===r.bin_id&&s.pathway===r.pathway)));
  current.operons=D.operons.filter(r=>current.neighborhoods.some(n=>n.bin_id===r.bin_id&&n.pathway===r.pathway&&n.neighborhood_id===r.neighborhood_id));
  var body=E('pathway-neighborhood-body');body.textContent='';current.neighborhoods.forEach(function(row){var tr=text('tr'),ops=current.operons.filter(o=>o.bin_id===row.bin_id&&o.neighborhood_id===row.neighborhood_id&&o.pathway===row.pathway);cell(tr,row.bin_id);cell(tr,row.pathway);cell(tr,button(row.neighborhood_id,function(){neighborhoodFigure(row);}));cell(tr,row.contig);cell(tr,number(row.start,0));cell(tr,number(row.end,0));cell(tr,row.protein_count);cell(tr,row.observed_order,'l wrap-text');cell(tr,ops.map(o=>o.template_id).join('; '));cell(tr,ops.length?ops.map(o=>number(o.score)).join('; '):'NA');body.appendChild(tr);});
  if(current.neighborhoods.length)neighborhoodFigure(current.neighborhoods[0]);else{E('pathway-neighborhood-figure').textContent='';E('pathway-neighborhood-details').textContent='No matching neighborhoods in these outputs.';}
  var qc=E('pathway-qc-body');qc.textContent='';bins.filter(b=>bin.value==='all'||b.bin_id===bin.value).forEach(function(b){var tr=text('tr');[b.bin_id,badge(b.status),badge(b.source_state),number(b.contigs,0),number(b.bases,0),number(b.orfs,0),b.accepted_models,b.accepted_proteins,b.review_hits,b.error].forEach(v=>cell(tr,v));qc.appendChild(tr);});
  var annotationKey=D.rank+'/'+bin.value;E('pathway-bakta').disabled=bin.value==='all'||!Array.from(document.querySelectorAll('[data-annotation-bin]')).some(b=>b.getAttribute('data-annotation-bin')===annotationKey);
 }
 function csvDownload(rows,filename){var keys=Array.from(new Set(rows.flatMap(r=>Object.keys(r).filter(k=>k!=='orfs'))));function escape(v){v=String(v??'');return /[",\r\n]/.test(v)?'"'+v.replace(/"/g,'""')+'"':v;}var lines=[keys.map(escape).join(',')].concat(rows.map(r=>keys.map(k=>escape(r[k])).join(',')));download(lines.join('\n')+'\n',filename,'text/csv;charset=utf-8');}
 function download(content,name,mime){var url=URL.createObjectURL(new Blob([content],{type:mime||'text/plain;charset=utf-8'})),a=text('a');a.href=url;a.download=name;document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);}
 function fasta(b,p,field){var g=D.proteins[key(b,p)];return g&&g[field]?'>'+b+'|'+p+' '+(g.original_contig||'')+'\n'+g[field].match(/.{1,70}/g).join('\n')+'\n':'';}
 function files(){var host=E('pathway-file-links');host.textContent='';if(D.preview){E('pathway-file-status').textContent='This standalone preview embeds the evidence and selected sequences. Filtered CSV/FASTA downloads above work here. Original files are in the supplied PathwayGenie output folder and become sample downloads after report generation/publication.';return;}var q=E('pathway-file-search').value.trim().toLowerCase(),all=D.files.filter(f=>f.label.toLowerCase().includes(q)&&(bin.value==='all'||!f.label.startsWith('bins/')||f.label.includes('/'+bin.value.replace(/[^A-Za-z0-9_.-]+/g,'_').slice(0,70)+'_')));all.slice(0,60).forEach(f=>host.appendChild(fileLink(f)));E('pathway-file-status').textContent=all.length+' matching files'+(all.length>60?'; first 60 shown. Refine the file filter to find another file.':'.');}
 ['summary','component','hit'].forEach(function(name){E('pathway-'+name+'-prev').addEventListener('click',()=>{pages[name]--;refresh();});E('pathway-'+name+'-next').addEventListener('click',()=>{pages[name]++;refresh();});});
 [bin,path,hide,E('pathway-hit-filter')].forEach(el=>el.addEventListener('change',()=>{resetPages();refresh();files();}));search.addEventListener('input',()=>{resetPages();refresh();});E('pathway-file-search').addEventListener('input',files);
 E('pathway-reset').addEventListener('click',()=>{bin.value='all';path.value='all';search.value='';hide.checked=true;E('pathway-hit-filter').value='all';resetPages();sorts={};refresh();files();});
 E('pathway-bakta').addEventListener('click',()=>document.dispatchEvent(new CustomEvent('metahopper:annotation-bin',{detail:{key:D.rank+'/'+bin.value}})));
 [['summary','pathways'],['component','components'],['hit','hits']].forEach(function(pair){var id=pair[0]==='hit'?'hits':pair[0]==='component'?'components':'summary';E('pathway-export-'+id).addEventListener('click',()=>csvDownload(current[pair[0]],pair[1]+'.filtered.csv'));});
 E('pathway-export-neighborhoods').addEventListener('click',()=>csvDownload(current.neighborhoods,'neighborhoods.filtered.csv'));E('pathway-export-operons').addEventListener('click',()=>csvDownload(current.operons,'operons.filtered.csv'));
 E('pathway-export-proteins').addEventListener('click',()=>{var seen=new Set(),content='';current.hit.forEach(h=>{var k=key(h.bin_id,h.protein);if(!seen.has(k)){content+=fasta(h.bin_id,h.protein,'aa');seen.add(k);}});download(content,'pathwaygenie.filtered_proteins.faa');});
 E('pathway-protein-close').addEventListener('click',()=>E('pathway-protein-dialog').close());
 [['aa','faa'],['nt','fna']].forEach(pair=>E('pathway-protein-'+pair[0]+'-download').addEventListener('click',()=>{if(selectedProtein)download(fasta(selectedProtein.bin,selectedProtein.protein,pair[0]),selectedProtein.bin.replace(/[^A-Za-z0-9_.-]/g,'_')+'.'+selectedProtein.protein+'.'+pair[1]);}));
 E('pathway-locate-contig').addEventListener('click',()=>{if(!selectedProtein)return;var g=D.proteins[key(selectedProtein.bin,selectedProtein.protein)]||{},h=D.hits.find(h=>h.bin_id===selectedProtein.bin&&h.protein===selectedProtein.protein),contig=g.original_contig||(h&&h.contig);var q=document.getElementById('q'),tab=document.querySelector('.tab[data-id="binarena"]');if(q&&tab&&contig){E('pathway-protein-dialog').close();tab.click();q.value=contig;q.dispatchEvent(new Event('input',{bubbles:true}));q.scrollIntoView({block:'center'});}});
 E('pathway-run-details').textContent='PathwayGenie '+D.version+' · rank '+D.rank+' · finished '+D.finished_utc+' · '+D.directory+'.';
 Object.keys(D.parameters||{}).forEach(k=>{E('pathway-parameters').appendChild(text('dt',k.replace(/_/g,' ')));E('pathway-parameters').appendChild(text('dd',D.parameters[k]===null?'Per-template default':String(D.parameters[k])));});
 bin.value='all';path.value='all';refresh();files();
})();
