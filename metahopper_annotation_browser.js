/* Embedded by metahopper_annotation_view.py; no network dependencies. */
(function(){
 'use strict';
 var dataNode=document.getElementById('annotation-feature-data');if(!dataNode)return;
 var data=JSON.parse(dataNode.textContent),bins=data.bins||[],fields=data.fields||[],index=Object.create(null),all=[],matches=[],pool=[],page=0,pageSize=100,sortField='',sortDirection=1,selectedRecord=null;
 fields.forEach(function(f,i){index[f]=i;});
 var el=function(id){return document.getElementById(id);},select=el('annotation-bin'),search=el('annotation-search'),type=el('annotation-type'),category=el('annotation-category');
 var typeNames={cds:'CDS',sorf:'sORF',trna:'tRNA',tmrna:'tmRNA',rrna:'rRNA',ncrna:'ncRNA','ncrna-region':'ncRNA region',crispr:'CRISPR',gap:'Gap',oric:'oriC',oriv:'oriV',orit:'oriT'};
 function protein(r){var kind=String(r.row[index.type]||'').toLowerCase();return kind==='cds'||kind==='sorf';}
 function typeName(kind){return typeNames[kind]||kind||'Unknown';}
 function hyp(r){if(r.detail.hypothetical!==null&&r.detail.hypothetical!==undefined)return r.detail.hypothetical===true;var p=String(r.row[index.product]||'');return protein(r)&&p?p.toLowerCase().indexOf('hypothetical protein')===0:null;}
 function pseudo(r){return r.detail.pseudogene===true?true:r.detail.pseudogene===false?false:null;}
 function categoryName(r){return !protein(r)?'Not a protein feature':hyp(r)===true?'Hypothetical':hyp(r)===false?'Non-hypothetical':'Unassigned product';}
 function pseudoName(r){return pseudo(r)===true?'Flagged':pseudo(r)===false?'Not flagged':protein(r)?'Unavailable':'—';}
 function rebuild(){all=[];bins.forEach(function(bin){(bin.features||[]).forEach(function(row,i){all.push({row:row,detail:(bin.details||[])[i]||{},bin:bin,uid:all.length});});});}
 function option(node,value,text){var n=document.createElement('option');n.value=value;n.textContent=text;node.appendChild(n);}
 option(select,'all','All annotated bins');bins.forEach(function(bin){option(select,bin.key,bin.rank+' / '+bin.bin.replace(/_/g,' '));});
 function updateTypes(){var old=type.value;type.textContent='';option(type,'all','All feature types');option(type,'proteins','Proteins (CDS + sORF)');Array.from(new Set(pool.map(function(r){return r.row[index.type].toLowerCase();}))).sort().forEach(function(kind){option(type,kind,typeName(kind));});type.value=Array.from(type.options).some(function(o){return o.value===old;})?old:'all';}
 function matchesQuery(r,query,featureType,featureCategory){
  var kind=String(r.row[index.type]||'').toLowerCase();
  if(featureType==='proteins'&&!protein(r))return false;
  if(featureType!=='all'&&featureType!=='proteins'&&kind!==featureType)return false;
  if(featureCategory==='hypothetical'&&(!protein(r)||hyp(r)!==true))return false;
  if(featureCategory==='named'&&(!protein(r)||hyp(r)!==false))return false;
  if(featureCategory==='unassigned'&&(!protein(r)||hyp(r)!==null))return false;
  if(featureCategory==='pseudogene'&&pseudo(r)!==true)return false;
  if(featureCategory==='pseudo-unknown'&&(!protein(r)||pseudo(r)!==null))return false;
  return !query||r.row.concat([r.bin.bin,r.bin.rank,r.detail.pseudogene_note||'']).some(function(v){return String(v).toLocaleLowerCase().indexOf(query)>=0;});
 }
 function link(id,path){var n=el(id);n.hidden=!path;if(path){n.setAttribute('data-sample-asset',path);n.setAttribute('href',path.split('/').map(encodeURIComponent).join('/'));}}
 function chooseBin(){pool=all.filter(function(r){return select.value==='all'||r.bin.key===select.value;});page=0;updateTypes();render(true);}
 function chip(label,value,filter){var button=document.createElement('button');button.type='button';button.className=category.value===filter?'active':'';button.textContent=label+' · '+value.toLocaleString();button.addEventListener('click',function(){category.value=filter;type.value='all';page=0;render(false);});return button;}
 function charts(){if(!window.MHCharts)return;var types={};pool.forEach(function(r){var kind=r.row[index.type].toLowerCase();types[kind]=(types[kind]||0)+1;});
  window.MHCharts.bar(el('annotation-types-figure'),Object.keys(types).sort(function(a,b){return types[b]-types[a];}).map(function(kind){return {label:typeName(kind),kind:kind,values:[types[kind]]};}),{title:'Feature types',unit:'Annotated features',empty:'Gene records are not embedded in this report yet.',onPick:function(row){type.value=row.kind;category.value='all';page=0;render(false);el('annotation-browser').scrollIntoView({block:'start'});}});
  var named=pool.filter(function(r){return protein(r)&&hyp(r)===false;}).length,hypothetical=pool.filter(function(r){return protein(r)&&hyp(r)===true;}).length,unknown=pool.filter(function(r){return protein(r)&&hyp(r)===null;}).length;
  var rows=[{label:'Non-hypothetical',kind:'named',values:[named],colors:['#207b86']},{label:'Hypothetical',kind:'hypothetical',values:[hypothetical],colors:['#dc9747']},{label:'Unassigned product',kind:'unassigned',values:[unknown],colors:['#9baab2']}].filter(function(r){return r.values[0]>0;});
  window.MHCharts.bar(el('annotation-products-figure'),rows,{title:'Protein product categories',unit:'Predicted protein features',empty:'No protein feature records are loaded.',onPick:function(row){category.value=row.kind;type.value='proteins';page=0;render(false);el('annotation-browser').scrollIntoView({block:'start'});}});
 }
 function sortValue(r,field){if(field==='bin')return r.bin.key;if(field==='pseudogene')return pseudoName(r);if(field==='aa_length')return r.detail.aa_length;if(field==='start'||field==='end')return Number(r.row[index[field]])||0;return r.row[index[field]]||'';}
 function textCell(tr,text,cls){var td=document.createElement('td');td.textContent=text===undefined||text===null||text===''?'—':String(text);if(cls)td.className=cls;tr.appendChild(td);return td;}
 function render(redrawCharts){
  var query=search.value.trim().toLocaleLowerCase();matches=pool.filter(function(r){return matchesQuery(r,query,type.value,category.value);});
  if(sortField)matches.sort(function(a,b){var x=sortValue(a,sortField),y=sortValue(b,sortField);return (typeof x==='number'&&typeof y==='number'?x-y:String(x??'').localeCompare(String(y??''),undefined,{numeric:true}))*sortDirection;});
  var body=el('annotation-feature-body');body.textContent='';var pages=Math.max(1,Math.ceil(matches.length/pageSize));page=Math.max(0,Math.min(page,pages-1));
  matches.slice(page*pageSize,(page+1)*pageSize).forEach(function(r){var row=r.row,tr=document.createElement('tr');
   textCell(tr,r.bin.rank+' / '+r.bin.bin,'l');var locus=textCell(tr,'','l'),button=document.createElement('button');button.type='button';button.className='annotation-feature-link';button.textContent=row[index.locus_tag]||row[index.type]+' '+row[index.start];button.addEventListener('click',function(){openFeature(r);});locus.textContent='';locus.appendChild(button);
   textCell(tr,row[index.gene],'l');textCell(tr,row[index.product],'l product-cell');textCell(tr,typeName(row[index.type]),'l');var cat=textCell(tr,'','l'),badge=document.createElement('span');badge.className='feature-tag '+(hyp(r)===true?'hypothetical':protein(r)?'named':'');badge.textContent=categoryName(r);cat.textContent='';cat.appendChild(badge);
   var ps=textCell(tr,pseudoName(r),'l');if(pseudo(r)===true)ps.className+=' feature-pseudogene';textCell(tr,row[index.contig],'l');textCell(tr,row[index.start]);textCell(tr,row[index.end]);textCell(tr,row[index.strand],'l');textCell(tr,r.detail.aa_length);textCell(tr,row[index.db_xrefs],'l reference-cell');body.appendChild(tr);
  });
  var selected=bins.find(function(b){return b.key===select.value;});
  el('annotation-browser-title').textContent='Gene and protein explorer'+(selected?' — '+selected.bin.replace(/_/g,' '):' — all bins');
  var missing=bins.filter(function(b){return (select.value==='all'||b===selected)&&b.error;});
  var message=matches.length.toLocaleString()+' matching features of '+pool.length.toLocaleString()+' loaded features.';
  if(!pool.length)message='No gene records are loaded for this selection. Regenerate this report from the existing Bakta outputs, or select a bin and open its Bakta JSON below.';
  if(selected&&selected.source)message+=' Source: '+selected.source+'.';
  if(missing.length)message+=' Gene records unavailable for '+missing.map(function(b){return b.bin+': '+b.error;}).join('; ')+'.';
  el('annotation-feature-status').textContent=message;el('annotation-page').textContent='Page '+(page+1)+' of '+pages+' · '+pageSize+' rows per page';el('annotation-prev').disabled=page===0;el('annotation-next').disabled=page+1>=pages;
  el('annotation-export').disabled=!matches.length;var aa=matches.filter(function(r){return protein(r)&&r.detail.aa;}).length,nt=matches.filter(function(r){return r.detail.nt;}).length;
  el('annotation-export-proteins').disabled=!aa;el('annotation-export-proteins').textContent='Filtered proteins FASTA ('+aa.toLocaleString()+')';el('annotation-export-genes').disabled=!nt;el('annotation-export-genes').textContent='Filtered genes FASTA ('+nt.toLocaleString()+')';
  link('annotation-feature-download',selected&&selected.features_path);link('annotation-protein-download',selected&&selected.protein_path);link('annotation-fasta-download',selected&&selected.fasta_path);
  var chips=el('annotation-category-chips');chips.textContent='';chips.appendChild(chip('All features',pool.length,'all'));chips.appendChild(chip('Hypothetical proteins',pool.filter(function(r){return protein(r)&&hyp(r)===true;}).length,'hypothetical'));chips.appendChild(chip('Pseudogenes',pool.filter(function(r){return pseudo(r)===true;}).length,'pseudogene'));
  var unknown=pool.filter(function(r){return protein(r)&&pseudo(r)===null;}).length;if(unknown)chips.appendChild(chip('Pseudogene status unavailable',unknown,'pseudo-unknown'));
  if(redrawCharts)charts();
 }
 function download(text,name,mime){var url=URL.createObjectURL(new Blob([text],{type:mime||'text/plain;charset=utf-8'})),a=document.createElement('a');a.href=url;a.download=name;document.body.appendChild(a);a.click();a.remove();setTimeout(function(){URL.revokeObjectURL(url);},1000);}
 function safeName(s){return String(s).replace(/[^A-Za-z0-9_.-]/g,'_');}
 function fasta(r,field){var seq=r.detail[field];if(!seq)return '';var row=r.row,head=row[index.locus_tag]||row[index.type]+'_'+row[index.start];return '>'+head+' '+r.bin.bin+' '+String(row[index.product]||'').replace(/[\r\n]/g,' ')+'\n'+seq.match(/.{1,70}/g).join('\n')+'\n';}
 function openFeature(r){selectedRecord=r;var row=r.row;el('annotation-detail-title').textContent=row[index.locus_tag]||'Annotated feature';el('annotation-detail-product').textContent=row[index.product]||'No product assigned';var dl=el('annotation-detail-fields');dl.textContent='';
  [['Bin',r.bin.key],['Gene',row[index.gene]],['Feature type',typeName(row[index.type])],['Product category',categoryName(r)],['Pseudogene',pseudoName(r)],['Pseudogene note',r.detail.pseudogene_note],['Contig',row[index.contig]],['Coordinates',row[index.start]+'–'+row[index.end]+' ('+row[index.strand]+')'],['Protein length',r.detail.aa_length===null||r.detail.aa_length===undefined?'Unavailable':r.detail.aa_length+' aa'],['Truncated end',r.detail.truncated],['Database references',row[index.db_xrefs]]].forEach(function(pair){var dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=pair[0];dd.textContent=pair[1]||'—';dl.appendChild(dt);dl.appendChild(dd);});
  el('annotation-detail-aa').textContent=r.detail.aa?r.detail.aa.match(/.{1,70}/g).join('\n'):'No protein sequence recorded for this feature.';el('annotation-detail-nt').textContent=r.detail.nt?r.detail.nt.match(/.{1,70}/g).join('\n'):'No nucleotide sequence is embedded for this feature.';
  el('annotation-detail-protein').disabled=!r.detail.aa;el('annotation-detail-gene').disabled=!r.detail.nt;el('annotation-feature-dialog').showModal();
 }
 function openBin(key,focus,cat){if(!bins.some(function(b){return b.key===key;}))return false;select.value=key;search.value='';category.value=cat||'all';type.value='all';chooseBin();if(focus){var tab=document.querySelector('.tab[data-id="annotation"]');if(tab)tab.click();el('annotation-browser').scrollIntoView({block:'start'});}return true;}
 select.addEventListener('change',chooseBin);search.addEventListener('input',function(){page=0;render(false);});type.addEventListener('change',function(){page=0;render(false);});category.addEventListener('change',function(){page=0;render(false);});
 el('annotation-clear').addEventListener('click',function(){search.value='';type.value='all';category.value='all';page=0;render(false);});
 el('annotation-prev').addEventListener('click',function(){page--;render(false);});el('annotation-next').addEventListener('click',function(){page++;render(false);});
 document.querySelectorAll('[data-feature-sort]').forEach(function(button){button.addEventListener('click',function(){var field=button.getAttribute('data-feature-sort');sortDirection=sortField===field?-sortDirection:1;sortField=field;page=0;render(false);});});
 el('annotation-export').addEventListener('click',function(){if(!matches.length)return;function cell(v){v=String(v??'');return /[\t\r\n"]/.test(v)?'"'+v.replace(/"/g,'""')+'"':v;}var headers=['rank','bin'].concat(fields,['hypothetical','pseudogene','aa_length','nt_length','pseudogene_note','truncated']);var lines=[headers.map(cell).join('\t')];matches.forEach(function(r){lines.push([r.bin.rank,r.bin.bin].concat(r.row,[hyp(r)===null?'unavailable':hyp(r),pseudo(r)===null?'unavailable':pseudo(r),r.detail.aa_length,r.detail.nt_length,r.detail.pseudogene_note,r.detail.truncated]).map(cell).join('\t'));});download(lines.join('\n')+'\n',safeName(select.value)+'.filtered_features.tsv','text/tab-separated-values;charset=utf-8');});
 el('annotation-export-proteins').addEventListener('click',function(){download(matches.filter(protein).map(function(r){return fasta(r,'aa');}).join(''),safeName(select.value)+'.filtered_proteins.faa');});
 el('annotation-export-genes').addEventListener('click',function(){download(matches.map(function(r){return fasta(r,'nt');}).join(''),safeName(select.value)+'.filtered_genes.ffn');});
 el('annotation-detail-close').addEventListener('click',function(){el('annotation-feature-dialog').close();});
 el('annotation-detail-protein').addEventListener('click',function(){if(selectedRecord)download(fasta(selectedRecord,'aa'),safeName(selectedRecord.row[index.locus_tag]||'feature')+'.faa');});
 el('annotation-detail-gene').addEventListener('click',function(){if(selectedRecord)download(fasta(selectedRecord,'nt'),safeName(selectedRecord.row[index.locus_tag]||'feature')+'.ffn');});
 document.addEventListener('click',function(e){var target=e.target.closest&&e.target.closest('[data-annotation-bin]');if(target){e.preventDefault();openBin(target.getAttribute('data-annotation-bin'),true);}});
 document.addEventListener('metahopper:annotation-bin',function(e){openBin(e.detail.key,true,e.detail.category);});
 el('annotation-load-json').addEventListener('change',async function(){var status=el('annotation-import-status'),file=this.files&&this.files[0],bin=bins.find(function(b){return b.key===select.value;});if(!file)return;if(!bin){status.textContent='Choose a specific bin before opening its Bakta JSON.';this.value='';return;}
  try{var imported=JSON.parse(await file.text());if(!Array.isArray(imported.features))throw Error('The file has no Bakta features array.');var rows=[],details=[],allowed=bin.expected_contigs&&bin.expected_contigs.length?new Set(bin.expected_contigs):null;
   function flag(v){if(typeof v==='string')return !['','false','no','0','none','null'].includes(v.trim().toLowerCase());if(v&&typeof v==='object'&&!Array.isArray(v))return true;return !!v;}
   imported.features.forEach(function(f){if(!f||typeof f!=='object'||Array.isArray(f))throw Error('A Bakta feature is not an object.');if((f.aa!=null&&typeof f.aa!=='string')||(f.nt!=null&&typeof f.nt!=='string'))throw Error('Bakta feature sequences must be strings.');var kind=String(f.type||'').toLowerCase(),contig=f.sequence||f.contig||f.sequence_id||'';if(contig&&typeof contig==='object')contig=contig.id||'';contig=String(contig);if(allowed&&!allowed.has(contig))throw Error('A feature belongs to a contig outside the selected bin: '+contig);var product=String(f.product||''),refs=f.db_xrefs||f.dbxref||[],aa=typeof f.aa==='string'?f.aa:'',nt=typeof f.nt==='string'?f.nt:'',pk='pseudogene' in f?'pseudogene':'pseudo' in f?'pseudo':null,pv=pk?f[pk]:null;
    rows.push([String(f.locus||f.locus_tag||f.id||''),String(f.gene||''),product,kind,contig,String(f.start||''),String(f.stop??f.end??''),String(f.strand||''),Array.isArray(refs)?refs.join(';'):String(refs)]);
    details.push({hypothetical:['cds','sorf'].includes(kind)?('hypothetical' in f&&f.hypothetical!=null?flag(f.hypothetical):product?product.toLowerCase().startsWith('hypothetical protein'):null):null,pseudogene:pk?flag(pv):false,pseudogene_note:pv&&typeof pv==='object'?(pv.description||pv.reason||JSON.stringify(pv)):typeof pv==='string'&&!['true','false','yes','no','0','1'].includes(pv.toLowerCase())?pv:'',truncated:String(f.truncated||''),aa:aa,nt:nt,aa_length:aa?aa.replace(/\*+$/,'').length:null,nt_length:nt?nt.length:null});
   });bin.features=rows;bin.details=details;bin.error='';bin.source='Bakta JSON opened locally for this selected bin';rebuild();chooseBin();status.textContent='Loaded '+rows.length.toLocaleString()+' real feature records into this browser session.';
  }catch(error){status.textContent='Could not open annotation: '+error.message;}this.value='';
 });
 rebuild();select.value='all';type.value='all';category.value='all';chooseBin();
})();
