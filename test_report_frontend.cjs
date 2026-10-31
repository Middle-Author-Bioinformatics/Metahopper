/* DOM/canvas logic regressions; no external browser or packages required. */
const fs=require('fs'),vm=require('vm'),assert=require('assert');
const reportPath=process.argv[2];
if(!reportPath)throw Error('Pass an updated HTML report path');
const html=fs.readFileSync(reportPath,'utf8');
assert(!html.includes('pane-expansion'),'Expansion must be integrated into Bins, not a separate pane');
for(const pane of ['overview','inspection','files'])assert(!html.includes('id="pane-'+pane+'"'),'unexpected pane '+pane);
assert(html.includes('id="pane-bins" class="pane on"'),'Bins must be the initial pane');
assert(!html.includes('>all ranks</option>'),'all ranks must not be selectable');
const script=html.match(/<script>([\s\S]*?)<\/script>/)[1];
const initial=JSON.parse(script.match(/var DATA = (.*?);\s*\/\*/s)[1]);
const elements=new Map(),events={};
class Element{
 constructor(tag){this.tagName=tag;this.children=[];this.style={};this.dataset={};this.attributes={};this.listeners={};this.value='';this.textContent='';this.disabled=false;this.files=[];this.clientWidth=1000;this.classList={toggle(){}};this._innerHTML='';}
 set innerHTML(value){this._innerHTML=value;this.children=[];}
 get innerHTML(){return this._innerHTML;}
 get options(){return this.children;}
 appendChild(child){this.children.push(child);child.parentNode=this;if(this.tagName==='select'&&(child.selected||this.children.length===1))this.value=child.value;return child;}
 remove(){if(this.parentNode)this.parentNode.children=this.parentNode.children.filter(c=>c!==this);}
 setAttribute(k,v){this.attributes[k]=String(v);if(k==='id')elements.set(v,this);}
 getAttribute(k){return this.attributes[k];}
 addEventListener(k,fn){this.listeners[k]=fn;}
 getBoundingClientRect(){return {left:0,top:0,width:1000,height:540};}
 setPointerCapture(){}
 focus(){this.focused=true;}
 showModal(){this.open=true;}
 close(){this.open=false;}
 click(){if(this.onclick)return this.onclick({target:this});}
 getContext(){return context2d;}
}
const canvasFills=[];
const context2d={globalAlpha:1,fillStyle:'',clearRect(){canvasFills.length=0;},arc(x,y,r){this.lastArc={x,y,r};},fill(){canvasFills.push({colour:this.fillStyle,alpha:this.globalAlpha,point:this.lastArc});},setTransform(){},beginPath(){},moveTo(){},lineTo(){},stroke(){},fillText(){},save(){},translate(){},rotate(){},restore(){},closePath(){},setLineDash(){}};
for(const match of html.matchAll(/<(\w+)\b([^>]*\bid="([^"]+)"[^>]*)>/g)){
 const element=new Element(match[1]);elements.set(match[3],element);
 const value=match[2].match(/\bvalue="([^"]*)"/);if(value)element.value=value[1];
 element.disabled=/\bdisabled\b/.test(match[2]);
}
const rankOptions=html.match(/<select id="rank-filter"[^>]*>([\s\S]*?)<\/select>/)[1];
for(const match of rankOptions.matchAll(/<option value="([^"]+)"([^>]*)>(.*?)<\/option>/g)){
 const option=new Element('option');option.value=match[1];option.selected=/\bselected\b/.test(match[2]);
 option.textContent=match[3];elements.get('rank-filter').appendChild(option);
}
function decodeHtml(text){return text.replace(/&quot;/g,'"').replace(/&#x27;/g,"'").replace(/&lt;/g,'<').replace(/&gt;/g,'>').replace(/&amp;/g,'&');}
const rankRows=Array.from(html.matchAll(/<tr data-rank="([^"]+)"([^>]*)>([\s\S]*?)<\/tr>/g),match=>{
 const row=new Element('tr');row.dataset.rank=match[1];
 row.dataset.binName=decodeHtml(match[2].match(/data-bin-name="([^"]*)"/)[1]);row.innerHTML=match[3];return row;
});
for(const id of ['x','y','xs','ys','c','bf','graph-hops','graph-limit','graph-labels'])elements.get(id).tagName='select';
elements.get('graph-hops').value='1';elements.get('graph-limit').value='500';elements.get('graph-labels').value='yes';
const document={getElementById:id=>elements.get(id)||null,createElement:tag=>new Element(tag),createElementNS:(ns,tag)=>new Element(tag),createTextNode:text=>({textContent:text}),querySelectorAll(query){return query==='#bin-table tbody tr'?rankRows:[];},body:new Element('body')};
const window={...(process.env.TEST_POINTER==='1'?{PointerEvent:function(){}}:{}),devicePixelRatio:1,addEventListener:(name,fn)=>{events[name]=fn;}};
let downloaded=0;
const sandbox={document,window,console,Blob,Set,Map,Math,JSON,Number,String,Array,Object,Infinity,parseFloat,parseInt,
 setTimeout(fn){fn();return 1;},clearTimeout(){},alert(message){throw Error(message);},
 URL:{createObjectURL(){downloaded++;return 'blob:test';},revokeObjectURL(){}},
 XMLSerializer:class{serializeToString(){return '<svg xmlns="http://www.w3.org/2000/svg"/>';}}};
vm.runInNewContext(script,sandbox,{timeout:10000});
const el=id=>elements.get(id);
function change(id,value){el(id).value=value;el(id).onchange();}
assert.deepEqual(el('tabs').children.map(button=>button.textContent),['Bins','BinaRena','Annotation']);
assert.equal(el('rank-filter').value,'phylum');
assert(rankRows.some(row=>row.dataset.rank==='phylum'&&row.style.display===''));
assert(rankRows.every(row=>(row.style.display==='')===(row.dataset.rank==='phylum')));
change('rank-filter','genus');
assert(rankRows.some(row=>row.dataset.rank==='genus'&&row.style.display===''));
assert(rankRows.every(row=>(row.style.display==='')===(row.dataset.rank==='genus')));
const originalMetrics=rankRows.map(row=>row.innerHTML);
function searchBins(query){el('bin-search').value=query;el('bin-search').listeners.input();return rankRows.filter(row=>row.style.display==='');}
let binMatches=searchBins('sOdAlIs');
assert(binMatches.length>0);
assert(binMatches.every(row=>row.dataset.binName.toLowerCase().includes('sodalis')));
assert(new Set(binMatches.map(row=>row.dataset.rank)).size>=2,'search remained constrained to the selected rank');
assert.equal(el('rank-filter').disabled,true);assert.equal(el('rank-filter').value,'genus');
assert(el('bin-search-status').textContent.includes(binMatches.length+' matching bins'));
const underscored=rankRows.find(row=>row.dataset.binName.includes('_'));assert(underscored);
binMatches=searchBins('  '+underscored.dataset.binName.replace(/_/g,'  ').toUpperCase()+'  ');
assert(binMatches.includes(underscored),'spaces/underscores/case were not normalized');
assert.equal(searchBins('this_taxon_does_not_exist_93264').length,0);
assert.equal(el('bin-search-empty').hidden,false);assert(el('bin-search-status').textContent.startsWith('0 matching bins'));
assert.equal(searchBins('<img src=x onerror=alert(1)>').length,0);
el('bin-search-clear').click();assert.equal(el('bin-search').value,'');assert.equal(el('rank-filter').disabled,false);
assert.equal(el('bin-search-empty').hidden,true);assert(el('bin-search').focused);
assert(rankRows.every(row=>(row.style.display==='')===(row.dataset.rank==='genus')));
searchBins('Sodalis');el('bin-search').listeners.keydown({key:'Escape'});
assert.equal(el('bin-search').value,'');assert.equal(el('rank-filter').disabled,false);
searchBins('   ');assert.equal(el('rank-filter').disabled,false);
assert(rankRows.every(row=>(row.style.display==='')===(row.dataset.rank==='genus')));
searchBins('Sodalis');searchBins('');
assert(rankRows.every(row=>(row.style.display==='')===(row.dataset.rank==='genus')));
assert.deepEqual(rankRows.map(row=>row.innerHTML),originalMetrics,'search changed bin metrics');
assert(html.match(/Classified bins<\/div><div class="v">[0-9,]+<\/div>/));
assert(html.match(/Assembly size<\/div><div class="v">[^<]+<\/div>/));
assert(html.includes('Contigs plotted</div><div class="v">'+initial.scatter.nShown.toLocaleString()+'</div>'));
function lassoAll(){
 const cv=el('cv');cv.listeners[process.env.TEST_POINTER==='1'?'pointerdown':'mousedown']({clientX:0,clientY:0});
 for(const [x,y] of [[1000,0],[1000,540],[0,540],[0,0]])cv.listeners[process.env.TEST_POINTER==='1'?'pointermove':'mousemove']({clientX:x,clientY:y});
 events[process.env.TEST_POINTER==='1'?'pointerup':'mouseup']({shiftKey:false});
}
function assertColours(field){
 change('c',field);
 const button=el('lg').children.find(x=>x.tagName==='button'&&x.textContent.startsWith('Show all'));if(button)button.click();
 const swatches=new Map(el('lg').children.filter(x=>x.dataset&&x.dataset.label).map(x=>[x.dataset.label,x.dataset.colour]));
 const cat=initial.scatter.categorical[field];
 assert.equal(canvasFills.length,initial.scatter.nShown);
 canvasFills.forEach((fill,i)=>{const label=cat.levels[cat.codes[i]];assert.equal(fill.colour,swatches.get(label),'dot/legend mismatch for '+label);assert.equal(fill.alpha,1);});
 if(field==='taxon_genus')assert.equal(new Set(swatches.values()).size,swatches.size,'repeated genus palette colours');
 return swatches;
}
const genus=assertColours('taxon_genus');
assertColours('bin');assertColours('taxon_species');
// A CSS-resized canvas must use the same point coordinates for draw and hover.
const hoverCanvas=el('cv');
const originalBounds=hoverCanvas.getBoundingClientRect;
hoverCanvas.getBoundingClientRect=()=>({left:100,top:40,width:500,height:260});
const dot=canvasFills[0].point;
hoverCanvas.listeners[process.env.TEST_POINTER==='1'?'pointermove':'mousemove']({clientX:100+dot.x/2,clientY:40+dot.y/2});
assert.equal(el('tip').style.display,'block');
assert(el('tip').innerHTML.includes(initial.scatter.ids[0]));
assert(el('tip').innerHTML.includes('coverage'));
hoverCanvas.listeners[process.env.TEST_POINTER==='1'?'pointerdown':'mousedown']({clientX:100,clientY:40});
hoverCanvas.listeners[process.env.TEST_POINTER==='1'?'lostpointercapture':'mouseleave']();
hoverCanvas.listeners[process.env.TEST_POINTER==='1'?'pointermove':'mousemove']({clientX:100+dot.x/2,clientY:40+dot.y/2});
assert.equal(el('tip').style.display,'block','leaving the iframe left lasso active and suppressed hover');
hoverCanvas.getBoundingClientRect=originalBounds;
change('c','taxon_genus');
const first=initial.scatter.ids[0],length=initial.scatter.numeric.length[0];
el('q').value=first;el('q').listeners.input();
assert(canvasFills.length>0);
lassoAll();assert.equal(el('graph-open').disabled,false);
const before=canvasFills[0].colour;
el('lbl').value='Manual test';el('assign').click();
assert.equal(el('c').value,'taxon_genus','manual assignment silently changed the colour mode');
assert.equal(canvasFills[0].colour,before,'manual assignment silently overrode taxonomy colour');
change('c','Manual sub-bin');
assert.equal(canvasFills[0].colour,el('lg').children[0].dataset.colour);
change('c','taxon_genus');lassoAll();el('graph-open').click();
assert(el('graph-status').textContent.includes(initial.assemblyGraphs.contigs[first]?'source graph(s)':'No source graph mapping'));
const catalogue={schema:'metahopper-graphs-1',stage:'final',run_name:initial.runName,lengths:{[first]:length},unmapped:{},warnings:[],
 graphs:[{id:'test/one.gfa',nodes:[{id:'1',length,depth:2},{id:'2',length:100,depth:1}],links:[{source:'1',source_orient:'+',target:'2',target_orient:'-',overlap:'0M',kind:'L'}]},
 {id:'test/two.gfa',nodes:[{id:'1',length:90,depth:1}],links:[]}],
 contigs:{[first]:[{graph:'test/one.gfa',steps:[['1','+']],method:'test_fixture'}]}};
async function finish(){
 el('graph-file').files=[{size:2000,text:async()=>JSON.stringify(catalogue)}];await el('graph-file').onchange.call(el('graph-file'));
 assert(el('graph-status').textContent.includes('1 source graph(s); 2 segments; 1 links.'));
 let svg=el('graph-stage').children[0];assert(svg,'graph SVG absent');
 const group=svg.children.find(x=>x.tagName==='g');
 const nodes=group.children.filter(x=>x.tagName==='g');assert.equal(nodes.length,2);
 assert(nodes.every(n=>n.attributes['data-source']==='test/one.gfa'),'sources were merged by numeric segment ID');
 const selected=nodes.find(n=>n.attributes['data-segment']==='1');assert.equal(selected.children[0].attributes.fill,before);
 change('graph-hops','0');assert(el('graph-status').textContent.includes('1 segments; 0 links.'));
 el('graph-download').click();assert.equal(downloaded,1);
 catalogue.run_name='wrong-run';await el('graph-file').onchange.call(el('graph-file'));
 assert(el('graph-status').textContent.includes('does not match this report'));
 el('graph-close').click();assert.equal(el('graph-dialog').open,false);
 console.log(JSON.stringify({contigs:initial.scatter.nShown,genus_colours:genus.size,checks:['exactly three tabs','Bins summary cards','Phylum default and rank switching','taxonomic-name search across ranks','case/space/underscore matching','no-match and literal-input handling','clear/Escape restore selected rank','bin metrics preserved by search','dot/legend agreement','opaque palette','unique genus colours','explicit manual colour mode','graph button selection','graph neighbourhood and strand links','source graph isolation','SVG download','wrong-run catalogue rejection']},null,2));
}
finish().catch(error=>{console.error(error);process.exitCode=1;});

