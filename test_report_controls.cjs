/* Exercise real report controls against a small DOM, without network or a browser. */
const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm');
const script=fs.readFileSync(__dirname+'/metahopper_report_controls.js','utf8');
class Event {constructor(type,options={}){this.type=type;this.bubbles=!!options.bubbles;this.detail=options.detail;}preventDefault(){this.defaultPrevented=true;}}
let document;
function match(el,selector){
 const tag=selector.match(/^[a-z][a-z0-9-]*/i);if(tag&&el.tagName!==tag[0].toUpperCase())return false;
 const id=selector.match(/#([\w-]+)/);if(id&&el.id!==id[1])return false;
 for(const m of selector.matchAll(/\.([\w-]+)/g))if(!el.className.split(/\s+/).includes(m[1]))return false;
 for(const m of selector.matchAll(/\[([\w-]+)(?:(\^?=)"([^"]*)")?\]/g)){const value=el.getAttribute(m[1]);if(value===null)return false;if(m[2]==='='&&value!==m[3])return false;if(m[2]==='^='&&!value.startsWith(m[3]))return false;}
 return true;
}
class Element {
 constructor(tag){this.tagName=tag.toUpperCase();this.nodeType=1;this.children=[];this.attributes={};this.dataset={};this.listeners={};this.style={};this._text='';this.value='';this.className='';this.disabled=false;this.hidden=false;this.offsetWidth=280;this.offsetHeight=210;}
 get parentElement(){return this.parentNode||null;}get firstChild(){return this.children[0]||null;}
 get textContent(){return this._text+this.children.map(c=>c.textContent).join('');}set textContent(v){this._text=String(v);this.children=[];}
 get classList(){return {add:(s)=>{if(!this.className.split(/\s+/).includes(s))this.className=(this.className+' '+s).trim();},remove:(s)=>{this.className=this.className.split(/\s+/).filter(c=>c!==s).join(' ');},toggle:(s,v)=>{if(v)this.classList.add(s);else this.classList.remove(s);},contains:(s)=>this.className.split(/\s+/).includes(s)};}
 appendChild(c){c.parentNode=this;this.children.push(c);return c;}
 insertBefore(c,before){c.parentNode=this;const i=this.children.indexOf(before);if(i<0)this.children.push(c);else this.children.splice(i,0,c);return c;}
 remove(){if(this.parentNode)this.parentNode.children=this.parentNode.children.filter(c=>c!==this);this.parentNode=null;}
 setAttribute(k,v){if(k==='id')this.id=String(v);else if(k==='class')this.className=String(v);else if(k==='hidden')this.hidden=true;else this.attributes[k]=String(v);}
 getAttribute(k){if(k==='id')return this.id??null;if(k==='class')return this.className;if(k.startsWith('data-')){const key=k.slice(5).replace(/-([a-z])/g,(_,c)=>c.toUpperCase());if(key in this.dataset)return this.dataset[key];}return this.attributes[k]??null;}
 querySelectorAll(selector){const nodes=[];function walk(n){n.children.forEach(c=>{nodes.push(c);walk(c);});}walk(this);return nodes.filter(el=>selector.split(',').some(s=>{const parts=s.trim().split(/\s+/);if(!match(el,parts.pop()))return false;let parent=el.parentElement;while(parts.length){const part=parts.pop();while(parent&&!match(parent,part))parent=parent.parentElement;if(!parent)return false;parent=parent.parentElement;}return true;}));}
 querySelector(s){return this.querySelectorAll(s)[0]||null;}
 closest(s){for(let n=this;n;n=n.parentElement)if(match(n,s))return n;return null;}
 addEventListener(t,cb){(this.listeners[t]||=[]).push(cb);}
 dispatchEvent(e){e.target=this;document.emit(e,true);(this.listeners[e.type]||[]).forEach(cb=>cb(e));this['on'+e.type]?.(e);if(e.bubbles)document.emit(e,false);return true;}
 click(){if(!this.disabled)this.dispatchEvent(new Event('click',{bubbles:true}));}
 getBoundingClientRect(){return {left:100,top:190,right:112,bottom:202};}
 scrollIntoView(){this.scrolled=true;}
}
function fixture(support){
 const body=new Element('body'),listeners={},microtasks=[];
 document={body,readyState:'complete',createElement:tag=>new Element(tag),querySelector:s=>body.querySelector(s),querySelectorAll:s=>body.querySelectorAll(s),getElementById:id=>body.querySelector('#'+id),addEventListener:(t,cb,capture=false)=>(listeners[t]||=[]).push({cb,capture}),emit(e,capture){(listeners[e.type]||[]).filter(x=>x.capture===capture).forEach(x=>x.cb(e));},dispatchEvent(e){e.target=this;this.emit(e,true);this.emit(e,false);}};
 const add=(parent,tag,id,text)=>{const el=new Element(tag);el.id=id;el.textContent=text||'';parent.appendChild(el);return el;};
 const host=add(body,'div');host.className='wrap';const tabs=add(host,'div','tabs');
 const panes={};for(const name of ['bins','binarena','annotation']){const pane=add(host,'div','pane-'+name);pane.className='pane'+(name==='bins'?' on':'');panes[name]=pane;const b=add(tabs,'button',null,name);b.className='tab'+(name==='bins'?' on':'');b.dataset.id=name;b.onclick=()=>{tabs.children.forEach(t=>t.classList.toggle('on',t===b));Object.entries(panes).forEach(([id,p])=>p.classList.toggle('on',id===name));};}
 const controls={};for(const [id,value]of Object.entries({'rank-filter':'genus','bin-search':'','bf':'*','q':'','annotation-bin':'all','annotation-search':'','annotation-type':'all','annotation-category':'all','pathway-bin':'all','pathway-pathway':'all','pathway-search':'','pathway-hit-filter':'overview'})){const el=add(panes.bins,'input',id);el.value=value;controls[id]=el;}
 const hide=add(panes.annotation,'input','pathway-hide-zero');hide.type='checkbox';hide.checked=false;
 const binTable=add(panes.bins,'table','bin-table');let visible=['Taxon_A','Taxon_B'];controls['bin-search'].addEventListener('input',()=>{visible=['Taxon_A','Taxon_B'].filter(n=>n.includes(controls['bin-search'].value));});
 const genes=add(panes.annotation,'section','annotation-genes-section');const genesTab=add(panes.annotation,'button','annotation-genes-tab');genesTab.setAttribute('role','tab');genesTab.setAttribute('aria-selected','true');genesTab.onclick=()=>{genes.hidden=false;genesTab.setAttribute('aria-selected','true');};
 const originalReset=add(genes,'button','annotation-clear');originalReset.onclick=()=>{controls['annotation-search'].value='';controls['annotation-type'].value='all';controls['annotation-category'].value='all';};
 const pathwayReset=add(panes.annotation,'button','pathway-reset');pathwayReset.onclick=()=>{controls['pathway-bin'].value='all';controls['pathway-pathway'].value='all';controls['pathway-search'].value='';controls['pathway-hit-filter'].value='overview';hide.checked=false;};
 let exports=0;add(genes,'button','annotation-export').onclick=()=>exports++;
 const figure=add(panes.bins,'div','bin-contiguity-figure'),svg=add(figure,'svg');
 const data=[{rank:'genus',bin:'Taxon_A',total_length_bp:1200000,N50:15000,num_contigs:42,GC_percent:0,completeness_percent:97,contamination_percent:0},{rank:'genus',bin:'Taxon_B',total_length_bp:400000,N50:2000,num_contigs:400}];
 const points=data.map(r=>{const p=add(svg,'circle');add(p,'title',null,r.rank+' / '+r.bin.replace(/_/g,' ')+'\nBin size');return p;});
 add(body,'script','bin-figure-data',JSON.stringify(data));
 add(body,'script','metahopper-assets',JSON.stringify({files:{'report_support/run_info.json':{},'report_support/full_run.log':{}}}));
 if(support!==undefined)add(body,'script','metahopper-run-downloads',JSON.stringify(support));
 let selected=['c1'],manual={c1:'Kept manual bin'},resetEvents=0;
 const win={innerWidth:500,innerHeight:500,scrollY:120,scrollTo:({top})=>win.scrollY=top,MHBinarenaView:{capture:()=>selected.slice(),restore:ids=>selected=ids.slice(),reset:()=>selected=[]}};
 document.addEventListener('metahopper:reset-view',()=>resetEvents++);
 const context=vm.createContext({document,window:win,Event,CustomEvent:Event,queueMicrotask:cb=>microtasks.push(cb),Date});vm.runInContext(script,context);
 function flush(){let n=0;while(microtasks.length){assert(n++<60,'event loop failed to settle');microtasks.shift()();}}
 function change(id,value,type='input'){controls[id].value=value;controls[id].dispatchEvent(new Event(type,{bubbles:true}));flush();}
 return {body,host,controls,points,context,win,flush,change,visible:()=>visible.slice(),selected:()=>selected,manual,resetEvents:()=>resetEvents,exports:()=>exports,E:id=>document.getElementById(id)};
}
{
 const f=fixture();assert(f.E('mh-report-toolbar'));assert.equal(f.E('mh-back-view').disabled,true);
 f.points[0].dispatchEvent(new Event('pointerenter'));const tip=f.E('mh-bin-tooltip');assert.equal(tip.hidden,false);
 for(const text of ['Taxon A','1.2 Mb','15 kb','42','97%','0%'])assert(tip.textContent.includes(text),text);
 assert.equal(f.points[0].getAttribute('tabindex'),'0');assert.equal(f.points[0].querySelector('title'),null);
 f.points[0].dispatchEvent(new Event('pointerleave'));assert.equal(tip.hidden,true);
 f.points[0].dispatchEvent(new Event('focus'));assert.equal(tip.hidden,false);
 f.points[0].click();f.flush();assert.deepEqual(f.visible(),['Taxon_A']);assert.equal(f.E('mh-back-view').disabled,false);assert(f.E('bin-table').scrolled);
 f.win.MHReportControls.back();f.flush();assert.deepEqual(f.visible(),['Taxon_A','Taxon_B']);assert.equal(f.E('mh-back-view').disabled,true);
 vm.runInContext(script,f.context);assert.equal(f.body.querySelectorAll('#mh-report-toolbar').length,1,'viewer plus standalone scripts duplicate controls');
 const keys=['bin-search','q','annotation-bin','annotation-search','annotation-type','annotation-category','pathway-bin','pathway-pathway','pathway-search','pathway-hit-filter'];
 const values=['Taxon_A','c1','genus/Taxon_A','dnaA','cds','pseudogene','Taxon_A','biotin','bioB','missing'];
 keys.forEach((id,i)=>f.change(id,values[i],id.includes('search')||id==='q'?'input':'change'));
 f.E('pathway-hide-zero').checked=true;f.E('pathway-hide-zero').dispatchEvent(new Event('change',{bubbles:true}));f.flush();
 f.win.MHReportControls.reset();f.flush();assert.equal(f.controls['rank-filter'].value,'genus');assert.deepEqual(f.visible(),['Taxon_A','Taxon_B']);assert.equal(f.controls['annotation-bin'].value,'all');assert.equal(f.controls['pathway-bin'].value,'all');assert.equal(f.E('pathway-hide-zero').checked,false);assert.equal(f.selected().length,0);assert.equal(f.manual.c1,'Kept manual bin');assert.equal(f.resetEvents(),1);
 f.win.MHReportControls.back();f.flush();assert.equal(f.controls['annotation-bin'].value,'genus/Taxon_A');assert.equal(f.controls['pathway-bin'].value,'Taxon_A');assert.equal(f.controls['pathway-hit-filter'].value,'missing');assert.deepEqual(f.selected(),['c1']);assert.equal(f.manual.c1,'Kept manual bin');
 assert.equal(f.E('mh-report-toolbar').querySelectorAll('a[data-sample-asset]').length,2,'run downloads missing');
}
{
 const f=fixture({info:'report_support/run_info.json',log:null});
 const links=f.E('mh-report-toolbar').querySelectorAll('a[data-sample-asset]');
 assert.equal(links.length,1,'missing saved logs must not expose a Full log link');
 assert.equal(links[0].textContent,'Run info');
}
console.log('Report controls passed: hover/focus facts, point selection and Back, full filter reset, rank and manual assignments preserved, restored filters and selection, idempotent setup and scoped run downloads.');
