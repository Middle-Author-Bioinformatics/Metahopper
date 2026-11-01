// Executes the actual report browser script with deterministic feature records.
const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const html=process.argv[2]?fs.readFileSync(process.argv[2],'utf8'):require('node:child_process').execFileSync('python3',[require('node:path').join(__dirname,'annotation_browser_fixture.py')],{encoding:'utf8'});
const data=JSON.parse(html.match(/id="annotation-feature-data">([\s\S]*?)<\/script>/)[1]);
const code=Array.from(html.matchAll(/<script>([\s\S]*?)<\/script>/g),m=>m[1]).find(s=>s.includes('function matchesQuery('));
const elements=new Map(),listeners={};let selectedTab=false,blob;
class Element{
 constructor(tag){this.tagName=tag;this.children=[];this.attributes={};this.listeners={};this.value='';this._text='';this.disabled=false;}
 set textContent(value){this._text=String(value);this.children=[];}
 get textContent(){return this._text;}
 appendChild(child){this.children.push(child);child.parentNode=this;return child;}
 setAttribute(key,value){this.attributes[key]=String(value);}
 getAttribute(key){return this.attributes[key];}
 addEventListener(key,callback){this.listeners[key]=callback;}
 click(){this.listeners.click?.({target:this});}
 remove(){}
 focus(){this.focused=true;}
 scrollIntoView(){this.scrolled=true;}
}
for(const m of html.matchAll(/<(\w+)[^>]*\bid="([^"]+)"[^>]*>/g))elements.set(m[2],new Element(m[1]));
elements.get('annotation-feature-data').textContent=JSON.stringify(data);
elements.get('annotation-type').value='proteins';
const document={getElementById:id=>elements.get(id),createElement:tag=>new Element(tag),
 addEventListener:(name,fn)=>{listeners[name]=fn;},querySelector:()=>({click(){selectedTab=true;}}),body:new Element('body')};
vm.runInNewContext(code,{document,console,Blob,setTimeout:fn=>fn(),URL:{createObjectURL(value){blob=value;return 'blob:test';},revokeObjectURL(){}}});
const el=id=>elements.get(id), trigger=(id,name)=>el(id).listeners[name]();
assert.equal(el('annotation-feature-body').children.length,100);
assert(el('annotation-feature-status').textContent.includes('205 matching features'));
assert.equal(el('annotation-next').disabled,false);trigger('annotation-next','click');
assert(el('annotation-page').textContent.startsWith('Page 2 of 3'));
trigger('annotation-next','click');assert.equal(el('annotation-feature-body').children.length,5);
el('annotation-search').value='chromosomal replication';trigger('annotation-search','input');
assert.equal(el('annotation-feature-body').children.length,1);
assert.equal(el('annotation-feature-body').children[0].children[1].textContent,'dnaA');
assert(el('annotation-page').textContent.startsWith('Page 1 of 1'));
el('annotation-search').value='COG:0593';trigger('annotation-search','input');assert.equal(el('annotation-feature-body').children.length,1);
trigger('annotation-clear','click');el('annotation-type').value='trna';trigger('annotation-type','change');assert.equal(el('annotation-feature-body').children.length,1);
el('annotation-type').value='all';trigger('annotation-type','change');assert(el('annotation-feature-status').textContent.includes('206 matching features'));
el('annotation-search').value='this does not exist';trigger('annotation-search','input');assert.equal(el('annotation-export').disabled,true);
el('annotation-search').value='<img src=x';trigger('annotation-search','input');
assert.equal(el('annotation-feature-body').children.length,1);
assert(el('annotation-feature-body').children[0].children[2].textContent.startsWith('<img src=x'));
trigger('annotation-clear','click');el('annotation-search').value='dnaA';trigger('annotation-search','input');trigger('annotation-export','click');
assert(blob instanceof Blob);
listeners.click({target:{closest(){return {getAttribute:()=>data.bins[1].key};}},preventDefault(){}});
assert(selectedTab);assert.equal(el('annotation-bin').value,data.bins[1].key);
assert.equal(el('annotation-feature-body').children.length,1);assert(el('annotation-browser-title').textContent.includes('Second'));
assert.equal(el('annotation-protein-download').attributes['data-sample-asset'],data.bins[1].protein_path);
blob.text().then(text=>{assert(text.includes('dnaA'));assert(!text.includes('tRNA-Leu'));console.log('Gene browser: bin links, search, feature filters, pagination, TSV export, safe text and sample downloads passed.');}).catch(e=>{console.error(e);process.exitCode=1;});
