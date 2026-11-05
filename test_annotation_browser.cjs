// Execute the embedded explorer against explicitly synthetic feature records.
const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const html=process.argv[2]?fs.readFileSync(process.argv[2],'utf8'):require('node:child_process').execFileSync('python3',[require('node:path').join(__dirname,'annotation_browser_fixture.py')],{encoding:'utf8'});
const data=JSON.parse(html.match(/id="annotation-feature-data">([\s\S]*?)<\/script>/)[1]);
const code=Array.from(html.matchAll(/<script>([\s\S]*?)<\/script>/g),m=>m[1]).find(s=>s.includes('function matchesQuery('));
const elements=new Map(),listeners={};let selectedTab=false,blob;
class Element{
 constructor(tag){this.tagName=tag;this.children=[];this.attributes={};this.listeners={};this.value='';this._text='';this.disabled=false;}
 set textContent(value){this._text=String(value);this.children=[];}
 get textContent(){return this._text;}
 get options(){return this.children;}
 appendChild(child){this.children.push(child);child.parentNode=this;if(this.tagName==='select'&&this.children.length===1)this.value=child.value;return child;}
 setAttribute(key,value){this.attributes[key]=String(value);}
 getAttribute(key){return this.attributes[key];}
 addEventListener(key,callback){this.listeners[key]=callback;}
 click(){this.listeners.click?.({target:this});}
 remove(){}
 scrollIntoView(){}
 showModal(){this.open=true;}
 close(){this.open=false;}
}
for(const m of html.matchAll(/<(\w+)[^>]*\bid="([^"]+)"[^>]*>/g))elements.set(m[2],new Element(m[1]));
elements.get('annotation-feature-data').textContent=JSON.stringify(data);
const document={getElementById:id=>elements.get(id),createElement:tag=>new Element(tag),querySelectorAll:()=>[],
 addEventListener:(name,fn)=>{listeners[name]=fn;},querySelector:()=>({click(){selectedTab=true;}}),body:new Element('body')};
vm.runInNewContext(code,{document,window:{MHCharts:{bar(){}}},console,Blob,setTimeout:fn=>fn(),URL:{createObjectURL(value){blob=value;return 'blob:test';},revokeObjectURL(){}}});
const el=id=>elements.get(id),trigger=(id,name)=>el(id).listeners[name].call(el(id)),rows=()=>el('annotation-feature-body').children;
(async()=>{
 assert.equal(rows().length,100);assert(el('annotation-feature-status').textContent.includes('207 matching features'));
 el('annotation-bin').value=data.bins[0].key;trigger('annotation-bin','change');
 assert(el('annotation-feature-status').textContent.includes('206 matching features'));
 trigger('annotation-next','click');assert(el('annotation-page').textContent.startsWith('Page 2 of 3'));
 trigger('annotation-next','click');assert.equal(rows().length,6);
 el('annotation-search').value='chromosomal replication';trigger('annotation-search','input');
 assert.equal(rows().length,1);assert.equal(rows()[0].children[2].textContent,'dnaA');
 rows()[0].children[1].children[0].click();assert(el('annotation-feature-dialog').open);assert.equal(el('annotation-detail-aa').textContent,'MK*');
 trigger('annotation-detail-protein','click');assert((await blob.text()).includes('>MH_1 Test Chromosomal replication initiator\nMK*\n'));
 trigger('annotation-detail-close','click');assert(!el('annotation-feature-dialog').open);
 el('annotation-search').value='COG:0593';trigger('annotation-search','input');assert.equal(rows().length,1);
 trigger('annotation-clear','click');el('annotation-type').value='trna';trigger('annotation-type','change');assert.equal(rows().length,1);
 el('annotation-type').value='cds';el('annotation-category').value='pseudogene';trigger('annotation-category','change');assert.equal(rows().length,2);
 el('annotation-search').value='hypothetical';trigger('annotation-search','input');assert.equal(rows().length,1);assert.equal(rows()[0].children[5].children[0].textContent,'Hypothetical');
 trigger('annotation-clear','click');el('annotation-category').value='pseudo-unknown';trigger('annotation-category','change');assert.equal(rows().length,1);assert.equal(rows()[0].children[6].textContent,'Unavailable');
 trigger('annotation-export','click');assert((await blob.text()).includes('unavailable'));
 trigger('annotation-clear','click');el('annotation-search').value='this does not exist';trigger('annotation-search','input');assert.equal(el('annotation-export').disabled,true);
 el('annotation-search').value='<img src=x';trigger('annotation-search','input');assert.equal(rows().length,1);assert(rows()[0].children[3].textContent.startsWith('<img src=x'));
 trigger('annotation-clear','click');el('annotation-search').value='dnaA';trigger('annotation-search','input');trigger('annotation-export','click');
 const tsv=await blob.text();assert(tsv.includes('dnaA'));assert(!tsv.includes('tRNA-Leu'));
 trigger('annotation-export-genes','click');assert((await blob.text()).includes('ATGAAATAA'));
 listeners.click({target:{closest(){return {getAttribute:()=>data.bins[1].key};}},preventDefault(){}});
 assert(selectedTab);assert.equal(el('annotation-bin').value,data.bins[1].key);assert.equal(rows().length,1);
 assert.equal(el('annotation-protein-download').attributes['data-sample-asset'],data.bins[1].protein_path);
 el('annotation-load-json').files=[{text:async()=>JSON.stringify({features:[{locus:'LOCAL',type:'cds',sequence:'contig_1',product:'hypothetical protein',pseudogene:{},aa:'MQQ',nt:'ATG'}]})}];
 await trigger('annotation-load-json','change');assert(el('annotation-import-status').textContent.startsWith('Loaded 1 real'));
 assert.equal(rows()[0].children[6].textContent,'Flagged');assert.equal(rows()[0].children[5].children[0].textContent,'Hypothetical');
 el('annotation-load-json').files=[{text:async()=>JSON.stringify({features:[{type:'cds',sequence:'foreign_contig'}]})}];
 await trigger('annotation-load-json','change');assert(el('annotation-import-status').textContent.includes('outside the selected bin'));assert.equal(rows()[0].children[1].children[0].textContent,'LOCAL');
 console.log('Explorer passed: bin links, independent categories, unknown flags, safe text, search, pagination, sequences, filtered exports, local JSON and contig validation.');
})().catch(e=>{console.error(e);process.exitCode=1;});
