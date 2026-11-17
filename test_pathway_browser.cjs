// Run the generated pathway UI with explicitly synthetic results and no browser.
const assert=require('node:assert/strict'),vm=require('node:vm');
const html=require('node:child_process').execFileSync('python3',['-c',String.raw`
import tempfile
from pathlib import Path
from pathway_test_fixture import fixture
from metahopper_pathway_view import read_results, render_pathways
with tempfile.TemporaryDirectory() as tmp:
    root=Path(tmp)/'Sample'
    print(render_pathways(read_results(fixture(root),root)))
`],{cwd:__dirname,encoding:'utf8'});
assert(!html.includes('id="pathway-component-body"'));
assert(!html.includes('Gene neighborhoods and reference-operon support'));
assert(!html.includes('Filter file names'));
for(const header of ['Assay scope','Scored panel','Selected route','Evidence'])assert(!html.includes('>'+header+'</th>'));
const data=JSON.parse(html.match(/id="metahopper-pathways-data">([\s\S]*?)<\/script>/)[1]);
const code=Array.from(html.matchAll(/<script>([\s\S]*?)<\/script>/g),m=>m[1])[0];
const elements=new Map();let downloaded,filename;
class Element {
  constructor(tag){this.tagName=tag;this.children=[];this.attributes={};this.listeners={};this.dataset={};this.style={};this.value='';this._text='';this.nodeType=1;}
  set textContent(value){this._text=String(value);this.children=[];}
  get textContent(){return this._text+this.children.map(c=>c.textContent).join('');}
  appendChild(child){child.parent=this;this.children.push(child);if(this.tagName==='select'&&this.children.length===1)this.value=child.value;return child;}
  get options(){return this.children;}
  setAttribute(key,value){this.attributes[key]=String(value);}
  getAttribute(key){return this.attributes[key];}
  addEventListener(type,callback){(this.listeners[type]||=[]).push(callback);}
  click(){if(this.tagName==='a')filename=this.download;for(const callback of this.listeners.click||[])callback();}
  closest(tag){return this.tagName===tag?this:this.parent?.closest(tag);}
  querySelectorAll(){return this.children.flatMap(child=>child.tagName==='th'?[child]:child.querySelectorAll());}
  showModal(){this.open=true;}
  close(){this.open=false;}
  scrollIntoView(){this.scrolled=true;}
  remove(){}
}
for(const m of html.matchAll(/<(\w+)\b([^>]*\bid="([^"]+)"[^>]*)>/g)){
  const node=new Element(m[1]);node.checked=/\bchecked\b/.test(m[2]);elements.set(m[3],node);
}
elements.get('metahopper-pathways-data').textContent=JSON.stringify(data);
for(const m of html.matchAll(/<table\b[^>]*>([\s\S]*?)<tbody id="([^"]+)"/g)){
  const table=new Element('table'),head=new Element('thead');
  for(const th of m[1].matchAll(/<th\b[^>]*>(.*?)<\/th>/g)){const node=new Element('th');node.textContent=th[1];head.appendChild(node);}
  table.appendChild(head);table.appendChild(elements.get(m[2]));
}
const document={getElementById:id=>elements.get(id),createElement:tag=>new Element(tag),
  createTextNode:value=>{const node=new Element('#text');node.textContent=value;return node;},
  querySelectorAll:()=>[],body:new Element('body')};
vm.runInNewContext(code,{document,window:{},console,Blob,URL:{createObjectURL(blob){downloaded=blob;return 'blob:test';},revokeObjectURL(){}},setTimeout:fn=>fn()});
const el=id=>elements.get(id),rows=id=>el('pathway-'+id+'-body').children;
const trigger=(id,event='change')=>{for(const callback of el(id).listeners[event]||[])callback();};
(async()=>{
  assert.equal(rows('summary').length,3); // The all-zero bin is hidden; failed/stale remain explicit.
  assert.equal(rows('hit').length,2);assert.equal(rows('hit')[0].children.length,8);
  assert(rows('hit')[1].children[5].textContent.includes('Shared family'));
  assert(!rows('hit').some(row=>row.children[0].textContent==='Stale'));
  rows('summary')[0].children[5].children[0].click();
  assert(el('pathway-panel-dialog').open);assert.equal(rows('panel').length,3);
  assert(rows('panel')[0].children[1].textContent.includes('SYNTHETIC_c1'));
  assert(rows('panel')[2].textContent.includes('Not measured'));
  el('pathway-panel-close').click();assert(!el('pathway-panel-dialog').open);
  rows('hit')[0].children[1].children[0].children[0].click();assert(el('pathway-panel-dialog').open);
  rows('hit')[0].children[3].children[0].click();assert(el('pathway-protein-dialog').open);
  assert.equal(el('pathway-protein-aa').textContent,'MK*');
  el('pathway-protein-aa-download').click();assert((await downloaded.text()).includes('MK*'));
  el('pathway-protein-close').click();
  el('pathway-search').value='g1';trigger('pathway-search','input');assert.equal(rows('hit').length,1);
  el('pathway-export-hits').click();assert.equal(filename,'hmm_hits.tsv');
  const tsv=await downloaded.text();assert(tsv.startsWith('bin\tpathway\tgene\tproduct\tprotein\thmm\tmatch\tevalue\thmm_coverage'));
  assert.equal(tsv.trim().split('\n').length,2);assert(tsv.includes('SYNTHETIC_c1'));assert(!tsv.includes('SYNTHETIC_c2'));
  el('pathway-export-summary').click();assert.equal(filename,'pathways.tsv');assert((await downloaded.text()).includes('coverage_percent'));
  el('pathway-export-panels').click();assert.equal(filename,'pathway_panels.tsv');assert((await downloaded.text()).includes('SYNTHETIC_c1'));
  el('pathway-reset').click();el('pathway-bin').value='Zero';trigger('pathway-bin');
  assert.equal(rows('summary')[0].children[2].textContent,'0%');assert.equal(rows('hit').length,0);
  el('pathway-hit-filter').value='review';trigger('pathway-hit-filter');assert.equal(rows('hit').length,1);
  assert(rows('hit')[0].children[5].textContent.includes('Context candidate'));
  el('pathway-bin').value='Stale';trigger('pathway-bin');assert.equal(rows('hit').length,0);
  assert(rows('summary')[0].children[4].textContent.includes('Stale input'));
  el('pathway-reset').click();el('pathway-search').value='Synthetic panel';trigger('pathway-search','input');
  assert.equal(rows('hit').length,2,'pathway-name searches must retain their gene hits');
  console.log('Pathway view passed: one simple hit table, panel links, three TSV exports, protein details, filters and stale/failed/zero handling.');
})().catch(error=>{console.error(error);process.exitCode=1;});
