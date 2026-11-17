// Exercise the generated charts and real table filter together, without a browser.
const assert = require('node:assert/strict');
const vm = require('node:vm');
const {execFileSync} = require('node:child_process');

function fixture(ranks) {
  return JSON.parse(execFileSync('python3', ['-c', String.raw`
import json, re, sys
import metahopper_report as report
rows, partitions = [], {}
for i, rank in enumerate(json.loads(sys.argv[1])):
    bins = [dict(rank=rank, unclassified=False, bin='Shared_taxon' if j == 0 else rank + '_bin_' + str(j),
                 total_length_bp=(j + 1) * 1000, num_contigs=j + 1, N50=(j + 1) * 100)
            for j in range(i + 2)]
    rows.extend(bins)
    partitions[rank] = dict(parts=[dict(bin=b['bin'], bp=b['total_length_bp'],
                                       contigs=b['num_contigs'], unbinned=False) for b in bins]
                           + [dict(bin='Unbinned', bp=100000-sum(b['total_length_bp'] for b in bins),
                                   contigs=10, unbinned=True)])
html = report.render_bins(rows, report.rank_totals(rows),
                         partition=dict(total_bp=100000, source='synthetic.fasta', ranks=partitions))
table_js = re.search(r'/\* ---------------- bins tab rank filter.*?\*/(.*?)/\* ---------------- BinaRena',
                     report.HTML, re.S).group(1)
print(json.dumps(dict(html=html, table_js=table_js)))
`, JSON.stringify(ranks)], {cwd:__dirname, encoding:'utf8'}));
}

class Event {
  constructor(type, options={}) { this.type=type; Object.assign(this, options); }
  preventDefault() {}
}
class Element {
  constructor(tag) {
    this.tagName=tag; this.children=[]; this.attributes={}; this.listeners={};
    this.style={}; this.dataset={}; this.value=''; this._text=''; this.classList={add(){}};
  }
  set textContent(value) { this._text=String(value); this.children=[]; }
  get textContent() { return this._text+this.children.map(c=>c.textContent).join(''); }
  get options() { return this.children; }
  appendChild(child) {
    this.children.push(child);
    if(this.tagName==='select' && (child.selected || this.children.length===1))this.value=child.value;
    return child;
  }
  setAttribute(key,value) { this.attributes[key]=String(value); }
  addEventListener(type,callback) { (this.listeners[type] ||= []).push(callback); }
  dispatchEvent(event) {
    event.target=this;
    this['on'+event.type]?.(event);
    for(const callback of this.listeners[event.type] || [])callback(event);
  }
  focus() { this.focused=true; }
  scrollIntoView() { this.scrolled=true; }
}

function reportDOM(ranks) {
  const {html,table_js}=fixture(ranks), elements=new Map(), doc=new Element('document');
  for(const m of html.matchAll(/<(\w+)\b[^>]*\bid="([^"]+)"[^>]*>/g))elements.set(m[2],new Element(m[1]));
  for(const m of html.matchAll(/<script type="application\/json" id="([^"]+)">([\s\S]*?)<\/script>/g))elements.get(m[1]).textContent=m[2];
  const rankOptions=html.match(/<select id="rank-filter"[^>]*>([\s\S]*?)<\/select>/)[1];
  for(const m of rankOptions.matchAll(/<option value="([^"]+)"([^>]*)>(.*?)<\/option>/g)) {
    const option=new Element('option');option.value=m[1];option.selected=m[2].includes('selected');
    option.textContent=m[3];elements.get('rank-filter').appendChild(option);
  }
  const tableRows=Array.from(html.matchAll(/<tr data-rank="([^"]+)" data-bin-name="([^"]+)"/g),m=>{
    const row=new Element('tr');row.dataset={rank:m[1],binName:m[2]};return row;
  });
  Object.assign(doc, {
    getElementById:id=>elements.get(id), createElement:tag=>new Element(tag),
    createElementNS:(ns,tag)=>new Element(tag), querySelectorAll:()=>tableRows,
  });
  const context=vm.createContext({document:doc,window:{},Event,CustomEvent:Event});
  for(const m of html.matchAll(/<script>([\s\S]*?)<\/script>/g))vm.runInContext(m[1],context);
  vm.runInContext(table_js,context);
  const el=id=>elements.get(id);
  const change=(id,value)=>{el(id).value=value;el(id).dispatchEvent(new Event('change'));};
  const search=value=>{el('bin-search').value=value;el('bin-search').dispatchEvent(new Event('input'));};
  const visibleRows=()=>tableRows.filter(r=>r.style.display==='');
  const svg=id=>el(id).children[0];
  const circles=()=>svg('bin-contiguity-figure').children.filter(n=>n.tagName==='circle');
  const flows=()=>svg('bin-partition-figure').children.filter(n=>n.attributes.class==='flow');
  return {html,el,change,search,visibleRows,svg,circles,flows};
}

const ranks=['phylum','class','order','family','genus','species'];
const report=reportDOM(ranks);
assert(!report.html.includes('bin-figure-rank'),'plots must not introduce another rank selector');
assert.deepEqual(report.el('rank-filter').options.map(o=>o.value),ranks);
assert.equal(report.el('rank-filter').value,'phylum');
for(const [i,rank] of ranks.entries()) {
  // One rank control updates the table, both titles and actual SVG data.
  report.change('rank-filter',rank);
  assert.equal(report.el('rank-filter').value,rank);
  assert.equal(report.visibleRows().length,i+2);
  assert(report.visibleRows().every(r=>r.dataset.rank===rank));
  assert.equal(report.circles().length,i+2);
  assert.equal(report.flows().length,i+3);
  assert(report.svg('bin-partition-figure').attributes['aria-label'].includes(rank+' rank'));
  for(const id of ['bin-partition-title','bin-contiguity-title'])assert(report.el(id).textContent.endsWith(rank[0].toUpperCase()+rank.slice(1)));
  assert.equal(report.el('bin-partition-total').textContent,'100 kb');
  const svg=report.svg('bin-partition-figure'), height=Number(svg.attributes.viewBox.split(' ')[3]);
  const texts=report.flows().flatMap(g=>g.children.filter(c=>c.tagName==='text'));
  assert(texts.every(n=>Number(n.attributes.y)+4<height),'SVG viewBox clips a destination label');
  const black=report.flows().at(-1);
  assert.equal(black.children.find(n=>n.tagName==='path').attributes.fill,'#000000');
  const heights=report.flows().map(g=>Number(g.children.find(n=>n.tagName==='rect').attributes.height));
  assert(Math.abs(heights.reduce((a,b)=>a+b,0)-280)<1e-8,'ribbons must conserve the assembly at every rank');
}

report.change('rank-filter','genus');
report.search('sHaReD tAxOn');
assert.equal(report.visibleRows().length,1,'table search must use the shared rank');
assert.equal(report.circles().length,1,'scatter search must stay at the selected rank');
assert.equal(report.flows().length,7,'search must retain the complete partition');
report.change('rank-filter','species');
assert.equal(report.flows().length,8,'plot rank must remain selectable during search');
report.search('no such taxon');
assert(report.el('bin-contiguity-figure').textContent.includes('No paired values'));
assert.equal(report.flows().length,8);
report.el('bin-search').dispatchEvent(new Event('keydown',{key:'Escape'}));
assert.equal(report.circles().length,7,'Escape should refresh the figures');
report.search('Shared');
report.el('bin-search-clear').dispatchEvent(new Event('click'));
assert.equal(report.circles().length,7,'Clear should refresh the figures');

// The dashboard focuses a bin by assigning the rank, then dispatching search input.
report.el('rank-filter').value='genus';report.search('Shared');
assert.equal(report.el('rank-filter').value,'genus');
assert(report.el('bin-partition-title').textContent.endsWith('Genus'));
report.flows()[1].dispatchEvent(new Event('click'));
assert.equal(report.visibleRows().length,1);
assert.equal(report.circles().length,1);
assert(report.el('bin-table').scrolled);
assert.equal(report.flows().length,7);

// No assumptions about Phylum being present in a saved report.
const genusOnly=reportDOM(['genus']);
assert.equal(genusOnly.el('rank-filter').value,'genus');
assert.equal(genusOnly.circles().length,2);
assert(report.html.includes('#bin-figures{grid-template-columns:minmax(0,1fr)}'));
assert(report.html.includes('#bin-figures .figure-host{max-height:none;overflow:visible}'));
console.log('Bin figures passed: six ranks, synchronized controls/table/SVGs, search, Clear/Escape, dashboard focus, flow clicks, complete partitions and unclipped labels.');
