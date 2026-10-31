/* Exercise the real report's lazy graph and FASTA paths with a parent asset bridge. */
const fs=require('fs'),vm=require('vm'),path=require('path');
const harness=fs.readFileSync(path.join(__dirname,'test_report_frontend.cjs'),'utf8');
const setup=harness.slice(0,harness.indexOf('vm.runInNewContext(script,sandbox'));
const body=String.raw`
const first=initial.scatter.ids[0],length=initial.scatter.numeric.length[0];
const catalog={schema:'metahopper-graphs-1',stage:'final',run_name:initial.runName,lengths:{[first]:length},unmapped:{},warnings:[],
 graphs:[{id:'megahit/assembly.gfa',nodes:[{id:'k3_1',length},{id:'k3_2',length:10}],links:[{source:'k3_1',source_orient:'+',target:'k3_2',target_orient:'-',kind:'L',overlap:'3M'}]}],
 contigs:{[first]:[{graph:'megahit/assembly.gfa',steps:[['k3_1','+']],method:'source_segment_name'}]}};
const manifest=new Element('script');manifest.textContent=JSON.stringify({files:{'assembly_graphs.final.json':{},'final/assembly/consolidated_contigs.fasta':{}}});elements.set('metahopper-assets',manifest);
const requests=[],downloads=[];
window.metahopperLoadAsset=async name=>{requests.push(name);return new TextEncoder().encode(name.endsWith('.json')?JSON.stringify(catalog):'>'+first+'\n'+'A'.repeat(length)+'\n').buffer;};
sandbox.TextDecoder=TextDecoder;sandbox.File=File;
sandbox.FileReader=class{readAsText(blob){blob.text().then(text=>{this.result=text;this.onload();}).catch(error=>this.onerror(error));}};
sandbox.URL.createObjectURL=blob=>{downloads.push(blob);return 'blob:test';};
vm.runInNewContext(script,sandbox,{timeout:10000});
const el=id=>elements.get(id);
el('q').value=first;el('q').listeners.input();
const canvas=el('cv');
function select(){canvas.listeners.mousedown({clientX:0,clientY:0});for(const [x,y]of[[1000,0],[1000,540],[0,540],[0,0]])canvas.listeners.mousemove({clientX:x,clientY:y});events.mouseup({shiftKey:false});}
(async()=>{
 select();await el('graph-open').click();
 assert(el('graph-status').textContent.includes('2 segments; 1 links.'));
 el('graph-close').click();
 el('lbl').value='Cloud test';el('assign').click();
 assert.equal(el('expfa').disabled,false,'cloud FASTA should be available without a local file');
 await el('expfa').click();
 for(let n=0;n<100&&!downloads.length;n++)await new Promise(resolve=>setImmediate(resolve));
 assert.equal(downloads.length,1);
 const lines=(await downloads[0].text()).trim().split('\n');
 assert.equal(lines.shift(),'>'+first);const sequence=lines.join('');assert.equal(sequence.length,length);assert(/^A+$/.test(sequence));
 assert.deepEqual(requests,['assembly_graphs.final.json','final/assembly/consolidated_contigs.fasta']);
 console.log('Cloud graph and FASTA report paths passed for '+initial.runName);
})().catch(error=>{console.error(error);process.exitCode=1;});
`;
vm.runInNewContext(setup+body,{require,process,console,Blob,File,TextEncoder,TextDecoder,setImmediate},{timeout:10000,filename:__filename});
