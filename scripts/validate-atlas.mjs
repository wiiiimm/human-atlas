import fs from 'node:fs';
import assert from 'node:assert/strict';
const filename=process.argv[2]??'atlas.json',female=filename.startsWith('atlas-female'),reconstructed=filename==='atlas-female-reconstructed.json';
const base=new URL('../public/models/',import.meta.url),atlas=JSON.parse(fs.readFileSync(new URL(filename,base)));
assert.equal(atlas.parts.length,reconstructed?2243:female?888:2234);assert.equal(atlas.concepts.length,reconstructed?4246:female?1073:3432);
const ids=new Set(atlas.parts.map(p=>p.id));assert.equal(ids.size,reconstructed?2243:female?888:2234);
const files=atlas.chunks.map(c=>{const b=fs.readFileSync(new URL(c.url.split('/').pop(),base));assert.equal(b.length,c.bytes);return b;});
if(female){assert.equal(atlas.sex,'female');assert.equal(atlas.parts.filter(p=>p.system==='pregnancy').length,reconstructed?0:8);for(const label of ['uterus','ovary','vagina'])assert.ok(atlas.parts.some(p=>p.name.toLowerCase().includes(label)));assert.ok(!atlas.parts.some(p=>/prostate|testis|penis/i.test(p.name)));}
if(reconstructed){
 const source=JSON.parse(fs.readFileSync(new URL('atlas.json',base)));
 const femaleSource=JSON.parse(fs.readFileSync(new URL('atlas-female.json',base)));
 const sourceParts=new Map(source.parts.map(p=>[p.id,p]));
 const femaleParts=new Map(femaleSource.parts.map(p=>[p.id,p]));
 const report=JSON.parse(fs.readFileSync(new URL('female-fit-report.json',base)));
 const sourceFiles=source.chunks.map(c=>fs.readFileSync(new URL(c.url.split('/').pop(),base)));
 const femaleFiles=femaleSource.chunks.map(c=>fs.readFileSync(new URL(c.url.split('/').pop(),base)));
 // Re-apply the recorded whole-body morph and breast drape (same definitions as the builder).
 const morph=report.morph,drape=report.drape;
 const smoothstep=(a,b,t)=>{const u=Math.min(1,Math.max(0,(t-a)/(b-a)));return u*u*(3-2*u);};
 const lateral=y=>{const k=morph.lateralKnots;if(y<k[0][0])return k[0][1];for(let i=0;i<k.length-1;i++){const [y0,s0]=k[i],[y1,s1]=k[i+1];if(y>=y0&&y<y1)return s0+(s1-s0)*smoothstep(y0,y1,y);}return k[k.length-1][1];};
 const apply=(x,y,z)=>{
  const R=morph.lateralRadius,t=morph.thoraxDepth,h=morph.head;
  const dx=Math.sign(x)*(lateral(y)-1)*R*Math.tanh(Math.abs(x)/R);
  const w=smoothstep(t.ramp[0],t.ramp[1],y)*(1-smoothstep(t.ramp[2],t.ramp[3],y));
  let qx=x+dx,qy=y,qz=z+z*(t.scale-1)*w;
  const wh=smoothstep(h.ramp[0],h.ramp[1],y)*(1-h.scale);
  qx-=wh*(qx-h.center[0]);qy-=wh*(qy-h.center[1]);qz-=wh*(qz-h.center[2]);
  return [qx*morph.stature,qy*morph.stature,qz*morph.stature];
 };
 const drapeShift=(x,y)=>{
  const fx=Math.min(Math.max((x-drape.origin[0])/drape.cell,0),drape.columns-1.000001),fy=Math.min(Math.max((y-drape.origin[1])/drape.cell,0),drape.rows-1.000001);
  const x0=Math.floor(fx),y0=Math.floor(fy),tx=fx-x0,ty=fy-y0,g=drape.values;
  return g[y0][x0]*(1-tx)*(1-ty)+g[y0][x0+1]*tx*(1-ty)+g[y0+1][x0]*(1-tx)*ty+g[y0+1][x0+1]*tx*ty;
 };
 assert.equal(drape.values.length,drape.rows);assert.ok(drape.values.every(r=>r.length===drape.columns));
 let retained=0,added=0,maxError=0;
 const check=(expected,actual,label)=>{for(let i=0;i<expected.length;i++){const e=Math.abs(expected[i]-actual[i]);maxError=Math.max(maxError,e);assert.ok(e<2e-4,`${label}: vertex mismatch ${e}`);}};
 for(const p of atlas.parts){
  const positions=new Float32Array(files[p.chunk].buffer,files[p.chunk].byteOffset+p.positions,p.vertexCount*3);
  const normals=new Int16Array(files[p.chunk].buffer,files[p.chunk].byteOffset+p.normals,p.vertexCount*3);
  for(let i=0;i<normals.length;i+=3){const l=Math.hypot(normals[i],normals[i+1],normals[i+2])/32767;assert.ok(l>.98&&l<1.02,`${p.id}: unnormalized normal`);}
  if(sourceParts.has(p.id)){
   const ref=sourceParts.get(p.id);
   assert.equal(p.provenance.source,'BodyParts3D 4.0');assert.equal(p.name,ref.name);assert.equal(p.system,ref.system);assert.equal(p.conceptId,ref.conceptId);
   assert.equal(p.vertexCount,ref.vertexCount);assert.equal(p.indexCount,ref.indexCount);
   const original=new Float32Array(sourceFiles[ref.chunk].buffer,sourceFiles[ref.chunk].byteOffset+ref.positions,ref.vertexCount*3);
   const originalIndices=new Uint32Array(sourceFiles[ref.chunk].buffer,sourceFiles[ref.chunk].byteOffset+ref.indices,ref.indexCount);
   const indices=new Uint32Array(files[p.chunk].buffer,files[p.chunk].byteOffset+p.indices,p.indexCount);
   assert.deepEqual(indices,originalIndices,`${p.id}: topology changed`);
   const expected=new Float32Array(original.length);
   for(let i=0;i<original.length;i+=3){const q=apply(original[i],original[i+1],original[i+2]);expected[i]=q[0];expected[i+1]=q[1];expected[i+2]=q[2];}
   check(expected,positions,p.id);retained++;
  } else {
   const reference=femaleParts.get(p.id);assert.ok(reference);
   assert.equal(p.provenance.source,'HRA united-female v1.5');
   const entry=report.added.find(a=>a.id===p.id);assert.ok(entry,`${p.id}: not in fit report`);assert.equal(entry.system,p.system);
   assert.ok(['reproductive','mammary','skeletal'].includes(p.system));
   const transform=report.transforms[entry.transform];assert.ok(transform);
   if(report.regenerated.includes(p.id)){
    // Regenerated breast fat: a closed heightfield body, every edge shared by exactly two triangles.
    const indices=new Uint32Array(files[p.chunk].buffer,files[p.chunk].byteOffset+p.indices,p.indexCount),edges=new Map();
    for(let t=0;t<indices.length;t+=3)for(let k=0;k<3;k++){const a=indices[t+k],b=indices[t+(k+1)%3],key=a<b?`${a}:${b}`:`${b}:${a}`;edges.set(key,(edges.get(key)??0)+1);}
    for(const count of edges.values())assert.equal(count,2,`${p.id}: breast fat body is not a closed surface`);
    assert.ok(p.bounds[0][1]>1.0&&p.bounds[1][1]<1.4&&Math.abs(p.bounds[0][2])<.2,`${p.id}: breast fat outside the chest`);
   } else {
    assert.equal(p.vertexCount,reference.vertexCount);assert.equal(p.indexCount,reference.indexCount);
    const original=new Float32Array(femaleFiles[reference.chunk].buffer,femaleFiles[reference.chunk].byteOffset+reference.positions,reference.vertexCount*3);
    const expected=new Float32Array(original.length);
    for(let i=0;i<original.length;i+=3){
     let x=original[i]*transform.scale[0]+transform.offset[0],y=original[i+1]*transform.scale[1]+transform.offset[1],z=original[i+2]*transform.scale[2]+transform.offset[2];
     if(p.system==='mammary')z+=drapeShift(x,y);
     const q=apply(x,y,z);expected[i]=q[0];expected[i+1]=q[1];expected[i+2]=q[2];
    }
    check(expected,positions,p.id);
   }
   added++;
  }
  for(let axis=0;axis<3;axis++){let lo=Infinity,hi=-Infinity;for(let i=axis;i<positions.length;i+=3){lo=Math.min(lo,positions[i]);hi=Math.max(hi,positions[i]);}assert.ok(Math.abs(lo-p.bounds[0][axis])<1e-6&&Math.abs(hi-p.bounds[1][axis])<1e-6,`${p.id}: stale bounds`);}
 }
 assert.equal(retained,2181);assert.equal(added,62);
 for(const [maleId,replacementIds] of Object.entries(report.replacements)){assert.ok(!ids.has(maleId),'Replaced male pelvis bone still present');for(const id of replacementIds)assert.ok(ids.has(id),`Missing pelvis replacement ${id}`);}
 for(const name of ['Left ilium','Right ilium','Left ischium','Right ischium','Left pubis','Right pubis','Sacrum','Coccyx'])assert.ok(atlas.parts.some(p=>p.name===name&&p.system==='skeletal'),`Missing ${name}`);
 assert.equal(report.retained.length,retained);assert.equal(report.added.length,added);
 for(const id of report.retained)assert.ok(ids.has(id));
 for(const p of report.excluded)assert.ok(!ids.has(p.id));
 assert.ok(!atlas.parts.some(p=>/penis|testicular|testis|prostat|seminal|deferent|epididym|spermatic|perineal|levator ani/i.test(p.name)));
 assert.ok(!ids.has('FJ2810')&&!ids.has('FJ3148'),'Male skin or urethra retained');
 for(const name of ['Left humerus','Right humerus','Left radius','Right radius','Left ulna','Right ulna','Left femur','Right femur','Left tibia','Right tibia','Left fibula','Right fibula']){
  const p=source.parts.find(p=>p.name===name);assert.ok(ids.has(p.id),'Missing original limb bone');
 }
 for(const p of atlas.parts)assert.ok(atlas.concepts.some(c=>c.elements.includes(p.id)),'Unsearchable structure');
 assert.equal(new Set(atlas.concepts.map(c=>c.id)).size,atlas.concepts.length);
 assert.ok(report.checks.minimumTransformDeterminant>0);assert.ok(report.checks.minimumMorphJacobian>0,'Morph folds space');
 assert.ok(report.checks.breastWallMaxResidualM<.002,'Breast tissue floats off the chest wall');
 const l=report.landmarks;assert.ok(l.stature.after<l.stature.before&&l.biacromialWidth.after<l.biacromialWidth.before&&l.biIliacWidth.after>l.biacromialWidth.after*.9,'Female proportions not applied');
 console.log(`Every retained mesh keeps its source topology and follows the recorded female morph (max deviation ${(maxError*1000).toFixed(3)} mm); fitted female meshes match their transforms and drape; male-specific anatomy is excluded.`);
}
let tris=0;
for(const p of atlas.parts){assert.ok(p.name.trim()&&p.name!=='-'&&!p.name.includes('Bounds('));assert.ok(p.conceptId!=='-');assert.ok(p.system);const b=files[p.chunk];assert.ok(p.indices+p.indexCount*4<=b.length);const pos=new Float32Array(b.buffer,b.byteOffset+p.positions,p.vertexCount*3),indices=new Uint32Array(b.buffer,b.byteOffset+p.indices,p.indexCount);assert.ok(indices.length>=3);for(const i of indices)assert.ok(i<p.vertexCount,`${p.id}: invalid vertex`);for(const value of pos)assert.ok(Number.isFinite(value));tris+=p.indexCount/3;}
for(const c of atlas.concepts){assert.ok(c.elements.length);for(const id of c.elements)assert.ok(ids.has(id),`${c.id}: missing ${id}`);}
assert.equal(tris,atlas.triangles);
console.log(`Verified ${ids.size} individually indexed meshes, ${atlas.concepts.length} complete concept mappings, ${tris.toLocaleString()} triangles, and every binary buffer.`);
