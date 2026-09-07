import fs from 'node:fs';
import assert from 'node:assert/strict';
import path from 'node:path';
import {fileURLToPath,pathToFileURL} from 'node:url';
import {createHash} from 'node:crypto';
const filename=process.argv[2]??'atlas.json',female=filename.startsWith('atlas-female'),reconstructed=filename==='atlas-female-reconstructed.json';
const option=name=>{const i=process.argv.indexOf(name);return i<0?undefined:process.argv[i+1];};
const defaultModels=fileURLToPath(new URL('../public/models/',import.meta.url));
const base=pathToFileURL(path.resolve(option('--models-dir')??defaultModels)+path.sep),sourceBase=pathToFileURL(path.resolve(option('--source-dir')??defaultModels)+path.sep);
const atlas=JSON.parse(fs.readFileSync(new URL(filename,base)));
const sha=bytes=>createHash('sha256').update(bytes).digest('hex');
assert.equal(atlas.parts.length,reconstructed?2243:female?888:2234);assert.equal(atlas.concepts.length,reconstructed?4246:female?1073:3432);
const ids=new Set(atlas.parts.map(p=>p.id));assert.equal(ids.size,reconstructed?2243:female?888:2234);
const files=atlas.chunks.map(c=>{const b=fs.readFileSync(new URL(c.url.split('/').pop(),base));assert.equal(b.length,c.bytes);return b;});
if(female){assert.equal(atlas.sex,'female');assert.equal(atlas.parts.filter(p=>p.system==='pregnancy').length,reconstructed?0:8);for(const label of ['uterus','ovary','vagina'])assert.ok(atlas.parts.some(p=>p.name.toLowerCase().includes(label)));assert.ok(!atlas.parts.some(p=>/prostate|testis|penis/i.test(p.name)));}
if(reconstructed){
 const source=JSON.parse(fs.readFileSync(new URL('atlas.json',sourceBase)));
 const femaleSource=JSON.parse(fs.readFileSync(new URL('atlas-female.json',sourceBase)));
 const sourceParts=new Map(source.parts.map(p=>[p.id,p]));
 const femaleParts=new Map(femaleSource.parts.map(p=>[p.id,p]));
 const report=JSON.parse(fs.readFileSync(new URL('female-fit-report.json',base)));
 const sourceFiles=source.chunks.map(c=>fs.readFileSync(new URL(c.url.split('/').pop(),sourceBase)));
 const femaleFiles=femaleSource.chunks.map(c=>fs.readFileSync(new URL(c.url.split('/').pop(),sourceBase)));
 // Independently reproduce the recorded non-breast morph and final-frame HRA breast contour.
 const morph=report.morph,breast=report.breastContour;
 const smoothstep=(a,b,t)=>{const u=Math.min(1,Math.max(0,(t-a)/(b-a)));return u*u*(3-2*u);};
 const lateral=(y,k=morph.lateralKnots)=>{if(y<k[0][0])return k[0][1];for(let i=0;i<k.length-1;i++){const [y0,s0]=k[i],[y1,s1]=k[i+1];if(y>=y0&&y<y1)return s0+(s1-s0)*smoothstep(y0,y1,y);}return k[k.length-1][1];};
 const apply=(x,y,z,partId)=>{
  const R=morph.lateralRadius,t=morph.thoraxDepth,h=morph.head;
  let dx=Math.sign(x)*(lateral(y)-1)*R*Math.tanh(Math.abs(x)/R);
  const blend=morph.lateralBlend;
  if(blend){
   const reference=lateral(y,blend.referenceKnots),factor=lateral(y),weight=smoothstep(blend.innerRadius,blend.outerRadius,Math.abs(x)),r=blend.heightRamp;
   const outer=Math.sign(x)*blend.outerTranslation*smoothstep(r[0],r[1],y)*(1-smoothstep(r[2],r[3],y));
   dx=(reference-1)*R*Math.tanh(x/R)+(1-weight)*(factor-reference)*R*Math.tanh(x/R)+weight*outer;
  }
  const waist=morph.waistRefinement;
  if(waist){const r=waist.heightRamp;const weight=smoothstep(r[0],r[1],y)*(1-smoothstep(r[1],r[2],y))*(1-smoothstep(...waist.radialRamp,Math.abs(x)));dx+=waist.delta*R*Math.tanh(x/R)*weight;}
  const w=smoothstep(t.ramp[0],t.ramp[1],y)*(1-smoothstep(t.ramp[2],t.ramp[3],y));
  let qx=x+dx,qy=y,qz=z+z*(t.scale-1)*w;
  const g=morph.gluteProjection;
  if(g?.partIds?.includes(partId)){
   const r2=((Math.abs(x)-g.centerX)/g.radiusX)**2+((y-g.centerY)/g.radiusY)**2;
   let posterior=1-smoothstep(...g.depthRamp,z);
   if(g.lowerDepth){let low=(1-smoothstep(...g.lowerDepth.heightFade,y))*(1-smoothstep(...g.lowerDepth.depthRamp,z));if(g.lowerDepth.heightRise)low*=smoothstep(...g.lowerDepth.heightRise,y);posterior=posterior+low-posterior*low;}
   const gate=(1-smoothstep(0,1,r2))*smoothstep(...g.midlineRamp,Math.abs(x))*posterior;qz-=g.amplitude*gate;
   if(g.inferior){const i=g.inferior,r=i.heightRamp;let weight=smoothstep(...i.midlineRamp,Math.abs(x))*(1-smoothstep(...i.lateralFade,Math.abs(x)))*smoothstep(r[0],r[1],y)*(1-smoothstep(r[2],r[3],y))*(1-smoothstep(...i.depthRamp,z));if(i.posteriorFade)weight*=smoothstep(...i.posteriorFade,z);qy-=i.amplitude*weight;}
  }
  const n=morph.nose;if(n){const r2=((x-n.center[0])/n.radius[0])**2+((y-n.center[1])/n.radius[1])**2,wn=(1-smoothstep(0,1,r2))*smoothstep(n.plane-n.rampDepth,n.plane+n.rampDepth,z);qz-=wn*(1-n.scale)*(z-n.plane);}
  const wh=smoothstep(h.ramp[0],h.ramp[1],y)*(1-h.scale);
  qx-=wh*(qx-h.center[0]);qy-=wh*(qy-h.center[1]);qz-=wh*(qz-h.center[2]);
  return [qx*morph.stature,qy*morph.stature,qz*morph.stature];
 };
 assert.ok(breast,'Missing source-derived breast contour report');
 assert.deepEqual(report.regenerated,[],'The canonical breast assembly must retain HRA topology');
 assert.equal(breast.bodyMorphApplied,false,'Breast geometry is already in the final atlas frame');
 assert.equal(breast.schemaVersion,1);
 assert.equal(breast.sourceManifestSha256,sha(fs.readFileSync(new URL('atlas-female.json',sourceBase))));
 const breastIds=['L','R'].flatMap(side=>['fat','areolar_tubercles','nipple','areola','mammary_lobes','main_lactiferous_ducts','main_lactiferous_sinuses','suspensory_ligaments'].map(name=>`VH_F_${name}_${side}`));
 assert.deepEqual([...breast.ids].sort(),[...breastIds].sort());
 assert.deepEqual([...report.contoured].sort(),[...breastIds].sort());
 const breastChunks=new Set(breastIds.map(id=>femaleParts.get(id).chunk));
 assert.equal(Object.keys(breast.sourceChunks).length,breastChunks.size);
 for(const i of breastChunks)assert.equal(breast.sourceChunks[femaleSource.chunks[i].url],sha(femaleFiles[i]),'Stale HRA breast source chunk hash');
 assert.deepEqual(breast.frontDepthRamp,[.35,.9]);assert.deepEqual(breast.lowerHeightRamp,[0,.15,.45,.7]);assert.equal(breast.normalDifferenceStepM,1e-5);
 for(const side of ['L','R']){
  const t=breast.transforms[side],ref=femaleParts.get(`VH_F_fat_${side}`);
  assert.deepEqual(t.sourceBounds,ref.bounds);assert.deepEqual(t.center,ref.bounds[0].map((v,i)=>(v+ref.bounds[1][i])/2));
  assert.equal(t.scale,1);assert.deepEqual(t.target,[side==='L'?.082:-.082,1.21,.088]);assert.equal(t.projectionReduction,.15);assert.equal(t.lowerLiftM,.01);
 }
 const contour=(point,t)=>{
  const [x,y,z]=point,[lo,hi]=t.sourceBounds;
  const depth=(z-lo[2])/(hi[2]-lo[2]),height=(y-lo[1])/(hi[1]-lo[1]);
  const front=smoothstep(.35,.9,depth),lower=smoothstep(0,.15,height)*(1-smoothstep(.45,.7,height));
  return [x,y+t.lowerLiftM*front*lower,z-t.projectionReduction*(z-lo[2])*front];
 };
 let minimumBreastJacobian=Infinity,maxBreastError=0,maxBreastNormalError=0,contoured=0;
 const contourNormal=(point,normal,t)=>{
  const h=1e-5,yp=[...point],ym=[...point],zp=[...point],zm=[...point];yp[1]+=h;ym[1]-=h;zp[2]+=h;zm[2]-=h;
  const a=(contour(yp,t)[1]-contour(ym,t)[1])/(2*h),b=(contour(zp,t)[1]-contour(zm,t)[1])/(2*h),c=(contour(zp,t)[2]-contour(zm,t)[2])/(2*h);
  assert.ok(Number.isFinite(a*c)&&a*c>0,'Invalid contour Jacobian at a breast source vertex');minimumBreastJacobian=Math.min(minimumBreastJacobian,a*c);
  const n=[normal[0],normal[1]/a,(normal[2]-b*normal[1]/a)/c],length=Math.hypot(...n);assert.ok(length>0&&Number.isFinite(length));
  return n.map(v=>Math.round(v/length*32767));
 };
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
   for(let i=0;i<original.length;i+=3){const q=apply(original[i],original[i+1],original[i+2],p.id);expected[i]=q[0];expected[i+1]=q[1];expected[i+2]=q[2];}
   check(expected,positions,p.id);retained++;
  } else {
   const reference=femaleParts.get(p.id);assert.ok(reference);
   assert.equal(p.provenance.source,'HRA united-female v1.5');
   const entry=report.added.find(a=>a.id===p.id);assert.ok(entry,`${p.id}: not in fit report`);assert.equal(entry.system,p.system);
   assert.ok(['reproductive','mammary','skeletal','integumentary'].includes(p.system));
   assert.equal(p.vertexCount,reference.vertexCount);assert.equal(p.indexCount,reference.indexCount);
   const original=new Float32Array(femaleFiles[reference.chunk].buffer,femaleFiles[reference.chunk].byteOffset+reference.positions,reference.vertexCount*3);
   const originalIndices=new Uint32Array(femaleFiles[reference.chunk].buffer,femaleFiles[reference.chunk].byteOffset+reference.indices,reference.indexCount);
   const indices=new Uint32Array(files[p.chunk].buffer,files[p.chunk].byteOffset+p.indices,p.indexCount);
   assert.deepEqual(indices,originalIndices,`${p.id}: HRA source topology changed`);
   if(breastIds.includes(p.id)){
    assert.equal(entry.transform,'breastContour');
    const t=breast.transforms[p.id.slice(-1)];
    const originalNormals=new Int16Array(femaleFiles[reference.chunk].buffer,femaleFiles[reference.chunk].byteOffset+reference.normals,reference.vertexCount*3);
    for(let i=0;i<original.length;i+=3){
     const point=[original[i],original[i+1],original[i+2]],changed=contour(point,t);
     const normal=contourNormal(point,[originalNormals[i]/32767,originalNormals[i+1]/32767,originalNormals[i+2]/32767],t);
     for(let axis=0;axis<3;axis++){
      const expected=Math.fround((changed[axis]-t.center[axis])*t.scale+t.target[axis]),error=Math.abs(expected-positions[i+axis]);
      maxBreastError=Math.max(maxBreastError,error);assert.ok(error<=1.5e-7,`${p.id}: final-frame breast contour mismatch ${error}`);
      const normalError=Math.abs(normal[axis]-normals[i+axis]);maxBreastNormalError=Math.max(maxBreastNormalError,normalError);
      assert.ok(normalError<=1,`${p.id}: inverse-transpose breast normal mismatch ${normalError}`);
     }
    }
    contoured++;
   } else {
    const transform=report.transforms[entry.transform];assert.ok(transform,`${p.id}: missing affine placement`);
    const expected=new Float32Array(original.length);
    for(let i=0;i<original.length;i+=3){
     const x=original[i]*transform.scale[0]+transform.offset[0],y=original[i+1]*transform.scale[1]+transform.offset[1],z=original[i+2]*transform.scale[2]+transform.offset[2];
     const q=apply(x,y,z,p.id);expected[i]=q[0];expected[i+1]=q[1];expected[i+2]=q[2];
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
 assert.equal(contoured,16);assert.ok(minimumBreastJacobian>0);
 assert.ok(Math.abs(minimumBreastJacobian-breast.minimumVertexJacobian)<1e-8,'Stale breast Jacobian screening result');
 assert.ok(Math.abs(minimumBreastJacobian-report.checks.minimumBreastVertexJacobian)<1e-8);
 assert.equal(report.checks.breastSourceTopologyPreserved,true);
 for(const key of ['drape','breastProfile','lobules','tissueInsets'])assert.ok(!(key in report),`Obsolete breast report field ${key}`);
 for(const key of Object.keys(report.checks))assert.ok(!key.startsWith('breastWall'),'A contour reproduction check is not a chest-wall attachment measurement');
 for(const p of atlas.parts.filter(p=>/VH_F_(nipple|areola)/.test(p.id)))assert.equal(p.system,'integumentary',`${p.id}: nipple and areola stay in the optional Body surface layer`);
 assert.equal(atlas.parts.filter(p=>p.system==='integumentary').length,6,'Expected six optional breast surface structures');
 const l=report.landmarks;assert.ok(l.stature.after<l.stature.before&&l.biacromialWidth.after<l.biacromialWidth.before&&l.biIliacWidth.after>l.biacromialWidth.after*.9,'Female proportions not applied');
 console.log(`Every retained mesh keeps its source topology and follows the recorded female morph (max deviation ${(maxError*1000).toFixed(3)} mm); fitted non-breast meshes match their transforms; ${contoured} HRA breast meshes preserve topology and reproduce the final-frame contour (max position error ${(maxBreastError*1000).toFixed(6)} mm, normal error ${maxBreastNormalError} signed-short units). These integrity checks do not validate breast placement or attachments.`);
}
let tris=0;
for(const p of atlas.parts){assert.ok(p.name.trim()&&p.name!=='-'&&!p.name.includes('Bounds('));assert.ok(p.conceptId!=='-');assert.ok(p.system);const b=files[p.chunk];assert.ok(p.indices+p.indexCount*4<=b.length);const pos=new Float32Array(b.buffer,b.byteOffset+p.positions,p.vertexCount*3),indices=new Uint32Array(b.buffer,b.byteOffset+p.indices,p.indexCount);assert.ok(indices.length>=3);for(const i of indices)assert.ok(i<p.vertexCount,`${p.id}: invalid vertex`);for(const value of pos)assert.ok(Number.isFinite(value));tris+=p.indexCount/3;}
for(const c of atlas.concepts){assert.ok(c.elements.length);for(const id of c.elements)assert.ok(ids.has(id),`${c.id}: missing ${id}`);}
assert.equal(tris,atlas.triangles);
console.log(`Verified ${ids.size} individually indexed meshes, ${atlas.concepts.length} complete concept mappings, ${tris.toLocaleString()} triangles, and every binary buffer.`);
