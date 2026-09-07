import assert from 'node:assert/strict';
import fs from 'node:fs';
import {createHash} from 'node:crypto';
import {gunzipSync} from 'node:zlib';
export const JOINT_ADDITION_IDS = ['VH_F_anterolateral_ligament_of_knee_L', 'VH_F_anterolateral_ligament_of_knee_R'];
const hash = b => createHash('sha256').update(b).digest('hex');
export function validateJointAdditions(atlas, base, sourceBase) {
 const report = JSON.parse(fs.readFileSync(new URL('female-fit-report.json', base)));
 const joint = report.jointAdditions;
 if (!joint) {
  assert.ok(!atlas.jointAdditions, 'Atlas additions lack a fitting record');
  assert.ok(!atlas.parts.some(p => JOINT_ADDITION_IDS.includes(p.id)), 'Unrecorded joint addition');
  return new Set();
 }
 assert.equal(joint.schemaVersion, 1);
 assert.deepEqual(joint.ids, JOINT_ADDITION_IDS);
 assert.deepEqual(atlas.jointAdditions.ids, JOINT_ADDITION_IDS);
 assert.equal(joint.baseManifestSha256, '4134384790ef6be9526105466ab946f4291afdc0e8237145dd0d64cd5809ffa8');
 assert.equal(joint.baseFitReportSha256, '594580f8c99323092efb49121ae7e0a2ef154b73e9e4ca8f6820db74e2e77c76');
 assert.deepEqual(joint.baseCounts, {parts:2243, concepts:4246, triangles:2436412, chunks:16});
 assert.equal(atlas.chunks.length, 17);
 for (const chunk of atlas.chunks.slice(0,16)) {
  assert.equal(hash(fs.readFileSync(new URL(chunk.url.split('/').pop(), base))), joint.baseChunkHashes[chunk.url]);
 }
 const rawSource=fs.readFileSync(new URL('atlas-female.json',sourceBase));
 assert.equal(hash(rawSource),joint.sourceManifestSha256);
 const source=JSON.parse(rawSource);
 const addedChunk=atlas.chunks[16], raw=fs.readFileSync(new URL(addedChunk.url.split('/').pop(),base));
 const gz=fs.readFileSync(new URL(addedChunk.gzip.split('/').pop(),base));
 assert.equal(joint.chunk.url,addedChunk.url);
 assert.equal(hash(raw),joint.chunk.sha256);assert.equal(hash(gz),joint.chunk.gzipSha256);
 assert.deepEqual(gunzipSync(gz),raw);
 assert.equal(raw.length,addedChunk.bytes);assert.equal(gz.length,addedChunk.gzipBytes);
 assert.equal(joint.parts.length,2);
 for (const id of JOINT_ADDITION_IDS) {
  const part=atlas.parts.find(p=>p.id===id), ref=source.parts.find(p=>p.id===id), record=joint.parts.find(p=>p.id===id);
  assert.ok(part&&ref&&record,id);assert.equal(part.system,'connective');assert.equal(part.chunk,16);
  assert.equal(part.provenance.source,'HRA united-female v1.5');assert.equal(part.provenance.sourceId,id);
  assert.equal(part.name,`${id.endsWith('_L')?'Left':'Right'} ${ref.name}`);
  assert.equal(part.vertexCount,ref.vertexCount);assert.equal(part.indexCount,ref.indexCount);
  const sourceChunk=source.chunks[ref.chunk], sourceBytes=fs.readFileSync(new URL(sourceChunk.url.split('/').pop(),sourceBase));
  assert.equal(hash(sourceBytes),joint.sourceChunkHashes[sourceChunk.url]);
  const arrays={positions:raw.subarray(part.positions,part.positions+part.vertexCount*12),normals:raw.subarray(part.normals,part.normals+part.vertexCount*6),indices:raw.subarray(part.indices,part.indices+part.indexCount*4)};
  for(const [field,bytes] of Object.entries(arrays))assert.equal(hash(bytes),record.finalArrayHashes[field],`${id} ${field}`);
  assert.deepEqual(arrays.indices,sourceBytes.subarray(ref.indices,ref.indices+ref.indexCount*4));
  const lo=[Infinity,Infinity,Infinity],hi=[-Infinity,-Infinity,-Infinity];
  for(let i=0;i<part.vertexCount;i++) {
   const n=[0,1,2].map(a=>arrays.normals.readInt16LE(i*6+a*2));
   assert.ok(Math.abs(Math.hypot(...n)-32767)<2,`${id}: normals`);
   for(let a=0;a<3;a++){const v=arrays.positions.readFloatLE(i*12+a*4);assert.ok(Number.isFinite(v));lo[a]=Math.min(lo[a],v);hi[a]=Math.max(hi[a],v);}
  }
  assert.deepEqual(part.bounds,[lo,hi]);
  assert.deepEqual(atlas.concepts.find(c=>c.id===`PART_${id}`)?.elements,[id]);
 }
 assert.deepEqual(atlas.parts.slice(2243).map(p=>p.id),JOINT_ADDITION_IDS);
 assert.deepEqual(atlas.concepts.slice(4246).map(c=>c.id),JOINT_ADDITION_IDS.map(id=>`PART_${id}`));
 assert.equal(hash(JSON.stringify(atlas.concepts.slice(0,4246))),'268333fda59e06ff28644cb9d509da2cd63611a6a500366766d990af690678cd');
 assert.equal(report.parts,atlas.parts.length);assert.equal(report.concepts,atlas.concepts.length);assert.equal(report.triangles,atlas.triangles);
 return new Set(JOINT_ADDITION_IDS);
}

// Reproduce the source-space joint field independently of the Python builder.
export function reproduceJointPositions(id, original, joint, morph) {
 const side=joint.config.sides[id.slice(-1)], fit=side.fit, field=side.refinement;
 assert.equal(field.ridge,.08);assert.equal(field.kernelRadiusM,.025);
 assert.equal(field.controlSourcePointsM.length,field.coefficients.length);
 const expected=new Float32Array(original.length);
 for(let i=0;i<original.length;i+=3) {
  const p=[original[i],original[i+1],original[i+2]];
  const q=[0,1,2].map(a=>fit.translationM[a]+fit.scale*p.reduce((sum,v,b)=>sum+v*fit.rotationRowVectors[b][a],0));
  field.controlSourcePointsM.forEach((center,j)=>{
   const d=p.reduce((sum,v,a)=>sum+(v-center[a])**2,0), w=Math.exp(-d/(2*field.kernelRadiusM**2));
   for(let a=0;a<3;a++)q[a]+=w*field.coefficients[j][a];
  });
  expected.set(morph(...q,id),i);
 }
 return expected;
}
