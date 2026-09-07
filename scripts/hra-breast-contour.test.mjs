// Run: node --test scripts/hra-breast-contour.test.mjs
// Tests the single canonical female manifest; no comparison model is required.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {createHash} from 'node:crypto';
import {gunzipSync} from 'node:zlib';
import {measureAdiposePectoralisClearance, measureAtlasBreastChestClearance} from './breast-chest-clearance.mjs';
const root = new URL('../public/models/', import.meta.url);
const read = name => fs.readFileSync(new URL(name, root));
const json = name => JSON.parse(read(name));
const hash = value => createHash('sha256').update(value).digest('hex');
const atlas = json('atlas-female-reconstructed.json');
const source = json('atlas-female.json');
const fit = json('female-fit-report.json');
const fixture = JSON.parse(fs.readFileSync(new URL('../data/anatomy/female-category-expectations.json', import.meta.url)));
const ids = [...fixture.expected.mammary, ...fixture.expected.integumentary];
const breastIds = new Set(ids);
const parts = new Map(atlas.parts.map(p => [p.id, p]));
const sourceParts = new Map(source.parts.map(p => [p.id, p]));
const buffers = new Map();
const chunk = (a, p) => {
 const url = a.chunks[p.chunk].url;
 if (!buffers.has(url)) buffers.set(url, read(url.split('/').pop()));
 return buffers.get(url);
};
const meshBytes = (a, p) => {
 const b = chunk(a, p);
 return [[p.positions, p.vertexCount*12], [p.normals, p.vertexCount*6], [p.indices, p.indexCount*4]].map(([offset, length]) => b.subarray(offset, offset+length));
};
const smooth = (lo, hi, value) => {const t = Math.max(0, Math.min(1, (value-lo)/(hi-lo))); return t*t*(3-2*t);};

test('canonical report identifies all 16 contoured HRA meshes without regenerated or trial assets', () => {
 assert.deepEqual(fit.regenerated, []);
 assert.deepEqual([...fit.contoured].sort(), [...ids].sort());
 const c = fit.breastContour;
 assert.equal(c.schemaVersion, 1);
 assert.equal(c.source, 'HRA united-female v1.5');
 assert.equal(c.sourceManifestSha256, hash(read('atlas-female.json')));
 assert.deepEqual([...c.ids].sort(), [...ids].sort());
 assert.deepEqual(c.frontDepthRamp, [.35, .9]);
 assert.deepEqual(c.lowerHeightRamp, [0, .15, .45, .7]);
 assert.equal(c.bodyMorphApplied, false);
 assert.ok(Number.isFinite(c.minimumVertexJacobian) && c.minimumVertexJacobian > 0);
 const urls = new Set(ids.map(id => source.chunks[sourceParts.get(id).chunk].url));
 for (const url of urls) assert.equal(c.sourceChunks[url], hash(read(url.split('/').pop())), url);
 assert.ok(!atlas.breastPresentation && !atlas.reconstructionReport?.includes('comparison'));
 assert.ok(atlas.chunks.every(c => !/comparison|hra-breast/.test(c.url)));
 for (const id of ids) {
  const p = parts.get(id), original = sourceParts.get(id);
  assert.equal(p.system, fixture.expected.mammary.includes(id) ? 'mammary' : 'integumentary');
  assert.equal(p.provenance.sourceId, id);
  assert.equal(p.provenance.source, 'HRA united-female v1.5');
  assert.equal(p.vertexCount, original.vertexCount);
  assert.equal(p.indexCount, original.indexCount);
  assert.deepEqual(meshBytes(atlas, p)[2], meshBytes(source, original)[2], `${id}: source triangle indices changed`);
 }
});

test('source base footprint stays fixed while only anterior projection and lower front are contoured', () => {
 let baseVertices = 0, lifted = 0, reduced = 0;
 for (const id of ids) {
  const p = parts.get(id), src = sourceParts.get(id), b = chunk(atlas, p), original = chunk(source, src);
  const t = fit.breastContour.transforms[id.slice(-1)];
  assert.equal(t.scale, 1); assert.equal(t.projectionReduction, .15); assert.equal(t.lowerLiftM, .01);
  assert.deepEqual(t.target, [id.endsWith('_L') ? .082 : -.082, 1.21, .088]);
  assert.deepEqual(t.sourceBounds, sourceParts.get(`VH_F_fat_${id.slice(-1)}`).bounds);
  const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  const sourceLo = [Infinity, Infinity, Infinity], sourceHi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < p.vertexCount; i++) {
   const point = [0, 1, 2].map(axis => original.readFloatLE(src.positions+i*12+axis*4));
   point.forEach((value, axis) => {sourceLo[axis] = Math.min(sourceLo[axis], value); sourceHi[axis] = Math.max(sourceHi[axis], value);});
   const actual = [0, 1, 2].map(axis => b.readFloatLE(p.positions+i*12+axis*4));
   const translated = point.map((value, axis) => Math.fround(value-t.center[axis]+t.target[axis]));
   const depth = (point[2]-t.sourceBounds[0][2])/(t.sourceBounds[1][2]-t.sourceBounds[0][2]);
   const height = (point[1]-t.sourceBounds[0][1])/(t.sourceBounds[1][1]-t.sourceBounds[0][1]);
   const front = smooth(.35, .9, depth), lower = smooth(0, .15, height)*(1-smooth(.45, .7, height));
   const expected = [point[0], point[1]+.01*front*lower, point[2]-.15*(point[2]-t.sourceBounds[0][2])*front].map((value, axis) => Math.fround(value-t.center[axis]+t.target[axis]));
   for (let axis = 0; axis < 3; axis++) {
    assert.equal(actual[axis], expected[axis], `${id}: unexpected contour at vertex ${i}, axis ${axis}`);
    lo[axis] = Math.min(lo[axis], actual[axis]); hi[axis] = Math.max(hi[axis], actual[axis]);
   }
   assert.equal(actual[0], translated[0]);
   assert.ok(actual[1] >= translated[1]-2e-7 && actual[1] <= translated[1]+.010001);
   assert.ok(actual[2] <= translated[2]+2e-7);
   if (depth <= .35) {assert.deepEqual(actual, translated); baseVertices++;}
   if (actual[1] > translated[1]+.001) lifted++;
   if (actual[2] < translated[2]-.001) reduced++;
   const normal = [0, 1, 2].map(axis => b.readInt16LE(p.normals+i*6+axis*2));
   assert.ok(Math.abs(Math.hypot(...normal)-32767) < 2, `${id}: invalid unit normal`);
  }
  assert.deepEqual(p.bounds, [lo, hi]);
  if (fixture.adiposeEnvelopes.includes(id)) for (const bound of [0, 1]) for (const axis of [0, 1]) {
   assert.ok(Math.abs(p.bounds[bound][axis]-([sourceLo, sourceHi][bound][axis]-t.center[axis]+t.target[axis])) < 2e-6, `${id}: base footprint changed`);
  }
 }
 assert.ok(baseVertices > 100 && lifted > 100 && reduced > 100, 'Exercise both fixed and contoured source regions.');
});

test('all breast position, normal and index bytes match the reviewed contour trial', () => {
 // Golden bytes read with git show from saved trial d0b747d; runtime tests do not
 // require that development branch, an extra manifest, or any old binary asset.
 const digest = createHash('sha256');
 for (const id of ids) {digest.update(id+'\n'); for (const bytes of meshBytes(atlas, parts.get(id))) digest.update(bytes);}
 assert.equal(digest.digest('hex'), '530070aae6d938f08babff6a3db5d234911d30e35cd2fc3bc924aea7f1d96ba2');
});

test('all 2227 nonbreast meshes and original concepts match the main baseline exactly', () => {
 // Main baseline dc553d59dc66ad4792f2890b0d645fc671a7012c, independently read
 // through git show before this promotion. Offsets/chunk packing may change;
 // every position, normal, triangle index and semantic part field must not.
 const digest = createHash('sha256'); let count = 0;
 for (const p of atlas.parts) {
  if (breastIds.has(p.id)) continue;
  count++;
  digest.update(JSON.stringify({id:p.id, name:p.name, conceptId:p.conceptId, system:p.system, vertexCount:p.vertexCount, indexCount:p.indexCount, bounds:p.bounds, provenance:p.provenance})+'\n');
  for (const bytes of meshBytes(atlas, p)) digest.update(bytes);
 }
 assert.equal(count, 2227);
 assert.equal(digest.digest('hex'), '212aa73e0094e24e55511fc49dcc6105182aa566d588c90c0da7fa520d68f9ff');
 assert.equal(hash(JSON.stringify(atlas.concepts)), '268333fda59e06ff28644cb9d509da2cd63611a6a500366766d990af690678cd');
});

test('adipose posterior stays in a geometric pectoralis-core band', () => {
 const files = atlas.chunks.map(c => read(c.url.split('/').pop()));
 const measured = measureAtlasBreastChestClearance(atlas, files);
 assert.equal(measured.coreCoveredCells, fit.checks.breastChestCoreCoveredCells);
 assert.equal(measured.coreGapMedianM, fit.checks.breastChestCoreGapMedianM);
 assert.ok(measured.sides.L.coveredCells >= 350 && measured.sides.R.coveredCells >= 350);
 assert.ok(measured.coreGapMedianM > -0.04 && measured.coreGapMedianM < -0.01);
 assert.ok(measured.coreGapMinM > -0.05 && measured.coreGapMaxM < 0.03);
 const quad = (id, z, x0, x1, y0, y1) => ({
  id, bounds: [[x0, y0, z], [x1, y1, z]],
  positions: new Float32Array([x0, y0, z, x1, y0, z, x1, y1, z, x0, y1, z]),
  indices: new Uint32Array([0, 1, 2, 0, 2, 3]),
 });
 const chest = [quad('FJ1464', .08, .02, .12, 1.16, 1.28), quad('FJ1464M', .08, -.12, -.02, 1.16, 1.28)];
 const seated = measureAdiposePectoralisClearance(chest, [quad('VH_F_fat_L', .055, .02, .12, 1.16, 1.28), quad('VH_F_fat_R', .055, -.12, -.02, 1.16, 1.28)]);
 const detached = measureAdiposePectoralisClearance(chest, [quad('VH_F_fat_L', .14, .02, .12, 1.16, 1.28), quad('VH_F_fat_R', .14, -.12, -.02, 1.16, 1.28)]);
 assert.ok(seated.coreCoveredCells > 100 && seated.coreGapMedianM < -0.01 && seated.coreGapMedianM > -0.04);
 assert.ok(detached.coreGapMedianM > 0.03, 'A forward-shifted envelope must fail the geometric front-gap band');
});

test('canonical compressed chunks reproduce raw payloads and source-based triangle totals', () => {
 for (const c of atlas.chunks) {
  const raw = read(c.url.split('/').pop()), compressed = read(c.gzip.split('/').pop());
  assert.equal(raw.length, c.bytes); assert.equal(compressed.length, c.gzipBytes);
  assert.deepEqual(gunzipSync(compressed), raw);
 }
 assert.equal(atlas.parts.length, 2243);
 assert.equal(atlas.triangles, 2436412);
 assert.equal(atlas.triangles, atlas.parts.reduce((sum, part) => sum+part.indexCount/3, 0));
});
