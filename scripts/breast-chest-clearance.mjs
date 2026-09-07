// Geometric xy-grid screen of adipose posterior z versus pectoralis anterior z.
// Positive gap is a front clearance; negative gap is intersection. Not attachment proof.
import fs from 'node:fs';
import path from 'node:path';
import {fileURLToPath} from 'node:url';

export const CELL = .004;
export const PAD = .01;
export const CORE_R2 = .36;

export const meshBuffers = (files, p) => ({
 positions: new Float32Array(files[p.chunk].buffer, files[p.chunk].byteOffset + p.positions, p.vertexCount * 3),
 indices: new Uint32Array(files[p.chunk].buffer, files[p.chunk].byteOffset + p.indices, p.indexCount),
 bounds: p.bounds,
});

export const pectoralisParts = parts => parts.filter(p => p.system === 'muscular' && /pectoralis/i.test(p.name));

const fmax = (a, b) => Number.isNaN(a) ? b : Number.isNaN(b) ? a : Math.max(a, b);
const fmin = (a, b) => Number.isNaN(a) ? b : Number.isNaN(b) ? a : Math.min(a, b);
const median = values => {
 const sorted = Float64Array.from(values).sort();
 const n = sorted.length, mid = n >> 1;
 return n % 2 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
};

export function projectZ(meshes, gx0, gx1, gy0, gy1, cell, reduce) {
 const nx = Math.round((gx1 - gx0) / cell) + 1, ny = Math.round((gy1 - gy0) / cell) + 1;
 const grid = new Float64Array(ny * nx).fill(NaN);
 const pick = reduce === 'max' ? fmax : fmin;
 const at = (i, j, z) => { const k = j * nx + i; grid[k] = pick(grid[k], z); };
 for (const {positions, indices} of meshes) {
  for (let t = 0; t < indices.length; t += 3) {
   const ia = indices[t] * 3, ib = indices[t + 1] * 3, ic = indices[t + 2] * 3;
   const ax = positions[ia], ay = positions[ia + 1], az = positions[ia + 2];
   const bx = positions[ib], by = positions[ib + 1], bz = positions[ib + 2];
   const cx = positions[ic], cy = positions[ic + 1], cz = positions[ic + 2];
   if (Math.max(ay, by, cy) < gy0 || Math.min(ay, by, cy) > gy1 || Math.max(ax, bx, cx) < gx0 || Math.min(ax, bx, cx) > gx1) continue;
   const i0 = Math.max(0, Math.ceil((Math.min(ax, bx, cx) - gx0) / cell));
   const i1 = Math.min(nx - 1, Math.floor((Math.max(ax, bx, cx) - gx0) / cell));
   const j0 = Math.max(0, Math.ceil((Math.min(ay, by, cy) - gy0) / cell));
   const j1 = Math.min(ny - 1, Math.floor((Math.max(ay, by, cy) - gy0) / cell));
   if (i1 < i0 || j1 < j0) continue;
   const det = (bx - ax) * (cy - ay) - (cx - ax) * (by - ay);
   if (Math.abs(det) < 1e-12) continue;
   for (let j = j0; j <= j1; j++) for (let i = i0; i <= i1; i++) {
    const x = gx0 + i * cell, y = gy0 + j * cell;
    const u = ((x - ax) * (cy - ay) - (cx - ax) * (y - ay)) / det;
    const v = ((bx - ax) * (y - ay) - (x - ax) * (by - ay)) / det;
    if (u < -1e-9 || v < -1e-9 || u + v > 1 + 1e-9) continue;
    at(i, j, az + u * (bz - az) + v * (cz - az));
   }
  }
  for (let v = 0; v < positions.length; v += 3) {
   const x = positions[v], y = positions[v + 1];
   if (x < gx0 || x > gx1 || y < gy0 || y > gy1) continue;
   at(Math.max(0, Math.min(nx - 1, Math.round((x - gx0) / cell))), Math.max(0, Math.min(ny - 1, Math.round((y - gy0) / cell))), positions[v + 2]);
  }
 }
 return {grid, nx, ny};
}

export function measureAdiposePectoralisClearance(pectoralis, adipose) {
 const lo = [Infinity, Infinity], hi = [-Infinity, -Infinity];
 for (const mesh of adipose) for (const axis of [0, 1]) {
  lo[axis] = Math.min(lo[axis], mesh.bounds[0][axis]);
  hi[axis] = Math.max(hi[axis], mesh.bounds[1][axis]);
 }
 const gx0 = lo[0] - PAD, gx1 = hi[0] + PAD, gy0 = lo[1] - PAD, gy1 = hi[1] + PAD;
 const wall = projectZ(pectoralis, gx0, gx1, gy0, gy1, CELL, 'max');
 const back = projectZ(adipose, gx0, gx1, gy0, gy1, CELL, 'min');
 const coreGaps = [], allGaps = [];
 const sides = {};
 for (const mesh of adipose) {
  const side = mesh.id.endsWith('_L') ? 'L' : 'R';
  const [b0, b1] = mesh.bounds, cx = (b0[0] + b1[0]) / 2, cy = (b0[1] + b1[1]) / 2, rx = (b1[0] - b0[0]) / 2, ry = (b1[1] - b0[1]) / 2;
  const gaps = [];
  for (let j = 0; j < wall.ny; j++) for (let i = 0; i < wall.nx; i++) {
   const w = wall.grid[j * wall.nx + i], b = back.grid[j * wall.nx + i];
   if (Number.isNaN(w) || Number.isNaN(b)) continue;
   const gap = b - w, x = gx0 + i * CELL, y = gy0 + j * CELL;
   const r2 = ((x - cx) / rx) ** 2 + ((y - cy) / ry) ** 2;
   if (r2 > 1) continue;
   allGaps.push(gap);
   if (r2 <= CORE_R2) {gaps.push(gap); coreGaps.push(gap);}
  }
  sides[side] = {coveredCells: gaps.length, intersectingCells: gaps.filter(g => g < 0).length, gapMinM: Math.min(...gaps), gapMedianM: median(gaps), gapMaxM: Math.max(...gaps)};
 }
 return {
  cellM: CELL, padM: PAD, coreR2: CORE_R2, grid: {origin: [gx0, gy0], extent: [gx1, gy1]},
  coveredCells: allGaps.length, intersectingCells: allGaps.filter(g => g < 0).length,
  gapMinM: Math.min(...allGaps), gapMedianM: median(allGaps), gapMaxM: Math.max(...allGaps),
  coreCoveredCells: coreGaps.length, coreIntersectingCells: coreGaps.filter(g => g < 0).length,
  coreGapMinM: Math.min(...coreGaps), coreGapMedianM: median(coreGaps), coreGapMaxM: Math.max(...coreGaps),
  sides,
 };
}

export function measureAtlasBreastChestClearance(atlas, files) {
 const pectoralis = pectoralisParts(atlas.parts).map(p => ({id: p.id, ...meshBuffers(files, p)}));
 const adipose = ['VH_F_fat_L', 'VH_F_fat_R'].map(id => {const p = atlas.parts.find(part => part.id === id); return {id, ...meshBuffers(files, p)};});
 return measureAdiposePectoralisClearance(pectoralis, adipose);
}

if (process.argv[1] && fileURLToPath(import.meta.url) === path.resolve(process.argv[1])) {
 const dir = process.argv[2];
 const atlas = JSON.parse(fs.readFileSync(path.join(dir, 'atlas-female-reconstructed.json')));
 const files = atlas.chunks.map(c => fs.readFileSync(path.join(dir, c.url.split('/').pop())));
 process.stdout.write(JSON.stringify(measureAtlasBreastChestClearance(atlas, files)));
}
