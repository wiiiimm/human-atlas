import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {buildSearchConcepts, matchesAnatomySearch} from '../app/anatomy-search.ts';

const read = name => JSON.parse(fs.readFileSync(new URL(`../public/models/${name}`, import.meta.url)));
const male = read('atlas.json');
const female = read('atlas-female-reconstructed.json');
const expected = {
 'ATLAS:group:hamstrings': ['FJ1395', 'FJ1395M', 'FJ1444', 'FJ1444M', 'FJ1435', 'FJ1435M', 'FJ1436', 'FJ1436M'],
 'ATLAS:group:quadriceps': ['FJ1433', 'FJ1433M', 'FJ1441', 'FJ1441M', 'FJ1442', 'FJ1442M', 'FJ1443', 'FJ1443M'],
 'ATLAS:group:calf': ['FJ1394', 'FJ1394M', 'FJ1397', 'FJ1397M', 'FJ1437', 'FJ1437M', 'FJ1429', 'FJ1429M'],
 'ATLAS:group:gluteal': ['FJ1418', 'FJ1418M', 'FJ1419', 'FJ1419M', 'FJ1420', 'FJ1420M'],
};

test('male and female shortcuts select the exact intended bilateral muscles, preserving source concepts', () => {
 for (const atlas of [male, female]) {
  const original = JSON.stringify(atlas);
  const concepts = buildSearchConcepts(atlas);
  assert.equal(concepts.length, atlas.concepts.length + 4);
  for (let i = 0; i < atlas.concepts.length; i++) assert.equal(concepts[i], atlas.concepts[i]);
  for (const [id, elements] of Object.entries(expected)) {
   const group = concepts.find(concept => concept.id === id);
   assert.deepEqual(group.elements, elements);
   assert.deepEqual(group.searchGroup.missingPartIds, []);
   assert.equal(new Set(group.elements).size, group.elements.length);
   assert.ok(group.elements.every(id => atlas.parts.some(part => part.id === id && part.system === 'muscular')));
  }
  assert.equal(JSON.stringify(atlas), original, 'Search shortcuts must not modify source atlas data.');
  assert.ok(concepts.find(c => c.id === 'FMA22429'), 'Original quadriceps source concept remains searchable.');
 }
});

test('familiar aliases work alongside original anatomical names and source IDs', () => {
 const concepts = buildSearchConcepts(female);
 for (const [query, id] of [['HAMSTRING', 'hamstrings'], [' quads ', 'quadriceps'], ['calves', 'calf'], ['glutes', 'gluteal']]) {
  assert.ok(concepts.filter(concept => matchesAnatomySearch(concept, query)).some(c => c.id === `ATLAS:group:${id}`));
 }
 const source = female.concepts.find(c => c.id === 'FMA22429');
 assert.equal(matchesAnatomySearch(source, source.name), true);
 assert.equal(matchesAnatomySearch(source, 'fma22429'), true);
 assert.equal(matchesAnatomySearch(source, ' '), false);
 assert.equal(matchesAnatomySearch(source, 'not a structure'), false);
});

test('missing or mismatched parts produce an explicit available subset without dangling IDs', () => {
 const partial = {...female, parts: female.parts.filter(p => p.id !== 'FJ1433').map(p => p.id === 'FJ1441' ? {...p, name: 'Unrelated replacement'} : p)};
 const group = buildSearchConcepts(partial).find(c => c.id === 'ATLAS:group:quadriceps');
 assert.equal(group.name, 'Quadriceps (available muscles)');
 assert.match(group.searchGroup.description, /^Available subset: 6 of 8/);
 assert.deepEqual(group.searchGroup.missingPartIds, ['FJ1433', 'FJ1441']);
 assert.ok(!group.elements.includes('FJ1433') && !group.elements.includes('FJ1441'));
 assert.ok(group.elements.every(id => partial.parts.some(p => p.id === id)));
});

test('partial HRA source does not receive unsupported complete groups or empty selections', () => {
 const source = read('atlas-female.json');
 const concepts = buildSearchConcepts(source);
 assert.deepEqual(concepts, source.concepts);
 assert.deepEqual(buildSearchConcepts({parts: [], concepts: []}), []);
});

test('building an augmented catalogue twice preserves unique group IDs and original entries', () => {
 const once = buildSearchConcepts(male);
 const twice = buildSearchConcepts({...male, concepts: once});
 assert.deepEqual(twice, once);
 assert.equal(new Set(twice.map(c => c.id)).size, twice.length);
 const reserved = {id: 'ATLAS:group:quadriceps', name: 'Existing source entry', elements: ['FJ1433']};
 const concepts = buildSearchConcepts({...male, concepts: [...male.concepts, reserved]});
 assert.equal(concepts.find(c => c.id === reserved.id), reserved);
});
