import type {Atlas, Concept} from './anatomy';

export interface SearchConcept extends Concept {
 searchAliases?: readonly string[];
 searchGroup?: {
  description: string;
  expectedPartIds: readonly string[];
  missingPartIds: string[];
 };
}

interface MuscleGroup {
 id: string;
 name: string;
 aliases: readonly string[];
 description: string;
 members: readonly (readonly [id: string, name: string])[];
}

// Search shortcuts, not new source ontology or muscle-activation groups.
// Membership is pinned to actual BodyParts3D IDs AND names in both bundled atlases.
// Anatomical naming reference: UAMS Department of Neuroscience, lower-limb table:
// https://medicine.uams.edu/neuroscience/education/medical-school-courses/human-structure-module/anatomy-tables/muscle-tables/muscles-of-the-lower-limb/
const GROUPS: readonly MuscleGroup[] = [
 {
  id: 'ATLAS:group:hamstrings', name: 'Hamstrings', aliases: ['hamstring', 'posterior thigh muscles'],
  description: 'Biceps femoris (both heads), semitendinosus and semimembranosus on both sides. The short head of biceps femoris crosses only the knee. This search group does not imply that every selected muscle has the same action.',
  members: [
   ['FJ1395', 'Long head of right biceps femoris'], ['FJ1395M', 'Long head of left biceps femoris'],
   ['FJ1444', 'Short head of right biceps femoris'], ['FJ1444M', 'Short head of left biceps femoris'],
   ['FJ1435', 'Right semimembranosus'], ['FJ1435M', 'Left semimembranosus'],
   ['FJ1436', 'Right semitendinosus'], ['FJ1436M', 'Left semitendinosus'],
  ],
 },
 {
  id: 'ATLAS:group:quadriceps', name: 'Quadriceps', aliases: ['quads', 'quadriceps femoris', 'anterior thigh muscles'],
  description: 'Rectus femoris, vastus intermedius, vastus lateralis and vastus medialis on both sides. This selection groups the represented muscle meshes; it does not measure activation.',
  members: [
   ['FJ1433', 'Right rectus femoris'], ['FJ1433M', 'Left rectus femoris'],
   ['FJ1441', 'Right vastus intermedius'], ['FJ1441M', 'Left vastus intermedius'],
   ['FJ1442', 'Right vastus lateralis'], ['FJ1442M', 'Left vastus lateralis'],
   ['FJ1443', 'Right vastus medialis'], ['FJ1443M', 'Left vastus medialis'],
  ],
 },
 {
  id: 'ATLAS:group:calf', name: 'Calf muscles', aliases: ['calf', 'calves', 'superficial posterior leg muscles'],
  description: 'Both gastrocnemius heads, soleus and plantaris on both sides. This shortcut covers the represented superficial calf muscles, not every muscle in the lower leg.',
  members: [
   ['FJ1394', 'Lateral head of right gastrocnemius'], ['FJ1394M', 'Lateral head of left gastrocnemius'],
   ['FJ1397', 'Medial head of right gastrocnemius'], ['FJ1397M', 'Medial head of left gastrocnemius'],
   ['FJ1437', 'Right soleus'], ['FJ1437M', 'Left soleus'],
   ['FJ1429', 'Right plantaris'], ['FJ1429M', 'Left plantaris'],
  ],
 },
 {
  id: 'ATLAS:group:gluteal', name: 'Gluteal muscles', aliases: ['glutes', 'gluteus muscles'],
  description: 'Gluteus maximus, medius and minimus on both sides. This shortcut does not include every muscle of the gluteal region or imply a common activation pattern.',
  members: [
   ['FJ1418', 'Right gluteus maximus'], ['FJ1418M', 'Left gluteus maximus'],
   ['FJ1419', 'Right gluteus medius'], ['FJ1419M', 'Left gluteus medius'],
   ['FJ1420', 'Right gluteus minimus'], ['FJ1420M', 'Left gluteus minimus'],
  ],
 },
];

/** Preserve source concepts; add selectable groups containing verified present meshes only. */
export function buildSearchConcepts(atlas: Pick<Atlas, 'parts'|'concepts'>): SearchConcept[] {
 const parts = new Map(atlas.parts.map(part => [part.id, part]));
 const ids = new Set(atlas.concepts.map(concept => concept.id));
 const result: SearchConcept[] = [...atlas.concepts];
 for (const group of GROUPS) {
  if (ids.has(group.id)) continue;
  const elements = group.members.filter(([id, name]) => {
   const part = parts.get(id);
   return part?.system === 'muscular' && part.name === name;
  }).map(([id]) => id);
  if (!elements.length) continue;
  const present = new Set(elements);
  const expectedPartIds = group.members.map(([id]) => id);
  const missingPartIds = expectedPartIds.filter(id => !present.has(id));
  result.push({
   id: group.id,
   name: group.name + (missingPartIds.length ? ' (available muscles)' : ''),
   elements,
   searchAliases: group.aliases,
   searchGroup: {
    description: missingPartIds.length
     ? `Available subset: ${elements.length} of ${expectedPartIds.length} expected muscle meshes. ${group.description}`
     : group.description,
    expectedPartIds,
    missingPartIds,
   },
  });
  ids.add(group.id);
 }
 return result;
}

/** Familiar aliases augment the existing name and source-ID search. */
export function matchesAnatomySearch(concept: SearchConcept, query: string): boolean {
 const term = query.trim().toLowerCase();
 return term.length > 0 && [concept.name, concept.id, ...(concept.searchAliases ?? [])]
  .some(value => value.toLowerCase().includes(term));
}
