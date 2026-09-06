# Female study model

Female anatomy is built from the BodyParts3D male reference model and the HRA female organ set, then reshaped toward estimated female proportions. The male atlas remains unchanged. Female · source only still shows the original HRA assembly for comparison.

## Current assembly

- 2,181 BodyParts3D meshes retained with their source topology, reshaped by the shared female body morph.
- 8 HRA female pelvis meshes (ilium, ischium, pubis on each side, sacrum, coccyx) fitted to the male hip bone envelope, replacing the male hip bones and sacrum. The acetabula meet the retained femoral heads within 2 mm and the sacrum meets the fifth lumbar vertebra.
- 53 source meshes omitted or replaced: male reproductive structures and associated vessels, the male urethra, the male skin/hair envelope, and selected pelvic-floor structures requiring separate redesign.
- 38 HRA female reproductive meshes fitted into the pelvis.
- 16 HRA breast meshes fitted over the chest and draped onto the pectoral wall, in a separate Breast tissue layer. The two fat bodies are regenerated as feathered heightfields (HRA thickness profile, edges tapering into the chest wall) so the breast rises from the chest instead of sitting on it.
- 2,243 selectable meshes, 4,246 searchable concepts, and 2,352,682 triangles in total.

## Build pipeline

```sh
python3 scripts/build-female-reconstruction.py
```

Requires Python and NumPy. The script reads only the bundled atlases and writes a separate manifest, fit report, and `female-base-*.bin` chunks (about 35 MB compressed; the female model no longer shares chunks with the male model).

1. **Exclusions and replacements.** Male-specific BodyParts3D structures are dropped from the manifest, and the male pelvis is swapped for the HRA female pelvis. Male concepts that named the pelvis keep working through the replacement mapping. Source identities and reasons are listed in `public/models/female-fit-report.json`.
2. **Affine fits.** All female reproductive meshes receive one positive diagonal affine transform that aligns the HRA bladder bounds with the BodyParts3D bladder bounds, preserving the assembly's internal alignment. The breast assembly keeps its HRA size and is placed so the nipple line sits at the fourth intercostal level of the chest.
3. **Breast drape.** A depth map of the anterior chest wall (front-most surface of chest muscles and bones on a 4 mm grid) is compared with the posterior envelope of the breast assembly. Each breast column is shifted in depth so the assembly rests on the pectoral wall instead of floating in front of it, with a 2 mm standoff so the muscle never pokes through the thin upper breast. The shift grid is stored in the fit report.
4. **Whole-body morph.** One smooth displacement field is applied to every mesh, retained and fitted alike: stature scaled to about 1.62 m, shoulders narrowed, thorax and waist narrowed, pelvis widened at the iliac crests, and the skull scaled down about its centre. Lateral changes saturate beyond the torso's half-width so the arms translate with the shoulders and hips rather than being squeezed. The field is continuous with a positive Jacobian everywhere (sampled minimum recorded in the report), so bones, muscles, vessels, and organs deform together and stay attached. Normals use the inverse-transpose of the local Jacobian.

Every mesh carries source attribution and an adaptation note for the inspector. The morph parameters, drape grid, transforms, exclusions, and a before/after landmark table (stature, biacromial width, bi-iliac width, head width, chest depth) are recorded in the fit report.

## Scope

This is a female study prototype, not a validated female anatomical atlas. The proportions are estimates chosen to read as female at atlas scale, not measurements from a female subject. The musculature is the male reference musculature, reshaped but not redesigned. Female organ placement has not been reviewed anatomically. The model does not yet contain a female external skin envelope, external genital assembly, or a validated female urethra/pelvic-floor replacement; the original male urethra and genital skin are omitted rather than presented as female structures.

Further work should focus on a female body surface, local pelvic morphology (subpubic angle, pelvic inlet shape), and soft-tissue distribution. Those changes should be independently reviewed before claiming anatomical accuracy.

## Validation

```sh
npm run check
node scripts/validate-atlas.mjs
node scripts/validate-atlas.mjs atlas-female.json
node scripts/validate-atlas.mjs atlas-female-reconstructed.json
node scripts/validate-interactions.mjs
npm run build
```

The reconstruction validator re-applies the recorded morph, transforms, and drape grid to every source vertex and compares the result with the bundled geometry, checks that retained meshes keep their source topology and normals stay unit length, verifies male-specific exclusions, bounds, search membership, source IDs, a positive morph Jacobian, and that the breast tissue rests on the chest wall. All original limb bones must remain present. Buffer, index, concept, and exploded-layout checks cover all three atlases. These verify software/data integrity, not anatomical accuracy.
