# Female study model

Female anatomy is a derived study model, not a scanned reference like the male atlas. It is built from the BodyParts3D male reference model and the HRA female organ set, then reshaped toward estimated female proportions guided by ecorché illustrations. The male atlas remains unchanged. The HRA source atlas (`atlas-female.json`) is kept for rebuilds and validation but is no longer offered in the viewer.

## Current assembly

- 2,181 BodyParts3D meshes retained with their source topology, reshaped by the shared female body morph.
- 8 HRA female pelvis meshes (ilium, ischium, pubis on each side, sacrum, coccyx) fitted to the male hip bone envelope, replacing the male hip bones and sacrum. The acetabula meet the retained femoral heads within 2 mm and the sacrum meets the fifth lumbar vertebra.
- 53 source meshes omitted or replaced: male reproductive structures and associated vessels, the male urethra, the male skin/hair envelope, and selected pelvic-floor structures requiring separate redesign.
- 38 HRA female reproductive meshes fitted into the pelvis.
- 16 HRA breast meshes fitted over the chest and draped onto the pectoral wall, across the Breast tissue and optional Body surface layers. The six nipple, areola, and areolar-tubercle meshes sit in the optional Body surface layer, hidden by default. The two breast bodies are regenerated as smooth closed concentric-ring surfaces guided by yoga and pilates ecorché illustrations: centred on the sternocostal pectoralis (top at mid-pectoral, bottom at the pectoral's lower border, lateral edge toward the axilla) with a full, rounded projection. They are artistic estimates, not source-measured anatomy.
- 2,243 selectable meshes, 4,246 searchable concepts, and 2,352,682 triangles in total.

## Build pipeline

```sh
python3 scripts/build-female-reconstruction.py
```

Requires Python and NumPy. The script reads only the bundled atlases and writes a separate manifest, fit report, and `female-base-*.bin` chunks (about 35 MB compressed; the female model no longer shares chunks with the male model).

1. **Exclusions and replacements.** Male-specific BodyParts3D structures are dropped from the manifest, and the male pelvis is swapped for the HRA female pelvis. Male concepts that named the pelvis keep working through the replacement mapping. Source identities and reasons are listed in `public/models/female-fit-report.json`.
2. **Affine fits.** All female reproductive meshes receive one positive diagonal affine transform that aligns the HRA bladder bounds with the BodyParts3D bladder bounds, preserving the assembly's internal alignment. The breast assembly is compressed in height and depth to fit the reference-guided contour.
3. **Breast contour.** A chest depth map supports a smooth elliptical mound on each side. The shell back embeds into the chest wall; internal HRA structures shift in depth beneath the new surface, with an extra inset for ducts, lobes and ligaments so they do not protrude. Nipples remain slightly proud. Surface normals are recomputed after fitting. The contour parameters and tissue displacement grid are recorded in the fit report. The stock reference images are visual guides only and are not bundled as assets.
4. **Whole-body morph.** One smooth displacement field is applied to every mesh, retained and fitted alike: stature scaled to about 1.62 m, shoulders narrowed, thorax and waist narrowed, pelvis widened at the iliac crests, the skull scaled down about its centre, and the nasal bridge and projection reduced (depth beyond the face plane is compressed inside a smooth window around the nasal bones and cartilages, for a flatter East Asian profile). Lateral changes saturate beyond the torso's half-width so the arms translate with the shoulders and hips rather than being squeezed. The field is continuous with a positive Jacobian everywhere (sampled minimum recorded in the report), so bones, muscles, vessels, and organs deform together and stay attached. Normals use the inverse-transpose of the local Jacobian.

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

The breast body renders in a muscle-family salmon tone with the same material response as muscle, with fine fibre striations radiating from the apex drawn in the shader, matching the illustration style. It is an illustrative estimate, not measured gland anatomy.
