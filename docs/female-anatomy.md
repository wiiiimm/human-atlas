# Female study model

Female anatomy is a derived study model, not a scanned reference like the male atlas. It is built from the BodyParts3D male reference model and the HRA female organ set, then reshaped toward estimated female proportions guided by ecorché illustrations. The male atlas remains unchanged. The HRA source atlas (`atlas-female.json`) is kept for rebuilds and validation but is no longer offered in the viewer.

## Current assembly

- 2,181 BodyParts3D meshes retained with their source topology, reshaped by the shared female body morph.
- 8 HRA female pelvis meshes (ilium, ischium, pubis on each side, sacrum, coccyx) fitted to the male hip bone envelope, replacing the male hip bones and sacrum. The bounding-box fit does not establish acetabular congruence, femoral-head clearance or L5–sacrum alignment; those relationships require landmark and attachment measurements (SWR-513).
- 53 source meshes omitted or replaced: male reproductive structures and associated vessels, the male urethra, the male skin/hair envelope, and selected pelvic-floor structures requiring separate redesign.
- 38 HRA female reproductive meshes fitted into the pelvis.
- 16 HRA breast meshes retain their source vertex counts and triangle topology under a shared per-side contour field and translation onto the chest. The six nipple, areola, and areolar-tubercle meshes remain in the optional Body surface layer, hidden by default. The two adipose envelopes use the original HRA base geometry instead of generated concentric-ring shells. Their source base footprint is unscaled; anterior projection is reduced and the lower anterior tissue lifted. These adjustments and placement are artistic estimates, not measurements of a female subject.
- 2 HRA anterolateral knee ligaments fitted to retained joint surfaces, with source topology preserved and Connective tissue controls.
- 2,245 selectable meshes, 4,248 searchable concepts, and 2,437,148 triangles in total.

## Build pipeline

```sh
python3 scripts/build-female-reconstruction.py
```

Requires Python and NumPy. The script reads only the bundled atlases and writes a separate manifest, fit report, and `female-base-*.bin` chunks (about 35 MB compressed; the female model no longer shares chunks with the male model).

1. **Exclusions and replacements.** Selected BodyParts3D structures are dropped from the manifest, including shared pelvic-floor anatomy that still needs a female replacement, and the male pelvis is swapped for the HRA female pelvis. Male concepts that named the pelvis keep working through the replacement mapping. Source identities and reasons are listed in `public/models/female-fit-report.json`.
2. **Affine fits.** All female reproductive meshes receive one positive diagonal affine transform that aligns the HRA bladder bounds with the BodyParts3D bladder bounds, preserving relative positions under that affine transform. This does not establish correct relationships with the separately fitted pelvis or retained muscles. The breast assembly uses the separate source-based contour placement described below; it does not receive the former compressed affine fit or chest drape.
3. **Source-based breast contour.** All 16 breast meshes use the original HRA source topology. For each side, the adipose envelope bounds define a common center and normalized depth/height frame. Source coordinates are translated at scale 1 to the target center `[±0.082, 1.21, 0.088]` m (positive x for left). A smooth front mask rises from normalized depth 0.35 to 0.90; the posterior base remains unchanged relative to that translation. The forward depth is reduced by 0.15 times source depth above the posterior bound, weighted by the front mask. The lower anterior tissue is lifted by at most 0.010 m with smooth height gates at 0, 0.15, 0.45 and 0.70. The same field applies to fat, glands, ducts, supports and optional nipple surfaces on each side, preserving coincident source points. Normals follow the field's inverse-transpose Jacobian. These meshes are placed directly in the final model frame, without a second whole-body morph, drape or tissue-specific inset. `breastContour` records source hashes and all field parameters; `contoured` lists the 16 IDs and `regenerated` is empty. No source triangles are replaced. See [the current breast presentation](breast-presentation.md).
4. **Body morph, localized waist refinement, and scoped glute contour.** A shared whole-body field applies estimated stature, lateral torso/limb proportions, thorax depth, head, and nasal adjustments to retained BodyParts3D meshes and the fitted pelvis/reproductive meshes. The breast meshes use their separate final-frame placement above. The current additional waist field narrows the measured external-oblique band from 245.9 to **229.4 mm (−6.70%)**. It rises and falls smoothly across source heights 1.03–1.12–1.21 m and fades out between |x| = 0.14 and 0.17 m, before reaching the arm bones. It changes lateral coordinates only: every y/z coordinate, 68 selected upper-limb bone meshes including normals, the checked glute meshes, and the measured hip/shoulder envelope widths are unchanged in this waist pass. Of 72 checked arm-muscle meshes, 70 are identical; two vertices in each long triceps head change by less than 0.001 mm. The preceding universal torso/limb field and 2.85 mm nearly rigid inward outer-arm translation remain unchanged. See [the current waist review](female-waist-slimmer-review.md) and its measured scope report.

   Separate posterior and inferior contour components still apply only to the explicitly listed BodyParts3D bilateral gluteus maximus and inferior gluteal vein IDs. Their parameters are unchanged by the waist refinement; bones, pelvic organs, uterine ligaments, lumbar muscles, and unknown IDs receive no such contour. The earlier contour pass broadened the concentrated posterior bump and lowered selected inferior glute vertices by up to 16.77 mm with a tapered return toward the hamstring. These are illustration-guided artistic choices, not subject measurements. Both applicable field variants have a positive sampled minimum Jacobian of 0.3771; normals use the full xyz inverse-transpose Jacobian. These checks do not prove physiological mechanics, correct attachments, or freedom from intersections. The scoped contour exceptions mean all shape changes do not share one universal spatial field.

Every mesh carries source attribution and an adaptation note for the inspector. The body-morph parameters, breast contour, pelvis/reproductive transforms and exclusions are recorded in the fit report. Its legacy `landmarks` field contains **bounding-box proxies**, not anatomically identified landmark measurements:

| Report key | Actual calculation | Female result |
|---|---|---:|
| `stature` | Highest vertex-bound Y coordinate relative to the model origin; assumes the foot baseline is zero | 1.6236 m |
| `biacromialWidth` | Combined lateral bounds of the scapulae, not acromion-to-acromion points | 0.2871 m |
| `biIliacWidth` | Combined lateral bounds of the ilia after fitting; original comparison uses complete hip bones | 0.3087 m |
| `headWidth` | Combined lateral bounds of the parietal bones, not an external head measurement | 0.1363 m |
| `chestDepth` | Anteroposterior extent of the sternum body alone, not thorax depth | 0.0443 m |

These values describe the current assets and are not target measurements or validation tolerances. The report retains its old keys for compatibility; their names must not be interpreted as measured anthropometry.

## Proportion assumptions

All shape adjustments are artistic estimates: the stature factor 0.95, lateral torso height-dependent factors 0.88–1.15, the recorded torso/limb blending field and the compact waist-refinement delta −0.095 (effective central-waist factor 0.860 at source y = 1.12 m), thorax depth factor 0.96, head factor 0.95, nasal projection factor 0.52, pelvis and reproductive bounding-box fits, the source-based breast contour and translation, and the explicitly scoped glute fields (source posterior amplitude 0.020 m and inferior amplitude 0.022 m, further limited by recorded smooth gates). The source data and transforms are traceable, but no declared female subject supplies these target proportions. Stature scaling combines with the head adjustment, so the final extent is not simply 95% of the initial height.

Women vary in body proportions, muscle mass and distribution; no single male-to-female scale factor represents that variation. A large MRI study reports differences associated with sex, age and body size, but it does not validate this model's chosen factors ([Janssen et al., 2000](https://pubmed.ncbi.nlm.nih.gov/10904038/)). SWR-517 now records an illustration-guided target and measured geometric changes; calibrated subject landmarks and attachment validation remain open. Bone shape, muscle volume and soft-tissue distribution cannot be validated by silhouette alone.

## Scope

This is a female study prototype, not a validated female anatomical atlas. The proportions are estimates chosen to read as female at atlas scale, not measurements from a female subject. The musculature is the male reference musculature, reshaped but not redesigned. Female organ placement has not been reviewed anatomically. The model does not yet contain a female external skin envelope, external genital assembly, or a validated female urethra/pelvic-floor replacement; the original male urethra and genital skin are omitted rather than presented as female structures.

The coverage baseline identifies unresolved core-muscle discoverability and explicitly excluded shared pelvic-floor structures. Consult [the current coverage report](anatomy-coverage.md); absence of a separately named mesh is not proof that tissue is geometrically absent. Further work includes core layers, pelvic-floor support, female urethra and organ registration, joint and muscle attachments, tissue presentation and a documented reference build.

The viewer is static: orbit, isolation and exploded layout do not model joint motion or muscle action. Yoga/Pilates mechanics require a separate pose system and reviewed attachment paths and joint relationships (SWR-518). For comparison, [OpenSim's scaling workflow](https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53089158/How%2BScaling%2BWorks) handles anatomical markers, joint frames, muscle attachment points and length-dependent parameters; this application's surface morph does not provide those validations.

## Validation

```sh
npm run check
npm run audit:coverage
npm run test:coverage
npm run test:female-readiness
node scripts/validate-atlas.mjs
node scripts/validate-atlas.mjs atlas-female.json
node scripts/validate-atlas.mjs atlas-female-reconstructed.json
node scripts/validate-interactions.mjs
node --test scripts/hra-breast-contour.test.mjs
node --experimental-strip-types --test scripts/female-category-toggle.test.mjs scripts/chest-visibility.test.mjs
npm run build
```

The reconstruction validator re-applies the recorded body morph and pelvis/reproductive transforms, and the separate HRA breast contour, to source vertices. It checks source topology, unit normals, male-specific exclusions, bounds, search membership, source IDs and sampled Jacobians. All original limb bones remain present. The breast regression checks independently pin the reviewed contour's position/normal/index bytes, verify all 16 breast meshes retain HRA source indices, and compare all 2,227 nonbreast meshes with the main-branch baseline before promotion. An xy-grid geometric screen compares each adipose posterior against the pectoralis anterior on the inner 60% ellipse; it catches a detached or further-intersecting contour and is not attachment validation. Buffer, index, concept and exploded-layout checks cover the canonical and source atlases. These verify software/data integrity, not anatomical accuracy.

All 16 breast structures use one plain reproductive-tissue material, distinct from red skeletal muscle. The source mesh geometry supplies the surface form; the previous procedural lobule/bump texture is not used. The inspector identifies the HRA source base and the estimated contour/placement. The Breast tissue layer can be hidden to reveal pectorals, while nipple and areolar structures remain in the separate optional Body surface layer.

The [previous breast-containment report](breast-containment-audit.md) and its associated candidate, inset and production-comparison reports measure the superseded generated-envelope model. They are preserved as historical evidence and do not validate the current HRA contour. Anatomical placement, containment and chest-wall relationships of the current breast assembly remain unvalidated.
The tissue distinction follows the [US National Cancer Institute SEER anatomy module](https://training.seer.cancer.gov/anatomy/reproductive/female/glands.html), which describes mammary glands, ducts, fat and connective tissue over the pectoral region. That teaching distinction does not validate our tissue shapes or placement; nipple smooth muscle is also distinct from skeletal-muscle mechanics.

## Teaching release readiness

Run `npm run validate:female-readiness` separately from a development build. **The current model is expected to fail this gate.** It requires passing atlas integrity checks, current named-target coverage or independently reviewed teaching-scope exclusions, and revision-bound evidence for anatomical landmarks, attachments, sectional relationships, tissue presentation and representative poses. Final scope needs independent anatomist and movement-educator reviews.

The [review checklist](../data/anatomy/female-review-checklist.json) records unresolved requirements. [Evidence instructions](../data/anatomy/reviews/README.md) describe reviewer roles, artifact hashes and the model fingerprint. Geometry or presentation changes invalidate an older fingerprint. Automated checks verify evidence completeness and freshness, not the truth of a professional judgment. Passing unit tests or building the viewer cannot mark SWR-519 or the female milestone complete; independent review remains outstanding.

The enhancement branch adds the left and right HRA anterolateral knee ligaments to the same `/female` atlas under **Connective tissue**. They use a recorded local bone-constrained fit followed by the existing body morph. Existing body/breast geometry is unchanged. These are estimated static placements; shallow insertion overlap is documented in the joint fitting report. Other knee structures, spinal cord and kidney internals remain offline candidates.
