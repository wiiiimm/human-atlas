# Female pelvis assembly: independent surface audit

[SWR-513](https://linear.app/stealth-company/issue/SWR-513) remains **In Progress**. This work adds repeatable geometric evidence and tests an isolated registration proposal. It does not certify the pelvis, change production geometry, or finish the anatomical correction.

## Reproduce

Requires Python 3 and NumPy, already used by the reconstruction builder. Run from the repository root:

```bash
python3 scripts/pelvis-surface-audit-test.py
python3 scripts/pelvis-surface-audit.py
python3 scripts/pelvis-registration-experiment.py
python3 scripts/pelvis-surface-views.py
```

Outputs are [surface measurements](../data/anatomy/pelvis-surface-audit.json), [registration experiment](../data/anatomy/pelvis-registration-experiment.json), and [assembly views](pelvis-surface-views.svg). The surface audit records manifest, fit-report, and used binary-chunk SHA-256 hashes. Regenerate after changing any mesh or transform; old measurements do not describe new geometry.

## What is measured

The 39 screens cover both femur/hip assemblies, 32 muscle/pelvis relationships, L5 and its disc against the sacrum, both ilium/sacrum assemblies, and the two pubis meshes. Every query vertex is measured against actual target **triangle surfaces** using an AABB hierarchy and exact point-to-triangle distance; this is not a comparison of bounding boxes or nearest vertices.

The control is the original male geometry transformed with the stored female body morph. For each retained femur, muscle, vertebra, and disc, the checker verifies identical triangle indices and reproduces current vertices within 0.0002 mm. This isolates the effect of replacing the bones from the shared whole-body warp, without treating the original assembly as validated ground truth.

A **control proximity patch** contains the query vertices within 3 mm of the control bones. Their indices are recorded, and those same vertices are measured against the replacement female bones. The 3 mm patch threshold and 5 mm review threshold are engineering screening values, not physiological tolerances. A proximity patch is **not** an annotated muscle origin, insertion, articular surface, or femoral-head landmark.

## Findings

These values are the 95th percentile of female distances for each *control proximity patch*, not distances across a whole muscle:

| Screen | Right p95 | Left p95 | Interpretation |
| --- | ---: | ---: | --- |
| Femur / hip envelope | 9.07 mm | 8.89 mm | A tiny overall minimum does not demonstrate congruence across the nearby femur surface. |
| Gluteus maximus / hip and sacrum | 25.10 mm | 25.09 mm | Large displacement of formerly nearby surface regions requires local review. |
| Piriformis / hip and sacrum | 13.21 mm | 12.98 mm | Whole-bone replacement has changed nearby relationships. |
| Semitendinosus / hip | 13.14 mm | 14.06 mm | Review the relevant attachment regions before teaching movement. |
| Long head of biceps femoris / hip | 14.19 mm | 15.52 mm | On the left, even the smallest sampled whole-mesh distance is 5.62 mm. |

For the L5 disc, 462 of the 620 control-patch vertices are now more than 5 mm from the female sacrum; patch p95 is 12.88 mm. Its overall minimum is only 0.018 mm. This directly illustrates why a single closest-point distance can conceal poor assembly correspondence.

The female ilium/sacrum and pubis/pubis minimum sample distances are below 0.3 mm. These minima **do not establish healthy SI joint or symphysis spacing**. The closest locations might be wrong-facing surfaces, segment boundaries, or intersections.

## A small registration proposal was tested and rejected

The separate experiment computes nearest-surface correspondence displacements on the right femur, left femur, and L5-disc control proximity patches. Each region has equal weight so the more densely sampled disc does not dominate. A translation cap of 5 mm prevents an unbounded proposal. No production asset is written.

The candidate translation is **(+0.40, −0.87, +0.42) mm** in model x/y/z. It improves 18 of 36 patch p95 values and worsens 18:

| Screen | Current p95 | Candidate p95 |
| --- | ---: | ---: |
| Right femur | 9.07 mm | 9.02 mm |
| Left femur | 8.89 mm | 9.16 mm |
| L5 disc | 12.88 mm | 12.07 mm |
| Left long-head biceps femoris | 15.52 mm | 15.81 mm |

The geometric correspondence suggestions conflict: femur patches suggest roughly 0.8–1.15 mm superior movement; the disc patch suggests 4.56 mm inferior movement. A single translation cannot follow both suggestions. The experiment uses unannotated nearest points and does not justify rotating, deforming, or reattaching anatomy. **The candidate is rejected for production.** Even uniformly better scores would not establish anatomical accuracy.

## Evidence required before changing registration

No reviewed anatomical surface annotations or attachment coordinates were found in the bundled manifests or fit report. The next correction needs annotations in both source frames, preserving their provenance and triangle/barycentric locations:

1. Bilateral femoral-head articular surface regions and corresponding acetabular surface/rim regions, with an independently reviewed definition of the fitted center and orientation. Record fit residuals and excluded neck/trochanter vertices; do not substitute femur bounding-box corners or whole-hip closest points.
2. Sacral superior and L5/disc opposing surface regions, distinguished from posterior processes. Review their relative orientation and spacing together with both hips.
3. Opposing SI joint and symphysis surface regions, distinguished from neighboring bone segmentation boundaries. Account for the connective structures represented or omitted in the source.
4. Source-supported attachment regions for gluteals, iliacus, adductors, hamstrings, and deep rotators, including the femoral insertions and any tendon/ligament mediation. Whole-muscle proximity to an entire hip bone is insufficient.
5. A reviewed target subject/reference build and correspondence between donors. Use those constraints to compare rigid registration, local adjustments, attachment relocation, and body proportion effects; retain any rejected candidate results.

This follows the distinction in [OpenSim's scaling documentation](https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53089158/How%2BScaling%2BWorks): anatomical marker correspondence, joint frames, muscle attachment positions, and wrapping geometry are explicit parts of scaling. Our meshes currently provide no equivalent validated annotation or kinematic model.

## Visual and numerical limits

![Current pelvis assembly in anterior, lateral, and transverse section views](pelvis-surface-views.svg)

The SVG is generated directly from current indexed triangles. The transverse section uses exact triangle intersections with the numerical plane y = 0.860 m. That plane is not labeled as a clinical/anatomical landmark. Femora are cropped for the assembly view. Colors only distinguish structures.

Distances are unsigned: they cannot diagnose penetration, and intersecting or nested surfaces can still yield small values. Source vertices are not uniform area samples; minimum sample distance is an upper bound on continuous surface separation. This is not a collision test, continuous Hausdorff measurement, cartilage-gap assessment, or physiological validation. The five analytic tests verify triangle interiors, edge/vertex cases, degeneracy, hierarchy agreement with an exhaustive oracle, and the limitation of unsigned samples across an intersection.

The bone sources are BodyParts3D/DBCLS and HRA/Visible Human female. Their source attribution and licenses remain in [ATTRIBUTION](../public/ATTRIBUTION.md). [NIH's HRA female pelvis record](https://3d.nih.gov/entries/3DPX-020984) documents the Visible Human origin; it does not validate this project's cross-source registration.
