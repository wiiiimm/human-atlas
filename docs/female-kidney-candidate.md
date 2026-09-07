# Female kidney internal candidate

The registration experiments, including the new constrained surface fit, remain disabled. No viewer mesh is added or replaced. The script exports eight temporary assemblies and measured fit evidence; none establishes a coherent fit to the retained kidneys and ureters.

```sh
PYTHONDONTWRITEBYTECODE=1 python3 scripts/female_kidney_candidate.py --output-dir /tmp/female-kidney-surface-refinement
```

NumPy is the only third-party dependency. The output directory must resolve outside this repository, including through symlinks. `report.json` records `enabled: false`, input hashes, the complete inventory, transforms, measurements, and artifact hashes. Each `{L,R}-{translation,similarity_icp,affine_bounds_diagnostic,surface_affine_constrained}.npz` contains transformed internals and capsule plus unchanged target kidney and ureter context. These are review artifacts, not a viewer manifest.

## Source and inventory

The source is Kristen Browne and Heidi Schlehlein, Human Reference Atlas / HuBMAP, [3D Reference Organ Set for Female v1.5 (2023)](https://doi.org/10.48539/HBM352.BTSQ.586), CC BY 4.0. The bundled atlas and chunk are pinned before geometry is read:

| Input | SHA-256 |
| --- | --- |
| `public/models/atlas-female.json` | `1525c07d2ed46263c086d1c6b2e52eb9b8f8985746e9a8259036bd6e1d9a1ced` |
| `public/models/female-6.bin` | `8b9d1c038bf4b7951bc19a0e30dd3c7c65c29af3f2bd9f605b8048204c66ee40` |
| Exact inventory and geometry digest | `e3715d7db0bf5b713afc688593842a35d035da09ca9b68a2cdc98ee66b503cf4` |

The exact 75 internal meshes comprise 39 left and 36 right parts. Each side has an outer cortex, renal column, and renal pelvis. Left/right counts are 11/10 papillae, 11/10 pyramids, 10/10 minor calyces, and 4/3 major calyces. Capsules, hila, and source ureters supply registration context and are excluded from the 75. Explicit IDs, expected names, urinary category, index ranges, finite positions, and source geometry hashes are checked. The target kidney IDs are `FJ3145`/`FJ3147`; retained ureters are `FJ3144`/`FJ3146`.

## Baseline measurements

These results use canonical manifest `4134384790ef6be9526105466ab946f4291afdc0e8237145dd0d64cd5809ffa8`. Surface p95 is the greater of the two directional 95th-percentile distances, evaluated at every actual capsule/target vertex against nearest triangles. It is a vertex-sampled surface screen, not a continuous Hausdorff distance or anatomical acceptance threshold.

| Side | Whole-assembly transform | Surface p95, mm | Source volume multiplier | Ureter patch p95, mm |
| --- | --- | ---: | ---: | ---: |
| Left | Translation | 19.82 | 1.000 | 7.19 |
| Left | Similarity ICP | 8.49 | 0.320 | 2.60 |
| Left | Bounds affine diagnostic | 6.79 | 0.273 | 1.10 |
| Right | Translation | 18.50 | 1.000 | 4.42 |
| Right | Similarity ICP | 16.88 | 0.361 | 4.45 |
| Right | Bounds affine diagnostic | 10.57 | 0.470 | 4.70 |

Similarity ICP proposes a rotation and uniform scale using deterministic vertex correspondences; triangle distances independently score the proposal. Bounds matching requires left x/y/z scales of 0.708/0.806/0.478 and right scales of 0.738/1.013/0.629. These substantial, unequal compressions still leave measurable surface disagreement. They are recorded as diagnostics, not promoted fits.

Each side uses one shared affine transform for all internals and the reference capsule. Triangle indices and coincident source points remain preserved, normals use the inverse transpose, and the constant Jacobian is positive. Translation and similarity preserve angles; the bounds affine changes internal proportions. Positive Jacobians and preserved topology do not establish anatomical compatibility.

The ureter screen follows the same 80 left and 64 right renal-pelvis vertices originally within 2 mm of source ureter triangles, then measures their transformed distance to the retained ureter. This engineering proximity patch is not a reviewed ureteropelvic-junction landmark. Small patch distances do not establish connection, orientation, or vascular compatibility.

## Containment is unavailable

The source capsules are closed, consistently oriented, single-component tissue meshes. Their signed material volumes are 12.54/9.96 ml, versus bounding-box volumes of 706.95/499.29 ml. Outer cortex meshes are also single connected tissue components, with material volumes of 130.74/83.81 ml. No separate component is justified as a complete organ enclosure. The retained target kidneys have 10/20 boundary edges, so their signed triangle integrals are not valid enclosed volumes.

The report retains three-direction ray-parity diagnostics for reproducibility but explicitly marks them invalid as organ-envelope containment. In particular, an internal point outside capsule tissue is not evidence that the source anatomy is displaced. Unsigned internal sample-to-surface distances are also recorded, without turning them into inside/outside claims. A reviewed enclosure and renal landmark/interface correspondence are needed before further fitting can justify promotion.


## Constrained surface refinement

`scripts/female_kidney_registration.py` extracts a diagnostic capsule exterior patch: an outward-facing triangle is retained when its centroid is the farthest intersection along a ray from the capsule bounds center. The patch consists entirely of existing source triangles. It may omit concavities and the hilum, and it is explicitly not an organ enclosure. Its indices in each NPZ refer to the reference capsule positions. No gaps are filled.

The added candidate uses bidirectional nearest-triangle correspondences, with five deterministic orientation seeds from principal axes and bounds scaling. Principal axes provide geometric long-axis proposals, not annotated anatomical axes. Each iteration solves one affine transform for the entire side, including all its internal parts. The original renal-pelvis proximity patch is constrained toward retained ureter triangles with relative weight 0.25. This discourages an exterior-only fit that rotates the collecting system away from its interface; it cannot establish the missing anatomical correspondence.

Singular scales are constrained to 0.55–1.10 and volume multiplier to at least 0.35. These are engineering limits on this experiment, not physiologically approved tolerances. The fit uses up to 192 vertices per surface direction and 12 iterations per seed. Final comparisons use every capsule and target vertex against actual triangles, matching the baseline metric. The report also includes both directions of exterior-patch distances, transformed hilum-to-ureter distances, every seed score, and iteration residuals.

The analytic regression tests cover rejection of inner capsule material faces, translation invariance of exterior extraction, positive bounded affine determinants, exact shared-transform recovery and normal tangency, and recovery of an exact synthetic similarity through triangle-surface fitting:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 scripts/female-kidney-registration.test.py
```


The refined candidates produce the following measured tradeoffs against bounds-only fitting (left/right):

| Measurement | Bounds-only baseline | Surface refinement |
| --- | --- | --- |
| Full capsule/target surface p95, mm | 6.79 / 10.57 | 5.91 / 7.10 |
| Volume multiplier | 0.273 / 0.470 | 0.350 / 0.434 |
| Ureter proximity patch p95, mm | 1.10 / 4.70 | 1.53 / 1.50 |
| Exterior patch/target surface p95, mm | See report | 6.14 / 7.15 |

The left result reduces compression and surface residual while slightly worsening ureter proximity. The right improves surface residual and ureter proximity with slightly greater volume compression. The singular scales are 0.871/0.707/0.568 on the left and 0.949/0.832/0.550 on the right. Thus both remain substantially adapted assemblies, not approved source-preserving placements. Polar rotations are about 34.6° and 22.8°; the report records geometric long-axis alignment separately.

Exterior patches retain 1,473/1,623 source triangles and have 63/43 boundary edges, confirming that the extraction supplies open fitting patches. Containment remains unavailable. These are practical improvements to an experimental fit, with clear residuals and tradeoffs; the report continues to set `enabled: false`.
