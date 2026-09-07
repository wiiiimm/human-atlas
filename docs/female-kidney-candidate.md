# Female kidney internal candidate

The first registration experiments remain disabled. No viewer mesh is added or replaced. The script exports six temporary assemblies and measured fit evidence; none establishes a coherent fit to the retained kidneys and ureters.

```sh
PYTHONDONTWRITEBYTECODE=1 python3 scripts/female_kidney_candidate.py --output-dir /tmp/female-kidney-candidate
```

NumPy is the only third-party dependency. The output directory must resolve outside this repository, including through symlinks. `report.json` records `enabled: false`, input hashes, the complete inventory, transforms, measurements, and artifact hashes. Each `{L,R}-{translation,similarity_icp,affine_bounds_diagnostic}.npz` contains transformed internals and capsule plus unchanged target kidney and ureter context. These are review artifacts, not a viewer manifest.

## Source and inventory

The source is Kristen Browne and Heidi Schlehlein, Human Reference Atlas / HuBMAP, [3D Reference Organ Set for Female v1.5 (2023)](https://doi.org/10.48539/HBM352.BTSQ.586), CC BY 4.0. The bundled atlas and chunk are pinned before geometry is read:

| Input | SHA-256 |
| --- | --- |
| `public/models/atlas-female.json` | `1525c07d2ed46263c086d1c6b2e52eb9b8f8985746e9a8259036bd6e1d9a1ced` |
| `public/models/female-6.bin` | `8b9d1c038bf4b7951bc19a0e30dd3c7c65c29af3f2bd9f605b8048204c66ee40` |
| Exact inventory and geometry digest | `e3715d7db0bf5b713afc688593842a35d035da09ca9b68a2cdc98ee66b503cf4` |

The exact 75 internal meshes comprise 39 left and 36 right parts. Each side has an outer cortex, renal column, and renal pelvis. Left/right counts are 11/10 papillae, 11/10 pyramids, 10/10 minor calyces, and 4/3 major calyces. Capsules, hila, and source ureters supply registration context and are excluded from the 75. Explicit IDs, expected names, urinary category, index ranges, finite positions, and source geometry hashes are checked. The target kidney IDs are `FJ3145`/`FJ3147`; retained ureters are `FJ3144`/`FJ3146`.

## Measured first pass

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
