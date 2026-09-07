# Female spinal cord candidate

The bundled HRA assembly supplies **29 spinal cord segment meshes and 63,534 triangles**. `scripts/female_spinal_candidate.py` now exports those actual meshes into an isolated registration candidate, with source hashes, normals, unchanged source triangle indices, continuity measurements and a front/side geometry overlay. It does not change either displayed atlas. The report and candidate manifest both use `enabled: false` because the first measured fit has unresolved cervical interference.

## Source and inventory

The inputs are the checked-in HRA united-female v1.5 assembly, `public/models/atlas-female.json`, and the current `atlas-female-reconstructed.json`. The HRA [official reference-library page](https://hubmapconsortium.github.io/ccf/pages/ccf-3d-reference-library.html) describes the Visible Human female source, lists a female spinal cord object, and supplies the reference-object licensing context. Preserve **Human Reference Atlas / Visible Human Female, CC BY 4.0** attribution with exported geometry. The existing versioned [united-female digital-object link](https://lod.humanatlas.io/ref-organ/united-female/v1.5) returned 404 during this review; this experiment reads the bundled, hashed assets and does not depend on that link remaining available or claim to use the latest release.

| Reviewed input | SHA-256 |
|---|---|
| HRA source manifest | `1525c07d2ed46263c086d1c6b2e52eb9b8f8985746e9a8259036bd6e1d9a1ced` |
| Current female manifest | `4134384790ef6be9526105466ab946f4291afdc0e8237145dd0d64cd5809ffa8` |

The exact inventory is C1–C8, T1–T12, L1–L5 and S1–S4. The script preserves source IDs, including `VH_F_eigth_thoracic_spinal_cord_segment`. S5 and a coccygeal cord segment are not separately represented in this bundled inventory; no replacement geometry is invented. The report lists every selected name, ID, source system, actual vertex bounds, vertex/triangle counts, topology and raw geometry hash. It also hashes every source and target chunk read and both analysis scripts used.

All 29 source segment meshes have open boundary edges after exact-coordinate seam welding: **98–273 boundary edges per segment**, with zero nonmanifold edges and zero zero-area triangles at the stated threshold. Small nearest-surface distances between neighboring ends therefore do not establish a single closed or watertight cord. The source assembly's vertebra labels include six lumbar vertebrae; this experiment uses only L1 at the lower end and does not silently assume the extra lumbar level maps to a conventional five-level target.

## Registration and measured result

The fit uses corresponding **20 actual vertebral surfaces**: C1–C7, T1–T12 and L1. Their triangle-area-weighted surface centroids define a continuous displacement field with piecewise linear height mapping and fixed transverse scale. Cord segment names are never treated as same-numbered vertebral levels. Both ends use recorded linear extrapolation; **2,713 source vertices** lie beyond the fitted centroid-height interval. Inverse-transpose Jacobians transform normals. Candidate measurements use the same float32 positions exported to the binary.

This is an exploratory registration using reproducible surface descriptors, not identified canal landmarks. Its Jacobian determinant remains positive, **0.887–1.195**, so the field preserves local orientation. Mean sampled source-bone-to-target-triangle residuals range from **2.51 to 4.75 mm** across the 20 bones. Those residuals measure surface agreement, not spinal cord attachment or anatomical correctness.

Adjacent segment end strips are sampled in both directions against actual opposite triangles. Each interface's closest sampled distance is below **0.001 mm** before and after fitting. The full report also includes signed axial-bound gaps, median and maximum sampled distances. These checks establish near-contact at sampled locations; they do not prove complete seam alignment, rule out overlap, or repair the open segments.

The canal screen intersects actual cord and vertebral triangles with a horizontal plane at each cord segment's midpoint. It compares the outermost cord intersection against the first bone boundary on 64 radial directions. The first candidate has:

| Sampled cord level | Worst radial clearance | Directions with negative clearance |
|---|---:|---:|
| C5 | approximately −1.92 mm | 21 / 64 |
| C6 | approximately −3.27 mm | 23 / 64 |

A negative clearance means the sampled cord outline passes the nearest bone boundary along that ray. It is a fit warning, **not a signed three-dimensional penetration-depth measurement**. No candidate cross-section center produced odd bone-intersection parity, but **27 of 29** sections have at least one escaping ray. Horizontal planes can pass through spaces between tilted vertebral arches; these incomplete rings cannot establish canal containment and do not by themselves prove displacement outside the canal. The source assembly also has three sections with negative radial clearance and 22 incomplete ray fans, so its own assembly must not be treated as a validated collision-free baseline.

`overlay.svg` draws all actual candidate cord and target vertebral triangles from the front and side. The broad alignment is plausible, while the measured cervical conflict remains unresolved. The next fitting step is to identify and review actual canal reference positions around the cervical conflict, then check more planes or full triangle intersections. A production import should also resolve segment-end boundaries and explain the absent distal labels. The current candidate remains isolated.

## Reproduce

Run from a checkout with Python 3 and NumPy; no Blender, network download or optional scientific library is needed:

```bash
python3 scripts/female_spinal_candidate.py --output-dir /tmp/human-atlas-spinal-candidate
```

The output directory must resolve outside the repository. Direct repository paths, directory symlink aliases into the repository and symlinked output files are rejected. `--models-dir` can select another directory containing both required manifests and their bundled chunk files, so a reviewed candidate framework can be compared without replacing production.

The output contains:

- `report.json`: deterministic source inventory, hashes, field anchors, Jacobians, bone residuals, both continuity analyses and both cross-section screens; top-level `enabled: false`.
- `candidate.json` and `candidate.bin`: 29 original source segments in the target frame. Byte offsets describe little-endian float32 positions/normals and uint32 indices. Source topology is retained. The manifest records the binary hash and CC BY 4.0 attribution.
- `overlay.svg`: standalone two-view visual review using actual triangles, viewable directly in a browser.

For an optional offline screenshot with locally installed Chrome:

```bash
google-chrome --headless --no-sandbox --disable-gpu --hide-scrollbars \
  --screenshot=/tmp/human-atlas-spinal-candidate/overlay.png --window-size=1000,950 \
  file:///tmp/human-atlas-spinal-candidate/overlay.svg
```

The numerical helper for exact point-to-triangle distances is the existing `scripts/pelvis-surface-audit.py`; its file hash is part of the report. Measurements are deterministic for the checked-in geometry and numerical runtime. The focused verification checks the analytical normal transform against finite differences, a ray/segment intersection with known distance, output-path rejection and preservation of the production model files. This experiment does not establish anatomy approval or publish a viewer option.
