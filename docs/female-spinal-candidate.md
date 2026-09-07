# Female spinal cord candidate

The bundled HRA assembly supplies **29 spinal cord segment meshes and 63,534 triangles**. `scripts/female_spinal_candidate.py` exports those actual meshes into an isolated candidate, retaining the source topology and recording the entire transformation. Its new continuous correction uses actual vertebral cross-sections: the original C5/C6 interference is removed on both fitting and held-out sections.

The candidate remains **unpublished** (`enabled: false`). The remaining concrete work is the overlapping upper cord–medulla junction and the source's open segment boundaries. These limitations are separate from the improved canal fit; they are not an external approval requirement.
## Source and inventory

The inputs are the checked-in HRA united-female v1.5 assembly, `public/models/atlas-female.json`, and the current `atlas-female-reconstructed.json`. The HRA [official reference-library page](https://hubmapconsortium.github.io/ccf/pages/ccf-3d-reference-library.html) describes the Visible Human female source, lists a female spinal cord object, and supplies the reference-object licensing context. Preserve **Human Reference Atlas / Visible Human Female, CC BY 4.0** attribution with exported geometry. The existing versioned [united-female digital-object link](https://lod.humanatlas.io/ref-organ/united-female/v1.5) returned 404 during this review; this experiment reads the bundled, hashed assets and does not depend on that link remaining available or claim to use the latest release.

| Reviewed input | SHA-256 |
|---|---|
| HRA source manifest | `1525c07d2ed46263c086d1c6b2e52eb9b8f8985746e9a8259036bd6e1d9a1ced` |
| Current female manifest | `4134384790ef6be9526105466ab946f4291afdc0e8237145dd0d64cd5809ffa8` |

The exact inventory is C1–C8, T1–T12, L1–L5 and S1–S4. The script preserves source IDs, including `VH_F_eigth_thoracic_spinal_cord_segment`. S5 and a coccygeal cord segment are not separately represented in this bundled inventory; no replacement geometry is invented. The report lists every selected name, ID, source system, actual vertex bounds, vertex/triangle counts, topology and raw geometry hash. It also hashes every source and target chunk read and both analysis scripts used.

All 29 source segment meshes have open boundary edges after exact-coordinate seam welding: **98–273 boundary edges per segment**, with zero nonmanifold edges and zero zero-area triangles at the stated threshold. Small nearest-surface distances between neighboring ends therefore do not establish a single closed or watertight cord. The source assembly's vertebra labels include six lumbar vertebrae; this experiment uses only L1 at the lower end and does not silently assume the extra lumbar level maps to a conventional five-level target.

## Registration and measured result

The initial registration uses corresponding **20 actual vertebral surfaces**: C1–C7, T1–T12 and L1. Their triangle-area-weighted surface centroids define a continuous displacement field with piecewise linear height mapping and fixed transverse scale. Cord segment names are never treated as same-numbered vertebral levels. Both ends use recorded linear extrapolation; **2,713 source vertices** lie beyond the centroid-height interval. The initial field has a positive Jacobian determinant, **0.887–1.195**. Mean sampled source-bone-to-target-triangle residuals range from **2.51 to 4.75 mm**. These are surface descriptors and residuals, not canal landmarks or anatomical approval.

The refinement adds five compact, twice continuously differentiable Wendland functions of height, translating only x and depth. Their coefficients are fitted directly against actual cord and vertebral triangle/plane intersections on **50 sections** across the affected C1–T2 cord segments. The deterministic search targets 0.5 mm radial clearance, penalizes bone-interior parity and loss of measured boundary coverage, limits sampled translation to 6 mm, and regularizes coefficient size. It does not replace contours, scale the cord's transverse dimensions, independently move individual segments or change vertex/triangle counts.

The fitted correction moves vertices by at most approximately **4.010 mm**. Its analytical Jacobian determinant is exactly 1; its inverse-transpose updates normals coherently with the existing registration. Height and total vertical length are unchanged. The report contains all basis centers, radii, coefficients, objective parameters and iterations. Measurements use the actual exported float32 positions.

| Actual section screen | Before correction | After correction |
|---|---:|---:|
| C5 midpoint minimum radial clearance | −1.917 mm | +0.781 mm |
| C6 midpoint minimum radial clearance | −3.271 mm | +0.486 mm |
| Fitting sections with negative clearance | 16 / 50 | 0 / 50 |
| Worst fitting-section radial clearance | −3.835 mm | +0.486 mm |
| Held-out sections with negative clearance | 17 / 60 | 0 / 60 |
| Worst held-out radial clearance | −3.436 mm | +0.275 mm |
| Held-out sections with odd bone parity | 3 / 60 | 0 / 60 |

The **60 held-out planes were not used to fit the correction**. They include near-end slices as well as intermediate slices. Every plane compares the actual outermost cord intersection against the nearest bone boundary along 64 directions. The report now distinguishes rays missing a bone intersection from those missing a cord intersection; open source contours must not be mistaken for complete boundaries. At the original 29 midpoints, 27 still have an escaping ray, and none has negative clearance or odd bone parity. Incomplete rings can occur where a horizontal plane passes between tilted vertebral arches. Neither these partial rings nor a finite set of clear samples establish complete canal containment. Negative radial clearance is a section-fit metric, not a three-dimensional penetration depth.

## Segment ends and upper medulla junction

All **28 adjacent segment pairs** retain overlapping axial bounds; there is no positive axial gap between their bounds. Each pair has a sampled nearest end-to-surface distance below **0.001 mm** before and after correction. Maximum distances in the 2 mm end strips reach approximately 2.05 mm because the strips include points away from the exact interface; they must not be presented as seam gaps. Open boundary edges remain as documented in the source inventory. The correction changes neither source topology nor height, and does not weld the segments or reconstruct missing distal labels.

The superior junction is a separate issue from cervical canal clearance. Actual triangle distances and three ray parities per closed medulla mesh establish the following:

| Upper 1 mm of first cord segment | Bundled HRA source | Corrected candidate against current medulla |
|---|---:|---:|
| Sampled strip vertices | 478 | 475 |
| Inside medulla union | 0 | 359 |
| Outside medulla union | 478 | 116 |
| Ambiguous classifications | 0 | 0 |
| Minimum / median / maximum surface distance | 2.141 / 3.090 / 4.022 mm | 0.028 / 0.923 / 3.550 mm |
| Signed axial bounds gap | +1.839 mm | −15.219 mm |

The bundled source has a real positive separation from its Allen medulla assembly. The registered candidate reaches the current BodyParts3D medulla and overlaps it, but **does not form a clean matched junction**: some superior cord vertices remain outside the medulla. The negative axial value is only overlapping height range, and the listed unsigned surface distances are not penetration depths. No cap, fake connection or source-seam repair was added to hide this mismatch.

`overlay.svg` and `baseline-overlay.svg` show actual cord triangles with the current vertebrae **and medulla** from front and side. `cervical-sections.svg` directly compares old and new C5/C6 outlines against the actual vertebral sections. The broad improved alignment supports a useful experimental visual review, but the current iteration deliberately leaves the cord out of the canonical model while the upper junction and open segmentation are refined.
## Reproduce

Run from a checkout with Python 3 and NumPy; no Blender, network download or optional scientific library is needed:

```bash
python3 scripts/female_spinal_candidate.py --output-dir /tmp/human-atlas-spinal-refined
```

The output directory must resolve outside the repository. Direct repository paths, directory symlink aliases into the repository and symlinked output files are rejected. `--models-dir` can select another directory containing both required manifests and their bundled chunk files, so a reviewed candidate framework can be compared without replacing production.

The output contains:

- `report.json`: deterministic source inventory, hashes, both transformation fields, source/baseline/corrected continuity, fitting and held-out section measurements, and source/candidate medulla interfaces; top-level `enabled: false`. `previewReadiness` records the sampled fit result separately from publication.
- `candidate.json` and `candidate.bin`: 29 original source segments in the target frame. Byte offsets describe little-endian float32 positions/normals and uint32 indices. Source topology is retained. The manifest records the binary hash and CC BY 4.0 attribution.
- `overlay.svg`, `baseline-overlay.svg` and `cervical-sections.svg`: standalone geometry reviews using actual triangles, viewable directly in a browser.

For an optional offline screenshot with locally installed Chrome:

```bash
google-chrome --headless --no-sandbox --disable-gpu --hide-scrollbars \
  --screenshot=/tmp/human-atlas-spinal-refined/overlay.png --window-size=1000,950 \
  file:///tmp/human-atlas-spinal-refined/overlay.svg
```

The numerical helper for exact point-to-triangle distances is the existing `scripts/pelvis-surface-audit.py`; its file hash is part of the report. Measurements are deterministic for the checked-in geometry and numerical runtime. Focused verification checks analytical normals against finite differences, compact support, positive Jacobians, ray parity including a duplicated seam, and all 29 exported position/normal/index buffers against the recorded field and original source topology:

```bash
python3 scripts/female-spinal-candidate.test.py --candidate-dir /tmp/human-atlas-spinal-refined
```

The aggregate runner separately verifies that production model files remain unchanged. This experiment does not establish anatomy approval or publish a viewer option.
