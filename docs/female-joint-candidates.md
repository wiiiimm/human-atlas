# Female knee and quadriceps tendon candidates

The ten-part candidate package remains **disabled as a whole** and is generated separately from the atlas. The bilateral anterolateral ligaments were subsequently selected for an additive experimental study preview, as documented below; the other eight parts remain off. A bone-constrained smooth deformation substantially improves the initial fit, but meniscal/attachment gaps and deformation strain remain. Existing body shape, bone positions, muscles, indices and normals remain unchanged. The script verifies SHA-256 identity of all three atlas manifests, their binary chunks and the female fit report before and after generation.

Run from the repository:

```sh
python3 scripts/female_joint_candidates.py --self-test
python3 scripts/female_joint_candidates.py --output-dir /tmp/female-joint-candidates
```

The output directory must resolve outside the repository, including through symlinks. It contains `report.json`, ten NPZ meshes (positions, triangle indices, recomputed normals), ten OBJ meshes, and `joint-diagnostic.svg` with frontal and sagittal projections against retained bones. The plots are marked unapproved. Raw position/index hashes, input hashes, exact IDs, transforms and measurements are in the report. NPZ archive bytes are not the reproducibility contract; their raw mesh arrays are.

## Inventory and fitting surfaces

The eight knee additions are `VH_F_anterior_cruciate_ligament_of_knee_L/R`, `VH_F_posterior_cruciate_ligament_of_knee_L/R`, `VH_F_anterolateral_ligament_of_knee_L/R` and `VH_F_meniscus_L/R`. The two tendons are `VH_F_tendon_of_quadriceps_femoris_L/R`. Each meniscus mesh has two connected components after welding; these remain together under the source ID. Component count does not establish medial/lateral labels.

The HRA source is already in meters, Y-up, with +Z anterior. The converter applies a vertical translation; an extra knee axis flip would be incorrect. These HRA knees are approximately reflected assembly counterparts, not established independent bilateral scans. The HRA assembly is not a single-person scan.

The main HRA femur leaf omits separately segmented joint regions. Fitting only that leaf falsely loses the cruciate attachment surfaces. This script uses the 15 non-cartilage leaves of each `HRA:VH_F_femur_L/R` compound, including the named condylar, intercondylar, trochlear and enthesis regions. Their exact IDs are recorded. Cartilage is excluded from the bone surface. None of these fitting bones or regions is added to the rendered candidate.

| Target | Left | Right |
|---|---|---|
| Femur | FJ3259 | FJ3365 |
| Tibia | FJ3282 | FJ3387 |
| Patella | FJ3275 | FJ3381 |

The fitting regions are the lower 90 mm of the femur, upper 80 mm of the tibia, and whole patella. Only triangles fully within a height crop are included. **These are geometric regions, not reviewed anatomical landmarks.** A deterministic translation initializes the already aligned coordinate axes. Eighteen similarity ICP iterations match sampled points to exact target triangles, with uniform scale constrained to 0.85–1.20 and proper rotation. Each bone contributes up to 72 fitting points; up to 96 disjoint odd-index points provide a sampled holdout screen. This is not complete collision validation or a global-optimum fit.

All ten structures share their side's continuous displacement field, initialized by the similarity described above. Registration occurs in the pre-morph BodyParts3D frame; the recorded female morph is then applied exactly once. The script retains source triangle topology and recomputes area-weighted normals for exported candidate surfaces. No existing morph parameter changes.

## Historical similarity baseline

Measured against reconstructed female manifest SHA-256 `4134384790ef6be9526105466ab946f4291afdc0e8237145dd0d64cd5809ffa8`. Both fits reached the **0.85 lower scale bound**, indicating that this simple coherent fit is under strain. Sampled fitting RMS decreased from 9.64/9.75 mm to approximately 5.75/5.76 mm (left/right), but local interfaces remained inadequate.

| Sampled holdout surface-distance p95 | Left | Right |
|---|---:|---:|
| Distal femur | 7.19 mm | 6.13 mm |
| Proximal tibia | 12.15 mm | 11.56 mm |
| Patella | 12.74 mm | 13.44 mm |

Attachment screens use **the same source vertices** lying within 3 mm of each source bone, then measure their distance to the retained female bone. These patches are proximity-derived, not certified entheses. All vertices in each identified patch are measured against exact triangles.

| Candidate patch-distance p95 | Left | Right |
|---|---:|---:|
| ACL → femur | 8.26 mm | 7.99 mm |
| ACL → tibia | 5.14 mm | 4.88 mm |
| PCL → femur | 6.15 mm | 5.86 mm |
| PCL → tibia | 6.04 mm | 5.94 mm |
| Anterolateral ligament → femur | 3.32 mm | 4.04 mm |
| Anterolateral ligament → tibia | 4.18 mm | 3.77 mm |
| Menisci → femur | 7.82 mm | 8.37 mm |
| Menisci → tibia | 7.01 mm | 6.87 mm |
| Quadriceps tendon → patella | **13.79 mm** | **14.43 mm** |

The source tendon patellar patches had p95 distance 2.71 mm, with 59 vertices per side. Their large candidate gaps are a concrete blocker; they must not be enabled based on a plausible overall knee location.

For the proximal tendon screen, the source-near-rectus patch has 28/27 vertices. Its distance to the union of retained rectus femoris and three vasti has candidate p95 2.32/2.17 mm. However, 80 of 174 tendon vertices per side are within 3 mm of that muscle union. This is only a proximity/possible-overlap screen: it cannot distinguish intended contact from duplicate distal tendinous geometry embedded in the existing muscle meshes. Unsigned surface distances do not prove containment or absence of intersections.

## Current bone-constrained refinement

The current script adds 16 independent local similarity ICP iterations for each bone, initialized by the coherent baseline. Each bone contributes up to 72 source surface controls. Controls are paired with exact nearest target triangle points following the local fit. One Gaussian radial-basis displacement field fits these correspondences, with 25 mm kernel radius and 0.015 diagonal regularization. It is added to the baseline mapping and applied identically to all knee/tendon meshes on that side. It fades toward the baseline away from its surface controls; there are no per-ligament translations or endpoint vertex clipping.

All control source coordinates, requested displacements and fitted coefficients are in `report.json`. Soft-tissue vertices and attachment-proximity patches are **not** fit inputs. The disjoint source bone holdout samples and the same soft-tissue source patches therefore provide independent evaluation of the automatic fitting constraints, though none is a clinically reviewed landmark set.

The default CLI writes this refined candidate to `/tmp/female-joint-candidates`, matching the command above. That report also recomputes and retains the similarity baseline under `sides.*.fit` for direct comparison. Historical similarity-only numbers in the previous section are recorded measurements, not a second output directory. The stronger-regularization preview uses `/tmp/female-joint-preview-candidates`.

| Sampled holdout p95 | Left baseline → refined | Right baseline → refined |
|---|---:|---:|
| Distal femur | 7.19 → 2.11 mm | 6.13 → 2.64 mm |
| Proximal tibia | 12.15 → 2.88 mm | 11.56 → 2.99 mm |
| Patella | 12.74 → 2.03 mm | 13.44 → 1.70 mm |

Every reported source-proximity attachment patch p95 improved:

| Patch p95 | Left baseline → refined | Right baseline → refined |
|---|---:|---:|
| ACL → femur | 8.26 → 3.94 mm | 7.99 → 4.47 mm |
| ACL → tibia | 5.14 → 4.43 mm | 4.88 → 3.96 mm |
| PCL → femur | 6.15 → 2.29 mm | 5.86 → 1.69 mm |
| PCL → tibia | 6.04 → 3.26 mm | 5.94 → 3.11 mm |
| Anterolateral ligament → femur | 3.32 → 1.84 mm | 4.04 → 2.09 mm |
| Anterolateral ligament → tibia | 4.18 → 2.75 mm | 3.77 → 2.65 mm |
| Menisci → femur | 7.82 → 6.39 mm | 8.37 → 5.78 mm |
| Menisci → tibia | 7.01 → 3.07 mm | 6.87 → 2.90 mm |
| Quadriceps tendon → patella | 13.79 → 4.40 mm | 14.43 → 3.98 mm |
| Tendon source-near-rectus patch → target quadriceps | 2.32 → 2.07 mm | 2.17 → 2.11 mm |

This is a measurable fitting improvement, not anatomical approval. The refinement's analytical Jacobian was evaluated at every candidate vertex and triangle centroid: no nonpositive determinant was found, and the minimum determinant was 0.357. Local singular values across these samples ranged from 0.311 to 2.047, indicating substantial compression/stretch despite positive orientation. These are pre-morph field measurements. The report separately records actual final triangle-edge length ratios against both source and the similarity candidate, including the unchanged female morph. Sampling does not prove global injectivity or rule out intersections between nonadjacent surfaces.

Source topology, all original atlas bytes and body shape remain unchanged. Candidate normals are recomputed from the deformed triangles. The self-test additionally checks the analytical RBF Jacobian against finite differences on a synthetic three-bone fixture and verifies identical mapping for coincident points.

The remaining meniscal femoral gaps, roughly 4 mm tendon patellar gaps, high local strain and unresolved tendon duplication semantics keep both candidates disabled. The next iteration should refine surface correspondences and attachment constraints, not relax the reported distances or hide them with independent mesh offsets.

## Required next fit work

The coherent refinement demonstrates improved joint-local placement while confirming that donor-shape and interface differences still require reviewed constraints. Do not lower the scale bound or independently move ligaments merely to conceal residuals.

1. Review corresponding distal-condylar, intercondylar, tibial attachment and superior-patellar patches on actual source and target surfaces. Record triangle/barycentric annotations rather than substituting bbox corners for anatomical landmarks.
2. Refine the current constrained deformation using reviewed bone-surface correspondences and held-out attachment regions. Preserve the source assembly's ligament course while assessing both ends. Report distortion and the effect on menisci.
3. Refit the quadriceps tendon with the accepted joint frame, reviewing its patellar footprint and the existing distal quadriceps geometry before deciding whether the added tendon duplicates tissue already represented.
4. Review section views and actual triangle intersections/containment where valid. Current projections and unsigned distances are geometric diagnostics, not anatomical approval or movement validation.

The analytic self-check recovers a known proper similarity and verifies exact triangle-face/edge distances. Full candidate generation verifies source topology and unchanged input atlas hashes. Neither check substitutes for the unresolved registration and attachment work above. HRA source licensing remains CC BY 4.0 as documented in `public/ATTRIBUTION.md`.


## Selected experimental study preview: bilateral ALL only

A bounded comparison retained the same bone correspondences and 25 mm Gaussian radius while increasing ridge regularization from 0.015 to **0.08**. A separate 35 mm-radius trial did not offer a better overall tradeoff for the selected parts. The default CLI still reproduces the historical 0.015 full candidate; the selected preview is explicit:

```sh
python3 scripts/female_joint_candidates.py --regularization .08 --output-dir /tmp/female-joint-preview-candidates
python3 scripts/female_joint_preview_review.py --self-test
python3 scripts/female_joint_preview_review.py --candidate-dir /tmp/female-joint-preview-candidates
```

Only `VH_F_anterolateral_ligament_of_knee_L` and `VH_F_anterolateral_ligament_of_knee_R` were selected for the canonical study preview. Candidate generation itself never publishes them. Its full-package `enabled:false` remains intentional. The append-only integration must use this explicit two-ID allowlist, preserving all existing meshes.

The selected fit report SHA-256 is `dd1e16dcbdb49d497128c2db7bad915cba2726fa3362bcc9bf8d6d11cb0ffdd1`; the independent `all-preview-review.json` SHA-256 is `21144342ed06133ef9b634f47e7182aa1f273478a2da56257596cc64aa89cca4`. The latter pins the fit report, both exported NPZ files and original atlas inputs, and verifies exact source triangle topology, finite positions/normals, unit normals consistent with the deformed surface, reported bounds and unchanged original atlas bytes.

| Selected ALL measurement | Left | Right |
|---|---:|---:|
| Source femoral proximity patch p95 | 2.82 mm | 2.88 mm |
| Selected femoral proximity patch p95 | 2.18 mm | 2.47 mm |
| Source tibial proximity patch p95 | 2.79 mm | 2.84 mm |
| Selected tibial proximity patch p95 | 2.54 mm | 2.40 mm |
| Sampled minimum Jacobian determinant | 0.523 | 0.480 |
| Sampled local singular-value range | 0.548–1.166 | 0.461–1.209 |

Source-versus-fitted frontal and sagittal inspection found continuous lateral ribbon shapes with recognizable source contours. The source comparison PNG is `/tmp/female-joint-preview-candidates/pcl-all-source-comparison.png`; the reproducible candidate SVG remains `joint-diagnostic.svg`. This is a static anatomy study interpretation, not an attachment-exact or kinematically validated ligament model.

The explicit intersection review does **not** claim intersection-free placement. It tests noncoplanar triangle crossings in both edge/face directions and classifies every vertex and triangle centroid against closed target bones using multiple ray directions. Target femur/tibia closure is checked after position welding. Coplanar/endpoint-only contacts are outside the triangle-crossing count; sampled penetration maxima are not continuous global maxima.

| ALL intersection screen | Left | Right |
|---|---:|---:|
| Candidate femoral crossing triangles | 27 | 22 |
| Source femoral crossing triangles | 8 | 8 |
| Candidate sampled maximum femoral penetration | 1.31 mm | 0.92 mm |
| Candidate tibial crossing triangles | 4 | 21 |
| Source tibial crossing triangles | 0 | 21 |
| Candidate sampled maximum tibial penetration | No inside samples | 0.54 mm |
| Source sampled maximum tibial penetration | No inside samples | 0.63 mm |

In particular, the four left tibial crossing triangles mean that zero interior vertex/centroid samples is **not** evidence of zero surface intersections. Source femoral inside depth is unavailable because its segmented non-cartilage surface has open boundaries; only source crossing counts and unsigned surface distances are reported there. The shallow contact-region overlap remains a documented approximation of the study preview.

### Why the other structures remain off

- **PCL:** coherent visual shape and mostly source-level attachment gaps, but additional tibial embedding remains. The source already has sampled maximum tibial penetration 3.48/3.45 mm; the smoother candidate has 4.84/4.34 mm. These are not comparisons to an arbitrary zero target. The additional approximately 1.36/0.89 mm and increased embedded sample count were sufficient to defer this pair while enabling the stronger ALL result.
- **ACL:** larger local stretching and residual femoral/tibial patch gaps remain than for the selected ALL pair; it was not selected in this bounded refinement.
- **Menisci:** the 5.78–6.39 mm refined femoral proximity-patch p95 is compared with **the same source patches at 2.88–2.89 mm**, rather than requiring the entire meniscus to touch bone. The entire source meniscus has median femoral distance about 6.50 mm and p95 about 18.59 mm; those are different populations. Source meniscus-to-cartilage all-vertex p95 is also about 18.52 mm, so a single uniform 'expected cartilage gap' cannot explain or excuse the added near-patch residual. Correct component placement and donor joint-shape fitting remain unresolved.
- **Quadriceps tendons:** distal patch gaps remain larger than source, and proximal contact cannot distinguish joining from duplication of tendinous geometry already embedded in existing muscle meshes. Their larger local strain further argues against selecting them now.

The selection is based on source-relative geometric fit, actual shape inspection and explicitly measured limitations for the current study model. It does not require or imply external clinical approval.
