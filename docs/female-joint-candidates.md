# Female knee and quadriceps tendon candidates

These are **disabled candidates**, generated separately from the atlas. Neither candidate is ready to enable: a shared similarity fit leaves material joint and attachment gaps. Existing body shape, bone positions, muscles, indices and normals remain unchanged. The script verifies SHA-256 identity of all three atlas manifests, their binary chunks and the female fit report before and after generation.

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

All ten structures share their side's transform. Registration occurs in the pre-morph BodyParts3D frame; the recorded female morph is then applied exactly once. The script retains source triangle topology and recomputes area-weighted normals for exported candidate surfaces. No existing morph parameter changes.

## Measured first-pass result

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

## Required next fit work

The candidate demonstrates that joint-local placement is possible but that a single similarity transform is insufficient for these source/target joint shapes. Do not lower the scale bound or independently move ligaments merely to conceal residuals.

1. Review corresponding distal-condylar, intercondylar, tibial attachment and superior-patellar patches on actual source and target surfaces. Record triangle/barycentric annotations rather than substituting bbox corners for anatomical landmarks.
2. Fit a constrained joint deformation with independent bone-surface constraints and held-out attachment regions. Preserve the source assembly's ligament course while assessing both ends. Report distortion and the effect on menisci.
3. Refit the quadriceps tendon with the accepted joint frame, reviewing its patellar footprint and the existing distal quadriceps geometry before deciding whether the added tendon duplicates tissue already represented.
4. Review section views and actual triangle intersections/containment where valid. Current projections and unsigned distances are geometric diagnostics, not anatomical approval or movement validation.

The analytic self-check recovers a known proper similarity and verifies exact triangle-face/edge distances. Full candidate generation verifies source topology and unchanged input atlas hashes. Neither check substitutes for the unresolved registration and attachment work above. HRA source licensing remains CC BY 4.0 as documented in `public/ATTRIBUTION.md`.
