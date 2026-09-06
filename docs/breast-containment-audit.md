# Breast assembly: containment audit and partial inset correction

[SWR-516](https://linear.app/stealth-company/issue/SWR-516) remains **In Progress**. A geometric placement defect has been partially corrected: the six internal lobe, main duct, and sinus meshes used an extra **−8 mm** posterior inset before the shared body morph. They now use **−3 mm**. Source vertex counts and triangle indices are unchanged; nipple/areola, suspensory supports, and regenerated outer envelopes are unchanged.

This is an improvement to an estimated assembly, not a completed or anatomically validated breast model. There are still **3,805 confidently outside internal vertices and 7,748 outside triangle centroids**, so further registration work is required.

## Interpretation and measurements

The regenerated `VH_F_fat_L/R` meshes serve as the **outer illustrative breast envelope** in this viewer. They are not occupied-volume segmentations of interlobar adipose tissue. Screening whether glands and ducts lie inside that outer presentation envelope does **not** imply that glands should be histologically contained inside adipose tissue.

The audit measures every internal vertex and every triangle centroid against actual envelope triangles using Blender's BVH:

- Nearest-triangle distances identify a 0.2 mm boundary uncertainty band.
- Five oblique ray directions classify interior/exterior. At least three non-ambiguous rays must agree; edge/vertex hits, nearly tangential hits, insufficient votes, and disagreements are uncertain.
- Envelope boundary edges, non-manifold edges, inconsistent winding, and degenerate triangles are checked. Invalid topology prevents confident containment classification.
- Suspensory supports are excluded because they connect beyond the envelope toward the chest. External nipple, areola, and areolar tubercle structures are excluded from the internal-tissue test.

Closed, consistently oriented topology alone does not prove an envelope has no self-intersections. Vertex/centroid samples do not guarantee containment of every triangle interior or volume. Thresholds are engineering tolerances, not clinical criteria. Duct endpoints also require anatomical interface review.

## Before and after

All 9,402 confidently outside baseline vertices were nearest the verified **back cap** of the procedural envelope; none were nearest the front or rim. Some lobe samples were more than 10 mm behind that envelope surface. This localized the defect to posterior placement rather than an exposed breast contour problem.

| Internal structure | Baseline outside vertices | Current outside vertices | Baseline outside centroids | Current outside centroids |
| --- | ---: | ---: | ---: | ---: |
| Left lobes | 4,244 | 2,022 | 8,730 | 4,139 |
| Left main ducts | 488 | 79 | 940 | 156 |
| Left sinuses | 0 | 0 | 0 | 0 |
| Right lobes | 4,200 | 1,638 | 8,579 | 3,333 |
| Right main ducts | 470 | 66 | 936 | 120 |
| Right sinuses | 0 | 0 | 0 | 0 |
| **Total** | **9,402** | **3,805** | **19,185** | **7,748** |

This reduces confidently outside vertex and centroid counts by about **60%**. Both sinus meshes have every current vertex and centroid confidently inside the envelope. Remaining outside samples in the offline candidate were all nearest the back cap; there were no new front or rim exits.

The production candidate was compared against **all 2,243 parts**. Exactly the six approved internal meshes changed position, all triangle indices are identical, every other part's positions are exactly unchanged, and the whole-body morph is unchanged. The largest position difference from the offline candidate is **0.00000626 mm**. One left-lobe vertex that was outside in the offline calculation becomes ray-uncertain after float32 serialization; it is retained as uncertain, producing 3,805 rather than 3,806 confidently outside vertices. All outside centroid counts match the offline prediction.

## Why −3 mm was selected

An isolated sweep reduced the extra inset in 1 mm steps. It preserved each source mesh and applied the change before the recorded whole-body morph; no clipping or per-vertex squeezing was used.

| Extra inset before morph | Outside vertices | Vertices posterior to projected chest surface |
| --- | ---: | ---: |
| −8 mm (baseline) | 9,402 | 14,034 |
| −7 mm | 8,248 | 12,854 |
| −6 mm | 7,072 | 11,747 |
| −5 mm | 5,847 | 10,717 |
| −4 mm | 4,835 | 9,618 |
| **−3 mm (selected candidate)** | **3,806** | **8,436** |
| −2 mm | 3,023 | 7,381 |
| −1 mm | 3,385 | 6,321 |
| 0 mm | 4,400 | 5,300 |

These are **offline candidate** vertex measurements. The chest screen casts a ray through each unchanged x/y position against the actual frontmost selected chest musculoskeletal triangles, using the same source-part selection as the builder. It does not use the builder's blurred or nearest-filled grid. Being posterior to that surface is a projected geometric discrepancy, not a closed-volume muscle-penetration diagnosis; the selected meshes do not establish complete chest anatomy.

Removing the inset entirely improves the posterior position but produces **1,932 outside sinus vertices**, where the baseline had none. The −2 mm candidate has fewer aggregate outside samples, but its duct counts worsen relative to −3 mm and 190 sinus vertices enter the boundary band. The −3 mm candidate improves each lobe and duct group, keeps both sinuses confidently inside, and introduces no front/rim escapes in the sampled candidate geometry.

For ducts and sinuses, the nearest 5% of baseline vertices to each unchanged nipple and areola were also tracked. Every patch's median and p95 surface distance improves with the selected candidate; nipple patch medians fall from roughly 3.5–3.9 mm to 1.3–1.5 mm. These are reproducible geometric patches, **not annotations of actual duct openings**. Better proximity does not prove continuity or absence of overlap.

## Reproduce and preserve evidence

Use Blender 4.0.2 as provisioned by [setup-blender-review.sh](../scripts/setup-blender-review.sh), or a compatible Blender installation with NumPy and `mathutils.bvhtree`:

```bash
blender --background --factory-startup --disable-autoexec --python-exit-code 1 \
  --python scripts/breast-containment-test.py
blender --background --factory-startup --disable-autoexec --python-exit-code 1 \
  --python scripts/breast-containment-audit.py -- \
  --output data/anatomy/breast-containment-audit.json
```

Nine analytic tests cover known inside/outside/boundary points, ambiguous vertex/edge rays, broken topology, inconsistent orientation, ray disagreement/insufficient votes, separate uncertainty counts, and exact equivalence between an offline inset adjustment and applying that adjustment before the shared morph.

Evidence artifacts:

- [Immutable baseline](../data/anatomy/breast-containment-baseline.json): captured before the correction; do not replace it with a new audit. Its source manifests/binary hashes match the pre-correction model assets at commit `81e8651`.
- [Current audit](../data/anatomy/breast-containment-audit.json): regenerate after future geometry changes.
- [Verified post-correction audit](../data/anatomy/breast-containment-after.json): actual candidate output measured before publication.
- [Full inset sweep](../data/anatomy/breast-inset-experiment.json): comparisons against the preserved baseline.
- [Conservative candidate validation](../data/anatomy/breast-candidate-validation.json): every vertex and centroid, boundary localization, and external interface distances.
- [All-part production comparison](../data/anatomy/breast-production-compare.json): verifies only the six intended meshes changed and the generated positions match the offline candidate.

To audit a separately generated candidate, pass `--root /path/to/candidate` and an explicit `--output`; that root must contain `public/models/atlas-female-reconstructed.json` and its referenced binary chunks. Offline sweep and candidate validation also accept `--root`, `--baseline`, and `--output`. Their root must contain the **preserved −8 mm baseline assets**, not the corrected model; input hashes reject stale or mismatched baselines. Production comparison accepts `--before-root`, `--after-root`, and `--output`.

## Still required

Review the remaining posterior lobe/duct placements, the relationship to actual chest tissues, and the envelope dimensions together. Establish source-supported duct/nipple interfaces and suspensory attachments rather than optimizing containment alone. Preserve the source anatomy and document any future local reshaping or envelope changes. This issue cannot be completed from improved containment counts or breast appearance alone.
