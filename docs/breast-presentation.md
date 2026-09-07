# Source-based breast tissue views

The single canonical female model at `/female` uses all 16 original HRA breast meshes, including both adipose envelopes. There is no alternate breast-model manifest, selector or query parameter. This replaces the previous generated concentric-ring envelopes.

## Geometry and scope

Every breast mesh retains its HRA vertex count and triangle indices. Each side uses the source adipose envelope's bounds to define one shared contour field. The base footprint is kept at source scale and translated to the chest; forward projection is reduced and the lower anterior portion lifted. Fat, glands, ducts, sinuses, supports and optional nipple surfaces all use that same per-side field. The breast assembly is placed directly in the final frame, without the old chest drape, tissue-specific insets or another application of the body morph.

The exact source hashes and settings are in `female-fit-report.json` under `breastContour`: scale 1, target centers `[±0.082, 1.21, 0.088]` m, projection-reduction factor 0.15 and maximum lower-front lift 0.010 m, moderated by the recorded smooth gates. `contoured` lists all 16 breast IDs and `regenerated` is empty. Normals follow the inverse-transpose Jacobian. Source topology and a positive sampled Jacobian do not validate anatomical placement, tissue volume, chest-wall contact or containment. The reconstruction validator also screens adipose posterior versus pectoralis anterior on a 4 mm xy grid inside each fat envelope's inner 60% ellipse. That measurement is a geometric regression bound, not attachment or containment proof.

The approach was informed by [Mahendra Beniwal's Female Atlas](https://github.com/HiMahendraBeniwal/female-atlas/tree/b0a20fa3a009ca5f6bfa07b7479312e04c862af2). This implementation uses the project's already bundled CC BY 4.0 HRA source assets; it does not copy that project's code or distribute reference-image pixels. The model remains an estimated study model.

## Chest controls

- **Tissue** shows the ten internal breast pieces, including the two contoured HRA adipose envelopes, with plain tissue materials and source-derived surface geometry. The previous procedural lobule/bump texture is removed.
- **Glands** hides the two adipose envelopes and shows the eight lobe, duct, sinus and support pieces. This is a visibility preset, not a geometric cut plane.
- **Pectorals** hides all 16 breast pieces and enables Muscles, exposing the underlying chest muscles.

The six nipple, areola and areolar-tubercle meshes remain in the separate optional **Body surface** layer. Search selection can reveal a hidden structure; isolation shows only selected structures. Explicit breast-layer controls restore Tissue mode. Reset restores the default Tissue view. Male geometry and visibility behavior are unchanged.

Breast glandular, adipose and connective tissues overlie the pectoral region; they are distinct from the skeletal muscles that move the shoulder. The [NCI SEER mammary gland overview](https://training.seer.cancer.gov/anatomy/reproductive/female/glands.html) supports this tissue distinction, not the model's estimated contour or placement.

## Validation and historical evidence

```bash
node --test scripts/hra-breast-contour.test.mjs
node --experimental-strip-types --test scripts/female-category-toggle.test.mjs scripts/chest-visibility.test.mjs
```

The geometry tests compare all 16 source triangle-index buffers, verify the unchanged posterior base and lateral coordinates relative to source translation, check the anterior/lower contour and normals, pin the reviewed trial's geometry bytes, and require all 2,227 nonbreast meshes to match the pre-promotion main baseline. The category tests verify the same 62 inserted IDs and their controls, including all 16 source-topology-preserving breast meshes.

The [older containment audit](breast-containment-audit.md), inset experiments and generated-envelope comparisons remain historical records. Their measurements describe the previous geometry and must not be used as current containment or chest-wall validation. The current HRA contour still needs anatomical review.
