# Inserted female structures: category-toggle audit

All **62 inserted HRA-derived structures** have an existing, active control-panel category and respond to that category's visibility toggle. No missing toggle or concrete wrong primary-category assignment was found. This is a category/visibility audit, not validation of the model's geometry or completeness.

All 16 contoured breast meshes preserve their HRA source topology and are included within those 62 structures. No breast envelope is regenerated. Exact independent expectations are recorded in [female-category-expectations.json](../data/anatomy/female-category-expectations.json).

| Control-panel toggle | Inserted pieces | Scope |
| --- | ---: | --- |
| Skeleton | 8 | Bilateral ilium, ischium, pubis; sacrum; coccyx |
| Reproductive | 38 | Ovaries, uterine tubes, uterus/cervix, vagina, and associated supports/folds/recess |
| Breast tissue | 10 | Two contoured HRA adipose envelopes, two lobe assemblies, two main-duct assemblies, two sinus assemblies, two suspensory-support assemblies |
| Body surface | 6 | Bilateral nipples, areolae, and areolar tubercles |

## Category interpretation

The application uses **one primary category per mesh**. These are anatomical assembly groupings, not multiple tissue/function tags. Ovaries remain under Reproductive despite their endocrine role, and reproductive ligaments/folds remain with their organ assembly. NCI identifies the ovaries, uterine tubes, uterus and vagina as female reproductive anatomy; its ovarian description also includes the associated ligamentous attachments. [NCI female reproductive system](https://training.seer.cancer.gov/anatomy/reproductive/female/), [NCI ovaries](https://training.seer.cancer.gov/anatomy/reproductive/female/ovaries.html).

Breast glands, ducts, fat and suspensory supports belong to the Breast tissue assembly; they are not assigned to the skeletal Muscles toggle. Nipple/areola structures are optional external surfaces. NCI describes internal glandular, fatty and connective components and the external nipple/areola, supporting this distinction. The split into these two UI groups is the viewer's presentation policy rather than a claim that there is only one anatomical taxonomy. [NCI mammary glands](https://training.seer.cancer.gov/anatomy/reproductive/female/glands.html).

Three details are intentional, not missing assignments:

- All 16 breast structures, including suspensory supports, use the plain reproductive-tissue rendering material. They remain controlled by **Breast tissue** or **Body surface**. Render material grouping in `scene.tsx` does not replace `part.system` for visibility.
- Ovaries do not also appear when only **Endocrine** is enabled; reproductive and breast ligaments do not also appear when only **Connective tissue** is enabled. Supporting multiple simultaneous category memberships would be a separate feature, not an existing behavior.
- The uterovesical pouch is retained with the female reproductive assembly; its name does not make it a urinary organ. The classification is an organ-associated grouping, not a claim that the recess is glandular tissue. Source HRA assignments are preserved for all 38 reproductive entries.

No inserted female pelvic-floor muscles or female urethra exist in this model yet. A working Muscles or Urinary toggle therefore does not establish that these pending structures have been added.

## Chest controls and default behavior

| Control | Inserted breast pieces shown | Optional surface pieces | Other enabled anatomy |
| --- | ---: | ---: | --- |
| Tissue | 10 | 0 | Retained; Muscles enabled |
| Glands | 8 (only the two outer envelopes removed) | 0 | Retained; Muscles enabled |
| Pectorals | 0 | 0 | Retained; Muscles enabled |
| Breast tissue category alone | 10 | 0 | Hidden unless explicitly selected |
| Body surface category alone | 0 internal pieces | 6 | Hidden unless explicitly selected |

The Glands view deliberately includes duct, sinus and support pieces, not only gland lobes. Tissue/Glands/Pectorals are chest presets, not substitutes for the system list. Switching to Tissue or Glands enables Breast tissue and Muscles and disables Body surface; switching to Pectorals disables both breast categories and enables Muscles. Other system switches retain their state.

By default, all 56 non-surface inserted structures are enabled and the six nipple/areola surface structures are optional. Clicking the Body surface category or enabling its switch returns the chest preset to Tissue so these six can be seen. Selecting a hidden structure through search, or isolating it, intentionally reveals it even if its category is off; clearing that selection restores category and chest filtering. The system handlers clear existing selection/isolation so disabling a category works as expected.

## Validation

Reviewed the actual source and generated manifests, `female-fit-report.json`, `app/anatomy.ts`, `app/page.tsx`, `app/scene.tsx`, and the existing chest tests. The renderer and UI visible count both call `partIsVisible`; tests compare its direct and cached-set lookup paths.

```bash
node --experimental-strip-types --test scripts/female-category-toggle.test.mjs scripts/chest-visibility.test.mjs
node --test scripts/hra-breast-contour.test.mjs
```

**All 19 category, chest-visibility and breast-contour tests pass against the promoted canonical assets.** The category suite contains six focused cases plus seven chest-visibility cases. Coverage includes the exact 62-ID inventory, all 16 source-topology-preserving breast meshes, all reviewed categories and active controls, every inserted piece against every system toggle, defaults, all three chest modes, every inserted piece under selection/isolation, and unchanged male behavior. The expected category fixture is independent of the generated atlas's actual `system` field, so an erroneous reassignment cannot pass by reproducing itself in the test. Six additional breast-contour tests verify source indices, the recorded shape field, normals, reviewed geometry bytes, unchanged nonbreast meshes, canonical compressed payloads, and the adipose/pectoralis geometric clearance screen.

The previous browser audit used manifest `9b6f0fc90c80e4cd50a7ec65fa9d43ddeb0d6e5d79b6d1b9a78105443610730b`, before source-based breast promotion. Its category counts were Reproductive 38, Breast tissue 10 and Body surface 6, and each toggled to zero and back; Skeleton included all eight inserted pelvic pieces. Those observations are historical UI evidence. Current geometry and category-test results must be checked against the promoted manifest rather than inferred from that earlier browser session.
Current automated-audit manifest SHA-256: `4134384790ef6be9526105466ab946f4291afdc0e8237145dd0d64cd5809ffa8`.

## Every inserted ID

| ID | Source/display name | Primary control-panel toggle |
| --- | --- | --- |
| `VH_F_sacrum` | Sacrum | Skeleton |
| `VH_F_coccyx` | Coccyx | Skeleton |
| `VH_F_pubis_compact_bone_R` | Right pubis | Skeleton |
| `VH_F_pubis_compact_bone_L` | Left pubis | Skeleton |
| `VH_F_ilium_compact_bone_R` | Right ilium | Skeleton |
| `VH_F_ilium_compact_bone_L` | Left ilium | Skeleton |
| `VH_F_ischium_compact_bone_R` | Right ischium | Skeleton |
| `VH_F_ischium_compact_bone_L` | Left ischium | Skeleton |
| `VH_F_fat_L` | Adipose tissue of left breast | Breast tissue |
| `VH_F_fat_R` | Adipose tissue of right breast | Breast tissue |
| `VH_F_mammary_lobes_L` | mammary lobe | Breast tissue |
| `VH_F_mammary_lobes_R` | mammary lobe | Breast tissue |
| `VH_F_main_lactiferous_ducts_L` | main lactiferous duct | Breast tissue |
| `VH_F_main_lactiferous_ducts_R` | main lactiferous duct | Breast tissue |
| `VH_F_main_lactiferous_sinuses_L` | Lactiferous sinus | Breast tissue |
| `VH_F_main_lactiferous_sinuses_R` | Lactiferous sinus | Breast tissue |
| `VH_F_suspensory_ligaments_L` | Suspensory ligament of left breast | Breast tissue |
| `VH_F_suspensory_ligaments_R` | Suspensory ligament of right breast | Breast tissue |
| `VH_F_areolar_tubercles_L` | areolar tubercle | Body surface |
| `VH_F_areolar_tubercles_R` | areolar tubercle | Body surface |
| `VH_F_nipple_L` | left nipple | Body surface |
| `VH_F_nipple_R` | right nipple | Body surface |
| `VH_F_areola_L` | Left female areola | Body surface |
| `VH_F_areola_R` | Right female areola | Body surface |
| `VH_F_vagina` | vagina | Reproductive |
| `VH_F_cervicovaginal_junction` | cervicovaginal junction | Reproductive |
| `VH_F_ampulla_of_uterine_tube_R` | Ampulla of right uterine tube | Reproductive |
| `VH_F_isthmus_of_fallopian_tube_R` | Isthmus of right uterine tube | Reproductive |
| `VH_F_fibria_of_uterine_tube_R` | Fimbria of right uterine tube | Reproductive |
| `VH_F_uterine_tube_infundibulum_R` | Infundibulum of right uterine tube | Reproductive |
| `VH_F_uterine_tube_infundibulum_L` | Infundibulum of left uterine tube | Reproductive |
| `VH_F_fibria_of_uterine_tube_L` | Fimbria of left uterine tube | Reproductive |
| `VH_F_isthmus_of_fallopian_tube_L` | Isthmus of left uterine tube | Reproductive |
| `VH_F_ampulla_of_uterine_tube_L` | Ampulla of left uterine tube | Reproductive |
| `VH_F_uterovesical_pouch` | Uterovesical pouch | Reproductive |
| `VH_F_right_round_ligament_of_uterus` | Right round ligament of uterus | Reproductive |
| `VH_F_left_round_ligament_of_uterus` | Left round ligament of uterus | Reproductive |
| `VH_F_right_uterosacral_ligament` | Right uterosacral ligament | Reproductive |
| `VH_F_left_uterosacral_ligament` | Left uterosacral ligament | Reproductive |
| `VH_F_right_cardinal_ligament_of_uterus` | Cardinal ligament | Reproductive |
| `VH_F_left_cardinal_ligament_of_uterus` | Cardinal ligament | Reproductive |
| `VH_F_suspensory_ligament_of_ovary_R` | Suspensory ligament of right ovary | Reproductive |
| `VH_F_suspensory_ligament_of_ovary_L` | Suspensory ligament of left ovary | Reproductive |
| `VH_F_ovarian_ligament_L` | Left ovarian ligament | Reproductive |
| `VH_F_ovarian_ligament_R` | Right ovarian ligament | Reproductive |
| `VH_F_broad_ligament` | broad ligament of uterus | Reproductive |
| `VH_F_mesosalpinx_L` | Left mesosalpinx | Reproductive |
| `VH_F_mesosalpinx_R` | Right mesosalpinx | Reproductive |
| `VH_F_mesovarium_R` | Right mesovarium | Reproductive |
| `VH_F_mesovarium_L` | Left mesovarium | Reproductive |
| `VH_F_left_ovary` | Left ovary | Reproductive |
| `VH_F_right_ovary` | Right ovary | Reproductive |
| `VH_F_abdominal_ostium_of_uterine_tube` | Abdominal ostium of uterine tube | Reproductive |
| `VH_F_body_of_uterus` | body of uterus | Reproductive |
| `VH_F_fundus_of_uterus` | Fundus of uterus | Reproductive |
| `VH_F_cornua` | cornua | Reproductive |
| `VH_F_lower_uterine_segment` | Lower uterine segment | Reproductive |
| `VH_F_posterior_wall_of_uterus` | Posterior wall of uterus | Reproductive |
| `VH_F_anterior_wall_of_uterus` | Anterior wall of uterus | Reproductive |
| `VH_F_cervix` | uterine cervix | Reproductive |
| `VH_F_internal_cervical_os` | internal cervical os | Reproductive |
| `VH_F_external_cervical_os` | external cervical os | Reproductive |
