# Female anatomy enhancement experiments

This work follows the canonical breast update in [PR #4](https://github.com/wiiiimm/human-atlas/pull/4). It develops four additions on `feature/female-anatomy-enhancements` without creating another female viewer or reshaping the accepted body.

| Candidate | Source meshes | Fitting work |
| --- | ---: | --- |
| Spinal cord | 29 | Continuous source assembly and retained spinal canal/brainstem |
| Knee ligaments and menisci | 8 | Bilateral local femur–tibia–patella registration |
| Quadriceps tendons | 2 | Shared knee registration, patellar interface and existing muscle overlap |
| Kidney internals | 75 | Coherent per-kidney assembly, capsule fit and retained ureter interface |

Run all experiments with:

```sh
npm run audit:female-enhancements -- --output-dir /tmp/female-enhancements
```

Python dependencies are recorded by the individual scripts. The runner writes a summary and per-experiment reports outside the repository. It hashes all public model files before and after execution and fails if an experiment changes them. Candidate artifacts are research outputs; the runner does not publish them into `/female`. A successful command means the experiment completed, not that the fit was accepted.

Individual methods and results:

- [Spinal cord](female-spinal-candidate.md)
- [Knees and quadriceps tendons](female-joint-candidates.md)
- [Kidney internals](female-kidney-candidate.md)

The HRA source already bundled in `atlas-female.json` is Kristen Browne and Heidi Schlehlein's *3D Reference Organ Set for Female v1.5* (2023), CC BY 4.0. See [source attribution](../public/ATTRIBUTION.md). The candidate selection was inspired by [Mahendra Beniwal's Female Atlas](https://github.com/HiMahendraBeniwal/female-atlas/tree/b0a20fa3a009ca5f6bfa07b7479312e04c862af2); these experiments use our bundled HRA assets and original fitting code.

The source spinal cord inventory has no S5 or coccygeal segment and is not a complete nervous system. Knee meniscus IDs each contain both source components. Repeated kidney names represent distinct meshes; the left papilla/pyramid count differs from its minor-calyx count. These identities are retained rather than silently renamed or completed.

## First fitting pass

The first pass generates all 114 candidate meshes, but none is enabled in the viewer. The spinal experiment has two sampled cervical sections with negative radial clearance. The joint fit leaves held-out bone gaps up to 13.44 mm (p95), and the tendon patellar patches reach 14.43 mm p95 versus approximately 2.71 mm in the source. Kidney capsule fitting exposes donor-shape mismatch; shell topology also limits organ-containment claims. The detailed reports distinguish measured proximity from unresolved anatomical attachment and enclosure questions.

These results identify the next implementation work: canal-constrained spinal registration, joint deformation constrained by corresponding attachment surfaces, a tendon fit that resolves the retained quadriceps interface, and a coherent renal assembly/enclosure fit. Existing source or target anatomy has not been edited to conceal the residuals.

Runner isolation checks: `python3 scripts/female-enhancement-runner.test.py` verifies repository/symlink output rejection, detection of changed model assets, and propagation of failed experiments. The individual scripts also perform their geometry-specific checks.

The checked-in [first-pass evidence summary](../data/anatomy/enhancements/summary.json) links the measured reports and records all 86 unchanged public model file hashes. Temporary candidate geometry and SVGs can be regenerated with the command above.
