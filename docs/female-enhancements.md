# Female anatomy enhancements

This branch includes the breast update merged in PR #4 and subsequent main changes. There is one female viewer at `/female`.

## Available in the viewer

The female model now includes the left and right **anterolateral knee ligaments** under **Connective tissue**. Search either full name to inspect or isolate it. The two source meshes use a recorded local bone-constrained fit followed by the existing female morph. They retain source triangle topology; all 2,243 previously displayed pieces retain their geometry and metadata. The added source-proximity patches have shallow bone overlap, documented in the [joint review](female-joint-candidates.md); these are estimated static study structures, not a validated joint-motion model.

The same male/female viewers support cursor-centred wheel zoom and bilateral muscle search shortcuts: hamstrings (8 meshes), quadriceps/quads (8), calf/calves (8), and gluteal muscles/glutes (6). Groups preserve the original source concepts and list their actual members. Missing or mismatched source meshes produce an explicitly partial group. Zoom respects hidden layers, isolation, view offsets and distance limits, with a target-plane fallback over empty space; touch controls remain available.

## Continuing geometry work

| Candidate | Latest measured progress | Viewer status |
| --- | --- | --- |
| Spinal cord, 29 segments | Held-out cervical negative-clearance sections fell from 17/60 to 0/60; medulla junction still needs refinement | Offline candidate |
| Knee ligaments and menisci, 8 meshes | Refined bone holdout p95 gaps about 1.7–3.0 mm; stronger regularization selected for ALL | ALL pair enabled; ACL, PCL and menisci deferred |
| Quadriceps tendons, 2 meshes | Patellar patch p95 improved from about 14 mm to about 4 mm; distortion/interface work remains | Offline candidate |
| Kidney internals, 75 meshes | Constrained fit improves full surface p95 to about 5.9/7.1 mm; compression and enclosure problems remain | Offline candidate |

Detailed methods and evidence: [spinal cord](female-spinal-candidate.md), [knees and tendons](female-joint-candidates.md), [kidneys](female-kidney-candidate.md). The [refined experiment summary](../data/anatomy/enhancements/summary.json) records the input model hashes for that run. Those experiments preceded the separate two-piece publication, so their baseline hashes intentionally identify the pre-addition model.

## Reproduction

```sh
# Rebuild the single female model, including its pinned ALL additions.
python3 scripts/build-female-reconstruction.py --output-dir /tmp/female-rebuilt

# Rerun all four fitting experiments without changing public assets.
npm run audit:female-enhancements -- --output-dir /tmp/female-enhancements

# Generate the selected stronger-regularization joint candidate.
python3 scripts/female_joint_candidates.py --regularization .08 --output-dir /tmp/female-joint-preview

npm run test:viewer-enhancements
npm run test:female-enhancements
```

The canonical builder first reconstructs the unchanged body/breast baseline, then appends the two selected ligament meshes using [the pinned recipe](../data/anatomy/joint-additions/report.json). The publication helper rejects changed source/report/mesh pins and independently reproduces the recorded transform. It writes a separate chunk and preserves every base chunk. The JavaScript atlas validator also reproduces the joint field independently, checks original source indices, and verifies category, normal, bound and compressed-buffer integrity.

Experiments require Python with NumPy. The aggregate runner writes reports and diagnostic geometry outside the repository, hashes all public model files before and after, and fails if an experiment changes them. A successful experiment means its computation completed, not that its output has been published. Further candidates must be deliberately selected and integrated; the runner never modifies the viewer.

HRA geometry is Kristen Browne and Heidi Schlehlein's *3D Reference Organ Set for Female v1.5* (2023), CC BY 4.0. See [source attribution](../public/ATTRIBUTION.md). Candidate selection and viewer ideas were inspired by [Mahendra Beniwal's Female Atlas](https://github.com/HiMahendraBeniwal/female-atlas/tree/b0a20fa3a009ca5f6bfa07b7479312e04c862af2); this work uses our bundled HRA assets and original fitting/interaction code.
