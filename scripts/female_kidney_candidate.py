#!/usr/bin/env python3
"""Read-only HRA kidney-internal registration experiment; NumPy only.

Writes review artifacts to --output-dir. It never modifies a viewer atlas.
All distances are engineering surface/proximity screens, not clinical landmarks.
"""
import argparse
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import zipfile
import numpy as np
from female_kidney_registration import radial_outer_patch, fit_surface_affine

ROOT = Path(__file__).resolve().parents[1]
SOURCE_SHA = '1525c07d2ed46263c086d1c6b2e52eb9b8f8985746e9a8259036bd6e1d9a1ced'
SOURCE_CHUNKS = {'/models/female-6.bin': '8b9d1c038bf4b7951bc19a0e30dd3c7c65c29af3f2bd9f605b8048204c66ee40'}
INVENTORY_SHA = 'e3715d7db0bf5b713afc688593842a35d035da09ca9b68a2cdc98ee66b503cf4'
spec = importlib.util.spec_from_file_location('surface_audit', ROOT/'scripts/pelvis-surface-audit.py')
audit = importlib.util.module_from_spec(spec); spec.loader.exec_module(audit)


def inventory():
    entries = []
    for side, papillae, pyramids, minors, majors in [('L', 11, 11, 10, 4), ('R', 10, 10, 10, 3)]:
        for stem, name in [('outer_cortex_of_kidney', 'outer cortex of kidney'), ('renal_column', 'renal column'), ('renal_pelvis', 'renal pelvis')]:
            entries.append((f'VH_F_{stem}_{side}', name, side))
        for stem, name, count in [('renal_papilla', 'renal papilla', papillae), ('renal_pyramid', 'renal pyramid', pyramids), ('minor_calyx', 'minor calyx', minors), ('major_calyx', 'major calyx', majors)]:
            entries += [(f'VH_F_{stem}_{side}_{chr(97+i)}', name, side) for i in range(count)]
    assert len(entries) == 75 and len({e[0] for e in entries}) == 75
    return entries


def quantiles(values):
    values = np.asarray(values)
    if not len(values): return None
    return dict(zip(['minMm', 'medianMm', 'p95Mm', 'maxMm'], [float(v)*1000 for v in np.quantile(values, [0, .5, .95, 1])]))


def topology(vertices, faces):
    directed = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    edges, inverse, counts = np.unique(np.sort(directed, axis=1), axis=0, return_inverse=True, return_counts=True)
    orientations = np.bincount(inverse, weights=np.where(directed[:, 0] < directed[:, 1], 1, -1))
    parent = np.arange(len(vertices))
    def find(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]; index = parent[index]
        return index
    for face in faces:
        for index in [1, 2]: parent[find(face[index])] = find(face[0])
    components = len({find(int(index)) for index in np.unique(faces)})
    triangles = vertices[faces]
    volume = float(np.sum(np.einsum('ij,ij->i', triangles[:, 0], np.cross(triangles[:, 1], triangles[:, 2])))/6)
    closed = bool(np.all(counts == 2) and np.all(orientations == 0))
    return dict(triangles=len(faces), connectedComponents=components, boundaryEdges=int(np.count_nonzero(counts == 1)), nonmanifoldEdges=int(np.count_nonzero(counts > 2)), orientationImbalanceEdges=int(np.count_nonzero(orientations)), closedOrientedManifoldEdges=closed, signedTriangleVolumeIntegralMl=volume*1e6, volumeIntegralDescribesClosedMaterial=closed, boundingBoxVolumeMl=float(np.prod(np.ptp(vertices, axis=0)))*1e6, warning='A closed tissue-material shell need not enclose the full organ. The signed integral is not a valid enclosed volume when boundary edges exist.')


def samples(vertices, maximum):
    return np.unique(np.linspace(0, len(vertices)-1, min(maximum, len(vertices)), dtype=int))


def transform(vertices, matrix, offset):
    return vertices @ matrix.T + offset


def nearest_vertex_indices(points, target):
    result = []
    for start in range(0, len(points), 128):
        distance = np.sum((points[start:start+128, None]-target[None])**2, axis=2)
        result.extend(np.argmin(distance, axis=1))
    return np.array(result)


def similarity(source, target, iterations=25):
    """Deterministic similarity ICP; vertices generate proposals, triangles score them."""
    source_center = (source.min(0)+source.max(0))/2
    target_center = (target.min(0)+target.max(0))/2
    scale = float(np.linalg.norm(np.ptp(target, axis=0))/np.linalg.norm(np.ptp(source, axis=0)))
    matrix = np.eye(3)*scale; offset = target_center-source_center*scale
    sample = source[samples(source, 512)]
    for _ in range(iterations):
        current = transform(sample, matrix, offset)
        paired = target[nearest_vertex_indices(current, target)]
        x = sample-sample.mean(0); y = paired-paired.mean(0)
        u, singular, vt = np.linalg.svd(x.T @ y)
        sign = np.eye(3); sign[2, 2] = np.linalg.det(vt.T @ u.T)
        rotation = vt.T @ sign @ u.T
        scale = float(np.sum(singular*np.diag(sign))/np.sum(x*x))
        matrix = scale*rotation; offset = paired.mean(0)-sample.mean(0) @ matrix.T
    return matrix, offset


def containment_votes(points, vertices, faces):
    """Three non-axis ray-parity screens; inconsistent votes remain uncertain.

    Duplicate hits on shared triangle edges are merged, not double-counted.
    These votes are only heuristic for open shells and nonmanifold source meshes.
    """
    triangles = vertices[faces]; a = triangles[:, 0]; e1 = triangles[:, 1]-a; e2 = triangles[:, 2]-a
    votes = np.zeros(len(points), dtype=int)
    for direction in [np.array([1., .371, .193]), np.array([.237, 1., .419]), np.array([.317, .173, 1.])]:
        direction /= np.linalg.norm(direction)
        h = np.cross(direction, e2); determinant = np.sum(e1*h, axis=1)
        valid = abs(determinant) > 1e-14
        inv = np.divide(1, determinant, out=np.zeros_like(determinant), where=valid)
        for i, point in enumerate(points):
            s = point-a; u = inv*np.sum(s*h, axis=1); q = np.cross(s, e1)
            v = inv*(q @ direction); distance = inv*np.sum(e2*q, axis=1)
            hits = np.sort(distance[valid & (u >= -1e-10) & (v >= -1e-10) & (u+v <= 1+1e-10) & (distance > 1e-9)])
            unique = int(bool(len(hits))) + int(np.count_nonzero(np.diff(hits) > 1e-7))
            votes[i] += unique % 2
    return votes


def require_output_directory(output):
    output = Path(output)
    if not output.is_absolute():
        output = Path.cwd() / output
    current = Path(output.anchor)
    for part in output.parts[1:]:
        current = current / part
        if current.is_symlink():
            raise ValueError('Output directory must not contain symbolic-link path components.')
    output = current.resolve()
    if ROOT == output or ROOT in output.parents:
        raise ValueError('Review artifacts must be outside the repository; use a /tmp output directory.')
    return output


def write_npz(path, arrays):
    # Fixed ZIP timestamps and sorted names make the review geometry reproducible.
    with zipfile.ZipFile(path, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for name, array in sorted(arrays.items()):
            buf = io.BytesIO(); np.save(buf, array, allow_pickle=False)
            info = zipfile.ZipInfo(name+'.npy', date_time=(1980, 1, 1, 0, 0, 0)); info.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(info, buf.getvalue())


def run(output):
    output = require_output_directory(output)
    output.mkdir(parents=True, exist_ok=True)
    filenames = ['report.json'] + [f'{side}-{fit}.npz' for side in ['L', 'R'] for fit in ['translation', 'similarity_icp', 'affine_bounds_diagnostic', 'surface_affine_constrained']]
    if any((output/name).is_symlink() for name in filenames):
        raise ValueError('Review artifact destinations must not be symbolic links.')
    source = audit.Atlas(ROOT, 'atlas-female.json'); target = audit.Atlas(ROOT, 'atlas-female-reconstructed.json')
    if audit.sha(source.path) != SOURCE_SHA: raise ValueError('HRA source manifest changed; re-review the pinned inventory before proceeding.')
    for url, expected in SOURCE_CHUNKS.items():
        if audit.sha(ROOT/'public'/url.lstrip('/')) != expected: raise ValueError(f'HRA source chunk changed: {url}')
    entries = inventory(); ids = {e[0] for e in entries}
    stems = ['VH_F_outer_cortex_of_kidney_', 'VH_F_renal_column_', 'VH_F_renal_papilla_', 'VH_F_renal_pyramid_', 'VH_F_minor_calyx_', 'VH_F_major_calyx_', 'VH_F_renal_pelvis_']
    assert {key for key in source.parts if any(key.startswith(stem) for stem in stems)} == ids
    meshes = {}; source_inventory = []; inventory_hash = hashlib.sha256()
    for key, name, side in entries:
        p = source.parts[key]
        assert p['name'] == name and p['system'] == 'urinary', key
        assert key not in target.parts, f'{key}: candidate ID already present in canonical model'
        vertices, faces = source.mesh([key]); assert np.isfinite(vertices).all()
        assert faces.min() >= 0 and faces.max() < len(vertices)
        normals = np.frombuffer(source.buffers[p['chunk']], '<i2', p['vertexCount']*3, p['normals']).reshape(-1, 3).astype(float)/32767
        geometry_hash = hashlib.sha256(vertices.astype('<f4').tobytes()+normals.astype('<f8').tobytes()+faces.astype('<u4').tobytes()).hexdigest()
        row = dict(id=key, name=name, side=side, vertices=len(vertices), triangles=len(faces), geometrySha256=geometry_hash)
        inventory_hash.update(json.dumps(row, sort_keys=True).encode())
        source_inventory.append(row); meshes[key] = (vertices, faces, normals)
    assert inventory_hash.hexdigest() == INVENTORY_SHA, 'Pinned75-part inventory/geometry hash changed'
    original_target_sha = audit.sha(target.path)
    side_results = []
    for side, target_id, ureter_id, source_ureter in [('L', 'FJ3145', 'FJ3144', 'VH_F_left_ureter'), ('R', 'FJ3147', 'FJ3146', 'VH_F_right_ureter')]:
        capsule_id = f'VH_F_kidney_capsule_{side}'
        capsule, capsule_faces = source.mesh([capsule_id]); kidney, kidney_faces = target.mesh([target_id])
        cortex, cortex_faces = source.mesh([f'VH_F_outer_cortex_of_kidney_{side}'])
        hilum, _ = source.mesh([f'VH_F_hilum_of_kidney_{side}'])
        ureter, ureter_faces = target.mesh([ureter_id]); src_ureter, src_ureter_faces = source.mesh([source_ureter])
        pelvis, pelvis_faces, _ = meshes[f'VH_F_renal_pelvis_{side}']
        distances = audit.Surface(src_ureter, src_ureter_faces).distances(pelvis)
        patch = distances <= .002
        # Preserve exact source-selected indices. If disconnected, report that rather
        # than invent an annotated endpoint or substitute a convenient arbitrary patch.
        patch_indices = np.flatnonzero(patch)
        target_surface = audit.Surface(kidney, kidney_faces); ureter_surface = audit.Surface(ureter, ureter_faces)
        query = []; query_ids = []
        for key, _, current_side in entries:
            if current_side == side:
                vertices = meshes[key][0]; index = samples(vertices, 24)
                query.extend(vertices[index]); query_ids.extend([(key, int(i)) for i in index])
        query = np.array(query); source_votes = containment_votes(query, capsule, capsule_faces)
        cortex_votes = containment_votes(query, cortex, cortex_faces)
        source_center = (capsule.min(0)+capsule.max(0))/2; target_center = (kidney.min(0)+kidney.max(0))/2
        diagonal = np.ptp(kidney, axis=0)/np.ptp(capsule, axis=0)
        matrix, offset = similarity(capsule, kidney)
        candidates = [
            ('translation', np.eye(3), target_center-source_center),
            ('similarity_icp', matrix, offset),
            ('affine_bounds_diagnostic', np.diag(diagonal), target_center-source_center*diagonal),
        ]
        exterior_ids = radial_outer_patch(capsule, capsule_faces)
        exterior_faces = capsule_faces[exterior_ids]
        print(side, 'radial exterior patch', len(exterior_ids), '/', len(capsule_faces), 'triangles; fitting constrained surface candidate', flush=True)
        surface_matrix, surface_offset, registration = fit_surface_affine(capsule, exterior_faces, kidney, kidney_faces, pelvis[patch], ureter_surface, audit.Surface)
        candidates.append(('surface_affine_constrained', surface_matrix, surface_offset))
        fits = []
        for name, matrix, offset in candidates:
            moved_capsule = transform(capsule, matrix, offset)
            cap_to_target = target_surface.distances(moved_capsule)
            target_to_cap = audit.Surface(moved_capsule, capsule_faces).distances(kidney)
            moved_query = transform(query, matrix, offset)
            votes = containment_votes(moved_query, kidney, kidney_faces)
            # Neither the thin source capsule material nor an open target kidney
            # establishes an organ enclosure. Do not convert parity into leakage.
            sample_surface_distance = target_surface.distances(moved_query)
            interface = ureter_surface.distances(transform(pelvis[patch], matrix, offset)) if patch.any() else np.array([])
            singular = np.linalg.svd(matrix, compute_uv=False)
            determinant = float(np.linalg.det(matrix))
            assert determinant > 0
            arrays = {}
            for key, _, current_side in entries:
                if current_side != side: continue
                vertices, faces, normals = meshes[key]
                transformed = transform(vertices, matrix, offset)
                transformed_normals = normals @ np.linalg.inv(matrix)
                transformed_normals /= np.maximum(np.linalg.norm(transformed_normals, axis=1, keepdims=True), 1e-20)
                arrays[key+'__positions'] = transformed.astype('<f4')
                arrays[key+'__indices'] = faces.astype('<u4')
                arrays[key+'__normals'] = np.rint(transformed_normals*32767).astype('<i2')
                assert np.isfinite(transformed).all() and np.isfinite(transformed_normals).all()
            arrays['diagnostic_outer_patch__indices'] = exterior_faces.astype('<u4')
            arrays['reference_capsule__positions'] = moved_capsule.astype('<f4'); arrays['reference_capsule__indices'] = capsule_faces.astype('<u4')
            arrays['target_kidney__positions'] = kidney.astype('<f4'); arrays['target_kidney__indices'] = kidney_faces.astype('<u4')
            arrays['retained_ureter__positions'] = ureter.astype('<f4'); arrays['retained_ureter__indices'] = ureter_faces.astype('<u4')
            artifact = output/f'{side}-{name}.npz'; write_npz(artifact, arrays)
            fits.append(dict(name=name, matrix=matrix.tolist(), offsetM=offset.tolist(), singularValues=singular.tolist(), maximumToMinimumScaleRatio=float(singular.max()/singular.min()), volumeScale=determinant, capsuleToTargetSurface=quantiles(cap_to_target), targetToCapsuleSurface=quantiles(target_to_cap), symmetricSampledP95Mm=max(quantiles(cap_to_target)['p95Mm'], quantiles(target_to_cap)['p95Mm']), sourceAssemblyInvariants=dict(sharedAffineTransformForEveryInternalMesh=True, sourceTriangleIndicesPreserved=True, positiveConstantJacobian=determinant, coincidentSourcePointsRemainCoincident=True, uniformPairwiseDistanceScale=float(singular.mean()) if np.allclose(singular, singular[0], atol=1e-10) else None, anglesPreserved=bool(np.allclose(singular, singular[0], atol=1e-10))), containmentScreen=dict(validAsOrganEnvelopeContainment=False, reason='Source capsule is thin closed tissue material; outer cortex is also segmented tissue material. Target kidney has open boundary edges. No verified organ-enclosure surface is available.', sampleCount=len(query), targetParityVotesDiagnostic=dict(threeInside=int((votes == 3).sum()), threeOutside=int((votes == 0).sum()), inconsistent=int(((votes != 0) & (votes != 3)).sum())), unsignedInternalSampleToTargetSurface=quantiles(sample_surface_distance)), ureterInterfaceScreen=dict(sourcePelvisToSourceUreterThresholdMm=2, sourcePelvisVertexIndices=patch_indices.tolist(), sourcePatchDistance=quantiles(distances[patch]), transformedPatchToRetainedUreter=quantiles(interface)), artifact=dict(path=artifact.name, sha256=audit.sha(artifact))))
            fits[-1]['registration'] = registration if name == 'surface_affine_constrained' else None
            fits[-1]['outerPatchToTargetSurface'] = quantiles(target_surface.distances(moved_capsule[np.unique(exterior_faces)]))
            fits[-1]['targetToOuterPatchSurface'] = quantiles(audit.Surface(moved_capsule, exterior_faces).distances(kidney))
            fits[-1]['transformedHilumToRetainedUreter'] = quantiles(ureter_surface.distances(transform(hilum, matrix, offset)))
            print(side, name, 'capsule p95 mm', round(fits[-1]['symmetricSampledP95Mm'], 3), 'volume scale', round(determinant, 3), 'ureter p95', None if not len(interface) else round(float(np.quantile(interface, .95))*1000, 3), flush=True)
        side_results.append(dict(side=side, capsuleId=capsule_id, targetKidneyId=target_id, retainedUreterId=ureter_id, internalParts=sum(1 for _, _, s in entries if s == side), sourceCapsuleActualExtentsMm=(np.ptp(capsule, axis=0)*1000).tolist(), targetKidneyActualExtentsMm=(np.ptp(kidney, axis=0)*1000).tolist(), sourceCapsuleTopology=topology(capsule, capsule_faces), diagnosticOuterPatch=dict(method='Outward-facing capsule triangles whose centroid is the farthest ray hit from the bounds center. No triangles are filled or invented; concavities and hilum may be omitted.', sourceTriangleIndices=exterior_ids.tolist(), topology=topology(capsule, exterior_faces), validAsOrganEnclosure=False), sourceOuterCortexTopology=topology(cortex, cortex_faces), targetKidneyTopology=topology(kidney, kidney_faces), sourceHilumBoundsM=[hilum.min(0).tolist(), hilum.max(0).tolist()], sourceMaterialParityDiagnostic=dict(validAsOrganEnvelopeContainment=False, sampleCount=len(query), capsuleThreeInside=int((source_votes == 3).sum()), capsuleThreeOutside=int((source_votes == 0).sum()), capsuleInconsistentVotes=int(((source_votes != 0) & (source_votes != 3)).sum()), outerCortexThreeInside=int((cortex_votes == 3).sum()), outerCortexThreeOutside=int((cortex_votes == 0).sum()), outerCortexInconsistentVotes=int(((cortex_votes != 0) & (cortex_votes != 3)).sum()), sampleVertices=query_ids), fits=fits))
    assert audit.sha(target.path) == original_target_sha
    source_chunks = {source.data['chunks'][i]['url']: hashlib.sha256(b).hexdigest() for i, b in source.buffers.items()}
    target_chunks = {target.data['chunks'][i]['url']: hashlib.sha256(b).hexdigest() for i, b in target.buffers.items()}
    for atlas in [source, target]:
        for i, original in atlas.buffers.items():
            path = ROOT/'public'/atlas.data['chunks'][i]['url'].lstrip('/')
            assert hashlib.sha256(path.read_bytes()).digest() == hashlib.sha256(original).digest()
    report = dict(schemaVersion=1, enabled=False, status='Review-only candidate; kidney registration and source assembly compatibility remain unvalidated', source='HRA united-female v1.5', license='CC BY 4.0', inventoryCount=len(entries), inventorySha256=inventory_hash.hexdigest(), inventory=source_inventory, inputHashes=dict(sourceManifestSha256=SOURCE_SHA, targetManifestSha256=original_target_sha, sourceChunks=source_chunks, targetChunks=target_chunks), canonicalAssetsUnchanged=True, method='Translate or fit each whole source-side assembly using the actual capsule, never independently squeeze individual internals. Similarity ICP uses deterministic vertex correspondences only to propose a fit; final bidirectional surface distances use exact nearest target triangles at every capsule/target vertex. A bounds-derived anisotropic candidate is diagnostic only. A second bounded affine proposal fits exact triangle surfaces bidirectionally over a radial exterior patch, with PCA seeds and a source ureter-proximity constraint; it remains disabled.', sides=side_results, limitations=['No reviewed renal hilum, pole, artery/vein or ureteropelvic-junction landmarks are available.', 'A capsule fit is not evidence that the separately sourced renal pelvis joins the existing ureter or vasculature.', 'No organ-envelope containment result is valid: the source capsule/outer cortex represent tissue material, and target kidneys have open boundaries. Parity diagnostics must not be interpreted as displaced source internals or containment approval.', 'The containment sample uses up to 24 deterministic vertices per internal mesh, not area weights or exhaustive collision checking.', 'The 2 mm source ureter-proximity cutoff is an engineering patch definition, not an anatomical endpoint annotation.', 'Affine bounds fitting can compress anatomy substantially and is not promoted merely because its envelope matches.', 'All output geometry is a temporary review artifact; no current viewer mesh is replaced or added.'])
    (output/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--output-dir', type=Path, default=Path('/tmp/female-kidney-candidate'))
    args = parser.parse_args(); run(args.output_dir)
