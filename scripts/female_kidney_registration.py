"""Deterministic geometric proposals, not anatomically annotated registration."""
import numpy as np


def radial_outer_patch(vertices, faces):
    """Keep outward-facing triangles visible from infinity along their center ray.

    This extracts a patch, not an enclosure: concavities and the hilum may be
    excluded. No faces are added and no internal component is called an organ.
    """
    center = (vertices.min(0)+vertices.max(0))/2
    tri = vertices[faces]; a = tri[:, 0]; e1 = tri[:, 1]-a; e2 = tri[:, 2]-a
    normal = np.cross(e1, e2); mid = tri.mean(1)
    candidates = np.flatnonzero(np.einsum('ij,ij->i', normal, mid-center) > 0)
    selected = []
    for index in candidates:
        direction = mid[index]-center; distance_to_face = np.linalg.norm(direction)
        direction /= distance_to_face
        h = np.cross(direction, e2); det = np.einsum('ij,ij->i', e1, h)
        valid = abs(det) > 1e-14
        inv = np.divide(1., det, out=np.zeros_like(det), where=valid)
        s = center-a; u = inv*np.einsum('ij,ij->i', s, h)
        q = np.cross(s, e1); v = inv*(q @ direction)
        distance = inv*np.einsum('ij,ij->i', e2, q)
        hits = distance[valid & (u >= -1e-9) & (v >= -1e-9) & (u+v <= 1+1e-9) & (distance > 0)]
        if len(hits) and hits.max() <= distance_to_face+1e-7:
            selected.append(index)
    if not selected: raise ValueError('No radial exterior patch found')
    return np.asarray(selected, dtype=int)


def principal_frame(vertices):
    center = vertices.mean(0)
    _, _, vt = np.linalg.svd(vertices-center, full_matrices=False)
    frame = vt.T
    if np.linalg.det(frame) < 0: frame[:, -1] *= -1
    return center, frame


def constrain_matrix(matrix, minimum_scale=.55, maximum_scale=1.10, minimum_volume=.35):
    """Engineering distortion bounds, not accepted anatomical scaling limits."""
    u, singular, vt = np.linalg.svd(matrix)
    if np.linalg.det(u @ vt) < 0: u[:, -1] *= -1
    singular = np.clip(singular, minimum_scale, maximum_scale)
    if minimum_volume > maximum_scale**3: raise ValueError('Infeasible volume bound')
    for _ in range(3):
        if np.prod(singular) >= minimum_volume-1e-12: break
        free = singular < maximum_scale
        singular[free] *= (minimum_volume/np.prod(singular))**(1/int(free.sum()))
        singular = np.minimum(maximum_scale, singular)
    return (u*singular) @ vt


def weighted_affine(source, target, weights, prior, regularization=.035):
    """Solve centered weighted point pairs with a scale-normalized linear prior."""
    weights = weights/weights.sum()
    center = np.sum(source*weights[:, None], axis=0)
    target_center = np.sum(target*weights[:, None], axis=0)
    x = source-center; y = target-target_center
    covariance = x.T @ (x*weights[:, None])
    strength = np.trace(covariance)*regularization/3
    matrix = np.linalg.solve(covariance+np.eye(3)*strength,
                             x.T @ (y*weights[:, None])+strength*prior.T).T
    matrix = constrain_matrix(matrix)
    return matrix, target_center-matrix @ center


def fit_surface_affine(source, source_faces, target, target_faces, patch, ureter_surface, surface_type, maximum_samples=192, iterations=12):
    """Bidirectional surface ICP with one shared, bounded affine transform.

    PCA supplies geometric long-axis seeds. Source interface points are projected
    onto actual retained ureter triangles; this remains a proximity constraint.
    """
    used = np.unique(source_faces)
    source_sample = source[used[np.linspace(0, len(used)-1, min(maximum_samples, len(used)), dtype=int)]]
    target_sample = target[np.linspace(0, len(target)-1, min(maximum_samples, len(target)), dtype=int)]
    source_center, source_frame = principal_frame(source[used])
    target_center, target_frame = principal_frame(target)
    target_surface = surface_type(target, target_faces)
    scale = np.linalg.norm(np.ptp(target, axis=0))/np.linalg.norm(np.ptp(source, axis=0))
    starts = [('bounds', np.diag(np.ptp(target, axis=0)/np.ptp(source, axis=0)))]
    for index, signs in enumerate([(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]):
        starts.append((f'principal_axes_{index}', scale*target_frame @ np.diag(signs) @ source_frame.T))
    records = []
    for name, initial in starts:
        matrix = constrain_matrix(initial); offset = target_center-matrix @ source_center
        history = []
        for _ in range(iterations):
            moved = source_sample @ matrix.T+offset
            forward = np.array([target_surface.nearest(p)[1] for p in moved])
            moved_surface = surface_type(source @ matrix.T+offset, source_faces)
            reverse = np.array([moved_surface.nearest(p)[1] for p in target_sample])
            reverse_source = (reverse-offset) @ np.linalg.inv(matrix).T
            groups_source = [source_sample, reverse_source]
            groups_target = [forward, target_sample]
            groups_weights = [np.full(len(source_sample), .5/len(source_sample)), np.full(len(target_sample), .5/len(target_sample))]
            if len(patch):
                moved_patch = patch @ matrix.T+offset
                paired_patch = np.array([ureter_surface.nearest(p)[1] for p in moved_patch])
                groups_source.append(patch); groups_target.append(paired_patch)
                groups_weights.append(np.full(len(patch), .25/len(patch)))
            matrix, offset = weighted_affine(np.concatenate(groups_source), np.concatenate(groups_target), np.concatenate(groups_weights), matrix)
            residual = np.linalg.norm(np.concatenate(groups_source) @ matrix.T+offset-np.concatenate(groups_target), axis=1)
            history.append(float(np.sqrt(np.average(residual**2, weights=np.concatenate(groups_weights))))*1000)
        # Score proposals on surface samples, with interface cost to reject flips
        # which match the exterior while moving the source hilum away.
        fdist = target_surface.distances(source_sample @ matrix.T+offset)
        rdist = surface_type(source @ matrix.T+offset, source_faces).distances(target_sample)
        interface = ureter_surface.distances(patch @ matrix.T+offset) if len(patch) else np.array([0.])
        score = float(np.sqrt((np.mean(fdist**2)+np.mean(rdist**2))/2+.25*np.mean(interface**2)))
        records.append(dict(seed=name, matrix=matrix, offset=offset, score=score, iterationRmsMm=history))
    best = min(records, key=lambda row: row['score'])
    u, _, vt = np.linalg.svd(best['matrix'])
    rotation_degrees = float(np.degrees(np.arccos(np.clip((np.trace(u @ vt)-1)/2, -1, 1))))
    moved_axis = best['matrix'] @ source_frame[:, 0]
    moved_axis /= np.linalg.norm(moved_axis)
    axis_degrees = float(np.degrees(np.arccos(np.clip(abs(moved_axis @ target_frame[:, 0]), 0, 1))))
    return best['matrix'], best['offset'], dict(method='Bidirectional nearest-triangle affine ICP; PCA orientation seeds and source ureter proximity constraint; one shared transform per side', engineeringBounds=dict(minimumSingularScale=.55, maximumSingularScale=1.10, minimumVolumeScale=.35), interfaceRelativeWeight=.25, samplesPerDirection=maximum_samples, iterationsPerSeed=iterations, selectedSeed=best['seed'], polarRotationDegrees=rotation_degrees, geometricLongAxisAngleDegrees=axis_degrees, seeds=[dict(seed=row['seed'], selectionScoreMm=row['score']*1000, iterationRmsMm=row['iterationRmsMm']) for row in records])
