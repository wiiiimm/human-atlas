"""Contour the original HRA breast assembly in its own source frame.

Preserve the source footprint and topology, soften anterior projection, lift the
lower anterior tissue, and translate the entire assembly into the final atlas.
This is a visual contour/placement, not an anatomical registration or approval.
"""
import copy
import numpy as np

BREAST_IDS = tuple(f'VH_F_{name}_{side}' for side in ('L', 'R') for name in (
    'fat', 'areolar_tubercles', 'nipple', 'areola', 'mammary_lobes',
    'main_lactiferous_ducts', 'main_lactiferous_sinuses', 'suspensory_ligaments'))
FRONT_DEPTH_RAMP = (.35, .9)
LOWER_HEIGHT_RAMP = (0., .15, .45, .7)
NORMAL_DIFFERENCE_STEP_M = 1e-5


def smoothstep(a, b, values):
    t = np.clip((values-a)/(b-a), 0., 1.)
    return t*t*(3.-2.*t)


def describe(source_parts):
    transforms = {}
    for side in ('L', 'R'):
        bounds = copy.deepcopy(source_parts[f'VH_F_fat_{side}']['bounds'])
        transforms[side] = dict(scale=1., center=((np.array(bounds[0])+bounds[1])/2).tolist(),
                                target=[.082 if side == 'L' else -.082, 1.21, .088],
                                sourceBounds=bounds, projectionReduction=.15, lowerLiftM=.01)
    return dict(schemaVersion=1, source='HRA united-female v1.5', license='CC BY 4.0', ids=list(BREAST_IDS),
                transforms=transforms, frontDepthRamp=list(FRONT_DEPTH_RAMP),
                lowerHeightRamp=list(LOWER_HEIGHT_RAMP), normalDifferenceStepM=NORMAL_DIFFERENCE_STEP_M,
                bodyMorphApplied=False,
                method='Original HRA breast meshes: retain source base footprint, reduce anterior projection, lift lower anterior tissue, then translate into final female coordinates.',
                normalMethod='Inverse-transpose Jacobian of the contour field, central finite differences in source meters; uniform scale and translation follow.',
                reviewStatus='Illustrative source-derived contour and experimental placement; no attachment, containment, or anatomical approval.')


def contour_points(points, transform):
    points = np.asarray(points, dtype=np.float64)
    lo, hi = np.array(transform['sourceBounds'], dtype=np.float64)
    if np.any(hi <= lo):
        raise ValueError('Breast source bounds must have positive extent')
    depth = (points[:, 2]-lo[2])/(hi[2]-lo[2])
    height = (points[:, 1]-lo[1])/(hi[1]-lo[1])
    front = smoothstep(*FRONT_DEPTH_RAMP, depth)
    lower = smoothstep(*LOWER_HEIGHT_RAMP[:2], height)*(1.-smoothstep(*LOWER_HEIGHT_RAMP[2:], height))
    result = points.copy()
    result[:, 1] += transform['lowerLiftM']*front*lower
    result[:, 2] -= transform['projectionReduction']*(points[:, 2]-lo[2])*front
    return result


def contour_normals(points, normals, transform):
    """Use the same triangular Jacobian structure as the contour, preserving x."""
    h = NORMAL_DIFFERENCE_STEP_M
    yp, ym, zp, zm = [np.array(points, dtype=np.float64, copy=True) for _ in range(4)]
    yp[:, 1] += h; ym[:, 1] -= h; zp[:, 2] += h; zm[:, 2] -= h
    plus_y, minus_y = contour_points(yp, transform), contour_points(ym, transform)
    plus_z, minus_z = contour_points(zp, transform), contour_points(zm, transform)
    a = (plus_y[:, 1]-minus_y[:, 1])/(2*h)
    b = (plus_z[:, 1]-minus_z[:, 1])/(2*h)
    c = (plus_z[:, 2]-minus_z[:, 2])/(2*h)
    determinant = a*c
    if not np.isfinite(determinant).all() or np.any(determinant <= 0):
        raise ValueError('Breast contour has an invalid Jacobian at a source vertex')
    result = np.array(normals, dtype=np.float64, copy=True)
    result[:, 1] = normals[:, 1]/a
    result[:, 2] = (normals[:, 2]-b*normals[:, 1]/a)/c
    lengths = np.linalg.norm(result, axis=1, keepdims=True)
    if not np.isfinite(lengths).all() or np.any(lengths <= 0):
        raise ValueError('Breast source has an invalid normal')
    return result/lengths, float(determinant.min())


def apply(points, normals, transform):
    if not np.isfinite(transform['scale']) or transform['scale'] <= 0:
        raise ValueError('Breast uniform scale must be positive')
    mapped = (contour_points(points, transform)-np.array(transform['center']))*transform['scale']+np.array(transform['target'])
    normal, minimum = contour_normals(points, normals, transform)
    if not np.isfinite(mapped).all():
        raise ValueError('Breast contour contains non-finite positions')
    # Match the approved JavaScript trial's Math.round convention, including ties.
    encoded = np.floor(normal*32767.+.5).astype('<i2')
    return mapped, encoded, minimum
