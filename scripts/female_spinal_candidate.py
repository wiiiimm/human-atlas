#!/usr/bin/env python3
"""Build an isolated HRA spinal-cord registration experiment (numpy only).

python3 scripts/female_spinal_candidate.py --output-dir /tmp/human-atlas-spinal-candidate
Never writes production models. An enabled:false report is a measured research
result, not approval to integrate the candidate into an anatomy viewer.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ORDINALS = 'first second third fourth fifth sixth seventh eigth ninth tenth eleventh twelfth'.split()
CORD_IDS = ([f'VH_F_C{i}_segment_of_cervical_spinal_cord' for i in range(1, 9)] +
            [f'VH_F_{word}_thoracic_spinal_cord_segment' for word in ORDINALS] +
            [f'VH_F_{word}_lumbar_spinal_cord_segment' for word in ORDINALS[:5]] +
            [f'VH_F_{word}_sacral_spinal_cord_segment' for word in ORDINALS[:4]])
TARGET_IDS = (['FJ3176', 'FJ3177', 'FJ3161', 'FJ3164', 'FJ3167', 'FJ3170', 'FJ3172'] +
              ['FJ3158', 'FJ3160', 'FJ3163', 'FJ3166', 'FJ3169', 'FJ3171', 'FJ3173',
               'FJ3174', 'FJ3175', 'FJ3154', 'FJ3155', 'FJ3156', 'FJ3157'])
BONE_PAIRS = list(zip([f'VH_F_cervical_vertebra_{i}' for i in range(1, 8)] +
                     [f'VH_F_thoracic_vertebra_{i}' for i in range(1, 13)] +
                     ['VH_F_lumbar_vertebra_1'], TARGET_IDS))


def sha(data):
    return hashlib.sha256(data).hexdigest()


class Atlas:
    def __init__(self, models, filename):
        self.path = models / filename
        self.raw = self.path.read_bytes()
        self.manifest = json.loads(self.raw)
        self.parts = {p['id']: p for p in self.manifest['parts']}
        self.models = models
        self.buffers = {}

    def mesh(self, key):
        p = self.parts[key]
        url = self.manifest['chunks'][p['chunk']]['url']
        if url not in self.buffers:
            self.buffers[url] = (self.models / Path(url).name).read_bytes()
        b = self.buffers[url]
        v = np.frombuffer(b, '<f4', p['vertexCount'] * 3, p['positions']).reshape(-1, 3).astype(float)
        n = np.frombuffer(b, '<i2', p['vertexCount'] * 3, p['normals']).reshape(-1, 3).astype(float) / 32767
        f = np.frombuffer(b, '<u4', p['indexCount'], p['indices']).reshape(-1, 3).copy()
        return v, n, f

    def evidence(self):
        return dict(manifest=self.path.name, manifestSha256=sha(self.raw),
                    chunks={k: sha(v) for k, v in sorted(self.buffers.items())})


def centroid(v, f):
    t = v[f]
    area2 = np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1)
    return np.average(t.mean(axis=1), axis=0, weights=area2)


def topology(v, f):
    unique, inverse = np.unique(v, axis=0, return_inverse=True)
    welded = inverse[f]
    edges = np.sort(np.concatenate([welded[:, [0, 1]], welded[:, [1, 2]], welded[:, [2, 0]]]), axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    area2 = np.linalg.norm(np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]]), axis=1)
    return dict(vertices=len(v), weldedVertices=len(unique), triangles=len(f),
                boundaryEdges=int(np.sum(counts == 1)), nonmanifoldEdges=int(np.sum(counts > 2)),
                zeroAreaTriangles=int(np.sum(area2 < 1e-14)))


class Registration:
    """One continuous field for all segments; never match cord labels to bone levels."""
    def __init__(self, source, target):
        self.anchors = []
        for sid, tid in BONE_PAIRS:
            sv, _, sf = source.mesh(sid)
            tv, _, tf = target.mesh(tid)
            self.anchors.append(dict(sourceId=sid, targetId=tid,
                                     sourceCentroid=centroid(sv, sf).tolist(),
                                     targetCentroid=centroid(tv, tf).tolist()))
        self.anchors.sort(key=lambda a: a['sourceCentroid'][1])
        self.s = np.array([a['sourceCentroid'] for a in self.anchors])
        self.t = np.array([a['targetCentroid'] for a in self.anchors])
        self.slopes = np.diff(self.t - self.s, axis=0) / np.diff(self.s[:, 1])[:, None]
        if not (np.diff(self.s[:, 1]) > 0).all() or not (np.diff(self.t[:, 1]) > 0).all():
            raise ValueError('Corresponding vertebral heights are not ordered')

    def map(self, v, n=None):
        i = np.clip(np.searchsorted(self.s[:, 1], v[:, 1], side='right') - 1, 0, len(self.s) - 2)
        slope = self.slopes[i]
        mapped = v + self.t[i] - self.s[i] + (v[:, 1] - self.s[i, 1])[:, None] * slope
        if n is None:
            return mapped
        normal = n.copy()
        normal[:, 1] = (n[:, 1] - slope[:, 0] * n[:, 0] - slope[:, 2] * n[:, 2]) / (1 + slope[:, 1])
        normal /= np.linalg.norm(normal, axis=1)[:, None]
        return mapped, normal


def surface_helper():
    spec = importlib.util.spec_from_file_location('pelvis_surface_audit', ROOT / 'scripts/pelvis-surface-audit.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.Surface


def continuity(meshes, Surface):
    result = []
    for upper, lower in zip(CORD_IDS, CORD_IDS[1:]):
        uv, _, uf = meshes[upper]
        lv, _, lf = meshes[lower]
        # Deterministic full end strips (2 mm) evaluated against actual triangles.
        up = uv[uv[:, 1] <= uv[:, 1].min() + .002]
        lp = lv[lv[:, 1] >= lv[:, 1].max() - .002]
        # Fixed 32 evenly indexed vertices per strip bounds runtime, not an exhaustive gap proof.
        up = up[np.linspace(0, len(up) - 1, min(32, len(up)), dtype=int)]
        lp = lp[np.linspace(0, len(lp) - 1, min(32, len(lp)), dtype=int)]
        distances = np.r_[Surface(lv, lf).distances(up), Surface(uv, uf).distances(lp)]
        result.append(dict(upper=upper, lower=lower, samples=len(distances),
                           axialBoundsGapMm=float((uv[:, 1].min() - lv[:, 1].max()) * 1000),
                           sampledEndToSurfaceMinMm=float(distances.min() * 1000),
                           sampledEndToSurfaceMedianMm=float(np.median(distances) * 1000),
                           sampledEndToSurfaceMaxMm=float(distances.max() * 1000)))
    return result


def section(triangles, height):
    """Triangle/plane intersections in the atlas x/z plane; skip coplanar triangles."""
    low, high = triangles[:, :, 1].min(axis=1), triangles[:, :, 1].max(axis=1)
    tri = triangles[(low < height) & (high > height)]
    segments = []
    for t in tri:
        points = []
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            if (a[1] <= height < b[1]) or (b[1] <= height < a[1]):
                q = a + (height - a[1]) / (b[1] - a[1]) * (b - a)
                points.append(q[[0, 2]])
        if len(points) == 2 and np.linalg.norm(points[0] - points[1]) > 1e-10:
            segments.append(points)
    return np.array(segments).reshape(-1, 2, 2)


def ray_hits(segments, center, direction):
    if not len(segments):
        return np.array([])
    a = segments[:, 0] - center
    edge = segments[:, 1] - segments[:, 0]
    cross = lambda p, q: p[..., 0] * q[..., 1] - p[..., 1] * q[..., 0]
    den = cross(direction, edge)
    valid = abs(den) > 1e-12
    r = np.zeros(len(segments)); u = r.copy()
    r[valid] = cross(a[valid], edge[valid]) / den[valid]
    u[valid] = cross(a[valid], direction) / den[valid]
    hits = np.sort(r[valid & (r > 1e-8) & (u >= 0) & (u < 1)])
    return hits[np.r_[True, np.diff(hits) > 1e-7]] if len(hits) else hits


def radial_boundaries(segments, center, directions):
    """Vectorized first boundary and parity; endpoint half-open rule avoids doubles."""
    if not len(segments):
        return np.full(len(directions), np.inf), np.zeros(len(directions), dtype=int)
    a = segments[:, 0] - center
    edge = segments[:, 1] - segments[:, 0]
    cross = lambda p, q: p[..., 0] * q[..., 1] - p[..., 1] * q[..., 0]
    denominator = cross(directions[:, None], edge[None])
    valid = abs(denominator) > 1e-12
    denominator = np.where(valid, denominator, 1)
    ray = cross(a, edge)[None] / denominator
    fraction = cross(a[None], directions[:, None]) / denominator
    hits = np.where(valid & (ray > 1e-8) & (fraction >= 0) & (fraction < 1), ray, np.inf)
    # Match the measured screen's 0.1 micron hit deduplication, including seams.
    hits.sort(axis=1)
    finite = np.isfinite(hits)
    differences = np.zeros_like(hits[:, 1:])
    np.subtract(hits[:, 1:], hits[:, :-1], out=differences, where=finite[:, 1:])
    distinct = finite.copy()
    distinct[:, 1:] &= differences > 1e-7
    return hits[:, 0], distinct.sum(axis=1) % 2


class CanalCorrection:
    """Compact C2 transverse translation fitted to actual cervical bone sections."""
    FIT_FRACTIONS = (.1, .3, .5, .7, .9)
    HOLDOUT_FRACTIONS = (.025, .2, .4, .6, .8, .975)

    def __init__(self, cords, bones):
        self.radius = .035
        mids = {key: float((v[:, 1].min() + v[:, 1].max()) / 2) for key, (v, _, _) in cords.items()}
        self.centers = np.linspace(mids[CORD_IDS[7]], mids[CORD_IDS[2]], 5)
        low, high = self.centers[0] - self.radius, self.centers[-1] + self.radius
        self.ids = [key for key in CORD_IDS if cords[key][0][:, 1].max() >= low and cords[key][0][:, 1].min() <= high]
        triangles = np.concatenate([v[f] for v, _, f in bones.values()])
        theta = np.linspace(.013, 2 * np.pi + .013, 64, endpoint=False)
        self.directions = np.column_stack([np.cos(theta), np.sin(theta)])
        self.sections = []
        for key in self.ids:
            v, _, f = cords[key]
            for fraction in self.FIT_FRACTIONS:
                height = float(v[:, 1].min() + fraction * np.ptp(v[:, 1]))
                cord = section(v[f], height)
                bone = section(triangles, height)
                center = (cord.min(axis=(0, 1)) + cord.max(axis=(0, 1))) / 2
                hits = [ray_hits(cord, center, direction) for direction in self.directions]
                radii = np.array([hit[-1] if len(hit) else np.inf for hit in hits])
                valid = np.isfinite(radii)
                first, _ = radial_boundaries(bone, center, self.directions)
                self.sections.append(dict(id=key, fraction=fraction, height=height, bone=bone,
                                          center=center, radii=radii, valid=valid,
                                          initialCoverage=int(np.sum(np.isfinite(first) & valid))))
        self.matrix = self.basis(np.array([s['height'] for s in self.sections]))
        self.coefficients = np.zeros((5, 2))
        self.iterations = 0
        self.initial_loss = self.loss(self.coefficients)
        # Deterministic coordinate descent directly over the smooth field, not
        # independent per-segment translations or fitted annotations.
        for step in (.002, .001, .0005, .00025, .0001):
            for _ in range(20):
                changed = False
                for i in range(5):
                    for axis in range(2):
                        best = self.coefficients
                        value = self.loss(best)
                        for delta in (-step, step):
                            trial = self.coefficients.copy()
                            trial[i, axis] += delta
                            candidate_loss = self.loss(trial)
                            if candidate_loss < value:
                                best, value = trial, candidate_loss
                        if best is not self.coefficients:
                            self.coefficients = best
                            changed = True
                self.iterations += 1
                if not changed:
                    break
        self.final_loss = self.loss(self.coefficients)

    def basis(self, heights, derivative=False):
        delta = np.asarray(heights)[:, None] - self.centers
        r = abs(delta) / self.radius
        inside = np.maximum(1 - r, 0)
        return -20 * delta / self.radius**2 * inside**3 if derivative else inside**4 * (4 * r + 1)

    def loss(self, coefficients):
        shifts = self.matrix @ coefficients
        # Bounds prevent solving interference by moving outside the local canal.
        if np.linalg.norm(shifts, axis=1).max() > .006:
            return 1e12
        loss = .001 * float(np.sum((coefficients * 1000)**2))
        for shift, s in zip(shifts, self.sections):
            first, parity = radial_boundaries(s['bone'], s['center'] + shift, self.directions)
            clearances = (first[s['valid']] - s['radii'][s['valid']]) * 1000
            coverage = int(np.sum(np.isfinite(first) & s['valid']))
            loss += float(np.sum(np.maximum(.5 - clearances, 0)**2))
            loss += 100 * int(parity.sum()) + max(0, s['initialCoverage'] - coverage - 4)**2
        return loss

    def map(self, points, normals=None):
        shift = self.basis(points[:, 1]) @ self.coefficients
        mapped = points.copy()
        mapped[:, [0, 2]] += shift
        if normals is None:
            return mapped
        slope = self.basis(points[:, 1], derivative=True) @ self.coefficients
        result = normals.copy()
        result[:, 1] -= slope[:, 0] * normals[:, 0] + slope[:, 1] * normals[:, 2]
        result /= np.linalg.norm(result, axis=1)[:, None]
        inactive = (shift == 0).all(axis=1) & (slope == 0).all(axis=1)
        result[inactive] = normals[inactive]
        return mapped, result

    def evidence(self):
        heights = np.linspace(self.centers[0] - self.radius, self.centers[-1] + self.radius, 2001)
        return dict(method='Five compact C2 Wendland functions translate the cord in x/depth as a continuous function of height. Coefficients minimize actual bone/cord section clearance violations; no source triangles or transverse sizes are replaced.',
                    basisCentersM=self.centers.tolist(), basisRadiusM=self.radius,
                    coefficientsXZ_M=self.coefficients.tolist(), fitIds=self.ids,
                    fitFractions=list(self.FIT_FRACTIONS), holdoutFractions=list(self.HOLDOUT_FRACTIONS),
                    fitSectionCount=len(self.sections), iterations=self.iterations,
                    objective=dict(targetRadialClearanceMm=.5, maximumFitSectionShiftMm=6,
                                   coefficientPenalty=.001, oddBoneParityPenalty=100,
                                   allowedLostCoveredRays=4, initialLoss=self.initial_loss, finalLoss=self.final_loss),
                    maximumSampledShiftMm=float(np.linalg.norm(self.basis(heights) @ self.coefficients, axis=1).max() * 1000),
                    maximumSampledSlope=float(np.linalg.norm(self.basis(heights, derivative=True) @ self.coefficients, axis=1).max()),
                    jacobianDeterminant=1,
                    normalMethod='Analytical inverse-transpose of the transverse-translation Jacobian composed with the initial registration normal transform.',
                    limitation='A deterministic local optimum on sampled open cross-sections is not a complete canal containment proof.')


def section_comparison_svg(path, baseline, candidate, bones):
    triangles = np.concatenate([v[f] for v, _, f in bones.values()])
    lines = ['<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="700" viewBox="0 0 1100 700">',
             '<rect width="1100" height="700" fill="#101826"/>',
             '<text x="30" y="35" fill="white" font-size="20">Actual cervical cross-sections: before and after canal-derived correction</text>']
    for column, key in enumerate(CORD_IDS[4:6]):
        v, _, f = baseline[key]
        height = (v[:, 1].min() + v[:, 1].max()) / 2
        bone = section(triangles, height)
        old = section(v[f], height)
        cv, _, cf = candidate[key]
        new = section(cv[cf], height)
        center = (old.min(axis=(0, 1)) + old.max(axis=(0, 1))) / 2
        x0 = 275 + 550 * column
        lines.append(f'<text x="{x0-220}" y="80" fill="white" font-size="18">C{column+5} cord segment, height {height:.5f} m</text>')
        for segments, color, width in [(bone, '#c0cddd', 1.5), (old, '#ff8276', 2.2), (new, '#65e6c2', 2.2)]:
            projected = (segments - center) * np.array([6500, -6500]) + [x0, 370]
            d = ' '.join(f'M{a[0]:.2f},{a[1]:.2f} L{b[0]:.2f},{b[1]:.2f}' for a, b in projected)
            lines.append(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}"/>')
    lines += ['<text x="30" y="660" fill="white" font-size="16">Grey: vertebral triangles | Red: initial cord | Green: corrected cord | Open contours are retained.</text>', '</svg>']
    path.write_text('\n'.join(lines) + '\n')


def canal_screen(cords, bones, fractions=(.5,), ids=None):
    bone_triangles = np.concatenate([v[f] for v, _, f in bones.values()])
    result = []
    for key, fraction in [(key, fraction) for key in (ids or CORD_IDS) for fraction in fractions]:
        v, _, f = cords[key]
        height = float(v[:, 1].min() + fraction * np.ptp(v[:, 1]))
        cord = section(v[f], height)
        if not len(cord):
            raise ValueError(f'No cord cross-section for {key}')
        center = (cord.min(axis=(0, 1)) + cord.max(axis=(0, 1))) / 2
        bone = section(bone_triangles, height)
        clearances = []; escaping = 0; bone_parities = []; bone_escaping = 0; cord_escaping = 0
        for theta in np.linspace(.013, 2 * np.pi + .013, 64, endpoint=False):
            direction = np.array([np.cos(theta), np.sin(theta)])
            bh = ray_hits(bone, center, direction)
            ch = ray_hits(cord, center, direction)
            bone_parities.append(len(bh) % 2)
            bone_escaping += int(not len(bh))
            cord_escaping += int(not len(ch))
            if not len(bh) or not len(ch):
                escaping += 1
            else:
                clearances.append(float((bh[0] - ch[-1]) * 1000))
        result.append(dict(id=key, sectionFraction=fraction, heightM=height, centerXZ=center.tolist(),
                           boneCrossSectionSegments=len(bone), cordCrossSectionSegments=len(cord),
                           radialDirections=64, escapingDirections=escaping,
                           boneEscapingDirections=bone_escaping, cordEscapingDirections=cord_escaping,
                           oddBoneParityDirections=sum(bone_parities),
                           minimumRadialClearanceMm=min(clearances) if clearances else None,
                           medianRadialClearanceMm=float(np.median(clearances)) if clearances else None,
                           negativeRadialClearanceDirections=sum(c < 0 for c in clearances)))
    return result


def mesh_point_parity(point, triangles, direction):
    """Actual-triangle ray parity with ambiguous edge/tangent hits rejected."""
    edge1 = triangles[:, 1] - triangles[:, 0]
    edge2 = triangles[:, 2] - triangles[:, 0]
    h = np.cross(direction, edge2)
    det = np.sum(edge1 * h, axis=1)
    valid = abs(det) > 1e-12
    inv = np.zeros_like(det)
    inv[valid] = 1 / det[valid]
    delta = point - triangles[:, 0]
    u = inv * np.sum(delta * h, axis=1)
    q = np.cross(delta, edge1)
    v = inv * np.sum(direction * q, axis=1)
    distance = inv * np.sum(edge2 * q, axis=1)
    hit = valid & (u >= 0) & (v >= 0) & (u + v <= 1) & (distance > 1e-8)
    if np.any(hit & ((u < 1e-6) | (v < 1e-6) | (1-u-v < 1e-6))):
        return None
    distances = np.sort(distance[hit])
    return int(np.sum(np.r_[True, np.diff(distances) > 1e-7]) % 2) if len(distances) else 0


def medulla_interface(cords, medulla, Surface):
    points = cords[CORD_IDS[0]][0]
    top = points[points[:, 1] >= points[:, 1].max() - .001]
    directions = np.array([[1, .371, .137], [-.219, 1, .413], [.173, -.281, 1.]])
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    combined_v = np.concatenate([v for v, _, _ in medulla.values()])
    surfaces = [(Surface(v, f), v[f], topology(v, f)) for v, _, f in medulla.values()]
    distances = np.min(np.array([surface.distances(top) for surface, _, _ in surfaces]), axis=0)
    states = []
    for point in top:
        classifications = []
        for surface, triangles, stats in surfaces:
            if stats['boundaryEdges'] or stats['nonmanifoldEdges']:
                classifications.append('uncertain')
                continue
            votes = [mesh_point_parity(point, triangles, direction) for direction in directions]
            classifications.append('uncertain' if None in votes or len(set(votes)) != 1 else ('inside' if votes[0] else 'outside'))
        states.append('inside' if 'inside' in classifications else ('uncertain' if 'uncertain' in classifications else 'outside'))
    return dict(method='All first-segment vertices in its upper 1 mm strip measured against actual medulla triangles; three non-ambiguous ray parities per closed medulla mesh classify points in their union.',
                cordId=CORD_IDS[0], medullaIds=list(medulla), sampledUpperStripVertices=len(top),
                upperStripCounts={state: states.count(state) for state in ('inside', 'outside', 'uncertain')},
                minimumSurfaceDistanceMm=float(distances.min() * 1000),
                medianSurfaceDistanceMm=float(np.median(distances) * 1000),
                maximumSurfaceDistanceMm=float(distances.max() * 1000),
                cordBounds=[points.min(axis=0).tolist(), points.max(axis=0).tolist()],
                medullaBounds=[combined_v.min(axis=0).tolist(), combined_v.max(axis=0).tolist()],
                axialBoundsGapMm=float((combined_v[:, 1].min() - points[:, 1].max()) * 1000),
                interpretation='Inside samples indicate geometric overlap, not a welded anatomical junction; unsigned surface distances are not penetration depths. Positive axial bounds gap guarantees vertical separation; negative values alone only indicate overlapping height ranges.')


def svg_overlay(path, cords, bones):
    lines = ['<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="950" viewBox="0 0 1000 950">',
             '<rect width="1000" height="950" fill="#101826"/>',
             '<text x="30" y="30" fill="white" font-size="18">Spinal registration candidate — unvalidated, source topology retained</text>']
    allv = np.concatenate([v for v, _, _ in bones.values()])
    low, high = allv[:, 1].min() - .015, allv[:, 1].max() + .015
    scale = 830 / (high - low)
    for axis, x0, title in [(0, 250, 'Front: x / height'), (2, 730, 'Side: depth / height')]:
        center = (allv[:, axis].min() + allv[:, axis].max()) / 2
        lines.append(f'<text x="{x0-150}" y="65" fill="white" font-size="16">{title}</text>')
        for group, color, opacity in [(bones, '#8b9eb7', '.22'), (cords, '#ffd66e', '.48')]:
            for v, _, f in group.values():
                projected = np.column_stack([x0 + (v[:, axis] - center) * scale, 85 + (high - v[:, 1]) * scale])
                # All triangles, no invented outline or body envelope.
                d = ' '.join('M' + ' L'.join(f'{x:.2f},{y:.2f}' for x, y in triangle) + 'Z' for triangle in projected[f])
                lines.append(f'<path d="{d}" fill="none" stroke="{color}" stroke-width=".35" opacity="{opacity}"/>')
    lines += ['<text x="30" y="930" fill="white" font-size="13">Grey: current vertebrae and medulla; gold: candidate HRA cord. Sampled fitting does not prove complete containment.</text>', '</svg>']
    path.write_text('\n'.join(lines) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--models-dir', type=Path, default=ROOT / 'public/models')
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output == ROOT.resolve() or ROOT.resolve() in output.parents:
        raise ValueError('Candidate output must resolve outside the repository, including symlink aliases')
    output.mkdir(parents=True, exist_ok=True)
    for name in ('report.json', 'candidate.json', 'candidate.bin', 'overlay.svg', 'baseline-overlay.svg', 'cervical-sections.svg'):
        if (output / name).is_symlink():
            raise ValueError(f'Refusing symlinked candidate output: {output / name}')
    source = Atlas(args.models_dir, 'atlas-female.json')
    target = Atlas(args.models_dir, 'atlas-female-reconstructed.json')
    discovered = sorted(p['id'] for p in source.manifest['parts'] if p['id'].startswith('VH_F_') and 'spinal_cord' in p['id'])
    if discovered != sorted(CORD_IDS):
        raise ValueError('Source spinal-cord inventory changed; review the explicit allowlist')
    original = {key: source.mesh(key) for key in CORD_IDS}
    registration = Registration(source, target)
    fitted = {}
    for key, (v, n, f) in original.items():
        p, normal = registration.map(v, n)
        if not np.isfinite(p).all() or not np.isfinite(normal).all():
            raise ValueError(f'Non-finite candidate geometry: {key}')
        # Audit precisely the float32 coordinates written to the candidate binary.
        fitted[key] = (p.astype('<f4').astype(float), normal.astype('<f4').astype(float), f)
    source_bones = {sid: source.mesh(sid) for sid, _ in BONE_PAIRS}
    target_bones = {tid: target.mesh(tid) for _, tid in BONE_PAIRS}
    baseline = fitted
    correction = CanalCorrection(baseline, target_bones)
    fitted = {}
    for key, (v, normal, faces) in baseline.items():
        p, n = correction.map(v, normal)
        fitted[key] = (p.astype('<f4').astype(float), n.astype('<f4').astype(float), faces)
    Surface = surface_helper()
    inventory = []
    blob = bytearray(); records = []
    for key in CORD_IDS:
        v, n, f = original[key]
        inventory.append(dict(id=key, name=source.parts[key]['name'], system=source.parts[key]['system'],
                              bounds=[v.min(axis=0).tolist(), v.max(axis=0).tolist()],
                              rawGeometrySha256=sha(v.astype('<f4').tobytes() + (np.rint(n * 32767).astype('<i2')).tobytes() + f.astype('<u4').tobytes()),
                              topology=topology(v, f)))
        p, normal, faces = fitted[key]
        record = dict(id=key, name=source.parts[key]['name'], system='nervous', vertexCount=len(p), indexCount=faces.size,
                      bounds=[p.min(axis=0).tolist(), p.max(axis=0).tolist()])
        for field, values in [('positions', p.astype('<f4')), ('normals', normal.astype('<f4')), ('indices', faces.astype('<u4'))]:
            record[field] = len(blob); blob.extend(values.tobytes())
        records.append(record)
    metrics = []
    for sid, tid in BONE_PAIRS:
        sv, _, sf = source_bones[sid]; tv, _, tf = target_bones[tid]
        mapped = registration.map(sv)
        sample = mapped[np.linspace(0, len(mapped)-1, min(64, len(mapped)), dtype=int)]
        distances = Surface(tv, tf).distances(sample) * 1000
        metrics.append(dict(sourceId=sid, targetId=tid, samples=len(sample),
                            meanSurfaceDistanceMm=float(np.mean(distances)),
                            p95SurfaceDistanceMm=float(np.quantile(distances, .95)),
                            maxSurfaceDistanceMm=float(distances.max())))
    target_medulla = {key: target.mesh(key) for key in ('FJ1769', 'FJ1831')}
    source_medulla = {key: source.mesh(key) for key in source.parts if key.startswith('Allen_') and 'medulla_oblongata' in key and 'central_canal' not in key}
    interface_source = medulla_interface(original, source_medulla, Surface)
    interface_candidate = medulla_interface(fitted, target_medulla, Surface)
    source_sections = canal_screen(original, source_bones)
    candidate_sections = canal_screen(fitted, target_bones)
    report = dict(schemaVersion=2, candidate='female-spinal-cord', enabled=False,
                  status='Canal-corrected study-model preview candidate; not installed in either production atlas.',
                  inputs=dict(source=source.evidence(), target=target.evidence(), scriptSha256=sha(Path(__file__).read_bytes()),
                              surfaceHelperSha256=sha((ROOT / 'scripts/pelvis-surface-audit.py').read_bytes())),
                  source=dict(dataset='HRA united-female v1.5 / Visible Human Female', license='CC BY 4.0',
                              sourceUrl='https://lod.humanatlas.io/ref-organ/united-female/v1.5',
                              segmentCount=len(CORD_IDS), triangles=sum(len(f) for _, _, f in original.values()),
                              inventory=inventory, missingFromBundledSource=['S5 spinal cord segment', 'Coccygeal spinal cord segment'],
                              note='Inventory labels do not imply a complete spinal cord or nerve-root model. The source also labels six lumbar vertebrae; only its L1 is used for this fit.'),
                  registration=dict(method='Continuous piecewise affine displacement through corresponding vertebral area-weighted surface centroids; fixed transverse scale, positive monotone height mapping, linear endpoint extrapolation.',
                                    anchors=registration.anchors,
                                    minimumJacobianDeterminant=float((1 + registration.slopes[:, 1]).min()),
                                    maximumJacobianDeterminant=float((1 + registration.slopes[:, 1]).max()),
                                    sourceAnchorHeightRangeM=[float(registration.s[0, 1]), float(registration.s[-1, 1])],
                                    extrapolatedVertices=sum(int(((v[:, 1] < registration.s[0, 1]) | (v[:, 1] > registration.s[-1, 1])).sum()) for v, _, _ in original.values()),
                                    surfaceResiduals=metrics,
                                    normalMethod='Inverse-transpose of the local piecewise affine Jacobian, normalized.'),
                  canalCorrection=correction.evidence(),
                  medullaInterface=dict(source=interface_source, candidate=interface_candidate),
                  continuity=dict(method='Adjacent segment axial bounds plus up to 32 sampled end-strip vertices in each direction measured against actual opposite triangles; unsigned proximity cannot prove welded continuity or absence of overlap.',
                                  source=continuity(original, Surface), baseline=continuity(baseline, Surface), candidate=continuity(fitted, Surface)),
                  canalScreen=dict(method='At each segment mid-height intersect actual triangles with a transverse plane; compare outer cord ray extent with first bone hit on 64 rays. Escape counts and odd bone-intersection parity expose incomplete rings or a center in bone. No invented canal mesh.',
                                   source=source_sections, baseline=canal_screen(baseline, target_bones), candidate=candidate_sections,
                                   fitBefore=canal_screen(baseline, target_bones, correction.FIT_FRACTIONS, correction.ids),
                                   fitAfter=canal_screen(fitted, target_bones, correction.FIT_FRACTIONS, correction.ids),
                                   holdoutBefore=canal_screen(baseline, target_bones, correction.HOLDOUT_FRACTIONS, correction.ids),
                                   holdoutAfter=canal_screen(fitted, target_bones, correction.HOLDOUT_FRACTIONS, correction.ids),
                                   limitations=['Twenty matched vertebrae only; discs, ligaments, dura and nerve roots are not modeled by this test.',
                                                'Separate fitting and held-out cross-sections use 64 rays each; crossings between samples remain possible.',
                                                'Bone centroids are reproducible surface descriptors, not anatomical canal landmarks.',
                                                'Source topology is retained, including separate capped or open segment boundaries; none are welded or repaired.',
                                                'The short existing central-canal fragment is not a continuous target canal and is not used as a fit guide.']),
                  blockers=['Incomplete vertebral-section rings and open source boundaries prevent a complete containment claim; these are limitations of this experimental preview.',
                            'S5 and coccygeal segment meshes are absent from the bundled cord allowlist; no synthetic replacement was created.',
                            'Source segment continuity is sampled, not a connected watertight cord reconstruction.'],
                  artifacts=dict(geometry='candidate.bin', manifest='candidate.json', overlay='overlay.svg', baselineOverlay='baseline-overlay.svg', cervicalSections='cervical-sections.svg'))
    report['summary'] = dict(candidateSectionsWithNegativeRadialClearance=sum(r['negativeRadialClearanceDirections'] > 0 for r in candidate_sections),
                            candidateSectionsWithEscapingRays=sum(r['escapingDirections'] > 0 for r in candidate_sections),
                            candidateSectionsWithOddBoneParity=sum(r['oddBoneParityDirections'] > 0 for r in candidate_sections),
                            sourceSectionsWithNegativeRadialClearance=sum(r['negativeRadialClearanceDirections'] > 0 for r in source_sections),
                            sourceSectionsWithEscapingRays=sum(r['escapingDirections'] > 0 for r in source_sections))
    for name in ('fitBefore', 'fitAfter', 'holdoutBefore', 'holdoutAfter'):
        rows = report['canalScreen'][name]
        report['summary'][name] = dict(sections=len(rows), negativeClearanceSections=sum(r['negativeRadialClearanceDirections'] > 0 for r in rows),
                                      minimumRadialClearanceMm=min(r['minimumRadialClearanceMm'] for r in rows if r['minimumRadialClearanceMm'] is not None),
                                      oddBoneParitySections=sum(r['oddBoneParityDirections'] > 0 for r in rows))
    refined_screens = [report['summary'][name] for name in ('fitAfter', 'holdoutAfter')]
    report['previewReadiness'] = dict(eligibleForExperimentalPreview=all(r['negativeClearanceSections'] == 0 and r['oddBoneParitySections'] == 0 for r in refined_screens),
                                    rationale='Fitting and held-out sampled sections clear measured bone boundaries; original topology and gaps remain visible. enabled:false records that this script does not publish the candidate.',
                                    completeCanalContainmentEstablished=False)
    (output / 'candidate.bin').write_bytes(blob)
    (output / 'candidate.json').write_text(json.dumps(dict(enabled=False, coordinateSystem='meters, Y up; current female atlas frame',
                                                          source='HRA united-female v1.5', license='CC BY 4.0', normalType='float32', indexType='uint32',
                                                          binary='candidate.bin', binarySha256=sha(blob), parts=records), indent=2) + '\n')
    svg_overlay(output / 'overlay.svg', fitted, {**target_bones, **target_medulla})
    svg_overlay(output / 'baseline-overlay.svg', baseline, {**target_bones, **target_medulla})
    section_comparison_svg(output / 'cervical-sections.svg', baseline, fitted, target_bones)
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(dict(output=str(output), enabled=False, **report['summary'])))


if __name__ == '__main__':
    main()
