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


def canal_screen(cords, bones):
    bone_triangles = np.concatenate([v[f] for v, _, f in bones.values()])
    result = []
    for key in CORD_IDS:
        v, _, f = cords[key]
        height = float((v[:, 1].min() + v[:, 1].max()) / 2)
        cord = section(v[f], height)
        if not len(cord):
            raise ValueError(f'No cord cross-section for {key}')
        center = (cord.min(axis=(0, 1)) + cord.max(axis=(0, 1))) / 2
        bone = section(bone_triangles, height)
        clearances = []; escaping = 0; bone_parities = []
        for theta in np.linspace(.013, 2 * np.pi + .013, 64, endpoint=False):
            direction = np.array([np.cos(theta), np.sin(theta)])
            bh = ray_hits(bone, center, direction)
            ch = ray_hits(cord, center, direction)
            bone_parities.append(len(bh) % 2)
            if not len(bh) or not len(ch):
                escaping += 1
            else:
                clearances.append(float((bh[0] - ch[-1]) * 1000))
        result.append(dict(id=key, heightM=height, centerXZ=center.tolist(),
                           boneCrossSectionSegments=len(bone), cordCrossSectionSegments=len(cord),
                           radialDirections=64, escapingDirections=escaping,
                           oddBoneParityDirections=sum(bone_parities),
                           minimumRadialClearanceMm=min(clearances) if clearances else None,
                           medianRadialClearanceMm=float(np.median(clearances)) if clearances else None,
                           negativeRadialClearanceDirections=sum(c < 0 for c in clearances)))
    return result


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
    lines += ['<text x="30" y="930" fill="white" font-size="13">Grey: current vertebrae; gold: registered HRA cord. Bone-centroid field is not a canal landmark fit.</text>', '</svg>']
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
    for name in ('report.json', 'candidate.json', 'candidate.bin', 'overlay.svg'):
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
    source_sections = canal_screen(original, source_bones)
    candidate_sections = canal_screen(fitted, target_bones)
    report = dict(schemaVersion=1, candidate='female-spinal-cord', enabled=False,
                  status='Isolated measured registration candidate; not installed in either production atlas.',
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
                  continuity=dict(method='Adjacent segment axial bounds plus up to 32 sampled end-strip vertices in each direction measured against actual opposite triangles; unsigned proximity cannot prove welded continuity or absence of overlap.',
                                  source=continuity(original, Surface), candidate=continuity(fitted, Surface)),
                  canalScreen=dict(method='At each segment mid-height intersect actual triangles with a transverse plane; compare outer cord ray extent with first bone hit on 64 rays. Escape counts and odd bone-intersection parity expose incomplete rings or a center in bone. No invented canal mesh.',
                                   source=source_sections, candidate=candidate_sections,
                                   limitations=['Twenty matched vertebrae only; discs, ligaments, dura and nerve roots are not modeled by this test.',
                                                'One cross-section per segment and 64 rays cannot exclude crossings between samples.',
                                                'Bone centroids are reproducible surface descriptors, not anatomical canal landmarks.',
                                                'Source topology is retained, including separate capped or open segment boundaries; none are welded or repaired.',
                                                'The short existing central-canal fragment is not a continuous target canal and is not used as a fit guide.']),
                  blockers=['Surface registration and sampled canal measurements need review before production integration.',
                            'S5 and coccygeal segment meshes are absent from the bundled cord allowlist; no synthetic replacement was created.',
                            'Source segment continuity is sampled, not a connected watertight cord reconstruction.'],
                  artifacts=dict(geometry='candidate.bin', manifest='candidate.json', overlay='overlay.svg'))
    report['summary'] = dict(candidateSectionsWithNegativeRadialClearance=sum(r['negativeRadialClearanceDirections'] > 0 for r in candidate_sections),
                            candidateSectionsWithEscapingRays=sum(r['escapingDirections'] > 0 for r in candidate_sections),
                            candidateSectionsWithOddBoneParity=sum(r['oddBoneParityDirections'] > 0 for r in candidate_sections),
                            sourceSectionsWithNegativeRadialClearance=sum(r['negativeRadialClearanceDirections'] > 0 for r in source_sections),
                            sourceSectionsWithEscapingRays=sum(r['escapingDirections'] > 0 for r in source_sections))
    (output / 'candidate.bin').write_bytes(blob)
    (output / 'candidate.json').write_text(json.dumps(dict(enabled=False, coordinateSystem='meters, Y up; current female atlas frame',
                                                          source='HRA united-female v1.5', license='CC BY 4.0', normalType='float32', indexType='uint32',
                                                          binary='candidate.bin', binarySha256=sha(blob), parts=records), indent=2) + '\n')
    svg_overlay(output / 'overlay.svg', fitted, target_bones)
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(dict(output=str(output), enabled=False, **report['summary'])))


if __name__ == '__main__':
    main()
