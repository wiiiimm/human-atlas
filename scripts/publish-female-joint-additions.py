#!/usr/bin/env python3
"""Build a temporary append-only bilateral ALL static preview; never publish it."""
import argparse
import copy
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
IDS = ('VH_F_anterolateral_ligament_of_knee_L', 'VH_F_anterolateral_ligament_of_knee_R')
APPROVED_REPORT_SHA = 'dd1e16dcbdb49d497128c2db7bad915cba2726fa3362bcc9bf8d6d11cb0ffdd1'
APPROVED_NPZ_SHA = dict(zip(IDS, ['f849d3123f45b96f4f1670d3851bb920d071ea0ffc2297b87c83b47e7e2900af', '1ef1a1b402b6d9526a70bd43db5ee75848b2251b23d15af41cce0a3645bfeeb8']))
STATUS = 'Experimental static study preview; joint attachments and movement are not validated'
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location('joint_append_surface', ROOT/'scripts/pelvis-surface-audit.py')
a = importlib.util.module_from_spec(spec); spec.loader.exec_module(a)


def sha(data): return hashlib.sha256(data).hexdigest()
def encoded(value): return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()
def require(condition, message):
    if not condition: raise ValueError(message)


def model_path(models_dir, url):
    require(url.startswith('/models/') and Path(url).name == url[len('/models/'):], 'Unsafe model URL')
    return models_dir/Path(url).name


def arrays(models_dir, manifest, part):
    data = model_path(models_dir, manifest['chunks'][part['chunk']]['url']).read_bytes()
    n = part['vertexCount']; count = part['indexCount']
    return (np.frombuffer(data, '<f4', n*3, part['positions']).reshape(-1, 3).copy(),
            np.frombuffer(data, '<i2', n*3, part['normals']).reshape(-1, 3).copy(),
            np.frombuffer(data, '<u4', count, part['indices']).reshape(-1, 3).copy())


def area_normals(positions, faces):
    normals = np.zeros_like(positions, dtype=float)
    triangles = positions.astype(float)[faces]
    cross = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
    for axis in range(3): np.add.at(normals, faces[:, axis], cross)
    lengths = np.linalg.norm(normals, axis=1)
    require(np.all(lengths > 1e-20), 'Degenerate vertex normals')
    return normals/lengths[:, None]


def reproduce(source, side, morph, key):
    fit = side['fit']; refinement = side['refinement']
    require(refinement['ridge'] == .08 and refinement['kernelRadiusM'] == .025, 'Only reviewed .08 ridge / .025 m radius configuration is allowed')
    points = source.astype(float)
    centers = np.asarray(refinement['controlSourcePointsM'], dtype=float)
    coefficients = np.asarray(refinement['coefficients'], dtype=float)
    require(centers.ndim == 2 and centers.shape[1] == 3 and centers.shape == coefficients.shape, 'Invalid RBF controls')
    require(np.isfinite(centers).all() and np.isfinite(coefficients).all(), 'Nonfinite RBF controls')
    moved = fit['scale']*points @ np.asarray(fit['rotationRowVectors'])+np.asarray(fit['translationM'])
    weights = np.exp(-np.sum((points[:, None]-centers[None])**2, axis=2)/(2*refinement['kernelRadiusM']**2))
    return a.morph(moved+weights @ coefficients, morph, key).astype('<f4')


def run(candidate_dir, output_dir, root=ROOT, base_models_dir=None, source_models_dir=None):
    root = root.resolve(); candidate_dir = candidate_dir.resolve(); output_dir = output_dir.resolve()
    for protected in [ROOT, root, candidate_dir]:
        require(output_dir != protected and protected not in output_dir.parents, 'Output must be outside the repository and candidate input directory')
    base_models_dir = (base_models_dir or root/'public/models').resolve()
    source_models_dir = (source_models_dir or root/'public/models').resolve()
    require(output_dir != source_models_dir or output_dir == base_models_dir, 'Cannot write source model directory')
    manifest_path = base_models_dir/'atlas-female-reconstructed.json'
    fit_path = base_models_dir/'female-fit-report.json'
    source_path = source_models_dir/'atlas-female.json'
    manifest_bytes = manifest_path.read_bytes(); fit_bytes = fit_path.read_bytes(); source_bytes = source_path.read_bytes()
    base = json.loads(manifest_bytes); fit = json.loads(fit_bytes); source = json.loads(source_bytes)
    report_bytes = (candidate_dir/'report.json').read_bytes(); report = json.loads(report_bytes)
    require(sha(report_bytes) == APPROVED_REPORT_SHA, 'Candidate report is not the pinned reviewed artifact')
    require(report.get('existingAtlasByteIdentityVerified') is True, 'Candidate did not verify original atlas identity')
    required_hashes = {'public/models/atlas-female-reconstructed.json': sha(manifest_bytes), 'public/models/female-fit-report.json': sha(fit_bytes), 'public/models/atlas-female.json': sha(source_bytes)}
    for path, digest in required_hashes.items(): require(report['inputHashes'].get(path) == digest, f'Stale required input: {path}')
    base_names = {manifest_path.name, fit_path.name} | {Path(chunk[key]).name for chunk in base['chunks'] for key in ['url', 'gzip']}
    for relative, digest in report['inputHashes'].items():
        require(relative.startswith('public/models/') and Path(relative).name == relative[len('public/models/'):], 'Unsafe input hash path')
        path = (base_models_dir if Path(relative).name in base_names else source_models_dir)/Path(relative).name
        require(sha(path.read_bytes()) == digest, f'Stale candidate input: {relative}')
    require(not any(p['id'] in IDS for p in base['parts']), 'ALL additions already present')
    require('jointAdditions' not in fit and 'jointAdditions' not in base, 'Existing joint additions require explicit migration')
    concept_ids = {c['id'] for c in base['concepts']}
    require(not any('PART_'+key in concept_ids for key in IDS), 'Addition concept already present')
    source_parts = {p['id']:p for p in source['parts']}
    config = dict(schemaVersion=1, ids=list(IDS), method=report['method'], morph=fit['morph'],
                  sides={side:{key:report['sides'][side][key] for key in ['fit', 'refinement']} for side in ['L', 'R']},
                  normalQuantization='Normalize candidate float32 normals in float64, round to nearest int16 using scale 32767')
    base_files = {}; source_chunks = {}; payload = bytearray(); new_parts = []; records = []
    for chunk in base['chunks']:
        for field in ['url', 'gzip']:
            url = chunk[field]; base_files[url] = model_path(base_models_dir, url).read_bytes()
    base_hashes = {url:sha(data) for url, data in base_files.items()}
    def pack(data):
        payload.extend(b'\0'*((-len(payload)) % 4)); offset = len(payload); payload.extend(data); return offset
    for key in IDS:
        side = key[-1]; original = source_parts[key]
        require(original['name'] == 'anterolateral ligament of knee', 'Unexpected source identity')
        rows = [row for row in report['sides'][side]['parts'] if row['id'] == key]
        require(len(rows) == 1, 'Missing or duplicate selected candidate')
        row = rows[0]; export = row['export']
        require(export['npz'] == key+'.npz', 'Unexpected candidate asset path')
        require(row.get('topologyPreserved') is True, 'Candidate topology flag missing')
        require(row['distortion']['nonpositiveJacobianSamples'] == 0 and row['distortion']['minimumJacobianDeterminant'] > 0, 'Candidate contains nonpositive sampled Jacobian')
        candidate_bytes = (candidate_dir/export['npz']).read_bytes()
        require(sha(candidate_bytes) == APPROVED_NPZ_SHA[key], 'Candidate NPZ is not the pinned reviewed artifact')
        with np.load(io.BytesIO(candidate_bytes), allow_pickle=False) as data:
            require(set(data.files) == {'positions', 'normals', 'indices'}, 'Unexpected candidate arrays')
            positions = data['positions']; normals = data['normals']; faces = data['indices']
        source_positions, source_normals, source_faces = arrays(source_models_dir, source, original)
        require(positions.dtype == np.dtype('<f4') and normals.dtype == np.dtype('<f4') and faces.dtype == np.dtype('<u4'), 'Unexpected candidate array dtype')
        require(positions.shape == source_positions.shape and normals.shape == positions.shape and faces.shape == source_faces.shape, 'Source array shapes changed')
        require(np.array_equal(faces, source_faces), 'Source triangle indices changed')
        require(np.isfinite(positions).all() and np.isfinite(normals).all(), 'Nonfinite candidate arrays')
        require(faces.size > 0 and int(faces.max()) < len(positions), 'Invalid index range')
        require(export['vertices'] == len(positions) and export['triangles'] == len(faces), 'Candidate counts changed')
        require(sha(positions.tobytes()) == export['positionsSha256'] and sha(faces.tobytes()) == export['indicesSha256'], 'Candidate raw array digest mismatch')
        require(np.array_equal(positions, reproduce(source_positions, report['sides'][side], fit['morph'], key)), 'Candidate does not reproduce the recorded source transform')
        bounds = np.stack([positions.min(0), positions.max(0)])
        require(np.allclose(bounds, export['boundsM'], rtol=0, atol=1e-7), 'Candidate bounds mismatch')
        lengths = np.linalg.norm(normals.astype(float), axis=1)
        require(np.all(abs(lengths-1) < 1e-4), 'Candidate normals are not unit length')
        normals_unit = normals.astype(float)/lengths[:, None]
        require(np.all(np.sum(normals_unit*area_normals(positions, faces), axis=1) > .995), 'Candidate normal direction does not match final geometry')
        quantized = np.rint(normals_unit*32767).astype('<i2')
        require(np.all(abs(np.linalg.norm(quantized.astype(float)/32767, axis=1)-1) < 5e-5), 'Quantized normals are not unit length')
        po = pack(positions.tobytes()); no = pack(quantized.tobytes()); ix = pack(faces.tobytes())
        name = ('Left' if side == 'L' else 'Right')+' '+original['name']
        provenance = dict(source='HRA united-female v1.5', sourceId=key, sourceConceptId=original.get('conceptId'), license='CC BY 4.0', adaptation=STATUS+'; shared bone-derived Gaussian RBF fit followed by the existing female morph once', reportSection='jointAdditions')
        new_parts.append(dict(id=key, name=name, conceptId='PART_'+key, system='connective', chunk=len(base['chunks']), positions=po, normals=no, indices=ix, vertexCount=len(positions), indexCount=faces.size, bounds=bounds.tolist(), provenance=provenance))
        source_url = source['chunks'][original['chunk']]['url']; source_chunks[source_url] = sha(model_path(source_models_dir, source_url).read_bytes())
        records.append(dict(id=key, sourceId=key, candidateNpzSha256=sha(candidate_bytes),
                            sourceArrayHashes=dict(positions=sha(source_positions.tobytes()), normals=sha(source_normals.tobytes()), indices=sha(source_faces.tobytes())),
                            finalArrayHashes=dict(positions=sha(positions.tobytes()), normals=sha(quantized.tobytes()), indices=sha(faces.tobytes())),
                            vertices=len(positions), triangles=len(faces), bounds=bounds.tolist(),
                            checks=dict(sourceTopologyIdentical=True, recordedSourceTransformReproduced=True, finitePositions=True, unitNormalTolerance=1e-4, quantizedUnitNormalTolerance=5e-5)))
    chunk_bytes = bytes(payload); chunk_name = 'female-joint-additions-'+sha(chunk_bytes)[:16]+'.bin'
    chunk_url = '/models/'+chunk_name
    require(chunk_url not in base_files, 'New chunk collides with existing URL')
    compressed = io.BytesIO()
    with gzip.GzipFile(fileobj=compressed, mode='wb', filename='', mtime=0, compresslevel=9) as archive: archive.write(chunk_bytes)
    gzip_bytes = compressed.getvalue()
    manifest = copy.deepcopy(base); result_fit = copy.deepcopy(fit)
    manifest['chunks'].append(dict(url=chunk_url, bytes=len(chunk_bytes), gzip=chunk_url+'.gz', gzipBytes=len(gzip_bytes)))
    manifest['parts'].extend(new_parts)
    manifest['concepts'].extend(dict(id='PART_'+p['id'], name=p['name'], elements=[p['id']]) for p in new_parts)
    manifest['triangles'] += sum(p['indexCount']//3 for p in new_parts)
    manifest['optimized']['preservedMeshes'] = len(manifest['parts'])
    manifest['jointAdditions'] = dict(ids=list(IDS), report='/models/female-fit-report.json', reviewStatus=STATUS)
    result_fit.update(parts=len(manifest['parts']), concepts=len(manifest['concepts']), triangles=manifest['triangles'])
    result_fit['jointAdditions'] = dict(schemaVersion=1, ids=list(IDS), reviewStatus=STATUS,
        baseManifestSha256=sha(manifest_bytes), baseFitReportSha256=sha(fit_bytes),
        baseCounts=dict(parts=len(base['parts']), concepts=len(base['concepts']), triangles=base['triangles'], chunks=len(base['chunks'])),
        baseChunkHashes=base_hashes, sourceManifestSha256=sha(source_bytes), sourceChunkHashes=source_chunks,
        candidateReportSha256=sha(report_bytes), candidateConfigSha256=sha(encoded(config)), config=config, parts=records,
        chunk=dict(url=chunk_url, sha256=sha(chunk_bytes), gzipSha256=sha(gzip_bytes)))
    files = {Path(url).name:data for url, data in base_files.items()}
    files.update({chunk_name:chunk_bytes, chunk_name+'.gz':gzip_bytes, 'atlas-female-reconstructed.json':(json.dumps(manifest, separators=(',', ':'), allow_nan=False)+'\n').encode(), 'female-fit-report.json':json.dumps(result_fit, indent=2, allow_nan=False).encode()})
    output_dir.mkdir(parents=True, exist_ok=True)
    require(all(not (output_dir/name).is_symlink() for name in files), 'Output file symlinks are not allowed')
    # All candidate validation finishes before the first artifact is written.
    for name, data in files.items(): (output_dir/name).write_bytes(data)
    require(manifest['parts'][:-2] == base['parts'] and manifest['concepts'][:-2] == base['concepts'] and manifest['chunks'][:-1] == base['chunks'], 'Append-only invariant failed')
    if output_dir != base_models_dir:
        require(manifest_path.read_bytes() == manifest_bytes and fit_path.read_bytes() == fit_bytes, 'Canonical metadata changed during build')
    for url, expected in base_hashes.items(): require(sha(model_path(base_models_dir, url).read_bytes()) == expected, 'Canonical geometry changed during build')
    return manifest, result_fit


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--base-models-dir', type=Path, default=ROOT/'public/models')
    parser.add_argument('--source-models-dir', type=Path, default=ROOT/'public/models')
    args = parser.parse_args()
    manifest, _ = run(args.candidate_dir, args.output_dir, base_models_dir=args.base_models_dir, source_models_dir=args.source_models_dir)
    print(json.dumps(dict(parts=len(manifest['parts']), concepts=len(manifest['concepts']), triangles=manifest['triangles'], added=list(IDS), publication='Temporary build only')))
