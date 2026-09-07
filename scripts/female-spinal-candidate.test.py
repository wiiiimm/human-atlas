#!/usr/bin/env python3
"""Focused correction geometry tests; optional verification of actual exported data."""
import argparse
import importlib.util
import json
from pathlib import Path
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('spinal', ROOT / 'scripts/female_spinal_candidate.py')
spinal = importlib.util.module_from_spec(spec)
spec.loader.exec_module(spinal)


class CorrectionTests(unittest.TestCase):
    def correction(self):
        c = spinal.CanalCorrection.__new__(spinal.CanalCorrection)
        c.radius = .035
        c.centers = np.linspace(1.39, 1.45, 5)
        c.coefficients = np.array([[-.0026, .0021], [-.0011, .0019], [-.002, .002], [-.0013, -.0009], [0, 0]])
        return c

    def test_smooth_compact_field_and_normal_jacobian(self):
        c = self.correction()
        points = np.array([[.003, h, -.035] for h in [1.2, 1.354, 1.39, 1.4, 1.417, 1.43, 1.485, 1.6]])
        normals = np.tile(np.array([.3, .4, .5]) / np.sqrt(.5), (len(points), 1))
        mapped, result = c.map(points, normals)
        step = 1e-7
        jacobian = np.stack([(c.map(points + np.eye(3)[axis] * step) - c.map(points - np.eye(3)[axis] * step)) / (2 * step) for axis in range(3)], axis=2)
        expected = np.linalg.solve(jacobian.transpose(0, 2, 1), normals[..., None])[..., 0]
        expected /= np.linalg.norm(expected, axis=1)[:, None]
        np.testing.assert_allclose(result, expected, atol=1e-8)
        np.testing.assert_allclose(np.linalg.det(jacobian), 1, atol=2e-9)
        np.testing.assert_array_equal(mapped[[0, 1, 6, 7]], points[[0, 1, 6, 7]])
        np.testing.assert_array_equal(result[[0, 1, 6, 7]], normals[[0, 1, 6, 7]])
        for edge in [c.centers[0] - c.radius, c.centers[-1] + c.radius]:
            self.assertLess(float(abs(c.basis(np.array([edge]), derivative=True)).max()), 1e-10)

    def test_radial_search_matches_independent_screen_with_duplicate_seam(self):
        square = np.array([[[-1., -1.], [1., -1.]], [[1., -1.], [1., 1.]], [[1., 1.], [-1., 1.]], [[-1., 1.], [-1., -1.]]])
        # Duplicate geometric seam must not flip odd/even bone parity.
        square = np.concatenate([square, square[1:2]])
        angles = np.linspace(.013, 2*np.pi + .013, 64, endpoint=False)
        directions = np.column_stack([np.cos(angles), np.sin(angles)])
        for center in [np.zeros(2), np.array([2., .3])]:
            first, parity = spinal.radial_boundaries(square, center, directions)
            for i, direction in enumerate(directions):
                hits = spinal.ray_hits(square, center, direction)
                self.assertEqual(parity[i], len(hits) % 2)
                self.assertAlmostEqual(first[i], hits[0] if len(hits) else np.inf)


def validate_candidate(directory):
    manifest = json.loads((directory / 'candidate.json').read_text())
    report = json.loads((directory / 'report.json').read_text())
    blob = (directory / 'candidate.bin').read_bytes()
    assert spinal.sha(blob) == manifest['binarySha256']
    source = spinal.Atlas(ROOT / 'public/models', 'atlas-female.json')
    target = spinal.Atlas(ROOT / 'public/models', 'atlas-female-reconstructed.json')
    registration = spinal.Registration(source, target)
    parameters = report['canalCorrection']
    correction = spinal.CanalCorrection.__new__(spinal.CanalCorrection)
    correction.radius = parameters['basisRadiusM']
    correction.centers = np.array(parameters['basisCentersM'])
    correction.coefficients = np.array(parameters['coefficientsXZ_M'])
    assert {p['id'] for p in manifest['parts']} == set(spinal.CORD_IDS)
    for part in manifest['parts']:
        v, n, f = source.mesh(part['id'])
        positions = np.frombuffer(blob, '<f4', len(v)*3, part['positions']).reshape(-1, 3)
        normals = np.frombuffer(blob, '<f4', len(v)*3, part['normals']).reshape(-1, 3)
        faces = np.frombuffer(blob, '<u4', f.size, part['indices']).reshape(-1, 3)
        np.testing.assert_array_equal(faces, f)
        baseline, normal = registration.map(v, n)
        baseline = baseline.astype('<f4').astype(float)
        normal = normal.astype('<f4').astype(float)
        expected, expected_normal = correction.map(baseline, normal)
        np.testing.assert_array_equal(positions, expected.astype('<f4'))
        np.testing.assert_array_equal(normals, expected_normal.astype('<f4'))
        assert np.isfinite(positions).all() and np.isfinite(normals).all()
        np.testing.assert_allclose(np.linalg.norm(normals, axis=1), 1, atol=2e-7)
    for name in ('fitAfter', 'holdoutAfter'):
        assert report['summary'][name]['negativeClearanceSections'] == 0
        assert report['summary'][name]['oddBoneParitySections'] == 0
    assert report['previewReadiness']['eligibleForExperimentalPreview']
    print('All 29 exported meshes retain source topology and match the recorded smooth field; fitting and holdout checks pass.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--candidate-dir', type=Path)
    args, remaining = parser.parse_known_args()
    result = unittest.main(argv=[__file__, *remaining], exit=False)
    if not result.result.wasSuccessful():
        raise SystemExit(1)
    if args.candidate_dir:
        validate_candidate(args.candidate_dir)
