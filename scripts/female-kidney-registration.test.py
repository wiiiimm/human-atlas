#!/usr/bin/env python3
"""Analytic checks for kidney surface proposals; no atlas writes."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np
from female_kidney_candidate import require_output_directory, run as kidney_run
from female_kidney_registration import radial_outer_patch, constrain_matrix, weighted_affine, fit_surface_affine

spec = importlib.util.spec_from_file_location('surface_audit', Path(__file__).with_name('pelvis-surface-audit.py'))
audit = importlib.util.module_from_spec(spec); spec.loader.exec_module(audit)


def box(extents=(1., 1., 1.)):
    vertices = np.array([[-1,-1,-1], [1,-1,-1], [1,1,-1], [-1,1,-1], [-1,-1,1], [1,-1,1], [1,1,1], [-1,1,1]], dtype=float)*np.array(extents)
    faces = np.array([[0,2,1], [0,3,2], [4,5,6], [4,6,7], [0,1,5], [0,5,4], [3,7,6], [3,6,2], [0,4,7], [0,7,3], [1,2,6], [1,6,5]])
    return vertices, faces


class KidneyRegistrationTests(unittest.TestCase):
    def test_radial_patch_excludes_inner_material_surface_without_filling(self):
        vertices, faces = box()
        shell_vertices = np.concatenate([vertices, vertices*.9])
        shell_faces = np.concatenate([faces, faces[:, ::-1]+8])
        chosen = radial_outer_patch(shell_vertices, shell_faces)
        np.testing.assert_array_equal(chosen, np.arange(12))
        np.testing.assert_array_equal(radial_outer_patch(shell_vertices+np.array([.123, 1.234, -.15]), shell_faces), chosen)

    def test_bounded_affine_cannot_reflect_or_collapse_assembly(self):
        for proposed in [np.diag([.1, .2, 2]), np.diag([-.8, .7, .9]), np.ones((3, 3))]:
            bounded = constrain_matrix(proposed)
            singular = np.linalg.svd(bounded, compute_uv=False)
            self.assertGreaterEqual(float(singular.min()), .55-1e-10)
            self.assertLessEqual(float(singular.max()), 1.10+1e-10)
            self.assertGreaterEqual(np.linalg.det(bounded), .35-1e-10)

    def test_shared_affine_recovers_correspondences_and_preserves_normal_tangency(self):
        source, _ = box()
        matrix = np.array([[.8, .06, 0], [0, .9, .02], [.01, 0, .85]])
        offset = np.array([.1, .2, -.04])
        target = source @ matrix.T+offset
        result, shift = weighted_affine(source, target, np.ones(len(source)), matrix)
        np.testing.assert_allclose(result, matrix, atol=1e-12)
        np.testing.assert_allclose(shift, offset, atol=1e-12)
        normal = np.cross(source[1]-source[0], source[2]-source[0]) @ np.linalg.inv(result)
        for tangent in [target[1]-target[0], target[2]-target[0]]:
            self.assertAlmostEqual(float(normal @ tangent), 0., places=12)

    def test_triangle_surface_icp_recovers_exact_similarity_without_interface(self):
        source, faces = box((.02, .05, .01))
        target = source*.8+np.array([.08, 1.12, -.06])
        matrix, offset, report = fit_surface_affine(source, faces, target, faces, np.empty((0, 3)), None, audit.Surface, maximum_samples=8, iterations=3)
        # Box symmetry allows several rotations; actual surface agreement is the invariant.
        moved = source @ matrix.T+offset
        self.assertLess(float(audit.Surface(target, faces).distances(moved).max()), 1e-10)
        self.assertLess(float(audit.Surface(moved, faces).distances(target).max()), 1e-10)
        self.assertGreater(np.linalg.det(matrix), .35)
        self.assertEqual(len(report['seeds']), 5)


class KidneyOutputDirectoryTests(unittest.TestCase):
    def test_symbolic_link_output_directory_is_rejected_before_resolution(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            real = tmp / 'real'
            real.mkdir()
            link = tmp / 'link'
            link.symlink_to(real, target_is_directory=True)
            with self.assertRaisesRegex(ValueError, 'symbolic-link path components'):
                require_output_directory(link)
            with self.assertRaisesRegex(ValueError, 'symbolic-link path components'):
                kidney_run(link)
            self.assertFalse((real / 'report.json').exists())
            self.assertEqual(list(real.iterdir()), [])

    def test_ordinary_external_directory_still_resolves(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / 'candidate'
            self.assertEqual(require_output_directory(output), output.resolve())


if __name__ == '__main__': unittest.main()
