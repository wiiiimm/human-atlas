"""Verify experimental fitting cannot silently change canonical model files."""
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('enhancement_audit', Path(__file__).with_name('audit-female-enhancements.py'))
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class RunnerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / 'repo'
        self.models = self.root / 'public/models'
        self.models.mkdir(parents=True)
        self.asset = self.models / 'geometry.bin'
        self.asset.write_bytes(b'original')
        self.output = Path(self.temp.name) / 'candidate'
        self.root_patch = patch.object(audit, 'ROOT', self.root)
        self.root_patch.start()
        self.addCleanup(self.root_patch.stop)

    def test_repository_and_symlink_output_rejected_before_execution(self):
        link = Path(self.temp.name) / 'alias'
        link.symlink_to(self.models, target_is_directory=True)
        with patch.object(audit.subprocess, 'run') as child:
            for output in (self.root, self.models, link / 'candidate'):
                with self.assertRaises(ValueError):
                    audit.run(output)
            child.assert_not_called()

    def test_modified_model_is_reported_even_when_child_succeeds(self):
        def child(args, **kwargs):
            output = Path(args[-1]); output.mkdir(parents=True)
            (output / 'report.json').write_text('{"enabled":false}')
            self.asset.write_bytes(b'changed')
            return SimpleNamespace(returncode=0)
        with patch.object(audit.subprocess, 'run', side_effect=child):
            with self.assertRaisesRegex(RuntimeError, 'changed public model assets'):
                audit.run(self.output)
        self.assertFalse(json.loads((self.output / 'summary.json').read_text())['canonicalAssetsUnchanged'])

    def test_failed_experiment_does_not_report_overall_success(self):
        with patch.object(audit.subprocess, 'run', return_value=SimpleNamespace(returncode=2)):
            self.assertEqual(audit.run(self.output), 1)
        report = json.loads((self.output / 'summary.json').read_text())
        self.assertTrue(report['canonicalAssetsUnchanged'])
        self.assertTrue(all(not item['completed'] for item in report['experiments'].values()))


if __name__ == '__main__':
    unittest.main()
