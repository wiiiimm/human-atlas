"""Run the four HRA enhancement experiments without publishing their geometry."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
EXPERIMENTS = {
    'spinal-cord': 'female_spinal_candidate.py',
    'knee-and-quadriceps': 'female_joint_candidates.py',
    'kidney-internals': 'female_kidney_candidate.py',
}


def fingerprint():
    """Protect both viewer atlases, rebuild inputs, geometry and fit provenance."""
    files = sorted((ROOT / 'public/models').glob('*'))
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in files if p.is_file()}


def run(output):
    output = output.resolve()
    # Candidate payloads must never overwrite source or canonical viewer files.
    if output == ROOT or ROOT in output.parents:
        raise ValueError('Use an output directory outside the repository, such as /tmp/female-enhancements')
    output.mkdir(parents=True, exist_ok=True)
    before = fingerprint()
    results = {}
    try:
        for name, script in EXPERIMENTS.items():
            target = output / name
            result = subprocess.run([sys.executable, str(ROOT / 'scripts' / script),
                                     '--output-dir', str(target)], cwd=ROOT)
            if result.returncode:
                results[name] = {'completed': False, 'exitCode': result.returncode}
                continue
            report = json.loads((target / 'report.json').read_text())
            results[name] = {'completed': True, 'report': str(target / 'report.json'),
                             'reportSha256': hashlib.sha256((target / 'report.json').read_bytes()).hexdigest(),
                             'enabled': report.get('enabled', False)}
    finally:
        unchanged = fingerprint() == before
        summary = {'schemaVersion': 1, 'canonicalAssetsUnchanged': unchanged,
                   'canonicalAssetHashes': before, 'experiments': results,
                   'publication': 'No viewer publication is performed by this runner.'}
        (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    if not unchanged:
        raise RuntimeError('An experiment changed public model assets; inspect changes before proceeding')
    if len(results) != len(EXPERIMENTS) or any(not r['completed'] for r in results.values()):
        return 1
    print(json.dumps({'canonicalAssetsUnchanged': True, 'experiments': results}, indent=2))
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    sys.exit(run(args.output_dir))
