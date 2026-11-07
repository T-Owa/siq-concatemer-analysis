#!/usr/bin/env python3
"""Reproduce SIQ post-processing from the supplied portable event tables."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--outdir', type=Path, default=Path('results'))
    p.add_argument('--dataset', choices=['all', 'modify', 'inhibitor'], default='all')
    p.add_argument('--data-dir', type=Path, default=ROOT / 'data')
    p.add_argument('--config', type=Path, default=ROOT / 'config/junction.json')
    p.add_argument('--catalog', type=Path, default=ROOT / 'config/pattern_catalog.csv')
    p.add_argument('--figures', action='store_true', help='Also make static SVG/PDF/PNG figures (matplotlib required)')
    a = p.parse_args()
    out = a.outdir.resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        p.error('Choose a new or empty --outdir; existing results are not overwritten.')
    if sys.version_info < (3, 10):
        p.error('Python 3.10 or newer is required.')
    datasets = ['modify', 'inhibitor'] if a.dataset == 'all' else [a.dataset]
    for dataset in datasets:
        for required in [a.data_dir / (dataset + '_events.csv.gz'), a.data_dir / 'samples.csv', a.config, a.catalog]:
            if not required.is_file():
                p.error('Missing input: ' + str(required))
    out.mkdir(parents=True, exist_ok=True)
    commands = []

    def run(script, args):
        cmd = [sys.executable, str(ROOT / 'scripts' / script), *map(str, args)]
        print('Running ' + script, flush=True)
        subprocess.run(cmd, check=True)
        # Record relative scientific inputs, without local user paths.
        commands.append(script)

    for dataset in datasets:
        common = ['--dataset', dataset, '--samples', a.data_dir / 'samples.csv', '--config', a.config]
        dest = out / dataset
        run('annotate_events.py', ['--input', a.data_dir / (dataset + '_events.csv.gz'), *common, '--outdir', dest / 'mh'])
        run('summarize_events.py', ['--events', dest / 'mh/events.csv', *common, '--outdir', dest / 'summary'])
        run('aggregate_patterns.py', ['--events', dest / 'mh/events.csv', *common, '--catalog', a.catalog, '--outdir', dest / 'patterns'])
    if a.figures:
        run('plot_figures.py', ['--results', out, '--outdir', out / 'figures'])
    manifest = {
        'package_version': '0.1.0', 'python_version': platform.python_version(),
        'datasets': datasets, 'scripts_executed': commands,
        'config_sha256': hashlib.sha256(a.config.read_bytes()).hexdigest(),
        'catalog_sha256': hashlib.sha256(a.catalog.read_bytes()).hexdigest(),
        'input_sha256': {x.name: hashlib.sha256(x.read_bytes()).hexdigest()
                         for x in [a.data_dir / 'samples.csv', *[a.data_dir / (d + '_events.csv.gz') for d in datasets]]},
        'scope': 'Post-processing of SIQ-reported events; does not rerun read merging or SIQ.'
    }
    (out / 'run_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    print('Completed: ' + str(out))


if __name__ == '__main__':
    main()
