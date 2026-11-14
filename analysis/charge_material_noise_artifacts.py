#!/usr/bin/env python3
"""Freeze native material-noise source proof and measured worker resources."""
import argparse
import json
from pathlib import Path
import re
import subprocess

from slope_tables import write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, output = args.input, args.output
    manifest = json.loads((root/'manifest.json').read_text())
    output.mkdir(exist_ok=False)
    jobs = []
    for point in manifest['points']:
        folder = Path(point['output'])
        assert (folder/'exit_code.txt').read_text().strip() == '0'
        resource = (folder/'resource.txt').read_text()
        fields = {}
        for label, key in [('User time (seconds)', 'user_cpu_s'), ('System time (seconds)', 'system_cpu_s'),
                           ('Maximum resident set size (kbytes)', 'max_RSS_KiB')]:
            fields[key] = float(re.search(re.escape(label)+r':\s*([\d.]+)', resource)[1])
        jobs.append(dict(output=str(folder), variant=point['variant'], **fields))
    unit = 'Tracking/Acts/FaserActsKalmanFilter/src/CalypsoSlopeDiagnosticsAlg.cxx'
    source = Path(manifest['source'])
    assert (source/unit).read_bytes() == (root/Path(unit).name).read_bytes()
    assert (root/'software_probe'/Path(unit).name).read_bytes() == (source/unit).read_bytes()
    assert (root/'CalypsoSlopeDiagnosticsAlg_before.cxx').read_bytes() == (
        Path(manifest['source_base'])/Path(unit).name).read_bytes()
    subprocess.run(['git', '-C', str(source), 'apply', '--reverse', '--check', manifest['increment']], check=True)
    acts_header = Path('/cvmfs/atlas.cern.ch/repo/sw/software/24.0/AthenaExternals/24.0.41/InstallArea/x86_64-el9-gcc13-opt/include/Acts/ActsVersion.hpp')
    write_json(output/'source_validation.json', dict(actual_source=str(source/unit),
        compiled_and_frozen_source_exact=True, previous_source_exact=True, increment_reverse_check=True,
        ACTS_version_header=str(acts_header), ACTS_version_header_contents=acts_header.read_text(),
        source_reference='https://raw.githubusercontent.com/acts-project/acts/v32.0.2/Core/src/Propagator/detail/PointwiseMaterialInteraction.cpp',
        source_reference_scope='Version-matched public implementation; runtime exact ACTS commit is not proven by the installed header'))
    write_json(output/'resource_summary.json', dict(jobs=jobs,
        worker_cpu_seconds=sum(j['user_cpu_s']+j['system_cpu_s'] for j in jobs),
        max_RSS_KiB=max(j['max_RSS_KiB'] for j in jobs),
        native_Athena_CPU_separately_measured=False))


if __name__ == '__main__':
    main()
