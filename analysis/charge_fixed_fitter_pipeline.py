#!/usr/bin/env python3
"""Run dependent native-output gates and reductions in one Condor worker."""
import argparse
from pathlib import Path

from charge_posterior import gate, analyze
from charge_fixed_fitter import analyze as compare_solver_coverage


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--analysis-version', default='v1')
    args = parser.parse_args()
    root = args.input
    gate_output = root/('default_gate_'+args.analysis_version)
    analysis_output = root/('analysis_'+args.analysis_version)
    gate(root, gate_output)
    analyze(root, analysis_output, gate_output/'validation.json')
    compare_solver_coverage(root, analysis_output)


if __name__ == '__main__':
    main()
