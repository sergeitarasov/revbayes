#!/usr/bin/env python3
"""Compile-only FBDSP const-virtual regression against the actual core header.

Usage: python3 tests/test_FBDSP_dispatch/check_override.py --boost-include PATH
Requires a C++23 compiler (CXX, default clang++); no RevBayes link is needed.
"""
import argparse
import os
from pathlib import Path
import shlex
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--boost-include', type=Path)
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
include_dirs = sorted({p.parent for p in (root / 'src').rglob('*.h')})
command = shlex.split(os.environ.get('CXX', 'clang++')) + ['-std=c++23', '-fsyntax-only']
command += ['-I' + str(p) for p in include_dirs]
if args.boost_include:
    command += ['-I' + str(args.boost_include.resolve())]
command += [str(Path(__file__).with_suffix('.cpp'))]
subprocess.run(command, check=True)
print('PASS: FBDSP overrides the const divergence-times likelihood virtual.')
