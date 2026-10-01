#!/usr/bin/env python3
"""Link and run the FBDSP dispatch test against built RevBayes core libraries."""
import argparse
import os
from pathlib import Path
import shlex
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--build-dir', type=Path, required=True)
parser.add_argument('--boost-root', type=Path, required=True)
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
include_dirs = sorted({p.parent for p in (root / 'src').rglob('*.h')})
exe = args.build_dir.resolve() / 'test_fbdsp_dispatch'
command = shlex.split(os.environ.get('CXX', 'clang++')) + ['-std=c++23']
command += ['-I' + str(p) for p in include_dirs]
command += ['-I' + str(args.boost_root.resolve() / 'include')]
command += [str(Path(__file__).with_name('test_dispatch.cpp'))]
command += [str(args.build_dir.resolve() / lib) for lib in ['librb-core.a', 'librb-revlanguage.a', 'librb-lib.a']]
command += [str(args.boost_root.resolve() / 'lib' / ('libboost_' + lib + '.a'))
            for lib in ['thread', 'serialization', 'date_time']]
command += ['-o', str(exe)]
subprocess.run(command, check=True)
subprocess.run([str(exe)], check=True)
