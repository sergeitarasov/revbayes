#!/usr/bin/env python3
"""Build the shared numerical kernel, compare unmodified DDD, then exercise Rev."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser()
parser.add_argument('--rb', type=Path, default=ROOT/'.local-build/revbayes-build/rb')
args = parser.parse_args()
os.chdir(ROOT)
out = ROOT/'tests/test_DDD_benchmark/results'
out.mkdir(exist_ok=True)
(ROOT/'.local-build').mkdir(exist_ok=True)
compiler = os.environ.get('CXX','c++')
kernel = 'src/core/functions/phylogenetics/DiversityDependentFbdLikelihood.cpp'
for source, target in [('driver.cpp','ddd-kernel'),('test_extant.cpp','ddd-extant-tests')]:
    subprocess.run([compiler,'-std=c++17','-O2','-I','src/core/functions/phylogenetics',
        'tests/test_DDD_benchmark/'+source,kernel,'-o','.local-build/'+target],check=True)
checks = subprocess.run(['.local-build/ddd-extant-tests'],text=True,capture_output=True,check=True)
(out/'extant_checks.txt').write_text(checks.stdout)
subprocess.run(['Rscript','tests/test_DDD_benchmark/benchmark.R'],check=True)
rev = subprocess.run([str(args.rb.resolve()),'tests/test_DDD_benchmark/results/compare.Rev'],
                     text=True,capture_output=True,timeout=300)
(out/'Rev_checks.txt').write_text(rev.stdout+rev.stderr)
if rev.returncode or 'DDD_REV_CHECKS_PASSED' not in rev.stdout or re.search(r'\b(Error|Exception):',rev.stdout+rev.stderr):
    raise RuntimeError('Rev benchmark failed; see results/Rev_checks.txt')
subprocess.run(['python3','tests/test_DDFBD/numerical_reference.py',
                '--output',str(out/'fbd_regression.json')],check=True,stdout=subprocess.DEVNULL)
report = dict(status='passed', kernel_sha256=hashlib.sha256(Path(kernel).read_bytes()).hexdigest(),
    rb_sha256=hashlib.sha256(args.rb.resolve().read_bytes()).hexdigest(),
    extant_checks=checks.stdout.strip(), rev_checks=rev.stdout.strip(),
    scope='Complete extant trees; age, survival or count conditioning; fixed-mu/K one-parameter fit; not statistical recovery')
(out/'run.json').write_text(json.dumps(report,indent=2)+'\n')
print(checks.stdout.strip())
print(rev.stdout.strip())
