#!/usr/bin/env python3
"""Exercise actual Rev wiring and MCMC; does not measure statistical calibration."""
import argparse
import csv
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('rb', type=Path)
    parser.add_argument('--generations', type=int, default=200)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    binary = args.rb.resolve()
    checks = []
    with tempfile.TemporaryDirectory(prefix='ddfbd-integration-') as directory:
        work = Path(directory)
        shutil.copytree(ROOT/'examples/DDFBD/data', work/'data')
        shutil.copy(ROOT/'examples/DDFBD/model.Rev', work/'model.Rev')
        def run(script, name):
            path = work/(name+'.Rev')
            path.write_text(script)
            result = subprocess.run([str(binary), str(path)], cwd=work,
                                    text=True, capture_output=True, timeout=600)
            (work/(name+'.stdout')).write_text(result.stdout+result.stderr)
            if result.returncode or re.search(r'\b(Error|Exception):',result.stdout+result.stderr):
                raise RuntimeError(name+' failed:\n'+(result.stdout+result.stderr)[-12000:])
            return result.stdout+result.stderr
        script = (ROOT/'tests/test_DDFBD/scripts/integration.Rev').read_text()
        text = run(script,'density')
        match = re.search(r'DDFBD_INITIAL_LOG_DENSITY=\s*([-+.\deE]+)',text)
        if not match: raise RuntimeError('Density marker missing:\n'+text)
        actual = float(match[1])
        expected = -23.414925199726387
        # The Rev script checks 1e-8 internally; string conversion can round.
        assert abs(actual-expected) < 1e-4, (actual, expected)
        checks.append(dict(check='Rev_absolute_labelled_density', value=actual, expected=expected))
        for mode,morph,dna in [('combined',True,True),('morphology',True,False),('dna',False,True)]:
            script = ('setOption("debugMCMC", "1", save=false)\n'
                      f'use_morphology = {str(morph).lower()}\nuse_dna = {str(dna).lower()}\n'
                      f'generations = {args.generations}\nprefix = "{mode}"\n'
                      'source("model.Rev")\nprint("DDFBD_MCMC_COMPLETE")\nq()\n')
            text = run(script,mode)
            assert 'DDFBD_MCMC_COMPLETE' in text
            log = work/'output'/(mode+'.log')
            lines = [line for line in log.read_text().splitlines() if not line.startswith('#')]
            rows = list(csv.DictReader(lines, delimiter='\t'))
            assert len(rows) >= args.generations//10
            for row in rows:
                for key in ['Posterior','Likelihood','Prior','lambda0','alpha','mu','psi']:
                    if key in row: assert math.isfinite(float(row[key])), (mode,key,row[key])
            for key in ['lambda0','alpha','mu','psi']:
                assert key in rows[0] and len({row[key] for row in rows}) > 1, (mode,key)
            trees = work/'output'/(mode+'.species.trees')
            assert trees.exists() and trees.stat().st_size > 100
            checks.append(dict(check=mode+'_joint_MCMC', samples=len(rows), generations=args.generations))
        for condition in ['sampling', 'sampledExtant']:
            (work/'conditioned_model.Rev').write_text((work/'model.Rev').read_text().replace(
                'condition="none"', 'condition="'+condition+'"'))
            script = ('setOption("debugMCMC", "1", save=false)\n'
                'use_morphology = true\nuse_dna = true\ngenerations = 100\n'
                f'prefix = "{condition}"\nsource("conditioned_model.Rev")\n'
                'print("DDFBD_CONDITIONED_COMPLETE")\nq()\n')
            text = run(script,condition)
            assert 'DDFBD_CONDITIONED_COMPLETE' in text
            checks.append(dict(check=condition+'_conditioned_joint_MCMC', generations=100))
        checkpoint_script = ('setOption("debugMCMC", "1", save=false)\n'
            'use_morphology = true\nuse_dna = true\ngenerations = 10\nprefix = "checkpoint"\n'
            'source("model.Rev")\n'
            'analysis.run(generations=20, checkpointFile="output/ddfbd.state", checkpointInterval=10)\n'
            'analysis.initializeFromCheckpoint("output/ddfbd.state")\n'
            'analysis.run(generations=20)\nprint("DDFBD_CHECKPOINT_COMPLETE")\nq()\n')
        text = run(checkpoint_script, 'checkpoint')
        assert 'DDFBD_CHECKPOINT_COMPLETE' in text
        checks.append(dict(check='checkpoint_restore_and_MCMC', status='passed'))
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            artifacts = args.output.with_suffix('')
            artifacts.mkdir(exist_ok=True)
            for file in work.glob('*.stdout'): shutil.copy(file,artifacts/file.name)
            shutil.copytree(work/'output',artifacts/'output',dirs_exist_ok=True)
    result = dict(status='passed', scope='runtime integration, not SBC', checks=checks)
    if args.output: args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__ == '__main__': main()
