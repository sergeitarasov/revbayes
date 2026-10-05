#!/usr/bin/env python3
"""Native BDSTP comparisons and morphology-only MCMC in an isolated directory.
Run: python3 tests/test_DDFBDP/run_integration.py path/to/rb --output result.json
"""
import argparse
import csv
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
ROOT=Path(__file__).resolve().parents[2]

def main():
 parser=argparse.ArgumentParser()
 parser.add_argument('rb',type=Path)
 parser.add_argument('--generations',type=int,default=200)
 parser.add_argument('--output',type=Path)
 args=parser.parse_args()
 binary=args.rb.resolve()
 checks=[]
 with tempfile.TemporaryDirectory(prefix='ddfbdp-integration-') as tmp:
  work=Path(tmp)
  def run(script,name,marker):
   for prefix in ('examples/DDFBDP/data/','tests/test_DDFBDP/data/'):
    script=script.replace('"'+prefix,'"'+str(ROOT/prefix)+'/')
   path=work/(name+'.Rev'); path.write_text(script)
   result=subprocess.run([str(binary),str(path)],cwd=work,text=True,capture_output=True,timeout=1800)
   text=result.stdout+result.stderr
   published_text=text.replace(str(ROOT),'<repository>').replace(str(work),'<temporary-directory>')
   (work/(name+'.stdout')).write_text(published_text)
   if args.output:
    artifacts=args.output.resolve().with_suffix(''); artifacts.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(work/(name+'.stdout'),artifacts/(name+'.stdout'))
   if result.returncode or marker not in text or re.search(r'\b(Error|Exception):',text):
    raise RuntimeError(name+' failed:\n'+text[-16000:])
   return text
  text=run((ROOT/'tests/test_DDFBDP/native_constant.Rev').read_text(),'native','DDFBDP_NATIVE_PASS')
  lines=[line for line in text.splitlines() if line.strip().startswith('NATIVE_COMPARE')]
  assert len(lines)==46,len(lines)
  errors=[float(line.split()[-1]) for line in lines]
  assert all(math.isfinite(float(x)) for line in lines for x in line.split()[-3:])
  assert max(errors)<1e-8
  checks.append(dict(check='native_BDSTP_alpha_zero',comparisons=len(lines),direct_matches=39,survival_rho_corrections=7,max_log_density_error=max(errors),note='Native survival divides by P(sampled extant)/rho; new prior divides by P(sampled extant). Compare native minus ln(rho) for survival.'))
  script=(ROOT/'examples/DDFBDP/morphology.Rev').read_text().replace('generations=1000',f'generations={args.generations}')
  script='setOption("debugMCMC","1",save=false)\n'+script
  text=run(script,'morphology','DDFBDP_MORPHOLOGY_PASS')
  lines=[x for x in (work/'output/DDFBDP/morphology.log').read_text().splitlines() if x and not x.startswith('#')]
  rows=list(csv.DictReader(lines,delimiter='\t'))
  assert len(rows)>=args.generations//10
  keys=['Posterior','Likelihood','Prior','lambda0','alpha','mu','psi','origin','shape','clock','numAncestors']
  for row in rows:
   for key in keys: assert math.isfinite(float(row[key])),(key,row[key])
  for key in keys:
   assert len({row[key] for row in rows})>1,('unchanged',key)
  tree_lines=[x for x in (work/'output/DDFBDP/morphology.trees').read_text().splitlines() if x and not x.startswith('#')]
  trees=list(csv.DictReader(tree_lines,delimiter='\t'))
  if args.output:
   artifacts=args.output.resolve().with_suffix('')
   for filename in ('morphology.log','morphology.trees'):
    shutil.copyfile(work/'output/DDFBDP'/filename,artifacts/filename)
  assert len(trees)==len(rows)
  assert len({x['tree'] for x in trees})>1
  checks.append(dict(check='morphology_joint_MCMC_debug_cache_checks',generations=args.generations,samples=len(rows),ancestor_counts=sorted({int(float(x['numAncestors'])) for x in rows})))
  # Re-evaluate sampled states with larger hidden-count cutoffs.
  cap_script=['taxa <- readTaxonData("examples/DDFBDP/data/bears_taxa.tsv")']
  for j in sorted({0,len(rows)//4,len(rows)//2,3*len(rows)//4,len(rows)-1}):
   row=rows[j]
   treepath=work/f'state{j}.tre'; treepath.write_text(trees[j]['tree']+'\n')
   cap_script += [f'initial <- readTrees("{treepath}")[1]']
   for cap in (64,128,256):
    cap_script += [f't{cap} ~ dnDDFBDP(originAge={row["origin"]}, lambda0={row["lambda0"]}, alpha={row["alpha"]}, mu={row["mu"]}, psi={row["psi"]}, rho=1, taxa=taxa, initialTree=initial, condition="sampling", maxHiddenLineages={cap})']
   cap_script += ['delta <- abs(t128.lnProbability()-t256.lnProbability())',f'print("CAP_COMPARE {j} " + abs(t64.lnProbability()-t128.lnProbability()) + " " + delta)', 'if (delta>1e-7) { stop("Hidden-count cutoff sensitivity") }']
  cap_script+=['print("DDFBDP_CAP_PASS")','q()']
  captext=run('\n'.join(cap_script),'cutoff','DDFBDP_CAP_PASS')
  caplines=[x.strip() for x in captext.splitlines() if x.strip().startswith('CAP_COMPARE')]
  assert len(caplines)==5
  checks.append(dict(check='sampled_state_cutoff_sensitivity',comparisons=caplines))
  # Check constrained-topology delegation and checkpoint deserialization with SA trees.
  short=script.replace(f'generations={args.generations}','generations=20')
  short=short.replace('tree ~ dnDDFBDP(', 'basePrior = dnDDFBDP(')
  short=short.replace('taxa=taxa,condition="sampling",maxHiddenLineages=128)', 'taxa=taxa,condition="sampling",maxHiddenLineages=128)\ntree ~ dnConstrainedTopology(basePrior, constraints=v(clade("Ursus_arctos","Ursus_maritimus")))')
  short=short.replace('q()', 'analysis.run(generations=20, checkpointFile="output/DDFBDP/state", checkpointInterval=10)\nanalysis.initializeFromCheckpoint("output/DDFBDP/state")\nanalysis.run(generations=20)\nprint("DDFBDP_CHECKPOINT_PASS")\nq()')
  run(short,'checkpoint','DDFBDP_CHECKPOINT_PASS')
  checks.append(dict(check='constrained_topology_checkpoint_restore',status='passed'))
 result=dict(status='passed',scope='Numerical and runtime integration; not parameter recovery or posterior convergence',checks=checks)
 if args.output:
  args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__': main()
