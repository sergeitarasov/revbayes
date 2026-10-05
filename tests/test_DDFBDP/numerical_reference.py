#!/usr/bin/env python3
"""Standalone tests: analytic constant-rate FBD and independent DD RK4.
No third-party Python libraries required. Does not exercise Rev tree extraction.
"""
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('old_reference', ROOT/'tests/test_DDFBD/numerical_reference.py')
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)
CASES = {
 'ancestor': ([(2,'A')],1),
 'terminal': ([(2.5,'B'),(1,'T')],1),
 'two_ancestors': ([(2,'A'),(1,'A')],1),
 'mixed': ([(2.5,'B'),(2,'A'),(1,'T')],1),
 'all_fossil': ([(2.5,'B'),(1.5,'T'),(.5,'T')],0),
 'root_ancestor': ([(2,'A'),(1.5,'B'),(.5,'T')],1),
 'extant_only': ([(2.5,'B'),(1.25,'B')],3),
}
def analytic(events,n,p):
 l,a,m,s,r=p
 out=ref.constant_functions(3,p)[1]+ref.lp(r,n)
 for t,kind in events:
  e,q=ref.constant_functions(t,p)
  out+= math.log(l)+q if kind=='B' else math.log(s)+(math.log(e)-q if kind=='T' else 0)
 return out

def rk4(events,n,p,cap=48,dt=.0005):
 l,a,m,s,r=p
 v=[1.]+[0.]*cap
 k=1
 previous=3.
 def advance(duration):
  nonlocal v
  lam=[l*math.exp(-a*(k+h-1)) for h in range(cap+1)]
  diag=[-(k+h)*(lam[h]+m+s) for h in range(cap+1)]
  up=[(h+2*k)*lam[h] for h in range(cap+1)]
  down=[h*m for h in range(cap+1)]
  def deriv(w):
   return [diag[h]*w[h]+(up[h-1]*w[h-1] if h else 0)+(down[h+1]*w[h+1] if h<cap else 0) for h in range(cap+1)]
  steps=max(1,math.ceil(duration/dt)); step=duration/steps
  for _ in range(steps):
   d1=deriv(v); d2=deriv([x+step*y/2 for x,y in zip(v,d1)])
   d3=deriv([x+step*y/2 for x,y in zip(v,d2)])
   d4=deriv([x+step*y for x,y in zip(v,d3)])
   v=[x+step*(b+2*c+2*d+e)/6 for x,b,c,d,e in zip(v,d1,d2,d3,d4)]
 for t,kind in events:
  advance(previous-t); previous=t
  if kind=='B':
   v=[x*l*math.exp(-a*(k+h-1)) for h,x in enumerate(v)]; k+=1
  elif kind=='A': v=[x*s for x in v]
  else: v=[0.]+[s*x for x in v[:-1]]; k-=1
 advance(previous)
 assert k==n
 return math.log(sum(x*(1-r)**h for h,x in enumerate(v)))+ref.lp(r,n)

def main():
 results=[]
 with tempfile.TemporaryDirectory(prefix='ddfbdp-kernel-') as tmp:
  exe=Path(tmp)/'kernel'
  kernel=ROOT/'src/core/functions/phylogenetics'
  subprocess.run([os.environ.get('CXX','clang++'),'-std=c++17','-O2','-Wall','-Wextra','-I'+str(kernel),str(ROOT/'tests/test_DDFBDP/kernel_driver.cpp'),str(kernel/'DiversityDependentFbdLikelihood.cpp'),'-o',str(exe)],check=True)
  def run(p,events,n,cap=128):
   cmd=[str(exe),*map(str,p),'3',str(cap),'1e-12',str(n)]
   for age,kind in events: cmd += [str(age),kind]
   return float(subprocess.check_output(cmd,text=True))
  for rho in (.4,1.):
   p=(.4,0.,.15,.2,rho)
   for name,(events,n) in CASES.items():
    got=run(p,events,n); expected=analytic(events,n,p)
    error=abs(got-expected)
    assert error<2e-9,(name,got,expected)
    results.append(dict(test='constant_'+name,rho=rho,error=error))
  for name,(events,n) in CASES.items():
   p=(.4,.12,.15,.2,.7)
   got=run(p,events,n,64); expected=rk4(events,n,p)
   error=abs(got-expected)
   assert error<2e-8,(name,got,expected)
   assert abs(run(p,events,n,128)-got)<2e-9
   results.append(dict(test='dd_'+name,error=error))
  # With complete living sampling and no extinction, a terminal fossil is impossible.
  assert run((.4,.1,0.,.2,1.),CASES['terminal'][0],1)==-math.inf
  # Fossil samples can still be ancestors in the same no-extinction process.
  assert math.isfinite(run((.4,.1,0.,.2,1.),CASES['ancestor'][0],1))
  # No birth: a chain of ancestors has an elementary joint density.
  p=(0.,.2,.15,.2,.7)
  got=run(p,CASES['two_ancestors'][0],1)
  expected=-3*(.15+.2)+2*math.log(.2)+math.log(.7)
  assert abs(got-expected)<2e-9
  # Complete pure birth: there can be no hidden lineage at the present.
  p=(.4,.2,0.,0.,1.)
  got=run(p,CASES['extant_only'][0],3)
  rates=[.4*math.exp(-.2*i) for i in range(3)]
  expected=math.log(rates[0])+math.log(rates[1])-.5*rates[0]-1.25*2*rates[1]-1.25*3*rates[2]
  assert abs(got-expected)<2e-9
  # Positive but extremely rare samples must retain a finite log density.
  got=run((0.,.2,.15,1e-250,1e-200),CASES['two_ancestors'][0],1)
  expected=-3*.15+2*math.log(1e-250)+math.log(1e-200)
  assert abs(got-expected)<2e-9
 print(json.dumps(dict(passed=len(results)+5,comparisons=results),indent=2))
if __name__=='__main__': main()
