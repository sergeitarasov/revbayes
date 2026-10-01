"""Summarize this single-dataset recovery pilot; never equate coverage with calibration."""
import csv,json,math,statistics
from pathlib import Path
BASE=Path(__file__).resolve().parent
PARAMS=['lambda0','alpha','mu','psi']
def quantile(x,p):
 x=sorted(x);a=(len(x)-1)*p;i=int(a);return x[i]+(a-i)*(x[min(i+1,len(x)-1)]-x[i])
def ess_batch(x):
 # Consistent batch-means estimate; approximate, not rank-normalized diagnostics.
 n=len(x);b=int(math.sqrt(n));a=n//b
 if a<3:return None
 v=statistics.variance(x)
 batches=[statistics.mean(x[i*b:(i+1)*b]) for i in range(a)]
 long_var=b*statistics.variance(batches)
 return min(n,n*v/long_var) if long_var>0 else None
def summarize():
 truth=json.loads((BASE/'data/manifest.json').read_text())['parameters'];chains={};vectors={}
 for name in ['fixed_tree','joint_1','joint_2']:
  file=BASE/'output'/(name+'.log')
  if not file.exists():continue
  try: rows=list(csv.DictReader(file.read_text().splitlines(),delimiter='\t'))
  except (ValueError,csv.Error):continue
  valid=[r for r in rows if all(r.get(p) is not None for p in PARAMS)]
  if len(valid)<4:continue
  # Burn-in was run separately, but conservatively discard another 20% of samples.
  valid=valid[len(valid)//5:];vectors[name]={p:[float(r[p]) for r in valid] for p in PARAMS}
  chains[name]={p:dict(truth=truth[p],mean=statistics.mean(x),median=statistics.median(x),ci95=[quantile(x,.025),quantile(x,.975)],truth_in_ci95=quantile(x,.025)<=truth[p]<=quantile(x,.975),approx_batch_ess=ess_batch(x)) for p,x in vectors[name].items()}
  chains[name]['retained_samples']=len(valid)
 rhat={}
 if all(k in vectors for k in ['joint_1','joint_2']):
  for p in PARAMS:
   n=min(len(vectors[k][p])//2 for k in ['joint_1','joint_2'])
   if n<3:continue
   splits=[part for k in ['joint_1','joint_2'] for part in [vectors[k][p][:n],vectors[k][p][-n:]]]
   W=statistics.mean(statistics.variance(x) for x in splits);B=n*statistics.variance(statistics.mean(x) for x in splits)
   rhat[p]=math.sqrt(((n-1)*W/n+B/n)/W) if W>0 else None
 result=dict(scope='Single size-selected simulated dataset; known ages/status/origin/rho/clock; not SBC',chains=chains,split_Rhat_joint=rhat,diagnostics='Classical split Rhat and approximate batch-means ESS; assess before interpreting recovery',additional_discard_fraction=.2)
 temp=BASE/'recovery_summary.tmp';temp.write_text(json.dumps(result,indent=2)+'\n');temp.replace(BASE/'recovery_summary.json')
 return result
if __name__=='__main__':print(json.dumps(summarize(),indent=2))
