#!/usr/bin/env python3
import datetime,hashlib,json,os,signal,subprocess,time
from pathlib import Path
from summarize import summarize
BASE=Path(__file__).resolve().parent
RB=BASE.parents[2]/'.local-build/revbayes-build/rb'
os.chdir(BASE)
children={};state=dict(supervisor_pid=os.getpid(),started=datetime.datetime.now().astimezone().isoformat(),status='running',binary=str(RB),binary_sha256=hashlib.sha256(RB.read_bytes()).hexdigest(),jobs={})
state['input_sha256']={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [BASE/'model.Rev',BASE/'fixed_tree.Rev',BASE/'joint_1.Rev',BASE/'joint_2.Rev',BASE/'data/manifest.json',BASE/'data/morphology.nex']}
def save():
 state['updated']=datetime.datetime.now().astimezone().isoformat()
 for name,p in children.items():
  row=state['jobs'][name];row['returncode']=p.poll()
  path=BASE/'output'/(name+'.log')
  if path.exists():
   lines=path.read_text().splitlines()
   if len(lines)>1:row['last_saved_generation']=lines[-1].split('\t')[0]
 tmp=BASE/'status.tmp';tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(BASE/'status.json')
def stop(sig,frame):
 for p in children.values():
  if p.poll() is None:p.terminate()
 state['status']='stopped';save();raise SystemExit(0)
signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
try:
 for name in ['fixed_tree','joint_1','joint_2']:
  with open(BASE/(name+'.stdout'),'w') as log:
   p=subprocess.Popen([str(RB),name+'.Rev'],stdout=log,stderr=subprocess.STDOUT,cwd=BASE)
  children[name]=p;state['jobs'][name]=dict(pid=p.pid,script=name+'.Rev',stdout=name+'.stdout')
 save()
 while any(p.poll() is None for p in children.values()):
  time.sleep(10);save();summarize()
 failures=[]
 for name,p in children.items():
  text=(BASE/(name+'.stdout')).read_text()
  if p.returncode or 'RECOVERY_CHAIN_COMPLETE' not in text:failures.append(name)
 state['status']='failed' if failures else 'complete';state['failures']=failures
 summarize();save()
except Exception as error:
 state['status']='supervisor_error';state['error']=repr(error);save()
 raise
