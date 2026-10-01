#!/usr/bin/env python3
"""Audit saved synthetic MCMC draws and recompute tree density at larger cutoffs."""
import argparse, csv, json, math, re, subprocess, tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]

def read_table(path):
    return list(csv.DictReader((x for x in path.read_text().splitlines() if not x.startswith('#')),delimiter='\t'))

def parse_tree(text):
    text=re.sub(r'\[[^\]]*\]','',text)
    tokens=re.findall(r'[(),:;]|[^(),:;\s]+',text)
    pos=0
    def node():
        nonlocal pos
        children=[]; name=''
        if tokens[pos]=='(':
            pos+=1; children.append(node())
            while tokens[pos]==',': pos+=1; children.append(node())
            assert tokens[pos]==')'; pos+=1
        else: name=tokens[pos]; pos+=1
        length=0.
        if tokens[pos]==':': pos+=1; length=float(tokens[pos]); pos+=1
        return dict(name=name,length=length,children=children)
    root=node()
    tips={}
    def distances(n,d):
        n['distance']=d
        if n['children']:
            for c in n['children']: distances(c,d+c['length'])
        else: tips[n['name']]=n
    distances(root,0.)
    # A and C are observed extant species; B is extinct.
    assert set(tips)=={'A','B','C'}
    assert abs(tips['A']['distance']-tips['C']['distance'])<5e-6
    height=(tips['A']['distance']+tips['C']['distance'])/2
    def ages(n):
        n['age']=max(0.,height-n['distance'])
        for c in n['children']: ages(c)
    ages(root)
    def topology(n): return '('+','.join(map(topology,n['children']))+')' if n['children'] else n['name']
    return root,tips,topology(root)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('artifacts',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    reports=[]
    with tempfile.TemporaryDirectory(prefix='ddfbd-draw-audit-') as tmp:
        driver=Path(tmp)/'kernel'
        subprocess.run(['c++','-std=c++17','-O2','-I'+str(ROOT/'src/core/functions/phylogenetics'),
            str(ROOT/'tests/test_DDFBD/kernel_driver.cpp'),str(ROOT/'src/core/functions/phylogenetics/DiversityDependentFbdLikelihood.cpp'),'-o',str(driver)],check=True)
        for mode in ['combined','morphology','dna']:
            rows=read_table(args.artifacts/'output'/(mode+'.log'))
            trees=read_table(args.artifacts/'output'/(mode+'.species.trees'))
            assert len(rows)==len(trees)
            moving=['lambda0','alpha','mu','psi','origin']+['fossil_age[%d]'%i for i in range(1,6)]
            if mode!='dna': moving+=['morph_clock']
            if mode!='morphology': moving+=['dna_clock']
            assert all(len({r[k] for r in rows})>1 for k in moving)
            shapes=set(); deaths=[]; errors=[]
            for i,(row,record) in enumerate(zip(rows,trees)):
                assert row['Iteration']==record['Iteration']
                root,tips,shape=parse_tree(record['species_tree']); shapes.add(shape)
                death=tips['B']['age']; deaths.append(death)
                fossil=[float(row['fossil_age[%d]'%j]) for j in range(1,6)]
                assert 0<death<min(fossil[2:4])+5e-6
                assert root['age']<float(row['origin'])+5e-6
                for j in range(5): assert abs(float(row['sample_ages[%d]'%(j+1)])-fossil[j])<1e-5
                if i % max(1,len(rows)//20): continue
                events=[]; attach={}
                def visit(n):
                    if not n['children']: return n['name']
                    assert len(n['children'])==2
                    a=visit(n['children'][0]); d=visit(n['children'][1])
                    attach[d]=n['age']; events.append((n['age'],'B')); return a
                attach[visit(root)]=float(row['origin'])
                oldest={'A':max(fossil[:2]),'B':max(fossil[2:4]),'C':fossil[4]}
                for name in tips:
                    assert attach[name]+5e-6>=oldest[name]>=tips[name]['age']-5e-6
                    events.append((oldest[name],'P'))
                events.append((death,'D')); events.sort(key=lambda x:(-x[0],'BPD'.index(x[1])))
                ll=[]
                for cap in [64,128,256]:
                    command=[str(driver),'ll',row['lambda0'],row['alpha'],row['mu'],row['psi'],'0.7',row['origin'],str(cap),'1e-13','5','2','0',str(len(events))]
                    for age,event in events: command.extend([str(age),event])
                    ll.append(float(subprocess.check_output(command,text=True)))
                errors.append(dict(iteration=int(row['Iteration']),log_density_raw=ll,
                    delta_64_256=ll[2]-ll[0],delta_128_256=ll[2]-ll[1]))
            assert len(shapes)>1 and len(set(deaths))>1
            maximum=max(abs(x['delta_64_256']) for x in errors)
            assert maximum<1e-6, (mode,maximum)
            reports.append(dict(mode=mode,draws=len(rows),distinct_ordered_topologies=len(shapes),
                extinction_age_range=[min(deaths),max(deaths)],moving_parameters=moving,
                max_cutoff_difference=maximum,cutoff_draws=errors))
    result=dict(status='passed',scope='Synthetic runtime/support and cutoff checks, not convergence or SBC',checks=reports)
    args.output.parent.mkdir(exist_ok=True,parents=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({r['mode']:dict(topologies=r['distinct_ordered_topologies'],max_cutoff_difference=r['max_cutoff_difference']) for r in reports},indent=2))
if __name__=='__main__': main()
