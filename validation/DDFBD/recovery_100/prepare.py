from pathlib import Path
import importlib.util,sys,json,random,math,csv
out=Path(__file__).resolve().parent
if (out/'status.json').exists(): raise RuntimeError('Existing run detected; prepare a separate directory instead of overwriting it')
spec=importlib.util.spec_from_file_location('sim',out.parent/'simulate_history.py'); sim=importlib.util.module_from_spec(spec);sys.modules['sim']=sim;spec.loader.exec_module(sim)
data=out/'data';p=sim.Parameters(.8,.015,.15,.3,.8);attempts=[]
for seed in range(1001,5001):
 h=sim.simulate(p,12,seed); e=sim.prune_extended(h); n=len(e['species']); attempts.append(dict(seed=seed,observed=n))
 if n==100:break
else: raise RuntimeError('No 100-species draw found')
manifest=sim.write_dataset(h,data)
manifest.update(selection='First seed >=1001 yielding exactly 100 observed named species; size-selected pilot, not SBC',attempts=attempts,observed_extinct=sum(s['extinct'] for s in e['species']),sampled_extant=e['sampledLiving'],unsampled_living=e['unsampledLiving'],character_model='100 independent binary symmetric sites, rate 0.2, all sites retained',fixed_in_inference=['origin=12','rho=0.8','character_clock=0.2','exact occurrence ages','living/extinct status'],inferred=['lambda0','alpha','mu','psi','extended tree and extinction endpoints in joint chains'])
(data/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
rng=random.Random(919191);states={1:[rng.randrange(2) for _ in range(100)]};last={1:0.};sequences={};ages={};species={};count=0
# Simulate characters independently on the complete history, including hidden species.
def evolve(i,t):
 prob=-math.expm1(-2*.2*(t-last[i]))/2
 states[i]=[1-x if rng.random()<prob else x for x in states[i]];last[i]=t
for event in h['events']:
 i=event['species'];t=event['time'];evolve(i,t)
 if event['type']=='birth':
  j=event['child'];states[j]=states[i].copy();last[j]=t
 elif event['type']=='fossil':
  count+=1;name='occ'+str(count);sequences[name]=''.join(map(str,states[i]));ages[name]=12-t;species[name]='sp'+str(i)
for i in h['sampled_present']:
 evolve(i,12);name='sp'+str(i)+'_now';sequences[name]=''.join(map(str,states[i]));ages[name]=0;species[name]='sp'+str(i)
(data/'morphology.nex').write_text('#NEXUS\nbegin data;\ndimensions ntax=%d nchar=100;\nformat datatype=standard symbols="01" missing=? gap=-;\nmatrix\n'%len(sequences)+'\n'.join(k+' '+v for k,v in sequences.items())+'\n;\nend;\n')
v=lambda xs:'v('+','.join(json.dumps(x) if isinstance(x,str) else format(x,'.17g') for x in xs)+')'
occ=list(csv.DictReader((data/'occurrences.tsv').open(),delimiter='\t'))
common='''# Synthetic recovery pilot: known origin, rho, character clock, ages and status.
initial_tree <- readTrees("data/extended_tree.nwk")[1]
occ_species <- %s
occ_ages <- %s
sample_names <- %s
sample_species <- %s
sample_ages <- %s
extant <- %s
origin <- 12.0
rho <- 0.8
clock <- 0.2
moves = VectorMoves()
lambda0 ~ dnExponential(1.0)
alpha ~ dnExponential(50.0)
mu ~ dnExponential(3.0)
psi ~ dnExponential(2.0)
lambda0.setValue(start_lambda)
alpha.setValue(start_alpha)
mu.setValue(start_mu)
psi.setValue(start_psi)
moves.append(mvScale(lambda0,lambda=0.5,weight=1))
moves.append(mvScale(alpha,lambda=0.5,weight=1))
moves.append(mvScale(mu,lambda=0.5,weight=1))
moves.append(mvScale(psi,lambda=0.5,weight=1))
species_tree ~ dnDDFBD(originAge=origin,lambda0=lambda0,alpha=alpha,mu=mu,psi=psi,rho=rho,
    occurrenceSpecies=occ_species,occurrenceAges=occ_ages,sampledExtant=extant,
    initialTree=initial_tree,condition="none",maxHiddenLineages=128,numericalTolerance=1e-11)
observation_tree := fnSpeciesObservationTree(species_tree,sampleNames=sample_names,speciesNames=sample_species,ages=sample_ages)
morph <- readDiscreteCharacterData("data/morphology.nex")
characters ~ dnPhyloCTMC(tree=observation_tree,Q=fnJC(2),branchRates=clock,type="Standard",coding="all")
characters.clamp(morph)
endpoint_min <- 0.0
'''%(v([x['species_id'] for x in occ]),v([float(x['age_min']) for x in occ]),v(sequences),v(species.values()),v(ages.values()),v(['sp'+str(i) for i in h['sampled_present']]))
common+='if (joint) {\n    moves.append(mvOrderedTreeSwap(species_tree,weight=2))\n    moves.append(mvNodeTimeSlideUniform(species_tree,weight=3))\n    moves.append(mvRootTimeSlideUniform(species_tree,origin=origin,weight=1))\n}\n'
bounds=[]
for s in e['species']:
 if not s['extinct']:continue
 name=s['species_id'];young=min(float(o['age_min']) for o in occ if o['species_id']==name);bound='youngest_'+name;bounds.append(bound)
 common+=f'{bound} <- {young:.17g}\nif (joint) {{ moves.append(mvTipTimeSlideUniform(species_tree,tip="{name}",min=endpoint_min,max={bound},delta=0.5,weight=0.1)) }}\n'
common+='joint_model = model(species_tree,endpoint_min,'+','.join(bounds)+')\n'
common+='''monitors = VectorMonitors()
monitors.append(mnModel(filename="output/"+prefix+".log",printgen=10))
monitors.append(mnFile(filename="output/"+prefix+".trees",printgen=100,species_tree))
analysis = mcmc(joint_model,monitors,moves)
if (burnin_generations > 0) { analysis.burnin(generations=burnin_generations,tuningInterval=200) }
analysis.run(generations=generations,checkpointFile="output/"+prefix+".state",checkpointInterval=1000)
analysis.operatorSummary()
print("RECOVERY_CHAIN_COMPLETE")
'''
(out/'model.Rev').write_text(common)
for name,joint,seed,starts,burn,gen in [('pilot',True,32001,(.55,.008,.1,.2),0,10),('fixed_tree',False,32002,(.55,.008,.1,.2),2000,20000),('joint_1',True,32003,(.55,.008,.1,.2),5000,30000),('joint_2',True,32004,(1.2,.025,.3,.5),5000,30000)]:
 script=f'seed({seed})\njoint = {str(joint).lower()}\nprefix = "{name}"\nburnin_generations = {burn}\ngenerations = {gen}\n'+''.join(f'start_{k} = {val}\n' for k,val in zip(['lambda','alpha','mu','psi'],starts))+'source("model.Rev")\nq()\n'
 (out/(name+'.Rev')).write_text(script)
print(json.dumps({k:v for k,v in manifest.items() if k!='attempts'},indent=2));print('Attempts:',len(attempts),'character samples',len(sequences))
