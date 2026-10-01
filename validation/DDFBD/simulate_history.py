#!/usr/bin/env python3
"""Independent complete-history budding DD-FBD simulator and pruning oracle.

Only the Python standard library is required. Simulation uses Gillespie total
population hazards, not the numerical likelihood filter. Every species keeps
its identity across births. Hidden species and all extinct endpoints remain in
history.json; pruning retains endpoints exactly for sampled named species.

Example:
 python3 validation/DDFBD/simulate_history.py --seed 7 --output /tmp/ddfbd-sim
 python3 validation/DDFBD/simulate_history.py --self-test

This is a simulator/verification artifact, not an MCMC implementation or SBC.
"""
import argparse
import csv
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import random


@dataclass(frozen=True)
class Parameters:
    lambda0: float = 0.8
    alpha: float = 0.15
    mu: float = 0.25
    psi: float = 0.7
    rho: float = 0.8

    def check(self):
        values=asdict(self)
        if any(not math.isfinite(x) or x < 0 for x in values.values()) or self.rho > 1:
            raise ValueError('rates must be finite nonnegative; rho must lie in [0,1]')

    def birth(self, n):
        return self.lambda0*math.exp(-self.alpha*(n-1)) if n else 0.0


def simulate(p, origin, seed, max_events=100000):
    p.check()
    if not math.isfinite(origin) or origin <= 0:
        raise ValueError('origin must be finite and positive')
    rng=random.Random(seed)
    alive=[1]
    next_id=2
    events=[]
    elapsed=0.0
    while alive:
        n=len(alive)
        per_species=p.birth(n)+p.mu+p.psi
        if per_species == 0:
            break
        event_time=elapsed+rng.expovariate(n*per_species)
        if event_time >= origin:
            break
        species=rng.choice(alive)
        choice=rng.random()*per_species
        if choice < p.birth(n):
            event=dict(time=event_time,type='birth',species=species,child=next_id)
            alive.append(next_id)
            next_id += 1
        elif choice < p.birth(n)+p.mu:
            event=dict(time=event_time,type='death',species=species)
            alive.remove(species)
        else:
            event=dict(time=event_time,type='fossil',species=species)
        event['total_before']=n
        events.append(event)
        elapsed=event_time
        if len(events) >= max_events:
            raise RuntimeError('simulation event safety limit reached; history not returned')
    sampled=[species for species in alive if rng.random() < p.rho]
    history=dict(parameters=asdict(p),origin=origin,seed=seed,events=events,
                 sampled_present=sampled,simulation_condition='unconditional')
    history['complete_history_log_density']=history_log_density(history)
    return history


def log_probability(probability, count=1):
    if not count:
        return 0.0
    return count*math.log(probability) if probability else -math.inf


def history_log_density(history):
    """Identity-specific event measure: no N factors in event-density terms.

    Reconstruct the living identities rather than trusting stored N or lifetimes.
    Fossil records are ordered events of the observed Poisson process.
    """
    p=Parameters(**history['parameters'])
    p.check()
    origin=history['origin']
    alive={1}
    ever={1}
    previous=0.0
    result=0.0
    for event in history['events']:
        t=event['time']
        if not previous <= t < origin or event['species'] not in alive:
            raise ValueError('incompatible event time or species identity')
        n=len(alive)
        if 'total_before' in event and event['total_before'] != n:
            raise ValueError('stored population does not match reconstructed population')
        result -= (t-previous)*n*(p.birth(n)+p.mu+p.psi)
        if event['type'] == 'birth':
            if event['child'] in ever:
                raise ValueError('birth reuses a species identity')
            result += log_probability(p.lambda0)-p.alpha*(n-1)
            ever.add(event['child'])
            alive.add(event['child'])
        elif event['type'] == 'death':
            result += log_probability(p.mu)
            alive.remove(event['species'])
        elif event['type'] == 'fossil':
            result += log_probability(p.psi)
        else:
            raise ValueError('unknown history event')
        previous=t
    n=len(alive)
    result -= (origin-previous)*n*(p.birth(n)+p.mu+p.psi)
    sampled=set(history['sampled_present'])
    if len(sampled) != len(history['sampled_present']) or not sampled <= alive:
        raise ValueError('present sample must contain unique living species')
    result += log_probability(p.rho,len(sampled))+log_probability(1-p.rho,n-len(sampled))
    return result


def species_histories(history):
    origin=history['origin']
    species={1:dict(birth=0.0,parent=None,death=origin,extinct=False,fossils=[],buds=[])}
    for event in history['events']:
        record=species[event['species']]
        if event['type'] == 'birth':
            child=event['child']
            record['buds'].append((event['time'],child))
            species[child]=dict(birth=event['time'],parent=event['species'],death=origin,
                                extinct=False,fossils=[],buds=[])
        elif event['type'] == 'death':
            record['death']=event['time']
            record['extinct']=True
        else:
            record['fossils'].append(event['time'])
    sampled=set(history['sampled_present'])
    for identity, record in species.items():
        record['sampled_present']=identity in sampled
        record['observed']=bool(record['fossils']) or identity in sampled
    return species


def prune_extended(history):
    """Return oriented extended tree; A child is ALWAYS first, D child second.

    Unary suppression may follow an unobserved new species before the first
    observed occurrence. It never discards an observed species' true endpoint.
    Node ages and structural attachment ages are returned explicitly.
    """
    history_log_density(history)  # Validate identities independently of pruning.
    species=species_histories(history)
    origin=history['origin']
    def lineage(identity, position=0):
        record=species[identity]
        if position == len(record['buds']):
            if not record['observed']:
                return None
            return dict(age=origin-record['death'],species=identity)
        time,child=record['buds'][position]
        a=lineage(identity,position+1)
        d=lineage(child)
        if a is None:
            return d
        if d is None:
            return a
        return dict(age=origin-time,children=[a,d])
    tree=lineage(1)
    observed={identity:record for identity,record in species.items() if record['observed']}
    attachments={}
    birth_events=[]
    def assign(node, attachment):
        if 'species' in node:
            attachments[node['species']]=attachment
            return
        birth_events.append((node['age'],'Birth'))
        assign(node['children'][0],attachment)
        assign(node['children'][1],node['age'])
    if tree:
        assign(tree,origin)
    rows=[]
    events=list(birth_events)
    for identity,record in observed.items():
        oldest=origin-min(record['fossils']) if record['fossils'] else 0.0
        endpoint=origin-record['death']
        row=dict(species_id='sp'+str(identity),species_number=identity,
                 structural_attachment_age=attachments[identity],oldest_observation_age=oldest,
                 endpoint_age=endpoint,true_species_birth_age=origin-record['birth'],
                 extinct=record['extinct'],sampled_present=record['sampled_present'],
                 fossil_count=len(record['fossils']))
        if not endpoint <= oldest <= attachments[identity]:
            raise AssertionError('pruned identity interval incompatible with retained path')
        rows.append(row)
        events.append((oldest,'Protect'))
        if record['extinct']:
            events.append((endpoint,'Death'))
    order={'Birth':0,'Protect':1,'Death':2}
    events.sort(key=lambda e:(-e[0],order[e[1]]))
    return dict(tree=tree,species=rows,kernel_events=events,
                fossilRecords=sum(row['fossil_count'] for row in rows),
                sampledLiving=sum(row['sampled_present'] for row in rows),
                unsampledLiving=sum(not row['extinct'] and not row['sampled_present'] for row in rows))


def newick(tree, origin, annotated=True):
    """Include an origin stem length; JSON node ages remain authoritative."""
    if tree is None:
        return None
    def encode(node,parent_age):
        if 'species' in node:
            body='sp'+str(node['species'])
        else:
            body='('+','.join(encode(child,node['age']) for child in node['children'])+')'
        if annotated:
            body+='[&age='+format(node['age'],'.17g')+']'
        return body+':'+format(parent_age-node['age'],'.17g')
    return encode(tree,origin)+';'


def write_dataset(history, output):
    output.mkdir(parents=True,exist_ok=True)
    projected=prune_extended(history)
    (output/'history.json').write_text(json.dumps(history,indent=2,allow_nan=False)+'\n')
    (output/'extended_tree.json').write_text(json.dumps(projected,indent=2)+'\n')
    for name,annotated in [('extended_tree.nwk',False),('extended_tree_annotated.nwk',True)]:
        text=newick(projected['tree'],history['origin'],annotated)
        if text:
            (output/name).write_text(text+'\n')
    occurrence_rows=[]
    for number,event in enumerate((e for e in history['events'] if e['type']=='fossil'),1):
        age=history['origin']-event['time']
        occurrence_rows.append(dict(occurrence_id='occ'+str(number),species_id='sp'+str(event['species']),
                                    age_min=age,age_max=age,count=1))
    for name,rows,fields in [
        ('occurrences.tsv',occurrence_rows,['occurrence_id','species_id','age_min','age_max','count']),
        ('species.tsv',projected['species'],['species_id','species_number','structural_attachment_age',
          'oldest_observation_age','endpoint_age','true_species_birth_age','extinct','sampled_present','fossil_count'])]:
        with (output/name).open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
    manifest=dict(seed=history.get('seed'),origin=history['origin'],parameters=history['parameters'],
                  selection=history.get('simulation_condition','explicit fixture'),
                  observed_species=len(projected['species']),fossils=projected['fossilRecords'],
                  complete_history_log_density=history_log_density(history),
                  no_observations=projected['tree'] is None,
                  newick_notes='A child first, D child second. Root branch is the origin stem. Use JSON node ages/endpoint ages; a generic fossil-only Newick importer may translate ages.')
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest


def self_test():
    p=Parameters(0.5,0.4,0.2,0.3,0.7)
    base=dict(parameters=asdict(p),origin=2.0,sampled_present=[1],events=[])
    assert abs(history_log_density(base)-(-2+math.log(.7))) < 1e-12
    extinct=dict(parameters=asdict(p),origin=2.0,sampled_present=[],events=[
        dict(time=.2,type='fossil',species=1),dict(time=.9,type='fossil',species=1),
        dict(time=1.5,type='death',species=1)])
    exact=-1.5+2*math.log(.3)+math.log(.2)
    assert abs(history_log_density(extinct)-exact) < 1e-12
    e=prune_extended(extinct)
    assert e['tree']==dict(age=.5,species=1) and e['fossilRecords']==2
    assert e['kernel_events']==[(1.8,'Protect'),(.5,'Death')]
    history=dict(parameters=asdict(p),origin=2.0,sampled_present=[2],events=[
        dict(time=.2,type='fossil',species=1),dict(time=.5,type='birth',species=1,child=2),
        dict(time=.9,type='fossil',species=1),dict(time=1.2,type='death',species=1)])
    exact=-(.5*1+.7*2*(p.birth(2)+p.mu+p.psi)+.8*1)+math.log(.5)+2*math.log(.3)+math.log(.2)+math.log(.7)
    assert abs(history_log_density(history)-exact) < 1e-12
    e=prune_extended(history)
    assert e['tree']['children'][0]['species']==1 and e['tree']['children'][1]['species']==2
    assert e['tree']['children'][0]['age']==.8
    assert e['species'][0]['oldest_observation_age']==1.8
    # The same genealogy with no parent fossils must suppress its endpoint and
    # the corresponding split while retaining the child species endpoint.
    unobserved_parent={**history,'events':[event for event in history['events'] if event['type']!='fossil']}
    pruned=prune_extended(unobserved_parent)
    assert pruned['tree']==dict(age=0.0,species=2)
    assert pruned['species'][0]['structural_attachment_age']==2.0
    assert pruned['species'][0]['true_species_birth_age']==1.5
    checks=['no_event_density','repeated_fossils_extinction_density','budding_identity_density',
            'ancestral_species_endpoint_preserved','unobserved_parent_suppression']
    # Across independent realizations, validate pruning and exact species
    # assignments. This is a simulator invariant test, not statistical recovery.
    for seed in range(30):
        draw=simulate(Parameters(),3.0,seed)
        result=prune_extended(draw)
        expected={species for species,record in species_histories(draw).items() if record['observed']}
        assert {r['species_number'] for r in result['species']}==expected
        assert result['fossilRecords']==sum(e['type']=='fossil' for e in draw['events'])
    checks.append('30_gillespie_pruning_invariants')
    return dict(status='passed',checks=checks,scope='simulator and pruning invariants, not SBC')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed',type=int,default=7)
    parser.add_argument('--origin',type=float,default=3)
    parser.add_argument('--lambda0',type=float,default=.8)
    parser.add_argument('--alpha',type=float,default=.15)
    parser.add_argument('--mu',type=float,default=.25)
    parser.add_argument('--psi',type=float,default=.7)
    parser.add_argument('--rho',type=float,default=.8)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--self-test',action='store_true')
    args=parser.parse_args()
    if args.self_test:
        print(json.dumps(self_test(),indent=2))
        return
    if args.output is None:
        parser.error('--output is required unless --self-test is used')
    history=simulate(Parameters(args.lambda0,args.alpha,args.mu,args.psi,args.rho),args.origin,args.seed)
    print(json.dumps(write_dataset(history,args.output),indent=2))

if __name__=='__main__':
    main()
