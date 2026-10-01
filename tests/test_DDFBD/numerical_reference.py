#!/usr/bin/env python3
"""Independent numerical regression of the standalone production DDFBD kernel.

Uses only Python's standard library and clang++ (or CXX). The constant-rate
reference is the unconditioned oriented extended-tree formula of Stadler et al.
2018, Theorem 2. A separately implemented RK4 forward solver checks DD cases.
Neither successful test implies that the Rev wrapper or topology moves work.
Run from anywhere: python3 tests/test_DDFBD/numerical_reference.py
"""
import argparse
import json
import importlib.util
import sys
import math
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
KERNEL = ROOT / 'src/core/functions/phylogenetics'


def schedule(species):
    # Each species has (structural attachment, oldest observation, endpoint,
    # sampled at present). The first species begins at the process origin.
    events = [(b, 'B') for b, o, d, sampled in species[1:]]
    events += [(o, 'P') for b, o, d, sampled in species]
    events += [(d, 'D') for b, o, d, sampled in species if d > 0]
    return sorted(events, key=lambda e: (-e[0], 'BPD'.index(e[1])))


def counts(species):
    sampled = sum(d == 0 and sampled for b, o, d, sampled in species)
    unsampled = sum(d == 0 and not sampled for b, o, d, sampled in species)
    return sampled, unsampled


def constant_functions(t, p):
    """Closed-form p(t), log q(t); handles linear and critical limits."""
    lam, alpha, mu, psi, rho = p
    assert alpha == 0
    gamma = lam+mu+psi
    if lam == 0:
        if gamma == 0:
            return 1-rho, 0.0
        e = math.exp(-gamma*t)
        return mu/gamma + (1-rho-mu/gamma)*e, -gamma*t
    delta = math.sqrt((lam-mu-psi)**2 + 4*lam*psi)
    if delta == 0:
        # Necessarily lambda=mu and psi=0 for nonnegative rates.
        denominator = 1+lam*rho*t
        return 1-rho/denominator, -2*math.log(denominator)
    c2 = (mu+psi-lam+2*lam*rho)/delta
    exponential = math.exp(-delta*t)
    denominator = (1+c2)+(1-c2)*exponential
    if denominator <= 0:
        raise ArithmeticError('singular analytic fixture')
    no_sample = (gamma + delta*((1-c2)*exponential-(1+c2))/denominator)/(2*lam)
    log_q = math.log(4)-delta*t-2*math.log(denominator)
    return no_sample, log_q


def lp(value, n):
    return 0.0 if n == 0 else (-math.inf if value == 0 else n*math.log(value))


def theorem2(species, fossils, p):
    lam, alpha, mu, psi, rho = p
    sampled, unsampled = counts(species)
    result = (lp(lam, len(species)-1)+lp(mu, sum(d > 0 for _, _, d, _ in species))
              +lp(psi, fossils)+lp(rho, sampled)+lp(1-rho, unsampled))
    if result == -math.inf:
        return result
    gamma = lam+mu+psi
    for b, o, d, _ in species:
        qb = constant_functions(b, p)[1]
        qo = constant_functions(o, p)[1]
        qd = constant_functions(d, p)[1]
        result += qb-qo + 0.5*(qo-qd-gamma*(o-d))
    return result


def independent_rk4(origin, events, fossils, sampled, unsampled, p,
                     cap=48, step=0.0005, ignore_protection=False):
    """Direct RK4, not uniformization; suitable only for these small fixtures."""
    lam, alpha, mu, psi, rho = p
    v = [1.0]+[0.0]*cap
    k = 1
    protected = 0
    current = origin
    logscale = lp(psi, fossils)
    def birth(n):
        return lam*math.exp(-alpha*(n-1)) if n else 0.0
    def propagate(duration):
        nonlocal v, logscale
        diag = [-(k+h)*(birth(k+h)+mu+psi) for h in range(cap+1)]
        up = [(h+2*k-(0 if ignore_protection else protected))*birth(k+h)
              for h in range(cap+1)]
        down = [h*mu for h in range(cap+1)]
        def rhs(x):
            y = [d*xx for d, xx in zip(diag, x)]
            for h in range(cap):
                y[h+1] += up[h]*x[h]
                y[h] += down[h+1]*x[h+1]
            return y
        steps = max(1, math.ceil(duration/step))
        dt = duration/steps
        for unused in range(steps):
            a = rhs(v)
            b = rhs([x+dt*z/2 for x,z in zip(v,a)])
            c = rhs([x+dt*z/2 for x,z in zip(v,b)])
            d = rhs([x+dt*z for x,z in zip(v,c)])
            v = [x+dt*(aa+2*bb+2*cc+dd)/6 for x,aa,bb,cc,dd in zip(v,a,b,c,d)]
        scale = max(v)
        assert min(v) >= -1e-14 and scale > 0
        v = [x/scale for x in v]
        logscale += math.log(scale)
    for age, event in events:
        propagate(current-age)
        if event == 'B':
            v = [x*birth(k+h) for h,x in enumerate(v)]
            k += 1
        elif event == 'P':
            protected += 1
        elif event == 'D':
            logscale += lp(mu, 1)
            k -= 1
            protected -= 1
        current = age
    propagate(current)
    return logscale+math.log(sum(x*(1-rho)**h for h,x in enumerate(v)))+lp(rho,sampled)+lp(1-rho,unsampled)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--driver', type=Path, help='Use an already compiled standalone driver')
    parser.add_argument('--output', type=Path, help='Optional JSON report path')
    args=parser.parse_args()
    reports=[]
    with tempfile.TemporaryDirectory(prefix='ddfbd-regression-') as temporary:
        driver=args.driver or Path(temporary)/'kernel_driver'
        if args.driver is None:
            subprocess.run([os.environ.get('CXX', 'clang++'), '-std=c++17', '-O2', '-Wall', '-Wextra',
                            '-I'+str(KERNEL), str(Path(__file__).with_name('kernel_driver.cpp')),
                            str(KERNEL/'DiversityDependentFbdLikelihood.cpp'), '-o', str(driver)], check=True)
        def call(mode, p, origin, cap=64, tol=1e-12, fossils=0, species=None, events=None, sampled=None, unsampled=None):
            argv=[str(driver), mode, *map(str,p), str(origin), str(cap), str(tol)]
            if mode == 'll':
                if species is not None:
                    events=schedule(species)
                    sampled, unsampled=counts(species)
                argv += list(map(str,[fossils,sampled,unsampled,len(events)]))
                for age, event in events:
                    argv += [str(age),event]
            return subprocess.run(argv, capture_output=True, text=True)
        def value(*a, **kw):
            result=call(*a, **kw)
            if result.returncode:
                raise AssertionError(result.stderr)
            return float(result.stdout)
        def agree(name, actual, expected, tolerance=2e-9):
            error=0 if actual == expected else abs(actual-expected)
            assert error < tolerance, (name,actual,expected,error)
            reports.append(dict(check=name, actual=actual, expected=expected, absolute_error=error))

        fixtures={
            'founder_protected': ([(3,2.4,0,True)],2),
            'fossil_extinct': ([(3,2.5,0.4,False)],3),
            'fossil_living_unsampled': ([(3,2.5,0,False)],3),
            'parent_protected_across_birth': ([(3,2.7,0.3,False),(2.1,1.6,0,True)],5),
            'unprotected_split_three_species': ([(3,1.2,0.5,False),(2.4,1.8,0,True),(1.9,1.1,0.2,False)],6),
            'two_buds_same_named_parent': ([(3,2.9,0,True),(2.4,1.8,0.5,False),(1.4,0.8,0,False)],7),
            'all_extinct': ([(3,2.7,0.8,False),(2.3,2.1,0.4,False)],4),
            'extant_only': ([(3,0,0,True),(1.8,0,0,True)],0),
        }
        p=(0.65,0.0,0.22,0.31,0.63)
        for name,(species,fossils) in fixtures.items():
            agree('constant_theorem2_'+name, value('ll',p,3,species=species,fossils=fossils),theorem2(species,fossils,p))
        simulator_path=ROOT/'validation/DDFBD/simulate_history.py'
        if simulator_path.exists():
            spec=importlib.util.spec_from_file_location('ddfbd_independent_simulator',simulator_path)
            simulator=importlib.util.module_from_spec(spec)
            sys.modules[spec.name]=simulator
            spec.loader.exec_module(simulator)
            for seed in range(8):
                hp=simulator.Parameters(lambda0=p[0],alpha=0,mu=p[2],psi=p[3],rho=p[4])
                history=simulator.simulate(hp,3,seed)
                pruned=simulator.prune_extended(history)
                if pruned['tree'] is None:
                    continue
                simulated_species=sorted([(row['structural_attachment_age'],row['oldest_observation_age'],
                    row['endpoint_age'],row['sampled_present']) for row in pruned['species']],reverse=True)
                fossils=pruned['fossilRecords']
                agree('simulated_extended_theorem2_seed'+str(seed),
                      value('ll',p,3,species=simulated_species,fossils=fossils),
                      theorem2(simulated_species,fossils,p))
        for mode in ['none','noextant']:
            reference_p=p if mode == 'none' else (*p[:3],0.0,p[4])
            agree('constant_'+mode, value(mode,p,3),constant_functions(3,reference_p)[0])

        for mode in ['logsampling','logextant']:
            reference_p=p if mode == 'logsampling' else (*p[:3],0.0,p[4])
            agree('direct_'+mode,value(mode,p,3,cap=256),math.log1p(-constant_functions(3,reference_p)[0]))
        agree('rare_fossil_sampling',value('logsampling',(0,0,0,1e-18,0),3),math.log(-math.expm1(-3e-18)))
        agree('rare_extant_sampling',value('logextant',(0,0,0,0,1e-18),3),math.log(1e-18))
        agree('no_birth_extant_log_underflow',value('logextant',(0,0,1,0,1),1000),-1000.)
        agree('no_birth_rare_extant_log_underflow',value('logextant',(0,0,1,0,1e-18),1000),-1000.+math.log(1e-18))
        agree('no_birth_tiny_fossil_exposure',value('logsampling',(0,0,0,1e-200,0),1e-200),2*math.log(1e-200))
        agree('no_birth_fossil_or_extant',value('logsampling',(0,0,.3,.4,.7),3),
              math.log(.4/.7*(-math.expm1(-.7*3))+.7*math.exp(-.7*3)))
        # To first order at 1e-18, P(fossil)=psi*integrated expected N, while
        # P(extant sample)=rho*expected N(T); second-order terms are negligible.
        agree('rare_fossils_with_birth_death',value('logsampling',(0.5,0,0.2,1e-18,0),3,cap=128),
              math.log(1e-18*math.expm1(0.3*3)/0.3),2e-8)
        agree('rare_extant_with_birth_death',value('logextant',(0.5,0,0.2,0,1e-18),3,cap=128),
              math.log(1e-18)+0.3*3,2e-8)
        for cap in [4,64]:
            agree('no_sampling_not_cap_leakage_'+str(cap),value('logsampling',(0.5,0,0.2,0,0),3,cap=cap),-math.inf)
        direct_caps=[value('logsampling',(0.8,0.03,0.3,0.2,0.4),3,cap=h) for h in [4,8,16,32,64]]
        assert all(b >= a-1e-10 for a,b in zip(direct_caps,direct_caps[1:]))
        assert abs(direct_caps[-1]-direct_caps[-2]) < 2e-7
        reports.append(dict(check='direct_conditioning_cap_convergence',caps=[4,8,16,32,64],
                            log_probabilities=direct_caps,final_difference=direct_caps[-1]-direct_caps[-2]))

        dd=(0.8,0.13,0.27,0.25,0.55)
        for name in ['parent_protected_across_birth','two_buds_same_named_parent','all_extinct']:
            species,fossils=fixtures[name]
            sampled,unsampled=counts(species)
            reference=independent_rk4(3,schedule(species),fossils,sampled,unsampled,dd)
            agree('DD_independent_RK4_'+name,value('ll',dd,3,species=species,fossils=fossils,cap=48),reference)
        species,fossils=fixtures['founder_protected']
        correct=value('ll',dd,3,species=species,fossils=fossils,cap=48)
        wrong=independent_rk4(3,schedule(species),fossils,1,0,dd,ignore_protection=True)
        assert wrong-correct > 0.1
        reports.append(dict(check='identity_constraint_changes_density', correct=correct, ignored_identity=wrong, log_difference=wrong-correct))

        pure=(0.8,0.35,0.0,0.0,1.0)
        species=[(3,0,0,True),(2.0,0,0,True),(0.7,0,0,True)]
        birth=lambda n:pure[0]*math.exp(-pure[1]*(n-1))
        exact=-birth(1)*1-2*birth(2)*1.3-3*birth(3)*0.7+math.log(birth(1))+math.log(birth(2))
        agree('DD_complete_pure_birth',value('ll',pure,3,species=species),exact)

        extreme=(0.5,1000.0,0.0,0.0,1.0)
        extreme_species=[(3,0,0,True),(2,0,0,True),(1,0,0,True)]
        extreme_exact=-0.5+2*math.log(0.5)-1000
        agree('DD_observed_birth_log_underflow',value('ll',extreme,3,species=extreme_species,cap=8),extreme_exact)
        agree('DD_extreme_complete_one_lineage',value('ll',(1,1000,0,0,1),1000,
              species=[(1000,0,0,True)],cap=8),-1000.0)
        species,fossils=fixtures['parent_protected_across_birth']
        near_constant=(p[0],1e-9,*p[2:])
        agree('DD_alpha_zero_limit',value('ll',near_constant,3,species=species,fossils=fossils),theorem2(species,fossils,p),2e-7)

        cap_p=(1.1,0.02,0.35,0.12,0.2)
        species,fossils=fixtures['two_buds_same_named_parent']
        caps=[4,8,16,32,64,128]
        capped=[value('ll',cap_p,3,species=species,fossils=fossils,cap=h) for h in caps]
        assert all(b >= a-1e-10 for a,b in zip(capped,capped[1:]))
        assert abs(capped[-1]-capped[-2]) < 2e-7
        reports.append(dict(check='DD_cap_convergence',caps=caps,log_likelihoods=capped,final_difference=capped[-1]-capped[-2]))
        loose=value('ll',dd,3,species=species,fossils=fossils,tol=1e-8)
        tight=value('ll',dd,3,species=species,fossils=fossils,tol=1e-13)
        agree('DD_tolerance_convergence',loose,tight,2e-7)

        boundary_cases=[
            ('zero_birth_fossil_extinct',(0,0,0.3,0.4,0.7),[(3,2,0.5,False)],2),
            ('zero_birth_present',(0,0,0.3,0.4,0.7),[(3,2,0,True)],2),
            ('zero_death_impossible',(0.4,0,0,0.3,0.7),[(3,2,0.5,False)],2),
            ('zero_fossils_possible',(0.4,0,0.3,0,0.7),[(3,0,0,True)],0),
            ('zero_psi_impossible',(0.4,0,0.3,0,0.7),[(3,2,0,True)],2),
            ('zero_rho_fossils',(0.4,0,0.3,0.2,0),[(3,2,0.5,False)],2),
            ('zero_rho_extant_impossible',(0.4,0,0.3,0.2,0),[(3,0,0,True)],0),
            ('complete_rho_unsampled_impossible',(0.4,0,0.3,0.2,1),[(3,2,0,False)],2),
            ('critical_birth_death',(0.4,0,0.4,0,0.7),[(3,0,0,True)],0),
            ('all_rates_zero',(0,0,0,0,1),[(3,0,0,True)],0),
        ]
        for name,pb,sp,f in boundary_cases:
            agree(name,value('ll',pb,3,species=sp,fossils=f),theorem2(sp,f,pb))
        agree('all_rates_zero_no_observation',value('none',(0,0,0,0,0.4),3),0.6)
        agree('no_sampling_no_observation',value('none',(0.4,0.2,0.3,0,0),3,cap=128),1.0)
        bad=call('ll',p,3,events=[(2,'P'),(1,'P')],fossils=2,sampled=1,unsampled=0)
        assert bad.returncode != 0 and 'protected' in bad.stderr
        reports.append(dict(check='invalid_protection_rejected',message=bad.stderr.strip()))

    # Strict JSON: infinite log densities are represented explicitly as strings.
    def safe(obj):
        if isinstance(obj,float) and not math.isfinite(obj):
            return '-Infinity' if obj < 0 else 'Infinity'
        if isinstance(obj,dict): return {k:safe(v) for k,v in obj.items()}
        if isinstance(obj,list): return [safe(v) for v in obj]
        return obj
    report=safe(dict(status='passed',scope='standalone extended-tree kernel only',checks=reports))
    text=json.dumps(report,indent=2,allow_nan=False)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text,end='')

if __name__ == '__main__':
    main()
