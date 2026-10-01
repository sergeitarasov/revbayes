"""Planning prototype, NOT a RevBayes distribution or production ODE solver.

Python standard library only. Fixed small RK4 steps are adequate for the supplied
fixtures, not arbitrary rate/age/cap combinations. Computes the raw oriented-tree
density with nonremoving sampling, before tree-label or conditioning corrections.
"""
import json
import math
from dataclasses import dataclass, replace


@dataclass(frozen=True)
class Parameters:
    lam: float = 0.5
    alpha: float = 0.2
    mu: float = 0.2
    psi: float = 0.3
    omega: float = 0.0
    rho: float = 0.7

    def birth(self, n):
        return self.lam * math.exp(-self.alpha * (n - 1)) if n else 0.0


def propagate(v, duration, k, p, step, transpose=False, protected=0):
    """exp(A_k duration) v, retaining original diagonal at the hidden cap."""
    hmax = len(v) - 1
    diag = [-(k+h)*(p.birth(k+h)+p.mu+p.psi+p.omega)
            for h in range(hmax+1)]
    up = [(h+2*k-protected)*p.birth(k+h) for h in range(hmax+1)]
    down = [h*p.mu for h in range(hmax+1)]

    def rhs(x):
        y = [diag[h]*x[h] for h in range(hmax+1)]
        for h in range(hmax):
            if transpose:
                y[h] += up[h]*x[h+1]
                y[h+1] += down[h+1]*x[h]
            else:
                y[h+1] += up[h]*x[h]
                y[h] += down[h+1]*x[h+1]
        return y

    nsteps = max(1, math.ceil(duration/step))
    dt = duration/nsteps
    for _ in range(nsteps):
        a = rhs(v)
        b = rhs([x+dt*y/2 for x,y in zip(v,a)])
        c = rhs([x+dt*y/2 for x,y in zip(v,b)])
        d = rhs([x+dt*y for x,y in zip(v,c)])
        v = [x+dt*(aa+2*bb+2*cc+dd)/6
             for x,aa,bb,cc,dd in zip(v,a,b,c,d)]
    if min(v) < -1e-14:
        raise ArithmeticError("RK4 lost positivity: reduce the step size")
    return v


def event(v, kind, k, p, transpose=False):
    if kind == "branch":
        return [x*p.birth(k+h) for h,x in enumerate(v)]
    if kind == "ancestor":
        return [x*p.psi for x in v]
    if kind == "fossil_tip":
        return ([p.psi*x for x in v[1:]] + [0.0] if transpose
                else [0.0] + [p.psi*x for x in v[:-1]])
    if kind == "occurrence":
        return [x*p.omega*(k+h) for h,x in enumerate(v)]
    raise ValueError(kind)


def normalize(v):
    scale = max(v)
    if scale <= 0:
        raise ArithmeticError("Zero density fixture")
    return [x/scale for x in v], math.log(scale)


def likelihood(origin, events, p, hmax=40, step=0.002, backward=False):
    """events are (age before present, kind), oldest first, excluding present."""
    k, age = 1, origin
    operations = []
    for next_age, kind in events:
        if not 0 <= next_age <= age or k < 0:
            raise ValueError("Invalid fixture event order")
        operations.append(("interval", age-next_age, k))
        operations.append((kind, 0, k))
        k += {"branch": 1, "fossil_tip": -1}.get(kind, 0)
        age = next_age
    operations.append(("interval", age, k))
    terminal = [p.rho**k * (1-p.rho)**h for h in range(hmax+1)]
    v = terminal if backward else [1.0] + [0.0]*hmax
    logscale = 0.0
    for kind, duration, active in (reversed(operations) if backward else operations):
        if kind == "interval":
            v = propagate(v, duration, active, p, step, backward)
        else:
            v = event(v, kind, active, p, backward)
        v, shift = normalize(v)
        logscale += shift
    result = v[0] if backward else sum(x*y for x,y in zip(v,terminal))
    return math.log(result) + logscale


def homogeneous_reference(origin, events, p):
    """Independent scalar constant-rate FBD solution (omega=alpha=0)."""
    assert p.alpha == p.omega == 0
    gamma = p.lam+p.mu+p.psi
    delta = math.sqrt(gamma*gamma-4*p.lam*p.mu)
    low, high = (gamma-delta)/(2*p.lam), (gamma+delta)/(2*p.lam)
    ratio0 = (1-p.rho-low)/(1-p.rho-high)

    def e(age):
        ratio = ratio0*math.exp(-delta*age)
        return (low-high*ratio)/(1-ratio)

    def log_d(age):
        return (-delta*age+2*math.log(abs(1-ratio0))
                -2*math.log(abs(1-ratio0*math.exp(-delta*age))))

    k, age, result = 1, origin, 0.0
    for next_age, kind in events:
        result += k*(log_d(age)-log_d(next_age))
        if kind == "branch":
            result += math.log(p.lam)
            k += 1
        elif kind == "ancestor":
            result += math.log(p.psi)
        elif kind == "fossil_tip":
            result += math.log(p.psi*e(next_age))
            k -= 1
        else:
            raise ValueError(kind)
        age = next_age
    return result + k*log_d(age) + k*math.log(p.rho)


def run_checks():
    reports = []
    fixtures = {
        "one_extant": [],
        "two_extant": [(1.4,"branch")],
        "sampled_ancestor": [(1.4,"branch"),(0.8,"ancestor")],
        "terminal_fossil": [(1.4,"branch"),(0.8,"fossil_tip")],
        "fossil_only": [(0.8,"fossil_tip")],
        "mixed": [(1.7,"branch"),(1.2,"branch"),
                  (1.0,"ancestor"),(0.7,"fossil_tip")],
    }
    p = Parameters(alpha=0)
    for name, events in fixtures.items():
        value = likelihood(2.0, events, p)
        reference = homogeneous_reference(2.0, events, p)
        error = abs(value-reference)
        assert error < 1e-8, (name,error)
        reports.append(dict(check="constant_FBD_"+name, log_likelihood=value,
                            reference=reference, absolute_error=error))

    p = Parameters(omega=0.15)
    events = fixtures["mixed"] + [(0.3,"occurrence")]
    forward = likelihood(2.0, events, p)
    backward = likelihood(2.0, events, p, backward=True)
    finer = likelihood(2.0, events, p, step=0.001)
    assert abs(forward-backward) < 1e-9
    assert abs(forward-finer) < 1e-8
    reports.append(dict(check="DD_forward_adjoint_and_step_halving",
                        forward=forward, backward=backward,
                        adjoint_error=abs(forward-backward),
                        step_halving_error=abs(forward-finer)))

    p = Parameters(alpha=0.4, mu=0, psi=0, rho=1)
    events = [(1.5,"branch"),(0.6,"branch")]
    exact = (-0.5*p.birth(1)-0.9*2*p.birth(2)-0.6*3*p.birth(3)
             +math.log(p.birth(1))+math.log(p.birth(2)))
    value = likelihood(2, events, p)
    assert abs(value-exact) < 1e-8
    reports.append(dict(check="DD_complete_pure_birth", log_likelihood=value,
                        reference=exact, absolute_error=abs(value-exact)))

    p = Parameters(lam=0.9, alpha=0.025, mu=0.35, psi=0.12, rho=0.25)
    values = [likelihood(3,fixtures["mixed"],p,hmax=h) for h in [4,8,16,32,64]]
    assert all(b >= a-1e-10 for a,b in zip(values,values[1:]))
    assert abs(values[-1]-values[-2]) < 1e-7
    reports.append(dict(check="DD_hidden_cap_convergence", caps=[4,8,16,32,64],
                        log_likelihoods=values, final_difference=values[-1]-values[-2]))

    p = Parameters(alpha=0)
    limiting = likelihood(2,fixtures["mixed"],replace(p,alpha=1e-8))
    baseline = homogeneous_reference(2,fixtures["mixed"],p)
    assert abs(limiting-baseline) < 1e-6
    reports.append(dict(check="DD_alpha_zero_limit", absolute_error=abs(limiting-baseline)))

    # A restricted named-species benchmark, NOT the general range likelihood:
    # the founder retains its identity for the entire interval, has two fossil
    # records, and is sampled alive at present. Buds can be hidden, but the
    # observed path cannot switch to a bud. Thus its coefficient is h+1.
    p = Parameters(alpha=0)
    v = propagate([1.0]+[0.0]*40, 2, 1, p, 0.002, protected=1)
    locked = math.log(sum(x*p.rho*(1-p.rho)**h for h,x in enumerate(v)))
    locked += 2*math.log(p.psi)
    log_d = homogeneous_reference(2, [], p)-math.log(p.rho)
    exact = math.log(p.rho)+2*math.log(p.psi)+0.5*log_d-(p.lam+p.mu+p.psi)
    assert abs(locked-exact) < 1e-8
    reports.append(dict(check="named_species_locked_founder_restricted_case",
                        log_likelihood=locked, reference=exact,
                        absolute_error=abs(locked-exact)))
    return dict(status="passed", scope="planning prototype only; raw oriented density",
                checks=reports)


if __name__ == "__main__":
    print(json.dumps(run_checks(), indent=2))
