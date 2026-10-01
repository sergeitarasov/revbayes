# Extant-tree benchmark against DDD

This side experiment implements a complete-extant-tree likelihood in RevBayes using the **same interval propagator** as the new fossil model. It integrates over unobserved extinct lineages: rates depend on total living diversity, not just the reconstructed tree. It is a callable likelihood function, not yet a new stochastic tree distribution.

From the repository root, after rebuilding RevBayes:

```sh
sh projects/meson/generate_sources.sh core
sh projects/meson/generate_sources.sh revlanguage
.local-build/run-build.sh
python3 tests/test_DDD_benchmark/run.py
```

The driver requires a C++ compiler, Python 3, R, and installed R packages DDD and ape. It does not install or modify DDD. The build helper and binary are local ignored artifacts from this development workspace; another checkout needs its normal RevBayes build, passed as `--rb /absolute/path/to/rb`.

## Rev usage

```rev
# Run from the repository root.
tree = readTrees("tests/test_DDD_benchmark/results/simulated.tre")[1]
lnL := fnDiversityDependentLogLikelihood(tree, lambda0=0.8, mu=0.2,
           K=8, rateModel="DDDpower", start="crown", condition="survival")
print(lnL)
```

For a stem analysis set `start="stem", originAge=4.0` (at least the root age). For the originally proposed decay use `rateModel="exponential", alpha=0.2`. No fossils or morphological/molecular data enter this function. It returns a number directly, so a `.lnProbability()` call is unnecessary. Assigning it with `:=` keeps it linked to its parameters; simply adding this deterministic number to a model does **not** add a likelihood term to MCMC.

## Matching conventions

| Setting | RevBayes | DDD |
|---|---|---|
| Linear birth | `DDDlinear` | `ddmodel=1` |
| Power birth | `DDDpower` | `ddmodel=2` |
| Age only | `condition="none"` | `cond=0` |
| Survival | `condition="survival"` | `cond=1` |
| Extant count | `condition="nTaxa"` | `cond=2` |
| Crown / stem | `start="crown"` / `"stem"` | `soc=2` / `1` |
| Phylogeny convention | `density="DDDphylogeny"` | `btorph=1` |
| Branching-time convention | `density="branchingTimes"` | `btorph=0` |
| Complete extant sampling | required | `missnumspec=0` |
| Hidden states 0,...,H | `maxHiddenLineages=H` | `res=H+1` before DDD's internal linear-model adjustment |

DDD calls model 2 “exponential”, but its implementation uses **a power law**:

- Our original exponential: lambda(N) = lambda0 exp[-alpha(N-1)].
- DDD model 1: lambda(N) = max(0, lambda0 - (lambda0-mu) N/K).
- DDD model 2: lambda(N) = lambda0 (N+1)^[-log(lambda0/mu)/log(K+1)].

The two nonconstant DDD forms are compared without modifying the package. Our exponential form has independent pure-birth and fossil-kernel reduction tests; its alpha=0 limit is also compared to DDD's constant-rate model. Nonzero exponential alpha is **not** falsely matched to DDD model 2.

## Shared kernel and normalization

During each inter-node interval, k retained lineages and h hidden lineages give N=k+h. There are no observations until the present, so the fossil kernel's protected count is zero. Its forward weighted embedding operator is

- A[h,h] = -(k+h)(lambda(k+h)+mu),
- A[h+1,h] = (h+2k)lambda(k+h),
- A[h-1,h] = h mu.

An internal node multiplies component h by lambda(k+h), then increments k. The endpoint selects h=0 because all living species are sampled. A crown start uses k=2 and omits the root birth factor; a stem start uses k=1 and includes all nodes. Branching-time density adds log Gamma(n); the default matches DDD's phylogeny measure, with no further FBD label factorial.

For DDD's conditioning denominator, propagate the same operator with k fixed at its initial value for the full start age. Divide component h by h+1 for a stem start, or (h+2)(h+3)/6 for a crown start. Sum these weights for survival, or select h=n-k for conditioning on extant count. This accounts for both crown sides surviving under the shared diversity process; it is not the product of two independent single-lineage survival probabilities.

Our upper boundary kills overflow. DDD 5.2.5's matrix implementation instead removes birth loss from its top diagonal. Low-cutoff values therefore need not match. `cutoff_convergence.tsv` deliberately exposes this difference and verifies convergence at larger caps. Propagation tolerance and cutoff error are separate controls.

## Checks and saved evidence

- 144 likelihood comparisons: two trees, two nonconstant rate laws, three parameter sets, stem/crown, three conditioning conventions, and two cutoffs.
- Six additional constant-rate comparisons against DDD.
- A second DDD backend using adaptive ODE integration.
- A difficult cutoff sweep through 256 hidden lineages.
- Nineteen direct reduction, analytic pure-birth, normalization, and invalid-age checks.
- Actual compiled Rev evaluation of all 144 cases, the fitted point, branching-time normalization, and a live DAG parameter update.
- Existing fossil-kernel numerical regressions.

`benchmark.R` simulates a tree with `DDD::dd_sim`, seed 1439, parameters lambda0=0.8, mu=0.2, K=8, model 2, crown age 3. DDD's simulation conditions on both crown sides surviving. It then fits **lambda0 only**, keeping mu=0.2 and K=8 fixed, through `DDD::dd_ML` and independently through the C++ likelihood. The likelihood surface, optimizer settings, and result are saved. This is an implementation benchmark, not a parameter-recovery or identifiability study.

`results/` contains numerical tables, the two Newick trees, complete simulated history (RDS), generated executable Rev comparisons, fit estimates, source-function snapshots from the installed reference package, session versions, and the final run report. `run.json` records the exact tested kernel and Rev binary hashes.

Reference: [DDD package and source](https://rsetienne.github.io/DDD/), [dd_loglik documentation](https://search.r-project.org/CRAN/refmans/DDD/html/dd_loglik.html), and Etienne et al. (2012), [doi:10.1098/rspb.2011.1439](https://doi.org/10.1098/rspb.2011.1439).
