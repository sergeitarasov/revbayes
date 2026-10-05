# Specimen diversity-dependent FBD prior

`dnDDFBDP` (full name `dnDiversityDependentSpecimenFBD`) is a new sampled-tree
prior, separate from the earlier named-species `dnDDFBD` model.

The per-lineage birth rate is

    lambda(N) = lambda0 * exp(-alpha * (N-1)), alpha >= 0.

N is the total number of living lineages, including unobserved lineages.
Extinction mu and preservation psi are constant. Sampling does not remove a
lineage. Each data label identifies one specimen or living observation. A fossil
can be a terminal observation or a sampled ancestor. Only biological splits are
speciation events; the artificial binary node representing a sampled ancestor
is not a split. Hidden splits and post-fossil continuations are marginalized.

```rev
tree ~ dnDDFBDP(originAge=origin, lambda0=lambda0, alpha=alpha,
                mu=mu, psi=psi, rho=1, taxa=taxa, condition="sampling")
morphology ~ dnPhyloCTMC(tree=tree, Q=fnJC(2), branchRates=clock,
                        type="Standard", coding="variable")
morphology.clamp(morpho)
moves.append(mvCollapseExpandFossilBranch(tree, origin, weight=6))
```

## Run the morphology example

See the [build and model guide](../../DDFBDP.md) to build this branch.
From the repository root, using a binary built with the new source files:

```sh
build/rb examples/DDFBDP/morphology.Rev
```

The script runs 1,000 generations and writes `output/DDFBDP/morphology.log` and
`morphology.trees`. This is an integration smoke test, not a converged empirical
analysis. It estimates lambda0, alpha, mu, psi, origin, morphology clock and gamma
shape, jointly with topology, node times, fossil ages and ancestor status. It uses
22 labels (8 living, 14 fossil) and 62 binary morphological characters; four labels
lack morphology and are added as missing data. There is no molecular partition.

The morphology and fossil age bounds come from the RevBayes
[specimen tutorial](https://revbayes.github.io/tutorials/fbd/fbd_specimen), by
Tracy Heath, Walker Pett and April Wright. The morphology matrix is copied
unchanged. The eight living taxa have age bounds set to exactly zero; the 14
fossil bounds are unchanged and represent uncertainty in individual occurrence
ages. Frozen original files and download hashes are in
`doc/diversity-dependent-fbd/specimen-tutorial-reference/`.

## Tests

```sh
python3 tests/test_DDFBDP/numerical_reference.py
python3 tests/test_DDFBDP/run_integration.py build/rb \
    --output tests/test_DDFBDP/results/integration.json
```

The first command compares the standalone kernel to analytic constant-rate FBD
and an independent RK4 implementation of the DD count equations. The second
uses a 200-generation diagnostic chain by default and compares actual Rev tree densities to `dnBDSTP`, runs morphology-only MCMC with
cache checking enabled, checks saved states at larger hidden-count cutoffs, and
exercises constrained topology and checkpoint restoration. Results are stored
in `tests/test_DDFBDP/results/`.

## Scope and numerical controls

- `condition="time"`: unconditioned density from a single lineage at originAge.
- `condition="sampling"`: condition on at least one fossil or living observation.
- `condition="survival"`: condition on at least one sampled living lineage.
- `maxHiddenLineages=128`: killed numerical boundary, not a biological carrying
  capacity. Increase it to assess density and posterior sensitivity.
- `numericalTolerance=1e-12`: propagation tolerance, not a bound on cutoff error.
- `initialTree`: optional ordinary sampled-ancestor TimeTree matching taxa and ages.

There are no skylines, species ranges, lineage-specific speciation modes, removal,
root-age starts or count conditioning. At least two observation labels are required.
General prior simulation is not yet implemented; MCMC initialization uses a valid
starting-tree generator and is explicitly not a prior draw. The derivation and
remaining simulation-recovery work are in
`doc/diversity-dependent-fbd/SPECIMEN_PLAN.md`. Validation of the earlier
named-species model does not establish parameter recovery for this new model.

Native compatibility detail: at alpha=0, `time` and `sampling` match `dnBDSTP`
directly. For `survival`, this prior divides by the actual probability of at
least one sampled extant descendant. Current native BDSTP divides that
probability by rho, so its log density is `ln(rho)` lower when rho<1. The
regression checks this difference explicitly. Both use BDSTP's legacy
log-factorial approximation for the fixed label-normalization constant.

## Validation completed on 4 October 2026

- Full RevBayes build succeeded.
- 26 standalone numerical checks passed; maximum independent comparison error
  was 2.4e-12 in log density.
- 46 native BDSTP comparisons passed: 39 direct, 7 with the documented survival
  rho correction. Maximum comparison error was 2.1e-12.
- The previous named-species kernel's 52 regression checks still passed.
- A 200-generation morphology-only diagnostic chain completed with debug cache
  checks enabled: 21 saved states, ancestor counts from 0 to 10, 212 accepted
  collapse/expand proposals out of 1,159, and changes in all estimated parameters.
- Five saved states were checked with hidden caps 64, 128 and 256. The largest
  absolute log-density difference was 5.4e-10 for 64 versus 128 and 1.3e-10 for
  128 versus 256. These are checks on these states, not a global error bound.
- A separate constrained-topology chain successfully saved a checkpoint,
  restored it and continued sampling (60 generations across the three segments).

Reports and the 200-generation diagnostic trace are in
`tests/test_DDFBDP/results/integration.json` and the adjacent `integration/`
directory. These tests establish numerical agreement and working MCMC integration;
parameter recovery and converged empirical inference remain untested.
