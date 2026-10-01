# Fossil MCMC integration and numerical validation

Run from the repository root, using the locally compiled binary:

```sh
mkdir -p tests/test_DDFBD/results
python3 tests/test_DDFBD/numerical_reference.py --output tests/test_DDFBD/results/numerical_reference.json
python3 tests/test_DDFBD/run_integration.py .local-build/revbayes-build/rb --generations 2000 --output tests/test_DDFBD/results/integration.json
python3 tests/test_DDFBD/audit_samples.py tests/test_DDFBD/results/integration --output tests/test_DDFBD/results/sample_audit.json
```

`run_integration.py` copies the example into a temporary directory. It checks the absolute labelled tree density against the independent constant-rate formula, then runs morphology, DNA and combined-data MCMC. It enables RevBayes's `debugMCMC=1` cache checks, exercises both sampling-conditioning options, and writes/restores/resumes a checkpoint. Successful output files and console logs are copied beside the requested JSON report. The current script includes the conditioning tests; the earlier saved 2,000-generation run predates their addition, so those results are separately recorded in `conditioned_integration.json`.

`audit_samples.py` reads all saved tree/parameter pairs, preserving child order. It checks extant versus extinct endpoints, fossil placement, synchronized fossil ages in character observations, movement of each parameter/clock/age and ordered topology. It recomputes raw fossil-tree densities at hidden caps 64, 128 and 256 for about 20 draws per mode. Density differences compare the same rounded saved state at each cutoff; these are not comparisons against full-precision in-memory posterior values.

## Verified results, 1 October 2026

- 52 independent numerical checks passed, including the analytical constant-rate limit, independent RK4 diversity-dependent calculations, simulated extended histories, parameter boundaries, cutoff/tolerance convergence and rare-observation log probabilities.
- Three 2,000-generation chains passed (201 saved states per mode).
- A further 200-generation chain in each mode passed with the final numerical-boundary build.
- Two additional 100-generation combined chains passed: any-sampling and sampled-extant conditioning.
- Checkpoint write, restore, and resumed MCMC passed.
- All 603 saved states passed the support/shared-age audit. The combined, morphological, and DNA chains visited 12, 12, and 11 distinct ordered topologies, respectively.
- Maximum raw log-density difference between caps 64 and 256 across 63 selected states: 2.92e-13.
- The existing DDD benchmark passed again after the shared kernel changes.

Results and traces are in `results/`. This verifies executable joint inference on the supplied synthetic example; these short runs do not establish mixing, statistical calibration, or parameter identifiability. Living/extinct species status remains fixed from the initial tree. All current examples use exponential decline in total living diversity with constant extinction and fossil-sampling rates.

## Bugs corrected during integration

1. Unnamed tip-move bounds prevented MCMC from reconnecting the move to its cloned DAG. Named bounds are included in `model(...)`.
2. The initial example's fixed extinction bound of 0.8 excluded valid endpoints up to the youngest fossil's actual age. The bound is now linked to that age.
3. Conditioning could report negative infinity for positive but underflowed observation probability. The no-birth boundary now uses an exact log-space expression, including log probabilities below -1000; other positive-mass underflow raises a numerical error.
