# 100-species DD-FBD recovery pilot

Prepared 1 October 2026. Dataset: 100 observed named species, 105 recovered occurrences, 64 sampled extant species, 28 observed extinct species, and 8 fossil-bearing species alive but unsampled at present. There are 169 character observations, each with 100 binary sites. The full unobserved species history is preserved in `data/history.json`.

Generating parameters: lambda0=0.8, alpha=0.015, mu=0.15, psi=0.3, rho=0.8, origin=12. Birth rate per species is lambda0 exp[-alpha(N-1)], where N includes all living species. Characters were simulated independently on the complete history with a symmetric two-state CTMC, rate 0.2, including constant sites. Inference uses `coding="all"`.

The simulation used seed 1069, the first seed from 1001 onward giving exactly 100 observed species. All attempted sizes are recorded. This deliberately size-selected single-dataset experiment is a recovery pilot, not a simulation-based calibration study. Inference uses the ordinary unconditioned likelihood; no correction for the dataset-selection rule is claimed.

## Jobs

- `fixed_tree`: 2,000 tuning/burn-in generations + 20,000 sampling generations. Known extended species tree; estimates lambda0, alpha, mu and psi.
- `joint_1` and `joint_2`: each 5,000 tuning/burn-in + 30,000 sampling generations. Infer the ordered extended tree, node times, extinct endpoints and the four rates. Scalar starting values differ between chains; both initialize at the simulated tree.
- Origin, rho, character clock, occurrence ages and living/extinct status are known/fixed. This first experiment isolates recovery of the diversification/fossil-sampling parameters; it does not test recovery of every possible unknown.
- Priors: lambda0 ~ Exp(1), alpha ~ Exp(50), mu ~ Exp(3), psi ~ Exp(2), using rate parameterization. Priors do not fix generating parameter values.
- Parameters are saved every 10 sampling generations, species trees every 100, checkpoints every 1,000. Hidden cap 128, propagation tolerance 1e-11. At the true values and both scalar starts, cap 128 versus 256 agrees within 1e-7 (`initial_cutoff_check.json`). Posterior-draw cutoff sensitivity should also be checked before scientific interpretation.

## Background execution and checking results

`background_job.py` supervises three concurrent RevBayes processes. `status.json` records supervisor/child PIDs, input hashes, last saved generations, and failures/completion. `launch.json` records how the detached job was started. Each chain has its own `.stdout`; posterior logs and checkpoints live in `output/`.

`recovery_summary.json` updates while chains sample and when they finish. It compares posterior mean, median and 95% equal-tail interval to the generating values, reports approximate batch-means ESS and classical split Rhat for the two joint chains, and discards the first 20% of saved samples in addition to explicit burn-in. Low ESS or poor Rhat means the recovery comparison is not yet trustworthy; truth lying inside an interval on one selected dataset is not a calibration result.

From this directory:

```sh
cat status.json
python3 summarize.py
# Inspect one chain's console log:
tail -n 30 joint_1.stdout
```

To stop this specific job, send SIGTERM to `supervisor_pid` in `status.json`; the supervisor terminates its own three children. Do not rerun preparation or launch duplicate jobs into this directory. Checkpoints support manual continuation via RevBayes `initializeFromCheckpoint` using the identical model.

A ten-generation pilot passed before the long jobs were launched. The fixed-tree control and two joint chains run concurrently; their results should be compared separately. macOS sleep inhibition, if available, is scoped to the lifetime of the supervisor process, not a persistent system setting.
