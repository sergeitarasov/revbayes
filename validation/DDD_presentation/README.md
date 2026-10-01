# Presentation report: diversity-dependent tree likelihoods

The finished report is `output/pdf/ddd_simulation_report.pdf` in the repository root. It is a landscape PDF suitable for presenting directly. Six independent figures are supplied as editable vector SVG and high-resolution PNG in `figures/`, and collected in `ddd_presentation_figures.zip`. The teaching scripts are in `examples/DDD/`.

## Experiment

We simulate **100 trees** with the unmodified DDD 5.2.5 package:

| Setting | Value |
|---|---|
| Rate models | DDD model 1 (linear), model 2 (power law) |
| Generating parameters | lambda0=0.8, mu=0.2, K=15 |
| Crown ages | 3 and 8 time units |
| Replicates | 25 per model-age combination |
| Sampling | All living species sampled |
| Ascertainment | Both crown sides survive, as in DDD::dd_sim |
| Additional selection | None; all accepted trees retained |
| Fitting | Only lambda0; mu and K fixed at generating values |
| Scalar search range | [0.20001,8], including endpoints |
| Hidden cutoff | H=128; fitted likelihood rechecked at H=256 |
| Propagation tolerance | 1e-13 |

Each tree is evaluated at three parameter points. These are lambda0=0.4,0.6,0.8 for the linear model and 0.6,0.8,1.2 for the power model. Using lower linear-model rates avoids parameter combinations that assign exactly zero likelihood to trees exceeding the model's attainable diversity. Fits still explore the stated full range and retain valid zero-likelihood regions in the objective.

The same bounded scalar search is applied to the two independent likelihood implementations: a 13-point logarithmic scan, local R `optimize` search around the best grid point, and explicit endpoint comparisons. This is not a certificate of global optimization under arbitrary likelihood shapes. Every simulation and bound hit is preserved.

The actual Rev executable checks 300 parameter points, 100 fitted points and an 81-point profile: **481 evaluations**. Scalar fitting uses the standalone C++ driver linked to the same production kernel; the final parameter points are also checked through the Rev wrapper.

## Results

- Maximum discrepancy in the 300 likelihood comparisons: approximately 2.21e-12.
- Maximum discrepancy in the 481 actual Rev checks: approximately 5.17e-12.
- Maximum fitted log-likelihood difference: approximately 7.93e-12.
- Maximum difference between the two lambda0 estimates: approximately 2.12e-6.
- Largest H=128 to H=256 fitted-likelihood change: approximately 2.06e-13.
- Bound hits: 11 lower, 3 upper, out of 100 fits.

The exact values, source and binary hashes are in `results/manifest.json`. `group_summary.tsv` reports the distribution of estimates. The numerical agreement establishes matched implementation behavior for the tested cases. It does not establish accurate joint recovery of lambda0, mu and K, calibrated confidence intervals, or fossil-model MCMC validity. The recovery experiment is a pilot with 25 trees per scenario and known mu/K. Bound hits are restricted-search results, not claims of unconstrained optima.

The teaching example is a separate six-tip tree with **K=8**, seed 1439. It is retained for continuity with the original benchmark; the presentation simulation study uses **K=15**.

## Reproduce

From `rb-fbdp/revbayes`, first build the local Rev executable and the standalone driver if needed:

```sh
.local-build/run-build.sh
c++ -std=c++17 -O2 -I src/core/functions/phylogenetics \
  tests/test_DDD_benchmark/driver.cpp \
  src/core/functions/phylogenetics/DiversityDependentFbdLikelihood.cpp \
  -o .local-build/ddd-kernel
```

Then:

```sh
Rscript validation/DDD_presentation/simulate_and_compare.R
.local-build/revbayes-build/rb validation/DDD_presentation/results/verify_in_Rev.Rev \
  > validation/DDD_presentation/results/Rev_verification.txt
.local-build/ddd-kernel DDDpower 1.2 0.8 0 30 1 1 512 10 10 6 3 \
  > validation/DDD_presentation/results/stress_H512.txt
Rscript validation/DDD_presentation/make_figures.R
python3 validation/DDD_presentation/build_report.py
```

R dependencies: `DDD`, `ape`, `ggplot2`, `cowplot`, `svglite`. PDF dependencies: Python `reportlab`, `Pillow`; DejaVu fonts. Set `DDD_FONT_DIR` to your DejaVu font directory if not using the bundled local runtime. PDF validation also uses `pypdf` and `pdftoppm`. On this workspace, the report was built with:

```
/Users/taravser/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3
```

`simulate_and_compare.R` accepts an optional replicate count for experiments, but the saved figure/report builders intentionally assert the 100-tree design. Update those checks and report text for a different study. Simulated histories and completed fits are reused on rerun; for a changed fitting protocol, use a fresh results directory or move the previous one aside first. Do not silently reuse estimates under a changed model or search range.

The cutoff figure uses the difficult case saved by `tests/test_DDD_benchmark/benchmark.R`, rather than pooling it with the 100 new simulations. Rev uses a killed upper boundary, while DDD's matrix backend modifies its top diagonal. Low-cutoff values differ; both require convergence. In the stress case the Rev H=256 versus H=512 change is about 4.4e-14.

## Evidence and figure guide

- `results/trees/*.rds`: full simulated histories, parameters, age and seed; one file per replicate.
- `results/trees/*.tre`: reconstructed extant trees actually passed to Rev.
- `simulations.tsv`: species counts and metadata, including extinct species.
- `likelihood_comparisons.tsv`, `fits.tsv`: unrounded numerical results.
- `verify_in_Rev.Rev`, `Rev_verification.txt`: executable checks and actual run output.
- `likelihood_profile.tsv`: 81 matched points for one illustrative tree.
- `illustrative_metadata.tsv`: selected power-law age-8 tree closest to its group's median tip count; no simulations excluded.
- `illustrative_diversity.tsv`: true diversity versus reconstructed lineages. A 1e-10 event-time tolerance reconciles roundoff between the full history and reconstructed tree; it is only for plotting event coincidences.
- `session.txt`: R and package versions; `manifest.json`: test metrics and hashes.
- `figures/01`: compare the rate laws; `02`: observed tree and hidden history; `03`: numerical agreement; `04`: likelihood profile; `05`: estimation variation; `06`: truncation convergence.

Source and theory references: [DDD package](https://rsetienne.github.io/DDD/), [dd_sim documentation](https://search.r-project.org/CRAN/refmans/DDD/html/dd_sim.html), Etienne et al. (2012), [doi:10.1098/rspb.2011.1439](https://doi.org/10.1098/rspb.2011.1439). Equations and implementation conventions are also recorded in `tests/test_DDD_benchmark/README.md`.
