# Implementation and resumption checkpoint

Updated 1 October 2026. Main fossil-model work resumed after the successful extant-tree/DDD benchmark. End-to-end synthetic fossil MCMC integration now passes; statistical calibration remains future work. The model and supporting materials are on branch `ddfbd-model`, based on upstream `development` commit `5bd797b1c0d3af4c903343af856d7931b289024f`.

For an accessible explanation of the model, data mapping and code, see the [biological model guide](MODEL_GUIDE.md).

## Completed side experiment: extant DD versus DDD

See [benchmark instructions and evidence](../../tests/test_DDD_benchmark/README.md) and its `results/run.json`. `fnDiversityDependentLogLikelihood` evaluates a complete extant time tree with rates depending on all living lineages. It shares the fossil kernel's propagator and matches DDD's linear and power rate laws, stem/crown convention, and conditions 0/1/2. The original exponential rate law is also available. DDD model 2 is a power law and must not be equated to exponential decay in N.

The benchmark includes actual Rev calls, an independent unmodified DDD installation, simulated and fixed trees, a matched one-parameter maximum-likelihood fit, cutoff convergence, and reduction to the fossil kernel. It is not an MCMC distribution or a statistical recovery study. Full results and exact build hashes are saved by `python3 tests/test_DDD_benchmark/run.py`.

## Presentation report and simulations (1 October 2026)

The user requested a runnable example and presentation material before resuming the main fossil task. `examples/DDD/likelihood.Rev` is now self-contained, with an included six-tip tree, cutoff check, three rate laws and a likelihood profile; `compare_in_R.R` is the paired DDD reference. Both were run successfully.

`validation/DDD_presentation/README.md` documents 100 new DDD simulations (two rate models, two crown ages, 25 replicates each), 300 paired likelihood points, 100 one-parameter fits, and 481 checks in the actual Rev executable. The maximum Rev-versus-DDD discrepancy is about 5.17e-12; 11 lower-bound and 3 upper-bound fits are preserved. This is a known-mu/K pilot, not joint parameter recovery. The landscape report is `output/pdf/ddd_simulation_report.pdf`; six SVG/PNG figures are also packaged in `validation/DDD_presentation/ddd_presentation_figures.zip`. Full histories, seeds, tables, reference session versions and hashes are saved. The fossil MCMC integration status is summarized separately below.

## Main fossil implementation already present

- `DiversityDependentFbdLikelihood.{h,cpp}`: shared positive-uniformization hidden-count kernel, killed upper boundary, log scaling, direct observation-probability absorber.
- `DiversityDependentFossilizedBirthDeathProcess.{h,cpp}` and `Dist_diversityDependentFBD.{h,cpp}`: `dnDiversityDependentFBD` / `dnDDFBD` DAG distribution; named repeated occurrences, sampled extants, initial species tree, origin and rates.
- `SpeciesObservationTreeFunction` / `fnSpeciesObservationTree`: projection onto fossil/extant character observations at their actual ages; copies preserve listeners and notify CTMC caches.
- `OrderedTreeSwapProposal` / `mvOrderedTreeSwap`: reversible swaps preserve child slots and support orientation changes.
- Examples in `examples/DDFBD`: morphology, DNA, combined data; five uncertain fossil occurrence ages assigned to three named species.
- `validation/DDFBD/simulate_history.py`: independent full-history simulation and pruning, with saved seed-7 history. Not SBC.
- Existing `FossilizedBirthDeathSpeciationProcess` const virtual override bug was corrected and demonstrated by an original-versus-fixed dispatch regression.

The extended budding-tree representation has one endpoint per observed named species. Child 0 continues the ancestor, child 1 is the daughter. Endpoints are true extinction ages or the present. Species identity is protected from its oldest occurrence to its endpoint. This makes the protection count deterministic and allows total diversity N=k+h without additional identity flags. See [the derivation](EXTENDED_TREE_DERIVATION.md); the early planning/review files preserve preliminary concerns superseded by that derivation.

## Main fossil validation and unresolved work

**Joint fossil MCMC now runs end to end.** Saved evidence lives in [tests/test_DDFBD/results](../../tests/test_DDFBD/results), with commands and scope in [the test README](../../tests/test_DDFBD/README.md).

- Morphology-only, DNA-only, and combined chains each ran 2,000 generations with `debugMCMC=1` (cache recomputation checks), producing 201 saved samples each.
- Additional 100-generation combined-data chains passed with `condition="sampling"` and `condition="sampledExtant"`.
- Checkpoint write/read/resume passed. A further 200-generation run of each data mode passed using the final numerical-boundary build.
- All inferred scalar parameters, character clocks and five fossil ages moved; B's extinction endpoint moved; the three chains visited 12, 12 and 11 distinct ordered topologies.
- All 603 saved draws passed age/support and shared-fossil-age checks. For 63 selected draws, cutoff 64 versus 256 changed raw log density by at most 2.92e-13.
- The independent numerical suite now passes 52 checks, including exact no-birth sampling probabilities well below floating-point probability range.
- The DDD benchmark was rerun successfully after these changes.

The initialization error was caused by unnamed constant min/max nodes used by the tip-time proposal. `model.Rev` now supplies named bounds and explicitly includes them in `model(...)`. The upper bound is linked to B's youngest fossil age, rather than fixed at 0.8: the fixed bound excluded otherwise valid extinction endpoints. Existing core projection, ordered-swap and FBDSP dispatch regressions had already passed; focused upstream FBD/projection results remain in `.local-build/focused-regressions/summary.json`.

`logObservationProbability` now computes the no-birth case in log space and reports numerical underflow for positive observation mass instead of silently returning a biologically impossible zero.

These are implementation/integration checks on synthetic data, not convergence evidence or SBC.

Other pending limitations:

- Alive/extinct status is fixed from the initial extended tree; it is currently a known-status/conditional analysis, not status inference.
- Prior redraw with fixed named occurrence data is deliberately unsupported. MCMC initialization restores the supplied compatible initial tree.
- Singleton projected character trees are incompatible with the existing CTMC root handling.
- Hidden-count convergence is separate from propagation tolerance, including the conditioning denominator.
- Recovery/SBC, replicated-chain convergence, larger-data performance and identification of parameters remain unfinished.

## Build and next commands

Local binary: `.local-build/revbayes-build/rb`. Dependencies are contained in ignored `.local-build`; no system installation was performed. Build uses Apple Clang, Meson, Ninja, Boost 1.88, and `-DMETA_NO_STD_FORWARD_DECLARATIONS` for bundled range-v3/libc++ compatibility.

```sh
# From the nested revbayes repository; regenerate lists if adding C++ files.
sh projects/meson/generate_sources.sh core
sh projects/meson/generate_sources.sh revlanguage
.local-build/run-build.sh
python3 tests/test_DDD_benchmark/run.py
# Full synthetic integration:
python3 tests/test_DDFBD/run_integration.py .local-build/revbayes-build/rb --generations 2000 --output tests/test_DDFBD/results/integration.json
python3 tests/test_DDFBD/audit_samples.py tests/test_DDFBD/results/integration --output tests/test_DDFBD/results/sample_audit.json
```

The CLI takes scripts positionally (`rb script.Rev`); it has no `--batch` option. Avoid running upstream `tests/run_integration_tests.sh` in a directory with valuable `data` or `output`: it deletes those directories. Use isolated test copies as in the saved focused-regression harness.

Original user-supplied supplement remains at `../papers/rspb20111439supp.pdf`. Planning sources: the RevBayes specimen tutorial, Stadler et al. (2018) species/range FBD paper, Etienne et al. (2012) and its supplement. Preserve the distinction between DDD's phylogeny-density convention and the labelled oriented fossil-tree density.
