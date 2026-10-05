# Diversity-dependent fossilized birth-death models

**New specimen implementation:** see [DDFBDP.md](DDFBDP.md) for
`dnDiversityDependentSpecimenFBD` / `dnDDFBDP`, the ordinary sampled-ancestor tree
prior with exponential decline in total diversity. It includes a morphology-only
example, standard fossil collapse/expand moves, build instructions and validation.

## Earlier named-species implementation


This research branch adds `dnDiversityDependentFBD` (`dnDDFBD`) to RevBayes, with speciation depending on total living diversity, including unobserved lineages. It supports repeated fossil occurrences assigned to named species and joint tree/parameter inference with morphological and molecular data.

- [Biological model and implementation guide](doc/diversity-dependent-fbd/MODEL_GUIDE.md)
- [Implementation status, validation and limitations](doc/diversity-dependent-fbd/IMPLEMENTATION.md)
- [Fossil-model examples](examples/DDFBD)
- [Extant-tree likelihood example and DDD comparison](examples/DDD/README.md)
- [Simulation report (PDF)](output/pdf/ddd_simulation_report.pdf)
- [Presentation simulations, figures and reproducibility](validation/DDD_presentation/README.md)
- [100-species recovery pilot](validation/DDFBD/recovery_100/README.md)

The branch starts from upstream `development` commit `5bd797b1c0d3af4c903343af856d7931b289024f`. Build RevBayes from this branch using the repository's build instructions; the model is not available in an unmodified release. Local build binaries and dependencies are excluded. Saved validation outputs document their original execution environment.

This remains a research implementation. Living/extinct species status is fixed, and successful likelihood and MCMC integration checks do not establish statistical calibration. See the model guide and validation documents before interpreting empirical results.
