## name
dnDiversityDependentFBD
## title
Diversity-dependent fossilized birth-death with named species occurrences
## description
A budding-speciation process on an oriented extended species tree. The per-species
speciation rate is `lambda0 * exp(-alpha * (N - 1))`, where N counts all living
species, including unsampled species. Extinction mu and fossil recovery psi are
constant. Hidden species histories are numerically marginalized.
## details
The alias is `dnDDFBD`. `occurrenceSpecies` and `occurrenceAges` contain one entry
per recovered fossil record, with repeated species names allowed. Ages can be
stochastic, for example with geological dating likelihoods or uniform bounds.
The supplied records must be the complete observed point pattern. Supplying only
first/last appearances while ignoring known interior counts is a different model
and is not supported by this first implementation.

`initialTree` is a binary extended tree with one tip per observed species. Child
zero is the ancestral continuation at every budding event. A positive tip age is
the species' latent extinction time, not its youngest fossil age. A zero tip age
means the species is alive. Living/extinct status is fixed by the initial tree.
`sampledExtant` lists species sampled alive at present; living fossil-bearing
species absent from that list receive a factor (1-rho). Every retained species
must have a fossil record or an extant sample. The species tree may include an
initial stem above its root, bounded by `originAge`.

Use `fnSpeciesObservationTree` to attach character data at specimen ages. Do not
attach fossil characters directly to latent extinction endpoints. Include each
character-bearing fossil exactly once among the occurrence records. DNA and
morphology are optional separate `dnPhyloCTMC` components.

Use `mvOrderedTreeSwap` for changes to topology and budding orientation. It
preserves child slots and reverses exactly. Other generic topology moves have not
been audited for the oriented state space. Standard node/root time moves and
bounded fossil-tip moves can change node ages and extinction endpoints.

The density retains orientations explicitly and includes the species-label
factor 1/n!, without an additional factor 2^(n-1). `condition="none"` gives the
unconditioned origin-started density. `"sampling"` conditions on at least one
fossil or extant sample; `"sampledExtant"` conditions on at least one extant
sample. Crown conditioning and conditioning on observed species number are not
implemented.

`maxHiddenLineages` is a numerical cutoff, not a biological carrying capacity.
Increase it to check likelihood and posterior sensitivity. Outward transitions
are killed while preserving the original diagonal hazard. Conditioning uses a
separate total-count cutoff equal to maxHiddenLineages plus the number of retained
species. `numericalTolerance` controls the positive matrix-exponential series,
not truncation in hidden diversity. Rates and all event times are validated at
each evaluation; numerical convergence failures are reported as errors.

MCMC initialization uses the supplied tree. General prior redraw conditional on
the specified named records is not implemented and throws rather than pretending
that a starting tree is a random draw. The independent complete-history simulator
is in `validation/DDFBD`.

Current limitations: budding only; constant mu and psi; fixed living/extinct
status; exact individual occurrence ages, optionally sampled in the DAG; no
range-only ascertainment, taxonomic assignment uncertainty, or interval-count
likelihood. Broad simulation-based calibration remains necessary before
scientific use.
## example
initial <- readTrees("extended_species.tre")[1]
tree ~ dnDDFBD(originAge=7, lambda0=0.5, alpha=0.2, mu=0.2,
    psi=0.3, rho=0.7,
    occurrenceSpecies=v("A","A","B","B","C"),
    occurrenceAges=v(4.0,2.0,2.0,1.0,1.0),
    sampledExtant=v("A","C"), initialTree=initial,
    condition="none", maxHiddenLineages=128)
## see_also
fnSpeciesObservationTree
mvOrderedTreeSwap
dnFossilizedBirthDeathSpeciation
dnOccurrenceBirthDeath
## references
- citation: Stadler et al. (2018). The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes.
  doi: 10.1016/j.jtbi.2018.03.005
- citation: Etienne et al. (2012). Diversity-dependence brings molecular phylogenies closer to agreement with the fossil record.
  doi: 10.1098/rspb.2011.1439
- citation: Andreoletti et al. (2022). The Occurrence Birth-Death Process for Combined-Evidence Analysis in Macroevolution and Epidemiology.
  doi: 10.1093/sysbio/syac037
