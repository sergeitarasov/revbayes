# Named-species diversity-dependent FBD examples

Run with a RevBayes binary built from this checkout, from this directory:

```sh
../../.local-build/revbayes-build/rb combined.Rev
../../.local-build/revbayes-build/rb morphology.Rev
../../.local-build/revbayes-build/rb dna.Rev
```

These use a tiny **synthetic** dataset, not empirical evidence. Each script runs
2,000 generations to demonstrate the data/model connection; this is not a
convergence recommendation. The common `model.Rev` infers the extended species
tree, budding orientation, node times, fossil B's extinction endpoint, fossil
ages, origin, lambda0, alpha, mu, psi, and character clocks. Rho is fixed at 0.7.
Molecular and morphological partitions can be included independently.

Five fossil records are assigned to species A, A, B, B, and C. They are five
samples of three species, not five independent species. Their geological age
uncertainties are independent uniform intervals in this example. The two living
samples are A_now and C_now. The initial species tree has B extinct at age 0.4;
its younger fossil lies near age 1.0. The extinction-time proposal uses the
current sampled age of that youngest fossil as its upper bound. Named proposal
bounds are included explicitly in the model so MCMC cloning preserves them.
The observation-tree function places its
characters at the fossil ages, not 0.4.

`data/occurrences.tsv` and `data/species.tsv` document these records and their
stable IDs. The small Rev example spells the same values out explicitly; it
does not read those two TSV files automatically.

Morphological data contain only variable binary characters, so the example uses
`coding="variable"`. Match ascertainment to your real matrix. DNA is present for
living samples only; fossil rows have missing states. Do not duplicate a
species-composite morphological row across occurrences. Associate character rows
with actual specimens and ages, or define a separate composite-data model.

The model uses budding speciation and constant mu/psi. Child zero retains the
parent identity. Use `mvOrderedTreeSwap` to explore this ordered state space;
generic topology moves that discard child order are not substitutes. The
species tree's tips are extinction/present endpoints, whereas the observation
tree's tips are actual fossil/extant samples.

Living/extinct status is fixed by the starting tree in this version. A fossil
species believed alive but unsampled at present can have a zero endpoint and be
omitted from `sampledExtant`. Inferring the discrete alive/dead status itself is
not implemented. If all retained species are extinct, generic Newick import may
translate tip ages; set the absolute node ages explicitly before constructing
the distribution (the simulator supplies authoritative age metadata).

The hidden-diversity cap is 64 for these small examples. Repeat with 128/256 and
compare likelihoods and posterior summaries. This cap never limits the
biological birth rate. The priors and origin bounds are illustrative and should
be assessed with the independent simulator in `validation/DDFBD`.

Output goes to `output/<mode>.log`, `.species.trees`, and `.observations.trees`.
The species trees retain child order, so keep that order when reading them.
Scientific use requires longer replicated chains, calibration, numerical
sensitivity checks, and data-appropriate ascertainment assumptions.

Runtime validation: the three 2,000-generation examples, both sampling-conditioning
options, and checkpoint restoration passed with debug cache checks. The saved
sample audit verifies moving trees/parameters, fossil-age consistency, and hidden
cutoff convergence. See [test instructions and results](../../tests/test_DDFBD/README.md).
