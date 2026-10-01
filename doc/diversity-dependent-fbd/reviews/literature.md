# DD-FBD literature and data-model planning notes

> Review snapshot. The final named-species scope is in ../PLAN.md. The later
> named-species addendum supersedes preliminary specimen-first recommendations.

## Scope and source access

Read the full RevBayes specimen tutorial, the accessible arXiv manuscript corresponding to the requested 2018 paper, and the methods/validation sections of the 2022 OBDP article (PMC and author PDF). ScienceDirect itself returned an internal/access error; its article identity is independently corroborated by arXiv. No claim of exhaustive literature coverage or novelty is justified. The Etienne 2012 supplement is assigned to another agent.

## Requested 2018 FBD paper

Stadler, Gavryushkina, Warnock, Drummond & Heath (2018), **The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes**, Journal of Theoretical Biology. Full manuscript: https://arxiv.org/html/1706.10106 ; record https://arxiv.org/abs/1706.10106 ; requested publisher https://www.sciencedirect.com/science/article/pii/S002251931830119X . The arXiv abstract title says “concepts”; current HTML says “modes.”

This is specifically a model linking multiple fossils to the same species through stratigraphic ranges, with budding, bifurcating, and anagenetic speciation. Budding preserves the ancestral species; bifurcation replaces it with two; anagenesis replaces it without increasing contemporaneous lineage count. It derives densities for sampled trees/ranges, extends some densities to unknown extinction times, and marginalizes unrecorded within-range fossils or interval-presence data. A range constrains continuity of species identity along a lineage; it is not an age-uncertainty interval for a single specimen. Closed-form branch factors here assume lineage-independent rates and cannot simply receive a plug-in reconstructed-tree diversity function.

## Requested specimen tutorial

https://revbayes.github.io/tutorials/fbd/fbd_specimen (Heath, Wright & Walker; updated March 26, 2025).

The tutorial composes a tree distribution with separate DNA, morphology, and fossil-age components. Its distribution is `dnBDSTP` with removal `r=0`, origin conditioning, speciation/extinction/fossil sampling rates, extant sampling probability and taxa. `addMissingTaxa` harmonizes partitions. Fossil-age intervals describe individual dated specimens; the tutorial explicitly warns against using species FAD–LAD ranges this way. It demonstrates sampled-ancestor proposals, fossil-time moves, topology/node-time moves, and clade constraints for poorly characterized fossils. Reusable module names: `mcmc_CEFBDP_Specimens.Rev`, `model_FBDP.Rev`, `model_UExp.Rev`, `model_GTRG.Rev`, `model_Morph.Rev`. Data: `bears_taxa.tsv`, `bears_cytb.nex`, `bears_morphology.nex`. Morphology uses variable-character ascertainment; preserve the actual dataset's ascertainment scheme. Existing substitution/clock modules can remain separate from the new diversification calculation.

## Closest existing framework: occurrence birth–death process

Andréoletti et al. (2022), **The Occurrence Birth–Death Process for Combined-Evidence Analysis in Macroevolution and Epidemiology**, Systematic Biology 71:1440–1452, https://doi.org/10.1093/sysbio/syac037 . Full text https://pmc.ncbi.nlm.nih.gov/articles/PMC9558841/ ; author PDF https://www.normalesup.org/~manceau/fr/files/AndreolettiZwaansEtAl_2022.pdf ; supplement/data https://doi.org/10.5061/dryad.p8cz8w9rq .

This is already implemented in RevBayes. It combines a sampled tree with occurrences lacking character data, using separate tree-sampling rate psi, occurrence rate omega, extant probability rho and removal probability r. Its forward/backward quantities index hidden lineage count i, with total count k(t)+i. Numerical master equations give the likelihood and smoothed diversity. Rates are piecewise constant through time, not diversity dependent. A finite lineage truncation approximates the infinite state space. Validation includes special-case likelihood comparisons and simulation-based calibration. This is the strongest existing implementation scaffold for our extension, especially if “fossil occurrence data” means unplaced events in addition to fossils represented in the tree.

## Other relevant primary sources

Etienne et al. (2012): https://pmc.ncbi.nlm.nih.gov/articles/PMC3282358/ . Hidden extinct/unsampled contemporaneous species affect rates; counts cannot be read from the reconstructed tree alone.

Etienne & Haegeman (2019), **Additional Analytical Support for a New Method to Compute the Likelihood of Diversification Models**, https://pmc.ncbi.nlm.nih.gov/articles/PMC6976549/ . Useful analytical audit of the 2012 count-filter framework and its missing-species interpretation; should inform independent derivation review.

Etienne, Pigot & Phillimore (2016), **How reliably can we infer diversity-dependent diversification from phylogenies?**, https://doi.org/10.1111/2041-210X.12565 . Evaluates estimator behavior under alternative conditioning. Supports a deliberate test matrix for conditioning, weak signal, extinction, and null false positives, rather than treating successful numerical optimization as identifiability.

## Proposed design decisions (our synthesis; not quoted literature results)

1. Define an explicit generative model before coding. Per-lineage lambda(N)=lambda0 exp[-alpha(N-1)], alpha>=0, constant mu, fossil sampling psi, Bernoulli extant sampling rho. N is total living diversity at the relevant instant, including hidden lineages. alpha=0 is the exact constant-rate boundary. alpha has inverse-lineage units; 1/alpha is a decay scale, not generally an equilibrium carrying capacity.
2. SUPERSEDED by the user's clarification: repeated occurrences assigned to named species are core scope. The required deliverable is a species-aware DD range/occurrence model, not a specimen-only or anonymous-occurrence MVP. See the concrete requirements below. An anonymous OBDP extension is useful only as a nested validation problem.
3. Partition fossil observation records into disjoint samples. A fossil represented in the tree cannot also be counted as an independent omega occurrence unless a deliberately modeled observation mechanism says so. Preserve occurrence IDs to detect duplicate use.
4. Distinguish three APIs/data meanings: uncertain specimen age; anonymous occurrence event (integrated over all potential lineage placements); repeated identified species occurrences/range. The third requires species-identity constraints and speciation-mode assumptions beyond an ordinary sampled-ancestor tree. It is now explicitly the core model, not an optional future extension or an implicit interpretation of an age interval.
5. A single clade-wide hidden count can couple all tree branches. Existing FBD per-branch multiplication relies on independence that DD removes. The numerical likelihood must propagate a global event schedule/count vector and integrate hidden diversity. Molecular/morphological pruning remains conditional on the sampled tree and needs no coupling to unsampled lineages under the stated independent character model.
6. Prefer an origin-conditioned first implementation. Crown conditioning couples the two descendant subclades through N; two independent clade survival factors are not automatically valid. Explicitly derive each ascertainment/conditioning normalizer and test it.
7. Give finite-state cutoff and ODE tolerances first-class diagnostics. Verify convergence over increasing cutoffs; never silently interpret a numerical cap as a biological maximum or normalize away escaped mass. Alpha=0 is a particularly demanding truncation test.
8. Reuse existing sampled-ancestor compatible tree moves and ordinary parameter moves where structurally valid. For alpha near zero, use an additive/nonnegative-compatible move and a prior with scientifically chosen scale; multiplicative moves initialized exactly at zero cannot explore positive values.
9. Acceptance gates: independent forward/backward agreement; constant-rate FBD/OBDP reduction with matched constants/conditioning; no-fossil DD reduction; pure-birth exact checks; full-data simulation recovery; simulation-based calibration; cutoff/tolerance robustness; joint tree/parameter MCMC and fossil-age/sampled-ancestor mixing. Compare null alpha=0 and DD models with appropriate boundary-aware model comparison rather than a naive interior chi-square assumption.
10. Empirical demonstrations must show uncertainty in N(t), lambda0, alpha, mu, psi (and omega if present), tree times and topology. Sparse fossils, unknown sampling, and uncertain origin can confound inferred feedback. A fitted diversity dependence is not proof of competition.

## Important current access limitation

The main 2022 OBDP paper identifies its master-equation derivation in Dryad supplementary appendices. I inspected its definitions and numerical implementation, but did not retrieve those appendices. The other agents' direct C++ audit and DD-supplement derivation should determine the exact event factors and indexing. Do not claim those factors have been fully independently verified from the OBDP supplement on this literature pass.

## Core scope after user clarification: repeated occurrences of named species

### Default biological model

Recommend **budding speciation**, with no anagenesis and no symmetric speciation in the first complete implementation. This is an explicit modeling assumption, not something inferable just from occurrence names. At birth, one daughter retains the parent's identity; the other receives a new identity. The parent may persist alongside descendants. Each contemporaneous species is counted once in N, regardless of the number of sampled specimens or subdivisions of its path in the tree. Fossil collection does not terminate the species. Constant extinction acts independently on each extant species. The requested exponential function applies to the total birth rate per species.

Mixed speciation requires extra beta and anagenesis parameters and event types. Symmetric birth terminates the parent species and produces two new identities while changing N by +1. Anagenesis changes identity but leaves N unchanged. Both are prohibited inside an interval known to belong continuously to one species; budding is allowed there only if the ancestral identity continues along the protected path. Defer mixed modes explicitly, while completing named-species budding inference.

### Exact full-history target (proposed mathematical specification)

Let H be the complete **oriented, species-labeled** history from one origin lineage at forward time 0 to present T. Every species has a birth, extinction-or-present time, and a parent identity; the initial species is a special root. H includes entirely unobserved species. Let N_H(t) be the number alive. Under a fixed labeled/oriented reference measure, before conditioning and any label-symmetry conversion,

    p(H, fossil times, extant sampling | theta)
      = exp[-integral_0^T N_H(t){lambda(N_H(t))+mu+psi} dt]
        * product_births lambda(N_H(t-))
        * mu^(number of extinctions)
        * psi^(number of fossil observations)
        * rho^(number of sampled extant species)
        * (1-rho)^(number of unsampled extant species).

There is no factor N on an individual labeled event when the affected lineage is specified. N lambda(N) appears when events are represented only as unmarked count changes and lineage choices are summed. The correct combinatorial conversion must be fixed once and tested against the existing constant-rate range likelihood; do not silently interchange these reference measures.

Multiply this density by: (a) compatibility of every named occurrence with the lifetime of that exact species; (b) dated fossil observation densities; (c) molecular/morphological likelihood on the sampled character genealogy induced by H; (d) priors. Integrate over all unobserved histories, latent fossil ages, species-birth/extinction times and orientations not retained in the MCMC state. Normalize for exactly the chosen sampling/survival condition. This defines an exact target even if a practical marginalization algorithm is initially incomplete.

For species s with fossil ages a_sj, require every a_sj inside its lifetime. Its FAD and LAD bound observed persistence, not true origination/extinction. A named extant species persists through the present and contributes rho once. All its fossils refer to the same species path. Species identities cannot disappear and then reappear. A single occurrence gives a zero-length observed range but does not imply an instantaneous species lifetime. Distinct species names must not be treated as interchangeable marks on one uninterrupted identity path.

### Can a simple marked-count extension work?

**A scalar hidden count plus only the number of active observed ranges is not generally sufficient.** Two configurations can have identical total count and identical number of protected paths, but differ in which identities must connect and whether a required identity change has occurred. Future named observations distinguish them, so they cannot be merged without proof of lumpability.

A tractable extension may nevertheless be possible for a fixed augmented species-aware skeleton. The state would be `(hidden_count, identity_state)` rather than only `hidden_count`. The finite identity state must encode sampled-path orientations and outstanding continuity/distinctness obligations. Use the following event rules as derivation requirements:

- Within a named observed range, the visible path must retain its identity. Hidden budding is allowed if the parent continues visibly and the new species is unobserved. A hidden event in which the visible path switches to the new bud is forbidden there.
- Outside protected ranges, both parent-continuation and new-bud-continuation may be possible. They have different identity-state effects and must not be combined if a future range can distinguish them.
- When two distinct named ranges occur serially on one sampled path, at least one identity-changing hidden budding event must occur in the gap. A binary “change has occurred” flag for that gap is the minimal additional memory in a simple case. Multiple simultaneous unresolved paths require a vector or other joint identity state, potentially growing exponentially.
- At a visible branching event crossed by a named range, that species must continue on exactly one daughter, which fixes its A/D orientation. Where neither daughter is identity-constrained, sample or sum the compatible orientations with the correct measure factors.
- Multiple named ranges on separate paths with overlapping ages are distinct contemporaneous species. Repeated observations on one protected path are not additional lineages.
- Fossil observations on a specified named path receive a per-lineage sampling factor psi; replacing this by N psi would instead integrate over anonymous lineage identity.

An alternative is explicit augmentation of required identity-changing birth times and the necessary oriented ancestral skeleton, then integrate only genuinely exchangeable hidden side histories by count. Another exact reference approach augments the full history and uses transdimensional MCMC. Neither should be advertised as a finished practical solution before event-kernel derivation, normalization, and reversible-move validation.

The 2018 paper provides a decisive constant-rate audit: its unrestricted branch equation has two possible hidden-birth continuations while its within-species equation has one, and its sampled-tree construction includes additional latent births between distinct serial species ranges. Under diversity dependence those hidden births affect all concurrent rates through N. Consequently the paper's final product of independent q-functions is not a valid plug-in shortcut. The count-plus-identity algorithm must reproduce its range likelihood at alpha=0 under matching orientation and conditioning conventions.

### Required implementation artifacts and acceptance gates

- Typed input records: immutable occurrence ID, species ID, age observation/interval, sampling provenance, character association. Store species metadata and living/dead status separately. Fossil age uncertainty is not the species FAD–LAD span.
- Explicit species-to-skeleton mapping, observed range endpoints, allowed orientations, and latent identity-transition constraints; validity checks should explain violations.
- A character-data policy: if morphology is a species composite, specify its attachment time/model; if multiple specimens have characters, retain their serial observations on the same species lineage. Do not copy a species composite onto every occurrence, which would replicate evidence.
- Tree/range moves that preserve or jointly change species mappings, ancestral status, orientation, latent origination/extinction times, and fossil ages. Standard sampled-ancestor moves alone do not explore this richer state correctly.
- Independent complete-history simulator preserving names and repeated occurrences, plus a tiny-history enumerator or explicit-history likelihood oracle.
- Tests specifically separating same-species serial fossils from different-species serial fossils; budding within a persisting ancestral range; overlapping ancestral and descendant ranges; extant species with multiple historical samples; zero-length singleton ranges; impossible recurrence; different but count-identical identity states.
- Constant-rate range-likelihood reduction, collapse to ordinary FBD when identity information is deliberately discarded with the correct marginalization, and count-cutoff convergence. Only after those pass proceed to repeated-species total-evidence MCMC calibration.
