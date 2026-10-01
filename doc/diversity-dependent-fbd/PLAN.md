# Diversity-dependent fossilized birth–death in RevBayes

Planning checkpoint: 30 September 2026. Upstream checkout: `development`, commit
`5bd797b1c0d3af4c903343af856d7931b289024f`.

Implementation follow-up: the extended-tree proof is now in
[EXTENDED_TREE_DERIVATION.md](EXTENDED_TREE_DERIVATION.md). Retaining all named
species endpoints makes protection a deterministic schedule and removes the
need for the general identity flags considered below. This document preserves
the planning checkpoint; implementation status and remaining validation are
recorded separately in `IMPLEMENTATION.md`.

## Recommendation and scope

Implement a **species-aware diversity-dependent FBD tree distribution**, coupled
to existing molecular and morphological character likelihoods. The user clarified
that fossil data contain **repeated occurrences assigned to named species**.
Species identity is therefore part of the model, not a later optional feature.

Use the existing FBD speciation/range machinery as the data and extended-tree
reference, and the occurrence birth–death process (OBDP) as the numerical
hidden-lineage reference. Neither existing likelihood can be modified by simply
replacing its scalar birth rate. Develop a new likelihood kernel, with a separate
complete-history simulator/evaluator for verification.

The proposed primary algorithm integrates unobserved diversity using a global
hidden-count filter augmented with species-identity constraints. Prove that the
chosen augmented state is sufficient before optimizing or integrating it into
RevBayes. A scalar hidden-count algorithm has already been checked for the
simpler specimen case, but does **not** complete the requested species model.

This checkpoint delivers planning, source audits, and a small numerical reference.
It does not deliver a registered Rev distribution, a compiled binary, or a
validated end-to-end inference method. Production source files are unchanged.

## 1. Generative model to implement

Start one species at origin age T. Traverse forward from the origin to the
present. Let N(t) count **all living species in the modeled clade**, including
species with no fossils, no character observations, and no sampled living
descendants. It is not cumulative species richness or the reconstructed-tree
lineage count. Interactions with species outside the modeled clade are excluded.

For N >= 1, define the per-species rates

\[
\lambda_N=\lambda_0\exp[-\alpha(N-1)],\qquad \mu_N=\mu,
\]

with lambda0 > 0, alpha >= 0, mu >= 0. Total birth and death hazards are
N lambda_N and N mu. Evaluate lambda at the population size immediately **before**
a birth. Zero is absorbing. Alpha=0 recovers constant-rate diversification.

Use **budding speciation** initially: a parent retains its identity and produces
one new species. This is a proposed biological assumption, not something implied
by diversity dependence. It allows the same species to occur before and after a
speciation event. Symmetric replacement and anagenesis require additional event
types and will be separate extensions. Constant mu is the background extinction
hazard; budding itself does not terminate the parent.

Every living species produces observed fossil occurrences at rate psi, initially
constant and nonremoving. Living species are sampled independently at present
with probability rho. A fossil from an extant species and its extant sample share
one identity. The youngest fossil is not an extinction event. Initial scope uses
one occurrence stream; a distinct anonymous-omega stream is not needed for the
user's stated data.

Alpha has inverse-diversity units; 1/alpha is a decay scale, not a hard capacity.
When lambda0 > mu > 0 and alpha > 0, the continuous approximation to zero net
growth is N*=1+log(lambda0/mu)/alpha. The stochastic process still permits
excursions and eventual extinction; N* is neither a truncation limit nor a stable
population guaranteed by the model.

## 2. Data contract and character observation times

Keep a lossless occurrence table:

| Field | Meaning |
|---|---|
| occurrence_id | Unique specimen/record ID, used to prevent double counting |
| species_id | Stable named-species assignment, shared across repeated records |
| age_min, age_max | Younger/older bounds in one declared time unit |
| count | Multiplicity only for an explicitly aggregated record |
| dating_group | Optional shared geological-age uncertainty group |
| character_sample_id | Links a character row to this specimen, when available |
| source, included | Provenance and analysis inclusion |

Maintain a separate species table with extant status, sampled-at-present status,
and extant character-sample links. Fossil-only species need not be declared
extinct solely because no extant sample is included when rho < 1. Known extinction
or comprehensive extant census information is an additional constraint.

Distinguish three observation designs:

1. All **observed/recovered** occurrence records are supplied: model their full
   Poisson point pattern. This does not assert that all fossils were recovered.
2. Counts in age intervals are supplied: use the corresponding interval-count
   likelihood or a correctly normalized ordered-age augmentation. Do not assign
   all counts to an interval midpoint. Include factorials appropriate to the
   count measure; explicitly labelled specimens use a different measure.
3. Only first/last occurrences or interval presence are supplied: marginalize
   omitted records under that ascertainment scheme. Do not treat omitted records
   as observed absence. Support this after the full-record implementation.

Uncertain occurrence ages become latent variables with the declared dating
likelihood. Shared stratum ages can be jointly latent; independent uniforms are
not an automatic assumption. Preserve internal occurrences even where some ages
are ancillary under constant psi: their number still informs recovery.

**Character observations occur at specimen biological/geological sampling ages,
not species extinction times.** An extended species tree may end a fossil species at its latent death
time. Passing that tree directly to dnPhyloCTMC would place fossil characters at
the wrong time. Build a deterministic observation tree containing the actual
character-bearing specimens and extant samples, with sampled-ancestor encoding
where needed. Characters can be attached to multiple observations of one species
without duplicating its diversification identity.

For a species-composite morphological row with no specimen/time attribution,
require an explicit observation-time/aggregation model. Do not replicate a
composite row at every fossil occurrence. Missing molecular or morphological
partitions use missing states; optional partitions contribute a factor of one.

## 3. Exact likelihood target, independent of numerical algorithm

Let H be the complete budding history: every species birth, death, parent,
identity, and lifetime. Let a be true occurrence ages, F observed occurrence
records, and S present-day sampling indicators. Under a history measure that
identifies the parent and dying species at each event,

use forward time t=0 at the origin and t=T at the present in the following
equation; occurrence ages elsewhere are measured backward from the present.

\[
 p(H,a,S\mid\theta,T)=
 \exp\left\{-\int_0^T N(t)[\lambda_{N(t)}+\mu+\psi]dt\right\}
 \prod_{b\in B(H)}\lambda_{N(b^-)}\;
 \mu^{|D(H)|}\psi^{|F|}
 \rho^{n_s}(1-\rho)^{N(T)-n_s}.
\]

This expression assumes all observed occurrence times are explicitly represented
and labelled; compatible occurrence-to-species assignments are enforced. Births
and deaths are individual events here. A process on counts alone instead has
N lambda_N and N mu event rates and a different reference measure. Hidden-label
permutation factors must be derived when summing histories, not inserted twice.

Multiply by the indicator that each occurrence belongs to a living instance of
its assigned species at that time, and by its dating observation likelihood.
Two observations of the same species require uninterrupted identity, not merely
some genealogical connection. Different named species cannot be the same
identity even when they occur serially along one sampled ancestral path.

Let E denote an oriented extended sampled tree with latent species endpoints and
identity variables. Its likelihood is the pushforward

\[
 p(E,F\mid\theta,T)=
\sum_{\text{compatible hidden structures}}\int_{\{H,a:\Pi(H,a)=E\}}
 p(H,a,S\mid\theta,T)\,
 \mathbf1\{\operatorname{species}(a)=F\}
 p(F_{\rm dating}\mid a)\,dH_{\rm hidden}\,da.
\]

This is an integral over the compatible fiber: retained event times are fixed
coordinates of E, and only hidden coordinates (and unretained uncertain sample
ages) are integrated. It is not an unrestricted continuous-history probability
of exact equality to one tree. Specify the induced reference measure when
deriving the labelled/oriented tree density.

If ages are sampled explicitly in MCMC, omit their integration from a single
likelihood evaluation. Likewise, latent species extinction times and orientation
can be retained in E and integrated by MCMC. The projection and density measure
must be documented precisely before implementing topology moves.

The posterior target is

\[
 p(E,a,\theta,\phi_m,\phi_d\mid X_m,X_d,F)
 \propto L_m(X_m\mid\mathcal T_{obs}(E,a),\phi_m)
 L_d(X_d\mid\mathcal T_{obs}(E,a),\phi_d)
 p(E,a,F\mid\theta)\,p(\theta,\phi_m,\phi_d).
\]

Use the fossil sampling likelihood once. Character models remain ordinary
conditional CTMCs, assuming character evolution does not itself depend on N or
drive diversification, and hidden speciation does not reset character states or
introduce unrepresented clock-rate changes. Unobserved character histories then
integrate out. A relaxed-clock model may instead be explicitly defined on the
observation tree; it must not be described as marginalizing a complete-history
clock with unmodeled resets at hidden births.

## 4. Verified simpler kernel: specimen tree without species constraints

This provides an implementation building block and independent reduction test.
It is not the final species-aware model.

At time t, let k(t) be the number of retained sampled-tree branches, h the hidden
living count, and q_h the unnormalized density of compatible embeddings. Between
observed events, for forward elapsed time,

\[
\dot q_h=(h-1+2k)\lambda_{k+h-1}q_{h-1}
 +(h+1)\mu q_{h+1}
 -(k+h)(\lambda_{k+h}+\mu+\psi+\omega)q_h.
\]

Omega is zero for the primary model; it denotes a separate anonymous observation
stream only in the OBDP reduction. A hidden birth has h possible parents. A birth
on a represented branch has two choices for which descendant retains the
represented continuation, producing the additional 2k. This is an embedding
density operator, **not a conservative or necessarily substochastic CTMC
generator**. Its column sums can be positive.

The raw oriented-tree event maps are:

| Event | State operation | k update |
|---|---|---|
| Particular observed bifurcation | q_h <- lambda_(k+h) q_h | +1, after rate evaluation |
| Sampled-ancestor fossil | q_h <- psi q_h | none |
| Terminal fossil specimen | q_(h+1) <- psi q_h | -1 |
| Separate anonymous occurrence | q_h <- omega(k+h)q_h | none |

A terminal fossil does not kill its species: its unobserved continuation becomes
hidden. Conversely, an **explicit species extinction endpoint** in an extended
tree contributes mu, reduces k by one, and leaves h unchanged. These are distinct
events and must not share one implementation flag.

Initialize k=1, q_0=1 at the origin. At present,

\[
 L_{raw}=\sum_{h\ge0}q_h\rho^{k}(1-\rho)^h.
\]

No binomial sampling coefficient is added to this embedded-tree boundary.
Implement rho=0/1 and fossil-only cases explicitly. Density conversion to the
chosen labelled-tree measure and any ascertainment normalizer come separately.

## 5. Species-aware filter: required derivation

For a fixed orientation and identity configuration, carry q_(h,z), where z
encodes the relevant identity state on retained branches and unresolved
identity-change requirements. States differing in future compatibility cannot
be merged just because they have the same N.

Enumerate elementary budding/death events on a representative complete labelled
history and aggregate only transitions with identical future constraints:

- Hidden species birth: h choices, rate lambda_N each; h increases.
- Retained species buds, and the parent continues the retained path: rate
  lambda_N; the new hidden daughter increases h and the retained identity stays.
- Retained species buds, and the daughter continues the retained path: rate
  lambda_N; the hidden parent increases h and the retained identity changes.
  Permit this only when compatible with z and update z.
- Hidden death: h choices at mu; h decreases.
- Retained death: forbidden except at its explicit extinction endpoint.
- Assigned occurrence: multiply by psi on its designated living species and
  apply the identity constraint. There is no N multiplier.
- Observed budding with two retained descendants: multiply by lambda_N, carry
  the parent identity down its A descendant and create a new identity on D.
  Sum or sample alternative orientations rather than assuming an exchangeable
  child order.

Keep the full no-event diagonal -N(lambda_N+mu+psi). A forbidden positive event
does not remove its hazard from that diagonal.

Within a same-species protected interval, only the parent-continuation route is
allowed. If s of k paths are protected and all other identity conditions are
absent, the hidden-birth coefficient becomes (h+2k-s)lambda_N. This local rule
alone is insufficient for a whole tree: serial **different** species also need
a birth that changes the retained identity in the separating gap.

For an extended tree, protection continues from the oldest assigned occurrence
through the species' latent extinction endpoint, including the segment after its
youngest fossil. If the species is alive at present, protect its identity to the
present instead. The species-aware present boundary weights every retained and
hidden living species: rho for each sampled living species and (1-rho) for each
unsampled living species, including retained fossil-bearing species alive but
unsampled at present. It is not automatically the specimen kernel's rho^k.

A useful exact local construction for one such gap is a two-state flag z=0/1
(identity unchanged/changed). With k=1:

- Within z=0: upward coefficient (h+1)lambda_(h+1).
- From (h,0) to (h+1,1): coefficient lambda_(h+1).
- Within z=1: upward coefficient (h+2)lambda_(h+1).
- Both blocks have hidden death h mu and the original no-event diagonal.

Start in flag zero and keep only flag one at the gap end if distinct species are
required. Summing the flags recovers the unrestricted kernel; flag zero alone
recovers the protected kernel. This avoids cancellation from subtracting two
similar likelihoods. Generalizing overlapping obligations and orientation
updates is the main outstanding mathematical task.

**Proof gate:** show that histories merged into (h,z) have the same aggregate
transition rates and future observation compatibility. Derive initialization,
event maps, and terminal identities for arbitrary valid oriented extended trees.
If a proposed compact z fails this test, expand it or explicitly augment the
unobserved identity-changing births. Do not ship a count-only approximation.
Worst-case state size may grow with simultaneous identity obligations; benchmark
this before promising scalability.

Full hidden-tree MCMC is a fallback, not the initial production choice. It has a
simple history density but needs reversible-jump hidden-subtree proposals,
Jacobian/Hastings proofs, consistent labels, and substantial mixing validation.
An unbiased particle likelihood is another research route; a noisy plug-in
likelihood does not automatically define valid MCMC.

## 6. Tree measure and conditioning

The Etienne supplement computes branching-time densities with its own sampling
and crown conventions. Its k lambda event factor and fixed-size missing-species
normalization cannot be copied into a likelihood for a particular labelled
tree with Bernoulli rho sampling.

Current ordinary FBD applies `(tips-SA-1)log(2)-log(tips!)` to its raw kernel.
Current OBDP inherits a different tree-shape factor. Named-species budding trees
add biologically meaningful orientations and identity restrictions, so neither
factor may be adopted without derivation. Factors varying with ancestor status
or orientation affect MCMC probabilities even when parameter independent.

First implement an origin-started unconditioned generative likelihood with
`condition="none"`. Then implement explicit, separately named ascertainment
events. Starting at an origin is not itself conditioning on survival.

For at least one sampled extant species, propagate the ordinary count CTMC with
birth N lambda_N, death N mu, and no fossil killing. If p_n(T) is its endpoint
law, divide by `1-sum_n p_n(T)(1-rho)^n`. For at least one observation of any
kind, add killing N psi to the no-observation count process before applying the
same endpoint weight. These conditions are different and must reflect dataset
selection.

Crown-started inference requires its own derivation. For both initial daughter
groups to have sampled extant descendants, a two-color population (n1,n2) evolves
with birth ni lambda_(n1+n2), death ni mu. Its endpoint weight is
`[1-(1-rho)^n1][1-(1-rho)^n2]`. The two groups are coupled; squaring one-lineage
survival is incorrect. Defer crown-start and conditioning on observed species
number until derived and independently checked.

## 7. Numerical implementation

Use a sparse/block-sparse matrix-exponential action or a stable adaptive ODE
solver for fixed-event intervals. A tridiagonal application costs O(H); species
state adds a factor depending on z. Do not copy dense exponentiation as the
production scaling strategy without benchmarks. Rescale positive vectors and
accumulate log scales. Detect zero, negative, NaN, and overflow cases explicitly.

Truncate hidden diversity at H with the original diagonal retained and outward
transitions discarded. This removes complete-history paths; it is not a
biological ceiling or reflecting boundary. For nonnegative kernels, increasing
H should increase the raw path sum. Check likelihood differences at H, 2H, 4H
over representative parameter/posterior states and solver tolerances.

Forward boundary values, normalized vector sums, and final boundary mass alone
are not rigorous error certificates for this weighted operator. Validate both
numerator and any conditioning denominator. A finite-cap calculation remains an
approximation to the infinite-state likelihood. Numerical failures should be
reported distinctly from biological zero density, not silently rejected as if
they were model constraints.

During MCMC use a fixed cap shown adequate in pilots, or a deterministic
state-based convergence procedure. Do not let past chain states determine the
current numerical target. Recompute the full schedule initially: alpha changes
all rates, and age/topology changes may reorder global events. Add caches only
after rejection restore, clone, and checkpoint tests pass.

Treat artificial zero-length sampled-ancestor nodes as one sampling event, not
a birth. Define deterministic ordering for exact-age ties, while recognizing
that independent continuous-time event coincidences normally have zero measure.
Age-interval counts need interval likelihoods rather than arbitrary tie ordering.

Forward/backward messages allow posterior diversity reconstruction by
`P(h,z at t | data,theta) proportional to q_(h,z)(t)b_(h,z)(t)`;
sum over z and report N=k+h. This is a smoothed uncertainty distribution, not a
single plug-in diversity trajectory for calculating lambda.

## 8. RevBayes integration design

Use the following existing code as references, not as unexamined base likelihoods:

| Existing file under src/ | Reuse or audit |
|---|---|
| core/distributions/phylogenetics/tree/birthdeath/FossilizedBirthDeathSpeciationProcess.{h,cpp} | Species-aware extended tree and orientations |
| core/distributions/phylogenetics/AbstractFossilizedBirthDeathRangeProcess.{h,cpp} | Repeated occurrence intervals/counts, range bookkeeping |
| core/datatypes/phylogenetics/Taxon.h | Existing occurrence map; assess lossless record-ID adapter |
| core/functions/phylogenetics/ComputeLikelihoodsLtMt.{h,cpp} | Hidden-count propagation and event architecture |
| core/distributions/phylogenetics/tree/birthdeath/OccurrenceBirthDeathProcess.{h,cpp} | Distribution wiring, not species semantics |
| core/distributions/phylogenetics/tree/birthdeath/BirthDeathSamplingTreatmentProcess.cpp | FBD reduction and sampled-ancestor measure |
| revlanguage/workspace/RbRegister_Dist.cpp | Distribution registration |

Proposed new components (names are provisional):

1. `DiversityDependentFbdLikelihood.{h,cpp}`: pure numerical rate law, identity
   state, event schedule, propagation, boundaries, conditioning, diagnostics.
2. `DiversityDependentFossilizedBirthDeathProcess.{h,cpp}`: distribution adapter
   over an extended species tree plus registered occurrence-age/identity DAG
   inputs. Prefer existing TimeTree infrastructure where its semantics suffice;
   introduce a typed history value only if required state cannot be safely
   represented. Own the configuration rather than storing references to temporary
   strings.
3. Rev wrapper `Dist_diversityDependentFBD.{h,cpp}` and help page
   `help/md/dnDiversityDependentFBD.md`.
4. A typed deterministic observation-tree function, linked to extended-tree
   changes, occurrence ages, and character sample assignments. Validate branch
   times and molecular/morphological likelihood invariance to unobserved tips.
5. Species-orientation and identity-augmentation proposals where required;
   existing topology/time proposals only after auditing biological child order.
6. Independent complete-history Gillespie simulator and full-history density
   evaluator, plus small identity-labelled enumeration fixtures.

Provisional public arguments: originAge, lambda0, alpha, mu, psi, rho, taxa,
speciesOccurrences, occurrenceAges, initialExtendedTree, condition,
maxHiddenLineages, numericalTolerance. Start with budding only and complete
observed occurrence records. This is an interface sketch, not runnable Rev code.

Register every stochastic dependency using addParameter and correct pointer
swapping. Implement const-correct overrides, cloning, touch/keep/restore, and
initialization. MCMC initialization may use a compatible constructed tree; that
is not an exact DD prior draw. Redraw/simulation must either implement the stated
distribution or clearly reject unsupported simulation requests.

### Findings that need attention before reusing existing classes

- FBDSP's `computeLnProbabilityDivergenceTimes()` is non-const while the base
  virtual is const; its header even notes the failed override. Static inspection
  identifies a dispatch mismatch. Reproduce its effect in a focused compiled
  test, repair separately if confirmed, and do not use unverified current output
  as a gold standard for the species-range likelihood.
- Extended-tree child zero carries ancestral identity; automatic child sorting
  can change the model. Tip ages encode species extinction, not sampling age.
- Existing OBDP conditioning and analytic shortcuts assume linear birth–death.
- OBDP's occurrence-age dependency registration and cutoff diagnostics need
  independent review rather than literal copying.

Current build configuration requires C++23. Source lists are generated by the
CMake/Meson helper scripts; use the supported build workflow. Add focused
`tests/test_DDFBD` integration cases and run the existing FBD suite too. No build
has been attempted in this planning checkpoint.

## 9. MCMC plan

Sample the extended species topology, node ages, origin, latent species endpoints,
occurrence ages, required identity/orientation variables, diversification rates,
and character/clock parameters. Hidden anonymous lineages are marginalized in
the primary likelihood; they do not require an RJ move if the marked filter is
proved sufficient.

Use scale moves for positive rates and origin offsets; nonnegative additive or
transformed moves for alpha. A scale move started at alpha=0 remains at zero.
Compare alpha=0 as an explicitly fixed nested model or a carefully defined
mixture, rather than interpreting the continuous prior's boundary as a model
probability. Choose rate and alpha priors in the dataset's units and assess them
with prior predictive simulation.

Use constrained node/endpoint slides, topology changes, orientation flips,
occurrence-age updates within geological support, and any augmented identity
birth updates. Ensure they connect alternative valid species placements;
ordinary specimen collapse/expand moves cannot be assumed to cover all range
states. Validate each move's reverse support and Hastings ratio.

If living/extinct status is unknown for a fossil-bearing species, support the
point mass at a present-day living endpoint as well as continuous extinction
times. This requires a dedicated status-changing proposal; a continuous endpoint
slide alone cannot cross that measure boundary. A first dataset with externally
known living/extinct status can fix those indicators explicitly.

Attach standard molecular and morphology CTMCs to the observation tree. Preserve
the actual morphology ascertainment correction and partition clocks. Begin with
small fixed-tree parameter inference, then free times, then full topology and
fossil-age inference. Monitor separate likelihood components and numerical
diagnostics along with rates, species lifetimes, tree summaries, and N(t).

## 10. Ordered work packages and acceptance gates

| Stage | Deliverable | Gate before proceeding |
|---|---|---|
| A. Model contract | Budding data schema, ascertainment, exact history measure | Same/different species constraints explicit; specimen times mapped |
| B. Identity likelihood | Proof and reference q_(h,z), or explicit identity augmentation | Tiny labelled-history sums and identity counterexamples agree |
| C. Numerical engine | Sparse deterministic likelihood, diagnostics, simulator | Constant-rate species FBD, DD reductions, cap/tolerance checks |
| D. Rev integration | Distribution, observation-tree function, moves, help | Builds; DAG dependency, rejection, clone, checkpoint tests |
| E. Total-evidence inference | Executable morphology-only, DNA-only, combined examples with repeated species records | Free-tree and parameter MCMC runs; age/orientation mixing checked |
| F. Scientific validation | Independent simulation study and documentation | Calibration, null false positives, recovery, numerical robustness |

Stages A–B are the next production work. Do not claim completion at a
specimen-only distribution or a fixed-tree optimizer.

Required checks include:

- Alpha=0 against the 2018 species-range likelihood under matched budding,
  sampling, orientation, and conditioning conventions; audit its implementation
  before trusting it as a reference.
- Remove identity constraints to recover ordinary FBD/OBDP; include fossil-tip
  versus sampled-ancestor states and absolute density constants.
- No fossil sampling to recover the extant DD process under matching sampling
  and origin/crown conventions, not the supplement's unmatched normalizer.
- Mu=0, rho=1, no fossils: analytic complete DD pure-birth density.
- Multiple occurrences of one species: one identity, no duplicated lineage
  count; reject histories replacing that species between its samples.
- Different serial species: zero likelihood when no compatible identity-changing
  event is possible. Orientation flips and gaps must affect the right terms.
- Fossil sample time versus extinction time: moving a death endpoint must not
  move a character observation or change its CTMC age spuriously.
- Explicitly enumerate tiny complete labelled histories; compare aggregated
  transitions with the proposed marked-state recursion.
- Independent Gillespie simulation with total hazard N[lambda_N+mu+psi], then
  sample/prune/map identities. Use this to test inference, not a simulator built
  from the same filter equations alone.
- Simulation-based calibration with adequate effective sample sizes, repeated
  chains, topological/orientation checks, and a designed grid of alpha, mu,
  recovery, rho, ages, fossil counts, and time uncertainty.
- Stress alpha near zero, weak fossil recovery, incomplete extant sampling,
  long origins, repeated species ranges, and many identity obligations.

An identifiable effective density dependence is not automatically evidence for
competition. Assess confounding among alpha, lambda0, mu, psi, rho, origin, and
sampling heterogeneity; use recovery and sensitivity results to bound claims.

## 11. Work completed at this checkpoint

Three parallel reviews covered likelihood theory, RevBayes architecture, and
literature/data semantics. Their notes are in `reviews/`; their later
named-species addenda supersede preliminary specimen-first recommendations.
This consolidated plan is the operative scope.

`reference_count_likelihood.py` is a standard-library-only planning prototype.
It uses small fixed-step RK4 for tiny examples, not a production numerical
solver. Reproduce its output from the repository root:

```sh
python3 doc/diversity-dependent-fbd/reference_count_likelihood.py
```

The saved `reference_checks.json` reports 11 passing checks:

- Six constant-rate FBD fixtures, maximum absolute log-density discrepancy
  approximately 5.2e-12.
- DD forward/adjoint agreement approximately 1.8e-15; step-halving discrepancy
  approximately 6.7e-12 on the same fixture.
- Complete DD pure birth discrepancy approximately 3.2e-14.
- Hidden caps 4, 8, 16, 32, 64: final log-density difference approximately 2.1e-8.
- Alpha approaching zero: discrepancy approximately 2.0e-8 at alpha=1e-8.
- A restricted named-species case with the founder's identity protected for the
  full interval: discrepancy approximately 1.1e-13 against a separate scalar
  constant-rate formula.

These are local algebra/numerical checks. They do not validate general named
species filtering, labelled-tree normalization, RevBayes MCMC, performance, or
scientific identifiability. Those are explicit acceptance gates above.

## Sources and reading record

1. [RevBayes specimen total-evidence tutorial](https://revbayes.github.io/tutorials/fbd/fbd_specimen).
   Read the workflow and its distinction between specimen dating and ranges.
2. Stadler et al. (2018), [The fossilized birth-death model for the analysis of
   stratigraphic range data under different speciation modes](https://doi.org/10.1016/j.jtbi.2018.03.005).
   Publisher access failed; read the [author manuscript](https://arxiv.org/html/1706.10106).
3. Etienne et al. (2012), [Diversity-dependence brings molecular phylogenies closer
   to agreement with the fossil record](https://pmc.ncbi.nlm.nih.gov/articles/PMC3282358/),
   doi:10.1098/rspb.2011.1439. Direct PMC access was intermittent; the theory
   review obtained readable article content and inspected the local supplement.
4. User-supplied [supplement](../../../papers/rspb20111439supp.pdf) is preserved
   at the workspace's `papers/rspb20111439supp.pdf`.
   The actual file supplied is named without “2”, is 10 pages, and contains
   Appendices S1–S2. Text and equation pages were inspected. SHA256:
   `66d0168d1e37b751e18ca82389308b8c1f59fc7b6a3e772cf4930613e26e8cb4`.
5. Andréoletti et al. (2022), [The Occurrence Birth–Death Process for Combined-Evidence
   Analysis in Macroevolution and Epidemiology](https://doi.org/10.1093/sysbio/syac037),
   [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC9558841/).
   Read methods and validation and audited current implementation; its Dryad
   supplementary derivation was not retrieved in this checkpoint.
6. [RevBayes source](https://github.com/revbayes/revbayes/tree/5bd797b1c0d3af4c903343af856d7931b289024f),
   inspected locally at the recorded commit. No novelty claim is made from this
   targeted literature review.
