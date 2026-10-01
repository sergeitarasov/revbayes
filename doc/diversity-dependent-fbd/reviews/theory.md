# Independent likelihood review: diversity-dependent fossilized birth-death

> Review snapshot. The final named-species scope is in ../PLAN.md. The later
> species-aware addendum supersedes the initial anonymous-specimen scope.

## Status and scope

This is a mathematical implementation plan, not a validated new likelihood implementation. I independently derived the count recursion below, then checked its constant-rate operators against current RevBayes `ComputeLikelihoodsLtMt.cpp` (forward OBDP traversal). They agree. The diversity-dependent extension still needs numerical and measure-normalization tests.

Assume a single exchangeable clade founded by one lineage at origin time, per-lineage birth rate `lambda(N)=lambda0 exp[-alpha(N-1)]`, constant death rate mu, nonremoving specimen sampling rate psi, present-day independent sampling rho. N counts all living lineages, including ones absent from the reconstructed tree. Alpha must be nonnegative. Optional independent unplaced occurrence sampling has rate omega per lineage. Rates and diversity refer to the modeled clade, not arbitrary global biodiversity.

The model's exact infinite-state likelihood is tractable through a one-dimensional count state under these assumptions. Tree branches are not probabilistically independent. Ordinary FBD pruning formulas cannot be modified merely by replacing lambda with the observed lineage count.

## Sources read and checked

- Etienne et al., *Diversity-dependence brings molecular phylogenies closer to agreement with the fossil record*: https://pmc.ncbi.nlm.nih.gov/articles/PMC3282358/ . Web search returned full readable article text; direct open sometimes returned a bot challenge.
- User-provided supplement `/Users/taravser/GitHub/rb-fbdp/papers/rspb20111439supp.pdf`, Appendix S1 pp. 1–5, S2 combinatorial construction. Text extracted; page 4 containing backward equation S12 visually inspected. The supplement explicitly distinguishes branching-time density from a particular topology's density; its matrix products include k lambda factors for branching times, and its missing-species sampling is fixed-size uniform sampling, not independent Bernoulli rho sampling. These conventions must not be mixed.
- Current checked-out RevBayes `src/core/functions/phylogenetics/ComputeLikelihoodsLtMt.cpp`: forward traversal has the same constant-rate hidden-birth, hidden-death, terminal-fossil, sampled-ancestor and present-boundary operators derived below, plus an independent occurrence channel. Its observed split multiplier is lambda.
- Current `AbstractBirthDeathProcess.cpp` documents an oriented unlabeled to labeled nonoriented tree factor `2^(n-1)/n!`. This generic formula alone is insufficient evidence for sampled-ancestor trees; check the specialized override and match its exact convention.

## Count-state derivation

Traverse forward in chronological time, oldest to youngest. Let k(t) be the number of branches of the observed sampled tree crossing t, and h be the number of additional living lineages with no represented continuation at this crossing. Set N=k+h. Define q_h(t) as the joint density of the observed past skeleton, summed over compatible hidden histories and its embeddings. q is not an ordinary probability vector: multiple ways to embed the observed skeleton in a complete genealogy are counted.

Between observation events, with k fixed:

`dq_h/dt = (h-1+2k) lambda(k+h-1) q_(h-1) + (h+1) mu q_(h+1) - (h+k)[lambda(h+k)+mu+psi+omega] q_h`.

Terms with negative indices vanish. This is a tridiagonal positive linear system `dq/dt=A_k q`.

Derivation by event enumeration:

1. A hidden lineage births: h choices, each rate lambda(N); h increases by one.
2. A represented lineage births with exactly one daughter continuing its represented branch: k choices and two daughter choices; h increases by one. Combined permitted birth rate is `(h+2k)lambda(N)`.
3. A hidden lineage dies: h choices at mu; h decreases by one.
4. A represented lineage dies between prescribed tree events: incompatible and has no positive transition.
5. Any psi or omega sampling event between prescribed observations is incompatible and has no positive transition. Their hazards remain in the diagonal.
6. The complete-history no-event hazard is N(lambda(N)+mu+psi+omega). Do not replace the diagonal by minus the sum of the allowed off-diagonal coefficients.

All lambda evaluations use the source population size before birth. In matrix notation `A[h+1,h]=(h+2k)lambda(k+h)`, while `A[h,h-1]=(h-1+2k)lambda(k+h-1)`. An off-by-one evaluation here changes the model.

Because the hidden birth coefficient includes embedding multiplicity, A is Metzler but generally not a substochastic CTMC generator. A forward sum greater than one is not itself an error; q contains a density and embedding weights.

## Event operators

The raw kernel below uses the oriented-tree convention with a lambda multiplier at an observed bifurcation, as does RevBayes OBDP. Apply the correct RevBayes topology/label conversion separately.

- Observed bifurcation on a particular branch: `q_h^+ = lambda(k^-+h) q_h^-`, followed by `k^+=k^-+1`. Thus h is unchanged and total N increases by one. Do not multiply by k for a particular topology; that corresponds to a different branching-time measure.
- Sampled ancestor: `q_h^+=psi q_h^-`, with k and h unchanged. There is no true speciation event and no birth-rate multiplier. RevBayes's zero-length side tip is a representation of this event.
- Terminal fossil specimen: `q_(h+1)^+=psi q_h^-`, with `k^+=k^--1`. Total N is unchanged. The lineage is biologically alive immediately after its last observation, so its continuation becomes hidden and must subsequently be integrated out. Multiplying by mu or removing the lineage would be wrong under nonremoving sampling.
- Unplaced occurrence from a *separate, disjoint observation channel*: `q_h^+=omega(k+h)q_h^-`, with k,h unchanged. Its lineage is not specified, hence N possible choices. The occurrence can occur on a hidden lineage. Assigning all occurrences to known branches would silently change the likelihood.

If removal r is later supported, terminal fossil update is the sum of a removed case `psi*r*q_h^-` and a continuing case `psi*(1-r)*q_(h-1)^-`, while k decreases. SA receives `psi*(1-r)`. Start with r=0 to avoid unnecessary scope.

Separate psi and omega only if data ascertainment actually defines two disjoint streams, or if both derive from one sampling rate with a specified probabilistic classification scheme. A specimen used as a phylogenetic terminal must not also be counted as an unplaced occurrence. Species-assigned or taxonomically constrained occurrences generally need explicit placements or extra marked state; omega*N handles wholly unassigned clade occurrences only.

## Boundary and full density

Origin start: `k=1, q_0=1`, all other entries zero. The origin age must precede the oldest observed event; the observed crown bifurcation, when present, is handled as a normal split event with a lambda(N) multiplier.

At present, with n represented extant tips:

`L_kernel = sum_(h>=0) q_h(0) rho^n (1-rho)^h`.

No binomial choose factor belongs here: represented identities have already been tracked through the skeleton embedding. Check boundary cases rho=0, rho=1 and n=0 explicitly rather than relying on floating-point 0^0 conventions.

`L_tree = C(tree, label convention) * L_kernel / P(conditioning event)`.

C must include the appropriate oriented/nonoriented and label normalization. Constants that depend on tree state cannot be discarded in topology or sampled-ancestor MCMC. In particular, changing a fossil from a tip to a sampled ancestor changes the number of real bifurcations; a missing factor of two per bifurcation then changes posterior odds. The extant-only generic factor is a useful test, not a substitute for auditing the FBD/OBDP tree-shape method.

Starting at a fixed crown with two initial lineages is a different distribution: initialize `k=2,q_0=1`, omit the forced root birth density, and specify what survival/sampling condition each initial lineage must satisfy. It is not an interchangeable shortcut for an origin-conditioned process.

## Conditioning must be recalculated under diversity dependence

Recommended first implementation: explicit `condition="none"` and an origin age. If survival conditioning is included, name its event unambiguously.

1. At least one sampled extant lineage: evolve an ordinary birth–death count CTMC, `G[n+1,n]=n lambda(n)`, `G[n-1,n]=n mu`, `G[n,n]=-n(lambda(n)+mu)`, from one lineage. If p_n is its endpoint law, denominator is `1-sum_n p_n(1-rho)^n`. Fossil rates are omitted because fossil observations are not part of this particular event.
2. At least one observation of any kind: use the same count CTMC with additional killing `(psi+omega)n` and endpoint `(1-rho)^n`; one minus the resulting no-observation probability is the denominator.
3. Both crown daughter groups have sampled extant descendants: use a two-type count process (n1,n2), birth rates `ni lambda(n1+n2)`, death rates `ni mu`, and endpoint weight `[1-(1-rho)^n1][1-(1-rho)^n2]`. The two groups are coupled through N, so squaring a one-lineage survival probability is wrong.
4. At least two sampled extant descendants of an origin process is different from both crown daughter groups surviving. Do not use these events interchangeably.

A conditional event should represent actual dataset ascertainment. Existing closed-form constant-rate `pSurvival` or OBDP `GetFunctionUandP` cannot be reused by evaluating an averaged lambda. A no-conditioning kernel should be validated first.

## Numerical scheme and inference

Use a sparse matrix-exponential action, or a stable adaptive ODE method, for each constant-k interval; multiply event operators; rescale q and accumulate log scales. Parameter likelihood updates affect all intervals. Tree moves can alter event ordering and every subsequent k; a local-branch cache is not generally sufficient. Distinct fossil and branch events must be ordered deterministically; sampled-ancestor zero-length representation must become one sampling event, not a sampling event plus an artificial birth.

The infinite-state formula is exact; a finite cap is an approximation. A killed boundary retaining all original diagonal hazards, with outbound transitions omitted, gives a lower raw nonnegative path sum as the cap increases. Do not reflect births or modify lambda to force a ceiling. Check log likelihood convergence over increasing caps at initial values and posterior draws. Small q at the final boundary or at only event times is not a rigorous error certificate: the bridge can visit large hidden counts in between, and downstream weights matter. Conditioning numerator and denominator truncation errors must both be controlled. Cap rules used inside MCMC must be deterministic for each state; path-dependent adjustments can break the intended Markov target.

Forward–backward products yield `P(H(t)=h|tree,parameters) proportional to q_h(t)b_h(t)`. This gives inferred total diversity N(t)=k(t)+h with uncertainty and offers a stringent forward/backward consistency check.

Sequence and morphological likelihoods remain ordinary tree likelihoods conditional on the sampled tree, multiplied by this coupled tree prior. Hidden unsampled branches with no characters integrate out of standard Markov character models; this assumes character evolution does not itself depend on total diversity or drive birth/death. Unknown fossil placements and ages remain tree/age MCMC variables. Fossils without character scores are represented by missing character rows if treated as tree specimens; this does not itself solve taxonomic assignment or ascertainment.

## Verification gates before calling this an implemented method

- Alpha=0 agrees with existing constant-rate OBDP/FBD, including absolute density and SA topology changes, not merely parameter-dependent terms.
- Psi=omega=0 agrees with extant diversity-dependent HMM at matched root/sampling/conditioning conventions.
- Mu=0,rho=1,no fossils: all hidden paths are eliminated by the endpoint. The surviving density has interval hazard `-k lambda(k)` and observed birth factors lambda(k).
- Lambda=0, one origin lineage: serial fossil ancestors followed by extant sampling have density `exp[-(mu+psi)T] psi^m rho`; omega=0. A terminal fossil at time s followed by no later observation integrates the lineage's unobserved survival/death rather than deleting it at s.
- Independently enumerate a tiny complete-history state space with lineage identities; prune to observed skeletons and compare count marginalization, including fossil-tip/SA ratios.
- Independently simulate complete histories under the actual N-dependent hazards, sample specimens/extant tips, prune, and perform prior-predictive checks and small simulation-based calibration for tree/parameter inference.
- Compare forward and backward evaluations; test positivity, no-event cases, impossible trees, chronological ties, rho boundaries, and cap convergence.

## Augmentation comparison

Full hidden-tree MCMC evaluates an elementary complete-history density (event-rate products and integrated total hazards) and avoids a count cap, but adds a variable-dimensional hidden genealogy. Its difficult birth/death/subtree proposals, label combinatorics, mixing and pruning consistency make it a poor first production path. It is valuable as an independent small-case validator and becomes necessary if lineage-specific ecological competition or occurrence identities destroy exchangeability. Pseudo-marginal particle methods need unbiased likelihood estimates and variance control; they are not obtained by simply plugging a stochastic count estimate into ordinary MCMC.

Conclusion: implement the deterministic count-marginalized likelihood first, adapting OBDP's event architecture while replacing rate evaluation, conditioning, and numerical verification. The central math is tractable; normalization and data-observation semantics are the main design decisions that must be explicit before production code.

# Addendum: core requirement is repeated occurrences assigned to named species

The anonymous specimen derivation above is a restricted validation submodel, not a solution to the user's named-species requirement. After clarification I read the 2018 FBD range paper in [arXiv HTML](https://arxiv.org/html/1706.10106), especially its speciation-mode definition and oriented sampled versus extended sampled trees. Its distinction is essential: several observations sharing a species name constrain biological species persistence, not merely a path's ancestry. Budding retains the parent species and creates one novel species; the two daughter roles are biologically different. The paper explicitly warns that its labeled-tree conversion is available for extended sampled trees but not a general sampled-tree conversion. Therefore neither an automatic factor of two nor treating a range endpoint as an ordinary fossil tip is defensible without a specified latent representation.

## Defensible exact target before an optimized algorithm

Let H be a complete budding genealogy with a unique species identity for every species lifetime, all hidden birth and death events, parent identities at births, occurrence-to-species assignments, and fossil/extant observations. A birth preserves its parent identity and introduces a new identity. For a fully specified identity genealogy and exact event times, define

`p(H,occurrences,extant-sampling|theta) = C(H) * [product over births lambda(N_before)] * mu^(number deaths) * psi^(number observed occurrences) * rho^(number sampled extant species) * (1-rho)^(number unsampled extant species) * exp{-integral N(t)[lambda(N(t))+mu+psi] dt}`.

C encodes only a declared labeling/reference-measure convention. There is no N factor for a birth or death whose lineage identity is specified; N appears after summing event identity possibilities. Identity-compatible observation constraints are indicators. Additional Poisson streams or uncertain ages require their observation factors, not duplicate counting. If occurrence records are a subsample of actual preserved fossils, psi denotes the actual observation process or an explicit detection model is required.

The target marginal likelihood is the integral/sum of this density over all hidden histories H compatible with the observed named occurrences and the inferred character tree, multiplied by the character likelihood conditional on that history's sampled subtree. This path-integral definition is a mathematically sound baseline and guides an exact simulator and small-history enumerator. Full-history reversible-jump MCMC is one possible computational realization, although mixing and correct proposals must be demonstrated. A total-evidence target can retain the relevant species-orientation/lifetime variables while integrating only exchangeable hidden histories.

## Candidate optimized route: count plus species-identity state

Use a sampled skeleton with all occurrence locations, latent budding orientations, and an identity-constraint state z. At each crossing z must determine which represented branches continue which constrained species identity and which identity changes have occurred in gaps between differently named observations. It must also encode any sampled-species lifetime or continuation variables retained beyond the last occurrence. Additional living lineages with no future represented sample are exchangeable and may still be counted by h. Thus the candidate messages are `q_(h,z)`, not q_h alone.

Conditional on a correctly specified z and k represented branches, births are enumerated as:

- A hidden lineage buds: h choices, h increases by one, z unchanged.
- A represented parent buds an unrepresented novel species, retaining the represented parent: one allowed continuation per branch, h increases by one, z stays on that identity.
- A represented branch follows a novel budding daughter while the parent becomes unrepresented: at most one additional continuation per branch, h increases by one, and z changes to the new species identity. This is forbidden across a segment anchored to one known species by observations at both ends.

If a of k branches are identity-anchored, and all other replacement choices are unconstrained and equivalent, the aggregate permitted hidden-birth coefficient is `(h+2k-a)lambda(N)`. It is `(h+k)lambda(N)` when every represented branch is identity-persistent. But `(h+2k-a)` is only a local check, not the entire species-aware algorithm: different allowed replacement events can lead to different z states and must not be collapsed prematurely.

The diagonal remains `-N[lambda(N)+mu+psi]`. At an observed split, sum over or sample only admissible parent/daughter role assignments, each with a lambda(N) factor. At an assigned occurrence, the named species must occupy that represented location; use psi once, with no N factor. A final observation is not necessarily the last represented continuation of its species: its species may persist and bud another sampled species later. Last-observation/lifetime/branch-continuation handling depends on the augmented tree representation.

A useful finite-automaton interpretation: along a represented path, same-name observations prohibit identity replacement between them; different-name observations demand at least one species-creating replacement unless an intervening observed budding node supplies it. Simply reducing 2 to 1 inside observed ranges does not enforce the latter condition. Known occurrences from one species cannot occupy two contemporaneous daughter paths. Names are never reused after species death or replacement. These requirements should be encoded and proven before implementing q_(h,z).

Minimality of z is not established here. Its size may grow exponentially with simultaneously unresolved assignments. A rigorous milestone is to derive a lumpability proof for a chosen oriented/extended skeleton; absent that proof, full-history augmentation remains the reference target. First validate the marked-count method for one species with repeated observations, two species separated by an unobserved budding event, and a species whose observed range spans a budding split. Then compare alpha=0 against the range FBD formulas at exactly matched orientation/lifetime conventions.

## Practical representation recommendation

Plan the production model around species data, with occurrence rows holding stable species IDs and age uncertainty, character observations linked to appropriate species/time positions, and latent species birth/death and budding orientation information. Keep ordinary specimen trees as a deliberate special case. Do not publish an anonymous q_h likelihood under a name implying it handles assigned occurrences.

A staged implementation can still use the anonymous scalar recursion as a verified numerical engine, but the essential scientific milestone is species-aware marginalization or augmentation and MCMC proposals that preserve named-species continuity. Independent simulation must produce species identities at birth and retain them across parent continuations, so that the central constraint is validated rather than merely checked on anonymous trees.
