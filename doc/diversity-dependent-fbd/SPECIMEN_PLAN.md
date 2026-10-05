# Specimen-based diversity-dependent FBD tree prior: implementation plan

Design begun 1 October 2026; implementation updated 4 October 2026. This is the replacement design requested after the named-species extended-tree prototype. Preserve the earlier prototype and its results as a separate model; its recovery experiments do not validate this specimen model.

## Updated implementation request — 4 October 2026

The user now requests a separate new tree-prior function, rather than adding an argument to `dnBDSTP`, and morphology-only integration testing. Use the name `dnDiversityDependentSpecimenFBD` (short alias `dnDDFBDP`) to distinguish it from the earlier named-species `dnDDFBD`. Retain the BDSTP sampled-ancestor representation and standard collapse/expand move. The public API recommendation below is superseded by this separate-function requirement.

The user selected the original exponential decline: lambda(N) = lambda0 exp[-alpha (N-1)], alpha >= 0. Use constant mu and psi, no skyline, and no molecular partition in the requested integration test.

## 1. Target model and user interface

Use the ordinary sampled-ancestor `TimeTree` used by `dnBDSTP(r=0, ...)`, directly as the tree argument to `dnPhyloCTMC`. Each data label identifies one observation: a fossil specimen with an age, or a sampled living taxon. A fossil can be terminal or ancestral to later observations; MCMC infers that status and placement. There are no species ranges, named-species persistence constraints, explicit true extinction endpoints, protected-lineage states, or ancestral/daughter orientation in this model.

All biological splits use the same speciation process. Unobserved splits and extinct/unsampled continuations are integrated out. Sampling does not kill a lineage (`r=0`). Constant per-lineage extinction and preservation/recovery rates are mu and psi; present sampling is independent Bernoulli rho. Speciation per lineage is

    lambda(N) = lambda0 exp[-alpha (N-1)],   N >= 1, alpha >= 0.

N is the total number alive, including hidden lineages. The total biological birth hazard is N lambda(N), not lambda(N). Lambda0 is the rate at N=1. Alpha=0 reproduces constant-rate `dnBDSTP` with time/sampling conditioning. The implementation audit below records a native survival-normalization difference at rho<1.

Interpret “single occurrence” as one occurrence per input label, following the requested BDSTP specimen formulation. Different observations may be linked as sampled ancestors on a continuing lineage; they are not preassigned to a persistent named species. We do not impose a prohibition on a lineage ever being sampled twice: that would change the ordinary Poisson sampling process and would not be the requested BDSTP(r=0) model.

Implemented API: a separate `dnDiversityDependentSpecimenFBD`, alias `dnDDFBDP`. The existing `dnBDSTP` is unchanged. Runnable morphology-only example: `examples/DDFBDP/morphology.Rev`.

```rev
alpha ~ dnExponential(10)
fbd_tree ~ dnDDFBDP(originAge=origin_time, lambda0=lambda0,
                    alpha=alpha, mu=mu, psi=psi, rho=rho,
                    taxa=taxa, condition="sampling")
phyMorpho ~ dnPhyloCTMC(tree=fbd_tree, Q=Q_morpho,
                       branchRates=clock_morpho, type="Standard", coding="variable")
phyMorpho.clamp(morpho)
moves.append(mvCollapseExpandFossilBranch(fbd_tree, origin_time, weight=6))
```

The character partitions are optional and retain their normal substitution models and ascertainment rules. The graph has one stochastic sampled tree shared by both CTMCs. There is no observation-tree projection and no extra DD likelihood added elsewhere in the graph:

    posterior ∝ L_morph × L_DNA × L_fossil_age × p_DD(sampled_tree | parameters) × parameter_priors.

First release: origin initialization, scalar rates, r fixed at zero, and conditions `time`, `sampling`, `survival`. Reject skyline vectors/timelines, burst events, diversified sampling, removal, and root-age initialization in the DD path with explicit messages. The origin is estimated if assigned a prior; “condition on origin” does not mean it must be fixed. Root-age conditioning is deferred because survival of two initial subtrees is coupled by total diversity; squaring a single-lineage probability is invalid for alpha>0.

## 2. Tutorial and source reading

Read the [specimen tutorial](https://revbayes.github.io/tutorials/fbd/fbd_specimen), all six linked Rev scripts, and the three data files. Exact downloaded copies and SHA256 hashes are in [specimen-tutorial-reference](specimen-tutorial-reference/sources.json).

The tutorial supplies a sampled tree to molecular and morphological CTMCs. Its `r=0` tree prior supports terminal and ancestral fossils, clade constraints, and specimen-age uncertainty. Its collapse/expand move changes fossil status while ordinary topology and age moves explore the rest of the tree.

File-level findings:

| File | Role / consequence for this implementation |
|---|---|
| `mcmc_CEFBDP_Specimens.Rev` | Loads both matrices, adds missing taxa, connects the modules and runs MCMC. Preserve that graph. |
| `model_FBDP.Rev` | Constructs `dnBDSTP`, wraps a topology constraint, adds fossil collapse/expand, FNPR, node/root/fossil-age moves and ancestor-count monitoring. Replace the prior's calculation, not this tree representation. |
| `model_GTRG.Rev` | GTR plus gamma molecular likelihood on `fbd_tree`; no DD-specific change. |
| `model_Morph.Rev` | Binary Mk plus gamma, a strict morphological clock and variable-character correction. Preserve the actual script's gamma component. |
| `model_UExp.Rev` | Per-branch molecular clocks indexed on the existing binary representation; zero-length sampled-ancestor branches need no extra evolutionary time. |
| `summarize_CEFBD.Rev` | Summarizes a saved tree trace; retain the full sampled tree as well as any pruned summaries when validating ancestors. |
| `bears_taxa.tsv` | 22 observation labels, eight with min age zero and 14 fossils; use actual age bounds and the existing Taxon contract. |
| `bears_cytb.nex` | Ten rows, 1,000 molecular sites; remaining taxa are missing data. |
| `bears_morphology.nex` | Eighteen rows, 62 binary morphological characters; remaining taxa are missing data. |

Some script comments call the age intervals “stratigraphic ranges”. For this design they are uncertainty bounds on a single occurrence, not durations occupied by a species. Preserve the downloaded files verbatim; document any adaptation rather than silently repairing the references. For example, the downloaded Ballusia bounds differ from the page's table, so tutorial regression must specify which input version it uses.

Relevant checked source:

- `src/revlanguage/distributions/phylogenetics/tree/Dist_BDSTP.{h,cpp}`: Rev parameters, `psi`/`phi` aliases, origin/root handling, optional initial tree; constructs an `AbstractBirthDeathProcess*`.
- `src/core/distributions/phylogenetics/tree/birthdeath/BirthDeathSamplingTreatmentProcess.{h,cpp}`: current constant/skyline likelihood, event classification, sampled-ancestor shape factor, initialization and validation simulation.
- `src/core/distributions/phylogenetics/tree/AbstractRootedTreeDistribution.cpp`: zero-length ancestor and age checks; adds the tree-shape term once.
- `src/core/distributions/phylogenetics/tree/birthdeath/AbstractBirthDeathProcess.cpp`: default survival handling assumes a root start and independent subtrees. Override the divergence-time calculation for the DD origin model rather than inheriting its automatic squared survival factor.
- `src/core/moves/proposal/tree/CollapseExpandFossilBranchProposal.cpp`: checks `allowsSA()`, changes the parent age and fossil ancestor flag, and supplies its own proposal-density terms. No species identities or DD rates are embedded in this proposal.
- `TopologyConstrainedTreeDistribution.h`: delegates `allowsSA()` to its base distribution, allowing the tutorial wrapper to remain usable.

## 3. Proposed hidden-count likelihood

The following defines the numerical kernel; validation results are recorded separately in `examples/DDFBDP/README.md`. Work forward from the origin to the present, processing the sampled tree chronologically.

Between observed events, let k be the number of branches of the observed skeleton currently continuing toward their next observation/split, h the number of additional living lineages with no future observations, and N=k+h. Hidden lineages include continuations after terminal fossils. Given exchangeable rates depending only on N and no species-identity constraints, h is sufficient; the observed topology supplies which retained branch each event belongs to.

Let v_h be a forward density/embedding weight. Initially k=1 and v_0=1. With column-vector convention, dv/dt=A_k v, where

    A_k[h,h]     = -(k+h) [lambda(k+h) + mu + psi]
    A_k[h+1,h]   = (h+2k) lambda(k+h)
    A_k[h-1,h]   = h mu.

Set the zero-population birth hazard to zero. The lower transition is omitted at h=0. The positive operator is not an ordinary probability-conserving generator: the 2k term sums the two ways one child of an unobserved split can carry the retained future. It must not be normalized as a stochastic matrix.

Interpretation:

- Hidden birth contributes h lambda(N).
- A birth on an observed branch with one unobserved child contributes 2k lambda(N).
- A hidden extinction contributes h mu.
- Extinction of a retained branch before its required next event, or any additional unrecorded sampling event, is incompatible; their hazards remain in the negative diagonal.
- A biological split leaving two observed descendants is handled at its explicit observed time, not by the between-event transition.

Observed-event operators:

| Event | Update | Biological interpretation |
|---|---|---|
| Genuine bifurcation | v_h <- lambda(k+h) v_h, then k <- k+1 | One new lineage; h unchanged. There is no k multiplier because the observed topology specifies the branch. |
| Sampled ancestor | v_h <- psi v_h; k,h unchanged | Sampling along a continuing observed lineage; no birth. |
| Terminal fossil | new_v[h+1] += psi v_h, then k <- k-1 | The sampled lineage becomes hidden after its last observation. N is unchanged at sampling. |

For a terminal fossil, do **not** multiply by mu or terminate the biological lineage. Sampling is not extinction when r=0. Also do not multiply by an independent scalar probability of no future samples: the shifted hidden state already integrates that continuation jointly with all other lineages and the shared N-dependent rates.

At the present require k=n_extant, then

    L_raw = rho^n_extant × sum_h v_h (1-rho)^h.

This integrates hidden living lineages that were not sampled at present. If rho=1 only h=0 contributes. After the last observation in an all-fossil tree, propagate with k=0 through to the present; do not stop at the youngest fossil.

Sampling factors occur exactly once through these event operators. There is no separate repeated-occurrence psi factor or species-range likelihood.

## 4. Sampled-ancestor representation and density measure

RevBayes represents an ancestral sample as a fossil leaf on a zero-length edge attached to a parent at the same age, with a continuing sibling. That binary container node is a sampling event, not a biological bifurcation.

Build a fresh event schedule from the current tree using the sampled-ancestor flags and parent-child relations:

1. A marked sampled-ancestor tip contributes one ancestral sampling event.
2. Its artificial parent contributes no speciation event.
3. Every other valid internal split contributes one birth event.
4. Every positive-age unmarked terminal observation contributes a terminal-fossil event.
5. Present tips supply the final rho factor.

Reject two sampled-ancestor children of the same parent, ancestor flags on extant tips, malformed zero edges, negative branch lengths, fossil ages outside their support, and nodes older than the origin. Allow an artificial sampled-ancestor root for an origin-conditioned tree. Equal fossil ages on different branches are valid; do not globally reject tied observations. Use topology-aware tie handling and test the legal exact-age cases.

Match BDSTP's labelled, non-oriented tree measure. The existing source uses

    log C = (n - a - 1) log 2 - log(n!),

where n is the number of observation labels including ancestral fossils, and a is the number of sampled ancestors. This equals one factor two per genuine bifurcation. Apply it once through `lnProbTreeShape()`. The collapse/expand move changes a, so this term must be recalculated and cannot be discarded as a constant. Do not carry over the old oriented-species model's label-only factor or its ordered topology move.

## 5. Conditioning

Use BDSTP's option names and meanings:

- `time`: origin specified, no extra sampling/survival division.
- `sampling`: divide by P(at least one fossil or extant sample).
- `survival`: divide by P(at least one sampled extant descendant), irrespective of fossils, for r=0.

Compute those denominators using the unobserved **total-count** birth-death process, started at N=1. Its birth and death rates are N lambda(N) and N mu. Any-sampling conditioning additionally absorbs at N psi; extant sampling is evaluated at the present with weight 1-(1-rho)^N. Existing direct-observation absorption machinery is reusable after API review; compute small probabilities stably in log space.

Do not accidentally apply the base class's root-conditioned factor as well. Conditional denominators need their own count-cutoff convergence checks. `condition="survival"` is not simply “something was ever fossilized”.

## 6. C++ and Rev implementation structure

Use the separate `dnDDFBDP` public entry point and a compact new core class derived from `AbstractBirthDeathProcess`. Do not extend the existing skyline's independent-lineage `E(t)` formulas by substituting a tree-based lambda: total diversity couples the lineages and invalidates that factorization.

Proposed files/work:

1. `DiversityDependentFbdLikelihood.{h,cpp}`: separate specimen entry point for the standalone event/count calculation with `Birth`, `AncestorSample`, `TerminalSample`, initialization, endpoint weights and conditioning. Reuse the existing tested positive uniformization primitive, keeping the old named-species event logic out of the new entry point.
2. `DiversityDependentSpecimenBirthDeathProcess.{h,cpp}`: scalar DAG parameters; sampled-tree event extraction; `allowsSA()=true`; cloning, swapping and dirty-state handling; origin-specific probability method and BDSTP shape factor.
3. Add `Dist_DDFBDP.{h,cpp}`, registered separately with numerical controls (`maxHiddenLineages`, `numericalTolerance`). Expose only scalar constant rates and an origin start; leave the existing `Dist_BDSTP` unchanged.
4. Reuse `TreeUtilities::startingTreeInitializer` and the BDSTP starting-tree approach for a valid MCMC initial state; support optional `initialTree`. A valid initializer is not a random prior draw. Implement a separate genuine full DD-FBD simulator for validation, or explicitly reject unsupported validation redraws until it exists.
5. Add help and a specimen example that follows the downloaded tutorial, changing the tree-prior call and adding alpha's prior/move. Preserve ordinary character models and optional clade constraints.

No custom CTMC, new fossil-branch proposal, observation-tree function, or latent extinction endpoints should be required.

Numerics: fixed hidden cap during a chain; kill overflow at the cap, including terminal-fossil shifts; log-scale intermediate vectors; use positive uniformization with nonstochastic norm bounds. Compare increasing caps and tighter tolerances outside the chain. Event intervals arise from observed times, not a skyline model. Underflow with positive density is a numerical failure, not automatic zero support. Start with recomputing the event schedule on every evaluation; optimize only after move/rollback tests establish correctness.

## 7. MCMC and fossil-age uncertainty

Reuse `mvFNPR`, `mvNodeTimeSlideUniform`, `mvRootTimeSlideUniform`, `mvFossilTimeSlideUniform`, and especially `mvCollapseExpandFossilBranch`. Sample lambda0, alpha, mu, psi and origin with ordinary scalar moves. The existing collapse/expand Hastings terms remain in the proposal; the new prior supplies the appropriate changed-state density.

Keep fossil ages on the same sampled tree supplied to the CTMCs. For an uncertain single occurrence [a_i,b_i], retain the tutorial's clamped uniform age-observation construction using `tmrca(fbd_tree, clade(fossil_i))`. Do not interpret a_i/b_i as first/last occurrences. For exact ages, fix the observation age through the taxon/tree support and omit a zero-width uniform distribution. Do not put the same age information into an additional independent density twice.

Always monitor the full tree and number of sampled ancestors for validation. CTMC likelihood and tree-prior caches must both update after fossil-age, topology and collapse/expand proposals, and both must restore after rejection.

## 8. Implementation and validation order

1. **Lock the specimen contract and likelihood.** Hand-check a fossil ancestor plus extant descendant, a terminal fossil plus an extant lineage, multiple ancestor samples, and an all-fossil tree. Prove/test terminal-sample demotion and the shape factor before MCMC.
2. **Standalone kernel.** Implement the interval and three event operators. Independently integrate small finite matrices and use pure-death/no-birth boundary formulas. Check rho=0/1, psi=0, alpha=0, zero-probability cases, tied observations and cutoff/tolerance convergence.
3. **Absolute constant-rate benchmark against native BDSTP.** Evaluate exactly the same tree and rates through both distributions, with alpha=0 explicitly selecting DD. Include terminal fossils, sampled ancestors, an ancestral artificial root, all-fossil trees, incomplete rho and all three origin-conditioning options. Require agreement around 1e-8 after cutoff convergence. Compare collapse/expand log-density differences, not just densities within one fixed ancestor configuration.
4. **Rev graph integration.** Exercise the standard fossil collapse/expand and age moves, with and without `dnConstrainedTopology`, under `debugMCMC=1`. Test rejected proposals, DAG clone/swap, checkpoint round trips, and changing ancestor counts. Couple to morphology alone for the requested implementation test. Other data partitions and large simulation recovery remain subsequent validation work.
5. **Extant-only regression.** Set psi=0 and remove fossil observations. Match the existing DD extant benchmark under the same origin, rho and density convention. DDD's power-law model is not exponential decay in N; compare rate laws exactly rather than calling different functions the same model.
6. **Independent simulation and pruning.** Simulate the full count-dependent birth-death process with non-destructive serial sampling, then prune to observations. A fossil with later sampled descendants becomes an ancestor; otherwise it is terminal even if its underlying lineage lived longer. Remove the previous named-species endpoint pruning. Preserve true hidden counts for debugging only; inference sees the sampled tree/data.
7. **Recovery.** First a small known-tree rate recovery test; then joint inference of topology, ancestor status, dates and rates. Reproduce the tutorial's graph on the frozen bear files for a smoke test. Only then launch a 100-observation simulated experiment and replicated recovery/calibration. Specify simulation selection/conditioning in advance; do not search for exactly 100 observations and then treat that as an unconditional SBC experiment. Assess tree/ancestor mixing as well as scalar diagnostics.

Acceptance criterion: the simpler model is one BDSTP-compatible sampled-tree prior for ordinary CTMC likelihoods, reduces to native BDSTP at alpha=0, changes ancestral status using the existing move, and accounts for total rather than reconstructed diversity. Earlier extended-species recovery results are not acceptance evidence for this changed observation model.

## Implementation audit: native normalization (4 October 2026)

Native BDSTP uses `RbMath::lnFactorial(n)` (a truncated Stirling formula); the new
wrapper uses the same fixed label constant for direct density comparisons.
The standalone kernel and analytic tests exclude this shape factor.

A substantive conditioning convention differs: native
`BirthDeathSamplingTreatmentProcess::pSurvival(origin,0)` computes
P(at least one sampled extant descendant)/rho. Our `survival` denominator is
P(at least one sampled extant descendant) itself. Therefore at alpha=0,
lnP_native - lnP_new = ln(rho). This known factor is explicitly removed in the
native regression for incomplete extant sampling. We retain the normalized
sampled-extant conditioning in the new model; it matters if rho is estimated.
The morphology-only example uses `sampling` and rho=1, where both agree directly.
