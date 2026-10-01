# Exact hidden-count likelihood for an oriented extended budding tree

## Decision

**The scalar hidden count is sufficient for the specifically defined oriented extended-tree representation below.** The previous need for a general marked state applied to ordinary sampled trees without retained species endpoints. It does not prohibit this narrower, useful representation. This is an analytical derivation; numerical and integration tests remain necessary before calling its implementation validated.

Primary-source check: [Stadler et al. (2018), Theorem 2 and its proof](https://arxiv.org/html/1706.10106). Its unrestricted q factor applies before first observation; its identity-preserving q-tilde factor applies from first observation to the retained death/present endpoint. The proof's difference of one versus two hidden birth routes is exactly the distinction used below. Its displayed likelihood is conditioned on obtaining at least one sample. The derivation here initially omits that conditioning.

## 1. Specify the represented object exactly

Let E be a rooted oriented binary extended tree with:

- A stem from the process origin T.
- Exactly one terminal endpoint for every observed named species i. A terminal at age d_i>0 is that species' actual extinction; d_i=0 means it survives to the present.
- At every retained split, an A child (continuation of parent identity) and D child (new daughter identity). These roles are part of the latent state.
- Every occurrence assigned to its named species and placed on its path. All occurrences of species i lie on the A-continuation path ending at endpoint i.
- The oldest observation age o_i, including the present sample if there are no fossils. The continuous path from o_i to d_i is protected: it must retain identity i, including across retained A continuations.
- The structural attachment age b_i where the straight A-continuation path leading to endpoint i starts (at a retained D split or at the process origin).

Require d_i <= all occurrence ages <= o_i <= b_i, and validate the topological A-path assignments, not merely these scalar inequalities. Extant samples have age zero and their species endpoint is zero. Distinct observed species have distinct terminal endpoints. An endpoint without either a fossil occurrence or a sampled present observation is not an observed species endpoint in this representation; it belongs in hidden history instead.

**b_i is not necessarily the biological origination time of species i.** Hidden daughter-following births before o_i can replace the identity on this path. Thus the filter integrates the true species origin in that segment. If the interface exposes b_i as an actual species birth age or requires an actual birth age, this derivation is not that distribution; further augmentation or protection would be required.

An ancestral sampled species has its own retained endpoint after its last observation, even if it produces observed descendant species. Consequently it cannot disappear into a serial sampled-ancestor path as it can in an ordinary sampled tree.

## 2. Deterministic schedule

Traverse forward with elapsed time u=0 at the origin and u=T at present. In an event-free interval let:

- k = number of retained E branches crossing the interval;
- s = number of these branches covered by protected species intervals;
- h = additional living species not on retained branches;
- N = k+h.

For a valid E, s can equivalently be counted as the number of species whose extended intervals (d_i,o_i) cross the time slice. Necessarily 0 <= s <= k. Which retained branch carries which protected identity is known from E; it must be retained for input validation and observation-tree construction, but it need not index the likelihood message.

Rates are homogeneous over species except for their common dependence on N: per-species lambda_N=lambda0 exp[-alpha(N-1)], death mu, nonremoving occurrence sampling psi. No character-dependent diversification, anagenesis, age-dependent hazards, species-specific recovery, or external competing clades.

## 3. Weighted lumping proof

Consider compatible complete histories with an embedding of the fixed retained skeleton E. At a time slice retain the actual labelled identities temporarily, then classify each additional lineage as hidden. Every hidden lineage must have no future retained observation or endpoint; its future genealogy is entirely marginalized. Protected retained lineages have specified identity; unprotected retained lineages have an exchangeable unknown identity until their first named observation.

For any two compatible detailed configurations with the same h at the same scheduled slice:

1. Their total population N is equal, so every elementary species event has the same rate law.
2. A hidden species buds: h interchangeable parent choices, each rate lambda_N. Both species remain hidden, so h increases by one.
3. A retained species buds an unretained daughter and the parent continues the retained branch: k choices, each rate lambda_N. This is compatible on protected and unprotected branches and increases h by one.
4. A retained species buds a retained daughter while the parent becomes hidden: k-s choices, each rate lambda_N. It is incompatible on a protected branch because it changes the identity. On an unprotected branch it produces an admissible unknown identity and increases h by one.
5. A hidden species dies: h interchangeable choices at mu, decreasing h by one.
6. A retained species death is incompatible except at its scheduled true extinction endpoint.
7. Any fossil occurrence outside the scheduled observations is incompatible, because the dataset is the complete set of observed occurrences from the modeled observation process.

All compatible detailed configurations therefore have identical aggregate weighted transition coefficients to h+1 and h-1, and the same total no-event hazard. Their future compatibility depends only on the fixed schedule plus h: there is no unresolved named identity among hidden species, and an unknown identity on an unprotected branch carries no future requirement to equal an earlier observed identity. Its later named endpoint is first anchored at o_i. Distinct endpoint paths have separated through a retained D birth, so their identities cannot accidentally be the same species.

This yields an exact aggregation of **embedded-history weights**, not the lumping of an ordinary conservative CTMC. The two retained-continuation possibilities count alternative skeleton embeddings. It is therefore expected that the resulting matrix can have positive column sums.

The argument applies to every finite compatible sequence of elementary events. Summing those sequences and integrating their times produces the semigroup of the tridiagonal system below. With alpha>=0 and finite rates, the birth process is dominated by a linear Yule process, so there is no finite-time explosion. The infinite nonnegative path sum defines the exact likelihood; a finite hidden cap remains a numerical approximation.

## 4. Between-event equation

Let v_h be the forward unnormalized embedded-history density. Then

`dv_h/du = (h-1+2k-s) lambda_(k+h-1) v_(h-1) + (h+1) mu v_(h+1) - (k+h)[lambda_(k+h)+mu+psi] v_h`.

Terms outside h>=0 vanish. In source-column notation:

- `A[h+1,h]=(h+2k-s)lambda_(k+h)`;
- `A[h-1,h]=h mu`;
- `A[h,h]=-(k+h)(lambda_(k+h)+mu+psi)`.

The birth rate is evaluated at the source count before birth. Do not replace the diagonal by minus the sum of retained transitions. At k=s=0 this is the ordinary count process killed by unscheduled fossil sampling.

## 5. Scheduled event maps

Initialize v_0=1, v_h=0 otherwise, k=1, s=0, except for deliberately specified observations at the origin (normally a zero-probability coincidence).

- Retained oriented split: multiply v_h by lambda_(k+h), then k increases by one. A protected parent remains protected on A; D is unprotected until its own oldest observation. Thus s does not change at a generic split.
- Occurrence: multiply v_h by psi. At that species' oldest observation, s additionally increases by one for the following younger interval. Other occurrences change neither k nor s.
- True species extinction endpoint at positive age: multiply v_h by mu; k and s each decrease by one. h is unchanged. The endpoint kills a retained living species; it is not a fossil-tip-to-hidden transition.
- Present: let l retained species be sampled at present, and r retained species survive but are not sampled at present. Then k=l+r and

  `L_raw(E,F)=sum_h v_h rho^l (1-rho)^(r+h)`.

A fossil species surviving unsampled to present receives (1-rho), not rho. A species observed only at present has a zero-length protected interval and receives rho. At rho=1, fossil species known to survive without present sampling are impossible.

No separate sampling multiplier is attached to a last fossil beyond its ordinary psi. The protected interval extends past that fossil to the true endpoint.

## 6. Why different-species gap flags are unnecessary here

In a nonextended sampled tree, an earlier named species and a later named species can lie serially on one retained straight path after their unsampled continuations have been deleted. Their difference requires at least one hidden identity-changing birth in the intervening gap.

In E, the earlier named species' actual endpoint remains in the tree. To connect the later different species there must therefore be a retained split separating its D ancestry from the earlier species' protected continuation. The requisite species-separating event is explicit in E. Hidden identity switches before each path's first occurrence are unrestricted, while switches after first occurrence are forbidden. There is no serial gap obligation left to encode in a filter flag.

A purported E with two distinct named species on the same uninterrupted A path ending in one endpoint violates the input representation. Do not silently feed it to the scalar filter.

## 7. Constant-rate theorem reduction

For alpha=0 define Gamma=lambda+mu+psi. Let p(t) solve `p'=mu-Gamma p+lambda p^2`, p(0)=1-rho, in age before present. Define normalized propagators q(0)=qtilde(0)=1 by

`d log q/dt = -Gamma+2lambda p(t)`;

`d log qtilde/dt = -Gamma+lambda p(t)`.

Thus `qtilde(t)=exp(-Gamma*t/2)*sqrt(q(t))`.

Because constant-rate hidden descendant lineages evolve independently, integrating them gives q factors on each unprotected portion and qtilde on each protected portion. For n retained species endpoints, m positive-age extinction endpoints, l sampled present endpoints, r unsampled present endpoints, and F fossil records, the raw oriented density must equal

`psi^F mu^m rho^l (1-rho)^r lambda^(n-1) product_i [q(b_i)/q(o_i)] [qtilde(o_i)/qtilde(d_i)]`.

This is the unconditioned version of the paper's Theorem 2. To compare its reported conditioned density, divide this raw result by `1-p(T)`, the probability of any fossil or present sample. Do not use extant-survival conditioning for that comparison.

This supplies a stringent independent reference for arbitrary small extended trees, including first observations on opposite sides of retained splits, fossils followed by extinction, and retained extant species without present samples. Agreement with only a founder-protection fixture is insufficient.

## 8. Measure, orientations, and conditioning

The raw formula is for an oriented embedded tree. If labels are uniformly assigned to n observed species, retaining orientations gives a 1/n! conversion under that convention. If labels are treated as a fixed conditioning convention, this fixed-n constant may instead be excluded consistently. Do not add a factor of two at every retained split while also sampling its A/D orientation.

For a labelled nonoriented tree, sum over compatible orientations. In this exchangeable model, any orientations producing the same k(t), s(t), observation and endpoint schedules have the same raw likelihood. A factor equal to their number is therefore permissible after verifying compatibility; if the data representation changes those schedules, perform the sum rather than assume it. The 2018 constant-rate corollary uses 2^(n-v-1)/n!, with v forced orientations. Use it as a matched reduction test and document the correspondence between its v and the actual representation.

Under diversity dependence, the probability of obtaining any sample is calculated from the ordinary total-count birth-death process killed at N psi, with terminal no-present-sample weight (1-rho)^N. One minus that no-sample probability is the conditioning denominator. This process does not use the embedding coefficients 2k-s. Start with an unconditioned origin model to isolate likelihood verification.

## 9. Counterexamples / scope boundaries

No counterexample was found under the precise E definition and exchangeable budding assumptions. The scalar closure fails, or represents a different model, if:

- Species endpoints are discarded but serial differently named sampled ancestors remain: needs marked gap state or explicit hidden identity events.
- More than one distinct species is assigned to one protected A path, or a species is assigned to two contemporaneous paths: incompatible E, likelihood zero.
- Named identity persists only until the last fossil rather than its true endpoint: permitting later D replacements would misrepresent the endpoint as another species' extinction.
- True species origination times are separately constrained before o_i: unrestricted replacement before o_i is no longer acceptable.
- A hidden lineage carries a future named observation/end-of-lifetime obligation: it was misclassified as hidden and must be retained/marked.
- Species have different birth, death, fossil recovery, or age-dependent hazards: h alone generally loses necessary state.
- Hidden speciation changes character states or branch-clock regimes not captured by the observation tree: the tree prior may still aggregate, but the joint character likelihood may require more latent state.

## 10. Production-gate tests

1. Arbitrary small constant-rate oriented extended trees: compare the full scalar filter with the q/qtilde product, including absolute raw density and optional matched any-sample conditioning.
2. n=1 fossil species with latent positive death d and two fossils: require mu endpoint, protection o→d, and continued h-only evolution after d. Vary d while fossil times stay fixed.
3. n=1 fossil species survives unobserved to present: use 1-rho. Contrast the same species sampled at present using rho; for fixed history their ratio must be rho/(1-rho).
4. n=2 parent observed before the split, child after it: parent remains protected through A; s unchanged at split and increases at child's first observation.
5. n=2 oldest observations after the split: both daughter paths initially unprotected; changing each first-occurrence age changes only the associated protection boundary.
6. n=3 with two budding events within one species' protected lifetime: its identity is the same through both A continuations; no duplicate species count.
7. Same data represented as an ordinary sampled tree with an ancestral sampled species: explicitly retain its endpoint and the species-separating split in E; demonstrate why dropping them requires an extra gap correction.
8. lambda=0 single species: density is `psi^F exp[-(mu+psi)(T-d)] mu` for an extinct endpoint, or `psi^F exp[-(mu+psi)T] rho` / `(1-rho)` for sampled/unsampled present endpoints.
9. mu=0,rho=1 guarantees no hidden lineages at endpoint; recover direct full-observation pure-birth event density at DD rates.
10. Tiny complete identity histories: enumerate all single and double hidden-event routes and their embeddings, verify upward multiplicities h+2k-s and death multiplicity h, including switches immediately before versus after o_i.
11. Paired orientations: where unconstrained both should have equal likelihood; orientation forcing a protected species down D must be rejected. Validate any summed-orientation factor against explicit enumeration.
12. DD forward/adjoint equality, step/tolerance convergence, hidden-cap convergence, alpha→0 limit, and input rejection for impossible identity placements.

The analytical sufficiency gate is satisfied for this narrowly defined oriented extended-tree model. The independent numerical constant-rate general-tree tests and full identity-history event enumeration should precede production claims.
