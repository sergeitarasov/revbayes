## name
dnDiversityDependentSpecimenFBD
## title
Diversity-dependent specimen fossilized birth-death tree prior
## description
A sampled-ancestor TimeTree prior with per-lineage speciation rate
`lambda0 * exp(-alpha * (N - 1))`, alpha >= 0, and constant extinction mu and
fossil sampling psi. N includes every living lineage, observed or hidden.
Hidden histories are numerically marginalized. The short alias is `dnDDFBDP`.
## details
Each taxon label identifies one fossil occurrence or one living sample. A fossil
may be a terminal tip or a zero-length sampled-ancestor tip. Sampling never kills
a lineage. All splits use the same speciation mechanism. Attach this tree
directly to dnPhyloCTMC and use mvCollapseExpandFossilBranch to change ancestor
status. No observation-tree projection or additional DD likelihood is needed.

Supply taxa with fossil occurrence-age bounds and present-day ages of zero.
Bounds describe uncertainty in one occurrence, not a species' stratigraphic
range. Repeated fossil specimens have distinct observation labels; they are not
preassigned to persistent species identities. There must be at least two labels.
An optional initialTree must match the taxa and ages. With no initialTree, MCMC
uses a starting-tree initializer; this is not a draw from the DD prior.
Unconditional prior redraw is not implemented and reports an error.

originAge is the age of the single-lineage origin, above the root. condition="time"
uses its unconditioned origin-started density; "sampling" conditions on at least
one fossil or living observation; "survival" conditions on at least one sampled
living lineage. There is no root-age, count, or skyline conditioning. rho is the
independent probability of observing each lineage alive at present.

At alpha=0 the time/sampling densities match constant-rate dnBDSTP with r=0.
For survival conditioning, native BDSTP divides the sampled-extant probability
by rho; this prior divides by the actual probability. Thus the native log density
is ln(rho) lower when rho<1. The density uses the same labelled, non-oriented
sampled-tree shape factor, including the legacy log-factorial approximation. A terminal fossil
marks the end of observation, not true extinction; its continuing lineage becomes
hidden and is integrated out. Increasing alpha strengthens exponential decline;
there is no hard carrying capacity and increasing birth rates are not supported.

maxHiddenLineages (default 128) is a numerical cutoff, not a biological maximum.
Outward transitions are killed, with the original diagonal hazard retained.
Increase the cutoff and compare densities and posterior results. The conditioning
calculation uses a total-count cutoff of maxHiddenLineages plus the number of
observation labels. numericalTolerance (default 1e-12) controls propagation,
not hidden-count truncation error. Large cutoffs or old trees may be expensive.

A morphology-only bear example is in examples/DDFBDP/morphology.Rev. Short MCMC
integration tests do not establish posterior convergence or parameter recovery.
## example
tree ~ dnDDFBDP(originAge=origin, lambda0=lambda0, alpha=alpha,
                mu=mu, psi=psi, rho=1, taxa=taxa, condition="sampling")
morphology ~ dnPhyloCTMC(tree=tree, Q=fnJC(2), branchRates=clock,
                        type="Standard", coding="variable")
morphology.clamp(morpho)
moves.append(mvCollapseExpandFossilBranch(tree, origin, weight=6))
## see_also
dnBDSTP
dnPhyloCTMC
mvCollapseExpandFossilBranch
dnDiversityDependentFBD
