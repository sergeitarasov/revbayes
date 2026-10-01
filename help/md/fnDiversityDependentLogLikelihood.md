## name
fnDiversityDependentLogLikelihood

## title
Diversity-dependent likelihood of a complete extant tree

## description
Returns the log likelihood of a rooted binary time tree whose tips are all at the present. Speciation depends on total diversity, including extinct lineages integrated out by a hidden-count calculation; extinction is constant. This is a likelihood evaluation function, not a tree-generating distribution.

## details
The default rate model is `exponential`: lambda(N) = lambda0 * exp(-alpha * (N-1)). `DDDlinear` uses max(0, lambda0 - (lambda0-mu)*N/K); `DDDpower` uses lambda0 * (N+1)^(-log(lambda0/mu)/log(K+1)). The latter two match DDD models 1 and 2, respectively, and require lambda0 > mu > 0 and K > 0.

`start="crown"` starts with two lineages at the root and omits the root speciation factor. `start="stem"` starts with one lineage at the supplied `originAge` and includes every internal speciation event. `condition` is `none`, `survival`, or `nTaxa`, matching DDD conditions 0, 1, or 2. Crown survival requires both initial sides to survive.

`density="DDDphylogeny"` matches DDD's btorph=1 convention. `branchingTimes` adds log((n-1)!). This DDD convention must not be silently substituted for another labelled-tree prior convention.

All extant species must be sampled. There are no fossil observations or missing extant species in this function. The hidden-count boundary is killed; increase `maxHiddenLineages` until the result converges. DDD's matrix backend uses a different top-boundary diagonal, so comparisons require convergence of both calculations. `numericalTolerance` controls propagation accuracy, not truncation error. A deterministic value assigned with `:=` updates with its DAG parameters but does not by itself contribute to an MCMC target density.

## example
    tree = readTrees("tree.tre")[1]
    lnL := fnDiversityDependentLogLikelihood(tree, lambda0=0.8, mu=0.2,
               K=8, rateModel="DDDpower", start="crown", condition="survival")
    print(lnL)

## see_also
dnDiversityDependentFBD
