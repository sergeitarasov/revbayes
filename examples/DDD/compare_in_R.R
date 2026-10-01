#!/usr/bin/env Rscript
# Run from repository root: Rscript examples/DDD/compare_in_R.R
# Optional independent reference for the Rev example. Requires DDD and ape.
suppressPackageStartupMessages(library(DDD))
tree <- ape::read.tree("examples/DDD/data/simulated.tre")
ages <- sort(ape::branching.times(tree),decreasing=TRUE)
# pars2 = hidden states, model, conditioning, density, verbose, crown/stem.
lnL <- dd_loglik(pars1=c(.8,.2,8),pars2=c(129,2,1,1,0,2),brts=ages,missnumspec=0)
cat(sprintf("DDD %s: lnL = %.14f\n",packageVersion("DDD"),lnL))
stopifnot(abs(lnL-(-7.92638283015368))<1e-7)
# Optional direct check of the Rev-generated profile (file output is rounded).
path <- "examples/DDD/output/profile.tsv"
if(file.exists(path)) {
  profile <- read.delim(path)
  profile$DDD <- vapply(profile$lambda0,function(l)
    dd_loglik(c(l,.2,8),c(129,2,1,1,0,2),ages,0),numeric(1))
  error <- max(abs(profile$lnL-profile$DDD))
  cat(sprintf("Maximum discrepancy from rounded Rev profile: %.6g\n",error))
  stopifnot(error<1e-4)
}
