#!/usr/bin/env Rscript
# From repository root. All simulated trees are retained, including two-tip trees.
suppressPackageStartupMessages(library(DDD))
suppressPackageStartupMessages(library(ape))
args <- commandArgs(TRUE)
reps <- if(length(args)) as.integer(args[1]) else 25L
out <- "validation/DDD_presentation/results"
dir.create(file.path(out,"trees"),recursive=TRUE,showWarnings=FALSE)
stopifnot(file.exists(".local-build/ddd-kernel"),reps>0)
options(digits=17)
fmt <- function(x) format(x,digits=17,scientific=FALSE,trim=TRUE)
write_tsv <- function(x,name) write.table(x,file.path(out,name),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(c(capture.output(sessionInfo()),capture.output(packageDescription("DDD"))),file.path(out,"session.txt"))
truth <- c(lambda0=.8,mu=.2,K=15)
bounds <- c(.20001,8)
# Same scalar optimizer, search bounds and coarse bracketing for both likelihoods.
# Endpoints are explicit candidates, so boundary estimates are preserved.
fit_one <- function(fn) {
  grid <- seq(log(bounds[1]),log(bounds[2]),length.out=13)
  y <- vapply(grid,function(x)fn(exp(x)),numeric(1))
  j <- which.max(y)
  fit <- optimize(function(x)-fn(exp(x)),c(grid[max(1,j-1)],grid[min(length(grid),j+1)]),tol=1e-7)
  candidates <- c(grid[j],fit$minimum,grid[1],tail(grid,1))
  values <- vapply(candidates,function(x)fn(exp(x)),numeric(1))
  winner <- which.max(values)
  estimate <- exp(candidates[winner])
  list(lambda=estimate,lnL=values[winner],boundary=if(estimate<=bounds[1]*1.00001)"lower" else if(estimate>=bounds[2]/1.00001)"upper" else "interior")
}
kernel <- function(ages,model,lambda,cap=128) {
  value <- system2(".local-build/ddd-kernel",c(model,fmt(lambda),fmt(truth['mu']),0,fmt(truth['K']),1,1,cap,fmt(max(ages)),fmt(ages)),stdout=TRUE)
  if(!is.null(attr(value,"status")))stop(paste(value,collapse="\n"))
  ans <- as.numeric(value); stopifnot(!is.na(ans),ans < Inf); ans
}
reference <- function(ages,model,lambda,cap=128) {
  ans <- DDD::dd_loglik(c(lambda,truth['mu'],truth['K']),c(cap+1,model,1,1,0,2),ages,0)
  stopifnot(!is.na(ans),ans < Inf); ans
}
cached_fits <- if(file.exists(file.path(out,"fits.tsv"))) read.delim(file.path(out,"fits.tsv")) else data.frame(id=character())
records <- comparisons <- fits <- list()
for(model in 1:2) for(age in c(3,8)) for(rep in seq_len(reps)) {
  name <- if(model==1)"DDDlinear" else "DDDpower"
  id <- sprintf("%s_age%d_rep%02d",name,age,rep)
  seed <- 20261001L+model*10000L+age*100L+rep
  saved <- file.path(out,"trees",paste0(id,".rds"))
  if(file.exists(saved)) {
    bundle <- readRDS(saved)
    stopifnot(identical(bundle$seed,seed),identical(bundle$truth,truth),bundle$model==model,bundle$age==age)
    sim <- bundle$simulation
  } else {
    set.seed(seed)
    sim <- DDD::dd_sim(truth,age,ddmodel=model)
    saveRDS(list(seed=seed,truth=truth,age=age,model=model,simulation=sim),saved)
  }
  ape::write.tree(sim$tes,file.path(out,"trees",paste0(id,".tre")),digits=16)
  ages <- sort(ape::branching.times(sim$tes),decreasing=TRUE)
  records[[length(records)+1]] <- data.frame(id=id,model=name,age=age,replicate=rep,seed=seed,
    extant_tips=length(sim$tes$tip.label),total_species=nrow(sim$L),extinct_species=sum(sim$L[,4]>=0))
  for(lambda in if(model==1)c(.4,.6,.8) else c(.6,.8,1.2)) {
    ddd <- reference(ages,model,lambda)
    rb <- kernel(ages,name,lambda)
    comparisons[[length(comparisons)+1]] <- data.frame(id=id,model=name,age=age,lambda0=lambda,DDD=ddd,Rev_kernel=rb,error=rb-ddd)
  }
  if(id %in% cached_fits$id) {
    f <- cached_fits[cached_fits$id==id,]
  } else {
  fddd <- fit_one(function(l)reference(ages,model,l))
  frb <- fit_one(function(l)kernel(ages,name,l))
  cutoff_delta <- kernel(ages,name,frb$lambda,256)-frb$lnL
  f <- data.frame(id=id,model=name,age=age,replicate=rep,extant_tips=length(sim$tes$tip.label),
    DDD_lambda=fddd$lambda,Rev_lambda=frb$lambda,DDD_lnL=fddd$lnL,Rev_lnL=frb$lnL,
    lambda_difference=frb$lambda-fddd$lambda,lnL_difference=frb$lnL-fddd$lnL,
    boundary=frb$boundary,DDD_boundary=fddd$boundary,cutoff_delta=cutoff_delta)
  }
  fits[[length(fits)+1]] <- f
  # Save completed numerical work incrementally.
  write_tsv(do.call(rbind,records),"simulations.tsv")
  write_tsv(do.call(rbind,comparisons),"likelihood_comparisons.tsv")
  write_tsv(do.call(rbind,fits),"fits.tsv")
  cat(sprintf("%3d/%d %s: n=%d, lambda_hat=%.6f (%s), deltaLnL=%.2g\n",length(records),4*reps,id,length(sim$tes$tip.label),f$Rev_lambda,f$boundary,f$lnL_difference)); flush.console()
}
records <- do.call(rbind,records); comparisons <- do.call(rbind,comparisons); fits <- do.call(rbind,fits)
stopifnot(nrow(records)==4*reps,max(abs(comparisons$error))<1e-7,max(abs(fits$lnL_difference))<1e-7,max(abs(fits$cutoff_delta))<1e-7)
# Illustrative tree chosen reproducibly as closest to the median number of tips
# in the power-law, age-8 group. No simulation is removed from any analysis.
pool <- subset(records,model=="DDDpower" & age==8)
chosen <- pool[which.min(abs(pool$extant_tips-median(pool$extant_tips))),]
bundle <- readRDS(file.path(out,"trees",paste0(chosen$id,".rds")))
sim <- bundle$simulation
file.copy(file.path(out,"trees",paste0(chosen$id,".tre")),file.path(out,"illustrative.tre"),overwrite=TRUE)
write_tsv(chosen,"illustrative_metadata.tsv")
ages <- sort(ape::branching.times(sim$tes),decreasing=TRUE)
profile <- do.call(rbind,lapply(exp(seq(log(bounds[1]),log(bounds[2]),length.out=81)),function(lambda) {
 data.frame(lambda0=lambda,DDD=reference(ages,2,lambda),Rev_kernel=kernel(ages,"DDDpower",lambda))
}))
write_tsv(profile,"likelihood_profile.tsv")
# Full realized diversity versus the reconstructed number of ancestral lineages.
time <- sort(unique(c(seq(0,8,length.out=401),8-sim$L[,1],8-sim$L[sim$L[,4]>=0,4],8-ages)))
ltt <- do.call(rbind,lapply(time,function(t) {
  ago <- 8-t
  full <- sum(sim$L[,1]>=ago-1e-10 & (sim$L[,4]<0 | sim$L[,4]<ago-1e-10))
  observed <- 2+sum(ages[-1]>=ago-1e-10)
  data.frame(time=t,total_diversity=full,reconstructed_lineages=observed)
}))
write_tsv(ltt,"illustrative_diversity.tsv")
# Actual Rev checks at every comparison point, fitted point and profile point.
rev <- c('max_error = 0.0','check_count = 0')
add_check <- function(lambda,expected,model) {
 c(sprintf('lnL := fnDiversityDependentLogLikelihood(tree,lambda0=%s,mu=0.2,K=15,rateModel="%s",start="crown",condition="survival",maxHiddenLineages=128,numericalTolerance=1e-13)',fmt(lambda),model),
 sprintf('err = abs(lnL-(%s))',fmt(expected)),
 'if (err > 1e-7) { stop("Simulation comparison failed") }',
 'if (err > max_error) { max_error = err }','check_count = check_count + 1')
}
for(i in seq_len(nrow(records))) {
 id <- records$id[i]
 rev <- c(rev,sprintf('tree = readTrees("%s/trees/%s.tre")[1]',out,id))
 cases <- subset(comparisons,id==records$id[i])
 for(j in seq_len(nrow(cases))) rev <- c(rev,add_check(cases$lambda0[j],cases$DDD[j],cases$model[j]))
 fit <- subset(fits,id==records$id[i])
 # At the Rev optimum, compare against the independent DDD value at the same parameter.
 a <- sort(ape::branching.times(readRDS(file.path(out,"trees",paste0(id,".rds")))$simulation$tes),decreasing=TRUE)
 expected <- reference(a,if(records$model[i]=="DDDlinear")1 else 2,fit$Rev_lambda)
 rev <- c(rev,add_check(fit$Rev_lambda,expected,records$model[i]))
}
rev <- c(rev,sprintf('tree = readTrees("%s/illustrative.tre")[1]',out))
for(i in seq_len(nrow(profile))) rev <- c(rev,add_check(profile$lambda0[i],profile$DDD[i],"DDDpower"))
rev <- c(rev,'print("SIMULATION_REV_CHECKS_PASSED count=" + check_count + " max_error=" + max_error)','q()')
writeLines(rev,file.path(out,"verify_in_Rev.Rev"))
writeLines(c(sprintf("replicates_per_scenario=%d",reps),"lambda0=0.8; mu=0.2; K=15; crown_age=3 or 8", "rate_models=DDDlinear, DDDpower; crown survival conditioning; complete extant sampling", "lambda_search=[0.20001,8]; mu and K fixed; all accepted dd_sim trees retained",sprintf("reference_DDD_version=%s",packageVersion("DDD"))),file.path(out,"design.txt"))
cat("SIMULATION_STUDY_COMPLETE\n")
