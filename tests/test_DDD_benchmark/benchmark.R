#!/usr/bin/env Rscript
# Run from repository root. Uses the installed, unmodified DDD package.
args <- commandArgs(TRUE)
out <- if(length(args)) args[1] else "tests/test_DDD_benchmark/results"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
stopifnot(requireNamespace("DDD",quietly=TRUE),requireNamespace("ape",quietly=TRUE))
writeLines(c(capture.output(sessionInfo()),capture.output(packageDescription("DDD"))),file.path(out,"session.txt"))
writeLines(capture.output(lapply(c("lambdamu","dd_loglik2","dd_loglik_M_aux"),function(x)get(x,asNamespace("DDD")))),file.path(out,"DDD_reference_functions.txt"))
set.seed(1439)
simulation <- DDD::dd_sim(c(0.8,0.2,8),age=3,ddmodel=2)
saveRDS(simulation,file.path(out,"simulation.rds"))
ape::write.tree(simulation$tes,file.path(out,"simulated.tre"),digits=16)
fixed <- ape::read.tree(text="((A:1,B:1):2,(C:2,D:2):1);")
ape::write.tree(fixed,file.path(out,"fixed.tre"),digits=16)
trees <- list(fixed=fixed,simulated=simulation$tes)
fmt <- function(x) format(x,digits=17,scientific=FALSE,trim=TRUE)
kernel <- function(model,lambda,mu,K,crown,condition,cap,origin,ages,alpha=0) {
  a <- c(model,fmt(lambda),fmt(mu),fmt(alpha),fmt(K),as.integer(crown),condition,cap,fmt(origin),fmt(ages))
  value <- system2(".local-build/ddd-kernel",a,stdout=TRUE)
  if(!is.null(attr(value,"status"))) stop(paste(value,collapse="\n"))
  as.numeric(value)
}
rows <- list(); rev <- c('max_error = 0.0')
for(tree_name in names(trees)) {
  ages <- sort(ape::branching.times(trees[[tree_name]]),decreasing=TRUE)
  rev <- c(rev,sprintf('tree = readTrees("%s/%s.tre")[1]',out,tree_name))
  for(model in c("DDDlinear","DDDpower")) for(par in list(c(.8,.2,8),c(.5,.35,5),c(1.2,.1,15)))
    for(crown in c(TRUE,FALSE)) for(condition in 0:2) for(cap in c(64,128)) {
      origin <- max(ages)+if(crown) 0 else 1
      brts <- if(crown) ages else c(origin,ages)
      ddd <- DDD::dd_loglik(par,c(cap+1,if(model=="DDDlinear")1 else 2,condition,1,0,if(crown)2 else 1),brts,0)
      rb <- kernel(model,par[1],par[2],par[3],crown,condition,cap,origin,ages)
      rows[[length(rows)+1]] <- data.frame(tree=tree_name,model=model,lambda0=par[1],mu=par[2],K=par[3],
        start=if(crown)"crown" else "stem",condition=condition,cap=cap,DDD=ddd,kernel=rb,error=rb-ddd)
      rev <- c(rev,sprintf('ll := fnDiversityDependentLogLikelihood(tree,lambda0=%s,mu=%s,K=%s,rateModel="%s",start="%s",originAge=%s,condition="%s",maxHiddenLineages=%d,numericalTolerance=1e-13)',
        fmt(par[1]),fmt(par[2]),fmt(par[3]),model,if(crown)"crown" else "stem",fmt(origin),c("none","survival","nTaxa")[condition+1],cap),
        sprintf('err = abs(ll - (%s))',fmt(ddd)),
        'if (err > 1e-7) { stop("DDD benchmark mismatch") }',
        'if (err > max_error) { max_error = err }')
    }
}
results <- do.call(rbind,rows)
write.table(results,file.path(out,"likelihoods.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
stopifnot(all(is.finite(results$error)),max(abs(results$error))<1e-7)
# A second independent DDD numerical backend (ODE rather than matrix exponentiation).
ages <- sort(ape::branching.times(simulation$tes),decreasing=TRUE)
p <- c(.8,.2,8)
ode <- DDD::dd_loglik(p,c(129,2,1,1,0,2,abstolint=1e-12,reltolint=1e-12),ages,0,methode="odeint::runge_kutta_cash_karp54")
ours <- kernel("DDDpower",p[1],p[2],p[3],TRUE,1,128,max(ages),ages)
stopifnot(abs(ode-ours)<1e-7)
# Match a one-parameter ML fit: mu and K fixed. Not a parameter-recovery study.
fitDDD <- DDD::dd_ML(ages,initparsopt=.8,idparsopt=1,idparsfix=c(2,3),parsfix=c(.2,8),
  res=129,ddmodel=2,cond=1,btorph=1,soc=2,tol=c(1e-8,1e-8,1e-8),verbose=FALSE)
fitKernel <- optimize(function(l) -kernel("DDDpower",l,.2,8,TRUE,1,128,max(ages),ages),
  interval=c(.200001,5),tol=1e-8)
write.table(fitDDD,file.path(out,"DDD_fit.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
fit <- data.frame(DDD_lambda=fitDDD$lambda, kernel_lambda=fitKernel$minimum,
                 DDD_loglik=fitDDD$loglik, kernel_loglik=-fitKernel$objective)
write.table(fit,file.path(out,"fit_comparison.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
stopifnot(abs(fit$DDD_lambda-fit$kernel_lambda)<1e-4,abs(fit$DDD_loglik-fit$kernel_loglik)<1e-7)
# Verify the branching-time normalization and live parameter updates in Rev.
rev <- c(rev,
  sprintf('bt := fnDiversityDependentLogLikelihood(tree,lambda0=0.8,mu=0.2,K=8,rateModel="DDDpower",density="branchingTimes")'),
  sprintf('ph := fnDiversityDependentLogLikelihood(tree,lambda0=0.8,mu=0.2,K=8,rateModel="DDDpower")'),
  sprintf('if (abs(bt-ph-(%s)) > 1e-8) { stop("Density convention mismatch") }',fmt(lgamma(length(ages)+1))),
  'lambda ~ dnUniform(0.21,2.0)', 'lambda.setValue(0.8)',
  'dd_live := fnDiversityDependentLogLikelihood(tree,lambda0=lambda,mu=0.2,K=8,rateModel="DDDpower")',
  'before = dd_live', 'lambda.setValue(1.1)',
  'after := fnDiversityDependentLogLikelihood(tree,lambda0=1.1,mu=0.2,K=8,rateModel="DDDpower")',
  'if (abs(dd_live-after) > 1e-8) { stop("DAG update mismatch") }',
  'if (abs(dd_live-before) < 1e-6) { stop("DAG parameter did not update") }',
  'print("DDD_REV_CHECKS_PASSED max_error=" + max_error)', 'q()')
# Evaluate the fitted optimum through the actual Rev function as well.
rev <- append(rev,c(sprintf('fitted := fnDiversityDependentLogLikelihood(tree,lambda0=%s,mu=0.2,K=8,rateModel="DDDpower")',fmt(fitKernel$minimum)),
  sprintf('if (abs(fitted-(%s)) > 1e-7) { stop("Fitted Rev likelihood mismatch") }',fmt(fitDDD$loglik))),after=length(rev)-2)
writeLines(rev,file.path(out,"compare.Rev"))
# DDD uses a modified top diagonal in its matrix backend; ours kills overflow.
# A deliberately difficult case demonstrates convergence, not equality at low caps.
sweep <- do.call(rbind,lapply(c(8,16,32,64,128,256),function(cap) {
  rb <- kernel("DDDpower",1.2,.8,30,TRUE,1,cap,10,c(10,6,3))
  ddd <- DDD::dd_loglik(c(1.2,.8,30),c(cap+1,2,1,1,0,2),c(10,6,3),0)
  data.frame(cap=cap,kernel=rb,DDD=ddd,error=rb-ddd)
}))
write.table(sweep,file.path(out,"cutoff_convergence.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
stopifnot(abs(tail(sweep$error,1))<1e-7,abs(diff(tail(sweep$kernel,2)))<1e-7)
constant <- do.call(rbind,lapply(c(TRUE,FALSE),function(crown) do.call(rbind,lapply(0:2,function(cond) {
  origin <- if(crown)3 else 4
  ddd <- DDD::dd_loglik(c(.8,.2,Inf),c(129,1,cond,1,0,if(crown)2 else 1),if(crown)c(3,2,1) else c(4,3,2,1),0)
  rb <- kernel("exponential",.8,.2,8,crown,cond,512,origin,c(3,2,1),alpha=0)
  data.frame(crown=crown,condition=cond,DDD=ddd,kernel=rb,error=rb-ddd)
}))))
write.table(constant,file.path(out,"constant_rate.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
stopifnot(max(abs(constant$error))<1e-7)
summary <- c(sprintf("DDD %s; %d matched likelihood comparisons",packageVersion("DDD"),nrow(results)),
  sprintf("Maximum absolute log-likelihood discrepancy: %.15g",max(abs(results$error))),
  sprintf("Independent ODE-backend discrepancy: %.15g",abs(ode-ours)),
  sprintf("Constant-rate maximum discrepancy (H=512): %.15g",max(abs(constant$error))),
  sprintf("Stress-case DDD discrepancy (H=256): %.15g",abs(tail(sweep$error,1))),
  sprintf("Stress-case kernel change H=128 to 256: %.15g",abs(diff(tail(sweep$kernel,2)))),
  sprintf("Simulated tree: %d extant tips, seed 1439",length(simulation$tes$tip.label)),capture.output(fit))
writeLines(summary,file.path(out,"summary.txt")); cat(summary,sep="\n")
