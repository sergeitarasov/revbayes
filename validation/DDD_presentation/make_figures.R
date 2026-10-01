#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ape))
root <- "validation/DDD_presentation"
out <- file.path(root,"results"); figs <- file.path(root,"figures")
dir.create(figs,recursive=TRUE,showWarnings=FALSE)
read <- function(name) read.delim(file.path(out,name),check.names=FALSE)
sims <- read("simulations.tsv"); comp <- read("likelihood_comparisons.tsv"); fits <- read("fits.tsv")
meta <- read("illustrative_metadata.tsv"); profile <- read("likelihood_profile.tsv"); ltt <- read("illustrative_diversity.tsv")
stopifnot(nrow(sims)==100,nrow(comp)==300,nrow(fits)==100,all(ltt$total_diversity>=ltt$reconstructed_lineages))
blue <- "#126A91"; orange <- "#C56724"; green <- "#298168"; ink <- "#183047"
cols <- c(DDDlinear=orange,DDDpower=blue)
labels <- c(DDDlinear="Linear (DDD 1)",DDDpower="Power law (DDD 2)")
theme_set(theme_minimal(base_size=14,base_family="Helvetica")+theme(
  plot.title=element_text(face="bold",colour=ink,size=17),plot.subtitle=element_text(size=12,colour="#526779"),
  panel.grid.minor=element_blank(),panel.grid.major=element_line(colour="#E5EAEE",linewidth=.35),
  axis.title=element_text(colour=ink),axis.text=element_text(colour="#46596B"),
  legend.position="bottom",legend.title=element_blank(),plot.margin=margin(10,15,8,10)))
save_plot <- function(plot,name,w=12,h=5.3) {
  ggsave(file.path(figs,paste0(name,".png")),plot,width=w,height=h,dpi=260,bg="white")
  ggsave(file.path(figs,paste0(name,".svg")),plot,width=w,height=h,device=svglite::svglite,bg="white")
}
N <- seq(1,30,length.out=300)
rates <- rbind(data.frame(N=N,rate=pmax(0,.8-(.8-.2)*N/15),model="Linear (DDD 1)"),
 data.frame(N=N,rate=.8*(N+1)^(-log(.8/.2)/log(15+1)),model="Power law (DDD 2)"),
 data.frame(N=N,rate=.8*exp(-log(4)/(15-1)*(N-1)),model="Exponential in N"))
p <- ggplot(rates,aes(N,rate,colour=model,linetype=model))+geom_line(linewidth=1.2)+
 geom_hline(yintercept=.2,colour=ink,linetype=3)+geom_vline(xintercept=15,colour="#A7B3BE",linetype=3)+
 annotate("text",x=26,y=.235,label="constant extinction = 0.2",size=4,colour=ink)+
 scale_colour_manual(values=c("Linear (DDD 1)"=orange,"Power law (DDD 2)"=blue,"Exponential in N"="#788795"))+
 scale_linetype_manual(values=c("Linear (DDD 1)"=1,"Power law (DDD 2)"=1,"Exponential in N"=2))+
 labs(x="Total living diversity N (observed + hidden)",y="Per-lineage speciation rate",
 title="Different rate laws, shared dependence on total diversity",subtitle="Simulation models: lambda0 = 0.8, mu = 0.2, K = 15. Grey curve: alpha = log(4)/14; shown for contrast only.")
save_plot(p,"01_rate_laws")
# Reconstructed tree in forward time; tip ordering follows the ape tree.
tree <- ape::read.tree(file.path(out,"illustrative.tre")); ntip <- length(tree$tip.label)
x <- ape::node.depth.edgelength(tree); y <- rep(NA,length(x)); y[1:ntip] <- seq_len(ntip)
place <- function(node) { if(is.na(y[node])) y[node]<<-mean(vapply(tree$edge[tree$edge[,1]==node,2],place,numeric(1))); y[node] }
place(ntip+1)
edge <- data.frame(x=x[tree$edge[,1]],xend=x[tree$edge[,2]],y=y[tree$edge[,2]],yend=y[tree$edge[,2]])
internal <- unique(tree$edge[,1]); vertical <- do.call(rbind,lapply(internal,function(n){ kids<-tree$edge[tree$edge[,1]==n,2]; data.frame(x=x[n],xend=x[n],y=min(y[kids]),yend=max(y[kids])) }))
ptree <- ggplot()+geom_segment(data=edge,aes(x=x,xend=xend,y=y,yend=yend),colour=blue,linewidth=.65)+
 geom_segment(data=vertical,aes(x=x,xend=xend,y=y,yend=yend),colour=blue,linewidth=.65)+
 geom_point(data=data.frame(x=x[1:ntip],y=y[1:ntip]),aes(x,y),size=1.6,colour=blue)+
 scale_x_continuous(limits=c(0,8.2),breaks=c(0,2,4,6,8))+labs(x="Time since crown origin",y=NULL,title=sprintf("Observed tree: %d living species",ntip),subtitle=sprintf("Illustrative history: %s",meta$id))+
 theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),panel.grid.major.y=element_blank())
lttlong <- rbind(data.frame(time=ltt$time,N=ltt$total_diversity,series="All living lineages"),data.frame(time=ltt$time,N=ltt$reconstructed_lineages,series="Ancestral lineages of extant tips"))
pltt <- ggplot(lttlong,aes(time,N,colour=series))+geom_step(linewidth=1.05)+
 scale_colour_manual(values=c("All living lineages"=orange,"Ancestral lineages of extant tips"=blue))+
 scale_x_continuous(breaks=c(0,2,4,6,8))+labs(x="Time since crown origin",y="Number of living lineages",title="The likelihood integrates over hidden history",subtitle=sprintf("%d species originated; %d became extinct.",meta$total_species,meta$extinct_species))+
 theme(legend.text=element_text(size=10))
save_plot(plot_grid(ptree,pltt,nrow=1,rel_widths=c(1,1.1)),"02_tree_and_hidden_diversity",13,5.3)
# Independent numerical agreement.
pagree <- ggplot(comp,aes(DDD,Rev_kernel,colour=model,shape=factor(age)))+
 geom_abline(slope=1,intercept=0,colour="#ADB9C2",linewidth=.7)+geom_point(size=2.3,alpha=.72)+
 scale_colour_manual(values=cols,labels=labels)+scale_shape_manual(values=c(16,17),labels=c("Crown age 3","Crown age 8"))+
 labs(x="DDD log likelihood",y="RevBayes kernel log likelihood",title="300 matched parameter-point comparisons",subtitle="100 trees; three parameter points per tree.\nComplete sampling and crown survival conditioning.")
perror <- ggplot(comp,aes(DDD,pmax(abs(error),1e-16),colour=model,shape=factor(age)))+
 geom_hline(yintercept=1e-7,linetype=2,colour="#963E3E")+geom_point(size=2.3,alpha=.7)+
 scale_y_log10(limits=c(1e-16,1e-6),breaks=c(1e-16,1e-13,1e-10,1e-7))+
 scale_colour_manual(values=cols,labels=labels)+scale_shape_manual(values=c(16,17),labels=c("Crown age 3","Crown age 8"))+
 labs(x="DDD log likelihood",y="Absolute log-likelihood difference",title="Numerical differences remain small",subtitle=sprintf("Maximum = %.2g; dashed line: tolerance 1e-7.\nExact zero differences are plotted at 1e-16.",max(abs(comp$error))))
legend <- get_legend(pagree+theme(legend.position="bottom"))
save_plot(plot_grid(plot_grid(pagree+theme(legend.position="none"),perror+theme(legend.position="none"),nrow=1),legend,ncol=1,rel_heights=c(1,.12)),"03_likelihood_agreement",13,5.5)
selected <- fits[fits$id==meta$id,]
baseline <- max(selected$Rev_lnL,profile$Rev_kernel)
p <- ggplot(profile,aes(lambda0,Rev_kernel-baseline))+geom_line(colour=blue,linewidth=1.25)+
 geom_point(data=profile[seq(1,nrow(profile),5),],aes(y=DDD-baseline),shape=1,size=2.6,stroke=1.1,colour=orange)+
 geom_vline(xintercept=.8,linetype=2,colour=green)+geom_vline(xintercept=selected$Rev_lambda,linetype=3,colour=ink)+
 scale_x_log10(breaks=c(.2,.4,.8,1.5,3,6,8))+
 labs(x=expression(lambda[0]~"(logarithmic axis; mu and K fixed)"),y="Log likelihood relative to maximum",
 title="The two likelihood profiles overlap",subtitle=sprintf("Blue line: RevBayes; orange circles: DDD. Generating lambda0 = 0.8 (green); fitted lambda0 = %.3f (grey).",selected$Rev_lambda))
save_plot(p,"04_likelihood_profile",12,5.0)
# Pilot estimation experiment: keep boundary estimates visible.
fits$age <- factor(fits$age); fits$model <- factor(fits$model,levels=c("DDDlinear","DDDpower"))
set.seed(173)
precover <- ggplot(fits,aes(age,Rev_lambda,colour=model))+
 geom_hline(yintercept=.8,colour=green,linetype=2,linewidth=.75)+
 geom_boxplot(width=.48,outlier.shape=NA,fill="white",linewidth=.65)+
 geom_point(aes(shape=boundary),position=position_jitter(width=.14,height=0,seed=173),size=2.1,alpha=.8)+
 facet_wrap(~model,nrow=1,labeller=as_labeller(labels))+
 scale_colour_manual(values=cols,guide="none")+
 scale_shape_manual(values=c(interior=16,lower=25,upper=24),drop=FALSE)+
 scale_y_log10(breaks=c(.2,.4,.8,1.6,3.2,8),limits=c(.18,9))+
 labs(x="Crown age",y=expression("Fitted "*lambda[0]*" (logarithmic axis)"),title="Matching software does not remove sampling variation",subtitle="25 trees per group. Green: generating value 0.8. Triangles: search-bound estimates.")
psize <- ggplot(sims,aes(factor(age),extant_tips,colour=model))+
 geom_boxplot(width=.48,outlier.shape=NA,fill="white",linewidth=.65)+
 geom_point(position=position_jitter(width=.14,height=0,seed=173),size=2.1,alpha=.8)+
 facet_wrap(~model,nrow=1,labeller=as_labeller(labels))+
 scale_colour_manual(values=cols,guide="none")+
 labs(x="Crown age",y="Number of extant tips",title="Information varies with realized tree size",subtitle="All accepted simulated trees retained; no filtering by tip count.")
save_plot(plot_grid(precover,psize,nrow=1,rel_widths=c(1.2,1)),"05_estimation_variation",14,5.7)
cut <- read.delim("tests/test_DDD_benchmark/results/cutoff_convergence.tsv")
cutlong <- rbind(data.frame(cap=cut$cap,lnL=cut$kernel,method="RevBayes"),data.frame(cap=cut$cap,lnL=cut$DDD,method="DDD"))
p <- ggplot(cutlong,aes(cap,lnL,colour=method,shape=method))+geom_line(linewidth=.9)+geom_point(size=3)+
 scale_x_log10(breaks=cut$cap)+scale_colour_manual(values=c(DDD=orange,RevBayes=blue))+
 labs(x="Maximum hidden lineages H (logarithmic axis)",y="Log likelihood",title="Converge the hidden-lineage cutoff before interpreting results",
 subtitle="Stress case: power law, lambda0 = 1.2, mu = 0.8, K = 30; crown age 10; internal ages 10, 6, 3.")
save_plot(p,"06_cutoff_convergence",12,5.0)
# Summary quantities used verbatim by the report builder.
summary <- do.call(rbind,lapply(split(seq_len(nrow(fits)),interaction(fits$model,fits$age,drop=TRUE)),function(ii) {
 f<-fits[ii,]; data.frame(model=as.character(f$model[1]),age=as.numeric(as.character(f$age[1])),n=nrow(f),
 median_tips=median(f$extant_tips),min_tips=min(f$extant_tips),max_tips=max(f$extant_tips),
 median_lambda=median(f$Rev_lambda),mean_lambda=mean(f$Rev_lambda),q25=unname(quantile(f$Rev_lambda,.25)),q75=unname(quantile(f$Rev_lambda,.75)),
 rmse=sqrt(mean((f$Rev_lambda-.8)^2)),lower=sum(f$boundary=="lower"),upper=sum(f$boundary=="upper"))
}))
write.table(summary,file.path(out,"group_summary.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
cat("SIX_FIGURES_COMPLETE\n")
