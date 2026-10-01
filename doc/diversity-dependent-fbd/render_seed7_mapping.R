#!/usr/bin/env Rscript
# Run from the RevBayes repository root; no fitting is performed.
suppressPackageStartupMessages(library(ggplot2))
h <- jsonlite::fromJSON('validation/DDFBD/example_seed7/history.json')
e <- h$events; origin <- h$origin
sp <- data.frame(species=1:10,birth=3,endpoint=0,parent=NA_integer_)
for (i in seq_len(nrow(e))) {
 if(e$type[i]=='birth') {sp$birth[sp$species==e$child[i]]<-origin-e$time[i];sp$parent[sp$species==e$child[i]]<-e$species[i]}
 if(e$type[i]=='death') sp$endpoint[sp$species==e$species[i]]<-origin-e$time[i]
}
fossils <- data.frame(species=e$species[e$type=='fossil'],age=origin-e$time[e$type=='fossil'])
sp$observed <- sp$species %in% c(fossils$species,h$sampled_present)
sp$status <- ifelse(sp$observed,'Represented species','Completely unobserved species')
sp$panel <- 'Complete simulated history'
fossils <- rbind(transform(fossils,panel='Complete simulated history'),transform(fossils,panel='Data supplied to the analysis'))
sampled <- expand.grid(species=h$sampled_present,panel=c('Complete simulated history','Data supplied to the analysis'))
sampled$age <- 0
links <- sp[!is.na(sp$parent),]
blank <- expand.grid(species=1:10,panel=c('Complete simulated history','Data supplied to the analysis'))
p <- ggplot()+geom_blank(data=blank,aes(x=0,y=species))+
 geom_segment(data=links,aes(x=birth,xend=birth,y=parent,yend=species),colour='#C8D2DA',linetype=3,linewidth=.5)+
 geom_segment(data=sp,aes(x=birth,xend=endpoint,y=species,yend=species,colour=status),linewidth=1.2)+
 geom_point(data=sp[sp$endpoint>0,],aes(endpoint,species),shape=4,size=3,colour='#435667')+
 geom_point(data=sp[sp$endpoint==0,],aes(endpoint,species),shape=1,size=3,colour='#126A91')+
 geom_point(data=fossils,aes(age,species,shape='Fossil occurrence'),size=2.8,colour='#C56724')+
 geom_point(data=sampled,aes(age,species,shape='Sampled at the present'),size=3,colour='#126A91')+
 facet_wrap(~panel,nrow=1)+scale_x_reverse(limits=c(3.06,-.09),breaks=c(3,2,1,0))+
 scale_y_reverse(breaks=1:10,labels=paste0('sp',1:10))+
 scale_colour_manual(values=c('Represented species'='#126A91','Completely unobserved species'='#A7B4BD'))+
 scale_shape_manual(values=c('Fossil occurrence'=16,'Sampled at the present'=15))+
 labs(title='From a complete DD-FBD history to fossil and living samples',
 subtitle='Seed 7: 10 species existed; 6 appear in the data; 9 fossil records and 2 present-day samples.',
 x='Age before the present (arbitrary time units)',y=NULL,colour=NULL,shape=NULL,
 caption='Left: bars are true species lifetimes; crosses mark extinction; open circles mark living species.\nRight: only the sampled records remain. Species identity assignments are retained; hidden lifetimes are not observed.')+
 theme_minimal(base_size=13)+theme(panel.grid.minor=element_blank(),panel.grid.major.y=element_blank(),
 plot.title=element_text(face='bold',colour='#183047'),legend.position='bottom',legend.box='vertical',plot.margin=margin(12,15,10,10))
dir.create('doc/diversity-dependent-fbd/figures',showWarnings=FALSE)
ggsave('doc/diversity-dependent-fbd/figures/seed7_data_mapping.png',p,width=13,height=6.5,dpi=220,bg='white')
ggsave('doc/diversity-dependent-fbd/figures/seed7_data_mapping.svg',p,width=13,height=6.5,device=svglite::svglite,bg='white')
