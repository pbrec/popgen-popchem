##Custom R script for genome-wide diversity analysis performed in Brand et al. 
##Title: "The evolution of sexual signaling is linked to odorant receptor tuning in perfume-collecting orchid bees"

##script requires snp data in vcf format
##files can be computed by mapping whole-genome sequencing reads to a reference genome and calling SNPs afterwards
##request VCF files (very large) by emailing Philipp Brand (pbrand@ucdavis.edu) or Santiago Ramirez (sanram@ucdavis.edu)

library(PopGenome)

snp<-readData("popgenome-vcf",format = "VCF",include.unknown = T,gffpath = "GFF") 

get.sum.data(snp)
show.slots(snp)

snp@n.sites

#ADD POPULATIONS
#pops<-get.individuals(snp)[1]
#viridissima
pop1<-c("PB0679","PB0738","PB0891","PB0909","PB0952","SR3371","SR3390","S290","S289","S291")
#dilemma
pop2<-c("PB0206","PB0209","PB0210","PB0471","PB0482","PB0485","PB0490","PB0709","PB0718","PB0740","PB0749","PB0754","PB0856","PB0862","PB0864","PB0900","PB0921","S286","S288","S287")

#allopatric dilemma 
da<-c("PB0206","B0209","PB0210","PB0471","PB0482","PB0485","PB0490","S286","S288","S287")
#sympatric dilemma
ds<-c("PB0709","PB0718","PB0740","PB0749","PB0754","PB0856","PB0862","PB0864","PB0900","PB0921")
snp  <- set.populations(snp, list(pop1,pop2,da,ds), diploid=F) 
snp@populations

#CALCULATE FST STATS
snp <- F_ST.stats(snp,FAST = T) 
head(get.F_ST(snp))
snp@nucleotide.F_ST

#CALCULATE NEUTRALITY STATS - takes a long time
snp<-neutrality.stats(snp,FAST = T)
get.neutrality(snp)
#population-wise
head(get.neutrality(snp)[[1]])
colMeans(snp@Tajima.D) #mean Tajima's D 

# Print diversities
get.diversity(snp)
get.diversity(snp)[[1]] # dilemma
get.diversity(snp)[[2]] # viridissima
get.diversity(snp)[[3]] # da
get.diversity(snp)[[3]] # ds


#SLIDING WINDOW STATS
winsize=50000
win_snp <- sliding.window.transform(snp, 
                                    width=winsize, jump=winsize, 
                                    type =  2,
                                    whole.data=F) # used False for whole genome
#retrieve genome positions corresponding to windows

genome.pos <- sapply(win_snp@region.names, function(x){
  split <- strsplit(x," ")[[1]][c(1,2,4)]
  val   <- mean(as.numeric(split[2:3]))
  return(val)
 
  })

names(genome.pos)

# Measurements per window
win_snp <- F_ST.stats(win_snp,FAST = T)
win_snp<-neutrality.stats(win_snp,FAST = T)


get.neutrality(win_snp)[[1]]

#win_snp@nucleotide.F_ST
#win_snp@nuc.diversity.within
colMeans(win_snp@nuc.diversity.within)/winsize #mean π for each population
colMeans(win_snp@Tajima.D) #mean Tajima's D for each population

#calculate some parameters
pairwise.FST <- t(win_snp@nuc.F_ST.pairwise)

v_div<- win_snp@nuc.diversity.within[,1] 
d_div<- win_snp@nuc.diversity.within[,2] 
da_div<-win_snp@nuc.diversity.within[,3] 
ds_div<-win_snp@nuc.diversity.within[,4] 

dv_FST<-pairwise.FST[,1]
dd_FST<-pairwise.FST[,6]
dv_FST[dv_FST<0]<-0
dd_FST[dd_FST<0]<-0
deltaFST<-dv_FST-dd_FST
plot(deltaFST)
deltaFSTnon0<-deltaFST
deltaFST[deltaFST<0]<-0
dav_FST<-pairwise.FST[,2]
dsv_FST<-pairwise.FST[,3]
dav_FST[dav_FST<0]<-0
dsv_FST[dsv_FST<0]<-0

#z-transform
dav_FSTz<-scale(dav_FST,center=T,scale=T)
dsv_FSTz<-scale(dsv_FST,center=T,scale=T)
dd_FSTz<-scale(dd_FST,center=T,scale=T)
dv_FSTz<-scale(dv_FST,center=T,scale=T)
deltaFSTz<-dv_FSTz-dd_FSTz

#make colors - all black, 99 perc FST outliers red. Probably not all needed in the end
col<-rep(1,length(dd_FST))
col_delta<-col
col_delta[deltaFST>quantile(deltaFST,.99,na.rm=T)]<-2
col_dv<-col
col_dv[dv_FST>quantile(dv_FST,.99,na.rm=T)]<-2
col_dd<-col
col_dd[dd_FST>quantile(dd_FST,.99,na.rm=T)]<-2
col_dav<-col
col_dav[dav_FST>quantile(dav_FST,.99,na.rm=T)]<-2
col_dsv<-col
col_dsv[dsv_FST>quantile(dsv_FST,.99,na.rm=T)]<-2
col<-rep(1,length(dd_FSTz))
col_dvz<-col
col_dvz[dv_FSTz>quantile(dv_FSTz,.99,na.rm=T)]<-2
col_ddz<-col
col_ddz[dd_FSTz>quantile(dd_FSTz,.99,na.rm=T)]<-2

col_davz<-col
col_davz[dav_FSTz>quantile(dav_FSTz,.99,na.rm=T)]<-2
col_dsvz<-col
col_dsvz[dsv_FSTz>quantile(dsv_FSTz,.99,na.rm=T)]<-2

col_deltaz<-col
col_deltaz[deltaFSTz>quantile(deltaFSTz,.99,na.rm=T)]<-2
col_deltaz[deltaFSTz<0]<-"lightsteelblue2"


###plot preparation
###Functions allow plotting of whole genome with scaffolds shaded gray and white
###

scaf.name <- sapply(names(genome.pos), function(x){
  split <- strsplit(x," ")[[1]][1]
  
  return(split)
  
})

col<-1
c<-1
ax<-data.frame(scaf="scaffold_0",tick=1)
for (i in 2:length(scaf.name)){
  if(scaf.name[i] == scaf.name[i-1]){
    col[i]<-c
    
  }
  else{
    if(c==1){
      c<-2
    }
    else{
      c<-1
    }
    col[i]<-c
    tmp<-data.frame(scaf=scaf.name[i],tick=i)
    ax<-rbind(ax,tmp)
  }
}



backg<-function(s,e){
  
  plot.new()
  plot.window(xlim=c(1,length(pairwise.FST[,1])),ylim=c(s,e))
  axis(side=2)
  
  for (i in seq(from=1,to = length(ax[,1]),by=2)){
  rect(ax[i,2],-14,(ax[i+1,2]-1),10000,col=rgb(0.5,0.5,0.5,alpha=0.3),border=NA) #alpha in rgb makes color transparent}
  }
}
#z-transformed plotting

layout(matrix(c(1,2,3,4), 4, 1, byrow = TRUE),heights = c(0.5,0.5,0.5,1))
par(omi=c(0.5,0.8,0.1,0.1))
par(mai=c(0.1,0,0,0))

backg(0,10) #0,10 are ylim coordinates needed
points(1:length(pairwise.FST[,1]), dav_FSTz,col=col_davz,cex=0.5,pch=16)
backg(0,10)
points(1:length(pairwise.FST[,1]), dsv_FSTz,col=col_dsvz,cex=0.5,pch=16)
backg(0,10)
points(1:length(pairwise.FST[,1]), dd_FSTz,col=col_ddz,cex=0.5,pch=16)
#plot(1:length(pairwise.FST[,1]), dv_FSTz,ylim=c(0,10),col=col_dvz,cex=0.5,xaxt='n',xlab="",pch=16,frame.plot = F)
backg(-10,10)
points(1:length(pairwise.FST[,1]), deltaFSTz,col=col_deltaz,cex=0.5,pch=16)


