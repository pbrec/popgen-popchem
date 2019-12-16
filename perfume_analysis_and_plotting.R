##This codes allows for the analysis of the perfume data from Brand et al.
##Title: "The evolution of sexual signaling is linked to odorant receptor tuning in perfume-collecting orchid bees"
##Based on the files peak_area_matrix_Brandetal.txt and pop_ids_gcms.txt.
##a full matrix with peak areas and retention times for individual bees is available in the same github folder 
setwd("/Users/pbrand/science/research/phd/manuscripts/popgen-popchem/Second version 032519/natureComm/revision2/code")
##########
##########

rm(list=ls())
require(vegan)
require(ecodist)

#################
################

####All external compounds (minus labial gland compounds)-------------------------------------------

###prepare data
data <- read.table("peak_area_matrix_Brandetal.txt", header=T)
d<-data[,2:length(data[1,])]
rownames(d)<-data[,1]
d<- d[,which(!apply(d,2,FUN = function(x){all(x == 0)}))] 

dr<-round(d)

dr<-dr[(length(dr[,1])-colSums(dr==0))>=3] #require a compound to be part of at least 3 individuals perfumes to get included

length(dr[,1]) #number individuals

length(dr) #number compounds

ds<-dr #duplicate matrix to allow analysis of subsets of individuals (below)

####################
###matrix subsetting - choose a combination of below subsets to perform more detailed analyses
####################

##remove peaks with less than 0.1% of total area

rsums <- rowSums(ds,na.rm=T)
drPP <-ds[,1:length(ds)]/rsums
ds[drPP<0.0005]<-0


##remove individuals with less than 10 peaks

dis<-which(rownames(ds) %in% rownames(ds[rowSums(ds!=0)<10,])) #get individuals with less than X compounds (default: 10)

# only sympatric or allopatric individuals

    popids<-read.table("pop_ids_gcms.txt",header=T)

# select allopatric individuals by deleting sympatric individuals
    delsym<-which(popids$sympatric==1)
    dis<-c(dis,delsym)

# select sympatric individuals by deleting allopatric individuals
    delallo<-which(popids$sympatric==0)
    dis<-c(dis,delallo)

#make sure to delete selected individuals only once
    dis<-sort(unique(unlist(dis)))




#create subset by deleting individuals that don't match *REQUIRED!*
ds_f<-ds[-c(dis),] 


#40 most common compounds
#this can be done before or after removing individuals with <10 perfume peaks first
#some code
  ds.discrete <- ds #before removing individuals selected in variable 'dis'
  ds.discrete <- ds_f #after removing individuals selected in variable 'dis'
  ds.discrete[ds.discrete > 0] <- 1
  csums <- colSums(ds.discrete)
  ds.c0 <- ds
  ds.c1 <-rbind(ds.c0, csums)
  ds.c2 <- ds.c1[,order(unlist(ds.c1[length(ds.c1[,1]),]), decreasing=TRUE)]
  ds.c3<-ds.c2[-c(length(ds.c1[,1])),]
  
  ds<-ds.c3[,1:40]#set second value to number of x most common compounds in the dataset to analyze

##############
##Simple Stats
##############

#Number of peaks per individual-unfiltered dataset
stats<-rowSums(dr!=0)
max(stats)
min(stats)
mean(stats)
sd(stats)
length(stats)
head(dr)


#Number of peaks per individual per species
  popids<-read.table("pop_ids_gcms.txt",header=T)
  stats<-rowSums(dr[popids$virCombined==1,]!=0) #switch to choose species to analyze: 1 = viridissima; 0 = dilemma
  
  sample(head(stats[order(stats,decreasing = T)],n = 100),30)

  max(stats)
  min(stats)
  mean(stats)
  sd(stats)
  length(stats)


#Number of peaks per individual on filtered subset
  pids<-popids[-c(dis),] # 
  length(pids$virCombined) #number individuals 
  stats<-rowSums(ds_f!=0) #all
  stats<-rowSums(ds_f[pids$virCombined==1,]!=0) #1 = vir; 0 = dil Species switch to use if desired
  max(stats)
  min(stats)
  mean(stats)
  sd(stats)
  length(stats)

#compounds of >10% prevalence in one species but <10% in the other
  popids<-read.table("pop_ids_gcms.txt",header=T)
  popids<-popids[-c(dis),]
  stats<-colSums(ds_f!=0)
  
  vstats<-colSums(ds_f[popids$virCombined==1,]!=0)
  dstats<-colSums(ds_f[popids$virCombined==0,]!=0)

  length(ds_f[1,])-length(which(vstats>=(length(ds_f[popids$virCombined==1,1])*0.1))) #?
  length(ds_f[1,])-length(which(dstats>=(length(ds_f[popids$virCombined==0,1])*0.1))) #?
  
  length(ds_f[popids$virCombined==0,1]) #how many 1: vir 0: dil
  
  # compounds present in less than 10% of individuals
  length(ds_f[1,])-length(which(stats/length(ds_f[,1])*100>=10)) # all
  length(ds_f[1,])-length(which(dstats/length(ds_f[popids$virCombined==0,1])*100>=10)) # just dilemma 
  length(ds_f[1,])-length(which(vstats/length(ds_f[popids$virCombined==1,1])*100>=10)) #just viridissima

  length(intersect(names(which(vstats/length(ds_f[popids$virCombined==1,1])*100<10)),names(which(dstats/length(ds_f[popids$virCombined==0,1])*100<10)))) #How many are shared by dil and vir?

  #get names of chemicals <10% prevalence
  overall<-names(which(stats/length(ds_f[,1])*100<10)) #for all
  dilsi<-names(which(dstats/length(ds_f[popids$virCombined==0,1])*100<10)) #for dil
  virsi<-names(which(vstats/length(ds_f[popids$virCombined==1,1])*100<10)) #for vir
  which(popids$virCombined==1)
  
  
  length(overall)
  
  length(virsi) #<10% in viridissima
  length(intersect(virsi,overall)) #those overlapping in dilemma and viridissima
  
  length(dilsi)
  length(intersect(dilsi,overall))
  
  colSums(ds_f[popids$virCombined==1,setdiff(overall,virsi)]!=0)/length(which(popids$virCombined==1))*100 #prevalence of compounds <10% in overall dataset but not vir
  colSums(ds_f[popids$virCombined==0,setdiff(overall,dilsi)]!=0)/length(which(popids$virCombined==0))*100 #same for dil
  
  
  low_v_in_d<-dstats[setdiff(virsi,overall)]/length(ds_f[popids$virCombined==0,1])*100 # which ones are in <10% viridissima individuals but >=10% in dilemma
  low_d_in_v<-vstats[setdiff(dilsi,overall)]/length(ds_f[popids$virCombined==1,1])*100 #reciprocal of above
  
  low_d_in_v[order(low_d_in_v,decreasing = T)] #order by highest in v
  low_v_in_d[order(low_v_in_d,decreasing = T)] #order by highest in d
   
  both_c<-c(low_d_in_v,low_v_in_d)
  both_c<-cbind(dstats[names(both_c)],vstats[names(both_c)])

  #plot Supplementary Figure 7
  plot(both_c[,1]/length(which(popids$virCombined==0))*100,xaxt='n',xlab="",ylab="Percent individuals with compound") #dilemma values in white
  axis(side = 1,at = 1:length(both_c[,1]),labels = rownames(both_c),las=2)
  points(both_c[,2]/length(which(popids$virCombined==1))*100,pch=16) # adds viridissima values in black
  abline(v=4.5,lty=2,col='grey')
  abline(h=50,lty=2,col=2)


###################
###nMDS plot
###################

##Calculate relative abundances

rsums <- rowSums(ds_f,na.rm=T)
drPP <-ds_f[,1:length(ds_f)]/rsums
rowSums(drPP)

length(drPP) # number of compounds

#sqrt transform - use if desired
drPP<-sqrt(drPP)

#nMDS plot
PF.bc1 <- bcdist(drPP) #compute bray-curtis distance

popids<-read.table("pop_ids_gcms.txt",header=T)
popids<-popids[-c(dis),]

translation<-data.frame(org = 1:15,tra = c(8,8,1,9,2,10,3,4,11,13,14,12,6,15,7))
popids$Pop<-translation[popids$Pop,2]
colP<-c("royalblue","green4","grey","grey")
#makes popids compatible with GBS analysis

require(rgl)
require(plot3Drgl) #makes it turn
library(plot3D)

PF.bc1 <- bcdist(drPP)
PF.nmds2 <- nmds(PF.bc1, mindim=3, maxdim=3, nits=50) #compute nmds in 3D
PF.nmin2 <- nmds.min(PF.nmds2,dims = 3)

theta<-30 #modify to get best orientation of 3D cube. Use plotgrl() function to find best orientation
phi<-20 #same
scatter3D(x = PF.nmin2[,2], y = PF.nmin2[,1], z = PF.nmin2[,3], colvar = NULL,type = 'p',pch=popids[,2],col=colP[(popids[,6]+1)],bty='b2',theta =theta, phi = phi,ticktype = "detailed",cex = 1.5)
plotrgl() #explore 3D plot in real time

## plot point sizes according to rel. abundance of hndb and L97

tmp<-drPP
hndb.s<-drPP$X50.9*10
l97.s<-drPP$X55*10
l97.s[which(l97.s!=0 & l97.s <1)]<-1
hndb.s[which(hndb.s!=0 & hndb.s <1)]<-1

scatter3D(x = PF.nmin2[,2], y = PF.nmin2[,1], z = PF.nmin2[,3], colvar = NULL,type = 'p',pch=16,col=colP[(popids[,6]+1)],bty='b2',theta =0, phi = 160,ticktype = "detailed",cex = 0.8,main = "Individuals with HNDB or L97")

scatter3D(x = PF.nmin2[,2], y = PF.nmin2[,1], z = PF.nmin2[,3], colvar = NULL,type = 'p',pch=1,col="red",bty='b2',theta =theta, phi = phi,cex = hndb.s, add = T)
scatter3D(x = PF.nmin2[,2], y = PF.nmin2[,1], z = PF.nmin2[,3], colvar = NULL,type = 'p',pch=1,col="black",bty='b2',theta =theta, phi = phi,cex = l97.s, add = T)


###########
###ANOSIM and SIMPER analyses
###########

#ANOSIM code
#Analysis of similarities (ANOSIM) provides a way to test statistically whether there is a significant difference between two or more groups of sampling units.
#The method is philosophically allied with NMDS ordination, in that it uses only the rank order of dissimilarity values.

require(vegan)
require(ecodist)
require(plyr)

setwd("/Users/pbrand/science/research/phd/manuscripts/popgen-popchem/Second version 032519/natureComm/revision2/code")

rm(list=ls())

data <- read.table("peak_area_matrix_Brandetal.txt", header=T)
d<-data[,2:length(data[1,])]
rownames(d)<-data[,1]

dr<-round(d)
dr<-dr[(length(dr[,1])-colSums(dr==0))>=3] #require a compound to be part of at least 3 individuals perfumes to get included

ds<-dr
#50 most commmon compounds##
ds.discrete <- ds
ds.discrete[ds.discrete > 0] <- 1
csums <- colSums(ds.discrete)
ds.c0 <- ds
ds.c1 <-rbind(ds.c0, csums)
ds.c2 <- ds.c1[,order(unlist(ds.c1[length(ds.c1[,1]),]), decreasing=TRUE)]
ds.c3<-ds.c2[-c(length(ds.c1[,1])),]
ds<-ds.c3[,1:50]

##remove individuals with less than 10 peaks

dis<-which(rownames(ds) %in% rownames(ds[rowSums(ds!=0)<10,])) #get individuals with less than X compounds

dis<-sort(unique(unlist(dis)))
ds_f<-ds[-c(dis),]

##Calculate relative abundances

rsums <- rowSums(ds_f,na.rm=T)
drPP <-ds_f[,1:length(ds_f)]/rsums
rowSums(drPP)
length(drPP[,1])

#create a data frame with grouping information that will be attached to the bray-curtis disimilairity matrix)

popids<-read.table("pop_ids_gcms.txt",header=T)
popids<-popids[-c(dis),]

grouping<- lapply(popids, function(x) if(is.factor(x)) factor(x) else x)
grouping<-lapply(grouping,function(x) factor(x))

summary(grouping$virCombined)

bc.dist<-bcdist(drPP)
jar.dist<-vegdist(drPP)

##ANOSIM
attach(grouping)
ano<-anosim(bc.dist,virCombined,permutations=999)
detach(grouping)
summary(ano)
plot(ano)
hist(ano$perm) #null distribution of R statistic
summary(grouping$virPerf)
ano$dis.rank

##SIMPER

sim<-with(grouping, simper(drPP, virCombined))
sim.sum<-summary(sim)
head(sim.sum$"1_0",n=18) #1_0 specifies dilemma-viridissima comparison




