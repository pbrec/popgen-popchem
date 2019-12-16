##Custom R script to plot Or41 functional data from Brand et al. 
##Title: "The evolution of sexual signaling is linked to odorant receptor tuning in perfume-collecting orchid bees"
##

data<-read.table("data.txt",header = T) #read in data provided in Source data file for Figure 4 of the paper

means<-aggregate(data$delta_spike,by=list(data$receptor,data$scent),FUN=mean)
sds<-aggregate(data$delta_spike,by=list(data$receptor,data$scent),FUN=sd)
n<-aggregate(data$delta_spike,by=list(data$receptor,data$scent),FUN=length)
se <- function(x) sqrt(var(x)/length(x))
ses<-aggregate(data$delta_spike,by=list(data$receptor,data$scent),FUN=se)

select<-c("edil","evir","etri","efla","eimp","female","male","beeswax","nestwax","dentalwax","paraffinwax","eugenol","terpinen-4-ol","1.4-dimethoxybenzene","methylsalicilate","benzylbenzoate","hndb","L97","1-octadecanol","delta-tetradecalactone","dodecyl-acetate","dodecanoic_acid-3","myristic_acid-3","palmitic_acid-3","stearic_acid-3","etoh","hexane","oil","cVA")


dil<-match(select, means$Group.2)
vir<-dil+1

both<-c(sapply(seq_along(dil), function(i) append(dil[i], vir[i], i)))

#vertical barplot
par(mar=c(6.1, 4.1, 4.1, 2.1))
names<-gsub("Heat", "",means[,2])
barCenters<-barplot(means[both,3],ylim = c(-8,70),beside = T,ylab="Corrected Response (Spikes per sec)",col=c("royalblue","green4"))
arrows(barCenters,means[both,3]-ses[both,3],barCenters,means[both,3]+ses[both,3],code=3,angle=90,length=0.05)

#add axis
ticks<-vector()
for( i in seq(from=1,to=length(both)-1,by=2)){
  ticks<-rbind(ticks,(barCenters[i,1]+barCenters[i+1,1])/2)
}

axis(side = 1,at =  ticks,labels = names[dil],las=2,tick = F)

#add n to bars
text(barCenters,means[both,3]+ses[both,3]+2,labels = n[both,3])

######dose-response-curves

#HNDB

dr<-read.table("hndb_dilutions.txt",header=T) #read hndb dilution responses from Source Data file
layout(matrix(c(1)))
par(mar=c(6.1,4.1,4.1,2.1))

means.dr<-aggregate(dr$delta_spike,by=list(dr$receptor,dr$scent),FUN=mean)
sds.dr<-aggregate(dr$delta_spike,by=list(dr$receptor,dr$scent),FUN=sd)
n.dr<-aggregate(dr$delta_spike,by=list(dr$receptor,dr$scent),FUN=length)
se<- function(x) sqrt(var(x)/length(x))
ses.dr<-aggregate(dr$delta_spike,by=list(dr$receptor,dr$scent),FUN=se)

select.dr<-c("hndb-6","hndb-5","hndb-4","hndb-3","hndb-2","hndb-1")

dil.dr<-match(select.dr, means.dr$Group.2)
vir.dr<-dil.dr+2
orco<-dil.dr+1

plot(means.dr[dil.dr,3],col=c("royalblue"),ylim=c(-2,30),cex=1,pch=16,ylab="Corrected Response (Spikes per sec)",xlab="log[HNDB]",xaxt="n")
axis(1,at=1:6,labels=c(-6,-5,-4,-3,-2,-1))
arrows(c(1:6),means.dr[dil.dr,3]-ses.dr[dil.dr,3],lwd=2,c(1:6),means.dr[dil.dr,3]+ses.dr[dil.dr,3],code=3,angle=90,length=0.05,col=c("royalblue"))
#smooth line with loess
#lo<-loess(means.dr[dil.dr,3]~c(1:6),degree = 2)
#xl<-seq(min(means.dr[dil.dr,3]),max(means.dr[dil.dr,3]), (max(means.dr[dil.dr,3]) - min(means.dr[dil.dr,3]))/1000)
#lines(xl,predict(lo,xl),col="royalblue",lwd=2)
xspline(1:6,means.dr[dil.dr,3],shape = 1,border="royalblue",lwd=2)

points(means.dr[orco,3],col=c("black"),pch=16,cex=1)
arrows(c(1:6),means.dr[orco,3]-ses.dr[orco,3],lwd=2,c(1:6),means.dr[orco,3]+ses.dr[orco,3],code=3,angle=90,length=0.05,col=c("black"))
xspline(1:6,means.dr[orco,3],shape = 1,border="black",lwd=2)

points(means.dr[vir.dr,3],col=c("green4"),pch=16,cex=1)
arrows(c(1:6),means.dr[vir.dr,3]-ses.dr[vir.dr,3],lwd=2,c(1:6),means.dr[vir.dr,3]+ses.dr[vir.dr,3],code=3,angle=90,length=0.05,col=c("green4"))
xspline(1:6,means.dr[vir.dr,3],shape = 1,border="green4",lwd=2)



#PA

dr<-read.table("pa_dilutions.txt",header=T) #read palmitic acid dilution responses from Source Data file

layout(matrix(c(1)))
par(mar=c(6.1,4.1,4.1,2.1))

means.dr<-aggregate(dr$delta_spike,by=list(dr$receptor,dr$scent),FUN=mean)
sds.dr<-aggregate(dr$delta_spike,by=list(dr$receptor,dr$scent),FUN=sd)
n.dr<-aggregate(dr$delta_spike,by=list(dr$receptor,dr$scent),FUN=length)
se<- function(x) sqrt(var(x)/length(x))
ses.dr<-aggregate(dr$delta_spike,by=list(dr$receptor,dr$scent),FUN=se)

select.dr<-c("palmitic_acid-5dr","palmitic_acid-4dr","palmitic_acid-3dr","palmitic_acid-2dr")

dil.dr<-match(select.dr, means.dr$Group.2)
vir.dr<-dil.dr+2
orco<-dil.dr+1

plot(means.dr[dil.dr,3],col=c("royalblue"),ylim=c(-2,180),cex=1,pch=16,ylab="Corrected Response (Spikes per sec)",xlab="log[0.5M Palmitic acid]",xaxt="n")
axis(1,at=1:4,labels=c(-4,-3,-2,-1))
arrows(c(1:4),means.dr[dil.dr,3]-ses.dr[dil.dr,3],lwd=2,c(1:4),means.dr[dil.dr,3]+ses.dr[dil.dr,3],code=3,angle=90,length=0.05,col=c("royalblue"))
#smooth line with loess
#lo<-loess(means.dr[dil.dr,3]~c(1:6),degree = 2)
#xl<-seq(min(means.dr[dil.dr,3]),max(means.dr[dil.dr,3]), (max(means.dr[dil.dr,3]) - min(means.dr[dil.dr,3]))/1000)
#lines(xl,predict(lo,xl),col="royalblue",lwd=2)
xspline(1:4,means.dr[dil.dr,3],shape = 1,border="royalblue",lwd=2)

points(means.dr[vir.dr,3],col=c("green4"),pch=16,cex=1)
arrows(c(1:4),means.dr[vir.dr,3]-ses.dr[vir.dr,3],lwd=2,c(1:4),means.dr[vir.dr,3]+ses.dr[vir.dr,3],code=3,angle=90,length=0.05,col=c("green4"))
xspline(1:4,means.dr[vir.dr,3],shape = 1,border="green4",lwd=2)

points(means.dr[orco,3],col=c("black"),pch=16,cex=1)
arrows(c(1:4),means.dr[orco,3]-ses.dr[orco,3],lwd=2,c(1:4),means.dr[orco,3]+ses.dr[orco,3],code=3,angle=90,length=0.05,col=c("black"))
xspline(1:4,means.dr[orco,3],shape = 1,border="black",lwd=2)

