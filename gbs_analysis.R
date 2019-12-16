##Custom R script for PCA analysis of GBS data from Brand et al. 
##Title: "The evolution of sexual signaling is linked to odorant receptor tuning in perfume-collecting orchid bees"
##


library(SNPRelate)

file="GBS_SNPs_MAF0.05_SNPrelate_format.txt"
data<-read.table(paste(file,"txt",sep="."),header=T)

codes<-read.table("pop_ids_gbs.txt",header=T)

snpgdsClose(genofile) #use if doing multiple analyses in one session
snpmat<-as.matrix(data[,c(5:length(data[1,]))])
alleles<-as.character(data[,4])
gdsfile<-paste(file,"gds",sep="1.") #adjust sep if doing multiple analyses in one session
snpgdsCreateGeno(gdsfile,genmat=snpmat,sample.id=colnames(data[,c(5:length(data[1,]))]),snp.id=data[,3],snp.chromosome = data[,1],snp.position = data[,2],snp.allele = alleles,snpfirstdim = T)
genofile<-snpgdsOpen(gdsfile)

#check number of SNPs per genotyped individual
numSNPsperInd<-colSums(read.gdsn(index.gdsn(genofile, "genotype"))<3)

#select what sample to use for PCA analysis
sample<-unlist(codes[codes$colcode==1,1]) #only dilemma
sample<-codes$id #all

#PCA
pca<-snpgdsPCA(genofile,autosome.only = F,remove.monosnp = T,maf = 0.05,missing.rate = 0.5,sample.id = sample)
pc.percent<-pca$varprop*100
head(round(pc.percent,2))

sample.id<-read.gdsn(index.gdsn(genofile,"sample.id"))

pch<-codes$popcode
col<-codes$colcode

col[col==1]<-"blue";col[col==2]<-"green4";col[col==3]<-"red2";col[col==4]<-"purple"
#plot
pcalab1<-paste("PC1: ",head(round(pc.percent,2))[1],"%",sep = "")
pcalab2<-paste("PC2: ",head(round(pc.percent,2))[2],"%",sep = "")
pcalab3<-paste("PC3: ",head(round(pc.percent,2))[3],"%",sep = "")
pcalab4<-paste("PC4: ",head(round(pc.percent,2))[4],"%",sep = "")

par(mfrow=c(1,1))

plot(pca$eigenvect[,1],pca$eigenvect[,2],col=col,pch=pch,xlab=pcalab1,ylab=pcalab2,cex=1.5, main="E. dilemma PC1 PC2",cex.axis=1.5,cex.lab = 1.5,lwd = 2)
#PCA_MER_pop8_scale.pdf
points(pca$eigenvect[,c(1,2)],cex=(numSNPsperInd/max(numSNPsperInd))*5)
text(pca$eigenvect[,c(1,2)],labels = pca$sample.id,cex=0.6,pos = 1)

popnames=c("CR","TAP","COR","SAY","HON","MOR","OAX","GUA","LOS","PAR","MER","CAM","CHE","ZOH","FL")
legend("topright",legend=popnames,col="black",pch=c(1:15),cex=0.7,y.intersp = 0.7,bty = "n"  ,ncol = 3)

