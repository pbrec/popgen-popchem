data<-read.table('moments_results.txt',header = F,sep = "\t")

means<-aggregate(data$V5,by=list(data$V6,data$V2),FUN=mean)

npara<-data[,c(2,4)]
npara<-unique(npara)
colnames(npara)<-c("Group.2","nparams")
npl<-merge(means,npara,by="Group.2")

#npl<-npl[grepl(pattern = "[V,",x = npl$Group.1, fixed=T),]

#AIC

names(npl)<-c("model","spec","ll","npara")

npara<-npl$npara
likes<-npl$ll

aic=2*npara-2*likes
delta.aic= aic-min(aic)
delta.aic[delta.aic>100]=100

npl$aic<-aic
npl$delta.aic<-delta.aic

# weights of evidence
exps=exp(-delta.aic/2)
wts=exps/sum(exps)

npl$wts<-wts


npl=npl[order(npl$wts,npl$aic,decreasing=c(T,F)),]
npl

