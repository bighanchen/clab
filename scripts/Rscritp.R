***Analysis of liver combined microarrays***

------------------------------------------------------------------------------------------------------
rawsignal=read.csv(file="rawsignal.txt",sep="\t");
rawsignalmRNA=subset(rawsignal, Type=="mRNA");
plot(c(-2,20),c(0,0.35),type="n",xlab="X",ylab="Y",main="Probability Density of mRNA")
for (i in 2:65){
	 lines(density(log(rawsignalmRNA[,i],2)))}
postscript("ProbabilityDensityofmRNA.ps")
par(mfrow=c(3,1)) #Figure Split


source("http://www.bioconductor.org/biocLite.R")
biocLite("limma")
library(limma) 
normsignalmRNA=normalizeBetweenArrays(as.matrix(log(rawsignalmRNA[,2:65],2)),method="s");
boxplot(log(rawsignalmRNA[,2:65],2))
boxplot(log(rawsignaluce[,2:65],2),xlab="Samples",ylab="Expression",main="Distribution of uce")   


pvalue=rep(1,64)
for (i in 1:64){
	rst=wilcox.test(qnormsignal[,1],normsignal[,i],paired=FALSE,conf.level=0.95)
	pvalue[i]=rst$p.value}
sort(pvalue)

normsignalmRNAfilter=normsignalmRNA[,pvalue>9.593942e-19]

source("plot.R")
write.table(normsignalmRNAfilter,file="normsignalmRNAfilter.txt",sep="\t",row.name=F, quote=F)
Gene=normsignalmRNAfilter_PG$Gene
normGene=array(1:6021*64,dim=c(6021,64))
medd=median(as.matrix(temp[,1]))

j=1
for (i in Gene){
	temp=subset(normsignalmRNA_PG,Gene==i)
	for (k in 1:64){
	normGene[j,k]=median(as.matrix(temp[,k]))}
	j=j+1}
------------------------------------------------------------------------------------------------------

***Analysis of Stem Cell Combined Microarray***

------------------------------------------------------------------------------------------------------
pdf("BoxplotofRawncRNA.pdf")
boxplot(log(rawsignalncRNA[,1:8],2),xlab="Samples",ylab="Expression",main="Boxplot of Raw ncRNA")
dev.off()

pdf("BoxplotofqnormncRNA.pdf")
boxplot(qnormsignalncRNA[,1:8],xlab="Samples",ylab="Expression",main="Boxplot of ncRNA After Quantile Normalization")
dev.off()

pdf("ProbabilityDensityofRawncRNA.pdf")
plot(c(-2,20),c(0,0.35),type="n",xlab="X",ylab="Y",main="Probability Density of Raw ncRNA")
for(i in 1:8){
	lines(density(log(rawsignalncRNA[,i],2)),col=i)}
dev.off()

pdf("ProbabilityDensityofqnormncRNA.pdf")
plot(c(-2,20),c(0,0.35),type="n",xlab="X",ylab="Y",main="Probability Density of ncRNA After Quantile Normalization")
for(i in 1:8){
	lines(density(qnormsignalncRNA[,i]),col=i)}
dev.off()

rawsignalncRNAM=array(1:6597*6,dim=c(6597,6))
for(i in 1:6597){
	for(j in 1:3){ 
		rawsignalncRNAM[i,j]=log(rawsignalncRNA[i,j+1]/rawsignalncRNA[i,1],2)}
	for(j in 4:6){
		rawsignalncRNAM[i,j]=log(rawsignalncRNA[i,j+2]/rawsignalncRNA[i,5],2)}
	}

rawsignalncRNAA=array(1:6597*6,dim=c(6597,6))
for(i in 1:6597){
	for(j in 1:3){ 
		rawsignalncRNAA[i,j]=log(rawsignalncRNA[i,j+1]*rawsignalncRNA[i,1],2)/2}
	for(j in 4:6){
		rawsignalncRNAA[i,j]=log(rawsignalncRNA[i,j+2]*rawsignalncRNA[i,5],2)/2}
	}
rawsignalncRNAMA=data.frame(M1=rawsignalncRNAM[,1],M2=rawsignalncRNAM[,2],M3=rawsignalncRNAM[,3],A1=rawsignalncRNAA[,1],A2=rawsignalncRNAA[,2],A3=rawsignalncRNAA[,3])

qnormsignalncRNAM=array(1:6597*6,dim=c(6597,6))
for(i in 1:6597){
	for(j in 1:3){ 
		qnormsignalncRNAM[i,j]=log(qnormsignalncRNA[i,j+1]/qnormsignalncRNA[i,1],2)}
	for(j in 4:6){
		qnormsignalncRNAM[i,j]=log(qnormsignalncRNA[i,j+2]/qnormsignalncRNA[i,5],2)}
	}

qnormsignalncRNAA=array(1:6597*6,dim=c(6597,6))
for(i in 1:6597){
	for(j in 1:3){ 
		qnormsignalncRNAA[i,j]=log(qnormsignalncRNA[i,j+1]*qnormsignalncRNA[i,1],2)/2}
	for(j in 4:6){
		qnormsignalncRNAA[i,j]=log(qnormsignalncRNA[i,j+2]*qnormsignalncRNA[i,5],2)/2}
	}
qnormsignalncRNAMA=data.frame(M1=qnormsignalncRNAM[,1],M2=qnormsignalncRNAM[,2],M3=qnormsignalncRNAM[,3],M4=qnormsignalncRNAM[,4],M5=qnormsignalncRNAM[,5],M6=qnormsignalncRNAM[,6],A1=qnormsignalncRNAA[,1],A2=qnormsignalncRNAA[,2],A3=qnormsignalncRNAA[,3],A4=qnormsignalncRNAA[,4],A5=qnormsignalncRNAA[,5],A6=qnormsignalncRNAA[,6],Gene=rawsignalncRNA$Gene,Type=rawsignalncRNA$Type)

library(limma)
qnormsignalncRNA=normalizeBetweenArrays(as.matrix(rawsignalncRNA[,1:8]),method="q")
source("getqnormsignalncRNAMA.R")
######################################################################################################
pdf("MAplotofRawncRNA21.pdf")
plot(rawsignalncRNAMA$A1,rawsignalncRNAMA$M1,xlab="A",ylab="M",main="MAplot of Raw ncRNA of 3 day")
points(lowess(rawsignalncRNAMA$A1,rawsignalncRNAMA$M1),cex=.1,col=2)
dev.off()

pdf("MAplotofqnormncRNA21.pdf")
plot(qnormsignalncRNAMA$A1,qnormsignalncRNAMA$M1,xlab="A",ylab="M",main="MAplot of ncRNA of 3 day After Quantile Normalization")
points(lowess(qnormsignalncRNAMA$A1,qnormsignalncRNAMA$M1),cex=.1,col=2)
lines(c(0,max(qnormsignalncRNAMA$A1)),c(2,2),col=2)
lines(c(0,max(qnormsignalncRNAMA$A1)),c(-2,-2),col=2)
dev.off()

up_3days=subset(qnormsignalncRNAMA,M1>2 & A1>5,select=c(M1,Gene))
write.table(up_3days,file="up_3days.txt",sep="\t",row.names=F,quote=F)
down_3days=subset(qnormsignalncRNAMA,M1<(-2) & A1>5,select=c(M1,Gene))
write.table(down_3days,file="down_3days.txt",sep="\t",row.names=F,quote=F)
#####################################################################################################
pdf("MAplotofRawncRNA31.pdf")
plot(rawsignalncRNAMA$A2,rawsignalncRNAMA$M2,xlab="A",ylab="M",main="MAplot of Raw ncRNA of 6 days")
points(lowess(rawsignalncRNAMA$A2,rawsignalncRNAMA$M2),cex=.1,col=2)
dev.off()

pdf("MAplotofqnormncRNA31.pdf")
plot(qnormsignalncRNAMA$A2,qnormsignalncRNAMA$M2,xlab="A",ylab="M",main="MAplot of ncRNA of 6 days After Quantile Normalization")
points(lowess(qnormsignalncRNAMA$A2,qnormsignalncRNAMA$M2),cex=.1,col=2)
lines(c(0,max(qnormsignalncRNAMA$A2)),c(2,2),col=2)
lines(c(0,max(qnormsignalncRNAMA$A2)),c(-2,-2),col=2)
dev.off()

up_6days=subset(qnormsignalncRNAMA,M2>2 & A2>5,select=c(M2,Gene))
write.table(up_6days,file="up_6days.txt",sep="\t",row.names=F,quote=F)
down_6days=subset(qnormsignalncRNAMA,M2<(-2) & A2>5,select=c(M2,Gene))
write.table(down_6days,file="down_6days.txt",sep="\t",row.names=F,quote=F)
#####################################################################################################
pdf("MAplotofRawncRNA41.pdf")
plot(rawsignalncRNAMA$A3,rawsignalncRNAMA$M3,xlab="A",ylab="M",main="MAplot of Raw ncRNA of 9 day")
points(lowess(rawsignalncRNAMA$A3,rawsignalncRNAMA$M3),cex=.1,col=2)
dev.off()

pdf("MAplotofqnormncRNA41.pdf")
plot(qnormsignalncRNAMA$A3,qnormsignalncRNAMA$M3,xlab="A",ylab="M",main="MAplot of ncRNA of 9 day After Quantile Normalization")
points(lowess(qnormsignalncRNAMA$A3,qnormsignalncRNAMA$M3),cex=.1,col=2)
lines(c(0,max(qnormsignalncRNAMA$A3)),c(2,2),col=2)
lines(c(0,max(qnormsignalncRNAMA$A3)),c(-2,-2),col=2)
dev.off()

up_9days=subset(qnormsignalncRNAMA,M3>2 & A3>5,select=c(M3,Gene))
write.table(up_9days,file="up_9days.txt",sep="\t",row.names=F,quote=F)
down_9days=subset(qnormsignalncRNAMA,M3<(-2) & A3>5,select=c(M3,Gene))
write.table(down_9days,file="down_9days.txt",sep="\t",row.names=F,quote=F)
#####################################################################################################
pdf("MAplotofRawncRNA65.pdf")
plot(rawsignalncRNAMA$A4,rawsignalncRNAMA$M4,xlab="A",ylab="M",main="MAplot of Raw ncRNA of 3 days")
points(lowess(rawsignalncRNAMA$A4,rawsignalncRNAMA$M4),cex=.1,col=2)
dev.off()

pdf("MAplotofqnormncRNA65.pdf")
plot(qnormsignalncRNAMA$A4,qnormsignalncRNAMA$M4,xlab="A",ylab="M",main="MAplot of ncRNA of 3 days After Quantile Normalization")
points(lowess(qnormsignalncRNAMA$A4,qnormsignalncRNAMA$M4),cex=.1,col=2)
lines(c(0,max(qnormsignalncRNAMA$A4)),c(2,2),col=2)
lines(c(0,max(qnormsignalncRNAMA$A4)),c(-2,-2),col=2)
dev.off()

up_3days_t=subset(qnormsignalncRNAMA,M4>2 & A4>5,select=c(M4,Gene))
write.table(up_3days_t,file="up_3days_t.txt",sep="\t",row.names=F,quote=F)
down_3days_t=subset(qnormsignalncRNAMA,M4<(-2) & A4>5,select=c(M4,Gene))
write.table(down_3days_t,file="down_3days_t.txt",sep="\t",row.names=F,quote=F)
#####################################################################################################
pdf("MAplotofRawncRNA75.pdf")
plot(rawsignalncRNAMA$A5,rawsignalncRNAMA$M5,xlab="A",ylab="M",main="MAplot of Raw ncRNA of 6 days")
points(lowess(rawsignalncRNAMA$A5,rawsignalncRNAMA$M5),cex=.1,col=2)
dev.off()

pdf("MAplotofqnormncRNA75.pdf")
plot(qnormsignalncRNAMA$A5,qnormsignalncRNAMA$M5,xlab="A",ylab="M",main="MAplot of ncRNA of 6 days After Quantile Normalization")
points(lowess(qnormsignalncRNAMA$A5,qnormsignalncRNAMA$M5),cex=.1,col=2)
lines(c(0,max(qnormsignalncRNAMA$A5)),c(2,2),col=2)
lines(c(0,max(qnormsignalncRNAMA$A5)),c(-2,-2),col=2)
dev.off()

up_6days_t=subset(qnormsignalncRNAMA,M5>2 & A5>5,select=c(M5,Gene))
write.table(up_6days_t,file="up_6days_t.txt",sep="\t",row.names=F,quote=F)
down_6days_t=subset(qnormsignalncRNAMA,M5<(-2) & A5>5,select=c(M5,Gene))
write.table(down_6days_t,file="down_6days_t.txt",sep="\t",row.names=F,quote=F)
#####################################################################################################
pdf("MAplotofRawncRNA85.pdf")
plot(rawsignalncRNAMA$A6,rawsignalncRNAMA$M6,xlab="A",ylab="M",main="MAplot of Raw ncRNA of 9 days")
points(lowess(rawsignalncRNAMA$A6,rawsignalncRNAMA$M6),cex=.1,col=2)
dev.off()

pdf("MAplotofqnormncRNA85.pdf")
plot(qnormsignalncRNAMA$A6,qnormsignalncRNAMA$M6,xlab="A",ylab="M",main="MAplot of ncRNA of 9 days After Quantile Normalization")
points(lowess(qnormsignalncRNAMA$A6,qnormsignalncRNAMA$M6),cex=.1,col=2)
lines(c(0,max(qnormsignalncRNAMA$A6)),c(2,2),col=2)
lines(c(0,max(qnormsignalncRNAMA$A6)),c(-2,-2),col=2)
dev.off()

up_9days_t=subset(qnormsignalncRNAMA,M6>2 & A6>5,select=c(M6,Gene))
write.table(up_9days_t,file="up_9days_t.txt",sep="\t",row.names=F,quote=F)
down_9days_t=subset(qnormsignalncRNAMA,M6<(-2) & A6>5,select=c(M6,Gene))
write.table(down_9days_t,file="down_9days_t.txt",sep="\t",row.names=F,quote=F)
------------------------------------------------------------------------------------------------------

***Survival Analysis Using R***

------------------------------------------------------------------------------------------------------
names(data.frame) #show the colomn name (or variables) of a data.frame
rownames(data.frame) # show the row names of a data.frame

#select column from a data.frame according to a list
sample_survival=read.csv(file="samplelist_survival",sep="\t");
availablenames=rownames(sample_survival)
normsignal=read.csv(file="Expression.txt",sep="\t")
normsignalNames=names(normsignal)
F=rep(1,length(normsignalNames))>0
for (i in 1:length(normsignalNames)){
	for (j in 1:length(availablenames)){ 
		if(normsignalNames[i]==availablenames[j]) F[i]=TRUE}}
F[1]=TRUE
subset=normsignal[,F]

GSE3141=read.csv(file="GSE3141_KMcurve.txt",sep="\t")
library(splines)
library(survival)
pvalue=matrix(data=0,nr=11,nc=2)

for (i in 3:13){
	GSE3141_order=GSE3141[order(GSE3141[,i]),]
	input1=GSE3141_order[c(1:floor(nrow(GSE3141_order)/2),(floor(nrow(GSE3141_order)/2)+1):nrow(GSE3141_order)),1:2]
	input1$Group=c(rep(1,floor(nrow(GSE3141_order)/2)),rep(2,nrow(GSE3141_order)-floor(nrow(GSE3141_order)/2)))
	output1=survdiff(Surv(input1$Time, input1$Status) ~ input1$Group, rho=0)
	pvalue[i-2,1]=1-pchisq(output1$chisq,1)
	input2=GSE3141_order[c(1:floor(nrow(GSE3141_order)/4),(floor(nrow(GSE3141_order)*3/4)+1):nrow(GSE3141_order)),1:2]
	input2$Group=c(rep(1,floor(nrow(GSE3141_order)/4)),rep(2,nrow(GSE3141_order)-floor(nrow(GSE3141_order)*3/4)))
	output2=survdiff(Surv(input2$Time, input2$Status) ~ input2$Group, rho=0)
	pvalue[i-2,2]=1-pchisq(output2$chisq,1)
}
MYO9B=GPL96_norm[GPL96_norm$ID_REF %in% c("208452_x_at","214780_s_at","217297_s_at"),]
SLIT2=GPL96_norm[GPL96_norm$ID_REF %in% c("209897_s_at"),] 
TRIM2=GPL96_norm[GPL96_norm$ID_REF %in% c("202341_s_at","202342_s_at","214248_s_at","214249_at","215945_s_at"),]
TRIM3=GPL96_norm[GPL96_norm$ID_REF %in% c("204910_s_at","204911_s_at","213884_s_at","213885_at"),]
USP33=GPL96_norm[GPL96_norm$ID_REF %in% c("212513_s_at","214843_s_at","217441_at"),]
EPHA2=GPL96_norm[GPL96_norm$ID_REF %in% c("203499_at"),] 
MICAL1=GPL96_norm[GPL96_norm$ID_REF %in% c("218376_s_at"),]
MICAL2=GPL96_norm[GPL96_norm$ID_REF %in% c("206275_s_at","212472_at","212473_s_at"),]
MICAL3=GPL96_norm[GPL96_norm$ID_REF %in% c("212715_s_at"),]
ARRB1=GPL96_norm[GPL96_norm$ID_REF %in% c("218832_x_at"),]
ARRB2=GPL96_norm[GPL96_norm$ID_REF %in% c("203388_at"),]

MYO9B_mean=mean(MYO9B[,2:101])
SLIT2_mean=mean(SLIT2[,2:101])
TRIM2_mean=mean(TRIM2[,2:101])
TRIM3_mean=mean(TRIM3[,2:101])
USP33_mean=mean(USP33[,2:101])
EPHA2_mean=mean(EPHA2[,2:101])
MICAL1_mean=mean(MICAL1[,2:101])
MICAL2_mean=mean(MICAL2[,2:101])
MICAL3_mean=mean(MICAL3[,2:101])
ARRB1_mean=mean(ARRB1[,2:101])
ARRB2_mean=mean(ARRB2[,2:101])

survivalanalysis=data.frame(MYO9B_mean,SLIT2_mean,TRIM2_mean,TRIM3_mean,USP33_mean,EPHA2_mean,MICAL1_mean,MICAL2_mean,MICAL3_mean,ARRB1_mean,ARRB2_mean)
write.table(survivalanalysis,file="survivalanalysis_expression.txt",sep="\t",row.name=T, quote=F)
------------------------------------------------------------------------------------------------------

***Affymetrix CEL files process***

------------------------------------------------------------------------------------------------------
source("http://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("CLL")
data(CLLbatch)
data(disease) 
head(disease)
rownames(disease)=disease$SampleID

sampleNames(CLLbatch)=sub("\\.CEL$", "", sampleNames(CLLbatch)) #rename sampleNames

mt=match(rownames(disease),sampleNames(CLLbatch)) # find the position of first argument in second argument

vmd=data.frame(labelDescription=c("Sample ID", "Disease ststus: progressive or stable disease"))
phenoData(CLLbatch) = new("AnnotatedDataFrame",data=disease[mt,], varMetadata=vmd)
CLLbatch = CLLbatch[, !is.na(CLLbatch$Disease)] 
library(simpleaffy)
saqc=qc(CLLbatch) 
dd = dist2(log2(exprs(CLLbatch)))
badArray=match("CLL1",sampleNames(CLLbatch))
CLLB=CLLbatch[,-badArray]
CLLrma=rma(CLLB) 
source("http://bioconductor.org/biocLite.R")
biocLite("hgu95av2")
biocLite("hgu95av2.db") 
library(hgu95av2)
library(hgu95av2.db)
library(annotate) 
ll = getEG(geneNames(CLLf), "hgu95av2")

smoothScatter(log2(mms[,1]),log2(pms[,1]),xlab=expression(log[2] * "MM values"), ylab=expression(log[2] * "PM values"),asp=1)

#Find indices of each probeset
pns = probeNames(CLLB)
indices = split(seq(along=pns), pns)
indices[["189_s_at"]]


library(affy)
library(simpleaffy) 
library(vsn)
GSE31546_rawdata=ReadAffy() 
sampleNames(GSE31546_rawdata)=sub("\\.CEL$","",sampleNames(GSE31546_rawdata))
GSE31546_rawdata_vsn = vsnrma(GSE31546_rawdata)
GSE31546_rawdata_vsnf = nsFilter(GSE31546_rawdata_vsn, remove.dupEntrez=FALSE,var.filter=FALSE)$eset
write.exprs(GSE31546_rawdata_vsnf, file="GSE31546_rawdata_vsnf.txt")
------------------------------------------------------------------------------------------------------

***R plot-KM-curve***

------------------------------------------------------------------------------------------------------
GSE4573_SLIT2_25=read.csv(file="GSE4573_SLIT2_25.txt",sep="\t")
GSE4573_ARRB1_50=read.csv(file="GSE4573_ARRB1_50.txt",sep="\t")

library(survival)
GSE4573_SLIT2_surv <- Surv(GSE4573_SLIT2_25$time, GSE4573_SLIT2_25$Status)
GSE4573_SLIT2_survKMest <-survfit(GSE4573_SLIT2_surv~Group,data=GSE4573_SLIT2_25)
GSE4573_ARRB1_surv <- Surv(GSE4573_ARRB1_50$time, GSE4573_ARRB1_50$Status)
GSE4573_ARRB1_survKMest <-survfit(GSE4573_ARRB1_surv~Group,data=GSE4573_ARRB1_50)

tiff("OSofSLIT2.tif")
plot (GSE4573_SLIT2_survKMest, lty=c(1,2),col=c("blue","red"),xlab="Survival Time (Months)")
legend(x=10,y=0.2,legend=c("low expression n=32","high expression n=33"),col=c("blue","red"),lty=c(1,1))
text(x=10, y=0, "P=0.022", cex=1, adj = 0)
text(x=90, y=0.9, "OS of SLIT2", font =2, cex=2, adj = 0) 
dev.off()

tiff("OSofARRB1.tif")
plot (GSE4573_ARRB1_survKMest, lty=c(1,2),col=c("blue","red"),xlab="Survival Time (Months)")
legend(x=10,y=0.2,legend=c("low expression n=65","high expression n=65"),col=c("blue","red"),lty=c(1,1))
text(x=10, y=0, "P=0.055", cex=1, adj = 0)  
text(x=90, y=0.9, "OS of ARRB1", font =2, cex=2, adj = 0)
dev.off()
