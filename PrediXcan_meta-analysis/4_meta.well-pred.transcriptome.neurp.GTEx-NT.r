###originally run on genegate in /nas40t2/hwheeler/predict_EXP/NEURP_PrediXcan_20140804/
###performs PrediXcan on each cohort in dislist and combines the results in an inverse-variance meta-analysis
###make tissue.genelist from header of *profile, replace - and : with . 

rm(list=ls())

"%&%" = function(a,b) paste(a,b,sep="")

.libPaths("/userhome/hwheeler/R/x86_64-redhat-linux-gnu-library/2.15/")
source('/nas40t2/hwheeler/predict_EXP/qqunif.r')

tissue = 'GTEx-NT'
#tissue = 'GTEx-WB'
#tissue = 'LCL-GD-all'
#tis = 'LCL.GD'
tis = 'GTEx-NT'
qval = 0.05

genelistfile = 'qval.lt' %&% qval %&% '.' %&% tissue %&% '.genelist'
genelist = scan(genelistfile,'character')
genelist = unique(genelist)
ngen <- length(genelist)
#dislist = c('CALGB40101.G3','CALGB90401.G3','E5103.G3','GoKinD.NEURP','T13b.G3','COG0433.G3')
dislist = c('CALGB40101.G3','E5103.G3')
#dislist = c('T13b.G3','COG0433.G3')
ndis <- length(dislist)
resultsarray = array(0,c(ndis,4,ngen))
dimnames(resultsarray)[[3]] = genelist

for(i in 1:length(dislist))
{
  disease = dislist[i]
  #filename = '/nas40t2/hwheeler/predict_EXP/NEURP_PrediXcan_20140804/results/' %&% disease %&% '.dos.GENEX.' %&% tissue %&% '.hapmap2.SNPs.profile'
  filename = '/nas40t2/hwheeler/predict_EXP/NEURP_PrediXcan_20140804/results/' %&% disease %&% '.dos.GENEX.' %&% tis %&% '.hapmap2.SNPs.profile'
  tempo = read.table(filename,header=T,as.is=T,sep='\t')
  pheno <- tempo[,6]
  for(j in 1:length(genelist)){
	gene = genelist[j]
	resultsarray[i,,gene] <- coef(summary(lm(tempo[,gene]~pheno)))[2,]
  }
}

coe = coef(summary(lm(tempo[,gene]~pheno)))
dimnames(resultsarray)[[2]] = colnames(coe)
dimnames(resultsarray)[[1]] = dislist

metaarray <- array(0, c(ngen,3))
for(i in 1:length(genelist)){
	gene <- genelist[i]
	res <- resultsarray[,,gene]
	wi <- 1/(res[,2]^2)
	se <- sqrt(1/sum(wi))
	beta <- (sum(res[,1]*wi))/(sum(wi))
	z <- beta/se
	p <- 2*pnorm(-abs(z))
	meta <- cbind(gene,z,p)
	metaarray[i,] <- meta
}

metahead <- c("Gene","Z-score","P-value")
colnames(metaarray) <- metahead
finalmeta<-metaarray[order(as.numeric(metaarray[,3])),]

date <-Sys.Date()
#write.table(finalmeta,file='results/META.G3.NEURP.resid.' %&% tissue %&% '.qval.lt' %&% qval %&% '.' %&% ndis %&% 'VCRcohorts.hapmap2.lasso.' %&% date %&% '.csv',quote=F,row.names=F,sep=',')
write.table(finalmeta,file='results/META.G3.NEURP.resid.' %&% tissue %&% '.qval.lt' %&% qval %&% '.' %&% ndis %&% 'cohorts.hapmap2.lasso.' %&% date %&% '.csv',quote=F,row.names=F,sep=',')
metap <- as.numeric(finalmeta[,3])
metap<-metap[!is.nan(metap)]
#png(file='results/META.G3.NEURP.resid.' %&% tissue %&% '.qval.lt' %&% qval %&% '.' %&% ndis %&% 'VCRcohorts.hapmap2.lasso.qqplot.' %&% date %&% '.png')
png(file='results/META.G3.NEURP.resid.' %&% tissue %&% '.qval.lt' %&% qval %&% '.' %&% ndis %&% 'cohorts.hapmap2.lasso.qqplot.' %&% date %&% '.png')
#qqunif(metap,plot=T,FDRthres=0.25,title='META.G3.NEURP.resid.' %&% tissue %&% '.qval.lt' %&% qval %&% '.' %&% ndis %&% 'VCRcohorts.hapmap2.lasso')
qqunif(metap,plot=T,FDRthres=0.25,title='META.G3.NEURP.resid.' %&% tissue %&% '.qval.lt' %&% qval %&% '.' %&% ndis %&% 'cohorts.hapmap2.lasso')
dev.off()


tvals<-resultsarray[,3,]
tvals<-t(tvals)

pvals<-resultsarray[,4,]
pvals<-t(pvals)

date<-Sys.Date()

alltp <- cbind(metaarray,tvals,pvals)
#write.table(alltp,file='results/NEURP.G3.lm.resid.' %&% tissue %&% '.qval.lt' %&% qval %&% '.' %&% ndis %&% 'VCRcohorts.hapmap2.lasso.t-stats.p-values.' %&% date %&% '.csv',quote=F,row.names=F,sep=',')
write.table(alltp,file='results/NEURP.G3.lm.resid.' %&% tissue %&% '.qval.lt' %&% qval %&% '.' %&% ndis %&% 'cohorts.hapmap2.lasso.t-stats.p-values.' %&% date %&% '.csv',quote=F,row.names=F,sep=',')

for(i in 1:length(dislist)){
  disease = dislist[i]
  p <- pvals[,i]
  p<-p[!is.nan(p)]
  png(file='results/' %&% disease %&% '.' %&% tissue %&% '.qval.lt' %&% qval %&% '.hapmap2.lasso.qqplot.' %&% date %&% '.png')
  qqunif(p,plot=T,FDRthres=0.25,title=disease %&% '.' %&% tissue %&% '.qval.lt' %&% qval %&% '.hapmap2.lasso')
  dev.off()
}



