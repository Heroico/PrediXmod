qqunif = function(p,BH=T,CI=T,FDRthres=0.05,plot=F,title="",...)
{
    nn = length(p)
    xx =  -log10((1:nn)/(nn+1))
    dat<-cbind(sort(p),1:nn)
    q<-(nn*dat[,1])/dat[,2] # calculate q-values from p-values
    dat<-cbind(dat,q)
    if (min(dat[,3]>FDRthres)) {
        nsnps<-0
    } else {
        nsnps<-round((sum(p<=max(dat[dat[,3]<FDRthres,1]))/nn)*100,4)
    }
    if (plot) {
        plot( xx,  -sort(log10(p)),
        xlab=expression(Expected~~-log[10](italic(p))),
        ylab=expression(Observed~~-log[10](italic(p))),
        cex.lab=1,mgp=c(2,1,0),pch=20,main=title,cex.main=1,
        ... )
        if(CI)
        {
            ## create the confidence intervals
            c95 <- rep(0,nn)
            c05 <- rep(0,nn)
            ## the jth order statistic from a
            ## uniform(0,1) sample
            ## has a beta(j,n-j+1) distribution
            ## (Casella & Berger, 2002,
            ## 2nd edition, pg 230, Duxbury)
            ## this portion was posted by anonymous on
            ## http://gettinggeneticsdone.blogspot.com/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html
    
            for(i in 1:nn)
            {
             	c95[i] <- qbeta(0.95,i,nn-i+1)
                c05[i] <- qbeta(0.05,i,nn-i+1)
            }
            polygon(c(xx, rev(xx)), c(-log10(c95), rev(-log10(c05))),col = "grey", border = NA)
            #    lines(xx,-log10(c95),col='gray')
            #    lines(xx,-log10(c05),col='gray')
        }
        points(xx,  -sort(log10(p)),pch=20)
        #points(xx[dat[,3]<0.05],  -log10(dat[dat[,3]<0.05,1]),pch=20,col="pink")


        y<-max(-log10(p))
	siggenes <- sum(p<=max(dat[dat[,3]<FDRthres,1]))
	totgenes <- length(p)
        text(0,y,paste(nsnps,"% of genes (",siggenes," of ",totgenes,") have a q-value <= ",FDRthres,sep=""),pos=4)
    
        #  print(paste(nsnps,"% of genes have a q-value <= ",FDRthres,sep=""))
        abline(0,1,col='red',lwd=2)
        if(BH)
        {
            abline(-log10(0.05),1, col='black',lty=2,lwd=1.5)
            abline(-log10(0.10),1, col='black',lty=3,lwd=1.5)
            abline(-log10(0.25),1, col='black',lty=4,lwd=1.5)
            legend('bottomright', c("FDR = 0.05","FDR = 0.10","FDR = 0.25"),
            col=c('black','black','black'),lty=2:4, lwd=2, cex=1)
            abline(h=-log10(0.05/nn),col="blue",lwd=2) ## bonferroni
        }
    } else {
        return(nsnps)
    }
}

