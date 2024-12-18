########################################################################
# plot of the expected information for single tail selective genotyping
# for the exponential distribution
########################################################################

source("info-calc.R")

postscript(file="exponential-info-censor.eps",height=6,width=8,horiz=F)
par(mar=c(3,3,2,2),pty="s")
h <- 0.001
q <- seq(h,1-h,by=h)
plot(q,infoExp(q),type="n",ylim=c(0,1),xlim=c(0,1),axes=F,
     xlab="",ylab="",yaxs="i",xaxs="i")
lines(q,infoExp(q))
lines(q,infoNorm(q),col="grey70")
lines(q,infoExpCensor(q,0.1),lty="11")
lines(q,infoExpCensor(q,0.2),lty="22")
lines(q,infoExpCensor(q,0.3),lty="33")
ptiles <- seq(0,100,by=10)
axis(1, at = ptiles/100, labels=formatC(ptiles))
axis(2, at = ptiles/100, labels=formatC(ptiles))
axis(3, at = ptiles/100, labels=formatC(ptiles))
axis(4, at = ptiles/100, labels=formatC(ptiles))
mtext("Expected information",side=2,line=2)
mtext("Selection fraction (percent)",side=1,line=2)
dev.off()

#######################################################################
# plot of the optimal selection fraction as a function of cost for
# exponential distribution
#######################################################################

library(qtlDesign)

a <- seq(h,1-h,by=h)
z0 <- vector(mode="numeric",length=length(cG))
z1 <- vector(mode="numeric",length=length(cG))
for( i in 1:length(cG) )
  {
    z <- infoExp(a)/(1+cG[i]*a)
    w <- which(z==max(z),arr.ind=T)
    z0[i] <- w[1]*h
    z1[i] <- optselection(cG[i],cross="bc")
  }

postscript(file="opt-alpha0.eps",height=8,width=10,horiz=F)
par(mar=c(3,3,3,3),pty="s")
ccc <- seq(-10,10,by=2)
plot(x=log2(cG),y=z0,ylim=c(0,1),yaxs="i",type="n",axes=F)
lines(x=log2(cG),y=z0)
lines(x=log2(cG),y=z1,lty="33")
axis(1, at = ccc, labels=c("1/1024","1/256","1/64","1/16","1/4","1",
                    "4","16","64","256","1024"))
ptiles <- seq(0,100,by=10)
axis(2, at = ptiles/100, labels=formatC(ptiles))
axis(1, at = ccc, labels=c("1/1024","1/256","1/64","1/16","1/4","1",
                    "4","16","64","256","1024"))
mtext("Genotyping cost",side=1,line=2)
mtext("Selection fraction (percent)",side=2,line=2)
dev.off()
