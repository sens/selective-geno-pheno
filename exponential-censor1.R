
########################################################################
# infoRatio: function to calculate the information ratio in an
# exponential distribution with censoring
########################################################################

# alpha = selection fraction
# beta = censoring fraction
# cGeno = cost of genotyping
# cTime = cost of followup

source("info-calc.R")

infoRatio <- function(alpha,beta,cGeno,cTime)
  {
    infoExpCensor(alpha,beta)/(1+cGeno*alpha+cTime*(-log(beta)))
  }

optim(par=c(0.5,0.5),lower=c(h,h),upper=c(1-h,1-h),method="L-BFGS-B",
      fn=function(x) {-log(infoRatio(x[1],x[2],cGeno,cTime))})

h <- 0.001
a <- seq(h,1-h,by=h)
b <- seq(h,1-h,by=h)

z <- outer(a,b,FUN=infoRatio,cGeno=10,cTime=9)


cG <- 2^(seq(-10,10,0.1))
cT <- 2^(seq(-10,10,0.1))

zAlpha <- matrix(nrow=length(cG),ncol=length(cT))
zBeta <- matrix(nrow=length(cG),ncol=length(cT))

for( i in 1:length(cG) )
  {
    print(i)
    for( j in 1:length(cT) )
      {
        z <- outer(a,b,FUN=infoRatio,cGeno=cG[i],cTime=cT[j])
        w <- which(z==max(z),arr.ind=T)
        zAlpha[i,j] <- w[1,1]
        zBeta[i,j] <- w[1,2]
      }
  }
  

#####################################################################
# calculate the optimal selection fraction as a function of genotyping
# and followup cost
#####################################################################

postscript(file="opt-alpha.eps",height=8,width=10,horiz=F)
par(mar=c(3,3,3,3),pty="s")
ccc <- seq(-10,10,by=2)
image(x=log2(cG),y=log2(cT),z=zAlpha,axes=F,col=topo.colors(20))
contour(x=log2(cG),y=log2(cT),z=zAlpha,axes=F,nlevels=19,add=T)
axis(1, at = ccc, labels=c("1/1024","1/256","1/64","1/16","1/4","1",
                    "4","16","64","256","1024"))
axis(2, at = ccc, labels=c("1/1024","1/256","1/64","1/16","1/4","1",
                    "4","16","64","256","1024"))
mtext("Followup cost",side=2,line=2)
mtext("Genotyping cost",side=1,line=2)
mtext("Optimal selection fraction (percent)",side=3,line=2)
dev.off()


#####################################################################
# calculate the optimal followup fraction as a function of genotyping
# and followup cost
#####################################################################

postscript(file="opt-beta.eps",height=8,width=10,horiz=F)
par(mar=c(3,3,3,3),pty="s")
ccc <- seq(-10,10,by=2)
image(x=log2(cG),y=log2(cT),z=zBeta,col=topo.colors(20))
contour(x=log2(cG),y=log2(cT),z=zBeta,axes=F,nlevels=19,add=T)
axis(1, at = ccc, labels=c("1/1024","1/256","1/64","1/16","1/4","1",
                    "4","16","64","256","1024"))
axis(2, at = ccc, labels=c("1/1024","1/256","1/64","1/16","1/4","1",
                    "4","16","64","256","1024"))
mtext("Followup cost",side=2,line=2)
mtext("Genotyping cost",side=1,line=2)
mtext("Optimal censoring fraction (percent)",side=3,line=2)
dev.off()

