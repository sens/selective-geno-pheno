######################################################################
# segRatio: function to calculate the segregation ratio as a function
# of distribution, and the effect size for a number of distributions
######################################################################

# n = sample size
# m = number of bins for quantiles
# dist = distribution (normal="norm", Cauchy="cauchy", Logistic="logis",
#        Gamma="gamma", Weibull="weibull", exponential ="exp"

segRatio <- function(n=1000,m=n/100,dist="norm",delta=0,...)
  {
    g <- rbinom(n,1,1/2)

    if( dist=="norm" )
      {
        y <- rnorm(n)
        y <- y+g*delta
      }
    if( dist=="cauchy" )
      {
        y <- rcauchy(n)
        y <- y+g*delta
      }
    if( dist=="logis" )
      {
        y <- rlogis(n)
        y <- y+g*delta
      }
    if( dist=="gamma" )
      {
        y <- rgamma(n,...)
        y <- y*exp(delta*g)
      }
        
    if( dist=="weibull" )
      {
        y <- rweibull(n,...)
        y <- y*exp(delta*g)
      }
    if( dist=="exp" )
      {
        y <- rexp(n)
        y <- y*exp(delta*g)
      }

    idx <- order(y)

    g1 <- g[idx]

    ptile <- floor((0:(n-1))/m)

    z <- tapply(g1,ptile,mean)
    z
  }


######################################################################
# information gain function for symmetric distributions
######################################################################

postscript(file="non-normal-symmetric.eps",height=6,width=8,horiz=F)
q <- seq(0.01,0.99,by=0.01)
par(mar=c(4,4,0,0))
y <- qnorm(q)
plot(q,y^2,type="n",ylim=c(0,6),xlim=c(0,1),axes=F,xlab="",ylab="",yaxs="i")
lines(q,y^2,lty=1)
y <- qcauchy(q)
lines(q,16*y^2/(y^4+2*y^2+1),lty="22")
ptiles <- seq(0,100,by=10)
axis(1, at = ptiles/100, labels=formatC(ptiles))
axis(2, at = 0:6, labels=formatC(0:6))
mtext("Information gain function",side=2,line=2)
mtext("Phenotype percentile",side=1,line=2)
text(0.93,5,expression("Normal"))
text(0.70,4,expression("Cauchy"))
text(0.77,1.85,expression("Logistic"))
y <- qlogis(q)
lines(q,4*(exp(y)-1)^2/(exp(y)+1)^2,lty="44")
dev.off()


###################################################################
# segregation distortion plot for symmetric distributions
###################################################################

# we get the IQR for each distribution, and then shift the mean by 0.1
# times that number (to keep things comparable)

iqr.norm <- 2*qnorm(0.75)
z.norm <- segRatio(n=10^7,m=10^5,dist="norm",delta=0.1)

iqr.cauchy <- 2*qcauchy(0.75)
z.cauchy <- segRatio(n=10^7,m=10^5,dist="cauchy",delta=0.1*iqr.cauchy/iqr.norm)

iqr.logis <- 2*qlogis(0.75)
z.logis <- segRatio(n=10^7,m=10^5,dist="logis",delta=0.1*iqr.logis/iqr.norm)

postscript(file="non-normal-symmetric1.eps",height=6,width=8,horiz=F)
par(mar=c(4,4,0,0))
plot((z.norm-0.5)^2,type="n",yaxs="i",ylim=c(0,0.004),
     ylab="Squared segregation distortion",xlab="Phenotype percentile")
lines((z.norm-0.5)^2)
lines((z.cauchy-0.5)^2,lty="22")
lines((z.logis-0.5)^2,lty="44")
dev.off()


######################################################################
# information gain function for survival distributions
######################################################################

postscript(file="non-normal-lifetime.eps",height=6,width=8,horiz=F)
par(mar=c(4,4,0,0))
y <- qgamma(q,1)
plot(q,4*(y-1)^2,type="n",ylim=c(0,25),xlim=c(0,1),
     axes=F,xlab="",ylab="",yaxs="i")
lines(q,4*(y-1)^2,lty=1)
y <- qgamma(q,3)
lines(q,4*(y-3)^2/3,lty="22")
y <- qgamma(q,10)
lines(q,4*(y-10)^2/10,lty="44")

ptiles <- seq(0,100,by=10)
axis(1, at = ptiles/100, labels=formatC(ptiles))
axis(2, at = seq(0,25,by=5), labels=formatC(seq(0,25,by=5)))
mtext("Information gain function",side=2,line=2)
mtext("Phenotype percentile",side=1,line=2)
text(0.05,4*0.6,expression("Exponential"))
text(0.03,4*1.2,expression("Gamma(3)"))
text(0.10,4*3,expression("Gamma(10)"))
dev.off()


###################################################################
# segregation distortion plot for survival distributions
###################################################################

z.exp <- segRatio(n=10^7,m=10^5,dist="exp",delta=0.1)

z.gamma3 <- segRatio(n=10^7,m=10^5,dist="gamma",delta=0.1,shape=3)
z.gamma10 <- segRatio(n=10^7,m=10^5,dist="gamma",delta=0.1,shape=10)

z.weibull10 <- segRatio(n=10^7,m=10^5,dist="weibull",delta=0.1,shape=10)

postscript(file="non-normal-lifetime1.eps",height=6,width=8,horiz=F)
par(mar=c(4,4,0,0))
plot((z.exp-0.5)^2,type="n",yaxs="i",ylim=c(0,0.004),
          ylab="Squared segregation distortion",xlab="Phenotype percentile")
lines((z.exp-0.5)^2)
lines((z.gamma3-0.5)^2/3,lty="22")
lines((z.gamma10-0.5)^2/10,lty="44")
lines((z.weibull10-0.5)^2/100,lty="66")
dev.off()

