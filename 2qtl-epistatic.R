#########################################################################
# functions to calculate the CDF, PDF, and quantiles of a mixture of
# normal distributions
#########################################################################

########################################################################
# pmixnorm: cdf of normal mixtures
########################################################################
# x = real number
# mm = mean vector
# ss = standard deviation
# alpha = mixing proportions
# level = subtract number from CDF

pmixnorm <- function(x,mm,ss,alpha,level=0)
{
  sum( alpha * pnorm(x,mm,ss) ) - level
}

########################################################################
# dmixnorm: pdf of normal mixtures
########################################################################
# x = real number
# mm = mean vector
# ss = standard deviation
# alpha = mixing proportions

dmixnorm <- function(x,mm,ss,alpha)
{
  sum( alpha * dnorm(x,mm,ss) ) 
}

########################################################################
# qmixnorm: quantiles of normal mixtures
########################################################################
# q = fraction for which quantile is desired
# mm = mean vector
# ss = standard deviation
# alpha = mixing proportions

qmixnorm <- function(q,mm,ss,alpha)
  {
    limit <- uniroot(pmixnorm,interval=c(min(mm)-5,max(mm)+5),mm=mm,
                     ss=ss,alpha=alpha,level=q)$root
    limit
  }

#########################################################################
# Function to calculate the observed information matrix for the model
# y = b1*q1 + b2*q2 + b3*a1*q2 + e in a backcross.  Details of the
# calculations are in Maxima (see file 2qtl-epistatic.max).  The
# function integrates over the mixture distribution to obtain the
# expected information matrix components.
#########################################################################

# a = selection fraction
# b = strength of common main effect

Io2X <- function(a,b)
  {
    # mean vectors when both QTL have equal main effects
    mm <- c(-2*b,0,0,+2*b)
    # sds
    ss <- c(1,1,1,1)
    # mixing proportions in a backcross
    alpha <- c(1,1,1,1)/4

    # this is a small function to calculate the different components
    # of the observed information matrix

    # y = phenotype
    # cc = which entry of the information matrix is desired
    #    [ A B D ]
    #    [ B A D ]
    #    [ D D C ]
    
    fun <- function(y,cc)
      {
        # denominator
        den <- (1+exp(4*b*y)+2*exp(2*b*y+2*b^2))^2
        z <- y
        for( i in 1:length(y) )
          {
            x <- y[i]
            # exponentiate this to get polynomial
            w <- c( 0, 4*b*x, 8*b*x, 2*b*x+2*b^2, 4*b*x+4*b^2, 6*b*x+2*b^2 )
            # coefficients of the polynomials
            if(cc=="A")
              ccx <- c( 1, 2-4*x^2, 1, -4*x^2-8*b*x-8*b^2+4,
                       4-4*x^2, -4*x^2+8*b*x-8*b^2+4 )
            if(cc=="B")
              ccx <- c( 1, 2-4*x^2, 1, -8*b*x-8*b^2,
                       4*x^2-4, 8*b*x-8*b^2 )
            if(cc=="C")
              ccx <- c( 1, 2-16*b^2, 1, -8*x^2-16*b*x-8*b^2+4,
                       4, -8*x^2+16*b*x-8*b^2+4 )
            if(cc=="D")
              ccx <- c( -1, 8*b*x, 1, 4*x^2+12*b*x+8*b^2-2,
                       0, -4*x^2+12*b*x-8*b^2+2 )
            z[i] <- sum(ccx*exp(w))*dmixnorm(x,mm,ss,alpha)
          }
        z/den
      }

    # get lower and upper limits of integration from quantiles
    y.lo <- qmixnorm(a/2,mm,ss,alpha)
    y.hi <- -y.lo

    # do the integration and return unique entries of the matrix
    IA <- a + integrate( fun,y.lo,y.hi, cc="A" )$value
    IB <- integrate( fun,y.lo,y.hi, cc="B" )$value
    IC <- a + integrate( fun,y.lo,y.hi, cc="C" )$value
    ID <- integrate( fun,y.lo,y.hi, cc="D" )$value
    c(IA,IB,IC,ID)
  }

#######################################################################
# function calculates variance matrix from information matrix terms
#######################################################################

var2X <- function(i)
  {
    aa <- i[1]
    bb <- i[2]
    cc <- i[3]
    dd <- i[4]
    
    m <- matrix(c(aa,bb,dd, bb,aa,dd, dd,dd,cc),nrow=3)
    minv <- solve(m)
    c(minv[1,1],minv[1,2],minv[3,3],minv[1,3])
  }

#######################################################################
# function calculates the information for each parameter from the
# inverse of the variance matrix
#######################################################################
    
info2X <- function(b,h=0.001)
  {
    p <- seq(h,1-h,by=h)
    z <- matrix(nrow=length(p),ncol=4)
    for( i in 1:length(p) )
      z[i,] <- Io2X(p[i],b)
    v <- apply(z,1,var2X)
    t(1/v)
  }

########################################################################
# make the information calculations for selected strengths of main effects
########################################################################

# make calculations for a sequence 
b <- seq(0,3,by=0.1)
z <- vector(mode="list",length=length(b))
for( i in 1:length(b) )
  {
    print(i)
    z[[i]] <- info2X(b[i])
  }

# get the info for b1, and then the max and min
z1 <- sapply(z,function(x) x[,1])
min.z1 <- apply(z1,1,min)
max.z1 <- apply(z1,1,max)

# get the info for b3, and then the max and min
z3 <- sapply(z,function(x) x[,3])
min.z3 <- apply(z3,1,min)
max.z3 <- apply(z3,1,max)

# a select sequence 
ii00 <- info2X(0)
ii20 <- info2X(0.5/sqrt(2))
ii50 <- info2X(1/sqrt(2))
ii75 <- info2X(sqrt(3/2))
ii100 <- info2X(sqrt(10/2))

#######################################################################
# information for b3 when two qtl have equal main effects
#######################################################################

postscript(file="2qtl110-b3.eps",height=6,width=6,horiz=F)
h <- 0.001
alpha <- seq(h,1-h,by=h)
plot(100*alpha,100*alpha,type="n",xlab="",ylab="",ylim=c(0,100),
     xlim=c(0,100),axes=F,yaxs="i",xaxs="i")
polygon(x=c(0,100*alpha,100,rev(100*alpha),0),
        y=c(0,100*min.z3,100,rev(100*max.z3),0),
        col="gray95",border=NA)

lines(100*alpha,100*alpha)
lines(100*alpha,100*ii00[,3])
lines(100*alpha,100*ii20[,3],lty="11")
lines(100*alpha,100*ii50[,3],lty="33")
lines(100*alpha,100*ii75[,3],lty="55")
lines(100*alpha,100*ii100[,3],lty="77")

byten <- seq(0,100,by=10)
axis(1,at=byten,labels=as.character(byten))
axis(2,at=byten,labels=as.character(byten))
# axis(3,at=byten,labels=as.character(byten))
axis(4,at=byten,labels=as.character(byten))

mtext("Expected information (%)",side=2,line=2)
mtext(expression(paste("Selection fraction,",alpha," (%)")),side=1,line=2)
dev.off()

#######################################################################
# information for b1 when two qtl have equal main effects
#######################################################################

postscript(file="2qtl110-b1.eps",height=6,width=6,horiz=F)

h <- 0.001
alpha <- seq(h,1-h,by=h)
plot(100*alpha,100*alpha,type="n",xlab="",ylab="",ylim=c(0,100),
     xlim=c(0,100),axes=F,yaxs="i",xaxs="i")
polygon(x=c(0,100*alpha,100,rev(100*alpha),0),
        y=c(0,100*min.z1,100,rev(100*max.z1),0),
        col="gray95",border=NA)

lines(100*alpha,100*alpha)
lines(100*alpha,100*ii00[,1])
lines(100*alpha,100*ii20[,1],lty="11")
lines(100*alpha,100*ii50[,1],lty="33")
lines(100*alpha,100*ii75[,1],lty="55")
lines(100*alpha,100*ii100[,1],lty="77")

axis(1,at=byten,labels=as.character(byten))
axis(2,at=byten,labels=as.character(byten))
# axis(3,at=byten,labels=as.character(byten))
axis(4,at=byten,labels=as.character(byten))

mtext("Expected information (%)",side=2,line=2)
mtext(expression(paste("Selection fraction,",alpha," (%)")),side=1,line=2)
dev.off()


#########################################################################
# Function to calculate the observed information matrix for the model
# y = b1*q1 + b2*q2 + e in a backcross, where the two QTL are linked.
# We consider the case when b2 has negligible effect and we vary the
# effect of b1.  Details of the calculations are in Maxima (see file
# 2qtl-epistatic.max).  The function integrates over the mixture
# distribution to obtain the expected information matrix components.
#########################################################################

# a = selection fraction
# b = strength of common main effect

Io2x <- function(a,b)
  {
    # mean vectors when both QTL have equal main effects
    mm <- c(-b,+b)
    # sds
    ss <- c(1,1)
    # mixing proportions in a backcross
    alpha <- c(1,1)/2

    # this is a small function to calculate the different components
    # of the observed information matrix

    # y = phenotype
    # cc = which entry of the information matrix is desired
    #    [ A 0 0 ]
    #    [ 0 B C ]
    #    [ 0 C B ]
    
    fun <- function(y,cc)
      {
        # denominator
        den1 <- (1+exp(2*b*y))^2
        den2 <- (1+exp(2*b*y))
        z <- y
        for( i in 1:length(y) )
          {
            x <- y[i]
            # exponentiate this to get polynomial
            w1 <- c( 0, 2*b*x, 4*b*x )
            w2 <- c( 0, 2*b*x )
            # coefficients of the polynomials
            if(cc=="A")
              {
                ccx <- c( 1, 2-4*x^2, 1 )
                z[i] <- sum(ccx*exp(w1))*dmixnorm(x,mm,ss,alpha)
              }
            if(cc=="B")
              {
                ccx <- c( 1-(x+b)^2, 1-(x-b)^2 )
                z[i] <- sum(ccx*exp(w2))*dmixnorm(x,mm,ss,alpha)
              }
            if(cc=="C")
              {
                ccx <- c( -1+(x+b)^2, 1-(x-b)^2 )
                z[i] <- sum(ccx*exp(w2))*dmixnorm(x,mm,ss,alpha)
              }
          }

        if(cc=="A")
          z <- z/den1
        if(cc=="B")
          z <- z/den2
        if(cc=="C")
          z <- z/den2
        z
      }

    # get lower and upper limits of integration from quantiles
    y.lo <- qmixnorm(a/2,mm,ss,alpha)
    y.hi <- -y.lo

    # do the integration and return unique entries of the matrix
    IA <- a + integrate( fun,y.lo,y.hi, cc="A" )$value
    IB <- a + integrate( fun,y.lo,y.hi, cc="B" )$value
    IC <- integrate( fun,y.lo,y.hi, cc="C" )$value
    c(IA,IB,IC)
  }

#######################################################################
# function calculates variance matrix from information matrix terms
#######################################################################

var2x <- function(i)
  {
    aa <- i[1]
    bb <- i[2]
    cc <- i[3]
    
    m <- matrix(c(aa,0,0, 0,bb,cc, 0,cc,bb),nrow=3)
    minv <- solve(m)
    c(minv[1,1],minv[2,2],minv[2,3])
  }

#######################################################################
# function calculates the information for each parameter from the
# inverse of the variance matrix
#######################################################################
    
info2x <- function(b,h=0.001)
  {
    p <- seq(h,1-h,by=h)
    z <- matrix(nrow=length(p),ncol=3)
    for( i in 1:length(p) )
      z[i,] <- Io2x(p[i],b)
    v <- apply(z,1,var2x)
    t(1/v)
  }

########################################################################
# make the information calculations for selected strengths of main effects
########################################################################

# make calculations for a sequence 
b <- seq(0,3,by=0.1)
z <- vector(mode="list",length=length(b))
for( i in 1:length(b) )
  {
    print(i)
    z[[i]] <- info2x(b[i])
  }

# get the info for b1, and then the max and min
z1 <- sapply(z,function(x) x[,1])
min.z1 <- apply(z1,1,min)
max.z1 <- apply(z1,1,max)

# get the info for b2 and then the max and min
z2 <- sapply(z,function(x) x[,2])
min.z2 <- apply(z2,1,min)
max.z2 <- apply(z2,1,max)

# make info calculations for e sequence of main effect strengths
ii00 <- info2x(0)
ii20 <- info2x(0.5)
ii50 <- info2x(1)
ii75 <- info2x(sqrt(3))
ii100 <- info2x(sqrt(10))

#######################################################################
# information for b1 when only first QTL has an effect
#######################################################################

postscript(file="2qtl100-b1.eps",height=6,width=6,horiz=F)

h <- 0.001
alpha <- seq(h,1-h,by=h)
plot(100*alpha,100*alpha,type="n",xlab="",ylab="",ylim=c(0,100),
     xlim=c(0,100),axes=F,yaxs="i",xaxs="i")
polygon(x=c(0,100*alpha,100,rev(100*alpha),0),
        y=c(0,100*min.z1,100,rev(100*max.z1),0),
        col="gray95",border=NA)

lines(100*alpha,100*alpha)
lines(100*alpha,100*ii00[,1])
lines(100*alpha,100*ii20[,1],lty="11")
lines(100*alpha,100*ii50[,1],lty="33")
lines(100*alpha,100*ii75[,1],lty="55")
lines(100*alpha,100*ii100[,1],lty="77")

byten <- seq(0,100,by=10)
axis(1,at=byten,labels=as.character(byten))
axis(2,at=byten,labels=as.character(byten))
# axis(3,at=byten,labels=as.character(byten))
axis(4,at=byten,labels=as.character(byten))

mtext("Expected information (%)",side=2,line=2)
mtext(expression(paste("Selection fraction,",alpha," (%)")),side=1,line=2)
dev.off()


#######################################################################
# information for b2 when only first QTL has an effect
#######################################################################

postscript(file="2qtl100-b2.eps",height=6,width=6,horiz=F)
h <- 0.001
alpha <- seq(h,1-h,by=h)
plot(100*alpha,100*alpha,type="n",xlab="",ylab="",ylim=c(0,100),
     xlim=c(0,100),axes=F,yaxs="i",xaxs="i")
polygon(x=c(0,100*alpha,100,rev(100*alpha),0),
        y=c(0,100*min.z2,100,rev(100*max.z2),0),
        col="gray95",border=NA)

lines(100*alpha,100*alpha)
lines(100*alpha,100*ii00[,2])
lines(100*alpha,100*ii20[,2],lty="11")
lines(100*alpha,100*ii50[,2],lty="33")
lines(100*alpha,100*ii75[,2],lty="55")
lines(100*alpha,100*ii100[,2],lty="77")

byten <- seq(0,100,by=10)
axis(1,at=byten,labels=as.character(byten))
axis(2,at=byten,labels=as.character(byten))
# axis(3,at=byten,labels=as.character(byten))
axis(4,at=byten,labels=as.character(byten))

mtext("Expected information (%)",side=2,line=2)
mtext(expression(paste("Selection fraction,",alpha," (%)")),side=1,line=2)
dev.off()

