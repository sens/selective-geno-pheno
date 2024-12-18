##################################################################
# opt.wait: function to calculate the optimal followup time as a
# function of cost
##################################################################

opt.wait <- function(cost,prop=TRUE,interval=c(0,1000))
  {
    n <- length(cost)
    z <- rep(0,n)

    for( i in 1:n )
      {
        f <- function(x) exp(x)-1-x-1/cost[i]
        w <- uniroot(f, interval=interval)$root
        if(prop)
          z[i] <- 1-exp(-w)
        else
          z[i] <- w
      }
    z
  }


#####################################################################
# function to calculate the optimal uncensoring fraction
#####################################################################

postscript(file="opt-wait.eps",height=8,width=11,horiz=F,paper="special")
ccc <- 2^seq(-10,16,by=0.1)
w <- opt.wait(ccc)
plot(log2(ccc),w,pch=".",ylim=c(0,1),axes=FALSE,xlab="Cost of followup",
     ylab="Proportion uncensored",yaxs="i")
lines(log2(ccc),w)

axis(1,log2(c(1/1000,1/250,1/50,1/10,1/3,1,3,
            10,50,250,1000,3000,10000,50000)),
     label=c("1/1000","1/250","1/50","1/10",
       "1/3","1","3","10","50","250","1000",
       "3000","10000","50000"))
axis(2,seq(0,1,by=0.1))
dev.off()

