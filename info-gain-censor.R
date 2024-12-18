########################################################################
# information gain with censoring with an exponential phenotype
########################################################################

postscript(file="exponential-censor.eps",height=6,width=8,horiz=F)
q <- seq(0.01,0.99,by=0.01)
par(mar=c(3,3,1,1))
y <- qexp(q)
plot(q,4*y^2,type="n",ylim=c(0,60),xlim=c(0,1),axes=F,xlab="",ylab="",yaxs="i")
lines(q,4*y^2,lty="44")
lines(q,4*(y-1)^2,lty="22")
idx <- 1:85
points(q[idx],4*(y[idx]-1)^2,pch=20)
idx <- 86:99
points(q[idx],rep(4*y[86]^2,99-85),pch=20)
ptiles <- seq(0,100,by=10)
axis(1, at = ptiles/100, labels=formatC(ptiles))
axis(2, at = seq(0,60,by=10), labels=formatC(seq(0,60,by=10)))
mtext("Information gain function",side=2,line=2)
mtext("Percentile",side=1,line=2)
dev.off()


