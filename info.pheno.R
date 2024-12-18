#######################################################################
# info.pheno: function to calculate information content as a function of
# the selection fraction for selective phenotyping
#######################################################################

# n = number of (unlinked) loci used for selection
# p = fraction phenotyped

info.pheno <- function(n,p,gain=FALSE)
  {
    if(gain==FALSE)
      approx(infoPhenoPoints(n),xout=p)$y
    else
      approx(infoPhenoPoints(n),xout=p)$y/p
  }

#######################################################################
#infoPhenoPoints: function to calculate information content for a
#range of selection fractions for selective phenotyping
#######################################################################

infoPhenoPoints <- function(n)
  {

    x <- 0:n
    p <- c(0,pbinom(x,n,1/2))
    dp <- diff(p)
    y <- (2*(n-x)/n)*dp
    y <- c(0,cumsum(y))
    cbind(p,y)
  }

#################################################################
# circle: function to draw a circle using a polygon approximation
#################################################################

# x = x-co-ordinate
# y = y-co-ordinate
# r = radius
# lty = line type
# npts = number of points used to draw the circle

circle <- function(x,y,r,lty,npts=1000)
  {
    theta <- seq(0,2*pi,len=npts)
    x1 <- x + r*cos(theta)
    y1 <- y + r*sin(theta)
    polygon(x1,y1,lty=lty)
  }

#######################################################################
# figure showing the information content as a function of selection fraction
#######################################################################

ps.options(pointsize=12)
# information with full genotyping
p <- seq(0,1,by=0.01)
J <- p

postscript(file="selpheno.eps", height=8,width=8,horizontal=F,paper="special")
plot(p,J,type="n",axes=F,xlim=c(0,1),ylim=c(0,1),
        xlab="",ylab="",yaxs="i",xaxs="i")

lines(p,J)
w <- infoPhenoPoints(1)
lines(w,lty="22")
w <- infoPhenoPoints(2)
lines(w,lty="44")
w <- infoPhenoPoints(3)
lines(w,lty="66")
w <- infoPhenoPoints(10)
lines(w,lty="88")
#w <- infoPhenoPoints(100)
#lines(w,lty=5)

text(0.7,0.7,"Random",pos=4)
text(0.45,0.9,"1",pos=2)
text(0.65,0.9,"2",pos=2)
text(0.74,0.92,"3",pos=2)
text(0.7,0.82,"10",pos=1)

axis(1,at=seq(0,1,by=0.1),labels=as.character(seq(0,100,by=10)))
axis(2,at=seq(0,1,by=0.1),labels=as.character(seq(0,100,by=10)))
# axis(3,at=seq(0,1,by=0.2),labels=as.character(seq(0,100,by=20)))
# axis(4,at=seq(0,1,by=0.2),labels=as.character(seq(0,100,by=20)))
mtext(side=1,line=2,
  expression(paste("Selection fraction, ",alpha, ", in percent")))
mtext(side=2,line=2,
  expression(paste("Expected percent information")))

dev.off()

#######################################################################
# figure showing the relative information content as a function of
# selection fraction
#######################################################################

postscript(file="selpheno-relative.eps", height=8,width=8,
           horizontal=F,paper="special")
p <- (0:2^10)/2^10
J <- p
plot(p,J/p,type="n",axes=F,xlim=c(0,1),ylim=c(0,2),
        xlab="",ylab="",yaxs="i",xaxs="i")
lines(p,J/p)
w <- approx(infoPhenoPoints(1),xout=p)
lines(w$x,w$y/w$x,lty="22")
w <- approx(infoPhenoPoints(2),xout=p)
lines(w$x,w$y/w$x,lty="44")
w <- approx(infoPhenoPoints(3),xout=p)
lines(w$x,w$y/w$x,lty="66")
w <- approx(infoPhenoPoints(10),xout=p)
lines(w$x,w$y/w$x,lty="88")

text(0.5,1.0,"Random",pos=1)
text(0.5,1.2,"10",pos=1)
text(0.2,1.7,"3",pos=1)
text(0.3,1.8,"2",pos=4)
text(0.55,1.8,"1",pos=4)

axis(1,at=seq(0,1,by=0.1),labels=as.character(seq(0,100,by=10)))
axis(2,at=seq(0,2,by=0.2),labels=as.character(seq(0,200,by=20)))
# axis(3,at=seq(0,1,by=0.2),labels=as.character(seq(0,100,by=20)))
# axis(4,at=seq(0,1,by=0.2),labels=as.character(seq(0,100,by=20)))
mtext(side=1,line=2,
  expression(paste("Selection fraction, ",alpha, ", in percent")))
mtext(side=2,line=2,
  expression(paste("Relative efficiency in percent")))

dev.off()


#################################################################
# draw the distances of genotypes for two loci in an F2
#################################################################
postscript(file="geno-dist.eps",height=8,width=8,horizontal=FALSE,
           paper="special")
par(mar=c(0,0,0,0))
g1 <- rep(0:2,3)
g2 <- rep(0:2,c(3,3,3))
plot(g1,g2,pty="s",ylim=c(-0.5,2.5),xlim=c(-0.5,2.5),axes=F,xlab="",ylab="")
circle(1,1,1,2)
circle(1,1,sqrt(2),2)
text(x=0,y=0,label="(0,0)",pos=2,cex=1)
text(x=0,y=1,label="(0,1)",pos=4,cex=1)
text(x=0,y=2,label="(0,2)",pos=2,cex=1)
text(x=1,y=2,label="(1,2)",pos=3,cex=1)
text(x=1,y=1,label="(1,1)",pos=1,cex=1)
text(x=1,y=0,label="(1,0)",pos=1,cex=1)
text(x=2,y=0,label="(2,0)",pos=4,cex=1)
text(x=2,y=1,label="(2,1)",pos=2,cex=1)
text(x=2,y=2,label="(2,2)",pos=4,cex=1)
dev.off()
