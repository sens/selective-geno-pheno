#######################################################################
# functions to calculate the expected information as a function of
# selection fraction for selected distributions
#######################################################################

# Normal
# alpha = two-tailed selection fraction

infoNorm <- function(alpha)
  {
    x <- qnorm(alpha/2,lower.tail=FALSE)
    
    I <- alpha + 2*x*dnorm(x)
    I
  }

# Cauchy
# alpha = two-tailed selection fraction

infoCauchy <- function(alpha)
  {
    x <- qcauchy(alpha/2,lower.tail=FALSE)
    
    I <- (1/4) - atan(x)/(2*pi) + (x-x^3)/(2*pi*(x^2+1)^2)
    4*I
  }

# Logistic
# alpha = two-tailed selection fraction

infoLogis <- function(alpha)
  {
    x <- qlogis(alpha/2,lower=FALSE)
    I <- (3*exp(2*x) + 1)/(3*(exp(x)+1)^3)
    6*I
  }

# Exponential
# alpha = one-tailed selection fraction

infoExp <- function(alpha)
  {
    alpha + alpha*(log(alpha)^2)
  }

# Exponential with censoring
# alpha = one-tailed selection fraction
# beta = censoring fraction

infoExpCensor <- function(alpha,beta)
  {
    as.numeric(alpha>=beta)*(alpha - beta + alpha*(log(alpha)^2)) +
      (1-as.numeric(alpha>=beta))*(alpha*(log(beta))^2)
  }

# Weibull
# alpha = one-tailed selection fraction
# m = shape parameter

infoWeibull <- function(alpha,m=1)
  {
    x <- qweibull(alpha,shape=m,lower=FALSE)

    (m^2*(x^(2*m)+1)*exp(-x^m))/m^2
  }


# Gamma
# alpha = one-tailed selection fraction
# m = shape parameter

infoGamma <- function(alpha,m=1)
  {
    x <- qgamma(alpha,shape=m,lower=FALSE)

    A <- pgamma(x,shape=m+1,lower=FALSE)*gamma(m+1)/gamma(m)
    B <- (x-m)*(exp(-x)*x^m)/gamma(m)
    (A+B)*gamma(m)/gamma(m+1)
  }
