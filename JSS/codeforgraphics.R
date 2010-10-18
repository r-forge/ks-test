setwd('/home/tba3/Desktop/today/ks-test/JSS')

  cvm.pval <- function(STAT, lambda, good=TRUE) {

    x <- STAT

    theta <- function(u) {
      VAL <- 0
      for(i in 1:length(lambda)) {
         VAL <- VAL + 0.5*atan(lambda[i]*u)
      }
      return(VAL - 0.5*x*u)
    }

    rho <- function(u) {
      VAL <- 0
      for(i in 1:length(lambda)) {
        VAL <- VAL + log(1 + lambda[i]^2*u^2)     
      }
      VAL <- exp(VAL*0.25)
      return(VAL)
    }

    fun <- function(u) return(sin(theta(u))/(u*rho(u)))

    pval <- 0
    try(pval <- 0.5 + integrate(fun, 0, Inf, subdivisions=1e6)$value/pi,
        silent=TRUE)
    if(pval > 0.001 || good == FALSE) return(pval)
    if(pval <= 0.001 & good == TRUE) {
      df <- sum(lambda != 0)
      return(dchisq(STAT/max(lambda),df))
    }
  } # End cvm.pval()

lambda <- sort(runif(10))
cvm.pval(1, lambda)

I <- seq(0,35,length.out=1000)
y <- rep(NA, length(I))
for(i in 1:length(I)) {
  y[i] <- cvm.pval(I[i], lambda, good=TRUE)
}
z <- rep(NA, length(I))
for(i in 1:length(I)) {
  z[i] <- cvm.pval(I[i], lambda, good=FALSE)
}

#pdf('fig1.pdf', height=6, width=8)
plot(I[I > 18], z[I > 18], type='l', xlab='test statistic', ylab='p-value', col='salmon', lwd=2)
abline(h=0, lty=2)
lines(I[I > 18], y[I > 18], lty=2, col='blue', lwd=2)
lines(I[I > 18], z[I > 18], col='salmon', lwd=2)
#dev.off()





