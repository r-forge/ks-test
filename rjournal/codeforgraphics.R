#setwd('/home/tba3/Desktop/today/ks-test/JSS')

  cvm.stat.disc <- function(x,y, type=c("W2", "U2", "A2")) {
    type <- match.arg(type)
    I <- knots(y)
    N <- length(x)
    e <- diff(c(0,N*y(I)))
    obs <- rep(0, length(I))
    for(j in 1:length(I)) {
      obs[j] <- length(which(x == I[j]))
    }
    S <- cumsum(obs)
    T <- cumsum(e)
    H <- T/N
    p <- e/N
    t <- (p + p[c(2:length(p), 1)])/2
    Z <- S - T
    Zbar <- sum(Z*t)

    S0 <- diag(p) - p %*% t(p)
    A <- matrix(1, length(p), length(p))
    A <- apply(row(A) >= col(A),2, as.numeric)
    E <- diag(t)
    One <- rep(1, nrow(E))
    K <- diag(0, length(H))
    diag(K)[-length(H)] <- 1/(H[-length(H)]*(1-H[-length(H)]))
    Sy <- A %*% S0 %*% t(A)
    M <- switch(type, W2 = E,
                U2 = (diag(1, nrow(E)) - 
                      E%*%One%*%t(One))%*%E%*%(diag(1, nrow(E)) -
                      One%*%t(One)%*%E),
                A2 = E%*%K)
    lambda <- eigen(M%*%Sy)$values

    STAT <- switch(type, W2 = sum(Z^2*t)/N, U2 = sum((Z-Zbar)^2*t)/N,
                   A2 = sum((Z^2*t/(H*(1-H) ))[-length(I)])/N)

    return(c(STAT, lambda))
  } # End cvm.stat()


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

#lambda <- sort(runif(10))
#lambda <- sort(c(0.020837 , 0.342575 , 0.531838 , 0.604825 , 0.645453 , 0.706874 , 0.750545 , 0.776914 , 0.810034 , 0.880271))
set.seed(42)
y <- ecdf(1:3)
x <- sample(1:3, replace=TRUE)
lambda <- cvm.stat.disc(x,y,'W2')[-1]
cvm.pval(1, lambda)

I <- seq(0,5,length.out=500)
y <- rep(NA, length(I))
for(i in 1:length(I)) {
  y[i] <- cvm.pval(I[i], lambda, good=TRUE)
}
z <- rep(NA, length(I))
for(i in 1:length(I)) {
  z[i] <- cvm.pval(I[i], lambda, good=FALSE)
}
w <- rep(NA, length(I))
for(i in 1:length(I)) {
  STAT <- I[i]
  logf <- function(t) {
    ans <- -t*STAT
    ans <- ans - 0.5*sum( log(1-2*t*lambda) )
    return(ans)
  }
  est2 <- exp(nlm(logf, 1/(4*max(lambda)))$minimum)
  #est2 <- exp(optimize(logf, interval=c(0,1/(2*max(lambda))))[[1]])
  w[i] <- est2
}



#pdf('fig1.pdf', height=6, width=8)
coff <- 1.2
coff2 <- 1.3
plot(I[I > coff], z[I > coff], type='l', xlab='test statistic', ylab='p-value', col='salmon', lwd=2)
abline(h=0, lty=2)
lines(I[I > coff2], y[I > coff2], lty=2, col='blue', lwd=2)
lines(I[I > coff], w[I > coff], col='green', lwd=2)
#dev.off()





