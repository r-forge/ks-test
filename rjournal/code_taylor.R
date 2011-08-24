# Function to get lambdas:
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
  } # End cvm.stat.disc()


logf <- function(t, x, lambda) {
  ans <- -t*x
  ans <- ans - 0.5*sum( log(1-2*t*lambda) )
  return(ans)
}

x <- 15
lambda <- rep(1,3)
tvals <- seq(0.00001, 1/(2*max(lambda)), length.out=100)
vals <- exp(sapply(tvals, logf, x, lambda))
plot(tvals, vals)
ans1 <- min(vals[is.finite(vals)])
df <- sum(lambda != 0)
ans2 <- dchisq(x/max(lambda),df)

## Simulation I
n <- 10
N <- 1000

A <- matrix(NA, nrow=N, ncol=3)
for(i in 1:1000) {
  X <- sample(c(1,3,5),n, replace=TRUE)
  out <- cvm.stat.disc(X, ecdf(1:5))
  STAT <- out[1]
  lambda <- out[-1]

  x <- STAT
  lambda <- lambda
  tvals <- seq(0.00001, 1/(2*max(lambda)), length.out=100)
  vals <- exp(sapply(tvals, logf, x, lambda))
  #plot(tvals, vals)
  ans1 <- min(vals[is.finite(vals)])
  df <- sum(lambda != 0)
  ans2 <- dchisq(x/max(lambda),df)
  A[i,1] <- STAT
  A[i,2] <- ans1
  A[i,3] <- ans2
}

 plot(A[,3], A[,2], xlab='p-value us', ylab='p-value reviewer', pch=19, main='simulation from sampling')
 abline(0,1)


## Simulation II
n <- 10
N <- 1000
X <- c(1,1,1,2,2,4,5,6,7)
out <- cvm.stat.disc(X, ecdf(1:7))
STAT <- out[1]
lambda <- out[-1]
xvals <- seq(1.2,STAT*12, length.out=N)

A <- matrix(NA, nrow=N, ncol=3)
for(i in 1:1000) {
  x <- xvals[i]
  tvals <- seq(0.00001, 1/(2*max(lambda)), length.out=100)
  vals <- exp(sapply(tvals, logf, x, lambda))
  #plot(tvals, vals)
  ans1 <- min(vals[is.finite(vals)])
  df <- sum(lambda != 0)
  ans2 <- dchisq(x/max(lambda),df)
  A[i,1] <- x
  A[i,2] <- ans1
  A[i,3] <- ans2
}

 plot(A[,1], A[,2], ylim=c(0,0.02), type='l', col='blue', xlab='test statistic', ylab='p-value')
 lines(A[,1], A[,3], col='red')







