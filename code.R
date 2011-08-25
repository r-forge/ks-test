N <- seq(10, 500, by = 20)
P <- seq(10,200, by=10)

N <-    seq(60, 85, by = 2)
P <- seq(10,30, by=1)
P <- rep(25, 15)

Q <- matrix(NA, nrow=length(P), ncol=length(N))
for(i in 1:length(N)) {
  for(j in 1:length(P)) {
    n <- N[i]
    p <- P[j]
    y <- ecdf(1:p)
    x <- sample(1:p, n, replace=TRUE)
    Q[j,i] <- ks.test(x,y, alternative="less")$p.value
  }
}
rownames(Q) <- P
colnames(Q) <- N

plot(1, 0, xlim=c(0,500), ylim=c(-650,650), col="white", xlab="Sample Size", ylab="log(|P-value|) * sign(P-value)")

for(i in 1:nrow(Q)) {
  lines(N, log(abs(Q[i,]))*sign(Q[i,]))
}

Q <- matrix(NA, nrow=40, ncol=82) 
A <- rep(NA, 40)

for(i in 1:40 ) {
   x <- sample(1:25, 82, replace=TRUE)
   y <- ecdf(1:25)
   Q[i,] <- x
   A[i] <- ks.test(x,y, alternative="less")$p.value
}




# x = data points
# y = stepfunction of cdf
cvm.stat <- function(x,y, type=c(W2, U2, A2)) {
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

  # Test Statistics:
  W2 <- sum(Z^2*t)/N
  U2 <- sum((Z-Zbar)^2*t)/N
  A2 <- sum((Z^2*t/(H*(1-H) ))[-length(I)])/N
  return(c(W2, U2, A2))
}

cvm.dist <- function(STAT, x,y, type=c(W2, U2, A2)) {

}

N <- 10000
Q <- matrix(rep(NA, N*3), nrow=N, ncol=3)
y <- ecdf(1:10)
for(i in 1:N) {
  x <- sample(1:10, 50, replace=TRUE)
  Q[i,] <- cvm.stat(x,y)
}

par(mfcol=c(1,3))
hist(Q[,1], breaks=100, main="W2")
hist(Q[,2], breaks=100, main="U2")
hist(Q[,3], breaks=100, main="A2")

# Should be 10:
sum(Q[,1] > 0.348)/nrow(Q)
sum(Q[,2] > 0.154)/nrow(Q)
sum(Q[,3] > 1.832)/nrow(Q)

## For distributions:

lambda <- c(4,1,0.5)
x <- 0.4
theta <- function(u) return(0.5*sum(atan(lambda*u)) - 0.5*x*u )
rho <- function(u) return(exp(0.25*sum(log(1+lambda^2*u^2))  ) )
fun <- function(u) return(sin(theta(u))/u*rho(u))

integrate(fun, 0, 25)



















