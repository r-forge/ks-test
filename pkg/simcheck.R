
setwd("/home/jay/Desktop/SVN/ks-test/pkg")

library(ks.test)

# Default simulation:
#   Sample size 10
#   Null distribution uniform on 1:5
#   10,000 replications

thesim <- function(NN=10, vals = 1:5, N=10000) {

  res <- matrix(0, N, 3)
  f0 <- ecdf(vals)

  for (i in 1:N) {
    x <- sample(vals, NN, replace=TRUE)
    conover <- ks.test(x, f0)
    ks.0 <- ks.test(x, f0, exact=FALSE)
    res[i,1] <- conover$stat
    res[i,2] <- conover$p.value
    res[i,3] <- ks.0$p.value
  }

  res[,1] <- round(res[,1], 5)
  D <- sort(unique(res[,1]))
  cbind(D,
        y <- sapply(D, function(d, x) return(sum(x>=d)/length(x)),
                    x=res[,1]))

  return(list(D=D, y=y, res=res))
}

doplots <- function(ans) {

  print(lapply(ans, function(x) apply(x$res, 2, summary)))

  par(mfrow=c(2,2))
  for (i in 1:4) {

    D <- ans[[i]]$D
    y <- ans[[i]]$y
    res <- ans[[i]]$res

    plot(D, y, pch=19, cex=0.5, xlab="Statistic d", ylab="P(D>=d)")
    points(res[,1], res[,2], col="red")
    points(res[,1], res[,3], col="green")
    legend(max(res[,1]), 1, pch=c(19,1,1),
           col=c("black", "red", "green"), yjust=1, xjust=1, cex=0.6,
           legend=c("Simulated", "Conover (exact)", "Asymptotic (default)"))

  }

}

###############################################################################

system.time({
  ans.lapply <- lapply(c(35,45,55,65), thesim, vals=1:5, N=1000)
})

doplots(ans.lapply)

###############################################################################


library(foreach)
library(doMC)
registerDoMC(4)

library(doSNOW)
cl <- makeCluster(4)
registerDoSNOW(cl)


system.time({

  ans.1.5 <- foreach(i=c(30,35,40,45), .packages="ks.test") %dopar% {
    ans <- thesim(NN=i, vals=1:5, N=1000)
    ans
  }

})

doplots(ans.1.5)




library(doSNOW)
cl <- makeCluster(4)
registerDoSNOW(cl)

b <- 100
ans <- foreach(i=1:10, .combine=rbind) %dopar% {

  ans <- b + i^2
  c(ans, b, i)

}


ans <- foreach(i=1:10, .combine='+') %dopar% {

  i^2

}

stopCluster(cl)

##################################################################################

mysim <- function(n, NSIM, theta) {
  # The following is pretty neat, avoiding the need for looping
  # over the simulations.  Not always possible, but certainly is
  # in this case:
  x <- matrix(runif(n*NSIM, 0, theta), NSIM, n)
  tml <- apply(x, 1, max)
  talt <- tml * (n+1) / n
  diffsq.ml <- (tml-theta)^2
  diffsq.alt <- (talt-theta)^2
  mse.ml <- mean(diffsq.ml)
  mse.alt <- mean(diffsq.alt)
  diffsq <-  diffsq.ml - diffsq.alt
  ans <- t.test(diffsq)
  c(mse.ml, mse.alt, ans$stat, ans$p.value)
}

# An example:
mysim(10, 100, 3)

library(doMC)
registerDoMC(2)

system.time({

n <- 10
NSIM <- 10000
thetas <- seq(0.1, 3, by=0.1)
ans <- matrix(0, length(thetas), 4)

# In sequence (pretty damn good, 2.7 seconds total):
#for (i in 1:length(thetas)) {
#  ans[i,] <- mysim(n, NSIM, thetas[i])
#}
# In parallel over the range of values of theta (about 1.5 seconds):
ans <- foreach (theta=thetas, .combine=rbind) %dopar% {
  mysim(n, NSIM, theta)
}

plot(thetas, ans[,1], type="l", xlab="theta", ylab="MSE")
lines(thetas, ans[,2], col="green")
legend(0, 0.1, legend=c("MLE", "Alternative"), col=c(1,3), lwd=1)

})

# Could also explore a range of values of the sample sizes...
# Different issue: how to properly do the comparison?
# Answer: should be the paired comparison.






