

library(ks.test)

N <- 10000
res <- matrix(0, N, 4)

vals <- 1:2
f0 <- ecdf(vals)
NN <- 40

for (i in 1:N) {
  x <- sample(vals, NN, replace=TRUE)
  conover <- ks.test(x, f0)
  ks.0 <- ks.test(x, f0, exact=FALSE)
  res[i,1] <- conover$stat
  res[i,2] <- conover$p.value
  res[i,3] <- ks.0$p.value
}

res[,1] <- round(res[,1], 3)
D <- sort(unique(res[,1]))
cbind(D,
      y <- sapply(D, function(d, x) return(sum(x>=d)/length(x)),
             x=res[,1]))

plot(D, y, pch=19, cex=0.5, xlab="Statistic d", ylab="P(D>=d)")
points(res[,1], res[,2], col="red")
points(res[,1], res[,3], col="green")
legend(0.25, 0.8, pch=c(19,1,1),
       col=c("black", "red", "green"),
       legend=c("Simulated", "Conover", "Asymptotic (default)"))






