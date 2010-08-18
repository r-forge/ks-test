\name{cvm.test}
\alias{cvm.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Discrete Cramer von Mises Goodness-of-Fit Test
}
\description{
Computes the test statistics for doing one sample Cramer von Mises Goodness-of-Fit Test
and calculates asymptotic p-values. Package tries to follow the same guidelines as ks.test in the R-base.
}
\usage{
cvm.test(x, y, type = c("W2", "U2", "A2"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a numerical vecotr of data values.
 
}
  \item{y}{
   an ecdf or step-function for specifying the null model.
}
  \item{type}{
 Which variant of the Cramer von Mises test to use. W2 is the default and most common method, U2 is for cyclical data,
 and A2 is the Anderson-Darling alternative. For details see references.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%% \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Taylor B. Arnold and John W. Emerson
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x,y, type=c("W2", "U2", "A2")) {
  cvm.pval <- function(STAT, lambda) {
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
    return(0.5 + integrate(fun, 0, Inf, subdivisions=1e6)$value/pi)
  }


  cvm.stat <- function(x,y, type=c("W2", "U2", "A2")) {
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

    # Test Statistics:
    #W2 <- sum(Z^2*t)/N
    #U2 <- sum((Z-Zbar)^2*t)/N
    #A2 <- sum((Z^2*t/(H*(1-H) ))[-length(I)])/N

    S0 <- diag(p) - p \%*\% t(p)
    A <- matrix(1, length(p), length(p))
    A <- apply(row(A) >= col(A),2, as.numeric)
    E <- diag(t)
    One <- rep(1, nrow(E))
    K <- diag(0, length(H))
    diag(K)[-length(H)] <- 1/(H[-length(H)]*(1-H[-length(H)]))
    Sy <- A \%*\% S0 \%*\% t(A)
    M <- switch(type, W2 = E, U2 = (diag(1, nrow(E))-E\%*\%One\%*\%t(One))\%*\%E\%*\%(diag(1, nrow(E))-One\%*\%t(One)\%*\%E) , A2 = E\%*\%K)
    lambda <- eigen(M\%*\%Sy)$values

    STAT <- switch(type, W2 = sum(Z^2*t)/N, U2 = sum((Z-Zbar)^2*t)/N,
                        A2 = sum((Z^2*t/(H*(1-H) ))[-length(I)])/N)

    return(c(STAT, lambda))
  }


  type <- match.arg(type)
  DNAME <- deparse(substitute(x))
  if(length(setdiff(x, knots(y))) != 0) stop("Data are incompatable with null distribution")
  tempout <- cvm.stat(x,y,type=type)
  STAT <- tempout[1]
  lambda <- tempout[2:length(tempout)]
  PVAL <- cvm.pval(STAT, lambda)
  METHOD <- paste("Cramer von Mises -", type)
  names(STAT) <- as.character(type)
  RVAL <- list(statistic = STAT, p.value = PVAL, alternative = "Two.sided",
              method = METHOD, data.name=DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
