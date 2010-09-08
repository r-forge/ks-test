\name{cvm.test}
\alias{cvm.test}
\title{Discrete Cramer-von Mises Goodness-of-Fit Tests}
\description{
Computes the test statistics for doing one-sample Cramer-von Mises goodness-of-fit tests
and calculates asymptotic p-values; cvm.test() tries to follow the same guidelines as ks.test() in the R-base.
}
\usage{
cvm.test(x, y, type = c("W2", "U2", "A2"))
}
\arguments{
  \item{x}{a numerical vector of data values.}
  \item{y}{an \code{\link{ecdf}} or step-function (\code{\link{stepfun}}) for specifying the null model.}
  \item{type}{the variant of the Cramer-von Mises test; \code{"W2"} is the default and most common method, \code{"U2"} is for cyclical data, and \code{"A2"} is the Anderson-Darling alternative. For details see references.}
}
\details{
Details here.
}
\value{
An object of class \code{htest}.
%%  ~Describe the value returned
%%  If it is a LIST, use
%% \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
NEED REFERENCES.
}
\author{
Taylor B. Arnold and John W. Emerson

Maintainer: Taylor B. Arnold <taylor.arnold@yale.edu>
}
\note{
Additional notes?
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ks.test}}, \code{\link{ecdf}}, \code{\link{stepfun}}
}
\examples{

x3 <- sample(1:10, 25, replace=TRUE)

# Using ecdf() to specify a discrete distribution:
if (require(ks.test, quietly=TRUE)) {
  ks.test(x3, ecdf(1:10))
} else {
  cat(paste("See package ks.test for a revised ks.test()",
            "handling discrete distributions.\n")
}
cvm.test(x3, ecdf(1:10))

# Using step() to specify the same discrete distribution:
myfun <- stepfun(1:10, cumsum(c(0, rep(0.1, 10))))
if (require(ks.test, quietly=TRUE)) ks.test(x3, myfun)
cvm.test(x3, myfun)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
