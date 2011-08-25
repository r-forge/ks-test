\name{cvm.test}
\alias{cvm.test}
\title{Discrete Cramer-von Mises Goodness-of-Fit Tests}
\description{
Computes the test statistics for doing one-sample Cramer-von Mises goodness-of-fit tests
and calculates asymptotic p-values.
}
\usage{
cvm.test(x, y, type = c("W2", "U2", "A2"))
}
\arguments{
  \item{x}{a numerical vector of data values.}
  \item{y}{an \code{\link{ecdf}} or step-function (\code{\link{stepfun}}) for specifying the hypothesized model.}
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
V. Choulakian, R. A. Lockhart, and M. A. Stephens (1994).
\emph{Cramer-von Mises statistics for discrete distributions}.
The Canadian Journal of Statistics, 22(1): 125-137.

H. Cramer (1928).
\emph{On the composition of elementary errors.}
Skand. Akt., 11:141-180.

M. A. Stephens (1974).
\emph{Edf statistics for goodness of fit and some comparisons}.
Journal of the American Statistical Association, 69(347): 730-737.

R. E. von Mises (1928).
\emph{Wahrscheinlichkeit, Statistik und Wahrheit.}
Julius Springer, Vienna, Austria.
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
require(dgof)

x3 <- sample(1:10, 25, replace=TRUE)

# Using ecdf() to specify a discrete distribution:
ks.test(x3, ecdf(1:10))
cvm.test(x3, ecdf(1:10))

# Using step() to specify the same discrete distribution:
myfun <- stepfun(1:10, cumsum(c(0, rep(0.1, 10))))
ks.test(x3, myfun)
cvm.test(x3, myfun)

# Usage of U2 for cyclical distributions (note U2 unchanged, but W2 not)

set.seed(1)
y <- sample(1:4, 20, replace=TRUE)
cvm.test(y, ecdf(1:4), type='W2')
cvm.test(y, ecdf(1:4), type='U2')
z <- y%%4 + 1
cvm.test(z, ecdf(1:4), type='W2')
cvm.test(z, ecdf(1:4), type = 'U2')


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
