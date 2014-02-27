\name{DUEMBGENshape}
\alias{DUEMBGENshape}

\title{
Duembgen's Shape Matrix
}
\description{
Iterative algorithm to estimate Duembgen's shape matrix using a partial Newton-Raphson approach. 
}
\usage{
DUEMBGENshape(X, nmax = 500, eps = 1e-06, maxiter = 100, perm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{numeric data matrix or dataframe. Missing values are not allowed.}
  \item{nmax}{integer, if the sample size n (number of rows of \code{X}) is smaller than \code{nmax}, then all n(n-1)/2 pairwise differences will be computed
and used in the algorithm. If n is larger, then the algorithm avoids storing all the pairwise differences and is more memory efficient.}
  \item{eps}{convergence tolerance.}
  \item{maxiter}{maximum number of iterations.}
  \item{perm}{logical. If TRUE the rows of the \code{X} will be randomly permuted before starting the computations. See details.}
}
\details{
The estimate is based on the new fast algorithm described in Duembgen et al. (2014).
Note that Duembgen's shape matrix is standardized such that it has determinant 1.

The function does not check if there are several identical observations. In that case the function will fail.

To get a good initial value for the algorithm, the estimator is first computed based on the pairwise differences of 
successive observations. Therefore the order of the rows of \code{X} is supposed to be random. If this is not the case, the data
should be first permuted using the argument \code{perm}.

In case \code{maxiter} is reached before convergence, the estimate at that iteration is returned and a warning given.
}
\value{
  A list containing:
  \item{Sigma}{Estimated shape matrix.}
  \item{iter}{Number of iterations of the algorithm.}
}
\references{
\cite{Duembgen, L. (1998), On Tyler's M-functional of scatter in high dimension, \emph{Annals of Institute of Statistical Mathematics}, \bold{50},  471--491. }

\cite{Duembgen, L., Nordhausen, K. and Schumacher, H. (2014),  New algorithms for M-estimation of multivariate location and scatter, preprint available as	arXiv:1312.6489.}
}
\author{
Lutz Duembgen and Klaus Nordhausen
}

\seealso{
\code{\link[ICSNP]{tyler.shape}}
}
\examples{
DUEMBGENshape(longley)
DUEMBGENshape(longley, nmax=10)
# compare to
# library(ICSNP)
# duembgen.shape(longley)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }

