\name{MVTMLE}
\alias{MVTMLE}

\title{
M-estimator of Location and Scatter Using Weights Coming From the Multivariate t-distribution
}
\description{
The algorithm of this function is based on a partial Newton approach and should be faster than the traditional fixed-point algorithm.
If the data follows a multivariate t-distribution with the correctly specified degrees of freedom this function gives the maximum likelihood estimate
of location and scatter.
}
\usage{
MVTMLE(X, nu = 1, location = TRUE, eps = 1e-06, maxiter = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{numeric data matrix or dataframe. Missing values are not allowed.}
  \item{nu}{assumed degrees of freedom of the t-distribution. Default is '1' which corresponds to the Cauchy distribution.}
  \item{location}{logical or numeric. If FALSE, it is assumed that the scatter should be computed wrt to the origin. If TRUE the location will be estimated and if it is a numeric
  vector it will be computed wrt to this vector.}
  \item{eps}{convergence tolerance, which means that the algorithm stops when the Frobenius norm of the gradient is smaller than eps.}
  \item{maxiter}{maximum number of iterations.}
}
\details{
The assumed degree of freedom nu must be at least 1 when the location and scatter should be estimated. If only the scatter is to be estimated, then it needs to be larger than zero only. 

In case \code{maxiter} is reached before convergence, the estimate at that iteration is returned and a warning is given.
}
\value{
 A list containing:
  \item{mu}{Estimated location if \code{location=TRUE}, otherwise the user specified location.}
  \item{Sigma}{Estimated scatter matrix.}
  \item{iter}{Number of iterations of the algorithm.}
}
\references{
\cite{Kent, J.T., Tyler, D.E. and Vardi, Y. (1994),  A curious likelihood identity for the multivariate t-distribution, \emph{Communications in Statistics, Theory and Methods}, \bold{23},  441--453.}

\cite{Duembgen, L., Nordhausen, K. and Schuhmacher, H. (2016),  New algorithms for M-estimation of multivariate location and scatter, \emph{Journal of Multivariate Analysis}, \bold{144},  200--217. \doi{10.1016/j.jmva.2015.11.009}}
}
\author{
Lutz Duembgen and Klaus Nordhausen
}



\seealso{
\code{\link[MASS]{cov.trob}}, \code{\link[ICS]{tM}}, \code{\link[fastM]{MVTMLEsymm}}
}
\examples{
MVTMLE(longley)
# compare to
# library(ICS)
# tM(longley)
# library(MASS)
# cov.trob(longley, nu=1, tol = 1e-06, maxit = 100)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
