\name{TYLERshape}
\alias{TYLERshape}

\title{Tyler's Shape Matrix}

\description{Iterative algorithm to estimate Tyler's shape matrix using a partial Newton-Raphson approach.}

\usage{TYLERshape(X, location = TRUE, eps = 1e-06, maxiter = 100)}


\arguments{
  \item{X}{numeric data matrix or dataframe. Missing values are not allowed.}
  \item{location}{logical or numeric. If FALSE, it is assumed that the scatter should be computed wrt to the origin. If TRUE the location will be estimated as the mean vector and if it is a numeric
  vector it will be computed wrt to the given vector.}
  \item{eps}{convergence tolerance, which means that the algorithm stops when the Frobenius norm of the gradient is smaller than eps.}
  \item{maxiter}{maximum number of iterations.}
}

\details{
The estimate is based on the new fast algorithm described in Duembgen et al. (2016).
Note that Tyler's shape matrix is standardized such that it has determinant 1.

The function does not check if there are observations equal to the mean (if \code{location=TRUE}), to the provided location vector or to the origin (if \code{location=FALSE}).
In these cases the function will fail.

In case \code{maxiter} is reached before convergence, the estimate at that iteration is returned and a warning is given.
}
\value{
  A list containing:
  \item{mu}{Estimated location if \code{location=TRUE}, otherwise the user specified location.}
  \item{Sigma}{Estimated shape matrix.}
  \item{iter}{Number of iterations of the algorithm.}
}
\references{
\cite{Tyler, D.E. (1987), A distribution-free M-estimator of scatter, \emph{Annals of Statistics}, \bold{15},  234--251. } 

\cite{Duembgen, L., Nordhausen, K. and Schuhmacher, H. (2016),  New algorithms for M-estimation of multivariate location and scatter, \emph{Journal of Multivariate Analysis}, \bold{144},  200--217. \doi{10.1016/j.jmva.2015.11.009}}
}
\author{
Lutz Duembgen and Klaus Nordhausen
}

\seealso{
\code{\link[ICSNP]{tyler.shape}}
}
\examples{
TYLERshape(longley)
# compare to
# library(ICSNP)
# tyler.shape(longley)

TYLERshape(longley, location=FALSE)
# compare to
# library(ICSNP)
# tyler.shape(longley, location=0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
