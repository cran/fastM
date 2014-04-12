\name{MVTMLE0r}
\alias{MVTMLE0r}
\alias{MVTMLE0r_FP}
\alias{MVTMLE0r_FP0}
\alias{MVTMLE0r_G}
\alias{MVTMLE0r_CG}
\title{
Different Algorithms for M-estimation of Scatter Using Weights Coming From the Multivariate t-distribution
}
\description{
The functions below are only for comparison purposes and are all written in R. 
Each function corresponds to a different algorithm for the scatter only problem for M-estimation
using weights coming from the multivariate t-distribution.
}
\usage{
MVTMLE0r(X, nu = 0, delta = 1e-06, prewhitened = FALSE, steps = FALSE)
MVTMLE0r_FP(X, nu = 0, delta = 1e-06, steps = FALSE)
MVTMLE0r_FP0(X, nu = 0, delta = 1e-06, steps = FALSE)
MVTMLE0r_G(X, nu = 0, delta = 1e-06, steps = FALSE)
MVTMLE0r_CG(X, nu = 0, delta = 1e-06, steps = FALSE)
}

\arguments{
  \item{X}{numeric data matrix or dataframe. Missing values are not allowed.}
  \item{nu}{assumed degrees of freedom of the t-distribution. Must be 0 or larger. Default is '0' which corresponds to Tyler's shape matrix.}
  \item{delta}{convergence tolerance, which means that the algorithms stop when the Frobenius norm of the gradient is smaller than delta.}
  \item{prewhitened}{logical. Is the data prewhitened or not.}
  \item{steps}{logial. If TRUE intermediate results are printed on the console.}
}
\details{
All functions are implemented in R and their purpose is only for demonstration of the differences of the different algorithms.
The function \code{MVTMLE0r} uses the recommended partial Newton approach as implemented also in (\code{\link{MVTMLE}} and \code{\link{TYLERshape}}).
\code{MVTMLE0r_FP} and \code{MVTMLE0r_FP0} are fixed-point algorithms where \code{MVTMLE0r_FP} iterates the fixed point equation with
'iterative whitening' of the data. The function \code{MVTMLE0r_G} uses a gradient approach and \code{MVTMLE0r_CG} a conjugate gradient approach.
Note that \code{MVTMLE0r_CG} does not check if the 'next' step is really an improvement and that all functions compute the scatter wrt to the origin. 

All functions have a hard coded maximum number of iterations of 1000. If that is reached the functions returns the final estimate, however without a warning.

For general purposes we recommend the functions \code{\link{MVTMLE}} and \code{\link{TYLERshape}}.


}
\value{
 A list containing at least:
  \item{S}{Estimated scatter matrix (or shape matrix if \code{nu=0}).}
  \item{iter}{Number of iterations of the algorithm.}
}
\references{
\cite{Duembgen, L., Nordhausen, K. and Schuhmacher, H. (2014),  New algorithms for M-estimation of multivariate location and scatter, preprint available as arXiv:1312.6489.}
}
\author{
Lutz Duembgen and Klaus Nordhausen
}


\seealso{ \code{\link{MVTMLE}}, \code{\link{TYLERshape}}}

\examples{
MVTMLE0r(longley,nu=1)
MVTMLE0r_FP(longley,nu=1)
MVTMLE0r_FP0(longley,nu=1)
MVTMLE0r_G(longley,nu=1)
MVTMLE0r_CG(longley,nu=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate}

