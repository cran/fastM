\name{MVTMLEsymm}
\alias{MVTMLEsymm}
\alias{MVTMLE_symm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Symmetrized M-estimator of Scatter Using Weights Coming From the t-distribution
}
\description{
Based on a partial Newton-Raphson approach offers this function two ways to compute the symmetrized M-estimator of scatter. The user can choose
if all pairwise differences are choosen and stored in the memory or if the computation and storage of this large matrix 
is to be avoided.
}
\usage{
MVTMLEsymm(X, nu = 1, nmax = 500, eps = 1e-06, maxiter = 100, perm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X}{numeric data matrix or dataframe with more rows than columns. Missing values are not allowed.}
\item{nu}{assumed degrees of freedom of the t-distribution, must be larger than 0. Default is '1'.}
\item{nmax}{integer, if the sample size n (number of rows of \code{X}) is smaller than \code{nmax}, then all n(n-1)/2 pairwise differences will be computed
and used in the algorithm. If n is larger, then the algorithm avoids storing all the pairwise differences and is more memory efficient.}
\item{eps}{convergence tolerance, which means that the algorithm stops when the Frobenius norm of the gradient is smaller than eps.}
\item{maxiter}{maximum number of iterations.}
\item{perm}{logical. If TRUE the rows of \code{X} will be randomly permuted before starting the computations. See details.}
}
\details{
To get a good initial value for the algorithm, the estimator is first computed based on the pairwise differences of 
successive observations. Therefore the order of the rows of \code{X} is supposed to be random. If this is not the case, the data
should be first permuted using the argument \code{perm}.

In case \code{maxiter} is reached before convergence, the estimate at that iteration is returned and a warning is given.
}
\value{
 A list containing:
  \item{Sigma}{Estimated scatter matrix.}
  \item{iter}{Number of iterations of the algorithm.}
}
\references{
\cite{Duembgen, L., Nordhausen, K. and Schuhmacher, H. (2014),  New algorithms for M-estimation of multivariate location and scatter, preprint available as arXiv:1312.6489.}
}
\author{
Lutz Duembgen and Klaus Nordhausen
}


\seealso{
 \code{\link{MVTMLE}}
}
\examples{
MVTMLEsymm(longley)
MVTMLEsymm(longley, nmax=10)
}
\keyword{multivariate }

