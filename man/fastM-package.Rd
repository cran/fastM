\name{fastM-package}
\alias{fastM-package}
\alias{fastM}
\docType{package}
\title{
Fast Computation of Multivariate M-estimators
}
\description{
The package implements the new algorithm for fast computation of M-scatter matrices using a partial Newton-Raphson procedure for several estimators. 
}
\details{
\tabular{ll}{
Package: \tab fastM\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2014-02-27\cr
License: \tab GPL (>= 2)
}
Multivariate M-estimators are usually computed using a fixed-point algorithm. As recently shown in Duembgen et al. (2014) 
a partial Newton-Raphson procedure applied to the second order Taylor expansion of the target function can make the computation considerably faster.
We implement this new algorithm for the multivariate M-estimator of location and scatter using weights coming from the multivariate t-distribution
(Kent et al., 1994), its symmetrized version, Tyler's shape matrix (Tyler, 1987) and Duembgen's shape matrix (Duembgen, 1998).  
For the symmetrized M-estimators we work with incomplete U-statistics to accelerate our procedures initially. 
}


\author{
Lutz Duembgen, Klaus Nordhausen, Heike Schumacher

Maintainer: Klaus Nordhausen <klaus.nordhausen@utu.fi> 

}
\references{
\cite{Duembgen, L. (1998), On Tyler's M-functional of scatter in high dimension, \emph{Annals of Institute of Statistical Mathematics}, \bold{50},  471--491. }

\cite{Duembgen, L., Nordhausen, K. and Schumacher, H. (2014),  New algorithms for M-estimation of multivariate location and scatter, preprint available as	arXiv:1312.6489.}

\cite{Kent, J.T., Tyler, D.E. and Vardi, Y. (1994),  A curious likelihood identity for the multivariate t-distribution, \emph{Communications in Statistics, Theory and Methods}, \bold{23},  441--453.}

\cite{Tyler, D.E. (1987), A distribution-free M-estimator of scatter, \emph{Annals of Statistics}, \bold{15},  234--251. } 

}

\keyword{package}
\keyword{multivariate}
