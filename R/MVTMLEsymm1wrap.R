MVTMLEsymm1 <- function(X, nu = 1, eps = 1e-06, maxiter = 100)
    {
    .Call( "cMVTMLEsymm1", X, nu, eps, maxiter, PACKAGE = "fastM")
    }


