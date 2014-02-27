#include <RcppArmadillo.h>

RcppExport SEXP MVTMLE0(SEXP x, SEXP nu, SEXP prewhitened, SEXP delta, SEXP maxiter);
RcppExport SEXP MVTMLEsymm1(SEXP x, SEXP nu, SEXP delta, SEXP maxiter);
RcppExport SEXP MVTMLEsymm2(SEXP x, SEXP nu, SEXP delta, SEXP maxiter);

RcppExport SEXP Tyler0(SEXP x, SEXP prewhitened, SEXP delta, SEXP maxiter);
RcppExport SEXP Tylersymm1(SEXP x, SEXP delta, SEXP maxiter);
RcppExport SEXP Tylersymm2(SEXP x, SEXP delta, SEXP maxiter);
