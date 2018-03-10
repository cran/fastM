#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP cMVTMLE0(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cMVTMLEsymm1(SEXP, SEXP, SEXP, SEXP);
extern SEXP cMVTMLEsymm2(SEXP, SEXP, SEXP, SEXP);
extern SEXP Tyler0(SEXP, SEXP, SEXP, SEXP);
extern SEXP Tylersymm1(SEXP, SEXP, SEXP);
extern SEXP Tylersymm2(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"cMVTMLE0",     (DL_FUNC) &cMVTMLE0,     5},
    {"cMVTMLEsymm1", (DL_FUNC) &cMVTMLEsymm1, 4},
    {"cMVTMLEsymm2", (DL_FUNC) &cMVTMLEsymm2, 4},
    {"Tyler0",       (DL_FUNC) &Tyler0,       4},
    {"Tylersymm1",   (DL_FUNC) &Tylersymm1,   3},
    {"Tylersymm2",   (DL_FUNC) &Tylersymm2,   3},
    {NULL, NULL, 0}
};

void R_init_fastM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
