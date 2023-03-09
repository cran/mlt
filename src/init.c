
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>


/* .Call calls */
extern SEXP R_rowSums(SEXP, SEXP);
extern SEXP R_offrowSums(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"R_rowSums", (DL_FUNC) &R_rowSums, 2},
    {"R_offrowSums", (DL_FUNC) &R_offrowSums, 3},
    {NULL, NULL, 0}
};

void attribute_visible R_init_mlt(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

