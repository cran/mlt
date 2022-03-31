
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>

/* nrow of a matrix */

int NROW
(
    SEXP x
) {

    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(XLENGTH(x));
    if (TYPEOF(a) == REALSXP)
        return(REAL(a)[0]);
    return(INTEGER(a)[0]);
}

/* ncol of a matrix */

int NCOL
(
    SEXP x
) {

    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(1);
    if (TYPEOF(a) == REALSXP)
        return(REAL(a)[1]);
    return(INTEGER(a)[1]);
}

SEXP R_rowSums
(
    SEXP x,
    SEXP coef
) {

    SEXP ret;
    double *dx, *dc, *dr;
    int p, n, idx;

    if (!isReal(x))
        error("x is not real");
    if (!isReal(coef))
        error("coef is not real");    
    n = NROW(x);
    p = NCOL(x);
    if (n != NROW(coef))
        error("number of rows do not match");
    if (p != NCOL(coef))
        error("number of columns do not match");
    
    dx = REAL(x);
    dc = REAL(coef);
    
    PROTECT(ret = allocVector(REALSXP, n));
    dr = REAL(ret);
    
    for (int i = 0; i < n; i++) {
        dr[i] = 0.0;
        for (int j = 0; j < p; j++) {
            idx = j * n + i;
            dr[i] += dx[idx] * dc[idx];
        }
    }
    
    UNPROTECT(1);
    return(ret);
}

SEXP R_offrowSums
(
    SEXP offset,
    SEXP x,
    SEXP coef
) {

    SEXP ret;
    double *dx, *doff, *dc, *dr;
    int p, n, idx;
    
    if (!isReal(offset))
        error("offset is not real");    
    if (!isReal(x))
        error("x is not real");
    if (!isReal(coef))
        error("coef is not real");    
    n = NROW(x);
    p = NCOL(x);
    if (n != NROW(coef))
        error("number of rows do not match");
    if (p != NCOL(coef))
        error("number of columns do not match");
    if (LENGTH(offset) != n)
        error("length(offset) does not match");
    
    dx = REAL(x);
    doff = REAL(offset);
    dc = REAL(coef);
    
    PROTECT(ret = allocVector(REALSXP, n));
    dr = REAL(ret);
    
    for (int i = 0; i < n; i++) {
        dr[i] = doff[i];
        for (int j = 0; j < p; j++) {
            idx = j * n + i;
            dr[i] += dx[idx] * dc[idx];
        }
    }
    
    UNPROTECT(1);
    return(ret);
}
