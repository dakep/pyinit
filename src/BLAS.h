/*
 * BLAS.h
 * pyinit
 *
 * Created by David Kepplinger on 2016-22-01.
 * Copyright (c) 2016 David Kepplinger. All rights reserved.
 */

#ifndef pyinit_BLAS_h
#define pyinit_BLAS_h

#define BLAS_CHAR const char * const
#define BLAS_INT const int

#include <R_ext/RS.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#ifndef BLAS_extern
#define BLAS_extern extern
#endif

/* DGEMM - perform one of the matrix-matrix operations    */
/* C := alpha*op( A )*op( B ) + beta*C */
BLAS_extern void
F77_NAME(dgemm)(const char *transa, const char *transb, const int *m,
        const int *n, const int *k, const double *alpha,
        const double *a, const int *lda,
        const double *b, const int *ldb,
        const double *beta, double *c, const int *ldc);

BLAS_extern void
F77_NAME(dtrsv)(const char *uplo, const char *trans,
        const char *diag, const int *n,
        const double *a, const int *lda,
        double *x, const int *incx);

BLAS_extern void
F77_NAME(dtrsm)(const char *side, const char *uplo,
        const char *transa, const char *diag,
        const int *m, const int *n, const double *alpha,
        const double *a, const int *lda,
        double *b, const int *ldb);

#define BLAS_DGEMV(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)                          \
    F77_NAME(dgemv)(trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy)

#define BLAS_DSYMV(uplo, n, alpha, a, lda, x, incx, beta, y, incy)                              \
    F77_NAME(dsymv)(uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy)

#define BLAS_DTRSV(uplo, trans, diag, n, a, lda, x, incx)                                       \
    F77_NAME(dtrsv)(uplo, trans, diag, &n, a, &lda, x, &incx)

#define BLAS_DGEMM(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)                \
    F77_NAME(dgemm)(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc)

#define BLAS_DTRSM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)                       \
    F77_NAME(dtrsm)(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb)

#define LAPACK_DPOTRF(uplo, n, a, lda, info)                                                    \
    F77_NAME(dpotrf)(uplo, &n, a, &lda, &info)

/*
 * Note that parameters LI and LU in DSYEVR_ are not referenced, but nevertheless must be valid
 * references to avoid a segfault on Solaris. Using a reference to parameter `n` as placeholder.
 */
#define LAPACK_DSYEVR_RANGE(uplo, n, a, lda, rangeLower, rangeUpper, abstol, nevalues,          \
                            evalues, evectors, ldevectors, isuppevectors, work, lwork, iwork,   \
                            liwork, info)                                                       \
    F77_NAME(dsyevr)("V", "V", uplo, &n, a, &lda, &rangeLower, &rangeUpper, &n, &n,             \
                     &abstol, &nevalues, evalues, evectors, &ldevectors, isuppevectors, work,   \
                     &lwork, iwork, &liwork, &info)

#endif // pyinit_BLAS_h
