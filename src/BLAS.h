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

#ifndef La_extern
#define La_extern extern
#endif

#define BLAS_R_NAME(x) F77_NAME(x)
#define LAPACK_R_NAME(x) F77_NAME(x)

/* DGEMM - perform one of the matrix-matrix operations    */
/* C := alpha*op( A )*op( B ) + beta*C */
BLAS_extern void
F77_NAME(dgemm)(const char *transa, const char *transb, const int *m,
        const int *n, const int *k, const double *alpha,
        const double *a, const int *lda,
        const double *b, const int *ldb,
        const double *beta, double *c, const int *ldc);


/* DSYR - perform the symmetric rank 1 operation A := alpha*x*x' + A */
BLAS_extern void
F77_NAME(dsyr)(const char *uplo, const int *n, const double *alpha,
           const double *x, const int *incx,
           double *a, const int *lda);

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

#define BLAS_DGEMV(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)							\
    BLAS_R_NAME(dgemv)(trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy)

#define BLAS_DSYMV(uplo, n, alpha, a, lda, x, incx, beta, y, incy)								\
    BLAS_R_NAME(dsymv)(uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy)

#define BLAS_DTRSV(uplo, trans, diag, n, a, lda, x, incx)                                       \
    BLAS_R_NAME(dtrsv)(uplo, trans, diag, &n, a, &lda, x, &incx)

#define BLAS_DGEMM(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)				\
	BLAS_R_NAME(dgemm)(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc)

#define BLAS_DDOT(n, x, incx, y, incy)															\
    BLAS_R_NAME(ddot)(&n, x, &incx, y, &incy)

#define BLAS_DSCAL(n, alpha, x, incx)															\
    BLAS_R_NAME(dscal)(&n, &alpha, x, &incx)

#define BLAS_DSYR(uplo, n, alpha, x, incx, a, lda)												\
	BLAS_R_NAME(dsyr)(uplo, &n, &alpha, x, &incx, a, &lda)

#define BLAS_DTRSM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)                       \
    BLAS_R_NAME(dtrsm)(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb)

#define BLAS_DTRMV(uplo, trans, diag, n, a, lda, x, incx)                                       \
    BLAS_R_NAME(dtrmv)(uplo, trans, diag, &n, a, &lda, x, &incx)


#define LAPACK_DPOTRF(uplo, n, a, lda, info)                                                    \
    LAPACK_R_NAME(dpotrf)(uplo, &n, a, &lda, &info)

#define LAPACK_DSYEVR_ALL_VECTORS(uplo, n, a, lda, abstol, nevalues, evalues, evectors,         \
                                  ldevectors, isuppevectors, work, lwork, iwork, liwork, info)  \
    LAPACK_R_NAME(dsyevr)("V", "A", uplo, &n, a, &lda, NULL, NULL, NULL, NULL, &abstol,         \
                     &nevalues, evalues, evectors, &ldevectors, isuppevectors, work, &lwork,    \
                     iwork, &liwork, &info)


#define LAPACK_DSYEVR_RANGE(uplo, n, a, lda, rangeLower, rangeUpper, abstol, nevalues,          \
                            evalues, evectors, ldevectors, isuppevectors, work, lwork, iwork,   \
                            liwork, info)                                                       \
    LAPACK_R_NAME(dsyevr)("V", "V", uplo, &n, a, &lda, &rangeLower, &rangeUpper, NULL, NULL,    \
                     &abstol, &nevalues, evalues, evectors, &ldevectors, isuppevectors, work,   \
                     &lwork, iwork, &liwork, &info)

#define LAPACK_DGELS(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)                      \
    LAPACK_R_NAME(dgels)(trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info)

#endif
