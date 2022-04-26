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

#define USE_FC_LEN_T

#include <Rconfig.h>
#include <R_ext/RS.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#ifndef FCONE
# define FCONE
#endif

#ifndef BLAS_extern
#define BLAS_extern extern
#endif

#define BLAS_DGEMV(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)         \
    F77_CALL(dgemv)(trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy  \
                    FCONE)

#define BLAS_DSYMV(uplo, n, alpha, a, lda, x, incx, beta, y, incy)             \
    F77_CALL(dsymv)(uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy FCONE)

#define BLAS_DTRSV(uplo, trans, diag, n, a, lda, x, incx)                      \
    F77_CALL(dtrsv)(uplo, trans, diag, &n, a, &lda, x, &incx FCONE FCONE FCONE)

#define BLAS_DGEMM(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c,    \
                   ldc) \
    F77_CALL(dgemm)(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb,      \
                    &beta, c, &ldc FCONE FCONE)

#define BLAS_DTRSM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)      \
    F77_CALL(dtrsm)(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb \
                    FCONE FCONE FCONE FCONE)

#define LAPACK_DPOTRF(uplo, n, a, lda, info)                                   \
    F77_CALL(dpotrf)(uplo, &n, a, &lda, &info FCONE)

/*
 * Note that parameters LI and LU in DSYEVR_ are not referenced, but
 * nevertheless must be valid references to avoid a segfault on Solaris.
 * Using a reference to parameter `n` as placeholder.
 */
#define LAPACK_DSYEVR_RANGE(uplo, n, a, lda, rangeLower, rangeUpper, abstol,   \
                            nevalues,evalues, evectors, ldevectors,            \
                            isuppevectors, work, lwork, iwork,                 \
                            liwork, info)                                      \
    F77_CALL(dsyevr)("V", "V", uplo, &n, a, &lda, &rangeLower, &rangeUpper, &n,\
                     &n, &abstol, &nevalues, evalues, evectors, &ldevectors,   \
                     isuppevectors, work, &lwork, iwork, &liwork, &info FCONE  \
                     FCONE FCONE)

#endif // pyinit_BLAS_h
