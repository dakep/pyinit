/*
 *  olsreg.c
 *  pyinit
 *
 *  Created by David Kepplinger on 2016-02-13.
 *  Copyright © 2016 David Kepplinger. All rights reserved.
 */

#include <string.h>
#include <Rmath.h>
#include <stdlib.h>

#include "olsreg.h"
#include "BLAS.h"
#include "evd.h"

static BLAS_INT BLAS_1L = 1;
static const double BLAS_0F = 0.0;
static const double BLAS_1F = 1.0;
static const double BLAS_M1F = -1.0;

static BLAS_CHAR BLAS_UPLO_UPPER = "U";
static BLAS_CHAR BLAS_TRANS_NO = "N";
static BLAS_CHAR BLAS_TRANS_TRANS = "T";
static BLAS_CHAR BLAS_DIAG_NO = "N";

int computeOLSCoefs(const double *restrict Xtr, const double *restrict y,
                    const int nobs, const int nvar,
					double *restrict coefs, AuxMemory *auxmem)
{
    int lapackInfo, retVal = OLS_COEFFICIENTS_OKAY;

    /* coefEst = t(X) %*% y */
    BLAS_DGEMV(BLAS_TRANS_NO, nvar, nobs, BLAS_1F, Xtr, nvar, y, BLAS_1L, BLAS_0F, coefs,
               BLAS_1L);

    /* Xsqrt = t(X) %*% X */
    BLAS_DGEMM(BLAS_TRANS_NO, BLAS_TRANS_TRANS, nvar, nvar, nobs,
               BLAS_1F, Xtr, nvar, Xtr, nvar,
               BLAS_0F, auxmem->Xsqrt, nvar);

    /* Xsqrt = chol(Xsqrt) */
    LAPACK_DPOTRF(BLAS_UPLO_UPPER, nvar, auxmem->Xsqrt, nvar, lapackInfo);

    auxmem->isXsqrtInverted = 0;

    if (lapackInfo == 0) {
        /* coefEst = inv(t(chol(Xsqrt))) %*% coefEst */
        BLAS_DTRSV(BLAS_UPLO_UPPER, BLAS_TRANS_TRANS, BLAS_DIAG_NO, nvar, auxmem->Xsqrt, nvar, coefs,
                   BLAS_1L);

        /* coefEst = inv(chol(Xsqrt)) %*% coefEst */
        BLAS_DTRSV(BLAS_UPLO_UPPER, BLAS_TRANS_NO, BLAS_DIAG_NO, nvar, auxmem->Xsqrt, nvar, coefs,
                   BLAS_1L);
    } else if (lapackInfo > 0) {
        /* Handle rank-deficient model here: Compute inv(Xsqrt) by EVD */
        int nevalues = 0, i, j;
        double *restrict workPtr;
        double *restrict evecPtr;
        double *restrict tmpXsqrtInv = (double *) malloc(nvar * nvar * sizeof(double));
        lapackInfo = 0;

        /* tmpXsqrtInv = t(X) %*% X */
        BLAS_DGEMM(BLAS_TRANS_NO, BLAS_TRANS_TRANS, nvar, nvar, nobs,
                   BLAS_1F, Xtr, nvar, Xtr, nvar,
                   BLAS_0F, tmpXsqrtInv, nvar);

        nevalues = symEigenValueDecomposition(BLAS_UPLO_UPPER, tmpXsqrtInv, nvar, auxmem);

        retVal = OLS_COEFFICIENTS_PSEUDO;

        if (nevalues > 0) {
            /*
             * First compute dblWorkMem = Q . A, with the first nevalues eigenvalues/vectors
             * where Q is the matrix of eigenvectors and A is the diagonal matrix with
             * sqrt(1 / eigenvalues)
             */
            resizeDblWorkAuxMemory(auxmem, nvar * nevalues);

            workPtr = auxmem->dblWorkMem;
            evecPtr = auxmem->eigenvectors;
            for (i = 0; i < nevalues; ++i) {
                auxmem->evalues[i] = sqrt(1. / auxmem->evalues[i]);
                for (j = 0; j < nvar; ++j) {
                    *workPtr = (*evecPtr) * auxmem->evalues[i];
                    ++evecPtr;
                    ++workPtr;
                }
            }

            /* Now compute tmpXsqrtInv = dblWorkMem . t(Q) */
            BLAS_DGEMM(BLAS_TRANS_NO, BLAS_TRANS_TRANS, nvar, nvar, nevalues,
                       BLAS_1F, auxmem->dblWorkMem, nvar, auxmem->eigenvectors, nvar,
                       BLAS_0F, tmpXsqrtInv, nvar);


            /* coefEst = XsqrtInv . XsqrtInv . coefEst */
            BLAS_DSYMV(BLAS_UPLO_UPPER, nvar, BLAS_1F, tmpXsqrtInv, nvar,
                       coefs, BLAS_1L, BLAS_0F, auxmem->dblWorkMem, BLAS_1L);
            BLAS_DSYMV(BLAS_UPLO_UPPER, nvar, BLAS_1F, tmpXsqrtInv, nvar,
                       auxmem->dblWorkMem, BLAS_1L, BLAS_0F, coefs, BLAS_1L);

            memcpy(auxmem->Xsqrt, tmpXsqrtInv, nvar * nvar * sizeof(double));

            auxmem->isXsqrtInverted = 1;
            free(tmpXsqrtInv);
        } else {
            lapackInfo = auxmem->intWorkMem[0];
        }
    }

    if (lapackInfo != 0) {
        auxmem->intWorkMem[0] = lapackInfo;
        return OLS_COEFFICIENTS_ERROR;
    }

    auxmem->intWorkMem[0] = 0;
    return retVal;
}


void computeResiduals(const double *restrict Xtr, const double *restrict y,
					  const int nobs, const int nvar,
					  const double *restrict coefs, double *restrict residuals)
{
    memcpy(residuals, y, nobs * sizeof(double));

    BLAS_DGEMV(BLAS_TRANS_TRANS, nvar, nobs,
               BLAS_M1F, Xtr, nvar,
               coefs, BLAS_1L,
               BLAS_1F, residuals, BLAS_1L);
}

