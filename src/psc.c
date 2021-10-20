/*
 *  psc.c
 *  pyinit
 *
 *  Created by David Kepplinger on 2016-01-30.
 *  Copyright © 2016 David Kepplinger. All rights reserved.
 */

#include "BLAS.h"  /* must come first! */

#include <float.h>
#include <stdlib.h>

#include "olsreg.h"
#include "config.h"
#include "psc.h"
#include "evd.h"


static const double BLAS_0F = 0.0;
static const double BLAS_1F = 1.0;

static BLAS_CHAR BLAS_SIDE_LEFT = "L";
static BLAS_CHAR BLAS_UPLO_UPPER = "U";
static BLAS_CHAR BLAS_TRANS_NO = "N";
static BLAS_CHAR BLAS_TRANS_TRANS = "T";
static BLAS_CHAR BLAS_DIAG_NO = "N";

int calculatePSCs(double *restrict pscs, AuxMemory* auxmem,
                  const double *restrict Xtr, const double *restrict y,
                  const int nobs, const int nvar)
{
    int i, j;
    int nevalues = 0;
    /*
     * G can point to the same address as auxmem->XtXinvX, because we don't need
     * auxmem->XtXinvX anymore when G is computed.
     * The same is true for Q and H, as well as Z and XtXinvX
     */
    double *restrict XtXinvX = pscs;
    double *restrict G = pscs;
    const double *restrict iterX;
    const double *restrict iterXtXinvX;
    const double *restrict iterXsqrtInvX;
    double *iterG;
    double Wii;

    /**
     * Make sure we have enough space in the auxiliary memory
     */
    resizeAuxMemory(auxmem, nvar, nobs);

    if (auxmem->isXsqrtInverted == 1) {
        /* auxmem->Xsqrt is the inverse of (t(X) . X)^1/2 */

        /* XsqrtInvX = XsqrtInv . t(X) */
        BLAS_DGEMM(BLAS_TRANS_NO, BLAS_TRANS_NO, nvar, nobs, nvar,
            BLAS_1F, auxmem->Xsqrt, nvar, Xtr, nvar,
            BLAS_0F, auxmem->XsqrtInvX, nvar);

        /* XtXinvX = XsqrtInv . XsqrtInvX */
        BLAS_DGEMM(BLAS_TRANS_NO, BLAS_TRANS_NO, nvar, nobs, nvar,
            BLAS_1F, auxmem->Xsqrt, nvar, auxmem->XsqrtInvX, nvar,
            BLAS_0F, XtXinvX, nvar);
    } else {
        /* auxmem->Xsqrt is the Cholesky decomposition of t(X) . X */
        memcpy(auxmem->XsqrtInvX, Xtr, nobs * nvar * sizeof(double));

        /* t(XsqrtInvX) = t(inv(Xsqrt)) %*% t(X) */
        BLAS_DTRSM(BLAS_SIDE_LEFT, BLAS_UPLO_UPPER, BLAS_TRANS_TRANS, BLAS_DIAG_NO,
                   nvar, nobs, BLAS_1F, auxmem->Xsqrt, nvar, auxmem->XsqrtInvX, nvar);

        memcpy(XtXinvX, auxmem->XsqrtInvX, nobs * nvar * sizeof(double));

        /* XtXinvX = inv(Xsqrt) %*% t(XsqrtInvX) */
        BLAS_DTRSM(BLAS_SIDE_LEFT, BLAS_UPLO_UPPER, BLAS_TRANS_NO, BLAS_DIAG_NO,
                   nvar, nobs, BLAS_1F, auxmem->Xsqrt, nvar, XtXinvX, nvar);
    }

    /* t(G) = t(XsqrtInvX) %*% diag(W[i, i]) */
    iterG = G;
    iterX = Xtr;
    iterXtXinvX = XtXinvX;
    iterXsqrtInvX = auxmem->XsqrtInvX;

    for (i = 0; i < nobs; ++i) {
        /*
         * Set Wii to i-th diagonal element of H,
         * i.e., inner product of x_i and i-th column of XtXinvX
         */
        Wii = 0;
        for (j = 0; j < nvar; ++j, ++iterX, ++iterXtXinvX) {
            Wii += (*iterX) * (*iterXtXinvX);
        }

        /* W[i, i] = */
        Wii = auxmem->residuals[i] / (1 - Wii);

        /* G[ , i] = W[i, i] * XsqrtInvX[, i] */
        for (j = 0; j < nvar; ++j, ++iterG, ++iterXsqrtInvX) {
            (*iterG) = Wii * (*iterXsqrtInvX);
        }
    }


    /* Q = t(G) %*% G */
    BLAS_DGEMM(BLAS_TRANS_NO, BLAS_TRANS_TRANS, nvar, nvar, nobs,
               BLAS_1F, G, nvar, G, nvar,
               BLAS_0F, auxmem->Q, nvar);

    nevalues = symEigenValueDecomposition(BLAS_UPLO_UPPER, auxmem->Q, nvar, auxmem);

    if (nevalues > 0) {
        BLAS_DGEMM(BLAS_TRANS_TRANS, BLAS_TRANS_NO, nobs, nevalues, nvar,
                   BLAS_1F, auxmem->XsqrtInvX, nvar, auxmem->eigenvectors, nvar,
                   BLAS_0F, pscs, nobs);
    } else if (nevalues < 0) {
        auxmem->intWorkMem[0] = (int) auxmem->evalues[0];
    }

    return nevalues;
}
