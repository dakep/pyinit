/*
 *  psc.c
 *  pyinit
 *
 *  Created by David Kepplinger on 2016-01-30.
 *  Copyright © 2016 David Kepplinger. All rights reserved.
 */

#include <float.h>
#include <stdlib.h>

#include "BLAS.h"
#include "olsreg.h"
#include "config.h"
#include "evd.h"


static BLAS_INT BLAS_M1L = -1;

static const double LAPACK_EV_ABSTOL = NUMERIC_EPS;
static const double LAPACK_EV_RANGE_LOWER = LAPACK_EV_MIN;
static const double LAPACK_EV_RANGE_UPPER = DBL_MAX;

int symEigenValueDecomposition(const char *const uplo, double *restrict matrix, const int n,
                               AuxMemory *const auxmem)
{
    int nevalues, lapackInfo;

    /* Compute needed size of working space */
    LAPACK_DSYEVR_RANGE(uplo, n, matrix, n,
                        LAPACK_EV_RANGE_LOWER, LAPACK_EV_RANGE_UPPER, LAPACK_EV_ABSTOL,
                        nevalues, auxmem->evalues,
                        auxmem->eigenvectors, n, auxmem->evectorsSupport,
                        auxmem->dblWorkMem, BLAS_M1L,
                        auxmem->intWorkMem, BLAS_M1L, lapackInfo);

    if (lapackInfo != 0) {
        auxmem->intWorkMem[0] = lapackInfo;
        return -1;
    }

    resizeDblWorkAuxMemory(auxmem, (int) auxmem->dblWorkMem[0]);
    resizeIntWorkAuxMemory(auxmem, auxmem->intWorkMem[0]);

    /* Actually perform eigenvalue decomposition */
    LAPACK_DSYEVR_RANGE(uplo, n, matrix, n,
                        LAPACK_EV_RANGE_LOWER, LAPACK_EV_RANGE_UPPER, LAPACK_EV_ABSTOL,
                        nevalues, auxmem->evalues,
                        auxmem->eigenvectors, n, auxmem->evectorsSupport,
                        auxmem->dblWorkMem, auxmem->dblWorkMemSize,
                        auxmem->intWorkMem, auxmem->intWorkMemSize, lapackInfo);

    /* One of the arguments was wrong */
    if (lapackInfo < 0) {
        auxmem->intWorkMem[0] = lapackInfo;
        return -1;
    }

    /*
     * Positive error codes probably indicate that some eigenvectors did
     * not converge --> we ignore that error
     */

    auxmem->intWorkMem[0] = 0;
    return nevalues;
}
