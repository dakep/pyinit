/*
 *  evd.h
 *  pyinit
 *
 *  Created by David Kepplinger on 2017-06-25.
 *  Copyright © 2017 David Kepplinger. All rights reserved.
 */

#ifndef evd_h
#define evd_h

#include "AuxMemory.h"

/**
 * Compute the Eigenvalue Decomposition of a symmetric matrix
 *
 * @param uplo	 use the upper 'U' or lower 'L' triangular part of `A`
 * @param A      symmetric matrix
 * @param n      number of rows/columns in `A`
 * @param auxmem Auxilliary memory used during computation
 * @return Returns the number of eigenvectors/eigenvalues or values less than 0 if an error
 *		   occured. The error code is stored as the first element in auxmem->intWorkMem
 */
int symEigenValueDecomposition(const char *const uplo, double *restrict A, const int n,
                               AuxMemory *const auxmem);

#endif /* evd_h */
