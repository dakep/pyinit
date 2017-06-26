/*
 *  olsreg.h
 *  pyinit
 *
 *  Created by David Kepplinger on 2016-02-13.
 *  Copyright © 2016 David Kepplinger. All rights reserved.
 */

#ifndef olsreg_h
#define olsreg_h

#include "config.h"
#include "AuxMemory.h"

#ifdef __cplusplus
extern "C" {
#endif

#define OLS_COEFFICIENTS_OKAY 0
#define OLS_COEFFICIENTS_PSEUDO 1
#define OLS_COEFFICIENTS_ERROR 2

/**
 * Compute the OLS coefficients using the Cholesky decomposition which
 * will be stored in `Xsqrt`.
 * If the system is rank-deficient, a slow SVD is used instead and `Xsqrt` will
 * be used to store the Moore-Penrose pseudo-inverse of (t(X) . X)^1/2
 *
 * @return 0 ... coefficients have been computed. Cholesky decomposition is stored in Xsqrt
 *         1 ... coefficients have been computed. Xsqrt holds the Moore-Penrose pseudo inverse
 *               of (t(X) . X)^1/2
 *         2 ... LAPACK routine gave an error. Error code is stored in the first element
 *               of auxmem->intWorkMem.
 */
int computeOLSCoefs(const double *restrict Xtr, const double *restrict y,
					const int nobs, const int nvar,
					double *restrict coefs, AuxMemory *auxmem);

/**
 * Compute residuals for linear regression given the coefficients and the data
 */
void computeResiduals(const double *restrict Xtr, const double *restrict y,
					  const int nobs, const int nvar,
					  const double *restrict coefs, double *restrict residuals);

#ifdef __cplusplus
}
#endif

#endif /* olsreg_h */
