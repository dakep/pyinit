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

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Compute the OLS coefficients using the Cholesky decomposition which
 * will be stored in `Xsqrt`.
 *
 * @return A return value not equal 0 indicates an error.
 */
int computeOLSCoefs(const double *restrict Xtr, const double *restrict y,
					const int nobs, const int nvar,
					double *restrict coefs, double *restrict Xsqrt);

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
