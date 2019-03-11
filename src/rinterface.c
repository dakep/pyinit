//
//  rinterface.c
//  pyinit
//
//  Created by David Kepplinger on 2017-10-21.
//  Copyright © 2017 David Kepplinger. All rights reserved.
//
#include "config.h"

#include "mscale.h"
#include "fastPY.h"
#include "psc.h"

#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/**
 * Calculate the M-Scale of a vector of numbers
 *
 * @param Rvalues  numeric vector of REAL's
 * @param Rb       numeric target average (right side) in the M-scale equation
 * @param Rcc      numeric constant for the rho function
 * @param Rmaxit   integer maximum number of iterations
 * @param Reps     numeric relative tolerance for convergence
 * @param Rrho_fun integer rho function to use (see Control.h for possible values).
 *      If the selected rho function is unknown, the bisquare function will be used by default.
 *
 * @return numeric The M-scale
 */
SEXP C_mscale(
    SEXP Rvalues,
    SEXP Rb,
    SEXP Rcc,
    SEXP Rmaxit,
    SEXP Reps,
    SEXP Rrho_fun
);

/**
 * Calculate the Pena-Yohai initial estimator
 *
 * @param RXtr                 numeric The (nvar by nobs) transposed X matrix
 * @param Ry	               numeric The (nobs) y vector
 * @param Rnumit               integer The maximum number of iterations
 * @param Reps                 numeric The relative tolerance for convergence
 * @param RkeepResidThreshold  numeric threshold for residuals to keep
 * @param RkeepResidProportion numeric proportion of observations to keep (based on residuals)
 * @param RpscProportion       numeric proportion of observations to keep (based on PSCs)
 * @param RmscaleDelta         numeric
 * @param RmscaleCc            numeric
 * @param RmscaleMaxit         integer
 * @param RmscaleEps           numeric
 * @param RmscaleRhoFun        integer
 *
 * @return Returns a list with three elements:
 *      item 1: number of initial estimators returned
 *      item 2: matrix of initial estimators
 *      item 3: value of the objective function for each estimator
 */
SEXP C_initpy(
    SEXP RXtr,
    SEXP Ry,
    SEXP Rnumit,
    SEXP Reps,
    SEXP RkeepResidThreshold,
    SEXP RkeepResidProportion,
    SEXP RpscProportion,
    SEXP RmscaleDelta,
    SEXP RmscaleCc,
    SEXP RmscaleMaxit,
    SEXP RmscaleEps,
    SEXP RmscaleRhoFun
);

static void getMatDims(SEXP matrix, int* nrows, int* ncols);

/**
 * .C entry point definitions for R
 */
static const R_CallMethodDef exportedCallMethods[] = {
    {"C_initpy", (DL_FUNC) &C_initpy, 12},
    {"C_mscale", (DL_FUNC) &C_mscale, 6},
    {NULL, NULL, 0}
};


void R_init_pyinit(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, exportedCallMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

/**
 * Calculate the M-Scale of a vector of numbers
 */
SEXP C_mscale(SEXP Rvalues, SEXP Rb, SEXP Rcc, SEXP Rmaxit, SEXP Reps, SEXP Rrho_fun)
{
    R_len_t nvalues = Rf_length(Rvalues);
    SEXP Rscale = PROTECT(Rf_allocVector(REALSXP, 1));
    RhoFunction rho_fun = getRhoFunctionByName(*INTEGER(Rrho_fun));

    *REAL(Rscale) = mscale(
        REAL(Rvalues),
        (int) nvalues,
        *REAL(Rb),
        *REAL(Reps),
        *INTEGER(Rmaxit),
        rho_fun,
        *REAL(Rcc)
    );

    UNPROTECT(1);
    return Rscale;
}



/**
 * Calculate the Pena-Yohai initial estimator
 */
SEXP C_initpy(SEXP RXtr, SEXP Ry, SEXP Rnumit,
              SEXP Reps, SEXP RkeepResidThreshold, SEXP RkeepResidProportion,
              SEXP RpscProportion, SEXP RmscaleDelta, SEXP RmscaleCc,
              SEXP RmscaleMaxit, SEXP RmscaleEps, SEXP RmscaleRhoFun)
{
    Control ctrl = {
        .numit = *INTEGER(Rnumit),
        .eps = *REAL(Reps),
        .keepResidThreshold = *REAL(RkeepResidThreshold),
        .keepResidProportion = *REAL(RkeepResidProportion),
        .pscProportion = *REAL(RpscProportion),
        .mscaleDelta = *REAL(RmscaleDelta),
        .mscaleCc = *REAL(RmscaleCc),
        .mscaleMaxit = *INTEGER(RmscaleMaxit),
        .mscaleEps = *REAL(RmscaleEps),
        .mscaleRhoFun = *INTEGER(RmscaleRhoFun)
    };

    int nobs = 0, nvar = 0;
    getMatDims(RXtr, &nvar, &nobs);

    SEXP result = PROTECT(Rf_allocVector(VECSXP, 3));
    SEXP numCoefs = PROTECT(Rf_allocVector(INTSXP, 1));
    SEXP coefs = PROTECT(Rf_allocVector(REALSXP, MAX_NUM_PSCS(nvar) * nvar));
    SEXP objF = PROTECT(Rf_allocVector(REALSXP, MAX_NUM_PSCS(nvar)));

    *INTEGER(numCoefs) = computePYEstimator(
        REAL(RXtr),
        REAL(Ry),
        nobs,
        nvar,
        &ctrl,
        REAL(coefs),
        REAL(objF)
    );

    SET_VECTOR_ELT(result, 0, numCoefs);
    SET_VECTOR_ELT(result, 1, coefs);
    SET_VECTOR_ELT(result, 2, objF);

    UNPROTECT(4);
    return result;
}

static inline void getMatDims(SEXP matrix, int* nrows, int* ncols)
{
    SEXP Rdims;
    int* dims;
    PROTECT(Rdims = Rf_getAttrib(matrix, R_DimSymbol));
    dims = INTEGER(Rdims);
    *nrows = dims[0];
    *ncols = dims[1];
    UNPROTECT(1);
}


