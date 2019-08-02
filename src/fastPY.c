/*
 *  fastPY.c
 *  pyinit
 *
 *  Created by David Kepplinger on 2016-01-28.
 *  Copyright © 2016 David Kepplinger. All rights reserved.
 */

#include <Rmath.h>
#include <stdlib.h>
#include <float.h>

#include "fastPY.h"

#include "BLAS.h"
#include "Control.h"
#include "mscale.h"
#include "PartialSort.h"
#include "AuxMemory.h"
#include "olsreg.h"
#include "psc.h"

static int filterDataThreshold(const double *restrict X, const double *restrict y,
                                double *restrict newX, double *restrict newY,
                                const int nobs, const int nvar,
                                const double *restrict values, const double threshold,
                                CompareFunction compare);

/*
 * Compare functions
 */
static double lessThan(const double a, const double b);
static double greaterThan(const double a, const double b);
static double absoluteLessThan(const double a, const double b);


/***************************************************************************************************
 *
 * Main function to compute the initial estimator
 *
 **************************************************************************************************/
int computePYEstimator(const double *restrict Xtr, const double *restrict y,
                       const int nobs, const int nvar, const Control* ctrl,
                       double *restrict estimates, double *restrict objFunScores)
{
    AuxMemory auxmem;
    double *restrict bestCoefEst;
    double *restrict currentEst;
    double *restrict currentXtr = (double*) malloc(nobs * nvar * sizeof(double));
    double *restrict currentY = (double*) malloc(nobs * sizeof(double));
    double *restrict filteredXtr = (double*) malloc(nobs * nvar * sizeof(double));
    double *restrict filteredY = (double*) malloc(nobs * sizeof(double));
    double *restrict pscs = (double*) malloc(nobs * nvar * sizeof(double));
    double *const minObjective = objFunScores;
    double * tmpObjective;
    double scaledThreshold = 0;
    int currentNobs = nobs;
    int iter = 0;
    int j;
    int numPSCs = 0;
    int linalgError = 0;
    int compCoefsStatus = 0;
    double *restrict currentPSC;
    int filteredNobs;
    double diff, normPrevBest = 0, normBest = 0;
    RhoFunction rhoFun = getRhoFunctionByName(ctrl->mscaleRhoFun);

    initAuxMemory(&auxmem);
    resizeAuxMemory(&auxmem, nvar, nobs);

    memcpy(currentXtr, Xtr, nobs * nvar * sizeof(double));
    memcpy(currentY, y, nobs * sizeof(double));

    /*
     * At first, we will leave the front of the estimates-matrix empty -- this is the
     * place where the best estimate of the previous run will be stored
     */
    bestCoefEst = estimates;
    currentEst = estimates + nvar;
    *minObjective = DBL_MAX;

    while(1) {
        tmpObjective = objFunScores + 1;

        memcpy(filteredXtr, currentXtr, currentNobs * nvar * sizeof(double));
        memcpy(filteredY, currentY, currentNobs * sizeof(double));

        /* 1. Estimate coefficients for the residuales filtered data (currentXtr) */
        compCoefsStatus = computeOLSCoefs(filteredXtr, filteredY, currentNobs, nvar, currentEst,
                                          &auxmem);
        linalgError = auxmem.intWorkMem[0];
        if (linalgError != 0) {
            break;
        }

        computeResiduals(Xtr, y, nobs, nvar, currentEst, auxmem.residuals);
        *tmpObjective = mscale(auxmem.residuals, nobs, ctrl->mscaleDelta, ctrl->mscaleEps,
                               ctrl->mscaleMaxit, rhoFun, ctrl->mscaleCc);

        if (*tmpObjective < *minObjective) {
            bestCoefEst = currentEst;
            *minObjective = *tmpObjective;
        }
        ++tmpObjective;

        /* 2. Calculate PSC for current work data */
        computeResiduals(currentXtr, currentY, currentNobs, nvar, currentEst, auxmem.residuals);
        numPSCs = calculatePSCs(pscs, &auxmem, currentXtr, currentY, currentNobs, nvar);

        if (numPSCs < 0) {
            linalgError = -numPSCs;
            break;
        }

        for(j = 0; j < numPSCs; ++j) {
            currentPSC = pscs + currentNobs * j;
            /* 4.1.1 Thin out X and y based on large values of PSCs */
            scaledThreshold = getQuantile(currentPSC, currentNobs, ctrl->pscProportion,
                                          lessThan);

            filteredNobs = filterDataThreshold(currentXtr, currentY, filteredXtr, filteredY,
                                               currentNobs, nvar, currentPSC,
                                               scaledThreshold, lessThan);

            /* 4.1.2. Estimate coefficients */
            currentEst += nvar;
            if (filteredNobs > 0) {
                compCoefsStatus = computeOLSCoefs(filteredXtr, filteredY, filteredNobs, nvar,
                                                  currentEst, &auxmem);
            } else {
                compCoefsStatus = OLS_COEFFICIENTS_ERROR;
            }

            if (compCoefsStatus == OLS_COEFFICIENTS_OKAY) {
                computeResiduals(Xtr, y, nobs, nvar, currentEst, auxmem.residuals);
                *tmpObjective = mscale(auxmem.residuals, nobs, ctrl->mscaleDelta, ctrl->mscaleEps,
                                       ctrl->mscaleMaxit, rhoFun, ctrl->mscaleCc);

                if (*tmpObjective < *minObjective) {
                    *minObjective = *tmpObjective;
                    bestCoefEst = currentEst;
                }
            } else {
                memset(currentEst, 0, nvar * sizeof(double));
            }
            ++tmpObjective;

            /* 4.2.1. Thin out X and y based on large values of PSCs */
            scaledThreshold = getQuantile(currentPSC, currentNobs, ctrl->pscProportion,
                                          greaterThan);

            filteredNobs = filterDataThreshold(currentXtr, currentY, filteredXtr, filteredY,
                                               currentNobs, nvar, currentPSC,
                                               scaledThreshold, greaterThan);

            /* 4.2.2. Estimate coefficients */
            currentEst += nvar;
            if (filteredNobs > 0) {
                compCoefsStatus = computeOLSCoefs(filteredXtr, filteredY, filteredNobs, nvar,
                                                  currentEst, &auxmem);
            } else {
                compCoefsStatus = OLS_COEFFICIENTS_ERROR;
            }

            if (compCoefsStatus == OLS_COEFFICIENTS_OKAY) {
                computeResiduals(Xtr, y, nobs, nvar, currentEst, auxmem.residuals);
                *tmpObjective = mscale(auxmem.residuals, nobs, ctrl->mscaleDelta, ctrl->mscaleEps,
                                       ctrl->mscaleMaxit, rhoFun, ctrl->mscaleCc);

                if (*tmpObjective < *minObjective) {
                    *minObjective = *tmpObjective;
                    bestCoefEst = currentEst;
                }
            } else {
                memset(currentEst, 0, nvar * sizeof(double));
            }
            ++tmpObjective;


            /* 4.3.1. Thin out X and y based on large values of PSCs */
            scaledThreshold = getQuantile(currentPSC, currentNobs, ctrl->pscProportion,
                                          absoluteLessThan);

            filteredNobs = filterDataThreshold(currentXtr, currentY, filteredXtr, filteredY,
                                               currentNobs, nvar, currentPSC,
                                               scaledThreshold, absoluteLessThan);

            /* 4.3.2. Estimate coefficients */
            currentEst += nvar;
            if (filteredNobs > 0) {
                compCoefsStatus = computeOLSCoefs(filteredXtr, filteredY, filteredNobs, nvar,
                                                  currentEst, &auxmem);
            } else {
                compCoefsStatus = OLS_COEFFICIENTS_ERROR;
            }

            if (compCoefsStatus == OLS_COEFFICIENTS_OKAY) {
                computeResiduals(Xtr, y, nobs, nvar, currentEst, auxmem.residuals);
                *tmpObjective = mscale(auxmem.residuals, nobs, ctrl->mscaleDelta, ctrl->mscaleEps,
                                       ctrl->mscaleMaxit, rhoFun, ctrl->mscaleCc);

                if (*tmpObjective < *minObjective) {
                    *minObjective = *tmpObjective;
                    bestCoefEst = currentEst;
                }
            } else {
                memset(currentEst, 0, nvar * sizeof(double));
            }
            ++tmpObjective;
        }

        if (linalgError != 0) {
            break;
        }

        /*
         * If iter > 0, compute the difference between the new best coefficient estimate and the previously
         * best coefficient estimate.
         */
        if (iter > 0) {
            diff = 0;
            normBest = 0;
            for (j = 0; j < nvar; ++j) {
                diff += fabs(bestCoefEst[j] - estimates[j]);
                normBest += fabs(bestCoefEst[j]);
            }
        }

        /* 5. Store best estimate for later at the beginning of estimates */
        memcpy(estimates, bestCoefEst, nvar * sizeof(double));
        bestCoefEst = estimates;
        currentEst = estimates + nvar;

        /*
         * Check if the difference is negligible or we reached the maximum number of iterations.
         */
        if (iter++ > 0 && ((diff < ctrl->eps * normPrevBest) || (iter >= ctrl->numit))) {
            break;
        }

        normPrevBest = normBest;

        /* 6. Calculate residuals with best coefficient estimate again */
        computeResiduals(Xtr, y, nobs, nvar, bestCoefEst, auxmem.residuals);

        if (ctrl->keepResidThreshold <= 0) {
            scaledThreshold = getQuantile(auxmem.residuals, nobs, ctrl->keepResidProportion,
                                          absoluteLessThan);
        } else {
            scaledThreshold = ctrl->keepResidThreshold * (*minObjective);
        }

        currentNobs = filterDataThreshold(Xtr, y, currentXtr, currentY, nobs, nvar, auxmem.residuals,
                                          scaledThreshold, absoluteLessThan);
    }

    freeAuxMemory(&auxmem);

    free(pscs);
    free(currentXtr);
    free(currentY);
    free(filteredXtr);
    free(filteredY);

    if (linalgError != 0) {
        Rf_error("There was an error in one of the calls to LINPACK (%d)", linalgError);
    }

    return MAX_NUM_PSCS(numPSCs);
}


/***************************************************************************************************
 *
 * Helper functions to filter data based on a threshold
 *
 **************************************************************************************************/
static int filterDataThreshold(const double *restrict Xtr, const double *restrict y,
                                double *restrict newXtr, double *restrict newY,
                                const int nobs, const int nvar,
                                const double *restrict values, const double threshold,
                                CompareFunction compare)
{
    int i, toCount = 0;
    double *restrict toYIt = newY;
    const double *restrict fromYIt = y;

    int copyRows = 0;

    const double *restrict startFromX = Xtr;
    double *restrict startToX = newXtr;

    for (i = 0; i < nobs; ++i, ++fromYIt) {
        if (compare(values[i], threshold) < 0) {
            /* Copy value from y to newY */
            *toYIt = *fromYIt;
            ++toYIt;
            ++toCount;
            ++copyRows;
        } else {
            /*
             * This row (column in Xtr) should not be copied.
             * So copy everything so far
             */
            if (copyRows > 0) {
                memcpy(startToX, startFromX, copyRows * nvar * sizeof(double));
                startToX += copyRows * nvar;
            }

            startFromX += (copyRows + 1) * nvar;
            copyRows = 0;
        }
    }

    /*
     * Copy last chunk of data
     */
    if (copyRows > 0) {
        memcpy(startToX, startFromX, copyRows * nvar * sizeof(double));
        startToX += copyRows * nvar;
    }

    return toCount;
}

/***************************************************************************************************
 *
 * Static compare functions
 *
 **************************************************************************************************/
static double lessThan(const double a, const double b)
{
    return a - b;
}

static double greaterThan(const double a, const double b)
{
    return b - a;
}

static double absoluteLessThan(const double a, const double b)
{
    return fabs(a) - fabs(b);
}
