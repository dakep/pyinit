//
//  fastPY.h
//  pyinit-csrc
//
//  Created by David Kepplinger on 2017-10-21.
//  Copyright © 2017 David Kepplinger. All rights reserved.
//

#ifndef fastPY_h
#define fastPY_h

#include "Control.h"

/**
 * NOTE: `estimates` must be able to hold nvar * (3 * nvar + 1) + nobs elements.
 * The first nvar * (3 * nvar + 1) elements are the initial estimates, the remaining
 * nobs elements can be discared!
 */
int computePYEstimator(
    const double *restrict Xtr,
    const double *restrict y,
    const int nobs,
    const int nvar,
    const Control* ctrl,
    double *restrict estimates,
    double *restrict objFunScores
);


#endif /* fastPY_h */
