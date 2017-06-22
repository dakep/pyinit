/*
 *  Control.h
 *  pyinit
 *
 *  Created by David Kepplinger on 2016-01-30.
 *  Copyright © 2016 David Kepplinger. All rights reserved.
 */

#ifndef Control_h
#define Control_h

typedef enum RhoFunctionNameTag {
    BISQUARE = 0,
    HUBER = 1,
    GAUSS_WEIGHT = 2
} RhoFunctionName;

typedef struct ControlTag {
    const int numIt;
    const double eps;
    const double residThreshold;
    const double residProportion;
    const double pscProportion;

    const double mscaleB;
    const double mscaleCC;
    const int mscaleMaxIt;
    const double mscaleEPS;
    const RhoFunctionName mscaleRhoFun;
} Control;

#endif /* Control_h */
