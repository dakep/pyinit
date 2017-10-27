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
    HUBER = 0,
    BISQUARE = 1,
    GAUSS = 5
} RhoFunctionName;

typedef struct ControlTag {
    const int numit;
    const double eps;
    const double keepResidThreshold;
    const double keepResidProportion;
    const double pscProportion;

    const double mscaleDelta;
    const double mscaleCc;
    const int mscaleMaxit;
    const double mscaleEps;
    const RhoFunctionName mscaleRhoFun;
} Control;

#endif /* Control_h */
