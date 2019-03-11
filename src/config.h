/*
 * config.h
 * pyinit
 *
 * Created by David Kepplinger on 2016-01-22.
 * Copyright (c) 2015 David Kepplinger. All rights reserved.
 */

#ifndef pyinit_config_h
#define pyinit_config_h

#include "autoconfig.h"

#ifdef HAVE_INTTYPES_H
#   include <inttypes.h>
#endif

#ifdef HAVE_STDINT_H
#   include <stdint.h>
#endif

#define NUMERIC_EPS 1e-16
#define LAPACK_EV_MIN 1e-12

#endif


