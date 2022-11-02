/******************************************************************************

SI UNITS DECLARATIONS/POSTINGS

These are units conversion factors to SI units (mks).
We let SI and mks denote generic SI units.

All derived units are defined ultimately in terms of the base
units m, kg, s, and Kelvin (and radians).  Other derived units may be added,
but they must be expressed in terms of previously defined units.

Note that temperature is not described as an absolute quantity, but only as
a relative quantity (as in "per degC").

Other bases can be defined by:
  Redefining the SI base units appropriately in terms of the new
  base units (all other units will be automatically adjusted to
  their proper values).

******************************************************************************/

#ifndef _UNITS_H
#define _UNITS_H

#include "dimen.h"
#define POST_UNIT_NULL
#include "units.inc"
#undef POST_UNIT_NULL

/* OTHER UNITS */

#define UNIT_AREA        UNIT_LENGTH*UNIT_LENGTH
#define UNIT_VOLUME      UNIT_LENGTH*UNIT_AREA
#define UNIT_PRESSURE    UNIT_MASS/(UNIT_LENGTH*UNIT_TIME*UNIT_TIME)
#define UNIT_VISCOSITY   UNIT_MASS/(UNIT_LENGTH*UNIT_TIME)

/* DECLARATIONS */

#define POST_UNIT_NAMES
static char unitName[][MAX_UNIT_NAME_LENGTH] = {
#include "units.inc"
};
#undef POST_UNIT_NAMES

#define POST_UNIT_VALUES
static double unitValue[] = {
#include "units.inc"
};
#undef POST_UNIT_VALUES

#define POST_UNIT_DIMEN
static int unitDimen[][NUMBER_DIMEN] = {
#include "units.inc"
};
#undef POST_UNIT_DIMEN

#endif

