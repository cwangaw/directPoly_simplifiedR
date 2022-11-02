#ifndef _DIMEN_H
#define _DIMEN_H

/* ABSOLUTE TEMPERATURE */

#define MAX_DEGREE_NAME_LENGTH 5

#define DEG_K "degK"
#define DEG_C "degC"
#define DEG_F "degF"
#define INDEX_DEG_K 0
#define INDEX_DEG_C 1
#define INDEX_DEG_F 2

/* BASE SI UNITS */

#define NUMBER_DIMEN 4

#define UNIT_LENGTH 1
#define UNIT_MASS   1
#define UNIT_TIME   1
#define UNIT_TEMP   1
//#define UNIT_ANGLE  1

#define UNIT_ABSOLUTE_TEMP        DEG_K
#define UNIT_INDEX_ABSOLUTE_TEMP  INDEX_DEG_K

/* DIMENSIONS */

#define DIMEN_LENGTH  1,0,0,0
#define DIMEN_MASS    0,1,0,0
#define DIMEN_TIME    0,0,1,0
#define DIMEN_TEMP    0,0,0,1
//#define DIMEN_ANGLE   0,0,0,0,1

#define DIMEN_NUMBER  0,0,0,0

#define DIMEN_PER_LENGTH       -1,0,0,0
#define DIMEN_AREA             2,0,0,0
#define DIMEN_VOLUME           3,0,0,0
#define DIMEN_VELOCITY         1,0,-1,0
#define DIMEN_ACCEL            1,0,-2,0
#define DIMEN_MASS_CONC        -3,1,0,0
#define DIMEN_MOLAR_CONC       -3,0,0,0
#define DIMEN_FORCE            1,1,-2,0
#define DIMEN_ENERGY           2,1,-2,0
#define DIMEN_PRESSURE         -1,1,-2,0
#define DIMEN_PER_PRESSURE     1,-1,2,0
#define DIMEN_VISCOSITY        -1,1,-1,0
#define DIMEN_PER_TIME         0,0,-1,0
#define DIMEN_VOLUME_PER_TIME  3,0,-1,0
#define DIMEN_DIFFUSION        2,0,-1,0

#endif
