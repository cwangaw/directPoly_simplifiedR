/*
        THE CENTER FOR SUBSURFACE MODELING
      PARALLEL SUBSURFACE SIMULATOR (PARSSIM)

            UNITS DECLARATIONS/POSTINGS
*/

/*****************************************************************************/

#ifdef POST_UNIT_NAMES
#define POST_UNIT(name,value,length,mass,time,temp) name
#endif

#ifdef POST_UNIT_VALUES
#define POST_UNIT(name,value,length,mass,time,temp) value
#endif

#ifdef POST_UNIT_DIMEN
#define POST_UNIT(name,value,length,mass,time,temp) {length,mass,time,temp}
#endif

#define COMMA ,

#ifdef POST_UNIT_NULL
#undef COMMA
#define COMMA
#define POST_UNIT(name,value,length,mass,time,temp)
#endif

/*****************************************************************************/

#ifndef MAX_UNIT_NAME_LENGTH
#define MAX_UNIT_NAME_LENGTH 13
#endif

/* LENGTH */

  POST_UNIT("m",             1.*UNIT_LENGTH,              1,0,0,0) COMMA
  POST_UNIT("meters",        1.*UNIT_LENGTH,              1,0,0,0) COMMA

  POST_UNIT("cm",            .01*UNIT_LENGTH,             1,0,0,0) COMMA
  POST_UNIT("centimeters",   .01*UNIT_LENGTH,             1,0,0,0) COMMA

  POST_UNIT("km",            1000.*UNIT_LENGTH,           1,0,0,0) COMMA
  POST_UNIT("kilometers",    1000.*UNIT_LENGTH,           1,0,0,0) COMMA

  POST_UNIT("in",            .0254*UNIT_LENGTH,           1,0,0,0) COMMA
  POST_UNIT("inches",        .0254*UNIT_LENGTH,           1,0,0,0) COMMA

  POST_UNIT("ft",            12*.0254*UNIT_LENGTH,        1,0,0,0) COMMA
  POST_UNIT("feet",          12*.0254*UNIT_LENGTH,        1,0,0,0) COMMA

  POST_UNIT("mi",            63360*.0254*UNIT_LENGTH,     1,0,0,0) COMMA
  POST_UNIT("miles",         63360*.0254*UNIT_LENGTH,     1,0,0,0) COMMA

/* MASS */

  POST_UNIT("kg",            1.*UNIT_MASS,                0,1,0,0) COMMA
  POST_UNIT("kilograms",     1.*UNIT_MASS,                0,1,0,0) COMMA

  POST_UNIT("g",             1.e-3*UNIT_MASS,             0,1,0,0) COMMA
  POST_UNIT("grams",         1.e-3*UNIT_MASS,             0,1,0,0) COMMA

  POST_UNIT("mg",            1.e-6*UNIT_MASS,             0,1,0,0) COMMA
  POST_UNIT("milligrams",    1.e-6*UNIT_MASS,             0,1,0,0) COMMA

  POST_UNIT("lb",            .45359237*UNIT_MASS,         0,1,0,0) COMMA
  POST_UNIT("pounds",        .45359237*UNIT_MASS,         0,1,0,0) COMMA

/* TIME */

  POST_UNIT("s",             1.*UNIT_TIME,                0,0,1,0) COMMA
  POST_UNIT("sec",           1.*UNIT_TIME,                0,0,1,0) COMMA
  POST_UNIT("seconds",       1.*UNIT_TIME,                0,0,1,0) COMMA

  POST_UNIT("min",           60.*UNIT_TIME,               0,0,1,0) COMMA
  POST_UNIT("minutes",       60.*UNIT_TIME,               0,0,1,0) COMMA

  POST_UNIT("hr",            3600.*UNIT_TIME,             0,0,1,0) COMMA
  POST_UNIT("hours",         3600.*UNIT_TIME,             0,0,1,0) COMMA

  POST_UNIT("day",           86400.*UNIT_TIME,             0,0,1,0) COMMA
  POST_UNIT("days",          86400.*UNIT_TIME,             0,0,1,0) COMMA

  POST_UNIT("yr",            365.*86400.*UNIT_TIME,        0,0,1,0) COMMA
  POST_UNIT("years",         365.*86400.*UNIT_TIME,        0,0,1,0) COMMA

/* TEMP */

  POST_UNIT(DEG_K,           1.*UNIT_TEMP,                0,0,0,1) COMMA
  POST_UNIT(DEG_C,           1.*UNIT_TEMP,                0,0,0,1) COMMA
  POST_UNIT(DEG_F,           1.8*UNIT_TEMP,               0,0,0,1) COMMA

/* ANGLE */

  //  POST_UNIT("rad",           1.*UNIT_ANGLE,               0,0,0,0,1) COMMA
  //  POST_UNIT("radians",       1.*UNIT_ANGLE,               0,0,0,0,1) COMMA

  //  POST_UNIT("deg",           .01745329252*UNIT_ANGLE,     0,0,0,0,1) COMMA
  //  POST_UNIT("degrees",       .01745329252*UNIT_ANGLE,     0,0,0,0,1) COMMA

/* NUMBER */

  POST_UNIT("mol",           1.,                          0,0,0,0) COMMA
  POST_UNIT("moles",         1.,                          0,0,0,0) COMMA

  POST_UNIT("deci",          .1,                          0,0,0,0) COMMA
  POST_UNIT("centi",         .01 ,                        0,0,0,0) COMMA
  POST_UNIT("milli",         .001,                        0,0,0,0) COMMA
  POST_UNIT("micro",         .000001,                     0,0,0,0) COMMA
  POST_UNIT("deka",          10,                          0,0,0,0) COMMA
  POST_UNIT("hecto",         100,                         0,0,0,0) COMMA
  POST_UNIT("kilo",          1000,                        0,0,0,0) COMMA
  POST_UNIT("mega",          1000000,                     0,0,0,0) COMMA

  POST_UNIT("%",             .01,                         0,0,0,0) COMMA

/* COMBINED UNITS */

/* AREA/PERMEABILITY */

  POST_UNIT("acre",          4046.873*UNIT_AREA,          2,0,0,0) COMMA
  POST_UNIT("acres",         4046.873*UNIT_AREA,          2,0,0,0) COMMA

  POST_UNIT("darcy",         9.8697e-13*UNIT_AREA,        2,0,0,0) COMMA
  POST_UNIT("darcies",       9.8697e-13*UNIT_AREA,        2,0,0,0) COMMA

  POST_UNIT("md",            9.8697e-16*UNIT_AREA,        2,0,0,0) COMMA
  POST_UNIT("millidarcies",  9.8697e-16*UNIT_AREA,        2,0,0,0) COMMA

/* VOLUME */

  POST_UNIT("l",             1.e-3*UNIT_VOLUME,           3,0,0,0) COMMA
  POST_UNIT("liters",        1.e-3*UNIT_VOLUME,           3,0,0,0) COMMA

  POST_UNIT("ml",            1.e-6*UNIT_VOLUME,           3,0,0,0) COMMA
  POST_UNIT("milliliters",   1.e-6*UNIT_VOLUME,           3,0,0,0) COMMA
  POST_UNIT("cc",            1.e-6*UNIT_VOLUME,           3,0,0,0) COMMA

  POST_UNIT("USgal",         3.785e-3*UNIT_VOLUME,        3,0,0,0) COMMA
  POST_UNIT("USgallons",     3.785e-3*UNIT_VOLUME,        3,0,0,0) COMMA

  POST_UNIT("UKgal",         4.546e-3*UNIT_VOLUME,        3,0,0,0) COMMA
  POST_UNIT("UKgallons",     4.546e-3*UNIT_VOLUME,        3,0,0,0) COMMA

  POST_UNIT("bbl",           .1589873*UNIT_VOLUME,        3,0,0,0) COMMA

/* PRESSURE */

  POST_UNIT("Pa",            1.*UNIT_PRESSURE,            -1,1,-2,0) COMMA
  POST_UNIT("Pascals",       1.*UNIT_PRESSURE,            -1,1,-2,0) COMMA

  POST_UNIT("bar",           1.e+5*UNIT_PRESSURE,         -1,1,-2,0) COMMA
  POST_UNIT("bars",          1.e+5*UNIT_PRESSURE,         -1,1,-2,0) COMMA

  POST_UNIT("torr",          133.322*UNIT_PRESSURE,       -1,1,-2,0) COMMA
  POST_UNIT("torrs",         133.322*UNIT_PRESSURE,       -1,1,-2,0) COMMA

  POST_UNIT("mmHg",          133.322*UNIT_PRESSURE,       -1,1,-2,0) COMMA
  POST_UNIT("mmWater",       9.8*UNIT_PRESSURE,           -1,1,-2,0) COMMA

  POST_UNIT("psi",           6.894757e+3*UNIT_PRESSURE,   -1,1,-2,0) COMMA
  POST_UNIT("psia",          6.894757e+3*UNIT_PRESSURE,   -1,1,-2,0) COMMA

  POST_UNIT("atm",           1.01325e+5*UNIT_PRESSURE,    -1,1,-2,0) COMMA
  POST_UNIT("atmospheres",   1.01325e+5*UNIT_PRESSURE,    -1,1,-2,0) COMMA

  POST_UNIT("inHg",          3.386e+3*UNIT_PRESSURE,      -1,1,-2,0) COMMA
  POST_UNIT("inWater",       248.92*UNIT_PRESSURE,        -1,1,-2,0) COMMA

/* VISCOSITY */

  POST_UNIT("cp",            .001*UNIT_VISCOSITY,         -1,1,-1,0) COMMA
  POST_UNIT("centipoise",    .001*UNIT_VISCOSITY,         -1,1,-1,0) COMMA

  POST_UNIT("poise",         .1*UNIT_VISCOSITY,           -1,1,-1,0) COMMA

/* LIST TERMINATOR */

  POST_UNIT("END",           -1,                          0,0,0,0)

#undef COMMA
#undef POST_UNIT
