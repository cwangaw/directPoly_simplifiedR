#ifndef _EXPR_H
#define _EXPR_H

#include <iostream>

#define MAX_NUMBER_OPS  42

/* ERROR AND RETURN CODES */
#define EXPR_NEGATIVE               1 /* not fatal */
#define EXPR_ADDITION_FOUND         2 /* not fatal */
#define EXPR_DIMENS_FAILED          3 /* not fatal */
#define EXPR_TOO_BIG                4
#define EXPR_NONNESTED_PARENTHESES  5
#define EXPR_ARITHMETIC_OP          6
#define EXPR_IMPROPER               7
#define EXPR_UNKNOWN_UNIT           8

int expr_setBaseUnits(double length, double mass, double time,
		      double temp, int baseTempIndex);
void expr_listUnits(std::ostream& os);
int expr_unitsIterator(char** name, double& value, int* dimen, int resetOnly);

int expr_eval(char* expr, double& result,
	      int mode, int dimenCheck, int* dimensions);

int expr_convertTemp(double& temp, char* fromDeg, char* toDeg);
int expr_convertTempBase(double& temp, char* fromDeg);
int expr_convertBaseTemp(double& temp, char* toDeg);
int expr_baseTempIndex();
int expr_convertedFromTempIndex();
int expr_convertedToTempIndex();

#endif
