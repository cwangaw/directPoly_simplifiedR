/******************************************************************************
Evaluate expressions

This package evaluates expressions that may contain units information.
  It consists of the three files  expr.c , expr.h , and units.h .

Todd Arbogast
The University of Texas at Austin
July 1997
July 2005
******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cstring>
#include <math.h>
#include <stdlib.h>
#include <cerrno>
using namespace std;

#include "units.h"
#include "expr.h"

/* OP CODES */
#define OP_BEGIN    -7
#define OP_LPAREN   -6

#define OP_EXPONEN  -5
#define OP_MULT     -4
#define OP_DIV      -3
#define OP_ADD      -2
#define OP_SUB      -1

#define OP_NUMBER    1

#define OP_RPAREN    2
#define OP_END       3

/* DECLARATIONS AND STATIC VARS */
static char baseDegrees[] = UNIT_ABSOLUTE_TEMP;
static int baseTempIndex = UNIT_INDEX_ABSOLUTE_TEMP;
static int convertedToTempIndex = -1;
static int convertedFromTempIndex = -1;

/******************************************************************************
Title: expr_setBaseUnits
Description: Set or reset the base units.
*/
int expr_setBaseUnits(double length, double mass, double time,
		      double temp, int tempIndex) {
  int i = 0;
  int j,k;
  double baseUnit[5];

  baseUnit[0] = length;
  baseUnit[1] = mass;
  baseUnit[2] = time;
  baseUnit[3] = temp;
  //baseUnit[4] = angle;

  while(strcmp(unitName[i],"END")) {
    for(j=0; j<NUMBER_DIMEN; j++) {
      if(baseUnit[j] <= 0) return EXPR_NEGATIVE;
      if(unitDimen[i][j] < 0) {
	for(k=0; k<-unitDimen[i][j]; k++) unitValue[i] *= baseUnit[j];
      } else {
	for(k=0; k<unitDimen[i][j]; k++) unitValue[i] /= baseUnit[j];
      }
    }
    i++;
  }

  switch(tempIndex) {
  case INDEX_DEG_K:
    strcpy(baseDegrees,DEG_K);
    baseTempIndex = INDEX_DEG_K;
    break;
  case INDEX_DEG_C:
    strcpy(baseDegrees,DEG_C);
    baseTempIndex = INDEX_DEG_C;
    break;
  case INDEX_DEG_F:
    strcpy(baseDegrees,DEG_F);
    baseTempIndex = INDEX_DEG_F;
    break;
  default:
    return EXPR_UNKNOWN_UNIT;
  }
  return 0;
}

/******************************************************************************
Title: expr_unitsCount
Description: Count the number of units in the units list.
Return: The count.
*/
int expr_unitsCount() {
  int numberOfUnits=-1;
  do { numberOfUnits++; } while(strcmp(unitName[numberOfUnits],"END"));
  return numberOfUnits;
}

/******************************************************************************
Title: expr_numberValue
Description: Determine the value of a numeric string.
Input:  str           The string containing the number.
Output  value         Its parsed value.
        nRead         Length of the unit name.
Return: Error flag.
*/
int expr_numberValue(char* str, double* value, int* nRead) {
  char* endPointer[1];
  endPointer[0] = str + strlen(str);

  errno = 0;
  *value = strtod(str,endPointer);
  if(errno) return EXPR_IMPROPER;

  *nRead = (int) (endPointer[0] - str);
  return 0;
}

/******************************************************************************
Title: expr_unitNumber
Description: Compare the string to the unit name list and return the index
             of the longest matching unit name.
Input:  str           The string containing the unit name.
Output  unitNumber    The index of the longest matching unit name.
        nRead         Length of the unit name.
Return: Error flag.
*/
int expr_unitNumber(char* str, int* unitNumber, int* nRead) {
  int i;
  int len;
  int unitLen = 0;
  int unitNo = -1;

  i = 0;
  while(strcmp(unitName[i],"END")) {
    len = strlen(unitName[i]);
    if(!strncmp(str,unitName[i],len)) {
      if(len > unitLen) {
	unitLen = len;
	unitNo = i;
      }
    }
    i++;
  }
  if(unitNo == -1) return EXPR_UNKNOWN_UNIT;

  *nRead = unitLen;
  *unitNumber = unitNo;
  return 0;
}

/******************************************************************************
Title: expr_listUnits
Description: List the currently available units
Input:  os   Output stream
*/
void expr_listUnits(ostream& os) {
  int i,j;

  os << "DEFINED UNITS\n";

  i = 0;
  while(strcmp(unitName[i],"END")) {
    os << "  " << setw(15) << left << unitName[i] << " = ";
    os << setw(10) << setprecision(4) << left << unitValue[i] << " [ ";
    for(j=0; j<NUMBER_DIMEN-1; j++)
      os << setw(2) << right << unitDimen[i][j] << ", ";
    os << setw(2) << right << unitDimen[i][NUMBER_DIMEN-1] << " ]\n";
    i++;
  }
}

/******************************************************************************
Title: expr_unitsIterator
Description: List the currently available units
Input:  resetOnly  Reset the iterator (and return immediately)
Output: name       Pointer to the name of the unit
        value      The value of the unit
        dimen      The NUMBER_DIMEN dimensions of the unit
Return: Done flag.
*/
int expr_unitsIterator(char** name, double& value, int* dimen, int resetOnly) {
  static int i=0;
  int j;

  if(resetOnly) {
    i=0;
  } else {
    *name = unitName[i];
    value = unitValue[i];
    for(j=0; j<NUMBER_DIMEN-1; j++) dimen[j] = unitDimen[i][j];
    i++;
  }
  return (strcmp(unitName[i],"END")) ? 0 : 1;
}

/******************************************************************************
Title: expr_baseTempIndex()
Description: Return the index of the base units.
Return: index of the base units.
*/
int expr_baseTempIndex() {
  return baseTempIndex;
}

/******************************************************************************
Title: expr_convertedFromTempIndex()
Description: Return the index of the input temp of the last temperature
             conversion.
Return: index of the input temp.
*/
int expr_convertedFromTempIndex() {
  return convertedFromTempIndex;
}

/******************************************************************************
Title: expr_convertedToTempIndex()
Description: Return the index of the output temp of the last temperature
             conversion.
Return: index of the output temp.
*/
int expr_convertedToTempIndex() {
  return convertedToTempIndex;
}

/******************************************************************************
Title: expr_convertTempBase
Description: Converts one absolute temp to base units
Input:  temp     The temperature.
        fromDeg  The degree scale on input (see dimen.h).
Output: temp     The converted temperature.
Return: Error flag.
*/
int expr_convertTempBase(double& temp, char* fromDeg) {
  return expr_convertTemp(temp, fromDeg, baseDegrees);
}

/******************************************************************************
Title: expr_convertBaseTemp
Description: Converts one absolute temp in base units to other units
Input:  temp     The temperature.
	toDeg    The degree scale on output (see dimen.h).
Output: temp     The converted temperature.
Return: Error flag.
*/
int expr_convertBaseTemp(double& temp, char* toDeg) {
  return expr_convertTemp(temp, baseDegrees, toDeg);
}

/******************************************************************************
Title: expr_convertTemp
Description: Converts one absolute temp to another
Input:  temp     The temperature.
        fromDeg  The degree scale on input (see dimen.h).
	toDeg    The degree scale on output (see dimen.h).
Output: temp     The converted temperature.
Return: Error flag.
*/
int expr_convertTemp(double& temp, char* fromDeg, char* toDeg) {

  /* convert from "fromDeg" to degrees C */
  if(!strcmp(fromDeg,DEG_F)) {
    temp = 5*(temp - 32)/9;
    convertedFromTempIndex = INDEX_DEG_F;
  } else if(!strcmp(fromDeg,DEG_K)) {
    temp -= 273.15;
    convertedFromTempIndex = INDEX_DEG_K;
  } else if(!strcmp(fromDeg,DEG_C)) {
    convertedFromTempIndex = INDEX_DEG_C;
  } else {
    return EXPR_UNKNOWN_UNIT;
  }

  /* convert to "toDeg" from degrees C */
  if(!strcmp(toDeg,DEG_F)) {
    temp = 1.8*temp + 32;
    convertedToTempIndex = INDEX_DEG_F;
  } else if(!strcmp(toDeg,DEG_K)) {
    temp += 273.15;
    convertedToTempIndex = INDEX_DEG_K;
  } else if(!strcmp(toDeg,DEG_C)) {
    convertedToTempIndex = INDEX_DEG_C;
  } else {
    return EXPR_UNKNOWN_UNIT;
  }

  return 0;
}

/******************************************************************************
Title: expr_stringValue
Description: Evaluates a string's value
Input:  str      The string (number or unit) to evaluate.
Output  value    The string's value.
        dimen    The string's dimensions.  Units from the units table
                 have dimensions; numbers do not.
        nRead    The number of characters read.
Return: Error flag.
*/
int expr_stringValue(char* str, double* value, int* dimen, int* nRead) {
  int i;
  int error = 0;
  int unitNumber = 0;

  if( ('0' <= str[0] && str[0] <= '9') || str[0] == '.' ) {

    /* ordinary number */
    error = expr_numberValue(&str[0],value,nRead);
    for(i=0; i< NUMBER_DIMEN; i++) dimen[i] = 0;

  } else {

    /* unit */
    error = expr_unitNumber(&str[0],&unitNumber,nRead);
    *value = unitValue[unitNumber];
    for(i=0; i< NUMBER_DIMEN; i++) dimen[i] = unitDimen[unitNumber][i];
  }

  return error;
}

/******************************************************************************
Title: expr_evalSubOps
Description: Evaluates a single parenthetical set of ops.
INPUT/OUTPUT:
  op        The operation/operand stack.
  value     The operand values.
  dimen     The operand dimensions.
OUTPUT:
  dimensFailed   Flag set to 1 if dimensions failed to be consistent.
                 Note incoming value unchanged otherwise, so multiple
                 calls can be made.
  additionFound  Flag set to 1 if addition/subtraction was found.
                 Note incoming value unchanged otherwise, so multiple
                 calls can be made.
Return: Error flag.
*/
static int expr_evalSubOps(int* op, double* value, int dimen[][NUMBER_DIMEN],
			   int* dimensFailed, int* additionFound) {
  int i;
  int curr,prev,next;
  double sign = 1;
  double dval1,dval2;
  int ival2;
  int unary;

  /* EXPONEN */
  curr = 0;
  while(op[curr] < OP_RPAREN) {
    if(op[curr] == OP_EXPONEN) {
      next = curr + 1; while(!op[next]) next++;
      prev = curr - 1; while(!op[prev]) prev--;

      /* negation */
      sign = 1;
      if(op[next] == OP_ADD || op[next] == OP_SUB) {
        if(op[next] == OP_SUB) sign = -1;
        op[next] = 0;
        next = next + 1; while(!op[next]) next++;
      }

      /* errors and signs */
      if(op[next] != OP_NUMBER || op[prev] != OP_NUMBER)
        return EXPR_IMPROPER;
      dval2 = sign*value[next];
      ival2 = (int) dval2;
      if(ival2 != dval2) return EXPR_IMPROPER;

      /* apply operation */
      if(value[prev] == 0) {
        value[curr] = 0;
      } else {
        errno = 0;
        value[curr] = pow(value[prev],dval2);
        if(errno) return EXPR_ARITHMETIC_OP;
      }
      for(i=0; i<NUMBER_DIMEN; i++) {
        if(dimen[next][i] != 0) *dimensFailed = 1;
        dimen[curr][i] = dimen[prev][i]*ival2;
      }
      op[prev] = 0;
      op[next] = 0;
      op[curr] = OP_NUMBER;

      curr = next;
    }
    curr++;
  }

  /* MULT and DIV */
  curr = 0;
  while(op[curr] < OP_RPAREN) {
    if(op[curr] == OP_MULT || op[curr] == OP_DIV) {
      next = curr + 1; while(!op[next]) next++;
      prev = curr - 1; while(!op[prev]) prev--;

      /* negation */
      sign = 1;
      if(op[next] == OP_ADD || op[next] == OP_SUB) {
        if(op[next] == OP_SUB) sign = -1;
        op[next] = 0;
        next = next + 1; while(!op[next]) next++;
      }

      /* errors and signs */
      if(op[next] != OP_NUMBER || op[prev] != OP_NUMBER)
        return EXPR_IMPROPER;
      dval2 = sign*value[next];

      /* apply operation */
      if(op[curr] == OP_MULT) {
        value[curr] = value[prev]*dval2;
        for(i=0; i<NUMBER_DIMEN; i++)
          dimen[curr][i] = dimen[prev][i] + dimen[next][i];
      } else {
        if(dval2 == 0) return EXPR_ARITHMETIC_OP;
        value[curr] = value[prev]/dval2;
        for(i=0; i<NUMBER_DIMEN; i++)
          dimen[curr][i] = dimen[prev][i] - dimen[next][i];
      }
      op[prev] = 0;
      op[next] = 0;
      op[curr] = OP_NUMBER;

      curr = next;
    }
    curr++;
  }

  /* ADD and SUB */
  curr = 0;
  while(op[curr] < OP_RPAREN) {
   if(op[curr] == OP_ADD || op[curr] == OP_SUB) {
      next = curr + 1; while(!op[next]) next++;
      prev = curr - 1; while(!op[prev]) prev--;

      /* errors */
      unary = (op[prev] != OP_NUMBER);
      if(op[next] != OP_NUMBER) return EXPR_IMPROPER;

      /* apply operation */
      dval1 = unary ? 0 : value[prev];
      if(op[curr] == OP_ADD) {
        value[next] = dval1 + value[next];
      } else {
        value[next] = dval1 - value[next];
      }
      if(!unary) {
	*additionFound = 1;
	for(i=0; i<NUMBER_DIMEN; i++)
	  if(dimen[prev][i] != dimen[next][i]) *dimensFailed = 1;
      }
      op[prev] = 0;
      op[curr] = 0;

      curr = next;
    }
    curr++;
  }

  /* TEST */
  curr = 0;
  prev = 0;
  while(op[curr] < OP_RPAREN) {
    if( !(op[curr] == 0 || (!prev && op[curr] == OP_NUMBER)) ) {
      /*
      cerr << "ERROR IN EVALSUBOPS:\n";
      i = 0;
      while(op[i] < OP_RPAREN) { cerr << " " << op[i]; i++; }
      cerr << endl;
      */
      return EXPR_IMPROPER;
    }
    if(op[curr] == OP_NUMBER) prev = 1;
    curr++;
  }

  return 0;
}

/******************************************************************************
Title: expr_eval
Description: Evaluates a units expression.
Input:
  expr           The expression to evaluate.
                  We allow * / + - of numbers or unit names, as well as
                  ( ) groupings and ** or ^ INTEGER powers
  mode            Flag:
                   0 = ordinary expression
                   1 = units expression. The expression is the units of
                       some other unspecified number.  Leading '*' and
                       '/' are allowed.  Binary '+' and '-' are not
		       allowed (though unary '+' and '-' are allowed).
                       Nonpositive results are not allowed.
		   2 = ordinary expression with leading '*' and '/' allowed.
  dimenCheck      Flag:  0 = ignore dimensions (do not check or give).
                         1 = check the dimensions
                             (and do NOT overwrite dimension).
                        -1 = give the dimensions found.
  dimensions      The dimensions that should result (if dimenCheck = 1).
Output:
  result          The value of the expression.
  dimensions      The dimensions of the result (if dimenCheck = -1).
Return: Error flag.
Remark: op is either an operator or an operand.
*/

int expr_eval(char* expr, double& result,
	      int mode, int dimenCheck, int* dimensions) {
  int i,j;
  int curr,next;
  int nRead;
  int error = 0;

  int additionFound = 0;
  int dimensFailed = 0;

  int op[MAX_NUMBER_OPS];
  double value[MAX_NUMBER_OPS];
  int dimen[MAX_NUMBER_OPS][NUMBER_DIMEN];

  /* Parse expr to set up op stack */

  op[0] = OP_BEGIN;
  i = 0;
  curr = 1;
  if(mode == 1 || mode == 2) {
    while(expr[i] == ' ' || expr[i] == '\t' || expr[i] == '\n') i++;
    if(expr[i] == '*' || expr[i] == '/') {
      op[1] = OP_NUMBER;
      value[1] = 1;
      for(j=0; j<NUMBER_DIMEN; j++) dimen[1][j] = 0;
      op[2] = (expr[i] == '*') ? OP_MULT : OP_DIV;
      i++;
      curr = 3;
    }
  }
  for(; curr < MAX_NUMBER_OPS; curr++) {
    while(expr[i] == ' ' || expr[i] == '\t' || expr[i] == '\n') i++;
    switch(expr[i]) {
    case '\0':
      op[curr] = OP_END;
      break;
    case '(':
      op[curr] = OP_LPAREN;
      break;
    case ')':
      op[curr] = OP_RPAREN;
      break;
    case '^':
      op[curr] = OP_EXPONEN;
      break;
    case '*': 
      if(expr[i+1] == '*') {
        op[curr] = OP_EXPONEN;
        i++;
      } else {
        op[curr] = OP_MULT;
      }
      break;
    case '/':
      op[curr] = OP_DIV;
      break;
    case '+':
      op[curr] = OP_ADD;
      break;
    case '-':
      op[curr] = OP_SUB;
      break;
    default:
      op[curr] = OP_NUMBER;
      error = expr_stringValue(&expr[i],&value[curr],dimen[curr],&nRead);
      if(error) return error;
      i += nRead-1;
    }
    i++;
    if(op[curr] == OP_END) break;
  }
  if(curr >= MAX_NUMBER_OPS) return EXPR_TOO_BIG;

  /* Collapse op stack */

  curr = 0;
  while(op[curr] != OP_END) {
    if(op[curr] == OP_LPAREN) {
      next = curr + 1; while(!op[next]) next++;
      if(op[next] == OP_RPAREN) {
        op[curr] = 0;
        op[next] = 0;
        curr = next;
      }
    } else if(op[curr] == OP_ADD || op[curr] == OP_SUB) {
      next = curr + 1; while(!op[next]) next++;
      if(op[next] == OP_ADD || op[next] == OP_SUB) {
        if(op[curr] == OP_SUB)
          op[next] = (op[next] == OP_ADD) ? OP_SUB : OP_ADD;
        op[curr] = 0;
        curr = next-1;
      }
    }
    curr++;
  }

  /* Evaluate op stack */

  while(1) {
    curr = 0;
    while(op[curr] != OP_END && op[curr] != OP_LPAREN) curr++;
    if(op[curr] == OP_END) break;

    next = 0;
    while(op[next] != OP_END && op[next] != OP_RPAREN) {
      if(op[next] == OP_LPAREN) curr = next;
      next++;
    }
    if(op[next] == OP_END) return EXPR_NONNESTED_PARENTHESES;

    error = expr_evalSubOps(&op[curr+1],&value[curr+1],&dimen[curr+1],
                            &dimensFailed,&additionFound);
    if(error) return error;
    op[curr] = 0;
    op[next] = 0;
  }

  curr = 0;
  while(op[curr] != OP_END && op[curr] != OP_RPAREN) curr++;
  if(op[curr] != OP_END) return EXPR_NONNESTED_PARENTHESES;

  curr = 1;
  while(op[curr] != OP_END) {
    if(op[curr]) { curr = -1; break; }
    curr++;
  }
  if(curr != -1) return EXPR_IMPROPER;

  error = expr_evalSubOps(&op[1],&value[1],&dimen[1],
			  &dimensFailed,&additionFound);
  if(error) return error;

  curr = 0;
  while(op[curr] != OP_END && op[curr] != OP_NUMBER) curr++;
  result = value[curr];

  if(dimenCheck && dimensFailed) return EXPR_DIMENS_FAILED;
  if(dimenCheck > 0) {
    if(mode == 1 && result <= 0) return EXPR_NEGATIVE;
    for(i=0; i<NUMBER_DIMEN; i++)
      if(dimensions[i] != dimen[curr][i]) return EXPR_DIMENS_FAILED;
  } else if(dimenCheck < 0) {
    for(i=0; i<NUMBER_DIMEN; i++) dimensions[i] = dimen[curr][i];
  }
  if(mode == 1 && additionFound) return EXPR_ADDITION_FOUND;
  return 0;
}

/******************************************************************************
TEST PROGRAM
*/
/*
int main(int argc, char** argv) {
  double result;
  int dimensions[NUMBER_DIMEN] = {0,0,0,0,0};

  expr_listUnits(cout);

  for(int i = 1; i < argc; i++) {
    int error = expr_eval(argv[i],result,2,2,dimensions);

    cout << argv[i] << " = " << result << " (";
    for(int j=0; j<NUMBER_DIMEN; j++) cout << " " << dimensions[j];

    cout << " )" << endl;
    if(error) cout << "  (ERROR CODE: " << error << ")" << endl;
  }

  return 0;
}
*/
