#ifndef __fcns_h_included__
#define __fcns_h_included__

#include "Mesh/baseObjects.h"

////////////////////////////////////////////////////////////////////////////////
// SPATIALLY DEPENDENT PDE FUNCTIONS

double coefA(double x, double y);
void coefB(double x, double y, Tensor1& b);
void coefC(double x, double y, Tensor1& c);
void coefD(double x, double y, Tensor2& d);
void coefD_inv(double x, double y, Tensor2& d);

double sourceVal(double x, double y);
double bcVal(double x, double y);

////////////////////////////////////////////////////////////////////////////////
// TRUE SOLUTION (IF KNOWN)

bool trueSolnKnown();
double trueSoln(double x, double y);
Tensor1 trueGradSoln(double x, double y);
Tensor2 trueHessianSoln(double x, double y);
Tensor1 trueUSoln(double x, double y);
double trueDivUSoln(double x, double y);

#endif
