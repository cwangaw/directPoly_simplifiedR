#include <cmath>
using namespace std;

#include "debug.h"
#include "fcns.h"
#include <assert.h>

////////////////////////////////////////////////////////////////////////////////
// SPATIALLY DEPENDENT PDE FUNCTIONS

// Coefficients
static const double PI = 3.141592653589793238462643383279502884197169399375105820974944;

double coefA(double x, double y) {
  return 0;
}

void coefB(double x, double y, Tensor1& b) {
  b.set(0,0);
}

void coefC(double x, double y, Tensor1& c) {
  c.set(0,0);
}

void coefD(double x, double y, Tensor2& d) {
  //d.set(0,0,0,0);
  d.set(1,0,0,1);
}

void coefD_inv(double x, double y, Tensor2& d) {
  //d.set(0,0,0,0);
  d.set(1,0,0,1);
}

Tensor2 trueHessianSoln(double x, double y);
// Source f
double sourceVal(double x, double y) {
  //Id D=I, a=b=0
  Tensor2 h = trueHessianSoln(x,y);
  return - h(1,1) - h(2,2);

}

// BC values g
double bcVal(double x, double y) {
  return trueSoln(x,y);
}

////////////////////////////////////////////////////////////////////////////////
// TRUE SOLUTION (IF KNOWN)

bool trueSolnKnown() { return true; }

// Real solution
double trueSoln(double x, double y) {
  return sin(PI*x)*sin(PI*y);
}

Tensor1 trueGradSoln(double x, double y) {
  return Tensor1(PI*cos(PI*x)*sin(PI*y),PI*sin(PI*x)*cos(PI*y));
}

Tensor2 trueHessianSoln(double x, double y) {
  return Tensor2(-PI*PI*sin(PI*x)*sin(PI*y),PI*PI*cos(PI*x)*cos(PI*y),
    PI*PI*cos(PI*x)*cos(PI*y),-PI*PI*sin(PI*x)*sin(PI*y));
}

Tensor1 trueUSoln(double x, double y) {
  Tensor1 b;
  coefB(x,y,b);
  Tensor2 d;
  coefD(x,y,d);
  
  return -1*(d*trueGradSoln(x,y)) + trueSoln(x,y)*b;
}

double trueDivUSoln(double x, double y) {
  return sourceVal(x,y)-coefA(x,y)*trueSoln(x,y);
}