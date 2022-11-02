#ifndef __ellipticpde_h_included__
#define __ellipticpde_h_included__

//#include <string>
//#include <cmath>

#include "Mesh/baseObjects.h"
#include "Mesh/polyMesh.h"
#include "directSerendipity.h"
#include "parameterData.h"
#include "Utilities/monitor.h"
#include "debug.h"


//=============================================================================
//  EQUATIONS
//
//    - div D grad p + div(b p) + c.grad p + a p = f
//                       in Omega (bounded, connected, Lipschitz domain in R^2)
//    p = g              on the boundary of Omega
//    where D is a tensor (matrix), b and c are vectors.
//
//    Find p in DS_{r,g} such that
//      (D grad p , grad q) + (b p , grad q) + (c.grad p , q) + (a p , q)
//          = (f , q) for all q in DS_{r,0}
//    where DS_r are the direct serendipity spaces and DS_{r,g} has the
//      boundary nodes set to g. We use a nodal basis for the space.
//      
//  ELEMENTS E
//    E is a convex, nondegenerate polygon (of at least 3 sides)
//    Polynomials of degree r >= 1
//    The mesh is conforrming
//
//=============================================================================

class EllipticPDE {
private:
  ParameterData* my_param;


public:
  EllipticPDE(ParameterData& param_in) : my_param(&param_in) {};

  ParameterData* parameterDataPtr() const { return my_param; };
  
  int solve(Monitor& monitor);
};

#endif
