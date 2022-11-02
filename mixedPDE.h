#ifndef __mixedpde_h_included__
#define __mixedpde_h_included__

//#include <string>
//#include <cmath>

#include "Mesh/baseObjects.h"
#include "Mesh/polyMesh.h"
#include "directMixed.h"
#include "parameterData.h"
#include "Utilities/monitor.h"
#include "debug.h"
#include "lapacke.h"
#include <cblas.h>



//=============================================================================
//  EQUATIONS
//
//    - div D grad p + div(b p) + a p = f
//                       in Omega (bounded, connected, Lipschitz domain in R^2)
//    p = 0              on the boundary of Omega
//    where D is a tensor (matrix)
//
//    Setting u = - D grad p + b p,
//    find (u, p, l) in V^s_r x W_s x L_r such that
//
//      (D^{-1} u , v) - (p , div v) + (D^{-1} b p, v)
//                     + Sum_E (l , v.nu)_{delta E \ delta Omega} 
//                     + Sum_E (p_D , v.nu)_{delta E \Cap delta Omega} = 0 for all v in V^s_r
//
//      (div u, q) + (a p, q) = (f , q) for all q in W_s
//
//      Sum_E (u.nu , l)_{delta E} = 0 for all l in L_r
//
//    where 
//      V^s_r is the direct sum of V^s_r(E), which is 
//        a full / reduced direct mixed space on element E for s = r, r - 1
//      W_s is the direct sum of W_s(E) = \Po_s(E) restriced to each element
//      L_r is \Po_r(e) restriced to each interior edge e, 0 on boundary edges
//      
//  ELEMENTS E
//    E is a convex, nondegenerate polygon (of at least 3 sides)
//    Polynomials of degree r >= 0
//    The mesh is hybrid
//
//=============================================================================

class MixedPDE {
private:
  ParameterData* my_param;

public:
  MixedPDE(ParameterData& param_in) : my_param(&param_in) {};

  ParameterData* parameterDataPtr() const { return my_param; };
  
  int solve_conf(Monitor& monitor);

  int solve_hybrid(Monitor& monitor);

  int solve(Monitor& monitor) {
    if (my_param->conforming) {
      return solve_conf(monitor);
    } else {
      return solve_hybrid(monitor);
    }
  };
};

lapack_int mat_inv(double *A, int n);
lapack_int block_mat_inv(double *A, int n, int *m, int num);
#endif
