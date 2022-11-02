#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <assert.h>
#include "mixedPDE.h"
#include "fcns.h"
#include "Mesh/baseObjects.h"
#include "Mesh/polyMesh.h"
#include "directMixed.h"
#include "parameterData.h"
#include "polyQuadrature.h"
#include "Utilities/monitor.h"
#include "Utilities/debug.h"
#include <complex.h>
#include <algorithm>

using namespace std;
using namespace directserendipity;
using namespace polymesh;
using namespace polyquadrature;


int MixedPDE::solve_conf(Monitor& monitor) {
  monitor(0,"Polynomial Degree = ", parameterDataPtr()->dmSpace.degPolyn());

  ParameterData& param = *parameterDataPtr();
  DirectMixedConf dmSpace = DirectMixedConf(parameterDataPtr()->dmSpace);
  // TEST SINGLE ELEMENT MESH //////////////////////////////////////////////

  if(false) {
    monitor(0,"Test Single Element Mesh and Space");

    polymesh::PolyMesh s_mesh(parameterDataPtr()->mesh.elementPtr(3));
    std::string fileName = parameterDataPtr()->directory_name;
    fileName += "mesh_e4";
    s_mesh.write_matlab(fileName);

    DirectMixed s_dmSpace(8,&s_mesh,false);
    fileName = parameterDataPtr()->directory_name;
    fileName += "dmSpace_e4";
    s_dmSpace.write_matlab(fileName);
  }
  
  // TEST BASIS FUNCTIONS //////////////////////////////////////////////////

  if(true) {
    monitor(0,"\nTest basis functions\n");

    DirectMixedConfArray u(&dmSpace, 'f');
    DirectDGArray p(&(parameterDataPtr()->dmSpace), 'f');

    // Set the coefficients of all the mixed conforming basis functions to be zero
    for(int i=0; i<u.size(); i++) {
      u[i]=0;
    }
    // Set the coefficient of the last mixed conforming basis function to be one
    u[u.size()-1]=1;

    // Set the coefficients of all the discrete Galerkin basis functions to be zero
    for(int i=0; i<p.size(); i++) {
      p[i]=0;
    }
    // Set the coefficients of the first discrete Galerkin basis function to be one
    p[0]=1;

    monitor(1,"Write Array");

    std::string fileName = parameterDataPtr()->directory_name;
    fileName += "basis_mesh_mixed";
    u.write_matlab_mesh(fileName,51,51);

    fileName = parameterDataPtr()->directory_name;
    fileName += "basis_mesh_mixed.txt";
    dmSpace.write_raw(fileName);

    fileName = parameterDataPtr()->directory_name;
    fileName += "basis_mesh_DG";
    p.write_matlab_mesh(fileName,51,51);
  }

  // TEST QUADRATURE ///////////////////////////////////////////////////////

  if(false) {
    monitor(1,"Test Quadrature Rules");
    testPolyQuadrature(&(parameterDataPtr()->mesh),1e-6);
    return 0;
  }
  
  // SOLVE THE PDE ///////////////////////////////////////////////////////////
if (true) {
  monitor(0,"\nSolve the PDE\n");

  // Only when r > 0, we calculate reduced space

  DirectMixedConfArray solution_u_f(&(dmSpace),'f');
  DirectMixedConfArray solution_u_r(&(dmSpace),'r');
  DirectDGArray solution_p_f(&(dmSpace),'f');
  DirectDGArray solution_p_r(&(dmSpace),'r');

  // Initialize matrix A for both full and reduced space
  int dimAfull = dmSpace.nMixedDoFs('f');
  int dimAreduced = dmSpace.nMixedDoFs('r');

  std::vector<double> mat_A_full_vector(dimAfull*dimAfull,0);
  double* mat_A_full = mat_A_full_vector.data();
  std::vector<double> mat_A_reduced_vector(dimAreduced*dimAreduced,0);
  double* mat_A_reduced = mat_A_reduced_vector.data();

  // Initialize matrix B for both full and reduced space
  int rowBfull = dimAfull, rowBreduced = dimAreduced;
  int colBfull = dmSpace.nDGDoFs('f');
  int colBreduced = dmSpace.nDGDoFs('r');

  std::vector<double> mat_B_full_vector(rowBfull*colBfull,0);
  double* mat_B_full = mat_B_full_vector.data();
  std::vector<double> mat_B_reduced_vector(rowBreduced*colBreduced,0);
  double* mat_B_reduced = mat_B_reduced_vector.data();

  // Initialize right hand side (Dirichlet BC)
  std::vector<double> rhs_bc_full_vector(dimAfull,0);
  double* rhs_bc_full = rhs_bc_full_vector.data();   

  std::vector<double> rhs_bc_reduced_vector(dimAreduced,0);
  double* rhs_bc_reduced = rhs_bc_reduced_vector.data();   

  // Initialize right hand side (W_s coefficients only)
  std::vector<double> rhs_full_vector(colBfull,0);
  double* rhs_full = rhs_full_vector.data();

  std::vector<double> rhs_reduced_vector(colBreduced,0);
  double* rhs_reduced = rhs_reduced_vector.data();

  // quadrature points
  // Note that we update quadRule in each iElement loop
  // But we define quadEdgeRule for all edges at once

  polyquadrature::PolyQuadrature quadRule(13,param.refinement_level);
  std::vector<polyquadrature::PolyEdgeQuadrature> quadEdgeRule(dmSpace.nEdges());

  for (int i = 0; i < dmSpace.nEdges(); i++) {
    if (dmSpace.bcType(i) == EdgeBCType::boundary) {
      quadEdgeRule[i].set(13, dmSpace.edgePtr(i));
    }
  }

  monitor(1,"Matrix and RHS Assembly"); ////////////////////////////////////////


  int starting_colBfull = 0, starting_colBreduced = 0;

  int loc_dimAfull = 0, loc_dimAreduced = 0, loc_colBfull = 0, loc_colBreduced = 0;

  int curr_full_index, curr_reduced_index;
  double evaluation;
  int numEdges = 0, edge_loc_to_glob;


  for(int iElement=0; iElement<param.mesh.nElements(); iElement++) {
    DirectMixedConfFE* mePtr = dmSpace.MixedElementPtr(iElement);
    DirectDGFE* dgePtr = dmSpace.DGElementPtr(iElement);

    quadRule.setElement(mePtr->elementPtr());
    
    mePtr->initBasis(quadRule.pts(), quadRule.num());
    dgePtr->initBasis(quadRule.pts(), quadRule.num());

    loc_dimAfull = mePtr -> dimVFull();
    loc_dimAreduced = mePtr -> dimVReduced();
    loc_colBfull = dgePtr -> dimFull();
    loc_colBreduced = dgePtr -> dimReduced();

    int nGon = param.mesh.elementPtr(iElement)->nVertices();
    int local_poly_dofs_full = (param.polynomial_degree + 2) * (param.polynomial_degree + 1) / 2 - 1;
    int local_poly_dofs_reduced = (param.polynomial_degree + 1) * param.polynomial_degree / 2 - 1;

    // Matrix A, B and rhs assembly over elements
 
    for(int iPt=0; iPt<quadRule.num(); iPt++) {
      double x = quadRule.pt(iPt)[0];
      double y = quadRule.pt(iPt)[1];

      Tensor2 valD;
      coefD_inv(x,y,valD);
  
      // Assemble matrix A
      for (int j = 0; j < loc_dimAfull; j++) {
        Tensor1 v_j = mePtr -> basis(j,iPt);
        if (mePtr -> isVertexDoF(j)) { v_j *= dmSpace.vertex_loc_to_glob_coeff(iElement,mePtr->mapDoFToVertex(j)); }
        for (int i = 0; i < loc_dimAfull; i++) {
          Tensor1 u_i = mePtr -> basis(i,iPt);
          if (mePtr -> isVertexDoF(i)) { u_i *= dmSpace.vertex_loc_to_glob_coeff(iElement,mePtr->mapDoFToVertex(i)); }
          curr_full_index = dimAfull * dmSpace.loc_to_glob(iElement,j,'f') + dmSpace.loc_to_glob(iElement,i,'f');
          evaluation = ((valD*u_i)*v_j)*quadRule.wt(iPt);
          mat_A_full[curr_full_index] += evaluation;
          if (j < loc_dimAreduced && i < loc_dimAreduced) {
            curr_reduced_index = dimAreduced * dmSpace.loc_to_glob(iElement,j,'r') + dmSpace.loc_to_glob(iElement,i,'r');
            mat_A_reduced[curr_reduced_index] += evaluation;
          }
        }
      }

      // Assemble matrix B

      for (int i = 0; i < loc_colBfull; i++ ) {
        double p_i = dgePtr -> basis(i,iPt);

        // First deal with vertex dofs ( they have nonzero divergence part )
        for (int j = 0; j < nGon; j++) {
          double divv_j = dmSpace.vertex_loc_to_glob_coeff(iElement,j) * mePtr -> vertexBasisDiv(j,iPt);
          evaluation = divv_j * p_i * quadRule.wt(iPt);

          curr_full_index = colBfull * dmSpace.vertex_loc_to_glob(iElement,j) + (starting_colBfull + i);
          mat_B_full[curr_full_index] += evaluation;

          if (i < loc_colBreduced) {
            curr_reduced_index = colBreduced * dmSpace.vertex_loc_to_glob(iElement,j) + (starting_colBreduced + i);
            mat_B_reduced[curr_reduced_index] += evaluation;
          }
        }

        // Now deal with poly dofs ( they have nonzero divergence part )
        for (int j = 0; j < local_poly_dofs_full; j++) {
          double divv_j = mePtr -> polyBasisDiv(j,iPt);
          evaluation = divv_j * p_i * quadRule.wt(iPt);

          curr_full_index = colBfull * dmSpace.poly_loc_to_glob(iElement,j,'f') + (starting_colBfull + i);
          mat_B_full[curr_full_index] += evaluation;

          if (j < local_poly_dofs_reduced && i < loc_colBreduced) {
            curr_reduced_index = colBreduced * dmSpace.poly_loc_to_glob(iElement,j,'r') + (starting_colBreduced + i);
            mat_B_reduced[curr_reduced_index] += evaluation;
          }
        }

      }

      // Assemble RHS (W_s coefficients only)
      double f_wted = sourceVal(quadRule.pt(iPt)[0],quadRule.pt(iPt)[1]) * quadRule.wt(iPt);
      for (int j = 0; j < loc_colBfull; j++) {
        curr_full_index = starting_colBfull + j;
        evaluation = f_wted * dgePtr->basis(j,iPt);
        rhs_full[curr_full_index] += evaluation;
        if (j < loc_colBreduced) {
          curr_reduced_index = starting_colBreduced + j;
          rhs_reduced[curr_reduced_index] += evaluation;
        }
      }
    }

    // rhs (Direchlet BC) assembly over elements

    numEdges = param.mesh.elementPtr(iElement) -> nVertices();

    // Column indexing of L: (global edge indexed by interior)
    // {edge(0),Func(0)}, {edge(0),Func(1)}, ..., {edge(0),Func(eePtr[0]->dim()-1)},
    // {edge(1),Func(0)}, {edge(1),Func(1)}, ..., {edge(1),Func(eePtr[1]->dim()-1)},
    // ..., {edge(nInteriorEdge()-1),Func(eePtr[nInteriorEdge()-1]->dim()-1)} 

    for (int iEdge = 0; iEdge < numEdges; iEdge ++) {
      edge_loc_to_glob = dmSpace.globalEdgeIndex(iElement,iEdge);

      if (dmSpace.bcType(edge_loc_to_glob) == EdgeBCType::boundary) { 
      // If the edge is on boundary, we add eval to rhs (Dirichlet BC)
      mePtr->initBasis(quadEdgeRule[edge_loc_to_glob].pts(), quadEdgeRule[edge_loc_to_glob].num());

      for (int iPt = 0; iPt < quadEdgeRule[edge_loc_to_glob].num(); iPt++) {
        for (int j = 0; j < loc_dimAfull; j++) {
          double v_jdotNu = mePtr->basisDotNu(j,iEdge,iPt);
          if (mePtr -> isVertexDoF(j)) { v_jdotNu *= dmSpace.vertex_loc_to_glob_coeff(iElement,mePtr->mapDoFToVertex(j)); }
            double bc = bcVal(quadEdgeRule[edge_loc_to_glob].pt(iPt)[0],quadEdgeRule[edge_loc_to_glob].pt(iPt)[1]);
            curr_full_index = dmSpace.loc_to_glob(iElement,j,'f');
            evaluation = bc * v_jdotNu * quadEdgeRule[edge_loc_to_glob].wt(iPt);
            rhs_bc_full[curr_full_index] -= evaluation;
            if (j < loc_dimAreduced) {
              curr_reduced_index = dmSpace.loc_to_glob(iElement,j,'r');
              rhs_bc_reduced[curr_reduced_index] -= evaluation;
            }
          }
        }
        
      }
    }

    starting_colBreduced += loc_colBreduced;
    starting_colBfull += loc_colBfull;
  }


  // OUTPUT MATRICES AND RHS //////////////////////////////////////////////////

  if (true) {
  
  std::ofstream fout1("test/A_full.txt");
  for(int j = 0; j < dimAfull; j++) {
    for(int i = 0; i < dimAfull; i++) {
      fout1 << mat_A_full[i + dimAfull*j] << "\t";
    }
    if (j < dimAfull - 1) fout1 << "\n";
  }

  std::ofstream fout2("test/A_reduced.txt");
  for(int j = 0; j < dimAreduced; j++) {
    for(int i = 0; i < dimAreduced; i++) {
      fout2 << mat_A_reduced[i + dimAreduced*j] << "\t";
    }
    if (j < dimAreduced - 1) fout2 << "\n";
  }

  std::ofstream fout3("test/B_full.txt");
  for(int j = 0; j < rowBfull; j++) {
    for(int i = 0; i < colBfull; i++) {
      fout3 << mat_B_full[i + colBfull*j] << "\t";
    }
    if (j < rowBfull - 1) fout3 << "\n";
  }

  std::ofstream fout4("test/B_reduced.txt");
  for(int j = 0; j < rowBreduced; j++) {
    for(int i = 0; i < colBreduced; i++) {
      fout4 << mat_B_reduced[i + colBreduced*j] << "\t";
    }
    if (j < rowBreduced - 1) fout4 << "\n";
  }


  std::ofstream fout7("test/rhs_full.txt");
  for(int i = 0; i < colBfull; i++) {
    fout7 << rhs_full[i];
    if (i < colBfull - 1) fout7 << "\n";
  }

  std::ofstream fout8("test/rhs_reduced.txt");
  for(int i = 0; i < colBreduced; i++) {
    fout8 << rhs_reduced[i];
    if (i < colBreduced - 1) fout8 << "\n";
  }

  std::ofstream fout9("test/rhs_bc_full.txt");
  for(int i = 0; i < dimAfull; i++) {
    fout9 << rhs_bc_full[i];
    if (i < dimAfull - 1) fout9 << "\n";
  }

  std::ofstream fout10("test/rhs_bc_reduced.txt");
  for(int i = 0; i < dimAreduced; i++) {
    fout10 << rhs_bc_reduced[i];
    if (i < dimAreduced - 1) fout10 << "\n";
  }

  }

  // SOLVE LINEAR SYSTEM //////////////////////////////////////////////////

  monitor(1,"Solution of linear system"); ////////////////////////////////////////
 
  // Solve the matrix
  // L^{T} A^{-1} { B (B^{T} A^{-1} B)^{-1} B^{T} A^{-1} L - L } c = - L^{T} A^{-1} B (B^{T} A^{-1} B)^{-1} rhs


  // Initialize A^{-1} and evaluate as A

  std::vector<double> A_full_inv_vector(dimAfull*dimAfull,0);
  double* A_full_inv = A_full_inv_vector.data();
  std::vector<double> A_reduced_inv_vector(dimAreduced*dimAreduced,0);
  double* A_reduced_inv = A_reduced_inv_vector.data();

  for (int j = 0; j < dimAfull; j++) {
    for (int i = 0; i < dimAfull; i++) {
      A_full_inv[j * dimAfull + i] = mat_A_full[j * dimAfull + i];
      if (i < dimAreduced && j < dimAreduced) {
        A_reduced_inv[j * dimAreduced + i] = mat_A_reduced[j * dimAreduced + i];
      }
    }
  }

  // Update A^{-1} to be the inverse matrix A^{-1} 


  lapack_int ierr = mat_inv(A_full_inv, dimAfull);
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }

  ierr = mat_inv(A_reduced_inv, dimAreduced);
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }
  

  // Calculate m2 = B^{T} A^{-1} B

  std::vector<double> m2_full_vector(colBfull*colBfull,0);
  double* m2_full = m2_full_vector.data();
  std::vector<double> m2_reduced_vector(colBreduced*colBreduced,0);
  double* m2_reduced = m2_reduced_vector.data();

  // First calculate B^{T}A^{-1}
  std::vector<double> btai_full_vector(colBfull*dimAfull,0);
  double* btai_full = btai_full_vector.data();
  std::vector<double> btai_reduced_vector(colBreduced*dimAreduced,0);
  double* btai_reduced = btai_reduced_vector.data();

  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, colBfull, dimAfull,
              dimAfull, 1, mat_B_full, colBfull, A_full_inv, dimAfull,
              0, btai_full, dimAfull);

  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, colBreduced, dimAreduced,
              dimAreduced, 1, mat_B_reduced, colBreduced, A_reduced_inv, dimAreduced,
              0, btai_reduced, dimAreduced);

  // Times B^{T}A^{-1} with B and store in m2
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colBfull, colBfull,
              dimAfull, 1, btai_full, dimAfull, mat_B_full, colBfull,
              0, m2_full, colBfull);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colBreduced, colBreduced,
              dimAreduced, 1, btai_reduced, dimAreduced, mat_B_reduced, colBreduced,
              0, m2_reduced, colBreduced);
 


 
  // We get rhs1 = - B^T A^{-1} rhs_bc + rhs
  std::vector<double> rhs1_full_vector(colBfull,0); // row number = dimAfull, col number = 1
  double* rhs1_full = rhs1_full_vector.data();
  std::vector<double> rhs1_reduced_vector(colBreduced,0);
  double* rhs1_reduced = rhs1_reduced_vector.data();

  for (int j = 0; j < colBfull; j++) {
    rhs1_full[j] = rhs_full[j];
    if (j < colBreduced) rhs1_reduced[j] = rhs_reduced[j];
  }

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colBfull, 1,
              dimAfull, -1, btai_full, dimAfull, rhs_bc_full, 1,
              1, rhs1_full, 1);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colBreduced, 1,
              dimAreduced, -1, btai_reduced, dimAreduced, rhs_bc_reduced, 1,
              1, rhs1_reduced, 1);  



  // Now we initialize b with rhs1 and solve the linear system
  // Notice that m2*b=rhs1

  std::vector<double> b_full_vector(colBfull,0);
  double* b_full = b_full_vector.data();
  std::vector<double> b_reduced_vector(colBreduced,0);
  double* b_reduced = b_reduced_vector.data();

  for (int j = 0; j < colBfull; j++) {
    b_full[j] = rhs1_full[j];
    if (j < colBreduced) b_reduced[j] = rhs1_reduced[j];
  }


  lapack_int* ipiv_full; lapack_int* ipiv_reduced; char norm = 'I'; 
  ipiv_full = (lapack_int*)malloc(colBfull * sizeof(lapack_int));
  ipiv_reduced = (lapack_int*)malloc(colBreduced * sizeof(lapack_int));
  double anorm_full = LAPACKE_dlange(LAPACK_ROW_MAJOR, norm, colBfull, colBfull, m2_full, colBfull);
  double anorm_reduced = LAPACKE_dlange(LAPACK_ROW_MAJOR, norm, colBreduced, colBreduced, m2_reduced, colBreduced);
  ierr = LAPACKE_dgesv(LAPACK_ROW_MAJOR, colBfull, 1, m2_full, colBfull, ipiv_full, b_full, 1); //mat updated to be LU
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }
    ierr = LAPACKE_dgesv(LAPACK_ROW_MAJOR, colBreduced, 1, m2_reduced, colBreduced, ipiv_reduced, b_reduced, 1); //mat updated to be LU
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }

  double rcond_full = 0; double rcond_reduced = 0;
  ierr = LAPACKE_dgecon(LAPACK_ROW_MAJOR, norm, colBfull, m2_full, colBfull, anorm_full, &rcond_full);
  if(ierr) {
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }
  ierr = LAPACKE_dgecon(LAPACK_ROW_MAJOR, norm, colBreduced, m2_reduced, colBreduced, anorm_reduced, &rcond_reduced);
  if(ierr) {
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }
  rcond_full = 1/rcond_full;
  rcond_reduced = 1/rcond_reduced;

  //Calculate inf condition number
  std::cout << "\tNorm Format:\t" << norm << std::endl;
  std::cout << "\tNorm of mat:\t" << "full:\t" << anorm_full << "\treduced:\t"<< anorm_reduced << std::endl;
  std::cout << "\tCond number:\t" << "full:\t" << rcond_full << "\treduced:\t"<< rcond_reduced  << std::endl;

  std::ofstream fout17("test/b_vec_full.txt");
  for(int i = 0; i < colBfull; i++) {
    fout17 << b_full[i];
    if (i < colBfull - 1) fout17 << "\n";
  }

  std::ofstream fout18("test/b_vec_reduced.txt");
  for(int i = 0; i < colBreduced; i++) {
    fout18 << b_reduced[i];
    if (i < colBreduced - 1) fout18 << "\n";
  }



  // Calculate a = A^{-1} [Bb+rhs_bc]

  // We first initialize and calculate  a1 = rhs_bc, and then a1 = Bb + a1


  std::vector<double> a1_full_vector(dimAfull,0);
  double* a1_full = a1_full_vector.data();
  std::vector<double> a1_reduced_vector(dimAreduced,0);
  double* a1_reduced = a1_reduced_vector.data();

  std::vector<double> a_full_vector(dimAfull,0);
  double* a_full = a_full_vector.data();
  std::vector<double> a_reduced_vector(dimAreduced,0);
  double* a_reduced = a_reduced_vector.data();


  // a1 = rhs_bc
  for (int j = 0; j < dimAfull; j++) {
    a1_full[j] = rhs_bc_full[j];
    if (j < dimAreduced) a1_reduced[j] = rhs_bc_reduced[j];
  }

  // a1 = Bb + a1

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimAfull, 1,
              colBfull, 1, mat_B_full, colBfull, b_full, 1,
              1, a1_full, 1);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimAreduced, 1,
              colBreduced, 1, mat_B_reduced, colBreduced, b_reduced, 1,
              1, a1_reduced, 1);


  // a = A^{-1}(a1)
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimAfull, 1,
              dimAfull, 1, A_full_inv, dimAfull, a1_full, 1,
              0, a_full, 1);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimAreduced, 1,
              dimAreduced, 1, A_reduced_inv, dimAreduced, a1_reduced, 1,
              0, a_reduced, 1);


  std::ofstream fout21("test/a_vec_full.txt");
  for(int i = 0; i < dimAfull; i++) {
    fout21 << a_full[i];
    if (i < dimAfull - 1) fout21 << "\n";
  }

  std::ofstream fout22("test/a_vec_reduced.txt");
  for(int i = 0; i < dimAreduced; i++) {
    fout22 << a_reduced[i];
    if (i < dimAreduced - 1) fout22 << "\n";
  }
  

  for(int i=0; i<solution_p_f.size(); i++) {
    solution_p_f[i] = b_full[i];
  }

  for(int i=0; i<solution_p_r.size(); i++) {
    solution_p_r[i] = b_reduced[i];
  }

  for(int i=0; i<solution_u_f.size(); i++) {
    solution_u_f[i] = a_full[i];
  }

  for(int i=0; i<solution_u_r.size(); i++) {
    solution_u_r[i] = a_reduced[i];
  }

  if(trueSolnKnown()) {
    monitor(0,"\nError estimate\n"); ///////////////////////////////////////////////
  
    double h = param.dsSpace.mesh()->maxElementDiameter();
    double maxChunk = param.dsSpace.mesh()->maxChunkParam();
    double averageChunk = param.dsSpace.mesh()->averageChunkParam();
    
    double l2Error_f = 0, l2UError_f = 0, l2DivUError_f = 0, l2Norm_f = 0, l2UNorm_f = 0, l2DivUNorm_f = 0;
    double l2Error_r = 0, l2UError_r = 0, l2DivUError_r = 0, l2Norm_r = 0, l2UNorm_r = 0, l2DivUNorm_r = 0;

    solution_p_f.l2normError(l2Error_f, l2Norm_f, param.refinement_level*0, trueSoln);
    solution_p_r.l2normError(l2Error_r, l2Norm_r, param.refinement_level*0, trueSoln);
    solution_u_f.l2normError(l2UError_f, l2UNorm_f, param.refinement_level*0, trueUSoln);
    solution_u_r.l2normError(l2UError_r, l2UNorm_r, param.refinement_level*0, trueUSoln);
    solution_u_f.l2normError_div(l2DivUError_f, l2DivUNorm_f, param.refinement_level*0, trueDivUSoln);
    solution_u_r.l2normError_div(l2DivUError_r, l2DivUNorm_r, param.refinement_level*0, trueDivUSoln);
    
    std::cout << "  Max Element Diameter h:  " << h << std::endl;
    std::cout << "  Max Chunkiness Parameter:  " << maxChunk << std::endl;
    std::cout << "  Average Chunkiness Parameter:  " << averageChunk << std::endl;
    std::cout << "  === p ===  " << std::endl;
    std::cout << "  L_2 Error full:      " << l2Error_f << std::endl;
    std::cout << "  L_2 Error reduced:      " << l2Error_r << std::endl;
    std::cout << "  Relative L_2 Error full:      " << l2Error_f/l2Norm_f << std::endl;
    std::cout << "  Relative L_2 Error reduced:      " << l2Error_r/l2Norm_r << std::endl;
    std::cout << std::endl;
    std::cout << "  === u ===  " << std::endl;
    std::cout << "  L_2 Error full:      " << l2UError_f << std::endl;
    std::cout << "  L_2 Error reduced:      " << l2UError_r << std::endl;
    std::cout << "  Relative L_2 Error full:      " << l2UError_f/l2UNorm_f << std::endl;
    std::cout << "  Relative L_2 Error reduced:      " << l2UError_r/l2UNorm_r << std::endl;
    std::cout << std::endl;
    std::cout << "  === div u ===  " << std::endl;
    std::cout << "  L_2 Error full:      " << l2DivUError_f << std::endl;
    std::cout << "  L_2 Error reduced:      " << l2DivUError_r << std::endl;
    std::cout << "  Relative L_2 Error full:      " << l2DivUError_f/l2DivUNorm_f << std::endl;
    std::cout << "  Relative L_2 Error reduced:      " << l2DivUError_r/l2DivUNorm_r << std::endl;
    std::cout << std::endl;


  }

  if(param.output_soln_Mixed_format > 0) {
  
    monitor(1,"Write Solution"); //////////////////////////////////////////////////

    switch(param.output_soln_Mixed_format) {
    case 1: {
      std::string fileName(param.directory_name);
      fileName += "solution_raw";

      std::string fileName_p_f = fileName + "_p_f";
      solution_p_f.write_raw(fileName_p_f);

      std::string fileName_p_r = fileName + "_p_r";
      solution_p_r.write_raw(fileName_p_r);

       
      std::string fileName_u_f = fileName + "_u_f";
      solution_u_f.write_raw(fileName_u_f);

      std::string fileName_u_r = fileName + "_u_r";
      solution_u_r.write_raw(fileName_u_r);

      break;
    }
    case 2: {
      std::string fileName(param.directory_name);
      fileName += "solution_mesh";

      std::string fileName_p_f = fileName + "_p_f";
      solution_p_f.write_matlab_mesh(fileName_p_f, 
                    param.output_mesh_numPts_Mixed_x,param.output_mesh_numPts_Mixed_y);

      std::string fileName_p_r = fileName + "_p_r";
      solution_p_r.write_matlab_mesh(fileName_p_r, 
                    param.output_mesh_numPts_Mixed_x,param.output_mesh_numPts_Mixed_y);

       
      std::string fileName_u_f = fileName + "_u_f";
      solution_u_f.write_matlab_mesh(fileName_u_f, 
                    param.output_mesh_numPts_Mixed_x,param.output_mesh_numPts_Mixed_y);

      std::string fileName_u_r = fileName + "_u_r";
      solution_u_r.write_matlab_mesh(fileName_u_r, 
                    param.output_mesh_numPts_Mixed_x,param.output_mesh_numPts_Mixed_y);


      if (trueSolnKnown()) {
        fileName_p_f += "_error";
        solution_p_f.write_matlab_mesh_error(fileName_p_f, 
                      param.output_mesh_numPts_Mixed_x,param.output_mesh_numPts_Mixed_y, trueSoln);

        fileName_p_r += "_error";
        solution_p_r.write_matlab_mesh_error(fileName_p_r, 
                      param.output_mesh_numPts_Mixed_x,param.output_mesh_numPts_Mixed_y, trueSoln);

        
        fileName_u_f += "_error";
        solution_u_f.write_matlab_mesh_error(fileName_u_f, 
                      param.output_mesh_numPts_Mixed_x,param.output_mesh_numPts_Mixed_y, trueUSoln);

        fileName_u_r += "_error";
        solution_u_r.write_matlab_mesh_error(fileName_u_r, 
                      param.output_mesh_numPts_Mixed_x,param.output_mesh_numPts_Mixed_y, trueUSoln);

        
        fileName_u_f += "_div";
        solution_u_f.write_matlab_mesh_div_error(fileName_u_f, 
                      param.output_mesh_numPts_Mixed_x,param.output_mesh_numPts_Mixed_y, trueDivUSoln);

        fileName_u_r += "_div";
        solution_u_r.write_matlab_mesh_div_error(fileName_u_r, 
                      param.output_mesh_numPts_Mixed_x,param.output_mesh_numPts_Mixed_y, trueDivUSoln);   
      }    
      
      break;
    }
    }
  }
}
  return 0;
};