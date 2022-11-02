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


using namespace directserendipity;
using namespace polymesh;
using namespace polyquadrature;

lapack_int mat_inv(double *A, int n)
{
  // inplace inverse n x n matrix A.
  // matrix A is Column Major (i.e. firts line, second line ... *not* C[][] order)
  // returns:
  //   ret = 0 on success
  //   ret < 0 illegal argument value
  //   ret > 0 singular matrix
  int ipiv[n+1];
  lapack_int ret;

  ret =  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n,n,A,n,ipiv);

  if (ret !=0) return ret;
  ret = LAPACKE_dgetri(LAPACK_ROW_MAJOR,n,A,n,ipiv);
  return ret;
}


lapack_int block_mat_inv(double *A, int n, int *m, int num)
{
  // inplace inverse n x n matrix A.
  // A is a block diagonal matrix with num blocks, each with size m_i x m_i

  lapack_int ret = 0;
  int starting_index = 0;

  for (int i = 0; i < num; i++) {

    std::vector<double> small_mat_vector(m[i]*m[i],0);
    double* small_mat = small_mat_vector.data();

    // assign the i-th block to small_mat
    for (int row = 0; row < m[i]; row++) {
      for (int col = 0; col < m[i]; col++) {
        small_mat[row*m[i]+col] = A[n*(starting_index+row)+starting_index+col];
      }
    }

    ret = mat_inv(small_mat,m[i]);

    // assign the inverse small matrix to A
    for (int row = 0; row < m[i]; row++) {
      for (int col = 0; col < m[i]; col++) {
        A[n*(starting_index+row)+starting_index+col] = small_mat[row*m[i]+col];
      }
    }

    starting_index += m[i];    
  }

  return ret;
}


int MixedPDE::solve_hybrid(Monitor& monitor) {
  monitor(0,"Polynomial Degree = ", parameterDataPtr()->dmSpace.degPolyn());

  ParameterData& param = *parameterDataPtr();
  DirectMixedHybrid dmSpace = DirectMixedHybrid(parameterDataPtr()->dmSpace);
  // TEST SINGLE ELEMENT MESH //////////////////////////////////////////////

  if(false) {
    monitor(0,"Test Single Element Mesh and Space");

    polymesh::PolyMesh s_mesh(parameterDataPtr()->mesh.elementPtr(3));
    std::string fileName = parameterDataPtr()->directory_name;
    fileName += "mesh_e4";
    s_mesh.write_matlab(fileName);

    DirectMixedHybrid s_dmSpace(8,&s_mesh);
    fileName = parameterDataPtr()->directory_name;
    fileName += "dmSpace_e4";
    s_dmSpace.write_matlab(fileName);
  }
  
  // TEST BASIS FUNCTIONS //////////////////////////////////////////////////

  if(true) {
    monitor(0,"\nTest basis functions\n");

    
    DirectMixedHybridArray u(&dmSpace, 'f');
    DirectDGArray p(&(parameterDataPtr()->dmSpace), 'f');
    DirectEdgeDGArray l(&dmSpace);

    // Set the coefficients of all the mixed hybrid basis functions to be zero
    for(int i=0; i<u.size(); i++) {
      u[i]=0;
    }
    // Set the coefficients of the first mixed hybrid basis function to be one
    u[0]=1;

    // Set the coefficients of all the discrete Galerkin basis functions to be zero
    for(int i=0; i<p.size(); i++) {
      p[i]=0;
    }
    // Set the coefficients of the first discrete Galerkin basis function to be one
    p[0]=1;

    // Set the coefficients of all the discrete Galerkin basis functions on edges to be zero
    for(int i=0; i<l.size(); i++) {
      l[i]=0;
    }
    // Set the coefficients of the first discrete Galerkin basis function on edges to be one
    l[0]=1;

    monitor(1,"Write Array");

    std::string fileName = parameterDataPtr()->directory_name;
    fileName += "basis_mesh_mixed";

    u.write_matlab_mesh(fileName,51,51);

    fileName = parameterDataPtr()->directory_name;
    fileName += "basis_mesh_DG";
    p.write_matlab_mesh(fileName,51,51);

    fileName = parameterDataPtr()->directory_name;
    fileName += "basis_mesh_EDG";
    l.write_matlab_mesh(fileName,51,51);
  }

  // TEST QUADRATURE ///////////////////////////////////////////////////////

  if(false) {
    monitor(1,"Test Quadrature Rules");
    testPolyQuadrature(&(parameterDataPtr()->mesh),1e-6);
    return 0;
  }
  
  // SOLVE THE PDE ///////////////////////////////////////////////////////////

  monitor(0,"\nSolve the PDE\n");

  // Only when r > 0, we calculate reduced space

  DirectMixedHybridArray solution_u_f(&(dmSpace),'f');
  DirectMixedHybridArray solution_u_r(&(dmSpace),'r');
  DirectDGArray solution_p_f(&(dmSpace),'f');
  DirectDGArray solution_p_r(&(dmSpace),'r');
  DirectEdgeDGArray solution_l_f(&(dmSpace));
  DirectEdgeDGArray solution_l_r(&(dmSpace));


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

  // Initialize matrix L
  int rowLfull = dimAfull, rowLreduced = dimAreduced;
  int colL = dmSpace.nIntEdgeDGDoFs();

  std::vector<double> mat_L_full_vector(rowLfull*colL,0);
  double* mat_L_full = mat_L_full_vector.data();
  std::vector<double> mat_L_reduced_vector(rowLreduced*colL,0);
  double* mat_L_reduced = mat_L_reduced_vector.data();

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

  // Initialize an array storing local dimension of A
  std::vector<int> locdim_full_vector(param.mesh.nElements(),0);
  int* locdim_full = locdim_full_vector.data();

  std::vector<int> locdim_reduced_vector(param.mesh.nElements(),0);
  int* locdim_reduced = locdim_reduced_vector.data();

  // quadrature points
  // Note that we update quadRule in each iElement loop
  // But we define quadEdgeRule for all edges at once

  polyquadrature::PolyQuadrature quadRule(13,param.refinement_level);

  std::vector<polyquadrature::PolyEdgeQuadrature> quadEdgeRule(dmSpace.nEdges());

  for (int i = 0; i < dmSpace.nEdges(); i++) {
    quadEdgeRule[i].set(13, dmSpace.edgePtr(i));
  }

  monitor(1,"Matrix and RHS Assembly"); ////////////////////////////////////////

  int starting_Afull = 0, starting_Areduced = 0;
  int starting_colBfull = 0, starting_colBreduced = 0;

  int loc_dimAfull = 0, loc_dimAreduced = 0, loc_colBfull = 0, loc_colBreduced = 0;

  int curr_full_index, curr_reduced_index;
  double evaluation;
  int numEdges = 0, loc_to_int = 0, loc_to_glob;

  // We construct an array of pointer to all the DirectEdgeDGFE
  // and get them by interior edge indexing

  std::vector<DirectEdgeDGFE*> eePtr(dmSpace.nEdges());
  for (int i = 0; i < dmSpace.nEdges(); i++) {
    eePtr[i] = dmSpace.DGEdgePtr(i);
    eePtr[i] -> initBasis(quadEdgeRule[i].pts(), quadEdgeRule[i].num());
  }


  for(int iElement=0; iElement<param.mesh.nElements(); iElement++) {
    DirectMixedHybridFE* mePtr = dmSpace.MixedElementPtr(iElement);
    DirectDGFE* dgePtr = dmSpace.DGElementPtr(iElement);

    quadRule.setElement(mePtr->elementPtr());
    
    mePtr->initBasis(quadRule.pts(), quadRule.num());
    dgePtr->initBasis(quadRule.pts(), quadRule.num());

    loc_dimAfull = mePtr -> dimVFull();
    loc_dimAreduced = mePtr -> dimVReduced();
    loc_colBfull = dgePtr -> dimFull();
    loc_colBreduced = dgePtr -> dimReduced();

    locdim_full[iElement] = loc_dimAfull;
    locdim_reduced[iElement] = loc_dimAreduced;

    // Matrix A, B and rhs assembly over elements
 
    for(int iPt=0; iPt<quadRule.num(); iPt++) {
      double x = quadRule.pt(iPt)[0];
      double y = quadRule.pt(iPt)[1];

      Tensor2 valD;
      coefD_inv(x,y,valD);
  
      // Assemble matrix A
      for (int j = 0; j < loc_dimAfull; j++) {
        Tensor1 v_j = mePtr -> basis(j,iPt);
        for (int i = 0; i < loc_dimAfull; i++) {
          Tensor1 u_i = mePtr -> basis(i,iPt);
          curr_full_index = dimAfull * (starting_Afull + j) + (starting_Afull + i);
          evaluation = ((valD*u_i)*v_j)*quadRule.wt(iPt);
          mat_A_full[curr_full_index] += evaluation;
          if (j < loc_dimAreduced && i < loc_dimAreduced) {
            curr_reduced_index = dimAreduced * (starting_Areduced + j) + (starting_Areduced + i);
            mat_A_reduced[curr_reduced_index] += evaluation;
          }
        }
      }

      // Assemble matrix B
      // For first j = 0 -> dimCurlPart()-1 rows, div(v_j) = 0, 
      // so we only need to consider divXPo part

      for (int j = mePtr -> dimCurlPart(); j < loc_dimAfull; j++) {
        double divv_j = mePtr -> basisdivXPo(j - mePtr -> dimCurlPart(),iPt);
        for (int i = 0; i < loc_colBfull; i++ ) {
          double p_i = dgePtr -> basis(i,iPt);
          curr_full_index = colBfull * (starting_Afull + j) + (starting_colBfull + i);
          evaluation = divv_j * p_i * quadRule.wt(iPt);
          mat_B_full[curr_full_index] += evaluation;
          if (j < loc_dimAreduced && i < loc_colBreduced) {
            curr_reduced_index = colBreduced * (starting_Areduced + j) + (starting_colBreduced + i);
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

    // Matrix L and rhs (Direchlet BC) assembly over elements

    numEdges = param.mesh.elementPtr(iElement) -> nVertices();

    // Column indexing of L: (global edge indexed by interior)
    // {edge(0),Func(0)}, {edge(0),Func(1)}, ..., {edge(0),Func(eePtr[0]->dim()-1)},
    // {edge(1),Func(0)}, {edge(1),Func(1)}, ..., {edge(1),Func(eePtr[1]->dim()-1)},
    // ..., {edge(nInteriorEdge()-1),Func(eePtr[nInteriorEdge()-1]->dim()-1)} 

    for (int iEdge = 0; iEdge < numEdges; iEdge ++) {
      loc_to_glob = dmSpace.globalEdgeIndex(iElement,iEdge);
      loc_to_int = dmSpace.interiorEdgeIndex(iElement,iEdge);
      mePtr->initBasis(quadEdgeRule[loc_to_glob].pts(), quadEdgeRule[loc_to_glob].num());

      for (int iPt = 0; iPt < quadEdgeRule[loc_to_glob].num(); iPt++) {
        for (int j = 0; j < loc_dimAfull; j++) {
          double v_jdotNu = mePtr->basisDotNu(j,iEdge,iPt);
          if (loc_to_int == -1) { // If the edge is on boundary, we add eval to rhs (Dirichlet BC)
            double bc = bcVal(quadEdgeRule[loc_to_glob].pt(iPt)[0],quadEdgeRule[loc_to_glob].pt(iPt)[1]);
            curr_full_index = starting_Afull + j;
            evaluation = bc * v_jdotNu * quadEdgeRule[loc_to_glob].wt(iPt);
            rhs_bc_full[curr_full_index] -= evaluation;
            if (j < loc_dimAreduced) {
              curr_reduced_index = starting_Areduced + j;
              rhs_bc_reduced[curr_reduced_index] -= evaluation;
            }
          } else { // If the edge is interior, we add eval to matrix L
            for (int i = 0; i < eePtr[loc_to_glob]->dim(); i++){
              double l_i = eePtr[loc_to_glob]->basis(i,iPt);
              // Here we use the property that eePtr[loc_to_glob]->dim() is the same (degPolyn()+1)for every edge in our mesh
              curr_full_index = colL * (starting_Afull + j) + (loc_to_int*(dmSpace.degPolyn()+1)+i);
              evaluation = l_i * v_jdotNu * quadEdgeRule[loc_to_glob].wt(iPt);
              mat_L_full[curr_full_index] += evaluation;
              if (j < loc_dimAreduced) {
                curr_reduced_index = colL * (starting_Areduced + j) + (loc_to_int*(dmSpace.degPolyn()+1)+i);
                mat_L_reduced[curr_reduced_index] += evaluation;
              }
            }
          }
        }  
      }
    }

    starting_Afull += loc_dimAfull;
    starting_Areduced += loc_dimAreduced;
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

  std::ofstream fout5("test/L_full.txt");
  for(int j = 0; j < rowLfull; j++) {
    for(int i = 0; i < colL; i++) {
      fout5 << mat_L_full[i + colL*j] << "\t";
    }
    if (j < rowLfull - 1) fout5 << "\n";
  }

  std::ofstream fout6("test/L_reduced.txt");
  for(int j = 0; j < rowLreduced; j++) {
    for(int i = 0; i < colL; i++) {
      fout6 << mat_L_reduced[i + colL*j] << "\t";
    }
    if (j < rowLreduced - 1) fout6 << "\n";
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

  lapack_int ierr = block_mat_inv(A_full_inv, dimAfull, locdim_full, param.mesh.nElements());
  //lapack_int ierr = mat_inv(A_full_inv, dimAfull);
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }

  ierr = block_mat_inv(A_reduced_inv, dimAreduced, locdim_reduced, param.mesh.nElements());
  //ierr = mat_inv(A_reduced_inv, dimAreduced);
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }
  
  // Calculate m1 = L^{T} A^{-1}

  std::vector<double> m1_full_vector(colL*dimAfull,0);
  double* m1_full = m1_full_vector.data();
  std::vector<double> m1_reduced_vector(colL*dimAreduced,0);
  double* m1_reduced = m1_reduced_vector.data();


  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, colL, dimAfull,
              dimAfull, 1, mat_L_full, colL, A_full_inv, dimAfull,
              0, m1_full, dimAfull);

  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, colL, dimAreduced,
              dimAreduced, 1, mat_L_reduced, colL, A_reduced_inv, dimAreduced,
              0, m1_reduced, dimAreduced);


  // Calculate m2 = ( B^{T} A^{-1} B )^{-1}

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

  // Take inverse
  ierr = mat_inv(m2_full, colBfull);
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }

  ierr = mat_inv(m2_reduced, colBreduced);
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }

  // Calculate m3 = B^{T} A^{-1} L

  std::vector<double> m3_full_vector(colBfull*colL,0);
  double* m3_full = m3_full_vector.data();
  std::vector<double> m3_reduced_vector(colBreduced*colL,0);
  double* m3_reduced = m3_reduced_vector.data();

  // We already have btai, we just time it with L to get m3
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colBfull, colL,
              dimAfull, 1, btai_full, dimAfull, mat_L_full, colL,
              0, m3_full, colL);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colBreduced, colL,
              dimAreduced, 1, btai_reduced, dimAreduced, mat_L_reduced, colL,
              0, m3_reduced, colL);

  // Calculate and store B(m2)

  std::vector<double> bm2_full_vector(dimAfull*colBfull,0);
  double* bm2_full = bm2_full_vector.data();
  std::vector<double> bm2_reduced_vector(dimAreduced*colBreduced,0);
  double* bm2_reduced = bm2_reduced_vector.data();

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimAfull, colBfull,
              colBfull, 1, mat_B_full, colBfull, m2_full, colBfull,
              0, bm2_full, colBfull);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimAreduced, colBreduced,
              colBreduced, 1, mat_B_reduced, colBreduced, m2_reduced, colBreduced,
              0, bm2_reduced, colBreduced);

  // Calculate and store m4 = B m2 m3 - L

  std::vector<double> m4_full_vector(dimAfull*colL,0);
  double* m4_full = m4_full_vector.data();
  std::vector<double> m4_reduced_vector(dimAreduced*colL,0);
  double* m4_reduced = m4_reduced_vector.data();  

  // We first evaluate m4 as L
  for (int j = 0; j < dimAfull; j++) {
    for (int i = 0; i < colL; i++) {
      m4_full[j * colL + i] = mat_L_full[j * colL + i];
      if (j < dimAreduced) {
        m4_reduced[j * colL + i] = mat_L_reduced[j * colL + i];
      }
    }
  }

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimAfull, colL,
              colBfull, 1, bm2_full, colBfull, m3_full, colL,
              -1, m4_full, colL);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimAreduced, colL,
              colBreduced, 1, bm2_reduced, colBreduced, m3_reduced, colL,
              -1, m4_reduced, colL);

  // Finally, we have mat = m1 m4
  std::vector<double> mat_full_vector(colL*colL,0);
  double* mat_full = mat_full_vector.data();
  std::vector<double> mat_reduced_vector(colL*colL,0);
  double* mat_reduced = mat_reduced_vector.data();

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colL, colL,
              dimAfull, 1, m1_full, dimAfull, m4_full, colL,
              0, mat_full, colL);          

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colL, colL,
              dimAreduced, 1, m1_reduced, dimAreduced, m4_reduced, colL,
              0, mat_reduced, colL);  

  // We get rhs1 = B^T A^{-1} rhs_bc - rhs
  std::vector<double> rhs1_full_vector(colBfull,0); // row number = dimAfull, col number = 1
  double* rhs1_full = rhs1_full_vector.data();
  std::vector<double> rhs1_reduced_vector(colBreduced,0);
  double* rhs1_reduced = rhs1_reduced_vector.data();

  for (int j = 0; j < colBfull; j++) {
    rhs1_full[j] = rhs_full[j];
    if (j < colBreduced) rhs1_reduced[j] = rhs_reduced[j];
  }

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colBfull, 1,
              dimAfull, 1, btai_full, dimAfull, rhs_bc_full, 1,
              -1, rhs1_full, 1);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colBreduced, 1,
              dimAreduced, 1, btai_reduced, dimAreduced, rhs_bc_reduced, 1,
              -1, rhs1_reduced, 1);  


  // rhs2 = B (m2) rhs1 - rhs_bc
  std::vector<double> rhs2_full_vector(dimAfull,0); // row number = dimAfull, col number = 1
  double* rhs2_full = rhs2_full_vector.data();
  std::vector<double> rhs2_reduced_vector(dimAreduced,0);
  double* rhs2_reduced = rhs2_reduced_vector.data();

  for (int j = 0; j < dimAfull; j++) {
    rhs2_full[j] = rhs_bc_full[j];
    if (j < dimAreduced) rhs2_reduced[j] = rhs_bc_reduced[j];
  }

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimAfull, 1,
              colBfull, 1, bm2_full, colBfull, rhs1_full, 1,
              -1, rhs2_full, 1);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimAreduced, 1,
              colBreduced, 1, bm2_reduced, colBreduced, rhs1_reduced, 1,
              -1, rhs2_reduced, 1);

  // Now we initialize c with (m1)(rhs2) and solve the linear system
  // Notice that mat*c = (m1)(rhs2)

  std::vector<double> c_full_vector(colL,0);
  double* c_full = c_full_vector.data();
  std::vector<double> c_reduced_vector(colL,0);
  double* c_reduced = c_reduced_vector.data();

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colL, 1,
              dimAfull, 1, m1_full, dimAfull, rhs2_full, 1,
              0, c_full, 1);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colL, 1,
              dimAreduced, 1, m1_reduced, dimAreduced, rhs2_reduced, 1,
              0, c_reduced, 1);


  lapack_int* ipiv; char norm = 'I'; 
  ipiv = (lapack_int*)malloc(colL * sizeof(lapack_int));
  double anorm = LAPACKE_dlange(LAPACK_ROW_MAJOR, norm, colL, colL, mat_full, colL);
  ierr = LAPACKE_dgesv(LAPACK_ROW_MAJOR, colL, 1, mat_full, colL, ipiv, c_full, 1); //mat updated to be LU
  if(ierr) {
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }
  double rcond = 0;
  ierr = LAPACKE_dgecon(LAPACK_ROW_MAJOR, norm, colL, mat_full, colL, anorm, &rcond);
  if(ierr) {
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }
  rcond = 1/rcond;


  double anorm_r = LAPACKE_dlange(LAPACK_ROW_MAJOR, norm, colL, colL, mat_reduced, colL);
  ierr = LAPACKE_dgesv(LAPACK_ROW_MAJOR, colL, 1, mat_reduced, colL, ipiv, c_reduced, 1); //mat updated to be LU
  if(ierr) {
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }
  double rcond_r = 0;
  ierr = LAPACKE_dgecon(LAPACK_ROW_MAJOR, norm, colL, mat_reduced, colL, anorm_r, &rcond_r);
  if(ierr) {
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }
  rcond_r = 1/rcond_r;

  //Calculate inf condition number
  std::cout << "\tNorm Format:\t" << norm << std::endl;
  std::cout << "\tNorm of mat:\t" << "full:\t" << anorm << "\treduced:\t"<< anorm_r << std::endl;
  std::cout << "\tCond number:\t" << "full:\t" << rcond << "\treduced:\t"<< rcond_r  << std::endl;

  std::ofstream fout17("test/c_full.txt");
  for(int i = 0; i < colL; i++) {
    fout17 << c_full[i];
    if (i < colL - 1) fout17 << "\n";
  }

  std::ofstream fout18("test/c_reduced.txt");
  for(int i = 0; i < colL; i++) {
    fout18 << c_reduced[i];
    if (i < colL - 1) fout18 << "\n";
  }

  // Calculate b = (m2) [ (m3) c - B^T A^{-1} rhs_bc + rhs ]
  // Note that rhs1 = B^T A^{-1} rhs_bc - rhs
  // So b = (m2) [ (m3) c - rhs1 ]

  std::vector<double> b0_full_vector(colBfull,0);
  double* b0_full = b0_full_vector.data();
  std::vector<double> b0_reduced_vector(colBreduced,0);
  double* b0_reduced = b0_reduced_vector.data();

  std::vector<double> b_full_vector(colBfull,0);
  double* b_full = b_full_vector.data();
  std::vector<double> b_reduced_vector(colBreduced,0);
  double* b_reduced = b_reduced_vector.data();

  // We first initialize b0 as rhs1 and calculate b0 = (m3) c - rhs1
  for (int j = 0; j < colBfull; j++) {
    b0_full[j] = rhs1_full[j];
    if (j < colBreduced) b0_reduced[j] = rhs1_reduced[j];
  }

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colBfull, 1,
              colL, 1, m3_full, colL, c_full, 1,
              -1, b0_full, 1);  

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colBreduced, 1,
              colL, 1, m3_reduced, colL, c_reduced, 1,
              -1, b0_reduced, 1);

  // Calculate b = m2(b0)
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colBfull, 1,
              colBfull, 1, m2_full, colBfull, b0_full, 1,
              0, b_full, 1);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colBreduced, 1,
              colBreduced, 1, m2_reduced, colBreduced, b0_reduced, 1,
              0, b_reduced, 1);

  // Calculate a = A^{-1} [Bb-Lc+rhs_bc]

  // We first initialize and calculate a0 = Lc, and then a1 = Bb - a0, then a1 = a1 + rhs_bc

  std::vector<double> a0_full_vector(dimAfull,0);
  double* a0_full = a0_full_vector.data();
  std::vector<double> a0_reduced_vector(dimAreduced,0);
  double* a0_reduced = a0_reduced_vector.data();

  std::vector<double> a1_full_vector(dimAfull,0);
  double* a1_full = a1_full_vector.data();
  std::vector<double> a1_reduced_vector(dimAreduced,0);
  double* a1_reduced = a1_reduced_vector.data();

  std::vector<double> a_full_vector(dimAfull,0);
  double* a_full = a_full_vector.data();
  std::vector<double> a_reduced_vector(dimAreduced,0);
  double* a_reduced = a_reduced_vector.data();

  // a0 = Lc
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimAfull, 1,
              colL, 1, mat_L_full, colL, c_full, 1,
              0, a0_full, 1);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimAreduced, 1,
              colL, 1, mat_L_reduced, colL, c_reduced, 1,
              0, a0_reduced, 1);

  // a1 = Bb - a0
  for (int j = 0; j < dimAfull; j++) {
    a1_full[j] = a0_full[j];
    if (j < dimAreduced) a1_reduced[j] = a0_reduced[j];
  }

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimAfull, 1,
              colBfull, 1, mat_B_full, colBfull, b_full, 1,
              -1, a1_full, 1);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimAreduced, 1,
              colBreduced, 1, mat_B_reduced, colBreduced, b_reduced, 1,
              -1, a1_reduced, 1);

  // a1 = a1 + rhs_bc
  for (int j = 0; j < dimAfull; j++) {
    a1_full[j] += rhs_bc_full[j];
    if (j < dimAreduced) a1_reduced[j] += rhs_bc_reduced[j];
  }

  // a = A^{-1}(a1)
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimAfull, 1,
              dimAfull, 1, A_full_inv, dimAfull, a1_full, 1,
              0, a_full, 1);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimAreduced, 1,
              dimAreduced, 1, A_reduced_inv, dimAreduced, a1_reduced, 1,
              0, a_reduced, 1);

  std::ofstream fout19("test/b_vec_full.txt");
  for(int i = 0; i < colBfull; i++) {
    fout19 << b_full[i];
    if (i < colBfull - 1) fout19 << "\n";
  }

  std::ofstream fout20("test/b_vec_reduced.txt");
  for(int i = 0; i < colBreduced; i++) {
    fout20 << b_reduced[i];
    if (i < colBreduced - 1) fout20 << "\n";
  }

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

  
    // Calculate conf = L^{T}*u

    std::vector<double> conf_full_vector(colL,0);
    double* conf_full = conf_full_vector.data();
    std::vector<double> conf_reduced_vector(colL,0);
    double* conf_reduced = conf_reduced_vector.data();


    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, colL, 1,
                dimAfull, 1, mat_L_full, colL, a_full,1,
                0, conf_full,1);

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, colL, 1,
                dimAreduced, 1, mat_L_reduced, colL, a_reduced, 1,
                0, conf_reduced, 1);


    double conf_l1Error_f = 0, conf_l2Error_f = 0, conf_linfError_f = 0;
    double conf_l1Error_r = 0, conf_l2Error_r = 0, conf_linfError_r = 0;

    for(int i=0; i<colL; i++) {
      conf_l1Error_f += fabs(conf_full[i]);
      conf_l2Error_f += pow(conf_full[i],2);
      conf_linfError_f = (conf_linfError_f < conf_full[i]) ? conf_full[i] : conf_linfError_f;
    }

    conf_l2Error_f = sqrt(conf_l2Error_f);

    for(int i=0; i<colL; i++) {
      conf_l1Error_r += fabs(conf_reduced[i]);
      conf_l2Error_r += pow(conf_reduced[i],2);
      conf_linfError_r = (conf_linfError_r < conf_reduced[i]) ? conf_reduced[i] : conf_linfError_r;
    }

    conf_l2Error_r = sqrt(conf_l2Error_r);

    std::cout << "  === Conforming Errors ===  " << std::endl;
    std::cout << "  L_2 Error full:      " << conf_l2Error_f << std::endl;
    std::cout << "  L_2 Error reduced:      " << conf_l2Error_f << std::endl;
    std::cout << "  L_1 Error full:      " << conf_l1Error_f << std::endl;
    std::cout << "  L_1 Error reduced:      " << conf_l1Error_r << std::endl;
    std::cout << "  L_inf Error full:      " << conf_linfError_f << std::endl;
    std::cout << "  L_inf Error reduced:      " << conf_linfError_r << std::endl;
    std::cout << std::endl;

    std::ofstream fout23("test/conf_vec_full.txt");
    for(int i = 0; i < colL; i++) {
      fout23 << conf_full[i];
      if (i < colL - 1) fout23 << "\n";
    }

    std::ofstream fout24("test/conf_vec_reduced.txt");
    for(int i = 0; i < colL; i++) {
      fout24 << conf_reduced[i];
      if (i < colL - 1) fout24 << "\n";
    }

    // Original Error Estimates
  

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

      
      std::string fileName_l_f = fileName + "_l_f";
      solution_l_f.write_raw(fileName_l_f);

      std::string fileName_l_r = fileName + "_l_r";
      solution_l_r.write_raw(fileName_l_r);      
      

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

      
      std::string fileName_l_f = fileName + "_l_f";
      solution_l_f.write_matlab_mesh(fileName_l_f, 
                    param.output_mesh_numPts_Mixed_x,param.output_mesh_numPts_Mixed_y);

      std::string fileName_l_r = fileName + "_l_r";
      solution_l_r.write_matlab_mesh(fileName_l_r, 
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

  return 0;
}; 
