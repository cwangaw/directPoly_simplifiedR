#include <string>
#include <vector>
#include <cmath>
#include <iostream>

#include "ellipticPDE.h"
#include "fcns.h"
#include "Mesh/baseObjects.h"
#include "Mesh/polyMesh.h"
#include "directSerendipity.h"
#include "parameterData.h"
#include "polyQuadrature.h"
#include "Utilities/monitor.h"
#include "Utilities/debug.h"
#include <complex.h>
#include "lapacke.h"
#include <cblas.h>
using namespace directserendipity;
using namespace polymesh;
using namespace polyquadrature;

double innerProduct(double* a, double*b, int size) {
  double product = 0;
  for (int i = 0; i < size; i++) {
    product += a[i] * b[i];
  }
  return product;
}

double l2Error(double* a, double*b, int size) {
  double error = 0;
  for (int i = 0; i < size; i++) {
    error += (a[i]-b[i]) * (a[i]-b[i]);
  }
  error = sqrt(error);
  return error;
}

// M is the left preconditioner
// We are solving for M^{-1} Ax = M^{-1} b
// Here M^{-1} is acting on vector as a function

// numVertices = indices[0][0] contains the number of vertices
// indices[0][1] indices[0][2] ... indices[0][numVertices] contains the dimension
// of each subspace corresponding to vertices
int noPrec(double* A, double* b, std::vector<int>* indices, int size) { 
  return 0;
}
int addSchPrec(double* A, double* b, std::vector<int>* indices, int size) {

  std::vector<double> sol_vec(size,0);
  double* sol = sol_vec.data();

  int numVertices = indices[0][0];

  //for (int i = 1; i <= numVertices+1; i++) {
  for (int i = 1; i <= numVertices+1; i++) {
    int size_of_subspace = indices[0][i];
    if (size_of_subspace == 0) continue;

    // Assemble matrix A_i and b_i
    std::vector<double> A_i_vec(size_of_subspace*size_of_subspace);
    double* A_i = A_i_vec.data();

    std::vector<double> b_i_vec(size_of_subspace);
    double* b_i = b_i_vec.data();

    for (int iRow = 0; iRow < size_of_subspace; iRow++) {
      // Assemble b_i
      b_i[iRow] = b[indices[i][iRow]];

      for (int iCol = 0; iCol < size_of_subspace; iCol++) {

        // Assemble A_i
        A_i[iRow*size_of_subspace+iCol] = A[indices[i][iRow]*size+indices[i][iCol]];
      }
    }

    // Solve A_i*x_i = b_i and store the result in b_i
    // Solve the matrix, result would be stored in rhs
    lapack_int* ipiv;
    ipiv = (lapack_int*)malloc(size_of_subspace * sizeof(lapack_int));
    int ierr = LAPACKE_dgesv(LAPACK_ROW_MAJOR, size_of_subspace, 1, A_i, size_of_subspace, ipiv, b_i, 1); //A_i updated to be LU
    if(ierr) { // ?? what should we do ???
      std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
    }
    // Assemble b_i (holding the information of x_i) back to sol
    for (int iRow = 0; iRow < size_of_subspace; iRow++) {
      sol[indices[i][iRow]] += b_i[iRow];
    }

  }
  for (int i = 0; i < size; i++) {
    b[i] = sol[i];
  }
  return 0;
};

int jacobiPrec(double* A, double* b, std::vector<int>* indices, int size) {
  for (int i = 0; i < size; i++) {
    b[i] /= A[i*size+i];
  }
  return 0;
}


double leftPCG(double* A, double* b, std::vector<int>* indices, double* sol, int problem_size, double* x_0, int max_iter,
              double tol, int (*leftPrec)(double*, double*, std::vector<int>*, int) = nullptr) {
  /* res = b - A * x_0 */
  std::vector <double> res_vec(problem_size);
  double* res = res_vec.data();

  for (int j = 0; j < problem_size; j++) {
    res[j] = b[j];
  }

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, problem_size, 1,
              problem_size, -1, A, problem_size, x_0, 1,
              1, res, 1);

  /* z = M^{-1} * res */
  std::vector <double> z_vec(problem_size);
  double* z = z_vec.data();

  // We first initialize z with res
  for (int j = 0; j < problem_size; j++) {
    z[j] = res[j];
  }

  if (leftPrec != nullptr) {
    // If left preconditioning function is null, we let it be identity matrix
    // Otherwise, z = M^{-1} * res;
    leftPrec(A, z, indices, problem_size);
  }

  /* p = z */
  std::vector <double> p_vec(problem_size);
  double* p = p_vec.data();

  for (int j = 0; j < problem_size; j++) {
    p[j] = z[j];
  }

  /* sol = x_0 */

  for (int j = 0; j < problem_size; j++) {
    sol[j] = x_0[j];
  }

  /* initialize iteration and error */
  int iter = 0;

  double delta_0 = innerProduct(res,z,problem_size);

  /* initialize terms that would be updated in iterations */
  double delta = delta_0;
  double alpha, beta;

  // Store res and z of last iteration
  std::vector <double> res_old_vec(problem_size);
  double* res_old = res_old_vec.data();

  std::vector <double> z_old_vec(problem_size);
  double* z_old = z_old_vec.data();


  /* iteration */
  while (iter < max_iter && delta > tol * tol * delta_0)
  {
    /* res_old = r */
    for (int i = 0; i < problem_size; i++) {
      res_old[i] = res[i];
    }

    /* z_old = z */
    for (int i = 0; i < problem_size; i++) {
      z_old[i] = z[i];
    }

    /* alpha = dot(res,z) / dot(A*p,p) */
    // We first calculate A*p
    std::vector <double> product_Ap_vec(problem_size);
    double* product_Ap = product_Ap_vec.data();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, problem_size, 1,
              problem_size, 1, A, problem_size, p, 1,
              0, product_Ap, 1);

    alpha = innerProduct(res,z,problem_size) / innerProduct(product_Ap,p,problem_size);


    /* x = x + alpha * p; */
    for (int i = 0; i < problem_size; i++) {
      sol[i] += alpha * p[i];
    }

    /* res = res - alpha * A * p; */
    for (int i = 0; i < problem_size; i++) {
      res[i] -= alpha * product_Ap[i];
    }

    /* z = M^{-1} * res */
    // We first initialize z with res
    for (int i = 0; i < problem_size; i++) {
      z[i] = res[i];
    }
    if (leftPrec != nullptr) {
      // If left preconditioning function is null, we let it be identity matrix
      // Otherwise, z = M^{-1} * res;
      leftPrec(A, z, indices, problem_size);

      /* beta = dot(res,z)/dot(res_old,z_old) */
      beta = innerProduct(res,z,problem_size)/innerProduct(res_old,z_old,problem_size);

      /* p = z + beta * p */
      for (int i = 0; i < problem_size; i++) {
        p[i] = z[i] + beta * p[i];
      }
    } 

    /* delta = res^T * z = res^T * M^{-1} * res */
    delta = innerProduct(res,z,problem_size);

    iter++;
  }
  cout << "Number of iterations: " << iter << endl;
  return delta;
}


int EllipticPDE::solve(Monitor& monitor) {
  monitor(0,"Polynomial Degree = ", parameterDataPtr()->dsSpace.degPolyn());

  ParameterData& param = *parameterDataPtr();



  // TEST SINGLE ELEMENT MESH //////////////////////////////////////////////

  if(false) {
    monitor(0,"Test Single Element Mesh and Space");

    polymesh::PolyMesh s_mesh(parameterDataPtr()->mesh.elementPtr(3));
    std::string fileName = parameterDataPtr()->directory_name;
    fileName += "mesh_e4";
    s_mesh.write_matlab(fileName);

    DirectSerendipity s_dsSpace(8,&s_mesh);
    fileName = parameterDataPtr()->directory_name;
    fileName += "dsSpace_e4";
    s_dsSpace.write_matlab(fileName);
  }
  
  // TEST PROCONDITIONING FUNCTIONS ////////////////////////////////////////
  if(true) {
    double A[4] = {2, 1, 1, 2};
    double b[2] = {4, 5};
    double x_0[2] = {0, 0};
    double sol[2] = {0, 0};
    int problem_size = 2;
    int max_iter = 1000;
    double tol = pow(10,-6);


    double delta = leftPCG(A, b, nullptr, sol, problem_size, x_0, max_iter,
              tol, jacobiPrec);

    cout << "Error of PCG testing problem: " << delta << endl;
  }

  // TEST BASIS FUNCTIONS //////////////////////////////////////////////////

  if(true) {
    monitor(0,"\nTest basis functions\n");

    DirectSerendipityArray u(&(parameterDataPtr()->dsSpace));
    for(int i=0; i<u.size(); i++) {
      // Set the coefficient of i-th basis function to be zero
      u[i]=0;

      // Get the coordinate (x,y) of the i-th node
      // Note that the i-th basis function is 1 at this node, and preambly 0 
      // at all the other nodes, unless modifications are made to codes 
      double x = parameterDataPtr()->dsSpace.nodePtr(i)->val(0);
      double y = parameterDataPtr()->dsSpace.nodePtr(i)->val(1);

      // Change the coefficient of the basis function corresponding
      // to the node (1/3,0) to 1
      if (abs(x-(double)1/3)<1e-6 && abs(y-0)<1e-6) u[i]=1;
    }
    
    monitor(1,"Write Array");

    std::string fileName = parameterDataPtr()->directory_name;
    fileName += "basis_mesh";
    std::string fileName_grad = parameterDataPtr()->directory_name;
    fileName_grad += "basis_grad_mesh";
    u.write_matlab_mesh(fileName,fileName_grad,101,101);
  }

  // TEST QUADRATURE ///////////////////////////////////////////////////////

  if(false) {
    monitor(1,"Test Quadrature Rules");
    testPolyQuadrature(&(parameterDataPtr()->mesh),1e-6);
    return 0;
  }
  
  // SOLVE THE PDE ///////////////////////////////////////////////////////////

  monitor(0,"\nSolve the PDE\n");
  
  DirectSerendipityArray solution(&(param.dsSpace));

  // Correct for BCSs: If node i is an interior node, i-th element of index_correction 
  // would give its index in our vector without BC nodes. If node i is a BC node, the  
  // i-th element would give -1.
  std::vector<int> index_correction(param.dsSpace.nNodes());
  int nn = 0;
  for (int i = 0; i < param.dsSpace.nNodes(); i++) {
    if (param.dsSpace.bcType(i) == BCType::interior) {
      index_correction[i] = nn;
      nn++;
    } else {
      index_correction[i] = -1;
    }
  }

  std::vector<double> mat_vector(nn*nn,0);
  double* mat = mat_vector.data();
  std::vector<double> rhs_vector(nn,0);
  double* rhs = rhs_vector.data();

  // quadrature points
  polyquadrature::PolyQuadrature quadRule(13,param.refinement_level);

  monitor(1,"Matrix and RHS Assembly"); ////////////////////////////////////////

  for(int iElement=0; iElement<param.mesh.nElements(); iElement++) {
    DirectSerendipityFE* fePtr = param.dsSpace.finiteElementPtr(iElement);
    
    quadRule.setElement(fePtr->elementPtr());

    // test 

    fePtr->initBasis(quadRule.pts(), quadRule.num());

    // Local matrix and rhs
    int nn_loc = fePtr->nNodes();
    std::vector<double> mat_loc(nn_loc*nn_loc,0), rhs_loc(nn_loc,0);

    // Determine local to global map
    int node_loc_to_gbl[nn_loc];
    for(int i=0; i<fePtr->nVertexNodes(); i++) {
      node_loc_to_gbl[i] = index_correction[fePtr->vertexNodePtr(i)->nodeIndex()];
    }
    for(int i=0; i<fePtr->nEdgeNodes(); i++) {
      node_loc_to_gbl[i + fePtr->nVertexNodes()] = index_correction[fePtr->edgeNodePtr(i)->nodeIndex()]; 
    }
    for(int i=0; i<fePtr->nCellNodes(); i++) {
      node_loc_to_gbl[i + fePtr->nVertexNodes() + fePtr->nEdgeNodes()]
	      = index_correction[fePtr->cellNodePtr(i)->nodeIndex()];
    }

    // Matrix and rhs assembly over elements
    for(int iPt=0; iPt<quadRule.num(); iPt++) {
      double x = quadRule.pt(iPt)[0];
      double y = quadRule.pt(iPt)[1];

      double valA = coefA(x,y);
      Tensor1 valB, valC;
      Tensor2 valD;
      coefB(x,y,valB); coefC(x,y,valC); coefD(x,y,valD);

      // Local interactions	      
      for(int jNode=0; jNode<nn_loc; jNode++) {
        if (node_loc_to_gbl[jNode] == -1) continue;
	      double valj = fePtr->basis(jNode, iPt); 
        Tensor1 gradValj = fePtr->basisGrad(jNode, iPt); 

        for(int iNode=0; iNode<nn_loc; iNode++) {
          //In case "shape functions" are not delta_{i,j} for BC nodes one day
          //if on BC, we use nodal basis function, otherwise we use shape functions

          double vali = fePtr->basis(iNode, iPt);
          Tensor1 gradVali = fePtr->basisGrad(iNode, iPt);

          if (node_loc_to_gbl[iNode] != -1) {
            // +(a N_i,N_j);
	          mat_loc[iNode + nn_loc*jNode] += valA*vali*valj*quadRule.wt(iPt);
            // +(c dot grad(N_i), N_j)
            mat_loc[iNode + nn_loc*jNode] += (valC*gradVali)*valj*quadRule.wt(iPt);
            // +(b N_i, grad(N_j))
            mat_loc[iNode + nn_loc*jNode] += vali*(valB*gradValj)*quadRule.wt(iPt);
            // +(D grad(N_i), grad(N_j))
            mat_loc[iNode + nn_loc*jNode] += ((valD*gradVali)*gradValj)*quadRule.wt(iPt);
          }
          
          //rhs
          Node node_i = *(fePtr->nodePtr(iNode));
          if (node_loc_to_gbl[iNode] == -1) {
	          double bcVali = bcVal(node_i[0],node_i[1]);
            // -(a g(x_i)*N^{BC}_i, N_j)
            rhs_loc[jNode] -= valA*vali*valj*bcVali*quadRule.wt(iPt); 
            // -(c dot g(x_i)*grad(N^{BC}_i), N_j)
            rhs_loc[jNode] -= (valC*gradVali)*valj*bcVali*quadRule.wt(iPt);
            // -(b g(x_i)*N^{BC}_i, grad(N_j))
            rhs_loc[jNode] -= vali*(valB*gradValj)*bcVali*quadRule.wt(iPt);
            // -(D g(x_i)*grad(N^{BC}_i), grad(N_j))
            rhs_loc[jNode] -= ((valD*gradVali)*gradValj)*bcVali*quadRule.wt(iPt);
          }
	}
        // +(f, N_j)
	rhs_loc[jNode] += sourceVal(quadRule.pt(iPt)[0],quadRule.pt(iPt)[1])*valj*quadRule.wt(iPt);
      }
    }

    // Map local matrix and rhs to global
    
    for(int jNode=0; jNode<nn_loc; jNode++) {
      if (node_loc_to_gbl[jNode] != -1) {
        for(int iNode=0; iNode<nn_loc; iNode++) {
          if (node_loc_to_gbl[iNode] != -1) {
            mat[node_loc_to_gbl[iNode] + nn*node_loc_to_gbl[jNode]]
              += mat_loc[iNode + nn_loc*jNode];
          }
        }
      }
    }

    for(int iNode=0; iNode<nn_loc; iNode++) {
      if (node_loc_to_gbl[iNode] != -1) {
        rhs[node_loc_to_gbl[iNode]] += rhs_loc[iNode];
      }
    }
  }

  // Assemble  indices

  // numVertices = indices[0][0] contains the number of vertices
  // indices[0][1] indices[0][2] ... indices[0][numVertices] contains the dimension
  // of each subspace corresponding to vertices
  // indices[0][numVertices+1] contains the dimension of low-order matrix


  std::vector<int>* indices = new std::vector<int>[param.dsSpace.nVertexNodes()+2];

  

  indices[0].push_back(param.dsSpace.nVertexNodes());

  // We first assemble local matrices
  for (int iVertex = 0; iVertex < param.dsSpace.nVertexNodes(); iVertex++) {
    int index = iVertex + 1;
    int local_dim = 0;
    // First we add vertex node (if it is interior)
    if (!param.mesh.vertexPtr(iVertex)->isOnBoundary()) {
      indices[index].push_back(index_correction[param.dsSpace.meshVertexToNodeIndex(iVertex)]);
      local_dim++;
    }

    // Now we add edge nodes (if it is interior)
    std::vector<int> connected_edges;
    param.mesh.nbrEdgesOfVertex(iVertex,connected_edges);
    for (unsigned int iEdge = 0; iEdge < connected_edges.size(); iEdge++) {
      if (!param.mesh.edgePtr(connected_edges[iEdge])->isOnBoundary()) {
        // If it is an interior edge, we would count its nodes
        int starting_index = param.dsSpace.meshEdgeToFirstNodeIndex(connected_edges[iEdge]);
        for (int nEdgeNode = 0; nEdgeNode < param.dsSpace.degPolyn()-1; nEdgeNode++){
          indices[index].push_back(index_correction[starting_index]);
          local_dim++;
          starting_index++;
        }
      }
    }

    // Now we add interior nodes
    std::vector<int> connected_elements;
    param.mesh.nbrElementsOfVertex(iVertex,connected_elements);
    for (unsigned int iElement = 0; iElement < connected_elements.size(); iElement++) {
      int nVertices = param.mesh.elementPtr(connected_elements[iElement])->nVertices();
      if (param.dsSpace.degPolyn() >= nVertices) {
        int numInteriorNodes = (param.dsSpace.degPolyn()-nVertices+2) * (param.dsSpace.degPolyn()-nVertices+1) / 2;
        int starting_index = param.dsSpace.meshElementToFirstNodeIndex(connected_elements[iElement]);
        for (int nInteriorNode = 0; nInteriorNode < numInteriorNodes; nInteriorNode++) {
          indices[index].push_back(index_correction[starting_index]);
          local_dim++;
          starting_index++;
        }
      }
    }

    // Now we store local_dim
    indices[0].push_back(local_dim);
    
  }

  // We now assemble the low-order matrix
  int lom_dim = 0; // dimension for low order matrix

  // First we add vertex nodes (if it is interior)
  for (int iVertex = 0; iVertex < param.dsSpace.nVertexNodes(); iVertex++) {
    if (!param.mesh.vertexPtr(iVertex)->isOnBoundary()) {
      indices[param.dsSpace.nVertexNodes()+1].push_back(index_correction[param.dsSpace.meshVertexToNodeIndex(iVertex)]);
      lom_dim++;
    }
  }

  // Now we add interior nodes
  
  for (int iElement = 0; iElement < param.mesh.nElements(); iElement++) {
      int nVertices = param.mesh.elementPtr(iElement)->nVertices();

      // If no vertex is in the interior, we skip the element
      int numInteriorVertex = 0;
      for (int iVertex = 0; iVertex < nVertices; iVertex++) {
        if (!param.mesh.elementPtr(iElement)->vertexPtr(iVertex)->isOnBoundary())
        numInteriorVertex += 1;
      }
      if (numInteriorVertex == 0) continue;

      if (param.dsSpace.degPolyn() >= nVertices) {
        int numInteriorNodes = (param.dsSpace.degPolyn()-nVertices+2) * (param.dsSpace.degPolyn()-nVertices+1) / 2;
        int starting_index = param.dsSpace.meshElementToFirstNodeIndex(iElement);
        for (int nInteriorNode = 0; nInteriorNode < numInteriorNodes; nInteriorNode++) {
          indices[param.dsSpace.nVertexNodes()+1].push_back(index_correction[starting_index]);
          lom_dim++;
          starting_index++;
        }
      }
    }
    

  indices[0].push_back(lom_dim);

  // Direct Solver would update mat to be LU, and rhs to be the solution
  // We need to store another set of mat and rhs for iterative solver
  std::vector<double> mat_iter_vector(nn*nn,0);
  double* mat_iter = mat_iter_vector.data();
  std::vector<double> rhs_iter_vector(nn,0);
  double* rhs_iter = rhs_iter_vector.data();

  for(int j=0; j<nn; j++) {
    for(int i=0; i<nn; i++) {
      mat_iter[i + nn*j] = mat[i + nn*j];
    }
  }

  for(int i=0; i<nn; i++) {
    rhs_iter[i] = rhs[i];
  }
  // Test

  std::ofstream fout("test/matrix.txt");
  fout.precision(24);
  for(int j=0; j<nn; j++) {
    for(int i=0; i<nn; i++) {
      fout << mat[i + nn*j] << "\t";
    }
    if (j < nn - 1) fout << "\n";
  }


  std::ofstream rout("test/rhs.txt");
  rout.precision(24);
  for(int i=0; i<nn; i++) {
    rout << rhs[i];
    if (i < nn - 1) rout << "\n";
  }


  monitor(1,"===============================================");
  monitor(1,"===Solution of linear system (DIRECT SOLVER)==="); ////////////////////////////////////////
  monitor(1,"===============================================");
  
  //Solve the matrix, result would be stored in rhs
  lapack_int* ipiv; char norm = 'I'; 
  ipiv = (lapack_int*)malloc(nn * sizeof(lapack_int));
  double anorm = LAPACKE_dlange(LAPACK_ROW_MAJOR, norm, nn, nn, mat, nn);
  int ierr = LAPACKE_dgesv(LAPACK_ROW_MAJOR, nn, 1, mat, nn, ipiv, rhs, 1); //mat updated to be LU
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }
  double rcond = 0;
  ierr = LAPACKE_dgecon(LAPACK_ROW_MAJOR, norm, nn, mat, nn, anorm, &rcond);
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }
  rcond = 1/rcond;

  // Calculate the residual of direct solver
  std::vector<double> Ax_direct_vec(nn,0);
  double* Ax_direct = Ax_direct_vec.data();

  for (int j = 0; j < nn; j++) {
    Ax_direct[j] = 0;
  }

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nn, 1,
              nn, 1, mat_iter, nn, rhs, 1,
              0, Ax_direct, 1);


  cout << "l2 norm of residual: " << l2Error(Ax_direct, rhs_iter, nn) << endl;


  std::ofstream sout("test/solution.txt");
  sout.precision(6);
  for(int i=0; i<nn; i++) {
    sout << rhs[i];
    if (i < nn - 1) sout << "\n";
  }


  //Calculate inf condition number
 
  std::cout << "\tNorm Format:\t" << norm << std::endl;
  std::cout << "\tNorm of mat:\t" << anorm << std::endl;
  std::cout << "\tCond number:\t" << rcond << std::endl;


  for(int i=0; i<solution.size(); i++) {
    if(index_correction[i] == -1) {
      double x = param.dsSpace.nodePtr(i)->val(0);
      double y = param.dsSpace.nodePtr(i)->val(1);
      solution[i] = bcVal(x,y);
    } else {
      solution[i] = rhs[index_correction[i]];
    }
  }

  if(param.output_soln_DS_format > 0) {

    monitor(1,"Write Chunkiness Parameter"); //////////////////////////////////////////////////

    switch(param.output_soln_DS_format) {
    case 1: {
      break;
    }
    case 2: {
      std::string fileName(param.directory_name);
      fileName += "mesh_one_over_chunkiness_parameter";
      solution.write_matlab_mesh_one_over_chunk(fileName,
				 param.output_mesh_numPts_DS_x,param.output_mesh_numPts_DS_y);
      break;
    }
    }
  
    monitor(1,"Write Solution"); //////////////////////////////////////////////////

    switch(param.output_soln_DS_format) {
    case 1: {
      std::string fileName(param.directory_name);
      fileName += "solution_raw";
      solution.write_raw(fileName);
      break;
    }
    case 2: {
      std::string fileName(param.directory_name);
      fileName += "solution_mesh";
      std::string fileNameGrad(param.directory_name);
      fileNameGrad += "solution_grad_mesh";
      solution.write_matlab_mesh(fileName,fileNameGrad,
				 param.output_mesh_numPts_DS_x,param.output_mesh_numPts_DS_y);
      break;
    }
    }
  }

  if(trueSolnKnown()) {
    monitor(0,"\nError estimate\n"); ///////////////////////////////////////////////
  
    double h = param.dsSpace.mesh()->maxElementDiameter();
    double maxChunk = param.dsSpace.mesh()->maxChunkParam();
    double minChunk = param.dsSpace.mesh()->minChunkParam();
    double averageChunk = param.dsSpace.mesh()->averageChunkParam();
    
    
    double l2Error = 0, l2GradError = 0, l2Norm = 0, l2GradNorm = 0;
    solution.l2normError(l2Error, l2GradError, l2Norm, l2GradNorm, param.refinement_level*0, trueSoln, trueGradSoln);
    
    std::cout << "  Max Element Diameter h:  " << h << std::endl;
    std::cout << "  Size of Matrix: " << nn << std::endl;
    std::cout << "  Max Chunkiness Parameter:  " << maxChunk << std::endl;
    std::cout << "  Min Chunkiness Parameter:  " << minChunk << std::endl;
    std::cout << "  Average Chunkiness Parameter:  " << averageChunk << std::endl;
    std::cout << "  L_2 Error:      " << l2Error << std::endl;
    std::cout << "  L_2 Grad Error: " << l2GradError << std::endl;
    std::cout << std::endl;
    std::cout << "  Relative L_2 Error:      " << l2Error/l2Norm << std::endl;
    std::cout << "  Relative L_2 Grad Error: " << l2GradError/l2GradNorm << std::endl;
    std::cout << std::endl;

    if(param.output_soln_DS_format > 0) {
      monitor(1,"Write True Solution"); ////////////////////////////////////////////

      DirectSerendipityArray u(&(param.dsSpace));

      for(int i=0; i<u.size(); i++) {
        double x = param.dsSpace.nodePtr(i)->val(0);
        double y = param.dsSpace.nodePtr(i)->val(1);
        u[i] = trueSoln(x,y);
      }

      switch(param.output_soln_DS_format) {
      case 1: {
        std::string fileName(param.directory_name);
        fileName += "true_solution_raw";
        u.write_raw(fileName);
        break;
      }
      case 2: {
        std::string fileName(param.directory_name);
        fileName += "true_solution_mesh";
        std::string fileNameGrad(param.directory_name);
        fileNameGrad += "true_solution_grad_mesh";
        u.write_matlab_mesh(fileName,fileNameGrad,
        param.output_mesh_numPts_DS_x,param.output_mesh_numPts_DS_y);
        break;
      }
      }


      monitor(1,"Write Error"); ////////////////////////////////////////////

      switch(param.output_soln_DS_format) {
      case 1: {
        break;
      }
      case 2: {
        std::string fileName(param.directory_name);
        fileName += "solution_mesh_error";
        std::string fileNameGrad(param.directory_name);
        fileNameGrad += "solution_mesh_grad_error";
        solution.write_matlab_mesh_error(fileName,
        param.output_mesh_numPts_DS_x,param.output_mesh_numPts_DS_y, trueSoln);
        solution.write_matlab_mesh_grad_error(fileNameGrad,
        param.output_mesh_numPts_DS_x,param.output_mesh_numPts_DS_y, trueGradSoln);
        fileName += "_on_element";
        fileNameGrad += "_on_element";
        solution.write_matlab_mesh_error_on_element(fileName,
        param.output_mesh_numPts_DS_x,param.output_mesh_numPts_DS_y, param.refinement_level*0, trueSoln);
        solution.write_matlab_mesh_grad_error_on_element(fileNameGrad,
        param.output_mesh_numPts_DS_x,param.output_mesh_numPts_DS_y, param.refinement_level*0, trueGradSoln);
        break;
      }
      }
    }
  }  




  monitor(1,"==================================================");
  monitor(1,"===Solution of linear system (ITERATIVE SOLVER)==="); ////////////////////////////////////////
  monitor(1,"==================================================");
  
  std::vector<double> sol_iter_vector(nn,0);
  double* sol_iter = sol_iter_vector.data();

  int max_iter = 2*nn;
  double tol = 1e-20;

  //Solve the matrix, result would be stored in sol
  double delta = leftPCG(mat_iter, rhs_iter, indices, sol_iter, nn, sol_iter, max_iter,
              tol, addSchPrec);

  std::ofstream sitout("test/solution_iteration.txt");
  sitout.precision(6);
  for(int i=0; i<nn; i++) {
    sitout << sol_iter[i];
    if (i < nn - 1) sitout << "\n";
  }

  cout << "Error of Preconditioned Conjugate Gradient Iteration: " << delta << endl;

  // Calculate the residual of iterative solver
  std::vector<double> Ax_iter_vec(nn,0);
  double* Ax_iter = Ax_iter_vec.data();

  for (int j = 0; j < nn; j++) {
    Ax_iter[j] = 0;
  }

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nn, 1,
              nn, 1, mat_iter, nn, sol_iter, 1,
              0, Ax_iter, 1);


  cout << "l2 norm of residual: " << l2Error(Ax_iter, rhs_iter, nn) << endl;

  // Here we rewrite the DirectSerendipity array "solution" to be 
  // the result of preconditioned iterative method
  for(int i=0; i<solution.size(); i++) {
    if(index_correction[i] == -1) {
      double x = param.dsSpace.nodePtr(i)->val(0);
      double y = param.dsSpace.nodePtr(i)->val(1);
      solution[i] = bcVal(x,y);
    } else {
      solution[i] = sol_iter[index_correction[i]];
    }
  }

  if(param.output_soln_DS_format > 0) {
  
    monitor(1,"Write Solution"); //////////////////////////////////////////////////

    switch(param.output_soln_DS_format) {
    case 1: {
      std::string fileName(param.directory_name);
      fileName += "solution_iter_raw";
      solution.write_raw(fileName);
      break;
    }
    case 2: {
      std::string fileName(param.directory_name);
      fileName += "solution_iter_mesh";
      std::string fileNameGrad(param.directory_name);
      fileNameGrad += "solution_grad_iter_mesh";
      solution.write_matlab_mesh(fileName,fileNameGrad,
				 param.output_mesh_numPts_DS_x,param.output_mesh_numPts_DS_y);
      break;
    }
    }
  }

  if(trueSolnKnown()) {
    monitor(0,"\nError estimate\n"); ///////////////////////////////////////////////

    double l2Error = 0, l2GradError = 0, l2Norm = 0, l2GradNorm = 0;
    solution.l2normError(l2Error, l2GradError, l2Norm, l2GradNorm, param.refinement_level*0, trueSoln, trueGradSoln);

    std::cout << "  L_2 Error:      " << l2Error << std::endl;
    std::cout << "  L_2 Grad Error: " << l2GradError << std::endl;
    std::cout << std::endl;
    std::cout << "  Relative L_2 Error:      " << l2Error/l2Norm << std::endl;
    std::cout << "  Relative L_2 Grad Error: " << l2GradError/l2GradNorm << std::endl;
    std::cout << std::endl;

    if(param.output_soln_DS_format > 0) {
      monitor(1,"Write Error"); ////////////////////////////////////////////

      switch(param.output_soln_DS_format) {
      case 1: {
        break;
      }
      case 2: {
        std::string fileName(param.directory_name);
        fileName += "solution_iter_mesh_error";
        std::string fileNameGrad(param.directory_name);
        fileNameGrad += "solution_iter_mesh_grad_error";
        solution.write_matlab_mesh_error(fileName,
        param.output_mesh_numPts_DS_x,param.output_mesh_numPts_DS_y, trueSoln);
        solution.write_matlab_mesh_grad_error(fileNameGrad,
        param.output_mesh_numPts_DS_x,param.output_mesh_numPts_DS_y, trueGradSoln);
        fileName += "_on_element";
        fileNameGrad += "_on_element";
        solution.write_matlab_mesh_error_on_element(fileName,
        param.output_mesh_numPts_DS_x,param.output_mesh_numPts_DS_y, param.refinement_level*0, trueSoln);
        solution.write_matlab_mesh_grad_error_on_element(fileNameGrad,
        param.output_mesh_numPts_DS_x,param.output_mesh_numPts_DS_y, param.refinement_level*0, trueGradSoln);
        break;
      }
      }
    }
  }  

  delete[] indices;
  return 0;
} 
