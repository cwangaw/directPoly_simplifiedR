#include <cmath>
#include <iostream>
#include <vector>

#include <complex.h>
#include "lapacke.h"
#include <stdio.h>
#include <assert.h>

using namespace std;

#include "Utilities/debug.h"
#include "Mesh/polyMesh.h"
using namespace polymesh;
#include "directSerendipity.h"
using namespace directserendipity;

////////////////////////////////////////////////////////////////////////////////
// Class DirectSerendipityFE

void DirectSerendipityFE::
         set_directserendipityfe(DirectSerendipity* dsSpace, PolyElement* element,
				 const std::vector<Node*>& vertexNodes, const std::vector<Node*>& edgeNodes,
				 const std::vector<Node*>& cellNodes) {
      my_ds_space = dsSpace;
      my_poly_element = element;

      num_vertices = element->nVertices();
      polynomial_degree = my_ds_space->degPolyn();
      deg_cell_poly = polynomial_degree - num_vertices;
      num_cell_nodes = (deg_cell_poly >=0) ? (deg_cell_poly+2)*(deg_cell_poly+1)/2 : 0;
	    num_nodes = num_vertices*polynomial_degree + num_cell_nodes;

      the_vertex_nodes = vertexNodes;
      the_edge_nodes = edgeNodes;
      the_cell_nodes = cellNodes;

      // Set up higher order element (if nescessary)
      if(polynomial_degree<num_vertices-2) {
	      if(one_element_mesh) delete one_element_mesh;
	      one_element_mesh = new polymesh::PolyMesh(my_poly_element);
	      if(high_order_ds_space) delete high_order_ds_space;
	      high_order_ds_space = new DirectSerendipity(num_vertices-2,one_element_mesh);
      }

      /*set up basis*/
      const int matrix_dim = polynomial_degree-1;
      int dim_of_x = num_vertices * (polynomial_degree - 1) * matrix_dim;
      x_vector_for_all.resize(dim_of_x);
      double* x_vector_for_all_array = x_vector_for_all.data();

      std::vector<double> matrix_A_vector(matrix_dim * matrix_dim);
      double* matrix_A = matrix_A_vector.data();
      std::vector<double> b_vector(matrix_dim);
      double* b = b_vector.data();

      for (int nEdge = 0; nEdge < num_vertices; nEdge++){
        for (int jNode = 0; jNode < (polynomial_degree-1); jNode++){
          //Determine each element of A
          for (int row = 0; row < matrix_dim; row++) {
            for (int col = 0; col < matrix_dim; col++) {
	          if (col <= polynomial_degree - num_vertices + 1 ) {
	            //coefficients for alpha_col at \x_{e,nEdge,row}
	            double entry = pow(projToEdge(nEdge,*edgeNodePtr(nEdge,row)),col);
	            for (int m=0; m<num_vertices; m++) {
	              if (m == nEdge || m == (nEdge + 1) % num_vertices || m == (nEdge + num_vertices - 1) % num_vertices) { continue; } else {
	                entry *= my_poly_element -> edgePtr(m) -> lambda(*edgeNodePtr(nEdge,row));
	              }
	            }
	            matrix_A[row*matrix_dim+col] = entry;
	          } else {
	            //coefficients for beta_(nEdge+col-polynomial_degree+num_vertices) at \x_{e,nEdge,row}
	            int index_for_beta = (nEdge+col-polynomial_degree+num_vertices) % num_vertices;
	            double entry = pow(lambda_supp(nEdge,index_for_beta,*edgeNodePtr(nEdge,row)),
			                          polynomial_degree-num_vertices+2);
	            for (int m=0; m<num_vertices; m++) {
	              if (m == nEdge || m == (nEdge + 1) % num_vertices
		                || m == (nEdge + num_vertices - 1) % num_vertices || m == index_for_beta) { continue; } else {
	                entry *= my_poly_element -> edgePtr(m) -> lambda(*edgeNodePtr(nEdge,row));
	              }
	            }
	            matrix_A[row*matrix_dim + col] = entry;
	          }
            }
          }
    
          //Determine each element of b
          for (int row = 0; row < matrix_dim; row++) {
            if (row == jNode) {
	            b[row] = 1 / my_poly_element ->edgePtr((nEdge + 1) % num_vertices) ->lambda(*edgeNodePtr(nEdge,row));
	            b[row] /= my_poly_element ->edgePtr((nEdge + num_vertices - 1) % num_vertices) ->lambda(*edgeNodePtr(nEdge,row));
            } else { b[row] = 0; }
          }
    
          //It remains to solve Ax=b

          lapack_int* ipiv;
          ipiv = (lapack_int*)malloc(matrix_dim * sizeof(lapack_int));
          int ierr = LAPACKE_dgesv(LAPACK_ROW_MAJOR, matrix_dim, 1, matrix_A, matrix_dim, ipiv, b, 1);
          free(ipiv);
          if(ierr) {
	          std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
          }

          //Pass value of x (updated b) to proper location of x_vector_for_all_array
          int global_edge_node_index = nEdge * (polynomial_degree -1) + jNode;
          int begin_position = global_edge_node_index * matrix_dim;
          for (int i = begin_position; i < begin_position + matrix_dim; i++) {
            x_vector_for_all_array[i] = b[i-begin_position];
          }
        }
      }
      


};

DirectSerendipityFE::~DirectSerendipityFE() {
  if(high_order_ds_space) delete high_order_ds_space;
  if(one_element_mesh) delete one_element_mesh;
  if(value_n) delete[] value_n;
  if (gradvalue_n) delete[] gradvalue_n;
}

// Helper functions for DirectSerendipityFE //////////////////////////////////////////
Node* DirectSerendipityFE::nodePtr(int i) const {
  assert(i>=0 && i<num_nodes);
  if (i < nVertexNodes()) {
    return vertexNodePtr(i);
  } else if ( i >= nVertexNodes() && i < (num_vertices*polynomial_degree) ) {
    return edgeNodePtr(i - nVertexNodes());
  } else {
    return cellNodePtr(i - num_vertices*polynomial_degree);
  }
};

int DirectSerendipityFE::mapI(int k, int deg_cell_poly) const {
    int j = mapJ(k, deg_cell_poly);
    return k - mapK(0, j, deg_cell_poly);
};

int DirectSerendipityFE::mapJ(int k, int deg_cell_poly) const {
    return my_ds_space->map_j_array[deg_cell_poly][k];
}

int DirectSerendipityFE::mapK(int i, int j, int deg_cell_poly) const {
    return num_cell_nodes - (deg_cell_poly + 2 - j) * (deg_cell_poly + 1 - j) / 2 + i;
};

void DirectSerendipityFE::projToEdge(int iEdge, double x, double y, double& result, Tensor1& gradresult) const {
  polymesh::OrientedEdge* e = my_poly_element->edgePtr(iEdge);
  Tensor1 tau(e->tangent());
  result = (x - e->vertexPtr(0)->val(0))*tau[0] + (y - e->vertexPtr(0)->val(1))*tau[1];
  gradresult = tau;
};

double DirectSerendipityFE::projToEdge(int iEdge, double x, double y) const {
  double result;
  Tensor1 gradresult;
  projToEdge(iEdge, x, y, result, gradresult);
  return result;
}


void directserendipity::lambda_for_two_points(const Point& pt0, const Point& pt1, double x, double y, double& result, Tensor1& gradresult) {
  Tensor1 tangent(pt1 - pt0);
  double edge_length = tangent.norm();
  if( edge_length != 0 ) tangent /= edge_length;
  
  Tensor1 normal(tangent[1],-tangent[0]);
  result = (pt0[0]-x)*normal[0] + (pt0[1]-y)*normal[1];
  gradresult = -1*normal;
}

double directserendipity::lambda_for_two_points(const Point& pt0, const Point& pt1,
						  double x, double y) {
  double result;
  Tensor1 gradresult;
  lambda_for_two_points(pt0, pt1, x, y, result, gradresult);
  return result;
}


const double* DirectSerendipityFE::get_result_of_coefficients(int i, int j) const{
  int global_edge_node_index = i * (polynomial_degree -1) + j;
  int begin_position = global_edge_node_index * (polynomial_degree -1);
  const double* pointer_to_begin_position= &x_vector_for_all.data()[begin_position];
  return pointer_to_begin_position;
}

//Supplemental functions

void DirectSerendipityFE::r_supp(int i, int j, const Point& p, double& result, Tensor1& gradresult) const {
  OrientedEdge* e_i = my_poly_element->edgePtr(i);
  OrientedEdge* e_j = my_poly_element->edgePtr(j);

  double r_simple = (e_i -> lambda(p) - e_j -> lambda(p)) / (e_i -> lambda(p) + e_j -> lambda(p));
  result = 0.5 * (1-r_simple);
  
  double sum = e_i -> lambda(p) + e_j -> lambda(p);
  gradresult(0) = (e_j->dLambda(0) * sum - (e_i->dLambda(0) + e_j->dLambda(0)) * e_j -> lambda(p)) / (sum * sum);
  gradresult(1) = (e_j->dLambda(1) * sum - (e_i->dLambda(1) + e_j->dLambda(1)) * e_j -> lambda(p)) / (sum * sum);
}

double DirectSerendipityFE::r_supp(int i, int j, const Point& p) const {
  double result;
  Tensor1 gradresult;
  r_supp(i, j, p, result, gradresult);
  return result;
}


void DirectSerendipityFE::lambda_supp (int i, int j, const Point& p, double& result, Tensor1& gradresult) const {
  if (i > j) {
    int m = i;
    i = j;
    j = m;
  }

  Vertex v_j = Vertex(vertexNodePtr(j));
  Vertex v_im1 = Vertex(vertexNodePtr((i+num_vertices-1)%num_vertices));
  Vertex v_i = Vertex(vertexNodePtr(i));
  Vertex v_jm1 = Vertex(vertexNodePtr((j+num_vertices-1)%num_vertices));
  Edge e1(&v_j,&v_im1);
  Edge e2(&v_i,&v_jm1);
  Tensor1 normal_vec = e1.normal()-e2.normal();

  result = (e1.lambda(p)-e2.lambda(p))/normal_vec.norm();
  gradresult(0) = (e1.dLambda(0)-e2.dLambda(0))/normal_vec.norm();
  gradresult(1) = (e1.dLambda(1)-e2.dLambda(1))/normal_vec.norm();
}

double DirectSerendipityFE::lambda_supp (int i, int j, const Point& p) const{
  double result;
  Tensor1 gradresult;
  lambda_supp(i, j, p, result, gradresult);
  return result;
}


void DirectSerendipityFE::phi_k_l (int k, int l, const Point& p, double& result, Tensor1& gradresult) const {
  assert( fabs(k-l) >= 2 && fabs(k-l) <= num_vertices-2 );
  //Define and initialize variable for grad
  int num_term = num_vertices+2;
  std::vector<double> term_grad_coef_part(num_term,1);
  std::vector<Tensor1> term_grad(num_term, Tensor1(0,0));
  result = 1;
  gradresult = Tensor1(0,0);
  double new_part; Tensor1 new_grad;

  for (int m=0; m<num_vertices; m++) {
    if (m==k || m==l) { continue; } 
    else {
      new_part = my_poly_element->edgePtr(m)->lambda(p);
      result *= new_part;
      for (int n = 0; n < num_term; n++) { term_grad_coef_part[n] *= (n==m) ? 1 : new_part; }
      term_grad[m] = my_poly_element->edgePtr(m)->dLambda(); 
    }
  }
  double lambda_supp_result; Tensor1 dlambda_supp_result;
  lambda_supp(k,l,p,lambda_supp_result,dlambda_supp_result);
  new_part = pow(lambda_supp_result,polynomial_degree-num_vertices+2);
  result *= new_part;
  for (int n = 0; n < num_term; n++) { term_grad_coef_part[n] *= (n==num_vertices) ? 1 : new_part; }
  if ((polynomial_degree-num_vertices+2)!=0) {
    term_grad[num_vertices] = ((polynomial_degree-num_vertices+2) * pow(lambda_supp_result,polynomial_degree-num_vertices+1)) * dlambda_supp_result;
  }


  r_supp(k,l,p,new_part,new_grad);
  result *= new_part;
  for (int n = 0; n < num_term; n++) { term_grad_coef_part[n] *= (n==num_vertices+1) ? 1 : r_supp(k,l,p); }
  term_grad[num_vertices+1] = new_grad;

  for (int i=0; i<num_term; i++) {
    gradresult += term_grad_coef_part[i] * term_grad[i];
  }
}

double DirectSerendipityFE::phi_k_l (int k, int l, const Point& p) const {
  double result;
  Tensor1 gradresult;
  phi_k_l(k, l, p, result, gradresult);
  return result;
}


//Lagrange polynomial (of degree r) that is 1 on x_nEdge_jNode, 0 on other nodes and the two vertices on edge_nEdge
void DirectSerendipityFE::lagrange(int nEdge, int jNode, const Point& pt, double& result, Tensor1& gradresult) const {
  int num_term = polynomial_degree + 1;
  std::vector<double> term_grad_coef_part(num_term,1);
  std::vector<Tensor1> term_grad(num_term, Tensor1(0,0));
  result = 1;
  gradresult = Tensor1(0,0);
  double new_part;

  //result *= (pt-x_{e,n,k})/(x_{e,n,j}-x_{e,n,k}) for all k != j  
  for (int k=0; k<polynomial_degree-1; k++) {
	  if (k == jNode) continue;
	  else {
	    Point p(*edgeNodePtr(nEdge, jNode));
	    Point zero_node(*edgeNodePtr(nEdge,k));
      double proj_pt;
      Tensor1 dproj_pt;
      projToEdge(nEdge,pt,proj_pt,dproj_pt);
	    new_part =  ( proj_pt - projToEdge(nEdge,zero_node) )
	             / ( projToEdge(nEdge,p)  - projToEdge(nEdge,zero_node) );
      result *= new_part;
      for (int m = 0; m < num_term; m++) { term_grad_coef_part[m] *= (m==k) ? 1 : new_part; }
      term_grad[k] = dproj_pt / ( projToEdge(nEdge,p)  - projToEdge(nEdge,zero_node) );
	  }
  }
  //result *= (pt-x_{v,n})/(x_{e,n,j}-x_{v,n})
  for (int n=nEdge-1+num_vertices; n<=nEdge+num_vertices; n++) {
    Point p(*edgeNodePtr(nEdge, jNode)), zero_node(*vertexNodePtr(n%num_vertices));
    double proj_pt; Tensor1 dproj_pt;
    projToEdge(nEdge,pt,proj_pt,dproj_pt);
    new_part =  (proj_pt-projToEdge(nEdge,zero_node)) / (projToEdge(nEdge,p)-projToEdge(nEdge,zero_node));
    result *= new_part;
    for (int m = 0; m < num_term; m++) { term_grad_coef_part[m] *= (m==polynomial_degree+n-(nEdge+num_vertices)) ? 1 : new_part; }
    term_grad[polynomial_degree+n-(nEdge+num_vertices)] = dproj_pt / (projToEdge(nEdge,p)-projToEdge(nEdge,zero_node));
  }

  for (int i=0; i<num_term; i++) {
    gradresult += term_grad_coef_part[i] * term_grad[i];
  }
};

double DirectSerendipityFE::lagrange(int nEdge, int jNode, const Point& pt) const {
  double result;
  Tensor1 gradresult;
  lagrange(nEdge, jNode, pt, result, gradresult);
  return result;
}

//Lagrange polynomial that is 1 on x_{v,i}, 0 on other nodes including vertex on edge_nEdge
void DirectSerendipityFE::lagrange_v(int i, int nEdge, const Point& pt, double& result, Tensor1& gradresult) const {
  assert( nEdge == i || nEdge == (i+1)%num_vertices );
  int num_term = polynomial_degree;
  std::vector<double> term_grad_coef_part(num_term,1);
  std::vector<Tensor1> term_grad(num_term, Tensor1(0,0));
  result = 1;
  gradresult = Tensor1(0,0);
  double new_part;

  //result *= (pt-x_{e,n,k})/(x_{v,i}-x_{e,n,k}) for all k 
  for (int k=0; k<polynomial_degree-1; k++) {
	  Point p(*vertexNodePtr(i)), zero_node(*edgeNodePtr(nEdge, k));
    double proj_pt; Tensor1 dproj_pt;
    projToEdge(nEdge, pt, proj_pt, dproj_pt);
	  new_part =  (proj_pt-projToEdge(nEdge, zero_node)) / (projToEdge(nEdge, p)-projToEdge(nEdge, zero_node));
    result *= new_part;
    for (int n = 0; n < num_term; n++) { term_grad_coef_part[n] *= (n==k) ? 1 : new_part; }
    term_grad[k] = dproj_pt / (projToEdge(nEdge, p)-projToEdge(nEdge, zero_node));
  }
  //result *= (pt-x_{v,i-1})/(x_{v,i}-x_{v,i-1}) if n=i;
  //result *= (pt-x_{v,i+1})/(x_{v,i}-x_{v,i+1}) if n=i+1;
  Point p(*vertexNodePtr(i));
  Point zero_node;
  if (nEdge == i) {
    zero_node = Point(*vertexNodePtr((i-1+num_vertices)%num_vertices));
  } else { 
    zero_node = Point(*vertexNodePtr((i+1)%num_vertices));
  }
  double proj_pt; Tensor1 dproj_pt;
  projToEdge(nEdge, pt, proj_pt, dproj_pt);
  new_part = (proj_pt-projToEdge(nEdge, zero_node)) / (projToEdge(nEdge, p)-projToEdge(nEdge, zero_node));
  result *= new_part;
  for (int n = 0; n < num_term; n++) { term_grad_coef_part[n] *= (n==polynomial_degree-1) ? 1 : new_part; }
  term_grad[polynomial_degree-1] = dproj_pt / (projToEdge(nEdge, p)-projToEdge(nEdge, zero_node));

  for (int i=0; i<num_term; i++) {
    gradresult += term_grad_coef_part[i] * term_grad[i];
  }
};

double DirectSerendipityFE::lagrange_v(int i, int nEdge, const Point& pt) const {
  double result;
  Tensor1 gradresult;
  lagrange_v(i, nEdge, pt, result, gradresult);
  return result;
}


//varphi before normalization for small r and edge nodes in A_po
void DirectSerendipityFE::varphi_edge_nodes_A_Po(int nEdge, int jNode, const Point& pt, double& result, Tensor1& gradresult) const {
  int num_term = polynomial_degree+1;
  std::vector<double> term_grad_coef_part(num_term,1);
  std::vector<Tensor1> term_grad(num_term, Tensor1(0,0));
  result = 1;
  gradresult = Tensor1(0,0);
  double new_part; Tensor1 new_grad;

  if (nEdge <= polynomial_degree-2 && jNode <= nEdge) {   
    for (int node = 0; node <= nEdge; node++) {
      if (node == jNode) { continue; }
      lambda_for_two_points(*edgeNodePtr(0,0),*edgeNodePtr(nEdge,node),pt,new_part,new_grad);
      result *= new_part;
      for (int n = 0; n < num_term; n++) { term_grad_coef_part[n] *= (n==node) ? 1 : new_part; }
      term_grad[node] = new_grad;
    }
    for (int edge = nEdge+1; edge <= polynomial_degree; edge++) {
      new_part = my_poly_element -> edgePtr(edge) -> lambda(pt);
      result *= new_part;
      for (int n = 0; n < num_term; n++) { term_grad_coef_part[n] *= (n==edge) ? 1 : new_part; }
      term_grad[edge] = my_poly_element -> edgePtr(edge) -> dLambda();
    }

  } 
  else if (nEdge > polynomial_degree-2 && nEdge <= polynomial_degree) {
    for (int node = 0; node < polynomial_degree-1; node++) {
      if (node == jNode) { continue; }
      lambda_for_two_points(*edgeNodePtr(0,0),*edgeNodePtr(nEdge,node),pt,new_part,new_grad);
      result *= new_part;
      for (int n = 0; n < num_term; n++) { term_grad_coef_part[n] *= (n==node) ? 1 : new_part; }
      term_grad[node] = new_grad;
    }
    //Cope with vertices involved for e_{r-1} and e_r
    lambda_for_two_points(*edgeNodePtr(0,0),*vertexNodePtr((nEdge-1+num_vertices)%num_vertices),pt,new_part,new_grad);
    result *= new_part;
    for (int n = 0; n < num_term; n++) { term_grad_coef_part[n] *= (n==polynomial_degree-1) ? 1 : new_part; }
    term_grad[polynomial_degree-1] = new_grad;

    if (nEdge == polynomial_degree) {
      lambda_for_two_points(*edgeNodePtr(0,0),*vertexNodePtr(nEdge),pt,new_part,new_grad);
      term_grad[polynomial_degree] = new_grad;
    }
    if (nEdge == polynomial_degree - 1) {
      new_part = my_poly_element -> edgePtr(nEdge+1) -> lambda(pt);
      term_grad[polynomial_degree] = my_poly_element -> edgePtr(nEdge+1) -> dLambda();
    }
    result *= new_part;
    for (int n = 0; n < num_term; n++) { term_grad_coef_part[n] *= (n==polynomial_degree) ? 1 : new_part; }
} 
  else { cout << "Index out of scope." << endl; return; }

  for (int i=0; i<num_term; i++) {
    gradresult += term_grad_coef_part[i] * term_grad[i];
  }
}

double DirectSerendipityFE::varphi_edge_nodes_A_Po(int nEdge, int jNode, const Point& pt) const {
  double result;
  Tensor1 gradresult;
  varphi_edge_nodes_A_Po(nEdge, jNode, pt, result, gradresult);
  return result;
}

//varphi before normalization for small r and vertex nodes in A_po
void DirectSerendipityFE::varphi_vertex_nodes_A_Po(int n, const Point& pt, double& result, Tensor1& gradresult) const {
  int num_term = polynomial_degree;
  std::vector<double> term_grad_coef_part(num_term,1);
  std::vector<Tensor1> term_grad(num_term, Tensor1(0,0));
  result = 1;
  gradresult = Tensor1(0,0);
  double new_part;

  if (polynomial_degree == 1) {
    int start = polynomial_degree - 2 + num_vertices;
    int stop = polynomial_degree + num_vertices;
    Point p[2];
    int i = 0;
    for (int nVertex = start; nVertex <= stop; nVertex++) {
      if (n == nVertex % num_vertices) { continue; }
      p[i] = *vertexNodePtr(nVertex % num_vertices);
      i++;
    }
    double new_part; Tensor1 new_grad;
    lambda_for_two_points(p[0],p[1],pt,new_part,new_grad);
    result *= new_part;
    term_grad[0] = new_grad;
  }
  else { //Polynomial degree >= 2
    if (n == polynomial_degree - 2) {
    new_part = my_poly_element -> edgePtr(polynomial_degree) -> lambda(pt);
    result *= new_part;
    for (int n = 0; n < num_term; n++) { term_grad_coef_part[n] *= (n==0) ? 1 : new_part; }
    term_grad[0] = my_poly_element -> edgePtr(polynomial_degree) -> dLambda();

    for (int node=0; node < polynomial_degree-1; node++) {
      double new_part; Tensor1 new_grad;
      lambda_for_two_points(*edgeNodePtr(0,0),*edgeNodePtr(polynomial_degree-1,node),pt,new_part,new_grad);
      result *= new_part;
      for (int n = 0; n < num_term; n++) { term_grad_coef_part[n] *= (n==node+1) ? 1 : new_part; }
      term_grad[node+1] = new_grad;
    }
  } else if (n == polynomial_degree - 1) {
      for (int node=0; node < polynomial_degree-1; node++) {
        double new_part; Tensor1 new_grad;
        lambda_for_two_points(*edgeNodePtr(0,0),*edgeNodePtr(polynomial_degree,node),pt,new_part,new_grad);
        result *= new_part;
        for (int n = 0; n < num_term; n++) { term_grad_coef_part[n] *= (n==node) ? 1 : new_part; }
        term_grad[node] = new_grad;
      }
      double new_part; Tensor1 new_grad;
      lambda_for_two_points(*edgeNodePtr(0,0),*vertexNodePtr(polynomial_degree),pt,new_part,new_grad);
      result *= new_part;
      for (int n = 0; n < num_term; n++) { term_grad_coef_part[n] *= (n==polynomial_degree-1) ? 1 : new_part; }
      term_grad[polynomial_degree-1] = new_grad;
  } 
  else if (n == polynomial_degree) {
    for (int node=0; node < polynomial_degree-1; node++) {
      double new_part; Tensor1 new_grad;
      lambda_for_two_points(*edgeNodePtr(0,0),*edgeNodePtr(polynomial_degree,node),pt,new_part,new_grad);
      result *= new_part;
      for (int n = 0; n < num_term; n++) { term_grad_coef_part[n] *= (n==node) ? 1 : new_part; }
      term_grad[node] = new_grad;
    }
    double new_part; Tensor1 new_grad;
    lambda_for_two_points(*edgeNodePtr(0,0),*vertexNodePtr(polynomial_degree-1),pt,new_part,new_grad);
    result *= new_part;
    for (int n = 0; n < num_term; n++) { term_grad_coef_part[n] *= (n==polynomial_degree-1) ? 1 : new_part; }
    term_grad[polynomial_degree-1] = new_grad;
  } else { cout << "Index is out of scope!" << endl; return; }
  }
  
  for (int i=0; i<num_term; i++) {
    gradresult += term_grad_coef_part[i] * term_grad[i];
  }
}

double DirectSerendipityFE::varphi_vertex_nodes_A_Po(int n, const Point& pt) const {
  double result;
  Tensor1 gradresult;
  varphi_vertex_nodes_A_Po(n, pt, result, gradresult);
  return result;
}


// Basis functions for DirectSerendipityFE //////////////////////////////////////////


void DirectSerendipityFE::initBasis(const Point* pt, int num_pts) {
  if(num_pts <= 0) return;
  num_eval_pts = num_pts;


  // Allocate space for the resulting values
  int sizeOfArray = num_pts * num_nodes;

  if(value_n) delete[] value_n;
  value_n = new double[sizeOfArray];
  if (gradvalue_n) delete[] gradvalue_n;
  gradvalue_n = new Tensor1[sizeOfArray];

  if (polynomial_degree < num_vertices - 2) {
    high_order_ds_space->finiteElementPtr(0)->initBasis(pt, num_pts);
    
    // Update phi_{v,i}
    for (int i = 0; i < num_vertices; i++) {
      // Define the array storing coefficients got by lagrange_v
      // Now we directly take higher order vertex basis functions that are linear on each edge

      for (int pt_index = 0; pt_index < num_pts; pt_index++) {
        double phi_pt = high_order_ds_space->finiteElementPtr(0)->vertexBasis(i,pt_index);
        Tensor1 gradresult = high_order_ds_space->finiteElementPtr(0)->gradVertexBasis(i,pt_index);
        value_n[pt_index*num_nodes+i] = phi_pt;
        gradvalue_n[pt_index*num_nodes+i] = gradresult;
      }
    }

    // Update phi_{e,nEdge,jNode}
    for (int nEdge = 0; nEdge < num_vertices; nEdge++) {
      for (int jNode = 0; jNode < (polynomial_degree - 1); jNode++) {
        // Define the array storing coefficients got by lagrange
        std::vector<double> coef_e_vector(num_vertices - 3);
        double* coef_e = coef_e_vector.data();

        for (int pt_index = 0; pt_index < num_pts; pt_index++) {
          double phi_pt = 0;
          Tensor1 gradresult(0,0);
          for (int sNode = 0; sNode < num_vertices - 3; sNode++) {
            if (pt_index == 0) {
              Point sNodePosition(*high_order_ds_space->finiteElementPtr(0)->edgeNodePtr(nEdge,sNode));
              coef_e[sNode] = lagrange(nEdge,jNode,sNodePosition);
            }
            phi_pt += coef_e[sNode] * high_order_ds_space->finiteElementPtr(0)->edgeBasis(nEdge,sNode,pt_index);
            gradresult += coef_e[sNode] * high_order_ds_space->finiteElementPtr(0)->gradEdgeBasis(nEdge,sNode,pt_index);
          }
          value_n[pt_index*num_nodes + num_vertices + nEdge*(polynomial_degree-1) + jNode] = phi_pt;
          gradvalue_n[pt_index*num_nodes + num_vertices + nEdge*(polynomial_degree-1) + jNode] = gradresult;
        }
      }
    }  
  }
  else {
    //Cell Nodal Basis Functions
    for (int k = 0; k < nCellNodes(); k++) {
      double phi_c = 1; //Store phi_{E,k} evaluated at x_{E,k}, will be updated only when pt_index == 0
      int i = mapI(k,deg_cell_poly); int j = mapJ(k,deg_cell_poly);
      for (int pt_index = 0; pt_index < num_pts; pt_index++) {
        //Define and initialize variable for nodal basis function
        double phi_pt = 1;
        //Define and initialize variable for grad
        int num_term = num_vertices+deg_cell_poly;
        std::vector<double> term_grad_phi_pt(num_term,1);
        std::vector<Tensor1> term_grad(num_term);


        for (int edge_index = 0; edge_index < num_vertices; edge_index++) {
          if (pt_index == 0) { phi_c *= lambda(edge_index,*the_cell_nodes[k]); }
          double new_part_phi_pt = lambda(edge_index, pt[pt_index]);
          phi_pt *= new_part_phi_pt;
          //Update terms for grad
          term_grad[edge_index] = dLambda(edge_index);
          for (int m = 0; m < num_term; m++) { term_grad_phi_pt[m] *= (edge_index==m) ? 1 : new_part_phi_pt; }
        }

        for (int l=0; l<i; l++){
          if (pt_index == 0) { phi_c *= (sqrt(3)/2) * (cellNodePtr(i,j)->val(0)-cellNodePtr(l,j)->val(0)) 
                                  - 0.5 * (cellNodePtr(i,j)->val(1)-cellNodePtr(l,j)->val(1)); }
          double new_part_phi_pt = (sqrt(3)/2) * (pt[pt_index](0)-cellNodePtr(l,j)->val(0)) - 0.5 * (pt[pt_index](1)-cellNodePtr(l,j)->val(1));
          phi_pt *= new_part_phi_pt;
          //Update terms for grad
          term_grad[num_vertices+l] = Tensor1(sqrt(3)/2,-0.5);
          for (int m = 0; m < num_term; m++) { term_grad_phi_pt[m] *= (l==(m-num_vertices)) ? 1 : new_part_phi_pt; }
        }

        for (int l=i+1; l<=deg_cell_poly-j; l++){
          if (pt_index == 0) {  phi_c *= (sqrt(3)/2) * (cellNodePtr(l,j)->val(0)-cellNodePtr(i,j)->val(0))
                                  - 0.5 * (cellNodePtr(i,j)->val(1)-cellNodePtr(l,j)->val(1)); }
          double new_part_phi_pt = (sqrt(3)/2) * (cellNodePtr(l,j)->val(0)-pt[pt_index](0)) - 0.5 * (pt[pt_index](1)-cellNodePtr(l,j)->val(1));
          phi_pt *= new_part_phi_pt;
          //Update terms for grad
          term_grad[num_vertices+l-1] = Tensor1(-sqrt(3)/2,-0.5);
           for (int m = 0; m < num_term; m++) { term_grad_phi_pt[m] *= (l==(m-num_vertices+1)) ? 1 : new_part_phi_pt; }
        }

        for (int l=0; l<j; l++){
          //Define lambda for a horizontal line by myself
          if (pt_index == 0) { phi_c *= cellNodePtr(i,j)->val(1) - cellNodePtr(0,l)->val(1); }
          double new_part_phi_pt = pt[pt_index](1) - cellNodePtr(0,l)->val(1);
          phi_pt *= new_part_phi_pt;
          //Update terms for grad
          term_grad[num_vertices+deg_cell_poly-j+l] = Tensor1(0,1);
          for (int m = 0; m < num_term; m++) { term_grad_phi_pt[m] *= (l==(m-num_vertices-deg_cell_poly+j)) ? 1 : new_part_phi_pt; }
        }
        //k-th cell nodal basis function evaluated at pt-index-th point
        value_n[k + num_vertices*polynomial_degree + pt_index*num_nodes] = phi_pt / phi_c;
        //The grad of k-th cell nodal basis function evaluated at pt-index-th point
        Tensor1 gradresult(0,0);
        for (int m = 0; m < num_term; m++) { gradresult += term_grad[m] * term_grad_phi_pt[m]; }
        gradresult /= phi_c;
        gradvalue_n[k + num_vertices*polynomial_degree + pt_index*num_nodes] = gradresult;
      }
    }

    //Edge Nodal Basis Functions
    if (num_vertices == 3) {
        //Store evaluation of any phi_{e,nEdge,jNode} at each cell node
        std::vector<double> phi_e_at_c_vector(nCellNodes(),1);
        double* phi_e_at_c = phi_e_at_c_vector.data();
        for (int nEdge = 0; nEdge < 3; nEdge++) {
          for (int jNode = 0; jNode < polynomial_degree - 1; jNode++) {
            double phi = 1; // Evaluate phi_{e,nEdge,jNode} at x_{e,nEdge,jNode}, will be updated only when pt_index == 0

            //Renew phi_e_at_c for this nEdge and jNode
            for (int n=0; n<nCellNodes(); n++) {
              phi_e_at_c[n] = 1;
              //Evaluating Lagrange Basis function(nEdge, jNode) (of degree r-2) at pt[pt_index]   
              for (int k=0; k<polynomial_degree-1; k++) {
                if (k == jNode) { continue; }
                else {
                  Point p(*edgeNodePtr(nEdge, jNode)), zero_node(*edgeNodePtr(nEdge,k)); 
                  double proj_pt; Tensor1 dproj_pt;
                  projToEdge(nEdge, *cellNodePtr(n), proj_pt, dproj_pt);
                  phi_e_at_c[n] *= (proj_pt-projToEdge(nEdge,zero_node)) / (projToEdge(nEdge,p)-projToEdge(nEdge,zero_node));
                }
              }
              for (int i = 1; i <= 2; i++) {
                phi_e_at_c[n] *= lambda((nEdge+i) % 3, *cellNodePtr(n));
              }
            }

            for (int pt_index = 0; pt_index < num_pts; pt_index++){
              int global_index =  pt_index * num_nodes + num_vertices + nEdge * (polynomial_degree-1) + jNode;
              
              int num_term = polynomial_degree+1;
              std::vector<double> term_grad_coef_part(num_term,1);
              std::vector<Tensor1> term_grad(num_term, Tensor1(0,0));
              double new_part;
              double phi_pt = 1;
              Tensor1 gradresult(0,0);

              //Evaluating Lagrange Basis function(nEdge, jNode) (of degree r-2) at pt[pt_index]   
              for (int k=0; k<polynomial_degree-1; k++) {
                if (k == jNode) { continue; }
                else {
                  Point p(*edgeNodePtr(nEdge, jNode)), zero_node(*edgeNodePtr(nEdge,k)); 
                  double proj_pt; Tensor1 dproj_pt;
                  projToEdge(nEdge, pt[pt_index], proj_pt, dproj_pt);
                  new_part =  (proj_pt-projToEdge(nEdge,zero_node)) / (projToEdge(nEdge,p)-projToEdge(nEdge,zero_node));
                  phi_pt *= new_part;
                  for (int n = 0; n < num_term; n++) { term_grad_coef_part[n] *= (n==k) ? 1 : new_part; }
                  term_grad[k] = dproj_pt / (projToEdge(nEdge,p)-projToEdge(nEdge,zero_node));
                }
              }

              for (int i = 1; i <= 2; i++) {
                new_part = lambda((nEdge+i) % 3, pt[pt_index]);
                phi_pt *= new_part;
                for (int n = 0; n < num_term; n++) { term_grad_coef_part[n] *= (n==polynomial_degree-2+i) ? 1 : new_part; }
                term_grad[polynomial_degree-2+i] = dLambda((nEdge+i) % 3);
              }

              for (int i=0; i<num_term; i++) { gradresult += term_grad_coef_part[i] * term_grad[i]; }

              //Deduct value at interior nodes if needed
              
              for (int k=0; k<nCellNodes(); k++) {
                  phi_pt -= phi_e_at_c[k] * value_n[k + num_vertices*polynomial_degree + pt_index*num_nodes];
                  gradresult -= phi_e_at_c[k] * gradvalue_n[k + num_vertices*polynomial_degree + pt_index*num_nodes];
              }
              

              if (pt_index == 0) { phi *=  lambda((nEdge+1) % 3, *edgeNodePtr(nEdge, jNode)) * lambda((nEdge+2) % 3, *edgeNodePtr(nEdge, jNode)); }
              value_n[global_index] = phi_pt / phi;
              gradresult /= phi;
              gradvalue_n[global_index] = gradresult;
            }
          }
        }

      } 
      else {
        //Store evaluation of any phi_{e,nEdge,jNode} at each cell node
        std::vector<double> phi_e_at_c_vector(nCellNodes());
        double* phi_e_at_c = phi_e_at_c_vector.data();

        for (int nEdge = 0; nEdge < num_vertices; nEdge++) {
          for (int jNode = 0; jNode < polynomial_degree - 1; jNode++) {
            const int matrix_dim = polynomial_degree-1;
            const double* x = get_result_of_coefficients(nEdge, jNode);
            /* Evaluation of phi_{e,nEdge,jNode} at x_{e,nEdge,jNode} and some helper variable,
              will only be updated when pt_index == 0 */
            double phi = 0, phi_alpha_part = 1, phi_beta_part = 0, p = 0;

            //Renew phi_e_at_c for this nEdge and jNode
            for (int k=0; k<nCellNodes(); k++) {
              double phi_cellNode_alpha_part = 1;
              double phi_cellNode_beta_part = 0;

              for (int m=0; m < num_vertices; m++) {
                if (m == nEdge) { continue; } else {
                phi_cellNode_alpha_part *= lambda(m, *cellNodePtr(k));
                }
              }

              double p_cellNode = 0;
              for (int row = 0; row <= polynomial_degree - num_vertices + 1; row++) {
                p_cellNode += *(x+row) * pow(projToEdge(nEdge,*cellNodePtr(k)),row);
              }
              phi_cellNode_alpha_part *= p_cellNode;

              for (int row = polynomial_degree - num_vertices + 2; row < matrix_dim; row++) {
                int index_for_beta = (nEdge+row-polynomial_degree+num_vertices) % num_vertices;
                phi_cellNode_beta_part += *(x+row) * phi_k_l(nEdge,index_for_beta,*cellNodePtr(k));
              }
              phi_e_at_c[k] = phi_cellNode_alpha_part + phi_cellNode_beta_part;
            }

            //Evaluate phi_{e,nEdge,jNode} at each point in the pt array
            for (int pt_index = 0; pt_index < num_pts; pt_index++) {
              int global_index = pt_index * num_nodes + num_vertices + nEdge * (polynomial_degree-1) + jNode;
              double phi_pt = 0, phi_pt_alpha_part = 1, phi_pt_beta_part = 0, p_pt = 0;

              int num_term_alpha_part = num_vertices+1;
              std::vector<double> term_grad_coef_alpha_part(num_term_alpha_part,1);
              std::vector<Tensor1> term_grad_alpha_part(num_term_alpha_part, Tensor1(0,0));
              double new_part;
              Tensor1 gradresult(0,0), dp_pt(0,0);

              //ALPHA PART
              for (int m=0; m < num_vertices; m++) {
                if (m == nEdge) { continue; } else {                  
                  if (pt_index == 0) { phi_alpha_part *= lambda(m, *edgeNodePtr(nEdge,jNode)); }
                  new_part = lambda(m, pt[pt_index]);
                  phi_pt_alpha_part *= new_part;
                  for (int n = 0; n < num_term_alpha_part; n++) { term_grad_coef_alpha_part[n] *= (n==m) ? 1 : new_part; }
                  term_grad_alpha_part[m] = dLambda(m);
                }
              }

            //p, p_pt and dp_pt were defined and initialized to be 0, 0, and (0,0)
            for (int row = 0; row <= polynomial_degree - num_vertices + 1; row++) {              
              if (pt_index == 0) { p += *(x+row) * pow(projToEdge(nEdge,*edgeNodePtr(nEdge,jNode)),row); }
              double proj_pt; Tensor1 dproj_pt;
              projToEdge(nEdge, pt[pt_index], proj_pt, dproj_pt);
              p_pt += *(x+row) * pow(proj_pt,row);
              if (row!=0) { dp_pt += *(x+row) * row * pow(proj_pt,row-1) * dproj_pt; }
            }
            if (pt_index == 0) { phi_alpha_part *= p; }
            phi_pt_alpha_part *= p_pt;
            for (int n = 0; n < num_term_alpha_part; n++) { term_grad_coef_alpha_part[n] *= (n==num_vertices) ? 1 : p_pt; }
            term_grad_alpha_part[num_vertices] = dp_pt;
          
            for (int i=0; i<num_term_alpha_part; i++) { 
              gradresult += term_grad_coef_alpha_part[i] * term_grad_alpha_part[i];                 
            }
      
            //BETA PART
            for (int row = polynomial_degree - num_vertices + 2; row < matrix_dim; row++) {
              int index_for_beta = (nEdge+row-polynomial_degree+num_vertices) % num_vertices;
              if (pt_index == 0) { phi_beta_part += *(x+row) * phi_k_l(nEdge,index_for_beta,*edgeNodePtr(nEdge,jNode)); }
              double phi_k_l_result; Tensor1 dphi_k_l_result;
              phi_k_l(nEdge,index_for_beta,pt[pt_index],phi_k_l_result,dphi_k_l_result);
              phi_pt_beta_part += *(x+row) * phi_k_l_result;
              gradresult += *(x+row) * dphi_k_l_result;  
            }
    
            phi_pt = phi_pt_alpha_part + phi_pt_beta_part;
            if (pt_index == 0) {  phi = phi_alpha_part + phi_beta_part; }

            // Deduct value at interior nodes if needed
            
            for (int k=0; k<nCellNodes(); k++) {
                phi_pt -= phi_e_at_c[k] * value_n[k + num_vertices*polynomial_degree + pt_index*num_nodes];
                gradresult -= phi_e_at_c[k] * gradvalue_n[k + num_vertices*polynomial_degree + pt_index*num_nodes];
            }
            
            value_n[global_index] = phi_pt / phi;

            gradresult /= phi;
            gradvalue_n[global_index] = gradresult;
            }
          }
        }
      }

    //Vertex Nodal Basis Functions
    for (int i = 0; i<num_vertices; i++) {
      double phi = 1;

      //Evaluate phi_{v,i} at each edge node
      std::vector<double> phi_v_at_e_vector(2*(polynomial_degree-1));
      double* phi_v_at_e = phi_v_at_e_vector.data();
      
      for (int k=i; k<i+2; k++){
        for (int l=0; l<polynomial_degree-1; l++){
          int index_v_at_e = (k-i)*(polynomial_degree-1) + l;
          phi_v_at_e[index_v_at_e] = 1;
          for (int j=i+2; j<i+num_vertices; j++){
            phi_v_at_e[index_v_at_e] *= lambda(j % num_vertices, *edgeNodePtr(k % num_vertices, l));
          }
        }
      }
      
      //Evaluate phi_{v,i} at each cell node
      std::vector<double> phi_v_at_c_vector(nCellNodes());
      double* phi_v_at_c = phi_v_at_c_vector.data();
      for (int k=0; k<nCellNodes(); k++) {
        phi_v_at_c[k] = 1;
        for (int j=i+2; j<i+num_vertices; j++){
          phi_v_at_c[k] *= lambda(j % num_vertices, *cellNodePtr(k));
        }
      }

      for (int pt_index = 0; pt_index<num_pts; pt_index++){
        int num_term = num_vertices-2;
        std::vector<double> term_grad_coef_part(num_term,1);
        std::vector<Tensor1> term_grad(num_term, Tensor1(0,0));
        double new_part;
        double phi_pt = 1;
        Tensor1 gradresult(0,0);

        for (int j=i+2; j<i+num_vertices; j++) {
          if (pt_index == 0) { phi *= lambda(j % num_vertices, *vertexNodePtr(i)); }
          new_part = lambda(j % num_vertices, pt[pt_index]);
          phi_pt*= new_part;
          for (int n = 0; n < num_term; n++) { term_grad_coef_part[n] *= (n==j-i-2) ? 1 : new_part; }
          term_grad[j-i-2] = dLambda(j % num_vertices);
        }

        for (int i=0; i<num_term; i++) { gradresult += term_grad_coef_part[i] * term_grad[i]; }

      //Deduct value at edge nodes
      //And then add value that make vertex basis functions linear on each edge
      
      
      for (int k=i; k<i+2; k++){
        for (int l=0; l<polynomial_degree-1; l++){
          int index_v_at_e = (k-i)*(polynomial_degree-1) + l;
          int global_index = pt_index * num_nodes + num_vertices + (k % num_vertices) * (polynomial_degree-1) + l;
          phi_pt -= phi_v_at_e[index_v_at_e] * value_n[global_index];
          gradresult -= phi_v_at_e[index_v_at_e] * gradvalue_n[global_index];
 
          // Add value that make vertex basis functions linear on each edge
          double linear_correction = 0;

          // If you want vertex basis functions linear on each edge,
          // just uncomment the following block
          
          /*
          if (k == i) {
            linear_correction = lambda((i+num_vertices-1)%num_vertices,*edgeNodePtr(k % num_vertices,l))
                                  /lambda((i+num_vertices-1)%num_vertices,*vertexNodePtr(i % num_vertices));
          } else {
            linear_correction = lambda((i+2)%num_vertices,*edgeNodePtr(k % num_vertices,l))
                                  /lambda((i+2)%num_vertices,*vertexNodePtr(i % num_vertices));
          }  
          */
          phi_pt += linear_correction * phi *  value_n[global_index];
          gradresult += linear_correction * phi * gradvalue_n[global_index]; 
        }
      }
      
      //Deduct value at interior nodes
      
      for (int k=0; k<nCellNodes(); k++) {
        phi_pt -=  phi_v_at_c[k] * value_n[k + num_vertices*polynomial_degree + pt_index*num_nodes];
        gradresult -= phi_v_at_c[k] * gradvalue_n[k + num_vertices*polynomial_degree + pt_index*num_nodes];
      }
      

      value_n[i + pt_index * num_nodes] = phi_pt / phi;

      gradresult /= phi;
      gradvalue_n[i + pt_index * num_nodes] = gradresult;
    }
  }
  }
};

void DirectSerendipityFE::initBasisLowOrderQuad(const Point* pt, int num_pts) {
  assert( polynomial_degree <= 2 && num_vertices == 4 );

  if(num_pts <= 0) return;
  num_eval_pts = num_pts;

  // Allocate space for the resulting values
  int sizeOfArray = num_pts * num_nodes;

  if(value_n) delete[] value_n;
  value_n = new double[sizeOfArray];
  if (gradvalue_n) delete[] gradvalue_n;
  gradvalue_n = new Tensor1[sizeOfArray];

  if (polynomial_degree == 2) {
    // Evaluation of phi_{e,iEdge,0} at \x_{e,iEdge,0}
    // Used for normalizing edge nodal basis function
    std::vector<double> phi_e_normalize(4);

    // Evaluation of phi_{v,iEdge} at \x_{v,iEdge}
    // Used for normalizing vertex nodal basis function    
    std::vector<double> phi_v_normalize(4);

    for (int iEdge = 0 ; iEdge < num_vertices; iEdge++) { 
      phi_e_normalize[iEdge] = lambda(iEdge+1,*edgeNodePtr(iEdge,0)) * lambda(iEdge+num_vertices-1, *edgeNodePtr(iEdge,0)) 
                              * r_supp(iEdge, iEdge+2, *edgeNodePtr(iEdge,0));

      phi_v_normalize[iEdge] = lambda(iEdge+2,*vertexNodePtr(iEdge)) * lambda(iEdge+3,*vertexNodePtr(iEdge));

      for (int k = iEdge; k <= iEdge+1; k++) {
          double eval_phi_e_k = lambda(k+1,*vertexNodePtr(iEdge)) * lambda(k+num_vertices-1, *vertexNodePtr(iEdge)) 
                        * r_supp(k,k+2,*vertexNodePtr(iEdge)) 
                        / (lambda(k+1,*edgeNodePtr(k%num_vertices,0)) * lambda(k+num_vertices-1, *edgeNodePtr(k%num_vertices,0)) 
                        * r_supp(k,k+2,*edgeNodePtr(k%num_vertices,0)));
          phi_v_normalize[iEdge] -= lambda(iEdge+2,*edgeNodePtr(k%num_vertices,0)) * lambda(iEdge+3,*edgeNodePtr(k%num_vertices,0)) * eval_phi_e_k;
        }

    }

    for (int pt_index = 0; pt_index < num_pts; pt_index++) { 
      // We first assign all the evaluated edge nodal 
      // basis functions to the corresponding location in array
      for (int iEdge = 0 ; iEdge < num_vertices; iEdge++) { 
        int global_index = pt_index * num_nodes + num_vertices + iEdge;

        double rSupp = 0; Tensor1 dRSupp(0,0);
        r_supp(iEdge, iEdge+2, pt[pt_index], rSupp, dRSupp);

        double value = lambda(iEdge+1,pt[pt_index]) * lambda(iEdge+num_vertices-1, pt[pt_index]) * rSupp;
        Tensor1 gradvalue = dLambda(iEdge+1) * lambda(iEdge+num_vertices-1, pt[pt_index]) * rSupp 
                          + lambda(iEdge+1,pt[pt_index]) * dLambda(iEdge+num_vertices-1) * rSupp 
                          + lambda(iEdge+1,pt[pt_index]) * lambda(iEdge+num_vertices-1, pt[pt_index]) * dRSupp;

        value_n[global_index] = value / phi_e_normalize[iEdge];
        gradvalue_n[global_index] = gradvalue / phi_e_normalize[iEdge];
      }

      // Now we assign all the evaluated vertex nodal 
      // basis functions to the corresponding location in array
      for (int iVertex = 0; iVertex < num_vertices; iVertex++) {
        int global_index = pt_index * num_nodes + iVertex;

        double value = lambda(iVertex+2,pt[pt_index]) * lambda(iVertex+3,pt[pt_index]);
        Tensor1 gradvalue = dLambda(iVertex+2) * lambda(iVertex+3,pt[pt_index]) + lambda(iVertex+2,pt[pt_index]) * dLambda(iVertex+3);

        for (int k = iVertex; k <= iVertex+1; k++) {
          double eval_phi_e_k = value_n[pt_index * num_nodes + num_vertices + k %num_vertices];
          Tensor1 eval_grad_phi_e_k = gradvalue_n[pt_index * num_nodes + num_vertices + k%num_vertices];

          value -= lambda(iVertex+2,*edgeNodePtr(k%num_vertices,0)) * lambda(iVertex+3,*edgeNodePtr(k%num_vertices,0)) * eval_phi_e_k;
          gradvalue -= lambda(iVertex+2,*edgeNodePtr(k%num_vertices,0)) * lambda(iVertex+3,*edgeNodePtr(k%num_vertices,0)) * eval_grad_phi_e_k;
        }

        value_n[global_index] = value / phi_v_normalize[iVertex];
        gradvalue_n[global_index] = gradvalue / phi_v_normalize[iVertex];

      }

    }
  }

}

// Eval for DirectSerendipityFE

void DirectSerendipityFE::eval(const Point* pt, double* result, Tensor1* gradResult, int num_pts,
			       double* vertex_dofs, double* edge_dofs, double* cell_dofs) {
  initBasis(pt,num_pts);

  for(int n=0; n<num_pts; n++) {
    result[n] = 0;
    gradResult[n].set(0,0);
    Point p(pt[n]);

    for(int i=0; i<num_vertices; i++) {
      result[n] += vertex_dofs[i]*basis(i,n);
      gradResult[n] += vertex_dofs[i]*gradVertexBasis(i,n);
    } 

    for(int k=0; k<num_vertices*(polynomial_degree-1); k++) {
      result[n] += edge_dofs[k]*edgeBasis(k,n);
      gradResult[n] += edge_dofs[k]*gradEdgeBasis(k,n);
    }

    for(int k=0; k<nCellNodes(); k++) {
      result[n] += cell_dofs[k]*cellBasis(k,n);
      gradResult[n] += cell_dofs[k]*gradCellBasis(k,n);
    }
  }
}

void DirectSerendipityFE::eval(const Point* pt, double* result, int num_pts, 
				 double* vertex_dofs, double* edge_dofs, double* cell_dofs) {
  initBasis(pt,num_pts);

  for(int n=0; n<num_pts; n++) {
    result[n] = 0;

    for(int i=0; i<num_vertices; i++) {
      result[n] += vertex_dofs[i]*vertexBasis(i,n);
    } 

    for(int k=0; k<num_vertices*(polynomial_degree-1); k++) {
      result[n] += edge_dofs[k]*edgeBasis(k,n);
    }

    for(int k=0; k<nCellNodes(); k++) {
      result[n] += cell_dofs[k]*cellBasis(k,n);
    }
  }
}

// Output for DirectSerendipityFE

void DirectSerendipityFE::write_raw(std::ofstream& fout) const {
  fout << "    DIRECT SERENDIPITY FE\n";
  fout << "    my_ds_space       = " << my_ds_space << "\n";
  fout << "    my_poly_element   = " << my_poly_element << "\n";
  fout << "    num_vertices      = " << num_vertices << "\n";
  fout << "    polynomial_degree = " <<  polynomial_degree << "\n";

  fout << "    vertex nodes:";
  for(unsigned long int i=0; i<the_vertex_nodes.size(); i++) {
    fout << " " << *the_vertex_nodes[i]
	 << " (" << the_vertex_nodes[i]->nodeIndex() << ")";
  }
  fout << "\n";
  fout << "    edge nodes:";
  for(unsigned long int i=0; i<the_edge_nodes.size(); i++) {
    fout << " " << *the_edge_nodes[i]
	 << " (" << the_edge_nodes[i]->nodeIndex() << ")";
  }
  fout << "\n";
  fout << "    cell nodes:";
  for(unsigned long int i=0; i<the_cell_nodes.size(); i++) {
    fout << " " << *the_cell_nodes[i]
	 << " (" << the_cell_nodes[i]->nodeIndex() << ")";
  }
  fout << "\n";
}

int DirectSerendipityFE::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}
