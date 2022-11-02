#include <cmath>

#include "polyQuadrature.h"
#include "debug.h"
#include <assert.h>

using namespace polyquadrature;

////////////////////////////////////////////////////////////////////////////////
// class PolyQuadrature


void PolyQuadrature::set_rule(int desired_dop, int refinement_level) {
  my_desired_dop = desired_dop;
  my_refinement_level = (refinement_level < 0) ? 0 : refinement_level;

  int num_of_units = pow(2,my_refinement_level);
  double length_refined_edge = 1/(double)num_of_units;

  my_rule = ruleForTriangle.size()-1;
  for(unsigned long int i=0; i<ruleForTriangle.size(); i++) {
    if(ruleForTriangle[i].dop >= desired_dop) { my_rule = i; break; }
  }
  my_dop  = ruleForTriangle[my_rule].dop;

  // Refined triangle
  num_pts_ref = pow(4,refinement_level)*ruleForTriangle[my_rule].num;

  if(my_pts_ref) delete[] my_pts_ref;
  if(my_wts_ref) delete[] my_wts_ref;
  my_pts_ref = new Point[num_pts_ref];
  my_wts_ref = new double[num_pts_ref];

  // set them
  int loc_pt_index = 0;
  
  // Construction of points in refined reference triangle 


  //       v2           .
  //       | \          .      
  //      v0--v1        .
  for (int i_loc = 0; i_loc < num_of_units; i_loc++) {
    for (int j_loc = 0; j_loc < num_of_units - i_loc; j_loc++) {
      Point refv0_refined(i_loc*length_refined_edge, j_loc*length_refined_edge);
      for (int iPt = 0; iPt < ruleForTriangle[my_rule].num; iPt++) {
        my_pts_ref[loc_pt_index] = refv0_refined + ruleForTriangle[my_rule].pts[iPt] * length_refined_edge;
        my_wts_ref[loc_pt_index] = ruleForTriangle[my_rule].wts[iPt] * length_refined_edge * length_refined_edge / 2;
        loc_pt_index++;
      }
    }
  }
          
  //      v1--v0     .
  //        \ |      .
  //         v2      .
  for (int i_loc = 1; i_loc < num_of_units; i_loc++) {
    for (int j_loc = 1; j_loc <= num_of_units - i_loc; j_loc++) {
      Point refv0_refined(i_loc*length_refined_edge, j_loc*length_refined_edge);
      for (int iPt = 0; iPt < ruleForTriangle[my_rule].num; iPt++) {
        my_pts_ref[loc_pt_index] = refv0_refined - ruleForTriangle[my_rule].pts[iPt] * length_refined_edge;
        my_wts_ref[loc_pt_index] = ruleForTriangle[my_rule].wts[iPt] * length_refined_edge * length_refined_edge / 2;
        loc_pt_index++;
      }
    }
  }
}

void PolyQuadrature::set_element(polymesh::PolyElement* element) {
  my_element = element;
  if(!element) return;
  
  int num_triangles = my_element->nVertices();
  num_pts = num_triangles*num_pts_ref;

  // Quadrature points and weights

  if(my_pts) delete[] my_pts;
  my_pts = new Point[num_pts];

  if(my_wts) delete[] my_wts;
  my_wts = new double[num_pts];

  // Triangles from center of polygon
  Point center = my_element->center();

  for(int i=0; i<num_triangles; i++) {
    // v0 = (0,0)
    Point v1(*(my_element->vertexPtr(i))); v1 -= center;
    Point v2(*(my_element->vertexPtr((i+1) % num_triangles))); v2 -= center;

    Tensor2 mappingMatrix(v1[0], v2[0], v1[1], v2[1]);
    double jacobian = std::abs(mappingMatrix.determinant())/2;

    int kk = i*num_pts_ref;
    for(int j=0; j<num_pts_ref; j++) {
      mappingMatrix.mult(my_pts_ref[j], my_pts[kk+j]);
      my_pts[kk+j] += center;
      my_wts[kk+j] = my_wts_ref[j]*jacobian;
    }
  }
}

PolyQuadrature::~PolyQuadrature() {
  if(my_pts) delete[] my_pts;
  if(my_wts) delete[] my_wts;
  if(my_pts_ref) delete[] my_pts_ref;
  if(my_wts_ref) delete[] my_wts_ref;
}

// Test degree of precision on a mesh over the domain [0,10]^2
void polyquadrature::testPolyQuadrature(polymesh::PolyMesh* mesh, int refinement_level, double eps, int toDOP) {
  auto f = [](Point& x, int i, int j) { return std::pow(x[0],i)*std::pow(x[1],j); };
  auto trueIntegF = [](int i, int j) { return std::pow(10,i+1)/(i+1)*std::pow(10,j+1)/(j+1); };

  PolyQuadrature quadrature;

  std::cout << std::endl;
  for(int testDOP=2; testDOP<=toDOP; testDOP++) {
    quadrature.setRule(testDOP, refinement_level);
    if(testDOP != quadrature.dop() ) continue;
    std::cout << "DOP = " << testDOP << "\n";

    for(int i=0; i<=testDOP; i++) {
      for(int j=0; j<=testDOP-i; j++) {

	double full_integ = 0;
	for(int iElem=0; iElem<mesh->nElements(); iElem++) {
	  quadrature.setElement(mesh->elementPtr(iElem));

	  double integ = 0;
	  for(int k=0; k<quadrature.num(); k++) {
	    integ += f(quadrature.pts()[k],i,j)*quadrature.wt(k);
	  }
	  full_integ += integ;
	}
	double true_integ = trueIntegF(i,j);

	double err = std::abs(true_integ - full_integ);
	if(err > eps) std::cout << "  i = " << i << ", j = " << j << ", err = " << err << "\n";
      }
    }
  }
  std::cout << std::endl;
};

////////////////////////////////////////////////////////////////////////////////
// class PolyEdgeQuadrature

void PolyEdgeQuadrature::set_rule(int desired_dop) {
  my_desired_dop = desired_dop;

    my_rule = ruleForEdge.size()-1;

    for(unsigned long int i=0; i<ruleForEdge.size(); i++) {
      if(ruleForEdge[i].num >= desired_dop) { my_rule = i; break; }
    }

    my_dop  = ruleForEdge[my_rule].num;
    num_pts = ruleForEdge[my_rule].num;
    my_pts_ref  = ruleForEdge[my_rule].pts;
    my_wts_ref  = ruleForEdge[my_rule].wts;
}

void PolyEdgeQuadrature::set_edge(polymesh::Edge* edge) {
  my_edge = edge;
  if(!edge) return;

  // Quadrature points and weights

  if (my_pts) delete[] my_pts;
  my_pts = new Point[num_pts];

  if (my_wts) delete[] my_wts;
  my_wts = new double[num_pts];

  Point v0(*(my_edge->vertexPtr(0)));
  Point v1(*(my_edge->vertexPtr(1)));
  Tensor1 vec(v1-v0);

  for (int i = 0; i < num_pts; i++) {
    my_pts[i].set(v0);
    my_pts[i] += (my_pts_ref[i]+1)/2 * (v1 - v0);
    my_wts[i] = my_wts_ref[i] * vec.norm()/ 2;
  }

}

PolyEdgeQuadrature::~PolyEdgeQuadrature() {
  if (my_pts) delete[] my_pts;
  if (my_wts) delete[] my_wts;
}

