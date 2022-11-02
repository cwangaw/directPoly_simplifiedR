#include <cmath>
#include <iostream>
#include <vector>

#include <stdio.h>
#include <assert.h>


using namespace std;

#include "Utilities/debug.h"
#include "Mesh/polyMesh.h"
using namespace polymesh;
#include "directSerendipity.h"
using namespace directserendipity;
#include "polyQuadrature.h"
using namespace polyquadrature;
////////////////////////////////////////////////////////////////////////////////
// class Node

NodeType Node::nodeType() const {
  return my_ds_space->the_node_type[my_node_index];
};

BCType Node::bcType() const {
  return my_ds_space->the_bc_type[my_node_index];
};

void Node::write_raw(std::ofstream& fout) const {
  fout << "      DIRECT SERENDIPITY NODE (" << *this << ")\n";
  fout << "      my_ds_space   = " << my_ds_space << "\n";
  fout << "      my_node_index = " << my_node_index << "\n";
  fout << "      my_mesh_index = " << my_mesh_index << "\n";
};

int Node::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
};

////////////////////////////////////////////////////////////////////////////////
// class DirectSerendipityArray

void DirectSerendipityArray::set_directserendipityarray(DirectSerendipity* dsSpace) {
  my_ds_space = dsSpace;
  num_nodes = my_ds_space->nNodes();
  
  if(the_array) delete[] the_array;
  the_array = new double[num_nodes];
}

DirectSerendipityArray::~DirectSerendipityArray() {
  if(the_array) delete[] the_array;
}

void DirectSerendipityArray::eval(const Point* pts, double* result,
				  Tensor1* gradResult, int num_pts) const {
  bool ptEvaluated[num_pts];
  for(int i=0; i<num_pts; i++) ptEvaluated[i] = false;

  // Loop through the elements
  std::vector<Point> elementPts;
  std::vector<int> elementPtsIndex;
  for(int iElement=0; iElement < my_ds_space->my_mesh->nElements(); iElement++) {
    DirectSerendipityFE* finiteElement = &(my_ds_space->the_ds_elements[iElement]);
    PolyElement* element = finiteElement->elementPtr();

    // Set list of points in element
    elementPts.clear();
    elementPtsIndex.clear();
    for(int i=0; i<num_pts; i++) {
      if(ptEvaluated[i]) continue;
      
      if(element->isInElement(pts[i])) {
        elementPts.push_back(pts[i]);
        elementPtsIndex.push_back(i);
        ptEvaluated[i] = true;
      }
    }
    if(elementPts.size() == 0) continue;
    
    // SET DoFs for element
    const int nGon = finiteElement->nVertices();
    double vertex_dofs[nGon];
    double edge_dofs[nGon*(finiteElement->degPolyn() - 1)];
    double cell_dofs[finiteElement->nCellNodes()];
  
    for(int i=0; i<nGon; i++) {
      int iMesh = element->vertexPtr(i)->meshIndex();
      int iVertexDof = my_ds_space->mesh_vertex_to_node_index[iMesh];
      vertex_dofs[i] = the_array[iVertexDof];
    }
    
    for(int i=0; i<nGon; i++) {
      int iMesh = element->edgePtr(i)->meshIndex();
      int iEdgeDof = my_ds_space->mesh_edge_to_first_node_index[iMesh];
      for(int j=0; j<(finiteElement->degPolyn() - 1); j++) {
	      edge_dofs[j + i*(finiteElement->degPolyn() - 1)] = the_array[iEdgeDof+j];
      }
    }

    {
      int iMesh = element->meshIndex();
      int iCellDof = my_ds_space->mesh_element_to_first_node_index[iMesh];
      for(int j=0; j<finiteElement->nCellNodes(); j++) {
	      cell_dofs[j] = the_array[iCellDof+j];
      }
    }

    // Evaluate array at points on element
    double elementResult[elementPts.size()];
    Tensor1* elementGradResult = new Tensor1[elementPts.size()];
    finiteElement->eval(elementPts.data(), elementResult, elementGradResult, elementPts.size(),
			vertex_dofs, edge_dofs, cell_dofs);
    // Place results in global array
    for(unsigned long int i=0; i<elementPts.size(); i++) {
      result[elementPtsIndex[i]] = elementResult[i];
      gradResult[elementPtsIndex[i]] = elementGradResult[i];
    }

    delete[] elementGradResult;
  }

  // Attend to unset points (outside the mesh)
  for(int i=0; i<num_pts; i++) {
    if(!ptEvaluated[i]) {
      result[i] = 0;
      gradResult[i].set(0,0);
    }
  }  
}

void DirectSerendipityArray::eval(const Point& pt, double& result, Tensor1& gradResult) const {
  int iElement = my_ds_space->my_mesh->elementIndex(pt);
  if(iElement < 0) {
    result = 0;
    gradResult.set(0,0);
    return;
  }
  DirectSerendipityFE* elem = &(my_ds_space->the_ds_elements[iElement]);

  // SET DoFs for element
  const int nGon = elem->nVertices();
  double vertex_dofs[nGon];
  double edge_dofs[nGon*(elem->degPolyn() - 1)];
  double cell_dofs[elem->nCellNodes()];

  for(int i=0; i<nGon; i++) {
    int iMesh = elem->elementPtr()->vertexPtr(i)->meshIndex();
    int iVertexDof = my_ds_space->mesh_vertex_to_node_index[iMesh];
    vertex_dofs[i] = the_array[iVertexDof];
  }

  for(int i=0; i<nGon; i++) {
    int iMesh = elem->elementPtr()->edgePtr(i)->meshIndex();
    int iEdgeDof = my_ds_space->mesh_edge_to_first_node_index[iMesh];
    for(int j=0; j<(elem->degPolyn() - 1); j++) {
      edge_dofs[j + i*(elem->degPolyn() - 1)] = the_array[iEdgeDof+j];
    }
  }

  {
    int iMesh = elem->elementPtr()->meshIndex();
    int iCellDof = my_ds_space->mesh_element_to_first_node_index[iMesh];
    for(int j=0; j<elem->nCellNodes(); j++) {
      cell_dofs[j] = the_array[iCellDof+j];
    }
  }

  // Evaluate
  elem->eval(pt, result, gradResult, vertex_dofs, edge_dofs, cell_dofs);
};

double DirectSerendipityArray::eval(const Point& pt) const {
  int iElement = my_ds_space->my_mesh->elementIndex(pt);
  if(iElement < 0) return 0;
  DirectSerendipityFE* elem = &(my_ds_space->the_ds_elements[iElement]);

  // SET DoFs for element
  const int nGon = elem->nVertices();
  double vertex_dofs[nGon];
  double edge_dofs[nGon*(elem->degPolyn() - 1)];
  double cell_dofs[elem->nCellNodes()];
  
  for(int i=0; i<nGon; i++) {
    int iMesh = elem->elementPtr()->vertexPtr(i)->meshIndex();
    int iVertexDof = my_ds_space->mesh_vertex_to_node_index[iMesh];
    vertex_dofs[i] = the_array[iVertexDof];
  }

  for(int i=0; i<nGon; i++) {
    int iMesh = elem->elementPtr()->edgePtr(i)->meshIndex();
    int iEdgeDof = my_ds_space->mesh_edge_to_first_node_index[iMesh];
    for(int j=0; j<(elem->degPolyn() - 1); j++) {
      edge_dofs[j + i*(elem->degPolyn() - 1)] = the_array[iEdgeDof+j];
    }
  }

  {
    int iMesh = elem->elementPtr()->meshIndex();
    int iCellDof = my_ds_space->mesh_element_to_first_node_index[iMesh];
    for(int j=0; j<elem->nCellNodes(); j++) {
      cell_dofs[j] = the_array[iCellDof+j];
    }
  }

  // Evaluate
  return elem->eval(pt, vertex_dofs, edge_dofs, cell_dofs);
};

void DirectSerendipityArray::eval_chunk(const Point* pts, double* result, int num_pts) const {
  bool ptEvaluated[num_pts];
  for(int i=0; i<num_pts; i++) ptEvaluated[i] = false;

  // Loop through the elements
  std::vector<Point> elementPts;
  std::vector<int> elementPtsIndex;
  for(int iElement=0; iElement < my_ds_space->my_mesh->nElements(); iElement++) {
    DirectSerendipityFE* finiteElement = &(my_ds_space->the_ds_elements[iElement]);
    PolyElement* element = finiteElement->elementPtr();

    // Set list of points in element
    elementPts.clear();
    elementPtsIndex.clear();
    for(int i=0; i<num_pts; i++) {
      if(ptEvaluated[i]) continue;
      
      if(element->isInElement(pts[i])) {
        elementPts.push_back(pts[i]);
        elementPtsIndex.push_back(i);
        ptEvaluated[i] = true;
      }
    }
    if(elementPts.size() == 0) continue;
  
    // Place results in global array
    for(unsigned long int i=0; i<elementPts.size(); i++) {
      result[elementPtsIndex[i]] = element->chunkParam();;
    }
  }

  // Attend to unset points (outside the mesh)
  for(int i=0; i<num_pts; i++) {
    if(!ptEvaluated[i]) {
      result[i] = 0;
    }
  }  
}

void DirectSerendipityArray::eval_chunk(const Point& pt, double& result) const {
  int iElement = my_ds_space->my_mesh->elementIndex(pt);
  if(iElement < 0) {
    result = 0;
    return;
  }
  DirectSerendipityFE* elem = &(my_ds_space->the_ds_elements[iElement]);
  // Evaluate
  result = elem->elementPtr()->chunkParam();
};

double DirectSerendipityArray::eval_chunk(const Point& pt) const {
  int iElement = my_ds_space->my_mesh->elementIndex(pt);
  if(iElement < 0) {
    return 0;
  }
  DirectSerendipityFE* elem = &(my_ds_space->the_ds_elements[iElement]);
  // Evaluate
  return elem->elementPtr()->chunkParam();
};

void DirectSerendipityArray::eval_error_on_element(const Point* pts, int num_pts, double* l2Error, double* l2GradError, int refinement_level,
					 double (*referenceFcn)(double,double),
					 Tensor1 (*referenceGradFcn)(double,double)) const {
  bool ptEvaluated[num_pts];
  for(int i=0; i<num_pts; i++) ptEvaluated[i] = false;

  // Loop through the elements
  std::vector<Point> elementPts;
  std::vector<int> elementPtsIndex;
  for(int iElement=0; iElement < my_ds_space->my_mesh->nElements(); iElement++) {
    DirectSerendipityFE* finiteElement = &(my_ds_space->the_ds_elements[iElement]);
    PolyElement* element = finiteElement->elementPtr();

    // Set list of points in element
    elementPts.clear();
    elementPtsIndex.clear();
    for(int i=0; i<num_pts; i++) {
      if(ptEvaluated[i]) continue;
      
      if(element->isInElement(pts[i])) {
        elementPts.push_back(pts[i]);
        elementPtsIndex.push_back(i);
        ptEvaluated[i] = true;
      }
    }
    if(elementPts.size() == 0) continue;


    double l2Error_for_this_element = 0, l2GradError_for_this_element = 0;
    PolyQuadrature quadRule(13,refinement_level); 
    quadRule.setElement(finiteElement->elementPtr());
    
    for(int iPt=0; iPt<quadRule.num(); iPt++) {
      double x = quadRule.pt(iPt).val(0);
      double y = quadRule.pt(iPt).val(1);
      
      double result; Tensor1 gradResult;
      eval(quadRule.pt(iPt), result, gradResult);
      
      double diff = (referenceFcn == nullptr) ? result : (result - referenceFcn(x,y));
      Tensor1 diffGrad = (referenceGradFcn == nullptr) ? gradResult : (gradResult - referenceGradFcn(x,y));
      
      l2Error_for_this_element += pow(diff,2) * quadRule.wt(iPt);
      l2GradError_for_this_element += diffGrad * diffGrad * quadRule.wt(iPt);
    }

    l2Error_for_this_element = sqrt(l2Error_for_this_element);
    l2GradError_for_this_element = sqrt(l2GradError_for_this_element);

    // Place results in global array
    for(unsigned long int i=0; i<elementPts.size(); i++) {
      l2Error[elementPtsIndex[i]] = l2Error_for_this_element;
      l2GradError[elementPtsIndex[i]] =l2GradError_for_this_element;
    }
  }

  // Attend to unset points (outside the mesh)
  for(int i=0; i<num_pts; i++) {
    if(!ptEvaluated[i]) {
      l2Error[i] = 0;
      l2GradError[i] = 0;
    }
  }  
};


void DirectSerendipityArray::l2normError(double& l2Error, double& l2GradError, double& l2Norm, double& l2GradNorm, int refinement_level,
					 double (*referenceFcn)(double,double),
					 Tensor1 (*referenceGradFcn)(double,double)) {
  l2Error = 0, l2GradError = 0, l2Norm = 0, l2GradNorm = 0;
  PolyQuadrature quadRule(13,refinement_level);

  for(int iElement=0; iElement < my_ds_space->mesh()->nElements(); iElement++) {
    DirectSerendipityFE* fePtr = my_ds_space->finiteElementPtr(iElement);
    quadRule.setElement(fePtr->elementPtr());
    
    for(int iPt=0; iPt<quadRule.num(); iPt++) {
      double x = quadRule.pt(iPt).val(0);
      double y = quadRule.pt(iPt).val(1);
      
      double result; Tensor1 gradResult;
      eval(quadRule.pt(iPt), result, gradResult);
      
      double diff = (referenceFcn == nullptr) ? result : (result - referenceFcn(x,y));
      Tensor1 diffGrad = (referenceGradFcn == nullptr) ? gradResult : (gradResult - referenceGradFcn(x,y));
      
      l2Error += pow(diff,2) * quadRule.wt(iPt);
      l2GradError += diffGrad * diffGrad * quadRule.wt(iPt);

      l2Norm += pow(result,2) * quadRule.wt(iPt);
      l2GradNorm += gradResult * gradResult * quadRule.wt(iPt);
    }
  }
  
  l2Error = sqrt(l2Error);
  l2GradError = sqrt(l2GradError);
  l2Norm = sqrt(l2Norm);
  l2GradNorm = sqrt(l2GradNorm);
}

void DirectSerendipityArray::write_matlab_mesh(std::ofstream* fout, std::ofstream* fout_grad,
					       int num_pts_x, int num_pts_y) const {
  if(num_pts_x <= 1) num_pts_x = 2;
  if(num_pts_y <= 1) num_pts_y = 2;

  // Determine mesh of points
  double xMin = my_ds_space->my_mesh->minX();
  double xMax = my_ds_space->my_mesh->maxX();
  double yMin = my_ds_space->my_mesh->minY();
  double yMax = my_ds_space->my_mesh->maxY();

  double dx = (xMax - xMin)/(num_pts_x-1);
  double dy = (yMax - yMin)/(num_pts_y-1);

  Point* pts = new Point[num_pts_x*num_pts_y];
  
  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      pts[j + num_pts_y*i].set(xMin+i*dx, yMin+j*dy);
    }
  }

  // Evaluate
  double result[num_pts_x*num_pts_y];
  Tensor1* gradResult = new Tensor1[num_pts_x*num_pts_y];

  eval(pts, result, gradResult, num_pts_x*num_pts_y);

  // Write file  
  *fout << "mesh(" << xMin << ":" << dx << ":" << xMax << ","
	<< yMin << ":" << dy << ":" << yMax <<",[ ";
  if(fout_grad) {
    *fout_grad << "quiver(" << xMin << ":" << dx << ":" << xMax << ","
	       << yMin << ":" << dy << ":" << yMax <<",[ ";
  }

  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      *fout << result[i + num_pts_x*j] << " ";
      if(fout_grad) *fout_grad << gradResult[i + num_pts_x*j][0] << " ";
    }
    *fout << "; ";
    if(fout_grad) *fout_grad << "; ";
  }
  *fout << "]);\n";
  *fout << "xlabel('x'); ylabel('y');\n";
 
  if(fout_grad) {
    *fout_grad << "],[ \n";

    for(int i=0; i<num_pts_x; i++) {
      for(int j=0; j<num_pts_y; j++) {
	*fout_grad << gradResult[i + num_pts_x*j][1] << " ";
      }
      *fout_grad << "; ";
    }
    *fout_grad << "]);\n";
    *fout_grad << "xlabel('x'); ylabel('y');\n";
  }

  delete[] gradResult;
  delete[] pts;
};

int DirectSerendipityArray::write_matlab_mesh(std::string& filename, std::string& filename_grad,
					      int num_pts_x, int num_pts_y) const {
  std::ofstream fout(filename+".m");
  if( !fout ) return 1;
  std::ofstream fout_grad(filename_grad+".m");
  if( !fout_grad ) return 1;
  write_matlab_mesh(&fout, &fout_grad, num_pts_x, num_pts_y);
  return 0;
}

int DirectSerendipityArray::write_matlab_mesh(std::string& filename,
					      int num_pts_x, int num_pts_y) const {
  std::ofstream fout(filename+".m");
  if( !fout ) return 1;
  write_matlab_mesh(&fout, nullptr, num_pts_x, num_pts_y);
  return 0;
}

void DirectSerendipityArray::write_matlab_mesh_by_pt(std::ofstream& fout, std::ofstream& fout_grad,
						     int num_pts_x, int num_pts_y) const {
  if(num_pts_x <= 1) num_pts_x = 2;
  if(num_pts_y <= 1) num_pts_y = 2;
  
  double xMin = my_ds_space->my_mesh->minX();
  double xMax = my_ds_space->my_mesh->maxX();
  double yMin = my_ds_space->my_mesh->minY();
  double yMax = my_ds_space->my_mesh->maxY();

  double dx = (xMax - xMin)/(num_pts_x-1);
  double dy = (yMax - yMin)/(num_pts_y-1);

  double result;
  Tensor1 gradResult;
  double tensorY[num_pts_x][num_pts_y];
  
  fout << "mesh(" << xMin << ":" << dx << ":" << xMax << ","
       << yMin << ":" << dy << ":" << yMax <<",[ ";
  fout_grad << "quiver(" << xMin << ":" << dx << ":" << xMax << ","
	    << yMin << ":" << dy << ":" << yMax <<",[ ";

  for(int j=0; j<num_pts_y; j++) {
    for(int i=0; i<num_pts_x; i++) {
      eval(Point(xMin+i*dx,yMin+j*dy), result, gradResult);
      
      fout << result << " ";
      fout_grad << gradResult[0] << " ";
      tensorY[i][j] = gradResult[1];
    }
    fout << "; ";
    fout_grad << "; ";
  }
  fout << "]);\n";
  fout << "xlabel('x'); ylabel('y');\n";

  fout_grad << "],[ \n";

  for(int j=0; j<num_pts_y; j++) {
    for(int i=0; i<num_pts_x; i++) {
      fout_grad << tensorY[i][j] << " ";
    }
    fout_grad << "; ";
  }
  fout_grad << "]);\n";
  fout_grad << "xlabel('x'); ylabel('y');\n";
};

void DirectSerendipityArray::write_matlab_mesh_by_pt(std::ofstream& fout,
						     int num_pts_x, int num_pts_y) const {
  if(num_pts_x <= 1) num_pts_x = 2;
  if(num_pts_y <= 1) num_pts_y = 2;
  
  double xMin = my_ds_space->my_mesh->minX();
  double xMax = my_ds_space->my_mesh->maxX();
  double yMin = my_ds_space->my_mesh->minY();
  double yMax = my_ds_space->my_mesh->maxY();

  double dx = (xMax - xMin)/(num_pts_x-1);
  double dy = (yMax - yMin)/(num_pts_y-1);

  fout << "mesh(" << xMin << ":" << dx << ":" << xMax << ","
       << yMin << ":" << dy << ":" << yMax <<",[ ";

  for(int j=0; j<num_pts_y; j++) {
    for(int i=0; i<num_pts_x; i++) {
      fout << eval(Point(xMin+i*dx,yMin+j*dy)) << " ";
    }
    fout << "; ";
  }
  fout << "]);\n";
  fout << "xlabel('x'); ylabel('y');\n";
};

int DirectSerendipityArray::write_matlab_mesh_by_pt(std::string& filename, std::string& filename_grad,
						    int num_pts_x, int num_pts_y) const {
  std::ofstream fout(filename+".m");
  if( !fout ) return 1;
  std::ofstream fout_grad(filename_grad+".m");
  if( !fout_grad ) return 1;
  write_matlab_mesh_by_pt(fout, fout_grad, num_pts_x, num_pts_y);
  return 0;
}

int DirectSerendipityArray::write_matlab_mesh_by_pt(std::string& filename,
						    int num_pts_x, int num_pts_y) const {
  std::ofstream fout(filename+".m");
  if( !fout ) return 1;
  write_matlab_mesh_by_pt(fout, num_pts_x, num_pts_y);
  return 0;
}

void DirectSerendipityArray::write_raw(std::ofstream& fout) const {
  for(int i=0; i<num_nodes; i++) {
    fout << the_array[i] << "  ";
    if( !(i%10) ) fout << "\n";
  }
};

int DirectSerendipityArray::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}

void DirectSerendipityArray::write_matlab_mesh_error(std::ofstream* fout, int num_pts_x, int num_pts_y, double (*referenceFcn)(double,double)) const {
  if(num_pts_x <= 1) num_pts_x = 2;
  if(num_pts_y <= 1) num_pts_y = 2;

  // Determine mesh of points
  double xMin = my_ds_space->my_mesh->minX();
  double xMax = my_ds_space->my_mesh->maxX();
  double yMin = my_ds_space->my_mesh->minY();
  double yMax = my_ds_space->my_mesh->maxY();

  double dx = (xMax - xMin)/(num_pts_x-1);
  double dy = (yMax - yMin)/(num_pts_y-1);

  Point* pts = new Point[num_pts_x*num_pts_y];
  
  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      pts[j + num_pts_y*i].set(xMin+i*dx, yMin+j*dy);
    }
  }

  // Evaluate
  double result[num_pts_x*num_pts_y];
  Tensor1* gradResult = new Tensor1[num_pts_x*num_pts_y];
  eval(pts, result, gradResult, num_pts_x*num_pts_y);

  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      result[j + num_pts_y*i] -= referenceFcn(xMin+i*dx, yMin+j*dy);
    }
  }


  // Write file  
  *fout << "mesh(" << xMin << ":" << dx << ":" << xMax << ","
	<< yMin << ":" << dy << ":" << yMax <<",[ ";

  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      *fout << fabs(result[i + num_pts_x*j]) << " ";
    }
    *fout << "; ";
  }
  *fout << "]);\n";
  *fout << "xlabel('x'); ylabel('y');\n";

  delete[] gradResult;
  delete[] pts;  
}

int DirectSerendipityArray::write_matlab_mesh_error(std::string& filename, int num_pts_x, int num_pts_y, double (*referenceFcn)(double,double)) const {
  std::ofstream fout(filename+".m");
  if( !fout ) return 1;
  write_matlab_mesh_error(&fout, num_pts_x, num_pts_y, referenceFcn);
  return 0;
}


void DirectSerendipityArray::write_matlab_mesh_grad_error(std::ofstream* fout, int num_pts_x, int num_pts_y, 
                                                Tensor1 (*referenceFcn)(double,double)) const {
  if(num_pts_x <= 1) num_pts_x = 2;
  if(num_pts_y <= 1) num_pts_y = 2;

  // Determine mesh of points
  double xMin = my_ds_space->my_mesh->minX();
  double xMax = my_ds_space->my_mesh->maxX();
  double yMin = my_ds_space->my_mesh->minY();
  double yMax = my_ds_space->my_mesh->maxY();

  double dx = (xMax - xMin)/(num_pts_x-1);
  double dy = (yMax - yMin)/(num_pts_y-1);

  Point* pts = new Point[num_pts_x*num_pts_y];
  
  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      pts[j + num_pts_y*i].set(xMin+i*dx, yMin+j*dy);
    }
  }

  // Evaluate
  double result[num_pts_x*num_pts_y];
  Tensor1* gradResult = new Tensor1[num_pts_x*num_pts_y];
  eval(pts, result, gradResult, num_pts_x*num_pts_y);

  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      gradResult[j + num_pts_y*i] -= referenceFcn(xMin+i*dx, yMin+j*dy);
    }
  }

  // Write file  
  *fout << "mesh(" << xMin << ":" << dx << ":" << xMax << ","
	<< yMin << ":" << dy << ":" << yMax <<",[ ";

  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      *fout << gradResult[i + num_pts_x*j].norm() << " ";
    }
    *fout << "; ";
  }
  *fout << "]);\n";
  *fout << "xlabel('x'); ylabel('y');\n";


  delete[] gradResult;
  delete[] pts;  
}

int DirectSerendipityArray::write_matlab_mesh_grad_error(std::string& filename, int num_pts_x, int num_pts_y, Tensor1 (*referenceFcn)(double,double)) const {
  std::ofstream fout(filename+".m");
  if( !fout ) return 1;
  write_matlab_mesh_grad_error(&fout, num_pts_x, num_pts_y, referenceFcn);
  return 0;
}




void DirectSerendipityArray::write_matlab_mesh_error_on_element(std::ofstream* fout, int num_pts_x, int num_pts_y, int refinement_level, double (*referenceFcn)(double,double)) const {

   double* err_array = new double[dsSpace()->mesh()->nElements()];

  // Evaluate error on each element
  for(int i=0; i<dsSpace()->mesh()->nElements(); i++) { 
    PolyElement* elem = dsSpace()->mesh()->elementPtr(i);
    double x = (elem->edgePtr(0)->vertexPtr(1)->val(0)+elem->edgePtr(2)->vertexPtr(1)->val(0))/2;
    double y = (elem->edgePtr(0)->vertexPtr(1)->val(1)+elem->edgePtr(2)->vertexPtr(1)->val(1))/2;
    Point pt(x,y);
    double result; double gradResult;
    eval_error_on_element(&pt, 1, &result, &gradResult, refinement_level, referenceFcn, nullptr);
    err_array[i] = result;
  }


  double minVal = err_array[0];
  double maxVal = err_array[0];
  
  for(int i=1; i<dsSpace()->mesh()->nElements(); i++) {
    if(minVal > err_array[i]) { minVal = err_array[i]; continue; }
    if(maxVal < err_array[i]) maxVal = err_array[i];
  }
  
  *fout << "clf;\n";
  *fout << "hold on;\n";
  *fout << "colormap('turbo');\n";
  *fout << "caxis([" << minVal << "  " << maxVal << "]);\n";
 
  for(int iElement=0; iElement<dsSpace()->mesh()->nElements(); iElement++) {
    PolyElement* elem = dsSpace()->mesh()->elementPtr(iElement);
    const int nGon = elem->nVertices();

    // Plot patch outline
    *fout << "patch([";
    for(int j=0; j<nGon-1; j++) {
      *fout << elem->edgePtr(j)->vertexPtr(1)->val(0) <<",";
    }
    *fout << elem->edgePtr(nGon-1)->vertexPtr(1)->val(0) <<"],[";
    for(int j=0; j<nGon-1; j++) {
      *fout << elem->edgePtr(j)->vertexPtr(1)->val(1) <<",";
    }
    *fout << elem->edgePtr(nGon-1)->vertexPtr(1)->val(1) <<"]";

    // Color from the meshArray
    *fout << "," << err_array[iElement] << ");\n";
  }
  

  delete[] err_array;
}

int DirectSerendipityArray::write_matlab_mesh_error_on_element(std::string& filename, int num_pts_x, int num_pts_y, int refinement_level, double (*referenceFcn)(double,double)) const {
  std::ofstream fout(filename+".m");
  if( !fout ) return 1;
  write_matlab_mesh_error_on_element(&fout, num_pts_x, num_pts_y, refinement_level, referenceFcn);
  return 0;
}


void DirectSerendipityArray::write_matlab_mesh_grad_error_on_element(std::ofstream* fout, int num_pts_x, int num_pts_y, int refinement_level, 
                                                Tensor1 (*referenceFcn)(double,double)) const {

   double* err_array = new double[dsSpace()->mesh()->nElements()];

  // Evaluate error on each element
  for(int i=0; i<dsSpace()->mesh()->nElements(); i++) { 
    PolyElement* elem = dsSpace()->mesh()->elementPtr(i);
    double x = (elem->edgePtr(0)->vertexPtr(1)->val(0)+elem->edgePtr(2)->vertexPtr(1)->val(0))/2;
    double y = (elem->edgePtr(0)->vertexPtr(1)->val(1)+elem->edgePtr(2)->vertexPtr(1)->val(1))/2;
    Point pt(x,y);
    double result; double gradResult;
    eval_error_on_element(&pt, 1, &result, &gradResult, refinement_level, nullptr, referenceFcn);
    err_array[i] = gradResult;
  }


  double minVal = err_array[0];
  double maxVal = err_array[0];
  
  for(int i=0; i<dsSpace()->mesh()->nElements(); i++) {
    if(minVal > err_array[i]) { minVal = err_array[i]; continue; }
    if(maxVal < err_array[i]) maxVal = err_array[i];
  }
  
  *fout << "clf;\n";
  *fout << "hold on;\n";
  *fout << "colormap('turbo');\n";
  *fout << "caxis([" << minVal << "  " << maxVal << "]);\n";
 
  for(int iElement=0; iElement<dsSpace()->mesh()->nElements(); iElement++) {
    PolyElement* elem = dsSpace()->mesh()->elementPtr(iElement);
    const int nGon = elem->nVertices();

    // Plot patch outline
    *fout << "patch([";
    for(int j=0; j<nGon-1; j++) {
      *fout << elem->edgePtr(j)->vertexPtr(1)->val(0) <<",";
    }
    *fout << elem->edgePtr(nGon-1)->vertexPtr(1)->val(0) <<"],[";
    for(int j=0; j<nGon-1; j++) {
      *fout << elem->edgePtr(j)->vertexPtr(1)->val(1) <<",";
    }
    *fout << elem->edgePtr(nGon-1)->vertexPtr(1)->val(1) <<"]";

    // Color from the meshArray
    *fout << "," << err_array[iElement] << ");\n";
  }
  

  delete[] err_array;
}

int DirectSerendipityArray::write_matlab_mesh_grad_error_on_element(std::string& filename, int num_pts_x, int num_pts_y, int refinement_level, Tensor1 (*referenceFcn)(double,double)) const {
  std::ofstream fout(filename+".m");
  if( !fout ) return 1;
  write_matlab_mesh_grad_error_on_element(&fout, num_pts_x, num_pts_y, refinement_level, referenceFcn);
  return 0;
}


void DirectSerendipityArray::write_matlab_mesh_one_over_chunk(std::ofstream* fout, int num_pts_x, int num_pts_y) const {
  if(num_pts_x <= 1) num_pts_x = 2;
  if(num_pts_y <= 1) num_pts_y = 2;

  // Determine mesh of points
  double xMin = my_ds_space->my_mesh->minX();
  double xMax = my_ds_space->my_mesh->maxX();
  double yMin = my_ds_space->my_mesh->minY();
  double yMax = my_ds_space->my_mesh->maxY();

  double dx = (xMax - xMin)/(num_pts_x-1);
  double dy = (yMax - yMin)/(num_pts_y-1);

  Point* pts = new Point[num_pts_x*num_pts_y];
  
  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      pts[j + num_pts_y*i].set(xMin+i*dx, yMin+j*dy);
    }
  }

  // Evaluate
  double result[num_pts_x*num_pts_y];
  eval_chunk(pts, result, num_pts_x*num_pts_y);

  // Write file  
  *fout << "mesh(" << xMin << ":" << dx << ":" << xMax << ","
	<< yMin << ":" << dy << ":" << yMax <<",[ ";

  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      *fout << fabs(1/result[i + num_pts_x*j]) << " ";
    }
    *fout << "; ";
  }
  *fout << "]);\n";
  *fout << "xlabel('x'); ylabel('y');\n";

  delete[] pts;  
}

int DirectSerendipityArray::write_matlab_mesh_one_over_chunk(std::string& filename, int num_pts_x, int num_pts_y) const {
  std::ofstream fout(filename+".m");
  if( !fout ) return 1;
  write_matlab_mesh_one_over_chunk(&fout, num_pts_x, num_pts_y);
  return 0;
}


////////////////////////////////////////////////////////////////////////////////
// Class DirectSerendipity

void DirectSerendipity::set_directserendipity(int polyDeg, PolyMesh* mesh) {
  polynomial_degree = polyDeg;
  my_mesh = mesh;

  // ALLOCATE NODES AND ELEMENTS

  // Determine number of nodes
  
  num_nodes = my_mesh->nVertices() + my_mesh->nEdges()*(polynomial_degree-1);
  for(int i=0; i<my_mesh->nElements(); i++) {
    int degCellPoly = polynomial_degree - my_mesh->elementPtr(i)->nVertices();
    if(degCellPoly >= 0) {
      num_nodes += (degCellPoly+2)*(degCellPoly+1)/2;
    }
  }

  // Allocate
  
  if(the_ds_nodes) delete[] the_ds_nodes;
  the_ds_nodes = new Node[num_nodes];
  
  if(the_node_type) delete[] the_node_type;
  the_node_type = new NodeType[num_nodes];
  
  if(the_bc_type) delete[] the_bc_type;
  the_bc_type = new BCType[num_nodes];

  if(the_ds_elements) delete[] the_ds_elements;
  the_ds_elements = new DirectSerendipityFE[my_mesh->nElements()];

  if(mesh_vertex_to_node_index) delete[] mesh_vertex_to_node_index;
  mesh_vertex_to_node_index = new int[my_mesh->nVertices()];
  
  if(mesh_edge_to_first_node_index) delete[] mesh_edge_to_first_node_index;
  mesh_edge_to_first_node_index = new int[my_mesh->nEdges()];
  for(int i=0; i<my_mesh->nEdges(); i++) mesh_edge_to_first_node_index[i] = -1;
  
  if(mesh_element_to_first_node_index) delete[] mesh_element_to_first_node_index;
  mesh_element_to_first_node_index = new int[my_mesh->nElements()];
  for(int i=0; i<my_mesh->nElements(); i++) mesh_element_to_first_node_index[i] = -1;

  // DETERMINE NODES (ORDERING AND TYPES)
  
  // Loop through elements

  int node_num = -1;
  for(int iElement=0; iElement<my_mesh->nElements(); iElement++) {
    int nGon = my_mesh->elementPtr(iElement)->nVertices();

    // Loop through element vertices (order by lowest element number)
    
    for(int iVertex=0; iVertex<nGon; iVertex++) {
      Vertex* v = my_mesh->elementPtr(iElement)->vertexPtr(iVertex);
      std::vector<int> vertexNbrIndices;
      v->nbrElements(vertexNbrIndices);

      bool newVertex = true;
      for(unsigned long int j=0; j<vertexNbrIndices.size(); j++) {
	if(vertexNbrIndices[j] < iElement) { newVertex = false; break; }
      }

      // Add vertex nodes (if new)

      if(newVertex) {
	node_num++;
	the_ds_nodes[node_num].set(*v, v->meshIndex(), node_num, this);
	the_node_type[node_num] = NodeType::vertex;
	the_bc_type[node_num] = (v->isOnBoundary()) ? BCType::dirichlet : BCType::interior;
	
	mesh_vertex_to_node_index[v->meshIndex()] = node_num;
      }

      // Add edge nodes (if new, order by edge, not oriented edge)

      if(polynomial_degree > 1) {
	Edge* e = my_mesh->edgePtr( my_mesh->elementPtr(iElement)->edgePtr(iVertex)->meshIndex() );
	int edgeNbrIndices[2];
	e->nbrElements(edgeNbrIndices);

	if( (edgeNbrIndices[0] == -1 || edgeNbrIndices[1] == -1
	     || edgeNbrIndices[0] > iElement || edgeNbrIndices[1] > iElement) ) {
	  BCType bc = (e->isOnBoundary()) ? BCType::dirichlet : BCType::interior;
	  Point v0(*(e->vertexPtr(0)));
	  Point v1(*(e->vertexPtr(1)));

	  mesh_edge_to_first_node_index[e->meshIndex()] = node_num + 1;
	  for(int j=0; j<polynomial_degree-1; j++) {
	    double wt = (double)(j+1) / (double)(polynomial_degree);
	    Point node(v1); node -= v0; node *= wt; node += v0; // node = v0 + wt*(v1-v0)
	    
	    node_num++;
	    the_ds_nodes[node_num].set(node, e->meshIndex(), node_num, this);
	    the_node_type[node_num] = NodeType::edge;
	    the_bc_type[node_num] = bc;
	  }
	}
      }
    }

    // Add cell nodes
    
    int degCellPoly = polynomial_degree - nGon;
    if(degCellPoly >= 0) {
      mesh_element_to_first_node_index[iElement] = node_num + 1;
      
      Point center(my_mesh->elementPtr(iElement)->center());
      if(degCellPoly > 0) {// Need interior triangle

	// Find distance to nearest edge
	double dist = my_mesh->elementPtr(iElement)->edgePtr(0)->lambda(center);
	for(int i=1; i<nGon; i++) {
	  double new_dist = my_mesh->elementPtr(iElement)->edgePtr(i)->lambda(center);
	  if(dist > new_dist) dist = new_dist;
	}
	dist *= 0.75;
	
	// Equilateral triangle
	//       (0,sqrt(3)/2)
	//          /   \        Scale by dist
	//        /       \      Center at "center"
	//   (-1/2,0)---(1/2,0)
 	Point basePoint(-0.5*dist,-dist/3);
	basePoint += center;

	Point dx(dist/degCellPoly,0);
	Point dy(0.5*dist/degCellPoly,0.5*std::sqrt(3)*dist/degCellPoly);

	for(int i=0; i<=degCellPoly; i++) {
	  Point nodePoint(basePoint);
	  for(int j=0; j<=degCellPoly-i; j++) {
	    node_num++;
	    the_ds_nodes[node_num].set(nodePoint, iElement, node_num, this);
	    the_node_type[node_num] = NodeType::cell;
	    the_bc_type[node_num] = BCType::interior;
	    
	    nodePoint += dx;
	  }    
	  basePoint += dy;
	}

      } else { // Just use center
	node_num++;
	the_ds_nodes[node_num].set(center, iElement, node_num, this);
	the_node_type[node_num] = NodeType::cell;
	the_bc_type[node_num] = BCType::interior;
      }
    }
    num_nodes = node_num + 1;
  }

  // ALLOCATE FINITE ELEMENTS

  for(int iElement=0; iElement<my_mesh->nElements(); iElement++) {
    std::vector<Node*> elementVertexNodes;
    std::vector<Node*> elementEdgeNodes;
    std::vector<Node*> elementCellNodes;

    elementVertexNodes.clear();
    elementEdgeNodes.clear();
    elementCellNodes.clear();
    
    PolyElement* element = my_mesh->elementPtr(iElement);
    int nGon = element->nVertices();
      
    // Vertices
    for(int i=0; i<nGon; i++) {
      Vertex* v = element->edgePtr(i)->vertexPtr(1);
      int nodeIndex = mesh_vertex_to_node_index[v->meshIndex()];
      Node* node = &the_ds_nodes[nodeIndex];
	
      elementVertexNodes.push_back(node);
    }

    // Edges
    if(polynomial_degree > 1) {
      for(int i=0; i<nGon; i++) {
	OrientedEdge* edge = element->edgePtr(i);
	int nodeIndex = mesh_edge_to_first_node_index[edge->meshIndex()];
	
	for(int j=0; j<polynomial_degree - 1; j++) {
	  Node* node = &the_ds_nodes[nodeIndex+j];
	
	  elementEdgeNodes.push_back(node);
	}
      }
    }

    // Cells
    if(polynomial_degree >= nGon) {
      int nodeIndex = mesh_element_to_first_node_index[iElement];
      
      int degCellPoly = polynomial_degree - nGon;
      for(int j=0; j<(degCellPoly+2)*(degCellPoly+1)/2; j++) {
	Node* node = &the_ds_nodes[nodeIndex+j];
	
	elementCellNodes.push_back(node);
      }
    }
      
    the_ds_elements[iElement].set(this, element, elementVertexNodes,
				  elementEdgeNodes, elementCellNodes);
  }

  // Set local cell node indexing array
  if(map_j_array) delete[] map_j_array;
  int maxDegCellPoly = polynomial_degree - 3;
  if(maxDegCellPoly<0) maxDegCellPoly = -1;
  map_j_array = new std::vector<int>[maxDegCellPoly+1];

  for(int k=0; k<maxDegCellPoly+1; k++) {
    for(int i=0; i<=k; i++) {
      for(int j=k-i; j>=0; j--) {
	map_j_array[k].push_back(i);
      }
    }
  }
};

DirectSerendipity::~DirectSerendipity() {
  if(map_j_array) delete[] map_j_array;
  if(the_ds_elements) delete[] the_ds_elements;

  if(the_bc_type) delete[] the_bc_type;
  if(the_node_type) delete[] the_node_type;
  if(the_ds_nodes) delete[] the_ds_nodes;

  if(mesh_vertex_to_node_index) delete[] mesh_vertex_to_node_index;
  if(mesh_edge_to_first_node_index) delete[] mesh_edge_to_first_node_index;
  if(mesh_element_to_first_node_index) delete[] mesh_element_to_first_node_index;
};

void DirectSerendipity::write_raw(std::ofstream& fout) const {
  fout << "DIRECT SERENDIPITY SPACE\n";
  fout << "polynomial_degree = " <<  polynomial_degree << "\n";
  fout << "my_mesh           = " << my_mesh << "\n";
  fout << "num_nodes         = " << num_nodes << "\n";

  fout << "\ndsSpace nodes:\n";
  for(int i=0; i<num_nodes; i++) {
    fout << "  Node " << i << " (type ";
    switch(the_node_type[i]) {
    case NodeType::vertex: { fout << "vertex"; break; }
    case NodeType::edge:   { fout << "edge"; break; }
    case NodeType::cell:   { fout << "cell"; break; }
    default: { fout << "???"; break; }
    }
    fout << ", bc ";
    switch(the_bc_type[i]) {
    case BCType::interior:  { fout << "interior"; break; }
    case BCType::dirichlet: { fout << "dirichlet"; break; }
    case BCType::neumann:   { fout << "neumann"; break; }
    case BCType::robin:     { fout << "robin"; break; }
    default: { fout << "???"; break; }
   }
    fout << ")\n";
    the_ds_nodes[i].write_raw(fout);
  }

  fout << "\ndsSpace connectivity\n";
  fout << "\nmesh_vertex_to_node_index:\n";
  for(int i=0; i<my_mesh->nVertices(); i++) {
    fout << "  " << mesh_vertex_to_node_index[i];
  }
  fout << "\n";
  fout << "\nmesh_edge_to_first_node_index:\n";
  for(int i=0; i<my_mesh->nEdges(); i++) {
    fout << "  " << mesh_edge_to_first_node_index[i];
  }
  fout << "\n";
  fout << "\nmesh_element_to_first_node_index:\n";
  for(int i=0; i<my_mesh->nElements(); i++) {
    fout << "  " << mesh_element_to_first_node_index[i];
  }
  fout << "\n";

  fout << "\ndsSpace elements:\n";
  for(int i=0; i<my_mesh->nElements(); i++) {
    fout << "  Element " << i << "\n";
    the_ds_elements[i].write_raw(fout);
  }
}

int DirectSerendipity::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}

int DirectSerendipity::write_matlab(std::string& filename) const {
  std::ofstream fout(filename + ".m");
  if( !fout ) return 1;

  // MESH
    
  fout << "clf;\n";
  fout << "hold on;\n";
  for(int i=0; i<my_mesh->nElements(); i++) {
    int nGon = my_mesh->nVerticesOfElement(i);

    fout << "patch([";
    for(int j=0; j<nGon-1; j++) {
      fout << my_mesh->elementPtr(i)->edgePtr(j)->vertexPtr(1)->val(0) <<",";
    }
    fout << my_mesh->elementPtr(i)->edgePtr(nGon-1)->vertexPtr(1)->val(0) <<"],[";
    for(int j=0; j<nGon-1; j++) {
      fout << my_mesh->elementPtr(i)->edgePtr(j)->vertexPtr(1)->val(1) <<",";
    }
    fout << my_mesh->elementPtr(i)->edgePtr(nGon-1)->vertexPtr(1)->val(1) <<"],'w')\n";
  }
  fout << "\n";

  // NODES

  // Vertex nodes
  fout << "scatter([ ";
  for(int i=0; i<num_nodes; i++) {
    if(the_ds_nodes[i].isVertex()) { fout << the_ds_nodes[i][0] << " "; }
  }
  fout << "],[ ";
  for(int i=0; i<num_nodes; i++) {
    if(the_ds_nodes[i].isVertex()) { fout << the_ds_nodes[i][1] << " "; }
  }
  fout << "],'ok','filled')\n";

  // Edge nodes
  fout << "scatter([ ";
  for(int i=0; i<num_nodes; i++) {
    if(the_ds_nodes[i].isEdge()) { fout << the_ds_nodes[i][0] << " "; }
  }
  fout << "],[ ";
  for(int i=0; i<num_nodes; i++) {
    if(the_ds_nodes[i].isEdge()) { fout << the_ds_nodes[i][1] << " "; }
  }
  fout << "],'sb','filled')\n";

  // Cell nodes
  fout << "scatter([ ";
  for(int i=0; i<num_nodes; i++) {
    if(the_ds_nodes[i].isCell()) { fout << the_ds_nodes[i][0] << " "; }
  }
  fout << "],[ ";
  for(int i=0; i<num_nodes; i++) {
    if(the_ds_nodes[i].isCell()) { fout << the_ds_nodes[i][1] << " "; }
  }
  fout << "],'^r','filled')\n";

  // Centers
  fout << "scatter([ ";
  for(int i=0; i<my_mesh->nElements(); i++) {
    fout << my_mesh->elementPtr(i)->center()[0] << " ";
  }
  fout << "],[ ";
  for(int i=0; i<my_mesh->nElements(); i++) {
    fout << my_mesh->elementPtr(i)->center()[1] << " ";
  }
  fout << "],'*g')\n";

  return 0;
}

