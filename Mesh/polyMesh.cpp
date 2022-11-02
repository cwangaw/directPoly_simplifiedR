#include <cmath>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <assert.h>

#include "debug.h"
#include "polyMesh.h"
using namespace polymesh;


// Function to print the
// index of an element
int getIndex(std::vector<int> v, int K)
{
    auto it = find(v.begin(), v.end(), K);
 
    // If element was found
    if (it != v.end())
    {
        // calculating the index
        // of K
        return it - v.begin();
    }
    else {
        // If the element is not
        // present in the vector
        return -1;
    }
}

////////////////////////////////////////////////////////////////////////////////
// Class Vertex

void Vertex::nbrElements(std::vector<int>& theNbrIndices) const {
  theNbrIndices.clear();
  if(my_mesh) {
    my_mesh->nbrElementsOfVertex(my_mesh_index, theNbrIndices);
  }
}

void Vertex::nbrEdges(std::vector<int>& theNbrIndices) const {
  theNbrIndices.clear();
  if(my_mesh) {
    my_mesh->nbrEdgesOfVertex(my_mesh_index, theNbrIndices);
  }
};

bool Vertex::isOnBoundary() const {
  return (!my_mesh) ? false :
    my_mesh->nbr_elements_of_vertex[my_mesh_index].size()
    != my_mesh->nbr_edges_of_vertex[my_mesh_index].size();
}

void Vertex::write_raw(std::ofstream& fout) const {
  fout << "      VERTEX (" << val(0) << "," <<val(1) << ")\n";
  fout << "      my_mesh       = " << my_mesh << "\n";
  fout << "      my_mesh_index = " << my_mesh_index << "\n";
}

int Vertex::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Class Edge
static const double edge_eps = 1e-8;

void Edge::set_edge(Vertex* pt0, Vertex* pt1, int myIndex, PolyMesh* myMesh) {
  the_vertex[0] = pt0;
  the_vertex[1] = pt1;
  my_mesh = myMesh;
  my_mesh_index = myIndex;

  tangent_vector = *pt1;
  tangent_vector -= *pt0;
  edge_length = tangent_vector.norm();
  if( edge_length != 0 ) tangent_vector /= edge_length;

  normal_vector[0] = tangent_vector[1];
  normal_vector[1] = -tangent_vector[0];
}

double Edge::lambda(double x, double y) const {
  return ((*the_vertex[0])[0] - x)*normal_vector[0] + ((*the_vertex[0])[1] - y)*normal_vector[1];
}

void Edge::nbrElements(int theNbrIndices[2]) const {
  if(!my_mesh) {
    theNbrIndices[0] = -1;
    theNbrIndices[1] = -1;
  } else {
    my_mesh->nbrElementsOfEdge(my_mesh_index, theNbrIndices);
  }
}

bool Edge::isOnBoundary() const {
  return (!my_mesh) ? false :
    my_mesh->nbr_elements_of_edge[my_mesh_index][0] == -1
    || my_mesh->nbr_elements_of_edge[my_mesh_index][1] == -1;
}

bool Edge::isInEdge(const Point& pt) const {
  double distToV0 = Tensor1(*vertexPtr(0)-pt).norm();
  double distToV1 = Tensor1(*vertexPtr(1)-pt).norm();
  if (fabs(distToV0 + distToV1 - edge_length) < edge_eps) {
    return true;
  } else { return false; }
}

void Edge::write_raw(std::ofstream& fout) const {
  fout << "    EDGE\n";
  fout << "    vertex 0:\n"; the_vertex[0]->write_raw(fout);
  fout << "    vertex 1:\n"; the_vertex[1]->write_raw(fout);
  fout << "    normal_vector: (" << normal_vector << ")\n";
  fout << "    tangent_vector: (" << tangent_vector << ")\n";
  fout << "    edge_length   = " << edge_length << "\n";
  fout << "    my_mesh       = " << my_mesh << "\n";
  fout << "    my_mesh_index = " << my_mesh_index << "\n";
}

int Edge::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Class OrientedEdge

void OrientedEdge::set_oriented_edge(Edge* e, double orient) {
  the_edge = e;
  if(orient >= 0) {
    the_orientation = 1;
    iOffset = 0;
  } else {
    the_orientation = -1;
    iOffset = 1;
  }
}

void OrientedEdge::nbrElements(int theNbrIndices[2]) {
  the_edge->nbrElements(theNbrIndices);
  if(the_orientation < 0) {
    int ii = theNbrIndices[0];
    theNbrIndices[0] = theNbrIndices[1];
    theNbrIndices[1] = ii;
  }
}

void OrientedEdge::write_raw(std::ofstream& fout) const {
  fout << "    ORIENTEDEDGE\n";
  fout << "    orientation = " << the_orientation << "\n";
  the_edge->write_raw(fout);
}

int OrientedEdge::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Class PolyElement

static const double polyelement_eps = 1e-8;

void PolyElement::set_polyelement(int nGon, Edge** theEdges, int myIndex, PolyMesh* myMesh) {
  num_vertices = nGon;
  my_mesh = myMesh;
  my_mesh_index = myIndex;

  // Determin orientation
  double orientation[num_vertices];
  for(int i=0; i<num_vertices; i++) {
    if(theEdges[i]->lambda(*(theEdges[(i+1) % num_vertices]->vertexPtr(0))) > 0 ||
       theEdges[i]->lambda(*(theEdges[(i+1) % num_vertices]->vertexPtr(1))) > 0 ) {
      orientation[i] = 1;
    } else {
      orientation[i] = -1;
    }
  }

  // Set up OrientedEdges
  if(the_oriented_edge) delete[] the_oriented_edge;
  the_oriented_edge = new OrientedEdge[num_vertices];
  for(int i=0; i<num_vertices; i++) {
    the_oriented_edge[i].set(theEdges[i],orientation[i]);
  }

  //Find diameter
  my_diameter = 0;
  double h = 0;
  for (int i = 0; i < num_vertices; i++) {
    for (int j = i + 1; j < num_vertices; j++) {
      h = Tensor1(*vertexPtr(i)-*vertexPtr(j)).norm();
      if (my_diameter < h) my_diameter = h;
    }
  }
 
  // Misc
  good_element = (num_vertices >= 3) && isConvexCounterclockwise();

  if(!good_element) {
    the_area = 0;

    // take arithmetic average of vertices
    my_centroid = *(the_oriented_edge[0].vertexPtr(1));
    for(int i=1; i<num_vertices; i++) {
      my_centroid += *(the_oriented_edge[i].vertexPtr(1));
    }
    my_centroid /= num_vertices;
    my_center = my_centroid;
  }

  if(good_element) {
    // Compute centroid and area
    my_centroid.set(0,0);
    the_area = 0;
    for(int i=1; i<num_vertices-1; i++) {
      double area_triangle = computeAreaTriangle(the_oriented_edge[0].vertexPtr(0),
						 the_oriented_edge[i].vertexPtr(0),
						 the_oriented_edge[i].vertexPtr(1));
      the_area += area_triangle;

      Point centroid_triangle(*(the_oriented_edge[0].vertexPtr(0)));
      centroid_triangle += *(the_oriented_edge[i].vertexPtr(0));
      centroid_triangle += *(the_oriented_edge[i].vertexPtr(1));
      centroid_triangle /= 3;

      centroid_triangle *= area_triangle;
      my_centroid += centroid_triangle;
    }
    my_centroid /= the_area;

    // Compute center of largest inscribed circle
    max_radius = 0;
    std::vector<Point> center_candidates;
    center_candidates.clear();

    for(int i=0; i<num_vertices-2; i++) {
      for(int j=i+1; j<num_vertices-1; j++) {
	for(int k=j+1; k<num_vertices; k++) {
	  bool triangle_unique = true;
	  Point triangle_center;

	  computeCenterOfTriangle(&the_oriented_edge[i], &the_oriented_edge[j],
				  &the_oriented_edge[k], triangle_center, triangle_unique);

	  double new_radius = the_oriented_edge[i].lambda(triangle_center);
	  bool inside = true;
	  for(int n=0; n<num_vertices; n++) {
	    if(n==i || n==j || n==k) continue;
	    double dist = the_oriented_edge[n].lambda(triangle_center);
	    if(dist + polyelement_eps < new_radius) { inside = false; break; }
	  }

	  if(inside) {
	    if(max_radius < new_radius) {
	      max_radius = new_radius;
	      center_candidates.clear();
	      center_candidates.push_back(triangle_center);
	    } else if(max_radius == new_radius) {
	      center_candidates.push_back(triangle_center);
	    }
	  }
	}
      }
    }

    if(center_candidates.size() == 0) {
        my_center = my_centroid;
    } else {
      my_center.set(0,0);
      for(unsigned long int i=0; i<center_candidates.size(); i++) {
	my_center += center_candidates[i];
      }
      my_center /= center_candidates.size();	  
    }
  }
};

PolyElement::~PolyElement() {
  if(the_oriented_edge) delete[] the_oriented_edge;
};

bool PolyElement::isConnected() {
  for(int i=0; i<num_vertices; i++) {
    if( the_oriented_edge[i].vertexPtr(1) != the_oriented_edge[(i+1) % num_vertices].vertexPtr(0) )
      return false;
  }
  return true;
};

bool PolyElement::isConvexCounterclockwise() {
  if( !isConnected() ) return false;
  
  for(int i=0; i<num_vertices; i++) {
    double lam = the_oriented_edge[i].lambda(*(the_oriented_edge[(i+1) % num_vertices].vertexPtr(1)));
    if( lam <= polyelement_eps ) return false;
  }
  return true;
};

double PolyElement::computeAreaTriangle(const Point* p0, const Point* p1, const Point* p2) {
  return 0.5*std::abs( ((*p1)[0] - (*p0)[0])*((*p2)[1] - (*p0)[1])
		       - ((*p1)[1] - (*p0)[1])*((*p2)[0] - (*p0)[0]) );
}

double PolyElement::computeArea() {
  double sum_area = 0;
  for(int i=1; i<num_vertices-1; i++) {
    sum_area += computeAreaTriangle(the_oriented_edge[0].vertexPtr(num_vertices-1),
				    the_oriented_edge[i].vertexPtr(0), the_oriented_edge[i].vertexPtr(1));
  }
  return sum_area;
}

void PolyElement::computeCenterOfTriangle(const OrientedEdge* e0, const OrientedEdge* e1,
				   const OrientedEdge* e2, Point& center, bool& unique) {
  unique = true;

  // Reorder so e0 not parallel to e1 and e2
  if(std::abs(e0->tangent()*e1->normal()) <= polyelement_eps) { // e0 || e1
    const OrientedEdge* ee = e2;
    e2 = e0;
    e0 = ee;
    unique = false;
  } else if(std::abs(e0->tangent()*e2->normal()) <= polyelement_eps) { // e0 || e2
    const OrientedEdge* ee = e1;
    e1 = e0;
    e0 = ee;
    unique = false;
  } else if(std::abs(e1->tangent()*e2->normal()) <= polyelement_eps) { // e1 || e2
    unique = false;
  }
  
  // Find ends of triangle base (through e0), v1 and v2
  Point p0(*(e0->vertexPtr(0)));
  Point p1(*(e1->vertexPtr(0)));
  Point p2(*(e2->vertexPtr(0)));

  Tensor1 t0(e0->tangent());
  Tensor1 t1(e1->tangent());
  Tensor1 t2(e2->tangent());

  double s1 = ( (p1[1] - p0[1])*t1[0] - (p1[0] - p0[0])*t1[1] ) / ( t0[1]*t1[0] - t0[0]*t1[1] );
  Point v1(t0); v1 *= s1; v1 += p0; // v1 = p0 + s1*t0;

  double s2 = ( (p2[1] - p0[1])*t2[0] - (p2[0] - p0[0])*t2[1] ) / ( t0[1]*t2[0] - t0[0]*t2[1] );
  Point v2(t0); v2 *= s2; v2 += p0; // v2 = p0 + s2*t0;

  // Find center
  t1 -= t0;
  t2 -= t0;

  double sc = ( (v2[1] - v1[1])*t2[0] - (v2[0] - v1[0])*t2[1] ) / ( t1[1]*t2[0] - t1[0]*t2[1] );
  center = t1; center *= sc; center += v1; // center = v1 + sc*t1;
}

int PolyElement::iLongestEdge() {
  // Find the longest edge
  int iLongestEdge = 0;
  for (int i = 1; i < num_vertices; i++) {
    if (edgePtr(i)->length() > edgePtr(iLongestEdge)->length()) {
      iLongestEdge = i;
    }
  }
  return iLongestEdge;
}


int PolyElement::countShortEdges(double ratio) {
  int nShortEdges = 0;

  // Count short edges
  double shortEdgeLength = edgePtr(iLongestEdge())->length() * ratio;
  for (int i = 0; i < num_vertices; i++) {
    if ( edgePtr(i)->length() <= shortEdgeLength ) {
      nShortEdges += 1;
    }
  }

  return nShortEdges;
}

bool PolyElement::isInElement(const Point& pt) const {
  for(int i=0; i<num_vertices; i++) {
    if(the_oriented_edge[i].lambda(pt) < -polyelement_eps) return false;
  }
  return true;
}

bool PolyElement::isOnElementBoundary(const Point& pt) const {
  for(int i=0; i<num_vertices; i++) {
    double val = the_oriented_edge[i].lambda(pt);
    if(std::abs(val) < polyelement_eps) return true;
    if(val < -polyelement_eps) return false;
  }
  return false;
}

double PolyElement::chunkParam() {
  double rho = diameter();
  for (int i = 0; i < num_vertices; i++) {
    for (int j = i+1; j < num_vertices; j++) {
      for (int k = j+1; k < num_vertices; k++) {
        Edge e0(vertexPtr(i),vertexPtr(j));
        Edge e1(vertexPtr(j),vertexPtr(k));
        Edge e2(vertexPtr(k),vertexPtr(i));

        OrientedEdge oe0(&e0, 1);
        OrientedEdge oe1(&e1, 1);
        OrientedEdge oe2(&e2, 1);

        Point center(0,0);
        bool unique = 0;
        computeCenterOfTriangle(&oe0, &oe1, &oe2, center, unique);
        rho = (oe1.lambda(center) < rho)? oe1.lambda(center) : rho;
      }
    }
  }

  return 4*rho/diameter();
}

void PolyElement::write_raw(std::ofstream& fout) const {
  fout << "  POLYELEMENT\n";
  fout << "  num_vertices  = " << num_vertices << "\n";
  fout << "  good_element  = " << good_element << "\n";
  fout << "  the_area      = " << the_area << "\n";
  fout << "  my_center     = " << my_center << "\n";
  fout << "  my_centroid   = " << my_centroid << "\n";
  fout << "  my_mesh       = " << my_mesh << "\n";
  fout << "  my_mesh_index = " << my_mesh_index << "\n";

  for(int i=0; i<num_vertices; i++) {
    fout << "  the_oriented_edge "<< i << ":\n";
    the_oriented_edge[i].write_raw(fout);
    fout << "\n";
  }
}

int PolyElement::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Class PolyMesh

void PolyMesh::set_polymesh(int numVertices, double* pts, int numElements,
			    int* numEdgesOfElement, int** elementByVertexIndex) {
  int num_elements_old = num_elements;
  num_vertices = numVertices;
  num_elements = numElements;
  mesh_is_good = false;

// Vertices (ordered by input array)

  if(the_vertices) delete[] the_vertices;
  the_vertices = new Vertex[num_vertices];
  if(!the_vertices) return;

  min_x = pts[0];
  max_x = pts[0];
  min_y = pts[1];
  max_y = pts[1];
  
  for(int i=0; i<num_vertices; i++) {
    the_vertices[i].set(pts[2*i],pts[2*i+1],i,this);

    if(pts[2*i] < min_x) {
      min_x = pts[2*i];
    } else if(pts[2*i] > max_x) { max_x = pts[2*i]; }
    if(pts[2*i+1] < min_y) {
      min_y = pts[2*i+1];
    } else if(pts[2*i+1] > max_y) { max_y = pts[2*i+1]; }
  }

  // Number of vertices of element (ordered by input array)
  
  if(num_vertices_of_element) delete[] num_vertices_of_element;
  num_vertices_of_element = new int[num_elements];
  if(!num_vertices_of_element) return;
  int maxNum = 0;
  for(int i=0; i<num_elements; i++) {
    int num = numEdgesOfElement[i];
    if(num<3) return;
    num_vertices_of_element[i] = num;
    if(maxNum < num) maxNum = num;
  }
  max_edges_of_element = maxNum;

  // Save vertex representation
  
  if(element_by_vertex_index) {
    for(int i=0; i<num_elements_old; i++) {
      delete[] element_by_vertex_index[i];
    }
    delete[] element_by_vertex_index;
  }
  element_by_vertex_index = new int*[num_elements];
  if(!element_by_vertex_index) return;
  for(int i=0; i<num_elements; i++) {
    element_by_vertex_index[i] = new int[num_vertices_of_element[i]];
    if(!element_by_vertex_index[i]) return;
  }
  for(int i=0; i<num_elements; i++) {
    for(int j=0; j<num_vertices_of_element[i]; j++) {
      element_by_vertex_index[i][j] = elementByVertexIndex[i][j];
    }
  }

  // Edges (order from low to high vertex index)
  
  if(nbr_edges_of_vertex) delete[] nbr_edges_of_vertex;
  nbr_edges_of_vertex = new std::vector<int>[num_vertices];

  // Determine which edges exist
  bool* pattern = new bool[num_vertices*num_vertices];
  if(!pattern) return;
  for(int i=0; i<num_vertices; i++) {
    for(int j=0; j<num_vertices; j++) pattern[j + num_vertices*i] = false;
  }
  for(int i=0; i<num_elements; i++) {
    for(int j=0; j<num_vertices_of_element[i]; j++) {
      int jPlusOne = (j+1) % num_vertices_of_element[i];
      int index_ij = elementByVertexIndex[i][j];
      int index_ijPlusOne = elementByVertexIndex[i][jPlusOne];
      if(index_ijPlusOne < index_ij) {
	int ihold = index_ij;
	index_ij = index_ijPlusOne;
	index_ijPlusOne = ihold;
      }
      pattern[index_ijPlusOne + num_vertices*index_ij] = true;
    }
  }
  // Edge representation
  num_edges = 0;
  for(int i=0; i<num_vertices; i++) {
    for(int j=i+1; j<num_vertices; j++) if(pattern[j + num_vertices*i]) num_edges++;
  }
  int endVerticesOfEdge[num_edges][2];
  int edgeIndex = 0;
  for(int i=0; i<num_vertices; i++) {
    for(int j=i+1; j<num_vertices; j++) {
      if(pattern[j + num_vertices*i]) {
	endVerticesOfEdge[edgeIndex][0] = i;
	endVerticesOfEdge[edgeIndex][1] = j;
	edgeIndex++;
      }
    }
  }
  delete[] pattern;
  
  // Create Edges (orient from low to high vertex index)
  if(the_edges) delete[] the_edges;
  the_edges = new Edge[num_edges];
  if(!the_edges) return;
  for(int i=0; i<num_edges; i++) {
    the_edges[i].set(&the_vertices[endVerticesOfEdge[i][0]],
		     &the_vertices[endVerticesOfEdge[i][1]],i,this);
    // Mesh connectivity info
    nbr_edges_of_vertex[endVerticesOfEdge[i][0]].push_back(i);
    nbr_edges_of_vertex[endVerticesOfEdge[i][1]].push_back(i);
  }
  
  // Elements (order by input array)
  
  if(nbr_elements_of_vertex) delete[] nbr_elements_of_vertex;
  nbr_elements_of_vertex = new std::vector<int>[num_vertices];
  if(nbr_elements_of_edge) delete[] nbr_elements_of_edge;
  nbr_elements_of_edge = new int2[num_edges];
  for(int i=0; i<num_edges; i++) {
    nbr_elements_of_edge[i][0] = -1;
    nbr_elements_of_edge[i][1] = -1;
  }
  
  if(the_elements) delete[] the_elements;
  the_elements = new PolyElement[num_elements];
  if(!the_elements) return;
  for(int i=0; i<num_elements; i++) {
    int nGon = num_vertices_of_element[i];
    Edge* theEdges[nGon];
    
    // Determine edges
    for(int j=0; j<nGon; j++) {
      int v0 = element_by_vertex_index[i][j];
      int v1 = element_by_vertex_index[i][(j+1) % nGon];
      if( v0 > v1 ) {
	int vv = v0;
	v0 = v1;
	v1 = vv;
      }
      bool noMatch = true;
      for(int k=0; k<num_edges; k++) {
	if( endVerticesOfEdge[k][0] == v0 &&
	    endVerticesOfEdge[k][1] == v1 ) {
	  noMatch = false;
	  theEdges[j] = &the_edges[k];
	  break;
	}
      }
      if(noMatch) return;
    }
    // Set elements
    the_elements[i].set(nGon,theEdges,i,this);
    if(!the_elements[i].isGood()) return;
    // Mesh connectivity info
    for(int j=0; j<nGon; j++) {
      int iEdge = theEdges[j]->meshIndex();
      if(nbr_elements_of_edge[iEdge][0] == -1) {
	nbr_elements_of_edge[iEdge][0] = i;
      } else if(nbr_elements_of_edge[iEdge][1] == -1) {
	nbr_elements_of_edge[iEdge][1] = i;
      } else {
	return;
      }
      int v = (the_elements[i].the_oriented_edge[j].orientation() > 0) ? 1 : 0;
      nbr_elements_of_vertex[endVerticesOfEdge[iEdge][v]].push_back(i);
    }
  }
  
  // Orient nbr elements of edges
  
  for(int j=0; j<num_edges; j++) {
    bool reverse = false;
    int i0 = nbr_elements_of_edge[j][0];
    int i1 = nbr_elements_of_edge[j][1];
    
    if(i0 != -1) {
      if(the_edges[j].lambda(the_elements[i0].center()) < 0) reverse = true;
    } else if(i1 != -1) {
      if(the_edges[j].lambda(the_elements[i1].center()) > 0) reverse = true;
    }
    
    if(reverse) {
      nbr_elements_of_edge[j][0] = i1;
      nbr_elements_of_edge[j][1] = i0;
    }
  }

  //Maximun diameter among elements
  max_element_diameter = 0;
  for (int i = 0; i < nElements(); i++) {
    if (elementPtr(i)->diameter()>max_element_diameter) max_element_diameter = elementPtr(i)->diameter();
  }

  //Maximum chunkiness parameter among elements
  max_chunkiness_parameter = 0;
  for (int i = 0; i < nElements(); i++) {
    if (elementPtr(i)->chunkParam()>max_chunkiness_parameter) max_chunkiness_parameter = elementPtr(i)->chunkParam();
  }

  //Miniimum chunkiness parameter among elements
  min_chunkiness_parameter = max_chunkiness_parameter;
  for (int i = 0; i < nElements(); i++) {
    if (elementPtr(i)->chunkParam()<min_chunkiness_parameter) min_chunkiness_parameter = elementPtr(i)->chunkParam();
  }


  //Average chunkiness parameter of mesh
  average_chunkiness_parameter = 0;
  for (int i = 0; i < nElements(); i++) {
    average_chunkiness_parameter += elementPtr(i)->chunkParam();
  }
  average_chunkiness_parameter /= nElements();


  // Final settings

  num_boundary_vertices = 0;
  for(int i=0; i<num_vertices; i++) {
    if(isVertexOnBoundary(i)) num_boundary_vertices++;
  }
  num_boundary_edges = 0;
  for(int i=0; i<num_edges; i++) {
    if(isEdgeOnBoundary(i)) num_boundary_edges++;
  }
    
  mesh_is_good = true;  
};

// Create a one element mesh from a single element
PolyMesh::PolyMesh(PolyElement* single_element) {
  const int nGon = single_element->nVertices();
  PolyMesh* mesh = single_element->mesh();

  // Extract element defined by vertices array (from the external/global mesh)
  
  int elementByVertexIndex_global[nGon];
  int* elementByVertexIndex_local[1];
  elementByVertexIndex_local[0] = new int[nGon];
  bool vIsSet[nGon];

  for(int i=0; i<nGon; i++) {
    elementByVertexIndex_global[i] = mesh->element_by_vertex_index[single_element->meshIndex()][i];
    elementByVertexIndex_local[0][i] = elementByVertexIndex_global[i];
    vIsSet[i] = false;
  }
  
  // Reduce indices to local sequence
  for(int i=nGon-1; i>=0; i--) {
    int iMax = nGon;
    int maxVal = -1;
    for(int j=0; j<nGon; j++) {
      if(vIsSet[j]) continue;
      
      if(maxVal < elementByVertexIndex_local[0][j]) {
	maxVal = elementByVertexIndex_local[0][j];
	iMax = j;
      }
    }
    elementByVertexIndex_local[0][iMax] = i;
    vIsSet[iMax] = true;
  }

  // Set the points array
  
  double pts[2*nGon];
  for(int i=0; i<nGon; i++) {
    int k_local = elementByVertexIndex_local[0][i];
    int k_global = elementByVertexIndex_global[i];
    pts[2*k_local] = mesh->vertexPtr(k_global)->val(0);
    pts[2*k_local+1] = mesh->vertexPtr(k_global)->val(1);
  }

  // Set local mesh
  
  int numEdgesOfElement[1];
  numEdgesOfElement[0] = nGon;

  set_polymesh(nGon, pts, 1, numEdgesOfElement, elementByVertexIndex_local);

  delete[] elementByVertexIndex_local[0];
}  

PolyMesh::~PolyMesh() {
  if(the_vertices) delete[] the_vertices;
  if(the_edges)    delete[] the_edges;
  if(the_elements) delete[] the_elements;
  
  if(num_vertices_of_element) delete[] num_vertices_of_element;

  if(element_by_vertex_index)  {
    for(int i=0; i<num_elements; i++) { if(element_by_vertex_index[i]) delete[] element_by_vertex_index[i]; }
    delete[] element_by_vertex_index;
  }
  
  if(nbr_edges_of_vertex)    delete[] nbr_edges_of_vertex;
  if(nbr_elements_of_vertex) delete[] nbr_elements_of_vertex;
  if(nbr_elements_of_edge)   delete[] nbr_elements_of_edge;
};

// Create a rectangularr mesh, and possibly distort
// meshTypeC: q=random quadrilateral, d=shape regular deviated, t=random cross-hatched triangular
// Error return 1=bad meshTypeC, 2=bad mesh
int PolyMesh::createMesh(char meshTypeC, int nx, int ny, double xMin, double xMax,
			 double yMin, double yMax, double distortionFactor) {
  int nVertices = (nx+1)*(ny+1);
  int nElements = ( meshTypeC == 't' ) ? 2*nx*ny : nx*ny;
  
  double vertices[2*nVertices];
  int nVerticesPerElement[nElements];
  int* elements[nElements];

  // Vertices

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(-1,1); //doubles from -1 to 1

  double h_x = (xMax-xMin)/nx;
  double h_y = (yMax-yMin)/ny;

  for (int row=0; row<ny+1; row++) {
    for (int col=0; col<nx+1; col++) {
      int iVertex = row*(nx+1)+col;
	
      vertices[2*iVertex]   = xMin + col*h_x; // x
      vertices[2*iVertex+1] = yMin + row*h_y; // y

      // distort (except boundary)
      if(0 < col && col < nx) {
	if(meshTypeC != 'd') {
	  double rand_x = distribution(generator);
	  vertices[2*iVertex]   += distortionFactor*rand_x*h_x;
	} else {
	  int parity = ( 2*(row % 2) - 1 ) * ( 2*(col % 2) - 1 );
	  vertices[2*iVertex]   += parity*distortionFactor*h_x;
	}
      }
      if(meshTypeC != 'd' && 0 < row && row < ny) {
	double rand_y = distribution(generator);
	vertices[2*iVertex+1] += distortionFactor*rand_y*h_y;
      }
    }
  }


  // Elements
  
  switch(meshTypeC) {
  case 'q': //Quadrilateral Mesh construction
  case 'd': { //Deviated Mesh construction
    for (int row=0; row<ny; row++) {
      for (int col=0; col<nx; col++) {
        int iElem = row*nx+col;
        nVerticesPerElement[iElem] = 4;
        elements[iElem] = new int[4];
	
        elements[iElem][0] = row*(nx+1) + col;
        elements[iElem][1] = row*(nx+1) + col + 1;
        elements[iElem][2] = (row+1)*(nx+1) + col + 1;
	elements[iElem][3] = (row+1)*(nx+1) + col;
      }
    }
    break;
  }

  case 't': { //Triangular Mesh construction
    for (int row=0; row<ny; row++) {
      for (int col=0; col<nx; col++) {
        int iRectangle = row*nx+col;
	
        int a = (row+1)*(nx+1) + col;
        int b = row*(nx+1) + col;
        int c = row*(nx+1) + col + 1;
        int d = (row+1)*(nx+1) + col + 1;
	
        //Left triangle
        int iElem = 2*iRectangle;
        nVerticesPerElement[iElem] = 3;
        elements[iElem] = new int[3];
	
        elements[iElem][0] = a;
        elements[iElem][1] = b;
        elements[iElem][2] = d;
	
        //Right triangle
        iElem = 2*iRectangle+1;
        nVerticesPerElement[iElem] = 3;
        elements[iElem] = new int[3];
	
        elements[iElem][0] = d;
        elements[iElem][1] = b;
        elements[iElem][2] = c;
      }
    }
    break;
  }

  default: return 1;
  }

  set(nVertices,vertices,nElements,nVerticesPerElement,elements);
  if(!isGood()) return 2;
  
  for(int i=0; i<nElements; i++) delete[] elements[i];
  return 0;
}



int PolyMesh::removeShortEdges(double ratio) {
  if (!nVertices()) return -1;
  if (ratio <= 0) return 0;

  // We first sort out global indices that are going to be removed
  // And they will first be considered as another index that we are going to keep 
  //(For example, if vertex 0 and vertex 1 form an short edge, 
  // we would first mark vertex 1 as vertex 0, and delete it later)

  std::vector<int> indexToBeRemoved;
  indexToBeRemoved.clear();
  std::vector<int> indexToBeRemoved_mapping;
  indexToBeRemoved_mapping.clear();

  // We loop through all the elements to find small edges
  // and look for indices to be removed, as well as the indices they would be mapped to
  for (int iElement = 0; iElement < nElements(); iElement++) { 
    PolyElement* thisElement = elementPtr(iElement);
    
    // Decide the criteria for short edges
    int iLongEdge = thisElement -> iLongestEdge();
    double shortEdgeLength = thisElement -> edgePtr(iLongEdge) -> length() * ratio;

    // Loop through all the edges of this element
    for (int iEdge = 0; iEdge < thisElement->nVertices(); iEdge++) {
      // If iEdge is an short edge in this element
    
      if ( thisElement -> edgePtr(iEdge) -> length() <= shortEdgeLength )  {
        // We consider the two ends of this edge,
        // If they are not on the list of indexToBeRemoved,
        // we mark both of them as that with smaller global index
        int v0_local = (iEdge+thisElement->nVertices()-1)%thisElement->nVertices();
        int v1_local = iEdge;

        int smaller_global_index = std::min(thisElement -> vertexPtr(v0_local) -> meshIndex(), thisElement -> vertexPtr(v1_local) -> meshIndex());
        int larger_global_index = std::max(thisElement -> vertexPtr(v0_local) -> meshIndex(), thisElement -> vertexPtr(v1_local) -> meshIndex());

        bool larger_edge_on_the_boundary = vertexPtr(larger_global_index)->isOnBoundary();

        // If they are not on the list, we add them to the list
        if (std::count(indexToBeRemoved.begin(), indexToBeRemoved.end(), larger_global_index) == 0 ) {
          if (larger_edge_on_the_boundary) {
            indexToBeRemoved.push_back(smaller_global_index);
            indexToBeRemoved_mapping.push_back(larger_global_index);
          } else {
            indexToBeRemoved.push_back(larger_global_index);
            indexToBeRemoved_mapping.push_back(smaller_global_index);
          }
        }
      }
    }
  }


  // We loop through all the elements to find small edges
  // and look for indices to be removed, as well as the indices they would be mapped to
  for (int iElement = 0; iElement < nElements(); iElement++) { 
    PolyElement* thisElement = elementPtr(iElement);
    
    // Decide the criteria for short edges
    int iLongEdge = thisElement -> iLongestEdge();
    double shortEdgeLength = thisElement -> edgePtr(iLongEdge) -> length() * ratio;

    // Loop through all the edges of this element
    for (int iEdge = 0; iEdge < thisElement->nVertices(); iEdge++) {
      // If iEdge is an short edge in this element

      int v0_local = (iEdge+thisElement->nVertices()-1)%thisElement->nVertices();
      int v1_local = iEdge;

      double v0_x = thisElement -> vertexPtr(v0_local) -> val(0);
      double v0_y = thisElement -> vertexPtr(v0_local) -> val(1);

      double v1_x = thisElement -> vertexPtr(v1_local) -> val(0);
      double v1_y = thisElement -> vertexPtr(v1_local) -> val(1);

      // If v0 is point A and v1 is point B, OR v0 is point B and v1 is point A

      if ( thisElement -> edgePtr(iEdge) -> length() <= shortEdgeLength ) {
 
        // We consider the two ends of this edge,
        // If they are not on the list of indexToBeRemoved,
        // we mark both of them as that with smaller global index

        int smaller_global_index = std::min(thisElement -> vertexPtr(v0_local) -> meshIndex(), thisElement -> vertexPtr(v1_local) -> meshIndex());
        int larger_global_index = std::max(thisElement -> vertexPtr(v0_local) -> meshIndex(), thisElement -> vertexPtr(v1_local) -> meshIndex());

        bool larger_edge_on_the_boundary = vertexPtr(larger_global_index)->isOnBoundary();

        // If they are not on the list, we add them to the list
        if (std::count(indexToBeRemoved.begin(), indexToBeRemoved.end(), larger_global_index) == 0 ) {
          if (larger_edge_on_the_boundary) {
            indexToBeRemoved.push_back(smaller_global_index);
            indexToBeRemoved_mapping.push_back(larger_global_index);
          } else {
            indexToBeRemoved.push_back(larger_global_index);
            indexToBeRemoved_mapping.push_back(smaller_global_index);
          }
        }
      }
    }
  }

  for (unsigned int i = 0; i < indexToBeRemoved.size(); i++) {
    std::cout << "Index to be removed: " << indexToBeRemoved[i] << std::endl;
    std::cout << "Coordinate of this index: " << vertexPtr(indexToBeRemoved[i]) -> val(0) << "," << vertexPtr(indexToBeRemoved[i]) -> val(1) << std::endl;
  } 
  // We loop through all the elements to do the mapping
  // And store them in the double array vector elementByVertexIndexBeforeReordering
  // We set up numEdgesOfElement at the same time by counting #unique indices in each element
  // And we can also get numOfEdgesToBeRemoved
  int numOfEdgesToBeRemoved = 0;
  std::vector<int> numEdgesOfElement(nElements());
  std::vector<int>* elementByVertexIndexBeforeReordering = new std::vector<int>[nElements()];

  for (int iElement = 0; iElement < nElements(); iElement++) {
    PolyElement* thisElement = elementPtr(iElement);

    for (int iVertex = 0; iVertex < thisElement -> nVertices(); iVertex++) {
      int vertex_index = thisElement -> vertexPtr(iVertex) -> meshIndex();
      int locate_vertex_index_on_the_list = getIndex(indexToBeRemoved, vertex_index);
      if ( locate_vertex_index_on_the_list == -1 ) {
        elementByVertexIndexBeforeReordering[iElement].push_back(vertex_index);
      } else {
        elementByVertexIndexBeforeReordering[iElement].push_back(indexToBeRemoved_mapping[locate_vertex_index_on_the_list]);
      }
    }

    // Now we count the new number of vertices;
    int newNumVertices = 0;

    for (int iVertex = 0; iVertex < thisElement->nVertices(); iVertex++) { 
      int repeat = 0;
      for (int i = 0; i < iVertex; i++) {
        if (elementByVertexIndexBeforeReordering[iElement][i] == elementByVertexIndexBeforeReordering[iElement][iVertex]) repeat++;
      }
      if (repeat == 0) newNumVertices++;
    }

    numEdgesOfElement[iElement] = newNumVertices;
    numOfEdgesToBeRemoved += thisElement->nVertices() - newNumVertices;
  }

  numOfEdgesToBeRemoved /= 2; //Each edge was counted twice in two elements sharing it


  // To set up elementByVertexIndex,
  // we first loop through elementByVertexIndexBeforeReordering
  // And construct an array mapping indices on the list to the new indices with order 1,2,3,4,5,...
  // We record their coordinates and find out new numVertices at the same times
  std::vector<int> newIndex(nVertices());
  std::vector<double> pts;
  pts.clear();
  int index = 0;

  for (int i = 0; i < nVertices(); i++) {
     // If this index of vertex is not on elementByVertexIndexBeforeReordering
     int count = 0;
     for (int iElement = 0; iElement < nElements(); iElement++) {
       count += std::count(elementByVertexIndexBeforeReordering[iElement].begin(), 
                           elementByVertexIndexBeforeReordering[iElement].end(), i);
     }
    if (count > 0) {
      newIndex[i] = index;
      // We add its x and y value to pts array
      pts.push_back(vertexPtr(i)->val(0));
      pts.push_back(vertexPtr(i)->val(1));   
      index++;   
    } else { newIndex[i] = -1; }
  }

  int numVertices = index;
  int numElements = nElements();

  int** elementByVertexIndex = new int*[numElements];
  for (int i = 0; i < numElements; i++) {
    elementByVertexIndex[i] = new int[numEdgesOfElement[i]];
  }

  for (int iElement = 0; iElement < numElements; iElement++) {
    int local_index = 0;
    for (int iVertex = 0; iVertex < nVerticesOfElement(iElement); iVertex++) {
      // Test if it repeats
      int repeat = 0;
      for (int i = 0; i < iVertex; i++) {
        if (elementByVertexIndexBeforeReordering[iElement][i] == elementByVertexIndexBeforeReordering[iElement][iVertex]) repeat++;
      }

      // We only log it if it is not a repeated index
      if (repeat == 0) {
        elementByVertexIndex[iElement][local_index] = newIndex[elementByVertexIndexBeforeReordering[iElement][iVertex]];
        local_index++;
      }
    }
    assert(local_index == numEdgesOfElement[iElement]);
  }

  set(numVertices, pts.data(), numElements, numEdgesOfElement.data(), elementByVertexIndex);

  for (int i = 0; i < nElements(); i++) {
    delete[] elementByVertexIndex[i];
  }
  delete[] elementByVertexIndex;

  delete[] elementByVertexIndexBeforeReordering;
  std::cout << "Number of short edges removed: " << numOfEdgesToBeRemoved << std::endl;
  return numOfEdgesToBeRemoved;

}

bool PolyMesh::isVertexOnBoundary(int i) const {
  return nbr_elements_of_vertex[i].size() != nbr_edges_of_vertex[i].size();
}

bool PolyMesh::isEdgeOnBoundary(int i) const {
  return nbr_elements_of_edge[i][0] == -1 || nbr_elements_of_edge[i][1] == -1;
}

int PolyMesh::elementIndex(const Point& pt) const {
  static int previousElementIndex = 0;

  for(int i=previousElementIndex; i<num_elements; i++) {
    if(the_elements[i].isInElement(pt)) {
      previousElementIndex = i;
      return i;
    }
  }
  for(int i=0; i<previousElementIndex; i++) {
    if(the_elements[i].isInElement(pt)) {
      previousElementIndex = i;
      return i;
    }
  }
  return -1;
}

// The input parameter ratio is the criteria of defining small edges
int PolyMesh::countShortEdges(double ratio) {
  int nShortEdges = 0;
  for (int iElement = 0; iElement < nElements(); iElement++) {
    nShortEdges += elementPtr(iElement) -> countShortEdges(ratio);
  }
  return nShortEdges;
}

void PolyMesh::write_raw(std::ofstream& fout) const {
  fout << "POLYMESH\n";
  fout << "num_vertices = " << num_vertices << "\n";
  fout << "num_edges    = " << num_edges << "\n";
  fout << "num_elements = " << num_elements << "\n";
  fout << "mesh_is_good = " <<  mesh_is_good << "\n";
  fout << "num_boundary_vertices = " << num_boundary_vertices << "\n";
  fout << "num_boundary_edges    = " << num_boundary_edges << "\n";
  fout << "max_edges_of_element  = " << max_edges_of_element << "\n";
  fout << "[min_x, max_x]x[min_y, max_y] = [" << min_x << "," << max_x
       << "]x[" << min_y <<"," << max_y << "]\n";
  
  fout << "\npolymesh vertices:\n";
  for(int i=0; i<num_vertices; i++) {
    fout << "  " << the_vertices[i] << "\n";
  }

  for(int i=0; i<num_edges; i++) {
    fout << "\npolymesh edges "<<i<<":\n";
    the_edges[i].write_raw(fout);
    fout << "\n";
  }

  for(int i=0; i<num_elements; i++) {
    fout << "\npolymesh elements "<<i<<":\n";
    the_elements[i].write_raw(fout);
    fout << "\n";
  }

  fout << "\nnum_vertices_of_element: ";
  for(int i=0; i<num_elements; i++) {
    fout << num_vertices_of_element[i] << " ";
  }
  fout << "\n";

  fout << "\nelement_by_vertex_index:\n";
  for(int i=0; i<num_elements; i++) {
    for(int j=0; j<num_vertices_of_element[i]; j++) {
      fout << element_by_vertex_index[i][j] << " ";
    }
    fout << "\n";
  }

  fout << "\nnbr_edges_of_vertex:\n";
  for(int i=0; i<num_vertices; i++) {
    fout << " vertex " << i << ":";
    for(unsigned long int j=0; j<nbr_edges_of_vertex[i].size(); j++) {
      fout << "  " << nbr_edges_of_vertex[i][j];
    }
    fout << "\n";
  }
  
  fout << "\nnbr_elements_of_vertex:\n";
  for(int i=0; i<num_vertices; i++) {
    fout << " vertex " << i << ":";
    for(unsigned long int j=0; j<nbr_elements_of_vertex[i].size(); j++) {
      fout << "  " << nbr_elements_of_vertex[i][j];
    }
    fout << "\n";
  }
  
  fout << "\nnbr_elements_of_edge:\n";
  for(int i=0; i<num_edges; i++) {
    fout << " edge " << i << ": " << nbr_elements_of_edge[i][0]
	 << "  " << nbr_elements_of_edge[i][1] << "\n";
  }
}

int PolyMesh::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}

int PolyMesh::write_matlab(std::string& filename) const {
  std::ofstream fout(filename + ".m");
  if( !fout ) return 1;

  fout << "clf;\n";
  fout << "hold on;\n";
  for(int i=0; i<num_elements; i++) {
    int nGon = num_vertices_of_element[i];

    fout << "patch([";
    for(int j=0; j<nGon-1; j++) {
      fout << the_elements[i].edgePtr(j)->vertexPtr(1)->val(0) <<",";
    }
    fout << the_elements[i].edgePtr(nGon-1)->vertexPtr(1)->val(0) <<"],[";
    for(int j=0; j<nGon-1; j++) {
      fout << the_elements[i].edgePtr(j)->vertexPtr(1)->val(1) <<",";
    }
    fout << the_elements[i].edgePtr(nGon-1)->vertexPtr(1)->val(1) <<"],'w')\n";
  }
  
  return 0;
}
