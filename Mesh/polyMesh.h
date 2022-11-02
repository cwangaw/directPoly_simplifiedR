#ifndef __polymesh_h_included__
#define __polymesh_h_included__

#include <string>
#include <fstream>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
// PolyMesh class
//   Includes classes for Edge and PolyElement
//   Vertices are base objects called Point
// Mesh of convex polygons in 2D
//
// Assume base objects: Point and Tensor1 (both two numbers)
////////////////////////////////////////////////////////////////////////////////

#include "baseObjects.h"

namespace polymesh {

  class PolyMesh;

  ////////////////////////////////////////////////////////////////////////////////
  // class Vertex
  ////////////////////////////////////////////////////////////////////////////////

  class Vertex : public Point
  {
  private:
    PolyMesh* my_mesh;
    int my_mesh_index;

    void set_vertex(double p0, double p1, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      the_point[0] = p0; the_point[1] = p1; my_mesh_index = myIndex; my_mesh = myMesh; };

  public:
    Vertex() { set_vertex(0,0); };
    Vertex(double p0, double p1, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      set_vertex(p0, p1, myIndex, myMesh); };
    Vertex(const Point& p, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      set_vertex(p[0], p[1], myIndex, myMesh); };
    Vertex(const Point* p, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      set_vertex((*p)[0], (*p)[1], myIndex, myMesh); };
    
    void set() { set_vertex(0,0); };
    void set(double p0=0, double p1=0, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      set_vertex(p0, p1, myIndex, myMesh); };
    void set(const Point& p, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      set_vertex(p[0], p[1], myIndex, myMesh); };
    void set(const Point* p, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      set_vertex((*p)[0], (*p)[1], myIndex, myMesh); };

    PolyMesh* mesh() const { return my_mesh; };
    int meshIndex() const { return my_mesh_index; };
     
    void nbrElements(std::vector<int>& theNbrIndices) const;
    void nbrEdges(std::vector<int>& theNbrIndices) const;
    bool isOnBoundary() const;

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };

  ////////////////////////////////////////////////////////////////////////////////
  // class Edge
  //      X ------------->------- X
  //     v[0]    |    tangent    v[1]
  //             V
  //           normal
  //    Linear function lambda for the edge
  //    Length of the edge
  //    Knows the two adjoining elements (or sets as outside)
  ////////////////////////////////////////////////////////////////////////////////

  class Edge
  {
  private:
    Vertex* the_vertex[2];
    Tensor1 normal_vector;
    Tensor1 tangent_vector;
    double edge_length;

    PolyMesh* my_mesh;
    int my_mesh_index;

    void set_edge(Vertex* pt0, Vertex* pt1, int myIndex=-1, PolyMesh* myMesh=nullptr);
    
  public:
    Edge(Vertex* pt0, Vertex* pt1, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      set_edge(pt0, pt1, myIndex, myMesh); };
    Edge() : normal_vector(), tangent_vector(), edge_length(0), my_mesh(nullptr), my_mesh_index(-1) {
      the_vertex[0] = nullptr; the_vertex[1] = nullptr; };
    
    void set(Vertex* pt0, Vertex* pt1, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      set_edge(pt0, pt1, myIndex, myMesh); };

    Vertex* vertexPtr(int i) const { return the_vertex[i % 2]; }
    double length() const { return edge_length; };

    PolyMesh* mesh() const { return my_mesh; };
    int meshIndex() const { return my_mesh_index; };

    void nbrElements(int theNbrIndices[2]) const;
    bool isOnBoundary() const;
    
    // Decide if a point is in this edge
    bool isInEdge(const Point& pt) const;

    Tensor1 normal() const { return normal_vector; };
    double normal(int i) const { return normal_vector[i % 2]; };
    Tensor1 tangent() const { return tangent_vector; };
    double tangent(int i) const { return tangent_vector[i % 2]; };
    
    double lambda(double x, double y) const;
    double lambda(const Point& p) const { return lambda(p[0], p[1]); };
    Tensor1 dLambda() const { return -1*normal_vector; };
    double dLambda(int i) const { return -normal_vector[i % 2]; };
    Tensor1 curlLambda() const { return tangent_vector; };
    double curlLambda(int i) const { return tangent_vector[i % 2]; };

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;

    friend class OrientedEdge;
  };

  ////////////////////////////////////////////////////////////////////////////////
  // class OrientedEdge
  //   An edge, or its reverse
  ////////////////////////////////////////////////////////////////////////////////

  class OrientedEdge
  {
  private:
    Edge* the_edge;
    double the_orientation; // 1 or -1
    int iOffset; // 0 or 1
    
    void set_oriented_edge(Edge* e, double orient);

  public:
    OrientedEdge(Edge* e, double orient) { set_oriented_edge(e,orient); };
    OrientedEdge() : the_edge(nullptr), the_orientation(0), iOffset(0) {};
    
    void set(Edge* e, double orient) { set_oriented_edge(e,orient); };

    Edge* edgePtr() const { return the_edge; };
    double orientation() const { return the_orientation; }

    Vertex* vertexPtr(int i) const { return the_edge->the_vertex[(i+iOffset) % 2]; }
    double length() const { return the_edge->edge_length; };

    PolyMesh* mesh() const { return the_edge->my_mesh; };
    int meshIndex() const { return the_edge->my_mesh_index; };
    
    void nbrElements(int theNbrIndices[2]);
    bool isOnBoundary() { return the_edge->isOnBoundary(); }

    Tensor1 normal() const { return the_orientation*the_edge->normal_vector; };
    double normal(int i) const { return the_orientation*the_edge->normal_vector[i % 2]; };
    Tensor1 tangent() const { return the_orientation*the_edge->tangent_vector; };
    double tangent(int i) const { return the_orientation*the_edge->tangent_vector[i % 2]; };
    
    double lambda(double x, double y) const { return the_orientation*the_edge->lambda(x,y); };
    double lambda(const Point& p) const { return the_orientation*the_edge->lambda(p[0], p[1]); };
    Tensor1 dLambda() const { return -1*(the_orientation*the_edge->normal_vector); };
    double dLambda(int i) const { return -the_orientation*the_edge->normal_vector[i % 2]; };
    Tensor1 curLambda() const { return the_orientation*the_edge->tangent_vector; };
    double curlLambda(int i) const { return the_orientation*the_edge->tangent_vector[i % 2]; };
    
    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };

  ////////////////////////////////////////////////////////////////////////////////
  // class PolyElement
  //       1      e1       5
  //        X ------------X      nVertices = number of vertices or edges
  //    e2 /              |      Composed of nVertices edges
  //     /                | e5   Convex, and ordered counterclockwise
  //  2 X                 |
  //   e3 \               |      The edges are stored as pointers to en external
  //       X ------------ X         data structure (the PolyMesh).
  //      3       e4       4
  //   Area of element
  ////////////////////////////////////////////////////////////////////////////////

  class PolyElement
  {
  private:
    int num_vertices;
    OrientedEdge* the_oriented_edge = nullptr;
    bool good_element;
    double the_area;
    Point my_center;  // of largest inscribed circle
    Point my_centroid; // center of mass
    double my_diameter;
    double max_radius;

    PolyMesh* my_mesh;
    int my_mesh_index;
    
    bool isConnected();
    bool isConvexCounterclockwise();
    double computeAreaTriangle(const Point* p0, const Point* p1, const Point* p2);
    double computeArea();
    void computeCenterOfTriangle(const OrientedEdge* e0, const OrientedEdge* e1,
				 const OrientedEdge* e2, Point& center, bool& unique);

    void set_polyelement(int nGon, Edge** theEdges, int myIndex=-1, PolyMesh* myMesh=nullptr);

  public:
    PolyElement(int nGon, Edge** theEdges, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      set_polyelement(nGon, theEdges, myIndex, myMesh); };
    PolyElement() : num_vertices(0), the_oriented_edge(nullptr), good_element(false),
		    the_area(0), my_center(), my_centroid(), my_diameter(0),
		    my_mesh(nullptr), my_mesh_index(-1) {};
    ~PolyElement();

    void set(int nGon, Edge** theEdges, int myIndex=-1, PolyMesh* myMesh=nullptr) {
      set_polyelement(nGon, theEdges, myIndex, myMesh); };
    
    bool isGood() const { return good_element; };
    int nVertices() const { return num_vertices; };

    int iLongestEdge(); // Return the index of longest edge 
    // The input parameter ratio is the criteria of defining small edges
    int countShortEdges(double ratio);
    
    OrientedEdge* edgePtr(int i) { return &the_oriented_edge[i % num_vertices]; };
    Vertex* vertexPtr(int i) { return the_oriented_edge[i % num_vertices].vertexPtr(1); };

    PolyMesh* mesh() const { return my_mesh; };
    int meshIndex() const { return my_mesh_index; };
    
    double area() const { return the_area; };
    Point center()   const { return my_center; };
    Point centroid() const { return my_centroid; };
    double maxRadius() const { return max_radius; }
    double diameter() const { return my_diameter; };
    double chunkParam();

    bool isInElement(const Point& pt) const;
    bool isOnElementBoundary(const Point& pt) const;

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;

    friend class PolyMesh;
  };

  ////////////////////////////////////////////////////////////////////////////////
  // class PolyMesh
  //   Lists of vertices, edges, and elements
  // Construction
  //    On input there are two arrays:
  //      pts: The points or vertices (x1, y1, x2, y2; ... ).
  //           Gives an indexing of the vertices
  //      elementByVerticesIndex: List of elements, by integer numbers for the vertex
  //             index (p1, p2, p3, p7; p1, p4, p5; ...)
  //           Gives an indexing of the elements
  //    Create the intermediate array
  //      edges_by_vertex_and_elements_index: List of edges
  //             (pt1-1, pt1-2, element1-1, element1-2; ...
  //           Gives an indexing of the edges
  //           Gives the neighboring elements (and whether on the boundary).
  ////////////////////////////////////////////////////////////////////////////////

  class PolyMesh
  {
  public:
    
  private:
    typedef int int2[2];

    int num_vertices;
    int num_edges;
    int num_elements;
    int max_edges_of_element;
    bool mesh_is_good;
    double min_x, max_x, min_y, max_y;
    int num_boundary_vertices;
    int num_boundary_edges;
    double max_element_diameter;
    double max_chunkiness_parameter;
    double min_chunkiness_parameter;
    double average_chunkiness_parameter;

    // Polymesh of elements
    Vertex* the_vertices = nullptr;
    Edge* the_edges = nullptr;
    PolyElement* the_elements = nullptr;
    
    // Original mesh
    int* num_vertices_of_element = nullptr;
    int** element_by_vertex_index = nullptr;
    
    // Mesh connectivity
    std::vector<int>* nbr_edges_of_vertex = nullptr;
    std::vector<int>* nbr_elements_of_vertex = nullptr;
    int2* nbr_elements_of_edge = nullptr;

    void set_polymesh(int numVertices, double* pts, int numElements,
	     int* numEdgesOfElement, int** elementByVertexIndex);

  public:
    PolyMesh(int numVertices, double* pts, int numElements,
	     int* numEdgesOfElement, int** elementByVertexIndex) {
      set_polymesh(numVertices, pts, numElements,
		   numEdgesOfElement, elementByVertexIndex); };
    PolyMesh() : num_vertices(0), num_edges(0), num_elements(0),
		 max_edges_of_element(0), mesh_is_good(false),
		 min_x(0), max_x(1), min_y(0), max_y(1),
		 num_boundary_vertices(0), num_boundary_edges(0), max_element_diameter(0),
		 the_vertices(nullptr), the_edges(nullptr), the_elements(nullptr),
		 num_vertices_of_element(nullptr), element_by_vertex_index(nullptr),
		 nbr_edges_of_vertex(nullptr), nbr_elements_of_vertex(nullptr),
		 nbr_elements_of_edge(nullptr) {};
    PolyMesh(PolyElement* single_element); // Create a one element mesh from an element
    ~PolyMesh();

    void set(int numVertices, double* pts, int numElements,
	     int* numEdgesOfElement, int** elementByVertexIndex) {
      set_polymesh(numVertices, pts, numElements,
		   numEdgesOfElement, elementByVertexIndex); };

    int createMesh(char meshType='q', int nx=1, int ny=1, double xMin=0, double xMax=1,
		   double yMin=0, double yMax=1, double distortionFactor=0); // Create simple mesh


    int removeShortEdges(double ratio); // Remove the short edges of the mesh

    // Access functions
    int nVertices() const { return num_vertices; };
    int nEdges() const { return num_edges; };
    int nElements() const { return num_elements; };
    int maxEdgesOfElement() const { return max_edges_of_element; };
    int maxVerticesOfElement() const { return max_edges_of_element; };
    double minX() const { return min_x; }
    double maxX() const { return max_x; }
    double minY() const { return min_y; }
    double maxY() const { return max_y; }
    
    double maxElementDiameter() const { return max_element_diameter; };
    double maxChunkParam() const { return max_chunkiness_parameter; }
    double minChunkParam() const { return min_chunkiness_parameter; }
    double averageChunkParam() const { return average_chunkiness_parameter; }
    
    bool isGood() const { return mesh_is_good; };
    
    PolyElement* elementPtr(int i) { return &the_elements[i % num_elements]; };
    Edge* edgePtr(int i) { return &the_edges[i % num_edges]; };
    Vertex* vertexPtr(int i) { return &the_vertices[i % num_vertices]; };
    
    // Mesh connectivity
    int nVerticesOfElement(int i) const { return num_vertices_of_element[i]; }
     
    void nbrEdgesOfVertex(int i, std::vector<int>& theNbrIndices) const {
      theNbrIndices = nbr_edges_of_vertex[i]; }
    void nbrElementsOfVertex(int i, std::vector<int>& theNbrIndices) const {
      theNbrIndices = nbr_elements_of_vertex[i]; }
    
    void nbrElementsOfEdge(int i, int theNbrIndices[2]) const {
      theNbrIndices[0] = nbr_elements_of_edge[i][0];
      theNbrIndices[1] = nbr_elements_of_edge[i][1]; }

    // The input parameter ratio is the criteria of defining small edges
    int countShortEdges(double ratio);
    
    bool isVertexOnBoundary(int i) const;
    bool isEdgeOnBoundary(int i) const;

    int elementIndex(const Point& pt) const;
    
    // Output
    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;

    int write_matlab(std::string& filename) const;

    friend class PolyElement;
    friend class OrientedEdge;
    friend class Edge;
    friend class Vertex;
  };

};

#endif
