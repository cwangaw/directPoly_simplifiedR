#ifndef __directmixed_h_included__
#define __directmixed_h_included__

#include "directSerendipity.h"
#include <algorithm>

using namespace std;

namespace directserendipity {

  enum class EdgeBCType { interior, boundary };
  enum class MixedDoFType { vertex, edge, cell, poly };
  enum class MixedDoFBCType { interior, boundary };
  class DirectSerendipity;
  class DirectMixed;
  class DirectMixedHybrid;
  class DirectMixedConf;
  
  ////////////////////////////////////////////////////////////////////////////////
  // class DirectMixedFE
  //    Defined on a poly-element (class PolyElement)
  //    
  //    Gives basis functions that
  //          1. initBasis()
  //              separated into curl(\cDS_{r+1}) part and \x\Po_s(E) part
  //          2. initConfBasis()
  //              has H(div)-conforming properties ordered as in paper
  //      First call initBasis to evaluate all basis functions at a given set of points
  //      Then access the basis functions by function index and pt number
  //
  // Store the result in such a way that the evaluations of all basis functions at
  // one point are put together. The basis functions are stored in the following order:
  //  curl(\cDS_{r+1}) -> \x\Po_{r-1}(E) -> \x\tilde\Po_r(E) ( polynomials of order r only )
  ////////////////////////////////////////////////////////////////////////////////

  class DirectMixedFE
  {
  private:
    DirectMixed* my_dm_space;
    polymesh::PolyElement* my_poly_element;
    bool my_conformity;

    int dim_v;
    int dim_curlpart; // H / C
    int dim_v_div; // H / C
    int dim_supp; // H / C

    int num_vertices; // Redundant with my_poly_element
    int polynomial_degree; // Redundant with my_ds_space
    Point ref_origin;

    // We take curl of r+1 degree basis functions
    polymesh::PolyMesh* one_element_mesh = nullptr;
    DirectSerendipity* high_order_ds_space = nullptr;
    
    // Pointer to Evaluation storage
    int num_eval_pts;
    Tensor1* v_value_n = nullptr;

    double* v_div_value_n = nullptr;
    double* v_edge_value_n = nullptr;

    void set_directmixedfe(DirectMixed* dmSpace, polymesh::PolyElement* element, bool conforming);

  protected:
    DirectMixedFE() : my_dm_space(nullptr), my_poly_element(nullptr), my_conformity(false) {};
    DirectMixedFE( DirectMixed* dmSpace, polymesh::PolyElement* element, bool conforming) { 
      set_directmixedfe(dmSpace, element, conforming); };
    
    ~DirectMixedFE();

    void set(DirectMixed* dmSpace, polymesh::PolyElement* element, bool conforming) {
      set_directmixedfe(dmSpace, element, conforming); };
  
  public:
 
    virtual void initBasis(const Point* pt, int num_pts) { };

    DirectMixed* spacePtr() const { return my_dm_space; };
    polymesh::PolyElement* elementPtr() const { return my_poly_element; };
    bool isConforming() const { return my_conformity; }

    // Access basis functions evaluated at pt[iPt]
    Tensor1 basis(int iFunc, int iPt) const {
      return v_value_n[iPt * dim_v + iFunc];
    };

    double basisDotNu(int iFunc, int nEdge, int iPt) const {
      return v_edge_value_n[iPt * dim_v * num_vertices + dim_v * nEdge + iFunc];
    }

    // If the functions are indexed including curl part, 
    // we need to pass iFunc = index of the function - dimCurlPart
    double basisdivXPo(int iFunc, int iPt) const {
      return v_div_value_n[iPt * dim_v_div + iFunc];
    }

    // Get dimensions of spaces
    int dimVFull() const { return dim_v; };
    int dimVReduced() const { return dim_v - polynomial_degree - 1; };

    // Move to Hybrid class later
    int dimCurlPart() const { return dim_curlpart; }
    int dimCurlPoly() const { return dim_curlpart - dim_supp;}
    int dimCurlSupp() const { return dim_supp; }

    int dimXPoFull() const { return (polynomial_degree+2) * (polynomial_degree+1) /2; };
    int dimXPoReduced() const { return polynomial_degree * (polynomial_degree + 1)/2; };

    // Evaluation \u in full or reduced space or both at a point
    void eval(const Point* pt, Tensor1* result, int num_pts, char type, double* dofs=nullptr);
    void eval(const Point& pt, Tensor1& result, char type, double* dofs=nullptr) {
      eval(&pt, &result, 1, type, dofs); };

    Tensor1 eval(const Point& pt, char type, double* dofs=nullptr) { 
      Tensor1 result;
      eval(&pt, &result, 1, type, dofs); 
      return result; 
    }

    void eval(const Point* pt, Tensor1* fullResult, Tensor1* reducedResult, int num_pts, 
              double* full_dofs=nullptr, double* reduced_dofs=nullptr);
    void eval(const Point& pt, Tensor1& fullResult, Tensor1& reducedResult, 
              double* full_dofs=nullptr, double* reduced_dofs=nullptr) {
      eval(&pt, &fullResult, &reducedResult, 1, full_dofs, reduced_dofs); };

    // Evaluation div \u in full or reduced space or both at a point
    virtual void eval_div(const Point* pt, double* result, int num_pts, char type, double* dofs=nullptr) { };
    void eval_div(const Point& pt, double& result, char type, double* dofs=nullptr) {
      eval_div(&pt, &result, 1, type, dofs); };

    double eval_div(const Point& pt, char type, double* dofs=nullptr) { 
      double result;
      eval_div(&pt, &result, 1, type, dofs); 
      return result; 
    }

    virtual void eval_div(const Point* pt, double* fullResult, double* reducedResult, int num_pts, 
              double* full_dofs=nullptr, double* reduced_dofs=nullptr) { };
    void eval_div(const Point& pt, double& fullResult, double& reducedResult, 
              double* full_dofs=nullptr, double* reduced_dofs=nullptr) {
      eval_div(&pt, &fullResult, &reducedResult, 1, full_dofs, reduced_dofs); };
    
    // Output functions
    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;

    friend class DirectMixedConfFE;
    friend class DirectMixedHybridFE;
  };

  class DirectMixedHybridFE : public DirectMixedFE
  {
    private:
    void set_directmixedhybridfe(DirectMixedHybrid* dmSpace, polymesh::PolyElement* element);
    public:
    DirectMixedHybridFE() { DirectMixedFE(); };
    DirectMixedHybridFE(DirectMixedHybrid* dmSpace, polymesh::PolyElement* element) {
      set_directmixedhybridfe(dmSpace, element);
    }

    void set() { DirectMixedFE(); };
    void set(const DirectMixedFE& fe) { set_directmixedfe(fe.spacePtr(),fe.elementPtr(),false); }
    void set(const DirectMixedHybridFE& fe) { set_directmixedfe(fe.spacePtr(),fe.elementPtr(),false); }
    void set(DirectMixedHybrid* dmSpace, polymesh::PolyElement* element) {
      set_directmixedhybridfe(dmSpace, element);
    }

    void initBasis(const Point* pt, int num_pts);

    void eval_div(const Point* pt, double* result, int num_pts, char type, double* dofs=nullptr);
    void eval_div(const Point* pt, double* fullResult, double* reducedResult, int num_pts, 
              double* full_dofs=nullptr, double* reduced_dofs=nullptr);

  };

  class DirectMixedConfFE : public DirectMixedFE
  {
    private:
    // Binomial coefficients
    int* binomCoeff = nullptr;

    void set_directmixedconffe(DirectMixedConf* dmSpace, polymesh::PolyElement* element);

    public:
    DirectMixedConfFE() { DirectMixedFE(); };
    DirectMixedConfFE(DirectMixedConf* dmSpace, polymesh::PolyElement* element) {
      set_directmixedconffe(dmSpace, element);
    };
    ~DirectMixedConfFE();

    void set() { DirectMixedFE(); };
    void set(DirectMixedConf* dmSpace, polymesh::PolyElement* element) {
      set_directmixedconffe(dmSpace, element);
    };

    // Get binomial coefficients from stored array
    int binomialCoefficient(int n, int k) const {
      return binomCoeff[n*(n+1)/2+k];
    };

    // Integration of x^m*y^n from 0 to t, 
    // on a line parametrized by \x = (1-t)*\v0 + t*\v1
    double integralPolyEdge(const Point v0, const Point v1, int m, int n, double t) const;

    void initBasis(const Point* pt, int num_pts);

    // eval_div for conforming element
    void eval_div(const Point* pt, double* result, int num_pts, char type, double* dofs=nullptr);
    void eval_div(const Point* pt, double* fullResult, double* reducedResult, int num_pts, 
              double* full_dofs=nullptr, double* reduced_dofs=nullptr);


    // Get dimension of cell and xPo part
    int dimCellBasis() { return max( 0, ((polynomial_degree - num_vertices + 3)
                            * (polynomial_degree - num_vertices + 2)) / 2); };
    int dimPolyBasis(char type) { return (type == 'f')? (dimXPoFull()-1):(dimXPoReduced()-1); };


    // return true if the dof is a vertex dof
    bool isVertexDoF(int i) {
      if ( (i >= dimCellBasis() + num_vertices * polynomial_degree) 
          && (i < dimCellBasis() + num_vertices * polynomial_degree + num_vertices)) {
            return true;
      } else { return false; }
    };
    // return the index of vertex corresponding to vertex dof
    int mapDoFToVertex(int i) {
      return (isVertexDoF(i))? (i - (dimCellBasis() + num_vertices * polynomial_degree)):-1;
    }

    // return true if the dof is an edge dof
    bool isEdgeDoF(int i) {
      if ( (i >= dimCellBasis()) 
          && (i < dimCellBasis() + num_vertices * polynomial_degree)) {
            return true;
      } else { return false; }
    };

    // return the corresponding EDGE INDEX if the dof is an edge dof
    int mapDoFToEdge(int i) {
      return (isEdgeDoF(i))? ((i - dimCellBasis())/polynomial_degree):-1;
    }

    // Access functions

    // It actually gives us the basis functions with outer normal fluxes on each edge
    // They have nonzero divergence
    Tensor1 vertexBasis(int iFunc, int iPt) {
      iFunc += dimCellBasis() + num_vertices * polynomial_degree;
      return v_value_n[iPt * dim_v + iFunc];
    };

    double vertexBasisDotNu(int iFunc, int nEdge, int iPt) {
      iFunc += dimCellBasis() + num_vertices * polynomial_degree;
      return v_edge_value_n[iPt * dim_v * num_vertices + dim_v * nEdge + iFunc];
    };

    double vertexBasisDiv(int iFunc, int iPt) {
      return v_div_value_n[iPt * dim_v_div + iFunc];
    };

    // It gives us curl of edge basis functions from \cDS^{r+1}
    // They have zero divergence
    Tensor1 edgeBasis(int iEdge, int jNode, int iPt) {
      int iFunc = dimCellBasis() + iEdge * polynomial_degree + jNode;
      return v_value_n[iPt * dim_v + iFunc];
    };

    double edgeBasisDotNu(int iEdge, int jNode, int nEdge, int iPt) {
      int iFunc = dimCellBasis() + iEdge * polynomial_degree + jNode;
      return v_edge_value_n[iPt * dim_v * num_vertices + dim_v * nEdge + iFunc];
    };

    double edgeBasisDiv(int iEdge, int jNode, int iPt) { return 0; };

    // It gives us curl of cell basis functions from \cDS^{r+1}
    // They have zero divergence
    Tensor1 cellBasis(int iFunc, int iPt) {
      return v_value_n[iPt * dim_v + iFunc];
    };

    double cellBasisDotNu(int iFunc, int nEdge, int iPt) {
      return v_edge_value_n[iPt * dim_v * num_vertices + dim_v * nEdge + iFunc];
    };

    double cellBasisDiv(int iFunc, int iPt) { return 0; };

    // It gives us bubble functions constructed with \x * p_i(\x)
    Tensor1 polyBasis(int iFunc, int iPt) {
      iFunc += dimCellBasis() + num_vertices * polynomial_degree + num_vertices;
      return v_value_n[iPt * dim_v + iFunc];
    };

    double polyBasisDotNu(int iFunc, int nEdge, int iPt) {
      iFunc += dimCellBasis() + num_vertices * polynomial_degree + num_vertices;
      return v_edge_value_n[iPt * dim_v * num_vertices + dim_v * nEdge + iFunc];
    };

    double polyBasisDiv(int iFunc, int iPt) {
      iFunc += num_vertices;
      return v_div_value_n[iPt * dim_v_div + iFunc];
    };
  };
  
  ////////////////////////////////////////////////////////////////////////////////
  // class DirectDGFE
  //    Defined on a poly-element (class PolyElement)
  //    
  //    Gives polynomials of order up to s on the element E, s = r-1, r
  //    Serve as paired space of mixed space that approximate scalar functions
  //      First call initBasis to evaluate all basis functions at a given set of points
  //      Then access the basis functions by function index and pt number
  //
  // Store the result in such a way that the evaluations of all basis functions at
  // one point are put together. The basis functions are stored in the following order:
  //  \Po_{r-1}(E) -> \tilde\Po_r(E) ( polynomials of order r only )
  ////////////////////////////////////////////////////////////////////////////////

  class DirectDGFE
  {
  private:
    DirectMixed* my_dm_space;
    polymesh::PolyElement* my_poly_element;

    int dim_w;

    int num_vertices; // Redundant with my_poly_element
    int polynomial_degree; // Redundant with my_ds_space
    Point ref_origin;
    
    // Pointer to Evaluation storage
    int num_eval_pts;
    double* value_n = nullptr;

    void set_directdgfe(DirectMixed* dmSpace, polymesh::PolyElement* element);

  public:
    DirectDGFE() : my_dm_space(nullptr), my_poly_element(nullptr) {};
    DirectDGFE( DirectMixed* dmSpace, polymesh::PolyElement* element ) { 
      set_directdgfe(dmSpace, element); };
    
    ~DirectDGFE();

    void set(DirectMixed* dmSpace, polymesh::PolyElement* element) {
      set_directdgfe(dmSpace, element); };

    polymesh::PolyElement* elementPtr() const { return my_poly_element; };

    void initBasis(const Point* pt, int num_pts); // (re)evaluate all basis fcns at the points

    // Access basis functions evaluated at pt[iPt]
    double basis(int iFunc, int iPt) const {
      return value_n[iPt * dim_w + iFunc];
    };

    // Get dimensions of spaces
    int dimFull() const { return dim_w; };
    int dimReduced() const { return (polynomial_degree+1) * polynomial_degree /2;}

    // Evaluation \p of full or reduced space at a point
    void eval(const Point* pt, double* result, int num_pts, char type, double* dofs=nullptr);
    void eval(const Point& pt, double& result,char type, double* dofs=nullptr) {
      eval(&pt, &result, 1, type, dofs); };

    double eval(const Point& pt, char type, double* dofs=nullptr) { 
      double result;
      eval(&pt, &result, 1, type, dofs); 
      return result; 
    }

    void eval(const Point* pt, double* fullResult, double* reducedResult, int num_pts, 
              double* full_dofs=nullptr, double* reduced_dofs=nullptr);
    void eval(const Point& pt, double& fullResult, double& reducedResult, 
              double* full_dofs=nullptr, double* reduced_dofs=nullptr) {
      eval(&pt, &fullResult, &reducedResult, 1, full_dofs, reduced_dofs); };
    
    // Output functions
    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };


  ////////////////////////////////////////////////////////////////////////////////
  // class DirectEdgeDGFE
  //
  //    Defined on an interior edge
  //    
  //    Gives polynomials of order up to r on each interior edge of the global mesh
  //
  //    Serve as Lagrange multipliers of mixed space
  //      First call initBasis to evaluate all basis functions at a given set of points
  //      Then access the basis functions by function index and pt number
  //
  // Store the result in such a way that the evaluations of all basis functions at
  // one point are put together. 
  ////////////////////////////////////////////////////////////////////////////////

  class DirectEdgeDGFE
  {
  private:
    DirectMixed* my_dm_space;
    polymesh::Edge* my_edge;

    int dim_l;

    int polynomial_degree; // Redundant with my_dm_space
    
    // Private helper function
    double projToEdge(const Point& p) const;

    // Pointer to Evaluation storage
    int num_eval_pts;
    double* value_n = nullptr;

    void set_directedgedgfe(DirectMixed* dmSpace, polymesh::Edge* edge);

  public:
    DirectEdgeDGFE() : my_dm_space(nullptr), my_edge(nullptr) {};
    DirectEdgeDGFE( DirectMixed* dmSpace, polymesh::Edge* edge ) { 
      set_directedgedgfe(dmSpace, edge); };
    
    ~DirectEdgeDGFE();

    void set(DirectMixed* dmSpace, polymesh::Edge* edge) { set_directedgedgfe(dmSpace, edge); };

    polymesh::Edge* edgePtr() const { return my_edge; };

    void initBasis(const Point* pt, int num_pts); // (re)evaluate all basis fcns at the points

    // Access basis functions evaluated at pt[iPt]
    double basis(int iFunc, int iPt) const {
      return value_n[iPt * dim_l + iFunc];
    };

    // Get dimensions of spaces
    int dim() { return dim_l; };

    // Evaluation \lambda at a point
    void eval(const Point* pt, double* result, int num_pts, double* dofs=nullptr);
    void eval(const Point& pt, double& result, double* dofs=nullptr) {
      eval(&pt, &result, 1, dofs); };

    double eval(const Point& pt, double* dofs=nullptr) { 
      double result;
      eval(&pt, &result, 1, dofs); 
      return result; 
    }
    
    // Output functions
    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };

  ////////////////////////////////////////////////////////////////////////////////
  // class DirectMixed
  //   Defined on a poly-mesh (class PolyMesh), consisting of direct mixed
  //     finite elements (class DirectMixedFE)
  //    
  //   Gives function evaluation from global coefficients for the basis functions
  //
  ////////////////////////////////////////////////////////////////////////////////

  class DirectMixed
  {
  private:
    int polynomial_degree;
    polymesh::PolyMesh* my_mesh;
    bool my_conformity;
    //DirectSerendipity my_high_order_ds_space;

    int num_edges; //redundant with my_mesh
    int num_interior_edges;

    polymesh::Edge* the_dm_edges = nullptr;
    EdgeBCType* the_bc_type = nullptr;
    int* interior_edge_indexing = nullptr;
    int* global_edge_indexing = nullptr;

    DirectDGFE* the_dg_elements = nullptr;

    int* dg_elem_first_to_global_dof_full = nullptr;
    int* dg_elem_first_to_global_dof_reduced = nullptr;

    int dg_dofs_full;
    int dg_dofs_reduced;

    int mixed_dofs_full = 0;
    int mixed_dofs_reduced = 0;

    void set_directmixed(int polyDeg, polymesh::PolyMesh* mesh, bool conforming);
    
  public:
    DirectMixed() : polynomial_degree(-1), my_mesh(nullptr), my_conformity(false), num_edges(0), 
                    num_interior_edges(0), the_dm_edges(nullptr), the_bc_type(nullptr), 
                    interior_edge_indexing(nullptr), global_edge_indexing(nullptr), 
                    the_dg_elements(nullptr),
                    dg_elem_first_to_global_dof_full(nullptr),
                    dg_elem_first_to_global_dof_reduced(nullptr), 
                    dg_dofs_full(0), dg_dofs_reduced(0) { };
    DirectMixed(int polyDeg, polymesh::PolyMesh* mesh, bool conforming) {
      set_directmixed(polyDeg, mesh, conforming); };
    ~DirectMixed();

    void set(int polyDeg, polymesh::PolyMesh* mesh, bool conforming) { set_directmixed(polyDeg, mesh, conforming); };
    
    int nEdges() const { return num_edges; };
    int nInteriorEdges() const { return num_interior_edges; };
    int degPolyn() const { return polynomial_degree; };
    polymesh::PolyMesh* mesh() const { return my_mesh; };
    polymesh::PolyElement* elementPtr(int i) const {return my_mesh->elementPtr(i); };

    bool isConforming() const { return my_conformity; }

    DirectDGFE* DGElementPtr(int i) const { return &the_dg_elements[i]; };
    EdgeBCType bcType(int i) const { return the_bc_type[i]; }; //Boundary type for each edge

    
    // Return the interior edge indexing from global edge indexing
    int glob_to_int(int i) const { return interior_edge_indexing[i]; };

    // Return the global edge indexing from interior edge indexing
    int int_to_glob(int i) const {return global_edge_indexing[i]; }; 

    // Return the interior edge indexing of the iEdge-th edge of the iElement-th element
    int interiorEdgeIndex(int iElement, int iEdge) const { 
      return interior_edge_indexing[my_mesh->elementPtr(iElement)->edgePtr(iEdge)->meshIndex()]; 
    };

    int globalEdgeIndex(int iElement, int iEdge) const {
      return my_mesh->elementPtr(iElement)->edgePtr(iEdge)->meshIndex();
    }

    // Return the number of global dofs of spaces
    // Gives number of global dofs
    int nMixedDoFs(char type) const {
      return (type == 'f')? mixed_dofs_full : mixed_dofs_reduced;
    };

    int nDGDoFs(char type) const;

    // Map the (first local dof of) element/edge to global dof
    int dg_Elem_First_To_Global_Dof(int i, char type) const;

    // Return the edge from interior edge index
    polymesh::Edge* edgeInteriorPtr(int i) const {
      return my_mesh -> edgePtr(int_to_glob(i));
    };

    polymesh::Edge* edgePtr(int i) const {
      return my_mesh -> edgePtr(i);
    };

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;

    int write_matlab(std::string& filename) const;


    friend class DirectMixedFE;
    friend class DirectDGFE;
    friend class DirectEdgeDGFE;
    friend class DirectMixedArray;
    friend class DirectDGArray;
    friend class DirectEdgeDGArray;

    friend class DirectMixedHybrid;
    friend class DirectMixedConf;

    friend class DirectMixedHybridArray;
    friend class DirectMixedConfArray;
  };

  class DirectMixedHybrid : public DirectMixed
  {
  private:
    DirectMixedHybridFE* the_dm_elements = nullptr;
    DirectEdgeDGFE* the_dg_edge_elements = nullptr;

    int* mixed_elem_first_to_global_dof_full = nullptr;
    int* mixed_elem_first_to_global_dof_reduced = nullptr;
    int* edge_elem_first_to_global_dof = nullptr;

    int dg_edge_dofs = 0; // Here we also include DoFs on the boundary
    int dg_int_edge_dofs = 0;

    void set_directmixedhybrid(int polyDeg, polymesh::PolyMesh* mesh);
    
  public:
    DirectMixedHybrid(const DirectMixed& dm) {
      set_directmixedhybrid(dm.degPolyn(), dm.mesh());
    }
    DirectMixedHybrid(int polyDeg, polymesh::PolyMesh* mesh) {
      set_directmixedhybrid(polyDeg, mesh); };
    ~DirectMixedHybrid();

    void set(int polyDeg, polymesh::PolyMesh* mesh) {
      set_directmixedhybrid(polyDeg, mesh); };

    DirectMixedHybridFE* MixedElementPtr(int i) const { return &the_dm_elements[i]; };
    DirectEdgeDGFE* DGEdgePtr(int i) const { return &the_dg_edge_elements[i]; };
    DirectEdgeDGFE* DGEdgeInteriorPtr(int i) const { return &the_dg_edge_elements[int_to_glob(i)]; };

    // Return the number of global dofs of spaces

    int nEdgeDGDoFs() const;
    int nIntEdgeDGDoFs() const;

    // Map the (first local dof of) element/edge to global dof
    int mixed_Elem_First_To_Global_Dof(int i, char type) const;
    int edge_Elem_First_To_Global_Dof(int i) const;

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };

  class DirectMixedConf : public DirectMixed
  {
    private:
    DirectMixedConfFE* the_dm_elements = nullptr;

    MixedDoFType* the_dof_type_full = nullptr;
    MixedDoFType* the_dof_type_reduced = nullptr;
    
    MixedDoFBCType* the_dof_bc_type_full = nullptr;
    MixedDoFBCType* the_dof_bc_type_reduced = nullptr;

    // Map local vertex basis functions to global dof
    int** vertex_loc_to_glob_full = nullptr;
    int** vertex_loc_to_glob_coeff_full = nullptr;

    // Map the global vertex index to corresponding global dof
    int* mesh_vertex_to_dof_index = nullptr;

    // Map local edge basis functions to global dof
    int** edge_loc_to_glob_full = nullptr;

    // Map local cell basis functions to global dof
    int** cell_loc_to_glob_full = nullptr;

    int** poly_loc_to_glob_full = nullptr;
    int** poly_loc_to_glob_reduced = nullptr;

    int** loc_to_glob_full = nullptr;
    int** loc_to_glob_reduced = nullptr;
    
    void set_directmixedconf(int polyDeg, polymesh::PolyMesh* mesh);
    
    public:
    DirectMixedConf(const DirectMixed& dm) {
      set_directmixedconf(dm.degPolyn(), dm.mesh());
    }
    DirectMixedConf(int polyDeg, polymesh::PolyMesh* mesh) {
      set_directmixedconf(polyDeg, mesh);
    }
    ~DirectMixedConf();

    void set(int polyDeg, polymesh::PolyMesh* mesh) {
      set_directmixedconf(polyDeg, mesh);
    }

    DirectMixedConfFE* MixedElementPtr(int i) const { return &the_dm_elements[i]; };

    // Return the type and BC type of a global dof with indexing
    MixedDoFType global_dof_type(int i, char type) const {
      return (type == 'f')? the_dof_type_full[i] : the_dof_type_reduced[i];
    };

    MixedDoFBCType global_dof_bc_type(int i, char type) const {
      return (type == 'f')? the_dof_bc_type_full[i] : the_dof_bc_type_reduced[i];
    };

    // Map local vertex basis functions 
    // (the basis functions with outer normal fluxes) 
    // to global dof index
    int vertex_loc_to_glob(int nElement, int iVertex) const {
      return vertex_loc_to_glob_full[nElement][iVertex];
    };

    // Gives this function counts for positive part or negative part
    // for this global dof
    // If local edge[iVertex] has the same orientation with the
    // corresponding global edge, then it gives 1, otherwise -1  
    int vertex_loc_to_glob_coeff(int nElement, int iVertex) const {
      return vertex_loc_to_glob_coeff_full[nElement][iVertex];
    };


    // Map local edge basis functions to global dof index
    int edge_loc_to_glob(int nElement, int iEdge, int jNode) const {
      return edge_loc_to_glob_full[nElement][iEdge*polynomial_degree+jNode];
    };
    
    // Map local cell bubble functions to global dof index
    int cell_loc_to_glob(int nElement, int i) const {
      return cell_loc_to_glob_full[nElement][i];
    };

    // Map local \x\Po bubble functions to global dof index
    int poly_loc_to_glob(int nElement, int i, char type) const {
      return (type == 'f')? poly_loc_to_glob_full[nElement][i] : poly_loc_to_glob_reduced[nElement][i];
    };

    // Map local dof of an element to global dof index
    int loc_to_glob(int nElement, int i, char type) const {
      return (type == 'f')? loc_to_glob_full[nElement][i] : loc_to_glob_reduced[nElement][i];
    };

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };

  ////////////////////////////////////////////////////////////////////////////////
  // class DirectMixedArray
  //    Gives values for each coefficient of basis function
  ////////////////////////////////////////////////////////////////////////////////

  class DirectMixedArray
  {
  private:

    int num_dofs = 0;
    int num_elements = 0; // Redundant with my_dm_space
    char space_type; // either 'f' or 'r'

    double* the_array = nullptr;
    polymesh::PolyMesh* my_mesh = nullptr; // Redundant with my_dm_space
    DirectMixed* my_dm_space = nullptr;

    void set_directmixedarray(DirectMixed* dmSpace, char spacetype);

  protected:

    DirectMixedArray() : space_type('f') {};
    DirectMixedArray(DirectMixed* dmSpace, char spacetype) {
      set_directmixedarray(dmSpace, spacetype); };
    ~DirectMixedArray();
    
    void set(DirectMixed* dmSpace, char spacetype) { set_directmixedarray(dmSpace, spacetype); };
  
  public:
    DirectMixed* dmSpace() const { return my_dm_space; }
    char spaceType() const { return space_type; }
    int size() const { return num_dofs; };

    double& operator() (int i)       { return the_array[i]; }
    double  operator() (int i) const { return the_array[i]; }
    double& operator[] (int i)       { return the_array[i]; }
    double  operator[] (int i) const { return the_array[i]; }

    virtual void eval(const Point* pts, Tensor1* result, int num_pts) const { };
    virtual void eval(const Point& pt, Tensor1& result) const { };
    Tensor1 eval(const Point& pt) const {
      Tensor1 result; eval(pt, result); return result;
    };

    virtual void eval_div(const Point* pts, double* result, int num_pts) const { };
    virtual void eval_div(const Point& pt, double& result) const { };
    double eval_div(const Point& pt) const {
      double result; eval_div(pt, result); return result;
    };

    void l2normError(double& l2Error, double& l2Norm, int refinement_level = 0, Tensor1 (*referenceFcn)(double,double) = nullptr);
    void l2norm(double& l2Norm, int refinement_level = 0) {
      double null; l2normError(l2Norm,null, refinement_level); };

    // Reference function should be
    void l2normError_div(double& l2Error, double& l2Norm, int refinement_level = 0, double (*referenceFcn)(double,double) = nullptr);
    void l2norm_div(double& l2Norm, int refinement_level = 0) {
      double null; l2normError_div(l2Norm,null, refinement_level); };

    void write_matlab_mesh(std::ofstream* fout, int num_pts_x, int num_pts_y) const;
    void write_matlab_mesh(std::ofstream& fout, int num_pts_x, int num_pts_y) const {
      write_matlab_mesh(&fout, num_pts_x, num_pts_y); };

    int write_matlab_mesh(std::string& filename, int num_pts_x, int num_pts_y) const;

    void write_matlab_mesh_by_pt(std::ofstream& fout, int num_pts_x, int num_pts_y) const;
    int write_matlab_mesh_by_pt(std::string& filename, int num_pts_x, int num_pts_y) const;

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;

    void write_matlab_mesh_div(std::ofstream* fout, int num_pts_x, int num_pts_y) const;
    void write_matlab_mesh_div(std::ofstream& fout, int num_pts_x, int num_pts_y) const {
      write_matlab_mesh_div(&fout, num_pts_x, num_pts_y); };

    int write_matlab_mesh_div(std::string& filename, int num_pts_x, int num_pts_y) const;

    void write_matlab_mesh_error(std::ofstream* fout, int num_pts_x, int num_pts_y, Tensor1 (*referenceFcn)(double,double) = nullptr) const;
    void write_matlab_mesh_error(std::ofstream& fout, int num_pts_x, int num_pts_y, Tensor1 (*referenceFcn)(double,double) = nullptr) const {
      write_matlab_mesh_error(&fout, num_pts_x, num_pts_y, referenceFcn); };

    int write_matlab_mesh_error(std::string& filename, int num_pts_x, int num_pts_y, Tensor1 (*referenceFcn)(double,double) = nullptr) const;

    void write_matlab_mesh_div_error(std::ofstream* fout, int num_pts_x, int num_pts_y, double (*referenceFcn)(double,double) = nullptr) const;
    void write_matlab_mesh_div_error(std::ofstream& fout, int num_pts_x, int num_pts_y, double (*referenceFcn)(double,double) = nullptr) const {
      write_matlab_mesh_div_error(&fout, num_pts_x, num_pts_y, referenceFcn); };

    int write_matlab_mesh_div_error(std::string& filename, int num_pts_x, int num_pts_y, double (*referenceFcn)(double,double) = nullptr) const;

    friend class DirectMixedHybridArray;
    friend class DirectMixedConfArray;
  };

  class DirectMixedHybridArray : public DirectMixedArray {
  private:
    DirectMixedHybrid* my_dm_space;

    void set_directmixedhybridarray(DirectMixedHybrid* dmSpace, char spacetype);
  public:

    DirectMixedHybridArray(DirectMixedHybrid* dmSpace, char spacetype) {
    set_directmixedhybridarray(dmSpace, spacetype); };

    void set(DirectMixedHybrid* dmSpace, char spacetype) {
    set_directmixedhybridarray(dmSpace, spacetype); };

    DirectMixedHybrid* dmSpace() const { return my_dm_space; };

    void eval(const Point* pts, Tensor1* result, int num_pts) const;
    void eval(const Point& pt, Tensor1& result) const;

    void eval_div(const Point* pts, double* result, int num_pts) const;
    void eval_div(const Point& pt, double& result) const;
  };




  class DirectMixedConfArray : public DirectMixedArray {
  private:
    DirectMixedConf* my_dm_space;
    void set_directmixedconfarray(DirectMixedConf* dmSpace, char spacetype);

  public:
    DirectMixedConfArray(DirectMixedConf* dmSpace, char spacetype) {
    set_directmixedconfarray(dmSpace, spacetype); };

    void set(DirectMixedConf* dmSpace, char spacetype) {
    set_directmixedconfarray(dmSpace, spacetype); };

    DirectMixedConf* dmSpace() const { return my_dm_space; };

    void eval(const Point* pts, Tensor1* result, int num_pts) const;
    void eval(const Point& pt, Tensor1& result) const;

    void eval_div(const Point* pts, double* result, int num_pts) const;
    void eval_div(const Point& pt, double& result) const;
  };

  ////////////////////////////////////////////////////////////////////////////////
  // class DirectDGArray
  //    Gives values for each coefficient of basis function
  ////////////////////////////////////////////////////////////////////////////////

  class DirectDGArray
  {
  private:
    int num_dofs;
    int num_elements; // Redundant with my_dm_space
    char space_type; // either 'f' or 'r'

    double* the_array = nullptr;

    DirectMixed* my_dm_space;

    void set_directdgarray(DirectMixed* dmSpace, char spacetype);

  public:
    DirectDGArray() : num_dofs(0), num_elements(0), space_type('f'), the_array(nullptr), my_dm_space(nullptr) {};
    DirectDGArray(DirectMixed* dmSpace, char spacetype) {
      set_directdgarray(dmSpace, spacetype); };
    DirectDGArray(const DirectDGArray& a) : the_array(nullptr) {
      set_directdgarray(a.dmSpace(), a.spaceType()); };
    ~DirectDGArray();
    
    void set(DirectMixed* dmSpace, char spacetype) { set_directdgarray(dmSpace, spacetype); };

    DirectMixed* dmSpace() const { return my_dm_space; };
    char spaceType() const { return space_type; }
    int size() const { return num_dofs; };

    double& operator() (int i)       { return the_array[i]; }
    double  operator() (int i) const { return the_array[i]; }
    double& operator[] (int i)       { return the_array[i]; }
    double  operator[] (int i) const { return the_array[i]; }

    void eval(const Point* pts, double* result, int num_pts) const;
    void eval(const Point& pt, double& result) const;
    double eval(const Point& pt) const {
      double result; eval(pt, result); return result;
    };

    void l2normError(double& l2Error, double& l2Norm, int refinement_level = 0, double (*referenceFcn)(double,double) = nullptr);
    void l2norm(double& l2Norm, int refinement_level = 0) {
      double null; l2normError(l2Norm,null,refinement_level); };

    void write_matlab_mesh(std::ofstream* fout, int num_pts_x, int num_pts_y) const;
    void write_matlab_mesh(std::ofstream& fout, int num_pts_x, int num_pts_y) const {
      write_matlab_mesh(&fout, num_pts_x, num_pts_y); };

    int write_matlab_mesh(std::string& filename, int num_pts_x, int num_pts_y) const;

    void write_matlab_mesh_by_pt(std::ofstream& fout, int num_pts_x, int num_pts_y) const;
    int write_matlab_mesh_by_pt(std::string& filename, int num_pts_x, int num_pts_y) const;

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;

    void write_matlab_mesh_error(std::ofstream* fout, int num_pts_x, int num_pts_y, double (*referenceFcn)(double,double) = nullptr) const;
    void write_matlab_mesh_error(std::ofstream& fout, int num_pts_x, int num_pts_y, double (*referenceFcn)(double,double) = nullptr) const {
      write_matlab_mesh_error(&fout, num_pts_x, num_pts_y, referenceFcn); };

    int write_matlab_mesh_error(std::string& filename, int num_pts_x, int num_pts_y, double (*referenceFcn)(double,double) = nullptr) const;
  };

  ////////////////////////////////////////////////////////////////////////////////
  // class DirectEdgeDGArray
  //    Gives values for each coefficient of basis function
  //    Here we do not distinguish boundary and
  ////////////////////////////////////////////////////////////////////////////////

  class DirectEdgeDGArray
  {
  private:
    int num_dofs;
    int num_edges; // Redundant with my_dm_space

    double* the_array = nullptr;

    DirectMixedHybrid* my_dm_space;

    void set_directedgedgarray(DirectMixedHybrid* dmSpace);

  public:
    DirectEdgeDGArray() : num_dofs(0), num_edges(0), the_array(nullptr), my_dm_space(nullptr) {};
    DirectEdgeDGArray(DirectMixedHybrid* dmSpace) {
      set_directedgedgarray(dmSpace); };
    DirectEdgeDGArray(const DirectEdgeDGArray& a) : the_array(nullptr) {
      set_directedgedgarray(a.dmSpace()); };
    ~DirectEdgeDGArray();
    
    void set(DirectMixedHybrid* dmSpace) { set_directedgedgarray(dmSpace); };

    DirectMixedHybrid* dmSpace() const { return my_dm_space; };
    int size() const { return num_dofs; };

    double& operator() (int i)       { return the_array[i]; }
    double  operator() (int i) const { return the_array[i]; }
    double& operator[] (int i)       { return the_array[i]; }
    double  operator[] (int i) const { return the_array[i]; }

    void eval(const Point* pts, double* result, int num_pts) const;
    void eval(const Point& pt, double& result) const;
    double eval(const Point& pt) const {
      double result; eval(pt, result); return result;
    };

    void l2normError(double& l2Error, double& l2Norm, double (*referenceFcn)(double,double) = nullptr);
    void l2norm(double& l2Norm) {
      double null; l2normError(l2Norm,null); };

    void write_matlab_mesh(std::ofstream* fout, int num_pts_x, int num_pts_y) const;
    void write_matlab_mesh(std::ofstream& fout, int num_pts_x, int num_pts_y) const {
      write_matlab_mesh(&fout, num_pts_x, num_pts_y); };

    int write_matlab_mesh(std::string& filename, int num_pts_x, int num_pts_y) const;

    void write_matlab_mesh_by_pt(std::ofstream& fout, int num_pts_x, int num_pts_y) const;
    int write_matlab_mesh_by_pt(std::string& filename, int num_pts_x, int num_pts_y) const;

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };
};
#endif