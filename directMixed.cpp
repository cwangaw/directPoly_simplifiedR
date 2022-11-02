#include <cmath>
#include <iostream>
#include <vector>
#include <assert.h>
#include <stdio.h>
#include <assert.h>
#include <algorithm>

using namespace std;

#include "Utilities/debug.h"
#include "Mesh/polyMesh.h"
using namespace polymesh;
#include "directMixed.h"
using namespace directserendipity;
#include "polyQuadrature.h"
using namespace polyquadrature;


////////////////////////////////////////////////////////////////////////////////
// Class DirectMixed

void DirectMixed::set_directmixed(int polyDeg, polymesh::PolyMesh* mesh, bool conforming) {
  polynomial_degree = polyDeg;
  my_mesh = mesh;
  my_conformity = conforming;
  num_edges = my_mesh -> nEdges();
  
  int index = 0;
  for (int i = 0; i < nEdges(); i++) {
    if (!my_mesh->edgePtr(i)->isOnBoundary()) index++;
  }
  num_interior_edges = index;

  // ALLOCATE ELEMENTS

  // Allocate
  
  if(the_dm_edges) delete[] the_dm_edges;
  the_dm_edges = new Edge[num_edges];
  
  if(the_bc_type) delete[] the_bc_type;
  the_bc_type = new EdgeBCType[num_edges];

  if(the_dg_elements) delete[] the_dg_elements;
  the_dg_elements = new DirectDGFE[my_mesh->nElements()];

  if(interior_edge_indexing) delete[] interior_edge_indexing;
  interior_edge_indexing = new int[num_edges];

  if(global_edge_indexing) delete[] global_edge_indexing;
  global_edge_indexing = new int[num_interior_edges];
  
  if(dg_elem_first_to_global_dof_full) delete[] dg_elem_first_to_global_dof_full;
  dg_elem_first_to_global_dof_full = new int[my_mesh->nElements()];

  if(dg_elem_first_to_global_dof_reduced) delete[] dg_elem_first_to_global_dof_reduced;
  dg_elem_first_to_global_dof_reduced = new int[my_mesh->nElements()];

  // DETERMINE EDGES (ORDERING AND TYPES)

  // interior_edge_indexing maps from global edge index to interior edge index
  // global_edge_indexing maps from interior edge index to global edge index
  index = 0;
  for (int i = 0; i < nEdges(); i++) {
    if (my_mesh->edgePtr(i)->isOnBoundary()) {
      the_bc_type[i] = EdgeBCType::boundary;
      interior_edge_indexing[i] = -1;
    } else {
      the_bc_type[i] = EdgeBCType::interior;
      interior_edge_indexing[i] = index;
      global_edge_indexing[index] = i;
      index++;
    }
  }

  // ALLOCATE FINITE ELEMENTS AND DECIDE GLOBAL DOFS
  dg_dofs_full = 0;
  dg_dofs_reduced = 0;

  for(int iElement=0; iElement<my_mesh->nElements(); iElement++) {
    dg_elem_first_to_global_dof_full[iElement] = dg_dofs_full;
    dg_elem_first_to_global_dof_reduced[iElement] = dg_dofs_reduced;

    PolyElement* element = my_mesh->elementPtr(iElement);
    the_dg_elements[iElement].set(this, element);

    dg_dofs_full += the_dg_elements[iElement].dimFull();
    dg_dofs_reduced += the_dg_elements[iElement].dimReduced();
  }
};

DirectMixed::~DirectMixed() {
    if (the_dm_edges) delete[] the_dm_edges;
    if (the_bc_type) delete[] the_bc_type;
    if (interior_edge_indexing) delete[] interior_edge_indexing;
    if (global_edge_indexing) delete[] global_edge_indexing;
    if (the_dg_elements) delete[] the_dg_elements;
    if (dg_elem_first_to_global_dof_full) delete[] dg_elem_first_to_global_dof_full;
    if (dg_elem_first_to_global_dof_reduced) delete[] dg_elem_first_to_global_dof_reduced;
}



int DirectMixed::nDGDoFs(char type) const {
  return (type == 'f') ?   dg_dofs_full : dg_dofs_reduced;
};


    
int DirectMixed::dg_Elem_First_To_Global_Dof(int i, char type) const {
  return (type == 'f') ? dg_elem_first_to_global_dof_full[i] : dg_elem_first_to_global_dof_reduced[i];
}



void DirectMixed::write_raw(std::ofstream& fout) const {
  fout << "DIRECT MIXED SPACE: PRESSURE SPACE\n";
  fout << "polynomial_degree      = " <<  polynomial_degree << "\n";
  fout << "my_mesh                = " << my_mesh << "\n";
  fout << "num_dofs of W_full     = " << dg_dofs_full << "\n";
  fout << "num_dofs of W_reduced  = " << dg_dofs_reduced << "\n";


  fout << "\ndmSpace edges:\n";
  for(int i=0; i<num_edges; i++) {
    fout << "  Edge " << i << " (type ";
    switch(the_bc_type[i]) {
    case EdgeBCType::interior: { fout << "interior"; break; }
    case EdgeBCType::boundary:   { fout << "boundary"; break; }
    default: { fout << "???"; break; }
    }
    fout << ")\n";
    
  }


  fout << "\n";
  fout << "\ndg_elem_first_to_global_dof_full:\n";
  for(int i=0; i<my_mesh->nElements(); i++) {
    fout << "  " << dg_elem_first_to_global_dof_full[i];
  }
  fout << "\n";
    fout << "\ndg_elem_first_to_global_dof_reduced:\n";
  for(int i=0; i<my_mesh->nElements(); i++) {
    fout << "  " << dg_elem_first_to_global_dof_reduced[i];
  }
  fout << "\n";
 

  fout << "\nDG elements:\n";
  for(int i=0; i<my_mesh->nElements(); i++) {
    fout << " Element " << i << "\n";
    the_dg_elements[i].write_raw(fout);
  }
  fout << "\n";


}

int DirectMixed::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}

int DirectMixed::write_matlab(std::string& filename) const {
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

////////////////////////////////////////////////////////////////////////////////
// Class DirectMixedHybrid

void DirectMixedHybrid::set_directmixedhybrid(int polyDeg, polymesh::PolyMesh* mesh) {
  set_directmixed(polyDeg, mesh, false);
  // ALLOCATE ELEMENTS


  // Allocate

  if(the_dm_elements) delete[] the_dm_elements;
  the_dm_elements = new DirectMixedHybridFE[my_mesh->nElements()];

  if(the_dg_edge_elements) delete[] the_dg_edge_elements;
  the_dg_edge_elements = new DirectEdgeDGFE[num_edges];

  if(mixed_elem_first_to_global_dof_full) delete[] mixed_elem_first_to_global_dof_full;
  mixed_elem_first_to_global_dof_full = new int[my_mesh->nElements()];

  if(mixed_elem_first_to_global_dof_reduced) delete[] mixed_elem_first_to_global_dof_reduced;
  mixed_elem_first_to_global_dof_reduced = new int[my_mesh->nElements()];
  
  if(edge_elem_first_to_global_dof) delete[] edge_elem_first_to_global_dof;
  edge_elem_first_to_global_dof = new int[num_edges];
  
  // ALLOCATE FINITE ELEMENTS AND DECIDE GLOBAL DOFS

  mixed_dofs_full = 0;
  mixed_dofs_reduced = 0;
  dg_edge_dofs = 0;

  for(int iElement=0; iElement<my_mesh->nElements(); iElement++) {
    mixed_elem_first_to_global_dof_full[iElement] = mixed_dofs_full;
    mixed_elem_first_to_global_dof_reduced[iElement] = mixed_dofs_reduced;

    PolyElement* element = my_mesh->elementPtr(iElement);

    the_dm_elements[iElement].set(this, element);

    mixed_dofs_full += the_dm_elements[iElement].dimVFull();
    mixed_dofs_reduced += the_dm_elements[iElement].dimVReduced();
  }
  
  for (int iEdge = 0; iEdge < num_edges; iEdge++ ) {
    edge_elem_first_to_global_dof[iEdge] = dg_edge_dofs;

    Edge* edge = my_mesh->edgePtr(iEdge);
    the_dg_edge_elements[iEdge].set(this,edge);
    dg_edge_dofs += the_dg_edge_elements[iEdge].dim();
    if (the_bc_type[iEdge] == EdgeBCType::interior) dg_int_edge_dofs += the_dg_edge_elements[iEdge].dim();
  }
};

DirectMixedHybrid::~DirectMixedHybrid() {
  if (the_dm_elements) delete[] the_dm_elements;
  if (the_dg_edge_elements) delete[] the_dg_edge_elements;
  if (mixed_elem_first_to_global_dof_full) delete[] mixed_elem_first_to_global_dof_full;
  if (mixed_elem_first_to_global_dof_reduced) delete[] mixed_elem_first_to_global_dof_reduced;
  if (edge_elem_first_to_global_dof) delete[] edge_elem_first_to_global_dof;
}

int DirectMixedHybrid::nEdgeDGDoFs() const {
  return dg_edge_dofs;
};

int DirectMixedHybrid::nIntEdgeDGDoFs() const {
  return dg_int_edge_dofs;
};

int DirectMixedHybrid::mixed_Elem_First_To_Global_Dof(int i, char type) const {
  return (type == 'f') ? mixed_elem_first_to_global_dof_full[i] : mixed_elem_first_to_global_dof_reduced[i];
}

int DirectMixedHybrid::edge_Elem_First_To_Global_Dof(int i) const {
  return edge_elem_first_to_global_dof[i];
}

void DirectMixedHybrid::write_raw(std::ofstream& fout) const { 
  DirectMixed::write_raw(fout);

  fout << "DIRECT MIXED HYBRID SPACE AND LAGRANGE MULTIPLYER:\n";
  fout << "num_dofs of V_full     = " << mixed_dofs_full << "\n";
  fout << "num_dofs of V_reduced  = " << mixed_dofs_reduced << "\n";
  fout << "num_dofs of Lambda     = " << dg_edge_dofs << "\n";

  fout << "\ndmSpace indexing\n";
  fout << "\nmixed_elem_first_to_global_dof_full:\n";
  for(int i=0; i<my_mesh->nElements(); i++) {
    fout << "  " << mixed_elem_first_to_global_dof_full[i];
  }
  fout << "\n";

  fout << "\nmixed_elem_first_to_global_dof_reduced:\n";
  for(int i=0; i<my_mesh->nElements(); i++) {
    fout << "  " << mixed_elem_first_to_global_dof_reduced[i];
  }
  fout << "\n";

  fout << "\nedge_elem_first_to_global_dof:\n";
  for(int i=0; i<my_mesh->nEdges(); i++) {
    fout << "  " << edge_elem_first_to_global_dof[i];
  }
  fout << "\n";

  fout << "\nMixed elements:\n";
  for(int i=0; i<my_mesh->nElements(); i++) {
    fout << " Element " << i << "\n";
    the_dm_elements[i].write_raw(fout);
  }
  fout << "\n";

  fout << "\nEdge DG elements:\n";
  for(int i=0; i<my_mesh->nEdges(); i++) {
    fout << " Edge " << i << "\n";
    the_dg_edge_elements[i].write_raw(fout);
  }
  fout << "\n";
};

int DirectMixedHybrid::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Class DirectMixedConf

void DirectMixedConf::set_directmixedconf(int polyDeg, polymesh::PolyMesh* mesh) {
  set_directmixed(polyDeg, mesh, true);

  int num_cell_dofs = 0;
  
  for (int iElement = 0; iElement < my_mesh->nElements(); iElement++) {
    num_cell_dofs += max( 0, (polynomial_degree - my_mesh->elementPtr(iElement)->nVertices() + 3)
                            * (polynomial_degree - my_mesh->elementPtr(iElement)->nVertices() + 2) / 2);
  }


  mixed_dofs_full = nEdges() * (polynomial_degree + 1) // number of vertex and edge basis functions
                  + num_cell_dofs 
                  + my_mesh->nElements() * ( (polynomial_degree + 2) * (polynomial_degree + 1) / 2 - 1 );

  mixed_dofs_reduced = nEdges() * (polynomial_degree + 1) // number of vertex and edge basis functions
                     + num_cell_dofs 
                     + my_mesh->nElements() * ( (polynomial_degree + 1) * polynomial_degree / 2 - 1 );
  
  // Allocate 

  if (the_dm_elements) delete[] the_dm_elements;
  the_dm_elements = new DirectMixedConfFE[my_mesh->nElements()];

  if (the_dof_type_full) delete[] the_dof_type_full;
  the_dof_type_full = new MixedDoFType[mixed_dofs_full];

  if (the_dof_type_reduced) delete[] the_dof_type_reduced;
  the_dof_type_reduced = new MixedDoFType[mixed_dofs_reduced];

  if (the_dof_bc_type_full) delete[] the_dof_bc_type_full;
  the_dof_bc_type_full = new MixedDoFBCType[mixed_dofs_full];

  if (the_dof_bc_type_reduced) delete[] the_dof_bc_type_reduced;
  the_dof_bc_type_reduced = new MixedDoFBCType[mixed_dofs_reduced];

  if (vertex_loc_to_glob_full) delete[] vertex_loc_to_glob_full;
  vertex_loc_to_glob_full = new int*[my_mesh->nElements()];

  if (vertex_loc_to_glob_coeff_full) delete[] vertex_loc_to_glob_coeff_full;
  vertex_loc_to_glob_coeff_full = new int*[my_mesh->nElements()];


  // Globally it is not really corresponding to vertices
  // But to edge indexes
  if (mesh_vertex_to_dof_index) delete[] mesh_vertex_to_dof_index;
  mesh_vertex_to_dof_index = new int[my_mesh->nEdges()];

  if (edge_loc_to_glob_full) delete[] edge_loc_to_glob_full;
  edge_loc_to_glob_full = new int*[my_mesh->nElements()];

  if (cell_loc_to_glob_full) delete[] cell_loc_to_glob_full;
  cell_loc_to_glob_full = new int*[my_mesh->nElements()];

  if (poly_loc_to_glob_full) delete[] poly_loc_to_glob_full;
  poly_loc_to_glob_full = new int*[my_mesh->nElements()];

  if (poly_loc_to_glob_reduced) delete[] poly_loc_to_glob_reduced;
  poly_loc_to_glob_reduced = new int*[my_mesh->nElements()];

  if (loc_to_glob_full) delete[] loc_to_glob_full;
  loc_to_glob_full = new int*[my_mesh->nElements()];

  if (loc_to_glob_reduced) delete[] loc_to_glob_reduced;
  loc_to_glob_reduced = new int*[my_mesh->nElements()];

  // SET UP THE ELEMENTS
  for(int iElement=0; iElement<my_mesh->nElements(); iElement++) {
    PolyElement* element = my_mesh->elementPtr(iElement);
    the_dm_elements[iElement].set(this, element);
    // Initialize loc_to_glob array
    loc_to_glob_full[iElement] = new int[the_dm_elements[iElement].dimVFull()];
    loc_to_glob_reduced[iElement] = new int[the_dm_elements[iElement].dimVReduced()];
  }

  // DETERMINE DOFS (ORDERING AND TYPES)
  
  // Global Dofs are stored in the following order: 
  // Edge paf dof for Global Edge 0 -> Edge dofs for Global Edge 0 
  // -> Edge paf dof for Global Edge 1 -> Edge dofs for Global Edge 1
  // -> ... -> Edge paf dof for Global Edge ( my_mesh -> nEdges()-1 ) 
  // -> Edge dofs for Global Edge ( my_mesh -> nEdges()-1 ) 
  // -> Cell dofs for Element 0 -> Cell dofs for Element 1
  // -> ... -> Cell dofs for Element ( my_mesh -> nElements()-1 )
  // -> Poly dofs for Element 0 -> Poly dofs for Element 1
  // -> ... -> Poly dofs for Element ( my_mesh -> nElements()-1 )

  // Loop through elements

  // Loop through element vertices to set up global edge paf and global edge dofs
  // (order by lowest element number)
  int dof_num = -1;

  for(int iElement = 0; iElement < my_mesh->nElements(); iElement++) {
    int nGon = my_mesh->elementPtr(iElement)->nVertices();
    int local_num_cell_dofs = max( 0, (polynomial_degree - nGon + 3)
                            * (polynomial_degree - nGon + 2) / 2);
    vertex_loc_to_glob_full[iElement] = new int[nGon];
    vertex_loc_to_glob_coeff_full[iElement] = new int[nGon];
    edge_loc_to_glob_full[iElement] = new int[nGon*polynomial_degree];

    for(int iVertex = 0; iVertex < nGon; iVertex++) {
      Edge* e = my_mesh->edgePtr( my_mesh->elementPtr(iElement)->edgePtr(iVertex)->meshIndex() );
	    int edgeNbrIndices[2];
	    e -> nbrElements(edgeNbrIndices);
      // If local edge has same orientation as corresponding global edge, it is 1; otherwise it is zero
      bool local_edge_orientation = (e->vertexPtr(1) == my_mesh->elementPtr(iElement)->vertexPtr(iVertex));

      // Add edge nodes (if new, order by edge, not oriented edge)
      bool newEdge = false;
      if( (edgeNbrIndices[0] == -1 || edgeNbrIndices[1] == -1
          || edgeNbrIndices[0] > iElement || edgeNbrIndices[1] > iElement) )
          { newEdge = true; }

      // Add vertex nodes (if new)

      if (newEdge) {
        dof_num++;
        the_dof_type_full[dof_num] = MixedDoFType::vertex;
        the_dof_bc_type_full[dof_num] = (e->isOnBoundary()) ? MixedDoFBCType::boundary : MixedDoFBCType::interior;
        the_dof_type_reduced[dof_num] = MixedDoFType::vertex;
        the_dof_bc_type_reduced[dof_num] = (e->isOnBoundary()) ? MixedDoFBCType::boundary : MixedDoFBCType::interior;
        mesh_vertex_to_dof_index[e->meshIndex()] = dof_num;
      }

      // Decide vertex_loc_to_glob
      int v_loc_to_glob = (newEdge) ? dof_num : mesh_vertex_to_dof_index[e->meshIndex()];
      vertex_loc_to_glob_full[iElement][iVertex] = v_loc_to_glob;
      loc_to_glob_full[iElement][iVertex + nGon * degPolyn() + local_num_cell_dofs] = v_loc_to_glob;
      loc_to_glob_reduced[iElement][iVertex + nGon * degPolyn() + local_num_cell_dofs] = v_loc_to_glob;

      // Decide vertex_loc_to_glob_coeff by the orientation of local edge(iVertex)
      vertex_loc_to_glob_coeff_full[iElement][iVertex] = (local_edge_orientation) ? 1 : -1;

      // Deal with edge
      for(int j = 0; j < polynomial_degree; j++) {
        // Add new edge dofs following the corresponding edge paf dof
        if (newEdge) {
          dof_num++;
          the_dof_type_full[dof_num] = MixedDoFType::edge;
          the_dof_bc_type_full[dof_num] = (e->isOnBoundary()) ? MixedDoFBCType::boundary : MixedDoFBCType::interior;
          the_dof_type_reduced[dof_num] = MixedDoFType::edge;
          the_dof_bc_type_reduced[dof_num] = (e->isOnBoundary()) ? MixedDoFBCType::boundary : MixedDoFBCType::interior;
        }

        // Decide edge_loc_to_glob
          edge_loc_to_glob_full[iElement][iVertex*polynomial_degree+j] = v_loc_to_glob + j + 1;
          loc_to_glob_full[iElement][local_num_cell_dofs+iVertex*polynomial_degree+j] = v_loc_to_glob + j + 1;
          loc_to_glob_reduced[iElement][local_num_cell_dofs+iVertex*polynomial_degree+j] = v_loc_to_glob + j + 1;
      }

    }
  }

  // Loop through elements to set up cell dofs

  if (num_cell_dofs > 0) {
    for(int iElement = 0; iElement < my_mesh->nElements(); iElement++) {
      int nGon = my_mesh->elementPtr(iElement)->nVertices();
      int local_num_cell_dofs = max( 0, (polynomial_degree - nGon + 3)
                            * (polynomial_degree - nGon + 2) / 2);
      cell_loc_to_glob_full[iElement] = new int[local_num_cell_dofs];
    
      for (int i = 0; i < local_num_cell_dofs; i++) {
        dof_num++;
        the_dof_type_full[dof_num] = MixedDoFType::cell;
        the_dof_bc_type_full[dof_num] = MixedDoFBCType::interior;
        the_dof_type_reduced[dof_num] = MixedDoFType::cell;
        the_dof_bc_type_reduced[dof_num] = MixedDoFBCType::interior;
        cell_loc_to_glob_full[iElement][i] = dof_num;
        loc_to_glob_full[iElement][i] = dof_num;
        loc_to_glob_reduced[iElement][i] = dof_num;
      }
    }
  }

  // Loop through elements to set up poly dofs
  int local_poly_dofs_full = (polynomial_degree + 2) * (polynomial_degree + 1) / 2 - 1;
  int local_poly_dofs_reduced = (polynomial_degree + 1) * polynomial_degree / 2 - 1;
  int dof_num_full = dof_num;
  int dof_num_reduced = dof_num;
  for(int iElement = 0; iElement < my_mesh->nElements(); iElement++) {
    poly_loc_to_glob_full[iElement] = new int[local_poly_dofs_full];
    poly_loc_to_glob_reduced[iElement] = new int[local_poly_dofs_reduced];

    int num_of_previous_dofs = my_mesh->elementPtr(iElement)->nVertices()*(polynomial_degree+1)
                              +max( 0, (polynomial_degree - my_mesh->elementPtr(iElement)->nVertices() + 3)
                              * (polynomial_degree - my_mesh->elementPtr(iElement)->nVertices() + 2) / 2);

    for (int i = 0; i < local_poly_dofs_reduced; i++) {
      dof_num_full++;
      dof_num_reduced++;
      the_dof_type_full[dof_num_full] = MixedDoFType::poly;
      the_dof_bc_type_full[dof_num_full] = MixedDoFBCType::interior;
      the_dof_type_reduced[dof_num_reduced] = MixedDoFType::poly;
      the_dof_bc_type_reduced[dof_num_reduced] = MixedDoFBCType::interior;
      poly_loc_to_glob_full[iElement][i] = dof_num_full;
      poly_loc_to_glob_reduced[iElement][i] = dof_num_reduced;
      loc_to_glob_full[iElement][num_of_previous_dofs+i]=dof_num_full;
      loc_to_glob_reduced[iElement][num_of_previous_dofs+i]=dof_num_reduced;
    }

    for (int i = local_poly_dofs_reduced; i < local_poly_dofs_full; i++) {
      dof_num_full++;
      the_dof_type_full[dof_num_full] = MixedDoFType::poly;
      the_dof_bc_type_full[dof_num_full] = MixedDoFBCType::interior;
      poly_loc_to_glob_full[iElement][i] = dof_num_full;
      loc_to_glob_full[iElement][num_of_previous_dofs+i]=dof_num_full;
    }
  }

  assert(dof_num_full == mixed_dofs_full - 1);
  assert(dof_num_reduced == mixed_dofs_reduced - 1);

}

DirectMixedConf::~DirectMixedConf() {
  if (the_dm_elements) delete[] the_dm_elements;
  if (the_dof_type_full) delete[] the_dof_type_full;
  if (the_dof_type_reduced) delete[] the_dof_type_reduced;
  if (the_dof_bc_type_full) delete[] the_dof_bc_type_full;
  if (the_dof_bc_type_reduced) delete[] the_dof_bc_type_reduced;
  if (vertex_loc_to_glob_full) delete[] vertex_loc_to_glob_full;
  if (vertex_loc_to_glob_coeff_full) delete[] vertex_loc_to_glob_coeff_full;
  if (mesh_vertex_to_dof_index) delete[] mesh_vertex_to_dof_index;
  if (edge_loc_to_glob_full) delete[] edge_loc_to_glob_full;
  if (cell_loc_to_glob_full) delete[] cell_loc_to_glob_full;
  if (poly_loc_to_glob_full) delete[] poly_loc_to_glob_full;
  if (poly_loc_to_glob_reduced) delete[] poly_loc_to_glob_reduced;
  if (loc_to_glob_full) delete[] loc_to_glob_full;
  if (loc_to_glob_reduced) delete[] loc_to_glob_reduced;
}


void DirectMixedConf::write_raw(std::ofstream& fout) const { 

  DirectMixed::write_raw(fout);

  fout << "DIRECT MIXED CONFORMING SPACE:\n";
  fout << "num_dofs of V_full     = " << mixed_dofs_full << "\n";
  fout << "num_dofs of V_reduced  = " << mixed_dofs_reduced << "\n";

  fout << "\nglobal dofs of V_full:\n";
  for(int i=0; i<mixed_dofs_full; i++) {
    fout << "  DoF " << i << " (type ";
    switch(the_dof_type_full[i]) {
    case MixedDoFType::vertex: { fout << "vertex"; break; }
    case MixedDoFType::edge:   { fout << "edge"; break; }
    case MixedDoFType::cell:   { fout << "cell"; break; }
    case MixedDoFType::poly:   { fout << "poly"; break; }
    default: { fout << "???"; break; }
    }
    fout << ", bc ";
    switch(the_dof_bc_type_full[i]) {
    case MixedDoFBCType::interior:  { fout << "interior"; break; }
    case MixedDoFBCType::boundary: { fout << "boundary"; break; }
    default: { fout << "???"; break; }
   }
    fout << ")\n";
  }
  
  fout << "\nglobal dofs of V_reduced:\n";
  for(int i=0; i<mixed_dofs_reduced; i++) {
    fout << "  DoF " << i << " (type ";
    switch(the_dof_type_reduced[i]) {
    case MixedDoFType::vertex: { fout << "vertex"; break; }
    case MixedDoFType::edge:   { fout << "edge"; break; }
    case MixedDoFType::cell:   { fout << "cell"; break; }
    case MixedDoFType::poly:   { fout << "poly"; break; }
    default: { fout << "???"; break; }
    }
    fout << ", bc ";
    switch(the_dof_bc_type_reduced[i]) {
    case MixedDoFBCType::interior:  { fout << "interior"; break; }
    case MixedDoFBCType::boundary: { fout << "boundary"; break; }
    default: { fout << "???"; break; }
   }
    fout << ")\n";
  }

  fout << "\ndmSpace connectivity\n";
  fout << "\nLocal edge indexing mapped to global dof:\n";
  for(int i=0; i<my_mesh->nElements(); i++) {
    fout << " Element " << i << " :\n";
    for (int iVertex = 0; iVertex < my_mesh -> elementPtr(i)->nVertices(); iVertex++) {
      fout << " Edge [" << iVertex << "] map to Global DoF [" << vertex_loc_to_glob(i,iVertex) << "]\n";
      fout << "\tOrientation Coefficient: "<< vertex_loc_to_glob_coeff(i, iVertex) << "\n";
    }
  }
  fout << "\n";

  fout << "\nMixed elements:\n";
  for(int i=0; i<my_mesh->nElements(); i++) {
    fout << " Element " << i << "\n";
    the_dm_elements[i].write_raw(fout);
  }
  fout << "\n";
};

int DirectMixedConf::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// class DirectMixedArray

void DirectMixedArray::set_directmixedarray(DirectMixed* dmSpace, char spacetype) {
  my_dm_space = dmSpace;
  num_dofs = my_dm_space -> nMixedDoFs(spacetype);
  num_elements = my_dm_space -> my_mesh -> nElements();
  my_mesh = my_dm_space -> my_mesh;

  space_type = spacetype;
  // We must make sure that space type is either 'f' or 'r'
  assert(space_type == 'f' || space_type == 'r');

  if (the_array) delete[] the_array;
  the_array = new double[num_dofs];
}

DirectMixedArray::~DirectMixedArray() {
  if(the_array) delete[] the_array;
}

void DirectMixedArray::l2normError(double& l2Error, double& l2Norm, int refinement_level, Tensor1 (*referenceFcn)(double,double)) {
  l2Error = 0, l2Norm = 0;
  PolyQuadrature quadRule(13,refinement_level);

  for(int iElement=0; iElement < num_elements; iElement++) {
    polymesh::PolyElement* mixedElementPtr = my_dm_space->elementPtr(iElement);
    quadRule.setElement(mixedElementPtr);
    
    for(int iPt=0; iPt<quadRule.num(); iPt++) {
      double x = quadRule.pt(iPt).val(0);
      double y = quadRule.pt(iPt).val(1);
      
      Tensor1 result;
      eval(quadRule.pt(iPt), result);

      Tensor1 diff = (referenceFcn == nullptr) ? result : (result - referenceFcn(x,y));
      
      l2Error += diff * diff * quadRule.wt(iPt);
      l2Norm += result * result * quadRule.wt(iPt);
    }
  }
  
  l2Error = sqrt(l2Error);
  l2Norm = sqrt(l2Norm);
}

void DirectMixedArray::l2normError_div(double& l2Error, double& l2Norm, int refinement_level, double (*referenceFcn)(double,double)) {
  l2Error = 0, l2Norm = 0;
  PolyQuadrature quadRule(13,refinement_level);

  for(int iElement=0; iElement < num_elements; iElement++) {
    polymesh::PolyElement* mixedElementPtr = my_dm_space->elementPtr(iElement);
    quadRule.setElement(mixedElementPtr);
    
    for(int iPt=0; iPt<quadRule.num(); iPt++) {
      double x = quadRule.pt(iPt).val(0);
      double y = quadRule.pt(iPt).val(1);
      
      double result;
      eval_div(quadRule.pt(iPt), result);

      double diff = (referenceFcn == nullptr) ? result : (result - referenceFcn(x,y));
      
      l2Error += diff * diff * quadRule.wt(iPt);
      l2Norm += result * result * quadRule.wt(iPt);
    }
  }
  
  l2Error = sqrt(l2Error);
  l2Norm = sqrt(l2Norm);
}


void DirectMixedArray::write_matlab_mesh(std::ofstream* fout, int num_pts_x, int num_pts_y) const {

  if(num_pts_x <= 1) num_pts_x = 2;
  if(num_pts_y <= 1) num_pts_y = 2;

  // Determine mesh of points
  double xMin = my_dm_space->my_mesh->minX();
  double xMax = my_dm_space->my_mesh->maxX();
  double yMin = my_dm_space->my_mesh->minY();
  double yMax = my_dm_space->my_mesh->maxY();

  double dx = (xMax - xMin)/(num_pts_x-1);
  double dy = (yMax - yMin)/(num_pts_y-1);

  Point* pts = new Point[num_pts_x*num_pts_y];

  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      pts[j + num_pts_y*i].set(xMin+i*dx, yMin+j*dy);
    }
  }

  // Evaluate
  Tensor1* result = new Tensor1[num_pts_x*num_pts_y];

  eval(pts, result, num_pts_x*num_pts_y);

  // Write file  
  *fout << "quiver(" << xMin << ":" << dx << ":" << xMax << ","
	<< yMin << ":" << dy << ":" << yMax <<",[ ";

  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      *fout << result[i + num_pts_x*j][0] << " ";
    }
    *fout << "; ";
    *fout << "; ";
  }
  *fout << "],[ \n";

  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      *fout << result[i + num_pts_x*j][1] << " ";
    }
    *fout << "; ";
  }
  *fout << "]);\n";
  *fout << "xlabel('x'); ylabel('y');\n";


  delete[] result;
  delete[] pts;
};

int DirectMixedArray::write_matlab_mesh(std::string& filename, int num_pts_x, int num_pts_y) const {
  std::ofstream fout(filename+".m");
  if( !fout ) return 1;
  write_matlab_mesh(&fout, num_pts_x, num_pts_y);
  return 0;
}

void DirectMixedArray::write_matlab_mesh_by_pt(std::ofstream& fout, int num_pts_x, int num_pts_y) const {
  if(num_pts_x <= 1) num_pts_x = 2;
  if(num_pts_y <= 1) num_pts_y = 2;
  
  double xMin = my_dm_space->my_mesh->minX();
  double xMax = my_dm_space->my_mesh->maxX();
  double yMin = my_dm_space->my_mesh->minY();
  double yMax = my_dm_space->my_mesh->maxY();

  double dx = (xMax - xMin)/(num_pts_x-1);
  double dy = (yMax - yMin)/(num_pts_y-1);

  Tensor1 result;
  double tensorY[num_pts_x][num_pts_y];
  
  fout << "quiver(" << xMin << ":" << dx << ":" << xMax << ","
	    << yMin << ":" << dy << ":" << yMax <<",[ ";

  for(int j=0; j<num_pts_y; j++) {
    for(int i=0; i<num_pts_x; i++) {
      eval(Point(xMin+i*dx,yMin+j*dy), result);
      fout << result[0] << " ";
      tensorY[i][j] = result[1];
    }
    fout << "; ";
  }
  fout << "],[ \n";

  for(int j=0; j<num_pts_y; j++) {
    for(int i=0; i<num_pts_x; i++) {
      fout << tensorY[i][j] << " ";
    }
    fout << "; ";
  }
  fout << "]);\n";
  fout << "xlabel('x'); ylabel('y');\n";
};


int DirectMixedArray::write_matlab_mesh_by_pt(std::string& filename, int num_pts_x, int num_pts_y) const {
  std::ofstream fout(filename+".m");
  if( !fout ) return 1;
  write_matlab_mesh_by_pt(fout, num_pts_x, num_pts_y);
  return 0;
}

void DirectMixedArray::write_raw(std::ofstream& fout) const {
  for(int i=0; i<num_dofs; i++) {
    fout << the_array[i] << "  ";
    if( !(i%10) ) fout << "\n";
  }
};

int DirectMixedArray::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}

void DirectMixedArray::write_matlab_mesh_div(std::ofstream* fout, int num_pts_x, int num_pts_y) const {
  if(num_pts_x <= 1) num_pts_x = 2;
  if(num_pts_y <= 1) num_pts_y = 2;

  // Determine mesh of points
  double xMin = my_dm_space->my_mesh->minX();
  double xMax = my_dm_space->my_mesh->maxX();
  double yMin = my_dm_space->my_mesh->minY();
  double yMax = my_dm_space->my_mesh->maxY();

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
  eval_div(pts, result, num_pts_x*num_pts_y);

  // Write file  
  *fout << "mesh(" << xMin << ":" << dx << ":" << xMax << ","
	<< yMin << ":" << dy << ":" << yMax <<",[ ";

  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      *fout << result[i + num_pts_x*j] << " ";
    }
    *fout << "; ";
  }
  *fout << "]);\n";
  *fout << "xlabel('x'); ylabel('y');\n";

  delete[] pts;
}

int DirectMixedArray::write_matlab_mesh_div(std::string& filename, int num_pts_x, int num_pts_y) const {
  std::ofstream fout(filename+".m");
  if( !fout ) return 1;
  write_matlab_mesh_div(&fout, num_pts_x, num_pts_y);
  return 0;
}

void DirectMixedArray::write_matlab_mesh_error(std::ofstream* fout, int num_pts_x, int num_pts_y, 
                                                Tensor1 (*referenceFcn)(double,double)) const {
  if(num_pts_x <= 1) num_pts_x = 2;
  if(num_pts_y <= 1) num_pts_y = 2;

  // Determine mesh of points
  double xMin = my_dm_space->my_mesh->minX();
  double xMax = my_dm_space->my_mesh->maxX();
  double yMin = my_dm_space->my_mesh->minY();
  double yMax = my_dm_space->my_mesh->maxY();

  double dx = (xMax - xMin)/(num_pts_x-1);
  double dy = (yMax - yMin)/(num_pts_y-1);

  Point* pts = new Point[num_pts_x*num_pts_y];
  
  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      pts[j + num_pts_y*i].set(xMin+i*dx, yMin+j*dy);
    }
  }

  // Evaluate
  Tensor1* result = new Tensor1[num_pts_x*num_pts_y];
  eval(pts, result, num_pts_x*num_pts_y);

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
      *fout << result[i + num_pts_x*j].norm() << " ";
    }
    *fout << "; ";
  }
  *fout << "]);\n";
  *fout << "xlabel('x'); ylabel('y');\n";


  delete[] result;
  delete[] pts;  
}

int DirectMixedArray::write_matlab_mesh_error(std::string& filename, int num_pts_x, int num_pts_y, Tensor1 (*referenceFcn)(double,double)) const {
  std::ofstream fout(filename+".m");
  if( !fout ) return 1;
  write_matlab_mesh_error(&fout, num_pts_x, num_pts_y, referenceFcn);
  return 0;
}

void DirectMixedArray::write_matlab_mesh_div_error(std::ofstream* fout, int num_pts_x, int num_pts_y, double (*referenceFcn)(double,double)) const {
  if(num_pts_x <= 1) num_pts_x = 2;
  if(num_pts_y <= 1) num_pts_y = 2;

  // Determine mesh of points
  double xMin = my_dm_space->my_mesh->minX();
  double xMax = my_dm_space->my_mesh->maxX();
  double yMin = my_dm_space->my_mesh->minY();
  double yMax = my_dm_space->my_mesh->maxY();

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
  eval_div(pts, result, num_pts_x*num_pts_y);

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

  delete[] pts;
}

int DirectMixedArray::write_matlab_mesh_div_error(std::string& filename, int num_pts_x, int num_pts_y, double (*referenceFcn)(double,double)) const {
  std::ofstream fout(filename+".m");
  if( !fout ) return 1;
  write_matlab_mesh_div_error(&fout, num_pts_x, num_pts_y, referenceFcn);
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// class DirectMixedHybridArray
void DirectMixedHybridArray::set_directmixedhybridarray(DirectMixedHybrid* dmSpace, char spacetype) {
  set_directmixedarray(dmSpace, spacetype);
  my_dm_space = dmSpace;
};

void DirectMixedHybridArray::eval(const Point* pts, 
	                          Tensor1* result, int num_pts) const {

  bool ptEvaluated[num_pts];
  for(int i=0; i<num_pts; i++) ptEvaluated[i] = false;

  // Loop through the elements
  std::vector<Point> elementPts;
  std::vector<int> elementPtsIndex;
  int global_index = 0;
  for(int iElement=0; iElement < num_elements; iElement++) {
    DirectMixedHybridFE* mixedElement = my_dm_space->MixedElementPtr(iElement);
    PolyElement* element = mixedElement->elementPtr();
    

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
    int nDoFs = ( space_type == 'f' ) ? mixedElement -> dimVFull() : mixedElement -> dimVReduced();

    double dofs[nDoFs];
    for(int i=0; i<nDoFs; i++) {
      dofs[i] = the_array[global_index];
      global_index ++;
    }
    // Evaluate array at points on element
    Tensor1* elementResult = new Tensor1[elementPts.size()];
    

    mixedElement->eval(elementPts.data(), elementResult, elementPts.size(), space_type, dofs);

    // Place results in global array
    for(unsigned long int i=0; i<elementPts.size(); i++) {
      result[elementPtsIndex[i]] = elementResult[i];
    }

    delete[] elementResult;
  }

  // Attend to unset points (outside the mesh)
  for(int i=0; i<num_pts; i++) {
    if(!ptEvaluated[i]) {
      result[i].set(0,0);
    }
  }  
}

void DirectMixedHybridArray::eval(const Point& pt, Tensor1& result) const {
  int iElement = my_dm_space->my_mesh->elementIndex(pt);
  if(iElement < 0) {
    result.set(0,0);
    return;
  }
  DirectMixedHybridFE* elem = my_dm_space->MixedElementPtr(iElement);

  // SET DoFs for element
  int nDoFs = ( space_type == 'f' ) ? elem -> dimVFull() : elem -> dimVReduced();

  // Find corresponding global index of first local dof
  int global_index = my_dm_space->mixed_Elem_First_To_Global_Dof(iElement, space_type);

  double dofs[nDoFs];
  for(int i=0; i<nDoFs; i++) {
    dofs[i] = the_array[global_index];
    global_index ++;
  }

  // Evaluate
  elem->eval(pt, result, space_type, dofs);
};

void DirectMixedHybridArray::eval_div(const Point* pts, 
	                          double* result, int num_pts) const {
  bool ptEvaluated[num_pts];
  for(int i=0; i<num_pts; i++) ptEvaluated[i] = false;

  // Loop through the elements
  std::vector<Point> elementPts;
  std::vector<int> elementPtsIndex;
  int global_index = 0;
  for(int iElement=0; iElement < num_elements; iElement++) {
    DirectMixedHybridFE* mixedElement = my_dm_space->MixedElementPtr(iElement);
    PolyElement* element = mixedElement->elementPtr();
    

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
    int nDoFs = ( space_type == 'f' ) ? mixedElement -> dimVFull() : mixedElement -> dimVReduced();

    double dofs[nDoFs];
    for(int i=0; i<nDoFs; i++) {
      dofs[i] = the_array[global_index];
      global_index ++;
    }

    // Evaluate array at points on element
    double* elementResult = new double[elementPts.size()];
    mixedElement->eval_div(elementPts.data(), elementResult, elementPts.size(), space_type, dofs);
    // Place results in global array
    for(unsigned long int i=0; i<elementPts.size(); i++) {
      result[elementPtsIndex[i]] = elementResult[i];
    }

    delete[] elementResult;
  }

  // Attend to unset points (outside the mesh)
  for(int i=0; i<num_pts; i++) {
    if(!ptEvaluated[i]) {
      result[i] = 0;
    }
  }  
}

void DirectMixedHybridArray::eval_div(const Point& pt, double& result) const {
  int iElement = my_dm_space->my_mesh->elementIndex(pt);
  if(iElement < 0) {
    result = 0;
    return;
  }
  DirectMixedHybridFE* elem = my_dm_space->MixedElementPtr(iElement);

  // SET DoFs for element
  int nDoFs = ( space_type == 'f' ) ? elem -> dimVFull() : elem -> dimVReduced();

  // Find corresponding global index of first local dof
  int global_index = my_dm_space->mixed_Elem_First_To_Global_Dof(iElement, space_type);

  double dofs[nDoFs];
  for(int i=0; i<nDoFs; i++) {
    dofs[i] = the_array[global_index];
    global_index ++;
  }

  // Evaluate
  elem->DirectMixedFE::eval_div(pt, result, space_type, dofs);
};

////////////////////////////////////////////////////////////////////////////////
// class DirectMixedConfArray

void DirectMixedConfArray::set_directmixedconfarray(DirectMixedConf* dmSpace, char spacetype) {
  set_directmixedarray(dmSpace, spacetype);
  my_dm_space = dmSpace;
};

void DirectMixedConfArray::eval(const Point* pts, Tensor1* result, int num_pts) const {

  bool ptEvaluated[num_pts];
  for(int i=0; i<num_pts; i++) ptEvaluated[i] = false;

  // Loop through the elements
  std::vector<Point> elementPts;
  std::vector<int> elementPtsIndex;

  for(int iElement=0; iElement < my_dm_space->my_mesh->nElements(); iElement++) {
    DirectMixedConfFE* mixedElement = my_dm_space->MixedElementPtr(iElement);
    PolyElement* element = mixedElement->elementPtr();

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
    int nGon = element -> nVertices();
    int nDoFs = ( space_type == 'f' ) ? mixedElement -> dimVFull() : mixedElement -> dimVReduced();

    double dofs[nDoFs];
    int local_index = 0;

 for (int i=0; i<mixedElement -> dimCellBasis(); i++) {
    int iCellDoF = my_dm_space -> cell_loc_to_glob(iElement,i);
    dofs[local_index] = the_array[iCellDoF];
    local_index++;
    
  }

  
  for(int i=0; i<nGon; i++) {
    for(int j=0; j<(my_dm_space->degPolyn()); j++) {
      int iEdgeDoF = my_dm_space -> edge_loc_to_glob(iElement,i,j);
      dofs[local_index] = the_array[iEdgeDoF];
      local_index++;
    }
  }
  for(int i=0; i<nGon; i++) {
    int iVertexOrientation = my_dm_space -> vertex_loc_to_glob_coeff(iElement,i);
    int iVertexDof =  my_dm_space -> vertex_loc_to_glob(iElement,i);
    dofs[local_index] =  iVertexOrientation * the_array[iVertexDof];
    local_index++;
  }

  int polyDoFs = mixedElement -> dimPolyBasis(space_type);
  for (int i=0; i<polyDoFs; i++) {
    int iPolyDoF = my_dm_space -> poly_loc_to_glob(iElement, i, space_type);
    dofs[local_index] = the_array[iPolyDoF];
    local_index++;
  }
    assert(local_index == nDoFs);

    // Evaluate array at points on element
    Tensor1* elementResult = new Tensor1[elementPts.size()];
    mixedElement->eval(elementPts.data(), elementResult, elementPts.size(), space_type, dofs);
    // Place results in global array
    for(unsigned long int i=0; i<elementPts.size(); i++) {
      result[elementPtsIndex[i]] = elementResult[i];
    }

    delete[] elementResult;
  }

  // Attend to unset points (outside the mesh)
  for(int i=0; i<num_pts; i++) {
    if(!ptEvaluated[i]) {
      result[i].set(0,0);
    }
  }  
}

void DirectMixedConfArray::eval(const Point& pt, Tensor1& result) const {
  int iElement = my_dm_space->my_mesh->elementIndex(pt);
  if(iElement < 0) {
    result.set(0,0);
    return;
  }
  DirectMixedConfFE* mixedElement = my_dm_space->MixedElementPtr(iElement);



  // SET DoFs for element
  int nGon = mixedElement->elementPtr() -> nVertices();
  int nDoFs = ( space_type == 'f' ) ? mixedElement -> dimVFull() : mixedElement -> dimVReduced();

  double dofs[nDoFs];
  int local_index = 0;

  for (int i=0; i<mixedElement -> dimCellBasis(); i++) {
    int iCellDoF = my_dm_space -> cell_loc_to_glob(iElement,i);
    dofs[local_index] = the_array[iCellDoF];
    local_index++;
    
  }

  
  for(int i=0; i<nGon; i++) {
    for(int j=0; j<(my_dm_space->degPolyn()); j++) {
      int iEdgeDoF = my_dm_space -> edge_loc_to_glob(iElement,i,j);
      dofs[local_index] = the_array[iEdgeDoF];
      local_index++;
    }
  }
  for(int i=0; i<nGon; i++) {
    int iVertexOrientation = my_dm_space -> vertex_loc_to_glob_coeff(iElement,i);
    int iVertexDof =  my_dm_space -> vertex_loc_to_glob(iElement,i);
    dofs[local_index] =  iVertexOrientation * the_array[iVertexDof];
    local_index++;
  }

  int polyDoFs = mixedElement -> dimPolyBasis(space_type);
  for (int i=0; i<polyDoFs; i++) {
    int iPolyDoF = my_dm_space -> poly_loc_to_glob(iElement, i, space_type);
    dofs[local_index] = the_array[iPolyDoF];
    local_index++;
  }

  assert(local_index == nDoFs);

  // Evaluate
  mixedElement->eval(pt, result, space_type, dofs);
}

void DirectMixedConfArray::eval_div(const Point* pts, double* result, int num_pts) const {
  bool ptEvaluated[num_pts];
  for(int i=0; i<num_pts; i++) ptEvaluated[i] = false;

  // Loop through the elements
  std::vector<Point> elementPts;
  std::vector<int> elementPtsIndex;

  for(int iElement=0; iElement < my_dm_space->my_mesh->nElements(); iElement++) {
    DirectMixedConfFE* mixedElement = my_dm_space->MixedElementPtr(iElement);
    PolyElement* element = mixedElement->elementPtr();

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
    int nGon = element -> nVertices();
    int nDoFs = ( space_type == 'f' ) ? mixedElement -> dimVFull() : mixedElement -> dimVReduced();

    double dofs[nDoFs];
    int local_index = 0;    
    
    for (int i=0; i<mixedElement -> dimCellBasis(); i++) {
      int iCellDoF = my_dm_space -> cell_loc_to_glob(iElement,i);
      dofs[local_index] = the_array[iCellDoF];
      local_index++;
    }

  for(int i=0; i<nGon; i++) {
    for(int j=0; j<(my_dm_space->degPolyn()); j++) {
      int iEdgeDoF = my_dm_space -> edge_loc_to_glob(iElement,i,j);
      dofs[local_index] = the_array[iEdgeDoF];
      local_index++;
    }
  }
    for(int i=0; i<nGon; i++) {
      int iVertexOrientation = my_dm_space -> vertex_loc_to_glob_coeff(iElement,i);
      int iVertexDof =  my_dm_space -> vertex_loc_to_glob(iElement,i);
      dofs[local_index] =  iVertexOrientation * the_array[iVertexDof];
      local_index++;
    }

    int polyDoFs = mixedElement -> dimPolyBasis(space_type);
    for (int i=0; i<polyDoFs; i++) {
      int iPolyDoF = my_dm_space -> poly_loc_to_glob(iElement, i, space_type);
      dofs[local_index] = the_array[iPolyDoF];
      local_index++;
    }

    assert(local_index == nDoFs);

    // Evaluate array at points on element
    double* elementResult = new double[elementPts.size()];
    mixedElement->eval_div(elementPts.data(), elementResult, elementPts.size(), space_type, dofs);
    // Place results in global array
    for(unsigned long int i=0; i<elementPts.size(); i++) {
      result[elementPtsIndex[i]] = elementResult[i];
    }

    delete[] elementResult;
  }

  // Attend to unset points (outside the mesh)
  for(int i=0; i<num_pts; i++) {
    if(!ptEvaluated[i]) {
      result[i] = 0;
    }
  }  
}

void DirectMixedConfArray::eval_div(const Point& pt, double& result) const {
  int iElement = my_dm_space->my_mesh->elementIndex(pt);
  if(iElement < 0) {
    result = 0;
  }
  DirectMixedConfFE* mixedElement = my_dm_space->MixedElementPtr(iElement);

  // SET DoFs for element
  int nGon = mixedElement->elementPtr() -> nVertices();
  int nDoFs = ( space_type == 'f' ) ? mixedElement -> dimVFull() : mixedElement -> dimVReduced();

  double dofs[nDoFs];
  int local_index = 0;

  for (int i=0; i<mixedElement -> dimCellBasis(); i++) {
    int iCellDoF = my_dm_space -> cell_loc_to_glob(iElement,i);
    dofs[local_index] = the_array[iCellDoF];
    local_index++;
  }

  for(int i=0; i<nGon; i++) {
    for(int j=0; j<(my_dm_space->degPolyn()); j++) {
      int iEdgeDoF = my_dm_space -> edge_loc_to_glob(iElement,i,j);
      dofs[local_index] = the_array[iEdgeDoF];
      local_index++;
    }
  }

  for(int i=0; i<nGon; i++) {
    int iVertexOrientation = my_dm_space -> vertex_loc_to_glob_coeff(iElement,i);
    int iVertexDof =  my_dm_space -> vertex_loc_to_glob(iElement,i);
    dofs[local_index] =  iVertexOrientation * the_array[iVertexDof];
    local_index++;
  }
  
  int polyDoFs = mixedElement -> dimPolyBasis(space_type);
  for (int i=0; i<polyDoFs; i++) {
    int iPolyDoF = my_dm_space -> poly_loc_to_glob(iElement, i, space_type);
    dofs[local_index] = the_array[iPolyDoF];
    local_index++;
  }

  assert(local_index == nDoFs);

  // Evaluate
  mixedElement->DirectMixedFE::eval_div(pt, result, space_type, dofs);
}

////////////////////////////////////////////////////////////////////////////////
// class DirectDGArray

void DirectDGArray::set_directdgarray(DirectMixed* dmSpace, char spacetype) {
  my_dm_space = dmSpace;
  num_elements = my_dm_space->my_mesh->nElements();
  space_type = spacetype;
  num_dofs = my_dm_space->nDGDoFs(space_type);

  // We must make sure that space type is either 'f' or 'r'
  assert(space_type == 'f' || space_type == 'r');

  if(the_array) delete[] the_array;
  the_array = new double[num_dofs];
}

DirectDGArray::~DirectDGArray() {
  if(the_array) delete[] the_array;
}

void DirectDGArray::eval(const Point* pts, double* result, int num_pts) const {
  bool ptEvaluated[num_pts];
  for(int i=0; i<num_pts; i++) ptEvaluated[i] = false;

  // Loop through the elements
  std::vector<Point> elementPts;
  std::vector<int> elementPtsIndex;
  int global_index = 0;
  for(int iElement=0; iElement < my_dm_space->my_mesh->nElements(); iElement++) {
    DirectDGFE* dgElement = my_dm_space->DGElementPtr(iElement);
    PolyElement* element = dgElement->elementPtr();
    

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
    int nDoFs = ( space_type == 'f' ) ? dgElement -> dimFull() : dgElement -> dimReduced();

    double dofs[nDoFs];
    for(int i=0; i<nDoFs; i++) {
      dofs[i] = the_array[global_index];
      global_index ++;
    }

    // Evaluate array at points on element
    double elementResult[elementPts.size()];
    dgElement->eval(elementPts.data(), elementResult, elementPts.size(), space_type, dofs);

    // Place results in global array
    for(unsigned long int i=0; i<elementPts.size(); i++) {
      result[elementPtsIndex[i]] = elementResult[i];
    }
  }

  // Attend to unset points (outside the mesh)
  for(int i=0; i<num_pts; i++) {
    if(!ptEvaluated[i]) {
      result[i] = 0;
    }
  }  
}

void DirectDGArray::eval(const Point& pt, double& result) const {
  int iElement = my_dm_space->my_mesh->elementIndex(pt);
  if(iElement < 0) {
    result = 0;
    return;
  }
  DirectDGFE* elem = my_dm_space->DGElementPtr(iElement);

  // SET DoFs for element
  int nDoFs = ( space_type == 'f' ) ? elem -> dimFull() : elem -> dimReduced();

  // Find corresponding global index of first local dof
  int global_index = my_dm_space->dg_Elem_First_To_Global_Dof(iElement, space_type);

  double dofs[nDoFs];
  for(int i=0; i<nDoFs; i++) {
    dofs[i] = the_array[global_index];
    global_index ++;
  }

  // Evaluate
  elem->eval(pt, result, space_type , dofs);
};



void DirectDGArray::l2normError(double& l2Error, double& l2Norm, int refinement_level, 
					 double (*referenceFcn)(double,double)) {
  l2Error = 0, l2Norm = 0;
  PolyQuadrature quadRule(13,refinement_level);

  for(int iElement=0; iElement < my_dm_space->mesh()->nElements(); iElement++) {
    DirectDGFE* dgPtr = my_dm_space->DGElementPtr(iElement);
    quadRule.setElement(dgPtr->elementPtr());
    
    for(int iPt=0; iPt<quadRule.num(); iPt++) {
      double x = quadRule.pt(iPt).val(0);
      double y = quadRule.pt(iPt).val(1);
      
      double result; 
      eval(quadRule.pt(iPt), result);
      
      double diff = (referenceFcn == nullptr) ? result : (result - referenceFcn(x,y));
      
      l2Error += pow(diff,2) * quadRule.wt(iPt);
      l2Norm += pow(result,2) * quadRule.wt(iPt);
    }
  }
  l2Error = sqrt(l2Error);
  l2Norm = sqrt(l2Norm);
}

void DirectDGArray::write_matlab_mesh(std::ofstream* fout, int num_pts_x, int num_pts_y) const {
  if(num_pts_x <= 1) num_pts_x = 2;
  if(num_pts_y <= 1) num_pts_y = 2;

  // Determine mesh of points
  double xMin = my_dm_space->my_mesh->minX();
  double xMax = my_dm_space->my_mesh->maxX();
  double yMin = my_dm_space->my_mesh->minY();
  double yMax = my_dm_space->my_mesh->maxY();

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
  eval(pts, result, num_pts_x*num_pts_y);

  // Write file  
  *fout << "mesh(" << xMin << ":" << dx << ":" << xMax << ","
	<< yMin << ":" << dy << ":" << yMax <<",[ ";

  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      *fout << result[i + num_pts_x*j] << " ";
    }
    *fout << "; ";
  }
  *fout << "]);\n";
  *fout << "xlabel('x'); ylabel('y');\n";

  delete[] pts;
};

int DirectDGArray::write_matlab_mesh(std::string& filename, int num_pts_x, int num_pts_y) const {
  std::ofstream fout(filename+".m");
  if( !fout ) return 1;
  write_matlab_mesh(&fout, num_pts_x, num_pts_y);
  return 0;
}


void DirectDGArray::write_matlab_mesh_by_pt(std::ofstream& fout, int num_pts_x, int num_pts_y) const {
  if(num_pts_x <= 1) num_pts_x = 2;
  if(num_pts_y <= 1) num_pts_y = 2;
  
  double xMin = my_dm_space->my_mesh->minX();
  double xMax = my_dm_space->my_mesh->maxX();
  double yMin = my_dm_space->my_mesh->minY();
  double yMax = my_dm_space->my_mesh->maxY();

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


int DirectDGArray::write_matlab_mesh_by_pt(std::string& filename, int num_pts_x, int num_pts_y) const {
  std::ofstream fout(filename+".m");
  if( !fout ) return 1;
  write_matlab_mesh_by_pt(fout, num_pts_x, num_pts_y);
  return 0;
}

void DirectDGArray::write_raw(std::ofstream& fout) const {
  for(int i=0; i<num_dofs; i++) {
    fout << the_array[i] << "  ";
    if( !(i%10) ) fout << "\n";
  }
};
int DirectDGArray::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}

void DirectDGArray::write_matlab_mesh_error(std::ofstream* fout, int num_pts_x, int num_pts_y, double (*referenceFcn)(double,double)) const {
  if(num_pts_x <= 1) num_pts_x = 2;
  if(num_pts_y <= 1) num_pts_y = 2;

  // Determine mesh of points
  double xMin = my_dm_space->my_mesh->minX();
  double xMax = my_dm_space->my_mesh->maxX();
  double yMin = my_dm_space->my_mesh->minY();
  double yMax = my_dm_space->my_mesh->maxY();

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
  eval(pts, result, num_pts_x*num_pts_y);

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

  delete[] pts;  
}

int DirectDGArray::write_matlab_mesh_error(std::string& filename, int num_pts_x, int num_pts_y, double (*referenceFcn)(double,double)) const {
  std::ofstream fout(filename+".m");
  if( !fout ) return 1;
  write_matlab_mesh_error(&fout, num_pts_x, num_pts_y, referenceFcn);
  return 0;
}
////////////////////////////////////////////////////////////////////////////////
// class DirectEdgeDGArray

void DirectEdgeDGArray::set_directedgedgarray(DirectMixedHybrid* dmSpace) {
  my_dm_space = dmSpace;
  num_edges = my_dm_space->my_mesh->nEdges();
  num_dofs = my_dm_space->nEdgeDGDoFs();

  if(the_array) delete[] the_array;
  the_array = new double[num_dofs];
}

DirectEdgeDGArray::~DirectEdgeDGArray() {
  if(the_array) delete[] the_array;
}

void DirectEdgeDGArray::eval(const Point* pts, double* result, int num_pts) const {
  bool ptEvaluated[num_pts];
  for(int i=0; i<num_pts; i++) ptEvaluated[i] = false;

  // Loop through the edges
  std::vector<Point> edgePts;
  std::vector<int> edgePtsIndex;
  int global_index = 0;
  for(int iEdge=0; iEdge < my_dm_space->my_mesh->nEdges(); iEdge++) {
    DirectEdgeDGFE* edgElement = my_dm_space->DGEdgePtr(iEdge);
    Edge* edge = edgElement->edgePtr();
    
    // Set list of points in element
    edgePts.clear();
    edgePtsIndex.clear();
    for(int i=0; i<num_pts; i++) {
      if(ptEvaluated[i]) continue;
      

      if(edge->isInEdge(pts[i])) {
        edgePts.push_back(pts[i]);
        edgePtsIndex.push_back(i);
        ptEvaluated[i] = true;
      }
    }
    if(edgePts.size() == 0) continue;
    
     // SET DoFs for edge
    int nDoFs = edgElement -> dim();

    double dofs[nDoFs];
    for(int i=0; i<nDoFs; i++) {
      dofs[i] = the_array[global_index];
      global_index++;
    }

    // Evaluate array at points on element
    double edgeResult[edgePts.size()];
    edgElement->eval(edgePts.data(), edgeResult, edgePts.size(), dofs);

    // Place results in global array
    for(unsigned long int i=0; i<edgePts.size(); i++) {
      result[edgePtsIndex[i]] = edgeResult[i];
    }
  }

  // Attend to unset points (outside the mesh)
  for(int i=0; i<num_pts; i++) {
    if(!ptEvaluated[i]) {
      result[i] = 0;
    }
  }  
}

void DirectEdgeDGArray::eval(const Point& pt, double& result) const {
  int iElement = my_dm_space->my_mesh->elementIndex(pt);
  if(iElement < 0) {
    result = 0;
    return;
  }

  // We test if the point is on edge
  polymesh::PolyElement* elemPtr = my_dm_space->my_mesh->elementPtr(iElement);
  if (!elemPtr->isOnElementBoundary(pt)) {
    result = 0;
    return;
  }

  int edge_index, iEdge = -1;
  for (int nEdge = 0; nEdge < elemPtr->nVertices(); nEdge++) {
    edge_index = elemPtr->edgePtr(nEdge)->meshIndex(); //corresponding global index
    if (my_dm_space->my_mesh->edgePtr(edge_index)->isInEdge(pt)) {
      // point is on edge with the global index iEdge
      iEdge = edge_index;
    }
  }
  DirectEdgeDGFE* edge = my_dm_space->DGEdgePtr(iEdge);

  // SET DoFs for element
  int nDoFs = edge -> dim();

  // Find corresponding global index of first local dof
  int global_index = my_dm_space->edge_Elem_First_To_Global_Dof(iElement);
  
  double dofs[nDoFs];
  for(int i=0; i<nDoFs; i++) {
    dofs[i] = the_array[global_index];
    global_index ++;
  }

  // Evaluate
  edge->eval(pt, result, dofs);
};



void DirectEdgeDGArray::l2normError(double& l2Error, double& l2Norm, 
					 double (*referenceFcn)(double,double)) {
  l2Error = 0, l2Norm = 0;
  PolyEdgeQuadrature quadRule(13);

  for(int iEdge=0; iEdge < my_dm_space->mesh()->nEdges(); iEdge++) {
    DirectEdgeDGFE* edgPtr = my_dm_space->DGEdgePtr(iEdge);
    quadRule.setEdge(edgPtr->edgePtr());
    
    for(int iPt=0; iPt<quadRule.num(); iPt++) {
      double x = quadRule.pt(iPt).val(0);
      double y = quadRule.pt(iPt).val(1);
      
      double result; 
      eval(quadRule.pt(iPt), result);
      
      double diff = (referenceFcn == nullptr) ? result : (result - referenceFcn(x,y));
      
      l2Error += pow(diff,2) * quadRule.wt(iPt);
      l2Norm += pow(result,2) * quadRule.wt(iPt);
    }
  }
  l2Error = sqrt(l2Error);
  l2Norm = sqrt(l2Norm);
}

void DirectEdgeDGArray::write_matlab_mesh(std::ofstream* fout, int num_pts_x, int num_pts_y) const {
  if(num_pts_x <= 1) num_pts_x = 2;
  if(num_pts_y <= 1) num_pts_y = 2;

  // Determine mesh of points
  double xMin = my_dm_space->my_mesh->minX();
  double xMax = my_dm_space->my_mesh->maxX();
  double yMin = my_dm_space->my_mesh->minY();
  double yMax = my_dm_space->my_mesh->maxY();

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
  eval(pts, result, num_pts_x*num_pts_y);

  // Write file  
  *fout << "mesh(" << xMin << ":" << dx << ":" << xMax << ","
	<< yMin << ":" << dy << ":" << yMax <<",[ ";

  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      *fout << result[i + num_pts_x*j] << " ";
    }
    *fout << "; ";
  }
  *fout << "]);\n";
  *fout << "xlabel('x'); ylabel('y');\n";

  if (pts) delete[] pts;
};

int DirectEdgeDGArray::write_matlab_mesh(std::string& filename, int num_pts_x, int num_pts_y) const {
  std::ofstream fout(filename+".m");
  if( !fout ) return 1;
  write_matlab_mesh(&fout, num_pts_x, num_pts_y);
  return 0;
}


void DirectEdgeDGArray::write_matlab_mesh_by_pt(std::ofstream& fout, int num_pts_x, int num_pts_y) const {
  if(num_pts_x <= 1) num_pts_x = 2;
  if(num_pts_y <= 1) num_pts_y = 2;
  
  double xMin = my_dm_space->my_mesh->minX();
  double xMax = my_dm_space->my_mesh->maxX();
  double yMin = my_dm_space->my_mesh->minY();
  double yMax = my_dm_space->my_mesh->maxY();

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


int DirectEdgeDGArray::write_matlab_mesh_by_pt(std::string& filename, int num_pts_x, int num_pts_y) const {
  std::ofstream fout(filename+".m");
  if( !fout ) return 1;
  write_matlab_mesh_by_pt(fout, num_pts_x, num_pts_y);
  return 0;
}

void DirectEdgeDGArray::write_raw(std::ofstream& fout) const {
  for(int i=0; i<num_dofs; i++) {
    fout << the_array[i] << "  ";
    if( !(i%10) ) fout << "\n";
  }
};


int DirectEdgeDGArray::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}
