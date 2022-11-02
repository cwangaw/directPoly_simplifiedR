#ifndef __parameterData_h_included__
#define __parameterData_h_included__

#include <string>
#include <cmath>

#include "Mesh/baseObjects.h"
#include "Mesh/polyMesh.h"
#include "directSerendipity.h"
#include "directMixed.h"
#include "debug.h"

class ParameterData {
public:
  ParameterData(); 
  ParameterData(const std::string& infileName);
  ParameterData(const std::string& infileName,
		const std::string& echofileName, bool echo0=1);
  ~ParameterData() { };

  int read();

  std::string infile_name;
  std::string echofile_name;
  bool echo;

  // PARAMETER CASE NUMBERS =====================================================

  void print_cases(std::ostream& = std::cout);

  enum case_soln_DS_output
    { case_soln_DS_output_omit=-1, case_soln_DS_output_none=0,
      case_soln_DS_output_raw, case_soln_DS_output_matlab,
      n_case_soln_DS_output };

  enum case_soln_Mixed_output
    { case_soln_Mixed_output_omit=-1, case_soln_Mixed_output_none=0,
      case_soln_Mixed_output_raw, case_soln_Mixed_output_matlab,
      n_case_soln_Mixed_output };

  enum case_mesh_output
    { case_mesh_output_none=0, case_mesh_output_raw, case_mesh_output_matlab,
      n_case_mesh_output };
  
  enum case_dsSpace_output
    { case_dsSpace_output_none=0, case_dsSpace_output_raw, case_dsSpace_output_matlab,
      n_case_dsSpace_output };
  
  enum case_dmSpace_output
    { case_dmSpace_output_none=0, case_dmSpace_output_raw, case_dmSpace_output_matlab,
      n_case_dmSpace_output };

  // DATA =======================================================================

  // OUTPUT DIRECTORY
  std::string directory_name; //output directory (preamble)

  // MESH PARAMETERS
  std::string mesh_type;
  int nVertices;
  int nElements;
  std::string mesh_vertices_fileName;
  std::string mesh_elements_fileName;
  double distortion_factor;
  polymesh::PolyMesh mesh;
  double mesh_small_edge_ratio;

  // FINITE ELEMENTS
  int polynomial_degree;
  int conforming;
  directserendipity::DirectSerendipity dsSpace;
  directserendipity::DirectMixed dmSpace;
  // ALGORITHM PARAMETERS
  int refinement_level;

  // OUTPUT PARAMETERS
  int output_soln_DS_format;
  int output_mesh_numPts_DS_x, output_mesh_numPts_DS_y;

  int output_soln_Mixed_format;
  int output_mesh_numPts_Mixed_x, output_mesh_numPts_Mixed_y;
  
  int output_mesh_format;
  int output_dsSpace_format;
  int output_dmSpace_format;

  int monitor_to_level; // level up to which to monitor code progress

 private:
  int readMesh();
  int writeMesh();
  int writeDSSpace();
  int writeDMSpace();
};

#endif
