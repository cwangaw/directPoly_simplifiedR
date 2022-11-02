///////////////////////////////////////////////////////////////////////////////
//
//        DIRECT SERENDIPITY AND MIXED FINITE ELEMENTS ON POLYGONS
//
//  Todd Arbogast <arbogast@ices.utexas.edu>
//  Chuning Wang <cwangaw@utexas.edu>
//
//  Center for Subsurface Modeling
//  Oden Institute for Computational Engineering and Sciences
//  The University of Texas at Austin
//  Austin, Texas  78712
//
//  Begun: July 2020
//
//=============================================================================
//  EQUATIONS
//
//    - div D grad p + div(b p) + c.grad p + a p = f
//                       in Omega (bounded, connected, Lipschitz domain in R^2)
//    p = g              on the boundary of Omega
//    where D is a tensor (matrix), b and c are vectors.
//
//    Find p in DS_{r,g} such that
//      (D grad p , grad q) + (b p , grad q) + (c.grad p , q) + (a p , q)
//          = (f , q) for all q in DS_{r,0}
//    where DS_r are the direct serendipity spaces and DS_{r,g} has the
//      boundary nodes set to g. We use a nodal basis for the space.
//
//    For mixed formulation, we ask c=0.    
//    
//    Find (u,p,lambda) in V_r^s X W_s X Lambda_r  such that
//    Setting u = - D grad p + b p,
//    find (u, p, l) in V^s_r x W_s x L_r (hybrid form) such that
//      (D^{-1} u , v) - (p , div v) + (D^{-1} b p, v)
//                  + Sum_E (l , v.nu)_{delta E} = 0 for all v in V^s_r
//      (div u, q) + (a p, q) = (f , q) for all q in W_s
//      Sum_E (u.nu , l)_{delta E} = 0 for all l in L_r
//    where 
//      V^s_r is the direct sum of V^s_r(E), which is 
//        a full / reduced direct mixed space on element E for s = r, r - 1
//      W_s is the direct sum of W_s(E) = \Po_s(E) restriced to each element
//      L_r is \Po_r(e) restriced to each interior edge e, 0 on boundary edges
//      
//  ELEMENTS E
//    E is a convex, nondegenerate polygon (of at least 3 sides)
//    Polynomials of degree r >= 1
//    The mesh is conforrming
//
//=============================================================================
//
//  ACKNOWLEDGMENTS
//    The development of this code has been supported by
//    the U.S. National Science Foundation.
//
//  COPYRIGHT
//    Copyright (C) 2022 Todd Arbogast and Chuning Wang
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    General Public License for more details (<https://www.gnu.org/licenses/>).
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
#include <sstream>

#include "debug.h"
#include "fcns.h"
#include "Utilities/monitor.h"
#include "parameterData.h"
#include "Utilities/version.h"
#include "Reader/expr.h"
#include "ellipticPDE.h"
#include "mixedPDE.h"
//#include "baseObjects.h"
//#include "polyMesh.h"
//#include "directSerendipity.h"

////////////////////////////////////////////////////////////////////////////////
// MAIN ROUTINES

// OUTPUT HEADER INFORMATION
static void write_header (const Version& code, const ParameterData& param)
{
  std::cout << "\n";
  std::cout << "  DIRECT SERENDIPITY AND MIXED FINITE ELEMENTS ON POLYGONS\n\n";

  std::cout << "           Name:       " << code.name << "\n";
  std::cout << "           Version:    " << code.version[0];
  for(int i=1; i<code.version_depth; i++) std::cout << "." << code.version[i];
  std::cout << "\n";
  std::cout << "           Date:       " << code.date << "\n";
  std::cout << code.description << "\n\n";

  std::cout << "           Input file: " << param.infile_name << "\n";
  if(param.echo)
  std::cout << "           Echo file:  " << param.echofile_name << "\n";
  std::cout << std::endl;
  return;
}

// PARSE COMMAND LINE ARGUMENTS
static int parse_flags (int argc, char** argv, ParameterData& param, Version& code)
{
  for(int i = 1; i < argc; i++) {
    if(argv[i][0] == '-')
      {
	switch(argv[i][1])
	  {
	  case 'h':
	    std::cout << "Usage: " << argv[0] << " [command-options]\n";
	    std::cout << "  Where the optional [command-options] are:\n";
	    std::cout << "    -i [<file name>]   The initial general input <file"
		      << " name>, or \"" << param.infile_name << "\" if omitted\n";
	    std::cout << "    -e [<file name>]   Echo input to <file name>, or \""
		      << param.echofile_name << "\" if omitted\n";
	    std::cout << "    -E                 Echo input to standard out\n";
	    std::cout << "    -c                 Just display input case information\n";
	    std::cout << "    -C                 Just display copyright information\n";
	    std::cout << "    -u                 Just display input physical units information\n";
	    std::cout << "    -v                 Just display version information\n";
	    std::cout << "    -h                 Print this usage information only";
	    std::cout << std::endl;
	    return 1;
	  case 'i':
	    i++;
	    if(i < argc && argv[i][0] != '-') param.infile_name = argv[i]; else i--;
	    break;
	  case 'e':
	    param.echo = 1;
	    i++;
	    if(i < argc && argv[i][0] != '-') param.echofile_name = argv[i]; else i--;
	    break;
	  case 'E':
	    param.echo = 1;
	    param.echofile_name = '\0';
	    break;
	  case 'c':
	    param.print_cases();
	    return 1;
	  case 'C':
	    std::cout << "Copyright (C) 2022 Todd Arbogast and Chuning Wang\n\n";
	    std::cout << "This program is free software: you can redistribute it and/or modify\n";
	    std::cout << "it under the terms of the GNU General Public License as published by\n";
	    std::cout << "the Free Software Foundation, either version 3 of the License, or\n";
	    std::cout << "any later version.\n\n";
	    std::cout << "This program is distributed in the hope that it will be useful,\n";
	    std::cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of\n";
	    std::cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU\n";
	    std::cout << "General Public License for more details (<https://www.gnu.org/licenses/>).\n";
	    return 1;
	  case 'u':
	    expr_listUnits(std::cout);
	    return 1;
	  case 'v':
	    write_header (code, param);
	    return 1;
	  default:
	    std::cerr << "Unrecognized command line option: " << argv[i] << std::endl;
	    return 1;
	  }
      }
    else
      {
	std::cerr << "Unrecognized item on command line: " << argv[i] << "\n";
	std::cerr << "For command line arguments, try: " << argv[0] << " -h"<< std::endl;
	return 1;
      }
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// MAIN

int main(int argc, char* argv[]) {
  // Initialize
  Version code("directpoly","May 2022",3);
  code.set_version(1,1,0);
  code.set_description("Solves a second order PDE using DS and hybrid mixed spaces on polygons\nAllows simple mesh construction");

  ParameterData param;
  Monitor monitor;
  
  try {
    if ( parse_flags(argc, argv, param, code) ) return 0;

    write_header (code, param);

    if ( param.read () ) return 1;
    monitor.reset (0, param.monitor_to_level);
    monitor (0, "THE DATA WAS READ");

    if(param.output_soln_DS_format != ParameterData::case_soln_DS_output_omit) {
      monitor (0, "\nBEGIN DS COMPUTATION\n");

      EllipticPDE ellipticPDE(param);
      ellipticPDE.solve(monitor);

      monitor(0,"\nEND DS COMPUTATION");
    }

    if(param.output_soln_Mixed_format != ParameterData::case_soln_Mixed_output_omit) {
      monitor(0,"\nBEGIN MIXED COMPUTATION");

      MixedPDE mixedPDE(param);
      mixedPDE.solve(monitor);

      monitor(0,"\nEND COMPUTATION");
    }
  }
  catch (std::exception &exc) {
    std::cerr << "\n\n----------------------------------------------------\n";
    std::cerr << "ERROR: " << exc.what() << std::endl;
    std::cerr << "----------------------------------------------------" << std::endl;
    return 1;
  }
  catch (int error) {
    std::cerr << "\n\n----------------------------------------------------\n";
    std::cerr << "ERROR, CODE: " << error << std::endl;
    std::cerr << "----------------------------------------------------" << std::endl;
    return 1;
  }
  catch (...) {
    std::cerr << "\n\n----------------------------------------------------\n";
    std::cerr << "ERROR, UNSPECIFIED" << std::endl;
    std::cerr << "----------------------------------------------------" << std::endl;
    return 1;
  }

  return 0;
}
