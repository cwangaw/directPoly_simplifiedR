# directPoly
This is a demonstration code for the paper Direct Serendipity and Mixed Finite Elements on Convex Polygons.

## Authors
- [Todd Arbogast](arbogast@ices.utexas.edu)
- [Chuning Wang](cwangaw@utexas.edu)


## Build and Run
You should make the files to build an executable file directpoly.$(TARG), where TARG is a variable containing the instruction set architecture of your processor. Note that you need CBLAS and LAPACKE library to successfully build the file.

```bash
make
./directpoly.$(TARG)
```

## Direct Serendipity Space
Nodal basis functions as defined in the paper on each elements are constructed in `directSerendipityFE.cpp` and assembled in `ellipticPDE.cpp` to serve as global basis functions. However, some modifications could be made easily in `directSerendipityFE.cpp`:

1. If you do not want to remove cell degrees of freedom, find the following lines at three places and comment them out.

```cpp
//Deduct value at interior nodes if needed

for (int k=0; k<nCellNodes(); k++) {
    phi_pt -= phi_e_at_c[k] * value_n[k + num_vertices*polynomial_degree + pt_index*num_nodes];
    gradresult -= phi_e_at_c[k] * gradvalue_n[k + num_vertices*polynomial_degree + pt_index*num_nodes];
}
```

2. If you want to make vertex nodal basis functions linear at each edge, find the following lines and uncomment the block.

```cpp
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
```  

## Direct Mixed Space
You can use either hybrid mixed method or H(div)-conforming mixed spaces by setting the corresponding input in `infile`.  

## Interpolation
You could interpolate functions by direct serendipity basis functions and graph the interpolated function as well as its gradient in the following part of `ellipticPDE.cpp`. You can also directly modify the coefficients of basis functions for mixed spaces, as well as discrete Galerkin spaces for elements and edges, then graph their sum in the following part of `mixedPDEConf.cpp` and `mixedPDEHybrid.cpp`.

```cpp
// TEST BASIS FUNCTIONS //////////////////////////////////////////////////
if(true) {
    ...
}
```

## Solve a PDE
Problem formulation is given in the heading comments of `main.cpp`. Coefficients a, b, c, D, and all the other related data could be modified in `fcns.cpp`. Note that the source function as well as Riemann boundary condition are defaultly calculated by true solution. If in your formulation, analytical solution is unknown, please rewrite these parts.

## Output
Basic results of running the codes are printed to the terminal. All the output files would be stored preambly in `test/` directory, which could be modified in the first line of `infile`.

## Acknowledgments
The development of this code has been supported by the U.S. National Science Foundation.

## Copyrights
Copyright (C) 2022 Todd Arbogast and Chuning Wang

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the [GNU General Public License](https://www.gnu.org/licenses/) for more details.