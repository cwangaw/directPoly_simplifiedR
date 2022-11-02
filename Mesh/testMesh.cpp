#include <iostream>
using namespace std;

#include "debug.h"
#include "baseObjects.h"
#include "polyMesh.h"
using namespace polymesh;

int main() {
  
  const int nVertices = 8;
  double pts[2*nVertices] = {0,0, 10,0, 10,5, 6,5, 5,8, 0,10, 5,3, 10,10};

  const int nElements = 5;
  int numEdgesPerElement[5] = {3,4,3,5,4};
  int* elementByVertexIndex[5];

  const int iOff = -1;
  
  elementByVertexIndex[0] = new int[3];
  elementByVertexIndex[0][0] = 1+iOff;
  elementByVertexIndex[0][1] = 2+iOff;
  elementByVertexIndex[0][2] = 7+iOff;

  elementByVertexIndex[1] = new int[4];
  elementByVertexIndex[1][0] = 3+iOff;
  elementByVertexIndex[1][1] = 4+iOff;
  elementByVertexIndex[1][2] = 7+iOff;
  elementByVertexIndex[1][3] = 2+iOff;

  elementByVertexIndex[2] = new int[3];
  elementByVertexIndex[2][0] = 5+iOff;
  elementByVertexIndex[2][1] = 8+iOff;
  elementByVertexIndex[2][2] = 6+iOff;

  elementByVertexIndex[3] = new int[5];
  elementByVertexIndex[3][0] = 5+iOff;
  elementByVertexIndex[3][1] = 6+iOff;
  elementByVertexIndex[3][2] = 1+iOff;
  elementByVertexIndex[3][3] = 7+iOff;
  elementByVertexIndex[3][4] = 4+iOff;

  elementByVertexIndex[4] = new int[4];
  elementByVertexIndex[4][0] = 8+iOff;
  elementByVertexIndex[4][1] = 5+iOff;
  elementByVertexIndex[4][2] = 4+iOff;
  elementByVertexIndex[4][3] = 3+iOff;

  PolyMesh mesh(nVertices, pts, nElements, numEdgesPerElement, elementByVertexIndex);

  std::string file("testMesh.txt");
  mesh.write_raw(file);

  std::string fileName("testMesh");
  mesh.write_matlab(fileName);
}
