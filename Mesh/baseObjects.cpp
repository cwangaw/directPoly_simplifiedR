#include "baseObjects.h"
#include "debug.h"

////////////////////////////////////////////////////////////////////////////////
// Class Point

Point operator+(const Point& p0, const Point& p1) {
  Point p(p0); p += p1; return p;
};

Point operator-(const Point& p0, const Point& p1) {
  Point p(p0); p -= p1; return p;
};

Point operator*(double scalar, const Point& p0) {
  Point p(p0); p *= scalar; return p;
};

Point operator*(const Point& p0, double scalar) {
  Point p(p0); p *= scalar; return p;
};

Point operator/(const Point& p0, double scalar) {
  Point p(p0); p /= scalar; return p;
};

////////////////////////////////////////////////////////////////////////////////
// Class Tensor1

static const double pi = 3.141592653589793238463; 

Tensor1& Tensor1::rotateCCW(double degrees) {
  double radians = pi*degrees/180;
  double c = cos(radians);
  double s = sin(radians);

  double x = c*the_point[0] - s*the_point[1];
  double y = s*the_point[0] + c*the_point[1];

  the_point[0] = x;
  the_point[1] = y;

  return *this;
};

////////////////////////////////////////////////////////////////////////////////
// Class Tensor2

void Tensor2::mult(double scalar, Tensor2& result) const {
  result.the_tensor[0][0] = scalar*the_tensor[0][0];
  result.the_tensor[0][1] = scalar*the_tensor[0][1];
  result.the_tensor[1][0] = scalar*the_tensor[1][0];
  result.the_tensor[1][1] = scalar*the_tensor[1][1];
}

void Tensor2::mult(const Tensor1& v, Tensor1& result) const {
  result[0] = the_tensor[0][0]*v[0] + the_tensor[0][1]*v[1];
  result[1] = the_tensor[1][0]*v[0] + the_tensor[1][1]*v[1];
}
  
void Tensor2::mult(const Point& v, Point& result) const {
  result[0] = the_tensor[0][0]*v[0] + the_tensor[0][1]*v[1];
  result[1] = the_tensor[1][0]*v[0] + the_tensor[1][1]*v[1];
}
  
void Tensor2::mult(const Tensor2& t, Tensor2& result) const {
  result.the_tensor[0][0] = the_tensor[0][0]*t.the_tensor[0][0] + the_tensor[0][1]*t.the_tensor[1][0];
  result.the_tensor[0][1] = the_tensor[0][0]*t.the_tensor[0][1] + the_tensor[0][1]*t.the_tensor[1][1];
  result.the_tensor[1][0] = the_tensor[1][0]*t.the_tensor[0][0] + the_tensor[1][1]*t.the_tensor[1][0];
  result.the_tensor[1][1] = the_tensor[1][0]*t.the_tensor[0][1] + the_tensor[1][1]*t.the_tensor[1][1];
}

Tensor2& Tensor2::operator*=(const Tensor2& t) {
  double t00 = the_tensor[0][0];
  the_tensor[0][0] = t00*t.the_tensor[0][0] + the_tensor[0][1]*t.the_tensor[1][0];
  the_tensor[0][1] = t00*t.the_tensor[0][1] + the_tensor[0][1]*t.the_tensor[1][1];
  double t10 = the_tensor[1][0];
  the_tensor[1][0] = t10*t.the_tensor[0][0] + the_tensor[1][1]*t.the_tensor[1][0];
  the_tensor[1][1] = t10*t.the_tensor[0][1] + the_tensor[1][1]*t.the_tensor[1][1];

  return *this;
}

bool Tensor2::inverse(Tensor2& invT) const {
  double det = determinant();
  if(det == 0) return false;

  invT(0,0) =  the_tensor[1][1]/det;
  invT(0,1) = -the_tensor[0][1]/det;
  invT(1,0) = -the_tensor[1][0]/det;
  invT(1,1) =  the_tensor[0][0]/det;

  return true;
}

bool Tensor2::realEigenvalues(double& lam0, double& lam1) const {
  double trace = the_tensor[0][0] + the_tensor[1][1];
  double det = determinant();

  if(det == 0) {
    lam0 = (trace < 0) ? trace : 0;
    lam1 = (trace < 0) ? 0 : trace;
    return true;
  }
  
  double discriminant = trace*trace - 4*det;
  if(discriminant < 0) return false;
  discriminant = std::sqrt( trace*trace - 4*det);
  lam0 = ( trace - discriminant ) / 2;
  lam1 = ( trace + discriminant ) / 2;
  return true;
}

bool Tensor2::realEigenvectors(double& lam0, Tensor1& v0, double& lam1, Tensor1& v1) const {
  bool isReal = realEigenvalues(lam0, lam1);
  if(!isReal) return false;

  v0.set(the_tensor[0][1], lam0 - the_tensor[0][0]);
  v1.set(the_tensor[0][1], lam1 - the_tensor[0][0]);
  return true;
}
