#include <iostream>
using namespace std;

#include "debug.h"
#include "baseObjects.h"

int main() {
  cout << "Point and Tensor1\n";
  
  Tensor1 p0(2,3);
  Tensor1 p1(4,-2);
  cout << "p0 =" << p0 << ", norm = " << p0.norm() << endl;
  cout << "p1 =" << p1 << endl;

  cout << "p0[4] = " << p0[4] << ", p0[3] = " << p0[3] << endl;
  
  Point p(p1+p0);
  cout << "p = p0 + p1 = " << p << endl;
  
  p = p1 - p0;
  cout << "p = p1 - p0 = " << p << endl;

  p = 2*p0 + p1/3;
  cout << " p = 2*p0 + p1/3 = " << p << endl;
				       
  cout << "\n\nTensor2\n";

  Tensor2 t0(2,3,-1,0);
  Tensor2 t1(4,-2,-2,1);
  cout << "t0 = " << t0 << "  det = " << t0.determinant() << endl;
  cout << "t1 = " << t1 << "  det = " << t1.determinant() << endl;
  cout << "t1(0,1) = " << t1(0,1) << ", t1(3,5) = " << t1(3,5) << endl;
  
  Tensor2 t(t0); t += t1;
  cout << "t = t0 + t1 = " << t << endl;
  t -= t1;
  cout << "t -= t1 = t0 = " << t << endl;
  t *= 3;
  cout << "t *= 3 = " << t << endl;
  t /= 2;
  cout << "t /= 2 = " << t << endl;
  t = t1;
  cout << "t = t1 = " << t << endl;

  cout << "(t1==t1) (t0==t1) (t1!=t1) (t0!=t1) = "
       << (t1==t1) << " " << (t0==t1) << " " << (t1!=t1) << " " << (t0!=t1) << endl;

  t0.mult(t1,t);
  cout << "t = t0*t1 = " << t << endl;

  Tensor1 pp;
  t0.mult(p0,pp);
  cout << "p = t0*p0 = " << pp << endl;

  t0.mult(5,t);
  cout << "t = 5*t0 = " << t << endl;

  t = t0;
  t *= t1;
  cout << "t = t0, t *= t1 = " << t << endl;

  double lam0, lam1;
  if(t0.realEigenvectors(lam0,p0,lam1,p1)) {
    cout << "t0: lam0  p0 = " << lam0 << "  " << p0 << "  lam1  p1 = " << lam1 << "  " << p1 << endl; 
  }
  if(t1.realEigenvectors(lam0,p0,lam1,p1)) {
    cout << "t1: lam0  p0 = " << lam0 << "  " << p0 << "  lam1  p1 = " << lam1 << "  " << p1 << endl; 
  }

 Tensor1 ttt(1,std::sqrt(3));
 cout << "ttt = " << ttt << ", rotated 90 CCW = " << ttt.rotateCCW()
       << ", rotated 60 CCW = " << ttt.rotateCCW(60) << "\n";  
 
  Tensor2 t2(3,1,2,4);
  if(t2.realEigenvectors(lam0,p0,lam1,p1)) {
    cout << "t2 = " << t2 << " : lam0  p0 = "
	 << lam0 << "  " << p0 << "  lam1  p1 = " << lam1 << "  " << p1 << endl; 
  }
  
  if(t0.inverse(t)) {
    cout << "t0^{-1} = " << t << endl;
  } else {
    cout << "t0 not invertible" << endl;
  }
  if(t1.inverse(t)) {
    cout << "t1^{-1} = " << t << endl;
  } else {
    cout << "t1 not invertible" << endl;
  }
}	  
