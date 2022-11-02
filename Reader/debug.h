#ifndef __DEBUG_H_INCLUDED__
#define __DEBUG_H_INCLUDED__

#include <iostream>

////////////////////////////////////////////////////////////////////////////////
//  Debugging support functions and macros

#define PLACE() {std::cerr<<"File "<<__FILE__<<", Line "<<__LINE__<<std::endl;}

inline void MY_DEBUG() PLACE()

template <typename TYPE>
inline void MY_DEBUG(const char* name, const TYPE& value) {
  std::cerr << name << ": " << value << std::endl;
}

template <typename TYPE>
inline void MY_DEBUG(const char* name, const TYPE* a, int n, int nPerLine=5) {
  std::cerr << name;
  if(a == 0)
    std::cerr << " is null";
  else {
    std::cerr << ":\n";
    for(int i=0; i<n; i++) {
      std::cerr << a[i];
      if((i+1)%nPerLine) {
	std::cerr << "  ";
      } else if(i+1<n) {
	std::cerr << "\n";
      }
    }
  }
  std::cerr << std::endl;
}

#endif
