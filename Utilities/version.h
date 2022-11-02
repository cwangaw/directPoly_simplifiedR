#ifndef __version_h_included__
#define __version_h_included__

#include <string>

class Version
{
public:
  Version(std::string name0, std::string date0, int version_depth0)
    : name(name0), date(date0), version_depth(version_depth0) {
    version = new int[version_depth]; }

  ~Version() { delete[] version; }

  void set_version(int i0, int i1=-1, int i2=-1, int i3=-1, int i4=-1,
		   int i5=-1, int i6=-1, int i7=-1, int i8=-1, int i9=-1) {
    int v[10]={i0,i1,i2,i3,i4,i5,i6,i7,i8,i9};
    for(int i=0; i<version_depth; i++) version[i]=v[i];
  }

  void set_description(std::string descr) { description = descr; };

  std::string name;
  std::string date;
  std::string description;
  int version_depth;
  int* version;
};

#endif
