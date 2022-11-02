#include <iostream>
using namespace std;

#include "monitor.h"

////////////////////////////////////////////////////////////////////////////////
// Monitor class

void Monitor::operator() (int level, std::string str) {
  level += base_level;
  if(level <= level_to_monitor) {
    for(int i = 0; i < level; i++) cout << "  ";
    cout << str << endl;
  }
}

void Monitor::operator() (int level, std::string str, double x) {
  level += base_level;
  if(level <= level_to_monitor) {
    for(int i = 0; i < level; i++) cout << "  ";
    cout << str << x << endl;
  }
}

void Monitor::operator() (int level, std::string str, int x) {
  level += base_level;
  if(level <= level_to_monitor) {
    for(int i = 0; i < level; i++) cout << "  ";
    cout << str << x << endl;
  }
}

void Monitor::operator() (int level, std::string str1, double x1, std::string str2, double x2) {
  level += base_level;
  if(level <= level_to_monitor) {
    for(int i = 0; i < level; i++) cout << "  ";
    cout << str1 << x1 << str2 << x2 << endl;
  }
}

void Monitor::operator() (int level, std::string str1, int    x1, std::string str2, double x2) {
  level += base_level;
  if(level <= level_to_monitor) {
    for(int i = 0; i < level; i++) cout << "  ";
    cout << str1 << x1 << str2 << x2 << endl;
  }
}

void Monitor::operator() (int level, std::string str1, double x1, std::string str2, int    x2) {
  level += base_level;
  if(level <= level_to_monitor) {
    for(int i = 0; i < level; i++) cout << "  ";
    cout << str1 << x1 << str2 << x2 << endl;
  }
}

void Monitor::operator() (int level, std::string str1, int    x1, std::string str2, int    x2) {
  level += base_level;
  if(level <= level_to_monitor) {
    for(int i = 0; i < level; i++) cout << "  ";
    cout << str1 << x1 << str2 << x2 << endl;
  }
}
