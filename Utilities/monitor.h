#ifndef __MONITOR_H_INCLUDED__
#define __MONITOR_H_INCLUDED__

#include <string>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////
//  Monitoring support class
//
// Monitor class
//
// Code monitoring output (to follow progress as the code runs).
// Print indented output if  level + base_level <= level_to_monitor .

class Monitor {
private:
  int base_level;       // base (offset) level;
  int level_to_monitor; // Levels that get monotoring output

public:
  Monitor(int base=0, int level_to=0) : base_level(base), level_to_monitor(level_to) {};
  Monitor(const Monitor& m, int change_base=0) : base_level(m.base_level+change_base),
						 level_to_monitor(m.level_to_monitor) {};

  void setBaseLevel(int base) { base_level = base; };
  void setLevelToMonitor(int level_to) { level_to_monitor = level_to; };
  void reset(int base=0, int level_to=0) { base_level = base; level_to_monitor = level_to; };
  void changeBaseLevel(int change=1) { base_level += change; };

  void operator() () { std::cout << std::endl; };
  void operator() (int level, std::string str);
  void operator() (int level, std::string str, double x);
  void operator() (int level, std::string str, int    x);
  void operator() (int level, std::string str1, double x1, std::string str2, double x2);
  void operator() (int level, std::string str1, int    x1, std::string str2, double x2);
  void operator() (int level, std::string str1, double x1, std::string str2, int    x2);
  void operator() (int level, std::string str1, int    x1, std::string str2, int    x2);

  int baseLevel() const { return base_level; };
  int levelToMonitor() const { return level_to_monitor; };
};

#endif
