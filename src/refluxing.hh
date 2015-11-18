#ifndef REFLUXING_HH
#define REFLUXING_HH

#include <cctk.h>

#include <vector>


namespace Refluxing {
  using namespace std;
  // List of variable to which refluxing should be applied
  extern vector<int> refluxing_vars;
} // namespace Refluxing


#endif  // #ifndef REFLUXING_HH
