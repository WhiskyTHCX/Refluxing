#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cstdlib>
#include <string>

#include <carpet.hh>

using namespace std;
using namespace Carpet;



struct args_t {
  cGH const* cctkGH;
};



static
void
output(const int vi,
       const char *const options,
       void *const args_)
{
  args_t& args = *static_cast<args_t*>(args_);
  const cGH *const cctkGH = args.cctkGH;
  char *const fullname = CCTK_FullName(vi);
  const string alias = string("Refluxing::") + string(fullname);
  CCTK_OutputVarAs(cctkGH, fullname, alias.c_str());
  free(fullname);
}



extern "C"
void Refluxing_Output(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  assert(is_level_mode());
  
  CCTK_INFO("Outputting refluxing debug variables");
  
  args_t args;
  args.cctkGH = cctkGH;
  CCTK_TraverseString(debug_variables, output, &args, CCTK_GROUP_OR_VAR);
}
