#include "refluxing.hh"

#include <carpet.hh>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <vector>



namespace Refluxing {
  
  using namespace std;
  
  
  
  // List of variable to which refluxing should be applied
  vector<int> refluxing_vars;
  
  
  
  // Find out whether a variable should be refluxed. We reflux those
  // variables that have storage.
  static
  void
  add_refluxing_var(const int vi, const char *const optstring,
                    void *const callback_arg)
  {
    const cGH *const cctkGH CCTK_ATTRIBUTE_UNUSED =
      static_cast<const cGH*>(callback_arg);
    assert(vi >= 0);
    bool active = true;
    if (optstring) {
      char *const parameter = strdup(optstring);
      char *const scope_pos = strstr(parameter, "::");
      assert(scope_pos);
      *scope_pos = '\0';
      const char *const thorn = parameter;
      const char *const name = scope_pos+2;
      int type;
      const void *const value_ptr = CCTK_ParameterGet(name, thorn, &type);
      assert(value_ptr);
      assert(type == PARAMETER_KEYWORD or type == PARAMETER_STRING);
      free(parameter);
      const char *const value =
        *static_cast<const char**>(const_cast<void*>(value_ptr));
      assert(value);
      active = strcmp(value, "none");
    }
    if (active) {
      char* const fullname = CCTK_FullName(vi);
      CCTK_VInfo(CCTK_THORNSTRING,
                 "   %d: %s", int(refluxing_vars.size()), fullname);
      free(fullname);
      refluxing_vars.push_back(vi);
    }
  }
  
  
  
  static
  void register_evolved(const int var, const int rhs)
  {
    DECLARE_CCTK_PARAMETERS;
    
    CCTK_INT (*register_var)(CCTK_INT EvolvedIndex, CCTK_INT RHSIndex);
    if (use_MoL_slow_multirate_sector) {
      register_var = MoLRegisterEvolvedSlow;
    } else {
      register_var = MoLRegisterEvolved;
    }
    
    int const ierr = register_var(var, rhs);
    if (ierr) {
      CCTK_ERROR("Could not register evolved variable");
    }
  }
  
  
  
  extern "C"
  void
  Refluxing_SetupVars(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    // Determine list of variables that should be refluxed
    CCTK_INFO("Determining refluxing variables:");
    assert(refluxing_vars.empty());
    const int iret =
      CCTK_TraverseString(refluxing_variables, add_refluxing_var,
                          const_cast<cGH*>(cctkGH), CCTK_GROUP_OR_VAR);
    assert(iret>=0);
    assert(nvars == int(refluxing_vars.size()));
    
    // Determine indices for GRHydro's variables (if any)
    *index_dens   = -1;
    *index_sx     = -1;
    *index_sy     = -1;
    *index_sz     = -1;
    *index_tau    = -1;
    *index_ye     = -1;
    *index_Bconsx = -1;
    *index_Bconsy = -1;
    *index_Bconsz = -1;
    for (int n=0; n<nvars; ++n) {
      const int vi = refluxing_vars.AT(n);
      if (vi == CCTK_VarIndex("GRHydro::dens"    )) *index_dens   = n;
      if (vi == CCTK_VarIndex("GRHydro::scon[0]" )) *index_sx     = n;
      if (vi == CCTK_VarIndex("GRHydro::scon[1]" )) *index_sy     = n;
      if (vi == CCTK_VarIndex("GRHydro::scon[2]" )) *index_sz     = n;
      if (vi == CCTK_VarIndex("GRHydro::tau"     )) *index_tau    = n;
      if (vi == CCTK_VarIndex("GRHydro::Y_e_con" )) *index_ye     = n;
      if (vi == CCTK_VarIndex("GRHydro::Bcons[0]")) *index_Bconsx = n;
      if (vi == CCTK_VarIndex("GRHydro::Bcons[1]")) *index_Bconsy = n;
      if (vi == CCTK_VarIndex("GRHydro::Bcons[2]")) *index_Bconsz = n;
    }
    
    // Register these variables with MoL
    const int vi_register_fine =
      CCTK_VarIndex("Refluxing::register_fine[0]");
    assert(vi_register_fine >= 0);
    const int vi_register_coarse =
      CCTK_VarIndex("Refluxing::register_coarse[0]");
    assert(vi_register_coarse >= 0);
    const int vi_flux =
      CCTK_VarIndex("Refluxing::flux[0]");
    assert(vi_flux >= 0);
    for (int n=0; n<nvars; ++n) {
      // some time updates are handled by thorns other than MoL. Only register
      // those fluxes with MoL whose variable is evolved by MoL.
      if (MoLQueryEvolvedRHS(refluxing_vars.AT(n)) >= 0) {
        for (int d=0; d<3; ++d) {
          register_evolved(vi_register_fine+3*n+d, vi_flux+3*n+d);
          register_evolved(vi_register_coarse+3*n+d, vi_flux+3*n+d);
        }
      }
    }
  }
  
} // namespace Refluxing
