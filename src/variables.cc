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
  vector<int> refluxing_vars;
}
using namespace Refluxing;

extern "C"
int
Refluxing_RegisterVariable(int const vi_var, int const use_slow_sector,
                           int * fluxid) {
  DECLARE_CCTK_PARAMETERS
  int ierr = 0;

  // Flux index
  *fluxid = refluxing_vars.size();
  assert(*fluxid < refluxing_nvars);

  // some time updates are handled by thorns other than MoL. Only register
  // those fluxes with MoL whose variable is evolved by MoL.
  if (MoLQueryEvolvedRHS(vi_var) >= 0) {
    const int vi_register_fine =
      CCTK_VarIndex("Refluxing::register_fine[0]");
    assert(vi_register_fine >= 0);
    const int vi_register_coarse =
      CCTK_VarIndex("Refluxing::register_coarse[0]");
    assert(vi_register_coarse >= 0);
    const int vi_flux =
      CCTK_VarIndex("Refluxing::flux[0]");
    assert(vi_flux >= 0);

    CCTK_INT (*register_var)(CCTK_INT EvolvedIndex, CCTK_INT RHSIndex);
    if (use_slow_sector != 0) {
      register_var = MoLRegisterEvolvedSlow;
    } else {
      register_var = MoLRegisterEvolved;
    }

    for(int d = 0; d < 3; ++d) {
      int const shift = 3*(*fluxid) + d;
      ierr |= register_var(vi_register_fine   + shift, vi_flux + shift);
      ierr |= register_var(vi_register_coarse + shift, vi_flux + shift);
    }
  }

  // Register variable that needs to be refluxed
  refluxing_vars.push_back(vi_var);

  // Return flux index
  return ierr;
}
