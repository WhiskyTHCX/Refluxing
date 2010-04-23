#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

static
void register_evolved (char const* const group, char const* const rhs)
{
  int const ierr =
    MoLRegisterEvolvedGroup (CCTK_GroupIndex(group),
                             CCTK_GroupIndex(rhs));
  if (ierr) {
    CCTK_WARN (CCTK_WARN_ABORT, "Could not register evolved groups");
  }
}



extern "C"
void Refluxing_Register (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  register_evolved ("refluxing::densflux_register_fine",
                    "refluxing::densflux_stored");
  register_evolved ("refluxing::sconflux_register_fine",
                    "refluxing::sconflux_stored");
  register_evolved ("refluxing::tauflux_register_fine",
                    "refluxing::tauflux_stored");
  
  register_evolved ("refluxing::densflux_register_coarse",
                    "refluxing::densflux_stored");
  register_evolved ("refluxing::sconflux_register_coarse",
                    "refluxing::sconflux_stored");
  register_evolved ("refluxing::tauflux_register_coarse",
                    "refluxing::tauflux_stored");
}
