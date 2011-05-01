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
  
  register_evolved ("Refluxing::densflux_register_fine",
                    "Refluxing::densflux_stored");
  register_evolved ("Refluxing::sconflux_register_fine",
                    "Refluxing::sconflux_stored");
  register_evolved ("Refluxing::tauflux_register_fine",
                    "Refluxing::tauflux_stored");
  if (CCTK_Equals(Y_e_evolution_method, "GRHydro")) {
    register_evolved ("Refluxing::yeflux_register_fine",
                      "Refluxing::yeflux_stored");
  }
  
  register_evolved ("Refluxing::densflux_register_coarse",
                    "Refluxing::densflux_stored");
  register_evolved ("Refluxing::sconflux_register_coarse",
                    "Refluxing::sconflux_stored");
  register_evolved ("Refluxing::tauflux_register_coarse",
                    "Refluxing::tauflux_stored");
  if (CCTK_Equals(Y_e_evolution_method, "GRHydro")) {
    register_evolved ("Refluxing::yeflux_register_coarse",
                      "Refluxing::yeflux_stored");
  }
}
