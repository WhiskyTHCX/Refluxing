#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine Refluxing_Init (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  densflux_register_fine = 0
  sxflux_register_fine   = 0
  syflux_register_fine   = 0
  szflux_register_fine   = 0
  tauflux_register_fine  = 0
  
  densflux_register_coarse = 0
  sxflux_register_coarse   = 0
  syflux_register_coarse   = 0
  szflux_register_coarse   = 0
  tauflux_register_coarse  = 0
  
  if (refluxing_debug_variables /= 0) then
     dens_correction_total = 0
     sx_correction_total   = 0
     sy_correction_total   = 0
     sz_correction_total   = 0
     tau_correction_total  = 0
  end if
end subroutine Refluxing_Init
