#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine Refluxing_CaptureFluxes (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  densflux_stored(:,:,:,flux_direction) = densflux
  sxflux_stored  (:,:,:,flux_direction) = sxflux
  syflux_stored  (:,:,:,flux_direction) = syflux
  szflux_stored  (:,:,:,flux_direction) = szflux
  tauflux_stored (:,:,:,flux_direction) = tauflux
end subroutine Refluxing_CaptureFluxes
