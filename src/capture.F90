#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine Refluxing_CaptureFluxes (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  ! Shift the fluxes, so that each cell contains both its
  ! cell-centered value and its left face flux. (GRHydro uses the
  ! opposite convention, where each cell contains its right face
  ! flux.) In other words, here we define that flux_i is the flux to
  ! the left of density_i, whereas GRHydro defines that flus_(i-1) is
  ! the flux to the left of density_i.
  integer, parameter :: grhydro_offset = +1

  integer :: dir
  
  dir = flux_direction
  densflux_stored(:,:,:,dir) = eoshift(densflux, -grhydro_offset, dim=dir)
  sxflux_stored  (:,:,:,dir) = eoshift(sxflux  , -grhydro_offset, dim=dir)
  syflux_stored  (:,:,:,dir) = eoshift(syflux  , -grhydro_offset, dim=dir)
  szflux_stored  (:,:,:,dir) = eoshift(szflux  , -grhydro_offset, dim=dir)
  tauflux_stored (:,:,:,dir) = eoshift(tauflux , -grhydro_offset, dim=dir)
end subroutine Refluxing_CaptureFluxes
