#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine Refluxing_Init (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  integer :: epoch, reflevel, reflevels
  character(1000) :: msg
  
  epoch     = GetRegriddingEpoch(cctkGH)
  reflevel  = GetRefinementLevel(cctkGH)
  reflevels = GetRefinementLevels(cctkGH)
  write (msg, '("Initialising refluxing information on level ",i2," of ",i2," (epoch ",i6,")")') reflevel, reflevels, epoch
  call CCTK_INFO (msg)
  
  ! TODO: Ensure that we are not called "in the middle of things",
  ! i.e. that we are not asked to reset a level for which we are
  ! currently integrating fluxes, i.e. which is in between coarse grid
  ! time steps.
  
  ! TODO: When initialising level L, initialise the coarse grid
  ! registers stored on level L and the fine grid registers stored on
  ! level L+1. The code below may be resetting a level that needs to
  ! be left alone.
  densflux_register_fine = 0
  sxflux_register_fine   = 0
  syflux_register_fine   = 0
  szflux_register_fine   = 0
  tauflux_register_fine  = 0
  if (CCTK_EQUALS(Y_e_evolution_method, "GRHydro")) then
     yeflux_register_fine  = 0
  end if
  if (CCTK_EQUALS(Bvec_evolution_method, "GRHydro")) then
    Bconsxflux_register_fine   = 0
    Bconsyflux_register_fine   = 0
    Bconszflux_register_fine   = 0
  end if
  
  densflux_register_coarse = 0
  sxflux_register_coarse   = 0
  syflux_register_coarse   = 0
  szflux_register_coarse   = 0
  tauflux_register_coarse  = 0
  if (CCTK_EQUALS(Y_e_evolution_method, "GRHydro")) then
     yeflux_register_coarse  = 0
  end if
  if (CCTK_EQUALS(Bvec_evolution_method, "GRHydro")) then
    Bconsxflux_register_coarse   = 0
    Bconsyflux_register_coarse   = 0
    Bconszflux_register_coarse   = 0
  end if
  
  if (refluxing_debug_variables /= 0) then
     dens_correction_total = 0
     sx_correction_total   = 0
     sy_correction_total   = 0
     sz_correction_total   = 0
     tau_correction_total  = 0
     if (CCTK_EQUALS(Y_e_evolution_method, "GRHydro")) then
        ye_correction_total  = 0
     end if
     if (CCTK_EQUALS(Bvec_evolution_method, "GRHydro")) then
        Bconsx_correction_total   = 0
        Bconsy_correction_total   = 0
        Bconsz_correction_total   = 0
     end if
  end if
  
  ! Initialise some variables so that they don't appear uninitialised
  ! to NaNChecker, although they are really unused at this point (and
  ! hence also uninitialised).
  
  densflux_correction = -1
  sxflux_correction   = -1
  syflux_correction   = -1
  szflux_correction   = -1
  tauflux_correction  = -1
  if (CCTK_EQUALS(Y_e_evolution_method, "GRHydro")) then
     yeflux_correction  = -1
  end if
  if (CCTK_EQUALS(Bvec_evolution_method, "GRHydro")) then
    Bconsxflux_correction   = -1
    Bconsyflux_correction   = -1
    Bconszflux_correction   = -1
  end if
  
  densflux_stored = -1
  sxflux_stored   = -1
  syflux_stored   = -1
  szflux_stored   = -1
  tauflux_stored  = -1
  if (CCTK_EQUALS(Y_e_evolution_method, "GRHydro")) then
     yeflux_stored  = -1
  end if
  if (CCTK_EQUALS(Bvec_evolution_method, "GRHydro")) then
    Bconsxflux_stored   = -1
    Bconsyflux_stored   = -1
    Bconszflux_stored   = -1
  end if
  
  if (refluxing_debug_variables /= 0) then
     flux_weight_fine   = -1
     flux_weight_coarse = -1
  end if
end subroutine Refluxing_Init
