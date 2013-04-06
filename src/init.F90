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
  call CCTK_INFO(msg)
  
  ! TODO: Ensure that we are not called "in the middle of things",
  ! i.e. that we are not asked to reset a level for which we are
  ! currently integrating fluxes, i.e. which is in between coarse grid
  ! time steps.
  
  ! TODO: When initialising level L, initialise the coarse grid
  ! registers stored on level L and the fine grid registers stored on
  ! level L+1. The code below may be resetting a level that needs to
  ! be left alone.
  
  if (suppress_refluxing_in_atmosphere.ne.0) then
     reflux_atmosphere_mask = 0.0d0
  end if
  
  register_fine   = 0
  register_coarse = 0
  if (delayed_refluxing /= 0) then
     delayed_correction = 0
  end if
  
  if (refluxing_debug_variables /= 0) then
     correction_total = 0
  end if
  
  ! Initialise some variables so that they don't appear uninitialised
  ! to NaNChecker, although they are really unused at this point (and
  ! hence also uninitialised).
  
  flux       = -1
  correction = -1
  
  if (refluxing_debug_variables /= 0) then
     flux_weight_fine   = -1
     flux_weight_coarse = -1
  end if
end subroutine Refluxing_Init
