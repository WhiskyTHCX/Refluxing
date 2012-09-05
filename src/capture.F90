#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine Refluxing_Sync(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
end subroutine Refluxing_Sync



subroutine Refluxing_CaptureFluxes (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_REAL, parameter :: poison = 42.0d+42
  
  ! Shift the fluxes, so that each cell contains both its
  ! cell-centered value and its left face flux. (GRHydro uses the
  ! opposite convention, where each cell contains its right face
  ! flux.) In other words, here we define that flux_i is the flux to
  ! the left of density_i, whereas GRHydro defines that flux_(i-1) is
  ! the flux to the left of density_i.
  integer, parameter :: grhydro_offset = -1
  
  integer :: imin(3), imax(3), ioff(3)
  
  ! Region with valid data
  imin(:) = 1           + cctk_nghostzones(:)
  imax(:) = cctk_lsh(:) - cctk_nghostzones(:)
  ! There is one more flux value in this direction
  imax(flux_direction) = imax(flux_direction) + 1
  
  ! Offset between source and destination
  ioff(:) = 0
  ioff(flux_direction) = grhydro_offset
  
  ! Add delayed correction terms to GRHydro fluxes
  ! (Yes, it is intended that this correction will then later be
  !  captured with the fluxes.)
  if (delayed_refluxing /= 0) then
     call correct (densflux, densflux_delayed_correction(:,:,:,flux_direction))
     call correct (sxflux  , sxflux_delayed_correction  (:,:,:,flux_direction))
     call correct (syflux  , syflux_delayed_correction  (:,:,:,flux_direction))
     call correct (szflux  , szflux_delayed_correction  (:,:,:,flux_direction))
     call correct (tauflux , tauflux_delayed_correction (:,:,:,flux_direction))
     if (CCTK_EQUALS(Y_e_evolution_method, "GRHydro")) then
        call correct &
             (Y_e_con_flux, yeflux_delayed_correction (:,:,:,flux_direction))
     end if
     if (CCTK_EQUALS(Bvec_evolution_method, "GRHydro")) then
        call correct &
             (Bconsxflux, Bconsxflux_delayed_correction (:,:,:,flux_direction))
        call correct &
             (Bconsyflux, Bconsyflux_delayed_correction (:,:,:,flux_direction))
        call correct &
             (Bconszflux, Bconszflux_delayed_correction (:,:,:,flux_direction))
     end if
  end if
  
  ! Capture fluxes from GRHydro
  call capture (densflux_stored(:,:,:,flux_direction), densflux)
  call capture (sxflux_stored  (:,:,:,flux_direction), sxflux  )
  call capture (syflux_stored  (:,:,:,flux_direction), syflux  )
  call capture (szflux_stored  (:,:,:,flux_direction), szflux  )
  call capture (tauflux_stored (:,:,:,flux_direction), tauflux )
  if (CCTK_EQUALS(Y_e_evolution_method, "GRHydro")) then
     call capture (yeflux_stored (:,:,:,flux_direction), Y_e_con_flux)
  end if
  if (CCTK_EQUALS(Bvec_evolution_method, "GRHydro")) then
     call capture (Bconsxflux_stored (:,:,:,flux_direction), Bconsxflux)
     call capture (Bconsyflux_stored (:,:,:,flux_direction), Bconsyflux)
     call capture (Bconszflux_stored (:,:,:,flux_direction), Bconszflux)
  end if
  
contains
  
  subroutine correct (grhydro, refluxing)
    CCTK_REAL, intent(inout) :: grhydro(:,:,:)
    CCTK_REAL, intent(in)    :: refluxing(:,:,:)
    integer, parameter :: rk = kind(grhydro)
    integer   :: di,dj,dk
    CCTK_REAL :: factor
    integer   :: i,j,k
    CCTK_REAL :: avg_alp
    
    ! Add fraction of Refluxing delayed correction variable to GRHydro
    ! flux variable
    
    factor = delayed_refluxing_fraction / CCTK_DELTA_TIME
    
    di = ioff(1)
    dj = ioff(2)
    dk = ioff(3)
    !$omp parallel do private(i,j,k, avg_alp)
    do k=imin(3),imax(3)
       do j=imin(2),imax(2)
          do i=imin(1),imax(1)
             avg_alp = 0.5_rk * (alp(i,j,k) + alp(i+di,j+dj,k+dk))
             grhydro(i+di,j+dj,k+dk) = grhydro(i+di,j+dj,k+dk) + &
                  factor * refluxing(i,j,k) / avg_alp
          end do
       end do
    end do
  end subroutine correct
  
  subroutine capture (refluxing, grhydro)
    CCTK_REAL, intent(out) :: refluxing(:,:,:)
    CCTK_REAL, intent(in)  :: grhydro(:,:,:)
    integer, parameter :: rk = kind(refluxing)
    integer   :: di,dj,dk
    integer   :: i,j,k
    CCTK_REAL :: avg_alp
    
    ! Poison Refluxing flux variable
    
    !$omp parallel do private(i,j,k)
    do k=1,size(refluxing,3)
       do j=1,size(refluxing,2)
          do i=1,size(refluxing,1)
             refluxing(i,j,k) = poison
          end do
       end do
    end do
    
    ! Copy interior of GRHydro flux variable to Refluxing flux variable
    
    di = ioff(1)
    dj = ioff(2)
    dk = ioff(3)
    !$omp parallel do private(i,j,k, avg_alp)
    do k=imin(3),imax(3)
       do j=imin(2),imax(2)
          do i=imin(1),imax(1)
             avg_alp = 0.5_rk * (alp(i,j,k) + alp(i+di,j+dj,k+dk))
             refluxing(i,j,k) = avg_alp * grhydro(i+di,j+dj,k+dk)
          end do
       end do
    end do
  end subroutine capture
  
end subroutine Refluxing_CaptureFluxes
