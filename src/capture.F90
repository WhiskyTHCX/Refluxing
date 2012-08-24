#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine Refluxing_Sync(CCTK_ARGUMENTS)

  implicit none
  DECLARE_CCTK_ARGUMENTS

  return

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
  
  subroutine capture (dest, source)
    CCTK_REAL, intent(out) :: dest  (:,:,:)
    CCTK_REAL, intent(in)  :: source(:,:,:)
    integer, parameter :: rk = kind(dest)
    integer   :: di,dj,dk
    integer   :: i,j,k
    CCTK_REAL :: avg_alp
    
    ! Poison destination array
    
    !$omp parallel do private(i,j,k)
    do k=1,size(dest,3)
       do j=1,size(dest,2)
          do i=1,size(dest,1)
             dest(i,j,k) = poison
          end do
       end do
    end do
    
    ! Copy interior of source to dest
    
    di = ioff(1)
    dj = ioff(2)
    dk = ioff(3)
    !$omp parallel do private(i,j,k, avg_alp)
    do k=imin(3),imax(3)
       do j=imin(2),imax(2)
          do i=imin(1),imax(1)
             avg_alp = 0.5_rk * (alp(i,j,k) + alp(i+di,j+dj,k+dk))
             dest(i,j,k) = avg_alp * source(i+di,j+dj,k+dk)
          end do
       end do
    end do
  end subroutine capture
  
end subroutine Refluxing_CaptureFluxes
