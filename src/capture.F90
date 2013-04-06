#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine Refluxing_Sync(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
end subroutine Refluxing_Sync



subroutine Refluxing_CaptureFluxes(CCTK_ARGUMENTS)
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
  if (delayed_refluxing/=0) then
     if (delayed_refluxing_sources==0) then
        ! Apply delayed correction to fluxes
        call correct_flux &
             (densflux, delayed_correction(:,:,:,3*index_dens+flux_direction))
        call correct_flux &
             (sxflux, delayed_correction(:,:,:,3*index_sx+flux_direction))
        call correct_flux &
             (syflux, delayed_correction(:,:,:,3*index_sy+flux_direction))
        call correct_flux &
             (szflux, delayed_correction(:,:,:,3*index_sz+flux_direction))
        call correct_flux &
             (tauflux, delayed_correction(:,:,:,3*index_tau+flux_direction))
        if (CCTK_EQUALS(Y_e_evolution_method, "GRHydro")) then
           call correct_flux &
                (Y_e_con_flux, &
                delayed_correction(:,:,:,3*index_ye+flux_direction))
        end if
        if (CCTK_EQUALS(Bvec_evolution_method, "GRHydro")) then
           call correct_flux &
                (Bconsxflux, &
                delayed_correction(:,:,:,3*index_Bconsx+flux_direction))
           call correct_flux &
                (Bconsyflux, &
                delayed_correction(:,:,:,3*index_Bconsy+flux_direction))
           call correct_flux &
                (Bconszflux, &
                delayed_correction(:,:,:,3*index_Bconsz+flux_direction))
        end if
     else
        ! Apply delayed correction to sources
        if (flux_direction == 1) then
           ! The really isn't a direction when correcting the sources
           call correct_source &
                (densrhs, delayed_correction(:,:,:,3*index_dens+1))
           call correct_source &
                (srhs(:,:,:,1), delayed_correction(:,:,:,3*index_sx+1))
           call correct_source &
                (srhs(:,:,:,2), delayed_correction(:,:,:,3*index_sy+1))
           call correct_source &
                (srhs(:,:,:,3), delayed_correction(:,:,:,3*index_sz+1))
           call correct_source &
                (taurhs, delayed_correction(:,:,:,3*index_tau+1))
           if (CCTK_EQUALS(Y_e_evolution_method, "GRHydro")) then
              call correct_source &
                   (Y_e_con_rhs, delayed_correction(:,:,:,3*index_ye+1))
           end if
           if (CCTK_EQUALS(Bvec_evolution_method, "GRHydro")) then
              call correct_source &
                   (Bconsrhs(:,:,:,1), &
                   delayed_correction(:,:,:,3*index_Bconsx+1))
              call correct_source &
                   (Bconsrhs(:,:,:,2), &
                   delayed_correction(:,:,:,3*index_Bconsy+1))
              call correct_source &
                   (Bconsrhs(:,:,:,3), &
                   delayed_correction(:,:,:,3*index_Bconsz+1))
           end if
        end if
     end if
  end if
  
  ! Capture fluxes from GRHydro
  call capture(flux(:,:,:,3*index_dens+flux_direction), densflux)
  call capture(flux(:,:,:,3*index_sx  +flux_direction), sxflux  )
  call capture(flux(:,:,:,3*index_sy  +flux_direction), syflux  )
  call capture(flux(:,:,:,3*index_sz  +flux_direction), szflux  )
  call capture(flux(:,:,:,3*index_tau +flux_direction), tauflux )
  if (CCTK_EQUALS(Y_e_evolution_method, "GRHydro")) then
     call capture(flux(:,:,:,3*index_ye+flux_direction), Y_e_con_flux)
  end if
  if (CCTK_EQUALS(Bvec_evolution_method, "GRHydro")) then
     call capture(flux(:,:,:,3*index_Bconsx+flux_direction), Bconsxflux)
     call capture(flux(:,:,:,3*index_Bconsy+flux_direction), Bconsyflux)
     call capture(flux(:,:,:,3*index_Bconsz+flux_direction), Bconszflux)
  end if
  
contains
  
  subroutine correct_flux(grhydro, refluxing)
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
  end subroutine correct_flux
  
  subroutine correct_source(grhydro, refluxing)
    CCTK_REAL, intent(inout) :: grhydro(:,:,:)
    CCTK_REAL, intent(in)    :: refluxing(:,:,:)
    CCTK_REAL :: factor
    integer   :: i,j,k
    
    ! Add fraction of Refluxing delayed correction variable to GRHydro
    ! source variable
    
    factor = delayed_refluxing_fraction / CCTK_DELTA_TIME
    
    !$omp parallel do private(i,j,k)
    do k=imin(3),imax(3)
       do j=imin(2),imax(2)
          do i=imin(1),imax(1)
             grhydro(i,j,k) = grhydro(i,j,k) + factor * refluxing(i,j,k)
          end do
       end do
    end do
  end subroutine correct_source
  
  subroutine capture(refluxing, grhydro)
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



subroutine Refluxing_DelayedCorrectionReduction(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  integer :: n
  integer :: imin(3), imax(3)

  ! Region with valid data
  imin(:) = 1           + cctk_nghostzones(:)
  imax(:) = cctk_lsh(:) - cctk_nghostzones(:)
  
  ! Reduce the amount of correction required after each time step
  
  do n=0,nvars-1
     call reduce_correction(delayed_correction(:,:,:,3*n+1))
  end do
  
contains
  
  subroutine reduce_correction(refluxing)
    CCTK_REAL, intent(inout) :: refluxing(:,:,:)
    CCTK_REAL :: factor
    integer   :: i,j,k
    
    factor = 1 - delayed_refluxing_fraction
    
    !$omp parallel do private(i,j,k)
    do k=imin(3),imax(3)
       do j=imin(2),imax(2)
          do i=imin(1),imax(1)
             refluxing(i,j,k) = factor * refluxing(i,j,k)
          end do
       end do
    end do
  end subroutine reduce_correction
  
end subroutine Refluxing_DelayedCorrectionReduction




subroutine Refluxing_ResetToAtmosphere(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_REAL :: rho_min
  integer :: i,j,k
  
  rho_min = GRHydro_rho_min*(1.0d0+GRHydro_atmo_tolerance)
  
  ! Loop over points and check if any point has been flagged for atmo
  ! reset or was reset to atmopshere.
  ! In that case, we set reflux_atmosphere_mask
  !$omp parallel do private(i,j,k)
  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
           
           if (atmosphere_mask(i,j,k) /= 0 .or. rho(i,j,k) <= rho_min) then
              reflux_atmosphere_mask(i,j,k) = 1.0d0
           end if
           
        end do
     end do
  end do
  
end subroutine Refluxing_ResetToAtmosphere
