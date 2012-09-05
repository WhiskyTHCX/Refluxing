#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

extern "C"
void Refluxing_DelayedCorrectionBoundaries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  bool const do_Ye   = CCTK_EQUALS(Y_e_evolution_method, "GRHydro");
  bool const do_Bvec = CCTK_EQUALS(Bvec_evolution_method, "GRHydro");
  
  CCTK_LOOP3_BND(delayed_correction, cctkGH, i,j,k, ni,nj,nk) {
    int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
    densflux_delayed_correction[ind3d] = 0.0;
    sxflux_delayed_correction[ind3d]   = 0.0;
    syflux_delayed_correction[ind3d]   = 0.0;
    szflux_delayed_correction[ind3d]   = 0.0;
    tauflux_delayed_correction[ind3d]  = 0.0;
    if (do_Ye) {
      yeflux_delayed_correction[ind3d] = 0.0;
    }
    if (do_Bvec) {
      Bconsxflux_delayed_correction[ind3d] = 0.0;
      Bconsyflux_delayed_correction[ind3d] = 0.0;
      Bconszflux_delayed_correction[ind3d] = 0.0;
    }
  } CCTK_ENDLOOP3STR_BND(delayed_correction);
}
