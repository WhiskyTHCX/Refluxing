#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

extern "C"
void Refluxing_DelayedCorrectionBoundaries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_LOOP3_BND(delayed_correction, cctkGH, i,j,k, ni,nj,nk) {
    for (int n=0; n<nvars; ++n) {
      for (int d=0; d<3; ++d) {
        int const ind = CCTK_VECTGFINDEX3D(cctkGH, i,j,k, 3*n+d);
        delayed_correction[ind] = 0.0;
      }
    }
  } CCTK_ENDLOOP3STR_BND(delayed_correction);
}
