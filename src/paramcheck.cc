#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>



extern "C"
void Refluxing_ParamCheck (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_IsThornActive ("CarpetRegrid2")) {
    CCTK_INT const snap_to_coarse =
      * static_cast<CCTK_INT const *>
      (CCTK_ParameterGet ("snap_to_coarse", "CarpetRegrid2", 0));
    if (not snap_to_coarse) {
      // Setting this parameter is not really required, but it can't
      // hurt
      CCTK_PARAMWARN ("Refluxing requires that the fine grid boundaries are aligned with coarse grid boundaries. Set CarpetRegrid2::snap_to_coarse to ensure this.");
    }
  }
}
