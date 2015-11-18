#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include <loopcontrol.h>

#include <carpet.hh>
#include <defs.hh>
#include <dh.hh>
#include <gh.hh>
#include <vect.hh>

#include <operator_prototypes_3d.hh>

#include "refluxing.hh"

using namespace std;
using namespace Carpet;
using namespace CarpetLib;
using namespace Refluxing;



static vector<CCTK_INT> saved_regridding_epochs;



#define SWITCH_TO_LEVEL(cctkGH, rl)             \
  do {                                          \
    bool switch_to_level_ = true;               \
    assert(is_level_mode());                    \
    const int rl_ = (rl);                       \
    BEGIN_GLOBAL_MODE (cctkGH) {                \
      ENTER_LEVEL_MODE (cctkGH, rl_) {
#define END_SWITCH_TO_LEVEL                     \
      } LEAVE_LEVEL_MODE;                       \
    } END_GLOBAL_MODE;                          \
    assert(switch_to_level_);                   \
    switch_to_level_ = false;                   \
  } while (false)



// Set weights on fine grid (for debugging only)
static
void flux_weight_fine_set(const cGH *restrict const cctkGH)
{
  DECLARE_CCTK_PARAMETERS;

  // Set weights on fine grid
  SWITCH_TO_LEVEL(cctkGH, reflevel+1) {

    if (veryverbose) {
      stringstream buf;
      buf << "Initialising fine grid on level " << reflevel << " to weight 0:";
      CCTK_INFO(buf.str().c_str());
    }

    BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
      BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        DECLARE_CCTK_ARGUMENTS;
        const dh& dd = *vdd.AT(Carpet::map);
        const dh::local_dboxes& local_box =
          dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component);



        // Set fine weight for restricted region

        // Initialise fine weight to zero
#pragma omp parallel
        CCTK_LOOP3_ALL(fine_restrict_init, cctkGH, i,j,k) {
          const int ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
          restrict_weight_fine[ind] = 0.0;
        } CCTK_ENDLOOP3_ALL(fine_restrict_init);

        // Set fine weight in restricted region to one
        const ibset& fine_restrict = local_box.restricted_region;
        LOOP_OVER_BSET(cctkGH, fine_restrict, box, imin, imax) {

          if (veryverbose) {
            stringstream buf;
            buf << "Setting fine grid restricted region on level " << reflevel << " map " << Carpet::map << " component " << component << " to weight 1: " << imin << ":" << imax-1;
            CCTK_INFO(buf.str().c_str());
          }

          assert(dim == 3);
#pragma omp parallel
          CCTK_LOOP3(fine_restrict,
                     i,j,k,
                     imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                     cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
          {
            const int ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
            restrict_weight_fine[ind] = 1.0;
          } CCTK_ENDLOOP3(fine_restrict);

        } END_LOOP_OVER_BSET;



        // Set fine weight for refluxing boundaries

        for (int dir=0; dir<3; ++dir) {
          // Unit vector
          const ivect idir = ivect::dir(dir);
          assert(dim == 3);

          // Initialise fine weight for refluxing boundaries to zero
#pragma omp parallel
          CCTK_LOOP3_ALL(fine_init, cctkGH, i,j,k) {
            const int ind = CCTK_VECTGFINDEX3D(cctkGH, i, j, k, dir);
            flux_weight_fine[ind] = 0.0;
          } CCTK_ENDLOOP3_ALL(fine_init);

          // Set fine weight for refluxing boundaries to one
          for (int face=0; face<2; ++face) {
            const ibset fine_boundary =
              local_box.fine_boundary[dir][face].shift(idir, 2);
            LOOP_OVER_BSET(cctkGH, fine_boundary, box, imin, imax) {

              if (veryverbose) {
                stringstream buf;
                buf << "Setting fine grid boundary on level " << reflevel << " map " << Carpet::map << " component " << component << " direction " << dir << " face " << face << " to weight 1: " << imin << ":" << imax-1;
                CCTK_INFO(buf.str().c_str());
              }

              assert(dim == 3);
#pragma omp parallel
              CCTK_LOOP3(fine_boundary,
                         i,j,k,
                         imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                         cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
              {
                int ind = CCTK_VECTGFINDEX3D(cctkGH, i, j, k, dir);
                flux_weight_fine[ind] = 1.0;
              } CCTK_ENDLOOP3(fine_boundary);

            } END_LOOP_OVER_BSET;
          } // for face

        } // for dir



      } END_LOCAL_COMPONENT_LOOP;
    } END_LOCAL_MAP_LOOP;
  } END_SWITCH_TO_LEVEL;
}



// Set weights on coarse grid (for debugging only)
static
void flux_weight_coarse_set(const cGH *restrict const cctkGH)
{
  DECLARE_CCTK_PARAMETERS;

  // Initialise coarse weight to zero everywhere
  if (veryverbose) {
    stringstream buf;
    buf << "Initialising coarse grid on level " << reflevel << " to weight 0:";
    CCTK_INFO(buf.str().c_str());
  }

  // Set weights on coarse grid
  BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      const dh& dd = *vdd.AT(Carpet::map);
      const dh::local_dboxes& local_box =
        dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component);



      // Set coarse weight for restricted region

      // Initialise coarse weight to zero
#pragma omp parallel
      CCTK_LOOP3_ALL(coarse_restrict_init, cctkGH, i,j,k) {
        int const ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
        restrict_weight_coarse[ind] = 0.0;
      } CCTK_ENDLOOP3_ALL(coarse_restrict_init);

      // Set coarse weight in restricted region to one
      const ibset& coarse_restrict = local_box.restricted_region;
      LOOP_OVER_BSET(cctkGH, coarse_restrict, box, imin, imax) {

        if (veryverbose) {
          stringstream buf;
          buf << "Setting coarse grid restricted region on level " << reflevel << " map " << Carpet::map << " component " << component << " to weight 1: " << imin << ":" << imax-1;
          CCTK_INFO(buf.str().c_str());
        }

        assert(dim == 3);
#pragma omp parallel
        CCTK_LOOP3(coarse_restrict,
                   i,j,k,
                   imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                   cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
        {
          int const ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
          restrict_weight_coarse[ind] = 1.0;
        } CCTK_ENDLOOP3(coarse_restrict);

      } END_LOOP_OVER_BSET;



      // Set coarse weight for refluxing boundaries

      for (int dir=0; dir<3; ++dir) {
        // Unit vector
        const ivect idir = ivect::dir(dir);
        assert(dim == 3);

        // Initialise coarse weight for refluxing boundaries to zero
#pragma omp parallel
        CCTK_LOOP3_ALL(coarse_init, cctkGH, i,j,k) {
          const int ind = CCTK_VECTGFINDEX3D(cctkGH, i, j, k, dir);
          flux_weight_coarse[ind] = 0.0;
        } CCTK_ENDLOOP3_ALL(coarse_init);

        // Set coarse weight for refluxing boundaries to one
        for (int face=0; face<2; ++face) {
          const ibset coarse_boundary =
            local_box.coarse_boundary[dir][face].shift(idir, 2);
          LOOP_OVER_BSET(cctkGH, coarse_boundary, box, imin, imax) {

            if (veryverbose) {
              stringstream buf;
              buf << "Setting coarse grid boundary on level " << reflevel << " map " << Carpet::map << " component " << component << " direction " << dir << " face " << face << " to weight 1: " << imin << ":" << imax-1;
              CCTK_INFO(buf.str().c_str());
            }

            assert(dim == 3);
#pragma omp parallel
            CCTK_LOOP3(coarse_boundary,
                       i,j,k,
                       imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                       cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
            {
              const int ind = CCTK_VECTGFINDEX3D(cctkGH, i, j, k, dir);
              flux_weight_coarse[ind] = 1.0;
            } CCTK_ENDLOOP3(coarse_boundary);

          } END_LOOP_OVER_BSET;
        } // for face

      } // for dir



    } END_LOCAL_COMPONENT_LOOP;
  } END_LOCAL_MAP_LOOP;
}



// Reflux
static
void reflux(const cGH *restrict const cctkGH)
{
  DECLARE_CCTK_PARAMETERS;

  assert(is_level_mode());



  if (veryverbose) {
    stringstream buf;
    buf << "Refluxing on level " << reflevel << ":";
    CCTK_INFO(buf.str().c_str());
  }

  // Check regridding epoch of current level
  const int levels = GetRefinementLevels(cctkGH);
  vector<CCTK_INT> epochs(levels);
  const int size = GetRegriddingEpochs(cctkGH, levels, &epochs[0]);
  assert(size == levels);

  const int level = GetRefinementLevel(cctkGH);
  assert(int(saved_regridding_epochs.size()) > level);
  if (saved_regridding_epochs.AT(level) != epochs.AT(level)) {
    CCTK_ERROR("Grid hierarchy changed while integrating fluxes for regridding");
  }


  // Initialise the coarse correction to zero everywhere
  BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;

      for (int n=0; n<refluxing_nvars; ++n) {
        for (int dir=0; dir<3; ++dir) {
          assert(dim == 3);
#pragma omp parallel
          CCTK_LOOP3_ALL(correction_coarse_init, cctkGH, i,j,k) {
            const int ind = CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 3*n+dir);
            correction[ind] = 0.0;
          } CCTK_ENDLOOP3_ALL(correction_coarse_init);
        } // for dir
      }   // for n


    } END_LOCAL_COMPONENT_LOOP;
  } END_LOCAL_MAP_LOOP;



  // Set the fine correction to the fine grid register
  SWITCH_TO_LEVEL(cctkGH, reflevel+1) {
    BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
      BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        DECLARE_CCTK_ARGUMENTS;

        for (int n=0; n<refluxing_nvars; ++n) {
          for (int dir=0; dir<3; ++dir) {

            assert(dim == 3);
#pragma omp parallel
            CCTK_LOOP3_ALL(correction_fine_init, cctkGH, i,j,k) {
              const int ind = CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 3*n+dir);
              correction[ind] = register_fine[ind];
              // assert(not isnan(flux_correction_ptr[ind+dir*np]));
            } CCTK_ENDLOOP3_ALL(correction_fine_init);
          } // for dir
        }   // for n


      } END_LOCAL_COMPONENT_LOOP;
    } END_LOCAL_MAP_LOOP;
  } END_SWITCH_TO_LEVEL;



  // Restrict the correction
  {
    const int gi = CCTK_GroupIndex("Refluxing::correction");
    const int tl = 0;

    // Restrict
    for (comm_state state; not state.done(); state.step()) {
      for (int m=0; m<maps; ++m) {
        // Fluxes
        for (int dir=0; dir<3; ++dir) {
          for (int face=0; face<2; ++face) {
            for (int n=0; n<refluxing_nvars; ++n) {
              ggf *const gv = arrdata.AT(gi).AT(m).data.AT(3*n+dir);
              gv->ref_reflux_all(state, tl, reflevel, mglevel, dir, face);
            }
          }
        } // for dir
      } // m
    } // for state
  }



  // Scale the restricted fine grid register (in the correction), and
  // subtract the coarse grid register; then update the state
  BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      const dh& dd = *vdd.AT(Carpet::map);
      const dh::local_dboxes& local_box =
        dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component);
      const ivect& lsh = ivect::ref(cctk_lsh);

      const rvect delta(CCTK_DELTA_SPACE(0),
                        CCTK_DELTA_SPACE(1),
                        CCTK_DELTA_SPACE(2));

      for (int dir=0; dir<3; ++dir) {
        // Unit vector
        const ivect idir = ivect::dir(dir);
        for (int face=0; face<2; ++face) {
          // Apply the correction to the cell to the left (one index
          // lower) if on the lower face, or the cell to the right
          // (same index) if on the upper face. The remainder of the
          // correction is applied to the other cell abutting this
          // face.
          const int idelta = index(lsh, idir);
          const int ioff = (face ? 0 : -1) * idelta;
          const int ifactor = face ? +1 : -1;
          const CCTK_REAL factor = ifactor;
          const ibset coarse_boundary =
            local_box.coarse_boundary[dir][face].shift(idir, 2);
          LOOP_OVER_BSET(cctkGH, coarse_boundary, box, imin, imax) {

            if (veryverbose) {
              stringstream buf;
              buf << "Refluxing on level " << reflevel << " map " << Carpet::map << " component " << component << " direction " << dir << " face " << face << ": " << imin << ":" << imax-1;
              CCTK_INFO(buf.str().c_str());
            }

            for (int n=0; n<refluxing_nvars; ++n) {
              const int vi = refluxing_vars.AT(n);
              CCTK_REAL *restrict const var =
                static_cast<CCTK_REAL*> (CCTK_VarDataPtrI(cctkGH, 0, vi));

              assert(dim == 3);
#pragma omp parallel
              CCTK_LOOP3(correction_calculate,
                         i,j,k,
                         imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                         cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
              {
                const int ind = CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 3*n+dir);
                const int state_ind = CCTK_GFINDEX3D(cctkGH, i, j, k);

                // Calculate flux difference
                correction[ind] -= register_coarse[ind];
                // assert(not isnan(correction[ind]));

                // Calculate correction
                CCTK_REAL difference = factor * correction[ind] / delta[dir];

                if (refluxing_debug_variables) {
                  // Keep a total of the refluxing changes
                  correction_total[ind+ioff] += difference;
                  // assert(not isnan(correction_total[ind]));
                }

                CCTK_REAL corr;
                corr = difference;
                var[state_ind+ioff] += corr;

              } CCTK_ENDLOOP3(correction_calculate);

            }

          } END_LOOP_OVER_BSET;
        } // for face
      } // for dir

    } END_LOCAL_COMPONENT_LOOP;
  } END_LOCAL_MAP_LOOP;
}



// Poison the boundary of the state vector to ensure boundary
// conditions are applied after refluxing
static
void poison_state_boundary(const cGH* restrict const cctkGH)
{
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL poison;
  memset(&poison, get_poison_value(), sizeof poison);

  BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;

      for (int n=0; n<refluxing_nvars; ++n) {
        const int vi = refluxing_vars.AT(n);
        CCTK_REAL *restrict const ptr =
          static_cast<CCTK_REAL*>(CCTK_VarDataPtrI(cctkGH, 0, vi));

        // Loop over all ghost faces
        for (int dir=0; dir<3; ++dir) {
          for (int face=0; face<2; ++face) {
            if (not cctk_bbox[2*dir+face]) {
              int imin[3], imax[3];
              for (int d=0; d<3; ++d) {
                imin[d] = 0;
                imax[d] = cctk_lsh[d];
              }
              if (face==0) {
                imax[dir] = imin[dir] + cctk_nghostzones[dir];
              } else {
                imin[dir] = imax[dir] - cctk_nghostzones[dir];
              }

              CCTK_LOOP3(poison_state_boundary,
                         i,j,k,
                         imin[0], imin[1], imin[2],
                         imax[0], imax[1], imax[2],
                         cctk_ash[0], cctk_ash[1], cctk_ash[2])
              {
                const int ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
                ptr[ind3d] = poison;
              } CCTK_ENDLOOP3(poison_state_boundary);

            } // if !bbox
          }   // face
        }     // dir

      } // n

    } END_LOCAL_COMPONENT_LOOP;
  } END_LOCAL_MAP_LOOP;
}



// Reset the fine grid flux register (on the next finer level)
static
void flux_register_fine_reset(const cGH *restrict const cctkGH)
{
  DECLARE_CCTK_PARAMETERS;

  SWITCH_TO_LEVEL(cctkGH, reflevel+1) {
    BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
      BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        DECLARE_CCTK_ARGUMENTS;

        for (int n=0; n<refluxing_nvars; ++n) {
          for (int dir=0; dir<3; ++dir) {
            assert(dim == 3);
#pragma omp parallel
            CCTK_LOOP3_ALL(fine_reset, cctkGH, i,j,k) {
              int const ind = CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 3*n+dir);
              register_fine[ind] = 0.0;
            } CCTK_ENDLOOP3_ALL(fine_reset);
          }
        }

      } END_LOCAL_COMPONENT_LOOP;
    } END_LOCAL_MAP_LOOP;
  } END_SWITCH_TO_LEVEL;
}



// Reset the coarse grid flux register (on this level)
static
void flux_register_coarse_reset(const cGH *restrict const cctkGH)
{
  DECLARE_CCTK_PARAMETERS;

  assert(is_level_mode());

  // Remember regridding epoch of current level
  int const levels = GetRefinementLevels(cctkGH);
  vector<CCTK_INT> epochs(levels);
  int const size = GetRegriddingEpochs(cctkGH, levels, &epochs[0]);
  assert(size == levels);

  int const level = GetRefinementLevel(cctkGH);
  if (int(saved_regridding_epochs.size()) <= level) {
    assert(int(saved_regridding_epochs.size()) == level);
    saved_regridding_epochs.resize(level+1);
  }
  saved_regridding_epochs.AT(level) = epochs.AT(level);



  BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;

      for (int n=0; n<refluxing_nvars; ++n) {
        for (int dir=0; dir<3; ++dir) {
          assert(dim == 3);
#pragma omp parallel
          CCTK_LOOP3_ALL(coarse_reset, cctkGH, i,j,k) {
            int const ind = CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 3*n+dir);
            register_coarse[ind] = 0.0;
          } CCTK_ENDLOOP3_ALL(coarse_reset);
        }
      }

    } END_LOCAL_COMPONENT_LOOP;
  } END_LOCAL_MAP_LOOP;
}



extern "C"
void Refluxing_CorrectState(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;

  // Sanity check
  assert(mglevel>=0 and reflevel>=0 and Carpet::map==-1 and component==-1);

  // We reflux on level L by taking a correction from level L+1 into
  // account. (This corresponds to restriction, where level L is
  // "corrected" from level L+1.)

  if (verbose or veryverbose) {
    CCTK_VInfo(CCTK_THORNSTRING,
               "Refluxing at iteration %d on level %d of %d",
               cctkGH->cctk_iteration, reflevel, reflevels);
  }

  // There is no finer level; do nothing
  if (reflevel+1 >= reflevels) return;



  // This works only with cell centred grids
  for (int m=0; m<maps; ++m) {
    gh const& hh = *vhh.AT(m);
    if (hh.refcent != cell_centered) {
      CCTK_ERROR("Refluxing requires cell-centred grids");
    }
  }



  if (refluxing_debug_variables) {
    flux_weight_fine_set(cctkGH);
    flux_weight_coarse_set(cctkGH);
  }

  reflux(cctkGH);

  poison_state_boundary(cctkGH);

  flux_register_fine_reset(cctkGH);
  flux_register_coarse_reset(cctkGH);
}



extern "C"
void Refluxing_Init2(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Remember regridding epochs of all levels
  const int levels = GetRefinementLevels(cctkGH);
  saved_regridding_epochs.resize(levels);
  const int size =
    GetRegriddingEpochs(cctkGH, levels, &saved_regridding_epochs[0]);
  assert(size == levels);
}



extern "C"
void Refluxing_Reset(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;

  // Sanity check
  assert(mglevel>=0 and reflevel>=0 and Carpet::map==-1 and component==-1);

  if (verbose or veryverbose) {
    CCTK_VInfo(CCTK_THORNSTRING,
               "Resetting refluxing information at iteration %d on level %d of %d",
               cctkGH->cctk_iteration, reflevel, reflevels);
  }

  assert(reflevel > 0);

  // This assumes that level L-1 is aligned
  const int do_every =
    ipow(mgfact, mglevel) *
    (maxtimereflevelfact / timereffacts.AT(reflevel - 1));
  if (not ((cctkGH->cctk_iteration - 1) % do_every == 0)) {
    CCTK_ERROR("Cannot regrid with refluxing when the parent levels are not aligned");
  }

#if 0
  // There is no finer level; do nothing
  if (reflevel+1 >= reflevels) return;
#endif

  // Initialise the coarse values twice, so that the coarse values are
  // also initialised on the finest level. (There should be a cleaner
  // way to do this.)
  for (int dr=-1; dr<=0; ++dr) {
  SWITCH_TO_LEVEL(cctkGH, reflevel+dr) {

    if (dr==-1) {
      flux_register_fine_reset(cctkGH);
    }
    flux_register_coarse_reset(cctkGH);

    if (refluxing_debug_variables) {
      BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
          DECLARE_CCTK_ARGUMENTS;

          for (int n=0; n<refluxing_nvars; ++n) {
            for (int dir=0; dir<3; ++dir) {
              assert(dim == 3);
#pragma omp parallel
              CCTK_LOOP3_ALL(correction_total_reset,
                             cctkGH, i,j,k)
              {
                const int ind = CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 3*n+dir);
                correction_total[ind] = 0.0;
              } CCTK_ENDLOOP3_ALL(correction_total_reset);
            }
          }

        } END_LOCAL_COMPONENT_LOOP;
      } END_LOCAL_MAP_LOOP;
    }

    BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
      BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        DECLARE_CCTK_ARGUMENTS;

        for (int n=0; n<refluxing_nvars; ++n) {
          for (int dir=0; dir<3; ++dir) {
            assert(dim == 3);
#pragma omp parallel
            CCTK_LOOP3_ALL(correction_reset, cctkGH, i,j,k) {
              const int ind = CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 3*n+dir);
              correction[ind] = -1.0;
            } CCTK_ENDLOOP3_ALL(correction_reset);
          }
        }

      } END_LOCAL_COMPONENT_LOOP;
    } END_LOCAL_MAP_LOOP;

    if (refluxing_debug_variables) {
      if (dr==-1) {
        flux_weight_fine_set(cctkGH);
      }
      flux_weight_coarse_set(cctkGH);
    }

  } END_SWITCH_TO_LEVEL;
  } // for dr
}
