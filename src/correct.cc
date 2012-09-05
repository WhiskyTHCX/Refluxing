#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <loopcontrol.h>

#include <carpet.hh>
#include <dh.hh>
#include <gh.hh>
#include <vect.hh>

#include <operator_prototypes_3d.hh>

using namespace std;
using namespace Carpet;
using namespace CarpetLib;
  
  
  
#define SWITCH_TO_LEVEL(cctkGH, rl)             \
  do {                                          \
    bool switch_to_level_ = true;               \
    assert (is_level_mode());                   \
    int const rl_ = (rl);                       \
    BEGIN_GLOBAL_MODE (cctkGH) {                \
      ENTER_LEVEL_MODE (cctkGH, rl_) {
#define END_SWITCH_TO_LEVEL                     \
      } LEAVE_LEVEL_MODE;                       \
    } END_GLOBAL_MODE;                          \
    assert (switch_to_level_);                  \
    switch_to_level_ = false;                   \
  } while (false)



// Find out whether a variable should be refluxed. We reflux those
// variables that have storage.
static
bool
reflux_var (cGH const * restrict const cctkGH,
            char const * restrict const name)
{
  int const vi = CCTK_VarIndex(name);
  assert (vi>=0);
  int const gi = CCTK_GroupIndexFromVarI (vi);
  assert (gi >= 0);
  int const istat = CCTK_QueryGroupStorageI (cctkGH, gi);
  assert (istat >= 0);
  return istat > 0;
}



static
CCTK_REAL *
get_varptr (cGH const * restrict const cctkGH,
            int const rl,
            char const * restrict const name)
{
  int const varindex = CCTK_VarIndex (name);
  assert (varindex >= 0);
  CCTK_REAL * restrict const ptr =
    static_cast<CCTK_REAL *>
    (VarDataPtrI(cctkGH, Carpet::map, rl, component, 0, varindex));
  assert (ptr);
  return ptr;
}

static
vector<CCTK_REAL *>
get_varptrs (cGH const * restrict const cctkGH,
             int const rl,
             char const * restrict const * restrict const names)
{
  int nvars = 0;
  for (int n = 0; names[n]; ++n) {
    if (reflux_var(cctkGH, names[n])) ++nvars;
  }
  vector<CCTK_REAL *> ptrs;
  ptrs.reserve(nvars);
  for (int n = 0; names[n]; ++n) {
    if (reflux_var(cctkGH, names[n])) {
      CCTK_REAL *const ptr = get_varptr (cctkGH, rl, names[n]);
      ptrs.push_back (ptr);
    }
  }
  return ptrs;
}



static
void
get_varind (cGH const * restrict const cctkGH,
            char const * restrict const name,
            int& gi, int& vi)
{
  int const varindex = CCTK_VarIndex (name);
  assert (varindex >= 0);
  gi = CCTK_GroupIndexFromVarI (varindex);
  assert (gi >= 0);
  int const varbase = CCTK_FirstVarIndexI (gi);
  assert (varbase >= 0);
  vi = varindex - varbase;
}

static
void
get_varinds (cGH const * restrict const cctkGH,
             char const * restrict const * restrict const names,
             vector<int>& gis, vector<int>& vis)
{
  int nvars = 0;
  for (int n = 0; names[n]; ++n) {
    if (reflux_var(cctkGH, names[n])) ++nvars;
  }
  assert (gis.empty());
  assert (vis.empty());
  gis.reserve(nvars);
  vis.reserve(nvars);
  for (int n = 0; names[n]; ++n) {
    if (reflux_var(cctkGH, names[n])) {
      int gi, vi;
      get_varind (cctkGH, names[n], gi, vi);
      gis.push_back (gi);
      vis.push_back (vi);
    }
  }
}



namespace variables {
    
  char const * restrict const flux_register_fine[] = {
    "Refluxing::densflux_register_fine[0]",
    "Refluxing::sxflux_register_fine[0]",
    "Refluxing::syflux_register_fine[0]",
    "Refluxing::szflux_register_fine[0]",
    "Refluxing::tauflux_register_fine[0]",
    "Refluxing::yeflux_register_fine[0]",
    "Refluxing::Bconsxflux_register_fine[0]",
    "Refluxing::Bconsyflux_register_fine[0]",
    "Refluxing::Bconszflux_register_fine[0]",
     NULL
  };
  
  char const * restrict const flux_register_coarse[] = {
    "Refluxing::densflux_register_coarse[0]",
    "Refluxing::sxflux_register_coarse[0]",
    "Refluxing::syflux_register_coarse[0]",
    "Refluxing::szflux_register_coarse[0]",
    "Refluxing::tauflux_register_coarse[0]",
    "Refluxing::yeflux_register_coarse[0]",
    "Refluxing::Bconsxflux_register_coarse[0]",
    "Refluxing::Bconsyflux_register_coarse[0]",
    "Refluxing::Bconszflux_register_coarse[0]",
    NULL
  };
  
  char const * restrict const flux_correction[] = {
    "Refluxing::densflux_correction[0]",
    "Refluxing::sxflux_correction[0]",
    "Refluxing::syflux_correction[0]",
    "Refluxing::szflux_correction[0]",
    "Refluxing::tauflux_correction[0]",
    "Refluxing::yeflux_correction[0]",
    "Refluxing::Bconsxflux_correction[0]",
    "Refluxing::Bconsyflux_correction[0]",
    "Refluxing::Bconszflux_correction[0]",
    NULL
  };
  
  char const * restrict const flux_delayed_correction[] = {
    "Refluxing::densflux_delayed_correction[0]",
    "Refluxing::sxflux_delayed_correction[0]",
    "Refluxing::syflux_delayed_correction[0]",
    "Refluxing::szflux_delayed_correction[0]",
    "Refluxing::tauflux_delayed_correction[0]",
    "Refluxing::yeflux_delayed_correction[0]",
    "Refluxing::Bconsxflux_delayed_correction[0]",
    "Refluxing::Bconsyflux_delayed_correction[0]",
    "Refluxing::Bconszflux_delayed_correction[0]",
    NULL
  };
  
  char const * restrict const correction_total[] = {
    "Refluxing::dens_correction_total[0]",
    "Refluxing::sx_correction_total[0]",
    "Refluxing::sy_correction_total[0]",
    "Refluxing::sz_correction_total[0]",
    "Refluxing::tau_correction_total[0]",
    "Refluxing::ye_correction_total[0]",
    "Refluxing::Bconsx_correction_total[0]",
    "Refluxing::Bconsy_correction_total[0]",
    "Refluxing::Bconsz_correction_total[0]",
    NULL
  };
  
  char const * restrict const var[] = {
    "GRHydro::dens",
    "GRHydro::scon[0]",
    "GRHydro::scon[1]",
    "GRHydro::scon[2]",
    "GRHydro::tau",
    "GRHydro::Y_e_con",
    "GRHydro::Bcons[0]",
    "GRHydro::Bcons[1]",
    "GRHydro::Bcons[2]",
    NULL
  };
  
  char const * restrict const none[] = {
    NULL
  };
  
} // namespace variables



// Set weights on fine grid (for debugging only)
static
void flux_weight_fine_set (cGH const * restrict const cctkGH)
{
  DECLARE_CCTK_PARAMETERS;
  
  assert (is_level_mode());
  
  // Set weights on fine grid
  SWITCH_TO_LEVEL (cctkGH, reflevel+1) {
    
    if (veryverbose) {
      stringstream buf;
      buf << "Initialising fine grid on level " << reflevel << " to weight 0:";
      CCTK_INFO (buf.str().c_str());
    }
    
    BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
      BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
        DECLARE_CCTK_ARGUMENTS;
        dh const& dd = *vdd.AT(Carpet::map);
        dh::local_dboxes const& local_box =
          dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component);
        
        ivect const& lsh = ivect::ref(cctk_lsh);
        int const np = prod(lsh);
        
        
        
        // Set fine weight for restricted region
        
        // Initialise fine weight to zero
#pragma omp parallel
        CCTK_LOOP3_ALL(GRHydro_Reflux_fine_restrict_init, cctkGH, i,j,k) {
          int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
          restrict_weight_fine[ind] = 0.0;
        } CCTK_ENDLOOP3_ALL(GRHydro_Reflux_fine_restrict_init);
        
        // Set fine weight in restricted region to one
        ibset const fine_restrict = local_box.restricted_region;
        LOOP_OVER_BSET(cctkGH, fine_restrict, box, imin, imax) {
          
          if (veryverbose) {
            stringstream buf;
            buf << "Setting fine grid restricted region on level " << reflevel << " map " << Carpet::map << " component " << component << " to weight 1: " << imin << ":" << imax-1;
            CCTK_INFO (buf.str().c_str());
          }
          
          assert (dim == 3);
#pragma omp parallel
          CCTK_LOOP3(GRHydro_Reflux_fine_restrict,
                     i,j,k,
                     imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                     cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
          {
            int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
            restrict_weight_fine[ind] = 1.0;
          } CCTK_ENDLOOP3(GRHydro_Reflux_fine_restrict);
          
        } END_LOOP_OVER_BSET;
        
        
        
        // Set fine weight for refluxing boundaries
        
        for (int dir=0; dir<3; ++dir) {
          // Unit vector
          ivect const idir = ivect::dir(dir);
          assert (dim == 3);
          
          // Initialise fine weight for refluxing boundaries to zero
#pragma omp parallel
          CCTK_LOOP3_ALL(GRHydro_Reflux_fine_init, cctkGH, i,j,k) {
            int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
            flux_weight_fine[ind+dir*np] = 0.0;
          } CCTK_ENDLOOP3_ALL(GRHydro_Reflux_fine_init);
          
          // Set fine weight for refluxing boundaries to one
          for (int face=0; face<2; ++face) {
            ibset const fine_boundary =
              local_box.fine_boundary[dir][face].shift(idir, 2);
            LOOP_OVER_BSET(cctkGH, fine_boundary, box, imin, imax) {
              
              if (veryverbose) {
                stringstream buf;
                buf << "Setting fine grid boundary on level " << reflevel << " map " << Carpet::map << " component " << component << " direction " << dir << " face " << face << " to weight 1: " << imin << ":" << imax-1;
                CCTK_INFO (buf.str().c_str());
              }
              
              assert (dim == 3);
#pragma omp parallel
              CCTK_LOOP3(GRHydro_Reflux_fine_boundary,
                         i,j,k,
                         imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                         cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
              {
                int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
                flux_weight_fine[ind+dir*np] = 1.0;
              } CCTK_ENDLOOP3(GRHydro_Reflux_fine_boundary);
              
            } END_LOOP_OVER_BSET;
          } // for face
          
        } // for dir
        
        
        
      } END_LOCAL_COMPONENT_LOOP;
    } END_LOCAL_MAP_LOOP;
  } END_SWITCH_TO_LEVEL;
}



// Set weights on coarse grid (for debugging only)
static
void flux_weight_coarse_set (cGH const * restrict const cctkGH)
{
  DECLARE_CCTK_PARAMETERS;
  
  assert (is_level_mode());
  
  // Initialise coarse weight to zero everywhere
  if (veryverbose) {
    stringstream buf;
    buf << "Initialising coarse grid on level " << reflevel << " to weight 0:";
    CCTK_INFO (buf.str().c_str());
  }
  
  // Set weights on coarse grid
  BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      dh const& dd = *vdd.AT(Carpet::map);
      dh::local_dboxes const& local_box =
        dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component);
      
      ivect const& lsh = ivect::ref(cctk_lsh);
      int const np = prod(lsh);
      
      
      
      // Set coarse weight for restricted region
      
      // Initialise coarse weight to zero
#pragma omp parallel
      CCTK_LOOP3_ALL(GRHydro_Reflux_coarse_restrict_init, cctkGH, i,j,k) {
        int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
        restrict_weight_coarse[ind] = 0.0;
      } CCTK_ENDLOOP3_ALL(GRHydro_Reflux_coarse_restrict_init);
      
      // Set coarse weight in restricted region to one
      ibset const coarse_restrict = local_box.restricted_region;
      LOOP_OVER_BSET(cctkGH, coarse_restrict, box, imin, imax) {
        
        if (veryverbose) {
          stringstream buf;
          buf << "Setting coarse grid restricted region on level " << reflevel << " map " << Carpet::map << " component " << component << " to weight 1: " << imin << ":" << imax-1;
          CCTK_INFO (buf.str().c_str());
        }
        
        assert (dim == 3);
#pragma omp parallel
        CCTK_LOOP3(GRHydro_Reflux_coarse_restrict,
                   i,j,k,
                   imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                   cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
        {
          int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
          restrict_weight_coarse[ind] = 1.0;
        } CCTK_ENDLOOP3(GRHydro_Reflux_coarse_restrict);
        
      } END_LOOP_OVER_BSET;
      
      
      
      // Set coarse weight for refluxing boundaries
      
      for (int dir=0; dir<3; ++dir) {
        // Unit vector
        ivect const idir = ivect::dir(dir);
        assert (dim == 3);
        
        // Initialise coarse weight for refluxing boundaries to zero
#pragma omp parallel
        CCTK_LOOP3_ALL(GRHydro_Reflux_coarse_init, cctkGH, i,j,k) {
          int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
          flux_weight_coarse[ind+dir*np] = 0.0;
        } CCTK_ENDLOOP3_ALL(GRHydro_Reflux_coarse_init);
        
        // Set coarse weight for refluxing boundaries to one
        for (int face=0; face<2; ++face) {
          ibset const coarse_boundary =
            local_box.coarse_boundary[dir][face].shift(idir, 2);
          LOOP_OVER_BSET(cctkGH, coarse_boundary, box, imin, imax) {
            
            if (veryverbose) {
              stringstream buf;
              buf << "Setting coarse grid boundary on level " << reflevel << " map " << Carpet::map << " component " << component << " direction " << dir << " face " << face << " to weight 1: " << imin << ":" << imax-1;
              CCTK_INFO (buf.str().c_str());
            }
            
            assert (dim == 3);
#pragma omp parallel
            CCTK_LOOP3(GRHydro_Reflux_coarse_boundary,
                       i,j,k,
                       imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                       cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
            {
              int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
              flux_weight_coarse[ind+dir*np] = 1.0;
            } CCTK_ENDLOOP3(GRHydro_Reflux_coarse_boundary);
            
          } END_LOOP_OVER_BSET;
        } // for face
        
      } // for dir
      
      
      
    } END_LOCAL_COMPONENT_LOOP;
  } END_LOCAL_MAP_LOOP;
}



// Reflux
static
void reflux (cGH const * restrict const cctkGH)
{
  DECLARE_CCTK_PARAMETERS;
  
  assert (is_level_mode());
  
  
  
  if (veryverbose) {
    stringstream buf;
    buf << "Refluxing on level " << reflevel << ":";
    CCTK_INFO (buf.str().c_str());
  }
  
  // Obtain atmosphere variable data
  bool const need_mask =
    suppress_refluxing_in_atmosphere or apply_limiter_atmo;
  
  int mask_fine_vi = -1;
  int mask_coarse_vi = -1, mask_coarse_gi = -1, mask_coarse_v0 = -1;
  
  if (need_mask) {
    mask_fine_vi = CCTK_VarIndex ("GRHydro::atmosphere_mask");
    assert (mask_fine_vi >= 0);
    mask_coarse_vi = CCTK_VarIndex ("Refluxing::restricted_atmosphere_mask");
    assert (mask_coarse_vi >= 0);
    mask_coarse_gi = CCTK_GroupIndexFromVarI (mask_coarse_vi);
    assert (mask_coarse_gi >= 0);
    mask_coarse_v0 = CCTK_FirstVarIndexI (mask_coarse_gi);
    assert (mask_coarse_v0 >= 0);
  }
  
  
  
  // Initialise the coarse correction to zero everywhere
  BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      ivect const& lsh = ivect::ref(cctk_lsh);
      int const np = prod(lsh);
      
      vector<CCTK_REAL *> flux_correction_ptrs =
        get_varptrs (cctkGH, reflevel, variables::flux_correction);
      int const nvars = flux_correction_ptrs.size();
      
      for (int n=0; n<nvars; ++n) {
        CCTK_REAL *restrict const flux_correction_ptr =
          flux_correction_ptrs.AT(n);
        assert (dim == 3);
#pragma omp parallel
        CCTK_LOOP3_ALL(GRHydro_Reflux_correction_coarse_init, cctkGH, i,j,k) {
          int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
          for (int dir=0; dir<3; ++dir) {
            flux_correction_ptr[ind+dir*np] = 0.0;
          }
        } CCTK_ENDLOOP3_ALL(GRHydro_Reflux_correction_coarse_init);
      } // for n
      
      if (need_mask) {
        // Initialise the atmosphere mask (in case it will not be set
        // by the restriction below)
        CCTK_REAL *restrict const mask_coarse_ptr =
          (CCTK_REAL *) CCTK_VarDataPtrI (cctkGH, 0, mask_coarse_vi);
        assert (mask_coarse_ptr);
        
#pragma omp parallel
        CCTK_LOOP3_ALL(GRHydro_Reflux_mask_init, cctkGH, i,j,k) {
          int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
          mask_coarse_ptr[ind] = 0.0;
        } CCTK_ENDLOOP3_ALL(GRHydro_Reflux_mask_init);
      }
      
    } END_LOCAL_COMPONENT_LOOP;
  } END_LOCAL_MAP_LOOP;
  
  
  
  // Set the fine correction to the fine grid register
  SWITCH_TO_LEVEL (cctkGH, reflevel+1) {
    BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
      BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
        DECLARE_CCTK_ARGUMENTS;
        ivect const& lsh = ivect::ref(cctk_lsh);
        int const np = prod(lsh);
        
        vector<CCTK_REAL *> flux_correction_ptrs =
          get_varptrs (cctkGH, reflevel, variables::flux_correction);
        vector<CCTK_REAL *> flux_register_fine_ptrs =
          get_varptrs (cctkGH, reflevel, variables::flux_register_fine);
        int const 
        nvars = flux_correction_ptrs.size();
                
        for (int n=0; n<nvars; ++n) {
          CCTK_REAL *restrict const flux_correction_ptr =
            flux_correction_ptrs.AT(n);
          CCTK_REAL const *restrict const flux_register_fine_ptr =
            flux_register_fine_ptrs.AT(n);
          assert (dim == 3);
          // TODO: loop only over a region?
#pragma omp parallel
          CCTK_LOOP3_ALL(GRHydro_Reflux_correction_fine_init, cctkGH, i,j,k) {
            int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
            for (int dir=0; dir<3; ++dir) {
              flux_correction_ptr[ind+dir*np] =
                flux_register_fine_ptr[ind+dir*np];
              // assert (not isnan(flux_correction_ptr[ind+dir*np]));
            }
          } CCTK_ENDLOOP3_ALL(GRHydro_Reflux_correction_fine_init);
        } // for n
        
        if (need_mask) {
          // Copy the atmosphere mask
          CCTK_INT const *restrict const mask_fine_ptr =
            (CCTK_INT const*) CCTK_VarDataPtrI (cctkGH, 0, mask_fine_vi);
          assert (mask_fine_ptr);
          CCTK_REAL *restrict const mask_coarse_ptr =
            (CCTK_REAL *) CCTK_VarDataPtrI (cctkGH, 0, mask_coarse_vi);
          assert (mask_coarse_ptr);
          // TODO: loop only over a region?
#pragma omp parallel
          CCTK_LOOP3_ALL(GRHydro_Reflux_mask_copy, cctkGH, i,j,k) {
            int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
            // Otherwise restricting may cancel atmosphere points:
            assert (mask_fine_ptr[ind] >= 0);
            mask_coarse_ptr[ind] = (CCTK_REAL) mask_fine_ptr[ind];
          } CCTK_ENDLOOP3_ALL(GRHydro_Reflux_mask_copy);
        }
        
      } END_LOCAL_COMPONENT_LOOP;
    } END_LOCAL_MAP_LOOP;
  } END_SWITCH_TO_LEVEL;
  
  
  
  // Restrict the correction
  {
    vector<int> gis, vis;
    get_varinds (cctkGH, variables::flux_correction, gis, vis);
    int const nvars = gis.size();
    int const tl = 0;
    
#if 0
    // Synchronise the fine correction, since the fine grid fluxes
    // have not been calculated on the ghost points, but they may be
    // required for restricting
    for (comm_state state; not state.done(); state.step()) {
      for (int m=0; m<maps; ++m) {
        for (int n=0; n<nvars; ++n) {
          for (int dir=0; dir<3; ++dir) {
            int const gi = gis.AT(n);
            int const vi = vis.AT(n) + dir;
            ggf *const gv = arrdata.AT(gi).AT(m).data.AT(vi);
            gv->sync_all (state, tl, reflevel+1, mglevel);
          }
        }
      }
    } // for state
#endif
    
    // Restrict
    for (comm_state state; not state.done(); state.step()) {
      for (int m=0; m<maps; ++m) {
        // Fluxes
        for (int n=0; n<nvars; ++n) {
          for (int dir=0; dir<3; ++dir) {
            for (int face=0; face<2; ++face) {
              int const gi = gis.AT(n);
              int const vi = vis.AT(n) + dir;
              ggf *const gv = arrdata.AT(gi).AT(m).data.AT(vi);
              gv->ref_reflux_all (state, tl, reflevel, mglevel, dir, face);
            }
          }
        }
        if (need_mask) {
          // Atmosphere mask
          ggf *const gv =
            arrdata.AT(mask_coarse_gi).AT(m).data.
            AT(mask_coarse_vi - mask_coarse_v0);
          gv->ref_restrict_all (state, tl, reflevel, mglevel);
        }
      }
    } // for state
  }
  
  
  
  // Scale the restricted fine grid register (in the correction), and
  // subtract the coarse grid register; then update the state
  BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      dh const& dd = *vdd.AT(Carpet::map);
      dh::local_dboxes const& local_box =
        dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component);
      ivect const& lsh = ivect::ref(cctk_lsh);
      int const np = prod(lsh);
      
      rvect const delta (CCTK_DELTA_SPACE(0),
                         CCTK_DELTA_SPACE(1),
                         CCTK_DELTA_SPACE(2));
      
      // Fluxes etc.
      vector<CCTK_REAL *> flux_correction_ptrs =
        get_varptrs (cctkGH, reflevel, variables::flux_correction);
      vector<CCTK_REAL *> flux_register_coarse_ptrs =
        get_varptrs (cctkGH, reflevel, variables::flux_register_coarse);
      vector<CCTK_REAL *> correction_total_ptrs =
        get_varptrs (cctkGH, reflevel, (refluxing_debug_variables ?
                                        variables::correction_total :
                                        variables::none));
      vector<CCTK_REAL *> var_ptrs =
        get_varptrs (cctkGH, reflevel, variables::var);
      vector<CCTK_REAL *> flux_delayed_correction_ptrs =
        get_varptrs (cctkGH, reflevel, (delayed_refluxing ?
                                        variables::flux_delayed_correction :
                                        variables::none));
      int const nvars = var_ptrs.size();
      
      // Atmosphere mask
      CCTK_INT const *restrict mask_fine_ptr = NULL;
      CCTK_REAL const *restrict mask_coarse_ptr = NULL;
      if (need_mask) {
        mask_fine_ptr =
          (CCTK_INT const*) CCTK_VarDataPtrI (cctkGH, 0, mask_fine_vi);
        assert (mask_fine_ptr);
        mask_coarse_ptr =
          (CCTK_REAL const*) CCTK_VarDataPtrI (cctkGH, 0, mask_coarse_vi);
        assert (mask_coarse_ptr);
      }
      
      for (int dir=0; dir<3; ++dir) {
        // Unit vector
        ivect const idir = ivect::dir(dir);
        for (int face=0; face<2; ++face) {
          // Apply the correction to the cell to the left (one index
          // lower) if on the lower face, or the cell to the right
          // (same index) if on the upper face. The remainder of the
          // correction is applied to the other cell abutting this
          // face.
          int const idelta = index (lsh, idir);
          int const ioff = (face ? 0 : -1) * idelta;
          int const ifactor = face ? +1 : -1;
          int const ioff_other = ioff - ifactor * idelta;
          CCTK_REAL const factor = ifactor;
          ibset const coarse_boundary =
            local_box.coarse_boundary[dir][face].shift(idir, 2);
          LOOP_OVER_BSET(cctkGH, coarse_boundary, box, imin, imax) {
            
            if (veryverbose) {
              stringstream buf;
              buf << "Refluxing on level " << reflevel << " map " << Carpet::map << " component " << component << " direction " << dir << " face " << face << ": " << imin << ":" << imax-1;
              CCTK_INFO (buf.str().c_str());
            }
            
            for (int n=0; n<nvars; ++n) {
              assert (dim == 3);
#pragma omp parallel
              CCTK_LOOP3(GRHydro_Reflux_correction_calculate,
                         i,j,k,
                         imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                         cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
              {
                int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
                
                // Check for atmosphere
                bool is_atmosphere = false;
                if (need_mask) {
                  is_atmosphere =
                    mask_fine_ptr  [ind-idelta] != 0 or
                    mask_fine_ptr  [ind       ] != 0 or
                    mask_coarse_ptr[ind-idelta] != 0.0 or
                    mask_coarse_ptr[ind       ] != 0.0;
                }
                
                // Calculate flux difference
                flux_correction_ptrs.AT(n)[ind+dir*np] -=
                  flux_register_coarse_ptrs.AT(n)[ind+dir*np];
                // assert (not isnan(flux_correction_ptrs.AT(n)[ind+dir*np]));
                
                // Calculate correction
                CCTK_REAL difference =
                  factor * flux_correction_ptrs.AT(n)[ind+dir*np] / delta[dir];
                if (suppress_refluxing_in_atmosphere and is_atmosphere) {
                  difference = 0.0;
                }
                if (apply_limiter or (apply_limiter_atmo and is_atmosphere)) {
                  CCTK_REAL const absval = fabs(var_ptrs.AT(n)[ind+ioff]);
                  CCTK_REAL const limfact =
                    apply_limiter_atmo and is_atmosphere
                    ? limiter_atmo_factor
                    : limiter_factor;
                  CCTK_REAL const maxabsdiff = absval * limfact;
                  if (fabs(difference) > maxabsdiff) {
                    difference = copysign(maxabsdiff, difference);
                  }
                }
                
                if (refluxing_debug_variables) {
                  // Keep a total of the refluxing changes
                  correction_total_ptrs.AT(n)[ind+dir*np+ioff] += difference;
                  // assert (not isnan(correction_total_ptrs.AT(n)[ind+dir*np+ioff]));
                }
                
                // Calculate what part of the difference to apply; the
                // remainder will be applied to the fine grid instead
                // TODO: Try different criteria
                // TODO: Do not use max for the momentum
                CCTK_REAL correction, remainder;
                if (not reflux_prolongate) {
                  correction = difference;
                } else {
                  correction = 0.5 * difference;
                  remainder = difference - correction;
                  flux_correction_ptrs.AT(n)[ind+dir*np] =
                    factor * remainder * delta[dir];
                }
                
                // Update the state
                if (not suppress_refluxing) {
                  if (not delayed_refluxing) {
                    var_ptrs.AT(n)[ind+ioff] += correction;
                    if (reflux_prolongate) {
                      var_ptrs.AT(n)[ind+ioff_other] += remainder;
                    }
                  } else {
                    assert(not reflux_prolongate); // not implemented
                    flux_delayed_correction_ptrs.AT(n)[ind+dir*np] +=
                      factor * correction * delta[dir];
                  }
                  // assert (not isnan(var_ptrs.AT(n)[ind+ioff]));
                }
                
              } CCTK_ENDLOOP3(GRHydro_Reflux_correction_calculate);
            }
            
          } END_LOOP_OVER_BSET;
        } // for face
      } // for dir
      
    } END_LOCAL_COMPONENT_LOOP;
  } END_LOCAL_MAP_LOOP;
  
  
  
  if (reflux_prolongate) {
  SWITCH_TO_LEVEL (cctkGH, reflevel+1) {
    
    // Prolongate the remainder of the correction back to the fine grid
    {
      vector<int> gis, vis;
      get_varinds (cctkGH, variables::flux_correction, gis, vis);
      int const nvars = gis.size();
      int const tl = 0;
      // Prolongate
      for (comm_state state; not state.done(); state.step()) {
        for (int m=0; m<maps; ++m) {
          // Fluxes
          for (int n=0; n<nvars; ++n) {
            for (int dir=0; dir<3; ++dir) {
              for (int face=0; face<2; ++face) {
                int const gi = gis.AT(n);
                int const vi = vis.AT(n) + dir;
                ggf *const gv = arrdata.AT(gi).AT(m).data.AT(vi);
                gv->ref_reflux_prolongate_all
                  (state, tl, reflevel, mglevel, dir, face);
              }
            }
          }
        }
      } // for state
    }
    
    
    
    // Updated the fine grid state according to the prolongated
    // remaining correction
    BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
      BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
        DECLARE_CCTK_ARGUMENTS;
        dh const& dd = *vdd.AT(Carpet::map);
        dh::local_dboxes const& local_box =
          dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component);
        ivect const& lsh = ivect::ref(cctk_lsh);
        int const np = prod(lsh);
        
        rvect const delta (CCTK_DELTA_SPACE(0),
                           CCTK_DELTA_SPACE(1),
                           CCTK_DELTA_SPACE(2));
        
        // Fluxes etc.
        vector<CCTK_REAL *> flux_correction_ptrs =
          get_varptrs (cctkGH, reflevel, variables::flux_correction);
        vector<CCTK_REAL *> var_ptrs =
          get_varptrs (cctkGH, reflevel, variables::var);
        int const nvars = var_ptrs.size();
        
        for (int dir=0; dir<3; ++dir) {
          // Unit vector
          ivect const idir = ivect::dir(dir);
          for (int face=0; face<2; ++face) {
            // Apply the correction to the cell to the left (one index
            // lower) if on the lower face, or the cell to the right
            // (same index) if on the upper face. The remainder of the
            // correction is applied to the other cell abutting this
            // face.
            int const idelta = index (lsh, idir);
            int const ioff = (face ? 0 : -1) * idelta;
            int const ifactor = face ? +1 : -1;
            int const ioff_other = ioff - ifactor * idelta;
            CCTK_REAL const factor = ifactor;
            ibset const fine_boundary =
              local_box.fine_boundary[dir][face].shift(idir, 2);
            LOOP_OVER_BSET(cctkGH, fine_boundary, box, imin, imax) {
              
              if (veryverbose) {
                stringstream buf;
                buf << "Applying remainder of correction on level " << reflevel << " map " << Carpet::map << " component " << component << " direction " << dir << " face " << face << ": " << imin << ":" << imax-1;
                CCTK_INFO (buf.str().c_str());
              }
              
              for (int n=0; n<nvars; ++n) {
                assert (dim == 3);
#pragma omp parallel
                CCTK_LOOP3(GRHydro_Reflux_remainder_calculate,
                           i,j,k,
                           imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                           cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
                {
                  int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
                  
                  // Calculate remaining correction
                  CCTK_REAL const remainder =
                    factor *
                    flux_correction_ptrs.AT(n)[ind+dir*np] / delta[dir];
                  
                  // Update the state
                  if (not suppress_refluxing) {
                    var_ptrs.AT(n)[ind+ioff_other] += remainder;
                  }
                  
                } CCTK_ENDLOOP3(GRHydro_Reflux_remainder_calculate);
              }
              
            } END_LOOP_OVER_BSET;
          } // for face
        } // for dir
        
      } END_LOCAL_COMPONENT_LOOP;
    } END_LOCAL_MAP_LOOP;
    
  } END_SWITCH_TO_LEVEL;
  } // if reflux_prolongate
  
}



// Reset the fine grid flux register (on the next finer level)
static
void flux_register_fine_reset (cGH const * restrict const cctkGH)
{
  assert (is_level_mode());
  
  SWITCH_TO_LEVEL (cctkGH, reflevel+1) {
    BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
      BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
        DECLARE_CCTK_ARGUMENTS;
        
        ivect const& lsh = ivect::ref(cctk_lsh);
        int const np = prod(lsh);
        
        vector<CCTK_REAL *> flux_register_fine_ptrs =
          get_varptrs (cctkGH, reflevel, variables::flux_register_fine);
        int const nvars = flux_register_fine_ptrs.size();
        
        for (int n=0; n<nvars; ++n) {
          for (int dir=0; dir<3; ++dir) {
            CCTK_REAL *restrict const flux_register_fine_ptr =
              flux_register_fine_ptrs.AT(n);
            assert (dim == 3);
#pragma omp parallel
            CCTK_LOOP3_ALL(GRHydro_Reflux_fine_reset, cctkGH, i,j,k) {
              int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
              flux_register_fine_ptr[ind+dir*np] = 0.0;
            } CCTK_ENDLOOP3_ALL(GRHydro_Reflux_fine_reset);
          }
        }
        
      } END_LOCAL_COMPONENT_LOOP;
    } END_LOCAL_MAP_LOOP;
  } END_SWITCH_TO_LEVEL;
}  



// Reset the coarse grid flux register (on this level)
static
void flux_register_coarse_reset (cGH const * restrict const cctkGH)
{
  assert (is_level_mode());
  
  BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      
      ivect const& lsh = ivect::ref(cctk_lsh);
      int const np = prod(lsh);
      
      vector<CCTK_REAL *> flux_register_coarse_ptrs =
        get_varptrs (cctkGH, reflevel, variables::flux_register_coarse);
      int const nvars = flux_register_coarse_ptrs.size();
      
      for (int n=0; n<nvars; ++n) {
        for (int dir=0; dir<3; ++dir) {
          CCTK_REAL *restrict const flux_register_coarse_ptr =
            flux_register_coarse_ptrs.AT(n);
          assert (dim == 3);
#pragma omp parallel
          CCTK_LOOP3_ALL(GRHydro_Reflux_coarse_reset, cctkGH, i,j,k) {
            int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
            flux_register_coarse_ptr[ind+dir*np] = 0.0;
          } CCTK_ENDLOOP3_ALL(GRHydro_Reflux_coarse_reset);
        }
      }
      
    } END_LOCAL_COMPONENT_LOOP;
  } END_LOCAL_MAP_LOOP;
}



extern "C"
void Refluxing_CorrectState (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  
  // Sanity check
  assert (mglevel>=0 and reflevel>=0 and Carpet::map==-1 and component==-1);
  
  // We reflux on level L by taking a correction from level L+1 into
  // account. (This corresponds to restriction, where level L is
  // "corrected" from level L+1.)
  
  if (verbose or veryverbose) {
    CCTK_VInfo (CCTK_THORNSTRING,
                "Refluxing at iteration %d on level %d of %d",
                cctkGH->cctk_iteration, reflevel, reflevels);
  }
  
  // There is no finer level; do nothing
  if (reflevel+1 >= reflevels) return;
  
  
  
  // This works only with cell centred grids
  for (int m=0; m<maps; ++m) {
    gh const& hh = *vhh.AT(m);
    if (hh.refcent != cell_centered) {
      CCTK_WARN (CCTK_WARN_ABORT, "Refluxing requires cell-centred grids");
    }
  }
  
  
  
  if (refluxing_debug_variables) {
    flux_weight_fine_set (cctkGH);
    flux_weight_coarse_set (cctkGH);
  }
  
  reflux (cctkGH);
  
  flux_register_fine_reset (cctkGH);
  flux_register_coarse_reset (cctkGH);
}



extern "C"
void Refluxing_Reset (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  
  // Sanity check
  assert (mglevel>=0 and reflevel>=0 and Carpet::map==-1 and component==-1);
  
  if (verbose or veryverbose) {
    CCTK_VInfo (CCTK_THORNSTRING,
                "Resetting refluxing information at iteration %d on level %d of %d",
                cctkGH->cctk_iteration, reflevel, reflevels);
  }
  
  assert (reflevel > 0);
  
  // This assumes that level L-1 is aligned
  int const do_every =
    ipow(mgfact, mglevel) *
    (maxtimereflevelfact / timereffacts.AT(reflevel - 1));
  if (not ((cctkGH->cctk_iteration - 1) % do_every == 0)) {
    cout << "iteration=" << cctkGH->cctk_iteration << "\n"
         << "do_every=" << do_every << "\n";
    CCTK_WARN (CCTK_WARN_ABORT,
               "Cannot regrid with refluxing when the parent levels are not aligned");
  }
  
#if 0
  // There is no finer level; do nothing
  if (reflevel+1 >= reflevels) return;
#endif
  
  // Initialise the coarse values twice, so that the coarse values are
  // also initialised on the finest level. (There should be a cleaner
  // way to do this.)
  for (int dr=-1; dr<=0; ++dr) {
  SWITCH_TO_LEVEL (cctkGH, reflevel+dr) {
    
    if (dr==-1) {
      flux_register_fine_reset (cctkGH);
    }
    flux_register_coarse_reset (cctkGH);
    
    if (refluxing_debug_variables) {
      BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
          DECLARE_CCTK_ARGUMENTS;
          
          ivect const& lsh = ivect::ref(cctk_lsh);
          int const np = prod(lsh);
          
          vector<CCTK_REAL *> flux_correction_total_ptrs =
            get_varptrs (cctkGH, reflevel, variables::correction_total);
          int const nvars = flux_correction_total_ptrs.size();
          
          for (int n=0; n<nvars; ++n) {
            CCTK_REAL *restrict const flux_correction_total_ptr =
              flux_correction_total_ptrs.AT(n);
            assert (dim == 3);
#pragma omp parallel
            CCTK_LOOP3_ALL(GRHydro_Reflux_correction_total_reset, cctkGH, i,j,k) {
              int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
              for (int dir=0; dir<3; ++dir) {
                flux_correction_total_ptr[ind+dir*np] = 0.0;
              }
            } CCTK_ENDLOOP3_ALL(GRHydro_Reflux_correction_total_reset);
          }
          
        } END_LOCAL_COMPONENT_LOOP;
      } END_LOCAL_MAP_LOOP;
    }
    
    BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
      BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
        DECLARE_CCTK_ARGUMENTS;
        
        ivect const& lsh = ivect::ref(cctk_lsh);
        int const np = prod(lsh);
        
        vector<CCTK_REAL *> flux_correction_ptrs =
          get_varptrs (cctkGH, reflevel, variables::flux_correction);
        int const nvars = flux_correction_ptrs.size();
        
        for (int n=0; n<nvars; ++n) {
          CCTK_REAL *restrict const flux_correction_ptr =
            flux_correction_ptrs.AT(n);
          assert (dim == 3);
#pragma omp parallel
          CCTK_LOOP3_ALL(GRHydro_Reflux_correction_reset, cctkGH, i,j,k) {
            int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
            for (int dir=0; dir<3; ++dir) {
              flux_correction_ptr[ind+dir*np] = -1.0;
            }
          } CCTK_ENDLOOP3_ALL(GRHydro_Reflux_correction_reset);
        }
        
      } END_LOCAL_COMPONENT_LOOP;
    } END_LOCAL_MAP_LOOP;
    
    if (refluxing_debug_variables) {
      if (dr==-1) {
        flux_weight_fine_set (cctkGH);
      }
      flux_weight_coarse_set (cctkGH);
    }
    
  } END_SWITCH_TO_LEVEL;
  } // for dr
}
