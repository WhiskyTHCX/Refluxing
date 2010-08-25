#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cstdlib>
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
  
  
  
#define SWITCH_TO_LEVEL(cctkGH, rl)                     \
  do {                                                  \
    bool switch_to_level_ = true;                       \
    assert (is_singlemap_mode());                       \
    int const rl_ = (rl);                               \
    int const m_ = Carpet::map;                         \
    BEGIN_GLOBAL_MODE (cctkGH) {                        \
      ENTER_LEVEL_MODE (cctkGH, rl_) {                  \
        ENTER_SINGLEMAP_MODE (cctkGH, m_, CCTK_GF) {
#define END_SWITCH_TO_LEVEL                     \
        } LEAVE_SINGLEMAP_MODE;                 \
      } LEAVE_LEVEL_MODE;                       \
    } END_GLOBAL_MODE;                          \
    assert (switch_to_level_);                  \
    switch_to_level_ = false;                   \
  } while (false)



static
CCTK_REAL * restrict
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
  while (names[nvars]) ++nvars;
  vector<CCTK_REAL *> ptrs(nvars);
  for (int n=0; n<nvars; ++n) {
    ptrs.AT(n) = get_varptr (cctkGH, rl, names[n]);
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
  while (names[nvars]) ++nvars;
  gis.resize(nvars);
  vis.resize(nvars);
  for (int n=0; n<nvars; ++n) {
    get_varind (cctkGH, names[n], gis.AT(n), vis.AT(n));
  }
}



namespace variables {
    
  char const * restrict const flux_register_fine[] = {
    "Refluxing::densflux_register_fine[0]",
    "Refluxing::sxflux_register_fine[0]",
    "Refluxing::syflux_register_fine[0]",
    "Refluxing::szflux_register_fine[0]",
    "Refluxing::tauflux_register_fine[0]",
    NULL
  };
  
  char const * restrict const flux_register_coarse[] = {
    "Refluxing::densflux_register_coarse[0]",
    "Refluxing::sxflux_register_coarse[0]",
    "Refluxing::syflux_register_coarse[0]",
    "Refluxing::szflux_register_coarse[0]",
    "Refluxing::tauflux_register_coarse[0]",
    NULL
  };
  
  char const * restrict const flux_correction[] = {
    "Refluxing::densflux_correction[0]",
    "Refluxing::sxflux_correction[0]",
    "Refluxing::syflux_correction[0]",
    "Refluxing::szflux_correction[0]",
    "Refluxing::tauflux_correction[0]",
    NULL
  };
  
  char const * restrict const correction_total[] = {
    "Refluxing::dens_correction_total[0]",
    "Refluxing::sx_correction_total[0]",
    "Refluxing::sy_correction_total[0]",
    "Refluxing::sz_correction_total[0]",
    "Refluxing::tau_correction_total[0]",
    NULL
  };
  
  char const * restrict const var[] = {
    "GRHydro::dens",
    "GRHydro::scon[0]",
    "GRHydro::scon[1]",
    "GRHydro::scon[2]",
    "GRHydro::tau",
    NULL
  };
  
} // namespace variables



// Set weights on fine grid (for debugging only)
static
void flux_weight_fine_set (cGH const * restrict const cctkGH)
{
  assert (is_singlemap_mode());
  
  // Set weights on fine grid
  SWITCH_TO_LEVEL (cctkGH, reflevel+1) {
    BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      dh const& dd = *vdd.AT(Carpet::map);
      dh::local_dboxes const& local_box =
        dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component);
      
      ivect const& lsh = ivect::ref(cctk_lsh);
      int const np = prod(lsh);
      
      for (int dir=0; dir<3; ++dir) {
        
        // Initialise fine weight to zero everywhere
        cout << "Initialising fine grid on level " << reflevel << " to weight 0:\n";
        assert (dim == 3);
#pragma omp parallel
        LC_LOOP3(GRHydro_Reflux_fine_init,
                 i,j,k,
                 0,0,0, cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
                 cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
        {
          int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
          flux_weight_fine[ind+dir*np] = 0.0;
        } LC_ENDLOOP3(GRHydro_Reflux_fine_init);
        
        // Set fine weight to one on boundary
        for (int face=0; face<2; ++face) {
          ibset const& fine_boundary = local_box.fine_boundary[dir][face];
          LOOP_OVER_BSET(cctkGH, fine_boundary, box, imin, imax) {
            
            cout << "Setting fine grid boundary on level " << reflevel << " direction " << dir << " face " << face << " to weight 1: " << imin << ":" << imax-1 << "\n";
            assert (dim == 3);
#pragma omp parallel
            LC_LOOP3(GRHydro_Reflux_fine_boundary,
                     i,j,k,
                     imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                     cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
            {
              int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
              flux_weight_fine[ind+dir*np] = 1.0;
            } LC_ENDLOOP3(GRHydro_Reflux_fine_boundary);
            
          } END_LOOP_OVER_BSET;
        } // for face
        
      } // for dir
      
    } END_LOCAL_COMPONENT_LOOP;
  } END_SWITCH_TO_LEVEL;
}



// Set weights on coarse grid (for debugging only)
static
void flux_weight_coarse_set (cGH const * restrict const cctkGH)
{
  assert (is_singlemap_mode());
  
  // Set weights on coarse grid
  BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
    DECLARE_CCTK_ARGUMENTS;
    dh const& dd = *vdd.AT(Carpet::map);
    dh::local_dboxes const& local_box =
      dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component);
    
    ivect const& lsh = ivect::ref(cctk_lsh);
    int const np = prod(lsh);
    
    for (int dir=0; dir<3; ++dir) {
      
      // Initialise coarse weight to zero everywhere
      cout << "Initialising coarse grid on level " << reflevel << " to weight 0:\n";
      assert (dim == 3);
#pragma omp parallel
      LC_LOOP3(GRHydro_Reflux_coarse_init,
               i,j,k,
               0,0,0, cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
               cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
      {
        int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
        flux_weight_coarse[ind+dir*np] = 0.0;
      } LC_ENDLOOP3(GRHydro_Reflux_coarse_init);
      
      // Set coarse weight to one on boundary
      for (int face=0; face<2; ++face) {
        ibset const& coarse_boundary = local_box.coarse_boundary[dir][face];
        LOOP_OVER_BSET(cctkGH, coarse_boundary, box, imin, imax) {
          
          cout << "Setting coarse grid boundary on level " << reflevel << " direction " << dir << " face " << face << " to weight 1: " << imin << ":" << imax-1 << "\n";
          // Set weight on coarse grid boundary to one
          assert (dim == 3);
#pragma omp parallel
          LC_LOOP3(GRHydro_Reflux_coarse_boundary,
                   i,j,k,
                   imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                   cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
          {
            int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
            flux_weight_coarse[ind+dir*np] = 1.0;
          } LC_ENDLOOP3(GRHydro_Reflux_coarse_boundary);
          
        } END_LOOP_OVER_BSET;
      } // for face
      
    } // for dir
    
  } END_LOCAL_COMPONENT_LOOP;
}



// Reflux
static
void reflux_old (cGH const * restrict const cctkGH)
{
  DECLARE_CCTK_PARAMETERS;
  
  assert (is_singlemap_mode());
  
  BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
    DECLARE_CCTK_ARGUMENTS;
    gh const& hh = *vhh.AT(Carpet::map);
    dh const& dd = *vdd.AT(Carpet::map);
    dh::light_dboxes const& light_box =
      dd.light_boxes.AT(mglevel).AT(reflevel).AT(component);
    dh::local_dboxes const& local_box =
      dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component);
    
    dh::light_dboxes const& fine_light_box =
      dd.light_boxes.AT(mglevel).AT(reflevel+1).AT(component);
    
    ibbox const& ext = light_box.exterior;
    ivect const& lsh = ivect::ref(cctk_lsh);
    int const np = prod(lsh);
    
    rvect const delta (CCTK_DELTA_SPACE(0),
                       CCTK_DELTA_SPACE(1),
                       CCTK_DELTA_SPACE(2));
    
    ivect const reffact = hh.reffacts.AT(reflevel+1) / hh.reffacts.AT(reflevel);
    
    ibbox const& fine_ext = fine_light_box.exterior;
    ivect const fine_lsh = fine_ext.shape() / fine_ext.stride();
    int const fine_np = prod(fine_lsh);
    
    // Obtain pointer to fine grid flux register
    // NOTE: This assumes that the fine grid has the same component
    // structure as the coarse grid
    vector<CCTK_REAL *> flux_register_fine_ptrs =
      get_varptrs (cctkGH, reflevel+1, variables::flux_register_fine);
    vector<CCTK_REAL *> flux_register_coarse_ptrs =
      get_varptrs (cctkGH, reflevel, variables::flux_register_coarse);
    vector<CCTK_REAL *> flux_correction_ptrs =
      get_varptrs (cctkGH, reflevel, variables::flux_correction);
    vector<CCTK_REAL *> correction_total_ptrs =
      get_varptrs (cctkGH, reflevel, variables::correction_total);
    vector<CCTK_REAL *> var_ptrs =
      get_varptrs (cctkGH, reflevel, variables::var);
    int const nvars = var_ptrs.size();
    
    for (int dir=0; dir<3; ++dir) {
      // Unit vector
      ivect const idir = ivect::dir(dir);
      
      // Modify the extents, since the flux grid function is not
      // really staggered, but the restriction operator expects a
      // truly staggered grid function
      assert (all (ext.stride() % 2 == 0));
      ibbox const mext (ext.lower() + idir * ext.stride()/2,
                        ext.upper() + idir * ext.stride()/2,
                        ext.stride());
      assert (all (fine_ext.stride() % 2 == 0));
      ibbox const fine_mext (fine_ext.lower() + idir * fine_ext.stride()/2,
                             fine_ext.upper() + idir * fine_ext.stride()/2,
                             fine_ext.stride());
      
      void (*restrict_ptr) (CCTK_REAL const * restrict const src,
                            ivect3 const & restrict srcext,
                            CCTK_REAL * restrict const dst,
                            ivect3 const & restrict dstext,
                            ibbox3 const & restrict srcbbox,
                            ibbox3 const & restrict dstbbox,
                            ibbox3 const & restrict regbbox);
      switch (dir) {
      case 0: restrict_ptr = restrict_3d_vc_rf2<CCTK_REAL,0,1,1>; break;
      case 1: restrict_ptr = restrict_3d_vc_rf2<CCTK_REAL,1,0,1>; break;
      case 2: restrict_ptr = restrict_3d_vc_rf2<CCTK_REAL,1,1,0>; break;
      default: assert(0);
      }
      
      // Initialise correction to zero everywhere
      for (int n=0; n<nvars; ++n) {
        assert (dim == 3);
#pragma omp parallel
        LC_LOOP3(GRHydro_Reflux_correction_init,
                 i,j,k,
                 0,0,0, cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
                 cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
        {
          int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
          flux_correction_ptrs.AT(n)[ind+dir*np] = 0.0;
        } LC_ENDLOOP3(GRHydro_Reflux_correction_init);
      }
      
      // Restrict the fine grid register into the correction
      for (int face=0; face<2; ++face) {
        ibset const& coarse_boundary = local_box.coarse_boundary[dir][face];
        LOOP_OVER_BSET(cctkGH, coarse_boundary, box, imin, imax) {
          
          cout << "Refluxing on level " << reflevel << " direction " << dir << " face " << face << ": " << imin << ":" << imax-1 << "\n";
          
#warning "TODO: stagger fluxes according to face here, and don't stagger regions according to face in dh.cc"          
          
          // Modify the box
          assert (all (box.stride() % 2 == 0));
          ibbox const mbox (box.lower() + idir * box.stride()/2,
                            box.upper() + idir * box.stride()/2,
                            box.stride());
          
          // Restrict
          for (int n=0; n<nvars; ++n) {
            (*restrict_ptr)
              (flux_register_fine_ptrs.AT(n) + dir*fine_np,
               fine_lsh,
               flux_correction_ptrs.AT(n) + dir*np,
               lsh,
               fine_mext,
               mext,
               mbox);
          }
          
          // Scale the fine grid register (in the correction), and
          // subtract the coarse grid register; then update the
          // density
          int const ioff = index (lsh, (face ? +1 : 0) * idir);
          CCTK_REAL const factor = face ? +1 : -1;
          for (int n=0; n<nvars; ++n) {
            assert (dim == 3);
#pragma omp parallel
            LC_LOOP3(GRHydro_Reflux_correction_calculate,
                     i,j,k,
                     imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                     cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
            {
              int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
              flux_correction_ptrs.AT(n)[ind+dir*np] -=
                flux_register_coarse_ptrs.AT(n)[ind+dir*np];
              // Update the density
              CCTK_REAL const difference =
                factor * flux_correction_ptrs.AT(n)[ind+dir*np] / delta[dir];
              // Choose the cell to the left if on the lower face, or
              // the cell to the right if on the upper face
              if (refluxing_debug_variables) {
                // Keep a total of the refluxing changes
                correction_total_ptrs.AT(n)[ind+dir*np+ioff] += difference;
              }
              var_ptrs.AT(n)[ind+ioff] += difference;
            } LC_ENDLOOP3(GRHydro_Reflux_correction_calculate);
          }
          
        } END_LOOP_OVER_BSET;
      } // for face
    
    } // for dir
    
  } END_LOCAL_COMPONENT_LOOP;
}



// Reflux
static
void reflux (cGH const * restrict const cctkGH)
{
  DECLARE_CCTK_PARAMETERS;
  
  assert (is_singlemap_mode());
  
  
  
  // Initialise the coarse correction to zero everywhere
  BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
    DECLARE_CCTK_ARGUMENTS;
    ivect const& lsh = ivect::ref(cctk_lsh);
    int const np = prod(lsh);

    vector<CCTK_REAL *> flux_correction_ptrs =
      get_varptrs (cctkGH, reflevel, variables::flux_correction);
    int const nvars = flux_correction_ptrs.size();
    
    for (int n=0; n<nvars; ++n) {
      assert (dim == 3);
#pragma omp parallel
      LC_LOOP3(GRHydro_Reflux_correction_coarse_init,
               i,j,k,
               0,0,0, cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
               cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
      {
        int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
        for (int dir=0; dir<3; ++dir) {
          flux_correction_ptrs.AT(n)[ind+dir*np] = 0.0;
        }
      } LC_ENDLOOP3(GRHydro_Reflux_correction_coarse_init);
    } // for n
  } END_LOCAL_COMPONENT_LOOP;
  
  
  
  // Set the fine correction to the fine grid register
  SWITCH_TO_LEVEL (cctkGH, reflevel+1) {
    BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      ivect const& lsh = ivect::ref(cctk_lsh);
      int const np = prod(lsh);
      
      vector<CCTK_REAL *> flux_correction_ptrs =
        get_varptrs (cctkGH, reflevel, variables::flux_correction);
      vector<CCTK_REAL *> flux_register_fine_ptrs =
        get_varptrs (cctkGH, reflevel, variables::flux_register_fine);
      int const nvars = flux_correction_ptrs.size();
      
      for (int n=0; n<nvars; ++n) {
        assert (dim == 3);
#pragma omp parallel
        LC_LOOP3(GRHydro_Reflux_correction_fine_init,
                 i,j,k,
                 0,0,0, cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
                 cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
        {
          int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
          for (int dir=0; dir<3; ++dir) {
            flux_correction_ptrs.AT(n)[ind+dir*np] =
              flux_register_fine_ptrs.AT(n)[ind+dir*np];
          }
        } LC_ENDLOOP3(GRHydro_Reflux_correction_fine_init);
      } // for n
    } END_LOCAL_COMPONENT_LOOP;
  } END_SWITCH_TO_LEVEL;
  
  
  
  // Restrict the correction
  {
    vector<int> gis, vis;
    get_varinds (cctkGH, variables::flux_correction, gis, vis);
    int const nvars = gis.size();
    int const tl = 0;
    for (comm_state state; not state.done(); state.step()) {
      for (int n=0; n<nvars; ++n) {
        for (int dir=0; dir<3; ++dir) {
          for (int face=0; face<2; ++face) {
            int const gi = gis.AT(n);
            int const vi = vis.AT(n) + dir;
            ggf *const gv = arrdata.AT(gi).AT(Carpet::map).data.AT(vi);
            gv->ref_reflux_all (state, tl, reflevel, mglevel, dir, face);
          }
        }
      }
    } // for state
  }
  
  
  
  // Scale the restricted fine grid register (in the correction), and
  // subtract the coarse grid register; then update the state
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
    
    vector<CCTK_REAL *> flux_correction_ptrs =
      get_varptrs (cctkGH, reflevel, variables::flux_correction);
    vector<CCTK_REAL *> flux_register_coarse_ptrs =
      get_varptrs (cctkGH, reflevel, variables::flux_register_coarse);
    vector<CCTK_REAL *> correction_total_ptrs =
      get_varptrs (cctkGH, reflevel, variables::correction_total);
    vector<CCTK_REAL *> var_ptrs =
      get_varptrs (cctkGH, reflevel, variables::var);
    int const nvars = var_ptrs.size();
    
    for (int dir=0; dir<3; ++dir) {
      // Unit vector
      ivect const idir = ivect::dir(dir);
      for (int face=0; face<2; ++face) {
        int const ioff = index (lsh, (face ? +1 : 0) * idir);
        CCTK_REAL const factor = face ? +1 : -1;
        ibset const& coarse_boundary = local_box.coarse_boundary[dir][face];
        LOOP_OVER_BSET(cctkGH, coarse_boundary, box, imin, imax) {
          
          cout << "Refluxing on level " << reflevel << " direction " << dir << " face " << face << ": " << imin << ":" << imax-1 << "\n";
          
          for (int n=0; n<nvars; ++n) {
            assert (dim == 3);
#pragma omp parallel
            LC_LOOP3(GRHydro_Reflux_correction_calculate,
                   i,j,k,
                   imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
                   cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
            {
              int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
              flux_correction_ptrs.AT(n)[ind+dir*np] -=
                flux_register_coarse_ptrs.AT(n)[ind+dir*np];
              // Update the state
              CCTK_REAL const difference =
                factor * flux_correction_ptrs.AT(n)[ind+dir*np] / delta[dir];
              // Choose the cell to the left if on the lower face, or
              // the cell to the right if on the upper face
              if (refluxing_debug_variables) {
                // Keep a total of the refluxing changes
                correction_total_ptrs.AT(n)[ind+dir*np+ioff] += difference;
              }
              var_ptrs.AT(n)[ind+ioff] += difference;
            } LC_ENDLOOP3(GRHydro_Reflux_correction_calculate);
          }
          
        } END_LOOP_OVER_BSET;
      } // for face
    } // for dir
    
  } END_LOCAL_COMPONENT_LOOP;
}



// Reset the fine grid flux register (on the next finer level)
static
void flux_register_fine_reset (cGH const * restrict const cctkGH)
{
  assert (is_singlemap_mode());
  
  SWITCH_TO_LEVEL (cctkGH, reflevel+1) {
    BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      
      ivect const& lsh = ivect::ref(cctk_lsh);
      int const np = prod(lsh);
      
      vector<CCTK_REAL *> flux_register_fine_ptrs =
        get_varptrs (cctkGH, reflevel, variables::flux_register_fine);
      int const nvars = flux_register_fine_ptrs.size();
      
      for (int n=0; n<nvars; ++n) {
        for (int dir=0; dir<3; ++dir) {
          assert (dim == 3);
#pragma omp parallel
          LC_LOOP3(GRHydro_Reflux_fine_reset,
                   i,j,k,
                   0,0,0, cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
                   cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
          {
            int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
            flux_register_fine_ptrs.AT(n)[ind+dir*np] = 0.0;
          } LC_ENDLOOP3(GRHydro_Reflux_fine_reset);
        }
      }
      
    } END_LOCAL_COMPONENT_LOOP;
  } END_SWITCH_TO_LEVEL;
}  



// Reset the coarse grid flux register (on this level)
static
void flux_register_coarse_reset (cGH const * restrict const cctkGH)
{
  assert (is_singlemap_mode());
  
  BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
    DECLARE_CCTK_ARGUMENTS;
    
    ivect const& lsh = ivect::ref(cctk_lsh);
    int const np = prod(lsh);
    
    vector<CCTK_REAL *> flux_register_coarse_ptrs =
      get_varptrs (cctkGH, reflevel, variables::flux_register_coarse);
    int const nvars = flux_register_coarse_ptrs.size();
    
    for (int n=0; n<nvars; ++n) {
      for (int dir=0; dir<3; ++dir) {
        assert (dim == 3);
#pragma omp parallel
        LC_LOOP3(GRHydro_Reflux_coarse_reset,
                 i,j,k,
                 0,0,0, cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
                 cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
        {
          int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
          flux_register_coarse_ptrs.AT(n)[ind+dir*np] = 0.0;
        } LC_ENDLOOP3(GRHydro_Reflux_coarse_reset);
      }
    }
    
  } END_LOCAL_COMPONENT_LOOP;
}



extern "C"
void Refluxing_CorrectState (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  
  // Sanity check
  assert (mglevel>=0 and reflevel>=0 and Carpet::map>=0 and component==-1);
  
  // We reflux on level L by taking a correction from level L+1 into
  // account.  (This corresponds to restriction, where level L is
  // "corrected" from level L+1.)
  
  CCTK_VInfo (CCTK_THORNSTRING,
              "Refluxing on patch #%d on level %d from level %d",
              Carpet::map, reflevel, reflevel+1);
  
  // There is no finer level; do nothing
  if (reflevel+1 >= reflevels) return;
  
  
  
  // This works only with cell centred grids
  gh const& hh = *vhh.AT(Carpet::map);
  if (hh.refcent != cell_centered) {
    CCTK_WARN (CCTK_WARN_ABORT, "Refluxing requires cell-centred grids");
  }
  
  
  
  if (refluxing_debug_variables) {
    flux_weight_fine_set (cctkGH);
    flux_weight_coarse_set (cctkGH);
  }
  
  reflux (cctkGH);
  
  flux_register_fine_reset (cctkGH);
  flux_register_coarse_reset (cctkGH);
}
