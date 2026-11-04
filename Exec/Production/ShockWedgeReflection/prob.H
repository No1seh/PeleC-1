#ifndef PROB_H
#define PROB_H

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

#include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"
#include "Constants.H"
#include "PelePhysics.H"
#include "Tagging.H"

#include "ProblemSpecificFunctions.H"
#include "prob_parm.H"

// Oblique shock-wedge reflection (regular / Mach reflection) problem,
// replicating Valencia et al., COB-2023-0803 (27th ABCM / COBEM 2023):
// a moving planar shock (region 2, post-shock) propagates in +x into a
// quiescent ambient (region 1) and interacts with a rigid wedge modeled
// as an embedded boundary (eb2.geom_type = plane) further downstream.

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];

  if (x <= prob_parm.split[0]) {
    // Region 2: post-shock state, flowing in +x behind the incident shock
    state(i, j, k, URHO) = prob_parm.rho_l;
    state(i, j, k, UMX) = prob_parm.rho_l * prob_parm.u_l;
    state(i, j, k, UMY) = 0.0;
    state(i, j, k, UMZ) = 0.0;
    state(i, j, k, UEDEN) =
      prob_parm.rhoe_l +
      0.5 * prob_parm.rho_l * prob_parm.u_l * prob_parm.u_l;
    state(i, j, k, UEINT) = prob_parm.rhoe_l;
    state(i, j, k, UTEMP) = prob_parm.T_l;
  } else {
    // Region 1: quiescent ambient ahead of the incident shock
    state(i, j, k, URHO) = prob_parm.rho_r;
    state(i, j, k, UMX) = prob_parm.rho_r * prob_parm.u_r;
    state(i, j, k, UMY) = 0.0;
    state(i, j, k, UMZ) = 0.0;
    state(i, j, k, UEDEN) =
      prob_parm.rhoe_r +
      0.5 * prob_parm.rho_r * prob_parm.u_r * prob_parm.u_r;
    state(i, j, k, UEINT) = prob_parm.rhoe_r;
    state(i, j, k, UTEMP) = prob_parm.T_r;
  }

  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = 0.0;
  }
  state(i, j, k, UFS) = state(i, j, k, URHO);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
  if ((idir == 0) && (sgn == -1)) {
    // Fixed inlet boundary at x = prob_lo[0]: continuously feed the
    // post-shock (region 2) state, as in Figure 2 of the reference paper.
    s_ext[URHO] = prob_parm.rho_l;
    s_ext[UMX] = prob_parm.rho_l * prob_parm.u_l;
    s_ext[UMY] = 0.0;
    s_ext[UMZ] = 0.0;
    s_ext[UEDEN] =
      prob_parm.rhoe_l +
      0.5 * prob_parm.rho_l * prob_parm.u_l * prob_parm.u_l;
    s_ext[UEINT] = prob_parm.rhoe_l;
    s_ext[UTEMP] = prob_parm.T_l;
    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[UFS + n] = 0.0;
    }
    s_ext[UFS] = s_ext[URHO];
  } else {
    for (int n = 0; n < NVAR; n++) {
      s_ext[n] = s_int[n];
    }
  }
}

void pc_prob_close();

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;

#endif
