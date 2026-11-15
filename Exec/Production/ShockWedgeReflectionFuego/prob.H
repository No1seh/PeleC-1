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

// Variant of Production/ShockWedgeReflection (Valencia et al.,
// COB-2023-0803) that swaps the constant-gamma gas (GammaLaw) for
// thermally perfect air (Fuego EOS, non-reacting O2/N2 "air" mechanism,
// temperature-dependent cp). Region 2 is obtained from region 1 and
// prob.Ms by solving the Rankine-Hugoniot jump with that same EOS.

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
    state(i, j, k, UFS + n) = state(i, j, k, URHO) * prob_parm.massfrac[n];
  }
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
      s_ext[UFS + n] = s_ext[URHO] * prob_parm.massfrac[n];
    }
  } else {
    for (int n = 0; n < NVAR; n++) {
      s_ext[n] = s_int[n];
    }
  }
}

void pc_prob_close();

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;

#endif
