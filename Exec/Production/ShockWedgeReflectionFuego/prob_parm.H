#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>
#include "mechanism.H"

struct ProbParmDevice
{
  amrex::Real Ms = -1.0;   // incident shock Mach number; if > 0, region 2 is
                           // computed from region 1 with the actual EOS
  amrex::Real p_l = 1.0;   // post-shock (region 2) pressure
  amrex::Real u_l = 0.0;   // post-shock (region 2) x-velocity
  amrex::Real rho_l = 0.0; // post-shock (region 2) density
  amrex::Real rhoe_l = 0.0;
  amrex::Real T_l = 300.0;
  amrex::Real p_r = 0.1;     // ambient (region 1) pressure
  amrex::Real u_r = 0.0;     // ambient (region 1) x-velocity
  amrex::Real rho_r = 0.125; // ambient (region 1) density
  amrex::Real rhoe_r = 0.0;
  amrex::Real T_r = 300.0;
  amrex::Real frac = 0.5; // fraction of the domain where the initial shock sits
  int idir = 1;
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> split;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {0.0};
};

struct ProbParmHost
{
  ProbParmHost() = default;
};
#endif
