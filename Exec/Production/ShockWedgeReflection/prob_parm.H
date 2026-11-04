#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p_l = 1.0;   // post-shock (region 2) pressure
  amrex::Real u_l = 0.0;   // post-shock (region 2) x-velocity
  amrex::Real rho_l = 0.0; // post-shock (region 2) density
  amrex::Real rhoe_l = 0.0;
  amrex::Real T_l = 1.0;
  amrex::Real p_r = 0.1;     // ambient (region 1) pressure
  amrex::Real u_r = 0.0;     // ambient (region 1) x-velocity
  amrex::Real rho_r = 0.125; // ambient (region 1) density
  amrex::Real rhoe_r = 0.0;
  amrex::Real T_r = 1.0;
  amrex::Real frac = 0.5; // fraction of the domain where the initial shock sits
  bool use_Tinit = false;
  int idir = 1; // direction across which the initial shock jump occurs
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> split;
};

struct ProbParmHost
{
  ProbParmHost() = default;
};
#endif
