#include "prob.H"

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex::Real* problo,
  const amrex::Real* probhi)
{
  auto& pp_d = *PeleC::h_prob_parm_device;
  {
    amrex::ParmParse pp("prob");
    pp.query("Ms", pp_d.Ms);
    pp.query("p_l", pp_d.p_l);
    pp.query("u_l", pp_d.u_l);
    pp.query("rho_l", pp_d.rho_l);
    pp.query("p_r", pp_d.p_r);
    pp.query("u_r", pp_d.u_r);
    pp.query("rho_r", pp_d.rho_r);
    pp.query("frac", pp_d.frac);
    pp.query("idir", pp_d.idir);
  }
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    pp_d.split[idir] = pp_d.frac * (problo[idir] + probhi[idir]);
  }

  // Air as a thermally perfect O2/N2 mixture (temperature-dependent cp)
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[O2_ID] = 0.233;
  massfrac[N2_ID] = 0.767;
  for (int n = 0; n < NUM_SPECIES; n++) {
    pp_d.massfrac[n] = massfrac[n];
  }

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real wbar = 0.0;
  eos.Y2WBAR(massfrac, wbar);
  const amrex::Real Rgas = pele::physics::Constants::RU / wbar;

  // Region 1 temperature from (rho, p)
  pp_d.T_r = pp_d.p_r / (pp_d.rho_r * Rgas);

  if (pp_d.Ms > 0.0) {
    // Rankine-Hugoniot jump for the given Ms solved with the actual EOS
    // (mass, momentum, energy with h = h(T)), by bisection on
    // eps = rho_1 / rho_2.
    amrex::Real cp = 0.0, cv = 0.0, h1 = 0.0;
    eos.TY2Cp(pp_d.T_r, massfrac, cp);
    eos.TY2Cv(pp_d.T_r, massfrac, cv);
    eos.TY2H(pp_d.T_r, massfrac, h1);
    const amrex::Real a1 = std::sqrt(cp / cv * Rgas * pp_d.T_r);
    const amrex::Real us = pp_d.Ms * a1;
    auto T2_of = [=](amrex::Real eps) {
      return eps * (Rgas * pp_d.T_r + us * us * (1.0 - eps)) / Rgas;
    };
    auto resid = [&](amrex::Real eps) {
      amrex::Real h2 = 0.0;
      eos.TY2H(T2_of(eps), massfrac, h2);
      return h2 + 0.5 * us * us * eps * eps - h1 - 0.5 * us * us;
    };
    amrex::Real lo = 0.05, hi = 0.6;
    amrex::Real flo = resid(lo);
    if (flo * resid(hi) > 0.0) {
      amrex::Abort("prob: could not bracket the Rankine-Hugoniot solution");
    }
    for (int it = 0; it < 200; it++) {
      const amrex::Real mid = 0.5 * (lo + hi);
      const amrex::Real fm = resid(mid);
      if (flo * fm <= 0.0) {
        hi = mid;
      } else {
        lo = mid;
        flo = fm;
      }
    }
    const amrex::Real eps = 0.5 * (lo + hi);
    pp_d.T_l = T2_of(eps);
    pp_d.rho_l = pp_d.rho_r / eps;
    pp_d.p_l = pp_d.rho_l * Rgas * pp_d.T_l;
    pp_d.u_l = us * (1.0 - eps);
    amrex::Print() << "prob: Ms=" << pp_d.Ms << " gamma1=" << cp / cv
                   << " a1=" << a1 << " T1=" << pp_d.T_r
                   << "\nprob: region 2 (variable-cp RH): rho2/rho1="
                   << 1.0 / eps << " rho2=" << pp_d.rho_l << " p2=" << pp_d.p_l
                   << " u2=" << pp_d.u_l << " T2=" << pp_d.T_l << std::endl;
  } else {
    pp_d.T_l = pp_d.p_l / (pp_d.rho_l * Rgas);
  }

  amrex::Real e_l = 0.0, e_r = 0.0;
  eos.RYP2E(pp_d.rho_l, massfrac, pp_d.p_l, e_l);
  pp_d.rhoe_l = pp_d.rho_l * e_l;
  eos.RYP2E(pp_d.rho_r, massfrac, pp_d.p_r, e_r);
  pp_d.rhoe_r = pp_d.rho_r * e_r;
}
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}
