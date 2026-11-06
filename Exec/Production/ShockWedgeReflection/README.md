# Oblique Shock-Wave Wedge Reflection (COB-2023-0803)

Replication of:

> S. Valencia, C. Celis, A. Mendiburu, L. Bravo, P. Khare, "Fluid Dynamic
> Simulations of Mach and Regular Reflections in Oblique Shock-Wave
> Configurations using Adaptive Mesh Refinement," 27th ABCM International
> Congress of Mechanical Engineering (COBEM 2023), COB-2023-0803.
> DOI: 10.48550/arXiv.2309.05882

## Relationship to PeleC's built-in EB-C14 test

PeleC already ships a generic "Mach reflection" EB verification case at
`Exec/RegTests/EB-C14` (shock colliding with a ramp, `eb2.geom_type = plane`,
multi-level AMR, PPM, GammaLaw) — the same ingredients used here. This case
reuses that exact approach, but differs in one deliberate way:

- **EB-C14** initializes a two-state Riemann problem with **zero velocity on
  both sides** and lets the shock self-form from the pressure/density jump.
  That is fine for a qualitative check that a Mach stem/triple point forms,
  but the resulting incident shock's Mach number is whatever the full
  Riemann solution produces — it is not set directly from a target `Ms`.
- **This case** needs the *exact* `Ms` and wedge angle from Table 1 of
  COB-2023-0803, so instead the post-shock (region 2) state is computed
  directly from the Rankine-Hugoniot relations (Eqs. 1-4, which assume a
  quiescent ambient and therefore a **nonzero** induced post-shock velocity)
  and is fed in continuously as a `Hard` inflow at x-lo — matching the
  paper's Figure 2 "fixed inlet boundary". Using EB-C14's zero-velocity
  trick here would not reproduce the paper's prescribed `Ms` values.

Everything else (EB plane wedge, PPM, GammaLaw EOS, AMR) follows EB-C14's
pattern directly. The smoke test for Case 3 confirmed this: simulated
density/pressure min-max matched the hand-derived Rankine-Hugoniot states
exactly.

## Problem setup

- 2D domain: x in [0, 4] m, y in [0, 0.75] m (400 x 75 cm in the CGS units
  PeleC/AMReX use internally).
- Base grid: 512 x 96 cells (dx = dy = 7.8125 mm), refined up to 3 AMR
  levels (ref_ratio = 2), giving a finest cell size of 0.977 mm, matching
  the mesh-independence study in Section 5.1 of the paper.
- A planar incident shock is initialized at x = 0.5 m, separating:
  - Region 2 (post-shock, x <= 0.5 m): computed from the ambient state and
    the incident shock Mach number Ms via the Rankine-Hugoniot relations
    (Eqs. 1-4 in the paper).
  - Region 1 (quiescent ambient, x > 0.5 m): taken directly from Table 1.
- A rigid wedge is modeled as an embedded boundary (`eb2.geom_type = plane`)
  with its apex at x = 3.0 m, inclined at the case's wedge angle theta_w.
- Governing equations: inviscid compressible Euler (`pelec.do_react = 0`,
  `diffuse_vel/temp/spec = 0`), GammaLaw EOS (ideal gas, gamma = 1.4,
  air molecular weight), discretized with the unsplit PPM scheme
  (`pelec.do_mol = 0`) and a CFL number of 0.3, as specified in the paper.
- Boundary conditions: `Hard` (fixed post-shock inflow) at x-lo,
  `FOExtrap` (outflow) at x-hi, and `SlipWall` (inviscid reflecting wall)
  at y-lo/y-hi and on the wedge EB — the inviscid-flow equivalent of the
  paper's "no-slip wall" boundaries, since no molecular diffusion is solved.
- AMR refinement: density-ratio based tagging (`tagging.denratio = 1.05`,
  scale-independent) plus static refinement of the wedge embedded boundary
  to the finest level, matching Figure 5 of the paper.

## Cases (Table 1 of the paper)

| Case | theta_w | Ms   | Regime | Input file |
|------|---------|------|--------|------------|
| 1    | 63.4°   | 2.05 | RR     | `case1.inp` |
| 2    | 60°     | 4.70 | RR/MR dual region | `case2.inp` |
| 3    | 27°     | 2.03 | SMR    | `case3.inp` |
| 4    | 20°     | 7.19 | TMR    | `case4.inp` |

All post-shock (region 2) states were derived from Table 1's ambient
conditions using the Rankine-Hugoniot relations (Eqs. 1-4), with
gamma = 1.4 and R_air = 287 J/(kg K), then converted to the CGS unit
system (dyn/cm^2, g/cm^3, cm/s) that PeleC/AMReX use internally.

## SUNDIALS

This build does not solve chemistry (`pelec.do_react = 0`), so SUNDIALS
(CVODE) is not exercised at runtime for this case, but it is still built
automatically as part of the PeleC/PelePhysics build (see
`CMake/SetSundialsOptions.cmake`) from the `Submodules/PelePhysics/Submodules/sundials`
submodule — no separate/system SUNDIALS installation is required or used.

## Building

```bash
cd PeleC
cmake -B build -DCMAKE_BUILD_TYPE=Release -DPELE_DIM=2 \
  -DPELE_ENABLE_MPI=OFF -DPELE_ENABLE_TESTS=OFF
cmake --build build --target PeleC-ShockWedgeReflection -j"$(nproc)"
```

## Running

```bash
cd Exec/Production/ShockWedgeReflection
../../../build/Exec/Production/ShockWedgeReflection/PeleC-ShockWedgeReflection case1.inp
```

Plotfiles (`plt1_*`, etc.) can be visualized with ParaView, VisIt, or
`yt`, and compared against Figures 6-10 of the paper (density ratio
contours and wedge-wall density profiles).
