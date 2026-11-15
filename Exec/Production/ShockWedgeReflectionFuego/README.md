# Shock-wedge reflection with temperature-dependent cp (Fuego, air)

Variant of `../ShockWedgeReflection` (replication of COB-2023-0803) that
replaces the constant-gamma gas (GammaLaw) with thermally perfect air: Fuego
EOS with the non-reacting two-species `air` mechanism (O2/N2, mass fractions
0.233/0.767). Geometry, AMR, numerics and boundary conditions are unchanged.
Only Cases 2 and 4 (the high Mach number ones) were run.

Two ways of setting the post-shock (region 2) inflow state:

- `caseN_fuego.inp`: `prob.Ms` is given and region 2 is obtained by solving
  the Rankine-Hugoniot jump with the same variable-cp EOS (done in
  `prob.cpp`).
- `caseN_fuego_gammaRH.inp`: region 2 is the constant-gamma (1.4)
  Rankine-Hugoniot state of Eqs. 1-4 of the paper, imposed on the
  variable-cp gas.

## Results (wall density ratio, paper values read by eye off Figure 10)

|                  | Paper      | GammaLaw | Fuego, consistent RH | Fuego, gamma RH inflow |
|------------------|------------|----------|----------------------|------------------------|
| Case 2 inflow    | ~4.9       | 4.88     | 5.40                 | 4.84                   |
| Case 2 plateau   | ~15.3      | 14.6     | 18.2                 | 17.9                   |
| Case 4 inflow    | ~6.5       | 5.47     | 6.56                 | 5.38                   |
| Case 4 plateau   | ~10.0-10.2 | 8.3      | 10.95                | 10.8                   |

None of the three configurations reproduces the paper for both cases: Case 2
is closest with GammaLaw (-5%), while the Case 4 inflow value only matches
with variable cp, which then overshoots the wedge plateau by ~8%.

## Building and running

```bash
cmake --build build --target PeleC-ShockWedgeReflectionFuego -j"$(nproc)"
cd Exec/Production/ShockWedgeReflectionFuego
mpirun -np 8 ../../../build/Exec/Production/ShockWedgeReflectionFuego/PeleC-ShockWedgeReflectionFuego case4_fuego.inp
```
