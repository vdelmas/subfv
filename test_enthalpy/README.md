# `test_enthalpy` — enthalpy-preserving Simple Riemann solver

Test cases and `three_wave` / `three_wave_enthalpy` comparisons for the enthalpy-preserving
variant of the Gallice 1D Simple Riemann solver, after

> L. Tallois, *Enthalpy preserving Simple Riemann solvers for steady hypersonic flows*,
> CEA-CESTA preprint, June 2026 (`papers/talois.pdf`), section 5.

## What the scheme is

For a steady Euler solution the total enthalpy `h = e + p/rho` is a Riemann invariant along
streamlines, so `h == h_inf` everywhere downstream of a uniform inflow, **shocks included**.
Godunov-type solvers (Rusanov, HLL, HLLC, Roe — and the `three_wave` already in this code,
which is exactly the paper's Gallice-1D solver) do not reproduce that at the discrete level:
they leave an enthalpy undershoot inside a stationary shock, which pollutes the energy
distribution of the whole subsonic pocket behind a bow shock, and with it any
Navier-Stokes boundary-layer criterion based on total enthalpy.

`three_wave_enthalpy` is the same 3-wave solver written in the variable set `(rho, rho*u, rho*h)`
instead of `(rho, rho*u, rho*e)`: same Lagrangian slopes, same `u*`, same `rho*`, true physical
fluxes untouched — only the energy slot of the four solver states (both outer states and both
star states) carries `rho*h`. With `h_l = h_r = h_inf` every energy component is then `h_inf`
times the matching density component and the flux assembly factors exactly as
`F(5) = h_inf * F(1)`, which is Hänel's condition.

The paper also gives a **multidimensional** version (MGallice-2D), built on the same Simple
solver but with the intermediate normal velocity taken from a NODAL velocity (the paper's
linear system (16)) instead of the face value. The code already has its baseline: `multi_point`
in `subfvns` **is** the paper's Gallice-2D scheme, non-classical consistency term
`-/+ 1/2 (p*_r - p*_l) (0, n, u*)` included. `multi_point_enthalpy` is its enthalpy-preserving
counterpart — same substitution of `rho*h` for `rho*e` in all four states, **plus one thing the
1D case does not have**: rewriting the solver in `rho*h` makes the energy component of that
non-classical term vanish identically, so the correction vector becomes `(0, n, 0)` and not
`(0, n, u*)`. Keeping the `u*` component would reintroduce a one-sided energy dissipation and
destroy the preservation.

Implemented in three places:

| solver | file | selected by |
|---|---|---|
| aho Euler (`subfveulerho`), 1D | `src/euler_ho/euler_ho_module.F90::three_wave_enthalpy_flux` | `flux_scheme='three_wave_enthalpy'` |
| Navier-Stokes (`subfvns`), 1D | `src/ns/ns_euler_rs_module.F90::three_wave_enthalpy` | `scheme='three_wave_enthalpy'` |
| Navier-Stokes (`subfvns`), multidimensional | `src/ns/ns_euler_rs_module.F90::multi_point_enthalpy` | `scheme='multi_point_enthalpy'` |

Properties verified numerically on the exact implemented expressions: consistency at a constant
state, integral-form (HLL) consistency `sum_k Lambda_k dU_k = F_r - F_l` (so the solver is a
genuine Godunov-type solver and the symmetric flux assembly is the right one), Hänel's condition
to machine precision, and reduction to the upwind flux when the face is fully supersonic.

**Caveat at order >= 2.** The paper additionally reconstructs `h` itself and feeds the
reconstructed `h` to the solver. Neither implementation here does that — `h` is rebuilt from the
reconstructed `(rho, u, p)`, so above first order the reconstruction errors no longer cancel
exactly and `h` is preserved only to truncation order. At first order it is exact.

**Entropy stability is not claimed** — as for flux vector splitting schemes, no slope condition is
known that guarantees it for this modification (paper, Remark 3).

## Cases

### `half_cylinder/` — inviscid, `subfveulerho`

Hypersonic M=20 flow over a half cylinder, non-dimensionalised as in the paper's section 6.2:
`(rho, ux, uy, p) = (1, M*sqrt(gamma), 0, 1)`, `gamma = 7/5`, hence **`h_inf = 283.5`**.
Left half-annulus, cylinder radius 1, outer (freestream) radius 2.5. Two meshes built from the
project's own `test/half_cylinder_{quad,tri}_euler` geometries, refined: `cyl_quad` 48x128 hexes,
`cyl_tri` unstructured prisms.

```
NPROCS=10 ./run.sh quad three_wave_enthalpy 1 8.0    # <quad|tri> <scheme> <order> <tmax>
pvpython enthalpy_stats.py 283.5 outputs/*/output_-1.pvtu
./make_figures.sh                                    # schlieren + H panels, shared colour range
```

### `half_cylinder_ns/` — viscous, `subfvns`

The project's existing `test/half_cylinder_{quad,tri}_ns` setup unchanged (M=17.6,
`rho=1e-3`, `u=5000`, `p=57.615576818289696`, isothermal wall at 500 K, Sutherland), hence
`h_inf = 1.2701654518864e7`. Only the Riemann solver differs between runs.

Here `h` is **not** expected to be uniform: viscous work and wall heat flux make it drop through
the boundary layer. That is the point — with `three_wave_enthalpy` the remaining deficit is
physical (boundary layer) instead of the numerical undershoot `three_wave` spreads from the bow
shock, which is what a total-enthalpy-based boundary-layer criterion needs.

Run on Curta with `job_<mesh>_<scheme>_o1.sbatch` (32 ranks, `imb`); locally via
`./run.sh quad three_wave_enthalpy`. Each job also exports `coeffs_export.csv` (wall Cp and
heat flux) with `export_coeffs.py`; `./plot_coeffs.sh <quad|tri>` then plots Cp and the heating
rate C_H for all four schemes against the LAURA reference (`fun3D_cp.csv` / `fun3D_st.csv`,
the same reference the project's own `test/half_cylinder_*_ns` cases use).

### `half_cylinder_mp/` — inviscid, `subfvns`, multidimensional solvers

`multi_point` vs `multi_point_enthalpy` (Gallice-2D vs MGallice-2D) on exactly the meshes and
freestream of `half_cylinder/`, so the panels are directly comparable to the
`three_wave` / `three_wave_enthalpy` ones. Run through `subfvns` with
`activate_diffusion = .false.`, because the node-based solvers exist only there — `euler_ho`
carries the face-based ones only.

```
NPROCS=6 ./run.sh quad multi_point_enthalpy
pvpython enthalpy_stats.py 283.5 outputs/*/output_-1.pvtu
pvpython steadiness.py outputs/*                      # is the run actually steady?
./make_figures.sh && python3 montage.py
```

## Result so far (inviscid, first order)

On the **quad** mesh the property is reproduced exactly: with `three_wave_enthalpy` the median
cell has `|H - h_inf|/h_inf = 1.3e-13` and the 90th percentile is `1.6e-07`, against `5.1e-03`
for `three_wave` — a factor 3e4. The residual max (9.3e-03) sits on a single cell.

On the **tri** mesh it does not show: q90 goes 7.7e-03 -> 8.5e-03, i.e. no gain. That is not a
defect of the solver. `steadiness.py` shows the tri runs are not converged —
`||d rho|| / ||rho||` between the last two dumps is 1.3e-02 / 2.1e-02 against 1.3e-05 / 1.5e-04
on the quad, three orders of magnitude larger and the same size as the measured enthalpy
departure. Enthalpy preservation is a STEADY-state property, so it simply does not apply there.
The cause is the carbuncle instability of the ONE-DIMENSIONAL solvers along the stagnation line,
plainly visible in the tri schlieren panel, and the paper says exactly this (sec. 6.2): the 1D
schemes fail on this case while the multidimensional ones are carbuncle-free. Hence
`half_cylinder_mp/`.

## Scripts

| script | what it does |
|---|---|
| `enthalpy_stats.py` | distribution of `\|H - h_inf\|/h_inf` over the cells, per run |
| `steadiness.py` | `\|\|d rho\|\|/\|\|rho\|\|` between the last two dumps — is the run steady at all? |
| `ranges.py` | shared colour ranges over a set of runs, so panels are comparable |
| `render.py` | one run -> schlieren PNG + total-enthalpy PNG (viridis), offscreen pvbatch |
| `make_figures.sh` | `ranges.py` + `render.py` over everything in `outputs/` |
| `montage.py` | assembles the eight panels into one `figures/comparison.png` |

`render.py`, `ranges.py` and `steadiness.py` auto-detect which solver wrote the file
(`grad_rho_cell`/`rho` for `euler_ho`, `Nodal_Grad_Density`/`Density` for `subfvns`).
