# subfv

Subface-based Finite Volume solver written in Fortran, built with CMake. Provides four executables covering compressible NS, shear flow, acoustics and Lagrange hydrodynamics.

## Requirements

- Fortran compiler (gfortran, ifort, ...)
- CMake >= 3.20
- MPI (OpenMPI or MPICH)
- BLAS & LAPACK
- Gmsh (to generate meshes from `.geo` files)
- ParaView (optional, for visualization)

Install all dependencies on Ubuntu/Debian:
```bash
sudo apt update && sudo apt install -y cmake gfortran openmpi-bin libopenmpi-dev libblas-dev liblapack-dev gmsh paraview
```

## Build

```bash
cmake -B build
make -C build -j
```

If MPI is not detected automatically:
```bash
cmake -B build -DMPI_Fortran_COMPILER=mpif90
```

To set the compiler explicitly:
```bash
cmake -B build -DCMAKE_Fortran_COMPILER=gfortran
```

## Install

```bash
cmake -B build -DCMAKE_INSTALL_PREFIX=$HOME/.local
make -C build -j
cmake --install build
```

This installs into `$HOME/.local/bin/`:
- `subfvns` — compressible Navier-Stokes / Euler solver
- `subfvshear` — shear flow solver
- `subfvacoustic` — acoustic solver
- `subfvlagrange` — Lagrange hydrodynamics solver
- `subfv-gmsh` — mesh generation tool (wrapper around Gmsh)

Add to your shell configuration:
```bash
export PATH=$HOME/.local/bin:$PATH
```

## Running a test

Meshes are generated from `.geo` files using `subfv-gmsh`. Example with the Gresho vortex:
```bash
cd test/low_mach/gresho_quad
subfv-gmsh gresho_quad.geo
bash run.sh two_wave .false. 0
```

## Running the test suite

From the build directory:
```bash
ctest --test-dir build -j4
```

Filter by category with labels:
```bash
ctest --test-dir build -L quick_test   # fast validation cases
ctest --test-dir build -L riemann_problem_1d
ctest --test-dir build -L sedov
ctest --test-dir build -L half_cylinder
ctest --test-dir build -L low_mach
ctest --test-dir build -L shear
```

## Project structure

```
src/
├── utils/          precision, sorting
├── mpi/            MPI abstraction
├── mesh/           mesh reading, geometry, connectivity
├── linear_algebra/ sparse solvers (CSR, BCSR)
├── io/             input/output utilities
├── ns/             compressible NS/Euler (subfvns, subfvshear)
├── shear/          shear flow (subfvshear)
├── acoustic/       acoustics (subfvacoustic)
└── lagrange/       Lagrange hydrodynamics (subfvlagrange)

test/
├── riemann_problem_1d/   1D Riemann problems (Toro tests)
├── sedov/                Sedov blast wave (tri, hex, tet)
├── half_cylinder/        Half-cylinder Euler and NS
├── low_mach/             Gresho, Kelvin-Helmholtz, potential flow
├── shear/                Shear flow convergence
├── scaling_mpi/          MPI scaling
└── test_error_heat_flux_hybrid/  Heat flux error study

tools/
└── subfv-gmsh/     mesh generation tool
```

## Visualization

Results are written as `.vtu` files (one per MPI subdomain). Open the `.pvtu` header file in ParaView. Apply the D3 filter (`Filters > Search > D3`) to correctly handle subdomain boundaries for contours and iso-surfaces.

## Author

Vincent Delmas, University of Bordeaux, Institut de Mathématiques de Bordeaux (IMB)
