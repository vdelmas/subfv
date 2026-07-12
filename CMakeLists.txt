cmake_minimum_required(VERSION 3.20)

project(subfv Fortran)

find_package(MPI REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# ── Core library (partagée par tous les targets) ───────────────────────────────
add_library(core STATIC
  src/utils/precision_module.F90
  src/utils/sort_module.F90
  src/mpi/mpi_module.F90
  src/mesh/mesh_module.F90
  src/mesh/mesh_reading_module.F90
  src/mesh/mesh_geometry_module.F90
  src/mesh/mesh_connectivity_module.F90
  src/linear_algebra/sparse_csr_linear_module.F90
  src/linear_algebra/sparse_bcsr_linear_module.F90
  src/linear_algebra/sparse_linear_module.F90
  src/linear_algebra/linear_solver_module.F90
  src/io/io_module.F90
)
target_link_libraries(core PUBLIC MPI::MPI_Fortran LAPACK::LAPACK BLAS::BLAS)
set_target_properties(core PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/core)
target_include_directories(core PUBLIC ${CMAKE_BINARY_DIR}/mod/core)

# ── NS library (partagée par subfvns et subfvshear) ───────────────────────────
add_library(ns STATIC
  src/ns/ns_global_data_module.F90
  src/ns/ns_io_module.F90
  src/ns/ns_euler_primitives_module.F90
  src/ns/ns_euler_recon_module.F90
  src/ns/ns_euler_zb_module.F90
  src/ns/ns_euler_exact_sol_module.F90
  src/ns/ns_mesh_metric_module.F90
  src/ns/ns_euler_rs_module.F90
  src/ns/ns_euler_module.F90
  src/ns/ns_diffusion_module.F90
  src/ns/ns_vectorial_diffusion_module.F90
)
target_link_libraries(ns PUBLIC core)
set_target_properties(ns PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/ns)
target_include_directories(ns PUBLIC ${CMAKE_BINARY_DIR}/mod/ns)

# ── Acoustic library ──────────────────────────────────────────────────────────
add_library(acoustic STATIC
  src/acoustic/acoustic_module.F90
)
target_link_libraries(acoustic PUBLIC core)
set_target_properties(acoustic PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/acoustic)
target_include_directories(acoustic PUBLIC ${CMAKE_BINARY_DIR}/mod/acoustic)

# ── Lagrange library ──────────────────────────────────────────────────────────
add_library(lagrange STATIC
  src/lagrange/lagrange_global_data_module.F90
  src/lagrange/lagrange_io_module.F90
  src/lagrange/lagrange_module.F90
)
target_link_libraries(lagrange PUBLIC core)
set_target_properties(lagrange PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/lagrange)
target_include_directories(lagrange PUBLIC ${CMAKE_BINARY_DIR}/mod/lagrange)

# ── Exécutables ───────────────────────────────────────────────────────────────
add_executable(subfvns      src/ns/ns_main.F90)
add_executable(subfvshear   src/shear/shear_main.F90)
add_executable(subfvacoustic src/acoustic/acoustic_main.F90)
add_executable(subfvlagrange src/lagrange/lagrange_main.F90)

target_link_libraries(subfvns      ns)
target_link_libraries(subfvshear   ns)
target_link_libraries(subfvacoustic acoustic)
target_link_libraries(subfvlagrange lagrange)

install(TARGETS subfvns subfvshear subfvacoustic subfvlagrange
  RUNTIME DESTINATION bin
)

install(PROGRAMS tools/subfv-gmsh/subfv-gmsh
  DESTINATION bin
)

enable_testing()
add_subdirectory(test)
