cmake_minimum_required(VERSION 3.20)

project(subfv Fortran)

find_package(MPI REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# ── Core library (partagée par tous les targets) ───────────────────────────────
add_library(core STATIC
  src/utils/precision_module.F90
  src/utils/sort_module.F90
  src/utils/quadrature_module.F90
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

# ── RBF library ────────────────────────────────────────────────────────────────
add_library(rbf STATIC
  src/rbf/rbf_module.F90
)
target_link_libraries(rbf PUBLIC core)
set_target_properties(rbf PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/rbf)
target_include_directories(rbf PUBLIC ${CMAKE_BINARY_DIR}/mod/rbf)

# ── ALE library ────────────────────────────────────────────────────────────────
# Self-contained: its own local multi_point-style nodal Riemann solver and
# mesh-motion helpers, independent of ns/lagrange (only links core), so it
# can never affect subfvns/subfvlagrange.
add_library(ale STATIC
  src/ale/ale_global_data_module.F90
  src/ale/ale_rbf_module.F90
  src/ale/ale_module.F90
)
target_link_libraries(ale PUBLIC core)
set_target_properties(ale PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/ale)
target_include_directories(ale PUBLIC ${CMAKE_BINARY_DIR}/mod/ale)

# ── High-order Euler library ──────────────────────────────────────────────────
add_library(euler_ho STATIC
  src/euler_ho/euler_ho_module.F90
)
target_link_libraries(euler_ho PUBLIC core)
set_target_properties(euler_ho PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/euler_ho)
target_include_directories(euler_ho PUBLIC ${CMAKE_BINARY_DIR}/mod/euler_ho)

# ── Arbitrary high-order derivative library ───────────────────────────────────
add_library(arbitrary_high_order STATIC
  src/arbitrary_high_order/arbitrary_high_order_module.F90
)
target_link_libraries(arbitrary_high_order PUBLIC core)
set_target_properties(arbitrary_high_order PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/arbitrary_high_order)
target_include_directories(arbitrary_high_order PUBLIC ${CMAKE_BINARY_DIR}/mod/arbitrary_high_order)

# euler_ho can optionally source its order>=2 reconstruction from
# arbitrary_high_order_module (use_aho_reconstruction), for a like-for-like
# comparison against its own ls_reconstruction.
target_link_libraries(euler_ho PUBLIC arbitrary_high_order)

# ── Exécutables ───────────────────────────────────────────────────────────────
add_executable(subfvns           src/ns/ns_main.F90)
add_executable(subfvshear        src/shear/shear_main.F90)
add_executable(subfvacoustic     src/acoustic/acoustic_main.F90)
add_executable(subfvlagrange     src/lagrange/lagrange_main.F90)
add_executable(subfvgreshoconv   src/ns/ns_gresho_convergence_main.F90)
add_executable(subfvrbf          src/rbf/rbf_main.F90)
add_executable(subfveulerho      src/euler_ho/euler_ho_main.F90)
add_executable(subfvahotest      src/arbitrary_high_order/arbitrary_high_order_main.F90)
add_executable(subfvale          src/ale/ale_main.F90)

target_link_libraries(subfvns          ns)
target_link_libraries(subfvshear       ns)
target_link_libraries(subfvacoustic    acoustic)
target_link_libraries(subfvlagrange    lagrange)
target_link_libraries(subfvgreshoconv  ns)
target_link_libraries(subfvrbf         rbf)
target_link_libraries(subfveulerho     euler_ho)
target_link_libraries(subfvahotest     arbitrary_high_order)
target_link_libraries(subfvale         ale)

install(TARGETS subfvns subfvshear subfvacoustic subfvlagrange subfvgreshoconv subfvrbf subfvale
  RUNTIME DESTINATION bin
)

install(PROGRAMS tools/subfv-gmsh/subfv-gmsh
  DESTINATION bin
)

enable_testing()
add_subdirectory(test)
