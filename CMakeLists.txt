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

# ── Parallel CGNS reading + in-code partitioning (src/partition/) ─────────────
# A .cgns meshfile is then read block-wise by every rank, partitioned
# (ParMETIS) and given its node-based ghost layer in the code itself, with
# no subfv-gmsh pre-partitioning. Needs CGNS built on parallel HDF5 and,
# optionally, ParMETIS (otherwise only SUBFV_PARTITIONER=block is available).
option(SUBFV_ENABLE_CGNS "Read .cgns meshes in parallel with in-code partitioning" ON)
if(SUBFV_ENABLE_CGNS)
  set(CGNS_ROOT "$ENV{HOME}/opt/cgns" CACHE PATH "CGNS install prefix")
  set(HDF5_ROOT "$ENV{HOME}/opt/hdf5" CACHE PATH "HDF5 (parallel) install prefix")
  set(PARMETIS_ROOT "$ENV{HOME}/opt/parmetis" CACHE PATH "ParMETIS/METIS/GKlib install prefix")
  find_library(CGNS_LIBRARY cgns HINTS ${CGNS_ROOT}/lib)
  find_path(CGNS_MOD_DIR cgns.mod HINTS ${CGNS_ROOT}/include)
  find_library(HDF5_C_LIBRARY hdf5 HINTS ${HDF5_ROOT}/lib NO_DEFAULT_PATH)
  if(NOT CGNS_LIBRARY OR NOT CGNS_MOD_DIR OR NOT HDF5_C_LIBRARY)
    message(WARNING "subfv: CGNS/HDF5 not found (CGNS_ROOT=${CGNS_ROOT}, HDF5_ROOT=${HDF5_ROOT}): "
      "building WITHOUT the parallel .cgns reader (.msh + subfv-gmsh only)")
    set(SUBFV_ENABLE_CGNS OFF)
  endif()
endif()
if(SUBFV_ENABLE_CGNS)
  target_sources(core PRIVATE
    src/partition/dist_utils_module.F90
    src/partition/cgns_dist_read_module.F90
    src/partition/dist_mesh_module.F90
  )
  target_include_directories(core PUBLIC ${CGNS_MOD_DIR})
  target_compile_definitions(core PRIVATE SUBFV_HAVE_CGNS)
  target_link_libraries(core PUBLIC ${CGNS_LIBRARY} ${HDF5_C_LIBRARY})

  find_library(PARMETIS_LIBRARY parmetis HINTS ${PARMETIS_ROOT}/lib NO_DEFAULT_PATH)
  find_library(METIS_LIBRARY metis HINTS ${PARMETIS_ROOT}/lib NO_DEFAULT_PATH)
  find_library(GKLIB_LIBRARY GKlib HINTS ${PARMETIS_ROOT}/lib NO_DEFAULT_PATH)
  # ParMETIS is called directly from Fortran (its own Fortran entry points), which
  # requires 32-bit idx_t/real_t to match integer(c_int)/real(c_float).
  find_file(METIS_HEADER metis.h HINTS ${PARMETIS_ROOT}/include NO_DEFAULT_PATH)
  if(METIS_HEADER)
    file(STRINGS ${METIS_HEADER} METIS_IDXW REGEX "^#define IDXTYPEWIDTH")
    file(STRINGS ${METIS_HEADER} METIS_REALW REGEX "^#define REALTYPEWIDTH")
    if(NOT METIS_IDXW MATCHES "32" OR NOT METIS_REALW MATCHES "32")
      message(WARNING "subfv: ParMETIS not built with 32-bit idx_t/real_t, ignoring it")
      set(PARMETIS_LIBRARY PARMETIS_LIBRARY-NOTFOUND)
    endif()
  endif()
  if(PARMETIS_LIBRARY AND METIS_LIBRARY)
    target_compile_definitions(core PRIVATE SUBFV_HAVE_PARMETIS)
    target_link_libraries(core PUBLIC ${PARMETIS_LIBRARY} ${METIS_LIBRARY})
    if(GKLIB_LIBRARY)
      target_link_libraries(core PUBLIC ${GKLIB_LIBRARY})
    endif()
    # the static ParMETIS archive needs the C MPI library, placed after it (--as-needed)
    target_link_libraries(core PUBLIC ${MPI_mpi_LIBRARY})
    message(STATUS "subfv: CGNS reader with ParMETIS partitioning")
  else()
    message(STATUS "subfv: CGNS reader without ParMETIS (block partition only)")
  endif()
endif()

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

# ── Shared node-movement algorithm ────────────────────────────────────────────
# The shock-alignment node movement itself, with no solver state in it: both
# the standalone shock_adapt solver and euler_ho drive it, so it lives on its
# own rather than being duplicated (the curvature fit and the MPI vertex
# exchange are subtle enough that two copies would drift apart).
add_library(shock_adapt_move STATIC
  src/shock_adapt/shock_adapt_move_module.F90
)
target_link_libraries(shock_adapt_move PUBLIC core)
set_target_properties(shock_adapt_move PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/shock_adapt_move)
target_include_directories(shock_adapt_move PUBLIC ${CMAKE_BINARY_DIR}/mod/shock_adapt_move)

# ── Shock-adaptive node-movement library ──────────────────────────────────────
# Self-contained (own first-order face-loop Euler solver + 2-state Riemann
# fluxes, own copy, no dependency on ns/euler_ho/ale): minimal, local node
# movement near a detected shock, X-Mesh-style. Only links core.
add_library(shock_adapt STATIC
  src/shock_adapt/shock_adapt_global_data_module.F90
  src/shock_adapt/shock_adapt_module.F90
)
target_link_libraries(shock_adapt PUBLIC core shock_adapt_move)
set_target_properties(shock_adapt PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/shock_adapt)
target_include_directories(shock_adapt PUBLIC ${CMAKE_BINARY_DIR}/mod/shock_adapt)

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

# euler_ho can fold shock-alignment node movement into its time loop.
target_link_libraries(euler_ho PUBLIC shock_adapt_move)

# shock_adapt reuses arbitrary_high_order_module's Green-Gauss reconstruction
# for its own optional order=2 mode (still only pulls in core transitively,
# so shock_adapt's independence from ns/lagrange/euler_ho is preserved).
target_link_libraries(shock_adapt PUBLIC arbitrary_high_order)

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
add_executable(subfvshockadapt   src/shock_adapt/shock_adapt_main.F90)

target_link_libraries(subfvns          ns)
target_link_libraries(subfvshear       ns)
target_link_libraries(subfvacoustic    acoustic)
target_link_libraries(subfvlagrange    lagrange)
target_link_libraries(subfvgreshoconv  ns)
target_link_libraries(subfvrbf         rbf)
target_link_libraries(subfveulerho     euler_ho)
target_link_libraries(subfvahotest     arbitrary_high_order)
target_link_libraries(subfvale         ale)
target_link_libraries(subfvshockadapt  shock_adapt)

install(TARGETS subfvns subfvshear subfvacoustic subfvlagrange subfvgreshoconv subfvrbf subfvale subfvshockadapt
  RUNTIME DESTINATION bin
)

install(PROGRAMS tools/subfv-gmsh/subfv-gmsh
  DESTINATION bin
)

enable_testing()
add_subdirectory(test)
