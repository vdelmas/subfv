! High-order explicit Euler solver (finite volume, SSP-RK3 or RK4 in time).
module euler_ho_module
  use precision_module
  use mesh_module
  use quadrature_module
  use arbitrary_high_order_module, only: compute_next_order_derivative, &
    compute_next_order_derivative_cweno, apply_grad_bias_correction, &
    aho_set_wall_mirror_data => set_wall_mirror_data, &
    aho_module_mirror_wall_vec_start => mirror_wall_vec_start, &
    aho_module_use_weno_blend => use_weno_blend, &
    aho_module_use_green_gauss => use_green_gauss, &
    aho_module_use_cweno_center => use_cweno_center, &
    aho_module_cweno_center_weight => cweno_center_weight, &
    aho_module_cweno_center_power => cweno_center_power, &
    aho_module_eps_weight_num => eps_weight_num, &
    aho_module_eps_weight_num_deep => eps_weight_num_deep, &
    aho_module_eps_weight_num_gg => eps_weight_num_gg, &
    aho_module_eps_weight_num_deep_gg => eps_weight_num_deep_gg, &
    aho_module_weno_power => weno_power, &
    aho_module_grad_norm_derate => grad_norm_derate, &
    aho_module_grad_norm_derate_gg => grad_norm_derate_gg, &
    aho_module_use_alt_gg_weight => use_alt_gg_weight
  use mpi_module, only: mpi_send_recv_type, mpi_memory_exchange
  implicit none
  private

  real(kind=DOUBLE), parameter :: PI = 4.0_DOUBLE * datan(1.0_DOUBLE)

  integer(kind=ENTIER), parameter :: MAX_BC = 10
  character(len=255), public :: meshfile_path = ''
  character(len=255), public :: meshfile      = ''
  real(kind=DOUBLE),  public :: gamma_gas     = 1.4_DOUBLE
  ! Per-cell gamma for multi-material cases; gamma_gas stays the uniform-fill default.
  real(kind=DOUBLE), dimension(:), allocatable, public :: gamma_arr
  real(kind=DOUBLE),  public :: cfl           = 0.4_DOUBLE
  real(kind=DOUBLE),  public :: tmax          = 1.0_DOUBLE
  integer(kind=ENTIER), public :: order       = 1
  logical, public :: boundary_2d              = .false.
  integer(kind=ENTIER), public :: n_sol_vtu  = 10
  logical, public :: compute_error            = .false.
  logical, public :: error_2d                 = .true.
  ! Selects the reconstruction method: 0 = no aho reconstruction (ls_reconstruction
  ! / order 1, the pre-aho default); 1 = aho_gg (Green-Gauss nodal fit), the only
  ! supported aho variant -- aho_ls, aho_cls, and the pure cell-to-cell CLS variant
  ! were removed 2026-09-22 (aho_ls/aho_cls never beat aho_gg; the cell-to-cell CLS
  ! path was confirmed unstable in the real solver, see arbitrary-high-order memory).
  integer(kind=ENTIER), public :: aho_method  = 0
  ! Positivity safety net: if the linear (grad-only) reconstruction predicts negative rho or p
  ! at any of a cell's own vertices, zero that cell's grad/hess/third entirely (forces a plain
  ! first-order/constant state for that cell instead of extrapolating an unphysical one). See
  ! apply_positivity_kill. Off by default -- opt in per case.
  logical, public :: kill_recons              = .false.

  ! Shock-alignment node movement folded into the time loop (see
  ! shock_adapt_move_module). Off by default: a run with n_adapt_cycles = 0 is
  ! bit-for-bit the old fixed-mesh behaviour. The cycles are spaced a fixed
  ! number of ITERATIONS apart rather than waiting for a steady state, because
  ! this case never reaches one -- the wake stays unsteady at these
  ! resolutions, so "converged, then adapt" has no meaning here.
  integer(kind=ENTIER), public :: n_adapt_cycles     = 0
  integer(kind=ENTIER), public :: adapt_start_iter   = 2000
  integer(kind=ENTIER), public :: adapt_interval_iter = 300
  real(kind=DOUBLE), public :: adapt_grad_threshold  = 0.3_DOUBLE
  real(kind=DOUBLE), public :: adapt_max_move_frac   = 0.4_DOUBLE
  real(kind=DOUBLE), public :: adapt_relax           = 0.5_DOUBLE

  logical, public :: use_weno_blend           = .true.
  logical, public :: use_cweno_center          = .false.
  real(kind=DOUBLE), public :: cweno_center_weight = 1000.0_DOUBLE
  integer(kind=ENTIER), public :: cweno_center_power = 4
  logical, public :: use_grad_bias_correction  = .true.
  logical, public :: use_rk4 = .false.
  logical, public :: use_max_lambda_dt = .false.
  real(kind=DOUBLE), public :: eps_weight_num      = 1.0e-2_DOUBLE
  real(kind=DOUBLE), public :: eps_weight_num_deep = 1.0_DOUBLE
  real(kind=DOUBLE), public :: eps_weight_num_gg = 1.0e-2_DOUBLE
  real(kind=DOUBLE), public :: eps_weight_num_deep_gg = 1.0_DOUBLE
  integer(kind=ENTIER), public :: weno_power       = 1
  real(kind=DOUBLE), public :: grad_norm_derate    = 1.0e4_DOUBLE
  real(kind=DOUBLE), public :: grad_norm_derate_gg = 1.0e4_DOUBLE
  logical, public :: use_alt_gg_weight = .false.
  character(len=32), public :: flux_scheme    = 'three_wave'
  integer(kind=ENTIER), parameter :: FLUX_RUSANOV = 0, FLUX_TWO_WAVE = 1, FLUX_THREE_WAVE = 2, &
    FLUX_MODIFIED_THREE_WAVE = 3
  integer(kind=ENTIER) :: flux_scheme_id = FLUX_RUSANOV

  ! Per-face quadrature point/weight cache for face_flux_loop's order>=3 multi-point rule: face
  ! geometry is static across the whole time-marching run, so face_quad_pts's own point/weight
  ! generation (and the face_coords gather feeding it) only needs to run once per face, not once
  ! per face per RK stage per iteration. CSR-style ragged storage (offset + flat arrays), keyed by
  ! (n_faces, order) so a changed mesh or order forces a rebuild.
  real(kind=DOUBLE), dimension(:, :), allocatable, save :: face_quad_pts_cache
  real(kind=DOUBLE), dimension(:), allocatable, save :: face_quad_wts_cache
  integer(kind=ENTIER), dimension(:), allocatable, save :: face_quad_offset_cache
  integer(kind=ENTIER), save :: face_quad_cache_n_faces = -1
  integer(kind=ENTIER), save :: face_quad_cache_order = -1

  ! Persistent buffers for aho_reconstruction's packed grad/hess/third arrays -- previously
  ! allocated and deallocated fresh every single call (every RK stage; third_flat alone is
  ! 135*n_elems*8 bytes, tens of MB on a real 3D mesh), a major and entirely avoidable source of
  ! minor page faults (profiled: ~46 million over one run, System time comparable to User time).
  ! order is a namelist constant for the whole run, so which of these are actually needed never
  ! changes mid-run -- guarded on n_elems alone, same pattern as the WENO buffers in
  ! arbitrary_high_order_module.
  real(kind=DOUBLE), dimension(:, :), allocatable, save :: grad_flat_buf, hess_flat_buf, third_flat_buf
  integer(kind=ENTIER), save :: aho_flat_buf_n_elems = -1

  integer(kind=ENTIER), public :: n_bc = 0
  character(len=255), dimension(MAX_BC), public :: bc_name = ''
  character(len=255), dimension(MAX_BC), public :: bc_type = ''
  real(kind=DOUBLE), dimension(5, MAX_BC), public :: bc_val  = 0.0_DOUBLE
  integer(kind=ENTIER), parameter :: BC_WALL = 0, BC_FREESTREAM = 1, &
    BC_OUTFLOW = 2, BC_DMR_TOP = 3
  integer(kind=ENTIER), dimension(MAX_BC) :: bc_type_id = BC_WALL

  ! init: 0=uniform, 1=Sod 1D, 2=isentropic vortex, 3=DMR, 4=Shu-Osher,
  ! 5=Woodward-Colella, 6=shock-bubble (Haas & Sturtevant 1987).
  real(kind=DOUBLE), public :: u_bg_vortex  = 0.0_DOUBLE
  real(kind=DOUBLE), public :: v_bg_vortex  = 0.0_DOUBLE
  integer(kind=ENTIER), public :: init         = 0
  real(kind=DOUBLE), dimension(5), public :: sol_uniform    = [1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE]
  real(kind=DOUBLE), public :: x1drp          = 0.5_DOUBLE
  real(kind=DOUBLE), dimension(5), public :: sol_w_1drp_l   = [1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE]
  real(kind=DOUBLE), dimension(5), public :: sol_w_1drp_r   = [0.125_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.1_DOUBLE]

  public :: read_params
  public :: init_sol
  public :: compute_prim
  public :: compute_rhs
  public :: compute_dt
  public :: sync_gamma_arr
  public :: compute_error_vortex
  public :: test_reconstruction_exactness
  public :: invalidate_face_quad_cache
  public :: compute_cell_moments
  public :: reconstruct, ls_reconstruction, aho_reconstruction
  public :: setup_wall_mirror

  real(kind=DOUBLE), dimension(:, :, :), allocatable :: cell_moment
  real(kind=DOUBLE), dimension(:, :, :, :), allocatable :: cell_moment3

contains

  ! Registers the aho module's wall-mirror ghost-cell fix (see arbitrary_high_order_module's own
  ! cache-block header comment): for every boundary vertex whose touching boundary faces are ALL
  ! wall-type (never mixed with freestream/outflow/etc.), average their face-area-weighted outward
  ! normal (mesh%face(iface)%norm already points outward from the real cell for a boundary face --
  ! same convention ghost_prim relies on) into a unit wall normal. Call once per mesh, after
  ! compute_geometry_mesh and before the first aho_reconstruction call.
  subroutine setup_wall_mirror(mesh)
    implicit none

    type(mesh_type), intent(in) :: mesh

    integer(kind=ENTIER) :: iface, ir, id_bc, iv, v
    real(kind=DOUBLE), dimension(:, :), allocatable :: norm_acc
    real(kind=DOUBLE), dimension(:, :), allocatable :: wall_norm_v
    logical, dimension(:), allocatable :: touches_wall, touches_nonwall, wall_valid_v
    real(kind=DOUBLE), dimension(3) :: n

    allocate(norm_acc(3, mesh%n_vert))
    allocate(touches_wall(mesh%n_vert), touches_nonwall(mesh%n_vert))
    norm_acc = 0.0_DOUBLE
    touches_wall = .false.
    touches_nonwall = .false.

    do iface = 1, mesh%n_faces
      ir = mesh%face(iface)%right_neigh
      if (ir >= 0) cycle
      id_bc = -ir
      if (id_bc < 1 .or. id_bc > n_bc) cycle
      if (.not. allocated(mesh%face(iface)%vert)) cycle
      ! z-faces on a boundary_2d extrusion carry no physics (same guard as face_flux_loop's is_zface).
      if (boundary_2d .and. abs(abs(mesh%face(iface)%norm(3)) - 1.0_DOUBLE) < 1.0e-6_DOUBLE) cycle
      do iv = 1, mesh%face(iface)%n_vert
        v = mesh%face(iface)%vert(iv)
        if (bc_type_id(id_bc) == BC_WALL) then
          norm_acc(:, v) = norm_acc(:, v) + mesh%face(iface)%area * mesh%face(iface)%norm
          touches_wall(v) = .true.
        else
          touches_nonwall(v) = .true.
        end if
      end do
    end do

    allocate(wall_norm_v(3, mesh%n_vert), wall_valid_v(mesh%n_vert))
    wall_norm_v = 0.0_DOUBLE
    wall_valid_v = .false.
    do v = 1, mesh%n_vert
      if (.not. mesh%vert(v)%is_bound) cycle
      if (.not. touches_wall(v) .or. touches_nonwall(v)) cycle
      n = norm_acc(:, v)
      if (dot_product(n, n) < 1.0e-24_DOUBLE) cycle
      wall_norm_v(:, v) = n / sqrt(dot_product(n, n))
      wall_valid_v(v) = .true.
    end do

    aho_module_mirror_wall_vec_start = 2_ENTIER
    call aho_set_wall_mirror_data(mesh, wall_norm_v, wall_valid_v)

    deallocate(norm_acc, touches_wall, touches_nonwall, wall_norm_v, wall_valid_v)
  end subroutine setup_wall_mirror

  subroutine read_params(filename)
    implicit none

    character(len=*), intent(in) :: filename

    integer(kind=ENTIER) :: funit, i_bc
    namelist /INPUT_PARAM/ &
      meshfile_path, meshfile, &
      n_bc, bc_name, bc_type, bc_val, &
      boundary_2d, &
      init, sol_uniform, x1drp, sol_w_1drp_l, sol_w_1drp_r, &
      u_bg_vortex, v_bg_vortex, &
      gamma_gas, order, cfl, tmax, n_sol_vtu, &
      compute_error, error_2d, aho_method, kill_recons, use_weno_blend, &
      use_cweno_center, cweno_center_weight, cweno_center_power, use_grad_bias_correction, use_rk4, &
      use_max_lambda_dt, flux_scheme, &
      eps_weight_num, eps_weight_num_deep, weno_power, grad_norm_derate, &
      grad_norm_derate_gg, use_alt_gg_weight, eps_weight_num_gg, eps_weight_num_deep_gg, &
      n_adapt_cycles, adapt_start_iter, adapt_interval_iter, &
      adapt_grad_threshold, adapt_max_move_frac, adapt_relax

    open(newunit=funit, file=trim(adjustl(filename)))
    read(nml=INPUT_PARAM, unit=funit)
    close(funit)

    select case (trim(adjustl(flux_scheme)))
    case ('three_wave')
      flux_scheme_id = FLUX_THREE_WAVE
    case ('modified_three_wave')
      flux_scheme_id = FLUX_MODIFIED_THREE_WAVE
    case ('two_wave')
      flux_scheme_id = FLUX_TWO_WAVE
    case default
      flux_scheme_id = FLUX_RUSANOV
    end select

    do i_bc = 1, n_bc
      select case (trim(adjustl(bc_type(i_bc))))
      case ('freestream')
        bc_type_id(i_bc) = BC_FREESTREAM
      case ('outflowsupersonic', 'outflow')
        bc_type_id(i_bc) = BC_OUTFLOW
      case ('dmr_top')
        bc_type_id(i_bc) = BC_DMR_TOP
      case default
        bc_type_id(i_bc) = BC_WALL
      end select
    end do
  end subroutine read_params

  subroutine init_sol(mesh, sol)
    implicit none

    type(mesh_type), intent(in)    :: mesh
    ! sol(6,:) = rho*Gamma, Gamma=1/(gamma-1), the Abgrall (1996) gamma-transport variable.
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(out) :: sol

    integer(kind=ENTIER) :: i
    real(kind=DOUBLE), dimension(5) :: w
    real(kind=DOUBLE), dimension(3) :: xc
    real(kind=DOUBLE) :: rb
    real(kind=DOUBLE), parameter :: xc_bub    = 0.350_DOUBLE
    real(kind=DOUBLE), parameter :: yc_bub    = 0.0445_DOUBLE
    real(kind=DOUBLE), parameter :: r_bub     = 0.025_DOUBLE
    real(kind=DOUBLE), parameter :: x_shock   = 0.400_DOUBLE
    real(kind=DOUBLE), parameter :: gamma_air = 1.4_DOUBLE
    real(kind=DOUBLE), parameter :: gamma_bub = 1.648_DOUBLE
    real(kind=DOUBLE), parameter :: rho_bub   = 0.287_DOUBLE / 1.578_DOUBLE
    real(kind=DOUBLE), parameter :: rho_post  = 1.376364_DOUBLE
    real(kind=DOUBLE), parameter :: p_post    = 1.569800_DOUBLE
    real(kind=DOUBLE), parameter :: vx_post   = -0.394729_DOUBLE
    ! Triple-point shock interaction (init=7), canonical three-material setup.
    real(kind=DOUBLE), parameter :: x_tp      = 1.0_DOUBLE
    real(kind=DOUBLE), parameter :: y_tp      = 1.5_DOUBLE
    real(kind=DOUBLE), parameter :: gamma_tp1 = 1.5_DOUBLE
    real(kind=DOUBLE), parameter :: gamma_tp2 = 1.4_DOUBLE
    real(kind=DOUBLE), parameter :: gamma_tp3 = 1.5_DOUBLE

    if (.not. allocated(gamma_arr)) allocate(gamma_arr(mesh%n_elems))
    gamma_arr = gamma_gas

    do i = 1, mesh%n_elems
      xc = mesh%elem(i)%coord
      select case (init)
      case (0)
        w = sol_uniform
      case (1)
        if (xc(1) < x1drp) then
          w = sol_w_1drp_l
        else
          w = sol_w_1drp_r
        end if
      case (3)
        w = dmr_state(xc(1), xc(2), 0.0_DOUBLE)
      case (4)
        w = shu_osher_state(xc(1))
      case (5)
        w = woodward_colella_state(xc(1))
      case (2)
        call vortex_cell_average(mesh, i, 0.0_DOUBLE, w)
      case (6)
        rb = sqrt((xc(1) - xc_bub)**2 + (xc(2) - yc_bub)**2)
        if (xc(1) > x_shock) then
          gamma_arr(i) = gamma_air
          w = [rho_post, vx_post, 0.0_DOUBLE, 0.0_DOUBLE, p_post]
        else if (rb <= r_bub) then
          gamma_arr(i) = gamma_bub
          w = [rho_bub, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE]
        else
          gamma_arr(i) = gamma_air
          w = [1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE]
        end if
      case (7)
        ! Triple-point shock interaction on [0,7]x[0,3], fluid at rest, reflecting
        ! walls all around (Galera/Maire/Breil, JCP 2010). Three materials:
        !   x <= 1          : rho=1,     p=1,   gamma=1.5  (driver)
        !   x > 1, y <= 1.5 : rho=1,     p=0.1, gamma=1.4  (dense, slow shock)
        !   x > 1, y > 1.5  : rho=0.125, p=0.1, gamma=1.5  (light, fast shock)
        ! The shock outruns itself across y=1.5, and the resulting shear layer
        ! rolls up into the vortex at the triple point.
        if (xc(1) <= x_tp) then
          gamma_arr(i) = gamma_tp1
          w = [1.0_DOUBLE,   0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE]
        else if (xc(2) <= y_tp) then
          gamma_arr(i) = gamma_tp2
          w = [1.0_DOUBLE,   0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.1_DOUBLE]
        else
          gamma_arr(i) = gamma_tp3
          w = [0.125_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.1_DOUBLE]
        end if
      case default
        w = sol_uniform
      end select
      sol(1:5, i) = primit_to_conserv(w, gamma_arr(i))
      sol(6, i)   = sol(1, i) / (gamma_arr(i) - 1.0_DOUBLE)
    end do
  end subroutine init_sol

  subroutine compute_prim(mesh, sol, prim)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(in)  :: sol
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(out) :: prim

    integer(kind=ENTIER) :: i

    do i = 1, mesh%n_elems
      prim(:, i) = conserv_to_primit(sol(1:5, i), gamma_arr(i))
    end do
  end subroutine compute_prim

  subroutine sync_gamma_arr(mesh, sol)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(in) :: sol

    integer(kind=ENTIER) :: i

    do i = 1, mesh%n_elems
      gamma_arr(i) = 1.0_DOUBLE + sol(1, i) / sol(6, i)
    end do
  end subroutine sync_gamma_arr

  subroutine compute_rhs(mesh, sol, prim, rhs, sum_lambda, t, num_procs, mpi_send_recv)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(in)  :: sol
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in)  :: prim
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(out) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems),    intent(out) :: sum_lambda
    real(kind=DOUBLE), intent(in) :: t
    integer(kind=ENTIER), intent(in), optional :: num_procs
    type(mpi_send_recv_type), intent(inout), optional :: mpi_send_recv

    real(kind=DOUBLE), allocatable :: grad(:, :, :)
    real(kind=DOUBLE), allocatable :: hess(:, :, :, :)
    real(kind=DOUBLE), allocatable :: third(:, :, :, :, :)
    integer(kind=ENTIER) :: i

    allocate(grad(5, 3, mesh%n_elems))
    grad = 0.0_DOUBLE

    if (order >= 3) then
      allocate(hess(5, 3, 3, mesh%n_elems))
      hess = 0.0_DOUBLE
    end if

    if (aho_method /= 0 .and. aho_method /= 1) then
      print *, 'FATAL: aho_method must be 0 (off, ls_reconstruction) or 1 (aho_gg) -- ', &
        'aho_ls/aho_cls/the cell-to-cell CLS variant have been removed'
      error stop 1
    end if

    if (order >= 4) then
      if (aho_method == 0) then
        print *, 'FATAL: order>=4 requires aho_method = 1 (aho_gg)'
        error stop 1
      end if
      allocate(third(5, 3, 3, 3, mesh%n_elems))
      third = 0.0_DOUBLE
    end if

    if (order >= 2) then
      if (aho_method /= 0) then
        call aho_reconstruction(mesh, prim, grad, hess, third, num_procs, mpi_send_recv)
      else
        call ls_reconstruction(mesh, prim, grad, hess)
      end if
    end if

    rhs        = 0.0_DOUBLE
    sum_lambda = 0.0_DOUBLE
    call face_flux_loop(mesh, sol, prim, grad, hess, rhs, sum_lambda, t, third)

    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) rhs(:, i) = rhs(:, i) / mesh%elem(i)%volume
    end do

    if (allocated(third)) deallocate(third)
    if (allocated(hess)) deallocate(hess)
    deallocate(grad)
  end subroutine compute_rhs

  function compute_dt(mesh, sum_lambda) result(dt)
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: sum_lambda

    real(kind=DOUBLE) :: dt, dt_i
    integer(kind=ENTIER) :: i
    integer :: mpi_ierr

    dt = huge(1.0_DOUBLE)
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost .and. sum_lambda(i) > 0.0_DOUBLE) then
        dt_i = cfl * mesh%elem(i)%volume / sum_lambda(i)
        dt = min(dt, dt_i)
      end if
    end do

    call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_DOUBLE, &
      MPI_MIN, MPI_COMM_WORLD, mpi_ierr)
  end function compute_dt

  ! Checks reconstruct()'s k-exactness on a synthetic degree-1/2 polynomial field;
  ! order_override lets the caller deliberately mismatch order vs. degree as a negative control.
  subroutine test_reconstruction_exactness(mesh, degree, dx, dy, half_width, order_override)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: degree
    real(kind=DOUBLE), intent(in) :: dx, dy, half_width
    integer(kind=ENTIER), intent(in), optional :: order_override

    real(kind=DOUBLE), dimension(5, mesh%n_elems) :: prim
    real(kind=DOUBLE), dimension(5, 3, mesh%n_elems) :: grad_ls, grad_aho
    real(kind=DOUBLE), dimension(:, :, :, :), allocatable :: hess_ls, hess_aho
    real(kind=DOUBLE), dimension(5) :: A, B, C, D, E, F
    real(kind=DOUBLE), dimension(3) :: xq, xc_vec
    real(kind=DOUBLE), dimension(5) :: wexact, w_ls, w_aho
    real(kind=DOUBLE) :: xc, yc, err_ls, err_aho
    integer(kind=ENTIER) :: i, v, s, saved_order
    real(kind=DOUBLE), dimension(4, 2) :: offsets

    saved_order = order
    if (present(order_override)) then
      order = order_override
    else
      order = degree + 1
    end if

    A = [100.0_DOUBLE, -0.7_DOUBLE, 2.3_DOUBLE, 0.4_DOUBLE, 100.0_DOUBLE]
    B = [0.07_DOUBLE, 0.13_DOUBLE, -0.05_DOUBLE, 0.09_DOUBLE, 0.02_DOUBLE]
    C = [-0.03_DOUBLE, 0.06_DOUBLE, 0.11_DOUBLE, -0.08_DOUBLE, 0.04_DOUBLE]
    D = [0.013_DOUBLE, -0.009_DOUBLE, 0.021_DOUBLE, 0.005_DOUBLE, -0.017_DOUBLE]
    E = [-0.021_DOUBLE, 0.011_DOUBLE, -0.008_DOUBLE, 0.019_DOUBLE, 0.006_DOUBLE]
    F = [0.017_DOUBLE, 0.008_DOUBLE, -0.014_DOUBLE, -0.006_DOUBLE, 0.012_DOUBLE]

    do i = 1, mesh%n_elems
      xc = mesh%elem(i)%coord(1)
      yc = mesh%elem(i)%coord(2)
      do v = 1, 5
        prim(v, i) = A(v) + B(v)*xc + C(v)*yc
        if (degree >= 2) then
          prim(v, i) = prim(v, i) + D(v)*(xc**2 + dx**2/12.0_DOUBLE) &
                                   + E(v)*(xc*yc) &
                                   + F(v)*(yc**2 + dy**2/12.0_DOUBLE)
        end if
      end do
    end do

    if (order >= 3) then
      allocate(hess_ls(5, 3, 3, mesh%n_elems), hess_aho(5, 3, 3, mesh%n_elems))
      hess_ls = 0.0_DOUBLE; hess_aho = 0.0_DOUBLE
    end if
    grad_ls = 0.0_DOUBLE; grad_aho = 0.0_DOUBLE

    call ls_reconstruction(mesh, prim, grad_ls, hess_ls)
    call aho_reconstruction(mesh, prim, grad_aho, hess_aho)

    offsets(1, :) = [ 0.30_DOUBLE,  0.20_DOUBLE]
    offsets(2, :) = [-0.40_DOUBLE,  0.10_DOUBLE]
    offsets(3, :) = [ 0.15_DOUBLE, -0.35_DOUBLE]
    offsets(4, :) = [-0.20_DOUBLE, -0.25_DOUBLE]

    err_ls = 0.0_DOUBLE; err_aho = 0.0_DOUBLE

    do i = 1, mesh%n_elems
      if (mesh%elem(i)%is_ghost) cycle
      xc = mesh%elem(i)%coord(1)
      yc = mesh%elem(i)%coord(2)
      if (abs(xc) > 0.6_DOUBLE*half_width .or. abs(yc) > 0.6_DOUBLE*half_width) cycle
      xc_vec = mesh%elem(i)%coord

      do s = 1, 4
        xq(1) = xc + offsets(s, 1)*dx
        xq(2) = yc + offsets(s, 2)*dy
        xq(3) = xc_vec(3)

        do v = 1, 5
          wexact(v) = A(v) + B(v)*xq(1) + C(v)*xq(2)
          if (degree >= 2) then
            wexact(v) = wexact(v) + D(v)*xq(1)**2 + E(v)*xq(1)*xq(2) + F(v)*xq(2)**2
          end if
        end do

        w_ls  = reconstruct(prim, grad_ls,  hess_ls,  i, xq, xc_vec, gamma_gas)
        w_aho = reconstruct(prim, grad_aho, hess_aho, i, xq, xc_vec, gamma_gas)

        err_ls  = max(err_ls,  maxval(abs(w_ls  - wexact)))
        err_aho = max(err_aho, maxval(abs(w_aho - wexact)))
      end do
    end do

    print *, 'RECON_EXACTNESS degree=', degree, ' order=', order, &
      ' Linf(ls)=', err_ls, ' Linf(aho)=', err_aho

    if (allocated(hess_ls))  deallocate(hess_ls)
    if (allocated(hess_aho)) deallocate(hess_aho)
    order = saved_order
  end subroutine test_reconstruction_exactness

  subroutine compute_error_vortex(mesh, sol, t, h, l2err)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), intent(in) :: t
    real(kind=DOUBLE), intent(out) :: h, l2err

    integer(kind=ENTIER) :: i, kf, iface_loc, n_fv, k, n_inner
    real(kind=DOUBLE) :: err2, area_tot, area_q, area_2d
    real(kind=DOUBLE), dimension(5) :: wexact, wsol
    real(kind=DOUBLE), dimension(3) :: xc
    real(kind=DOUBLE), dimension(:, :), allocatable :: qpts_v, fc
    real(kind=DOUBLE), dimension(:),    allocatable :: qwts_v
    logical :: found_zface

    err2     = 0.0_DOUBLE
    area_tot = 0.0_DOUBLE
    n_inner  = 0

    do i = 1, mesh%n_elems
      if (mesh%elem(i)%is_ghost) cycle
      xc = mesh%elem(i)%coord
      if (abs(xc(1)) > 5.0_DOUBLE .or. abs(xc(2)) > 5.0_DOUBLE) cycle

      call vortex_cell_average(mesh, i, t, wexact)

      found_zface = .false.
      do kf = 1, mesh%elem(i)%n_faces
        iface_loc = mesh%elem(i)%face(kf)
        if (abs(abs(mesh%face(iface_loc)%norm(3)) - 1.0_DOUBLE) < 1.0e-6_DOUBLE) then
          n_fv = mesh%face(iface_loc)%n_vert
          if (n_fv >= 3 .and. allocated(mesh%face(iface_loc)%vert)) then
            allocate(fc(3, n_fv))
            do k = 1, n_fv
              fc(:, k) = mesh%vert(mesh%face(iface_loc)%vert(k))%coord
            end do
            call face_quad_pts(n_fv, fc, max(order, 3), qpts_v, qwts_v)
            deallocate(fc)
            area_q = sum(qwts_v)
            deallocate(qpts_v, qwts_v)
            found_zface = .true.
            exit
          end if
        end if
      end do

      wsol = conserv_to_primit(sol(1:5, i), gamma_arr(i))
      if (found_zface .and. area_q > 0.0_DOUBLE) then
        area_2d = area_q
      else
        area_2d = mesh%elem(i)%volume**(2.0_DOUBLE/3.0_DOUBLE)
      end if
      err2     = err2     + area_2d * (wsol(1) - wexact(1))**2
      area_tot = area_tot + area_2d
      n_inner  = n_inner  + 1
    end do

    call mpi_sum_error_vortex(err2, area_tot, n_inner)

    if (n_inner > 0) then
      h = sqrt(area_tot / n_inner)
    else
      h = 0.0_DOUBLE
    end if
    l2err = sqrt(err2)
  end subroutine compute_error_vortex

  subroutine mpi_sum_error_vortex(err2, area_tot, n_inner)
    use mpi
    implicit none

    real(kind=DOUBLE), intent(inout) :: err2, area_tot
    integer(kind=ENTIER), intent(inout) :: n_inner

    real(kind=DOUBLE) :: buf_real(2)
    integer(kind=ENTIER) :: n_inner_sum
    integer :: mpi_ierr

    call MPI_ALLREDUCE((/err2, area_tot/), buf_real, 2, MPI_DOUBLE_PRECISION, &
      MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
    err2     = buf_real(1)
    area_tot = buf_real(2)

    call MPI_ALLREDUCE(n_inner, n_inner_sum, 1, MPI_INTEGER, MPI_SUM, &
      MPI_COMM_WORLD, mpi_ierr)
    n_inner = n_inner_sum
  end subroutine mpi_sum_error_vortex

  ! ================================================================
  ! Internal subroutines
  ! ================================================================

  ! Cell-centred weighted least-squares polynomial reconstruction over vertex-neighbours.
  ! Tops out at order 3 (no cubic basis); order>=4 requires aho_method /= 0.
  subroutine ls_reconstruction(mesh, prim, grad, hess)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: prim
    real(kind=DOUBLE), dimension(5, 3, mesh%n_elems), intent(inout) :: grad
    real(kind=DOUBLE), dimension(:, :, :, :), allocatable, intent(inout) :: hess

    integer(kind=ENTIER) :: i, iv, mf, iface, il, ir, jn, kv, nv
    integer(kind=ENTIER) :: n_coeff, n_neigh, n_needed
    integer(kind=ENTIER) :: ring_start, ring_end, ic, ivv
    integer(kind=ENTIER), dimension(500) :: neigh_list
    real(kind=DOUBLE), allocatable :: Amat(:, :)
    real(kind=DOUBLE), dimension(500) :: wt_arr
    real(kind=DOUBLE), dimension(9, 9) :: ATA
    real(kind=DOUBLE), dimension(9)    :: ATb, x_coeff
    real(kind=DOUBLE) :: ddx, ddy, ddz, dist, wt, Lref, sdx, sdy, sdz
    real(kind=DOUBLE), dimension(3) :: xci, xcj

    if (boundary_2d) then
      n_coeff = merge(5, 2, order >= 3)
    else
      n_coeff = merge(9, 3, order >= 3)
    end if

    do i = 1, mesh%n_elems
      if (mesh%elem(i)%is_ghost) cycle

      xci    = mesh%elem(i)%coord
      n_neigh = 0

      do nv = 1, mesh%elem(i)%n_vert
        iv = mesh%elem(i)%vert(nv)
        do mf = 1, mesh%vert(iv)%n_faces_neigh
          iface = mesh%vert(iv)%face_neigh(mf)
          call add_unique_neighbor(mesh%face(iface)%left_neigh)
          call add_unique_neighbor(mesh%face(iface)%right_neigh)
        end do
      end do

      ! Ring-expand until enough neighbours to determine the unknowns actually alive
      ! on this mesh (e.g. a single-row-in-y mesh has fewer live directions).
      ring_start = 1
      do
        n_needed = count_needed_coeffs()
        if (n_neigh >= n_needed .or. n_neigh >= 500) exit
        ring_end = n_neigh
        do ic = ring_start, ring_end
          do ivv = 1, mesh%elem(neigh_list(ic))%n_vert
            iv = mesh%elem(neigh_list(ic))%vert(ivv)
            do mf = 1, mesh%vert(iv)%n_faces_neigh
              iface = mesh%vert(iv)%face_neigh(mf)
              call add_unique_neighbor(mesh%face(iface)%left_neigh)
              call add_unique_neighbor(mesh%face(iface)%right_neigh)
            end do
          end do
        end do
        if (ring_end == n_neigh) exit
        ring_start = ring_end + 1
      end do

      if (n_neigh < n_needed) cycle

      allocate(Amat(n_neigh, n_coeff))

      ! Lref nondimensionalizes design-matrix columns so linear/quadratic terms share a
      ! scale -- needed for a well-conditioned solve on a fine mesh.
      Lref = max(mesh%elem(i)%volume, 1.0e-300_DOUBLE)**(1.0_DOUBLE/3.0_DOUBLE)

      do jn = 1, n_neigh
        xcj = mesh%elem(neigh_list(jn))%coord
        ddx = xcj(1) - xci(1)
        ddy = xcj(2) - xci(2)
        ddz = xcj(3) - xci(3)
        if (boundary_2d) then
          dist = sqrt(ddx**2 + ddy**2)
        else
          dist = sqrt(ddx**2 + ddy**2 + ddz**2)
        end if
        wt = 1.0_DOUBLE / sqrt(max(dist, 1.0e-14_DOUBLE))
        wt_arr(jn) = wt
        sdx = ddx / Lref
        sdy = ddy / Lref
        sdz = ddz / Lref

        if (boundary_2d) then
          Amat(jn, 1) = sdx * wt
          Amat(jn, 2) = sdy * wt
          if (order >= 3) then
            Amat(jn, 3) = sdx * sdx * 0.5_DOUBLE * wt
            Amat(jn, 4) = sdx * sdy             * wt
            Amat(jn, 5) = sdy * sdy * 0.5_DOUBLE * wt
          end if
        else
          Amat(jn, 1) = sdx * wt
          Amat(jn, 2) = sdy * wt
          Amat(jn, 3) = sdz * wt
          if (order >= 3) then
            Amat(jn, 4) = sdx * sdx * 0.5_DOUBLE * wt
            Amat(jn, 5) = sdx * sdy             * wt
            Amat(jn, 6) = sdx * sdz             * wt
            Amat(jn, 7) = sdy * sdy * 0.5_DOUBLE * wt
            Amat(jn, 8) = sdy * sdz             * wt
            Amat(jn, 9) = sdz * sdz * 0.5_DOUBLE * wt
          end if
        end if
      end do

      ATA(1:n_coeff, 1:n_coeff) = matmul(transpose(Amat), Amat)

      do kv = 1, 5
        ATb(1:n_coeff) = 0.0_DOUBLE
        do jn = 1, n_neigh
          ATb(1:n_coeff) = ATb(1:n_coeff) + Amat(jn, :) * wt_arr(jn) * &
            (prim(kv, neigh_list(jn)) - prim(kv, i))
        end do

        call solve_ls(ATA, ATb, x_coeff, n_coeff)

        grad(kv, 1, i) = x_coeff(1) / Lref
        grad(kv, 2, i) = x_coeff(2) / Lref
        if (boundary_2d) then
          grad(kv, 3, i) = 0.0_DOUBLE
        else
          grad(kv, 3, i) = x_coeff(3) / Lref
        end if

        if (order >= 3 .and. allocated(hess)) then
          hess(kv, :, :, i) = 0.0_DOUBLE
          if (boundary_2d) then
            hess(kv, 1, 1, i) = x_coeff(3) / Lref**2
            hess(kv, 1, 2, i) = x_coeff(4) / Lref**2
            hess(kv, 2, 1, i) = x_coeff(4) / Lref**2
            hess(kv, 2, 2, i) = x_coeff(5) / Lref**2
          else
            hess(kv, 1, 1, i) = x_coeff(4) / Lref**2
            hess(kv, 1, 2, i) = x_coeff(5) / Lref**2
            hess(kv, 2, 1, i) = x_coeff(5) / Lref**2
            hess(kv, 1, 3, i) = x_coeff(6) / Lref**2
            hess(kv, 3, 1, i) = x_coeff(6) / Lref**2
            hess(kv, 2, 2, i) = x_coeff(7) / Lref**2
            hess(kv, 2, 3, i) = x_coeff(8) / Lref**2
            hess(kv, 3, 2, i) = x_coeff(8) / Lref**2
            hess(kv, 3, 3, i) = x_coeff(9) / Lref**2
          end if
        end if

      end do

      deallocate(Amat)
    end do

  contains

    subroutine add_unique_neighbor(cand)
      implicit none

      integer(kind=ENTIER), intent(in) :: cand

      integer(kind=ENTIER) :: jn
      logical :: found

      if (cand <= 0) return
      if (cand == i .or. mesh%elem(cand)%is_ghost) return
      found = .false.
      do jn = 1, n_neigh
        if (neigh_list(jn) == cand) then
          found = .true.
          exit
        end if
      end do
      if (.not. found .and. n_neigh < 500) then
        n_neigh = n_neigh + 1
        neigh_list(n_neigh) = cand
      end if
    end subroutine add_unique_neighbor

    function count_needed_coeffs() result(n_needed)
      implicit none

      integer(kind=ENTIER) :: n_needed

      real(kind=DOUBLE) :: spread_x, spread_y, spread_z
      real(kind=DOUBLE) :: ddxj, ddyj, ddzj
      logical :: y_alive, z_alive
      integer(kind=ENTIER) :: jn2, n_lin

      spread_x = 0.0_DOUBLE
      spread_y = 0.0_DOUBLE
      spread_z = 0.0_DOUBLE
      do jn2 = 1, n_neigh
        ddxj = mesh%elem(neigh_list(jn2))%coord(1) - xci(1)
        ddyj = mesh%elem(neigh_list(jn2))%coord(2) - xci(2)
        ddzj = mesh%elem(neigh_list(jn2))%coord(3) - xci(3)
        spread_x = max(spread_x, abs(ddxj))
        spread_y = max(spread_y, abs(ddyj))
        spread_z = max(spread_z, abs(ddzj))
      end do

      y_alive = spread_y >= 1.0e-8_DOUBLE * max(spread_x, 1.0e-300_DOUBLE)
      z_alive = (.not. boundary_2d) .and. &
        (spread_z >= 1.0e-8_DOUBLE * max(spread_x, 1.0e-300_DOUBLE))

      n_lin = 1
      if (y_alive) n_lin = n_lin + 1
      if (z_alive) n_lin = n_lin + 1

      n_needed = n_lin
      if (order >= 3) n_needed = n_needed + n_lin * (n_lin + 1) / 2
    end function count_needed_coeffs

  end subroutine ls_reconstruction

  ! Drop-in alternative to ls_reconstruction: grad/hess/third of the 5 primitives from
  ! arbitrary_high_order_module's dual/primal hierarchy, batched as one nc_in=5 field.
  ! Exchanges grad/hess/third across MPI ranks after each stage since a ghost cell's own
  ! vertex-neighbour stencil can be incomplete one hop in from the partition interface.
  subroutine aho_reconstruction(mesh, prim, grad, hess, third, num_procs, mpi_send_recv)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: prim
    real(kind=DOUBLE), dimension(5, 3, mesh%n_elems), intent(inout) :: grad
    real(kind=DOUBLE), dimension(:, :, :, :), allocatable, intent(inout) :: hess
    real(kind=DOUBLE), dimension(:, :, :, :, :), allocatable, intent(inout), optional :: third
    integer(kind=ENTIER), intent(in), optional :: num_procs
    type(mpi_send_recv_type), intent(inout), optional :: mpi_send_recv

    integer(kind=ENTIER) :: e, v, dir1, dir2, dir3
    logical :: do_exchange
    real(kind=DOUBLE), dimension(:, :), allocatable :: hess_v, third_v
    logical, dimension(:), allocatable :: valid_hess_v, valid_third_v
    real(kind=DOUBLE), dimension(:), allocatable :: grad_oi_v
    real(kind=DOUBLE), dimension(:, :), allocatable :: hess_v_local, third_v_local
    logical, dimension(:), allocatable :: valid_hess_v_local, valid_third_v_local
    real(kind=DOUBLE), dimension(:, :), allocatable :: grad_v
    logical, dimension(:), allocatable :: valid_grad_v
    logical :: need_hess_v

    do_exchange = present(num_procs)
    if (do_exchange) do_exchange = num_procs > 1

    aho_module_use_weno_blend = use_weno_blend
    aho_module_use_green_gauss = .true.
    aho_module_use_cweno_center = use_cweno_center
    aho_module_cweno_center_weight = cweno_center_weight
    aho_module_cweno_center_power = cweno_center_power
    aho_module_eps_weight_num = eps_weight_num
    aho_module_eps_weight_num_deep = eps_weight_num_deep
    aho_module_weno_power = weno_power
    aho_module_grad_norm_derate = grad_norm_derate
    aho_module_grad_norm_derate_gg = grad_norm_derate_gg
    aho_module_use_alt_gg_weight = use_alt_gg_weight
    aho_module_eps_weight_num_gg = eps_weight_num_gg
    aho_module_eps_weight_num_deep_gg = eps_weight_num_deep_gg

    if (aho_flat_buf_n_elems /= mesh%n_elems) then
      if (allocated(grad_flat_buf)) deallocate(grad_flat_buf, hess_flat_buf, third_flat_buf)
      ! Always sized to their full (order-4) extent regardless of the run's own order -- a fixed,
      ! one-time ~50MB for a 32k-element mesh, trivial next to what it replaces (reallocating up
      ! to the 135-row buffer fresh every RK stage).
      allocate(grad_flat_buf(15, mesh%n_elems))
      allocate(hess_flat_buf(45, mesh%n_elems))
      allocate(third_flat_buf(135, mesh%n_elems))
      aho_flat_buf_n_elems = mesh%n_elems
    end if

    associate (grad_flat => grad_flat_buf, hess_flat => hess_flat_buf, third_flat => third_flat_buf)

    allocate(grad_oi_v(mesh%n_vert))
    grad_oi_v = 0.0_DOUBLE
    if (use_cweno_center) then
      call compute_next_order_derivative_cweno(mesh, 3_ENTIER, 5_ENTIER, boundary_2d, &
        prim, grad_flat)
    else if (order >= 4 .and. present(third) .and. use_grad_bias_correction) then
      call compute_next_order_derivative(mesh, 3_ENTIER, 5_ENTIER, boundary_2d, &
        prim, grad_flat, deriv_order=1_ENTIER, oi_v_out=grad_oi_v, &
        dphi_v_out=grad_v, valid_v_out=valid_grad_v)
    else
      call compute_next_order_derivative(mesh, 3_ENTIER, 5_ENTIER, boundary_2d, &
        prim, grad_flat, deriv_order=1_ENTIER, oi_v_out=grad_oi_v)
    end if

    if (do_exchange) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 15_ENTIER, grad_flat)

    ! grad_flat is unpacked into grad only after the order-4 bias correction below.
    if (order >= 3 .and. allocated(hess)) then
      if (use_cweno_center) then
        call compute_next_order_derivative_cweno(mesh, 3_ENTIER, 15_ENTIER, boundary_2d, &
          grad_flat, hess_flat)
      else
        need_hess_v = (order >= 4 .and. present(third) .and. use_grad_bias_correction)
        if (need_hess_v) then
          call compute_next_order_derivative(mesh, 3_ENTIER, 15_ENTIER, boundary_2d, &
            grad_flat, hess_flat, deriv_order=2_ENTIER, &
            dphi_v_out=hess_v_local, valid_v_out=valid_hess_v_local)
          if (allocated(hess_v)) deallocate(hess_v)
          if (allocated(valid_hess_v)) deallocate(valid_hess_v)
          allocate(hess_v(size(hess_v_local,1), size(hess_v_local,2)))
          allocate(valid_hess_v(size(valid_hess_v_local)))
          hess_v = hess_v_local
          valid_hess_v = valid_hess_v_local
        else
          call compute_next_order_derivative(mesh, 3_ENTIER, 15_ENTIER, boundary_2d, &
            grad_flat, hess_flat, deriv_order=2_ENTIER)
        end if
      end if
      if (do_exchange) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 45_ENTIER, hess_flat)

      do e = 1, mesh%n_elems
        do dir2 = 1, 3
          do dir1 = 1, 3
            do v = 1, 5
              hess(v, dir1, dir2, e) = hess_flat((dir2-1)*15 + (dir1-1)*5 + v, e)
            end do
          end do
        end do
      end do

      if (order >= 4 .and. present(third)) then
        if (allocated(third)) then
          if (use_cweno_center) then
            call compute_next_order_derivative_cweno(mesh, 3_ENTIER, 45_ENTIER, boundary_2d, &
              hess_flat, third_flat)
          else
            if (use_grad_bias_correction) then
              call compute_next_order_derivative(mesh, 3_ENTIER, 45_ENTIER, boundary_2d, &
                hess_flat, third_flat, deriv_order=3_ENTIER, &
                dphi_v_out=third_v_local, valid_v_out=valid_third_v_local)
              if (allocated(third_v)) deallocate(third_v)
              if (allocated(valid_third_v)) deallocate(valid_third_v)
              allocate(third_v(size(third_v_local,1), size(third_v_local,2)))
              allocate(valid_third_v(size(valid_third_v_local)))
              third_v = third_v_local
              valid_third_v = valid_third_v_local
              ! RESTORED 2026-09-22 per user report: this original, hand-derived correction
              ! (accounts for each neighbor's own cell-average-vs-point-value gap via its second
              ! moment, AND the cell-blend curvature bias from averaging several vertex samples of
              ! a curved gradient field) is the version that actually reached order 4 -- corrects
              ! ONLY grad_flat in place; hess/third are left at their raw (uncorrected) values.
              ! This session tried two generic eq.13-based alternatives (directly correcting hess
              ! via a z^(m)-contraction operator, and correcting grad then rebuilding hess/third
              ! from it) -- both compiled and looked locally plausible, but neither moved the
              ! observed vortex convergence rate past ~3, so they are reverted in favor of this.
              call apply_grad_bias_correction(mesh, 3_ENTIER, 5_ENTIER, boundary_2d, grad_flat, &
                hess_v, third_v, valid_hess_v, valid_third_v, grad_oi_v)
              if (do_exchange) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 15_ENTIER, grad_flat)
            else
              call compute_next_order_derivative(mesh, 3_ENTIER, 45_ENTIER, boundary_2d, &
                hess_flat, third_flat, deriv_order=3_ENTIER)
            end if
          end if
          if (do_exchange) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 135_ENTIER, third_flat)

          do e = 1, mesh%n_elems
            do dir3 = 1, 3
              do dir2 = 1, 3
                do dir1 = 1, 3
                  do v = 1, 5
                    third(v, dir1, dir2, dir3, e) = &
                      third_flat((dir3-1)*45 + (dir2-1)*15 + (dir1-1)*5 + v, e)
                  end do
                end do
              end do
            end do
          end do
        end if
      end if
    end if

    do e = 1, mesh%n_elems
      do dir1 = 1, 3
        do v = 1, 5
          grad(v, dir1, e) = grad_flat((dir1-1)*5 + v, e)
        end do
      end do
    end do

    if (kill_recons) call apply_positivity_kill(mesh, prim, grad, hess, third)

    if (allocated(grad_oi_v)) deallocate(grad_oi_v)
    if (allocated(hess_v)) deallocate(hess_v)
    if (allocated(valid_hess_v)) deallocate(valid_hess_v)
    if (allocated(third_v)) deallocate(third_v)
    if (allocated(valid_third_v)) deallocate(valid_third_v)
    end associate
  end subroutine aho_reconstruction

  ! kill_recons safety net: if the linear (grad-only) reconstruction predicts negative rho or p
  ! at any of a cell's own vertices, zero that cell's grad/hess/third entirely -- forces a
  ! first-order/constant reconstruction for that cell rather than extrapolating a
  ! positivity-violating state. See the wall-node-reconstruction-todo memory: this is the
  ! previously-proposed "expedient" fix for the missing cell-average positivity limiter, now
  ! generalized to any cell (not just wall-adjacent ones).
  subroutine apply_positivity_kill(mesh, prim, grad, hess, third)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: prim
    real(kind=DOUBLE), dimension(5, 3, mesh%n_elems), intent(inout) :: grad
    real(kind=DOUBLE), dimension(:, :, :, :), allocatable, intent(inout) :: hess
    real(kind=DOUBLE), dimension(:, :, :, :, :), allocatable, intent(inout), optional :: third

    integer(kind=ENTIER) :: i, k, id_vert
    real(kind=DOUBLE), dimension(3) :: dx
    real(kind=DOUBLE) :: rho_v, p_v
    logical :: bad

    do i = 1, mesh%n_elems
      if (mesh%elem(i)%is_ghost) cycle
      bad = .false.
      do k = 1, mesh%elem(i)%n_vert
        id_vert = mesh%elem(i)%vert(k)
        dx = mesh%vert(id_vert)%coord - mesh%elem(i)%coord
        rho_v = prim(1, i) + dot_product(grad(1, :, i), dx)
        p_v   = prim(5, i) + dot_product(grad(5, :, i), dx)
        if (rho_v <= 0.0_DOUBLE .or. p_v <= 0.0_DOUBLE) then
          bad = .true.
          exit
        end if
      end do
      if (bad) then
        grad(:, :, i) = 0.0_DOUBLE
        if (allocated(hess)) hess(:, :, :, i) = 0.0_DOUBLE
        if (present(third)) then
          if (allocated(third)) third(:, :, :, :, i) = 0.0_DOUBLE
        end if
      end if
    end do
  end subroutine apply_positivity_kill

  ! Solves ATA_in*x=rhs_in (a Gram/normal matrix), dropping directions whose entire
  ! row/column is zero (no information in the data) rather than giving up on the
  ! whole vector at the first degenerate pivot.
  subroutine solve_ls(ATA_in, rhs_in, x, n)
    implicit none

    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), dimension(9, 9), intent(in) :: ATA_in
    real(kind=DOUBLE), dimension(9),    intent(in) :: rhs_in
    real(kind=DOUBLE), dimension(9),    intent(out) :: x

    integer(kind=ENTIER), dimension(9) :: idx
    integer(kind=ENTIER) :: n_keep, ii, jj, fail_at
    real(kind=DOUBLE), dimension(9, 9) :: ATA_red
    real(kind=DOUBLE), dimension(9) :: rhs_red, x_red
    real(kind=DOUBLE) :: col_tol

    x = 0.0_DOUBLE

    col_tol = 1.0e-12_DOUBLE * maxval(abs(ATA_in(1:n, 1:n)))

    n_keep = 0
    do jj = 1, n
      if (any(abs(ATA_in(1:n, jj)) > col_tol)) then
        n_keep = n_keep + 1
        idx(n_keep) = jj
      end if
    end do

    do
      if (n_keep == 0) return

      do ii = 1, n_keep
        do jj = 1, n_keep
          ATA_red(ii, jj) = ATA_in(idx(ii), idx(jj))
        end do
        rhs_red(ii) = rhs_in(idx(ii))
      end do

      call gauss_solve(ATA_red, rhs_red, x_red, n_keep, fail_at)

      if (fail_at == 0) then
        do ii = 1, n_keep
          x(idx(ii)) = x_red(ii)
        end do
        return
      end if

      do ii = fail_at, n_keep - 1
        idx(ii) = idx(ii + 1)
      end do
      n_keep = n_keep - 1
    end do
  end subroutine solve_ls

  ! Gaussian elimination with partial pivoting; fail_at=jj on a degenerate pivot at
  ! step jj (x left at 0), fail_at=0 on a complete solve. See solve_ls's own retry loop.
  subroutine gauss_solve(ATA_in, rhs_in, x, n, fail_at)
    implicit none

    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), dimension(9, 9), intent(in) :: ATA_in
    real(kind=DOUBLE), dimension(9),    intent(in) :: rhs_in
    real(kind=DOUBLE), dimension(9),    intent(out) :: x
    integer(kind=ENTIER), intent(out) :: fail_at

    real(kind=DOUBLE), dimension(9, 10) :: aug
    integer(kind=ENTIER) :: ii, jj, kk, piv_row
    real(kind=DOUBLE) :: pivot_val, fac, tmp, piv_tol

    x = 0.0_DOUBLE
    fail_at = 0
    aug(:, 1:n)  = ATA_in(1:9, 1:n)
    aug(:, n+1)  = rhs_in

    piv_tol = 1.0e-12_DOUBLE * maxval(abs(ATA_in(1:n, 1:n)))

    do jj = 1, n
      piv_row   = jj
      pivot_val = abs(aug(jj, jj))
      do ii = jj + 1, n
        if (abs(aug(ii, jj)) > pivot_val) then
          pivot_val = abs(aug(ii, jj))
          piv_row   = ii
        end if
      end do

      if (piv_row /= jj) then
        do kk = jj, n + 1
          tmp            = aug(jj,      kk)
          aug(jj,      kk) = aug(piv_row, kk)
          aug(piv_row, kk) = tmp
        end do
      end if

      if (abs(aug(jj, jj)) < piv_tol) then
        fail_at = jj
        return
      end if

      do ii = jj + 1, n
        fac = aug(ii, jj) / aug(jj, jj)
        aug(ii, jj:n+1) = aug(ii, jj:n+1) - fac * aug(jj, jj:n+1)
      end do
    end do

    do ii = n, 1, -1
      x(ii) = aug(ii, n + 1)
      do kk = ii + 1, n
        x(ii) = x(ii) - aug(ii, kk) * x(kk)
      end do
      x(ii) = x(ii) / aug(ii, ii)
    end do
  end subroutine gauss_solve

  ! Forces the next ensure_face_quad_cache call to regenerate the quadrature
  ! points from the current face geometry.
  !
  ! The cache is keyed on (n_faces, order), which assumes what the comment
  ! below states: that face geometry never changes during a run. A solver that
  ! MOVES nodes breaks that assumption without changing either key, so the
  ! quadrature points stay frozen at their pre-move positions while the fluxes
  ! keep being integrated on them. That silently destroys the solution at
  ! order >= 2 (order 1 survives, since it never evaluates a reconstruction
  ! away from the cell average). Call this after move_mesh/compute_geometry_mesh.
  subroutine invalidate_face_quad_cache()
    implicit none

    face_quad_cache_n_faces = -1
    face_quad_cache_order = -1
  end subroutine invalidate_face_quad_cache

  ! Builds face_quad_pts_cache/face_quad_wts_cache/face_quad_offset_cache once per (mesh, order):
  ! face_flux_loop's order>=3 branch used to regenerate each face's quadrature rule from scratch
  ! on every call (every RK stage, every iteration) even though face geometry never changes during
  ! a run -- this precomputes it once and face_flux_loop just looks it up.
  subroutine ensure_face_quad_cache(mesh)
    implicit none

    type(mesh_type), intent(in) :: mesh

    integer(kind=ENTIER) :: iface, k, iv, n_fvert, n_qpts, total_qpts, pos
    real(kind=DOUBLE), dimension(:, :), allocatable :: face_coords, qpts_local
    real(kind=DOUBLE), dimension(:), allocatable :: qwts_local
    logical :: use_face_quad

    if (face_quad_cache_n_faces == mesh%n_faces .and. face_quad_cache_order == order) return

    if (allocated(face_quad_offset_cache)) deallocate(face_quad_offset_cache)
    allocate(face_quad_offset_cache(mesh%n_faces + 1))
    face_quad_offset_cache(1) = 1

    do iface = 1, mesh%n_faces
      n_fvert = mesh%face(iface)%n_vert
      use_face_quad = order >= 3 .and. n_fvert > 0 .and. allocated(mesh%face(iface)%vert)
      if (use_face_quad) then
        allocate(face_coords(3, n_fvert))
        do k = 1, n_fvert
          iv = mesh%face(iface)%vert(k)
          face_coords(:, k) = mesh%vert(iv)%coord
        end do
        call face_quad_pts(int(n_fvert, ENTIER), face_coords, int(order - 1, ENTIER), &
          qpts_local, qwts_local)
        deallocate(face_coords)
        n_qpts = size(qwts_local)
        deallocate(qpts_local, qwts_local)
      else
        n_qpts = 1
      end if
      face_quad_offset_cache(iface + 1) = face_quad_offset_cache(iface) + n_qpts
    end do

    total_qpts = face_quad_offset_cache(mesh%n_faces + 1) - 1
    if (allocated(face_quad_pts_cache)) deallocate(face_quad_pts_cache)
    if (allocated(face_quad_wts_cache)) deallocate(face_quad_wts_cache)
    allocate(face_quad_pts_cache(3, total_qpts), face_quad_wts_cache(total_qpts))

    do iface = 1, mesh%n_faces
      n_fvert = mesh%face(iface)%n_vert
      pos = face_quad_offset_cache(iface)
      use_face_quad = order >= 3 .and. n_fvert > 0 .and. allocated(mesh%face(iface)%vert)
      if (use_face_quad) then
        allocate(face_coords(3, n_fvert))
        do k = 1, n_fvert
          iv = mesh%face(iface)%vert(k)
          face_coords(:, k) = mesh%vert(iv)%coord
        end do
        call face_quad_pts(int(n_fvert, ENTIER), face_coords, int(order - 1, ENTIER), &
          qpts_local, qwts_local)
        deallocate(face_coords)
        n_qpts = size(qwts_local)
        face_quad_pts_cache(:, pos:pos+n_qpts-1) = qpts_local
        face_quad_wts_cache(pos:pos+n_qpts-1) = qwts_local
        deallocate(qpts_local, qwts_local)
      else
        face_quad_pts_cache(:, pos) = mesh%face(iface)%coord
        face_quad_wts_cache(pos) = mesh%face(iface)%area
      end if
    end do

    face_quad_cache_n_faces = mesh%n_faces
    face_quad_cache_order = order
  end subroutine ensure_face_quad_cache

  subroutine face_flux_loop(mesh, sol, prim, grad, hess, rhs, sum_lambda, t, third)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(in)  :: sol
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in)  :: prim
    real(kind=DOUBLE), dimension(5, 3, mesh%n_elems), intent(in) :: grad
    real(kind=DOUBLE), dimension(:, :, :, :), allocatable, intent(in) :: hess
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems),    intent(inout) :: sum_lambda
    real(kind=DOUBLE), intent(in) :: t
    real(kind=DOUBLE), dimension(:, :, :, :, :), allocatable, intent(in), optional :: third

    integer(kind=ENTIER) :: iface, il, ir, n_fvert, n_qpts, q, qpos
    real(kind=DOUBLE), dimension(3) :: norm, xface
    real(kind=DOUBLE) :: qwt
    real(kind=DOUBLE), dimension(5) :: wL, wR, flux
    real(kind=DOUBLE) :: lambda, gL, gR, face_lambda_accum
    real(kind=DOUBLE) :: flux_rgm1, Gl_rgm1, Gr_rgm1
    logical :: is_zface

    call ensure_face_quad_cache(mesh)

    do iface = 1, mesh%n_faces
      il   = mesh%face(iface)%left_neigh
      ir   = mesh%face(iface)%right_neigh
      norm = mesh%face(iface)%norm

      if (mesh%elem(il)%is_ghost) cycle

      ! Ghost/boundary cell has no gamma_arr of its own -- assume same material as il.
      gL = gamma_arr(il)
      if (ir > 0) then
        gR = gamma_arr(ir)
      else
        gR = gL
      end if

      n_fvert = mesh%face(iface)%n_vert

      ! z-faces on a boundary_2d extrusion carry no physics: skip before touching qpts.
      is_zface = boundary_2d .and. abs(abs(norm(3)) - 1.0_DOUBLE) < 1.0e-6_DOUBLE
      if (is_zface) cycle

      n_qpts = face_quad_offset_cache(iface + 1) - face_quad_offset_cache(iface)

      face_lambda_accum = 0.0_DOUBLE
      do q = 1, n_qpts
        qpos  = face_quad_offset_cache(iface) + q - 1
        xface = face_quad_pts_cache(:, qpos)
        qwt   = face_quad_wts_cache(qpos)

        wL = reconstruct(prim, grad, hess, il, xface, mesh%elem(il)%coord, gL, third)
        if (ir > 0) then
          wR = reconstruct(prim, grad, hess, ir, xface, mesh%elem(ir)%coord, gR, third)
        else
          wR = ghost_prim(xface, norm, wL, -ir, t)
        end if

        select case (flux_scheme_id)
        case (FLUX_THREE_WAVE)
          flux = three_wave_flux(wL, wR, norm, gL, gR) * qwt
        case (FLUX_MODIFIED_THREE_WAVE)
          flux = modified_three_wave_flux(wL, wR, norm, gL, gR) * qwt
        case (FLUX_TWO_WAVE)
          flux = two_wave_flux(wL, wR, norm, gL, gR) * qwt
        case default
          flux = rusanov(wL, wR, norm, gL, gR) * qwt
        end select
        lambda = max(abs(dot_product(wL(2:4), norm)) + cs(wL, gL), &
                     abs(dot_product(wR(2:4), norm)) + cs(wR, gR)) * qwt
        face_lambda_accum = face_lambda_accum + lambda

        ! Gamma-transport: upwind by the sign of the mass flux just computed.
        Gl_rgm1 = sol(6, il) / sol(1, il)
        if (ir > 0) then
          Gr_rgm1 = sol(6, ir) / sol(1, ir)
        else
          Gr_rgm1 = Gl_rgm1
        end if
        if (flux(1) >= 0.0_DOUBLE) then
          flux_rgm1 = flux(1) * Gl_rgm1
        else
          flux_rgm1 = flux(1) * Gr_rgm1
        end if

        rhs(1:5, il) = rhs(1:5, il) - flux
        rhs(6, il)   = rhs(6, il)   - flux_rgm1
        if (ir > 0) then
          if (.not. mesh%elem(ir)%is_ghost) then
            rhs(1:5, ir) = rhs(1:5, ir) + flux
            rhs(6, ir)   = rhs(6, ir)   + flux_rgm1
          end if
        end if
      end do

      if (use_max_lambda_dt) then
        sum_lambda(il) = max(sum_lambda(il), face_lambda_accum)
      else
        sum_lambda(il) = sum_lambda(il) + face_lambda_accum
      end if
      if (ir > 0) then
        if (.not. mesh%elem(ir)%is_ghost) then
          if (use_max_lambda_dt) then
            sum_lambda(ir) = max(sum_lambda(ir), face_lambda_accum)
          else
            sum_lambda(ir) = sum_lambda(ir) + face_lambda_accum
          end if
        end if
      end if
    end do
  end subroutine face_flux_loop

  ! Rejects a reconstructed state with negative rho/p, or velocity far above the
  ! reference speed+sound-speed scale.
  pure function physical_state(w_cand, w_ref, g) result(ok)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: w_cand, w_ref
    real(kind=DOUBLE), intent(in) :: g
    logical :: ok

    real(kind=DOUBLE) :: spd_ref, spd_cand, c_ref

    ok = .false.
    if (w_cand(1) <= 0.0_DOUBLE) return
    if (w_cand(5) <= 0.0_DOUBLE) return

    spd_ref  = w_ref(2)**2  + w_ref(3)**2  + w_ref(4)**2
    spd_cand = w_cand(2)**2 + w_cand(3)**2 + w_cand(4)**2
    c_ref    = g * w_ref(5) / max(w_ref(1), 1.0e-16_DOUBLE)
    if (spd_cand > 400.0_DOUBLE * (spd_ref + c_ref)) return

    ok = .true.
  end function physical_state

  ! Second (order>=4: also third) geometric moment of each cell about its own centroid,
  ! so reconstruct()'s Taylor terms have the correct cell AVERAGE, not just centroid value.
  subroutine compute_cell_moments(mesh)
    implicit none

    type(mesh_type), intent(in) :: mesh

    integer(kind=ENTIER) :: i, k, n_v, quad_order
    real(kind=DOUBLE), dimension(:, :), allocatable :: vcoords, pts
    real(kind=DOUBLE), dimension(:), allocatable :: wts
    real(kind=DOUBLE), dimension(3) :: xc
    integer(kind=ENTIER) :: jj, kk, ll
    logical :: need_third

    need_third = (order >= 4)
    quad_order = merge(3_ENTIER, 2_ENTIER, need_third)

    if (allocated(cell_moment)) deallocate(cell_moment)
    allocate(cell_moment(3, 3, mesh%n_elems))
    cell_moment = 0.0_DOUBLE

    if (allocated(cell_moment3)) deallocate(cell_moment3)
    if (need_third) then
      allocate(cell_moment3(3, 3, 3, mesh%n_elems))
      cell_moment3 = 0.0_DOUBLE
    end if

    do i = 1, mesh%n_elems
      n_v = mesh%elem(i)%n_vert
      if (n_v /= 4 .and. n_v /= 5 .and. n_v /= 6 .and. n_v /= 8) cycle
      allocate(vcoords(3, n_v))
      do k = 1, n_v
        vcoords(:, k) = mesh%vert(mesh%elem(i)%vert(k))%coord
      end do
      call volume_quad_pts(n_v, vcoords, quad_order, pts, wts)
      deallocate(vcoords)

      xc = mesh%elem(i)%coord
      do jj = 1, 3
        do kk = 1, 3
          cell_moment(jj, kk, i) = sum(wts * (pts(jj, :) - xc(jj)) * (pts(kk, :) - xc(kk))) &
            / max(mesh%elem(i)%volume, 1.0e-300_DOUBLE)
        end do
      end do
      if (need_third) then
        do jj = 1, 3
          do kk = 1, 3
            do ll = 1, 3
              cell_moment3(jj, kk, ll, i) = sum(wts * (pts(jj, :) - xc(jj)) &
                * (pts(kk, :) - xc(kk)) * (pts(ll, :) - xc(ll))) &
                / max(mesh%elem(i)%volume, 1.0e-300_DOUBLE)
            end do
          end do
        end do
      end if
      deallocate(pts, wts)
    end do
  end subroutine compute_cell_moments

  ! Polynomial reconstruction of w at point xq from cell i; hierarchical fallback to
  ! order-(p-1) if the order-p candidate is unphysical.
  function reconstruct(prim, grad, hess, i, xq, xc, g, third) result(w)
    implicit none

    real(kind=DOUBLE), dimension(:, :),          intent(in) :: prim
    real(kind=DOUBLE), dimension(:, :, :),       intent(in) :: grad
    real(kind=DOUBLE), dimension(:, :, :, :), allocatable, intent(in) :: hess
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE), dimension(3), intent(in) :: xq, xc
    real(kind=DOUBLE), intent(in) :: g
    real(kind=DOUBLE), dimension(:, :, :, :, :), allocatable, intent(in), optional :: third
    real(kind=DOUBLE), dimension(5) :: w

    real(kind=DOUBLE), dimension(3) :: dx
    real(kind=DOUBLE), dimension(5) :: w_try, w_prev, grad_term
    integer(kind=ENTIER) :: j, k, l

    dx = xq - xc
    w  = prim(:, i)
    grad_term = matmul(grad(:, :, i), dx)

    if (order >= 2) then
      w_try = w + grad_term
      if (physical_state(w_try, w, g)) w = w_try
    end if

    if (order >= 3 .and. allocated(hess)) then
      if (.not. allocated(cell_moment)) then
        print *, 'FATAL: reconstruct() called at order>=3 without compute_cell_moments.'
        error stop 1
      end if
      ! Schwarz symmetry: dx(j)*dx(k) and cell_moment(j,k,i) both commute exactly in j,k (the
      ! latter by construction -- see compute_cell_moments), so the off-diagonal (j,k) and (k,j)
      ! terms of the original 3x3 sum share the same (dx(j)*dx(k) - cell_moment(j,k,i)) factor.
      ! Summing hess(:,j,k,i)+hess(:,k,j,i) once and reusing that factor is EXACT algebra (not an
      ! approximation assuming hess itself is perfectly symmetric) -- 6 vector terms instead of 9.
      w_try = prim(:, i) + grad_term
      do j = 1, 3
        w_try = w_try + 0.5_DOUBLE * hess(:, j, j, i) * (dx(j) * dx(j) - cell_moment(j, j, i))
      end do
      do j = 1, 2
        do k = j + 1, 3
          w_try = w_try + 0.5_DOUBLE * (hess(:, j, k, i) + hess(:, k, j, i)) &
            * (dx(j) * dx(k) - cell_moment(j, k, i))
        end do
      end do
      if (physical_state(w_try, prim(:, i), g)) w = w_try
    end if

    if (order >= 4 .and. present(third)) then
      if (allocated(third)) then
        if (.not. allocated(cell_moment3)) then
          print *, 'FATAL: reconstruct() called at order>=4 without cell_moment3 filled.'
          error stop 1
        end if
        ! Same exact regrouping for the fully symmetric third-order term: dx(j)*dx(k)*dx(l) and
        ! cell_moment3(j,k,l,i) commute in all 3 indices (again by construction), so the 27-term
        ! sum collapses to 10 sorted-index groups, each summing third(:,.,.,.,i) over its distinct
        ! permutations (1 for jjj, 3 for jjk, 6 for jkl) against one shared dx-moment factor.
        w_prev = w
        w_try = w_prev
        do j = 1, 3
          w_try = w_try + (1.0_DOUBLE / 6.0_DOUBLE) * third(:, j, j, j, i) &
            * (dx(j) * dx(j) * dx(j) - cell_moment3(j, j, j, i))
        end do
        do j = 1, 3
          do k = 1, 3
            if (j == k) cycle
            ! (j,j,k) pattern: 3 permutations jjk, jkj, kjj.
            w_try = w_try + (1.0_DOUBLE / 6.0_DOUBLE) &
              * (third(:, j, j, k, i) + third(:, j, k, j, i) + third(:, k, j, j, i)) &
              * (dx(j) * dx(j) * dx(k) - cell_moment3(j, j, k, i))
          end do
        end do
        ! (1,2,3) pattern: all 6 permutations.
        w_try = w_try + (1.0_DOUBLE / 6.0_DOUBLE) &
          * (third(:, 1, 2, 3, i) + third(:, 1, 3, 2, i) + third(:, 2, 1, 3, i) &
             + third(:, 2, 3, 1, i) + third(:, 3, 1, 2, i) + third(:, 3, 2, 1, i)) &
          * (dx(1) * dx(2) * dx(3) - cell_moment3(1, 2, 3, i))
        if (physical_state(w_try, w_prev, g)) w = w_try
      end if
    end if

    w(1) = max(w(1), 1.0e-12_DOUBLE)
    w(5) = max(w(5), 1.0e-12_DOUBLE)
  end function reconstruct

  function ghost_prim(xf, norm, wL, id_bc, t) result(wR)
    implicit none

    real(kind=DOUBLE), dimension(3), intent(in) :: xf, norm
    real(kind=DOUBLE), dimension(5), intent(in) :: wL
    integer(kind=ENTIER), intent(in) :: id_bc
    real(kind=DOUBLE), intent(in) :: t
    real(kind=DOUBLE), dimension(5) :: wR

    real(kind=DOUBLE) :: vn

    if (id_bc < 1 .or. id_bc > n_bc) then
      wR    = wL
      vn    = dot_product(wL(2:4), norm)
      wR(2) = wL(2) - 2.0_DOUBLE * vn * norm(1)
      wR(3) = wL(3) - 2.0_DOUBLE * vn * norm(2)
      wR(4) = wL(4) - 2.0_DOUBLE * vn * norm(3)
      return
    end if

    select case (bc_type_id(id_bc))
    case (BC_FREESTREAM)
      wR = bc_val(:, id_bc)
    case (BC_OUTFLOW)
      wR = wL
    case (BC_DMR_TOP)
      wR = dmr_state(xf(1), xf(2), t)
    case default
      wR    = wL
      vn    = dot_product(wL(2:4), norm)
      wR(2) = wL(2) - 2.0_DOUBLE * vn * norm(1)
      wR(3) = wL(3) - 2.0_DOUBLE * vn * norm(2)
      wR(4) = wL(4) - 2.0_DOUBLE * vn * norm(3)
    end select
  end function ghost_prim

  ! ================================================================
  ! Euler physics
  ! ================================================================

  pure function conserv_to_primit(u, g) result(w)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: u
    real(kind=DOUBLE), intent(in) :: g
    real(kind=DOUBLE), dimension(5) :: w

    real(kind=DOUBLE) :: rho_inv, ke

    rho_inv = 1.0_DOUBLE / max(u(1), 1.0e-16_DOUBLE)
    w(1) = u(1)
    w(2) = u(2) * rho_inv
    w(3) = u(3) * rho_inv
    w(4) = u(4) * rho_inv
    ke   = 0.5_DOUBLE * (w(2)**2 + w(3)**2 + w(4)**2)
    w(5) = (g - 1.0_DOUBLE) * (u(5) - u(1) * ke)
  end function conserv_to_primit

  pure function primit_to_conserv(w, g) result(u)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE), intent(in) :: g
    real(kind=DOUBLE), dimension(5) :: u

    real(kind=DOUBLE) :: ke

    ke   = 0.5_DOUBLE * (w(2)**2 + w(3)**2 + w(4)**2)
    u(1) = w(1)
    u(2) = w(1) * w(2)
    u(3) = w(1) * w(3)
    u(4) = w(1) * w(4)
    u(5) = w(1) * (ke + w(5) / ((g - 1.0_DOUBLE) * w(1)))
  end function primit_to_conserv

  pure function cs(w, g) result(c)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE), intent(in) :: g
    real(kind=DOUBLE) :: c

    c = sqrt(max(g * w(5) / max(w(1), 1.0e-16_DOUBLE), 0.0_DOUBLE))
  end function cs

  pure function euler_flux(w, n, g) result(F)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), intent(in) :: g
    real(kind=DOUBLE), dimension(5) :: F

    real(kind=DOUBLE) :: vn, rho, p, e

    rho = w(1); p = w(5)
    vn  = w(2)*n(1) + w(3)*n(2) + w(4)*n(3)
    e   = p / ((g - 1.0_DOUBLE) * rho) + 0.5_DOUBLE*(w(2)**2+w(3)**2+w(4)**2)

    F(1) = rho * vn
    F(2) = rho * w(2) * vn + p * n(1)
    F(3) = rho * w(3) * vn + p * n(2)
    F(4) = rho * w(4) * vn + p * n(3)
    F(5) = rho * e * vn + p * vn
  end function euler_flux

  pure function rusanov(wL, wR, n, gL, gR) result(F)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: wL, wR
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), intent(in) :: gL, gR
    real(kind=DOUBLE), dimension(5) :: F

    real(kind=DOUBLE) :: lam, vnL, vnR

    vnL = wL(2)*n(1) + wL(3)*n(2) + wL(4)*n(3)
    vnR = wR(2)*n(1) + wR(3)*n(2) + wR(4)*n(3)
    lam = max(abs(vnL) + cs(wL, gL), abs(vnR) + cs(wR, gR))

    F = 0.5_DOUBLE * (euler_flux(wL, n, gL) + euler_flux(wR, n, gR)) &
      - 0.5_DOUBLE * lam * (primit_to_conserv(wR, gR) - primit_to_conserv(wL, gL))
  end function rusanov

  ! HLLC-family 3-wave approximate Riemann solver (Toro, ch. 10); own local copy of
  ! ns_euler_rs_module's three_wave so gamma_arr multi-material support stays self-contained.
  pure function three_wave_flux(wL, wR, n, gL, gR) result(F)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: wL, wR
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), intent(in) :: gL, gR
    real(kind=DOUBLE), dimension(5) :: F

    real(kind=DOUBLE), dimension(5) :: uL, uR, fL, fR, uL_star, uR_star
    real(kind=DOUBLE) :: rhoL, rhoR, vnL, vnR, pL, pR, eL, eR, aL, aR
    real(kind=DOUBLE) :: lambdaL, lambdaR, v_star
    real(kind=DOUBLE) :: rhoL_star, rhoR_star, pL_star, pR_star, sL_wave, sR_wave

    rhoL = wL(1); vnL = dot_product(wL(2:4), n); pL = wL(5)
    uL = primit_to_conserv(wL, gL); eL = uL(5)/rhoL; aL = cs(wL, gL)

    rhoR = wR(1); vnR = dot_product(wR(2:4), n); pR = wR(5)
    uR = primit_to_conserv(wR, gR); eR = uR(5)/rhoR; aR = cs(wR, gR)

    fL(1)   = vnL*uL(1)
    fL(2:4) = vnL*uL(2:4) + pL*n
    fL(5)   = (uL(5) + pL)*vnL

    fR(1)   = vnR*uR(1)
    fR(2:4) = vnR*uR(2:4) + pR*n
    fR(5)   = (uR(5) + pR)*vnR

    lambdaL = max(aL*rhoL, sqrt(rhoL*max(0.0_DOUBLE, pR - pL)), -rhoL*(vnR - vnL))
    lambdaR = max(aR*rhoR, sqrt(rhoR*max(0.0_DOUBLE, pL - pR)), -rhoR*(vnR - vnL))
    v_star  = (lambdaL*vnL + lambdaR*vnR - (pR - pL)) / (lambdaR + lambdaL)

    rhoL_star = 1.0_DOUBLE / (1.0_DOUBLE/rhoL + (v_star - vnL)/lambdaL)
    pL_star   = pL - lambdaL*(v_star - vnL)
    uL_star(1)   = rhoL_star
    uL_star(2:4) = rhoL_star*(wL(2:4) + (v_star - vnL)*n)
    uL_star(5)   = rhoL_star*(eL + (pL*vnL - pL_star*v_star)/lambdaL)

    rhoR_star = 1.0_DOUBLE / (1.0_DOUBLE/rhoR + (vnR - v_star)/lambdaR)
    pR_star   = pR + lambdaR*(v_star - vnR)
    uR_star(1)   = rhoR_star
    uR_star(2:4) = rhoR_star*(wR(2:4) + (v_star - vnR)*n)
    uR_star(5)   = rhoR_star*(eR + (pR_star*v_star - pR*vnR)/lambdaR)

    sL_wave = vnL - lambdaL/rhoL
    sR_wave = vnR + lambdaR/rhoR

    F = 0.5_DOUBLE*(fL + fR) - 0.5_DOUBLE * ( &
      abs(sL_wave)*(uL_star - uL) + &
      abs(v_star)*(uR_star - uL_star) + &
      abs(sR_wave)*(uR - uR_star))
  end function three_wave_flux

  ! Modified 3-wave solver (Toro-family HLLC variant): blends a SINGLE shared tangential
  ! velocity between both star states, instead of three_wave_flux's plain carry-through of
  ! each side's own unchanged tangential velocity, with a matching tang_energy correction so
  ! the energy equation stays consistent -- own local copy of ns_euler_rs_module's
  ! modified_three_wave, gamma_arr multi-material aware like this module's own three_wave_flux.
  pure function modified_three_wave_flux(wL, wR, n, gL, gR) result(F)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: wL, wR
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), intent(in) :: gL, gR
    real(kind=DOUBLE), dimension(5) :: F

    real(kind=DOUBLE), dimension(5) :: uL, uR, fL, fR, uL_star, uR_star
    real(kind=DOUBLE) :: rhoL, rhoR, vnL, vnR, pL, pR, eL, eR, aL, aR
    real(kind=DOUBLE) :: lambdaL, lambdaR, v_star, tang_coeff, tang_energy
    real(kind=DOUBLE) :: rhoL_star, rhoR_star, pL_star, pR_star, sL_wave, sR_wave
    real(kind=DOUBLE), dimension(3) :: vtL, vtR, vt_star, dvt

    rhoL = wL(1); vnL = dot_product(wL(2:4), n); pL = wL(5)
    uL = primit_to_conserv(wL, gL); eL = uL(5)/rhoL; aL = cs(wL, gL)

    rhoR = wR(1); vnR = dot_product(wR(2:4), n); pR = wR(5)
    uR = primit_to_conserv(wR, gR); eR = uR(5)/rhoR; aR = cs(wR, gR)

    fL(1)   = vnL*uL(1)
    fL(2:4) = vnL*uL(2:4) + pL*n
    fL(5)   = (uL(5) + pL)*vnL

    fR(1)   = vnR*uR(1)
    fR(2:4) = vnR*uR(2:4) + pR*n
    fR(5)   = (uR(5) + pR)*vnR

    lambdaL = max(aL*rhoL, sqrt(rhoL*max(0.0_DOUBLE, pR - pL)), -rhoL*(vnR - vnL))
    lambdaR = max(aR*rhoR, sqrt(rhoR*max(0.0_DOUBLE, pL - pR)), -rhoR*(vnR - vnL))
    v_star  = (lambdaL*vnL + lambdaR*vnR - (pR - pL)) / (lambdaR + lambdaL)

    vtL = wL(2:4) - vnL*n
    vtR = wR(2:4) - vnR*n
    vt_star = (lambdaL*vtL + lambdaR*vtR) / (lambdaL + lambdaR)
    dvt = vtR - vtL

    tang_coeff  = -lambdaL*lambdaR / (lambdaL + lambdaR)
    tang_energy = tang_coeff * dot_product(dvt, vt_star)

    rhoL_star = 1.0_DOUBLE / (1.0_DOUBLE/rhoL + (v_star - vnL)/lambdaL)
    pL_star   = pL - lambdaL*(v_star - vnL)
    uL_star(1)   = rhoL_star
    uL_star(2:4) = rhoL_star*(v_star*n + vt_star)
    uL_star(5)   = rhoL_star*(eL + (pL*vnL - pL_star*v_star - tang_energy)/lambdaL)

    rhoR_star = 1.0_DOUBLE / (1.0_DOUBLE/rhoR + (vnR - v_star)/lambdaR)
    pR_star   = pR + lambdaR*(v_star - vnR)
    uR_star(1)   = rhoR_star
    uR_star(2:4) = rhoR_star*(v_star*n + vt_star)
    uR_star(5)   = rhoR_star*(eR + (pR_star*v_star - pR*vnR + tang_energy)/lambdaR)

    sL_wave = vnL - lambdaL/rhoL
    sR_wave = vnR + lambdaR/rhoR

    F = 0.5_DOUBLE*(fL + fR) - 0.5_DOUBLE * ( &
      abs(sL_wave)*(uL_star - uL) + &
      abs(v_star)*(uR_star - uL_star) + &
      abs(sR_wave)*(uR - uR_star))
  end function modified_three_wave_flux

  ! HLL-family 2-wave approximate Riemann solver: single star state, no contact wave.
  pure function two_wave_flux(wL, wR, n, gL, gR) result(F)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: wL, wR
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), intent(in) :: gL, gR
    real(kind=DOUBLE), dimension(5) :: F

    real(kind=DOUBLE), dimension(5) :: uL, uR, fL, fR, u_star
    real(kind=DOUBLE) :: rhoL, rhoR, vnL, vnR, pL, pR, aL, aR
    real(kind=DOUBLE) :: lambdaL, lambdaR, sL_wave, sR_wave

    rhoL = wL(1); vnL = dot_product(wL(2:4), n); pL = wL(5)
    aL = cs(wL, gL); uL = primit_to_conserv(wL, gL)

    rhoR = wR(1); vnR = dot_product(wR(2:4), n); pR = wR(5)
    aR = cs(wR, gR); uR = primit_to_conserv(wR, gR)

    fL(1)   = vnL*uL(1)
    fL(2:4) = vnL*uL(2:4) + pL*n
    fL(5)   = (uL(5) + pL)*vnL

    fR(1)   = vnR*uR(1)
    fR(2:4) = vnR*uR(2:4) + pR*n
    fR(5)   = (uR(5) + pR)*vnR

    lambdaL = max(aL*rhoL, sqrt(rhoL*max(0.0_DOUBLE, pR - pL)), -rhoL*(vnR - vnL))
    lambdaR = max(aR*rhoR, sqrt(rhoR*max(0.0_DOUBLE, pL - pR)), -rhoR*(vnR - vnL))

    sL_wave = vnL - lambdaL/rhoL
    sR_wave = vnR + lambdaR/rhoR

    u_star = (sR_wave*uR - sL_wave*uL - (fR - fL)) / (sR_wave - sL_wave)

    F = 0.5_DOUBLE*(fL + fR) - 0.5_DOUBLE*( &
      abs(sL_wave)*(u_star - uL) + abs(sR_wave)*(uR - u_star))
  end function two_wave_flux

  ! Isentropic vortex primitive state at (x, t), Shu vortex with beta=5.
  pure subroutine vortex_prim(x, t, w)
    implicit none

    real(kind=DOUBLE), dimension(3), intent(in)  :: x
    real(kind=DOUBLE), intent(in) :: t
    real(kind=DOUBLE), dimension(5), intent(out) :: w

    real(kind=DOUBLE) :: beta, r2, tmp, dx, dy

    beta = 5.0_DOUBLE
    dx   = x(1) - u_bg_vortex * t
    dy   = x(2) - v_bg_vortex * t
    r2   = dx**2 + dy**2

    tmp  = (gamma_gas - 1.0_DOUBLE) * beta**2 / (8.0_DOUBLE * gamma_gas * PI**2) * exp(1.0_DOUBLE - r2)
    w(1) = (1.0_DOUBLE - tmp)**(1.0_DOUBLE / (gamma_gas - 1.0_DOUBLE))
    w(2) = u_bg_vortex - dy * beta / (2.0_DOUBLE * PI) * exp(0.5_DOUBLE * (1.0_DOUBLE - r2))
    w(3) = v_bg_vortex + dx * beta / (2.0_DOUBLE * PI) * exp(0.5_DOUBLE * (1.0_DOUBLE - r2))
    w(4) = 0.0_DOUBLE
    w(5) = w(1)**gamma_gas
  end subroutine vortex_prim

  ! High-order cell average of the exact vortex field via volume quadrature (works for
  ! any cell type; falls back to a centroid point value for unsupported cell kinds).
  subroutine vortex_cell_average(mesh, i, t, w)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE), intent(in) :: t
    real(kind=DOUBLE), dimension(5), intent(out) :: w

    integer(kind=ENTIER) :: n_v, k, q
    real(kind=DOUBLE), dimension(:, :), allocatable :: vcoords, pts
    real(kind=DOUBLE), dimension(:), allocatable :: wts
    real(kind=DOUBLE), dimension(5) :: w_sum, wq
    real(kind=DOUBLE) :: vol_q

    n_v = mesh%elem(i)%n_vert
    if (n_v /= 4 .and. n_v /= 5 .and. n_v /= 6 .and. n_v /= 8) then
      call vortex_prim(mesh%elem(i)%coord, t, w)
      return
    end if

    allocate(vcoords(3, n_v))
    do k = 1, n_v
      vcoords(:, k) = mesh%vert(mesh%elem(i)%vert(k))%coord
    end do
    call volume_quad_pts(n_v, vcoords, 3_ENTIER, pts, wts)
    deallocate(vcoords)

    w_sum = 0.0_DOUBLE; vol_q = 0.0_DOUBLE
    do q = 1, size(wts)
      call vortex_prim(pts(:, q), t, wq)
      w_sum = w_sum + wts(q) * wq
      vol_q = vol_q + wts(q)
    end do
    if (vol_q > 0.0_DOUBLE) then
      w = w_sum / vol_q
    else
      call vortex_prim(mesh%elem(i)%coord, t, w)
    end if
    deallocate(pts, wts)
  end subroutine vortex_cell_average

  ! Double Mach reflection (Woodward & Colella 1984): Mach-10 shock inclined 60deg,
  ! touching the x-axis at x=1/6 at t=0. Used for init=3 (t=0) and the 'dmr_top' BC.
  pure function dmr_state(x, y, t) result(w)
    implicit none

    real(kind=DOUBLE), intent(in) :: x, y, t
    real(kind=DOUBLE), dimension(5) :: w

    real(kind=DOUBLE), parameter :: shock_angle = PI/3.0_DOUBLE
    real(kind=DOUBLE), parameter :: mach_shock = 10.0_DOUBLE
    real(kind=DOUBLE), parameter :: x0 = 1.0_DOUBLE/6.0_DOUBLE
    real(kind=DOUBLE), parameter :: a1 = 1.0_DOUBLE
    real(kind=DOUBLE) :: shock_speed_x, x_shock

    shock_speed_x = mach_shock*a1/sin(shock_angle)
    x_shock = x0 + y/tan(shock_angle) + shock_speed_x*t

    if (x < x_shock) then
      w = (/ 8.0_DOUBLE, 8.25_DOUBLE*cos(PI/6.0_DOUBLE), -8.25_DOUBLE*sin(PI/6.0_DOUBLE), &
             0.0_DOUBLE, 116.5_DOUBLE /)
    else
      w = (/ 1.4_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE /)
    end if
  end function dmr_state

  ! Shu-Osher shock/entropy-wave interaction (Shu & Osher, JCP 1989): Mach-3 shock at
  ! x=-4 moving into a sinusoidal density field, domain [-5,5], run to t=1.8.
  pure function shu_osher_state(x) result(w)
    implicit none

    real(kind=DOUBLE), intent(in) :: x
    real(kind=DOUBLE), dimension(5) :: w

    if (x < -4.0_DOUBLE) then
      w = (/ 3.857143_DOUBLE, 2.629369_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 10.33333_DOUBLE /)
    else
      w = (/ 1.0_DOUBLE + 0.2_DOUBLE*sin(5.0_DOUBLE*x), 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE /)
    end if
  end function shu_osher_state

  ! Woodward-Colella interacting blast waves (JCP 1984): domain [0,1], three uniform
  ! pressure zones (1e5 ratio at the extremes), reflecting walls, run to t=0.038.
  pure function woodward_colella_state(x) result(w)
    implicit none

    real(kind=DOUBLE), intent(in) :: x
    real(kind=DOUBLE), dimension(5) :: w

    if (x < 0.1_DOUBLE) then
      w = (/ 1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1000.0_DOUBLE /)
    else if (x < 0.9_DOUBLE) then
      w = (/ 1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.01_DOUBLE /)
    else
      w = (/ 1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 100.0_DOUBLE /)
    end if
  end function woodward_colella_state

end module euler_ho_module
