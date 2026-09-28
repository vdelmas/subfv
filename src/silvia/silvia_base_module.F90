module silvia_base_module
  use precision_module
  use mesh_module
  use quadrature_module
  use arbitrary_high_order_module, only: compute_next_order_derivative
  implicit none
  
  integer, parameter :: mnbc = 10 !Maximum number of boundaries
  real(kind=DOUBLE), parameter :: PI = 4.0_DOUBLE * datan(1.0_DOUBLE)

  !Default values
  real(kind=DOUBLE), parameter :: gamma = 1.4_DOUBLE 
  real(kind=DOUBLE), parameter :: Prandtl = 0.71_DOUBLE !Used to deduce kappa
  real(kind=DOUBLE) :: Cv_p = 720.19471_DOUBLE

  !If use_sutherland == FALSE, the given mu_p is used,
  !otherwise mu_p is found through sutherland's law 
  !with the reference state given by mu0, T0
  real(kind=DOUBLE) :: mu_p = 1.329e-5_DOUBLE
  logical :: use_sutherland = .FALSE.
  real(kind=DOUBLE) :: mu0 = 1.716e-5_DOUBLE
  real(kind=DOUBLE) :: T0 = 273.15_DOUBLE
  real(kind=DOUBLE), parameter :: C0 = 110.4_DOUBLE

  !Timestep
  real(kind=DOUBLE) :: cfl=1e5
  real(kind=DOUBLE) :: tmax=1.0, t, dt

  !Mesh
  character(len=255) :: meshfile_path, meshfile
  logical :: rescale = .false.
  logical :: periodic_mesh = .FALSE.
  real(kind=DOUBLE) :: rescale_factor=1.0_QUAD

  !Cylinder mapping of the original mesh
  logical :: use_cylinder_map = .FALSE.
  integer(kind=ENTIER) :: dim_cylinder_map = 2

  !Boundary conditions
  integer(kind=ENTIER) :: n_bc = 0
  character(len=255), dimension(mnbc) :: bc_name
  ! bc_val for Euler (rho, u, v, w, p)
  character(len=255), dimension(mnbc) :: bc_type
  real(kind=DOUBLE), dimension(5, mnbc) :: bc_val
  ! bc_val_V for Viscous (u, v, w)
  character(len=255), dimension(mnbc) :: bc_type_V
  real(kind=DOUBLE), dimension(3, mnbc) :: bc_val_V
  ! bc_val_T for Heat (T)
  character(len=255), dimension(mnbc) :: bc_type_T
  real(kind=DOUBLE), dimension(mnbc) :: bc_val_T

  ! Precomputed BC type IDs — filled by init_flags after reading input
  integer, parameter :: BC_V_NEUMANN              = 1  ! "neumann" or ""
  integer, parameter :: BC_V_DIRICHLET            = 2  ! "dirichlet"
  integer, parameter :: BC_T_NEUMANN              = 1  ! "neumann" or ""
  integer, parameter :: BC_T_DIRICHLET            = 2  ! "dirichlet"
  integer, parameter :: BC_EULER_WALL              = 1  ! "wall" or ""
  integer, parameter :: BC_EULER_ADHERENCE_WALL    = 2  ! "adherence_wall"
  integer, parameter :: BC_EULER_FREESTREAM        = 3  ! "freestream"
  integer, parameter :: BC_EULER_OUTFLOWSUPERSONIC = 4  ! "outflowsupersonic"
  integer, parameter :: BC_EULER_INFLOW_POND            = 5  ! "inflow_pond"
  integer, parameter :: BC_EULER_INOUT_DOUBLE_MACH     = 6  ! "inout_double_mach"
  integer, parameter :: BC_EULER_DOUBLE_MACH_BOTTOM    = 7  ! "double_mach_bottom"
  integer, parameter :: BC_EULER_POTENTIAL_FLOW_2D    = 8  ! "potential_flow_2d"
  integer, parameter :: BC_EULER_POTENTIAL_FLOW_3D    = 9  ! "potential_flow_3d"
  integer, dimension(:), allocatable :: bc_V_id
  integer, dimension(:), allocatable :: bc_T_id
  integer, dimension(:), allocatable :: bc_euler_id

  integer(kind=ENTIER), dimension(mnbc) :: bc_type_id = BC_EULER_WALL
  ! Scheme integer IDs
  integer, parameter :: SCHEME_MULTI_POINT          = 1  ! "multi_point"
  integer, parameter :: SCHEME_MULTI_POINT_ISO      = 2  ! "multi_point_iso"
  integer, parameter :: SCHEME_THREE_WAVE           = 3  ! "three_wave"
  integer, parameter :: SCHEME_TWO_WAVE             = 4  ! "two_wave"
  integer, parameter :: SCHEME_MODIFIED_THREE_WAVE  = 5  ! "modified_three_wave"
  integer, parameter :: SCHEME_MULTI_POINT_PRESSURE    = 6  ! "multi_point_pressure"
  integer, parameter :: SCHEME_MULTI_POINT_PRESSURE_PH = 8  ! "multi_point_pressure_ph"
  integer, parameter :: SCHEME_WIP                     = 7  ! "WIP"
  integer, parameter :: SCHEME_USI3D                   = 9  ! "USI3D"
  integer, parameter :: SCHEME_MULTI_POINT_VILAR       = 10 ! "multi_point_vilar"
  integer, parameter :: SCHEME_ZB                   = 42 ! "ZB_*_*"
  ! ZB advection sub-scheme IDs
  integer, parameter :: SCHEME_ADV_AR1D     = 1  ! "AR1D"
  integer, parameter :: SCHEME_ADV_AM       = 2  ! "AM"
  integer, parameter :: SCHEME_ADV_AMISO    = 3  ! "AMISO"
  integer, parameter :: SCHEME_ADV_ARMD     = 4  ! "ARMD"
  integer, parameter :: SCHEME_ADV_ARMDU    = 5  ! "ARMDU"
  integer, parameter :: SCHEME_ADV_ARMDM    = 6  ! "ARMDM"
  integer, parameter :: SCHEME_ADV_ARMDMAT  = 7  ! "ARMDMAT"
  integer, parameter :: SCHEME_ADV_ARMDUMAT = 8  ! "ARMDUMAT"
  integer, parameter :: SCHEME_ADV_ARMDMMAT = 9  ! "ARMDMMAT"
  integer, parameter :: SCHEME_ADV_ARMDWIP     = 10  ! "ARMDWIP"
  ! ZB Lagrange sub-scheme IDs
  integer, parameter :: SCHEME_LAG_LS  = 1  ! "LS"
  integer, parameter :: SCHEME_LAG_LSU = 2  ! "LSU"
  integer, parameter :: SCHEME_LAG_LSM = 3  ! "LSM"
  integer, parameter :: SCHEME_LAG_LSWIP = 4  ! "LSWIP"
  integer, parameter :: SCHEME_LAG_LPP = 5  ! "LPP"
  integer, parameter :: SCHEME_LAG_LVPPP = 6  ! "LVPPP"
  integer, parameter :: SCHEME_LAG_LS1D = 7  ! "LS1D"
  integer, parameter :: SCHEME_LAG_LPF  = 8  ! "LPF"
  ! Runtime scheme selection (set by init_bc_flags)
  integer :: scheme_id = 0
  integer :: scheme_adv_id = 0
  integer :: scheme_lag_id = 0

  !Space
  integer(kind=ENTIER) :: order = 1
  logical :: use_aho_reconstruction = .FALSE.
  logical :: second_order = .FALSE.
  integer(kind=ENTIER) :: method = 4

  !Scheme
  character(len=255) :: scheme = ""
  logical :: exclude_bound_vert = .FALSE.
  integer(kind=ENTIER) :: bc_style = 1
  logical :: boundary_2d = .FALSE.
  logical :: activate_diffusion
  logical :: local_time_step = .TRUE.

  !Init
  logical :: init_uniform = .FALSE.
  real(kind=DOUBLE), dimension(5) :: sol_uniform
  logical :: init_isentropic_vortex = .FALSE.
  logical :: init_gresho = .FALSE.
  logical :: init_potential_flow_2d = .FALSE.
  logical :: init_potential_flow_3d = .FALSE.
  logical :: compute_error = .FALSE.
  logical :: error_2d = .FALSE.
  real(kind=DOUBLE) :: error_2d_h = 0.025_DOUBLE
  logical :: thermal_couette = .false. !Computes error
  logical :: init_1drp
  real(kind=DOUBLE) :: x1drp
  real(kind=DOUBLE), dimension(5) :: sol_w_1drp_l, sol_w_1drp_r

  !Restart
  logical :: init_restart = .FALSE.
  integer(kind=ENTIER) :: restart_iter=0
  real(kind=DOUBLE) :: restart_time=0.0, restart_cpu_time=0.0
  character(len=255) :: restart_file
  integer(kind=ENTIER) :: id_vtk_restart = 0

  logical :: init_sedov = .FALSE.
  real(kind=DOUBLE) :: r_sedov

  logical :: init_kelvin = .FALSE.
  logical :: init_double_mach = .FALSE.

  ! Triple-point shock interaction (single material, gamma = 1.4)
  logical :: init_triple_point = .FALSE.

  !Writes a file containing for each vertex &
  !the cell size associated (for salome adaptation)
  logical :: write_cell_size = .FALSE.
  real(kind=DOUBLE) :: delta_mach_imp=1.0

  !Output on surface tag coeffs_surf
  logical :: compute_coeffs = .FALSE.
  integer(kind=ENTIER) :: coeffs_surf = 0
  real(kind=DOUBLE) :: pinf = 0., rhoinf = 0., vinf = 0.

  !Raw plot of solution for each cell
  logical :: plot_solution_dat = .FALSE.
  real(kind=DOUBLE) :: xmin_dat, xmax_dat
  real(kind=DOUBLE) :: ymin_dat, ymax_dat
  real(kind=DOUBLE) :: zmin_dat, zmax_dat

  !Residual output
  logical :: write_residual = .TRUE.
  integer(kind=ENTIER) :: n_iter_residual=100

  !Print residual every n_iter_print
  integer(kind=ENTIER) :: n_iter_print = 1
  integer(kind=ENTIER) :: n_iter_write_sol = 100

  !Final residual
  real(kind=DOUBLE) :: final_res = 1e-4_DOUBLE
  integer(kind=ENTIER) :: n_max_iter = 10000

  public :: read_input_parameters
  public :: init_sol
  public :: compute_dt
  public :: compute_rhs
  ! public :: test_reconstruction_exactness
  ! public :: compute_cell_moments
  public :: primit_to_conserv
  public :: conserv_to_primit
  public :: print_bcs
  public :: count_elems
  public :: face_geometry
	public :: reconstruct, ls_reconstruction, solve_ls, gauss_solve

  ! Cache of each cell's own second geometric moment about its centroid,
  ! M_jk(i) = (1/V_i) * int_cell (x_j - xc_j)(x_k - xc_k) dV -- a pure
  ! mesh-geometry quantity, independent of order/field, computed once by
  ! compute_cell_moments (called from euler_ho_main right after
  ! compute_geometry_mesh) and read by reconstruct() at order 3. See the
  ! header comment above reconstruct() for why it is needed.
  real(kind=DOUBLE), dimension(:, :, :), allocatable :: cell_moment

contains

  subroutine read_input_parameters(filename)
    implicit none

    character(len=*), intent(in) :: filename

    integer(kind=ENTIER) :: funit

    namelist /INPUT_PARAM/ &
      meshfile_path, meshfile, &
      rescale, rescale_factor, &
      periodic_mesh, &
      use_cylinder_map, dim_cylinder_map, &
      n_bc, bc_name, bc_type, bc_val, &
      bc_type_V, bc_val_V, bc_type_T, bc_val_T, &
      n_iter_print, n_iter_write_sol, &
      order, use_aho_reconstruction, &
      second_order, method, &
      scheme, exclude_bound_vert, &
      activate_diffusion, &
      local_time_step, &
      bc_style, boundary_2d, &
      write_residual, n_iter_residual, &
      final_res, n_max_iter, &
      init_uniform, sol_uniform, &
      init_sedov, r_sedov, &
      init_isentropic_vortex, &
      init_gresho, &
      init_potential_flow_2d, &
      init_potential_flow_3d, &
      init_kelvin, &
      init_double_mach, &
      init_triple_point, &
      init_1drp, sol_w_1drp_l, sol_w_1drp_r, x1drp, &
      compute_error, error_2d, error_2d_h, &
      n_iter_write_sol, &
      init_restart, restart_file, id_vtk_restart, &
      compute_coeffs, pinf, rhoinf, vinf, coeffs_surf, &
      plot_solution_dat, &
      xmin_dat, xmax_dat, &
      ymin_dat, ymax_dat, &
      zmin_dat, zmax_dat, &
      write_cell_size, &
      delta_mach_imp, &
      mu_p, Cv_p, use_sutherland, &
      mu0, T0, &
      cfl, tmax, &
      thermal_couette

    open (newunit=funit, file=trim(adjustl(filename)))
    read (nml=INPUT_PARAM, unit=funit)
    close (unit=funit)
  end subroutine read_input_parameters

  pure function primit_to_conserv(w) result(u)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE), dimension(5) :: u

    u(1) = w(1)
    u(2:4) = w(2:4)*w(1)
    u(5) = w(5)/(gamma - 1) &
      + 0.5_DOUBLE*w(1)*(w(2)**2 + w(3)**2 + w(4)**2)
  end function primit_to_conserv

  pure function conserv_to_primit(u) result(w)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: u
    real(kind=DOUBLE), dimension(5) :: w

    w(1) = u(1)
    w(2:4) = u(2:4)/u(1)
    w(5) = (gamma - 1)*(u(5) &
      - 0.5_DOUBLE*w(1)*(w(2)**2 + w(3)**2 + w(4)**2))
  end function conserv_to_primit

  pure subroutine sol_gresho_mach(x, w, mach)
    implicit none

    real(kind=DOUBLE), dimension(3), intent(in) :: x
    real(kind=DOUBLE), dimension(5), intent(inout) :: w
    real(kind=DOUBLE), intent(in) :: mach

    real(kind=DOUBLE) :: w_gresho, p_gresho, r
    real(kind=DOUBLE), dimension(3) :: coord2

    w_gresho = 0.2_DOUBLE
    p_gresho = 1.0_DOUBLE/(gamma*mach**2)

    coord2(:) = x - (/0._DOUBLE, 0._DOUBLE, 0._DOUBLE/)
    r = norm2(coord2(:2))

    if (r < w_gresho) then
      w(1) = 1.0_DOUBLE
      w(2) = -5.0_DOUBLE*coord2(2)
      w(3) = 5.0_DOUBLE*coord2(1)
      w(4) = 0.0_DOUBLE
      w(5) = p_gresho + 12.5_DOUBLE*r**2
    else if (r > 1.0_DOUBLE*w_gresho .and. r < 2.0_DOUBLE*w_gresho) then
      w(1) = 1.0_DOUBLE
      w(2) = (5.0_DOUBLE - 2.0_DOUBLE/r)*coord2(2)
      w(3) = (-5.0_DOUBLE + 2.0_DOUBLE/r)*coord2(1)
      w(4) = 0.0_DOUBLE
      w(5) = p_gresho + 12.5_DOUBLE*r**2 + &
        4.0_DOUBLE - 20.0_DOUBLE*r + 4.0_DOUBLE*log(5.0_DOUBLE*r)
    else
      w(1) = 1.0_DOUBLE
      w(2) = 0.0_DOUBLE
      w(3) = 0.0_DOUBLE
      w(4) = 0.0_DOUBLE
      w(5) = p_gresho - 2.0_DOUBLE + 4.0_DOUBLE*log(2.0_DOUBLE)
    end if
  end subroutine sol_gresho_mach

  pure subroutine sol_isentropic_vortex(coord, w, t)
    implicit none

    real(kind=DOUBLE), intent(in) :: t
    real(kind=DOUBLE), dimension(3), intent(in) :: coord
    real(kind=DOUBLE), dimension(5), intent(inout) :: w

    real(kind=DOUBLE), dimension(3) :: center_coord, vel
    real(kind=DOUBLE) :: beta, r

    beta = 5.0_DOUBLE
    vel(:) = (/0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE/)
    center_coord(:) = (/0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE/) + vel*t
    r = (coord(1) - center_coord(1))**2 + (coord(2) - center_coord(2))**2
    w(1) = 1.0_DOUBLE*(1.0_DOUBLE - ((gamma - 1)*beta**2)/(8.0_DOUBLE*gamma*pi**2)* &
      exp(1.0_DOUBLE - r))**(1.0_DOUBLE/(gamma - 1.0_DOUBLE))
    w(2) = vel(1) &
      - (coord(2) - center_coord(2))*beta/(2.0_DOUBLE*pi)*exp(0.5_DOUBLE*(1.0_DOUBLE - r))
    w(3) = vel(2) &
      + (coord(1) - center_coord(1))*beta/(2.0_DOUBLE*pi)*exp(0.5_DOUBLE*(1.0_DOUBLE - r))
    w(4) = 0.0_DOUBLE
    w(5) = w(1)**gamma
  end subroutine sol_isentropic_vortex

  subroutine init_sol(mesh, sol)
    type(mesh_type), intent(in)    :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(out) :: sol

    integer(kind=ENTIER) :: i, kf, iface_loc, n_fv, k, q
    real(kind=DOUBLE), dimension(5) :: w, w_sum, wq
    real(kind=DOUBLE), dimension(3) :: xc
    real(kind=DOUBLE), dimension(:, :), allocatable :: qpts_v, fc
    real(kind=DOUBLE), dimension(:),    allocatable :: qwts_v
    real(kind=DOUBLE) :: area_q

    if (init_uniform) then
      do i = 1, mesh%n_elems
        sol(:, i) = primit_to_conserv(sol_uniform)
      end do
    else if (init_1drp) then
      do i = 1, mesh%n_elems
        if( mesh%elem(i)%coord(1) < x1drp ) then
          sol(:, i) = primit_to_conserv(sol_w_1drp_l)
        else 
          sol(:, i) = primit_to_conserv(sol_w_1drp_r)
        end if
      end do
    else if (init_isentropic_vortex) then
      do i = 1, mesh%n_elems
        call sol_isentropic_vortex(mesh%elem(i)%coord, w, 0.0_DOUBLE)
        sol(:, i) = primit_to_conserv(w)
      end do
    else if (init_gresho) then
      do i = 1, mesh%n_elems
        call sol_gresho_mach(mesh%elem(i)%coord, w, 1.0_DOUBLE)
        sol(:, i) = primit_to_conserv(w)
      end do
    end if
  end subroutine init_sol

  subroutine compute_dt(mesh, sum_lambda, dt)
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: sum_lambda
    real(kind=DOUBLE) :: dt, dt_i
    integer(kind=ENTIER) :: i

    dt = huge(1.0_DOUBLE)
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost .and. sum_lambda(i) > 0.0_DOUBLE) then 
        dt_i = cfl * mesh%elem(i)%volume / sum_lambda(i)
        ! print*, 'elem=', i, 'sum_lambda=', sum_lambda(i), 'dt_i=',dt_i
        dt = min(dt, dt_i)
      end if
    end do
  end subroutine compute_dt

  subroutine compute_rhs(mesh, sol, prim, rhs, sum_lambda, t)
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in)  :: sol
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in)  :: prim
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(out) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems),    intent(out) :: sum_lambda
    ! Current stage's time, for a time-dependent BC (e.g. 'dmr_top'). Not
    ! tracked per-RK-substage (all 3 SSP-RK3 stages of one step reuse the
    ! step's start time) -- a minor, deliberate simplification for a
    ! qualitative test, not exact substage timing.
    real(kind=DOUBLE), intent(in) :: t

    ! Per-cell gradients (5 primitives x 3 spatial dims)
    real(kind=DOUBLE), allocatable :: grad(:, :, :)  ! (5, 3, n_elems)
    real(kind=DOUBLE), allocatable :: hess(:, :, :, :) ! (5, 3, 3, n_elems)

    allocate(grad(5, 3, mesh%n_elems))
    grad = 0.0_DOUBLE

    if (order >= 3) then
      allocate(hess(5, 3, 3, mesh%n_elems))
      hess = 0.0_DOUBLE
    end if

    if (order >= 2) then
      ! if (use_aho_reconstruction) then
      !   call aho_reconstruction(mesh, prim, grad, hess)
      ! else
        call ls_reconstruction(mesh, prim, grad, hess)
      ! end if
    end if

    rhs        = 0.0_DOUBLE
    sum_lambda = 0.0_DOUBLE
    call face_flux_loop(mesh, sol, prim, grad, hess, rhs, sum_lambda, t)

    ! Divide by cell volume
    block
      integer(kind=ENTIER) :: i
      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) then
          rhs(:, i) = rhs(:, i) / mesh%elem(i)%volume
        end if
      end do
    end block

  end subroutine compute_rhs

  subroutine face_flux_loop(mesh, sol, prim, grad, hess, rhs, sum_lambda, t)
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in)  :: sol, prim
    real(kind=DOUBLE), dimension(5, 3, mesh%n_elems), intent(in) :: grad
    real(kind=DOUBLE), dimension(:, :, :, :), allocatable, intent(in) :: hess
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems),    intent(inout) :: sum_lambda
    real(kind=DOUBLE), intent(in) :: t

    integer(kind=ENTIER) :: iface, il, ir, iv, k, n_fvert, n_qpts, q
    integer(kind=ENTIER) :: iloop
    real(kind=DOUBLE), dimension(3) :: norm, xface
    real(kind=DOUBLE), dimension(:, :), allocatable :: face_coords, qpts
    real(kind=DOUBLE), dimension(:),    allocatable :: qwts
    real(kind=DOUBLE), dimension(5) :: wL, wR, flux
    real(kind=DOUBLE) :: lambda
    logical :: is_zface

    do iface = 1, mesh%n_faces
      il   = mesh%face(iface)%left_neigh
      ir   = mesh%face(iface)%right_neigh
      norm = mesh%face(iface)%norm

      ! Skip ghost-owned faces
      if (mesh%elem(il)%is_ghost) cycle

      n_fvert = mesh%face(iface)%n_vert

      ! For z-faces on 2D-extruded meshes (norm_z ≈ ±1, boundary_2d=T), the
      ! solution has no z-variation. A single centroid point is exact AND
      ! avoids spurious net z-flux from GMSH's vertex ordering (top/bottom
      ! faces get different (x,y) Gauss points due to opposite windings).
      is_zface = boundary_2d .and. abs(abs(norm(3)) - 1.0_DOUBLE) < 1.0e-6_DOUBLE
      if (order >= 2 .and. n_fvert > 0 .and. allocated(mesh%face(iface)%vert) &
          .and. .not. is_zface) then
        ! Quadrature points on face
        allocate(face_coords(3, n_fvert))
        do k = 1, n_fvert
          iv = mesh%face(iface)%vert(k)
          face_coords(:, k) = mesh%vert(iv)%coord
        end do
        call face_quad_pts(int(n_fvert, ENTIER), face_coords, &
          int(order, ENTIER), qpts, qwts)
        deallocate(face_coords)
      else
        ! Order 1 or z-face (2D extrusion): single point at face centroid
        n_qpts = 1
        allocate(qpts(3, 1), qwts(1))
        qpts(:, 1) = mesh%face(iface)%coord
        qwts(1)    = mesh%face(iface)%area
      end if

      ! For z-faces in 2D mode: zero net physics, skip entirely
      if (is_zface) then
        deallocate(qpts, qwts)
        cycle
      end if

      n_qpts = size(qwts)
      do q = 1, n_qpts
        xface = qpts(:, q)

        wL = reconstruct(prim, grad, hess, il, xface, mesh%elem(il)%coord)
        if (ir > 0) then
          wR = reconstruct(prim, grad, hess, ir, xface, mesh%elem(ir)%coord)
        else
          ! Boundary: reconstruct left then apply BC
          ! here is outflowsupersonic
          wR = wL
          !wR = ghost_prim(xface, norm, wL, -ir, t)
        end if

        ! Numerical flux (physical: per-area * area_weight_at_qpt).
        ! flux_scheme_id resolved once in read_params -- see its
        ! declaration for why this is an integer compare, not a string one.
        ! select case (flux_scheme_id)
        ! case (FLUX_THREE_WAVE)
        !   flux = three_wave_flux(wL, wR, norm) * qwts(q)
        ! case (FLUX_TWO_WAVE)
        ! flux = two_wave_flux(wL, wR, norm) * qwts(q)
        flux = rusanov(wL, wR, norm) * qwts(q)
        
        lambda = max(abs(dot_product(wL(2:4), norm)) + cs(wL), &
                     abs(dot_product(wR(2:4), norm)) + cs(wR)) * qwts(q)

        ! print*, 'q=', q, 'qwts(q)=', qwts(q), 'lam=', lambda
        ! print*, 'rho=', wL(1), 'u=', wL(2:4), 'p=', wL(5), 'cs=', cs(wL)

        rhs(:, il)  = rhs(:, il)  - flux
        sum_lambda(il) = sum_lambda(il) + lambda
        if (ir > 0) then
          rhs(:, ir)     = rhs(:, ir)     + flux
          sum_lambda(ir) = sum_lambda(ir) + lambda
        end if
      end do

      deallocate(qpts, qwts)
    end do
  end subroutine face_flux_loop

  subroutine face_geometry(mesh)
    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER) :: iface, il, ir
    real(kind=DOUBLE), dimension(3) :: norm
    real(kind=DOUBLE) :: area, nabs
  
    do iface = 1, mesh%n_faces
      il   = mesh%face(iface)%left_neigh
      ir   = mesh%face(iface)%right_neigh
      ! Skip ghost-owned faces
      if (mesh%elem(il)%is_ghost) cycle

      norm = mesh%face(iface)%norm
      nabs = sqrt(norm(1)**2+norm(2)**2+norm(3)**2)
      area = mesh%face(iface)%area
      print*, 'face area=', area, 'norm abs=', nabs
    end do
  end subroutine face_geometry

  ! Ghost cell primitive state for boundary condition
  ! function ghost_prim(xf, norm, wL, id_bc, t) result(wR)
  !   real(kind=DOUBLE), dimension(3), intent(in) :: xf, norm
  !   real(kind=DOUBLE), dimension(5), intent(in) :: wL
  !   integer(kind=ENTIER), intent(in) :: id_bc
  !   real(kind=DOUBLE), intent(in) :: t
  !   real(kind=DOUBLE), dimension(5) :: wR

  !   real(kind=DOUBLE) :: vn

  !   if (id_bc < 1 .or. id_bc > n_bc) then
  !     ! Default: slip wall (mirror normal velocity)
  !     wR    = wL
  !     vn    = dot_product(wL(2:4), norm)
  !     wR(2) = wL(2) - 2.0_DOUBLE * vn * norm(1)
  !     wR(3) = wL(3) - 2.0_DOUBLE * vn * norm(2)
  !     wR(4) = wL(4) - 2.0_DOUBLE * vn * norm(3)
  !     return
  !   end if

  !   ! bc_type_id resolved once in read_params -- see its declaration.
  !   select case (bc_type_id(id_bc))
  !   case (BC_EULER_FREESTREAM)
  !     wR = bc_val(:, id_bc)
  !   case (BC_EULER_OUTFLOWSUPERSONIC)
  !     ! Zero-gradient extrapolation: valid where every characteristic
  !     ! leaves the domain (locally supersonic outflow).
  !     wR = wL
  !   case default
  !     ! BC_WALL (also the fallback for blank/unrecognized bc_type,
  !     ! matching the old string select case's behavior)
  !     wR    = wL
  !     vn    = dot_product(wL(2:4), norm)
  !     wR(2) = wL(2) - 2.0_DOUBLE * vn * norm(1)
  !     wR(3) = wL(3) - 2.0_DOUBLE * vn * norm(2)
  !     wR(4) = wL(4) - 2.0_DOUBLE * vn * norm(3)
  !   end select
  ! end function ghost_prim

  pure function euler_flux(w, n) result(F)
    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), dimension(5) :: F
    real(kind=DOUBLE) :: vn, rho, p, e

    rho = w(1); p = w(5)
    vn  = w(2)*n(1) + w(3)*n(2) + w(4)*n(3)
    e   = p / ((gamma - 1.0_DOUBLE) * rho) + 0.5_DOUBLE*(w(2)**2+w(3)**2+w(4)**2)

    F(1) = rho * vn
    F(2) = rho * w(2) * vn + p * n(1)
    F(3) = rho * w(3) * vn + p * n(2)
    F(4) = rho * w(4) * vn + p * n(3)
    F(5) = rho * e * vn + p * vn
  end function euler_flux

  ! Rusanov (local Lax-Friedrichs) numerical flux per unit area
  pure function rusanov(wL, wR, n) result(F)
    real(kind=DOUBLE), dimension(5), intent(in) :: wL, wR
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), dimension(5) :: F
    real(kind=DOUBLE) :: lam, vnL, vnR

    vnL = wL(2)*n(1) + wL(3)*n(2) + wL(4)*n(3)
    vnR = wR(2)*n(1) + wR(3)*n(2) + wR(4)*n(3)
    lam = max(abs(vnL) + cs(wL), abs(vnR) + cs(wR))

    F = 0.5_DOUBLE * (euler_flux(wL, n) + euler_flux(wR, n)) &
      - 0.5_DOUBLE * lam * (primit_to_conserv(wR) - primit_to_conserv(wL))
  end function rusanov

  pure function two_wave_flux(wL, wR, n) result(F)
    real(kind=DOUBLE), dimension(5), intent(in) :: wL, wR
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), dimension(5) :: F

    real(kind=DOUBLE), dimension(5) :: uL, uR, fL, fR, u_star
    real(kind=DOUBLE) :: rhoL, rhoR, vnL, vnR, pL, pR, aL, aR
    real(kind=DOUBLE) :: lambdaL, lambdaR, sL_wave, sR_wave

    rhoL = wL(1); vnL = dot_product(wL(2:4), n); pL = wL(5)
    aL = cs(wL); uL = primit_to_conserv(wL)

    rhoR = wR(1); vnR = dot_product(wR(2:4), n); pR = wR(5)
    aR = cs(wR); uR = primit_to_conserv(wR)

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

  ! Sound speed from primitive state
  pure function cs(w) result(c)
    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE) :: c
    c = sqrt(max(gamma * w(5) / max(w(1), 1.0e-16_DOUBLE), 0.0_DOUBLE))
  end function cs

  pure function sound_speed_w(w) result(a)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE) :: a

    a = sqrt(gamma*w(5)/w(1))
  end function sound_speed_w

  subroutine print_bcs
    implicit none
    
    type(mesh_type) :: mesh
    integer(kind=ENTIER) :: i, le, re, i_bc

    do i=1, mesh%n_faces
    le = mesh%face(i)%left_neigh
    re = mesh%face(i)%right_neigh
      if(re <= 0) then
        if(re < 0) then
          print*, "boundary ", i, -re, trim(adjustl(bc_name(-re))), trim(adjustl(bc_type(-re)))
        end if
      end if
    end do
  end subroutine print_bcs

  subroutine count_elems(mesh,n_elems_loc,n_elems_ghost)
    implicit none
    
    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(inout) :: n_elems_loc,n_elems_ghost
    integer(kind=ENTIER) :: i

    n_elems_loc=0
    n_elems_ghost=0
    do i=1,mesh%n_elems
      if (mesh%elem(i)%is_ghost) then
        n_elems_ghost=n_elems_ghost+1
        cycle
      end if
      n_elems_loc=n_elems_loc+1
    end do
  end subroutine count_elems

  ! Polynomial reconstruction of primitive variable w at point xq from cell i.
  ! Hierarchical fallback: if an order-p reconstruction gives unphysical
  ! rho or p, it is replaced by the order-(p-1) result.
  function reconstruct(prim, grad, hess, i, xq, xc) result(w)
    real(kind=DOUBLE), dimension(:, :),          intent(in) :: prim    ! (5, n_elems)
    real(kind=DOUBLE), dimension(:, :, :),       intent(in) :: grad    ! (5, 3, n_elems)
    real(kind=DOUBLE), dimension(:, :, :, :), allocatable, intent(in) :: hess
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE), dimension(3), intent(in) :: xq, xc
    real(kind=DOUBLE), dimension(5) :: w

    real(kind=DOUBLE), dimension(3) :: dx
    real(kind=DOUBLE), dimension(5) :: w_try
    integer(kind=ENTIER) :: j, k

    dx = xq - xc
    w  = prim(:, i)

    if (order >= 2) then
      w_try = w + matmul(grad(:, :, i), dx)
      if (physical_state(w_try, w)) w = w_try
    end if

    if (order >= 3 .and. allocated(hess)) then
      if (.not. allocated(cell_moment)) then
        print *, 'FATAL: reconstruct() called at order>=3 but compute_cell_moments ', &
          'was never called -- the order-3 Taylor polynomial would silently carry the ', &
          'wrong cell average (see header comment above compute_cell_moments). Call ', &
          'compute_cell_moments(mesh) once, right after compute_geometry_mesh.'
        error stop 1
      end if
      ! Always build from scratch: cell value + gradient + Hessian, MINUS
      ! the cell's own curvature-average bias (see compute_cell_moments)
      ! so that this polynomial's cell average is prim(:,i), not its
      ! value at xc plus a spurious O(h^2) offset.
      w_try = prim(:, i) + matmul(grad(:, :, i), dx)
      do j = 1, 3
        do k = 1, 3
          w_try = w_try + 0.5_DOUBLE * hess(:, j, k, i) &
            * (dx(j) * dx(k) - cell_moment(j, k, i))
        end do
      end do
      if (physical_state(w_try, prim(:, i))) w = w_try
    end if

    ! Enforce positivity (safety net)
    w(1) = max(w(1), 1.0e-12_DOUBLE)
    w(5) = max(w(5), 1.0e-12_DOUBLE)
  end function reconstruct

	! Returns .true. if w_cand is a physically acceptable reconstruction
  ! relative to the reference state w_ref.
  ! Rejects: negative rho or p; velocity more than 20x the reference speed + c.
  pure function physical_state(w_cand, w_ref) result(ok)
    real(kind=DOUBLE), dimension(5), intent(in) :: w_cand, w_ref
    logical :: ok
    real(kind=DOUBLE) :: spd_ref, spd_cand, c_ref

    ok = .false.
    if (w_cand(1) <= 0.0_DOUBLE) return
    if (w_cand(5) <= 0.0_DOUBLE) return

    ! Velocity magnitude check
    spd_ref  = w_ref(2)**2  + w_ref(3)**2  + w_ref(4)**2
    spd_cand = w_cand(2)**2 + w_cand(3)**2 + w_cand(4)**2
    c_ref    = gamma * w_ref(5) / max(w_ref(1), 1.0e-16_DOUBLE)
    if (spd_cand > 400.0_DOUBLE * (spd_ref + c_ref)) return   ! 20x speed limit

    ok = .true.
  end function physical_state

	! ================================================================
  ! Internal subroutines
  ! ================================================================

  ! Cell-centred least-squares polynomial reconstruction.
  !
  ! For each non-ghost cell i, collects all vertex-neighbours (cells sharing
  ! at least one vertex), then fits a degree-(order-1) polynomial to their
  ! mean values by solving the normal equations (A^T A) x = A^T b.
  !
  ! Polynomial basis (boundary_2d=T, i.e. 2D):
  !   order 2: [dx, dy]                                    (2 coefficients)
  !   order 3: [dx, dy, dx²/2, dx·dy, dy²/2]              (5 coefficients)
  ! Full 3D:
  !   order 2: [dx, dy, dz]                                (3 coefficients)
  !   order 3: [dx, dy, dz, dx²/2, dx·dy, dx·dz,
  !              dy²/2, dy·dz, dz²/2]                      (9 coefficients)
  !
  ! Coefficients are stored in the existing grad/hess arrays so that the
  ! reconstruct() function needs no changes.
  subroutine ls_reconstruction(mesh, prim, grad, hess)
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: prim
    real(kind=DOUBLE), dimension(5, 3, mesh%n_elems), intent(inout) :: grad
    real(kind=DOUBLE), dimension(:, :, :, :), allocatable, intent(inout) :: hess

    integer(kind=ENTIER) :: i, iv, mf, iface, il, ir, jn, kv, nv
    integer(kind=ENTIER) :: n_coeff, n_neigh, n_needed
    integer(kind=ENTIER) :: ring_start, ring_end, ic, ivv
    integer(kind=ENTIER), dimension(500) :: neigh_list
    real(kind=DOUBLE), allocatable :: Amat(:, :)
    real(kind=DOUBLE), dimension(500) :: wt_arr   ! sqrt(1/dist) per neighbour
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

      ! Collect unique real vertex-neighbours of cell i (ring 1)
      do nv = 1, mesh%elem(i)%n_vert
        iv = mesh%elem(i)%vert(nv)
        do mf = 1, mesh%vert(iv)%n_faces_neigh
          iface = mesh%vert(iv)%face_neigh(mf)
          call add_unique_neighbor(mesh%face(iface)%left_neigh)
          call add_unique_neighbor(mesh%face(iface)%right_neigh)
        end do
      end do

      ! Ring-expand (neighbours of neighbours, etc.) when the immediate
      ! vertex-neighbour ring doesn't have enough cells to determine the
      ! unknowns actually alive on THIS mesh -- e.g. a nominal order-3 fit
      ! has 5 unknowns (boundary_2d: dx,dy,dxx,dxy,dyy), but on a
      ! single-row-in-y mesh (subfv's Sod tube convention) every neighbour
      ! shares the same y, so dy/dxy/dyy are identically zero columns and
      ! only 2 unknowns (dx,dxx) are really determinable -- n_needed
      ! (below) tracks that reduced count, recomputed as the ring grows,
      ! so the stencil stays as narrow (and thus as locally accurate) as
      ! the data actually require instead of always demanding the full,
      ! nominal n_coeff neighbours even when most of those unknowns are
      ! degenerate. Previously this used a fixed n_coeff target and simply
      ! gave up (`cycle`, leaving grad=hess at their initialized 0) when
      ! the immediate ring fell short of it, silently degrading every
      ! order-3 ls_reconstruction cell to order 1 on the whole Sod-tube
      ! mesh family. Mirrors arbitrary_high_order_module's
      ! gather_ls_neighbors ring growth; the degeneracy detection mirrors
      ! its compute_nodal_derivative_at_vertex direction-spread check.
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
        if (ring_end == n_neigh) exit   ! stencil can't grow further (tiny/disconnected mesh)
        ring_start = ring_end + 1
      end do

      if (n_neigh < n_needed) cycle   ! still under-determined after ring expansion: leave grad=0

      allocate(Amat(n_neigh, n_coeff))

      ! Characteristic length for THIS cell, used to nondimensionalize the
      ! design matrix columns before solving (see Lref note below) --
      ! cube root of the cell volume is a reasonable, cheap proxy for its
      ! own size regardless of cell shape.
      Lref = max(mesh%elem(i)%volume, 1.0e-300_DOUBLE)**(1.0_DOUBLE/3.0_DOUBLE)

      ! Build the design matrix A weighted by w_j = 1/dist_j (inverse-distance).
      ! Row j of Amat is scaled by sqrt(w_j) = 1/sqrt(dist_j) so that
      ! ATA = A^T W A  and  ATb = A^T W b  (with W = diag(w_j)). Columns are
      ! built from ddx/Lref etc (order-1 in the SCALED offset), not raw
      ! ddx -- found necessary 2026-09-14: with raw offsets, the linear-term
      ! columns (~dx*wt) and quadratic-term columns (~dx^2*wt) differ in
      ! scale by a factor of ~1/dx, which at a fine mesh (dx=0.01) inflates
      ! ATA's condition number enough to genuinely corrupt the solve (order-3
      ! ls_reconstruction was found NOT exact -- Linf 5e-4 instead of machine
      ! precision -- on a degree-2 test polynomial on a genuinely-3D
      ! pseudo-1D mesh; a plain absolute-vs-relative pivot-tolerance fix,
      ! tried first, made no difference, pointing at conditioning rather
      ! than a threshold). x_coeff comes back in these SCALED units and is
      ! divided back down by the matching power of Lref where grad/hess are
      ! assigned below.
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
        wt = 1.0_DOUBLE / sqrt(max(dist, 1.0e-14_DOUBLE))   ! sqrt(1/dist)
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

      ! Weighted normal matrix A^T W A  (n_coeff x n_coeff)
      ATA(1:n_coeff, 1:n_coeff) = matmul(transpose(Amat), Amat)

      ! Solve once per primitive variable (ATb = A^T W b)
      do kv = 1, 5
        ATb(1:n_coeff) = 0.0_DOUBLE
        do jn = 1, n_neigh
          ATb(1:n_coeff) = ATb(1:n_coeff) + Amat(jn, :) * wt_arr(jn) * &
            (prim(kv, neigh_list(jn)) - prim(kv, i))
        end do

        call solve_ls(ATA, ATb, x_coeff, n_coeff)

        ! x_coeff came back in Lref-scaled units (see Amat construction
        ! above): d/dx = (1/Lref) d/d(sdx), d^2/dx^2 = (1/Lref^2) d^2/d(sdx)^2.
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

    ! Add cell `cand` to neigh_list(1:n_neigh) (host-associated, declared
    ! above in ls_reconstruction) if it is a real, non-ghost cell other
    ! than i itself and not already present.
    subroutine add_unique_neighbor(cand)
      integer(kind=ENTIER), intent(in) :: cand
      integer(kind=ENTIER) :: jn
      logical :: found

      ! if (cand <= 0 .or. cand == i .or. mesh%elem(cand)%is_ghost) return
      if (cand <= 0 .or. cand == i) return
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

    ! Number of unknowns actually determinable from the CURRENT
    ! neigh_list(1:n_neigh) (host-associated), given which raw coordinate
    ! directions have any spread among these neighbours relative to i.
    ! x is assumed always alive (a mesh degenerate even in x is not
    ! handled specially). A degenerate y drops y itself and every
    ! quadratic term involving it (xy, yy); same for z. Reduces exactly
    ! to the nominal n_coeff (5 or 9 for order>=3, 2 or 3 otherwise) on
    ! any ordinary 2D/3D mesh where all directions are alive.
    function count_needed_coeffs() result(n_needed)
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

  ! Solves ATA_in * x = rhs_in (ATA_in = A^T W A, a Gram/normal matrix,
  ! always symmetric PSD). Drops directions with no information in the
  ! data -- e.g. the y-gradient unknown on a mesh that is a single row of
  ! cells in y (ls_reconstruction's Sod-tube mesh: boundary_2d, every
  ! neighbour at the same y, so the whole y-row/column of ATA is exactly
  ! zero) -- and solves only the reduced, well-posed subsystem for the
  ! rest, rather than giving up on the whole vector.
  !
  ! BUG FIXED (previously): a plain Gaussian elimination with partial
  ! pivoting returned x=0 for EVERY unknown, including well-determined
  ! ones, the instant it hit any near-zero pivot -- so on the Sod mesh
  ! above, the perfectly-determined x-gradient was ALSO silently zeroed
  ! alongside the genuinely-degenerate y one, making ls_reconstruction
  ! produce bit-identical output to plain order 1 at every order (2 and
  ! 3) on that whole mesh family. Because ATA is a Gram matrix, a
  ! genuinely-degenerate direction j has its ENTIRE row/column j equal to
  ! zero (not just a small pivot reached after row operations), so this
  ! is detected directly from the untouched input, up front, rather than
  ! discovered mid-elimination -- and only actually-degenerate directions
  ! are dropped.
  subroutine solve_ls(ATA_in, rhs_in, x, n)
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

    ! Relative to the matrix's own scale, same reasoning as gauss_solve's
    ! piv_tol: an absolute floor here misclassifies a genuinely-alive but
    ! naturally-small quadratic-term column (scale ~dx^2*wt, vs ~dx*wt for
    ! the linear-term columns) as degenerate purely from that scale gap.
    col_tol = 1.0e-12_DOUBLE * maxval(abs(ATA_in(1:n, 1:n)))

    n_keep = 0
    do jj = 1, n
      if (any(abs(ATA_in(1:n, jj)) > col_tol)) then
        n_keep = n_keep + 1
        idx(n_keep) = jj
      end if
    end do

    ! Iteratively solve the currently-active subset (idx(1:n_keep)); if
    ! gauss_solve hits a genuinely-degenerate pivot MID-elimination (a
    ! direction that looked alive column-wise above but turns out to be a
    ! linear combination of others once the others are eliminated against
    ! it -- e.g. a symmetric structured-mesh stencil where some quadratic
    ! cross-term direction is exactly indeterminate), drop just that one
    ! direction and retry with the rest, rather than discarding the whole
    ! vector. Found necessary 2026-09-14: on a genuinely-3D pseudo-1D mesh
    ! (regular hex cells), the old single-shot gauss_solve threw away
    ! EVERY coefficient, including perfectly well-determined ones like the
    ! plain x-gradient, the moment ANY one of the 9 order-3 directions hit
    ! a degenerate pivot -- caught because order-3 ls_reconstruction was
    ! found not exact (Linf 5e-4, not machine precision) on a degree-2
    ! test polynomial on that mesh.
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

      ! Drop the direction that failed and retry with the rest.
      do ii = fail_at, n_keep - 1
        idx(ii) = idx(ii + 1)
      end do
      n_keep = n_keep - 1
    end do
  end subroutine solve_ls

  ! Gaussian elimination with partial pivoting for small dense systems.
  ! Solves ATA_in * x = rhs_in. On a degenerate pivot at step jj, returns
  ! immediately with fail_at=jj (x left at 0) instead of guessing --
  ! solve_ls's caller loop drops direction jj and retries with the rest,
  ! rather than discarding every other, well-determined direction along
  ! with the one genuinely-degenerate one (see solve_ls's own comment).
  ! fail_at=0 signals a complete, successful solve.
  subroutine gauss_solve(ATA_in, rhs_in, x, n, fail_at)
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

    ! Degeneracy threshold RELATIVE to the matrix's own scale, not an
    ! absolute constant -- an absolute floor misclassifies a legitimately
    ! nonzero but naturally-small pivot (e.g. from a quadratic-term
    ! column, scale ~dx^2*wt, vs ~dx*wt for linear-term columns) as
    ! degenerate purely from that scale gap, not from true rank-deficiency.
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

end module silvia_base_module