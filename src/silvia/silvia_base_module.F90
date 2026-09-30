module silvia_base_module
  use precision_module
  use mesh_module
  use quadrature_module
  use mpi_module
  use arbitrary_high_order_module, only: compute_next_order_derivative, &
    aho_use_green_gauss => use_green_gauss
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
  character(len=255), dimension(mnbc) :: bc_type = ""
  real(kind=DOUBLE), dimension(5, mnbc) :: bc_val
  ! bc_val_V for Viscous (u, v, w)
  character(len=255), dimension(mnbc) :: bc_type_V = ""
  real(kind=DOUBLE), dimension(3, mnbc) :: bc_val_V
  ! bc_val_T for Heat (T)
  character(len=255), dimension(mnbc) :: bc_type_T = ""
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
  public :: compute_cell_moments
  public :: primit_to_conserv
  public :: conserv_to_primit
  public :: print_bcs
  public :: count_elems
  public :: face_geometry
	public :: reconstruct, aho_reconstruction

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

  subroutine init_flags()
    implicit none

    integer :: i, p1, p2
    character(len=255) :: t, t_adv, t_lag
    logical :: recognized

    allocate(bc_V_id(n_bc), bc_T_id(n_bc), bc_euler_id(n_bc))
    do i = 1, n_bc
      t = trim(adjustl(bc_type_V(i)))
      if (t == "dirichlet") then
        bc_V_id(i) = BC_V_DIRICHLET
      else if (t == "neumann" .or. t == "") then
        bc_V_id(i) = BC_V_NEUMANN
      else
        print *, "ERROR: bc_type_V unrecognized for BC ", i, ": '", trim(t), "'"
        error stop
      end if
      t = trim(adjustl(bc_type_T(i)))
      if (t == "dirichlet") then
        bc_T_id(i) = BC_T_DIRICHLET
      else if (t == "neumann" .or. t == "") then
        bc_T_id(i) = BC_T_NEUMANN
      else
        print *, "ERROR: bc_type_T unrecognized for BC ", i, ": '", trim(t), "'"
        error stop
      end if
      t = trim(adjustl(bc_type(i)))
      if (t == "wall" .or. t == "") then
        bc_euler_id(i) = BC_EULER_WALL
      else if (t == "adherence_wall") then
        bc_euler_id(i) = BC_EULER_ADHERENCE_WALL
      else if (t == "freestream") then
        bc_euler_id(i) = BC_EULER_FREESTREAM
      else if (t == "outflowsupersonic") then
        bc_euler_id(i) = BC_EULER_OUTFLOWSUPERSONIC
      else if (t == "inflow_pond") then
        bc_euler_id(i) = BC_EULER_INFLOW_POND
      else if (t == "inout_double_mach") then
        bc_euler_id(i) = BC_EULER_INOUT_DOUBLE_MACH
      else if (t == "double_mach_bottom") then
        bc_euler_id(i) = BC_EULER_DOUBLE_MACH_BOTTOM
      else if (t == "potential_flow_2d") then
        bc_euler_id(i) = BC_EULER_POTENTIAL_FLOW_2D
      else if (t == "potential_flow_3d") then
        bc_euler_id(i) = BC_EULER_POTENTIAL_FLOW_3D
      else
        print *, "ERROR: bc_type unrecognized for BC ", i, ": '", trim(t), "'"
        error stop
      end if
    end do
    ! Parse scheme string to integer ID
    t = trim(adjustl(scheme))
    if (t == "multi_point") then
      scheme_id = SCHEME_MULTI_POINT
    else if (t == "multi_point_iso") then
      scheme_id = SCHEME_MULTI_POINT_ISO
    else if (t == "three_wave") then
      scheme_id = SCHEME_THREE_WAVE
    else if (t == "two_wave") then
      scheme_id = SCHEME_TWO_WAVE
    else if (t == "modified_three_wave") then
      scheme_id = SCHEME_MODIFIED_THREE_WAVE
    else
      print *, "ERROR: scheme unrecognized: '", trim(t), "'"
      error stop
    end if
  end subroutine init_flags

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

  subroutine compute_rhs(mesh, sol, prim, grad, hess, rhs, sum_lambda, t, &
      num_procs, mpi_send_recv)
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in)  :: sol
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in)  :: prim
    real(kind=DOUBLE), dimension(5, 3, mesh%n_elems), intent(inout) :: grad
    real(kind=DOUBLE), dimension(:, :, :, :), allocatable, intent(inout) :: hess
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(out) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems),    intent(out) :: sum_lambda
    real(kind=DOUBLE), intent(in) :: t
    integer, intent(in) :: num_procs
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    integer(kind=ENTIER) :: i

    grad = 0.0_DOUBLE
    hess = 0.0_DOUBLE
    if (order >= 2) call aho_reconstruction(mesh, prim, grad, hess, num_procs, mpi_send_recv)

    !BJ

    rhs        = 0.0_DOUBLE
    sum_lambda = 0.0_DOUBLE
    call subface_flux_loop(mesh, sol, prim, grad, hess, rhs, sum_lambda, t)

    do i = 1, mesh%n_elems
      rhs(:, i) = rhs(:, i) / mesh%elem(i)%volume
    end do
  end subroutine compute_rhs

  subroutine subface_flux_loop(mesh, sol, prim, grad, hess, rhs, sum_lambda, t)
    use ns_euler_rs_module, only: compute_lambdas_and_solve_nodal_velocity_new, &
      multi_point
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in)  :: sol, prim
    real(kind=DOUBLE), dimension(5, 3, mesh%n_elems), intent(in) :: grad
    real(kind=DOUBLE), dimension(:, :, :, :), allocatable, intent(in) :: hess
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems),    intent(inout) :: sum_lambda
    real(kind=DOUBLE), intent(in) :: t

    integer(kind=ENTIER) :: iface, il, ir, iv, n_fvert, n_qpts, q, inode
    integer(kind=ENTIER) :: iloop
    integer(kind=ENTIER) :: i, j, k, id_sub_face, id_face
    real(kind=DOUBLE), dimension(3) :: norm, xface
    real(kind=DOUBLE), dimension(:, :), allocatable :: face_coords, qpts
    real(kind=DOUBLE), dimension(:),    allocatable :: qwts
    real(kind=DOUBLE), dimension(5) :: wL, wR, flux, sol_l, sol_r
    real(kind=DOUBLE), dimension(5,2) :: lr_flux
    real(kind=DOUBLE) :: lambda
    logical :: is_zsubface

    integer(kind=ENTIER) :: ng
    real(kind=DOUBLE), dimension(:, :), allocatable :: sol_w_l, sol_w_r, norm_list
    real(kind=DOUBLE), dimension(:), allocatable :: lambda_l, lambda_r, weight
    real(kind=DOUBLE), dimension(3) :: v_node
    real(kind=DOUBLE) :: vn_nodal, sl, sr, area


    do i = 1, mesh%n_vert

      !Compute LR for every subface
      ng = mesh%vert(i)%n_sub_faces_neigh
      allocate(sol_w_l(5, ng))
      allocate(sol_w_r(5, ng))
      allocate(lambda_l(ng))
      allocate(lambda_r(ng))
      allocate(weight(ng))
      allocate(norm_list(3, ng))

      do j = 1, mesh%vert(i)%n_sub_faces_neigh
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        il   = mesh%sub_face(id_sub_face)%left_elem_neigh
        ir   = mesh%sub_face(id_sub_face)%right_elem_neigh
        norm = mesh%sub_face(id_sub_face)%norm
        area = mesh%sub_face(id_sub_face)%area

        is_zsubface = boundary_2d .and. abs(abs(norm(3)) - 1.0_DOUBLE) < 1.0e-6_DOUBLE
        if (order>=2 .and. .not. is_zsubface) then
          ! Quadrature points on subface
          allocate(face_coords(3, 4)) ! subfaces are always quadrilaterals
          face_coords(:,1) = mesh%vert(i)%coord
          face_coords(:,2) = mesh%vert(i)%coord ! edge midpoint
          face_coords(:,3) = mesh%face(id_face)%coord
          face_coords(:,4) = mesh%vert(i)%coord ! edge midpoint
          call face_quad_pts(4, face_coords, int(order, ENTIER), qpts, qwts)
          deallocate(face_coords)
        else
          allocate(qpts(3, 1), qwts(1))
          qpts(:, 1) = mesh%vert(i)%coord
          qwts(1)    = area
        end if

        ! we use MP scheme only in one quadrature point per subface, for now...
        wL = reconstruct(prim, grad, hess, il, qpts(:,1), mesh%elem(il)%coord)
        if (ir > 0) then
          wR = reconstruct(prim, grad, hess, ir, qpts(:,1), mesh%elem(ir)%coord)
        else
          sol_l = primit_to_conserv(wL)
          call compute_right_state(mesh, id_sub_face, ir, sol_l, sol_r)
          wR = conserv_to_primit(sol_r)
        end if

        sol_w_l(:, j) = wL
        sol_w_r(:, j) = wR
        lambda_l(j) = 0.0_DOUBLE
        lambda_r(j) = 0.0_DOUBLE
        weight(j) = area
        norm_list(:, j) = norm

        deallocate(qpts, qwts)
      end do

      !Solve the nodal system 
      call compute_lambdas_and_solve_nodal_velocity_new(ng, weight,&
         norm_list, lambda_l, lambda_r, sol_w_l, sol_w_r, v_node)

      !Compute the multi_point flux
      do j = 1, mesh%vert(i)%n_sub_faces_neigh
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        il   = mesh%sub_face(id_sub_face)%left_elem_neigh
        ir   = mesh%sub_face(id_sub_face)%right_elem_neigh
        norm = mesh%sub_face(id_sub_face)%norm
        area = mesh%sub_face(id_sub_face)%area

        ! First quadrature point corresponding to MP scheme, attached to vertex
        select case (scheme_id)
          case (1)
            ! Same as ns: zero normal nodal velocity on wall sub-faces
            ! (re == 0 is the default wall in compute_right_state)
            vn_nodal = dot_product(v_node, norm)
            if (ir == 0) then
              vn_nodal = 0.0_DOUBLE
            else if (ir < 0) then
              if (bc_euler_id(-ir) == BC_EULER_WALL) vn_nodal = 0.0_DOUBLE
            end if
            call multi_point(sol_w_l(:, j), sol_w_r(:, j), norm, lr_flux, vn_nodal, &
              lambda_l(j), lambda_r(j), sl, sr)
            ! lr_flux is per unit area; lr_flux(:, 2) is already the right cell's
            ! outgoing flux (multi_point negates it)
            lambda = max(abs(dot_product(wL(2:4), norm)) + cs(wL), &
                  abs(dot_product(wR(2:4), norm)) + cs(wR)) * mesh%sub_face(id_sub_face)%area
            rhs(:, il)  = rhs(:, il)  - area*lr_flux(:, 1)
            sum_lambda(il) = sum_lambda(il) + max(0.0_DOUBLE, -sl)*area
            if (ir > 0) then
              rhs(:, ir)     = rhs(:, ir)     - area*lr_flux(:, 2)
              sum_lambda(ir) = sum_lambda(ir) + max(0.0_DOUBLE, sr)*area
            end if
          case (4)
            wL = sol_w_l(:, j)
            wR = sol_w_r(:, j)
            flux = rusanov(wL, wR, norm) * mesh%sub_face(id_sub_face)%area
            rhs(:, il)  = rhs(:, il) - flux
            lambda = max(abs(dot_product(wL(2:4), norm)) + cs(wL), &
                  abs(dot_product(wR(2:4), norm)) + cs(wR)) * mesh%sub_face(id_sub_face)%area
            sum_lambda(il) = sum_lambda(il) + lambda
            if (ir > 0) then
              rhs(:, ir)  = rhs(:, ir)  + flux
              sum_lambda(ir) = sum_lambda(ir) + lambda
            end if
          case default
            wL = sol_w_l(:, j)
            wR = sol_w_r(:, j)
            flux = rusanov(wL, wR, norm) * mesh%sub_face(id_sub_face)%area
            rhs(:, il)  = rhs(:, il) - flux
            lambda = max(abs(dot_product(wL(2:4), norm)) + cs(wL), &
                  abs(dot_product(wR(2:4), norm)) + cs(wR)) * mesh%sub_face(id_sub_face)%area
            sum_lambda(il) = sum_lambda(il) + lambda
            if (ir > 0) then
              rhs(:, ir)  = rhs(:, ir)  + flux
              sum_lambda(ir) = sum_lambda(ir) + lambda
            end if
        end select
      end do

      deallocate(sol_w_l, sol_w_r, lambda_l, lambda_r, weight, norm_list)
      ! deallocate(qpts, qwts)
    end do
  end subroutine subface_flux_loop

  subroutine compute_right_state(mesh, id_sub_face, re, sol_l, sol_r)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_sub_face, re
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_l
    real(kind=DOUBLE), dimension(5), intent(inout) :: sol_r

    integer(kind=ENTIER) :: id_face
    real(kind=DOUBLE), dimension(3) :: face_coord
    real(kind=DOUBLE), dimension(5) :: w_pf

    if (re > 0) then
      print*, "re > 0, bad use of compute right sate !"
      error stop
    else if (re == 0) then !Default option wall
      sol_r(:) = sol_l
      sol_r(2:4) = sol_r(2:4) &
        - 2.0_DOUBLE*dot_product(sol_r(2:4), mesh%sub_face(id_sub_face)%norm)*mesh%sub_face(id_sub_face)%norm
    else !!Boundary
      select case (bc_euler_id(-re))
      case (BC_EULER_OUTFLOWSUPERSONIC)
        sol_r = sol_l
      case (BC_EULER_FREESTREAM)
        sol_r = primit_to_conserv(bc_val(:, -re))
      case (BC_EULER_INFLOW_POND)
        sol_r = primit_to_conserv(bc_val(:, -re))
        if (dot_product(sol_r(2:4), mesh%sub_face(id_sub_face)%norm) > 0.0_DOUBLE) then
          sol_r = sol_l
        end if
      case (BC_EULER_WALL)
        sol_r(:) = sol_l
        sol_r(2:4) = sol_r(2:4) &
          - 2.0_DOUBLE*dot_product(sol_r(2:4), mesh%sub_face(id_sub_face)%norm)*mesh%sub_face(id_sub_face)%norm
      case (BC_EULER_ADHERENCE_WALL)
        sol_r(:) = sol_l
        sol_r(2:4) = 2.0_DOUBLE*sol_l(1)*bc_val(2:4, -re) - sol_l(2:4)
        sol_r(5) = sol_l(5) - 0.5_DOUBLE*sol_l(1)*norm2(sol_l(2:4)/sol_l(1))**2 &
          + 0.5_DOUBLE*sol_l(1)*norm2(sol_r(2:4)/sol_r(1))**2
      case (BC_EULER_INOUT_DOUBLE_MACH)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        face_coord = mesh%face(id_face)%coord
        if (face_coord(2) > 1.732_DOUBLE*(face_coord(1) - 0.1667_DOUBLE - 10.0_DOUBLE*t)) then
          sol_r = primit_to_conserv((/8.0_DOUBLE, 7.145_DOUBLE, -4.125_DOUBLE, 0.0_DOUBLE, 116.5_DOUBLE/))
        else
          sol_r = primit_to_conserv((/1.4_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE/))
        end if
      case (BC_EULER_DOUBLE_MACH_BOTTOM)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        face_coord = mesh%face(id_face)%coord
        if (face_coord(1) < 0.1667_DOUBLE) then
          sol_r = primit_to_conserv((/8.0_DOUBLE, 7.145_DOUBLE, -4.125_DOUBLE, 0.0_DOUBLE, 116.5_DOUBLE/))
        else
          sol_r(:) = sol_l
          sol_r(2:4) = sol_r(2:4) &
            - 2.0_DOUBLE*dot_product(sol_r(2:4), mesh%sub_face(id_sub_face)%norm)*mesh%sub_face(id_sub_face)%norm
        end if
      case default
        print *, "BC TYPE NOT RECOGNIZED !"
        error stop
      end select
    end if
  end subroutine compute_right_state

  subroutine face_flux_loop(mesh, sol, prim, grad, hess, rhs, sum_lambda, t)
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in)  :: sol, prim
    real(kind=DOUBLE), dimension(5, 3, mesh%n_elems), intent(in) :: grad
    real(kind=DOUBLE), dimension(:, :, :, :), allocatable, intent(in) :: hess
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems),    intent(inout) :: sum_lambda
    real(kind=DOUBLE), intent(in) :: t

    integer(kind=ENTIER) :: iface, il, ir, iv, k, n_fvert, n_qpts, q
    integer(kind=ENTIER) :: iloop, id_sub_face
    real(kind=DOUBLE), dimension(3) :: norm, xface
    real(kind=DOUBLE), dimension(:, :), allocatable :: face_coords, qpts
    real(kind=DOUBLE), dimension(:),    allocatable :: qwts
    real(kind=DOUBLE), dimension(5) :: wL, wR, flux, sol_l, sol_r
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
          !wR = wL
          !wR = ghost_prim(xface, norm, wL, -ir, t)
          id_sub_face = mesh%face(iface)%sub_face(1)
          sol_l = primit_to_conserv(wL)
          call compute_right_state(mesh, id_sub_face, ir, sol_l, sol_r)
          wR = conserv_to_primit(sol_r)
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

  ! Second geometric moment of each cell about its own centroid, so that
  ! reconstruct()'s order-3 Taylor polynomial has the correct cell AVERAGE,
  ! not just the correct centroid value (same as euler_ho's version).
  subroutine compute_cell_moments(mesh)
    type(mesh_type), intent(in) :: mesh

    integer(kind=ENTIER) :: i, k, n_v, jj, kk
    real(kind=DOUBLE), dimension(:, :), allocatable :: vcoords, pts
    real(kind=DOUBLE), dimension(:), allocatable :: wts
    real(kind=DOUBLE), dimension(3) :: xc

    if (allocated(cell_moment)) deallocate(cell_moment)
    allocate(cell_moment(3, 3, mesh%n_elems))
    cell_moment = 0.0_DOUBLE

    do i = 1, mesh%n_elems
      n_v = mesh%elem(i)%n_vert
      if (n_v /= 4 .and. n_v /= 5 .and. n_v /= 6 .and. n_v /= 8) cycle
      allocate(vcoords(3, n_v))
      do k = 1, n_v
        vcoords(:, k) = mesh%vert(mesh%elem(i)%vert(k))%coord
      end do
      call volume_quad_pts(n_v, vcoords, 2_ENTIER, pts, wts)
      deallocate(vcoords)

      xc = mesh%elem(i)%coord
      do jj = 1, 3
        do kk = 1, 3
          cell_moment(jj, kk, i) = sum(wts * (pts(jj, :) - xc(jj)) * (pts(kk, :) - xc(kk))) &
            / max(mesh%elem(i)%volume, 1.0e-300_DOUBLE)
        end do
      end do
      deallocate(pts, wts)
    end do
  end subroutine compute_cell_moments

  ! Gradient (and Hessian at order >= 3) of the 5 primitives from
  ! arbitrary_high_order_module, same settings as euler_ho's
  ! aho_reconstruction: Green-Gauss nodal derivatives blended to the cells.
  ! grad is exchanged across MPI ranks before being differentiated again,
  ! since a ghost cell's vertex stencil can be incomplete.
  subroutine aho_reconstruction(mesh, prim, grad, hess, num_procs, mpi_send_recv)
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: prim
    real(kind=DOUBLE), dimension(5, 3, mesh%n_elems), intent(inout) :: grad
    real(kind=DOUBLE), dimension(:, :, :, :), allocatable, intent(inout) :: hess
    integer, intent(in) :: num_procs
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    real(kind=DOUBLE), dimension(:, :), allocatable :: grad_flat, hess_flat
    integer(kind=ENTIER) :: e, v, dir1, dir2

    aho_use_green_gauss = .true.

    ! Flat layout of the aho module: component (dir-1)*nc_in + v
    allocate(grad_flat(15, mesh%n_elems))
    call compute_next_order_derivative(mesh, 3_ENTIER, 5_ENTIER, boundary_2d, &
      prim, grad_flat, deriv_order=1_ENTIER)
    if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 15_ENTIER, grad_flat)

    do e = 1, mesh%n_elems
      do dir1 = 1, 3
        do v = 1, 5
          grad(v, dir1, e) = grad_flat((dir1-1)*5 + v, e)
        end do
      end do
    end do

    if (order >= 3 .and. allocated(hess)) then
      allocate(hess_flat(45, mesh%n_elems))
      call compute_next_order_derivative(mesh, 3_ENTIER, 15_ENTIER, boundary_2d, &
        grad_flat, hess_flat, deriv_order=2_ENTIER)
      if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 45_ENTIER, hess_flat)

      do e = 1, mesh%n_elems
        do dir2 = 1, 3
          do dir1 = 1, 3
            do v = 1, 5
              hess(v, dir1, dir2, e) = hess_flat((dir2-1)*15 + (dir1-1)*5 + v, e)
            end do
          end do
        end do
      end do
      deallocate(hess_flat)
    end if

    deallocate(grad_flat)
  end subroutine aho_reconstruction

end module silvia_base_module