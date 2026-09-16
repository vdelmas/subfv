module silvia_module
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
  ! public :: compute_error_vortex
  ! public :: test_reconstruction_exactness
  ! public :: compute_cell_moments
  ! Exposed for the standalone quadrature+reconstruction consistency test
  ! (scratchpad/quad_recon_test) -- not otherwise called outside this module.
  !public :: reconstruct, ls_reconstruction, aho_reconstruction
  public :: primit_to_conserv
  public :: conserv_to_primit
  public :: print_bcs

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
      error_2d, error_2d_h, &
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

  ! subroutine compute_error_isentropic(mesh, sol, t)
  !   use ns_global_data_module, only: error_2d, error_2d_h
  !   use ns_euler_primitives_module, only: conserv_to_primit
  !   implicit none

  !   type(mesh_type), intent(in) :: mesh
  !   real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
  !   real(kind=DOUBLE), intent(in) :: t

  !   integer(kind=ENTIER) :: i
  !   real(kind=DOUBLE) :: error, volume
  !   real(kind=DOUBLE), dimension(5) :: wexact, wsol

  !   error = 0.0_DOUBLE
  !   volume = 0.0_DOUBLE
  !   do i = 1, mesh%n_elems
  !     if (abs(mesh%elem(i)%coord(1)) < 3.0_DOUBLE &
  !       .and. abs(mesh%elem(i)%coord(2)) < 3.0_DOUBLE &
  !       .and. abs(mesh%elem(i)%coord(3)) < 3.0_DOUBLE) then
  !       call sol_isentropic_vortex(mesh%elem(i)%coord, wexact, t)
  !       wsol = conserv_to_primit(sol(:, i))
  !       error = error + mesh%elem(i)%volume*(wsol(1) - wexact(1))**2
  !       volume = volume + mesh%elem(i)%volume
  !     end if
  !   end do

  !   if (error_2d) then
  !     print *, "Error Vortex: ", sqrt((volume/error_2d_h)/mesh%n_elems), sqrt(error)
  !   else
  !     print *, "Error Vortex: ", (volume/mesh%n_elems)**(1.0_DOUBLE/3.0_DOUBLE), sqrt(error)
  !   end if
  ! end subroutine compute_error_isentropic

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
        call sol_gresho_mach(mesh%elem(i)%coord, w, 1.0e-5_DOUBLE)
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

    rhs        = 0.0_DOUBLE
    sum_lambda = 0.0_DOUBLE
    call face_flux_loop(mesh, sol, prim, rhs, sum_lambda, t)

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

  subroutine face_flux_loop(mesh, sol, prim, rhs, sum_lambda, t)
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in)  :: sol, prim
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems),    intent(inout) :: sum_lambda
    real(kind=DOUBLE), intent(in) :: t

    integer(kind=ENTIER) :: iface, il, ir, iv, k, n_fvert, n_qpts, q
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

      ! Find z-faces
      is_zface = boundary_2d .and. abs(abs(norm(3)) - 1.0_DOUBLE) < 1.0e-6_DOUBLE
      n_qpts = 1
      allocate(qpts(3, 1), qwts(1))
      qpts(:, 1) = mesh%face(iface)%coord
      qwts(1)    = mesh%face(iface)%area

      ! For z-faces in 2D mode: zero net physics, skip entirely
      if (is_zface) then
        deallocate(qpts, qwts)
        cycle
      end if

      n_qpts = size(qwts)
      do q = 1, n_qpts
        xface = qpts(:, q)

        wL = prim(:, il)
        if (ir > 0) then
          wR = prim(:, ir)
        else
          ! Boundary: reconstruct left then apply BC
          wR = wL
          ! if (mesh%elem(ir)%is_ghost) then
          !   print*, 'found ghost!'
          ! end if
          !wR = ghost_prim(xface, norm, wL, -ir, t)
        end if

        ! Numerical flux (physical: per-area * area_weight_at_qpt).
        ! flux_scheme_id resolved once in read_params -- see its
        ! declaration for why this is an integer compare, not a string one.
        ! select case (flux_scheme_id)
        ! case (FLUX_THREE_WAVE)
        !   flux = three_wave_flux(wL, wR, norm) * qwts(q)
        ! case (FLUX_TWO_WAVE)
        flux = two_wave_flux(wL, wR, norm) * qwts(q)
        ! case default
        !   flux = rusanov(wL, wR, norm) * qwts(q)
        ! end select
        lambda = max(abs(dot_product(wL(2:4), norm)) + cs(wL), &
                     abs(dot_product(wR(2:4), norm)) + cs(wR)) * qwts(q)

        rhs(:, il)  = rhs(:, il)  - flux
        sum_lambda(il) = sum_lambda(il) + lambda
        !if (ir > 0 .and. .not. mesh%elem(ir)%is_ghost) then
          rhs(:, ir)     = rhs(:, ir)     + flux
          sum_lambda(ir) = sum_lambda(ir) + lambda
        !end if
      end do

      deallocate(qpts, qwts)
    end do
  end subroutine face_flux_loop

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
end module silvia_module