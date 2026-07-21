module ns_global_data_module
  use precision_module
  implicit none

  integer, parameter :: mnbc = 10 !Maximum number of boundaries
  real(kind=DOUBLE), parameter :: PI = 4.0_DOUBLE*datan(1.0_DOUBLE)

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
    else if (t == "multi_point_pressure") then
      scheme_id = SCHEME_MULTI_POINT_PRESSURE
    else if (t == "multi_point_pressure_ph") then
      scheme_id = SCHEME_MULTI_POINT_PRESSURE_PH
    else if (t == "WIP") then
      scheme_id = SCHEME_WIP
    else if (t(1:2) == "ZB") then
      scheme_id = SCHEME_ZB
      ! Parse ZB_adv_lag into sub-scheme strings
      p1 = index(t, '_')
      p2 = index(t(p1+1:), '_') + p1
      t_adv = trim(t(p1+1:p2-1))
      t_lag = trim(t(p2+1:))
      recognized = .false.
      if (t_adv == "AR1D")     then; scheme_adv_id = SCHEME_ADV_AR1D;     recognized = .true.; end if
      if (t_adv == "AM")       then; scheme_adv_id = SCHEME_ADV_AM;       recognized = .true.; end if
      if (t_adv == "AMISO")    then; scheme_adv_id = SCHEME_ADV_AMISO;    recognized = .true.; end if
      if (t_adv == "ARMD")     then; scheme_adv_id = SCHEME_ADV_ARMD;     recognized = .true.; end if
      if (t_adv == "ARMDWIP")  then; scheme_adv_id = SCHEME_ADV_ARMDWIP;  recognized = .true.; end if
      if (t_adv == "ARMDU")    then; scheme_adv_id = SCHEME_ADV_ARMDU;    recognized = .true.; end if
      if (t_adv == "ARMDM")    then; scheme_adv_id = SCHEME_ADV_ARMDM;    recognized = .true.; end if
      if (t_adv == "ARMDMAT")  then; scheme_adv_id = SCHEME_ADV_ARMDMAT;  recognized = .true.; end if
      if (t_adv == "ARMDUMAT") then; scheme_adv_id = SCHEME_ADV_ARMDUMAT; recognized = .true.; end if
      if (t_adv == "ARMDMMAT") then; scheme_adv_id = SCHEME_ADV_ARMDMMAT; recognized = .true.; end if
      if (.not. recognized) then
        print *, "ERROR: ZB advection unrecognized: '", trim(t_adv), "'"
        error stop
      end if
      recognized = .false.
      if (t_lag == "LS")  then; scheme_lag_id = SCHEME_LAG_LS;    recognized = .true.; end if
      if (t_lag == "LSWIP")  then; scheme_lag_id = SCHEME_LAG_LSWIP; recognized = .true.; end if
      if (t_lag == "LPP")  then; scheme_lag_id = SCHEME_LAG_LPP; recognized = .true.; end if
      if (t_lag == "LVPPP")  then; scheme_lag_id = SCHEME_LAG_LVPPP; recognized = .true.; end if
      if (t_lag == "LSU") then; scheme_lag_id = SCHEME_LAG_LSU;   recognized = .true.; end if
      if (t_lag == "LSM") then; scheme_lag_id = SCHEME_LAG_LSM;   recognized = .true.; end if
      if (t_lag == "LS1D") then; scheme_lag_id = SCHEME_LAG_LS1D; recognized = .true.; end if
      if (t_lag == "LPF")  then; scheme_lag_id = SCHEME_LAG_LPF;  recognized = .true.; end if
      if (.not. recognized) then
        print *, "ERROR: ZB Lagrange unrecognized: '", trim(t_lag), "'"
        error stop
      end if
    else
      print *, "ERROR: scheme unrecognized: '", trim(t), "'"
      error stop
    end if
  end subroutine init_flags
end module ns_global_data_module