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
  integer, parameter :: SCHEME_USI3D                   = 9  ! "USI3D"
  integer, parameter :: SCHEME_MULTI_POINT_VILAR       = 10 ! "multi_point_vilar"
  integer, parameter :: SCHEME_WIP2                    = 11 ! "WIP2"
  integer, parameter :: SCHEME_WIP2_NOLM               = 12 ! "WIP2_NOLM" (theta=1, reference)
  integer, parameter :: SCHEME_WIP2_AADV               = 13 ! "WIP2_AADV" (theta=1, AMISO advection)
  integer, parameter :: SCHEME_WIP2_HYB                = 14 ! "WIP2_HYB" (theta=1, shock-sensor blended advection)
  ! eps_p sensor variants of WIP2_NOLM (see tex/wip.tex's four candidates)
  integer, parameter :: SCHEME_WIP2_EJUMP              = 15 ! "WIP2_EJUMP" (d) normalised pressure jump
  integer, parameter :: SCHEME_WIP2_EDIV               = 16 ! "WIP2_EDIV"  (a) -div(v) only
  integer, parameter :: SCHEME_WIP2_EMAX               = 17 ! "WIP2_EMAX"  (b) max(p jump, -div v), Ma^2 gate
  integer, parameter :: SCHEME_WIP2_ENTH                = 18 ! "WIP2_ENTH" Haenel/MGallice enthalpy fix
  integer, parameter :: SCHEME_WIP2_ENTH_EJUMP          = 19 ! "WIP2_ENTH_EJUMP" enthalpy fix + (d) sensor
  integer, parameter :: SCHEME_WIP2_TP                  = 20 ! "WIP2_TP" Tallois multi-D low-Mach correction
  integer, parameter :: SCHEME_WIP2_TP_ENTH             = 21 ! "WIP2_TP_ENTH" + Haenel enthalpy fix
  integer, parameter :: SCHEME_WIP2_EPOS                = 22 ! "WIP2_EPOS" per-face positivity-driven eps, no sensor
  integer, parameter :: SCHEME_WIP2_EPOS2               = 23 ! "WIP2_EPOS2" same, impedance-scaled pressure term
  integer, parameter :: SCHEME_WIP2_ED1                 = 24 ! "WIP2_ED1" Ducros sensor at amplitude 1*a_p
  integer, parameter :: SCHEME_WIP2_ED2                 = 25 ! "WIP2_ED2" Ducros sensor at amplitude 2*a_p
  integer, parameter :: SCHEME_WIP2_MD                  = 26 ! "WIP2_MD" multi-d nodal advection, AMISO hybridisation
  integer, parameter :: SCHEME_WIP2_MD5                 = 27 ! "WIP2_MD5" same, ARMD hybridisation w=0.5
  integer, parameter :: SCHEME_WIP2_MDH                 = 28 ! "WIP2_MDH" multi-d advection, enthalpy-consistent gradient
  integer, parameter :: SCHEME_WIP2_PPEP                = 29 ! "WIP2_PPEP" Eulerian port of Lagrangian vp_pp_ep
  integer, parameter :: SCHEME_WIP2_PPEPI               = 30 ! "WIP2_PPEPI" same, nodal INTERNAL energy in the relaxation
  integer, parameter :: SCHEME_WIP2_PPEPV               = 31 ! "WIP2_PPEPV" mu = rho |(v - v*).n|
  integer, parameter :: SCHEME_WIP2_PPEPVI              = 32 ! "WIP2_PPEPVI" same with internal energy
  ! 60-62: kept clear of the 1-21 block, which the WIP2 family keeps growing into. These three
  ! were 15/16/17 and silently collided with SCHEME_WIP2_EJUMP/EDIV/EMAX: the WIP2 if-chain in
  ! ns_euler_module fires before the flux select case, so the enthalpy schemes ran a WIP2
  ! solver instead, with no error message. Renumber here, not there, if more are added.
  integer, parameter :: SCHEME_THREE_WAVE_ENTHALPY     = 60 ! "three_wave_enthalpy" (enthalpy-preserving three_wave)
  integer, parameter :: SCHEME_MULTI_POINT_ENTHALPY    = 61 ! "multi_point_enthalpy" (enthalpy-preserving multi_point)
  integer, parameter :: SCHEME_THREE_WAVE_ENTHALPY2    = 62 ! "three_wave_enthalpy2" (HLL-consistent realisation; see solver header)
  ! multi_point form of the Lagrangian vp_pp_ep (nodal pressure + nodal-energy energy flux)
  integer, parameter :: SCHEME_MULTI_POINT_ENER        = 63 ! "multi_point_ener"       mu = lambda
  integer, parameter :: SCHEME_MULTI_POINT_ENER_ISO    = 64 ! "multi_point_ener_iso"   + three_wave blend
  integer, parameter :: SCHEME_MULTI_POINT_ENER_V      = 65 ! "multi_point_ener_v"     mu = rho |v - v*|
  integer, parameter :: SCHEME_MULTI_POINT_ENER_V_ISO  = 66 ! "multi_point_ener_v_iso" + three_wave blend
  ! nodal INTERNAL energy eps_p: eps*_l = eps*_r = eps_p on every sub-face (mu = lambda)
  integer, parameter :: SCHEME_MULTI_POINT_EPS         = 67 ! "multi_point_eps"
  integer, parameter :: SCHEME_MULTI_POINT_EPS_ISO     = 68 ! "multi_point_eps_iso"    + three_wave blend
  ! two DIFFERENT star internal energies, no internal-energy exchange between the sides; one nodal
  ! scalar kappa restores conservation: eps*_s = eps_s + kappa (eps2) or eps_s (1 + kappa) (eps2m)
  integer, parameter :: SCHEME_MULTI_POINT_EPS2        = 69 ! "multi_point_eps2"
  integer, parameter :: SCHEME_MULTI_POINT_EPS2M       = 70 ! "multi_point_eps2m"
  ! (pv)*_s = p_p u*_s, u*_s the acoustic star velocity of side s: no nodal energy at all; nodal
  ! energy conservation comes from p_p being the 1/lambda mean of the pbar (interior nodes)
  integer, parameter :: SCHEME_MULTI_POINT_PS          = 71 ! "multi_point_ps"
  ! multi_point_ps with p_p averaged over the INTERIOR sub-faces only and the classic multi_point
  ! flux on boundary and boundary_2d z-face sub-faces (one-sided pbar_l, zero energy flux at a wall
  ! or symmetry plane): the interior telescoping then closes exactly at boundary nodes too
  integer, parameter :: SCHEME_MULTI_POINT_PSW         = 72 ! "multi_point_psw"
  ! nodal p_p in the MOMENTUM flux, the classic Pi_s (v*.n) in the ENERGY flux (Vincent): conservative
  ! because the nodal velocity solve is sum A (Pi_l - Pi_r) n = 0; boundary / z-faces as psw
  integer, parameter :: SCHEME_MULTI_POINT_PV          = 73 ! "multi_point_pv"
  integer, parameter :: SCHEME_ZB                   = 42 ! "ZB_*_*"
  ! ZB advection sub-scheme IDs
  integer, parameter :: SCHEME_ADV_AR1D     = 1  ! "AR1D"
  integer, parameter :: SCHEME_ADV_AM       = 2  ! "AM"
  integer, parameter :: SCHEME_ADV_AMISO    = 3  ! "AMISO"
  integer, parameter :: SCHEME_ADV_ARMD     = 4  ! "ARMD"
  integer, parameter :: SCHEME_ADV_ARMDU    = 5  ! "ARMDU"
  integer, parameter :: SCHEME_ADV_ARMDMAT  = 7  ! "ARMDMAT"
  integer, parameter :: SCHEME_ADV_ARMDUMAT = 8  ! "ARMDUMAT"
  ! SCHEME_ADV_ARMDM (6), SCHEME_ADV_ARMDMMAT (9) removed 2026-09-15: only
  ! active users of compute_ellip's mat_h_p output, and unreachable via any
  ! schemes.txt entry (advection sub-schemes there are AMISO/AR1D only)
  integer, parameter :: SCHEME_ADV_ARMDWIP     = 10  ! "ARMDWIP"
  ! ZB Lagrange sub-scheme IDs
  integer, parameter :: SCHEME_LAG_LS  = 1  ! "LS"
  integer, parameter :: SCHEME_LAG_LSU = 2  ! "LSU"
  ! SCHEME_LAG_LSM (3) removed 2026-09-15: same reason, unreachable via
  ! schemes.txt (lag sub-schemes there are LVPPP/LPP/LS1D/LPF only)
  integer, parameter :: SCHEME_LAG_LSWIP = 4  ! "LSWIP"
  integer, parameter :: SCHEME_LAG_LPP = 5  ! "LPP"
  integer, parameter :: SCHEME_LAG_LVPPP = 6  ! "LVPPP"
  integer, parameter :: SCHEME_LAG_LS1D = 7  ! "LS1D"
  integer, parameter :: SCHEME_LAG_LPF  = 8  ! "LPF"
  ! Runtime scheme selection (set by init_bc_flags)
  integer :: scheme_id = 0
  ! Pilotage par namelist du balayage ZB/WIP2. -1 = laisser le scheme_id decider.
  ! zb_lag_mode : forme du flux d energie de la partie Lagrange (0 a 3).
  ! zb_adv_mode : flux d advection (0 a 5). Les deux se croisent librement.
  integer :: zb_lag_mode = -1
  integer :: zb_adv_mode = -1
  ! zb_eps_mode : capteur de dissipation d advection (0 a 8). -1 = scheme_id.
  ! zb_eps_coef : amplitude du mode 8, en fraction de la vitesse du son nodale.
  !   Le mode 8 n a AUCUN capteur : la stabilisation de l advection est la
  !   VITESSE seule, |vstar|, plus eps = zb_eps_coef * a_p. zb_eps_coef = 0
  !   donne le decentrement pur a vstar, sans aucune dissipation ajoutee.
  integer :: zb_eps_mode = -1
  real(kind=DOUBLE) :: zb_eps_coef = 0.0_DOUBLE
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
  ! multi_point_ener*: raise lambda_s to the root lambda_req,s that makes the star internal
  ! energy non-negative (fixed point on the node, at most 5 passes, lambdas only increase)
  logical :: ener_lambda_pos = .FALSE.
  ! diagnostic: at the first flux call, write every interior sub-face side of the
  ! multi_point_ener* node passes to ener_sides_<rank>.txt (no effect on the solution)
  logical :: ener_dump_sides = .FALSE.
  ! CFL number of the convex-combination argument in LOCAL time stepping: the point-implicit
  ! update divides by V/dt + sum_lambda/ener_theta, so the weight kept on the cell's own
  ! state is alpha_c = 1 - ener_theta (cfl = 1e5 makes V/dt negligible). 1 = historical
  ! behaviour (alpha_c = 0).
  real(kind=DOUBLE) :: ener_theta = 1.0_DOUBLE
  ! multi_point_ener*: relax the specific INTERNAL energy eps = E - |u|^2/2 instead of the total
  ! energy E in the one-sided energy flux (pv)_s = p_s u_s -/+ mu_s (X_p - X_s) and in X_p.
  ! Conservation is unchanged (it comes from the telescoping of the jumps, whatever X is);
  ! relaxing E removes total energy from a fast cell without removing momentum, i.e. drains
  ! its internal energy.
  logical :: ener_internal = .FALSE.
  ! multi_point_ener*: pressure work at the STAR velocity, (pv)*_l = p_p u*_l - mu_l (X_p - X_l),
  ! (pv)*_r = p_p u*_r + mu_r (X_p - X_r), u*_l = u_l + (p_l - p_p)/lambda_l,
  ! u*_r = u_r + (p_p - p_r)/lambda_r (nodal lambdas). Makes eps* convex in p_l - p_p again
  ! (as the classic pbar v*); X_p from the same telescoping balance.
  logical :: ener_pstar = .FALSE.

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
      zb_lag_mode, zb_adv_mode, zb_eps_mode, zb_eps_coef, &
      local_time_step, ener_lambda_pos, ener_dump_sides, ener_theta, ener_internal, ener_pstar, &
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
    else if (t == "three_wave_enthalpy") then
      scheme_id = SCHEME_THREE_WAVE_ENTHALPY
    else if (t == "multi_point_enthalpy") then
      scheme_id = SCHEME_MULTI_POINT_ENTHALPY
    else if (t == "three_wave_enthalpy2") then
      scheme_id = SCHEME_THREE_WAVE_ENTHALPY2
    else if (t == "multi_point_ener") then
      scheme_id = SCHEME_MULTI_POINT_ENER
    else if (t == "multi_point_ener_iso") then
      scheme_id = SCHEME_MULTI_POINT_ENER_ISO
    else if (t == "multi_point_ener_v") then
      scheme_id = SCHEME_MULTI_POINT_ENER_V
    else if (t == "multi_point_ener_v_iso") then
      scheme_id = SCHEME_MULTI_POINT_ENER_V_ISO
    else if (t == "multi_point_eps") then
      scheme_id = SCHEME_MULTI_POINT_EPS
    else if (t == "multi_point_eps_iso") then
      scheme_id = SCHEME_MULTI_POINT_EPS_ISO
    else if (t == "multi_point_eps2") then
      scheme_id = SCHEME_MULTI_POINT_EPS2
    else if (t == "multi_point_eps2m") then
      scheme_id = SCHEME_MULTI_POINT_EPS2M
    else if (t == "multi_point_ps") then
      scheme_id = SCHEME_MULTI_POINT_PS
    else if (t == "multi_point_psw") then
      scheme_id = SCHEME_MULTI_POINT_PSW
    else if (t == "multi_point_pv") then
      scheme_id = SCHEME_MULTI_POINT_PV
    else if (t == "multi_point_pressure") then
      scheme_id = SCHEME_MULTI_POINT_PRESSURE
    else if (t == "multi_point_pressure_ph") then
      scheme_id = SCHEME_MULTI_POINT_PRESSURE_PH
    else if (t == "multi_point_vilar") then
      scheme_id = SCHEME_MULTI_POINT_VILAR
    else if (t == "WIP") then
      scheme_id = SCHEME_WIP
    else if (t == "WIP2") then
      scheme_id = SCHEME_WIP2
    else if (t == "WIP2_NOLM") then
      scheme_id = SCHEME_WIP2_NOLM
    else if (t == "WIP2_AADV") then
      scheme_id = SCHEME_WIP2_AADV
    else if (t == "WIP2_HYB") then
      scheme_id = SCHEME_WIP2_HYB
    else if (t == "WIP2_EJUMP") then
      scheme_id = SCHEME_WIP2_EJUMP
    else if (t == "WIP2_EDIV") then
      scheme_id = SCHEME_WIP2_EDIV
    else if (t == "WIP2_EMAX") then
      scheme_id = SCHEME_WIP2_EMAX
    else if (t == "WIP2_ENTH") then
      scheme_id = SCHEME_WIP2_ENTH
    else if (t == "WIP2_ENTH_EJUMP") then
      scheme_id = SCHEME_WIP2_ENTH_EJUMP
    else if (t == "WIP2_TP") then
      scheme_id = SCHEME_WIP2_TP
    else if (t == "WIP2_TP_ENTH") then
      scheme_id = SCHEME_WIP2_TP_ENTH
    else if (t == "WIP2_EPOS") then
      scheme_id = SCHEME_WIP2_EPOS
    else if (t == "WIP2_EPOS2") then
      scheme_id = SCHEME_WIP2_EPOS2
    else if (t == "WIP2_ED1") then
      scheme_id = SCHEME_WIP2_ED1
    else if (t == "WIP2_ED2") then
      scheme_id = SCHEME_WIP2_ED2
    else if (t == "WIP2_MD") then
      scheme_id = SCHEME_WIP2_MD
    else if (t == "WIP2_MD5") then
      scheme_id = SCHEME_WIP2_MD5
    else if (t == "WIP2_MDH") then
      scheme_id = SCHEME_WIP2_MDH
    else if (t == "WIP2_PPEP") then
      scheme_id = SCHEME_WIP2_PPEP
    else if (t == "WIP2_PPEPI") then
      scheme_id = SCHEME_WIP2_PPEPI
    else if (t == "WIP2_PPEPV") then
      scheme_id = SCHEME_WIP2_PPEPV
    else if (t == "WIP2_PPEPVI") then
      scheme_id = SCHEME_WIP2_PPEPVI
    else if (t == "USI3D") then
      scheme_id = SCHEME_USI3D
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
      if (t_adv == "ARMDMAT")  then; scheme_adv_id = SCHEME_ADV_ARMDMAT;  recognized = .true.; end if
      if (t_adv == "ARMDUMAT") then; scheme_adv_id = SCHEME_ADV_ARMDUMAT; recognized = .true.; end if
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