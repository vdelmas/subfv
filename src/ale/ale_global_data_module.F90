module ale_global_data_module
  use precision_module
  implicit none

  integer, parameter :: mnbc = 10
  real(kind=DOUBLE), parameter :: gamma = 1.4_DOUBLE

  !Mesh
  character(len=255) :: meshfile_path = "", meshfile
  logical :: boundary_2d = .FALSE.

  !Boundary conditions: 'wall' (default, mirror state), 'freestream' (fixed
  !bc_val state), and 'piston' (impermeable like a wall for the flux, but
  !its nodes are pinned to move at bc_val(2:4,i) -- see ale_main.F90).
  integer(kind=ENTIER) :: n_bc = 0
  character(len=255), dimension(mnbc) :: bc_name = ""
  character(len=255), dimension(mnbc) :: bc_type = ""
  real(kind=DOUBLE), dimension(5, mnbc) :: bc_val = 0.0_DOUBLE
  logical, dimension(:), allocatable :: bc_is_wall
  logical, dimension(:), allocatable :: bc_is_piston

  !Time
  real(kind=DOUBLE) :: cfl = 0.5_DOUBLE
  real(kind=DOUBLE) :: tmax = 1.0_DOUBLE
  real(kind=DOUBLE) :: t, dt

  !Init: uniform state, a 1D Riemann problem split at x=x1drp, or the
  !Haas-Sturtevant shock-bubble (Ms=1.22 in air hitting a He+28%air bubble,
  !same physical setup as test/lagrange_shock_bubble, init=9).
  logical :: init_uniform = .FALSE.
  real(kind=DOUBLE), dimension(5) :: sol_uniform = 0.0_DOUBLE
  logical :: init_1drp = .FALSE.
  real(kind=DOUBLE) :: x1drp = 0.5_DOUBLE
  real(kind=DOUBLE), dimension(5) :: sol_w_1drp_l = 0.0_DOUBLE
  real(kind=DOUBLE), dimension(5) :: sol_w_1drp_r = 0.0_DOUBLE
  logical :: init_shock_bubble = .FALSE.
  real(kind=DOUBLE) :: xc_bub = 0.350_DOUBLE
  real(kind=DOUBLE) :: yc_bub = 0.0445_DOUBLE
  real(kind=DOUBLE) :: r_bub = 0.025_DOUBLE
  real(kind=DOUBLE) :: gamma_air = 1.4_DOUBLE
  real(kind=DOUBLE) :: gamma_bub = 1.648_DOUBLE
  real(kind=DOUBLE) :: rho_bub = 0.287_DOUBLE / 1.578_DOUBLE ! He+28%air, Table 1

  !Grid velocity mode:
  !  'zero'       -> w_p=0            (Eulerian limit)
  !  'lagrangian' -> w_p=v_p          (Lagrangian limit, every node)
  !  'rbf_hybrid' -> w_p=v_p on the bubble interface and piston nodes only,
  !                  RBF-interpolated elsewhere (domain boundary nodes
  !                  slide tangentially along their wall normal, interior
  !                  nodes fully free) -- see ale_main.F90.
  character(len=32) :: ale_grid_velocity = 'zero'
  real(kind=DOUBLE) :: rbf_radius = 0.1_DOUBLE
  ! Only every bubble_iface_stride-th mesh vertex found on the material
  ! interface is kept as a "hard" RBF landmark (w_p=v_p exactly); the
  ! others are left free (w_p interpolated like any interior vertex).
  ! Needed because a compact-support RBF radius large enough to span the
  ! whole domain makes the interpolation matrix near-rank-deficient when
  ! fed *every* interface vertex at the mesh's local (much finer)
  ! resolution -- the CG solve then stalls well short of convergence.
  integer(kind=ENTIER) :: bubble_iface_stride = 10
  ! Same idea for the generic domain-wall sliding points (not the piston,
  ! whose imposed value is constant regardless of how densely it's
  ! sampled -- see ale_main.F90): a finer mesh means more wall vertices at
  ! the *same* domain-spanning radius, and the RBF matrix conditioning
  ! degrades sharply (not gradually) past some density, turning a
  ! CG solve that converges in tens of iterations into one that exhausts
  ! its iteration budget every single timestep.
  integer(kind=ENTIER) :: wall_stride = 10

  !Output
  integer(kind=ENTIER) :: n_iter_print = 50
  integer(kind=ENTIER) :: n_iter_write_sol = 0

contains
  subroutine read_input_parameters(filename)
    implicit none

    character(len=*), intent(in) :: filename

    integer(kind=ENTIER) :: funit

    namelist /INPUT_PARAM/ &
      meshfile_path, meshfile, boundary_2d, &
      n_bc, bc_name, bc_type, bc_val, &
      cfl, tmax, &
      init_uniform, sol_uniform, &
      init_1drp, x1drp, sol_w_1drp_l, sol_w_1drp_r, &
      init_shock_bubble, xc_bub, yc_bub, r_bub, gamma_air, gamma_bub, &
      n_iter_print, n_iter_write_sol

    namelist /ALE_PARAM/ ale_grid_velocity, rbf_radius, bubble_iface_stride, wall_stride

    open(newunit=funit, file=trim(adjustl(filename)))
    read(unit=funit, nml=INPUT_PARAM)
    close(funit)

    open(newunit=funit, file=trim(adjustl(filename)))
    read(unit=funit, nml=ALE_PARAM)
    close(funit)
  end subroutine read_input_parameters

  subroutine init_bc_flags()
    implicit none

    integer(kind=ENTIER) :: i
    character(len=255) :: t

    allocate(bc_is_wall(n_bc), bc_is_piston(n_bc))
    bc_is_piston = .FALSE.
    do i = 1, n_bc
      t = trim(adjustl(bc_type(i)))
      if (t == "wall" .or. t == "") then
        bc_is_wall(i) = .TRUE.
      else if (t == "freestream") then
        bc_is_wall(i) = .FALSE.
      else if (t == "piston") then
        bc_is_wall(i) = .TRUE.
        bc_is_piston(i) = .TRUE.
      else
        print*, "ERROR: bc_type unrecognized for BC ", i, ": '", trim(t), &
          "' (only 'wall'/'freestream'/'piston' supported)"
        error stop
      end if
    end do
  end subroutine init_bc_flags
end module ale_global_data_module
