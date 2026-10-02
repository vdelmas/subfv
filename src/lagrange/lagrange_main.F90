program main
  use mpi
  use precision_module
  use mpi_module
  use mesh_module
  use mesh_reading_module
  use mesh_geometry_module
  use mesh_connectivity_module
  use lagrange_module
  use lagrange_io_module
  implicit none

  integer(kind=ENTIER) :: fn
  integer :: mpi_ierr, me, num_procs
  type(mpi_send_recv_type) :: mpi_send_recv
  character(len=255) :: meshfile, meshfile_path

  character(len=255) :: fln

  integer(kind=ENTIER), parameter :: n_max_bc = 10
  integer(kind=ENTIER) :: n_bc, iter, i, j, method_length = 0
  integer(kind=ENTIER) :: id_sub_face, id_face, piston_bc_idx
  real(kind=DOUBLE) :: cfl, cfl_max = 0.8, b2d_h = 1.0
  real(kind=DOUBLE) :: h_extrude = 1.0
  integer(kind=ENTIER) :: init = 0
  ! Diagnostic: drop boundary vertices from the o2 nodal gradient fit, as the Euler
  ! code does, instead of mirroring. Default .false. = behaviour unchanged.
  logical :: grad_bound_skip = .false.
  character(len=255) :: scheme = ""
  character(len=255), dimension(n_max_bc) :: bc_name
  character(len=255), dimension(n_max_bc) :: bc_type
  real(kind=DOUBLE), dimension(5, n_max_bc) :: bc_val
  real(kind=DOUBLE), dimension(5) :: sol_uniform
  real(kind=DOUBLE), dimension(3) :: piston_vel
  type(mesh_type) :: mesh

  logical :: boundary_2d = .true.

  real(kind=DOUBLE) :: t, t_max, dt
  real(kind=DOUBLE), dimension(:,:), allocatable :: vp
  real(kind=DOUBLE), dimension(:), allocatable :: mass, pp, gamma_arr
  real(kind=DOUBLE), dimension(:,:), allocatable :: sol
  real(kind=DOUBLE), dimension(:,:), allocatable :: rhs
  real(kind=DOUBLE), dimension(:,:), allocatable :: new_sol
  logical, dimension(:), allocatable :: vp_is_imposed
  ! Second-order reconstruction: per-cell grad_v (3x3), grad_p (3), div_v (1)
  logical :: second_order = .false.
  real(kind=DOUBLE), dimension(:,:,:), allocatable :: grad_v
  real(kind=DOUBLE), dimension(:,:), allocatable :: grad_p
  real(kind=DOUBLE), dimension(:), allocatable :: div_v
  real(kind=DOUBLE), dimension(:), allocatable :: alpha_p_arr
  integer(kind=ENTIER), dimension(:), allocatable :: p_limited
  real(kind=DOUBLE), dimension(:), allocatable :: h_p_arr

  integer(kind=ENTIER) :: n_sol_vtu=2
  integer(kind=ENTIER) :: i_sol_vtu
  integer(kind=ENTIER) :: iso_weight_mode = 2   ! w_pcf recipe for classic_iso_vol
  ! Isentropic-vortex (init=2) convergence diagnostics
  logical :: compute_error = .false.
  logical :: grad_weno = .true.
  logical :: grad_gg_corrected = .false.
  real(kind=DOUBLE) :: grad_eps_weno = 0.0_DOUBLE   ! <=0 -> tiny(), the historical value
  real(kind=DOUBLE) :: h_err, l1err, l2err, linferr

  namelist /INPUT_PARAM/ &
    meshfile_path, meshfile, &
    t_max, cfl, &
    init, &
    grad_bound_skip, &
    lambda_acoustic_only, ppvp_jump_mode, ppvp_nodal_average, diag_entropy, vp_ep_wip_mode, vp_ep_wip_floor, vp_ep_wip_clip, &
    sedov_nodal_deposit, sedov_energy, &
    scheme, method_length, b2d_h, &
    n_bc, bc_name, bc_type, bc_val, &
    sol_uniform, boundary_2d, &
    n_sol_vtu, second_order, iso_weight_mode, compute_error, grad_weno, grad_gg_corrected, grad_eps_weno, &
    gresho_mach, gamma_uniform

  call MPI_INIT(mpi_ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

  open(newunit=fn, file="input_data.f")
  read(unit=fn, nml=INPUT_PARAM)
  close(fn)

  ! "classic_rhoa" is the EUCCLHYD-original variant of "classic": the sub-face
  ! swept-mass coefficient uses the pure acoustic impedance lambda = rho*a,
  ! without the pressure-/velocity-jump (Dukowicz-like) terms. Selecting it by
  ! scheme name rather than by the lambda_acoustic_only namelist flag alone
  ! keeps the two variants distinguishable in a campaign's inputs and logs.
  if (scheme == "classic_rhoa") lambda_acoustic_only = .true.
  if (scheme == "vp_ep_wip") vp_ep_wip = .true.
  if (me == 0) then
    if (lambda_acoustic_only) then
      print*, "lambda = rho*a (EUCCLHYD original)"
    else
      print*, "lambda = max(rho*a, pressure-jump, velocity-jump)"
    end if
  end if

  call read_mesh_msh(mesh, meshfile_path, meshfile, &
    n_bc, bc_name, me, num_procs, mpi_send_recv)
  call build_mesh(mesh, num_procs, mpi_send_recv, &
    .true., boundary_2d)

  ! Saltzman (init=5): skewed mesh by mapping X_sk = X + (0.1-Y)*sin(pi*X), Y_sk = Y
  ! (Maire 2007 / Vilar et al.). sin(pi*X)=0 at X=0,1 so left/right boundaries unaffected.
  if (init == 5) then
    block
      real(kind=DOUBLE), parameter :: pi_cst = acos(-1.0_DOUBLE)
      real(kind=DOUBLE) :: xv, yv
      do i = 1, mesh%n_vert
        xv = mesh%vert(i)%coord(1)
        yv = mesh%vert(i)%coord(2)
        mesh%vert(i)%coord(1) = xv + (0.1_DOUBLE - yv) * sin(pi_cst*xv)
      end do
    end block
  end if

  call compute_geometry_mesh(mesh, .true., boundary_2d)

  ! Extrusion height of a boundary_2d mesh (single z-layer): used by the
  ! classic_iso_vol weight to turn a prism sub-element volume back into an
  ! in-plane area scale.  The nodal solver keeps vp_z = 0 in b2d, so this is
  ! constant in time.
  h_extrude = 1.0_DOUBLE
  if (boundary_2d .and. mesh%n_vert > 0) then
    h_extrude = maxval(mesh%vert(1:mesh%n_vert)%coord(3)) &
              - minval(mesh%vert(1:mesh%n_vert)%coord(3))
    if (h_extrude <= 0.0_DOUBLE) h_extrude = 1.0_DOUBLE
  end if

  allocate(vp(3, mesh%n_vert))
  allocate(pp(mesh%n_vert))
  allocate(sol(5, mesh%n_elems))
  allocate(new_sol(5, mesh%n_elems))
  allocate(rhs(5, mesh%n_elems))
  allocate(gamma_arr(mesh%n_elems))
  allocate(vp_is_imposed(mesh%n_vert))

  ! gamma_uniform > 0 overrides the per-init default. Needed for the *classic*
  ! Noh problem, which is posed with gamma = 5/3 (post-shock rho = 16 in 2D
  ! cylindrical, shock at r = t/3) whereas init == 3 defaults to 1.4 here.
  if (gamma_uniform > 0.0_DOUBLE) then
    gamma_arr = gamma_uniform
  else if (init == 5) then
    gamma_arr = 5.0_DOUBLE / 3.0_DOUBLE
  else
    gamma_arr = 1.4_DOUBLE
  end if
  vp = 0.0_DOUBLE
  call init_sol(mesh, sol, sol_uniform, init, me, num_procs, gamma_arr, boundary_2d)
  call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol)

  ! Identify piston nodes from bc_type='piston'
  vp_is_imposed = .false.
  piston_bc_idx = 0
  piston_vel = 0.0_DOUBLE
  do i = 1, n_bc
    if (trim(bc_type(i)) == 'piston') then
      piston_bc_idx = i
      piston_vel = bc_val(2:4, i)
      exit
    end if
  end do
  if (piston_bc_idx > 0) then
    do i = 1, mesh%n_vert
      do j = 1, mesh%vert(i)%n_sub_faces_neigh
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        if (mesh%face(id_face)%right_neigh == -piston_bc_idx) then
          vp_is_imposed(i) = .true.
          exit
        end if
      end do
    end do
  end if

  allocate(grad_v(3, 3, mesh%n_elems))
  allocate(grad_p(3, mesh%n_elems))
  allocate(div_v(mesh%n_elems))
  allocate(alpha_p_arr(mesh%n_vert))
  allocate(p_limited(mesh%n_elems))
  allocate(h_p_arr(mesh%n_vert))
  grad_v = 0.0_DOUBLE
  grad_p = 0.0_DOUBLE
  div_v  = 0.0_DOUBLE
  alpha_p_arr = 1.0_DOUBLE
  p_limited   = 0
  h_p_arr     = 0.0_DOUBLE
  if (scheme == "sidil") then
    do i = 1, mesh%n_vert
      h_p_arr(i) = compute_length(mesh, i, method_length, boundary_2d, h_extrude)
    end do
  end if

  i_sol_vtu = 0
  write(fln, *) i_sol_vtu
  write(fln, *) "output_"//trim(adjustl(fln))
  call write_sol_lag(mesh, fln, sol, vp, pp, gamma_arr, grad_v, grad_p, div_v, alpha_p_arr, p_limited, h_p_arr)
  call write_sol_dat_lag(mesh, fln, sol, gamma_arr)
  i_sol_vtu = i_sol_vtu + 1

  allocate(mass(mesh%n_elems))
  mass = 0.0_DOUBLE
  do i=1, mesh%n_elems
    mass(i) = mesh%elem(i)%volume/sol(1, i)
  end do

  t = 0.0_DOUBLE
  ! dt is read by compute_gradients / compute_rhs_* (GRP half-step) BEFORE
  ! compute_dt assigns it on the first iteration -- it must start defined.
  dt = 0.0_DOUBLE
  iter = 1
  do while ( t < t_max )

    ! Set imposed velocities (piston BC) before the RHS computation
    if (piston_bc_idx > 0) then
      do i = 1, mesh%n_vert
        if (vp_is_imposed(i)) vp(:, i) = piston_vel
      end do
    end if

    if (second_order) then
      call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol)
      call compute_gradients(mesh, sol, gamma_arr, dt, grad_v, grad_p, div_v, p_limited, grad_weno, &
        grad_gg_corrected, grad_eps_weno, bound_skip=grad_bound_skip)
      ! Exchange ghost cell gradients so second-order RHS can reconstruct across MPI boundaries
      call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 3, grad_p)
      call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 9, grad_v)
      call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 1, div_v)
    end if

    if( scheme == "classic" .or. scheme == "classic_rhoa" ) then
      call compute_rhs_lagrange(mesh, sol, vp, &
        dt, rhs, n_bc, bc_type, bc_val, boundary_2d, mass, gamma_arr, vp_is_imposed, &
        second_order, grad_v, grad_p, div_v, alpha_p_arr)
    else if( scheme == "classic_iso" ) then
      call compute_rhs_lagrange_classic_iso(mesh, sol, vp, &
        dt, rhs, n_bc, bc_type, bc_val, boundary_2d, mass, gamma_arr, vp_is_imposed, &
        1, h_extrude)
    else if( scheme == "classic_iso_vol" ) then
      call compute_rhs_lagrange_classic_iso(mesh, sol, vp, &
        dt, rhs, n_bc, bc_type, bc_val, boundary_2d, mass, gamma_arr, vp_is_imposed, &
        iso_weight_mode, h_extrude)
    else  if( scheme == "sidil" ) then
      ! h_extrude (auto-detected from the mesh's own z-extent above), not
      ! the user-set b2d_h namelist value: relying on a manually-set
      ! parameter to match the mesh's actual extrusion thickness is a
      ! standing footgun (easy to forget/mistype per test case) now that
      ! compute_length's b2d branch genuinely depends on it being exact.
      call compute_rhs_lagrange_sidil(mesh, sol, vp, &
        dt, rhs, n_bc, bc_type, bc_val, boundary_2d, mass, method_length, h_extrude, &
        gamma_arr, vp_is_imposed, h_p_arr)
    else if( scheme == "vp_ep" .or. scheme == "vp_ep_wip" ) then
      ! vitesse nodale + energie nodale, puis solveur 1D par sous-face.
      call compute_rhs_lagrange_vp_ep(mesh, sol, vp, &
        dt, rhs, n_bc, bc_type, bc_val, boundary_2d, mass, gamma_arr, vp_is_imposed, &
        second_order, grad_v, grad_p, div_v)
    else if( scheme == "vp_pp" .or. scheme == "nodal_vp" .or. scheme == "vp_pp_ppvp" ) then
      ! vitesse nodale ET pression nodale, toutes deux issues de systemes nodaux :
      ! pas de h_p, donc aucune dependance a method_length / b2d_h.
      ! "nodal_vp" est l'ancien nom du schema, conserve en alias : beaucoup de
      ! input_data.f de runs deja effectues le portent encore.
      call compute_rhs_lagrange_vp_pp(mesh, sol, vp, &
        dt, rhs, n_bc, bc_type, bc_val, boundary_2d, mass, gamma_arr, vp_is_imposed, &
        second_order, grad_v, grad_p, div_v, scheme == "vp_pp_ppvp")
    else
      print*, "No scheme !"
      error stop
    end if

    ! Taylor-Green needs an energy source to close the total energy equation
    ! (Vilar eq. 4.103). Added here, on the PRE-move mesh, because that is the
    ! geometry rhs was just assembled on.
    if (init == 10) call add_taylor_green_source(mesh, mass, rhs)

    call compute_dt(mesh, sol, dt, cfl, vp, me, num_procs, gamma_arr)
    if( t + dt > t_max ) dt = t_max - t
    call move_mesh(mesh, vp, dt)

    call mpi_memory_exchange_vert(mesh, mpi_send_recv)

    call compute_geometry_mesh(mesh, .true., boundary_2d, do_check=.false.)
    new_sol = sol + dt * rhs
    ! Reset specific volume from geometry to enforce geometric consistency.
    ! This prevents tau from drifting away from the actual volume/mass ratio
    ! due to discretization errors, especially important for large deformations.
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) new_sol(1, i) = mesh%elem(i)%volume / mass(i)
    end do
    sol = new_sol
    call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol)
    t = t + dt
    if( mod(iter, 100) == 0 .and. me == 0 ) print*, t, dt
    if( diag_entropy .and. mod(iter, 100) == 0 ) &
      call print_entropy_diag(mesh, t, me, num_procs)

    if( t >= i_sol_vtu * t_max / real(n_sol_vtu - 1) ) then
      write(fln, *) i_sol_vtu
      write(fln, *) "output_"//trim(adjustl(fln))
      call write_sol_lag(mesh, fln, new_sol, vp, pp, gamma_arr, grad_v, grad_p, div_v, alpha_p_arr, p_limited, h_p_arr)
      call write_sol_dat_lag(mesh, fln, new_sol, gamma_arr)
      i_sol_vtu = i_sol_vtu + 1
    end if

    iter = iter + 1
  end do

  i_sol_vtu = -1
  write(fln, *) i_sol_vtu
  write(fln, *) "output_"//trim(adjustl(fln))
  call write_sol_lag(mesh, fln, new_sol, vp, pp, gamma_arr, grad_v, grad_p, div_v, alpha_p_arr, p_limited, h_p_arr)
  call write_sol_dat_lag(mesh, fln, new_sol, gamma_arr)

  if (compute_error) then
    if (init == 10) then
      ! Taylor-Green: the rate is measured on the pressure, rho being uniform here.
      call compute_error_taylor_green(mesh, new_sol, gamma_arr, h_extrude, h_err, &
        l1err, l2err, linferr)
    else
      call compute_error_vortex(mesh, new_sol, t, h_extrude, h_err, l1err, l2err, linferr)
    end if
    if (me == 0) print *, "FINAL t=", t, "h=", h_err, "L1=", l1err, "L2=", l2err, "Linf=", linferr
  end if

  call MPI_FINALIZE(mpi_ierr)
end program main
