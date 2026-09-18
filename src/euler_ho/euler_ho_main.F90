! High-order explicit Euler solver — main program.
!
! Time integration: SSP-RK3 (Shu-Osher 1988)
!   u1       = u0 + dt * L(u0)
!   u2       = 3/4 * u0 + 1/4 * (u1 + dt * L(u1))
!   u^{n+1}  = 1/3 * u0 + 2/3 * (u2 + dt * L(u2))
program euler_ho_main
  use mpi
  use precision_module
  use mpi_module
  use mesh_module
  use mesh_reading_module
  use mesh_geometry_module
  use mesh_connectivity_module
  use io_module
  use euler_ho_module
  implicit none

  integer :: mpi_ierr, me, num_procs
  type(mesh_type) :: mesh
  type(mpi_send_recv_type) :: mpi_send_recv

  ! sol/sol1/sol2/rhs are 6-wide: components 1:5 are the usual
  ! [rho,rho*u,rho*v,rho*w,rho*E], component 6 is the gamma-transport
  ! variable rho*Gamma (Gamma=1/(gamma-1), see euler_ho_module's gamma_arr
  ! declaration) -- it rides through the exact same RK3/MPI-exchange
  ! arithmetic as the 5 physical equations, no separate bookkeeping.
  real(kind=DOUBLE), allocatable :: sol(:, :), sol1(:, :), sol2(:, :), sol3(:, :)
  real(kind=DOUBLE), allocatable :: prim(:, :), rhs(:, :)
  real(kind=DOUBLE), allocatable :: k1(:, :), k2(:, :), k3(:, :), k4(:, :)
  real(kind=DOUBLE), allocatable :: sum_lambda(:)

  real(kind=DOUBLE) :: t, dt, h_err, l2err
  integer(kind=ENTIER) :: iter, i_sol_vtu, fn_vtu, fn_pvtu, i
  character(len=255) :: fln, fln_adj

  character(len=255) :: input_file

  call MPI_INIT(mpi_ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

  if (command_argument_count() >= 1) then
    call get_command_argument(1, input_file)
  else
    input_file = "input_data.f"
  end if
  call read_params(trim(input_file))

  call read_mesh_msh(mesh, meshfile_path, meshfile, &
    n_bc, bc_name, me, num_procs, mpi_send_recv)
  call build_mesh(mesh, num_procs, mpi_send_recv, .true., boundary_2d)
  call compute_geometry_mesh(mesh, .true., boundary_2d)
  if (order >= 3) call compute_cell_moments(mesh)

  allocate(sol(6, mesh%n_elems), sol1(6, mesh%n_elems), sol2(6, mesh%n_elems))
  allocate(prim(5, mesh%n_elems), rhs(6, mesh%n_elems))
  allocate(sum_lambda(mesh%n_elems))
  if (use_rk4) then
    allocate(sol3(6, mesh%n_elems))
    allocate(k1(6, mesh%n_elems), k2(6, mesh%n_elems), k3(6, mesh%n_elems), k4(6, mesh%n_elems))
  end if

  call init_sol(mesh, sol)
  ! Ghost cells get an analytically-correct sol from init_sol on every rank
  ! independently (each rank knows its own ghost cells' geometry), so this
  ! exchange is a no-op for t=0 -- it matters starting at RK stage 1 below,
  ! see the comment there.
  if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 6_ENTIER, sol)

  ! Initial output
  i_sol_vtu = 0
  call write_vtu(mesh, sol, me, i_sol_vtu)
  if (compute_error) then
    call compute_error_vortex(mesh, sol, 0.0_DOUBLE, h_err, l2err)
    if (me == 0) print *, "t=", 0.0_DOUBLE, "h=", h_err, "L2(rho)=", l2err
  end if

  ! Compute initial dt for RK3 stages
  call compute_prim(mesh, sol, prim)
  call compute_rhs(mesh, sol, prim, rhs, sum_lambda, 0.0_DOUBLE, num_procs, mpi_send_recv)
  dt = compute_dt(mesh, sum_lambda)

  t    = 0.0_DOUBLE
  iter = 1
  do while (t < tmax)
    if (t + dt > tmax) dt = tmax - t

    if (use_rk4) then
      ! Classical (non-SSP) 4-stage RK4, genuinely 4th-order in time --
      ! added 2026-09-18 per the user, to rule out SSP-RK3's own O(dt^3)
      ! temporal error as the reason order 4's spatial accuracy doesn't
      ! separate from order 3 on some benchmarks. Same time-dependent-BC
      ! simplification as SSP-RK3 below (every stage's compute_rhs uses
      ! the step's start time t, not each stage's own t/t+dt/2/t+dt).
      !   k1 = L(u0);          k2 = L(u0+dt/2*k1); k3 = L(u0+dt/2*k2)
      !   k4 = L(u0+dt*k3);    u^{n+1} = u0 + dt/6*(k1+2k2+2k3+k4)
      call compute_prim(mesh, sol, prim)
      call compute_rhs(mesh, sol, prim, rhs, sum_lambda, t, num_procs, mpi_send_recv)
      dt = min(dt, compute_dt(mesh, sum_lambda))
      k1 = rhs
      sol1 = sol + 0.5_DOUBLE * dt * k1
      if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 6_ENTIER, sol1)
      call sync_gamma_arr(mesh, sol1)

      call compute_prim(mesh, sol1, prim)
      call compute_rhs(mesh, sol1, prim, rhs, sum_lambda, t, num_procs, mpi_send_recv)
      k2 = rhs
      sol2 = sol + 0.5_DOUBLE * dt * k2
      if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 6_ENTIER, sol2)
      call sync_gamma_arr(mesh, sol2)

      call compute_prim(mesh, sol2, prim)
      call compute_rhs(mesh, sol2, prim, rhs, sum_lambda, t, num_procs, mpi_send_recv)
      k3 = rhs
      sol3 = sol + dt * k3
      if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 6_ENTIER, sol3)
      call sync_gamma_arr(mesh, sol3)

      call compute_prim(mesh, sol3, prim)
      call compute_rhs(mesh, sol3, prim, rhs, sum_lambda, t, num_procs, mpi_send_recv)
      k4 = rhs
      sol = sol + (dt / 6.0_DOUBLE) * (k1 + 2.0_DOUBLE * k2 + 2.0_DOUBLE * k3 + k4)
      if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 6_ENTIER, sol)
      call sync_gamma_arr(mesh, sol)
    else
    ! SSP-RK3 stage 1: u1 = u0 + dt * L(u0). A time-dependent BC uses the
    ! step's start time t for all 3 stages below (not exact substage
    ! timing t, t+dt, t+dt/2 -- see compute_rhs).
    call compute_prim(mesh, sol, prim)
    call compute_rhs(mesh, sol, prim, rhs, sum_lambda, t, num_procs, mpi_send_recv)
    dt   = min(dt, compute_dt(mesh, sum_lambda))
    sol1 = sol + dt * rhs
    ! Under MPI, sol1's ghost cells were just formed from THIS rank's own
    ! rhs/sol, not the neighbour rank's -- compute_rhs's own flux stencil
    ! reads sol1 at ghost cells next (stage 2), so they must be refreshed
    ! from the owning rank first. Found 2026-09-13: this exchange (and the
    ! two below) were entirely missing -- num_procs>1 silently ran every
    ! partition boundary as if it were a wall/extrapolated edge, using
    ! whatever stale sol value the ghost cell last held (its t=0 analytic
    ! IC, forever) instead of the evolving neighbour solution. Caught by
    ! comparing a 6-rank order-1 vortex run against its serial reference:
    ! L2(rho) differed by 14% at h=0.25, where an MPI-transparent scheme
    ! must match to roundoff.
    if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 6_ENTIER, sol1)
    call sync_gamma_arr(mesh, sol1)

    ! SSP-RK3 stage 2: u2 = 3/4*u0 + 1/4*(u1 + dt*L(u1))
    call compute_prim(mesh, sol1, prim)
    call compute_rhs(mesh, sol1, prim, rhs, sum_lambda, t, num_procs, mpi_send_recv)
    sol2 = 0.75_DOUBLE * sol + 0.25_DOUBLE * (sol1 + dt * rhs)
    if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 6_ENTIER, sol2)
    call sync_gamma_arr(mesh, sol2)

    ! SSP-RK3 stage 3: u^{n+1} = 1/3*u0 + 2/3*(u2 + dt*L(u2))
    call compute_prim(mesh, sol2, prim)
    call compute_rhs(mesh, sol2, prim, rhs, sum_lambda, t, num_procs, mpi_send_recv)
    sol  = (1.0_DOUBLE/3.0_DOUBLE) * sol &
         + (2.0_DOUBLE/3.0_DOUBLE) * (sol2 + dt * rhs)
    if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 6_ENTIER, sol)
    call sync_gamma_arr(mesh, sol)
    end if

    t    = t + dt
    iter = iter + 1

    ! New dt for next step
    call compute_prim(mesh, sol, prim)
    call compute_rhs(mesh, sol, prim, rhs, sum_lambda, t, num_procs, mpi_send_recv)
    dt = compute_dt(mesh, sum_lambda)

    if (me == 0 .and. mod(iter, 100) == 0) print *, "iter=", iter, "t=", t, "dt=", dt

    ! Periodic output
    if (n_sol_vtu > 1) then
      if (t >= real(i_sol_vtu, DOUBLE) * tmax / real(n_sol_vtu - 1, DOUBLE)) then
        call write_vtu(mesh, sol, me, i_sol_vtu)
        if (compute_error) then
          call compute_error_vortex(mesh, sol, t, h_err, l2err)
          if (me == 0) print *, "t=", t, "h=", h_err, "L2(rho)=", l2err
        end if
        i_sol_vtu = i_sol_vtu + 1
      end if
    end if
  end do

  ! Final output
  call write_vtu(mesh, sol, me, -1)
  if (compute_error) then
    call compute_error_vortex(mesh, sol, t, h_err, l2err)
    if (me == 0) print *, "FINAL t=", t, "h=", h_err, "L2(rho)=", l2err
  end if

  call MPI_FINALIZE(mpi_ierr)

contains

  subroutine write_vtu(mesh, sol, me, idx)
    use euler_ho_module, only: compute_prim, gamma_arr, boundary_2d
    use arbitrary_high_order_module, only: compute_next_order_derivative, &
      aho_module_use_green_gauss => use_green_gauss
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(in) :: sol
    integer, intent(in) :: me, idx
    integer(kind=ENTIER) :: fn_v, fn_pv
    character(len=255) :: fname
    real(kind=DOUBLE), allocatable :: prim_loc(:, :), rho(:), p(:), temp(:)
    real(kind=DOUBLE), allocatable :: centroid(:, :), velocity(:, :)
    real(kind=DOUBLE), parameter :: r_gas = 287.0_DOUBLE ! air, for T = p/(rho*R)
    integer(kind=ENTIER) :: i, iv
    ! Density-gradient diagnostic output (2026-09-18, per the user): lets
    ! a schlieren-style render just ColorBy this field instead of running
    ! ParaView's own Gradient filter, which was found very expensive at
    ! scale on these meshes. Computed fresh at output time with a direct
    ! call to compute_next_order_derivative, independent of whatever
    ! reconstruction the run itself used for its own RK stages -- always
    ! forced to the Green-Gauss nodal fit (cheap, no LS solve) regardless
    ! of this run's own use_green_gauss setting, since this is a one-shot
    ! diagnostic, not part of the flux computation. dphi_v_out/valid_v_out
    ! expose the same per-vertex pass at no extra cost.
    logical :: saved_use_gg
    real(kind=DOUBLE), allocatable :: dphi_cell(:, :), dphi_v(:, :)
    logical, allocatable :: valid_v(:)
    real(kind=DOUBLE), allocatable :: grad_rho_cell(:, :), grad_rho_vert(:, :)

    write(fname, '(a,i0)') 'output_', idx

    allocate(prim_loc(5, mesh%n_elems))
    call compute_prim(mesh, sol, prim_loc)

    allocate(rho(mesh%n_elems), p(mesh%n_elems), temp(mesh%n_elems))
    allocate(centroid(3, mesh%n_elems), velocity(3, mesh%n_elems))
    do i = 1, mesh%n_elems
      rho(i) = prim_loc(1, i)
      velocity(:, i) = prim_loc(2:4, i)
      p(i)   = prim_loc(5, i)
      temp(i) = p(i) / (max(rho(i), 1.0e-16_DOUBLE) * r_gas)
      centroid(:, i) = mesh%elem(i)%coord
    end do

    saved_use_gg = aho_module_use_green_gauss
    aho_module_use_green_gauss = .true.
    allocate(dphi_cell(15, mesh%n_elems))
    call compute_next_order_derivative(mesh, 3_ENTIER, 5_ENTIER, boundary_2d, &
      prim_loc, dphi_cell, deriv_order=1_ENTIER, &
      dphi_v_out=dphi_v, valid_v_out=valid_v)
    aho_module_use_green_gauss = saved_use_gg

    ! prim's variable 1 is rho; dphi(:,e) is laid out (dir-1)*5+v, so
    ! indices 1/6/11 are d(rho)/dx, d(rho)/dy, d(rho)/dz (same convention
    ! as euler_ho_module's own grad_flat).
    allocate(grad_rho_cell(3, mesh%n_elems))
    grad_rho_cell(1, :) = dphi_cell(1, :)
    grad_rho_cell(2, :) = dphi_cell(6, :)
    grad_rho_cell(3, :) = dphi_cell(11, :)

    ! Boundary vertices are skipped entirely by compute_next_order_
    ! derivative's own accumulation loop (mesh%vert%is_bound, see that
    ! subroutine), leaving dphi_v_cache/valid_cache at whatever
    ! uninitialised allocate() gave them there -- valid_v(iv) is NOT
    ! trustworthy on its own for those. Check is_bound directly first
    ! (always well-defined) and only then fall back to valid_v for a
    ! genuinely-interior-but-degenerate stencil.
    allocate(grad_rho_vert(3, mesh%n_vert))
    do iv = 1, mesh%n_vert
      if (mesh%vert(iv)%is_bound .or. .not. valid_v(iv)) then
        grad_rho_vert(:, iv) = 0.0_DOUBLE
      else
        grad_rho_vert(1, iv) = dphi_v(1, iv)
        grad_rho_vert(2, iv) = dphi_v(6, iv)
        grad_rho_vert(3, iv) = dphi_v(11, iv)
      end if
    end do

    call open_file_vtu(mesh, trim(adjustl(fname)), fn_v, fn_pv)
    call write_file_vtu_start_cell_data(mesh, trim(adjustl(fname)), fn_v, fn_pv)
    call write_file_vtu_cell_scalar(mesh, trim(adjustl(fname)), fn_v, fn_pv, rho, 'rho')
    call write_file_vtu_cell_vector(mesh, trim(adjustl(fname)), fn_v, fn_pv, velocity, 'velocity')
    call write_file_vtu_cell_scalar(mesh, trim(adjustl(fname)), fn_v, fn_pv, p,   'p')
    call write_file_vtu_cell_scalar(mesh, trim(adjustl(fname)), fn_v, fn_pv, temp, 'T')
    call write_file_vtu_cell_scalar(mesh, trim(adjustl(fname)), fn_v, fn_pv, gamma_arr, 'gamma')
    call write_file_vtu_cell_vector(mesh, trim(adjustl(fname)), fn_v, fn_pv, centroid, 'Centroid')
    call write_file_vtu_cell_vector(mesh, trim(adjustl(fname)), fn_v, fn_pv, grad_rho_cell, 'grad_rho_cell')
    call write_file_vtu_end_cell_data(mesh, trim(adjustl(fname)), fn_v, fn_pv)
    call write_file_vtu_start_vert_data(mesh, trim(adjustl(fname)), fn_v, fn_pv)
    call write_file_vtu_vert_vector(mesh, trim(adjustl(fname)), fn_v, fn_pv, grad_rho_vert, 'grad_rho_vert')
    call write_file_vtu_end_vert_data(mesh, trim(adjustl(fname)), fn_v, fn_pv)
    call close_file_vtu(mesh, trim(adjustl(fname)), fn_v, fn_pv)

    deallocate(prim_loc, rho, p, temp, centroid, velocity)
    deallocate(dphi_cell, dphi_v, valid_v, grad_rho_cell, grad_rho_vert)
  end subroutine write_vtu

end program euler_ho_main
