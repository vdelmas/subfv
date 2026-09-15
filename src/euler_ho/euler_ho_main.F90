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

  real(kind=DOUBLE), allocatable :: sol(:, :), sol1(:, :), sol2(:, :)
  real(kind=DOUBLE), allocatable :: prim(:, :), rhs(:, :)
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

  allocate(sol(5, mesh%n_elems), sol1(5, mesh%n_elems), sol2(5, mesh%n_elems))
  allocate(prim(5, mesh%n_elems), rhs(5, mesh%n_elems))
  allocate(sum_lambda(mesh%n_elems))

  call init_sol(mesh, sol)
  ! Ghost cells get an analytically-correct sol from init_sol on every rank
  ! independently (each rank knows its own ghost cells' geometry), so this
  ! exchange is a no-op for t=0 -- it matters starting at RK stage 1 below,
  ! see the comment there.
  if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5_ENTIER, sol)

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
    if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5_ENTIER, sol1)

    ! SSP-RK3 stage 2: u2 = 3/4*u0 + 1/4*(u1 + dt*L(u1))
    call compute_prim(mesh, sol1, prim)
    call compute_rhs(mesh, sol1, prim, rhs, sum_lambda, t, num_procs, mpi_send_recv)
    sol2 = 0.75_DOUBLE * sol + 0.25_DOUBLE * (sol1 + dt * rhs)
    if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5_ENTIER, sol2)

    ! SSP-RK3 stage 3: u^{n+1} = 1/3*u0 + 2/3*(u2 + dt*L(u2))
    call compute_prim(mesh, sol2, prim)
    call compute_rhs(mesh, sol2, prim, rhs, sum_lambda, t, num_procs, mpi_send_recv)
    sol  = (1.0_DOUBLE/3.0_DOUBLE) * sol &
         + (2.0_DOUBLE/3.0_DOUBLE) * (sol2 + dt * rhs)
    if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5_ENTIER, sol)

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
    use euler_ho_module, only: compute_prim
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    integer, intent(in) :: me, idx
    integer(kind=ENTIER) :: fn_v, fn_pv
    character(len=255) :: fname
    real(kind=DOUBLE), allocatable :: prim_loc(:, :), rho(:), ux(:), uy(:), uz(:), p(:), temp(:)
    real(kind=DOUBLE), allocatable :: centroid(:, :)
    real(kind=DOUBLE), parameter :: r_gas = 287.0_DOUBLE ! air, for T = p/(rho*R)
    integer(kind=ENTIER) :: i

    write(fname, '(a,i0)') 'output_', idx

    allocate(prim_loc(5, mesh%n_elems))
    call compute_prim(mesh, sol, prim_loc)

    allocate(rho(mesh%n_elems), ux(mesh%n_elems), uy(mesh%n_elems))
    allocate(uz(mesh%n_elems), p(mesh%n_elems), temp(mesh%n_elems))
    allocate(centroid(3, mesh%n_elems))
    do i = 1, mesh%n_elems
      rho(i) = prim_loc(1, i)
      ux(i)  = prim_loc(2, i)
      uy(i)  = prim_loc(3, i)
      uz(i)  = prim_loc(4, i)
      p(i)   = prim_loc(5, i)
      temp(i) = p(i) / (max(rho(i), 1.0e-16_DOUBLE) * r_gas)
      centroid(:, i) = mesh%elem(i)%coord
    end do

    call open_file_vtu(mesh, trim(adjustl(fname)), fn_v, fn_pv)
    call write_file_vtu_start_cell_data(mesh, trim(adjustl(fname)), fn_v, fn_pv)
    call write_file_vtu_cell_scalar(mesh, trim(adjustl(fname)), fn_v, fn_pv, rho, 'rho')
    call write_file_vtu_cell_scalar(mesh, trim(adjustl(fname)), fn_v, fn_pv, ux,  'u')
    call write_file_vtu_cell_scalar(mesh, trim(adjustl(fname)), fn_v, fn_pv, uy,  'v')
    call write_file_vtu_cell_scalar(mesh, trim(adjustl(fname)), fn_v, fn_pv, uz,  'w')
    call write_file_vtu_cell_scalar(mesh, trim(adjustl(fname)), fn_v, fn_pv, p,   'p')
    call write_file_vtu_cell_scalar(mesh, trim(adjustl(fname)), fn_v, fn_pv, temp, 'T')
    call write_file_vtu_cell_vector(mesh, trim(adjustl(fname)), fn_v, fn_pv, centroid, 'Centroid')
    call write_file_vtu_end_cell_data(mesh, trim(adjustl(fname)), fn_v, fn_pv)
    call close_file_vtu(mesh, trim(adjustl(fname)), fn_v, fn_pv)

    deallocate(prim_loc, rho, ux, uy, uz, p, temp, centroid)
  end subroutine write_vtu

end program euler_ho_main
