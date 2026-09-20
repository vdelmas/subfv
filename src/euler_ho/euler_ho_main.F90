! High-order explicit Euler solver: SSP-RK3 (Shu-Osher 1988) or classical RK4 in time.
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
  ! Component 6 of sol/sol1/sol2/sol3/rhs is the gamma-transport variable
  ! rho*Gamma (Gamma=1/(gamma-1)); it rides through the same RK/MPI arithmetic
  ! as the 5 physical equations.
  real(kind=DOUBLE), allocatable :: sol(:, :), sol1(:, :), sol2(:, :), sol3(:, :)
  real(kind=DOUBLE), allocatable :: prim(:, :), rhs(:, :)
  real(kind=DOUBLE), allocatable :: k1(:, :), k2(:, :), k3(:, :), k4(:, :)
  real(kind=DOUBLE), allocatable :: sum_lambda(:)
  real(kind=DOUBLE) :: t, dt, h_err, l2err
  integer(kind=ENTIER) :: iter, i_sol_vtu
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
  ! setup_wall_mirror (wall-tangent fit at boundary vertices) is implemented but NOT enabled yet:
  ! it still develops a slow instability at some wall-adjacent vertices on the cylinder-tunnel mesh
  ! (elevated but finite tangential gradient that compounds over ~100+ iterations into NaN) --
  ! reproduces even at single rank, so it's not an MPI/ghost issue, root cause still open. Leave
  ! disabled (falls through to the proven phantom-zero-gradient hack in
  ! arbitrary_high_order_module's compute_next_order_derivative) until that's understood.
  ! call setup_wall_mirror(mesh)
  if (order >= 3) call compute_cell_moments(mesh)

  allocate(sol(6, mesh%n_elems), sol1(6, mesh%n_elems), sol2(6, mesh%n_elems))
  allocate(prim(5, mesh%n_elems), rhs(6, mesh%n_elems))
  allocate(sum_lambda(mesh%n_elems))
  if (use_rk4) then
    allocate(sol3(6, mesh%n_elems))
    allocate(k1(6, mesh%n_elems), k2(6, mesh%n_elems), k3(6, mesh%n_elems), k4(6, mesh%n_elems))
  end if

  call init_sol(mesh, sol)
  if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 6_ENTIER, sol)

  i_sol_vtu = 0
  call write_vtu(mesh, sol, me, i_sol_vtu)
  if (compute_error) then
    call compute_error_vortex(mesh, sol, 0.0_DOUBLE, h_err, l2err)
    if (me == 0) print *, "t=", 0.0_DOUBLE, "h=", h_err, "L2(rho)=", l2err
  end if

  call compute_prim(mesh, sol, prim)
  call compute_rhs(mesh, sol, prim, rhs, sum_lambda, 0.0_DOUBLE, num_procs, mpi_send_recv)
  dt = compute_dt(mesh, sum_lambda)

  t    = 0.0_DOUBLE
  iter = 1
  do while (t < tmax)
    if (t + dt > tmax) dt = tmax - t

    if (use_rk4) then
      ! k1=L(u0); k2=L(u0+dt/2*k1); k3=L(u0+dt/2*k2); k4=L(u0+dt*k3);
      ! u^{n+1}=u0+dt/6*(k1+2k2+2k3+k4). Every stage's compute_rhs uses the
      ! step's start time t (not each stage's own substage time), as SSP-RK3 does too.
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
      ! SSP-RK3 (Shu-Osher 1988). Ghost cells are refreshed after every stage:
      ! compute_rhs's flux stencil reads the neighbour rank's sol at the next
      ! stage, and without this exchange a partition boundary silently behaves
      ! like a wall holding its t=0 value forever.
      call compute_prim(mesh, sol, prim)
      call compute_rhs(mesh, sol, prim, rhs, sum_lambda, t, num_procs, mpi_send_recv)
      dt   = min(dt, compute_dt(mesh, sum_lambda))
      sol1 = sol + dt * rhs
      if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 6_ENTIER, sol1)
      call sync_gamma_arr(mesh, sol1)

      call compute_prim(mesh, sol1, prim)
      call compute_rhs(mesh, sol1, prim, rhs, sum_lambda, t, num_procs, mpi_send_recv)
      sol2 = 0.75_DOUBLE * sol + 0.25_DOUBLE * (sol1 + dt * rhs)
      if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 6_ENTIER, sol2)
      call sync_gamma_arr(mesh, sol2)

      call compute_prim(mesh, sol2, prim)
      call compute_rhs(mesh, sol2, prim, rhs, sum_lambda, t, num_procs, mpi_send_recv)
      sol  = (1.0_DOUBLE/3.0_DOUBLE) * sol &
           + (2.0_DOUBLE/3.0_DOUBLE) * (sol2 + dt * rhs)
      if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 6_ENTIER, sol)
      call sync_gamma_arr(mesh, sol)
    end if

    t    = t + dt
    iter = iter + 1

    call compute_prim(mesh, sol, prim)
    call compute_rhs(mesh, sol, prim, rhs, sum_lambda, t, num_procs, mpi_send_recv)
    dt = compute_dt(mesh, sum_lambda)

    if (me == 0 .and. mod(iter, 100) == 0) print *, "iter=", iter, "t=", t, "dt=", dt

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
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(in) :: sol
    integer, intent(in) :: me, idx

    integer(kind=ENTIER) :: fn_v, fn_pv
    character(len=255) :: fname
    real(kind=DOUBLE), allocatable :: prim_loc(:, :), rho(:), p(:), temp(:)
    real(kind=DOUBLE), allocatable :: centroid(:, :), velocity(:, :)
    real(kind=DOUBLE), parameter :: r_gas = 287.0_DOUBLE
    integer(kind=ENTIER) :: i, iv
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

    ! Density-gradient diagnostic (for schlieren-style ColorBy): computed fresh
    ! at output time, always forced to the Green-Gauss nodal fit (cheap, no LS
    ! solve) regardless of this run's own use_green_gauss setting.
    saved_use_gg = aho_module_use_green_gauss
    aho_module_use_green_gauss = .true.
    allocate(dphi_cell(15, mesh%n_elems))
    call compute_next_order_derivative(mesh, 3_ENTIER, 5_ENTIER, boundary_2d, &
      prim_loc, dphi_cell, deriv_order=1_ENTIER, &
      dphi_v_out=dphi_v, valid_v_out=valid_v)
    aho_module_use_green_gauss = saved_use_gg

    ! dphi(:,e) is laid out (dir-1)*5+v; var 1 is rho, so indices 1/6/11 are
    ! d(rho)/dx, d(rho)/dy, d(rho)/dz.
    allocate(grad_rho_cell(3, mesh%n_elems))
    grad_rho_cell(1, :) = dphi_cell(1, :)
    grad_rho_cell(2, :) = dphi_cell(6, :)
    grad_rho_cell(3, :) = dphi_cell(11, :)

    ! Boundary vertices are skipped by the accumulation loop, so valid_v alone
    ! isn't trustworthy there -- check is_bound first.
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
