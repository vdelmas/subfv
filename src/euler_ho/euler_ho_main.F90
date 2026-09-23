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
  real(kind=DOUBLE), allocatable :: adapt_orig_scale(:), adapt_cum_disp(:, :)
  integer(kind=ENTIER) :: iter, i_sol_vtu, nan_local
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

    ! A diverged run used to march to tmax and exit 0, so a NaN solution looked
    ! like a success in the logs and only showed up as a blank figure.
    !
    ! Test the SOLUTION, not dt: compute_dt takes a min over cells, and
    ! gfortran's min returns the non-NaN operand, so NaN cells are silently
    ! skipped and dt stays finite all the way to tmax. Guarding on dt looked
    ! right and caught nothing.
    if (mod(iter, 50) == 0) then
      nan_local = 0
      if (any(sol /= sol)) nan_local = 1
      call MPI_ALLREDUCE(MPI_IN_PLACE, nan_local, 1, MPI_INTEGER8, MPI_MAX, &
        MPI_COMM_WORLD, mpi_ierr)
      if (nan_local /= 0) then
        if (me == 0) print *, "[-] solution contains NaN at iter", iter, " -- diverged"
        call MPI_ABORT(MPI_COMM_WORLD, 2, mpi_ierr)
      end if
    end if

    call maybe_adapt_mesh()

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

  ! One shock-alignment node-movement cycle, if this iteration is due for one.
  !
  ! Cycles are triggered on iteration count, not on a residual: this case has
  ! an unsteady wake and never settles, so there is no steady state to wait
  ! for. The flow only needs to be roughly established (adapt_start_iter)
  ! before the density gradient marks a meaningful front.
  subroutine maybe_adapt_mesh()
    use shock_adapt_move_module, only: compute_shock_sensor_grad_rho, &
      build_vert_adjacency, compute_node_displacement_curvature, move_mesh, &
      min_elem_volume, mpi_memory_exchange_vert, compute_local_scale
    use arbitrary_high_order_module, only: invalidate_geometry_caches
    implicit none

    integer(kind=ENTIER) :: n_done, n_moved, n_moved_tot, i
    real(kind=DOUBLE) :: max_disp, max_disp_tot, vmin
    real(kind=DOUBLE), allocatable :: rho(:), node_sensor(:), disp(:, :)
    logical, allocatable :: node_flagged(:)
    integer(kind=ENTIER), allocatable :: n_neigh(:), vneigh(:, :)

    if (n_adapt_cycles <= 0) return
    if (iter < adapt_start_iter) return
    n_done = (iter - adapt_start_iter)/adapt_interval_iter
    if (n_done >= n_adapt_cycles) return
    if (mod(iter - adapt_start_iter, adapt_interval_iter) /= 0) return

    allocate(rho(mesh%n_elems), node_sensor(mesh%n_vert))
    allocate(node_flagged(mesh%n_vert), disp(3, mesh%n_vert))
    allocate(n_neigh(mesh%n_vert), vneigh(16, mesh%n_vert))

    do i = 1, mesh%n_elems
      rho(i) = sol(1, i)
    end do

    call compute_shock_sensor_grad_rho(mesh, rho, adapt_grad_threshold, node_sensor, node_flagged)
    call build_vert_adjacency(mesh, n_neigh, vneigh)
    ! Freeze the pre-movement length scale on the first cycle: the cumulative
    ! cap must be measured against the mesh we started from, not the one the
    ! previous cycles already compressed.
    if (.not. allocated(adapt_orig_scale)) then
      allocate(adapt_orig_scale(mesh%n_vert), adapt_cum_disp(3, mesh%n_vert))
      call compute_local_scale(mesh, adapt_orig_scale)
      adapt_cum_disp = 0.0_DOUBLE
    end if

    call compute_node_displacement_curvature(mesh, node_flagged, n_neigh, vneigh, &
      adapt_max_move_frac, adapt_relax, disp, n_moved, max_disp, &
      adapt_orig_scale, adapt_cum_disp)
    call move_mesh(mesh, disp)

    ! A vertex on a partition cut is is_ghost=.false. on BOTH ranks sharing it,
    ! so each computed its own displacement from its own side's stencil; the
    ! exchange averages them so the two sides agree on one position.
    if (num_procs > 1) call mpi_memory_exchange_vert(mesh, mpi_send_recv)

    call compute_geometry_mesh(mesh, .true., boundary_2d)
    ! Every geometry cache is keyed on a size (n_vert/n_elems/n_faces) that a
    ! move does not change, so without these two the reconstruction and the
    ! face quadrature silently keep using the pre-move geometry.
    call invalidate_geometry_caches()
    call invalidate_face_quad_cache()
    if (order >= 3) call compute_cell_moments(mesh)

    vmin = min_elem_volume(mesh)
    call MPI_ALLREDUCE(MPI_IN_PLACE, vmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, mpi_ierr)
    n_moved_tot = n_moved
    call MPI_ALLREDUCE(MPI_IN_PLACE, n_moved_tot, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
    max_disp_tot = max_disp
    call MPI_ALLREDUCE(MPI_IN_PLACE, max_disp_tot, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, mpi_ierr)

    if (me == 0) print *, "[adapt] cycle", n_done + 1, "iter=", iter, "t=", t, &
      "moved=", n_moved_tot, "max_disp=", max_disp_tot, "min_vol=", vmin

    if (vmin <= 0.0_DOUBLE) then
      if (me == 0) print *, "[-] adapt produced a non-positive cell volume, stopping"
      call MPI_ABORT(MPI_COMM_WORLD, 1, mpi_ierr)
    end if

    deallocate(rho, node_sensor, node_flagged, disp, n_neigh, vneigh)
  end subroutine maybe_adapt_mesh

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
    real(kind=DOUBLE), allocatable :: prim_loc(:, :), rho(:), p(:), temp(:), h_tot(:)
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

    allocate(rho(mesh%n_elems), p(mesh%n_elems), temp(mesh%n_elems), h_tot(mesh%n_elems))
    allocate(centroid(3, mesh%n_elems), velocity(3, mesh%n_elems))
    do i = 1, mesh%n_elems
      rho(i) = prim_loc(1, i)
      velocity(:, i) = prim_loc(2:4, i)
      p(i)   = prim_loc(5, i)
      temp(i) = p(i) / (max(rho(i), 1.0e-16_DOUBLE) * r_gas)
      ! Total enthalpy h = e + p/rho = gamma/(gamma-1) p/rho + |u|^2/2. A Riemann invariant of
      ! the steady Euler equations along streamlines, hence uniformly h_inf for a steady flow fed
      ! by a uniform inflow -- the diagnostic that separates an enthalpy-preserving Riemann solver
      ! (flux_scheme='three_wave_enthalpy') from one that is not.
      h_tot(i) = gamma_arr(i) / (gamma_arr(i) - 1.0_DOUBLE) * p(i) / max(rho(i), 1.0e-16_DOUBLE) &
        + 0.5_DOUBLE * dot_product(prim_loc(2:4, i), prim_loc(2:4, i))
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
    call write_file_vtu_cell_scalar(mesh, trim(adjustl(fname)), fn_v, fn_pv, h_tot, 'H')
    call write_file_vtu_cell_scalar(mesh, trim(adjustl(fname)), fn_v, fn_pv, gamma_arr, 'gamma')
    call write_file_vtu_cell_vector(mesh, trim(adjustl(fname)), fn_v, fn_pv, centroid, 'Centroid')
    call write_file_vtu_cell_vector(mesh, trim(adjustl(fname)), fn_v, fn_pv, grad_rho_cell, 'grad_rho_cell')
    call write_file_vtu_end_cell_data(mesh, trim(adjustl(fname)), fn_v, fn_pv)
    call write_file_vtu_start_vert_data(mesh, trim(adjustl(fname)), fn_v, fn_pv)
    call write_file_vtu_vert_vector(mesh, trim(adjustl(fname)), fn_v, fn_pv, grad_rho_vert, 'grad_rho_vert')
    call write_file_vtu_end_vert_data(mesh, trim(adjustl(fname)), fn_v, fn_pv)
    call close_file_vtu(mesh, trim(adjustl(fname)), fn_v, fn_pv)

    deallocate(prim_loc, rho, p, temp, h_tot, centroid, velocity)
    deallocate(dphi_cell, dphi_v, valid_v, grad_rho_cell, grad_rho_vert)
  end subroutine write_vtu

end program euler_ho_main
