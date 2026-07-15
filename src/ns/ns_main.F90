program main
  use mpi

  use precision_module
  use mpi_module
  use mesh_module
  use mesh_geometry_module
  use mesh_reading_module
  use mesh_connectivity_module

  use ns_global_data_module
  use ns_io_module
  use ns_euler_module
  use ns_vectorial_diffusion_module
  use ns_diffusion_module

  implicit none

  integer(kind=ENTIER) :: me, num_procs, mpi_ierr
  integer(kind=ENTIER) :: fn_residual
  integer(kind=ENTIER) :: iter, i
  integer(kind=ENTIER) :: count_rate, count_max, start_count, end_count
  real(kind=DOUBLE) :: res
  real(kind=DOUBLE), dimension(5) :: res_vect
  integer(kind=ENTIER_D) :: teu_total
  real(kind=DOUBLE) :: wall_time

  type(mesh_type) :: mesh
  type(mpi_send_recv_type) :: mpi_send_recv
  real(kind=DOUBLE), dimension(:, :, :), allocatable :: mat_diag, cell_grad, &
    nodal_grad, mat_h_p
  real(kind=DOUBLE), dimension(:, :), allocatable    :: sol, rhs, vp, soln, delta_sol
  real(kind=DOUBLE), dimension(:), allocatable       :: sum_lambda, h_p

  call MPI_INIT(mpi_ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

  call read_input_parameters("input_data.f")
  call init_flags()

  call read_mesh_msh(mesh, meshfile_path, meshfile, &
    n_bc, bc_name, me, num_procs, mpi_send_recv)
  if (use_cylinder_map) call map_cylinder(mesh, dim_cylinder_map)
  call build_mesh(mesh, num_procs, mpi_send_recv, .TRUE., boundary_2d, periodic_mesh)
  call compute_geometry_mesh(mesh, .TRUE., boundary_2d, me, num_procs)
  if (rescale) then
    call rescale_mesh(mesh, rescale_factor)
    call compute_geometry_mesh(mesh, .TRUE., boundary_2d, me, num_procs)
  end if

  allocate(sol(5, mesh%n_elems))           ; sol       = 0.0_DOUBLE
  allocate(soln(5, mesh%n_elems))          ; soln      = 0.0_DOUBLE
  allocate(delta_sol(5, mesh%n_elems))     ; delta_sol = 0.0_DOUBLE
  allocate(rhs(5, mesh%n_elems))           ; rhs       = 0.0_DOUBLE
  allocate(sum_lambda(mesh%n_elems))       ; sum_lambda = 0.0_DOUBLE
  allocate(mat_diag(5, 5, mesh%n_elems))   ; mat_diag  = 0.0_DOUBLE
  allocate(cell_grad(3, 5, mesh%n_elems))  ; cell_grad  = 0.0_DOUBLE
  allocate(nodal_grad(3, 5, mesh%n_vert))  ; nodal_grad = 0.0_DOUBLE
  allocate(vp(3, mesh%n_vert))             ; vp        = 0.0_DOUBLE
  allocate(h_p(mesh%n_vert))               ; h_p       = 0.0_DOUBLE
  allocate(mat_h_p(3, 3, mesh%n_vert))     ; mat_h_p   = 0.0_DOUBLE
  if (scheme_id == SCHEME_ZB) then
    do i = 1, mesh%n_vert
      h_p(i) = compute_length(mesh, i)
      mat_h_p(:, :, i) = compute_ellip(mesh, i)
    end do
  end if

  call init_residual_file(me, fn_residual)

  t = 0.0_DOUBLE
  iter = 0
  res = final_res + 1.0_DOUBLE
  call init_sol(mesh, sol, me, num_procs, mpi_send_recv)
  if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol)
  if (n_iter_write_sol > 0) &
    call write_sol(mesh, sol, delta_sol, cell_grad, vp, h_p, mat_h_p, me, iter/n_iter_write_sol)

  call system_clock(count_rate=count_rate, count_max=count_max)
  call system_clock(start_count)

  do while (res > final_res .and. iter <= n_max_iter)
    call mpi_barrier(mpi_comm_world, mpi_ierr)
    call compute_rhs_ns(mesh, sol, mpi_send_recv, num_procs, &
      rhs, mat_diag, sum_lambda, cell_grad, nodal_grad, vp, h_p, mat_h_p)
    if (second_order) then
      soln = sol
      call add_timedisc_and_solve(mesh, sol, rhs, mat_diag, sum_lambda, &
        delta_sol, res, res_vect, num_procs, me, 2.0_DOUBLE, 1.0_DOUBLE, cfl, t, dt)
      if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol)
      call compute_rhs_ns(mesh, sol, mpi_send_recv, num_procs, &
        rhs, mat_diag, sum_lambda, cell_grad, nodal_grad, vp, h_p, mat_h_p)
      call add_timedisc_and_solve(mesh, soln, rhs, mat_diag, sum_lambda, &
        delta_sol, res, res_vect, num_procs, me, 1.0_DOUBLE, 1.0_DOUBLE, cfl, t, dt)
      sol = soln
    else
      call add_timedisc_and_solve(mesh, sol, rhs, mat_diag, sum_lambda, &
        delta_sol, res, res_vect, num_procs, me, 1.0_DOUBLE, 1.0_DOUBLE, cfl, t, dt)
    end if
    if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol)
    call write_iter_info(iter, res, res_vect, me, fn_residual, start_count, count_rate, &
      mesh, sol, delta_sol, cell_grad, nodal_grad, vp, h_p, mat_h_p)
    iter = iter + 1
    t = t + dt
    if (.not. local_time_step .and. t >= tmax) exit
  end do

  call system_clock(end_count)
  call write_sol(mesh, sol, delta_sol, cell_grad, vp, h_p, mat_h_p, me, -1)
  if (write_residual .and. me == 0) close(fn_residual)

  wall_time = real(end_count - start_count, DOUBLE) / real(count_rate, DOUBLE)
  teu_total = int(iter, ENTIER_D) * int(mesh%n_elems, ENTIER_D)
  if (num_procs > 1) then
    call MPI_ALLREDUCE(MPI_IN_PLACE, teu_total, 1, MPI_INTEGER8,         MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, end_count, 1, MPI_INT,              MPI_MAX, MPI_COMM_WORLD, mpi_ierr)
    wall_time = real(end_count - start_count, DOUBLE) / real(count_rate, DOUBLE)
  end if
  if (me == 0) then
    print *, "  --------------------------------------------------------"
    print *, "  Total element updates (TEU) =", teu_total
    print *, "  Total CPU time (wall)       =", wall_time
    print *, "  Total CPU time / TEU        =", wall_time / real(teu_total, DOUBLE)
    print *, "  --------------------------------------------------------"
  end if
  call MPI_FINALIZE(mpi_ierr)

contains
  subroutine compute_rhs_ns(mesh, sol, mpi_send_recv, num_procs, &
      rhs, mat_diag, sum_lambda, cell_grad, nodal_grad, vp, h_p, mat_h_p)
    use precision_module
    use mpi_module
    use mesh_module
    use ns_global_data_module, only: second_order, activate_diffusion, cfl, method
    use ns_euler_module, only: compute_euler_flux, compute_cell_grad_from_nodal_grad, &
      flux_limiting
    use ns_vectorial_diffusion_module, only: compute_mu, &
      build_momentum_diff_mat_into_global_mat_diag
    use ns_diffusion_module, only: build_diff_mat_into_global_mat_diag
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv
    integer(kind=ENTIER), intent(in) :: num_procs
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(out) :: rhs
    real(kind=DOUBLE), dimension(5, 5, mesh%n_elems), intent(out) :: mat_diag
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(out) :: sum_lambda
    real(kind=DOUBLE), dimension(3, 5, mesh%n_elems), intent(out) :: cell_grad
    real(kind=DOUBLE), dimension(3, 5, mesh%n_vert), intent(out) :: nodal_grad
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: h_p
    real(kind=DOUBLE), dimension(3, 3, mesh%n_vert), intent(in) :: mat_h_p

    real(kind=DOUBLE), allocatable :: rhs_o2(:,:), sum_lambda_o2(:), mu_arr(:)
    real(kind=DOUBLE) :: mu0corr

    rhs = 0.0_DOUBLE; sum_lambda = 0.0_DOUBLE
    mat_diag = 0.0_DOUBLE; cell_grad = 0.0_DOUBLE; nodal_grad = 0.0_DOUBLE

    call compute_euler_flux(mesh, sol, rhs, cell_grad, sum_lambda, &
      .false., vp, h_p, mat_h_p)
    if (second_order) then
      allocate(rhs_o2(5, mesh%n_elems), sum_lambda_o2(mesh%n_elems))
      rhs_o2 = 0.0_DOUBLE; sum_lambda_o2 = 0.0_DOUBLE
      call compute_cell_grad_from_nodal_grad(mesh, sol, cell_grad, nodal_grad, method)
      if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 3*5, cell_grad)
      call compute_euler_flux(mesh, sol, rhs_o2, cell_grad, sum_lambda_o2, &
        .true., vp, h_p, mat_h_p)
      call flux_limiting(mesh, sol, rhs, rhs_o2, sum_lambda, cfl)
    end if
    if (activate_diffusion) then
      allocate(mu_arr(mesh%n_elems))
      call compute_mu(mesh, sol, num_procs, mu_arr, mu0corr)
      call build_momentum_diff_mat_into_global_mat_diag(mesh, mat_diag, rhs, sol, mu_arr, mu0corr)
      call build_diff_mat_into_global_mat_diag(mesh, mat_diag, rhs, sol, mu_arr)
    end if
  end subroutine compute_rhs_ns
end program main