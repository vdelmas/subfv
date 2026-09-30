program main
    use mpi
    use precision_module
    use mpi_module
    use mesh_module
    use mesh_reading_module
    use mesh_geometry_module
    use mesh_connectivity_module
    use io_module
    use silvia_base_module
    use silvia_errors_module
    use ns_global_data_module, only: ns_bc_style => bc_style, ns_boundary_2d => boundary_2d
  
    implicit none

    integer :: funit
    integer :: mpi_ierr, me, num_procs
    type(mesh_type) :: mesh
    type(mpi_send_recv_type) :: mpi_send_recv

    integer(kind=ENTIER) :: i, le, re, i_bc
    integer(kind=ENTIER) :: iter, iter_write_sol
    integer(kind=ENTIER) :: n_elems_loc, n_elems_ghost, n_elems_tot
    
    real(kind=DOUBLE) :: area, h_err, l2err
    real(kind=DOUBLE), dimension(3) :: n
    real(kind=DOUBLE), dimension(:,:), allocatable :: sol, sol_w, sol1, sol2
    real(kind=DOUBLE), dimension(:,:), allocatable :: rhs
    real(kind=DOUBLE), allocatable :: sum_lambda(:)

    real(kind=DOUBLE), allocatable :: grad(:, :, :)  ! (5, 3, n_elems)
    real(kind=DOUBLE), allocatable :: hess(:, :, :, :) ! (5, 3, 3, n_elems)

    call MPI_INIT(mpi_ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    call read_input_parameters('input_data.f')
    ! The ns nodal solver reads these from its own global module
    ns_bc_style = bc_style
    ns_boundary_2d = boundary_2d
    call init_flags()

    if (me == 0) print *, "SUBFV METHOD ORDER ", order, "TYPE ", scheme_id

    call read_mesh_msh(mesh, meshfile_path, meshfile, &
        n_bc, bc_name, me, num_procs, mpi_send_recv)
    call build_mesh(mesh, num_procs, mpi_send_recv, .true., boundary_2d)
    call compute_geometry_mesh(mesh, .true., boundary_2d)
    if (order >= 3) call compute_cell_moments(mesh)

    allocate(sol(5,mesh%n_elems))
    allocate(grad(5,3,mesh%n_elems))
    allocate(hess(5,3,3,mesh%n_elems))
    allocate(sol_w(5,mesh%n_elems))
    allocate(rhs(5,mesh%n_elems))
    allocate(sum_lambda(mesh%n_elems))

    if (order>=2) then
      allocate(sol1(5,mesh%n_elems))
      if (order>=3) then
        allocate(sol2(5,mesh%n_elems))
      end if
    end if

    do i=1, mesh%n_elems
      sol(:,i) = (/1.0_DOUBLE,1.0_DOUBLE,1.0_DOUBLE,1.0_DOUBLE,100.0_DOUBLE/)
      sol_w(:,i) = conserv_to_primit(sol(:,i))
      rhs(:,i) = (/0.0_DOUBLE,0.0_DOUBLE,0.0_DOUBLE,0.0_DOUBLE,0.0_DOUBLE/)
      sum_lambda(i) = 0.0_DOUBLE
    end do

    call init_sol(mesh, sol)
    if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol)
    
    ! call face_geometry(mesh)

    t=0.0_DOUBLE
    dt=0.0_DOUBLE
    iter=0
    iter_write_sol = 0
    call write_vtu(mesh, sol, me, iter_write_sol)
    if (me==0) print*, 'output save at iter=', iter, ' dt=', dt, ' time=', t
    iter_write_sol = 1
    
    do while (t<tmax .AND. iter<n_max_iter)
    ! do while (t<tmax)  
      if (t + dt > tmax) dt = tmax - t

      call mpi_barrier(mpi_comm_world, mpi_ierr)

      ! SSP-RK3 stage 1: u1 = u0 + dt * L(u0). A time-dependent BC uses the
      do i=1, mesh%n_elems
        sol_w(:,i) = conserv_to_primit(sol(:,i))
      end do
      call compute_rhs(mesh, sol, sol_w, grad, hess, rhs, sum_lambda, t, num_procs, mpi_send_recv)
      call compute_dt(mesh, sum_lambda, dt)
      call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, mpi_ierr)

      select case (order)
        case (1)
          !
          sol = sol + dt * rhs
          if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol)  
          !
        case (2)
          !
          ! First step of SSP RK 
          sol1 = sol + dt * rhs
          if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol1)
          ! Second step of SSP RK 
          do i=1, mesh%n_elems
            sol_w(:,i) = conserv_to_primit(sol1(:,i))
          end do
          call compute_rhs(mesh, sol1, sol_w, grad, hess, rhs, sum_lambda, t, num_procs, mpi_send_recv)
          sol = 0.5_DOUBLE * sol + 0.5_DOUBLE * (sol1 + dt * rhs)
          if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol)
          ! 
        case (3)
          ! 
          ! First step of SSP RK
          sol1 = sol + dt * rhs
          if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol1)
          ! Second step of SSP RK 
          do i=1, mesh%n_elems
            sol_w(:,i) = conserv_to_primit(sol1(:,i))
          end do
          call compute_rhs(mesh, sol1, sol_w, grad, hess, rhs, sum_lambda, t, num_procs, mpi_send_recv)
          sol2 = 0.75_DOUBLE * sol + 0.25_DOUBLE * (sol1 + dt * rhs)
          if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol2)
          ! Third step of SSP RK 
          do i=1, mesh%n_elems
            sol_w(:,i) = conserv_to_primit(sol2(:,i))
          end do
          call compute_rhs(mesh, sol2, sol_w, grad, hess, rhs, sum_lambda, t, num_procs, mpi_send_recv)
          sol  = (1.0_DOUBLE/3.0_DOUBLE) * sol + (2.0_DOUBLE/3.0_DOUBLE) * (sol2 + dt * rhs)
          if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol)
          !
        case default
          !
          sol = sol + dt * rhs
          if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol)  
          ! 
        end select

      ! Periodic output
      if (n_iter_write_sol > 1) then
        if (t >= real(iter_write_sol, DOUBLE) * tmax / real(n_iter_write_sol - 1, DOUBLE)) then
          call write_vtu(mesh, sol, me, iter_write_sol)
          if (me==0) print*, 'output save at iter=', iter, ' dt=', dt, ' time=', t
          if (compute_error) then
            call compute_error_test(mesh, sol, t, 1.0_DOUBLE, h_err, l2err)
            call MPI_ALLREDUCE(MPI_IN_PLACE, h_err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE, l2err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
            if (me == 0) print *, "t=", t, "h=", h_err, "L2(rho)=", l2err
          end if
          iter_write_sol = iter_write_sol + 1
        end if
      end if
      t = t + dt
      iter = iter + 1
    end do

    ! Final output
    call write_vtu(mesh, sol, me, -1)
    if (me==0) print*, 'output save at iter=', iter, ' dt=', dt, ' time=', t
    if (compute_error) then
      call compute_error_test(mesh, sol, t, 1.0_DOUBLE, h_err, l2err)
      call MPI_ALLREDUCE(MPI_IN_PLACE, h_err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, l2err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
      if (me == 0) print *, "FINAL t=", t, "h=", h_err, "L2(rho)=", l2err
    end if

    call MPI_FINALIZE(mpi_ierr)
contains

subroutine write_vtu(mesh, sol, me, idx)
  implicit none 

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    integer, intent(in) :: me, idx
    integer(kind=ENTIER) :: fn_v, fn_pv
    character(len=255) :: fname
    real(kind=DOUBLE), allocatable :: prim_loc(:, :), rho(:), ux(:), uy(:), uz(:), p(:), temp(:)
    real(kind=DOUBLE), allocatable :: centroid(:, :)
    real(kind=DOUBLE), parameter :: r_gas = 287.0_DOUBLE ! air, for T = p/(rho*R)
    integer(kind=ENTIER) :: i, iv
    logical :: saved_use_gg
    real(kind=DOUBLE), allocatable :: dphi_cell(:, :), dphi_v(:, :)
    logical, allocatable :: valid_v(:)
    real(kind=DOUBLE), allocatable :: grad_rho_cell(:, :), grad_rho_vert(:, :)

    write(fname, '(a,i0)') 'output_', idx

    allocate(prim_loc(5, mesh%n_elems))

    do i=1, mesh%n_elems
        prim_loc(:, i) = conserv_to_primit(sol(:, i))
    end do

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

    ! Density-gradient diagnostic (for schlieren), same as euler_ho's write_vtu:
    ! computed at output time with the aho Green-Gauss nodal fit.
    saved_use_gg = aho_use_green_gauss
    aho_use_green_gauss = .true.
    allocate(dphi_cell(15, mesh%n_elems))
    call compute_next_order_derivative(mesh, 3_ENTIER, 5_ENTIER, boundary_2d, &
      prim_loc, dphi_cell, deriv_order=1_ENTIER, &
      dphi_v_out=dphi_v, valid_v_out=valid_v)
    aho_use_green_gauss = saved_use_gg

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
    call write_file_vtu_cell_scalar(mesh, trim(adjustl(fname)), fn_v, fn_pv, ux,  'u')
    call write_file_vtu_cell_scalar(mesh, trim(adjustl(fname)), fn_v, fn_pv, uy,  'v')
    call write_file_vtu_cell_scalar(mesh, trim(adjustl(fname)), fn_v, fn_pv, uz,  'w')
    call write_file_vtu_cell_scalar(mesh, trim(adjustl(fname)), fn_v, fn_pv, p,   'p')
    call write_file_vtu_cell_scalar(mesh, trim(adjustl(fname)), fn_v, fn_pv, temp, 'T')
    call write_file_vtu_cell_vector(mesh, trim(adjustl(fname)), fn_v, fn_pv, centroid, 'Centroid')
    call write_file_vtu_cell_vector(mesh, trim(adjustl(fname)), fn_v, fn_pv, grad_rho_cell, 'grad_rho_cell')
    call write_file_vtu_end_cell_data(mesh, trim(adjustl(fname)), fn_v, fn_pv)
    call write_file_vtu_start_vert_data(mesh, trim(adjustl(fname)), fn_v, fn_pv)
    call write_file_vtu_vert_vector(mesh, trim(adjustl(fname)), fn_v, fn_pv, grad_rho_vert, 'grad_rho_vert')
    call write_file_vtu_end_vert_data(mesh, trim(adjustl(fname)), fn_v, fn_pv)
    call close_file_vtu(mesh, trim(adjustl(fname)), fn_v, fn_pv)

    deallocate(prim_loc, rho, ux, uy, uz, p, temp, centroid)
    deallocate(dphi_cell, dphi_v, valid_v, grad_rho_cell, grad_rho_vert)
  end subroutine write_vtu
end program main


