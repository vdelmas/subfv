program main
    use mpi
    use precision_module
    use mpi_module
    use mesh_module
    use mesh_reading_module
    use mesh_geometry_module
    use mesh_connectivity_module
    use io_module
    use silvia_module
    implicit none

    integer :: funit
    integer :: mpi_ierr, me, num_procs
    type(mesh_type) :: mesh
    type(mpi_send_recv_type) :: mpi_send_recv

    integer(kind=ENTIER) :: i, le, re, i_bc
    integer(kind=ENTIER) :: iter
    
    real(kind=DOUBLE) :: area
    real(kind=DOUBLE), dimension(3) :: n
    real(kind=DOUBLE), dimension(:,:), allocatable :: sol
    real(kind=DOUBLE), dimension(:,:), allocatable :: rhs
    
    namelist /INPUT_PARAM/ &
        meshfile_path, meshfile, &
        boundary_2d, bc_style, &
        n_bc, bc_name, bc_type, bc_val, &
        bc_type_V, bc_val_V, bc_type_T, bc_val_T, &
        periodic_mesh, use_cylinder_map, &
        dim_cylinder_map, &
        n_iter_print, n_iter_write_sol, &
        n_max_iter, final_res, &
        second_order, activate_diffusion, &
        local_time_step, cfl, tmax, &
        init_gresho, &
        error_2d, error_2d_h, exclude_bound_vert, &
        scheme, write_residual

    call MPI_INIT(mpi_ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    open (newunit=funit, file=trim(adjustl("input_data.f")))
    read (nml=INPUT_PARAM, unit=funit)
    close (unit=funit)

    call read_mesh_msh(mesh, meshfile_path, meshfile, &
        n_bc, bc_name, me, num_procs, mpi_send_recv)
    call build_mesh(mesh, num_procs, mpi_send_recv, .true., boundary_2d)
    call compute_geometry_mesh(mesh, .true., boundary_2d)

    allocate(sol(5,mesh%n_elems))
    do i=1, mesh%n_elems
      sol(:, i) = (/1.0_DOUBLE,1.0_DOUBLE,1.0_DOUBLE,1.0_DOUBLE,100.0_DOUBLE/)
    end do

    call init_sol(mesh, sol)

    t=0.0_DOUBLE
    tmax=1.0_DOUBLE
    iter=0
    call write_vtu(mesh, sol, 0, iter)
    do while (t<tmax)
      dt=0.1_DOUBLE
      !call compute_dt(mesh, sol, dt)
      !call compute_rhs(mesh, sol, rhs)
      do i=1,mesh%n_elems
        !sol(:, i) = sol(:, i) - dt/mesh%elem(i)%volume * rhs(:, i)
        sol(:, i) = sol(:, i)
      end do
      t = t + dt
      iter = iter + 1
      call write_vtu(mesh, sol, 0, iter)
    end do

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
    integer(kind=ENTIER) :: i

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
end program main


