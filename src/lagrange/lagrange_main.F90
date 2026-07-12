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
  integer(kind=ENTIER) :: n_bc, iter, i, method_length = 0
  real(kind=DOUBLE) :: cfl, cfl_max = 0.8, b2d_h = 1.0
  integer(kind=ENTIER) :: init = 0
  character(len=255) :: scheme = ""
  character(len=255), dimension(n_max_bc) :: bc_name
  character(len=255), dimension(n_max_bc) :: bc_type
  real(kind=DOUBLE), dimension(5, n_max_bc) :: bc_val
  real(kind=DOUBLE), dimension(5) :: sol_uniform
  type(mesh_type) :: mesh

  logical :: boundary_2d = .true.

  real(kind=DOUBLE) :: t, t_max, dt
  real(kind=DOUBLE), dimension(:,:), allocatable :: vp
  real(kind=DOUBLE), dimension(:), allocatable :: mass, pp
  real(kind=DOUBLE), dimension(:,:), allocatable :: sol
  real(kind=DOUBLE), dimension(:,:), allocatable :: rhs
  real(kind=DOUBLE), dimension(:,:), allocatable :: new_sol

  integer(kind=ENTIER) :: n_sol_vtu=2
  integer(kind=ENTIER) :: i_sol_vtu

  namelist /INPUT_PARAM/ &
    meshfile_path, meshfile, &
    t_max, cfl, &
    init, &
    scheme, method_length, b2d_h, &
    n_bc, bc_name, bc_type, bc_val, &
    sol_uniform, boundary_2d, &
    n_sol_vtu

  call MPI_INIT(mpi_ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

  open(newunit=fn, file="input_data.f")
  read(unit=fn, nml=INPUT_PARAM)
  close(fn)

  call read_mesh_msh(mesh, meshfile_path, meshfile, &
    n_bc, bc_name, me, num_procs, mpi_send_recv)
  call build_mesh(mesh, num_procs, mpi_send_recv, &
    .true., boundary_2d)
  call compute_geometry_mesh(mesh, .true., boundary_2d)

  allocate(vp(3, mesh%n_vert))
  allocate(pp(mesh%n_vert))
  allocate(sol(5, mesh%n_elems))
  allocate(new_sol(5, mesh%n_elems))
  allocate(rhs(5, mesh%n_elems))

  vp = 0.0_DOUBLE
  call init_sol(mesh, sol, sol_uniform, init, me, num_procs)

  i_sol_vtu = 0
  write(fln, *) i_sol_vtu
  write(fln, *) "output_"//trim(adjustl(fln))
  call write_sol_lag(mesh, fln, sol, vp, pp)
  i_sol_vtu = i_sol_vtu + 1

  allocate(mass(mesh%n_elems))
  mass = 0.0_DOUBLE
  do i=1, mesh%n_elems
    mass(i) = mesh%elem(i)%volume/sol(1, i)
  end do

  t = 0.0_DOUBLE
  iter = 1
  do while ( t < t_max ) 

    if( scheme == "classic" ) then
      call compute_rhs_lagrange(mesh, sol, vp, &
        dt, rhs, n_bc, bc_type, bc_val, boundary_2d, mass)
    else  if( scheme == "sidil" ) then
      call compute_rhs_lagrange_sidil(mesh, sol, vp, &
        dt, rhs, n_bc, bc_type, bc_val, boundary_2d, mass, method_length, b2d_h)
    else 
      print*, "No scheme !"
      error stop
    end if

    call compute_dt(mesh, sol, dt, cfl, vp, me, num_procs)
    if( t + dt > t_max ) dt = t_max - t
    call move_mesh(mesh, vp, dt)

    call mpi_memory_exchange_vert(mesh, mpi_send_recv)
    call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol)

    call compute_geometry_mesh(mesh, .true., boundary_2d)
    new_sol = sol + dt * rhs
    !do i=1, mesh%n_elems
      !new_sol(1, i) = mesh%elem(i)%volume/mass(i)
      !mesh%elem(i)%volume = mass(i) * new_sol(1, i)
    !end do
    sol = new_sol
    t = t + dt
    if( mod(iter, 100) == 0 .and. me == 0 ) print*, t, dt

    if( t >= i_sol_vtu * t_max / real(n_sol_vtu - 1) ) then
      write(fln, *) i_sol_vtu
      write(fln, *) "output_"//trim(adjustl(fln))
      call write_sol_lag(mesh, fln, new_sol, vp, pp)
      i_sol_vtu = i_sol_vtu + 1
    end if

    iter = iter + 1
  end do

  i_sol_vtu = -1
  write(fln, *) i_sol_vtu
  write(fln, *) "output_"//trim(adjustl(fln))
  call write_sol_lag(mesh, fln, new_sol, vp, pp)
  call MPI_FINALIZE(mpi_ierr)
end program main
