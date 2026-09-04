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
  integer(kind=ENTIER) :: init = 0
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

  allocate(vp(3, mesh%n_vert))
  allocate(pp(mesh%n_vert))
  allocate(sol(5, mesh%n_elems))
  allocate(new_sol(5, mesh%n_elems))
  allocate(rhs(5, mesh%n_elems))
  allocate(gamma_arr(mesh%n_elems))
  allocate(vp_is_imposed(mesh%n_vert))

  if (init == 5) then
    gamma_arr = 5.0_DOUBLE / 3.0_DOUBLE
  else
    gamma_arr = 1.4_DOUBLE
  end if
  vp = 0.0_DOUBLE
  call init_sol(mesh, sol, sol_uniform, init, me, num_procs, gamma_arr)
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

  i_sol_vtu = 0
  write(fln, *) i_sol_vtu
  write(fln, *) "output_"//trim(adjustl(fln))
  call write_sol_lag(mesh, fln, sol, vp, pp, gamma_arr)
  call write_sol_dat_lag(mesh, fln, sol, gamma_arr)
  i_sol_vtu = i_sol_vtu + 1

  allocate(mass(mesh%n_elems))
  mass = 0.0_DOUBLE
  do i=1, mesh%n_elems
    mass(i) = mesh%elem(i)%volume/sol(1, i)
  end do

  t = 0.0_DOUBLE
  iter = 1
  do while ( t < t_max )

    ! Set imposed velocities (piston BC) before the RHS computation
    if (piston_bc_idx > 0) then
      do i = 1, mesh%n_vert
        if (vp_is_imposed(i)) vp(:, i) = piston_vel
      end do
    end if

    if( scheme == "classic" ) then
      call compute_rhs_lagrange(mesh, sol, vp, &
        dt, rhs, n_bc, bc_type, bc_val, boundary_2d, mass, gamma_arr, vp_is_imposed)
    else  if( scheme == "sidil" ) then
      call compute_rhs_lagrange_sidil(mesh, sol, vp, &
        dt, rhs, n_bc, bc_type, bc_val, boundary_2d, mass, method_length, b2d_h, &
        gamma_arr, vp_is_imposed)
    else
      print*, "No scheme !"
      error stop
    end if

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

    if( t >= i_sol_vtu * t_max / real(n_sol_vtu - 1) ) then
      write(fln, *) i_sol_vtu
      write(fln, *) "output_"//trim(adjustl(fln))
      call write_sol_lag(mesh, fln, new_sol, vp, pp, gamma_arr)
      call write_sol_dat_lag(mesh, fln, new_sol, gamma_arr)
      i_sol_vtu = i_sol_vtu + 1
    end if

    iter = iter + 1
  end do

  i_sol_vtu = -1
  write(fln, *) i_sol_vtu
  write(fln, *) "output_"//trim(adjustl(fln))
  call write_sol_lag(mesh, fln, new_sol, vp, pp, gamma_arr)
  call write_sol_dat_lag(mesh, fln, new_sol, gamma_arr)
  call MPI_FINALIZE(mpi_ierr)
end program main
