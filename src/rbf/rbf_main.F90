program main
  use mpi
  use precision_module
  use mpi_module
  use mesh_module
  use mesh_reading_module
  use mesh_geometry_module
  use mesh_connectivity_module
  use io_module
  use rbf_module

  implicit none

  integer(kind=ENTIER) :: n_bc=0, me, num_procs, fn, mpi_ierr
  character(len=255) :: meshfile_path, meshfile
  character(len=255), dimension(10) :: bc_name
  type(mesh_type) :: mesh
  type(mpi_send_recv_type) :: mpi_send_recv
  logical :: b2d

  real(kind=DOUBLE) :: tmax=2.d0, cpu_t1, cpu_t2
  integer(kind=ENTIER) :: id_surf_imp=1
  real(kind=DOUBLE) :: rad=1

  namelist /params/ meshfile, tmax, id_surf_imp, rad, n_bc, bc_name

  integer(kind=ENTIER) :: fn_vtu, fn_pvtu

  integer(kind=ENTIER) :: nimp, nslide, ndof, n, i, j, k, id_sub_face
  integer(kind=ENTIER), dimension(:), allocatable :: imp_vert
  real(kind=DOUBLE), dimension(:,:), allocatable :: x, no_arr, disp
  real(kind=DOUBLE), dimension(:), allocatable :: delta_imp, delta
  real(kind=DOUBLE), dimension(3) :: no

  real(kind=DOUBLE) :: t, dt
  integer(kind=ENTIER) :: iter

  call MPI_INIT(mpi_ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

  meshfile_path=""
  b2d=.true.

  open(newunit=fn, file="input_data.f", status="old")
  read(fn, nml=params)
  close(fn)

  call read_mesh_msh(mesh, meshfile_path, meshfile, &
    n_bc, bc_name, me, num_procs, mpi_send_recv)
  call build_mesh(mesh, num_procs, mpi_send_recv, .true., b2d)
  call compute_geometry_mesh(mesh, .true., b2d)

  ! RBF point set = two separate lists, concatenated as 1:nimp then nimp+1:n:
  !  - hard points (1:nimp): the imposed surface (id_surf_imp) vertices,
  !    delta_imp = imposed_vel*dt, full 3d dof each
  !  - sliding points (nimp+1:nimp+nslide): the outer far-field boundary
  !    vertices (is_bound, far from the origin), delta_imp=0, 1 dof each
  !    along the local normal (no_arr) -- so the far boundary can't be
  !    dragged in/out, and is left completely free tangentially (not
  !    constrained/solved for at all).
  ! delta_imp/delta are packed to their true size ndof=3*nimp+nslide (see
  ! dof_offset in rbf_module): 3 slots per hard point, 1 per sliding point.
  ! Membership is a mesh-topology property, fixed once; only x (and no_arr,
  ! since it depends on the surface normal at each point) change per
  ! iteration as the mesh moves, so imp_vert is built once here and reused
  ! every iteration with the *same* point ordering.
  nimp = 0
  nslide = 0
  do i=1, mesh%n_vert
    if( vert_belongs_to_surf(mesh, i, id_surf_imp) ) then
      nimp = nimp + 1
    else if( mesh%vert(i)%is_bound .and. norm2(mesh%vert(i)%coord) > 10.0_DOUBLE ) then
      nslide = nslide + 1
    end if
  end do
  n = nimp + nslide
  ndof = 3*nimp + nslide

  print*, mesh%n_vert, nimp, nslide

  allocate(imp_vert(n))
  allocate(x(3, n))
  allocate(no_arr(3, n))
  allocate(delta_imp(ndof))
  allocate(delta(ndof))
  allocate(disp(3, mesh%n_vert))
  delta_imp = 0.0_DOUBLE
  delta = 0.0_DOUBLE
  disp = 0.0_DOUBLE
  no_arr = 0.0_DOUBLE

  k = 0
  do i=1, mesh%n_vert
    if( vert_belongs_to_surf(mesh, i, id_surf_imp) ) then
      k = k + 1
      imp_vert(k) = i
      x(:, k) = mesh%vert(i)%coord
      delta_imp(dof_offset(k, nimp)+1:dof_offset(k, nimp)+3) = (/1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE/)
    end if
  end do
  do i=1, mesh%n_vert
    if( mesh%vert(i)%is_bound .and. norm2(mesh%vert(i)%coord) > 10.0_DOUBLE ) then
      k = k + 1
      imp_vert(k) = i
      x(:, k) = mesh%vert(i)%coord
      do j=1, mesh%vert(i)%n_sub_faces_neigh
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        if( mesh%sub_face(id_sub_face)%right_elem_neigh <= 0 ) then
          if( abs(mesh%sub_face(id_sub_face)%norm(3)) < 1e-8_DOUBLE ) then
            no = mesh%sub_face(id_sub_face)%norm
          end if
        end if
      end do
      no_arr(:, k) = no
      delta_imp(dof_offset(k, nimp)+1) = 0.0_DOUBLE
    end if
  end do

  t=0.d0
  dt=0.1d0

  id_surf_imp = 1
  iter = 0

  if( nimp < 1 ) then
    print*, "bad imp", id_surf_imp
  end if

  call cpu_time(cpu_t1)

  call compute_rbf_field(nimp, nslide, x, no_arr, delta_imp, rad, delta)
  do i=1, mesh%n_vert
    call eval_rbf_field(n, nimp, x, no_arr, rad, delta, mesh%vert(i)%coord, disp(:, i))
  end do

  call cpu_time(cpu_t2)
  print*, cpu_t2-cpu_t1

  call open_file_vtu(mesh, "test", fn_vtu, fn_pvtu)
  call write_file_vtu_start_vert_data(mesh, "test", fn_vtu, fn_pvtu)
  call write_file_vtu_vert_vector(mesh, "test", fn_vtu, fn_pvtu, disp, "dispv")
  call write_file_vtu_end_vert_data(mesh, "test", fn_vtu, fn_pvtu)
  call close_file_vtu(mesh, "test", fn_vtu, fn_pvtu)

  call MPI_FINALIZE(mpi_ierr)

contains
  function vert_belongs_to_surf(mesh, id_vert, id_surf)
    implicit none

    type(mesh_type) :: mesh
    integer(kind=4) :: id_vert, id_surf
    integer(kind=4) :: i, id_face
    logical :: vert_belongs_to_surf

    vert_belongs_to_surf = .false.
    do i=1, mesh%vert(id_vert)%n_faces_neigh
      id_face = mesh%vert(id_vert)%face_neigh(i)
      if( mesh%face(id_face)%right_neigh == - id_surf ) then
        vert_belongs_to_surf = .true.
        return
      end if
    end do
  end function vert_belongs_to_surf
end program main
