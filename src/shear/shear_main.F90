program main_shear
  use mpi
  use precision_module
  use mpi_module
  use mesh_module
  use mesh_geometry_module
  use mesh_reading_module
  use mesh_connectivity_module
  use ns_global_data_module
  use ns_euler_module
  use ns_mesh_metric_module
  use ns_io_module

  implicit none

  integer(kind=ENTIER) :: me, num_procs, mpi_ierr, i
  integer(kind=ENTIER), parameter :: nangle=180
  integer(kind=ENTIER) :: iangle, fn_error
  real(kind=DOUBLE) :: error, error_global, angle
  real(kind=DOUBLE), dimension(5) :: unp, wnp, wn
  real(kind=DOUBLE), parameter :: dt_shear = 1e-5

  type(mesh_type) :: mesh
  type(mpi_send_recv_type) :: mpi_send_recv

  real(kind=DOUBLE), dimension(:, :),    allocatable :: sol, rhs
  real(kind=DOUBLE), dimension(:, :, :), allocatable :: cell_grad, mat_h_p
  real(kind=DOUBLE), dimension(:, :),    allocatable :: vp
  real(kind=DOUBLE), dimension(:),       allocatable :: sum_lambda, h_p

  call MPI_INIT(mpi_ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

  call read_input_parameters("input_data.f")
  call init_flags()

  call read_mesh_msh(mesh, meshfile_path, meshfile, &
    n_bc, bc_name, me, num_procs, mpi_send_recv)
  call build_mesh(mesh, num_procs, mpi_send_recv, .TRUE., boundary_2d)
  call compute_geometry_mesh(mesh, .TRUE., boundary_2d, me, num_procs)

  allocate(sol(5, mesh%n_elems))          ; sol        = 0.0_DOUBLE
  allocate(rhs(5, mesh%n_elems))          ; rhs        = 0.0_DOUBLE
  allocate(sum_lambda(mesh%n_elems))      ; sum_lambda  = 0.0_DOUBLE
  allocate(cell_grad(3, 5, mesh%n_elems)) ; cell_grad   = 0.0_DOUBLE
  allocate(vp(3, mesh%n_vert))           ; vp          = 0.0_DOUBLE
  allocate(h_p(mesh%n_vert))              ; h_p         = 0.0_DOUBLE
  allocate(mat_h_p(3, 3, mesh%n_vert))   ; mat_h_p     = 0.0_DOUBLE

  if (scheme_id == SCHEME_ZB) then
    do i = 1, mesh%n_vert
      h_p(i)         = compute_length(mesh, i)
      mat_h_p(:,:,i) = compute_ellip(mesh, i)
    end do
  end if

  if (me == 0) open(newunit=fn_error, file="error_shear.dat", status="replace")

  do iangle=1, nangle+1
    angle = real(iangle-1)/real(nangle)*360.0_DOUBLE
    call init_shear(mesh, sol, angle)
    call compute_euler_flux(mesh, sol, rhs, cell_grad, sum_lambda, &
      .false., vp, h_p, mat_h_p)

    !call write_sol(mesh, sol, cell_grad, vp, h_p, me, iangle)

    error = 0.0_DOUBLE
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        if ( abs(mesh%elem(i)%coord(1)) < 0.8 .and. &
          abs(mesh%elem(i)%coord(2)) < 0.8 ) then
          unp = sol(:, i) + dt_shear/mesh%elem(i)%volume*rhs(:, i)
          wnp = conserv_to_primit(unp)
          wn = conserv_to_primit(sol(:, i))
          error = error + mesh%elem(i)%volume*sum(abs(wnp(2:4)-wn(2:4)))
          sol(:, i) = sol(:, i) + dt_shear/mesh%elem(i)%volume*rhs(:, i)
        end if
      end if
    end do

    if (num_procs > 1) then
      call MPI_ALLREDUCE(error, error_global, 1, MPI_DOUBLE_PRECISION, &
        MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
    else
      error_global = error
    end if
    if (me == 0) write(fn_error, *) angle, error_global
  end do

  if (me == 0) close(fn_error)

  call MPI_FINALIZE(mpi_ierr)
contains
  subroutine init_shear(mesh, sol, angle)
    use ns_euler_primitives_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: sol
    real(kind=DOUBLE) :: angle

    integer(kind=ENTIER) :: i
    real(kind=DOUBLE) :: angle_rad
    real(kind=DOUBLE), dimension(5) :: w

    angle_rad = 2*pi*angle/360.0_DOUBLE
    do i=1, mesh%n_elems
      if (cos(angle_rad)*mesh%elem(i)%coord(2) &
        - sin(angle_rad)*mesh%elem(i)%coord(1) > 0.0_DOUBLE) then
        w(1) = 1.0_DOUBLE
        w(2) = cos(angle_rad)
        w(3) = sin(angle_rad)
        w(4) = 0.0_DOUBLE
        w(5) = 1.0_DOUBLE
      else
        w(1) = 1.0_DOUBLE
        w(2) = -cos(angle_rad)
        w(3) = -sin(angle_rad)
        w(4) = 0.0_DOUBLE
        w(5) = 1.0_DOUBLE
      end if
      sol(:, i) = primit_to_conserv(w)
    end do
  end subroutine init_shear
end program main_shear
