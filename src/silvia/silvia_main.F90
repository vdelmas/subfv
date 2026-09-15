program main
    use mpi
    use precision_module
    use mpi_module
    use mesh_module
    use mesh_reading_module
    use mesh_geometry_module
    use mesh_connectivity_module
    use io_module
    implicit none

    real(kind=DOUBLE), parameter :: gamma = 1.4_DOUBLE

    integer :: funit
    integer :: mpi_ierr, me, num_procs
    type(mesh_type) :: mesh
    type(mpi_send_recv_type) :: mpi_send_recv

    character(len=255) :: meshfile_path, meshfile
    integer(kind=ENTIER) :: n_bc=0, i, le, re, i_bc
    integer(kind=ENTIER) :: iter
    character(len=255), dimension(10) :: bc_name, bc_type
    logical :: boundary_2d

    real(kind=DOUBLE) :: t, dt, t_max, area
    real(kind=DOUBLE), dimension(3) :: n
    real(kind=DOUBLE), dimension(:,:), allocatable :: sol
    real(kind=DOUBLE), dimension(:,:), allocatable :: rhs
    real(kind=DOUBLE), dimension(5) :: sol_l, sol_r, sol_w_l, sol_w_r

    namelist /INPUT_PARAM/ &
        meshfile_path, meshfile, &
        boundary_2d, &
        n_bc, bc_name, bc_type

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

    ! Same treatment for bc_type -- see bc_type_id's declaration. Same
    ! fallback as ghost_prim's old string select case: anything
    ! unrecognized (including blank) defaults to a slip wall.
    ! do i_bc = 1, n_bc
    !   select case (trim(adjustl(bc_type(i_bc))))
    !   case ('freestream')
    !     bc_type_id(i_bc) = BC_FREESTREAM
    !   case ('outflowsupersonic', 'outflow')
    !     bc_type_id(i_bc) = BC_OUTFLOW
    !   case ('dmr_top')
    !     bc_type_id(i_bc) = BC_DMR_TOP
    !   case default
    !     bc_type_id(i_bc) = BC_WALL
    !   end select
    ! end do

    do i=1, mesh%n_faces
      le = mesh%face(i)%left_neigh
      re = mesh%face(i)%right_neigh
      if(re <= 0) then
        if(re < 0) then
          print*, "boundary", i, -re, trim(adjustl(bc_name(-re))), trim(adjustl(bc_type(-re)))
        end if
      end if
    end do

    allocate(sol(5,mesh%n_elems))
    do i=1, mesh%n_elems
      sol(:, i) = (/1.0_DOUBLE,1.0_DOUBLE,1.0_DOUBLE,1.0_DOUBLE,100.0_DOUBLE/)
    end do

    !call init_sol(mesh, sol)

    t=0.0_DOUBLE
    t_max=1.0_DOUBLE
    iter=0
    call write_vtu(mesh, sol, 0, iter)
    do while (t<t_max)
      dt=0.1_DOUBLE
      !call compute_dt(mesh, sol, dt)
      !call compute_rhs(mesh, sol, rhs)
      do i=1,mesh%n_elems
        !sol(:, i) = sol(:, i) - dt/mesh%elem(i)%volume * rhs(:, i)
        sol(:, i) = sol(:, i) - dt
      end do
      t = t + dt
      iter = iter + 1
      call write_vtu(mesh, sol, 0, iter)
    end do

    call MPI_FINALIZE(mpi_ierr)
contains
subroutine compute_rhs(mesh, sol, rhs)
  implicit none

  type(mesh_type), intent(in) :: mesh
  real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
  real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs

  real(kind=DOUBLE), dimension(5, 2) :: lr_flux

  do i=1, mesh%n_faces
    le = mesh%face(i)%left_neigh
    sol_l = sol(:, le)
    sol_w_l = conserv_to_primit(sol_l)

    re = mesh%face(i)%right_neigh
    if(re > 0) then
      sol_r = sol(:, re)
      sol_w_r = conserv_to_primit(sol_r)
    else
      ! reconstruct left state and then apply BCs
      ! sol_w_r = ghost_prim(xface, norm, wL, -ir, t)
    end if
  
    area = mesh%face(i)%area
    n = mesh%face(i)%norm
    !call two_wave(sol_w_l, sol_w_r, n, lr_flux, sl, sr)

    rhs(:, le) = rhs(:, le) + area*lr_flux(:, 1)
    if(re > 0) then
      rhs(:, re) = rhs(:, re) + area*lr_flux(:, 2)
    end if

  end do
end subroutine compute_rhs

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

pure function primit_to_conserv(w) result(u)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE), dimension(5) :: u

    u(1) = w(1)
    u(2:4) = w(2:4)*w(1)
    u(5) = w(5)/(gamma - 1) &
      + 0.5_DOUBLE*w(1)*(w(2)**2 + w(3)**2 + w(4)**2)
  end function primit_to_conserv

  pure function conserv_to_primit(u) result(w)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: u
    real(kind=DOUBLE), dimension(5) :: w

    w(1) = u(1)
    w(2:4) = u(2:4)/u(1)
    w(5) = (gamma - 1)*(u(5) &
      - 0.5_DOUBLE*w(1)*(w(2)**2 + w(3)**2 + w(4)**2))
  end function conserv_to_primit

  ! Ghost cell primitive state for boundary condition
  ! pure function ghost_prim(xf, norm, wL, id_bc, t) result(wR)
  !   real(kind=DOUBLE), dimension(3), intent(in) :: xf, norm
  !   real(kind=DOUBLE), dimension(5), intent(in) :: wL
  !   integer(kind=ENTIER), intent(in) :: id_bc
  !   real(kind=DOUBLE), intent(in) :: t
  !   real(kind=DOUBLE), dimension(5) :: wR

  !   real(kind=DOUBLE) :: vn

  !   if (id_bc < 1 .or. id_bc > n_bc) then
  !     ! Default: slip wall (mirror normal velocity)
  !     wR    = wL
  !     vn    = dot_product(wL(2:4), norm)
  !     wR(2) = wL(2) - 2.0_DOUBLE * vn * norm(1)
  !     wR(3) = wL(3) - 2.0_DOUBLE * vn * norm(2)
  !     wR(4) = wL(4) - 2.0_DOUBLE * vn * norm(3)
  !     return
  !   end if

  !   ! bc_type_id resolved once in read_params -- see its declaration.
  !   select case (bc_type_id(id_bc))
  !   case (BC_FREESTREAM)
  !     wR = bc_val(:, id_bc)
  !   case (BC_OUTFLOW)
  !     ! Zero-gradient extrapolation: valid where every characteristic
  !     ! leaves the domain (locally supersonic outflow).
  !     wR = wL
  !   case (BC_DMR_TOP)
  !     ! Double Mach reflection: exact post/pre-shock state at (xf(1), t),
  !     ! following the shock's known straight-line motion. See dmr_state.
  !     wR = dmr_state(xf(1), xf(2), t)
  !   case default
  !     ! BC_WALL (also the fallback for blank/unrecognized bc_type,
  !     ! matching the old string select case's behavior)
  !     wR    = wL
  !     vn    = dot_product(wL(2:4), norm)
  !     wR(2) = wL(2) - 2.0_DOUBLE * vn * norm(1)
  !     wR(3) = wL(3) - 2.0_DOUBLE * vn * norm(2)
  !     wR(4) = wL(4) - 2.0_DOUBLE * vn * norm(3)
  !   end select
  ! end function ghost_prim

  subroutine two_wave(sol_w_l, sol_w_r, n, lr_flux, sl, sr)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), dimension(5) :: sol_l, sol_r
    real(kind=DOUBLE), dimension(5, 2), intent(inout) :: lr_flux
    real(kind=DOUBLE), intent(inout) :: sl, sr

    real(kind=DOUBLE) :: rhol, rhor, vn_l, vn_r, pl, pr, &
      al, ar, lambda_l, lambda_r
    real(kind=DOUBLE), dimension(5) :: fl, fr, sol_et

    rhol = sol_w_l(1)
    vn_l = dot_product(sol_w_l(2:4), n)
    pl = sol_w_l(5)
    al = sound_speed_w(sol_w_l)
    sol_l = primit_to_conserv(sol_w_l)

    rhor = sol_w_r(1)
    vn_r = dot_product(sol_w_r(2:4), n)
    pr = sol_w_r(5)
    ar = sound_speed_w(sol_w_r)
    sol_r = primit_to_conserv(sol_w_r)

    fl(1)   = vn_l*sol_l(1)
    fl(2:4) = vn_l*sol_l(2:4) + pl*n
    fl(5)   = (sol_l(5) + pl)*vn_l

    fr(1)   = vn_r*sol_r(1)
    fr(2:4) = vn_r*sol_r(2:4) + pr*n
    fr(5)   = (sol_r(5) + pr)*vn_r

    lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(vn_r - vn_l))
    lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(vn_r - vn_l))

    sl = vn_l - lambda_l/rhol
    sr = vn_r + lambda_r/rhor

    sol_et = (sr*sol_r - sl*sol_l - (fr - fl))/(sr - sl)

    lr_flux(:, 1) = 0.5_DOUBLE*(fl + fr) &
      - 0.5_DOUBLE*(abs(sl)*(sol_et - sol_l) + abs(sr)*(sol_r - sol_et))
    lr_flux(:, 2) = -lr_flux(:, 1)
  end subroutine two_wave

  pure function sound_speed_w(w) result(a)
    implicit none

    real(kind=DOUBLE), parameter :: gamma = 1.4_DOUBLE
    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE) :: a

    a = sqrt(gamma*w(5)/w(1))
  end function sound_speed_w
end program main


