module lagrange_module
  use mpi
  use precision_module
  use mpi_module
  use mesh_module
  use mesh_reading_module
  use mesh_geometry_module
  use mesh_connectivity_module

  implicit none

  real(kind=DOUBLE), parameter :: gamma = 7.0_DOUBLE/5.0_DOUBLE
contains
  subroutine compute_rhs_lagrange(mesh, sol, vp, dt, rhs, n_bc, bc_type, bc_val, b2d, mass)
    use linear_solver_module, only: lu_solve, print_mat, inverse_3_by_3
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), intent(in) :: dt
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mass
    integer(kind=ENTIER), intent(in) :: n_bc
    character(len=255), dimension(n_bc) :: bc_type
    real(kind=DOUBLE), dimension(5, n_bc) :: bc_val
    logical, intent(in) :: b2d

    integer(kind=ENTIER) :: i, j, id_elem, id_face, id_sub_face
    integer(kind=ENTIER) :: nsfn, idl, idr, var
    real(kind=DOUBLE), dimension(:, :), allocatable :: lambda
    real(kind=DOUBLE), dimension(:, :, :), allocatable :: sol_lr
    real(kind=DOUBLE), dimension(5) :: flux
    real(kind=DOUBLE), dimension(3, 3) :: Mp, Mp_inv
    real(kind=DOUBLE), dimension(3) :: Rp, Bp
    real(kind=DOUBLE) :: vl, vr, pl, pr, v_bar, pl_et, pr_et, v_et, Pp
    logical :: corner

    rhs = 0.0_DOUBLE
    vp = 0.0_DOUBLE
    do i=1, mesh%n_vert
      nsfn = mesh%vert(i)%n_sub_faces_neigh
      allocate(lambda(2, nsfn))
      lambda = 0.0_DOUBLE
      allocate(sol_lr(5, 2, nsfn))
      sol_lr = 0.0_DOUBLE

      !bluid nodal system
      Bp = 0.0_DOUBLE
      Mp = 0.0_DOUBLE
      Rp = 0.0_DOUBLE
      lambda(:, :) = 0.0_DOUBLE
      do j=1, nsfn
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        idl = mesh%face(id_face)%left_neigh
        idr = mesh%face(id_face)%right_neigh
        sol_lr(:, 1, j) = sol(:, idl)
        pl = pressure(sol_lr(:, 1, j))
        vl = dot_product(sol_lr(2:4, 1, j), mesh%face(id_face)%norm)


        if( idr > 0 ) then
          sol_lr(:, 2, j) = sol(:, idr)
          pr = pressure(sol_lr(:, 2, j))
          vr = dot_product(sol_lr(2:4, 2, j), mesh%face(id_face)%norm)

          lambda(1, j) = max(sqrt(gamma*pl*sol_lr(1, 1, j))/sol_lr(1, 1, j), &
            sqrt(max(0.0_DOUBLE, pr-pl)/sol_lr(1, 1, j)), &
            -(vr-vl)/sol_lr(1, 1, j))
          lambda(2, j) = max(sqrt(gamma*pr*sol_lr(1, 2, j))/sol_lr(1, 2, j), &
            sqrt(max(0.0_DOUBLE, pl-pr)/sol_lr(1, 2, j)), &
            -(vr-vl)/sol_lr(1, 2, j))

          v_bar = (lambda(1, j)*vl + lambda(2, j)*vr - (pr - pl))&
            /(lambda(2, j) + lambda(1, j))
          Mp = Mp + mesh%sub_face(id_sub_face)%area*(lambda(1,j)+lambda(2,j)) &
            * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
          Rp = Rp + mesh%sub_face(id_sub_face)%area*(lambda(1,j)+lambda(2,j)) &
            *mesh%face(id_face)%norm*v_bar

        else
          if( b2d .and. abs(mesh%face(id_face)%norm(3)) > 1e-8_DOUBLE ) then
            lambda(1, j) = sqrt(gamma*pl*sol_lr(1, 1, j))/sol_lr(1, 1, j)
            Mp = Mp + mesh%sub_face(id_sub_face)%area*2*lambda(1,j) &
              * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
          else
            lambda(1, j) = sqrt(gamma*pl*sol_lr(1, 1, j))/sol_lr(1, 1, j)
            Bp = Bp + mesh%sub_face(id_sub_face)%area*mesh%face(id_face)%norm
            Mp = Mp + mesh%sub_face(id_sub_face)%area*lambda(1,j) &
              * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
            Rp = Rp + mesh%sub_face(id_sub_face)%area*(pl + lambda(1,j)*vl)&
              *mesh%face(id_face)%norm
          end if
        end if
      end do

      if( maxval(abs(Bp)) > 1e-8_DOUBLE ) then
        call inverse_3_by_3(Mp, Mp_inv)
        Pp = dot_product(Rp, matmul(Mp_inv,Bp))/dot_product(Bp, matmul(Mp_inv, Bp))
        vp(:, i) = matmul(Mp_inv, Rp-Pp*BP)
      else
        !call lu_solve(3, Mp, vp(:, i), Rp)
        call inverse_3_by_3(Mp, Mp_inv)
        vp(:, i) = matmul(Mp_inv, Rp)
      end if

      !if( b2d ) vp(3, i) = 0.0_DOUBLE
      !call lu_solve(3, Mp, vp(:, i), Rp)

      do j=1, nsfn
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        idl = mesh%face(id_face)%left_neigh
        idr = mesh%face(id_face)%right_neigh

        v_et = dot_product(vp(:, i), mesh%face(id_face)%norm)
        vl = dot_product(sol_lr(2:4, 1, j), mesh%face(id_face)%norm)
        pl = pressure(sol_lr(:, 1, j))
        pl_et = pl - lambda(1, j)*(v_et - vl)

        flux(1)   = -v_et
        flux(2:4) = pl_et*mesh%face(id_face)%norm
        flux(5)   = pl_et*v_et

        do var=1, 5
          rhs(var, idl) = rhs(var, idl) &
            - mesh%sub_face(id_sub_face)%area/mass(idl) * flux(var)
        end do

        if( idr > 0 ) then
          vr = dot_product(sol_lr(2:4, 2, j), mesh%face(id_face)%norm)
          pr = pressure(sol_lr(:, 2, j))
          pr_et = pr + lambda(2, j)*(v_et - vr)
          flux(1)   = -v_et
          flux(2:4) = pr_et*mesh%face(id_face)%norm
          flux(5)   = pr_et*v_et
          do var=1, 5
            rhs(var, idr) = rhs(var, idr) &
              + mesh%sub_face(id_sub_face)%area/mass(idr) * flux(var)
          end do
        end if
      end do

      deallocate(lambda, sol_lr)
    end do
  end subroutine compute_rhs_lagrange

  subroutine compute_rhs_lagrange_sidil(mesh, sol, vp, dt, rhs, &
      n_bc, bc_type, bc_val, b2d, mass, method, b2d_h)
    use linear_solver_module, only: lu_solve, print_mat, inverse_3_by_3
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), intent(in) :: dt
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mass
    integer(kind=ENTIER), intent(in) :: n_bc, method
    character(len=255), dimension(n_bc) :: bc_type
    real(kind=DOUBLE), dimension(5, n_bc) :: bc_val
    logical, intent(in) :: b2d
    real(kind=DOUBLE), intent(in) :: b2d_h

    integer(kind=ENTIER) :: i, j, id_elem, id_face, id_sub_face
    integer(kind=ENTIER) :: id_sub_elem
    integer(kind=ENTIER) :: nsfn, idl, idr, var
    real(kind=DOUBLE), dimension(:, :), allocatable :: lambda
    real(kind=DOUBLE), dimension(:, :, :), allocatable :: sol_lr
    real(kind=DOUBLE), dimension(5) :: flux, sol_w
    real(kind=DOUBLE), dimension(3, 3) :: Mp, Mp_inv
    real(kind=DOUBLE), dimension(3) :: Rp, Bp
    real(kind=DOUBLE) :: vl, vr, pl, pr, v_bar, pl_et, pr_et, v_et, pp, a_p
    logical :: corner

    rhs = 0.0_DOUBLE
    vp = 0.0_DOUBLE
    do i=1, mesh%n_vert

      call compute_nodal_velocity_sidil(mesh, i, sol, vp(:, i), method, b2d, b2d_h)
      call compute_nodal_pressure_sidil(mesh, i, sol, pp, method, b2d, b2d_h)

      a_p = 0.0_DOUBLE
      do j=1, mesh%vert(i)%n_sub_elems_neigh
        id_sub_elem = mesh%vert(i)%sub_elem_neigh(j)
        id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
        sol_w = lag_to_primit(sol(:, id_elem))
        a_p = a_p + mesh%sub_elem(id_sub_elem)%volume &
          * sqrt(gamma*sol_w(5)/sol_w(1))
      end do
      a_p = a_p / mesh%vert(i)%volume

      nsfn = mesh%vert(i)%n_sub_faces_neigh
      do j=1, nsfn
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        idl = mesh%face(id_face)%left_neigh
        idr = mesh%face(id_face)%right_neigh

        v_et = dot_product(vp(:, i), mesh%face(id_face)%norm)

        flux(1)   = -v_et
        flux(2:4) = pp*mesh%face(id_face)%norm
        flux(5)   = pp*v_et
        do var=1, 5
          rhs(var, idl) = rhs(var, idl) &
            - mesh%sub_face(id_sub_face)%area/mass(idl) * flux(var)
        end do

        if( idr > 0 ) then
          flux(1)   = -v_et
          flux(2:4) = pp*mesh%face(id_face)%norm
          flux(5)   = pp*v_et
          do var=1, 5
            rhs(var, idr) = rhs(var, idr) &
              + mesh%sub_face(id_sub_face)%area/mass(idr) * flux(var)
          end do
        end if
      end do

    end do
  end subroutine compute_rhs_lagrange_sidil

  function tensor_product(a,b) result(c)
    implicit none

    real(kind=DOUBLE), dimension(:) :: a, b
    real(kind=DOUBLE), dimension(size(a), size(b)) :: c

    integer(kind=ENTIER) :: i, j

    do i=1, size(a)
      do j=1, size(b)
        c(i, j) = a(i)*b(j)
      end do
    end do
  end function tensor_product

  subroutine compute_dt(mesh, sol, dt, cfl, vp, me, num_procs)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: vp
    real(kind=DOUBLE), intent(in) :: cfl
    real(kind=DOUBLE), intent(inout) :: dt
    integer(kind=ENTIER) :: me, num_procs

    real(kind=DOUBLE), parameter :: CV = 0.8_DOUBLE
    integer(kind=ENTIER) :: i, j, k, mpi_ierr
    integer(kind=ENTIER) :: id_sub_face, id_face
    integer(kind=ENTIER) :: id_vert, id_sub_elem
    real(kind=DOUBLE) :: dx, c
    real(kind=DOUBLE), dimension(5) :: w
    real(kind=DOUBLE), dimension(3) :: norm

    dt = 1e6
    do i=1, mesh%n_elems
      c = sqrt(gamma*pressure(sol(:, i))*sol(1, i))
      do j=1, mesh%elem(i)%n_sub_elems
        id_sub_elem = mesh%elem(i)%sub_elem(j)
        do k=1, mesh%sub_elem(id_sub_elem)%n_sub_faces
          id_sub_face = mesh%sub_elem(id_sub_elem)%sub_face(k)
          id_face = mesh%sub_face(id_sub_face)%mesh_face
          id_vert = mesh%sub_face(id_sub_face)%mesh_vert
          if( mesh%face(id_face)%left_neigh == i ) then
            norm = mesh%face(id_face)%norm
          else
            norm = -mesh%face(id_face)%norm
          end if
          dt = min(dt, mesh%sub_elem(i)%volume/(mesh%sub_face(id_sub_face)%area*c))
          dt = min(dt, CV*mesh%sub_elem(id_sub_elem)%volume&
            /(1e-8_DOUBLE + abs(dot_product(vp(:, id_vert), norm))&
            *mesh%face(id_face)%area))
        end do
      end do
    end do
    dt = cfl * dt

    if( num_procs > 1 ) call MPI_ALLREDUCE(MPI_IN_PLACE, &
      dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, mpi_ierr)
  end subroutine compute_dt

  subroutine init_sol(mesh, sol, sol_uniform, init, me, num_procs)
    use mpi
    use mpi_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: sol
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_uniform
    integer(kind=ENTIER), intent(in) :: init
    integer(kind=ENTIER), intent(in) :: me, num_procs

    integer(kind=ENTIER) :: i, mpi_ierr
    real(kind=DOUBLE) :: tot_vol
    !real(kind=DOUBLE), parameter :: r_sedov = 2.4_DOUBLE/20.0_DOUBLE
    real(kind=DOUBLE), parameter :: r_sedov = 0.06_DOUBLE
    real(kind=DOUBLE), parameter :: r_sedov_3d = 0.12_DOUBLE

    if ( init == 0 ) then
      do i=1, mesh%n_elems
        if( mesh%elem(i)%coord(1) < 0.5_DOUBLE ) then
          sol(1, i) = 1.0_DOUBLE
          sol(2:4, i) = 0.0_DOUBLE
          sol(5, i) = 1.0_DOUBLE
        else
          sol(1, i) = 0.125_DOUBLE
          sol(2:4, i) = 0.0_DOUBLE
          sol(5, i) = 0.1_DOUBLE
        end if
        !sol(:, i) = sol_uniform
        sol(1, i) = 1.0_DOUBLE/sol(1, i)
        sol(5, i) = sol(5, i)/((gamma - 1.0_DOUBLE)/sol(1, i)) &
          + 0.5_DOUBLE*norm2(sol(2:4, i))**2
      end do
    else if( init == 1) then

      tot_vol = 0.0_DOUBLE
      do i = 1, mesh%n_elems
        if( .not. mesh%elem(i)%is_ghost ) then
          if (sqrt(mesh%elem(i)%coord(1)**2 + mesh%elem(i)%coord(2)**2) < r_sedov) then
          !if (norm2(mesh%elem(i)%coord) < r_sedov) then
            tot_vol = tot_vol + mesh%elem(i)%volume
          end if
        end if
      end do

      if( num_procs > 1 ) call MPI_ALLREDUCE(MPI_IN_PLACE, &
        tot_vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)

      if (tot_vol < 1e-10_DOUBLE) then
        print *, "[-] r_sedov too small, tot_vol_init:", tot_vol
        error stop
      end if

      do i = 1, mesh%n_elems
        !if (norm2(mesh%elem(i)%coord) < r_sedov) then
        if (sqrt(mesh%elem(i)%coord(1)**2 + mesh%elem(i)%coord(2)**2) < r_sedov) then
          sol(1, i) = 1.0_DOUBLE
          sol(2, i) = 0.0_DOUBLE
          sol(3, i) = 0.0_DOUBLE
          sol(4, i) = 0.0_DOUBLE
          !sol(5,i) = 0.311357_DOUBLE / tot_vol
          !sol(5, i) = 0.244816_DOUBLE/tot_vol
          !sol(5, i) = 0.851072_DOUBLE/tot_vol
          sol(5, i) = 0.983909_DOUBLE/tot_vol
        else
          sol(1, i) = 1.0_DOUBLE
          sol(2, i) = 0.0_DOUBLE
          sol(3, i) = 0.0_DOUBLE
          sol(4, i) = 0.0_DOUBLE
          sol(5, i) = 1e-8_DOUBLE/(gamma - 1.0_DOUBLE)
        end if
      end do
    else if( init == 13) then

      tot_vol = 0.0_DOUBLE
      do i = 1, mesh%n_elems
        if( .not. mesh%elem(i)%is_ghost ) then
          if (norm2(mesh%elem(i)%coord) < r_sedov_3d) then
            tot_vol = tot_vol + mesh%elem(i)%volume
          end if
        end if
      end do

      if( num_procs > 1 ) call MPI_ALLREDUCE(MPI_IN_PLACE, &
        tot_vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)

      if (tot_vol < 1e-10_DOUBLE) then
        print *, "[-] r_sedov_3d too small, tot_vol_init:", tot_vol
        error stop
      end if

      do i = 1, mesh%n_elems
        if (norm2(mesh%elem(i)%coord) < r_sedov_3d) then
          sol(1, i) = 1.0_DOUBLE
          sol(2, i) = 0.0_DOUBLE
          sol(3, i) = 0.0_DOUBLE
          sol(4, i) = 0.0_DOUBLE
          !sol(5,i) = 0.311357_DOUBLE / tot_vol
          !sol(5, i) = 0.244816_DOUBLE/tot_vol
          !sol(5, i) = 0.851072_DOUBLE/tot_vol
          sol(5, i) = 0.983909_DOUBLE/tot_vol
        else
          sol(1, i) = 1.0_DOUBLE
          sol(2, i) = 0.0_DOUBLE
          sol(3, i) = 0.0_DOUBLE
          sol(4, i) = 0.0_DOUBLE
          sol(5, i) = 1e-8_DOUBLE/(gamma - 1.0_DOUBLE)
        end if
      end do
    else if( init == 2) then
      do i=1, mesh%n_elems
        call sol_isentropic_vortex(mesh%elem(i)%coord, sol(:, i), 0.0_DOUBLE)
        sol(5, i) = sol(5, i)/(sol(1, i)*(gamma-1.0_DOUBLE)) &
          + 0.5_DOUBLE*norm2(sol(2:4, i))**2
        sol(1, i) = 1.0_DOUBLE/sol(1, i)
      end do
    else if( init == 42 ) then
      do i=1, mesh%n_elems
        sol(1, i) = 1.0_DOUBLE
        !sol(2:4, i) = 0.0_DOUBLE
        sol(2:4, i) = 0.0_DOUBLE
        sol(2, i) = mesh%elem(i)%coord(1)
        !sol(5, i) = (10+mesh%elem(i)%coord(1))/(gamma-1.0_DOUBLE)
        sol(5, i) = 10/(gamma-1.0_DOUBLE)
      end do
    else if( init == 4 ) then
      do i=1, mesh%n_elems
        if( norm2(mesh%elem(i)%coord(:2)) < 0.5_DOUBLE ) then
          sol(1, i) = 1.0_DOUBLE
          sol(2:4, i) = 0.0_DOUBLE
          sol(5, i) = 1.0_DOUBLE
        else
          sol(1, i) = 0.125_DOUBLE
          sol(2:4, i) = 0.0_DOUBLE
          sol(5, i) = 0.1_DOUBLE
        end if
        sol(1, i) = 1.0_DOUBLE/sol(1, i)
        sol(5, i) = sol(5, i)/((gamma - 1.0_DOUBLE)/sol(1, i)) &
          + 0.5_DOUBLE*norm2(sol(2:4, i))**2
      end do
    end if

  end subroutine init_sol

  function pressure(u)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: u
    real(kind=DOUBLE) :: pressure

    pressure = (gamma-1.0_DOUBLE)/u(1)*(u(5) - 0.5_DOUBLE*norm2(u(2:4))**2)
  end function pressure

  subroutine mpi_memory_exchange_vert(mesh, mpi_send_recv)
    use mpi
    implicit none

    type(mesh_type), intent(inout) :: mesh
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    integer :: mpi_ierr
    integer(kind=ENTIER) :: i, k, j
    integer(kind=ENTIER) :: id_elem, n_vert_tot, id_vert

    do i = 1, mpi_send_recv%n_mpi_send_neigh

      n_vert_tot = 0
      do k = 1, mpi_send_recv%mpi_send_neigh(i)%n_elems
        id_elem = mpi_send_recv%mpi_send_neigh(i)%elem_id(k)
        n_vert_tot = n_vert_tot + mesh%elem(id_elem)%n_vert
      end do

      allocate(mpi_send_recv%mpi_send_neigh(i)%sol(3, n_vert_tot))

      n_vert_tot = 1
      do k = 1, mpi_send_recv%mpi_send_neigh(i)%n_elems
        id_elem = mpi_send_recv%mpi_send_neigh(i)%elem_id(k)
        do j=1, mesh%elem(id_elem)%n_vert
          id_vert = mesh%elem(id_elem)%vert(j)
          mpi_send_recv%mpi_send_neigh(i)%sol(:, n_vert_tot) = mesh%vert(id_vert)%coord
          n_vert_tot = n_vert_tot + 1
        end do
      end do

      call mpi_isend(mpi_send_recv%mpi_send_neigh(i)%sol(1, 1), &
        3*(n_vert_tot-1), MPI_DOUBLE, &
        mpi_send_recv%mpi_send_neigh(i)%partition_id, &
        0, MPI_COMM_WORLD, mpi_send_recv%mpi_reqsend(i), mpi_ierr)
    end do

    do i = 1, mpi_send_recv%n_mpi_recv_neigh
      n_vert_tot = 0
      do k = 1, mpi_send_recv%mpi_recv_neigh(i)%n_elems
        id_elem = mpi_send_recv%mpi_recv_neigh(i)%elem_id(k)
        n_vert_tot = n_vert_tot + mesh%elem(id_elem)%n_vert
      end do
      allocate(mpi_send_recv%mpi_recv_neigh(i)%sol(3, n_vert_tot))
      call mpi_irecv(mpi_send_recv%mpi_recv_neigh(i)%sol(1, 1), &
        3*n_vert_tot, MPI_DOUBLE, &
        mpi_send_recv%mpi_recv_neigh(i)%partition_id, &
        MPI_ANY_TAG, MPI_COMM_WORLD, mpi_send_recv%mpi_reqrecv(i), mpi_ierr)
    end do

    call mpi_waitall(mpi_send_recv%n_mpi_send_neigh, &
      mpi_send_recv%mpi_reqsend, mpi_send_recv%mpi_sendstat, mpi_ierr)
    call mpi_waitall(mpi_send_recv%n_mpi_recv_neigh, &
      mpi_send_recv%mpi_reqrecv, mpi_send_recv%mpi_recvstat, mpi_ierr)

    do i = 1, mpi_send_recv%n_mpi_recv_neigh
      n_vert_tot = 1
      do k = 1, mpi_send_recv%mpi_recv_neigh(i)%n_elems
        id_elem = mpi_send_recv%mpi_recv_neigh(i)%elem_id(k)
        do j = 1, mesh%elem(id_elem)%n_vert
          id_vert = mesh%elem(id_elem)%vert(j)
          mesh%vert(id_vert)%coord = mpi_send_recv%mpi_recv_neigh(i)%sol(:, n_vert_tot)
          n_vert_tot = n_vert_tot + 1
        end do
      end do
    end do

    do i = 1, mpi_send_recv%n_mpi_send_neigh
      deallocate(mpi_send_recv%mpi_send_neigh(i)%sol)
    end do
    do i = 1, mpi_send_recv%n_mpi_recv_neigh
      deallocate(mpi_send_recv%mpi_recv_neigh(i)%sol)
    end do
  end subroutine mpi_memory_exchange_vert

  pure subroutine sol_isentropic_vortex(coord, w, t)
    use lagrange_global_data_module
    implicit none

    real(kind=DOUBLE), intent(in) :: t
    real(kind=DOUBLE), dimension(3), intent(in) :: coord
    real(kind=DOUBLE), dimension(5), intent(inout) :: w

    real(kind=DOUBLE), dimension(3) :: center_coord, vel
    real(kind=DOUBLE) :: beta, r

    beta = 5.0_DOUBLE
    vel(:) = (/0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE/)
    center_coord(:) = (/0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE/) + vel*t
    r = (coord(1) - center_coord(1))**2 + (coord(2) - center_coord(2))**2
    w(1) = 1.0_DOUBLE*(1.0_DOUBLE - ((gamma - 1)*beta**2)/(8.0_DOUBLE*gamma*pi**2)* &
      exp(1.0_DOUBLE - r))**(1.0_DOUBLE/(gamma - 1.0_DOUBLE))
    w(2) = vel(1) &
      - (coord(2) - center_coord(2))*beta/(2.0_DOUBLE*pi)*exp(0.5_DOUBLE*(1.0_DOUBLE - r))
    w(3) = vel(2) &
      + (coord(1) - center_coord(1))*beta/(2.0_DOUBLE*pi)*exp(0.5_DOUBLE*(1.0_DOUBLE - r))
    w(4) = 0.0_DOUBLE
    w(5) = w(1)**gamma
  end subroutine sol_isentropic_vortex

  subroutine compute_nodal_velocity_sidil(mesh, id_vert, sol, vp, method, b2d, b2d_h)
    use lagrange_global_data_module, only : gamma, boundary_2d
    use linear_solver_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert, method
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3), intent(inout) :: vp
    real(kind=DOUBLE), intent(in) :: b2d_h
    logical, intent(in) :: b2d

    integer(kind=ENTIER) :: j, k, le, re
    integer(kind=ENTIER) :: id_sub_elem, id_elem
    integer(kind=ENTIER) :: id_sub_face, id_face
    real(kind=DOUBLE) :: rho_p, a_p, h_p
    real(kind=DOUBLE), dimension(3) :: grad_p, Bp
    real(kind=DOUBLE), dimension(5) :: sol_w, sol_l, sol_r
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(3,3) :: mat

    vp = 0.0_DOUBLE
    rho_p = 0.0_DOUBLE
    a_p = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_w = lag_to_primit(sol(:, id_elem))
      if( sol_w(5) < 0.0_DOUBLE ) then
        print*,"Error pos", id_elem, mesh%elem(id_elem)%coord
      end if
      vp = vp &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(2:4)
      rho_p = rho_p &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(1)
      a_p = a_p &
        + mesh%sub_elem(id_sub_elem)%volume &
        * sqrt(gamma*sol_w(5)/sol_w(1))
    end do
    rho_p = rho_p / mesh%vert(id_vert)%volume
    a_p = a_p / mesh%vert(id_vert)%volume
    vp = vp / mesh%vert(id_vert)%volume

    grad_p = 0.0_DOUBLE
    mat = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      le = mesh%face(id_face)%left_neigh
      re = mesh%face(id_face)%right_neigh
      sol_l = sol(:, le)
      if( re > 0 ) then
        sol_r = sol(:, re)
      else
        sol_r = sol(:, le)
        sol_r(2:4) = sol_r(2:4) - 2.0_DOUBLE*dot_product(sol_r(2:4), &
          mesh%face(id_face)%norm)*mesh%face(id_face)%norm
      end if
      sol_w_l = lag_to_primit(sol_l)
      sol_w_r = lag_to_primit(sol_r)
      grad_p = grad_p + (sol_w_r(5) - sol_w_l(5)) &
        * mesh%sub_face(id_sub_face)%area&
        * mesh%face(id_face)%norm
    end do
    grad_p = grad_p / mesh%vert(id_vert)%volume

    if( mesh%vert(id_vert)%is_bound ) then
      Bp = boundary_normal(mesh, id_vert)
      Bp = Bp/norm2(Bp)
      grad_p = grad_p - dot_product(grad_p, Bp)*Bp
    end if

    h_p = compute_length(mesh, id_vert, method, b2d, b2d_h)
    vp = vp - 0.5_DOUBLE*h_p/(rho_p*a_p)*grad_p
  end subroutine compute_nodal_velocity_sidil

  subroutine compute_nodal_pressure_sidil(mesh, id_vert, sol, pp, method, b2d, b2d_h)
    use lagrange_global_data_module, only : gamma
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert, method
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), intent(inout) :: pp
    real(kind=DOUBLE), intent(in) :: b2d_h
    logical :: b2d

    integer(kind=ENTIER) :: j, id_elem, id_sub_elem, le, re
    integer(kind=ENTIER) :: id_face, id_sub_face
    real(kind=DOUBLE), dimension(5) :: sol_w
    real(kind=DOUBLE) :: rho_p, a_p, h_p, div_v
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, sol_l, sol_r

    pp = 0.0_DOUBLE
    rho_p = 0.0_DOUBLE
    a_p = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_w = lag_to_primit(sol(:, id_elem))
      pp = pp &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(5)
      rho_p = rho_p &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(1)
      a_p = a_p &
        + mesh%sub_elem(id_sub_elem)%volume &
        * sqrt(gamma*sol_w(5)/sol_w(1))
    end do
    pp = pp / mesh%vert(id_vert)%volume
    rho_p = rho_p / mesh%vert(id_vert)%volume
    a_p = a_p / mesh%vert(id_vert)%volume

    div_v = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      le = mesh%face(id_face)%left_neigh
      re = mesh%face(id_face)%right_neigh
      sol_l = sol(:, le)
      if( re > 0 ) then
        sol_r = sol(:, re)
      else
        sol_r = sol(:, le)
        sol_r(2:4) = sol_r(2:4) - 2.0_DOUBLE*dot_product(sol_r(2:4), &
          mesh%face(id_face)%norm)*mesh%face(id_face)%norm
      end if
      sol_w_l = lag_to_primit(sol_l)
      sol_w_r = lag_to_primit(sol_r)
      div_v = div_v &
        + dot_product(sol_w_r(2:4) - sol_w_l(2:4), &
        mesh%sub_face(id_sub_face)%area&
        * mesh%face(id_face)%norm)
    end do
    div_v = div_v / mesh%vert(id_vert)%volume
    h_p = compute_length(mesh, id_vert, method, b2d, b2d_h)
    pp = pp - 0.5_DOUBLE*h_p*rho_p*a_p*div_v
    if( pp < 0.0_DOUBLE ) then
      print*, "Neg prex in sidil nodal"
      error stop
    end if
  end subroutine compute_nodal_pressure_sidil

  function compute_length(mesh, id_vert, method, b2d, b2d_h) result(h_p)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert, method
    logical, intent(in) :: b2d
    real(kind=DOUBLE), intent(in) :: b2d_h
    real(kind=DOUBLE) :: h_p

    integer(kind=ENTIER) :: j, k
    integer(kind=ENTIER) :: id_sub_elem, id_sub_face
    integer(kind=ENTIER) :: id_elem, id_face
    real(kind=DOUBLE) :: area_sum

    area_sum = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      area_sum = area_sum + norm2(corner_normal(mesh, id_sub_elem))
    end do

    if( method == 0 ) then
      h_p = 0.0_DOUBLE*mesh%vert(id_vert)%volume/area_sum
    else if( method == 1) then
      h_p = 1.0_DOUBLE*mesh%vert(id_vert)%volume/area_sum
    else if( method == 2) then
      h_p = 2.0_DOUBLE*mesh%vert(id_vert)%volume/area_sum
    else if( method == 3) then
      h_p = 4.0_DOUBLE*mesh%vert(id_vert)%volume/area_sum
    else if( method == 4) then
      h_p = 8.0_DOUBLE*mesh%vert(id_vert)%volume/area_sum
    end if
  end function compute_length

  function boundary_normal(mesh, id_vert) result(Bp)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert

    real(kind=DOUBLE), dimension(3) :: Bp
    integer(kind=ENTIER) :: j
    integer(kind=ENTIER) :: id_sub_face, id_face

    Bp = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      if( mesh%face(id_face)%right_neigh <= 0 ) then
        Bp = Bp &
          + mesh%sub_face(id_sub_face)%area*mesh%face(id_face)%norm
      end if
    end do
  end function boundary_normal

  function corner_normal(mesh, id_sub_elem) result(norm)
    use lagrange_global_data_module, only: boundary_2d
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_sub_elem
    real(kind=DOUBLE), dimension(3) :: norm

    integer(kind=ENTIER) :: j, id_sub_face, id_face

    norm = 0.0_DOUBLE
    do j=1, mesh%sub_elem(id_sub_elem)%n_sub_faces
      id_sub_face = mesh%sub_elem(id_sub_elem)%sub_face(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      if( .not. boundary_2d .or. &
        ( boundary_2d .and. abs(mesh%face(id_face)%norm(3)) < 1e-12_DOUBLE) ) then
        if( mesh%face(id_face)%left_neigh  &
          == mesh%sub_elem(id_sub_elem)%mesh_elem ) then
          norm = norm + mesh%sub_face(id_sub_face)%area*mesh%face(id_face)%norm
        else
          norm = norm - mesh%sub_face(id_sub_face)%area*mesh%face(id_face)%norm
        end if
      end if
    end do
  end function corner_normal

  function lag_to_primit(u) result(w)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: u
    real(kind=DOUBLE), dimension(5) :: w

    w(1) = 1.0_DOUBLE/u(1)
    w(2:4) = u(2:4)
    w(5) = pressure(u)
  end function lag_to_primit

  subroutine move_mesh(mesh, vp, dt)
    implicit none

    type(mesh_type), intent(inout) :: mesh
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: vp
    real(kind=DOUBLE), intent(in) :: dt

    integer(kind=ENTIER) :: i

    do i=1, mesh%n_vert
      mesh%vert(i)%coord = mesh%vert(i)%coord + dt*vp(:, i)
    end do
  end subroutine move_mesh
end module lagrange_module
