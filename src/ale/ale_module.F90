! Self-contained ALE Euler solver: only depends on core (precision, mesh,
! mpi, linear_algebra) and ale_global_data_module. It does not reuse or
! touch ns/lagrange's shared modules, so that changes here (in particular
! the ALE-specific multi_point below) can never affect subfvns/subfvlagrange.
module ale_module
  use precision_module
  use mesh_module
  use ale_global_data_module, only: gamma, n_bc, bc_type, bc_val, bc_is_wall, boundary_2d
  implicit none

contains

  ! `g` (ratio of specific heats) is optional and defaults to the module's
  ! single-material `gamma` -- every single-material call site is
  ! unaffected. Pass it explicitly for a multi-material cell/state (the
  ! shock-bubble test's helium bubble vs. ambient air).
  pure function primit_to_conserv(w, g) result(u)
    implicit none
    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE), intent(in), optional :: g
    real(kind=DOUBLE), dimension(5) :: u
    real(kind=DOUBLE) :: gg

    gg = gamma
    if (present(g)) gg = g

    u(1) = w(1)
    u(2:4) = w(2:4)*w(1)
    u(5) = w(5)/(gg - 1.0_DOUBLE) + 0.5_DOUBLE*w(1)*(w(2)**2 + w(3)**2 + w(4)**2)
  end function primit_to_conserv

  pure function conserv_to_primit(u, g) result(w)
    implicit none
    real(kind=DOUBLE), dimension(5), intent(in) :: u
    real(kind=DOUBLE), intent(in), optional :: g
    real(kind=DOUBLE), dimension(5) :: w
    real(kind=DOUBLE) :: gg

    gg = gamma
    if (present(g)) gg = g

    w(1) = u(1)
    w(2:4) = u(2:4)/u(1)
    w(5) = (gg - 1.0_DOUBLE)*(u(5) - 0.5_DOUBLE*w(1)*(w(2)**2 + w(3)**2 + w(4)**2))
  end function conserv_to_primit

  pure function sound_speed_w(w, g) result(a)
    implicit none
    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE), intent(in), optional :: g
    real(kind=DOUBLE) :: a
    real(kind=DOUBLE) :: gg

    gg = gamma
    if (present(g)) gg = g

    a = sqrt(gg*w(5)/w(1))
  end function sound_speed_w

  pure function is_wall(re)
    implicit none
    integer(kind=ENTIER), intent(in) :: re
    logical :: is_wall

    if (re > 0) then
      is_wall = .FALSE.
    else if (re == 0) then
      is_wall = .TRUE.
    else
      is_wall = bc_is_wall(-re)
    end if
  end function is_wall

  ! Only 'wall'/'piston' (mirror the normal velocity) and 'freestream'
  ! (fixed bc_val state) boundaries are supported -- see
  ! ale_global_data_module. The ghost is the same material as the interior
  ! cell, so it's built/read back with the interior cell's own gamma `g`.
  subroutine compute_right_state(mesh, id_sub_face, re, sol_l, sol_r, g)
    implicit none
    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_sub_face, re
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_l
    real(kind=DOUBLE), dimension(5), intent(out) :: sol_r
    real(kind=DOUBLE), intent(in), optional :: g

    if (re == 0 .or. bc_is_wall(-re)) then
      sol_r = sol_l
      sol_r(2:4) = sol_r(2:4) &
        - 2.0_DOUBLE*dot_product(sol_r(2:4), mesh%sub_face(id_sub_face)%norm)*mesh%sub_face(id_sub_face)%norm
    else
      sol_r = primit_to_conserv(bc_val(:, -re), g)
    end if
  end subroutine compute_right_state

  ! First-order left/right primitive reconstruction around a vertex.
  ! gl/gr are the (possibly different, multi-material) gammas of the
  ! left/right cells; for a boundary sub-face (re<=0) gr is unused, the
  ! ghost reuses gl since it's the same material as the interior cell.
  subroutine reconstruct_lr_w(mesh, sol, id_vert, id_sub_face, le, re, gl, gr, sol_w_l, sol_w_r)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    integer(kind=ENTIER), intent(in) :: id_vert, id_sub_face, le, re
    real(kind=DOUBLE), intent(in) :: gl, gr
    real(kind=DOUBLE), dimension(5), intent(out) :: sol_w_l, sol_w_r

    real(kind=DOUBLE), dimension(5) :: sol_ghost

    sol_w_l = conserv_to_primit(sol(:, le), gl)
    if (re > 0) then
      sol_w_r = conserv_to_primit(sol(:, re), gr)
    else
      call compute_right_state(mesh, id_sub_face, re, primit_to_conserv(sol_w_l, gl), sol_ghost, gl)
      sol_w_r = conserv_to_primit(sol_ghost, gl)
    end if
  end subroutine reconstruct_lr_w

  ! ALE multi_point Riemann solver: same MPCC-type two-state flux as the
  ! fixed-mesh scheme, but the whole Riemann fan is solved relative to a
  ! face moving at normal velocity wn (=w_p.n): vn_l, vn_r and the contact
  ! velocity vn_nodal are all relative quantities (vn_nodal is expected to
  ! already be (v_p-w_p).n from the caller). Wave speeds/upwind weights
  ! (vn_l, vn_r, v_et, sl, sr) stay relative -- that's the correct
  ! criterion for whether a wave reaches the moving face. The conserved
  ! states sol_l/sol_r stay absolute (computed from the unshifted
  ! sol_w_l/sol_w_r), and the three energy-flux terms that mix an
  ! absolute state with a relative velocity ((E+p)*vn in fl/fr, and the
  ! pl*vn_l/pr_bar*v_et work terms inside the star states) get an explicit
  ! +wn*p correction so the returned flux matches F(U_abs) - wn*U_abs
  ! exactly for mass and momentum, and up to the same convention for
  ! energy (derived by requiring every occurrence of a relative velocity
  ! multiplying a one-sided pressure to use the absolute velocity instead).
  subroutine multi_point_ale(sol_w_l, sol_w_r, n, lr_flux, vn_nodal, wn, &
      lambda_l, lambda_r, sl, sr, gl, gr, dbg_id_vert, dbg_coord)
    implicit none

    real(kind=DOUBLE), intent(in) :: vn_nodal, wn
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), intent(in), optional :: gl, gr
    integer(kind=ENTIER), intent(in), optional :: dbg_id_vert
    real(kind=DOUBLE), dimension(3), intent(in), optional :: dbg_coord
    real(kind=DOUBLE), dimension(5) :: sol_l, sol_r
    real(kind=DOUBLE), intent(inout) :: lambda_l, lambda_r
    real(kind=DOUBLE), dimension(5, 2), intent(inout) :: lr_flux
    real(kind=DOUBLE), intent(inout) :: sl, sr

    real(kind=DOUBLE) :: rhol, rhor, vn_l, vn_r
    real(kind=DOUBLE) :: pl, pr, rhol_et, rhor_et
    real(kind=DOUBLE) :: v_et, el, er
    real(kind=DOUBLE) :: pl_bar, pr_bar
    real(kind=DOUBLE), dimension(5) :: fl, fr
    real(kind=DOUBLE), dimension(5) :: sol_l_et, sol_r_et

    rhol = sol_w_l(1)
    vn_l = dot_product(sol_w_l(2:4), n) - wn
    pl = sol_w_l(5)
    sol_l = primit_to_conserv(sol_w_l, gl)
    el = sol_l(5)/rhol

    rhor = sol_w_r(1)
    vn_r = dot_product(sol_w_r(2:4), n) - wn
    pr = sol_w_r(5)
    sol_r = primit_to_conserv(sol_w_r, gr)
    er = sol_r(5)/rhor

    fl(1)   = vn_l*sol_l(1)
    fl(2:4) = vn_l*sol_l(2:4) + pl*n
    fl(5)   = (sol_l(5) + pl)*vn_l + wn*pl

    fr(1)   = vn_r*sol_r(1)
    fr(2:4) = vn_r*sol_r(2:4) + pr*n
    fr(5)   = (sol_r(5) + pr)*vn_r + wn*pr

    v_et = vn_nodal

    rhol_et = 1.0_DOUBLE/(1.0_DOUBLE/rhol + (v_et - vn_l)/lambda_l)
    pl_bar = pl - lambda_l*(v_et - vn_l)

    sol_l_et(1)   = rhol_et
    sol_l_et(2:4) = rhol_et*(sol_w_l(2:4) + (v_et - vn_l)*n)
    sol_l_et(5)   = rhol_et*(el + (pl*vn_l - pl_bar*v_et + wn*(pl - pl_bar))/lambda_l)

    rhor_et = 1.0_DOUBLE/(1.0_DOUBLE/rhor + (vn_r - v_et)/lambda_r)
    pr_bar = pr + lambda_r*(v_et - vn_r)

    sol_r_et(1)   = rhor_et
    sol_r_et(2:4) = rhor_et*(sol_w_r(2:4) + (v_et - vn_r)*n)
    sol_r_et(5)   = rhor_et*(er + (pr_bar*v_et - pr*vn_r + wn*(pr_bar - pr))/lambda_r)

    if (rhol_et < 0.0_DOUBLE .or. rhor_et < 0.0_DOUBLE) then
      print *, "Negative specific volume MPCC (ale) !", rhol_et, rhor_et
      if (present(dbg_id_vert)) print *, "  id_vert =", dbg_id_vert
      if (present(dbg_coord)) print *, "  coord =", dbg_coord
      print *, "  rhol,pl,vn_l,lambda_l =", rhol, pl, vn_l, lambda_l
      print *, "  rhor,pr,vn_r,lambda_r =", rhor, pr, vn_r, lambda_r
      print *, "  v_et,wn,n =", v_et, wn, n
      error stop
    end if

    sl = vn_l - lambda_l/rhol
    sr = vn_r + lambda_r/rhor

    lr_flux(:, 1) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(v_et)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et)) &
      - 0.5_DOUBLE*(pr_bar - pl_bar)* &
      (/0.0_DOUBLE, n, v_et + wn/)

    lr_flux(:, 2) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(v_et)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et)) &
      + 0.5_DOUBLE*(pr_bar - pl_bar)* &
      (/0.0_DOUBLE, n, v_et + wn/)

    lr_flux(:, 2) = -lr_flux(:, 2)
  end subroutine multi_point_ale

  ! Nodal Lagrangian velocity solve (bc_style=1 convention only: wall
  ! sub-faces are simply excluded from the local pseudo-inverse system,
  ! the node's velocity is set by its non-wall neighbors alone).
  subroutine compute_lambdas_and_solve_nodal_velocity(mesh, id_vert, sol_w_lr, gamma_lr, lambda, v_bars, v_node)
    use linear_solver_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, 2, mesh%vert(id_vert)%n_sub_faces_neigh), intent(in) :: sol_w_lr
    real(kind=DOUBLE), dimension(2, mesh%vert(id_vert)%n_sub_faces_neigh), intent(in) :: gamma_lr
    real(kind=DOUBLE), dimension(2, mesh%vert(id_vert)%n_sub_faces_neigh), intent(inout) :: lambda
    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_faces_neigh), intent(inout) :: v_bars
    real(kind=DOUBLE), dimension(3), intent(inout) :: v_node

    integer(kind=ENTIER) :: iter
    integer(kind=ENTIER) :: j, id_sub_face, id_face, re
    real(kind=DOUBLE) :: rhol, vn_l, pl, al, lambda_l
    real(kind=DOUBLE) :: rhor, vn_r, pr, ar, lambda_r
    real(kind=DOUBLE) :: v_et
    real(kind=DOUBLE), dimension(3) :: Rp
    real(kind=DOUBLE), dimension(3, 3) :: mat
    ! ns's own copy of this iterative lambda ramp-up (ns_euler_rs_module's
    ! compute_lambdas_and_solve_nodal_velocity) uses 4 iterations; strong,
    ! still-forming shocks (piston-driven, before the flow has smoothed
    ! out) can need a couple more passes for lambda_l/lambda_r to grow
    ! enough to keep the multi_point star-state densities/pressures
    ! positive -- multi_point_ale reuses this same lambda unmodified.
    integer(kind=ENTIER), parameter :: n_lambda_iter = 6

    v_node = 0.0_DOUBLE
    iter = 0
    do while (iter < n_lambda_iter)
      iter = iter + 1

      mat(:, :) = 0.0_DOUBLE
      Rp(:) = 0.0_DOUBLE

      do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
        id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        re = mesh%face(id_face)%right_neigh

        rhol = sol_w_lr(1, 1, j)
        vn_l = dot_product(sol_w_lr(2:4, 1, j), mesh%sub_face(id_sub_face)%norm)
        pl = sol_w_lr(5, 1, j)
        al = sound_speed_w(sol_w_lr(:, 1, j), gamma_lr(1, j))

        rhor = sol_w_lr(1, 2, j)
        vn_r = dot_product(sol_w_lr(2:4, 2, j), mesh%sub_face(id_sub_face)%norm)
        pr = sol_w_lr(5, 2, j)
        ar = sound_speed_w(sol_w_lr(:, 2, j), gamma_lr(2, j))

        lambda_l = lambda(1, j)
        lambda_r = lambda(2, j)

        if (iter == 1) then
          lambda_l = max(lambda_l, al*rhol, &
            sqrt(rhol*max(0.0_DOUBLE, pr - pl)), &
            -rhol*(vn_r - vn_l))
          lambda_r = max(lambda_r, ar*rhor, &
            sqrt(rhor*max(0.0_DOUBLE, pl - pr)), &
            -rhor*(vn_r - vn_l))

          v_bars(j) = (lambda_l*vn_l + lambda_r*vn_r - (pr - pl))/(lambda_r + lambda_l)

          lambda_l = max(lambda_l, al*rhol*(1.0_DOUBLE &
            + 1.5_DOUBLE*max(0.0_DOUBLE, -(v_bars(j) - vn_l)/al)))
          lambda_r = max(lambda_r, ar*rhor*(1.0_DOUBLE &
            + 1.5_DOUBLE*max(0.0_DOUBLE, (v_bars(j) - vn_r)/ar)))
        else
          if (is_wall(re)) then
            lambda_l = max(lambda_l, al*rhol*(1.0_DOUBLE &
              + 1.5_DOUBLE*max(0.0_DOUBLE, -(v_bars(j) - vn_l)/al)))
            lambda_r = max(lambda_r, ar*rhor*(1.0_DOUBLE &
              + 1.5_DOUBLE*max(0.0_DOUBLE, (v_bars(j) - vn_r)/ar)))
          else
            v_et = dot_product(v_node, mesh%sub_face(id_sub_face)%norm)
            lambda_l = max(lambda_l, al*rhol*(1.0_DOUBLE &
              + 1.5_DOUBLE*max(0.0_DOUBLE, -(v_et - vn_l)/al)))
            lambda_r = max(lambda_r, ar*rhor*(1.0_DOUBLE &
              + 1.5_DOUBLE*max(0.0_DOUBLE, (v_et - vn_r)/ar)))
          end if
        end if

        lambda(1, j) = lambda_l
        lambda(2, j) = lambda_r

        v_bars(j) = (lambda_l*vn_l + lambda_r*vn_r - (pr - pl))/(lambda_r + lambda_l)

        ! bc_style=0 convention (ns_euler_rs_module): every sub-face, wall
        ! or not, contributes to mat *and* Rp with the same acoustic
        ! v_bars(j) term -- a wall's ghost state is already the mirrored
        ! state built by reconstruct_lr_w/compute_right_state, so v_bars(j)
        ! there naturally works out to the correct no-penetration value
        ! (e.g. exactly 0 for a still wall with matching pl=pr) without
        ! needing separate pressure-force bookkeeping. Excluding walls
        ! entirely (bc_style=1) instead leaves the 3x3 system under-
        ! determined at a vertex touching two walls (e.g. the piston/floor
        ! corner), which showed up as an erratic v_node there and,
        ! downstream, an insufficient lambda ramp-up and a negative star
        ! density in multi_point_ale.
        mat = mat + (lambda_r + lambda_l)*mesh%sub_face(id_sub_face)%area &
          *tensor_product_3(mesh%sub_face(id_sub_face)%norm, mesh%sub_face(id_sub_face)%norm)
        Rp = Rp + mesh%sub_face(id_sub_face)%area* &
          (lambda_r + lambda_l)*v_bars(j)*mesh%sub_face(id_sub_face)%norm
      end do

      call pseudo_inverse_inplace_lapack(3, mat)
      v_node = matmul(mat, Rp)
    end do
  end subroutine compute_lambdas_and_solve_nodal_velocity

  ! Prescribed arbitrary grid velocity w(x,t) at the mesh vertices. Zero
  ! for now: reduces the ALE scheme to the fixed-mesh multi_point scheme
  ! exactly (v_p-w_p=v_p), the reference case to validate against.
  subroutine compute_grid_velocity(mesh, t, wp)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), intent(in) :: t
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(out) :: wp

    wp = 0.0_DOUBLE
  end subroutine compute_grid_velocity

  ! Nodal Lagrangian fluid velocity only (reconstruction + the same nodal
  ! solve used inside compute_rhs_ale, but no flux/rhs). v_p does not
  ! depend on w_p, so this can be called on its own to prescribe w_p=v_p
  ! (the Lagrangian limit of the ALE scheme) before the real flux pass.
  subroutine compute_nodal_velocity_field(mesh, sol, gamma_arr, vp)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(out) :: vp

    integer(kind=ENTIER) :: id_vert, j, nsfn, le, re, id_sub_face, id_face
    real(kind=DOUBLE), dimension(3) :: v_vert
    real(kind=DOUBLE), allocatable :: sol_w_lr(:, :, :), gamma_lr(:, :), lambda(:, :), v_bars(:)

    do id_vert = 1, mesh%n_vert
      nsfn = mesh%vert(id_vert)%n_sub_faces_neigh
      allocate(sol_w_lr(5, 2, nsfn), gamma_lr(2, nsfn), lambda(2, nsfn), v_bars(nsfn))
      sol_w_lr = 0.0_DOUBLE
      lambda = 0.0_DOUBLE

      do j = 1, nsfn
        id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        le = mesh%face(id_face)%left_neigh
        re = mesh%face(id_face)%right_neigh
        gamma_lr(1, j) = gamma_arr(le)
        gamma_lr(2, j) = merge(gamma_arr(re), gamma_arr(le), re > 0)
        call reconstruct_lr_w(mesh, sol, id_vert, id_sub_face, le, re, &
          gamma_lr(1, j), gamma_lr(2, j), sol_w_lr(:, 1, j), sol_w_lr(:, 2, j))
      end do

      call compute_lambdas_and_solve_nodal_velocity(mesh, id_vert, sol_w_lr, gamma_lr, lambda, v_bars, v_vert)
      vp(:, id_vert) = v_vert

      deallocate(sol_w_lr, gamma_lr, lambda, v_bars)
    end do
  end subroutine compute_nodal_velocity_field

  ! Explicit global CFL time step, from the per-cell acoustic-impedance sum.
  subroutine compute_dt_ale(mesh, sum_lambda, cfl, dt)
    use mpi
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: sum_lambda
    real(kind=DOUBLE), intent(in) :: cfl
    real(kind=DOUBLE), intent(out) :: dt

    integer(kind=ENTIER) :: i, mpi_ierr

    dt = 1e10_DOUBLE
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) dt = min(dt, mesh%elem(i)%volume/sum_lambda(i))
    end do
    dt = cfl*dt

    call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, mpi_ierr)
  end subroutine compute_dt_ale

  subroutine compute_rhs_ale(mesh, sol, gamma_arr, wp, rhs, sum_lambda, vp)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: wp
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(out) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(out) :: sum_lambda
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(out) :: vp

    integer(kind=ENTIER) :: id_vert, j, idse, ide, nsen, max_nsen
    real(kind=DOUBLE), dimension(3) :: v_vert
    real(kind=DOUBLE), allocatable :: sum_lambda_vert(:), flux_sum_vert(:, :)

    rhs = 0.0_DOUBLE
    sum_lambda = 1e-12_DOUBLE
    vp = 0.0_DOUBLE

    max_nsen = 0
    do id_vert = 1, mesh%n_vert
      max_nsen = max(max_nsen, mesh%vert(id_vert)%n_sub_elems_neigh)
    end do
    allocate(sum_lambda_vert(max_nsen), flux_sum_vert(5, max_nsen))

    do id_vert = 1, mesh%n_vert
      nsen = mesh%vert(id_vert)%n_sub_elems_neigh
      sum_lambda_vert(1:nsen) = 0.0_DOUBLE
      flux_sum_vert(:, 1:nsen) = 0.0_DOUBLE

      call compute_rhs_around_vert_ale(mesh, sol, gamma_arr, wp, &
        nsen, sum_lambda_vert, flux_sum_vert, id_vert, v_vert)
      vp(:, id_vert) = v_vert

      do j = 1, mesh%vert(id_vert)%n_sub_elems_neigh
        idse = mesh%vert(id_vert)%sub_elem_neigh(j)
        ide = mesh%sub_elem(idse)%mesh_elem
        rhs(:, ide) = rhs(:, ide) - flux_sum_vert(:, j)
        sum_lambda(ide) = sum_lambda(ide) + sum_lambda_vert(j)
      end do
    end do

    deallocate(sum_lambda_vert, flux_sum_vert)
  end subroutine compute_rhs_ale

  subroutine compute_rhs_around_vert_ale(mesh, sol, gamma_arr, wp, &
      nsen, sum_lambda_vert, flux_sum_vert, id_vert, v_vert)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: wp
    integer(kind=ENTIER), intent(in) :: nsen, id_vert
    real(kind=DOUBLE), dimension(nsen), intent(inout) :: sum_lambda_vert
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: flux_sum_vert
    real(kind=DOUBLE), dimension(3), intent(out) :: v_vert

    integer(kind=ENTIER) :: nsfn, j, id_sub_face, id_face
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE) :: sl, sr, vn_nodal, wn
    real(kind=DOUBLE), allocatable :: sol_w_lr(:, :, :), gamma_lr(:, :), lr_flux(:, :, :)
    real(kind=DOUBLE), allocatable :: lambda(:, :), v_bars(:)

    nsfn = mesh%vert(id_vert)%n_sub_faces_neigh
    allocate(sol_w_lr(5, 2, nsfn), gamma_lr(2, nsfn), lr_flux(5, 2, nsfn), lambda(2, nsfn), v_bars(nsfn))
    rse_loc = 0
    sol_w_lr = 0.0_DOUBLE
    lambda = 0.0_DOUBLE

    do j = 1, nsfn
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      le = mesh%face(id_face)%left_neigh
      re = mesh%face(id_face)%right_neigh
      gamma_lr(1, j) = gamma_arr(le)
      gamma_lr(2, j) = merge(gamma_arr(re), gamma_arr(le), re > 0)
      call reconstruct_lr_w(mesh, sol, id_vert, id_sub_face, le, re, &
        gamma_lr(1, j), gamma_lr(2, j), sol_w_lr(:, 1, j), sol_w_lr(:, 2, j))
    end do

    call compute_lambdas_and_solve_nodal_velocity(mesh, id_vert, sol_w_lr, gamma_lr, lambda, v_bars, v_vert)

    do j = 1, nsfn
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      re = mesh%face(id_face)%right_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if (rse > 0) rse_loc = mesh%sub_elem(rse)%id_loc_around_node

      wn = dot_product(wp(:, id_vert), mesh%sub_face(id_sub_face)%norm)

      if (is_wall(re)) then
        ! No-penetration at the (possibly moving) wall: zero relative
        ! normal velocity between fluid and grid -- same as the fixed-mesh
        ! scheme's vn_nodal=0, which is already the wp=0 special case of it.
        vn_nodal = 0.0_DOUBLE
      else
        vn_nodal = dot_product(v_vert - wp(:, id_vert), mesh%sub_face(id_sub_face)%norm)
      end if

      call multi_point_ale(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
        mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), vn_nodal, wn, &
        lambda(1, j), lambda(2, j), sl, sr, gamma_lr(1, j), gamma_lr(2, j), &
        id_vert, mesh%vert(id_vert)%coord)

      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        sum_lambda_vert(lse_loc) = sum_lambda_vert(lse_loc) &
          + mesh%sub_face(id_sub_face)%area*max(0.0_DOUBLE, -sl)
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) &
          + mesh%sub_face(id_sub_face)%area*lr_flux(:, 1, j)

        if (rse > 0) then
          if (mesh%sub_elem(rse)%mesh_vert == id_vert) then
            sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
              + mesh%sub_face(id_sub_face)%area*max(0.0_DOUBLE, sr)
            flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) &
              + mesh%sub_face(id_sub_face)%area*lr_flux(:, 2, j)
          end if
        end if
      end if
    end do

    deallocate(sol_w_lr, gamma_lr, lr_flux, lambda, v_bars)
  end subroutine compute_rhs_around_vert_ale

  ! gamma_arr is set here too (uniform `gamma` unless init_shock_bubble),
  ! since it's the initial condition that decides material placement.
  subroutine init_sol(mesh, sol, gamma_arr)
    use ale_global_data_module, only: init_uniform, sol_uniform, &
      init_1drp, x1drp, sol_w_1drp_l, sol_w_1drp_r, &
      init_shock_bubble, xc_bub, yc_bub, r_bub, gamma_air, gamma_bub, rho_bub
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: sol
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(inout) :: gamma_arr

    integer(kind=ENTIER) :: i
    real(kind=DOUBLE) :: rb
    real(kind=DOUBLE), dimension(5) :: w

    gamma_arr = gamma

    if (init_uniform) then
      do i = 1, mesh%n_elems
        sol(:, i) = primit_to_conserv(sol_uniform)
      end do
    else if (init_1drp) then
      do i = 1, mesh%n_elems
        if (mesh%elem(i)%coord(1) < x1drp) then
          sol(:, i) = primit_to_conserv(sol_w_1drp_l)
        else
          sol(:, i) = primit_to_conserv(sol_w_1drp_r)
        end if
      end do
    else if (init_shock_bubble) then
      ! Same physical setup as test/lagrange_shock_bubble (init=9): quiescent
      ! air (rho=p=1) everywhere except the bubble, in mechanical/thermal
      ! equilibrium with the ambient (same p), He+28%air by mass so
      ! gamma=1.648 and rho=rho_bub. The Ms=1.22 shock is driven by the
      ! 'piston' BC at the right wall, not by an initial discontinuity here.
      w(2:4) = 0.0_DOUBLE
      do i = 1, mesh%n_elems
        rb = sqrt((mesh%elem(i)%coord(1) - xc_bub)**2 + (mesh%elem(i)%coord(2) - yc_bub)**2)
        if (rb <= r_bub) then
          gamma_arr(i) = gamma_bub
          w(1) = rho_bub
        else
          gamma_arr(i) = gamma_air
          w(1) = 1.0_DOUBLE
        end if
        w(5) = 1.0_DOUBLE
        sol(:, i) = primit_to_conserv(w, gamma_arr(i))
      end do
    else
      print*, "[-] No init chosen!"
      error stop
    end if
  end subroutine init_sol

  subroutine move_mesh(mesh, wp, dt)
    implicit none
    type(mesh_type), intent(inout) :: mesh
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: wp
    real(kind=DOUBLE), intent(in) :: dt

    integer(kind=ENTIER) :: i

    do i = 1, mesh%n_vert
      mesh%vert(i)%coord = mesh%vert(i)%coord + dt*wp(:, i)
    end do
  end subroutine move_mesh

  subroutine mpi_memory_exchange_vert(mesh, mpi_send_recv)
    use mpi
    use mpi_module
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
        do j = 1, mesh%elem(id_elem)%n_vert
          id_vert = mesh%elem(id_elem)%vert(j)
          mpi_send_recv%mpi_send_neigh(i)%sol(:, n_vert_tot) = mesh%vert(id_vert)%coord
          n_vert_tot = n_vert_tot + 1
        end do
      end do

      call mpi_isend(mpi_send_recv%mpi_send_neigh(i)%sol(1, 1), &
        3*(n_vert_tot - 1), MPI_DOUBLE, &
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

end module ale_module
