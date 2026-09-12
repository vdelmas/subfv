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
  subroutine compute_rhs_lagrange(mesh, sol, vp, dt, rhs, n_bc, bc_type, bc_val, b2d, mass, &
      gamma_arr, vp_is_imposed, second_order, grad_v, grad_p, div_v, alpha_p_arr)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), intent(in) :: dt
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mass
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    logical, dimension(mesh%n_vert), intent(in) :: vp_is_imposed
    integer(kind=ENTIER), intent(in) :: n_bc
    character(len=255), dimension(n_bc) :: bc_type
    real(kind=DOUBLE), dimension(5, n_bc) :: bc_val
    logical, intent(in) :: b2d
    logical, intent(in) :: second_order
    real(kind=DOUBLE), dimension(3, 3, mesh%n_elems), intent(in) :: grad_v
    real(kind=DOUBLE), dimension(3, mesh%n_elems), intent(in) :: grad_p
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: div_v
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(out), optional :: alpha_p_arr

    integer(kind=ENTIER) :: i, j, k, nsfn, nsen, idse, ide
    real(kind=DOUBLE), dimension(:,:), allocatable :: rhs_o1, rhs_o2
    real(kind=DOUBLE), dimension(3) :: vp_o1_i, vp_o2_i
    real(kind=DOUBLE) :: alpha_p

    rhs = 0.0_DOUBLE
    if (present(alpha_p_arr)) alpha_p_arr = 1.0_DOUBLE

    do i = 1, mesh%n_vert
      nsfn = mesh%vert(i)%n_sub_faces_neigh
      nsen = mesh%vert(i)%n_sub_elems_neigh
      allocate(rhs_o1(5, nsen))

      vp_o1_i = vp(:, i)
      call compute_rhs_around_node(mesh, sol, i, nsfn, nsen, vp_o1_i, vp_is_imposed(i), &
        dt, rhs_o1, n_bc, bc_type, bc_val, b2d, mass, gamma_arr, &
        .false., grad_v, grad_p, div_v)

      if (second_order) then
        allocate(rhs_o2(5, nsen))
        vp_o2_i = vp(:, i)
        call compute_rhs_around_node(mesh, sol, i, nsfn, nsen, vp_o2_i, vp_is_imposed(i), &
          dt, rhs_o2, n_bc, bc_type, bc_val, b2d, mass, gamma_arr, &
          .true., grad_v, grad_p, div_v)

        ! If rhs_o2 or vp_o2_i has any NaN (negative reconstructed pressure → sqrt(neg) → NaN
        ! propagating through the nodal solver), skip blending entirely for this node.
        if (any(rhs_o2 /= rhs_o2) .or. any(vp_o2_i /= vp_o2_i)) then
          alpha_p = -0.01_DOUBLE
        else
          alpha_p = 1.0_DOUBLE
          do k = 1, nsen
            idse = mesh%vert(i)%sub_elem_neigh(k)
            ide = mesh%sub_elem(idse)%mesh_elem
            if (.not. mesh%elem(ide)%is_ghost) then
              alpha_p = min(alpha_p, &
                omega_pos(sol(:, ide) + dt*rhs_o1(:, k), &
                          dt*(rhs_o2(:, k) - rhs_o1(:, k))))
            end if
          end do
        end if

        if (present(alpha_p_arr)) alpha_p_arr(i) = alpha_p

        if (alpha_p > 0.0_DOUBLE) then
          vp(:, i) = vp_o1_i + alpha_p*(vp_o2_i - vp_o1_i)
          do j = 1, nsen
            idse = mesh%vert(i)%sub_elem_neigh(j)
            ide = mesh%sub_elem(idse)%mesh_elem
            rhs(:, ide) = rhs(:, ide) + rhs_o1(:,j) + alpha_p*(rhs_o2(:,j) - rhs_o1(:,j))
          end do
        else
          vp(:, i) = vp_o1_i
          do j = 1, nsen
            idse = mesh%vert(i)%sub_elem_neigh(j)
            ide = mesh%sub_elem(idse)%mesh_elem
            rhs(:, ide) = rhs(:, ide) + rhs_o1(:,j)
          end do
        end if
        deallocate(rhs_o2)
      else
        vp(:, i) = vp_o1_i
        do j = 1, nsen
          idse = mesh%vert(i)%sub_elem_neigh(j)
          ide = mesh%sub_elem(idse)%mesh_elem
          rhs(:, ide) = rhs(:, ide) + rhs_o1(:,j)
        end do
      end if

      deallocate(rhs_o1)
    end do
  end subroutine compute_rhs_lagrange

  subroutine compute_rhs_around_node(mesh, sol, i_vert, nsfn, nsen, vp_i, vp_is_imposed_i, &
      dt, rhs, n_bc, bc_type, bc_val, b2d, mass, gamma_arr, &
      second_order, grad_v, grad_p, div_v)
    use linear_solver_module, only: lu_solve, print_mat, inverse_3_by_3
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    integer(kind=ENTIER), intent(in) :: i_vert, nsfn, nsen
    real(kind=DOUBLE), dimension(3), intent(inout) :: vp_i
    logical, intent(in) :: vp_is_imposed_i
    real(kind=DOUBLE), intent(in) :: dt
    real(kind=DOUBLE), dimension(5, nsen), intent(out) :: rhs
    integer(kind=ENTIER), intent(in) :: n_bc
    character(len=255), dimension(n_bc) :: bc_type
    real(kind=DOUBLE), dimension(5, n_bc) :: bc_val
    logical, intent(in) :: b2d
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mass
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    logical, intent(in) :: second_order
    real(kind=DOUBLE), dimension(3, 3, mesh%n_elems), intent(in) :: grad_v
    real(kind=DOUBLE), dimension(3, mesh%n_elems), intent(in) :: grad_p
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: div_v

    integer(kind=ENTIER) :: j, id_face, id_sub_face
    integer(kind=ENTIER) :: idl, idr, lse, lse_loc, rse, rse_loc
    real(kind=DOUBLE), dimension(2, nsfn) :: lambda
    real(kind=DOUBLE), dimension(5, 2, nsfn) :: sol_lr
    real(kind=DOUBLE), dimension(5) :: flux
    real(kind=DOUBLE), dimension(3, 3) :: Mp, Mp_inv
    real(kind=DOUBLE), dimension(3) :: Rp, Bp
    real(kind=DOUBLE) :: vl, vr, pl, pr, v_bar, pl_et, pr_et, v_et, Pp
    real(kind=DOUBLE) :: gl, gr, al, ar, taul, taur
    logical :: is_corner, have_first_wall_norm
    real(kind=DOUBLE), dimension(3) :: first_wall_norm, dx
    real(kind=DOUBLE) :: cross_z

    lambda   = 0.0_DOUBLE
    sol_lr   = 0.0_DOUBLE
    rhs      = 0.0_DOUBLE
    Bp = 0.0_DOUBLE
    Mp = 0.0_DOUBLE
    Rp = 0.0_DOUBLE
    is_corner = .false.
    have_first_wall_norm = .false.
    first_wall_norm = 0.0_DOUBLE

    do j = 1, nsfn
      id_sub_face = mesh%vert(i_vert)%sub_face_neigh(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      idl = mesh%face(id_face)%left_neigh
      idr = mesh%face(id_face)%right_neigh
      gl = gamma_arr(idl)
      sol_lr(:, 1, j) = sol(:, idl)
      pl = pressure(sol_lr(:, 1, j), gl)
      vl = dot_product(sol_lr(2:4, 1, j), mesh%face(id_face)%norm)

      if (second_order) then
        taul = sol_lr(1, 1, j)
        al   = sqrt(gl * pl * taul)
        dx = mesh%vert(i_vert)%coord - mesh%elem(idl)%coord
        pl = pl + dot_product(grad_p(:, idl), dx) - 0.5_DOUBLE*dt*al**2/taul*div_v(idl)
        vl = vl + dot_product(matmul(grad_v(:,:,idl), dx), mesh%face(id_face)%norm) &
                - 0.5_DOUBLE*dt*taul*dot_product(grad_p(:, idl), mesh%face(id_face)%norm)
      end if

      if (idr > 0) then
        gr = gamma_arr(idr)
        sol_lr(:, 2, j) = sol(:, idr)
        pr = pressure(sol_lr(:, 2, j), gr)
        vr = dot_product(sol_lr(2:4, 2, j), mesh%face(id_face)%norm)

        if (second_order) then
          taur = sol_lr(1, 2, j)
          ar   = sqrt(gr * pr * taur)
          dx = mesh%vert(i_vert)%coord - mesh%elem(idr)%coord
          pr = pr + dot_product(grad_p(:, idr), dx) - 0.5_DOUBLE*dt*ar**2/taur*div_v(idr)
          vr = vr + dot_product(matmul(grad_v(:,:,idr), dx), mesh%face(id_face)%norm) &
                  - 0.5_DOUBLE*dt*taur*dot_product(grad_p(:, idr), mesh%face(id_face)%norm)
        end if

        lambda(1, j) = max(sqrt(gl*pl*sol_lr(1,1,j))/sol_lr(1,1,j), &
          sqrt(max(0.0_DOUBLE, pr-pl)/sol_lr(1,1,j)), -(vr-vl)/sol_lr(1,1,j))
        lambda(2, j) = max(sqrt(gr*pr*sol_lr(1,2,j))/sol_lr(1,2,j), &
          sqrt(max(0.0_DOUBLE, pl-pr)/sol_lr(1,2,j)), -(vr-vl)/sol_lr(1,2,j))

        if (.not. vp_is_imposed_i) then
          v_bar = (lambda(1,j)*vl + lambda(2,j)*vr - (pr - pl)) / (lambda(2,j) + lambda(1,j))
          Mp = Mp + mesh%sub_face(id_sub_face)%area*(lambda(1,j)+lambda(2,j)) &
            * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
          Rp = Rp + mesh%sub_face(id_sub_face)%area*(lambda(1,j)+lambda(2,j)) &
            * mesh%face(id_face)%norm * v_bar
        end if
      else
        if (b2d .and. abs(mesh%face(id_face)%norm(3)) > 1e-8_DOUBLE) then
          lambda(1, j) = sqrt(gl*pl*sol_lr(1,1,j))/sol_lr(1,1,j)
          if (.not. vp_is_imposed_i) then
            Mp = Mp + mesh%sub_face(id_sub_face)%area*2*lambda(1,j) &
              * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
          end if
        else if (idr < 0 .and. trim(bc_type(-idr)) == 'pressure') then
          lambda(1, j) = sqrt(gl*pl*sol_lr(1,1,j))/sol_lr(1,1,j)
          if (.not. vp_is_imposed_i) then
            Mp = Mp + mesh%sub_face(id_sub_face)%area*lambda(1,j) &
              * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
            Rp = Rp + mesh%sub_face(id_sub_face)%area &
              * (bc_val(1,-idr) + lambda(1,j)*vl) * mesh%face(id_face)%norm
          end if
        else
          lambda(1, j) = sqrt(gl*pl*sol_lr(1,1,j))/sol_lr(1,1,j)
          if (.not. vp_is_imposed_i) then
            if (.not. have_first_wall_norm) then
              first_wall_norm = mesh%face(id_face)%norm
              have_first_wall_norm = .true.
            else
              cross_z = first_wall_norm(1)*mesh%face(id_face)%norm(2) &
                      - first_wall_norm(2)*mesh%face(id_face)%norm(1)
              if (abs(cross_z) > 1.0e-8_DOUBLE) is_corner = .true.
            end if
            Bp = Bp + mesh%sub_face(id_sub_face)%area * mesh%face(id_face)%norm
            Mp = Mp + mesh%sub_face(id_sub_face)%area*lambda(1,j) &
              * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
            Rp = Rp + mesh%sub_face(id_sub_face)%area*(pl + lambda(1,j)*vl) &
              * mesh%face(id_face)%norm
          end if
        end if
      end if
    end do

    if (.not. vp_is_imposed_i) then
      if (is_corner) then
        vp_i = 0.0_DOUBLE
      else if (maxval(abs(Bp)) > 1e-8_DOUBLE) then
        call inverse_3_by_3(Mp, Mp_inv)
        Pp = dot_product(Rp, matmul(Mp_inv, Bp)) / dot_product(Bp, matmul(Mp_inv, Bp))
        vp_i = matmul(Mp_inv, Rp - Pp*Bp)
      else
        call inverse_3_by_3(Mp, Mp_inv)
        vp_i = matmul(Mp_inv, Rp)
      end if
    end if

    do j = 1, nsfn
      id_sub_face = mesh%vert(i_vert)%sub_face_neigh(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      idl = mesh%face(id_face)%left_neigh
      idr = mesh%face(id_face)%right_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      gl = gamma_arr(idl)

      v_et = dot_product(vp_i, mesh%face(id_face)%norm)
      vl = dot_product(sol_lr(2:4, 1, j), mesh%face(id_face)%norm)
      pl = pressure(sol_lr(:, 1, j), gl)
      if (second_order) then
        taul = sol_lr(1, 1, j)
        al   = sqrt(gl * pl * taul)
        dx = mesh%vert(i_vert)%coord - mesh%elem(idl)%coord
        pl = pl + dot_product(grad_p(:, idl), dx) - 0.5_DOUBLE*dt*al**2/taul*div_v(idl)
        vl = vl + dot_product(matmul(grad_v(:,:,idl), dx), mesh%face(id_face)%norm) &
                - 0.5_DOUBLE*dt*taul*dot_product(grad_p(:, idl), mesh%face(id_face)%norm)
      end if
      pl_et = pl - lambda(1, j)*(v_et - vl)

      flux(1)   = -v_et
      flux(2:4) = pl_et * mesh%face(id_face)%norm
      flux(5)   = pl_et * v_et
      if (mesh%sub_elem(lse)%mesh_vert == i_vert) then
        rhs(:, lse_loc) = rhs(:, lse_loc) - mesh%sub_face(id_sub_face)%area/mass(idl) * flux
      end if

      if (idr > 0) then
        gr = gamma_arr(idr)
        vr = dot_product(sol_lr(2:4, 2, j), mesh%face(id_face)%norm)
        pr = pressure(sol_lr(:, 2, j), gr)
        if (second_order) then
          taur = sol_lr(1, 2, j)
          ar   = sqrt(gr * pr * taur)
          dx = mesh%vert(i_vert)%coord - mesh%elem(idr)%coord
          pr = pr + dot_product(grad_p(:, idr), dx) - 0.5_DOUBLE*dt*ar**2/taur*div_v(idr)
          vr = vr + dot_product(matmul(grad_v(:,:,idr), dx), mesh%face(id_face)%norm) &
                  - 0.5_DOUBLE*dt*taur*dot_product(grad_p(:, idr), mesh%face(id_face)%norm)
        end if
        pr_et = pr + lambda(2, j)*(v_et - vr)
        flux(1)   = -v_et
        flux(2:4) = pr_et * mesh%face(id_face)%norm
        flux(5)   = pr_et * v_et
        rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == i_vert) then
          rse_loc = mesh%sub_elem(rse)%id_loc_around_node
          rhs(:, rse_loc) = rhs(:, rse_loc) + mesh%sub_face(id_sub_face)%area/mass(idr) * flux
        end if
      end if
    end do
  end subroutine compute_rhs_around_node

  subroutine compute_gradients(mesh, sol, gamma_arr, dt, grad_v, grad_p, div_v, p_limited)
    use linear_solver_module, only: tensor_product_3, inverse_3_by_3
    ! Nodal WENO Green-Gauss gradients weighted by ||grad_p||.
    ! div_v = trace(grad_v).
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    real(kind=DOUBLE), intent(in) :: dt
    real(kind=DOUBLE), dimension(3, 3, mesh%n_elems), intent(out) :: grad_v
    real(kind=DOUBLE), dimension(3, mesh%n_elems), intent(out) :: grad_p
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(out) :: div_v
    integer(kind=ENTIER), dimension(mesh%n_elems), intent(out), optional :: p_limited

    integer(kind=ENTIER) :: i, j, id_sub_face, id_sub_elem, id_face, id_elem, idl, idr
    real(kind=DOUBLE) :: omega_p, pl, pr, area, w_p, norm_Sp, p0
    real(kind=DOUBLE), parameter :: eps_weno = tiny(1.0_DOUBLE)
    real(kind=DOUBLE), dimension(3) :: S_p, vl_vec, vr_vec, dv_vec, n, dx
    real(kind=DOUBLE), dimension(3, 3) :: S_p_gv, mat, mat_inv
    real(kind=DOUBLE), dimension(:), allocatable :: sum_omega_p, sum_omega_gv

    grad_v = 0.0_DOUBLE
    grad_p = 0.0_DOUBLE
    div_v  = 0.0_DOUBLE
    if (present(p_limited)) p_limited = 0

    allocate(sum_omega_p(mesh%n_elems), sum_omega_gv(mesh%n_elems))
    sum_omega_p  = 0.0_DOUBLE
    sum_omega_gv = 0.0_DOUBLE

    do i = 1, mesh%n_vert
      omega_p = mesh%vert(i)%volume
      S_p    = 0.0_DOUBLE
      S_p_gv = 0.0_DOUBLE
      mat    = 0.0_DOUBLE

      do j = 1, mesh%vert(i)%n_sub_faces_neigh
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        idl    = mesh%sub_face(id_sub_face)%left_elem_neigh
        idr    = mesh%sub_face(id_sub_face)%right_elem_neigh
        area   = mesh%sub_face(id_sub_face)%area
        n      = mesh%sub_face(id_sub_face)%norm
        pl     = pressure(sol(:, idl), gamma_arr(idl))
        vl_vec = sol(2:4, idl)

        if (idr > 0) then
          pr     = pressure(sol(:, idr), gamma_arr(idr))
          vr_vec = sol(2:4, idr)
          mat = mat + tensor_product_3(area * n, &
            mesh%elem(idr)%coord - mesh%elem(idl)%coord)
        else
          pr     = pressure(sol(:, idl), gamma_arr(idl))
          vr_vec = sol(2:4, idl) - 2.0_DOUBLE*dot_product(sol(2:4, idl), n)*n
          mat = mat + tensor_product_3(area * n, &
            2.0_DOUBLE*(mesh%face(id_face)%coord - mesh%elem(idl)%coord))
        end if

        dv_vec = vr_vec - vl_vec
        S_p    = S_p    + area * n * (pr - pl)
        S_p_gv = S_p_gv + area * tensor_product_3(dv_vec, n)
      end do

      S_p    = S_p    / omega_p
      S_p_gv = S_p_gv / omega_p
      norm_Sp = norm2(S_p)
      w_p = omega_p / (eps_weno + norm_Sp**4)
      !w_p = omega_p

      do j = 1, mesh%vert(i)%n_sub_elems_neigh
        id_sub_elem = mesh%vert(i)%sub_elem_neigh(j)
        id_elem     = mesh%sub_elem(id_sub_elem)%mesh_elem
        grad_p(:, id_elem)    = grad_p(:, id_elem)    + w_p * S_p
        grad_v(:, :, id_elem) = grad_v(:, :, id_elem) + w_p * S_p_gv
        sum_omega_p(id_elem)  = sum_omega_p(id_elem)  + w_p
        sum_omega_gv(id_elem) = sum_omega_gv(id_elem) + w_p
      end do
    end do

    do i = 1, mesh%n_elems
      if (sum_omega_p(i)  > 0.0_DOUBLE) grad_p(:, i)    = grad_p(:, i)    / sum_omega_p(i)
      if (sum_omega_gv(i) > 0.0_DOUBLE) grad_v(:, :, i) = grad_v(:, :, i) / sum_omega_gv(i)
      div_v(i) = grad_v(1, 1, i) + grad_v(2, 2, i) + grad_v(3, 3, i)
    end do

    deallocate(sum_omega_p, sum_omega_gv)

    ! Positivity limiter: if the full reconstruction of p at any vertex gives p ≤ 0,
    ! zero out the gradients so that reconstruction falls back to cell-centred p > 0.
    do i = 1, mesh%n_vert
      do j = 1, mesh%vert(i)%n_sub_elems_neigh
        id_sub_elem = mesh%vert(i)%sub_elem_neigh(j)
        id_elem     = mesh%sub_elem(id_sub_elem)%mesh_elem
        if (mesh%elem(id_elem)%is_ghost) cycle
        dx = mesh%vert(i)%coord - mesh%elem(id_elem)%coord
        p0 = pressure(sol(:, id_elem), gamma_arr(id_elem))
        ! Full reconstruction: p0 + grad_p·dx - 0.5*dt*(γ*p0)*div_v
        if (p0 + dot_product(grad_p(:, id_elem), dx) &
            - 0.5_DOUBLE*dt*gamma_arr(id_elem)*p0*div_v(id_elem) <= 0.0_DOUBLE) then
          grad_p(:, id_elem)    = 0.0_DOUBLE
          grad_v(:, :, id_elem) = 0.0_DOUBLE
          div_v(id_elem)        = 0.0_DOUBLE
          if (present(p_limited)) p_limited(id_elem) = 1
        end if
      end do
    end do
  end subroutine compute_gradients

  subroutine compute_rhs_lagrange_classic_iso(mesh, sol, vp, dt, rhs, n_bc, bc_type, &
      bc_val, b2d, mass, gamma_arr, vp_is_imposed, iso_weight_mode, h_extrude)
    ! "classic_iso" family: same node-based nodal-velocity solve as the classic
    ! scheme, but every INTERNAL in-plane sub-face flux is blended with a 1D
    ! acoustic face flux:
    !
    !     F = w_pcf * F_classic + (1 - w_pcf) * F_face
    !
    ! - w_pcf(j) = scale_a / area_j  (clipped to [0,1]) over the sub-faces around
    !   the node, in the spirit of the Euler scheme multi_point_iso.
    !     iso_weight_mode = 1  ("classic_iso")     : scale_a = min_k area_k
    !     iso_weight_mode = 2  ("classic_iso_vol") : scale_a = (min_k V_sub,k)^(2/3)
    !       with V_sub,k the sub-element volumes around the vertex.
    ! - F_face uses classic 1D Lagrangian star values (p_f, v_f) from an acoustic
    !   Riemann solver, as done for the Lagrange part of the LS1D schemes.
    ! Boundary sub-faces and, in boundary_2d, the top/bottom z-faces keep the pure
    ! classic flux (w_pcf = 1, unchanged BC handling). The nodal system uses the
    ! matching weight A_j * w_pcf(j) so it stays the discrete adjoint of the flux.
    use linear_solver_module, only: lu_solve, print_mat, inverse_3_by_3
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), intent(in) :: dt
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mass
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    logical, dimension(mesh%n_vert), intent(in) :: vp_is_imposed
    integer(kind=ENTIER), intent(in) :: n_bc
    character(len=255), dimension(n_bc) :: bc_type
    real(kind=DOUBLE), dimension(5, n_bc) :: bc_val
    logical, intent(in) :: b2d
    integer(kind=ENTIER), intent(in) :: iso_weight_mode
    real(kind=DOUBLE), intent(in) :: h_extrude   ! extrusion height (b2d meshes)

    integer(kind=ENTIER) :: i, j, id_face, id_sub_face
    integer(kind=ENTIER) :: nsfn, idl, idr, var
    real(kind=DOUBLE), dimension(:, :), allocatable :: lambda
    real(kind=DOUBLE), dimension(:, :, :), allocatable :: sol_lr
    real(kind=DOUBLE), dimension(:), allocatable :: w_pcf
    real(kind=DOUBLE), dimension(5) :: flux, flux_c, flux_f
    real(kind=DOUBLE), dimension(3, 3) :: Mp, Mp_inv
    real(kind=DOUBLE), dimension(3) :: Rp, Bp
    real(kind=DOUBLE) :: vl, vr, pl, pr, v_bar, pl_et, pr_et, v_et, Pp
    real(kind=DOUBLE) :: gl, gr
    real(kind=DOUBLE) :: tau_l, tau_r, a_l, a_r, rho_m, a_m, p_f, v_f, wj, scale_a, aw
    real(kind=DOUBLE) :: p_w
    integer(kind=ENTIER) :: n_ip
    logical :: is_corner, have_first_wall_norm
    real(kind=DOUBLE), dimension(3) :: first_wall_norm
    real(kind=DOUBLE) :: cross_z

    rhs = 0.0_DOUBLE
    do i=1, mesh%n_vert
      if (.not. vp_is_imposed(i)) vp(:, i) = 0.0_DOUBLE
    end do

    do i=1, mesh%n_vert
      nsfn = mesh%vert(i)%n_sub_faces_neigh
      allocate(lambda(2, nsfn))
      lambda = 0.0_DOUBLE
      allocate(sol_lr(5, 2, nsfn))
      sol_lr = 0.0_DOUBLE
      allocate(w_pcf(nsfn))

      ! ---- Sub-face weight w_pcf(j) = min(1, (scale_a / A_j)^p_w) --------------
      ! Recipes for the reference area scale scale_a (and exponent p_w), chosen
      ! by iso_weight_mode.  In boundary_2d the top/bottom z-faces are excluded
      ! from every statistic and given w_pcf = 1 (pure classic).
      !   1  min_k A_k                                (p=1)  "classic_iso"
      !   2  b2d: sqrt(A_2D_min)*h ; 3D: V_min^(2/3)   (p=1)  "classic_iso_vol"
      !   3  arithmetic mean of A_k                   (p=1)
      !   4  min_k A_k                                (p=1/2)
      !   5  (min_k A_k / max_k A_k)                  (p=1)   node aspect ratio
      !   6  harmonic mean of A_k                     (p=1)
      !   7  b2d: sqrt(A_2D_min)*h ; 3D: V_min^(2/3)   (p=1/2)
      p_w = 1.0_DOUBLE
      if (iso_weight_mode == 2 .or. iso_weight_mode == 7) then
        scale_a = huge(1.0_DOUBLE)
        do j = 1, mesh%vert(i)%n_sub_elems_neigh
          scale_a = min(scale_a, &
            abs(mesh%sub_elem(mesh%vert(i)%sub_elem_neigh(j))%volume))
        end do
        if (b2d) then
          scale_a = sqrt(scale_a / h_extrude) * h_extrude
        else
          scale_a = scale_a ** (2.0_DOUBLE/3.0_DOUBLE)
        end if
        if (iso_weight_mode == 7) p_w = 0.5_DOUBLE
      else
        w_pcf = huge(1.0_DOUBLE)
        n_ip = 0
        do j=1, nsfn
          id_sub_face = mesh%vert(i)%sub_face_neigh(j)
          id_face = mesh%sub_face(id_sub_face)%mesh_face
          if( b2d .and. abs(mesh%face(id_face)%norm(3)) > 1e-8_DOUBLE ) cycle
          w_pcf(j) = mesh%sub_face(id_sub_face)%area
          n_ip = n_ip + 1
        end do
        select case (iso_weight_mode)
        case (3)
          scale_a = sum(w_pcf, mask = (w_pcf < huge(1.0_DOUBLE))) / max(1, n_ip)
        case (5)
          scale_a = minval(w_pcf) &
            / maxval(w_pcf, mask = (w_pcf < huge(1.0_DOUBLE)))
        case (6)
          scale_a = real(max(1, n_ip), DOUBLE) &
            / sum(1.0_DOUBLE / w_pcf, mask = (w_pcf < huge(1.0_DOUBLE)))
        case (4)
          scale_a = minval(w_pcf)
          p_w = 0.5_DOUBLE
        case default   ! 1
          scale_a = minval(w_pcf)
        end select
      end if
      do j=1, nsfn
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        idr = mesh%face(id_face)%right_neigh
        if( idr > 0 .and. .not. (b2d .and. &
            abs(mesh%face(id_face)%norm(3)) > 1e-8_DOUBLE) ) then
          if (iso_weight_mode == 5) then
            w_pcf(j) = min(1.0_DOUBLE, scale_a)
          else
            w_pcf(j) = min(1.0_DOUBLE, &
              (scale_a / mesh%sub_face(id_sub_face)%area) ** p_w)
          end if
        else
          ! boundary faces and z-faces: pure classic flux, keep full area weight.
          w_pcf(j) = 1.0_DOUBLE
        end if
      end do

      Bp = 0.0_DOUBLE
      Mp = 0.0_DOUBLE
      Rp = 0.0_DOUBLE
      is_corner = .false.
      have_first_wall_norm = .false.
      first_wall_norm = 0.0_DOUBLE
      do j=1, nsfn
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        idl = mesh%face(id_face)%left_neigh
        idr = mesh%face(id_face)%right_neigh
        gl = gamma_arr(idl)
        sol_lr(:, 1, j) = sol(:, idl)
        pl = pressure(sol_lr(:, 1, j), gl)
        vl = dot_product(sol_lr(2:4, 1, j), mesh%face(id_face)%norm)
        ! The nodal system is the discrete adjoint of the ACTUAL sub-face flux.
        ! The classic (multipoint) part of that flux carries the factor
        ! A_j * w_pcf(j); the (1-w_pcf) face-flux part cancels between left and
        ! right and drops out. So every A_j below is replaced by A_j * w_pcf(j)
        ! -> on internal in-plane faces the area cancels (constant min_area),
        ! on boundary / z-faces w_pcf(j) = 1 and the area weight is unchanged.
        aw = mesh%sub_face(id_sub_face)%area * w_pcf(j)

        if( idr > 0 ) then
          gr = gamma_arr(idr)
          sol_lr(:, 2, j) = sol(:, idr)
          pr = pressure(sol_lr(:, 2, j), gr)
          vr = dot_product(sol_lr(2:4, 2, j), mesh%face(id_face)%norm)

          lambda(1, j) = max(sqrt(gl*pl*sol_lr(1, 1, j))/sol_lr(1, 1, j), &
            sqrt(max(0.0_DOUBLE, pr-pl)/sol_lr(1, 1, j)), &
            -(vr-vl)/sol_lr(1, 1, j))
          lambda(2, j) = max(sqrt(gr*pr*sol_lr(1, 2, j))/sol_lr(1, 2, j), &
            sqrt(max(0.0_DOUBLE, pl-pr)/sol_lr(1, 2, j)), &
            -(vr-vl)/sol_lr(1, 2, j))

          if (.not. vp_is_imposed(i)) then
            v_bar = (lambda(1, j)*vl + lambda(2, j)*vr - (pr - pl)) &
              /(lambda(2, j) + lambda(1, j))
            Mp = Mp + aw*(lambda(1,j)+lambda(2,j)) &
              * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
            Rp = Rp + aw*(lambda(1,j)+lambda(2,j)) &
              *mesh%face(id_face)%norm*v_bar
          end if

        else
          if( b2d .and. abs(mesh%face(id_face)%norm(3)) > 1e-8_DOUBLE ) then
            lambda(1, j) = sqrt(gl*pl*sol_lr(1, 1, j))/sol_lr(1, 1, j)
            if (.not. vp_is_imposed(i)) then
              Mp = Mp + aw*2*lambda(1,j) &
                * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
            end if
          else if (idr < 0 .and. trim(bc_type(-idr)) == 'pressure') then
            lambda(1, j) = sqrt(gl*pl*sol_lr(1, 1, j))/sol_lr(1, 1, j)
            if (.not. vp_is_imposed(i)) then
              Mp = Mp + aw*lambda(1,j) &
                * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
              Rp = Rp + aw &
                *(bc_val(1, -idr) + lambda(1,j)*vl) &
                *mesh%face(id_face)%norm
            end if
          else
            lambda(1, j) = sqrt(gl*pl*sol_lr(1, 1, j))/sol_lr(1, 1, j)
            if (.not. vp_is_imposed(i)) then
              if (.not. have_first_wall_norm) then
                first_wall_norm = mesh%face(id_face)%norm
                have_first_wall_norm = .true.
              else
                cross_z = first_wall_norm(1)*mesh%face(id_face)%norm(2) &
                        - first_wall_norm(2)*mesh%face(id_face)%norm(1)
                if (abs(cross_z) > 1.0e-8_DOUBLE) is_corner = .true.
              end if
              Bp = Bp + aw*mesh%face(id_face)%norm
              Mp = Mp + aw*lambda(1,j) &
                * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
              Rp = Rp + aw*(pl + lambda(1,j)*vl) &
                *mesh%face(id_face)%norm
            end if
          end if
        end if
      end do

      if (.not. vp_is_imposed(i)) then
        if (is_corner) then
          vp(:, i) = 0.0_DOUBLE
        else if( maxval(abs(Bp)) > 1e-8_DOUBLE ) then
          call inverse_3_by_3(Mp, Mp_inv)
          Pp = dot_product(Rp, matmul(Mp_inv,Bp))/dot_product(Bp, matmul(Mp_inv, Bp))
          vp(:, i) = matmul(Mp_inv, Rp-Pp*BP)
        else
          call inverse_3_by_3(Mp, Mp_inv)
          vp(:, i) = matmul(Mp_inv, Rp)
        end if
      end if

      do j=1, nsfn
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        idl = mesh%face(id_face)%left_neigh
        idr = mesh%face(id_face)%right_neigh
        gl = gamma_arr(idl)

        v_et = dot_product(vp(:, i), mesh%face(id_face)%norm)
        vl = dot_product(sol_lr(2:4, 1, j), mesh%face(id_face)%norm)
        pl = pressure(sol_lr(:, 1, j), gl)
        pl_et = pl - lambda(1, j)*(v_et - vl)

        ! Blend the 1D acoustic face flux only on INTERNAL IN-PLANE sub-faces.
        ! Boundary faces and, in boundary_2d, the top/bottom z-faces stay pure
        ! classic (wj = 1), so their handling is identical to compute_rhs_lagrange.
        wj = 1.0_DOUBLE
        flux_f = 0.0_DOUBLE
        if( idr > 0 .and. .not. (b2d .and. &
            abs(mesh%face(id_face)%norm(3)) > 1e-8_DOUBLE) ) then
          wj = w_pcf(j)
          gr = gamma_arr(idr)
          vr = dot_product(sol_lr(2:4, 2, j), mesh%face(id_face)%norm)
          pr = pressure(sol_lr(:, 2, j), gr)
          tau_l = sol_lr(1, 1, j)
          tau_r = sol_lr(1, 2, j)
          a_l = sqrt(max(0.0_DOUBLE, gl*pl*tau_l))
          a_r = sqrt(max(0.0_DOUBLE, gr*pr*tau_r))
          rho_m = 0.5_DOUBLE*(1.0_DOUBLE/tau_l + 1.0_DOUBLE/tau_r)
          a_m   = max(a_l, a_r, 1.0e-10_DOUBLE)
          v_f = 0.5_DOUBLE*(vl + vr) - (pr - pl)/(2.0_DOUBLE*rho_m*a_m)
          p_f = 0.5_DOUBLE*(pl + pr) - 0.5_DOUBLE*rho_m*a_m*(vr - vl)
          flux_f(1)   = -v_f
          flux_f(2:4) = p_f*mesh%face(id_face)%norm
          flux_f(5)   = p_f*v_f
        end if

        ! left cell
        flux_c(1)   = -v_et
        flux_c(2:4) = pl_et*mesh%face(id_face)%norm
        flux_c(5)   = pl_et*v_et
        flux = wj*flux_c + (1.0_DOUBLE - wj)*flux_f
        do var=1, 5
          rhs(var, idl) = rhs(var, idl) &
            - mesh%sub_face(id_sub_face)%area/mass(idl) * flux(var)
        end do

        if( idr > 0 ) then
          gr = gamma_arr(idr)
          vr = dot_product(sol_lr(2:4, 2, j), mesh%face(id_face)%norm)
          pr = pressure(sol_lr(:, 2, j), gr)
          pr_et = pr + lambda(2, j)*(v_et - vr)
          flux_c(1)   = -v_et
          flux_c(2:4) = pr_et*mesh%face(id_face)%norm
          flux_c(5)   = pr_et*v_et
          flux = wj*flux_c + (1.0_DOUBLE - wj)*flux_f
          do var=1, 5
            rhs(var, idr) = rhs(var, idr) &
              + mesh%sub_face(id_sub_face)%area/mass(idr) * flux(var)
          end do
        end if
      end do

      deallocate(lambda, sol_lr, w_pcf)
    end do
  end subroutine compute_rhs_lagrange_classic_iso

  subroutine compute_rhs_lagrange_sidil(mesh, sol, vp, dt, rhs, &
      n_bc, bc_type, bc_val, b2d, mass, method, b2d_h, gamma_arr, vp_is_imposed)
    use linear_solver_module, only: lu_solve, print_mat, inverse_3_by_3
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), intent(in) :: dt
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mass
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    logical, dimension(mesh%n_vert), intent(in) :: vp_is_imposed
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
    do i=1, mesh%n_vert
      if (.not. vp_is_imposed(i)) vp(:, i) = 0.0_DOUBLE
    end do

    do i=1, mesh%n_vert

      if (.not. vp_is_imposed(i)) then
        call compute_nodal_velocity_sidil(mesh, i, sol, vp(:, i), method, b2d, b2d_h, gamma_arr)
      end if
      call compute_nodal_pressure_sidil(mesh, i, sol, pp, method, b2d, b2d_h, gamma_arr)

      a_p = 0.0_DOUBLE
      do j=1, mesh%vert(i)%n_sub_elems_neigh
        id_sub_elem = mesh%vert(i)%sub_elem_neigh(j)
        id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
        sol_w = lag_to_primit(sol(:, id_elem), gamma_arr(id_elem))
        a_p = a_p + mesh%sub_elem(id_sub_elem)%volume &
          * sqrt(gamma_arr(id_elem)*sol_w(5)/sol_w(1))
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

  subroutine compute_dt(mesh, sol, dt, cfl, vp, me, num_procs, gamma_arr)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: vp
    real(kind=DOUBLE), intent(in) :: cfl
    real(kind=DOUBLE), intent(inout) :: dt
    integer(kind=ENTIER) :: me, num_procs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr

    real(kind=DOUBLE), parameter :: CV = 0.8_DOUBLE
    integer(kind=ENTIER) :: i, j, k, mpi_ierr
    integer(kind=ENTIER) :: id_sub_face, id_face
    integer(kind=ENTIER) :: id_vert, id_sub_elem
    real(kind=DOUBLE) :: dx, c, avg_sub_vol
    real(kind=DOUBLE), dimension(5) :: w
    real(kind=DOUBLE), dimension(3) :: norm

    dt = 1e6
    do i=1, mesh%n_elems
      c = sqrt(max(0.0_DOUBLE, gamma_arr(i)*pressure(sol(:, i), gamma_arr(i))*sol(1, i)))
      c = max(c, 1.0e-10_DOUBLE)
      ! Use average sub-element volume: this prevents dt→0 when individual sub-elems
      ! shrink due to cell distortion, while matching the reference scheme's dt scale.
      avg_sub_vol = abs(mesh%elem(i)%volume) / real(mesh%elem(i)%n_sub_elems, DOUBLE)
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
          dt = min(dt, avg_sub_vol/(mesh%sub_face(id_sub_face)%area*c))
          dt = min(dt, CV*avg_sub_vol &
            /(1e-8_DOUBLE + abs(dot_product(vp(:, id_vert), norm))&
            *mesh%face(id_face)%area))
        end do
      end do
    end do
    dt = cfl * dt

    if( num_procs > 1 ) call MPI_ALLREDUCE(MPI_IN_PLACE, &
      dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, mpi_ierr)
  end subroutine compute_dt

  subroutine init_sol(mesh, sol, sol_uniform, init, me, num_procs, gamma_arr, boundary_2d)
    use mpi
    use mpi_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: sol
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_uniform
    integer(kind=ENTIER), intent(in) :: init
    integer(kind=ENTIER), intent(in) :: me, num_procs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(inout) :: gamma_arr
    logical, intent(in) :: boundary_2d

    integer(kind=ENTIER) :: i, imin, mpi_ierr
    real(kind=DOUBLE) :: tot_vol, r, rmin, rmin_share
    !real(kind=DOUBLE), parameter :: r_sedov = 2.4_DOUBLE/20.0_DOUBLE
    real(kind=DOUBLE), parameter :: r_sedov = 0.011_DOUBLE
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
        sol(5, i) = sol(5, i)/((gamma_arr(i) - 1.0_DOUBLE)/sol(1, i)) &
          + 0.5_DOUBLE*norm2(sol(2:4, i))**2
      end do
    else if( init == 1) then
      ! Find the cell nearest to the origin (same as Euler init_sedov).
      rmin = huge(1.0_DOUBLE)
      imin = -1
      do i = 1, mesh%n_elems
        if (mesh%elem(i)%is_ghost) cycle
        if (boundary_2d) then
          if (norm2(mesh%elem(i)%coord(:2)) < rmin) then
            rmin = norm2(mesh%elem(i)%coord(:2))
            imin = i
          end if
        else
          if (norm2(mesh%elem(i)%coord) < rmin) then
            rmin = norm2(mesh%elem(i)%coord)
            imin = i
          end if
        end if
      end do


      if (num_procs > 1) then
        rmin_share = rmin
        call MPI_ALLREDUCE(MPI_IN_PLACE, rmin_share, 1, MPI_DOUBLE, &
          MPI_MIN, MPI_COMM_WORLD, mpi_ierr)
        if (abs(rmin - rmin_share) > 1e-14_DOUBLE) imin = -1
      end if

      do i = 1, mesh%n_elems
        sol(1, i) = 1.0_DOUBLE
        sol(2:4, i) = 0.0_DOUBLE
        sol(5, i) = 1e-12_DOUBLE/(gamma_arr(i) - 1.0_DOUBLE)
      end do

      if (imin > 0) then
        if (boundary_2d) then
          sol(5, imin) = 0.984042_DOUBLE/mesh%elem(imin)%volume
        else
          sol(5, imin) = 0.851072_DOUBLE/mesh%elem(imin)%volume
        end if
      end if
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
          sol(5, i) = 1e-8_DOUBLE/(gamma_arr(i) - 1.0_DOUBLE)
        end if
      end do
    else if( init == 2) then
      do i=1, mesh%n_elems
        call sol_isentropic_vortex(mesh%elem(i)%coord, sol(:, i), 0.0_DOUBLE)
        sol(5, i) = sol(5, i)/(sol(1, i)*(gamma_arr(i)-1.0_DOUBLE)) &
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
        sol(5, i) = 10/(gamma_arr(i)-1.0_DOUBLE)
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
        sol(5, i) = sol(5, i)/((gamma_arr(i) - 1.0_DOUBLE)/sol(1, i)) &
          + 0.5_DOUBLE*norm2(sol(2:4, i))**2
      end do
    else if( init == 3 ) then
      ! Noh 2D cylindrical: uniform density rho=1, radial inward unit velocity, p~0
      do i=1, mesh%n_elems
        sol(1, i) = 1.0_DOUBLE
        r = sqrt(mesh%elem(i)%coord(1)**2 + mesh%elem(i)%coord(2)**2)
        if (r < 1.0e-10_DOUBLE) then
          sol(2, i) = 0.0_DOUBLE
          sol(3, i) = 0.0_DOUBLE
        else
          sol(2, i) = -mesh%elem(i)%coord(1) / r
          sol(3, i) = -mesh%elem(i)%coord(2) / r
        end if
        sol(4, i) = 0.0_DOUBLE
        ! E = p/((gamma-1)*rho) + 0.5*|v|^2, with p=1e-6 and |v|=1 away from center
        sol(5, i) = 1.0e-6_DOUBLE / (gamma_arr(i) - 1.0_DOUBLE) &
          + 0.5_DOUBLE * (sol(2, i)**2 + sol(3, i)**2)
      end do
    else if( init == 5 ) then
      ! Saltzman piston: uniform gas at rest (rho=1, p~0); piston BC at left wall sets v=(1,0,0)
      do i=1, mesh%n_elems
        sol(1, i) = 1.0_DOUBLE
        sol(2:4, i) = 0.0_DOUBLE
        sol(5, i) = 1.0e-6_DOUBLE / (gamma_arr(i) - 1.0_DOUBLE)
      end do
    else if( init == 6 ) then
      ! Sod bi-materiaux: left (x<0.5) gamma=1.4, right (x>=0.5) gamma=1.5
      do i=1, mesh%n_elems
        if( mesh%elem(i)%coord(1) < 0.5_DOUBLE ) then
          gamma_arr(i) = 1.4_DOUBLE
          sol(1, i) = 1.0_DOUBLE
          sol(2:4, i) = 0.0_DOUBLE
          sol(5, i) = 1.0_DOUBLE
        else
          gamma_arr(i) = 1.5_DOUBLE
          sol(1, i) = 0.125_DOUBLE
          sol(2:4, i) = 0.0_DOUBLE
          sol(5, i) = 0.1_DOUBLE
        end if
        sol(1, i) = 1.0_DOUBLE/sol(1, i)
        sol(5, i) = sol(5, i)/((gamma_arr(i) - 1.0_DOUBLE)/sol(1, i)) &
          + 0.5_DOUBLE*norm2(sol(2:4, i))**2
      end do
    else if( init == 7 ) then
      ! Triple point 3D (LLNL version): three states, two materials
      ! Left  (x<1):                     rho=1,   p=1,   gamma=1.5
      ! Right, outside cyl (r_yz>1.5):  rho=0.1, p=0.1, gamma=1.5
      ! Right, inside  cyl (r_yz<=1.5): rho=1,   p=0.1, gamma=1.4
      ! where r_yz = sqrt(y^2 + z^2) is the radius from the x-axis
      do i=1, mesh%n_elems
        r = sqrt(mesh%elem(i)%coord(2)**2 + mesh%elem(i)%coord(3)**2)
        if (mesh%elem(i)%coord(1) < 1.0_DOUBLE) then
          gamma_arr(i) = 1.5_DOUBLE
          sol(1, i)    = 1.0_DOUBLE
          sol(5, i)    = 1.0_DOUBLE
        else if (r > 1.5_DOUBLE) then
          gamma_arr(i) = 1.5_DOUBLE
          sol(1, i)    = 0.1_DOUBLE
          sol(5, i)    = 0.1_DOUBLE
        else
          gamma_arr(i) = 1.4_DOUBLE
          sol(1, i)    = 1.0_DOUBLE
          sol(5, i)    = 0.1_DOUBLE
        end if
        sol(2:4, i) = 0.0_DOUBLE
        sol(1, i) = 1.0_DOUBLE / sol(1, i)
        sol(5, i) = sol(5, i) / ((gamma_arr(i) - 1.0_DOUBLE) / sol(1, i)) &
          + 0.5_DOUBLE * norm2(sol(2:4, i))**2
      end do
    else if( init == 8 ) then
      ! Gradient test: linear pressure p = 1 + x + y, zero velocity, rho = 1.
      ! Exact grad_p = [1, 1, 0]; exact grad_v = 0; exact div_v = 0.
      do i=1, mesh%n_elems
        sol(1, i) = 1.0_DOUBLE
        sol(2:4, i) = 0.0_DOUBLE
        sol(5, i) = (1.0_DOUBLE + mesh%elem(i)%coord(1) + mesh%elem(i)%coord(2)) &
                    / (gamma_arr(i) - 1.0_DOUBLE)
      end do
    end if

  end subroutine init_sol

  function pressure(u, g)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: u
    real(kind=DOUBLE), intent(in) :: g
    real(kind=DOUBLE) :: pressure

    pressure = (g-1.0_DOUBLE)/u(1)*(u(5) - 0.5_DOUBLE*norm2(u(2:4))**2)
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

  subroutine compute_nodal_velocity_sidil(mesh, id_vert, sol, vp, method, b2d, b2d_h, gamma_arr)
    use lagrange_global_data_module, only : boundary_2d
    use linear_solver_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert, method
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3), intent(inout) :: vp
    real(kind=DOUBLE), intent(in) :: b2d_h
    logical, intent(in) :: b2d
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr

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
      sol_w = lag_to_primit(sol(:, id_elem), gamma_arr(id_elem))
      if( sol_w(5) < 0.0_DOUBLE ) then
        print*,"Error pos", id_elem, mesh%elem(id_elem)%coord
      end if
      vp = vp &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(2:4)
      rho_p = rho_p &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(1)
      a_p = a_p &
        + mesh%sub_elem(id_sub_elem)%volume &
        * sqrt(gamma_arr(id_elem)*sol_w(5)/sol_w(1))
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
      sol_w_l = lag_to_primit(sol_l, gamma_arr(le))
      sol_w_r = lag_to_primit(sol_r, gamma_arr(le))
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

  subroutine compute_nodal_pressure_sidil(mesh, id_vert, sol, pp, method, b2d, b2d_h, gamma_arr)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert, method
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), intent(inout) :: pp
    real(kind=DOUBLE), intent(in) :: b2d_h
    logical :: b2d
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr

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
      sol_w = lag_to_primit(sol(:, id_elem), gamma_arr(id_elem))
      pp = pp &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(5)
      rho_p = rho_p &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(1)
      a_p = a_p &
        + mesh%sub_elem(id_sub_elem)%volume &
        * sqrt(gamma_arr(id_elem)*sol_w(5)/sol_w(1))
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
      sol_w_l = lag_to_primit(sol_l, gamma_arr(le))
      sol_w_r = lag_to_primit(sol_r, gamma_arr(le))
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

  function lag_to_primit(u, g) result(w)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: u
    real(kind=DOUBLE), intent(in) :: g
    real(kind=DOUBLE), dimension(5) :: w

    w(1) = 1.0_DOUBLE/u(1)
    w(2:4) = u(2:4)
    w(5) = pressure(u, g)
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
  
  function omega_pos(V, deltaV)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: V, deltaV
    real(kind=DOUBLE) :: omega_pos

    real(kind=DOUBLE) :: newEps
    real(kind=DOUBLE), dimension(5) :: newV

    if (any(deltaV /= deltaV)) then
      omega_pos = 0.0_DOUBLE
      return
    end if

    omega_pos = 1.0_DOUBLE
    newV = V + omega_pos * deltaV
    newEps = newV(5) - 0.5_DOUBLE*dot_product(newV(2:4),newV(2:4))
    do while (.not. (newEps >= 0.0_DOUBLE) .and. omega_pos > 1e-14_DOUBLE)
      omega_pos = 0.9_DOUBLE*omega_pos
      newV = V + omega_pos * deltaV
      newEps = newV(5) - 0.5_DOUBLE*dot_product(newV(2:4),newV(2:4))
    end do
    if (.not. (newEps >= 0.0_DOUBLE)) omega_pos = 0.0_DOUBLE
  end function omega_pos
end module lagrange_module
