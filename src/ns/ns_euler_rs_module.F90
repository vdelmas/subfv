module ns_euler_rs_module
  use precision_module
  use mesh_module
  use ns_euler_primitives_module
  use ns_euler_recon_module
  implicit none

contains
  subroutine three_wave(sol_w_l, sol_w_r, lr_flux, sl, sr)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(5) :: sol_l, sol_r
    real(kind=DOUBLE), dimension(5, 2), intent(inout) :: lr_flux
    real(kind=DOUBLE), intent(inout) :: sl, sr

    real(kind=DOUBLE) :: rhol, rhor, ul, vl, wl, ur, vr, wr
    real(kind=DOUBLE) :: pl, pr, rhol_et, rhor_et
    real(kind=DOUBLE) :: ul_et, vl_et, wl_et, ur_et, vr_et, wr_et, u_et
    real(kind=DOUBLE) :: el, er, al, ar
    real(kind=DOUBLE) :: u_bar, pl_bar, pr_bar
    real(kind=DOUBLE) :: lambda_l, lambda_r
    real(kind=DOUBLE), dimension(5) :: fl, fr
    real(kind=DOUBLE), dimension(5) :: sol_l_et, sol_r_et

    rhol = sol_w_l(1)
    ul = sol_w_l(2)
    vl = sol_w_l(3)
    wl = sol_w_l(4)
    pl = sol_w_l(5)
    sol_l = primit_to_conserv(sol_w_l)
    sol_r = primit_to_conserv(sol_w_r)
    el = sol_l(5)/rhol
    al = sound_speed_w(sol_w_l)

    rhor = sol_w_r(1)
    ur = sol_w_r(2)
    vr = sol_w_r(3)
    wr = sol_w_r(4)
    pr = sol_w_r(5)
    er = sol_r(5)/rhor
    ar = sound_speed_w(sol_w_r)

    fl(:) = ul*sol_l(:) + (/0.0_DOUBLE, pl, 0.0_DOUBLE, 0.0_DOUBLE, pl*ul/)
    fr(:) = ur*sol_r(:) + (/0.0_DOUBLE, pr, 0.0_DOUBLE, 0.0_DOUBLE, pr*ur/)

    lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(ur - ul))
    lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(ur - ul))
    u_bar = (lambda_l*ul + lambda_r*ur - (pr - pl))/(lambda_r + lambda_l)

    u_et = u_bar

    rhol_et = 1.0_DOUBLE/(1.0_DOUBLE/rhol + (u_et - ul)/lambda_l)
    ul_et = u_et
    vl_et = vl
    wl_et = wl
    pl_bar = pl - lambda_l*(u_et - ul)

    sol_l_et(1) = rhol_et
    sol_l_et(2) = rhol_et*ul_et
    sol_l_et(3) = rhol_et*vl_et
    sol_l_et(4) = rhol_et*wl_et
    sol_l_et(5) = rhol_et*(el + (pl*ul - pl_bar*ul_et)/lambda_l)

    rhor_et = 1.0_DOUBLE/(1.0_DOUBLE/rhor + (ur - u_et)/lambda_r)
    ur_et = u_et
    vr_et = vr
    wr_et = wr
    pr_bar = pr + lambda_r*(u_et - ur)

    sol_r_et(1) = rhor_et
    sol_r_et(2) = rhor_et*ur_et
    sol_r_et(3) = rhor_et*vr_et
    sol_r_et(4) = rhor_et*wr_et
    sol_r_et(5) = rhor_et*(er + (pr_bar*ur_et - pr*ur)/lambda_r)

    sl = ul - lambda_l/rhol
    sr = ur + lambda_r/rhor

    lr_flux(:, 1) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(u_et)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et))

    lr_flux(:, 2) = -lr_flux(:, 1)
  end subroutine three_wave

  subroutine modified_three_wave(sol_w_l, sol_w_r, lr_flux, sl, sr)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(5) :: sol_l, sol_r
    real(kind=DOUBLE), dimension(5, 2), intent(inout) :: lr_flux
    real(kind=DOUBLE), intent(inout) :: sl, sr

    real(kind=DOUBLE) :: rhol, rhor, ul, vl, wl, ur, vr, wr
    real(kind=DOUBLE) :: pl, pr, rhol_et, rhor_et
    real(kind=DOUBLE) :: ul_et, vl_et, wl_et, ur_et, vr_et, wr_et, u_et
    real(kind=DOUBLE) :: el, er, al, ar
    real(kind=DOUBLE) :: piv_et, v_et, piw_et, w_et
    real(kind=DOUBLE) :: u_bar, pl_bar, pr_bar
    real(kind=DOUBLE) :: lambda_l, lambda_r
    real(kind=DOUBLE), dimension(5) :: fl, fr
    real(kind=DOUBLE), dimension(5) :: sol_l_et, sol_r_et

    rhol = sol_w_l(1)
    ul = sol_w_l(2)
    vl = sol_w_l(3)
    wl = sol_w_l(4)
    pl = sol_w_l(5)
    sol_l = primit_to_conserv(sol_w_l)
    sol_r = primit_to_conserv(sol_w_r)
    el = sol_l(5)/rhol
    al = sound_speed_w(sol_w_l)

    rhor = sol_w_r(1)
    ur = sol_w_r(2)
    vr = sol_w_r(3)
    wr = sol_w_r(4)
    pr = sol_w_r(5)
    er = sol_r(5)/rhor
    ar = sound_speed_w(sol_w_r)

    fl(:) = ul*sol_l(:) + (/0.0_DOUBLE, pl, 0.0_DOUBLE, 0.0_DOUBLE, pl*ul/)
    fr(:) = ur*sol_r(:) + (/0.0_DOUBLE, pr, 0.0_DOUBLE, 0.0_DOUBLE, pr*ur/)

    lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(ur - ul))
    lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(ur - ul))
    u_bar = (lambda_l*ul + lambda_r*ur - (pr - pl))/(lambda_r + lambda_l)

    u_et = u_bar

    piv_et = -lambda_l*lambda_r/(lambda_l + lambda_r)*(vr - vl)
    piw_et = -lambda_l*lambda_r/(lambda_l + lambda_r)*(wr - wl)

    v_et = (lambda_l*vl + lambda_r*vr)/(lambda_l + lambda_r)
    w_et = (lambda_l*wl + lambda_r*wr)/(lambda_l + lambda_r)

    rhol_et = 1.0_DOUBLE/(1.0_DOUBLE/rhol + (u_et - ul)/lambda_l)
    ul_et = u_et
    vl_et = v_et
    wl_et = w_et
    pl_bar = pl - lambda_l*(u_et - ul)

    sol_l_et(1) = rhol_et
    sol_l_et(2) = rhol_et*ul_et
    sol_l_et(3) = rhol_et*vl_et
    sol_l_et(4) = rhol_et*wl_et
    sol_l_et(5) = rhol_et*(el + (pl*ul - pl_bar*ul_et &
      - piv_et*v_et - piw_et*w_et)/lambda_l)

    rhor_et = 1.0_DOUBLE/(1.0_DOUBLE/rhor + (ur - u_et)/lambda_r)
    ur_et = u_et
    vr_et = v_et
    wr_et = w_et
    pr_bar = pr + lambda_r*(u_et - ur)

    sol_r_et(1) = rhor_et
    sol_r_et(2) = rhor_et*ur_et
    sol_r_et(3) = rhor_et*vr_et
    sol_r_et(4) = rhor_et*wr_et
    sol_r_et(5) = rhor_et*(er + (pr_bar*ur_et - pr*ur &
      + piv_et*v_et + piw_et*w_et)/lambda_r)

    sl = ul - lambda_l/rhol
    sr = ur + lambda_r/rhor

    lr_flux(:, 1) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(u_et)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et))

    lr_flux(:, 2) = -lr_flux(:, 1)
  end subroutine modified_three_wave

  subroutine two_wave(sol_w_l, sol_w_r, lr_flux, sl, sr)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(5) :: sol_l, sol_r
    real(kind=DOUBLE), dimension(5, 2), intent(inout) :: lr_flux
    real(kind=DOUBLE), intent(inout) :: sl, sr

    real(kind=DOUBLE) :: rhol, rhor, ul, ur, pl, pr, &
      al, ar, lambda_l, lambda_r
    real(kind=DOUBLE), dimension(5) :: fl, fr, sol_et

    rhol = sol_w_l(1)
    ul = sol_w_l(2)
    pl = sol_w_l(5)
    al = sound_speed_w(sol_w_l)
    sol_l = primit_to_conserv(sol_w_l)

    rhor = sol_w_r(1)
    ur = sol_w_r(2)
    pr = sol_w_r(5)
    ar = sound_speed_w(sol_w_r)
    sol_r = primit_to_conserv(sol_w_r)

    fl(:) = ul*sol_l(:) + (/0.0_DOUBLE, pl, 0.0_DOUBLE, 0.0_DOUBLE, pl*ul/)
    fr(:) = ur*sol_r(:) + (/0.0_DOUBLE, pr, 0.0_DOUBLE, 0.0_DOUBLE, pr*ur/)

    lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(ur - ul))
    lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(ur - ul))

    sl = ul - lambda_l/rhol
    sr = ur + lambda_r/rhor

    sol_et = (sr*sol_r - sl*sol_l - (fr - fl))/(sr - sl)

    lr_flux(:, 1) = 0.5_DOUBLE*(fl + fr) &
      - 0.5_DOUBLE*(abs(sl)*(sol_et - sol_l) + abs(sr)*(sol_r - sol_et))
    lr_flux(:, 2) = -lr_flux(:, 1)
  end subroutine two_wave

  subroutine multi_point(sol_w_l, sol_w_r, lr_flux, u_nodal, &
      lambda_l, lambda_r, sl, sr)
    implicit none

    real(kind=DOUBLE), intent(in) :: u_nodal
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(5) :: sol_l, sol_r
    real(kind=DOUBLE), intent(inout) :: lambda_l, lambda_r
    real(kind=DOUBLE), dimension(5, 2), intent(inout) :: lr_flux
    real(kind=DOUBLE), intent(inout) :: sl, sr

    real(kind=DOUBLE) :: rhol, rhor, ul, vl, wl, ur, vr, wr
    real(kind=DOUBLE) :: pl, pr, rhol_et, rhor_et
    real(kind=DOUBLE) :: ul_et, vl_et, wl_et, ur_et, vr_et, wr_et, u_et
    real(kind=DOUBLE) :: el, er
    real(kind=DOUBLE) :: u_bar, pl_bar, pr_bar
    real(kind=DOUBLE), dimension(5) :: fl, fr
    real(kind=DOUBLE), dimension(5) :: sol_l_et, sol_r_et

    rhol = sol_w_l(1)
    ul = sol_w_l(2)
    vl = sol_w_l(3)
    wl = sol_w_l(4)
    pl = sol_w_l(5)
    sol_l = primit_to_conserv(sol_w_l)
    sol_r = primit_to_conserv(sol_w_r)
    el = sol_l(5)/rhol

    rhor = sol_w_r(1)
    ur = sol_w_r(2)
    vr = sol_w_r(3)
    wr = sol_w_r(4)
    pr = sol_w_r(5)
    er = sol_r(5)/rhor

    fl(:) = ul*sol_l(:) + (/0.0_DOUBLE, pl, 0.0_DOUBLE, 0.0_DOUBLE, pl*ul/)
    fr(:) = ur*sol_r(:) + (/0.0_DOUBLE, pr, 0.0_DOUBLE, 0.0_DOUBLE, pr*ur/)

    u_bar = (lambda_l*ul + lambda_r*ur)/(lambda_l + lambda_r) &
      - (pr - pl)/(lambda_r + lambda_l)

    u_et = u_nodal

    rhol_et = 1.0_DOUBLE/(1.0_DOUBLE/rhol + (u_et - ul)/lambda_l)
    ul_et = u_et
    vl_et = vl
    wl_et = wl
    pl_bar = pl - lambda_l*(u_et - ul)

    sol_l_et(1) = rhol_et
    sol_l_et(2) = rhol_et*ul_et
    sol_l_et(3) = rhol_et*vl_et
    sol_l_et(4) = rhol_et*wl_et
    sol_l_et(5) = rhol_et*(el + (pl*ul - pl_bar*u_et)/lambda_l)

    rhor_et = 1.0_DOUBLE/(1.0_DOUBLE/rhor + (ur - u_et)/lambda_r)
    ur_et = u_et
    vr_et = vr
    wr_et = wr
    pr_bar = pr + lambda_r*(u_et - ur)

    sol_r_et(1) = rhor_et
    sol_r_et(2) = rhor_et*ur_et
    sol_r_et(3) = rhor_et*vr_et
    sol_r_et(4) = rhor_et*wr_et
    sol_r_et(5) = rhor_et*(er + (pr_bar*u_et - pr*ur)/lambda_r)

    if (rhol_et < 0.0_DOUBLE .or. rhor_et < 0.0_DOUBLE) then
      print *, "Negative specific volume MPCC !", rhol_et, rhor_et
      print*, rhol, rhor
      print*, lambda_l, lambda_r
      print*, ul_et, vl_et, wl_et
      print*, ur_et, vr_et, wr_et
      print*,"R"
      print*, 1.0_DOUBLE/rhor, (ur - u_et)/lambda_r
      error stop
    end if

    sl = ul - lambda_l/rhol
    sr = ur + lambda_r/rhor

    lr_flux(:, 1) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(u_et)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et)) &
      - 0.5_DOUBLE*(pr_bar - pl_bar)* &
      (/0.0_DOUBLE, 1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, u_et/)

    lr_flux(:, 2) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(u_et)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et)) &
      + 0.5_DOUBLE*(pr_bar - pl_bar)* &
      (/0.0_DOUBLE, 1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, u_et/)

    lr_flux(:, 2) = -lr_flux(:, 2)
  end subroutine multi_point

  subroutine compute_lambdas_and_solve_nodal_velocity(mesh, id_vert, sol_w_lr, lambda, u_bars, u_node, p_bound)
    use ns_global_data_module, only: bc_style, boundary_2d
    use linear_solver_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, 2, mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(in) :: sol_w_lr
    real(kind=DOUBLE), dimension(2, mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(inout) :: lambda
    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(inout) :: u_bars
    real(kind=DOUBLE), dimension(3), intent(inout) :: u_node
    real(kind=DOUBLE), intent(inout) :: p_bound

    integer(kind=ENTIER) :: iter
    integer(kind=ENTIER) :: j, id_sub_face, id_face, le, re
    real(kind=DOUBLE) :: rhol, ul, pl, al, lambda_l
    real(kind=DOUBLE) :: rhor, ur, pr, ar, lambda_r
    real(kind=DOUBLE) :: u_et
    real(kind=DOUBLE), dimension(3) :: Rp, Bp
    real(kind=DOUBLE), dimension(3, 3) :: mat

    Bp = wall_normal(mesh, id_vert)
    u_node = 0.0_DOUBLE
    iter = 0
    do while (iter < 4)
      iter = iter + 1

      mat(:, :) = 0.0_DOUBLE
      Rp(:) = 0.0_DOUBLE

      do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
        id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face

        le = mesh%face(id_face)%left_neigh
        re = mesh%face(id_face)%right_neigh

        rhol = sol_w_lr(1, 1, j)
        ul = sol_w_lr(2, 1, j)
        pl = sol_w_lr(5, 1, j)
        al = sound_speed_w(sol_w_lr(:, 1, j))

        rhor = sol_w_lr(1, 2, j)
        ur = sol_w_lr(2, 2, j)
        pr = sol_w_lr(5, 2, j)
        ar = sound_speed_w(sol_w_lr(:, 2, j))

        lambda_l = lambda(1, j)
        lambda_r = lambda(2, j)

        if (iter == 1) then
          lambda_l = max(lambda_l, al*rhol, &
            sqrt(rhol*max(0.0_DOUBLE, pr - pl)), &
            -rhol*(ur - ul))
          lambda_r = max(lambda_r, ar*rhor, &
            sqrt(rhor*max(0.0_DOUBLE, pl - pr)), &
            -rhor*(ur - ul))

          u_bars(j) = (lambda_l*ul + lambda_r*ur - (pr - pl)) &
            /(lambda_r + lambda_l)

          lambda_l = max(lambda_l, al*rhol*(1.0_DOUBLE &
            + 1.5_DOUBLE*max(0.0_DOUBLE, -(u_bars(j) - ul)/al)))
          lambda_r = max(lambda_r, ar*rhor*(1.0_DOUBLE &
            + 1.5_DOUBLE*max(0.0_DOUBLE, (u_bars(j) - ur)/ar)))
        else
          if (is_wall(re) .and. bc_style == 1) then
            lambda_l = max(lambda_l, al*rhol*(1.0_DOUBLE &
              + 1.5_DOUBLE*max(0.0_DOUBLE, -(u_bars(j) - ul)/al)))
            lambda_r = max(lambda_r, ar*rhor*(1.0_DOUBLE &
              + 1.5_DOUBLE*max(0.0_DOUBLE, (u_bars(j) - ur)/ar)))
          else
            u_et = dot_product(u_node, mesh%sub_face(id_sub_face)%norm)
            lambda_l = max(lambda_l, al*rhol*(1.0_DOUBLE &
              + 1.5_DOUBLE*max(0.0_DOUBLE, -(u_et - ul)/al)))
            lambda_r = max(lambda_r, ar*rhor*(1.0_DOUBLE &
              + 1.5_DOUBLE*max(0.0_DOUBLE, (u_et - ur)/ar)))
          end if
        end if

        lambda(1, j) = lambda_l
        lambda(2, j) = lambda_r

        u_bars(j) = (lambda_l*ul + lambda_r*ur - (pr - pl)) &
          /(lambda_r + lambda_l)

        if (is_wall(re)) then
          if (bc_style == 0 .or. bc_style == 2) then
            mat = mat &
              + (lambda_r + lambda_l)*mesh%sub_face(id_sub_face)%area &
              *tensor_product_3(mesh%sub_face(id_sub_face)%norm, &
              mesh%sub_face(id_sub_face)%norm)
            if (bc_style == 0) then
              Rp = Rp + mesh%sub_face(id_sub_face)%area &
                *(lambda_r + lambda_l)*u_bars(j) &
                *mesh%sub_face(id_sub_face)%norm
            else if (bc_style == 2) then
              if (boundary_2d) then
                if (abs(mesh%sub_face(id_sub_face)%norm(3)) < 1e-8_DOUBLE) then
                  Rp = Rp + mesh%sub_face(id_sub_face)%area &
                    *pl*mesh%sub_face(id_sub_face)%norm
                end if
              else
                Rp = Rp + mesh%sub_face(id_sub_face)%area &
                  *pl*mesh%sub_face(id_sub_face)%norm
              end if
            end if
          end if
        else
          mat = mat + (lambda_r + lambda_l)*mesh%sub_face(id_sub_face)%area &
            *tensor_product_3(mesh%sub_face(id_sub_face)%norm, &
            mesh%sub_face(id_sub_face)%norm)
          Rp = Rp + mesh%sub_face(id_sub_face)%area* &
            (lambda_r + lambda_l)*u_bars(j)*mesh%sub_face(id_sub_face)%norm
        end if
      end do

      if (bc_style == 2) then
        if (maxval(abs(Bp)) > 1e-8_DOUBLE) then
          call pseudo_inverse_inplace_lapack(3, mat)
          p_bound = dot_product(matmul(mat, Rp), Bp)/dot_product(matmul(mat, Bp), Bp)
          u_node = matmul(mat, Rp - p_bound*Bp)
          u_node = u_node - dot_product(u_node, Bp)/dot_product(Bp, Bp)*Bp
        else
          call lu_solve_inplace_lapack(3, mat, u_node, Rp)
        end if
      else
        call pseudo_inverse_inplace_lapack(3, mat)
        u_node = matmul(mat, Rp)
      end if
    end do
  end subroutine compute_lambdas_and_solve_nodal_velocity

  subroutine compute_lambdas_and_solve_nodal_velocity_iso(mesh, id_vert, sol_w_lr, lambda, u_bars, u_node, p_bound)
    use ns_global_data_module, only: bc_style, boundary_2d
    use linear_solver_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, 2, mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(in) :: sol_w_lr
    real(kind=DOUBLE), dimension(2, mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(inout) :: lambda
    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(inout) :: u_bars
    real(kind=DOUBLE), dimension(3), intent(inout) :: u_node
    real(kind=DOUBLE), intent(inout) :: p_bound

    integer(kind=ENTIER) :: iter
    integer(kind=ENTIER) :: j, id_sub_face, id_face, le, re
    real(kind=DOUBLE) :: rhol, ul, pl, al, lambda_l
    real(kind=DOUBLE) :: rhor, ur, pr, ar, lambda_r
    real(kind=DOUBLE) :: u_et
    real(kind=DOUBLE), dimension(3) :: Rp, Bp
    real(kind=DOUBLE), dimension(3, 3) :: mat

    Bp = wall_normal(mesh, id_vert)
    u_node = 0.0_DOUBLE
    iter = 0
    do while (iter < 4)
      iter = iter + 1

      mat(:, :) = 0.0_DOUBLE
      Rp(:) = 0.0_DOUBLE

      do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
        id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face

        le = mesh%face(id_face)%left_neigh
        re = mesh%face(id_face)%right_neigh

        rhol = sol_w_lr(1, 1, j)
        ul = sol_w_lr(2, 1, j)
        pl = sol_w_lr(5, 1, j)
        al = sound_speed_w(sol_w_lr(:, 1, j))

        rhor = sol_w_lr(1, 2, j)
        ur = sol_w_lr(2, 2, j)
        pr = sol_w_lr(5, 2, j)
        ar = sound_speed_w(sol_w_lr(:, 2, j))

        lambda_l = lambda(1, j)
        lambda_r = lambda(2, j)

        if (iter == 1) then
          lambda_l = max(lambda_l, al*rhol, &
            sqrt(rhol*max(0.0_DOUBLE, pr - pl)), &
            -rhol*(ur - ul))
          lambda_r = max(lambda_r, ar*rhor, &
            sqrt(rhor*max(0.0_DOUBLE, pl - pr)), &
            -rhor*(ur - ul))

          u_bars(j) = (lambda_l*ul + lambda_r*ur - (pr - pl)) &
            /(lambda_r + lambda_l)

          lambda_l = max(lambda_l, al*rhol*(1.0_DOUBLE &
            + 1.5_DOUBLE*max(0.0_DOUBLE, -(u_bars(j) - ul)/al)))
          lambda_r = max(lambda_r, ar*rhor*(1.0_DOUBLE &
            + 1.5_DOUBLE*max(0.0_DOUBLE, (u_bars(j) - ur)/ar)))
        else
          if (is_wall(re) .and. bc_style == 1) then
            lambda_l = max(lambda_l, al*rhol*(1.0_DOUBLE &
              + 1.5_DOUBLE*max(0.0_DOUBLE, -(u_bars(j) - ul)/al)))
            lambda_r = max(lambda_r, ar*rhor*(1.0_DOUBLE &
              + 1.5_DOUBLE*max(0.0_DOUBLE, (u_bars(j) - ur)/ar)))
          else
            u_et = dot_product(u_node, mesh%sub_face(id_sub_face)%norm)
            lambda_l = max(lambda_l, al*rhol*(1.0_DOUBLE &
              + 1.5_DOUBLE*max(0.0_DOUBLE, -(u_et - ul)/al)))
            lambda_r = max(lambda_r, ar*rhor*(1.0_DOUBLE &
              + 1.5_DOUBLE*max(0.0_DOUBLE, (u_et - ur)/ar)))
          end if
        end if

        lambda(1, j) = lambda_l
        lambda(2, j) = lambda_r

        u_bars(j) = (lambda_l*ul + lambda_r*ur - (pr - pl)) &
          /(lambda_r + lambda_l)

        if (is_wall(re)) then
          if (bc_style == 0 .or. bc_style == 2) then
            mat = mat &
              + (lambda_r + lambda_l) &
              *tensor_product_3(mesh%sub_face(id_sub_face)%norm, &
              mesh%sub_face(id_sub_face)%norm)
            if (bc_style == 0) then
              Rp = Rp + (lambda_r + lambda_l)*u_bars(j) &
                *mesh%sub_face(id_sub_face)%norm
            else if (bc_style == 2) then
              if (boundary_2d) then
                if (abs(mesh%sub_face(id_sub_face)%norm(3)) < 1e-8_DOUBLE) then
                  Rp = Rp + pl*mesh%sub_face(id_sub_face)%norm
                end if
              else
                Rp = Rp + pl*mesh%sub_face(id_sub_face)%norm
              end if
            end if
          end if
        else
          mat = mat + (lambda_r + lambda_l) &
            *tensor_product_3(mesh%sub_face(id_sub_face)%norm, mesh%sub_face(id_sub_face)%norm)
          Rp = Rp + (lambda_r + lambda_l)*u_bars(j)*mesh%sub_face(id_sub_face)%norm
        end if
      end do

      if (bc_style == 2) then
        if (maxval(abs(Bp)) > 1e-8_DOUBLE) then
          call pseudo_inverse_inplace_lapack(3, mat)
          p_bound = dot_product(matmul(mat, Rp), Bp)/dot_product(matmul(mat, Bp), Bp)
          u_node = matmul(mat, Rp - p_bound*Bp)
          u_node = u_node - dot_product(u_node, Bp)/dot_product(Bp, Bp)*Bp
        else
          call lu_solve_inplace_lapack(3, mat, u_node, Rp)
        end if
      else
        call pseudo_inverse_inplace_lapack(3, mat)
        u_node = matmul(mat, Rp)
      end if
    end do
  end subroutine compute_lambdas_and_solve_nodal_velocity_iso

  subroutine compute_lambdas_and_solve_nodal_velocity_new(ng, weight, norm, &
      lambda_l, lambda_r, sol_w_l, sol_w_r, up)
    use ns_global_data_module, only: bc_style, boundary_2d
    use linear_solver_module
    implicit none

    integer(kind=ENTIER), intent(in) :: ng
    real(kind=DOUBLE), dimension(ng), intent(in) :: weight
    real(kind=DOUBLE), dimension(ng), intent(inout) :: lambda_l, lambda_r
    real(kind=DOUBLE), dimension(5, ng), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(3, ng), intent(in) :: norm
    real(kind=DOUBLE), dimension(3), intent(inout) :: up

    integer(kind=ENTIER) :: iter, j
    real(kind=DOUBLE) :: rhol, ul, pl, al
    real(kind=DOUBLE) :: rhor, ur, pr, ar
    real(kind=DOUBLE) :: u_bar, u_star
    real(kind=DOUBLE), dimension(3) :: rhs
    real(kind=DOUBLE), dimension(3, 3) :: mat

    !Init lambdas with two-point positivity conditions 
    !and Duckowitz correction
    mat(:, :) = 0.0_DOUBLE
    rhs(:) = 0.0_DOUBLE
    do j=1, ng
      rhol = sol_w_l(1, j)
      ul = sol_w_l(2, j)
      pl = sol_w_l(5, j)
      al = sound_speed_w(sol_w_l(:, j))

      rhor = sol_w_r(1, j)
      ur = sol_w_r(2, j)
      pr = sol_w_r(5, j)
      ar = sound_speed_w(sol_w_r(:, j))

      lambda_l(j) = max(al*rhol, &
        sqrt(rhol*max(0.0_DOUBLE, pr - pl)), &
        -rhol*(ur - ul))
      lambda_r(j) = max(ar*rhor, &
        sqrt(rhor*max(0.0_DOUBLE, pl - pr)), &
        -rhor*(ur - ul))
      u_bar = (lambda_l(j)*ul + lambda_r(j)*ur - (pr - pl))&
        /(lambda_r(j) + lambda_l(j))

      !Duckowitz correction
      lambda_l(j) = max(lambda_l(j), al*rhol*(1.0_DOUBLE &
        + 1.5_DOUBLE*max(0.0_DOUBLE, -(u_bar - ul)/al)))
      lambda_r(j) = max(lambda_r(j), ar*rhor*(1.0_DOUBLE &
        + 1.5_DOUBLE*max(0.0_DOUBLE, (u_bar - ur)/ar)))
      u_bar = (lambda_l(j)*ul + lambda_r(j)*ur - (pr - pl))&
        /(lambda_r(j) + lambda_l(j))

      mat = mat + weight(j)*(lambda_r(j)+lambda_l(j))&
        *tensor_product_3(norm(:, j), norm(:, j))
      rhs = rhs + weight(j)*(lambda_r(j)+lambda_l(j))&
        *u_bar*norm(:, j)
    end do
    call pseudo_inverse_inplace_lapack(3, mat)
    up = matmul(mat, rhs)

    iter = 0
    do while (iter < 2)
      iter = iter + 1

      mat(:, :) = 0.0_DOUBLE
      rhs(:) = 0.0_DOUBLE
      do j = 1, ng
        rhol = sol_w_l(1, j)
        ul = sol_w_l(2, j)
        pl = sol_w_l(5, j)
        al = sound_speed_w(sol_w_l(:, j))

        rhor = sol_w_r(1, j)
        ur = sol_w_r(2, j)
        pr = sol_w_r(5, j)
        ar = sound_speed_w(sol_w_r(:, j))

        u_star = dot_product(up, norm(:, j))

        lambda_l(j) = max(lambda_l(j), al*rhol*(1.0_DOUBLE &
          + 1.5_DOUBLE*max(0.0_DOUBLE, -(u_star - ul)/al)))
        lambda_r(j) = max(lambda_r(j), ar*rhor*(1.0_DOUBLE &
          + 1.5_DOUBLE*max(0.0_DOUBLE, (u_star - ur)/ar)))
        u_bar = (lambda_l(j)*ul + lambda_r(j)*ur - (pr - pl)) &
          /(lambda_r(j) + lambda_l(j))

        mat = mat + weight(j)*(lambda_r(j)+lambda_l(j))&
          *tensor_product_3(norm(:, j), norm(:, j))
        rhs = rhs + weight(j)*(lambda_r(j)+lambda_l(j))&
          *u_bar*norm(:, j)
      end do
      call pseudo_inverse_inplace_lapack(3, mat)
      up = matmul(mat, rhs)

    end do
  end subroutine compute_lambdas_and_solve_nodal_velocity_new
end module ns_euler_rs_module