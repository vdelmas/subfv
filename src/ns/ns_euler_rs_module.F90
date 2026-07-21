module ns_euler_rs_module
  use precision_module
  use mesh_module
  use ns_euler_primitives_module
  use ns_euler_recon_module
  implicit none

contains
  subroutine three_wave(sol_w_l, sol_w_r, n, lr_flux, sl, sr)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), dimension(5) :: sol_l, sol_r
    real(kind=DOUBLE), dimension(5, 2), intent(inout) :: lr_flux
    real(kind=DOUBLE), intent(inout) :: sl, sr

    real(kind=DOUBLE) :: rhol, rhor, vn_l, vn_r
    real(kind=DOUBLE) :: pl, pr, rhol_et, rhor_et
    real(kind=DOUBLE) :: v_et, el, er, al, ar
    real(kind=DOUBLE) :: v_bar, pl_bar, pr_bar
    real(kind=DOUBLE) :: lambda_l, lambda_r
    real(kind=DOUBLE), dimension(5) :: fl, fr, sol_l_et, sol_r_et

    rhol = sol_w_l(1)
    vn_l = dot_product(sol_w_l(2:4), n)
    pl = sol_w_l(5)
    sol_l = primit_to_conserv(sol_w_l)
    el = sol_l(5)/rhol
    al = sound_speed_w(sol_w_l)

    rhor = sol_w_r(1)
    vn_r = dot_product(sol_w_r(2:4), n)
    pr = sol_w_r(5)
    sol_r = primit_to_conserv(sol_w_r)
    er = sol_r(5)/rhor
    ar = sound_speed_w(sol_w_r)

    fl(1)   = vn_l*sol_l(1)
    fl(2:4) = vn_l*sol_l(2:4) + pl*n
    fl(5)   = (sol_l(5) + pl)*vn_l

    fr(1)   = vn_r*sol_r(1)
    fr(2:4) = vn_r*sol_r(2:4) + pr*n
    fr(5)   = (sol_r(5) + pr)*vn_r

    lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(vn_r - vn_l))
    lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(vn_r - vn_l))
    v_bar = (lambda_l*vn_l + lambda_r*vn_r - (pr - pl))/(lambda_r + lambda_l)
    v_et = v_bar

    rhol_et = 1.0_DOUBLE/(1.0_DOUBLE/rhol + (v_et - vn_l)/lambda_l)
    pl_bar = pl - lambda_l*(v_et - vn_l)

    sol_l_et(1)   = rhol_et
    sol_l_et(2:4) = rhol_et*(sol_w_l(2:4) + (v_et - vn_l)*n)
    sol_l_et(5)   = rhol_et*(el + (pl*vn_l - pl_bar*v_et)/lambda_l)

    rhor_et = 1.0_DOUBLE/(1.0_DOUBLE/rhor + (vn_r - v_et)/lambda_r)
    pr_bar = pr + lambda_r*(v_et - vn_r)

    sol_r_et(1)   = rhor_et
    sol_r_et(2:4) = rhor_et*(sol_w_r(2:4) + (v_et - vn_r)*n)
    sol_r_et(5)   = rhor_et*(er + (pr_bar*v_et - pr*vn_r)/lambda_r)

    sl = vn_l - lambda_l/rhol
    sr = vn_r + lambda_r/rhor

    lr_flux(:, 1) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(v_et)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et))

    lr_flux(:, 2) = -lr_flux(:, 1)
  end subroutine three_wave

  subroutine modified_three_wave(sol_w_l, sol_w_r, n, lr_flux, sl, sr)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), dimension(5) :: sol_l, sol_r
    real(kind=DOUBLE), dimension(5, 2), intent(inout) :: lr_flux
    real(kind=DOUBLE), intent(inout) :: sl, sr

    real(kind=DOUBLE) :: rhol, rhor, vn_l, vn_r
    real(kind=DOUBLE) :: pl, pr, rhol_et, rhor_et
    real(kind=DOUBLE) :: v_et, el, er, al, ar
    real(kind=DOUBLE) :: v_bar, pl_bar, pr_bar
    real(kind=DOUBLE) :: lambda_l, lambda_r, tang_coeff, tang_energy
    real(kind=DOUBLE), dimension(3) :: vt_l, vt_r, vt_et, dvt
    real(kind=DOUBLE), dimension(5) :: fl, fr, sol_l_et, sol_r_et

    rhol = sol_w_l(1)
    vn_l = dot_product(sol_w_l(2:4), n)
    pl = sol_w_l(5)
    sol_l = primit_to_conserv(sol_w_l)
    el = sol_l(5)/rhol
    al = sound_speed_w(sol_w_l)

    rhor = sol_w_r(1)
    vn_r = dot_product(sol_w_r(2:4), n)
    pr = sol_w_r(5)
    sol_r = primit_to_conserv(sol_w_r)
    er = sol_r(5)/rhor
    ar = sound_speed_w(sol_w_r)

    fl(1)   = vn_l*sol_l(1)
    fl(2:4) = vn_l*sol_l(2:4) + pl*n
    fl(5)   = (sol_l(5) + pl)*vn_l

    fr(1)   = vn_r*sol_r(1)
    fr(2:4) = vn_r*sol_r(2:4) + pr*n
    fr(5)   = (sol_r(5) + pr)*vn_r

    lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(vn_r - vn_l))
    lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(vn_r - vn_l))
    v_bar = (lambda_l*vn_l + lambda_r*vn_r - (pr - pl))/(lambda_r + lambda_l)
    v_et = v_bar

    vt_l = sol_w_l(2:4) - vn_l*n
    vt_r = sol_w_r(2:4) - vn_r*n
    vt_et = (lambda_l*vt_l + lambda_r*vt_r)/(lambda_l + lambda_r)
    dvt = vt_r - vt_l

    tang_coeff = -lambda_l*lambda_r/(lambda_l + lambda_r)
    tang_energy = tang_coeff*dot_product(dvt, vt_et)

    rhol_et = 1.0_DOUBLE/(1.0_DOUBLE/rhol + (v_et - vn_l)/lambda_l)
    pl_bar = pl - lambda_l*(v_et - vn_l)

    sol_l_et(1)   = rhol_et
    sol_l_et(2:4) = rhol_et*(v_et*n + vt_et)
    sol_l_et(5)   = rhol_et*(el + (pl*vn_l - pl_bar*v_et - tang_energy)/lambda_l)

    rhor_et = 1.0_DOUBLE/(1.0_DOUBLE/rhor + (vn_r - v_et)/lambda_r)
    pr_bar = pr + lambda_r*(v_et - vn_r)

    sol_r_et(1)   = rhor_et
    sol_r_et(2:4) = rhor_et*(v_et*n + vt_et)
    sol_r_et(5)   = rhor_et*(er + (pr_bar*v_et - pr*vn_r + tang_energy)/lambda_r)

    sl = vn_l - lambda_l/rhol
    sr = vn_r + lambda_r/rhor

    lr_flux(:, 1) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(v_et)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et))

    lr_flux(:, 2) = -lr_flux(:, 1)
  end subroutine modified_three_wave

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

  subroutine multi_point(sol_w_l, sol_w_r, n, lr_flux, vn_nodal, &
      lambda_l, lambda_r, sl, sr)
    implicit none

    real(kind=DOUBLE), intent(in) :: vn_nodal
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), dimension(5) :: sol_l, sol_r
    real(kind=DOUBLE), intent(inout) :: lambda_l, lambda_r
    real(kind=DOUBLE), dimension(5, 2), intent(inout) :: lr_flux
    real(kind=DOUBLE), intent(inout) :: sl, sr

    real(kind=DOUBLE) :: rhol, rhor, vn_l, vn_r
    real(kind=DOUBLE) :: pl, pr, rhol_et, rhor_et
    real(kind=DOUBLE) :: v_et, el, er
    real(kind=DOUBLE) :: v_bar, pl_bar, pr_bar
    real(kind=DOUBLE), dimension(5) :: fl, fr
    real(kind=DOUBLE), dimension(5) :: sol_l_et, sol_r_et

    rhol = sol_w_l(1)
    vn_l = dot_product(sol_w_l(2:4), n)
    pl = sol_w_l(5)
    sol_l = primit_to_conserv(sol_w_l)
    el = sol_l(5)/rhol

    rhor = sol_w_r(1)
    vn_r = dot_product(sol_w_r(2:4), n)
    pr = sol_w_r(5)
    sol_r = primit_to_conserv(sol_w_r)
    er = sol_r(5)/rhor

    fl(1)   = vn_l*sol_l(1)
    fl(2:4) = vn_l*sol_l(2:4) + pl*n
    fl(5)   = (sol_l(5) + pl)*vn_l

    fr(1)   = vn_r*sol_r(1)
    fr(2:4) = vn_r*sol_r(2:4) + pr*n
    fr(5)   = (sol_r(5) + pr)*vn_r

    v_bar = (lambda_l*vn_l + lambda_r*vn_r - (pr - pl))/(lambda_r + lambda_l)

    v_et = vn_nodal

    rhol_et = 1.0_DOUBLE/(1.0_DOUBLE/rhol + (v_et - vn_l)/lambda_l)
    pl_bar = pl - lambda_l*(v_et - vn_l)

    sol_l_et(1)   = rhol_et
    sol_l_et(2:4) = rhol_et*(sol_w_l(2:4) + (v_et - vn_l)*n)
    sol_l_et(5)   = rhol_et*(el + (pl*vn_l - pl_bar*v_et)/lambda_l)

    rhor_et = 1.0_DOUBLE/(1.0_DOUBLE/rhor + (vn_r - v_et)/lambda_r)
    pr_bar = pr + lambda_r*(v_et - vn_r)

    sol_r_et(1)   = rhor_et
    sol_r_et(2:4) = rhor_et*(sol_w_r(2:4) + (v_et - vn_r)*n)
    sol_r_et(5)   = rhor_et*(er + (pr_bar*v_et - pr*vn_r)/lambda_r)

    if (rhol_et < 0.0_DOUBLE .or. rhor_et < 0.0_DOUBLE) then
      print *, "Negative specific volume MPCC !", rhol_et, rhor_et
      print*, rhol, rhor
      print*, lambda_l, lambda_r
      print*, sol_l_et(2:4)/rhol_et
      print*, sol_r_et(2:4)/rhor_et
      print*,"R"
      print*, 1.0_DOUBLE/rhor, (vn_r - v_et)/lambda_r
      error stop
    end if

    sl = vn_l - lambda_l/rhol
    sr = vn_r + lambda_r/rhor

    lr_flux(:, 1) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(v_et)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et)) &
      - 0.5_DOUBLE*(pr_bar - pl_bar)* &
      (/0.0_DOUBLE, n, v_et/)

    lr_flux(:, 2) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(v_et)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et)) &
      + 0.5_DOUBLE*(pr_bar - pl_bar)* &
      (/0.0_DOUBLE, n, v_et/)

    lr_flux(:, 2) = -lr_flux(:, 2)
  end subroutine multi_point

  subroutine compute_lambdas_and_solve_nodal_velocity(mesh, id_vert, sol_w_lr, lambda, v_bars, v_node, p_bound)
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
      intent(inout) :: v_bars
    real(kind=DOUBLE), dimension(3), intent(inout) :: v_node
    real(kind=DOUBLE), intent(inout) :: p_bound

    integer(kind=ENTIER) :: iter
    integer(kind=ENTIER) :: j, id_sub_face, id_face, le, re
    real(kind=DOUBLE) :: rhol, vn_l, pl, al, lambda_l
    real(kind=DOUBLE) :: rhor, vn_r, pr, ar, lambda_r
    real(kind=DOUBLE) :: v_et
    real(kind=DOUBLE), dimension(3) :: Rp, Bp
    real(kind=DOUBLE), dimension(3, 3) :: mat

    Bp = wall_normal(mesh, id_vert)
    v_node = 0.0_DOUBLE
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
        vn_l = dot_product(sol_w_lr(2:4, 1, j), mesh%sub_face(id_sub_face)%norm)
        pl = sol_w_lr(5, 1, j)
        al = sound_speed_w(sol_w_lr(:, 1, j))

        rhor = sol_w_lr(1, 2, j)
        vn_r = dot_product(sol_w_lr(2:4, 2, j), mesh%sub_face(id_sub_face)%norm)
        pr = sol_w_lr(5, 2, j)
        ar = sound_speed_w(sol_w_lr(:, 2, j))

        lambda_l = lambda(1, j)
        lambda_r = lambda(2, j)

        if (iter == 1) then
          lambda_l = max(lambda_l, al*rhol, &
            sqrt(rhol*max(0.0_DOUBLE, pr - pl)), &
            -rhol*(vn_r - vn_l))
          lambda_r = max(lambda_r, ar*rhor, &
            sqrt(rhor*max(0.0_DOUBLE, pl - pr)), &
            -rhor*(vn_r - vn_l))

          v_bars(j) = (lambda_l*vn_l + lambda_r*vn_r - (pr - pl)) &
            /(lambda_r + lambda_l)

          lambda_l = max(lambda_l, al*rhol*(1.0_DOUBLE &
            + 1.5_DOUBLE*max(0.0_DOUBLE, -(v_bars(j) - vn_l)/al)))
          lambda_r = max(lambda_r, ar*rhor*(1.0_DOUBLE &
            + 1.5_DOUBLE*max(0.0_DOUBLE, (v_bars(j) - vn_r)/ar)))
        else
          if (is_wall(re) .and. bc_style == 1) then
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

        v_bars(j) = (lambda_l*vn_l + lambda_r*vn_r - (pr - pl)) &
          /(lambda_r + lambda_l)

        if (is_wall(re)) then
          if (bc_style == 0 .or. bc_style == 2) then
            mat = mat &
              + (lambda_r + lambda_l)*mesh%sub_face(id_sub_face)%area &
              *tensor_product_3(mesh%sub_face(id_sub_face)%norm, &
              mesh%sub_face(id_sub_face)%norm)
            if (bc_style == 0) then
              Rp = Rp + mesh%sub_face(id_sub_face)%area &
                *(lambda_r + lambda_l)*v_bars(j) &
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
            (lambda_r + lambda_l)*v_bars(j)*mesh%sub_face(id_sub_face)%norm
        end if
      end do

      if (bc_style == 2) then
        if (maxval(abs(Bp)) > 1e-8_DOUBLE) then
          call pseudo_inverse_inplace_lapack(3, mat)
          p_bound = dot_product(matmul(mat, Rp), Bp)/dot_product(matmul(mat, Bp), Bp)
          v_node = matmul(mat, Rp - p_bound*Bp)
          v_node = v_node - dot_product(v_node, Bp)/dot_product(Bp, Bp)*Bp
        else
          call lu_solve_inplace_lapack(3, mat, v_node, Rp)
        end if
      else
        call pseudo_inverse_inplace_lapack(3, mat)
        v_node = matmul(mat, Rp)
      end if
    end do
  end subroutine compute_lambdas_and_solve_nodal_velocity

  subroutine compute_lambdas_and_solve_nodal_velocity_iso(mesh, id_vert, sol_w_lr, lambda, v_bars, v_node, p_bound)
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
      intent(inout) :: v_bars
    real(kind=DOUBLE), dimension(3), intent(inout) :: v_node
    real(kind=DOUBLE), intent(inout) :: p_bound

    integer(kind=ENTIER) :: iter
    integer(kind=ENTIER) :: j, id_sub_face, id_face, le, re
    real(kind=DOUBLE) :: rhol, vn_l, pl, al, lambda_l
    real(kind=DOUBLE) :: rhor, vn_r, pr, ar, lambda_r
    real(kind=DOUBLE) :: v_et
    real(kind=DOUBLE), dimension(3) :: Rp, Bp
    real(kind=DOUBLE), dimension(3, 3) :: mat

    Bp = wall_normal(mesh, id_vert)
    v_node = 0.0_DOUBLE
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
        vn_l = dot_product(sol_w_lr(2:4, 1, j), mesh%sub_face(id_sub_face)%norm)
        pl = sol_w_lr(5, 1, j)
        al = sound_speed_w(sol_w_lr(:, 1, j))

        rhor = sol_w_lr(1, 2, j)
        vn_r = dot_product(sol_w_lr(2:4, 2, j), mesh%sub_face(id_sub_face)%norm)
        pr = sol_w_lr(5, 2, j)
        ar = sound_speed_w(sol_w_lr(:, 2, j))

        lambda_l = lambda(1, j)
        lambda_r = lambda(2, j)

        if (iter == 1) then
          lambda_l = max(lambda_l, al*rhol, &
            sqrt(rhol*max(0.0_DOUBLE, pr - pl)), &
            -rhol*(vn_r - vn_l))
          lambda_r = max(lambda_r, ar*rhor, &
            sqrt(rhor*max(0.0_DOUBLE, pl - pr)), &
            -rhor*(vn_r - vn_l))

          v_bars(j) = (lambda_l*vn_l + lambda_r*vn_r - (pr - pl)) &
            /(lambda_r + lambda_l)

          lambda_l = max(lambda_l, al*rhol*(1.0_DOUBLE &
            + 1.5_DOUBLE*max(0.0_DOUBLE, -(v_bars(j) - vn_l)/al)))
          lambda_r = max(lambda_r, ar*rhor*(1.0_DOUBLE &
            + 1.5_DOUBLE*max(0.0_DOUBLE, (v_bars(j) - vn_r)/ar)))
        else
          if (is_wall(re) .and. bc_style == 1) then
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

        v_bars(j) = (lambda_l*vn_l + lambda_r*vn_r - (pr - pl)) &
          /(lambda_r + lambda_l)

        if (is_wall(re)) then
          if (bc_style == 0 .or. bc_style == 2) then
            mat = mat &
              + (lambda_r + lambda_l) &
              *tensor_product_3(mesh%sub_face(id_sub_face)%norm, &
              mesh%sub_face(id_sub_face)%norm)
            if (bc_style == 0) then
              Rp = Rp + (lambda_r + lambda_l)*v_bars(j) &
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
          Rp = Rp + (lambda_r + lambda_l)*v_bars(j)*mesh%sub_face(id_sub_face)%norm
        end if
      end do

      if (bc_style == 2) then
        if (maxval(abs(Bp)) > 1e-8_DOUBLE) then
          call pseudo_inverse_inplace_lapack(3, mat)
          p_bound = dot_product(matmul(mat, Rp), Bp)/dot_product(matmul(mat, Bp), Bp)
          v_node = matmul(mat, Rp - p_bound*Bp)
          v_node = v_node - dot_product(v_node, Bp)/dot_product(Bp, Bp)*Bp
        else
          call lu_solve_inplace_lapack(3, mat, v_node, Rp)
        end if
      else
        call pseudo_inverse_inplace_lapack(3, mat)
        v_node = matmul(mat, Rp)
      end if
    end do
  end subroutine compute_lambdas_and_solve_nodal_velocity_iso

  subroutine multi_point_pressure(sol_w_l, sol_w_r, n, lr_flux, p_nodal, &
      lambda_l, lambda_r, sl, sr)
    implicit none

    real(kind=DOUBLE), intent(in) :: p_nodal
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), dimension(5) :: sol_l, sol_r
    real(kind=DOUBLE), intent(inout) :: lambda_l, lambda_r
    real(kind=DOUBLE), dimension(5, 2), intent(inout) :: lr_flux
    real(kind=DOUBLE), intent(inout) :: sl, sr

    real(kind=DOUBLE) :: rhol, rhor, vn_l, vn_r, pl, pr, el, er
    real(kind=DOUBLE) :: vn_bar, vn_l_et, vn_r_et, rhol_et, rhor_et
    real(kind=DOUBLE), dimension(3) :: vt_l, vt_r
    real(kind=DOUBLE), dimension(5) :: fl, fr, sol_l_et, sol_r_et

    rhol  = sol_w_l(1)
    vn_l  = dot_product(sol_w_l(2:4), n)
    pl    = sol_w_l(5)
    vt_l  = sol_w_l(2:4) - vn_l*n
    sol_l = primit_to_conserv(sol_w_l)
    el    = sol_l(5)/rhol

    rhor  = sol_w_r(1)
    vn_r  = dot_product(sol_w_r(2:4), n)
    pr    = sol_w_r(5)
    vt_r  = sol_w_r(2:4) - vn_r*n
    sol_r = primit_to_conserv(sol_w_r)
    er    = sol_r(5)/rhor

    fl(1)   = vn_l*sol_l(1)
    fl(2:4) = vn_l*sol_l(2:4) + pl*n
    fl(5)   = (sol_l(5) + pl)*vn_l

    fr(1)   = vn_r*sol_r(1)
    fr(2:4) = vn_r*sol_r(2:4) + pr*n
    fr(5)   = (sol_r(5) + pr)*vn_r

    vn_bar = (lambda_l*vn_l + lambda_r*vn_r - (pr - pl))/(lambda_r + lambda_l)

    vn_l_et  = vn_l - (p_nodal - pl)/lambda_l
    rhol_et  = 1.0_DOUBLE/(1.0_DOUBLE/rhol + (vn_bar - vn_l)/lambda_l)
    sol_l_et(1)   = rhol_et
    sol_l_et(2:4) = rhol_et*(vt_l + vn_l_et*n)
    sol_l_et(5)   = rhol_et*(el + (pl*vn_l - p_nodal*vn_bar)/lambda_l)

    vn_r_et  = vn_r + (p_nodal - pr)/lambda_r
    rhor_et  = 1.0_DOUBLE/(1.0_DOUBLE/rhor + (vn_r - vn_bar)/lambda_r)
    sol_r_et(1)   = rhor_et
    sol_r_et(2:4) = rhor_et*(vt_r + vn_r_et*n)
    sol_r_et(5)   = rhor_et*(er + (p_nodal*vn_bar - pr*vn_r)/lambda_r)

    sl = vn_l - lambda_l/rhol
    sr = vn_r + lambda_r/rhor

    lr_flux(:, 1) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(vn_bar)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et))
    lr_flux(:, 2) = -lr_flux(:, 1)
  end subroutine multi_point_pressure

  subroutine multi_point_pressure_ph(sol_w_l, sol_w_r, n, lr_flux, p_nodal, &
      lambda_l, lambda_r, sl, sr)
    implicit none

    real(kind=DOUBLE), intent(in) :: p_nodal
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), dimension(5) :: sol_l, sol_r
    real(kind=DOUBLE), intent(inout) :: lambda_l, lambda_r
    real(kind=DOUBLE), dimension(5, 2), intent(inout) :: lr_flux
    real(kind=DOUBLE), intent(inout) :: sl, sr

    real(kind=DOUBLE) :: rhol, rhor, vn_l, vn_r, pl, pr, el, er
    real(kind=DOUBLE) :: vn_bar, vn_l_et, vn_r_et, rhol_et, rhor_et
    real(kind=DOUBLE), dimension(3) :: vt_l, vt_r
    real(kind=DOUBLE), dimension(5) :: fl, fr, sol_l_et, sol_r_et

    rhol  = sol_w_l(1)
    vn_l  = dot_product(sol_w_l(2:4), n)
    pl    = sol_w_l(5)
    vt_l  = sol_w_l(2:4) - vn_l*n
    sol_l = primit_to_conserv(sol_w_l)
    el    = sol_l(5)/rhol

    rhor  = sol_w_r(1)
    vn_r  = dot_product(sol_w_r(2:4), n)
    pr    = sol_w_r(5)
    vt_r  = sol_w_r(2:4) - vn_r*n
    sol_r = primit_to_conserv(sol_w_r)
    er    = sol_r(5)/rhor

    fl(1)   = vn_l*sol_l(1)
    fl(2:4) = vn_l*sol_l(2:4) + pl*n
    fl(5)   = (sol_l(5) + pl)*vn_l

    fr(1)   = vn_r*sol_r(1)
    fr(2:4) = vn_r*sol_r(2:4) + pr*n
    fr(5)   = (sol_r(5) + pr)*vn_r

    vn_bar = (lambda_l*vn_l + lambda_r*vn_r - (pr - pl))/(lambda_r + lambda_l)

    vn_l_et  = vn_l - (p_nodal - pl)/lambda_l
    rhol_et  = 1.0_DOUBLE/(1.0_DOUBLE/rhol + (vn_bar - vn_l)/lambda_l)
    sol_l_et(1)   = rhol_et
    sol_l_et(2:4) = rhol_et*(vt_l + vn_l_et*n)
    sol_l_et(5)   = rhol_et*(el + (pl*vn_l - p_nodal*vn_l_et)/lambda_l)

    vn_r_et  = vn_r + (p_nodal - pr)/lambda_r
    rhor_et  = 1.0_DOUBLE/(1.0_DOUBLE/rhor + (vn_r - vn_bar)/lambda_r)
    sol_r_et(1)   = rhor_et
    sol_r_et(2:4) = rhor_et*(vt_r + vn_r_et*n)
    sol_r_et(5)   = rhor_et*(er + (p_nodal*vn_r_et - pr*vn_r)/lambda_r)

    sl = vn_l - lambda_l/rhol
    sr = vn_r + lambda_r/rhor

    lr_flux(:, 1) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(vn_bar)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et))
    lr_flux(:, 2) = -lr_flux(:, 1)
  end subroutine multi_point_pressure_ph

  subroutine compute_lambdas_and_solve_nodal_pressure(mesh, id_vert, sol_w_lr, lambda, p_bars, p_node)
    use ns_global_data_module, only: bc_style
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, 2, mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(in) :: sol_w_lr
    real(kind=DOUBLE), dimension(2, mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(inout) :: lambda
    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(inout) :: p_bars
    real(kind=DOUBLE), intent(inout) :: p_node

    integer(kind=ENTIER) :: j, id_sub_face, le, re, iter, n_iter
    real(kind=DOUBLE) :: rhol, vn_l, pl, al, lambda_l
    real(kind=DOUBLE) :: rhor, vn_r, pr, ar, lambda_r
    real(kind=DOUBLE) :: s1, s2

    n_iter = 10
    p_node = 0.0_DOUBLE
    do iter = 1, n_iter
      s1 = 0.0_DOUBLE
      s2 = 0.0_DOUBLE

      do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
        id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
        le = mesh%face(mesh%sub_face(id_sub_face)%mesh_face)%left_neigh
        re = mesh%face(mesh%sub_face(id_sub_face)%mesh_face)%right_neigh

        rhol    = sol_w_lr(1, 1, j)
        vn_l    = dot_product(sol_w_lr(2:4, 1, j), mesh%sub_face(id_sub_face)%norm)
        pl      = sol_w_lr(5, 1, j)
        al      = sound_speed_w(sol_w_lr(:, 1, j))

        rhor    = sol_w_lr(1, 2, j)
        vn_r    = dot_product(sol_w_lr(2:4, 2, j), mesh%sub_face(id_sub_face)%norm)
        pr      = sol_w_lr(5, 2, j)
        ar      = sound_speed_w(sol_w_lr(:, 2, j))

        lambda_l = lambda(1, j)
        lambda_r = lambda(2, j)

        if (iter == 1) then
          lambda_l = max(lambda_l, al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(vn_r - vn_l))
          lambda_r = max(lambda_r, ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(vn_r - vn_l))

          p_bars(j) = (lambda_r*pl + lambda_l*pr - lambda_r*lambda_l*(vn_r - vn_l))/(lambda_r + lambda_l)
          lambda_l = max(lambda_l, al*rhol + 1.5_DOUBLE*sqrt(max(0.0_DOUBLE, (p_bars(j) - pl)*rhol)))
          lambda_r = max(lambda_r, ar*rhor + 1.5_DOUBLE*sqrt(max(0.0_DOUBLE, (p_bars(j) - pr)*rhor)))
        else
          if (bc_style == 1 .and. is_wall(re)) then
            lambda_l = max(lambda_l, al*rhol + 1.5_DOUBLE*sqrt(max(0.0_DOUBLE, (p_bars(j) - pl)*rhol)))
            lambda_r = max(lambda_r, ar*rhor + 1.5_DOUBLE*sqrt(max(0.0_DOUBLE, (p_bars(j) - pr)*rhor)))
          else
            lambda_l = max(lambda_l, al*rhol + 1.5_DOUBLE*sqrt(max(0.0_DOUBLE, (p_node - pl)*rhol)))
            lambda_r = max(lambda_r, ar*rhor + 1.5_DOUBLE*sqrt(max(0.0_DOUBLE, (p_node - pr)*rhor)))
          end if
        end if

        p_bars(j) = (lambda_r*pl + lambda_l*pr - lambda_r*lambda_l*(vn_r - vn_l))/(lambda_r + lambda_l)

        if (bc_style == 1 .and. is_wall(re)) then
          s1 = s1 + (1.0_DOUBLE/lambda_r + 1.0_DOUBLE/lambda_l)*mesh%sub_face(id_sub_face)%area &
            *(lambda_r*pl + lambda_r*lambda_l*vn_l)/(lambda_r + lambda_l)
          s2 = s2 + (1.0_DOUBLE/lambda_l)*mesh%sub_face(id_sub_face)%area
        else
          s1 = s1 + (1.0_DOUBLE/lambda_r + 1.0_DOUBLE/lambda_l)*mesh%sub_face(id_sub_face)%area*p_bars(j)
          s2 = s2 + (1.0_DOUBLE/lambda_r + 1.0_DOUBLE/lambda_l)*mesh%sub_face(id_sub_face)%area
        end if

        lambda(1, j) = lambda_l
        lambda(2, j) = lambda_r
      end do

      p_node = s1/s2
    end do
  end subroutine compute_lambdas_and_solve_nodal_pressure

  subroutine compute_lambdas_and_solve_nodal_velocity_new(ng, weight, norm, &
      lambda_l, lambda_r, sol_w_l, sol_w_r, v_node)
    use ns_global_data_module, only: bc_style, boundary_2d
    use linear_solver_module
    implicit none

    integer(kind=ENTIER), intent(in) :: ng
    real(kind=DOUBLE), dimension(ng), intent(in) :: weight
    real(kind=DOUBLE), dimension(ng), intent(inout) :: lambda_l, lambda_r
    real(kind=DOUBLE), dimension(5, ng), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(3, ng), intent(in) :: norm
    real(kind=DOUBLE), dimension(3), intent(inout) :: v_node

    integer(kind=ENTIER) :: iter, j
    real(kind=DOUBLE) :: rhol, vn_l, pl, al
    real(kind=DOUBLE) :: rhor, vn_r, pr, ar
    real(kind=DOUBLE) :: v_bar, vn_star
    real(kind=DOUBLE), dimension(3) :: rhs
    real(kind=DOUBLE), dimension(3, 3) :: mat

    !Init lambdas with two-point positivity conditions 
    !and Duckowitz correction
    mat(:, :) = 0.0_DOUBLE
    rhs(:) = 0.0_DOUBLE
    do j=1, ng
      rhol = sol_w_l(1, j)
      vn_l = dot_product(sol_w_l(2:4, j), norm(:, j))
      pl = sol_w_l(5, j)
      al = sound_speed_w(sol_w_l(:, j))

      rhor = sol_w_r(1, j)
      vn_r = dot_product(sol_w_r(2:4, j), norm(:, j))
      pr = sol_w_r(5, j)
      ar = sound_speed_w(sol_w_r(:, j))

      lambda_l(j) = max(al*rhol, &
        sqrt(rhol*max(0.0_DOUBLE, pr - pl)), &
        -rhol*(vn_r - vn_l))
      lambda_r(j) = max(ar*rhor, &
        sqrt(rhor*max(0.0_DOUBLE, pl - pr)), &
        -rhor*(vn_r - vn_l))
      v_bar = (lambda_l(j)*vn_l + lambda_r(j)*vn_r - (pr - pl))&
        /(lambda_r(j) + lambda_l(j))

      !Duckowitz correction
      lambda_l(j) = max(lambda_l(j), al*rhol*(1.0_DOUBLE &
        + 1.5_DOUBLE*max(0.0_DOUBLE, -(v_bar - vn_l)/al)))
      lambda_r(j) = max(lambda_r(j), ar*rhor*(1.0_DOUBLE &
        + 1.5_DOUBLE*max(0.0_DOUBLE, (v_bar - vn_r)/ar)))
      v_bar = (lambda_l(j)*vn_l + lambda_r(j)*vn_r - (pr - pl))&
        /(lambda_r(j) + lambda_l(j))

      mat = mat + weight(j)*(lambda_r(j)+lambda_l(j))&
        *tensor_product_3(norm(:, j), norm(:, j))
      rhs = rhs + weight(j)*(lambda_r(j)+lambda_l(j))&
        *v_bar*norm(:, j)
    end do
    call pseudo_inverse_inplace_lapack(3, mat)
    v_node = matmul(mat, rhs)

    iter = 0
    do while (iter < 2)
      iter = iter + 1

      mat(:, :) = 0.0_DOUBLE
      rhs(:) = 0.0_DOUBLE
      do j = 1, ng
        rhol = sol_w_l(1, j)
        vn_l = dot_product(sol_w_l(2:4, j), norm(:, j))
        pl = sol_w_l(5, j)
        al = sound_speed_w(sol_w_l(:, j))

        rhor = sol_w_r(1, j)
        vn_r = dot_product(sol_w_r(2:4, j), norm(:, j))
        pr = sol_w_r(5, j)
        ar = sound_speed_w(sol_w_r(:, j))

        vn_star = dot_product(v_node, norm(:, j))

        lambda_l(j) = max(lambda_l(j), al*rhol*(1.0_DOUBLE &
          + 1.5_DOUBLE*max(0.0_DOUBLE, -(vn_star - vn_l)/al)))
        lambda_r(j) = max(lambda_r(j), ar*rhor*(1.0_DOUBLE &
          + 1.5_DOUBLE*max(0.0_DOUBLE, (vn_star - vn_r)/ar)))
        v_bar = (lambda_l(j)*vn_l + lambda_r(j)*vn_r - (pr - pl)) &
          /(lambda_r(j) + lambda_l(j))

        mat = mat + weight(j)*(lambda_r(j)+lambda_l(j))&
          *tensor_product_3(norm(:, j), norm(:, j))
        rhs = rhs + weight(j)*(lambda_r(j)+lambda_l(j))&
          *v_bar*norm(:, j)
      end do
      call pseudo_inverse_inplace_lapack(3, mat)
      v_node = matmul(mat, rhs)

    end do
  end subroutine compute_lambdas_and_solve_nodal_velocity_new
end module ns_euler_rs_module