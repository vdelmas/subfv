module ns_euler_zb_module
  use precision_module
  use mesh_module
  use ns_euler_primitives_module
  use ns_euler_recon_module
  implicit none

contains
  subroutine compute_nodal_velocity_LS(mesh, id_vert, sol, grad, vp, h_p, second_order)
    use ns_global_data_module, only : boundary_2d, scheme
    use linear_solver_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: h_p
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, le, re, k
    integer(kind=ENTIER) :: id_sub_elem, id_elem
    integer(kind=ENTIER) :: id_sub_face, id_face
    real(kind=DOUBLE) :: rho_p, a_p
    real(kind=DOUBLE), dimension(3) :: grad_p, Bp
    real(kind=DOUBLE), dimension(5) :: sol_w, sol_l, sol_r
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(3, 3) :: mat

    vp(:, id_vert) = 0.0_DOUBLE
    rho_p = 0.0_DOUBLE
    a_p = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_w = conserv_to_primit(sol(:, id_elem))
      if( second_order ) sol_w = sol_w + matmul(transpose(grad(:, :, id_elem)), &
        mesh%vert(id_vert)%coord - mesh%elem(id_elem)%coord)
      vp(:, id_vert) = vp(:, id_vert) &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(2:4)
      rho_p = rho_p &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(1)
      a_p = a_p &
        + mesh%sub_elem(id_sub_elem)%volume &
        * sound_speed_w(sol_w)
    end do
    rho_p = rho_p / mesh%vert(id_vert)%volume
    a_p = a_p / mesh%vert(id_vert)%volume
    vp(:, id_vert) = vp(:, id_vert) / mesh%vert(id_vert)%volume

    grad_p = 0.0_DOUBLE
    mat = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)
      if( re > 0 ) then
        grad_p = grad_p + (sol_w_r(5) - sol_w_l(5)) &
          * mesh%sub_face(id_sub_face)%area&
          * mesh%sub_face(id_sub_face)%norm
        mat = mat + mesh%sub_face(id_sub_face)%area&
          *tensor_product(mesh%sub_face(id_sub_face)%norm, &
          mesh%elem(re)%coord - mesh%elem(le)%coord)
      end if
    end do
    !grad_p = grad_p / mesh%vert(id_vert)%volume
    call pseudo_inverse_inplace_lapack(3, mat)
    grad_p = matmul(mat, mesh%vert(id_vert)%volume*grad_p)

    if( mesh%vert(id_vert)%is_bound ) then
      Bp = wall_normal(mesh, id_vert)
      if( norm2(Bp) > 1e-12_DOUBLE ) then
        Bp = Bp/norm2(Bp)
        grad_p = grad_p - dot_product(grad_p, Bp)*Bp
        vp(:, id_vert) = vp(:, id_vert) - dot_product(vp(:, id_vert), Bp)*Bp
      end if
    end if

    vp(:, id_vert) = vp(:, id_vert) - 0.5_DOUBLE*h_p(id_vert)/(rho_p*a_p)*grad_p
  end subroutine compute_nodal_velocity_LS

  ! removed 2026-09-15: compute_nodal_velocity_LSM (dead code, unreachable via schemes.txt, depended on compute_ellip)

  subroutine compute_nodal_velocity_LSU(mesh, id_vert, sol, grad, vp, second_order)
    use ns_global_data_module, only : boundary_2d, scheme
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, le, re
    integer(kind=ENTIER) :: id_sub_elem, id_elem
    integer(kind=ENTIER) :: id_sub_face, id_face
    real(kind=DOUBLE) :: rho_p, a_p
    real(kind=DOUBLE), dimension(3) :: grad_p, Bp
    real(kind=DOUBLE), dimension(5) :: sol_w, sol_l, sol_r
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r

    vp(:, id_vert) = 0.0_DOUBLE
    rho_p = 0.0_DOUBLE
    a_p = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_w = conserv_to_primit(sol(:, id_elem))
      if( second_order ) sol_w = sol_w + matmul(transpose(grad(:, :, id_elem)), &
        mesh%vert(id_vert)%coord - mesh%elem(id_elem)%coord)
      vp(:, id_vert) = vp(:, id_vert) &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(2:4)
      rho_p = rho_p &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(1)
      a_p = a_p &
        + mesh%sub_elem(id_sub_elem)%volume &
        * sound_speed_w(sol_w)
    end do
    rho_p = rho_p / mesh%vert(id_vert)%volume
    a_p = a_p / mesh%vert(id_vert)%volume
    vp(:, id_vert) = vp(:, id_vert) / mesh%vert(id_vert)%volume

    grad_p = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)
      grad_p = grad_p + (sol_w_r(5) - sol_w_l(5)) &
        * mesh%sub_face(id_sub_face)%norm
    end do

    if( mesh%vert(id_vert)%is_bound ) then
      Bp = wall_normal(mesh, id_vert)
      if( norm2(Bp) > 1e-12_DOUBLE ) then
        Bp = Bp/norm2(Bp)
        grad_p = grad_p - dot_product(grad_p, Bp)*Bp
        vp(:, id_vert) = vp(:, id_vert) - dot_product(vp(:, id_vert), Bp)*Bp
      end if
    end if

    vp(:, id_vert) = vp(:, id_vert) - 0.5_DOUBLE/(rho_p*a_p)*grad_p
  end subroutine compute_nodal_velocity_LSU

  subroutine compute_nodal_pressure_LS(mesh, id_vert, sol, grad, pp, h_p, second_order)
    use ns_global_data_module, only : scheme, boundary_2d
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: h_p
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    real(kind=DOUBLE), intent(inout) :: pp
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, id_elem, id_sub_elem, le, re
    integer(kind=ENTIER) :: id_face, id_sub_face
    real(kind=DOUBLE), dimension(5) :: sol_w
    real(kind=DOUBLE) :: rho_p, a_p, div_v
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, sol_l, sol_r

    pp = 0.0_DOUBLE
    rho_p = 0.0_DOUBLE
    a_p = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_w = conserv_to_primit(sol(:, id_elem))
      if( second_order ) sol_w = sol_w + matmul(transpose(grad(:, :, id_elem)), &
        mesh%vert(id_vert)%coord - mesh%elem(id_elem)%coord)
      pp = pp &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(5)
      rho_p = rho_p &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(1)
      a_p = a_p &
        + mesh%sub_elem(id_sub_elem)%volume &
        * sound_speed_w(sol_w)
    end do
    pp = pp / mesh%vert(id_vert)%volume
    rho_p = rho_p / mesh%vert(id_vert)%volume
    a_p = a_p / mesh%vert(id_vert)%volume

    div_v = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)
      div_v = div_v &
        + dot_product(sol_w_r(2:4) - sol_w_l(2:4), &
        mesh%sub_face(id_sub_face)%area&
        * h_p(id_vert) * mesh%sub_face(id_sub_face)%norm)
    end do
    div_v = div_v / mesh%vert(id_vert)%volume

    pp = pp - 0.5_DOUBLE*rho_p*a_p*div_v
  end subroutine compute_nodal_pressure_LS

  subroutine compute_nodal_pressure_LSU(mesh, id_vert, sol, grad, pp, second_order)
    use ns_global_data_module, only : scheme
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    real(kind=DOUBLE), intent(inout) :: pp
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, id_elem, id_sub_elem, le, re
    integer(kind=ENTIER) :: id_face, id_sub_face
    real(kind=DOUBLE), dimension(5) :: sol_w
    real(kind=DOUBLE) :: rho_p, a_p, div_v
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, sol_l, sol_r

    pp = 0.0_DOUBLE
    rho_p = 0.0_DOUBLE
    a_p = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_w = conserv_to_primit(sol(:, id_elem))
      if( second_order ) sol_w = sol_w + matmul(transpose(grad(:, :, id_elem)), &
        mesh%vert(id_vert)%coord - mesh%elem(id_elem)%coord)
      pp = pp &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(5)
      rho_p = rho_p &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(1)
      a_p = a_p &
        + mesh%sub_elem(id_sub_elem)%volume &
        * sound_speed_w(sol_w)
    end do
    pp = pp / mesh%vert(id_vert)%volume
    rho_p = rho_p / mesh%vert(id_vert)%volume
    a_p = a_p / mesh%vert(id_vert)%volume

    div_v = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)
      div_v = div_v &
        + dot_product(sol_w_r(2:4) - sol_w_l(2:4), &
        mesh%sub_face(id_sub_face)%norm)
    end do

    pp = pp - 0.5_DOUBLE*rho_p*a_p*div_v
  end subroutine compute_nodal_pressure_LSU

  ! removed 2026-09-15: compute_nodal_pressure_LSM (dead code, unreachable via schemes.txt, depended on compute_ellip)

  subroutine compute_rhs_around_vert_ARMDMAT(mesh, sol, grad, &
      nsen, sum_lambda_vert, flux_sum_vert, &
      second_order, id_vert, vp, h_p)
    use ns_global_data_module, only: bc_style, scheme, &
      exclude_bound_vert, boundary_2d
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: h_p
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    integer(kind=ENTIER), intent(in) :: nsen
    real(kind=DOUBLE), dimension(nsen), intent(inout) :: &
      sum_lambda_vert
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: &
      flux_sum_vert
    logical, intent(in) :: second_order
    integer(kind=ENTIER), intent(in) :: id_vert

    integer(kind=ENTIER) :: j, id_sub_face, id_face
    integer(kind=ENTIER) :: id_elem, id_sub_elem
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE), dimension(3) :: norm

    real(kind=DOUBLE) :: vnl, vnr, al, ar, lambda, wpcf
    real(kind=DOUBLE), dimension(5) :: sol_p, sol_l, sol_r, ff
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, fminus, fplus
    real(kind=DOUBLE), dimension(5,3) :: fp, grad_sol
    real(kind=DOUBLE), dimension(3,3) :: SMAX

    rse_loc = 0
    sol_p = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_p = sol_p + mesh%sub_elem(id_sub_elem)%volume * sol(:, id_elem)
    end do
    sol_p = sol_p / mesh%vert(id_vert)%volume

    grad_sol = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm
      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)
      grad_sol = grad_sol + tensor_product(primit_to_conserv(sol_w_r) - primit_to_conserv(sol_w_l), &
        mesh%sub_face(id_sub_face)%area*mesh%sub_face(id_sub_face)%norm)
    end do
    grad_sol = grad_sol / mesh%vert(id_vert)%volume

    !call compute_nodal_velocity_LSM(mesh, id_vert, sol, grad, vp(:, id_vert), mat_h_p, second_order)
    call compute_nodal_velocity_LS(mesh, id_vert, sol, grad, vp, h_p, second_order)

    SMAX = 0.0_DOUBLE
    SMAX(1, 1) = abs(vp(1, id_vert))
    SMAX(2, 2) = abs(vp(2, id_vert))
    SMAX(3, 3) = abs(vp(3, id_vert))
    fp = tensor_product(sol_p, vp(:, id_vert)) &
      - 0.5_DOUBLE*h_p(id_vert)*matmul(grad_sol, SMAX)

    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      vnl = dot_product(sol_w_l(2:4), norm)
      vnr = dot_product(sol_w_r(2:4), norm)

      al = sound_speed_w(sol_w_l)
      ar = sound_speed_w(sol_w_r)

      sol_r = primit_to_conserv(sol_w_r)
      sol_l = primit_to_conserv(sol_w_l)

      ff = 0.5_DOUBLE*(vnr*sol_r + vnl*sol_l) &
        - 0.5_DOUBLE*(sol_r - sol_l)*max(abs(vnr), abs(vnl))

      wpcf = 0.5_DOUBLE
      fminus = wpcf*matmul(fp, norm) + (1.0_DOUBLE-wpcf)*ff
      fplus = fminus

      !Used for local timestepping
      lambda = max(1e-8_DOUBLE, -vnl, vnr)+max(al,ar)

      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        sum_lambda_vert(lse_loc) = sum_lambda_vert(lse_loc) &
          + mesh%sub_face(id_sub_face)%area*lambda
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) + mesh%sub_face(id_sub_face)%area*fminus

        if (rse > 0) then
        if (mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_ARMDMAT

  subroutine compute_rhs_around_vert_ARMDUMAT(mesh, sol, grad, &
      nsen, sum_lambda_vert, flux_sum_vert, &
      second_order, id_vert, vp)
    use ns_global_data_module, only: bc_style, scheme, &
      exclude_bound_vert, boundary_2d
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    integer(kind=ENTIER), intent(in) :: nsen
    real(kind=DOUBLE), dimension(nsen), intent(inout) :: &
      sum_lambda_vert
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: &
      flux_sum_vert
    logical, intent(in) :: second_order
    integer(kind=ENTIER), intent(in) :: id_vert

    integer(kind=ENTIER) :: j, id_sub_face, id_face
    integer(kind=ENTIER) :: id_elem, id_sub_elem
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE), dimension(3) :: norm

    real(kind=DOUBLE) :: vnl, vnr, al, ar, lambda, wpcf
    real(kind=DOUBLE), dimension(5) :: sol_p, sol_l, sol_r, ff
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, fminus, fplus
    real(kind=DOUBLE), dimension(5,3) :: fp, grad_sol
    real(kind=DOUBLE), dimension(3,3) :: SMAX

    rse_loc = 0
    sol_p = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_p = sol_p + mesh%sub_elem(id_sub_elem)%volume * sol(:, id_elem)
    end do
    sol_p = sol_p / mesh%vert(id_vert)%volume

    grad_sol = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm
      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)
      grad_sol = grad_sol + tensor_product(primit_to_conserv(sol_w_r) - primit_to_conserv(sol_w_l), &
        mesh%sub_face(id_sub_face)%norm)
    end do

    call compute_nodal_velocity_LSU(mesh, id_vert, sol, grad, vp, second_order)

    SMAX = 0.0_DOUBLE
    SMAX(1, 1) = abs(vp(1, id_vert))
    SMAX(2, 2) = abs(vp(2, id_vert))
    SMAX(3, 3) = abs(vp(3, id_vert))
    fp = tensor_product(sol_p, vp(:, id_vert)) &
      - 0.5_DOUBLE*matmul(grad_sol,SMAX)

    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      vnl = dot_product(sol_w_l(2:4), norm)
      vnr = dot_product(sol_w_r(2:4), norm)

      al = sound_speed_w(sol_w_l)
      ar = sound_speed_w(sol_w_r)

      sol_r = primit_to_conserv(sol_w_r)
      sol_l = primit_to_conserv(sol_w_l)

      ff = 0.5_DOUBLE*(vnr*sol_r + vnl*sol_l) &
        - 0.5_DOUBLE*(sol_r - sol_l)*max(abs(vnr), abs(vnl))

      wpcf = 0.5_DOUBLE
      fminus = wpcf*matmul(fp, norm) + (1.0_DOUBLE-wpcf)*ff
      fplus = fminus

      !Used for local timestepping
      lambda = max(1e-8_DOUBLE, -vnl, vnr)+max(al,ar)

      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        sum_lambda_vert(lse_loc) = sum_lambda_vert(lse_loc) &
          + mesh%sub_face(id_sub_face)%area*lambda
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) + mesh%sub_face(id_sub_face)%area*fminus

        if (rse > 0) then
        if (mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_ARMDUMAT

  ! removed 2026-09-15: compute_rhs_around_vert_ARMDMMAT (dead code, unreachable via schemes.txt, depended on compute_ellip)

  subroutine compute_rhs_around_vert_ARMD(mesh, sol, grad, &
      nsen, sum_lambda_vert, flux_sum_vert, &
      second_order, id_vert, vp, h_p)
    use ns_global_data_module, only: bc_style, scheme, &
      exclude_bound_vert, boundary_2d
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: h_p
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    integer(kind=ENTIER), intent(in) :: nsen
    real(kind=DOUBLE), dimension(nsen), intent(inout) :: &
      sum_lambda_vert
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: &
      flux_sum_vert
    logical, intent(in) :: second_order
    integer(kind=ENTIER), intent(in) :: id_vert

    integer(kind=ENTIER) :: j, id_sub_face, id_face
    integer(kind=ENTIER) :: id_elem, id_sub_elem
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE), dimension(3) :: norm

    real(kind=DOUBLE) :: vnl, vnr, al, ar, lambda, wpcf
    real(kind=DOUBLE), dimension(5) :: sol_p, sol_l, sol_r, ff
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, fminus, fplus
    real(kind=DOUBLE), dimension(5,3) :: fp, grad_sol

    rse_loc = 0
    sol_p = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_p = sol_p + mesh%sub_elem(id_sub_elem)%volume * sol(:, id_elem)
    end do
    sol_p = sol_p / mesh%vert(id_vert)%volume

    grad_sol = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm
      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)
      grad_sol = grad_sol + tensor_product(primit_to_conserv(sol_w_r) - primit_to_conserv(sol_w_l), &
        mesh%sub_face(id_sub_face)%area*mesh%sub_face(id_sub_face)%norm)
    end do
    grad_sol = grad_sol / mesh%vert(id_vert)%volume

    call compute_nodal_velocity_LS(mesh, id_vert, sol, grad, vp, h_p, second_order)

    fp = tensor_product(sol_p, vp(:, id_vert)) &
      - 0.5_DOUBLE*norm2(vp(:, id_vert))*h_p(id_vert)*grad_sol

    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      vnl = dot_product(sol_w_l(2:4), norm)
      vnr = dot_product(sol_w_r(2:4), norm)

      al = sound_speed_w(sol_w_l)
      ar = sound_speed_w(sol_w_r)

      sol_r = primit_to_conserv(sol_w_r)
      sol_l = primit_to_conserv(sol_w_l)

      ff = 0.5_DOUBLE*(vnr*sol_r + vnl*sol_l) &
        - 0.5_DOUBLE*(sol_r - sol_l)*max(abs(vnr), abs(vnl))

      wpcf = 0.5_DOUBLE
      fminus = wpcf*matmul(fp, norm) + (1.0_DOUBLE-wpcf)*ff
      fplus = fminus

      !Used for local timestepping
      lambda = max(1e-8_DOUBLE, -vnl, vnr)+max(al,ar)

      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        sum_lambda_vert(lse_loc) = sum_lambda_vert(lse_loc) &
          + mesh%sub_face(id_sub_face)%area*lambda
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) + mesh%sub_face(id_sub_face)%area*fminus

        if (rse > 0) then
        if (mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_ARMD

  subroutine compute_rhs_around_vert_ARMDU(mesh, sol, grad, &
      nsen, sum_lambda_vert, flux_sum_vert, &
      second_order, id_vert, vp)
    use ns_global_data_module, only: bc_style, scheme, &
      exclude_bound_vert, boundary_2d
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    integer(kind=ENTIER), intent(in) :: nsen
    real(kind=DOUBLE), dimension(nsen), intent(inout) :: &
      sum_lambda_vert
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: &
      flux_sum_vert
    logical, intent(in) :: second_order
    integer(kind=ENTIER), intent(in) :: id_vert

    integer(kind=ENTIER) :: j, id_sub_face, id_face
    integer(kind=ENTIER) :: id_elem, id_sub_elem
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE), dimension(3) :: norm

    real(kind=DOUBLE) :: vnl, vnr, al, ar, lambda, wpcf
    real(kind=DOUBLE), dimension(5) :: sol_p, sol_l, sol_r, ff
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, fminus, fplus
    real(kind=DOUBLE), dimension(5,3) :: fp, grad_sol

    rse_loc = 0
    sol_p = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_p = sol_p + mesh%sub_elem(id_sub_elem)%volume * sol(:, id_elem)
    end do
    sol_p = sol_p / mesh%vert(id_vert)%volume

    grad_sol = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm
      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)
      grad_sol = grad_sol + tensor_product(primit_to_conserv(sol_w_r) - primit_to_conserv(sol_w_l), &
        mesh%sub_face(id_sub_face)%norm)
    end do

    call compute_nodal_velocity_LSU(mesh, id_vert, sol, grad, vp, second_order)
    fp = tensor_product(sol_p, vp(:, id_vert)) &
      - 0.5_DOUBLE*norm2(vp(:, id_vert))*grad_sol

    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      vnl = dot_product(sol_w_l(2:4), norm)
      vnr = dot_product(sol_w_r(2:4), norm)

      al = sound_speed_w(sol_w_l)
      ar = sound_speed_w(sol_w_r)

      sol_r = primit_to_conserv(sol_w_r)
      sol_l = primit_to_conserv(sol_w_l)

      ff = 0.5_DOUBLE*(vnr*sol_r + vnl*sol_l) &
        - 0.5_DOUBLE*(sol_r - sol_l)*max(abs(vnr), abs(vnl))

      wpcf = 0.5_DOUBLE
      fminus = wpcf*matmul(fp, norm) + (1.0_DOUBLE-wpcf)*ff
      fplus = fminus

      !Used for local timestepping
      lambda = max(1e-8_DOUBLE, -vnl, vnr)+max(al,ar)

      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        sum_lambda_vert(lse_loc) = sum_lambda_vert(lse_loc) &
          + mesh%sub_face(id_sub_face)%area*lambda
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) + mesh%sub_face(id_sub_face)%area*fminus

        if (rse > 0) then
        if (mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_ARMDU

  ! removed 2026-09-15: compute_rhs_around_vert_ARMDM (dead code, unreachable via schemes.txt, depended on compute_ellip)

  subroutine compute_rhs_around_vert_AR1D(mesh, sol, grad, &
      nsen, sum_lambda_vert, flux_sum_vert, &
      second_order, id_vert)
    use ns_global_data_module, only: bc_style, scheme, &
      exclude_bound_vert, boundary_2d
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    integer(kind=ENTIER), intent(in) :: nsen
    real(kind=DOUBLE), dimension(nsen), intent(inout) :: &
      sum_lambda_vert
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: &
      flux_sum_vert
    logical, intent(in) :: second_order
    integer(kind=ENTIER), intent(in) :: id_vert

    integer(kind=ENTIER) :: j, id_sub_face, id_face
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE), dimension(3) :: norm

    real(kind=DOUBLE) :: vnl, vnr, al, ar, lambda, corr, corr_check, am, vm, machm
    real(kind=DOUBLE), dimension(5) :: sol_l, sol_r
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, fminus, fplus

    rse_loc = 0
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      vnl = dot_product(sol_w_l(2:4), norm)
      vnr = dot_product(sol_w_r(2:4), norm)

      al = sound_speed_w(sol_w_l)
      ar = sound_speed_w(sol_w_r)

      am = 0.5_DOUBLE*(al+ar)
      vm = 0.5_DOUBLE*(norm2(sol_w_l(2:4))+norm2(sol_w_r(2:4)))

      machm = vm/am

      lambda = max(1e-8_DOUBLE, abs(vnl), abs(vnr))
      !lambda = max(1e-8_DOUBLE, abs(vnl), abs(vnr)) + am
      !lambda = max(1e-8_DOUBLE, abs(vnl), abs(vnr), abs(vm))
      sol_l = primit_to_conserv(sol_w_l)
      sol_r = primit_to_conserv(sol_w_r)
      fminus = 0.5_DOUBLE*(vnl * sol_l + vnr * sol_r) - 0.5*lambda*(sol_r - sol_l)
      fplus = fminus

      !Used for local timestepping
      lambda = max(1e-8_DOUBLE, -vnl, vnr)+max(al,ar)

      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        sum_lambda_vert(lse_loc) = sum_lambda_vert(lse_loc) &
          + mesh%sub_face(id_sub_face)%area*lambda
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) + mesh%sub_face(id_sub_face)%area*fminus

        if (rse > 0) then
        if (mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_AR1D

  subroutine compute_rhs_around_vert_LSU(mesh, sol, grad, &
      nsen, flux_sum_vert, &
      id_vert, vp, second_order)
    use ns_global_data_module, only: bc_style, scheme, &
      exclude_bound_vert, boundary_2d
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    integer(kind=ENTIER), intent(in) :: nsen
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: &
      flux_sum_vert
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, id_sub_face, id_face
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE), dimension(3) :: norm

    real(kind=DOUBLE) :: pp
    real(kind=DOUBLE), dimension(5) :: fminus, fplus
    real(kind=DOUBLE), dimension(5,3) :: fp

    rse_loc = 0
    call compute_nodal_velocity_LSU(mesh, id_vert, sol, grad, vp, second_order)
    call compute_nodal_pressure_LSU(mesh, id_vert, sol, grad, pp, second_order)

    fp(1, :) = 0.0_DOUBLE
    fp(2:4, :) = pp*eye3
    fp(5, :) = pp*vp(:, id_vert)

    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm

      fminus = matmul(fp, norm)
      fplus = fminus

      if( boundary_2d &
        .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1e-3 ) then
        fminus = 0.0_DOUBLE
        fplus = 0.0_DOUBLE
      end if

      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) + mesh%sub_face(id_sub_face)%area*fminus

        if (rse > 0) then
        if (mesh%sub_elem(rse)%mesh_vert == id_vert) then
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_LSU

  subroutine compute_rhs_around_vert_LS(mesh, sol, grad, &
      nsen, flux_sum_vert, &
      id_vert, vp, h_p, second_order)
    use ns_global_data_module, only: bc_style, scheme, &
      exclude_bound_vert, boundary_2d
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: h_p
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    integer(kind=ENTIER), intent(in) :: nsen
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: &
      flux_sum_vert
    integer(kind=ENTIER), intent(in) :: id_vert
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, id_sub_face, id_face
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE), dimension(3) :: norm

    real(kind=DOUBLE) :: pp
    real(kind=DOUBLE), dimension(5) :: fminus, fplus
    real(kind=DOUBLE), dimension(5,3) :: fp

    rse_loc = 0
    call compute_nodal_velocity_LS(mesh, id_vert, sol, grad, vp, h_p, second_order)
    call compute_nodal_pressure_LS(mesh, id_vert, sol, grad, pp, h_p, second_order)

    fp(1, :) = 0.0_DOUBLE
    fp(2:4, :) = pp*eye3
    fp(5, :) = pp*vp(:, id_vert)

    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm

      fminus = matmul(fp, norm)
      fplus = fminus

      if( boundary_2d &
        .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1e-3 ) then
        fminus = 0.0_DOUBLE
        fplus = 0.0_DOUBLE
      end if

      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) + mesh%sub_face(id_sub_face)%area*fminus

        if (rse > 0) then
        if (mesh%sub_elem(rse)%mesh_vert == id_vert) then
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_LS

  ! removed 2026-09-15: compute_rhs_around_vert_LSM (dead code, unreachable via schemes.txt, depended on compute_ellip)

  subroutine compute_rhs_around_vert_AM(mesh, sol, grad, &
      nsen, u_vert, sum_lambda_vert, flux_sum_vert, &
      second_order, id_vert, h_p)
    use ns_global_data_module, only: bc_style, scheme, &
      exclude_bound_vert, boundary_2d
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    integer(kind=ENTIER), intent(in) :: nsen
    real(kind=DOUBLE), dimension(3), intent(inout) :: u_vert
    real(kind=DOUBLE), dimension(nsen), intent(inout) :: &
      sum_lambda_vert
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: &
      flux_sum_vert
    logical, intent(in) :: second_order
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: h_p

    integer(kind=ENTIER) :: j, id_sub_face, id_face
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE), dimension(3) :: norm

    real(kind=DOUBLE) :: sum_area, vnl, vnr, al, ar, lambda, corr
    real(kind=DOUBLE), dimension(5) :: sol_p, sol_l, sol_r, u_bar
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, fminus, fplus, ff


    rse_loc = 0
    call compute_corr(mesh, id_vert, sol, grad, corr, h_p, second_order)
    u_vert = 0.0_DOUBLE

    sol_p = 0.0_DOUBLE
    sum_area = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      norm = mesh%sub_face(id_sub_face)%norm
      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      vnl = dot_product(sol_w_l(2:4), norm)
      vnr = dot_product(sol_w_r(2:4), norm)

      al = sound_speed_w(sol_w_l)
      ar = sound_speed_w(sol_w_r)

      lambda = max(1e-8_DOUBLE, -vnl, vnr)+corr
      !lambda = max(1e-8_DOUBLE, -vnl, vnr)

      sol_l = primit_to_conserv(sol_w_l)
      sol_r = primit_to_conserv(sol_w_r)
      u_bar = sol_l*0.5_DOUBLE*(1.0_DOUBLE + vnl/lambda) &
        + sol_r*0.5_DOUBLE*(1.0_DOUBLE - vnr/lambda)

      if( re > 0 ) then
        sol_p = sol_p + mesh%sub_face(id_sub_face)%area * lambda * u_bar
        sum_area = sum_area + mesh%sub_face(id_sub_face)%area * lambda
      end if
    end do
    sol_p = sol_p / sum_area

    !!Compute flux across each sub_face
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm
      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      vnl = dot_product(sol_w_l(2:4), norm)
      vnr = dot_product(sol_w_r(2:4), norm)

      al = sound_speed_w(sol_w_l)
      ar = sound_speed_w(sol_w_r)

      lambda = max(1e-8_DOUBLE, -vnl, vnr)+corr
      !lambda = max(1e-8_DOUBLE, -vnl, vnr)

      sol_l = primit_to_conserv(sol_w_l)
      sol_r = primit_to_conserv(sol_w_r)

      if( re > 0 ) then
        fminus = sol_l*vnl - lambda*(sol_p - sol_l)
        fplus = sol_r*vnr + lambda*(sol_p - sol_r)
      else
        ff = 0.5_DOUBLE*(vnr*sol_r + vnl*sol_l) &
          - 0.5_DOUBLE*(sol_r - sol_l)*max(abs(vnr), abs(vnl))
        fminus = ff
        fplus = fminus
      end if

      !Used for local timestepping
      lambda = max(1e-8_DOUBLE, -vnl, vnr)+max(al,ar)

      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        sum_lambda_vert(lse_loc) = sum_lambda_vert(lse_loc) &
          + mesh%sub_face(id_sub_face)%area*lambda
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) + mesh%sub_face(id_sub_face)%area*fminus

        if (rse > 0) then
        if (mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_AM

  subroutine compute_rhs_around_vert_AMISO(mesh, sol, grad, &
      nsen, u_vert, sum_lambda_vert, flux_sum_vert, &
      second_order, id_vert, h_p)
    use ns_global_data_module, only: bc_style, scheme, &
      exclude_bound_vert, boundary_2d
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    integer(kind=ENTIER), intent(in) :: nsen
    real(kind=DOUBLE), dimension(3), intent(inout) :: u_vert
    real(kind=DOUBLE), dimension(nsen), intent(inout) :: &
      sum_lambda_vert
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: &
      flux_sum_vert
    logical, intent(in) :: second_order
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: h_p

    integer(kind=ENTIER) :: j, id_sub_face, id_face
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE), dimension(3) :: norm

    real(kind=DOUBLE) :: sum_area, vnl, vnr, al, ar, lambda, corr
    real(kind=DOUBLE) :: wpcf, min_apf, pr, pl, rhol, rhor
    real(kind=DOUBLE), dimension(5) :: sol_p, sol_l, sol_r, u_bar, ff, sol_m
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, fminus, fplus
    real(kind=DOUBLE) :: vm, rhom, am


    rse_loc = 0
    min_apf = huge(1.0_DOUBLE)
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      if( re > 0 ) then
        ! Une sous-face z de maillage faux-2D est inerte : si elle fixait le minimum
        ! d'aire, wpcf vaudrait h/t au lieu de 1 et iso cesserait d'egaler multi_point.
        if (.not. (boundary_2d .and. &
            abs(mesh%sub_face(id_sub_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE)) &
          min_apf = min(min_apf, mesh%sub_face(id_sub_face)%area)
      end if
    end do

    u_vert = 0.0_DOUBLE

    !call compute_corr(mesh, id_vert, sol, grad, corr, h_p, second_order)
    call compute_corr2(mesh, id_vert, sol, grad, corr, second_order)

    sol_p = 0.0_DOUBLE
    sum_area = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      norm = mesh%sub_face(id_sub_face)%norm
      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      vnl = dot_product(sol_w_l(2:4), norm)
      vnr = dot_product(sol_w_r(2:4), norm)

      al = sound_speed_w(sol_w_l)
      ar = sound_speed_w(sol_w_r)

      lambda = max(1e-8_DOUBLE, -vnl, vnr)+corr
      !lambda = max(1e-8_DOUBLE, -vnl, vnr)

      sol_l = primit_to_conserv(sol_w_l)
      sol_r = primit_to_conserv(sol_w_r)
      u_bar = sol_l*0.5_DOUBLE*(1.0_DOUBLE + vnl/lambda) &
        + sol_r*0.5_DOUBLE*(1.0_DOUBLE - vnr/lambda)


      if( re > 0 ) then
        wpcf = min_apf/mesh%sub_face(id_sub_face)%area
        sol_p = sol_p + mesh%sub_face(id_sub_face)%area * wpcf * lambda * u_bar
        sum_area = sum_area + mesh%sub_face(id_sub_face)%area * wpcf * lambda
      end if
    end do
    sol_p = sol_p / sum_area

    !!Compute flux across each sub_face
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm
      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      rhol = sol_w_l(1)
      vnl = dot_product(sol_w_l(2:4), norm)
      pl = sol_w_l(5)
      al = sound_speed_w(sol_w_l)

      rhor = sol_w_r(1)
      vnr = dot_product(sol_w_r(2:4), norm)
      pr = sol_w_r(5)
      ar = sound_speed_w(sol_w_r)

      lambda = max(1e-8_DOUBLE, -vnl, vnr)+corr
      !lambda = max(1e-8_DOUBLE, -vnl, vnr)

      sol_l = primit_to_conserv(sol_w_l)
      sol_r = primit_to_conserv(sol_w_r)

      !ff = 0.5_DOUBLE*(vnr*sol_r + vnl*sol_l) &
      !  - 0.5_DOUBLE*(sol_r - sol_l)*max(abs(vnr), abs(vnl))

      rhom = 0.5_DOUBLE*(rhol+rhor)
      am = 0.5_DOUBLE*(al+ar)
      !vm = (lambda*(vnr+vnl) - (pr-pl))/(2.0_DOUBLE*lambda)
      vm = 0.5_DOUBLE*(vnr+vnl) - 0.5_DOUBLE/(rhom*am)*(pr-pl)
      sol_m = 0.5_DOUBLE*(sol_r+sol_l)
      ff = vm*sol_m - 0.5_DOUBLE*(abs(vm)+corr)*(sol_r-sol_l)

      if( re > 0 ) then
        wpcf = min_apf/mesh%sub_face(id_sub_face)%area * 1.0_DOUBLE/3.0_DOUBLE
        !wpcf = min_apf/mesh%sub_face(id_sub_face)%area
        !wpcf = 0.0_DOUBLE
        fminus = wpcf*(sol_l*vnl - lambda*(sol_p - sol_l)) &
          + (1.0_DOUBLE-wpcf)*ff
        fplus = wpcf*(sol_r*vnr + lambda*(sol_p - sol_r)) &
          + (1.0_DOUBLE-wpcf)*ff
      else
        fminus = ff
        fplus = fminus
      end if

      !Used for local timestepping
      lambda = max(1e-8_DOUBLE, -vnl, vnr)+max(al,ar)

      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        sum_lambda_vert(lse_loc) = sum_lambda_vert(lse_loc) &
          + mesh%sub_face(id_sub_face)%area*lambda
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) + mesh%sub_face(id_sub_face)%area*fminus

        if (rse > 0) then
        if (mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_AMISO

  function compute_length(mesh, id_vert) result(h_p)
    use ns_global_data_module, only: boundary_2d
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE) :: h_p

    integer(kind=ENTIER) :: j, id_sub_elem
    real(kind=DOUBLE) :: area_sum

    area_sum = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      area_sum = area_sum + norm2(corner_normal(mesh, id_sub_elem))
    end do
    h_p = 2*mesh%vert(id_vert)%volume/area_sum
    if( boundary_2d ) then
      h_p = sqrt(2.)*h_p
    else 
      h_p = sqrt(3.)*h_p
    end if
  end function compute_length

  function compute_length_diag(mesh, id_vert) result(h_p)
    use ns_global_data_module, only: boundary_2d
    use linear_solver_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(3, 3) :: h_p

    integer(kind=ENTIER) :: j, k, id_sub_face, id_face
    integer(kind=ENTIER) :: id_sub_elem, id_elem
    real(kind=DOUBLE) :: limp

    real(kind=DOUBLE), dimension(3) :: ce1, ce2, cp, v3

    real(kind=DOUBLE), dimension(6) :: Hp6, Rp6, v6
    real(kind=DOUBLE), dimension(6, 6) :: Mp6

    integer(kind=ENTIER) :: idv, idvm, idvp, kp, km
    real(kind=DOUBLE), dimension(3, 3) :: Q
    real(kind=DOUBLE), dimension(3) :: Lambda

    cp = mesh%vert(id_vert)%coord

    Mp6 = 0.0_DOUBLE
    Rp6 = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face

      idvp = 0
      idvm = 0
      do k=1, mesh%face(id_face)%n_vert
        idv = mesh%face(id_face)%vert(k)
        if( idv == id_vert ) then
          kp = 1+mod((k+1)-1, mesh%face(id_face)%n_vert)
          km = 1+mod((k-1)-1+mesh%face(id_face)%n_vert, mesh%face(id_face)%n_vert)
          idvm = mesh%face(id_face)%vert(km)
          idvp = mesh%face(id_face)%vert(kp)
          exit
        end if
      end do
      ce1 = 0.5_DOUBLE*(mesh%vert(idvp)%coord + cp)
      ce2 = 0.5_DOUBLE*(mesh%vert(idvm)%coord + cp)

      !ce1
      v3 = ce1 - cp
      limp = 2*norm2(v3)
      v3 = v3 / norm2(v3)
      v6(1) = v3(1)*v3(1)
      v6(2) = 2*v3(1)*v3(2)
      v6(3) = 2*v3(1)*v3(3)
      v6(4) = v3(2)*v3(2)
      v6(5) = 2*v3(2)*v3(3)
      v6(6) = v3(3)*v3(3)
      Mp6 = Mp6 + tensor_product(v6, v6)
      Rp6 = Rp6 + limp*v6

      !ce2
      v3 = ce2 - cp
      limp = 2*norm2(v3)
      v3 = v3 / norm2(v3)
      v6(1) = v3(1)*v3(1)
      v6(2) = 2*v3(1)*v3(2)
      v6(3) = 2*v3(1)*v3(3)
      v6(4) = v3(2)*v3(2)
      v6(5) = 2*v3(2)*v3(3)
      v6(6) = v3(3)*v3(3)
      Mp6 = Mp6 + tensor_product(v6, v6)
      Rp6 = Rp6 + limp*v6

      !ce2
      if( boundary_2d ) then
        v3 = mesh%face(id_face)%coord - cp
        v3(3) = 0.0_DOUBLE
      else
        v3 = mesh%face(id_face)%coord - cp
      end if
      limp = 2*norm2(v3)
      v3 = v3 / norm2(v3)
      v6(1) = v3(1)*v3(1)
      v6(2) = 2*v3(1)*v3(2)
      v6(3) = 2*v3(1)*v3(3)
      v6(4) = v3(2)*v3(2)
      v6(5) = 2*v3(2)*v3(3)
      v6(6) = v3(3)*v3(3)
      Mp6 = Mp6 + tensor_product(v6, v6)
      Rp6 = Rp6 + limp*v6
    end do

    call pseudo_inverse_inplace_lapack(6, Mp6)
    !call inv_lapack(6, Mp6)
    Hp6 = matmul(Mp6, Rp6)

    h_p = 0.0_DOUBLE
    h_p(1, 1) = Hp6(1)
    h_p(1, 2) = Hp6(2)
    h_p(1, 3) = Hp6(3)
    h_p(2, 2) = Hp6(4)
    h_p(2, 3) = Hp6(5)
    h_p(3, 3) = Hp6(6)

    h_p(2, 1) = h_p(1, 2)
    h_p(3, 1) = h_p(1, 3)
    h_p(3, 2) = h_p(2, 3)

    call spectral_decomposition(h_p, Q, Lambda)
    h_p = matmul(Q, matmul(diag(abs(Lambda)), transpose(Q)))
    !Lambda = 1.0_DOUBLE / sqrt(abs(Lambda))
    !h_p = matmul(Q, matmul(diag(Lambda), transpose(Q)))
  end function compute_length_diag

  ! removed 2026-09-15: compute_ellip (dead code, unreachable via schemes.txt, depended on compute_ellip)

  function corner_normal(mesh, id_sub_elem) result(norm)
    use ns_global_data_module, only: boundary_2d
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
        ( boundary_2d .and. abs(mesh%sub_face(id_sub_face)%norm(3)) < 1e-12_DOUBLE) ) then
        if( mesh%face(id_face)%left_neigh  &
          == mesh%sub_elem(id_sub_elem)%mesh_elem ) then
          norm = norm + mesh%sub_face(id_sub_face)%area*mesh%sub_face(id_sub_face)%norm
        else
          norm = norm - mesh%sub_face(id_sub_face)%area*mesh%sub_face(id_sub_face)%norm
        end if
      end if
    end do
  end function corner_normal

  subroutine compute_corr(mesh, id_vert, sol, grad, corr, h_p, second_order)
    use ns_global_data_module, only : scheme
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: h_p
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    real(kind=DOUBLE), intent(inout) :: corr
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, id_elem, id_sub_elem, le, re
    integer(kind=ENTIER) :: id_face, id_sub_face
    real(kind=DOUBLE), dimension(5) :: sol_w
    real(kind=DOUBLE) :: rho_p, a_p, div_v
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, sol_l, sol_r
    real(kind=DOUBLE), dimension(3) :: grad_p

    rho_p = 0.0_DOUBLE
    a_p = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_w = conserv_to_primit(sol(:, id_elem))
      if( second_order ) sol_w = sol_w + matmul(transpose(grad(:, :, id_elem)), &
        mesh%vert(id_vert)%coord - mesh%elem(id_elem)%coord)
      rho_p = rho_p &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(1)
      a_p = a_p &
        + mesh%sub_elem(id_sub_elem)%volume &
        * sound_speed_w(sol_w)
    end do
    rho_p = rho_p / mesh%vert(id_vert)%volume
    a_p = a_p / mesh%vert(id_vert)%volume

    div_v = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)
      div_v = div_v &
        + dot_product(sol_w_r(2:4) - sol_w_l(2:4), &
        mesh%sub_face(id_sub_face)%area&
        * mesh%sub_face(id_sub_face)%norm)
    end do
    div_v = div_v / mesh%vert(id_vert)%volume

    grad_p = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)
      if( re > 0 ) then
        grad_p = grad_p + (sol_w_r(5) - sol_w_l(5)) &
          * mesh%sub_face(id_sub_face)%area&
          * mesh%sub_face(id_sub_face)%norm
      end if
    end do
    grad_p = grad_p / mesh%vert(id_vert)%volume

    !corr = min(1.0_DOUBLE, max(-div_v/a_p, 0.0_DOUBLE))*a_p
    !corr = max(0.0_DOUBLE, min(abs(div_v)/a_p, 1.0_DOUBLE))*a_p
    corr = max(0.0_DOUBLE, &
      min(h_p(id_vert)*abs(div_v)/a_p+h_p(id_vert)*norm2(grad_p)/a_p**2, 1.0_DOUBLE))&
    *a_p
  end subroutine compute_corr

  subroutine compute_corrM(mesh, id_vert, sol, grad, corr, mat_h_p, second_order)
    use ns_global_data_module, only : scheme
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, 3, mesh%n_vert), intent(in) :: mat_h_p
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    real(kind=DOUBLE), intent(inout) :: corr
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, id_elem, id_sub_elem, le, re
    integer(kind=ENTIER) :: id_face, id_sub_face
    real(kind=DOUBLE), dimension(5) :: sol_w
    real(kind=DOUBLE) :: rho_p, a_p, div_v
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, sol_l, sol_r

    rho_p = 0.0_DOUBLE
    a_p = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_w = conserv_to_primit(sol(:, id_elem))
      if( second_order ) sol_w = sol_w + matmul(transpose(grad(:, :, id_elem)), &
        mesh%vert(id_vert)%coord - mesh%elem(id_elem)%coord)
      rho_p = rho_p &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(1)
      a_p = a_p &
        + mesh%sub_elem(id_sub_elem)%volume &
        * sound_speed_w(sol_w)
    end do
    rho_p = rho_p / mesh%vert(id_vert)%volume
    a_p = a_p / mesh%vert(id_vert)%volume

    div_v = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)
      div_v = div_v &
        + dot_product(sol_w_r(2:4) - sol_w_l(2:4), &
        mesh%sub_face(id_sub_face)%area&
        * matmul(mat_h_p(:, :, id_vert), mesh%sub_face(id_sub_face)%norm))
    end do
    div_v = div_v / mesh%vert(id_vert)%volume

    !corr = min(1.0_DOUBLE, max(-div_v/a_p, 0.0_DOUBLE))*a_p
    corr = max(0.0_DOUBLE, min(abs(div_v)/a_p, 1.0_DOUBLE))*a_p
  end subroutine compute_corrM

  subroutine compute_rhs_around_vert_ARMDWIP(mesh, sol, grad, &
      nsen, sum_lambda_vert, flux_sum_vert, &
      second_order, id_vert, vp, h_p, mat_h_p)
    use ns_global_data_module, only: bc_style, scheme, &
      exclude_bound_vert, boundary_2d
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: h_p
    real(kind=DOUBLE), dimension(3, 3, mesh%n_vert), intent(in) :: mat_h_p
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    integer(kind=ENTIER), intent(in) :: nsen
    real(kind=DOUBLE), dimension(nsen), intent(inout) :: &
      sum_lambda_vert
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: &
      flux_sum_vert
    logical, intent(in) :: second_order
    integer(kind=ENTIER), intent(in) :: id_vert

    integer(kind=ENTIER) :: j, id_sub_face, id_face
    integer(kind=ENTIER) :: id_elem, id_sub_elem
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE), dimension(3) :: norm

    real(kind=DOUBLE) :: vnl, vnr, al, ar, lambda, wpcf, corr, smax
    real(kind=DOUBLE), dimension(5) :: sol_p, sol_l, sol_r, ff
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, fminus, fplus
    real(kind=DOUBLE), dimension(5,3) :: fp, grad_sol
    real(kind=DOUBLE), dimension(3,3) :: smax_mat

    rse_loc = 0
    sol_p = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_p = sol_p + mesh%sub_elem(id_sub_elem)%volume * sol(:, id_elem)
    end do
    sol_p = sol_p / mesh%vert(id_vert)%volume

    grad_sol = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm
      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)
      grad_sol = grad_sol + tensor_product(primit_to_conserv(sol_w_r) - primit_to_conserv(sol_w_l), &
        mesh%sub_face(id_sub_face)%area*mesh%sub_face(id_sub_face)%norm)
    end do
    grad_sol = grad_sol / mesh%vert(id_vert)%volume

    call compute_corr(mesh, id_vert, sol, grad, corr, h_p, second_order)
    !call compute_corrM(mesh, id_vert, sol, grad, corr, mat_h_p, second_order)

    !call compute_nodal_velocity_LSM(mesh, id_vert, sol, grad, vp, mat_h_p, second_order)
    call compute_nodal_velocity_LS(mesh, id_vert, sol, grad, vp, h_p, second_order)

    smax_mat = diag(abs(vp(:, id_vert))) + corr*eye(3)
    fp = tensor_product(sol_p, vp(:, id_vert)) &
      - 0.5_DOUBLE*h_p(id_vert)*matmul(grad_sol, smax_mat)

    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      vnl = dot_product(sol_w_l(2:4), norm)
      vnr = dot_product(sol_w_r(2:4), norm)

      al = sound_speed_w(sol_w_l)
      ar = sound_speed_w(sol_w_r)

      sol_r = primit_to_conserv(sol_w_r)
      sol_l = primit_to_conserv(sol_w_l)

      smax = max(abs(vnr), abs(vnl)) + corr
      !smax = max(abs(vnr), abs(vnl))
      ff = 0.5_DOUBLE*(vnr*sol_r + vnl*sol_l) &
        - 0.5_DOUBLE*smax*(sol_r - sol_l)
      !ff =  0.5_DOUBLE*dot_product(vp(:, id_vert), norm)*(sol_r + sol_l) &
      !  - 0.5_DOUBLE*smax*(sol_r-sol_l)

      wpcf = 1.0_DOUBLE/3.0_DOUBLE
      !wpcf = 1.0_DOUBLE
      fminus = wpcf*matmul(fp, norm) + (1.0_DOUBLE-wpcf)*ff
      fplus = fminus

      !Used for local timestepping
      lambda = max(1e-8_DOUBLE, -vnl, vnr)+max(al,ar)

      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        sum_lambda_vert(lse_loc) = sum_lambda_vert(lse_loc) &
          + mesh%sub_face(id_sub_face)%area*lambda
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) + mesh%sub_face(id_sub_face)%area*fminus

        if (rse > 0) then
        if (mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_ARMDWIP

  subroutine compute_rhs_around_vert_LSWIP(mesh, sol, grad, &
      nsen, flux_sum_vert, &
      id_vert, vp, h_p, mat_h_p, second_order)
    use ns_global_data_module, only: bc_style, scheme, &
      exclude_bound_vert, boundary_2d
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: h_p
    real(kind=DOUBLE), dimension(3, 3, mesh%n_vert), intent(in) :: mat_h_p
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    integer(kind=ENTIER), intent(in) :: nsen
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: &
      flux_sum_vert
    integer(kind=ENTIER), intent(in) :: id_vert
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, id_sub_face, id_face
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE), dimension(3) :: norm

    real(kind=DOUBLE) :: pp
    real(kind=DOUBLE), dimension(5) :: fminus, fplus
    real(kind=DOUBLE), dimension(5,3) :: fp

    real(kind=DOUBLE) :: rho_avg, a_avg, pbar, vbar, pl, pr
    real(kind=DOUBLE) :: vnl, vnr, al, ar, lambda, wpcf
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, ff, sol_l, sol_r

    rse_loc = 0
    call compute_nodal_velocity_LS(mesh, id_vert, sol, grad, vp, h_p, second_order)
    call compute_nodal_pressure_LS(mesh, id_vert, sol, grad, pp, h_p, second_order)

    !call compute_nodal_velocity_LSM(mesh, id_vert, sol, grad, vp, mat_h_p, second_order)
    !call compute_nodal_pressure_LSM(mesh, id_vert, sol, grad, pp, mat_h_p, second_order)

    fp(1, :) = 0.0_DOUBLE
    fp(2:4, :) = pp*eye3
    fp(5, :) = pp*vp(:, id_vert)

    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      vnl = dot_product(sol_w_l(2:4), norm)
      vnr = dot_product(sol_w_r(2:4), norm)

      al = sound_speed_w(sol_w_l)
      ar = sound_speed_w(sol_w_r)

      pl = sol_w_l(5)
      pr = sol_w_r(5)

      sol_r = primit_to_conserv(sol_w_r)
      sol_l = primit_to_conserv(sol_w_l)

      rho_avg = 0.5_DOUBLE*(sol_l(1)+sol_r(1))
      a_avg = 0.5_DOUBLE*(al+ar)
      pbar = 0.5_DOUBLE*(pr+pl) &
        - 0.5_DOUBLE*rho_avg*a_avg*(vnr-vnl)
      vbar = 0.5_DOUBLE*(vnl+vnr) &
        - 0.5_DOUBLE*(pr-pl)/(rho_avg*a_avg)

      ff(1) = 0.0_DOUBLE
      ff(2:4) = pbar * norm
      ff(5) = vbar * pbar

      !wpcf = 1.0_DOUBLE/3.0_DOUBLE
      wpcf = 1.0_DOUBLE
      fminus = wpcf*matmul(fp, norm) + (1.0_DOUBLE-wpcf)*ff
      fplus = fminus

      if( boundary_2d &
        .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1e-3 ) then
        fminus = 0.0_DOUBLE
        fplus = 0.0_DOUBLE
      end if

      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) + mesh%sub_face(id_sub_face)%area*fminus

        if (rse > 0) then
        if (mesh%sub_elem(rse)%mesh_vert == id_vert) then
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_LSWIP

  subroutine compute_rhs_around_vert_LPP(mesh, sol, grad, &
      nsen, flux_sum_vert, &
      id_vert, second_order)
    use ns_global_data_module, only: bc_style, scheme, &
      exclude_bound_vert, boundary_2d
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    integer(kind=ENTIER), intent(in) :: nsen
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: &
      flux_sum_vert
    integer(kind=ENTIER), intent(in) :: id_vert
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, id_sub_face, id_face
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE), dimension(3) :: norm

    real(kind=DOUBLE) :: pp, lambda_l, lambda_r, rhol, rhor
    real(kind=DOUBLE), dimension(5) :: fminus, fplus
    real(kind=DOUBLE), dimension(5) :: fmp_l, fmp_r

    real(kind=DOUBLE) :: rho_avg, a_avg, pbar, vbar, pl, pr
    real(kind=DOUBLE) :: vnl, vnr, al, ar, wpcf
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, ff

    rse_loc = 0
    call compute_nodal_pressure_LPP(mesh, id_vert, sol, grad, pp, second_order)

    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      rhol = sol_w_l(1)
      vnl = dot_product(sol_w_l(2:4), norm)
      pl = sol_w_l(5)
      al = sound_speed_w(sol_w_l)

      rhor = sol_w_r(1)
      vnr = dot_product(sol_w_r(2:4), norm)
      pr = sol_w_r(5)
      ar = sound_speed_w(sol_w_r)

      lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(vnr - vnl))
      lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(vnr - vnl))

      fmp_l(1) = 0.0_DOUBLE
      fmp_l(2:4) = pp * norm
      fmp_l(5) = pp * (vnl - (pp - pl)/lambda_l)

      fmp_r(1) = 0.0_DOUBLE
      fmp_r(2:4) = pp * norm
      fmp_r(5) = pp * (vnr + (pp - pr)/lambda_r)

      rho_avg = 0.5_DOUBLE*(sol_w_l(1)+sol_w_r(1))
      a_avg = 0.5_DOUBLE*(al+ar)
      pbar = 0.5_DOUBLE*(pr+pl) &
        - 0.5_DOUBLE*rho_avg*a_avg*(vnr-vnl)
      vbar = 0.5_DOUBLE*(vnl+vnr) &
        - 0.5_DOUBLE*(pr-pl)/(rho_avg*a_avg)

      ff(1) = 0.0_DOUBLE
      ff(2:4) = pbar * norm
      ff(5) = vbar * pbar

      pbar = (pl/lambda_l + pr/lambda_r - (vnr-vnl))/&
        (1.0_DOUBLE/lambda_l + 1.0_DOUBLE/lambda_r)
      vbar = vnl - (pbar - pl)/lambda_l

      ff(1) = 0.0_DOUBLE
      ff(2:4) = pbar * norm
      ff(5) = vbar * pbar

      !if( (boundary_2d .and. abs(norm(3)) > 1e-8_DOUBLE )) then
      if( re <= 0 ) then
        fminus = ff
        fplus = fminus
      else
        wpcf = 1.0_DOUBLE
        fminus = wpcf*fmp_l + (1.0_DOUBLE-wpcf)*ff
        fplus = wpcf*fmp_r + (1.0_DOUBLE-wpcf)*ff
      end if

      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) + mesh%sub_face(id_sub_face)%area*fminus

        if (rse > 0) then
        if (mesh%sub_elem(rse)%mesh_vert == id_vert) then
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_LPP

  subroutine compute_nodal_pressure_LPP(mesh, id_vert, sol, grad, pp, second_order)
    use ns_global_data_module, only : scheme, boundary_2d
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    real(kind=DOUBLE), intent(inout) :: pp
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, le, re
    integer(kind=ENTIER) :: id_face, id_sub_face
    real(kind=DOUBLE) :: pbar, lambda_l, lambda_r, invlamb
    real(kind=DOUBLE) :: vnl, vnr, al, ar, rhol, rhor, pl, pr
    real(kind=DOUBLE) :: denomsum, vcorr
    real(kind=DOUBLE), dimension(3) :: norm, Bp, vl, vr
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r

    pp = 0.0_DOUBLE
    denomsum = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      norm = mesh%sub_face(id_sub_face)%norm

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)
      
      rhol = sol_w_l(1)
      vnl = dot_product(sol_w_l(2:4), norm)
      pl = sol_w_l(5)
      al = sound_speed_w(sol_w_l)

      rhor = sol_w_r(1)
      vnr = dot_product(sol_w_r(2:4), norm)
      pr = sol_w_r(5)
      ar = sound_speed_w(sol_w_r)

      lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(vnr - vnl))
      lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(vnr - vnl))

      invlamb = (1.0_DOUBLE/lambda_l + 1.0_DOUBLE/lambda_r)
      pbar = (pl/lambda_l+pr/lambda_r - (vnr - vnl))/invlamb
      if (re > 0) then
          pp = pp + mesh%sub_face(id_sub_face)%area*invlamb*pbar
          denomsum = denomsum + mesh%sub_face(id_sub_face)%area*invlamb
      else
          pp = pp + 0.5_DOUBLE*mesh%sub_face(id_sub_face)%area*invlamb*pbar
          denomsum = denomsum + 0.5_DOUBLE*mesh%sub_face(id_sub_face)%area*invlamb
      end if
    end do
    pp = pp / denomsum
  end subroutine compute_nodal_pressure_LPP

  subroutine compute_nodal_velocity_LVP(mesh, id_vert, sol, grad, vp, second_order)
    use ns_global_data_module, only : scheme, boundary_2d
    use linear_solver_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    real(kind=DOUBLE), dimension(3), intent(inout) :: vp
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, le, re
    integer(kind=ENTIER) :: id_face, id_sub_face
    real(kind=DOUBLE) :: vbar, lambda_l, lambda_r, invlamb
    real(kind=DOUBLE) :: vnl, vnr, al, ar, rhol, rhor, pl, pr
    real(kind=DOUBLE) :: denomsum, vcorr
    real(kind=DOUBLE), dimension(3) :: norm, Bp, vl, vr, rhs
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(3,3) :: mat

    mat = 0.0_DOUBLE
    rhs = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      norm = mesh%sub_face(id_sub_face)%norm

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)
      
      rhol = sol_w_l(1)
      vnl = dot_product(sol_w_l(2:4), norm)
      pl = sol_w_l(5)
      al = sound_speed_w(sol_w_l)

      rhor = sol_w_r(1)
      vnr = dot_product(sol_w_r(2:4), norm)
      pr = sol_w_r(5)
      ar = sound_speed_w(sol_w_r)

      lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(vnr - vnl))
      lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(vnr - vnl))

      mat = mat + mesh%sub_face(id_sub_face)%area&
      *(lambda_l + lambda_r)*tensor_product(norm, norm)

      vbar = (lambda_l*vnl + lambda_r*vnr - (pr -pl))/(lambda_l+lambda_r)
      rhs = rhs + mesh%sub_face(id_sub_face)%area&
      *(lambda_l + lambda_r)*vbar*norm
    end do

    call inv_lapack(3, mat)
    vp = matmul(mat, rhs)
  end subroutine compute_nodal_velocity_LVP

  subroutine compute_corr2(mesh, id_vert, sol, grad, corr, second_order)
    use ns_global_data_module, only : scheme, boundary_2d
    use linear_solver_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    real(kind=DOUBLE), intent(inout) :: corr
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, le, re, id_elem, id_sub_elem
    integer(kind=ENTIER) :: id_face, id_sub_face
    real(kind=DOUBLE) :: vbar, lambda_l, lambda_r
    real(kind=DOUBLE) :: corr_div_v, corr_grad_p, pbar
    real(kind=DOUBLE) :: vnl, vnr, al, ar, rhol, rhor, pl, pr, a_p
    real(kind=DOUBLE) :: denomsum, vcorr, div_v, invlamb, rho_p
    real(kind=DOUBLE), dimension(3) :: norm, Bp, vl, vr, rhs, grad_p
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, sol_w
    real(kind=DOUBLE), dimension(3,3) :: mat

    a_p = 0.0_DOUBLE
    rho_p = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_w = conserv_to_primit(sol(:, id_elem))
      if( second_order ) sol_w = sol_w + matmul(transpose(grad(:, :, id_elem)), &
        mesh%vert(id_vert)%coord - mesh%elem(id_elem)%coord)
      a_p = a_p &
        + mesh%sub_elem(id_sub_elem)%volume &
        * sound_speed_w(sol_w)
      rho_p = rho_p &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(1)
    end do
    rho_p = rho_p / mesh%vert(id_vert)%volume
    a_p = a_p / mesh%vert(id_vert)%volume

    mat = 0.0_DOUBLE
    rhs = 0.0_DOUBLE
    div_v = 0.0_DOUBLE
    denomsum = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      norm = mesh%sub_face(id_sub_face)%norm

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)
      
      rhol = sol_w_l(1)
      vnl = dot_product(sol_w_l(2:4), norm)
      pl = sol_w_l(5)
      al = sound_speed_w(sol_w_l)

      rhor = sol_w_r(1)
      vnr = dot_product(sol_w_r(2:4), norm)
      pr = sol_w_r(5)
      ar = sound_speed_w(sol_w_r)

      lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(vnr - vnl))
      lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(vnr - vnl))

      !Grad(p)
      mat = mat + mesh%sub_face(id_sub_face)%area&
      *(lambda_l + lambda_r)*tensor_product(norm, norm)
      vbar = (- (pr - pl))/(lambda_l+lambda_r)
      rhs = rhs + mesh%sub_face(id_sub_face)%area&
      *(lambda_l + lambda_r)*vbar*norm

      !Div(v)
      invlamb = (1.0_DOUBLE/lambda_l + 1.0_DOUBLE/lambda_r)
      pbar = ( (vnr - vnl))/invlamb
      div_v = div_v + mesh%sub_face(id_sub_face)%area*invlamb*pbar
      denomsum = denomsum + mesh%sub_face(id_sub_face)%area*invlamb
    end do

    call inv_lapack(3, mat)
    grad_p = matmul(mat, rhs) * rho_p * a_p

    div_v = div_v/(denomsum*a_p*rho_p)

    corr_div_v = min(1.0_DOUBLE, max(0.0_DOUBLE, -div_v/a_p))*a_p
    corr_grad_p = min(1.0_DOUBLE, max(0.0_DOUBLE, norm2(grad_p)/a_p))*a_p
    !corr = max(corr_div_v, corr_grad_p)
    corr = corr_div_v
  end subroutine compute_corr2

  subroutine compute_corr_ducros(mesh, id_vert, sol, grad, corr, second_order)
    use ns_global_data_module, only : boundary_2d
    use linear_solver_module, only: cross_product
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    real(kind=DOUBLE), intent(inout) :: corr
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, id_sub_elem, id_elem, id_sub_face, le, re
    real(kind=DOUBLE) :: a_p, div_v, f_ducros, ma_node
    real(kind=DOUBLE), dimension(3) :: norm, curl_v, vl, vr, v_p
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, sol_w

    ! nodal sound speed and velocity (volume-weighted)
    a_p = 0.0_DOUBLE
    v_p = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_w = conserv_to_primit(sol(:, id_elem))
      a_p = a_p + mesh%sub_elem(id_sub_elem)%volume * sound_speed_w(sol_w)
      v_p = v_p + mesh%sub_elem(id_sub_elem)%volume * sol_w(2:4)
    end do
    a_p = a_p / mesh%vert(id_vert)%volume
    v_p = v_p / mesh%vert(id_vert)%volume
    ma_node = norm2(v_p) / a_p

    ! div(v) and curl(v) via Gauss theorem on sub-faces around node
    ! using cell-centered values (no reconstruction) as is standard for shock sensors
    div_v  = 0.0_DOUBLE
    curl_v = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      if (re <= 0) cycle
      norm = mesh%sub_face(id_sub_face)%norm
      sol_w_l = conserv_to_primit(sol(:, le))
      sol_w_r = conserv_to_primit(sol(:, re))
      vl = sol_w_l(2:4)
      vr = sol_w_r(2:4)
      div_v  = div_v  + mesh%sub_face(id_sub_face)%area * dot_product(norm, vr - vl)
      curl_v = curl_v + mesh%sub_face(id_sub_face)%area * cross_product(norm, vr - vl)
    end do
    div_v  = div_v  / mesh%vert(id_vert)%volume
    curl_v = curl_v / mesh%vert(id_vert)%volume

    ! Ducros filter: ~1 at shocks (irrotational compression), ~0 in BL (shear-dominated)
    f_ducros = div_v**2 / (div_v**2 + norm2(curl_v)**2 + 1.0e-10_DOUBLE * a_p**2)

    ! activate only for compression; suppress near stagnation via Mach filter
    corr = f_ducros * min(1.0_DOUBLE, ma_node) &
         * min(1.0_DOUBLE, max(0.0_DOUBLE, -div_v / a_p)) * 4 * a_p

    if (boundary_2d .and. mesh%vert(id_vert)%is_bound) corr = 0.0_DOUBLE

  end subroutine compute_corr_ducros

  subroutine compute_corr_pressure(mesh, id_vert, sol, grad, corr, second_order)
    use ns_global_data_module, only : boundary_2d
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    real(kind=DOUBLE), intent(inout) :: corr
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, id_sub_elem, id_elem, id_sub_face, le, re
    real(kind=DOUBLE) :: a_p, rho_p, dp_avg, sum_area, pl, pr
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, sol_w

    ! nodal sound speed and density (volume-weighted)
    a_p   = 0.0_DOUBLE
    rho_p = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_w = conserv_to_primit(sol(:, id_elem))
      a_p   = a_p   + mesh%sub_elem(id_sub_elem)%volume * sound_speed_w(sol_w)
      rho_p = rho_p + mesh%sub_elem(id_sub_elem)%volume * sol_w(1)
    end do
    a_p   = a_p   / mesh%vert(id_vert)%volume
    rho_p = rho_p / mesh%vert(id_vert)%volume

    ! area-averaged pressure jump across sub-faces
    dp_avg   = 0.0_DOUBLE
    sum_area = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      if (re <= 0) cycle
      sol_w_l = conserv_to_primit(sol(:, le))
      sol_w_r = conserv_to_primit(sol(:, re))
      pl = sol_w_l(5)
      pr = sol_w_r(5)
      dp_avg   = dp_avg   + mesh%sub_face(id_sub_face)%area * abs(pr - pl)
      sum_area = sum_area + mesh%sub_face(id_sub_face)%area
    end do
    dp_avg = dp_avg / (sum_area + 1.0e-30_DOUBLE)

    ! normalize by rho*a^2 (acoustic pressure scale ~ gamma*p); clip to [0, a_p]
    corr = min(1.0_DOUBLE, dp_avg / (rho_p * a_p**2)) * 2 * a_p

    if (boundary_2d .and. mesh%vert(id_vert)%is_bound) corr = 0.0_DOUBLE

  end subroutine compute_corr_pressure

  subroutine compute_corr_pressure_div(mesh, id_vert, sol, grad, corr, second_order)
    use ns_global_data_module, only : boundary_2d
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    real(kind=DOUBLE), intent(inout) :: corr
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, id_sub_elem, id_elem, id_sub_face, le, re
    real(kind=DOUBLE) :: a_p, rho_p, p_p, dp_max, div_v, sum_area, pl, pr, h_p, ma_node
    real(kind=DOUBLE) :: corr_pressure, corr_div
    real(kind=DOUBLE), dimension(3) :: norm, v_p
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, sol_w

    ! nodal sound speed, density, pressure and velocity (volume-weighted)
    a_p   = 0.0_DOUBLE
    rho_p = 0.0_DOUBLE
    p_p   = 0.0_DOUBLE
    v_p   = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_w = conserv_to_primit(sol(:, id_elem))
      a_p   = a_p   + mesh%sub_elem(id_sub_elem)%volume * sound_speed_w(sol_w)
      rho_p = rho_p + mesh%sub_elem(id_sub_elem)%volume * sol_w(1)
      p_p   = p_p   + mesh%sub_elem(id_sub_elem)%volume * sol_w(5)
      v_p   = v_p   + mesh%sub_elem(id_sub_elem)%volume * sol_w(2:4)
    end do
    a_p   = a_p   / mesh%vert(id_vert)%volume
    rho_p = rho_p / mesh%vert(id_vert)%volume
    p_p   = p_p   / mesh%vert(id_vert)%volume
    v_p   = v_p   / mesh%vert(id_vert)%volume
    ma_node = norm2(v_p) / a_p

    ! max pressure jump and div(v) via sub-face loop
    ! dp_max (not averaged) avoids dilution when the shock crosses only a few faces
    dp_max   = 0.0_DOUBLE
    div_v    = 0.0_DOUBLE
    sum_area = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      if (re <= 0) cycle
      norm    = mesh%sub_face(id_sub_face)%norm
      sol_w_l = conserv_to_primit(sol(:, le))
      sol_w_r = conserv_to_primit(sol(:, re))
      pl = sol_w_l(5)
      pr = sol_w_r(5)
      dp_max   = max(dp_max, abs(pr - pl))
      div_v    = div_v    + mesh%sub_face(id_sub_face)%area &
                          * dot_product(norm, sol_w_r(2:4) - sol_w_l(2:4))
      sum_area = sum_area + mesh%sub_face(id_sub_face)%area
    end do
    div_v = div_v / mesh%vert(id_vert)%volume

    ! h_p ~ volume/surface as local length scale (makes div_v*h_p/a_p dimensionless)
    h_p = mesh%vert(id_vert)%volume / (sum_area + 1.0e-30_DOUBLE)

    ! Blazek-style relative pressure jump: dp/p (large at shocks, small elsewhere)
    ! Mach filter min(1,Ma) suppresses stagnation-line false positives
    !corr_pressure = min(1.0_DOUBLE, ma_node**2) * min(1.0_DOUBLE, dp_max / p_p) * a_p
    corr_pressure = min(1.0_DOUBLE, ma_node**2) * min(1.0_DOUBLE, dp_max / (rho_p * a_p**2) ) * a_p
    corr_div      = min(1.0_DOUBLE, ma_node**2) * min(1.0_DOUBLE, max(0.0_DOUBLE, -div_v * h_p / a_p)) * a_p
    corr = 4*max(corr_pressure, corr_div)

    if (boundary_2d .and. mesh%vert(id_vert)%is_bound) corr = 0.0_DOUBLE

  end subroutine compute_corr_pressure_div

  subroutine compute_rhs_around_vert_LVPPP(mesh, sol, grad, &
      nsen, flux_sum_vert, &
      id_vert, vp, h_p, second_order)
    use ns_global_data_module, only: bc_style, scheme, &
      exclude_bound_vert, boundary_2d
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: h_p
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    integer(kind=ENTIER), intent(in) :: nsen
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: &
      flux_sum_vert
    integer(kind=ENTIER), intent(in) :: id_vert
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, id_sub_face, id_face
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE), dimension(3) :: norm

    real(kind=DOUBLE) :: pp
    real(kind=DOUBLE), dimension(5) :: fminus, fplus
    real(kind=DOUBLE), dimension(5,3) :: fp

    real(kind=DOUBLE) :: pl, pr, rhol, rhor, wpcf
    real(kind=DOUBLE) :: vnl, vnr, al, ar
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, ff
    real(kind=DOUBLE) :: rhom, pm, um, up, am, vm, machm

    rse_loc = 0
    call compute_nodal_velocity_LVP(mesh, id_vert, sol, grad, vp(:,id_vert), second_order)
    call compute_nodal_pressure_LPP(mesh, id_vert, sol, grad, pp, second_order)

    !call compute_nodal_velocity_LS(mesh, id_vert, sol, grad, vp, h_p, second_order)
    !call compute_nodal_pressure_LS(mesh, id_vert, sol, grad, pp, h_p, second_order)


    fp(1, :) = 0.0_DOUBLE
    fp(2:4, :) = pp*eye3
    fp(5, :) = pp*vp(:, id_vert)

    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      rhol = sol_w_l(1)
      vnl = dot_product(sol_w_l(2:4), norm)
      pl = sol_w_l(5)
      al = sound_speed_w(sol_w_l)

      rhor = sol_w_r(1)
      vnr = dot_product(sol_w_r(2:4), norm)
      pr = sol_w_r(5)
      ar = sound_speed_w(sol_w_r)

      rhom = 0.5_DOUBLE*(rhor+rhol)
      pm = 0.5_DOUBLE*(pl+pr)
      um = 0.5_DOUBLE*(vnl+vnr)
      am = max(al, ar)

      vm = 0.5_DOUBLE*(norm2(sol_w_l(2:4)) + norm2(sol_w_r(2:4)))
      machm = vm/am

      up = um - 0.5_DOUBLE/(rhom*am)*(pr-pl)
      pp = pm - 0.5_DOUBLE*rhom*am*(vnr-vnl)

      ff(1) = 0.0_DOUBLE
      ff(2:4) = pp * norm
      ff(5) = pp * up

      !wpcf = 1.0_DOUBLE/3.0_DOUBLE
      wpcf = 1.0_DOUBLE

      fminus = wpcf*matmul(fp, norm)+(1.0_DOUBLE-wpcf)*ff
      fplus = fminus

      if( boundary_2d &
        .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1e-3 ) then
        fminus = 0.0_DOUBLE
        fplus = 0.0_DOUBLE
      end if

      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) + mesh%sub_face(id_sub_face)%area*fminus

        if (rse > 0) then
        if (mesh%sub_elem(rse)%mesh_vert == id_vert) then
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_LVPPP

  subroutine compute_rhs_around_vert_LS1D(mesh, sol, grad, &
      nsen, flux_sum_vert, &
      id_vert, second_order)
    use ns_global_data_module, only: bc_style, scheme, &
      exclude_bound_vert, boundary_2d
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    integer(kind=ENTIER), intent(in) :: nsen
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: &
      flux_sum_vert
    integer(kind=ENTIER), intent(in) :: id_vert
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, id_sub_face, id_face
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE), dimension(3) :: norm

    real(kind=DOUBLE) :: pp, rhol, rhor
    real(kind=DOUBLE), dimension(5) :: fminus, fplus

    real(kind=DOUBLE) :: pl, pr, vm, machm
    real(kind=DOUBLE) :: vnl, vnr, al, ar
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, ff

    real(kind=DOUBLE) :: rhom, pm, um, up, am

    rse_loc = 0
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      rhol = sol_w_l(1)
      vnl = dot_product(sol_w_l(2:4), norm)
      pl = sol_w_l(5)
      al = sound_speed_w(sol_w_l)

      rhor = sol_w_r(1)
      vnr = dot_product(sol_w_r(2:4), norm)
      pr = sol_w_r(5)
      ar = sound_speed_w(sol_w_r)

      rhom = 0.5_DOUBLE*(rhor+rhol)
      pm = 0.5_DOUBLE*(pl+pr)
      um = 0.5_DOUBLE*(vnl+vnr)
      am = max(al, ar)

      vm = 0.5_DOUBLE*(norm2(sol_w_l(2:4)) + norm2(sol_w_r(2:4)))
      machm = vm/am

      up = um - 0.5_DOUBLE/(rhom*am)*(pr-pl)
      pp = pm - 0.5_DOUBLE*rhom*am*(vnr-vnl)

      ff(1) = 0.0_DOUBLE
      ff(2:4) = pp * norm
      ff(5) = pp * up

      fminus = ff
      fplus = fminus

      if( boundary_2d &
        .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1e-3 ) then
        fminus = 0.0_DOUBLE
        fplus = 0.0_DOUBLE
      end if

      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) &
          + mesh%sub_face(id_sub_face)%area*fminus

        if (rse > 0) then
        if (mesh%sub_elem(rse)%mesh_vert == id_vert) then
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) &
            - mesh%sub_face(id_sub_face)%area*fplus
        end if
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_LS1D

  subroutine build_grad_nodal_system(mesh, id_vert, N, S, S_tilde_T, B)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_faces_neigh, mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(inout) :: N
    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_faces_neigh, mesh%vert(id_vert)%n_sub_elems_neigh), &
      intent(inout) :: S
    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_elems_neigh, mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(inout) :: S_tilde_T
    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_faces_neigh), intent(inout) :: B

    integer(kind=ENTIER) :: l, p
    integer(kind=ENTIER) :: id_sub_elem, id_elem, id_sub_elem_loc
    integer(kind=ENTIER) :: nsfn, nsen
    integer(kind=ENTIER) :: isfl, isfp, isfl_loc, isfp_loc, ifl, ifp
    real(kind=DOUBLE) :: blc_1
    real(kind=DOUBLE), dimension(3) :: norml, normp

    nsfn = mesh%vert(id_vert)%n_sub_faces_neigh
    nsen = mesh%vert(id_vert)%n_sub_elems_neigh

    N = 0.0_DOUBLE
    S = 0.0_DOUBLE
    S_tilde_T = 0.0_DOUBLE
    B = 0.0_DOUBLE

    !Build NT = ST + B
    do l = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      isfl = mesh%vert(id_vert)%sub_face_neigh(l)
      isfl_loc = mesh%sub_face(isfl)%id_loc_around_node
      ifl = mesh%sub_face(isfl)%mesh_face

      !Left sub elem
      id_sub_elem = mesh%sub_face(isfl)%left_sub_elem_neigh
      id_elem = mesh%face(ifl)%left_neigh
      id_sub_elem_loc = mesh%sub_elem(id_sub_elem)%id_loc_around_node
      norml = mesh%sub_face(isfl)%norm

      !Add terms in N, S, and S_tilde
      do p = 1, mesh%sub_elem(id_sub_elem)%n_sub_faces
        isfp = mesh%sub_elem(id_sub_elem)%sub_face(p)
        isfp_loc = mesh%sub_face(isfp)%id_loc_around_node
        ifp = mesh%sub_face(isfp)%mesh_face
        if (mesh%sub_face(isfp)%left_sub_elem_neigh == id_sub_elem) then
          normp = mesh%sub_face(isfp)%norm
        else
          normp = -mesh%sub_face(isfp)%norm
        end if
        blc_1 = mesh%sub_face(isfl)%area*mesh%sub_face(isfp)%area &
          *(1.0_DOUBLE/mesh%sub_elem(id_sub_elem)%volume) &
          *dot_product(normp, norml)
        N(isfl_loc, isfp_loc) = N(isfl_loc, isfp_loc) + blc_1
        S(isfl_loc, id_sub_elem_loc) = S(isfl_loc, id_sub_elem_loc) + blc_1
        S_tilde_T(id_sub_elem_loc, isfl_loc) = &
          S_tilde_T(id_sub_elem_loc, isfl_loc) + blc_1
      end do

      !Right sub elem
      norml = -mesh%sub_face(isfl)%norm
      id_sub_elem = mesh%sub_face(isfl)%right_sub_elem_neigh
      id_elem = mesh%face(ifl)%right_neigh
      if (id_sub_elem > 0) then
        id_sub_elem_loc = mesh%sub_elem(id_sub_elem)%id_loc_around_node

        !Add terms in N, S, and S_tilde
        do p = 1, mesh%sub_elem(id_sub_elem)%n_sub_faces
          isfp = mesh%sub_elem(id_sub_elem)%sub_face(p)
          isfp_loc = mesh%sub_face(isfp)%id_loc_around_node
          ifp = mesh%sub_face(isfp)%mesh_face
          if (mesh%sub_face(isfp)%left_sub_elem_neigh == id_sub_elem) then
            normp = mesh%sub_face(isfp)%norm
          else
            normp = -mesh%sub_face(isfp)%norm
          end if
          blc_1 = mesh%sub_face(isfl)%area*mesh%sub_face(isfp)%area &
            *(1.0_DOUBLE/mesh%sub_elem(id_sub_elem)%volume) &
            *dot_product(normp, norml)
          N(isfl_loc, isfp_loc) = N(isfl_loc, isfp_loc) + blc_1
          S(isfl_loc, id_sub_elem_loc) = S(isfl_loc, id_sub_elem_loc) + blc_1
          S_tilde_T(id_sub_elem_loc, isfl_loc) = &
            S_tilde_T(id_sub_elem_loc, isfl_loc) + blc_1
        end do
      end if
    end do
  end subroutine build_grad_nodal_system

  subroutine compute_rhs_around_vert_LPF(mesh, sol, grad, &
      nsen, flux_sum_vert, &
      id_vert, second_order)
    use ns_global_data_module, only: bc_style, scheme, &
      exclude_bound_vert, boundary_2d, gamma
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    integer(kind=ENTIER), intent(in) :: nsen
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: &
      flux_sum_vert
    integer(kind=ENTIER), intent(in) :: id_vert
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, id_sub_face, id_face
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE), dimension(3) :: norm

    real(kind=DOUBLE) :: lambda_l, lambda_r, rhol, rhor
    real(kind=DOUBLE), dimension(5) :: fminus, fplus
    real(kind=DOUBLE), dimension(5) :: fmp_l, fmp_r
    real(kind=DOUBLE), dimension(3) :: vp

    real(kind=DOUBLE) :: rho_avg, a_avg, pbar, vbar, pl, pr
    real(kind=DOUBLE) :: vnl, vnr, al, ar, wpcf, pp, vpcf
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, ff, sol_w

    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_elems_neigh) :: p_elem
    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_faces_neigh) :: p_face

    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_faces_neigh, &
      mesh%vert(id_vert)%n_sub_faces_neigh) :: N
    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_faces_neigh, &
      mesh%vert(id_vert)%n_sub_elems_neigh) :: S
    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_elems_neigh, &
      mesh%vert(id_vert)%n_sub_faces_neigh) :: S_tilde_T
    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_faces_neigh) :: B

    integer(kind=ENTIER) :: id_elem, id_sub_elem, k, id_sub_face_loc
    real(kind=DOUBLE) :: aelem
    real(kind=DOUBLE), dimension(3, mesh%vert(id_vert)%n_sub_elems_neigh) :: grad_p_tilde
    real(kind=DOUBLE), dimension(3,3, mesh%vert(id_vert)%n_sub_elems_neigh) :: mat

    rse_loc = 0
    call compute_nodal_pressure_LPP(mesh, id_vert, sol, grad, pp, second_order)

    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      if (mesh%sub_elem(lse)%mesh_vert /= id_vert) cycle
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      rhol = sol_w_l(1)
      vnl = dot_product(sol_w_l(2:4), norm)
      pl = sol_w_l(5)
      al = sound_speed_w(sol_w_l)

      rhor = sol_w_r(1)
      vnr = dot_product(sol_w_r(2:4), norm)
      pr = sol_w_r(5)
      ar = sound_speed_w(sol_w_r)

      lambda_l = rhol*al
      lambda_r = rhor*ar
      vpcf = (lambda_l*vnl + lambda_r*vnr - (pr-pl))/(lambda_l+lambda_r)

      ff(1) = 0.0_DOUBLE
      ff(2:4) = pp * norm
      ff(5) = vpcf * pp

      fminus = ff
      fplus = fminus

      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) + mesh%sub_face(id_sub_face)%area*fminus

        if (rse > 0) then
        if (mesh%sub_elem(rse)%mesh_vert == id_vert) then
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_LPF

  subroutine compute_rhs_around_vert_WIP(mesh, sol, grad, &
      nsen, flux_sum_vert, sum_lambda_vert, &
      id_vert, vp, h_p, second_order)
    use ns_global_data_module, only: bc_style, scheme, &
      exclude_bound_vert, boundary_2d, wip_diss, grad_filled
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: h_p
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    integer(kind=ENTIER), intent(in) :: nsen
    real(kind=DOUBLE), dimension(nsen), intent(inout) :: &
      sum_lambda_vert
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: &
      flux_sum_vert
    integer(kind=ENTIER), intent(in) :: id_vert
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, id_sub_face, id_face
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE), dimension(3) :: norm

    ! wip_diss : forme de la dissipation de pression (voir ns_global_data_module)
    real(kind=DOUBLE) :: zface, dlen, divf, divl, divr
    real(kind=DOUBLE) :: afa, dlen_l, dlen_r

    real(kind=DOUBLE) :: pp
    real(kind=DOUBLE), dimension(5) :: fminus, fplus
    real(kind=DOUBLE), dimension(5,3) :: fmp_lag, fmp_adv

    integer(kind=ENTIER) :: id_sub_elem, id_elem, k
    real(kind=DOUBLE) :: pl, pr, rhol, rhor, wpcf, vbar
    real(kind=DOUBLE) :: vnl, vnr, al, ar, lambda_lts, lambda_adv
    real(kind=DOUBLE) :: lambda_l, lambda_r, sum_area
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, ff_lag, ff_adv, sol_p
    real(kind=DOUBLE), dimension(5) :: sol_l, sol_r, sol_m
    real(kind=DOUBLE), dimension(5,3) :: grad_sol_p
    real(kind=DOUBLE), dimension(3) :: vpm
    real(kind=DOUBLE) :: rhom, pm, um, up, am, vm, machm, pbar, corr
    real(kind=DOUBLE) :: divv_p, sum_area_lambda, pbar2

    !MULTI POINT LAG
    rse_loc = 0
    call compute_nodal_velocity_LVP(mesh, id_vert, sol, grad, vp(:,id_vert), second_order)
    call compute_nodal_pressure_LPP(mesh, id_vert, sol, grad, pp, second_order)

    fmp_lag(1, :) = 0.0_DOUBLE
    fmp_lag(2:4, :) = pp*eye3
    fmp_lag(5, :) = pp*vp(:, id_vert)

    !MULTI POINT ADV
    grad_sol_p = 0.0_DOUBLE
    sum_area = 0.0_DOUBLE
    sum_area_lambda = 0.0_DOUBLE
    divv_p = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm
      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)
      sol_l = primit_to_conserv(sol_w_l)
      rhol = sol_w_l(1)
      vnl = dot_product(sol_w_l(2:4), norm)
      pl = sol_w_l(5)
      al = sound_speed_w(sol_w_l)
      sol_r = primit_to_conserv(sol_w_r)
      rhor = sol_w_r(1)
      vnr = dot_product(sol_w_r(2:4), norm)
      pr = sol_w_r(5)
      ar = sound_speed_w(sol_w_r)

      lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(vnr - vnl))
      lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(vnr - vnl))

      grad_sol_p = grad_sol_p + mesh%sub_face(id_sub_face)%area &
        *tensor_product(sol_r - sol_l, mesh%sub_face(id_sub_face)%norm)
      sum_area = sum_area + mesh%sub_face(id_sub_face)%area

      if( re > 0 ) then
        divv_p = divv_p &
          + mesh%sub_face(id_sub_face)%area*(vnr - vnl)
        sum_area_lambda = sum_area_lambda &
          + mesh%sub_face(id_sub_face)%area &
          *(1.0_DOUBLE/lambda_l + 1.0_DOUBLE/lambda_r)
      else
        divv_p = divv_p &
          + 0.5_DOUBLE*mesh%sub_face(id_sub_face)%area*(vnr - vnl)
        sum_area_lambda = sum_area_lambda &
          + 0.5_DOUBLE*mesh%sub_face(id_sub_face)%area&
          *(1.0_DOUBLE/lambda_l + 1.0_DOUBLE/lambda_r)
      end if
    end do
    grad_sol_p = grad_sol_p / sum_area
    divv_p = divv_p / sum_area_lambda

    sol_p = 0.0_DOUBLE
    vpm = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_p = sol_p + mesh%sub_elem(id_sub_elem)%volume * sol(:, id_elem)
      vpm = vpm + mesh%sub_elem(id_sub_elem)%volume * sol(2:4, id_elem)/sol(1, id_elem)
    end do
    sol_p = sol_p / mesh%vert(id_vert)%volume
    vpm = vpm / mesh%vert(id_vert)%volume

    !fmp_adv = tensor_product(sol_p, vp(:, id_vert)) &
    !  - 0.5_DOUBLE*diag(abs(vp(:, id_vert)))*grad_sol_p
    fmp_adv = tensor_product(sol_p, vp(:, id_vert)) &
      - 0.5_DOUBLE*norm2(vp(:, id_vert))*grad_sol_p

    !call compute_corr2(mesh, id_vert, sol, grad, corr, second_order)
    call compute_corr_ducros(mesh, id_vert, sol, grad, corr, second_order)
    !call compute_corr_pressure(mesh, id_vert, sol, grad, corr, second_order)
    !call compute_corr_pressure_div(mesh, id_vert, sol, grad, corr, second_order)

    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      sol_l = primit_to_conserv(sol_w_l)
      rhol = sol_w_l(1)
      vnl = dot_product(sol_w_l(2:4), norm)
      pl = sol_w_l(5)
      al = sound_speed_w(sol_w_l)

      sol_r = primit_to_conserv(sol_w_r)
      rhor = sol_w_r(1)
      vnr = dot_product(sol_w_r(2:4), norm)
      pr = sol_w_r(5)
      ar = sound_speed_w(sol_w_r)

      rhom = 0.5_DOUBLE*(rhol+rhor)
      am = 0.5_DOUBLE*(al+ar)
      !am = max(al, ar)
      vm = 0.5_DOUBLE*(norm2(sol_w_l(2:4)) + norm2(sol_w_r(2:4)))
      machm = vm/am
      sol_m = 0.5_DOUBLE*(sol_r + sol_l)

      !lambda_l = rhol*al
      !lambda_r = rhor*ar
      lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(vnr - vnl))
      lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(vnr - vnl))
      vbar = (lambda_l*vnl + lambda_r*vnr - (pr-pl))/(lambda_l+lambda_r)
      !pbar = (lambda_r * pl + lambda_l * pr - lambda_l*lambda_r*(vnr-vnl))/(lambda_l+lambda_r)
      !pbar = (lambda_r * pl + lambda_l * pr - lambda_l*lambda_r*(vnr-vnl))/(lambda_l+lambda_r)
      ! --- dissipation de pression : les trois cas sont exclusifs ---
      select case (wip_diss)
      case (1)
        ! saut normal par face (forme three_wave)
        zface = lambda_l*lambda_r/(lambda_l+lambda_r)
        pbar = (lambda_r*pl+lambda_l*pr)/(lambda_l+lambda_r) &
          - zface*(vnr - vnl)
      case (2)
        ! vraie divergence par face. Les gradients de maille valent d_i w_k,
        ! donc div v = grad(1,2) + grad(2,3) + grad(3,4). Le terme en z est
        ! exclu sur un maillage extrude (boundary_2d) : la divergence physique
        ! y est celle du plan.
        !
        ! LONGUEUR. Surtout PAS |x_r - x_l| : stitch_periodic_faces recolle une
        ! face periodique en posant face%right_neigh = la maille de l'AUTRE BOUT
        ! du domaine, sans translation de coordonnees. La distance des deux
        ! centroides y vaut donc la largeur du domaine (facteur ~49 sur le
        ! Gresho 50x50), et les mailles de coin en portent deux. On prend a la
        ! place l'extension normale INTRINSEQUE de chaque maille a cette face,
        ! d_c = V_c / A_face, qui ne fait intervenir aucune difference de
        ! coordonnees : exacte pour un cube (d = h), symetrique sous
        ! retournement de face, et insensible a la connectivite.
        zface = lambda_l*lambda_r/(lambda_l+lambda_r)
        ! id_face n'est pas affecte ailleurs dans cette routine : le faire ici.
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        afa = mesh%face(id_face)%area
        if( .not. (second_order .or. grad_filled) ) then
          ! gradients absents (passe d'ordre 1 sans grad_filled) : on retombe
          ! sur le saut. Pas compte : ce n'est pas un fait de maillage.
          pbar = (lambda_r*pl+lambda_l*pr)/(lambda_l+lambda_r) &
            - zface*(vnr - vnl)
        else
          divl = grad(1, 2, le) + grad(2, 3, le)
          if( .not. boundary_2d ) divl = divl + grad(3, 4, le)
          dlen_l = mesh%elem(le)%volume/max(afa, 1e-300_DOUBLE)
          if( re > 0 ) then
            divr = grad(1, 2, re) + grad(2, 3, re)
            if( .not. boundary_2d ) divr = divr + grad(3, 4, re)
            dlen_r = mesh%elem(re)%volume/max(afa, 1e-300_DOUBLE)
            divf = 0.5_DOUBLE*(divl + divr)
            dlen = 0.5_DOUBLE*(dlen_l + dlen_r)
          else
            ! BORD REEL (pas de maille droite). On n'utilise PAS le saut ici :
            ! melanger une dissipation en saut au bord et une dissipation en
            ! divergence a l'interieur laisse une couronne de forcage parasite
            ! tout autour du domaine, qui ne s'annule pas sur un champ
            ! solenoidal alors que l'interieur s'annule. On prend la divergence
            ! de la seule maille gauche, avec sa propre extension : la propriete
            ! << nul sur un champ a divergence nulle >> reste alors UNIFORME
            ! jusqu'au bord.
            divf = divl
            dlen = dlen_l
          end if
          pbar = (lambda_r*pl+lambda_l*pr)/(lambda_l+lambda_r) &
            - zface*dlen*divf
        end if
      case default
        ! 0 : dissipation nodale, comportement existant
        pbar = (lambda_r*pl+lambda_l*pr)/(lambda_l+lambda_r) - divv_p
      end select
      !if( abs(pbar - pbar2)/pbar > 0.1_DOUBLE) then
      !  print*, pbar, pbar2
      !end if

      !vbar = 0.5_DOUBLE*(vnl+vnr) - 0.5_DOUBLE/(rhom*am) * (pr - pl)
      !pbar = 0.5_DOUBLE*(pl+pr) - 0.5_DOUBLE*rhom*am*(vnr-vnl)
      !ff_adv = vbar*sol_m - 0.5_DOUBLE*max(abs(vnl),abs(vnr))*(sol_r-sol_l)
      !ff_adv = 0.5_DOUBLE*(vnl*sol_l+vnr*sol_r) &
      !  - 0.5_DOUBLE*max(abs(vnl),abs(vnr))*(sol_r-sol_l)

      !ff_adv = vbar*sol_m - 0.5_DOUBLE*(abs(vbar)+min(1.0_DOUBLE, machm)*am)*(sol_r-sol_l)
      !print*, corr

      !Not ok start carbuncle
      ff_adv = vbar*sol_m - 0.5_DOUBLE*(abs(vbar)+corr)*(sol_r-sol_l)

      !ff_adv = vbar*sol_m - 0.5_DOUBLE*(abs(vbar)+corr)*(sol_r-sol_l)
      !ff_adv = vbar*sol_m - 0.5_DOUBLE*abs(vbar)*(sol_r-sol_l)

      !OK ish no carbuncle
      !ff_adv = 0.5_DOUBLE*(vnl*sol_l + vnr*sol_r) &
        !- 0.5_DOUBLE*(max(abs(vnl),abs(vnr))+0.01_DOUBLE*corr)*(sol_r-sol_l)

      ff_lag(1) = 0.0_DOUBLE
      ff_lag(2:4) = pbar * norm
      ff_lag(5) = pbar * vbar
      !ff_lag(2:4) = pp * norm
      !ff_lag(5) = pp * vbar

      wpcf = 1.0_DOUBLE/3.0_DOUBLE
      !wpcf = 0.0_DOUBLE
      !fminus = wpcf*matmul(fmp_adv + fmp_lag, norm) &
      !  + (1.0_DOUBLE-wpcf)*(ff_adv+ff_lag)
      !fminus = wpcf*matmul(fmp_adv, norm) + (1.0_DOUBLE-wpcf)*ff_adv + ff_lag
      !fminus = wpcf*matmul(fmp_adv, norm) + (1.0_DOUBLE-wpcf)*ff_adv + matmul(fmp_lag, norm)
      !fminus = wpcf*matmul(fmp_lag, norm) + (1.0_DOUBLE-wpcf)*ff_lag + ff_adv
      !fminus = wpcf*matmul(fmp_lag, norm) + (1.0_DOUBLE-wpcf)*ff_lag + ff_adv
      !fminus = matmul(fmp_lag, norm) + ff_adv
      !fminus = matmul(fmp_adv, norm) + ff_lag
      !fminus = matmul(fmp_adv + fmp_lag, norm)
      !fminus = wpcf*matmul(fmp_lag, norm) + (1.0_DOUBLE-wpcf)*ff_lag + matmul(fmp_adv, norm)
      !fminus = wpcf*matmul(fmp_adv, norm) + (1.0_DOUBLE-wpcf)*ff_adv + matmul(fmp_lag, norm)
      fminus = ff_adv + ff_lag
      fplus = fminus

      if( boundary_2d &
        .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1e-3 ) then
        fminus = 0.0_DOUBLE
        fplus = 0.0_DOUBLE
      end if

      lambda_lts = max(abs(vnl), abs(vnr)) + max(al, ar)

      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        sum_lambda_vert(lse_loc) = sum_lambda_vert(lse_loc) &
          + mesh%sub_face(id_sub_face)%area*lambda_lts
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) + mesh%sub_face(id_sub_face)%area*fminus

        if (rse > 0) then
        if (mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda_lts
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_WIP

  ! WIP2: nodal pressure built from the per-sub-face acoustic solver (as in LPF,
  ! which measured ~450x lower low-Mach floor than the one-sided LPP form because
  ! its energy-flux velocity is symmetric, hence conservative), with the velocity
  ! jump scaled by theta = min(1, Ma_node).
  !
  ! The jump term -theta*(vnr-vnl) is the "divv" contribution: it is what couples
  ! the multi-D velocity field into the nodal pressure (and what keeps the Gresho
  ! vortex alive on quads), but unscaled it is O(Ma) relative to p, one power short
  ! of the incompressible limit -- the measured consequence is the L2_rho floor.
  ! Scaling it by Ma restores O(Ma^2) while leaving theta=1 (i.e. the unmodified
  ! acoustic solver) everywhere the flow is transonic or faster.
  ! The pressure jump inside vstar is deliberately NOT scaled: it carries the
  ! pressure/velocity coupling that the low-Mach limit needs.
  subroutine compute_rhs_around_vert_WIP2(mesh, sol, grad, &
      nsen, flux_sum_vert, sum_lambda_vert, &
      id_vert, second_order, low_mach, adv_mode, eps_mode, enth_fix, tp_fix, ep_mode, lag_mode)
    use linear_solver_module
    use ns_global_data_module, only: boundary_2d, zb_eps_coef, zb_pivp, zb_ppus, zb_tang, zb_lag, &
      grad_filled, zb_pdiv
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    integer(kind=ENTIER), intent(in) :: nsen
    real(kind=DOUBLE), dimension(nsen), intent(inout) :: sum_lambda_vert
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: flux_sum_vert
    integer(kind=ENTIER), intent(in) :: id_vert
    logical, intent(in) :: second_order
    logical, intent(in) :: low_mach
    ! 0 -> WIP's Rusanov-at-vstar advection with the Ducros sensor
    !      (best low-Mach slope, least accurate on Sedov)
    ! 1 -> AMISO advection: nodal upwind state blended with a pressure-corrected
    !      central flux (best on Sedov, slope drops below 1)
    ! 2 -> shock-sensor blend of the two: AMISO where the Ducros indicator says
    !      "shock", Rusanov-at-vstar where the flow is smooth/vortical
    ! 3 -> ARMD-style MULTI-D nodal advection hybridised with the 1D face flux:
    !        F = w * (sol_p (x) v_p - h_p/2 * grad_p(U).SMAX).n + (1-w) * F_1D
    !      with SMAX = diag(|v_p,k|). The ZB ARMD* schemes hardwire w = 0.5 and
    !      are catastrophic here (1e-2 to 3e-1 on the wall heat flux, plus Sedov
    !      crashes); AMISO, the one robust member of that family, instead uses
    !      w = (min_f A_f / A_f)/3, i.e. at most a third and weighted by the
    !      local sub-face area ratio. Mode 3 takes the multi-d flux with AMISO's
    !      hybridisation. h_p is taken as V_p / sum_f A_f, the volume-to-surface
    !      length, so no external h_p array is needed.
    ! 4 -> same, but w = 0.5 as in ARMD, to separate "multi-d is wrong" from
    !      "multi-d was hybridised too strongly".
    ! 5 -> mode 3 with an enthalpy-consistent gradient. Measured: mode 3 still
    !      floors at low Mach (1.26e-6 at Ma=1e-4 against 6.03e-8 without it).
    !      The reason is the energy row of grad_p(U): rho*E ~ p/(gamma-1) is
    !      O(1/Ma^2), so dissipating its gradient injects an error that dwarfs
    !      the O(Ma^2) density fluctuation being measured. Since rho*E = rho*h - p,
    !      replacing d(rho*E) by h_bar*d(rho) in the dissipated gradient is the
    !      same Haenel condition already used by enth_fix on the 1D jump, and it
    !      removes exactly that term.
    integer(kind=ENTIER), intent(in) :: adv_mode
    ! Which eps_p carbuncle sensor feeds the advection viscosity. These are the
    ! four candidates of tex/wip.tex; Vincent's note there is that they give very
    ! different wall heat fluxes, because the term must fire inside the shock and
    ! NOT in the boundary layer or at the stagnation point.
    ! 0 -> compute_corr_ducros       (c) Ducros filter          [4*a_p]
    ! 1 -> compute_corr_pressure     (d) normalised p jump      [2*a_p]
    ! 2 -> compute_corr2             (a) -div(v) only           [1*a_p]
    ! 3 -> compute_corr_pressure_div (b) max(p jump, -div v)    [4*a_p]
    ! 4 -> no sensor at all: a PER-FACE, positivity-driven eps built the same way
    !      the Lagrange slopes already are, lambda = max(rho*a, sqrt(rho*dp),
    !      -rho*dv). Transposed to a velocity scale that gives
    !          eps_f = max(0, -(v_r-v_l).n, sqrt(|p_r-p_l| / rho_m))
    !      which needs no threshold, no Mach gate and no shock detector, and
    !      vanishes in a boundary layer by construction: dp/dn ~ 0 across a
    !      boundary layer and the velocity variation there is tangential, so both
    !      terms are ~0, while a shock has strong compression and a large dp.
    !      Note this eps is per-face, unlike modes 0-3 whose sensor is nodal and
    !      is applied identically to every sub-face around the node.
    !   8 -> AUCUN capteur : vitesse seule plus eps = zb_eps_coef*a_p (namelist).
    integer(kind=ENTIER), intent(in) :: eps_mode
    ! Haenel / MGallice enthalpy preservation (Tallois, papers/talois.pdf, eq 22-23):
    ! the energy dissipation of the flux must equal the total enthalpy times the
    ! mass dissipation, else total enthalpy is not preserved across the bow shock.
    ! Since rho*E = rho*h - p, replacing d(rho*E) by h_bar*d(rho) in the dissipation
    ! is exactly that condition. Tallois shows the unmodified multidimensional
    ! (Gallice-2D) solver puts the density maximum off the stagnation point, which
    ! is what ruins a wall heat flux.
    logical, intent(in) :: enth_fix
    ! Tallois PhD section 6.3.3, eq (6.3.3.2): the MULTIDIMENSIONAL low-Mach
    ! correction. p_l^theta = theta_n*[p_l - lambda_l*(u_n - u_l).n] + (1-theta_n)*q_n
    ! with theta_n = min(1, |u_n|/a). Here q_n (eq 6.3.3.1) is exactly this code's
    ! nodal pressure pp, so the low-Mach end (theta_n->0) reproduces WIP2_NOLM
    ! bit-for-bit, and the shock end (theta_n->1) becomes the one-sided EUCCLHYD
    ! nodal pressure flux built on the nodal velocity. Crucially BOTH ends are
    ! nodal, so unlike a convex blend with the 1D flux the multidimensional
    ! character is never lost -- that is the thesis's stated reason for this form.
    logical, intent(in) :: tp_fix
    ! Eulerian port of the Lagrangian 'vp_pp_ep' scheme (src/lagrange, pi_pp_ep):
    ! momentum keeps the nodal pressure pp, but the Lagrange ENERGY flux of every
    ! interior sub-face becomes the one-sided flux built on a NODAL total energy
    !     (pv)*_l = p_l v_l - lambda_l (e_p - e_l)
    !     (pv)*_r = p_r v_r + lambda_r (e_p - e_r)
    ! with e = rho*E/rho the specific total energy and mu = lambda (phi = 1, no free
    ! coefficient). e_p is DEFINED by the telescoping sum_f A[(pv)*_r - (pv)*_l] = 0,
    !     e_p = sum A[lam_l e_l + lam_r e_r + (p_l v_l - p_r v_r)] / sum A (lam_l + lam_r)
    ! so the node conserves energy exactly; it involves no momentum flux, which is
    ! why swapping the one-sided Pi for pp costs nothing in energy conservation.
    !
    ! Departure from the Lagrangian code, deliberate: there the boundary sub-faces
    ! enter e_p through a mirror state, and their flux lambda (e_p - e) carries no
    ! velocity factor, so on a wall it is an INVISCID energy flux through an
    ! impermeable surface -- on the hypersonic cylinder that is a spurious wall heat
    ! flux, exactly the quantity being measured. Here boundary sub-faces (re <= 0)
    ! and the boundary_2d z-faces (zeroed below anyway) are excluded from e_p and
    ! keep the WIP2 energy flux pp*vstar, which vanishes at a wall. Telescoping over
    ! the interior sub-faces of the node is unaffected, hence still exact.
    !
    ! ep_mode: 0 -> off (plain WIP2 energy flux)
    !          1 -> WIP2_PPEP,  X = E   the specific TOTAL energy (the Lagrangian vp_pp_ep)
    !          2 -> WIP2_PPEPI, X = eps = E - |u|^2/2, the specific INTERNAL energy
    ! The telescoping only needs the form +-lambda (X_p - X_s) with one X_p shared by both
    ! sides, so total energy is conserved exactly for ANY per-cell scalar X. Mode 2 exists
    ! because mode 1 overheats the M17.6 cylinder wall 6x (wall cell 1027 K vs 626 K for
    ! WIP2_NOLM at equal Cp): across a boundary layer the jump in E is dominated by the
    ! kinetic part (~87x the internal one at M17.6), and relaxing it turns that shear into
    ! heat. Mode 2 removes exactly that part, hence tests the mechanism.
    !          3 -> WIP2_PPEPV, X = E,   mu_s = rho_s |(v_s - v*).n|
    !          4 -> WIP2_PPEPVI, X = eps, mu_s = rho_s |(v_s - v*).n|
    ! Modes 1/2 dissipate X with mu = lambda = rho*a, i.e. at ACOUSTIC speed. A thermal
    ! boundary layer is an entropy (contact) wave, which Godunov dissipates with the wave's
    ! own speed; mode 1/2 instead act as an artificial heat conduction rho*a*c_v
    ! (~3.6e4 W/m2/K on the M17.6 cylinder, ~7x the physical first-cell conductance) and
    ! overheat the wall cell. Modes 3/4 take the impedance of the cell's normal velocity
    ! RELATIVE to the interface, rho_s |(v_s - v*).n|: derived, no min/switch/coefficient,
    ! it vanishes at a wall (v* = 0, v.n -> 0) and across an exact Lagrangian contact.
    ! Degenerate node: sum A(mu_l+mu_r) = 0 means v_l.n = v_r.n = v* on every interior
    ! sub-face, which forces p_l = p_r there, so the one-sided p v fluxes already telescope;
    ! e_p is then irrelevant and the node falls back to the plain WIP2 energy flux, which
    ! is also conservative. The guard is "den > 0", no threshold. Note the mu -> 0 LIMIT is
    ! not that fallback: mu_l e_p -> (mu_l/sum mu) * sum A(p_l v_l - p_r v_r), a
    ! redistribution of the nodal imbalance of the one-sided p v fluxes.
    integer(kind=ENTIER), intent(in) :: ep_mode
    ! lag_mode : forme du flux d ENERGIE de la partie Lagrange. La partie
    ! QUANTITE DE MOUVEMENT reste pp*norm dans tous les cas ; seul le flux
    ! d energie change, et 0 reproduit l existant AU BIT PRES.
    !   0 : pp * vstar      -- l existant. vstar est la vitesse etoile 1D DE
    !       FACE, PAS la vitesse nodale : c est une difference structurelle
    !       avec le schema lagrangien, qui est bati sur des quantites nodales.
    !   1 : pp * (v_p.n)    -- tout nodal : la forme lagrangienne telle quelle.
    !   2 : pp * u*_s       -- par cote, u*_l = vnl + (pl-pp)/lambda_l et
    !       u*_r = vnr + (pp-pr)/lambda_r.
    !   3 : Pi_s * (v_p.n)  -- par cote, pression de sous-face CLASSIQUE.
    ! eps* de chacune, en variables specifiques, verifie a 2e-16 :
    !   2 -> eps* = eps + (pp^2 - p^2)/(2 lambda^2), SANS AUCUNE VITESSE, donc
    !        positivite automatique par Einfeldt-Menikoff-Plohr ;
    !   1 -> eps* = eps + pp (u - v_p.n)/lambda - ...   terme de vitesse, non borne ;
    !   3 -> eps* = eps + u (d_Pi - d) - ...            terme de vitesse ET non
    !        galileen, puisque la pression du flux d energie (Pi_s) differe de
    !        celle du flux de quantite de mouvement (pp).
    integer(kind=ENTIER), intent(in) :: lag_mode

    integer(kind=ENTIER) :: j, id_sub_face, id_sub_elem, id_elem
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE), dimension(3) :: norm, v_p

    real(kind=DOUBLE) :: pp, pbar_f, denomsum, weight, invlamb
    ! zb_lag 4 et 5 : pression de face et son argument dissipatif
    ! zb_pdiv : argument dissipatif de la pression etoile par sous-face
    real(kind=DOUBLE) :: pdis, pdfa, pddl, pddr, pdel, pder
    real(kind=DOUBLE) :: lambda_l, lambda_r, lambda_lts
    real(kind=DOUBLE) :: rhol, rhor, pl, pr, vnl, vnr, al, ar
    real(kind=DOUBLE) :: a_p, ma_node, theta, corr, vstar
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, sol_w
    real(kind=DOUBLE), dimension(5) :: sol_l, sol_r, sol_m
    real(kind=DOUBLE), dimension(5) :: ff_lag, ff_lag_l, ff_lag_r, ff_adv, fminus, fplus

    ! AMISO advection workspace
    real(kind=DOUBLE) :: min_apf, sum_area_a, lambda_a, wpcf, rhom, am, vm
    real(kind=DOUBLE), dimension(5) :: sol_p, u_bar, ff_a, fadv_minus, fadv_plus
    ! shock-sensor blend workspace
    logical :: need_amiso, need_ducros
    real(kind=DOUBLE) :: corr_a, w_shock, eps_face
    real(kind=DOUBLE), dimension(3) :: utl_t, utr_t, ust_t, ft_t
    real(kind=DOUBLE) :: dut_t, nul_t, nur_t
    real(kind=DOUBLE) :: min_apf_md, sum_area_md, h_node, w_md
    real(kind=DOUBLE), dimension(3) :: v_md
    real(kind=DOUBLE), dimension(5) :: sol_p_md
    real(kind=DOUBLE), dimension(5,3) :: grad_md, fp_md
    real(kind=DOUBLE), dimension(3,3) :: smax_md
    real(kind=DOUBLE) :: h_l, h_r, h_bar
    real(kind=DOUBLE) :: theta_n, pflux_l, pflux_r, vn_node
    real(kind=DOUBLE), dimension(3) :: u_node
    real(kind=DOUBLE), dimension(5) :: djump
    real(kind=DOUBLE) :: e_p, ep_num, ep_den, e_l, e_r, mu_l, mu_r, vs_ep
    logical :: ep_active

    rse_loc = 0

    ! nodal Mach number (volume-weighted), drives the low-Mach scaling
    a_p = 0.0_DOUBLE
    v_p = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_w = conserv_to_primit(sol(:, id_elem))
      a_p = a_p + mesh%sub_elem(id_sub_elem)%volume * sound_speed_w(sol_w)
      v_p = v_p + mesh%sub_elem(id_sub_elem)%volume * sol_w(2:4)
    end do
    a_p = a_p / mesh%vert(id_vert)%volume
    v_p = v_p / mesh%vert(id_vert)%volume
    ma_node = norm2(v_p) / a_p

    if (low_mach) then
      theta = min(1.0_DOUBLE, ma_node)
    else
      theta = 1.0_DOUBLE
    end if

    ! nodal velocity + its Mach, needed by the Tallois multidimensional correction
    theta_n = 0.0_DOUBLE
    u_node = 0.0_DOUBLE
    if (tp_fix .or. lag_mode == 1 .or. lag_mode == 3 .or. zb_pivp .or. zb_lag > 0) then
      call compute_nodal_velocity_LVP(mesh, id_vert, sol, grad, u_node, second_order)
      if (tp_fix) theta_n = min(1.0_DOUBLE, norm2(u_node)/a_p)
    end if

    ! carbuncle / shock sensor: each advection variant keeps its own, so that
    ! switching adv_mode reproduces that family's behaviour exactly
    need_ducros = (adv_mode == 0 .or. adv_mode == 2 .or. adv_mode == 3 .or. adv_mode == 4 .or. adv_mode == 5)

    corr = 0.0_DOUBLE
    corr_a = 0.0_DOUBLE
    if (need_ducros .and. eps_mode /= 4 .and. eps_mode /= 5 .and. eps_mode /= 8) then
      select case (eps_mode)
      case (1)
        call compute_corr_pressure(mesh, id_vert, sol, grad, corr, second_order)
      case (2)
        call compute_corr2(mesh, id_vert, sol, grad, corr, second_order)
      case (3)
        call compute_corr_pressure_div(mesh, id_vert, sol, grad, corr, second_order)
      case default
        call compute_corr_ducros(mesh, id_vert, sol, grad, corr, second_order)
      end select
      ! Modes 6/7 keep the Ducros sensor EXACTLY as is and only rescale it.
      ! Purpose: the heat-flux campaign left an ambiguity -- the nodal -div(v)
      ! sensor (amplitude 1*a_p) is the best on quad while Ducros (amplitude
      ! 4*a_p) is the best on tri by a factor 4.5. That gap can come from the
      ! sensor's SHAPE (the div^2/(div^2+curl^2) filter suppressing shear) or
      ! merely from its AMPLITUDE. compute_corr_ducros returns w*4*a_p, so
      ! multiplying by 1/4 and 1/2 isolates the amplitude at fixed shape.
      if (eps_mode == 6) corr = 0.25_DOUBLE*corr
      if (eps_mode == 7) corr = 0.5_DOUBLE*corr
    end if

    ! eps_mode 8 : AUCUN capteur. La stabilisation de l advection est la VITESSE
    ! SEULE, |vstar|, plus un petit coefficient constant eps = zb_eps_coef*a_p.
    ! Ce n est donc ni un detecteur de choc ni un filtre : pas de seuil, pas de
    ! porte en Mach, pas de rapport div/curl. zb_eps_coef = 0 donne exactement
    ! le decentrement a vstar sans aucune dissipation ajoutee, ce qui est la
    ! borne inferieure de la famille et sert de reference au balayage.
    if (eps_mode == 8) corr = zb_eps_coef*a_p

    ! Blend weight: compute_corr_ducros returns corr = w * 4 * a_p, where w is
    ! its dimensionless [0,1] shock indicator (Ducros filter gated by Mach and
    ! by compression). Recover w by dividing out the 4*a_p it multiplied in,
    ! so the blend follows the same sensor the dissipation already uses.
    ! NB: this divides out compute_corr_ducros's own 4*a_p, so it is only valid
    ! for eps_mode == 0 -- the only blended scheme (WIP2_HYB) is registered so.
    w_shock = 0.0_DOUBLE
    if (adv_mode == 2) then
      w_shock = min(1.0_DOUBLE, max(0.0_DOUBLE, corr/(4.0_DOUBLE*a_p)))
    end if

    ! With w_shock == 0 the AMISO contribution is multiplied by zero, so skip
    ! it entirely (exactly equivalent, and it is the common case: the sensor
    ! only fires at shocks, leaving most of a mesh on the Ducros branch alone).
    need_amiso = (adv_mode == 1) .or. (adv_mode == 2 .and. w_shock > 0.0_DOUBLE)
    ! corr_a est le SECOND capteur, celui d AMISO (-div v), calcule
    ! independamment de eps_mode. Le laisser actif sous eps_mode 8 garderait un
    ! detecteur dans l advection de adv_mode 1 et 2 alors que le mode 8 est
    ! justement defini comme "la vitesse seule" -- la contradiction serait
    ! silencieuse. Il suit donc la meme regle que corr.
    if (need_amiso .and. eps_mode /= 8) then
      call compute_corr2(mesh, id_vert, sol, grad, corr_a, second_order)
    else if (need_amiso) then
      corr_a = zb_eps_coef*a_p
    end if

    ! --- AMISO advection: nodal upwind state sol_p (first pass) ---
    if (need_amiso) then
      min_apf = huge(1.0_DOUBLE)
      do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
        id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
        re = mesh%sub_face(id_sub_face)%right_elem_neigh
        if (re > 0) then
          ! Une sous-face z de maillage faux-2D est inerte : si elle fixait le minimum
          ! d'aire, wpcf vaudrait h/t au lieu de 1 et iso cesserait d'egaler multi_point.
          if (.not. (boundary_2d .and. &
              abs(mesh%sub_face(id_sub_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE)) &
            min_apf = min(min_apf, mesh%sub_face(id_sub_face)%area)
        end if
      end do

      sol_p = 0.0_DOUBLE
      sum_area_a = 0.0_DOUBLE
      do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
        id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
        le = mesh%sub_face(id_sub_face)%left_elem_neigh
        re = mesh%sub_face(id_sub_face)%right_elem_neigh
        norm = mesh%sub_face(id_sub_face)%norm

        call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
          second_order, sol_w_l, sol_w_r)

        vnl = dot_product(sol_w_l(2:4), norm)
        vnr = dot_product(sol_w_r(2:4), norm)
        lambda_a = max(1e-8_DOUBLE, -vnl, vnr) + corr_a

        sol_l = primit_to_conserv(sol_w_l)
        sol_r = primit_to_conserv(sol_w_r)
        u_bar = sol_l*0.5_DOUBLE*(1.0_DOUBLE + vnl/lambda_a) &
              + sol_r*0.5_DOUBLE*(1.0_DOUBLE - vnr/lambda_a)

        if (re > 0) then
          wpcf = min_apf/mesh%sub_face(id_sub_face)%area
          sol_p = sol_p + mesh%sub_face(id_sub_face)%area*wpcf*lambda_a*u_bar
          sum_area_a = sum_area_a + mesh%sub_face(id_sub_face)%area*wpcf*lambda_a
        end if
      end do
      sol_p = sol_p / sum_area_a
    end if

    ! --- nodal pressure with the Mach-scaled velocity-jump (divv) term ---
    pp = 0.0_DOUBLE
    denomsum = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      norm = mesh%sub_face(id_sub_face)%norm

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      rhol = sol_w_l(1)
      vnl = dot_product(sol_w_l(2:4), norm)
      pl = sol_w_l(5)
      al = sound_speed_w(sol_w_l)

      rhor = sol_w_r(1)
      vnr = dot_product(sol_w_r(2:4), norm)
      pr = sol_w_r(5)
      ar = sound_speed_w(sol_w_r)

      lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(vnr - vnl))
      lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(vnr - vnl))

      invlamb = 1.0_DOUBLE/lambda_l + 1.0_DOUBLE/lambda_r
      ! zb_pdiv : argument dissipatif de la pression etoile PAR SOUS-FACE.
      !   0 : (v_nr - v_nl), le SAUT normal (existant, bit pour bit)
      !   1 : ell_pf * D_pf, une VRAIE DIVERGENCE
      ! Le saut vaut ell * d_n u_n : un seul terme de la trace. Il ne s'annule
      ! donc PAS sur un champ a divergence nulle, alors que la trace complete
      ! d_n u_n + d_t u_t s'annule. Seul cet argument change : la structure du
      ! flux, le poids en 1/lambda et la moyenne nodale sont intacts.
      ! ell = 1/2 (V_l + V_r)/A_face, extension normale intrinseque, JAMAIS une
      ! difference de centroides (le recollement periodique renvoie la maille de
      ! l'autre bout du domaine sans translation).
      pdis = vnr - vnl
      ! BUG EVITE : la boucle qui construit p_p n'exclut PAS les sous-faces z de
      ! l'extrusion (celle de e_p, elle, les exclut). Avec le saut leur
      ! contribution vaut (v_zr - v_zl) ~ 0 et c'est benin ; avec la divergence
      ! du PLAN elle ne l'est pas, et polluerait chaque pression nodale. On y
      ! garde donc le saut, de sorte que zb_pdiv ne differe que sur les faces
      ! DANS le plan.
      if (zb_pdiv == 1 .and. (second_order .or. grad_filled) &
          .and. .not. (boundary_2d .and. abs(norm(3)) > 1e-3)) then
        pdfa = mesh%face(mesh%sub_face(id_sub_face)%mesh_face)%area
        pddl = grad(1, 2, le) + grad(2, 3, le)
        if (.not. boundary_2d) pddl = pddl + grad(3, 4, le)
        pdel = mesh%elem(le)%volume/max(pdfa, 1e-300_DOUBLE)
        if (re > 0) then
          pddr = grad(1, 2, re) + grad(2, 3, re)
          if (.not. boundary_2d) pddr = pddr + grad(3, 4, re)
          pder = mesh%elem(re)%volume/max(pdfa, 1e-300_DOUBLE)
          pdis = 0.5_DOUBLE*(pdel + pder)*0.5_DOUBLE*(pddl + pddr)
        else
          ! bord reel : divergence de la seule maille gauche, pour que
          ! << nul sur un champ solenoidal >> reste UNIFORME jusqu'au bord
          pdis = pdel*pddl
        end if
      end if
      pbar_f = (pl/lambda_l + pr/lambda_r - theta*pdis)/invlamb

      weight = mesh%sub_face(id_sub_face)%area*invlamb
      if (re <= 0) weight = 0.5_DOUBLE*weight
      pp = pp + weight*pbar_f
      denomsum = denomsum + weight
    end do
    pp = pp / denomsum

    ! --- nodal total energy e_p ('vp_pp_ep' port), interior sub-faces only ---
    e_p = 0.0_DOUBLE
    if (ep_mode > 0) then
      ep_num = 0.0_DOUBLE
      ep_den = 0.0_DOUBLE
      do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
        id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
        le = mesh%sub_face(id_sub_face)%left_elem_neigh
        re = mesh%sub_face(id_sub_face)%right_elem_neigh
        norm = mesh%sub_face(id_sub_face)%norm
        if (re <= 0) cycle
        if (boundary_2d .and. abs(norm(3)) > 1e-3) cycle

        call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
          second_order, sol_w_l, sol_w_r)

        sol_l = primit_to_conserv(sol_w_l)
        rhol = sol_w_l(1)
        vnl = dot_product(sol_w_l(2:4), norm)
        pl = sol_w_l(5)
        al = sound_speed_w(sol_w_l)
        e_l = sol_l(5)/rhol
        if (ep_mode == 2 .or. ep_mode == 4) e_l = e_l - 0.5_DOUBLE*dot_product(sol_w_l(2:4), sol_w_l(2:4))

        sol_r = primit_to_conserv(sol_w_r)
        rhor = sol_w_r(1)
        vnr = dot_product(sol_w_r(2:4), norm)
        pr = sol_w_r(5)
        ar = sound_speed_w(sol_w_r)
        e_r = sol_r(5)/rhor
        if (ep_mode == 2 .or. ep_mode == 4) e_r = e_r - 0.5_DOUBLE*dot_product(sol_w_r(2:4), sol_w_r(2:4))

        lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(vnr - vnl))
        lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(vnr - vnl))

        if (ep_mode >= 3) then
          vs_ep = (lambda_l*vnl + lambda_r*vnr - (pr - pl))/(lambda_l + lambda_r)
          mu_l = rhol*abs(vnl - vs_ep)
          mu_r = rhor*abs(vnr - vs_ep)
        else
          mu_l = lambda_l
          mu_r = lambda_r
        end if

        if (zb_pivp .or. zb_ppus) then
          ! le travail de pression telescope SEUL : X_p devient une combinaison
          ! convexe, bornee quel que soit mu. Voir zb_pivp / zb_ppus.
          ep_num = ep_num + mesh%sub_face(id_sub_face)%area*(mu_l*e_l + mu_r*e_r)
        else
          ep_num = ep_num + mesh%sub_face(id_sub_face)%area &
            *(mu_l*e_l + mu_r*e_r + (pl*vnl - pr*vnr))
        end if
        ep_den = ep_den + mesh%sub_face(id_sub_face)%area*(mu_l + mu_r)
      end do
      if (ep_den > 0.0_DOUBLE) e_p = ep_num/ep_den
    end if
    ep_active = (ep_mode > 0)
    if (ep_mode >= 3) ep_active = (ep_den > 0.0_DOUBLE)

    ! --- multi-d nodal advection flux (modes 3/4) ---
    if (adv_mode == 3 .or. adv_mode == 4 .or. adv_mode == 5) then
      call compute_nodal_velocity_LVP(mesh, id_vert, sol, grad, v_md, second_order)

      sol_p_md = 0.0_DOUBLE
      do j = 1, mesh%vert(id_vert)%n_sub_elems_neigh
        id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
        id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
        sol_p_md = sol_p_md + mesh%sub_elem(id_sub_elem)%volume*sol(:, id_elem)
      end do
      sol_p_md = sol_p_md / mesh%vert(id_vert)%volume

      grad_md = 0.0_DOUBLE
      sum_area_md = 0.0_DOUBLE
      min_apf_md = huge(1.0_DOUBLE)
      do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
        id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
        le = mesh%sub_face(id_sub_face)%left_elem_neigh
        re = mesh%sub_face(id_sub_face)%right_elem_neigh
        norm = mesh%sub_face(id_sub_face)%norm
        call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
          second_order, sol_w_l, sol_w_r)
        grad_md = grad_md + tensor_product( &
          primit_to_conserv(sol_w_r) - primit_to_conserv(sol_w_l), &
          mesh%sub_face(id_sub_face)%area*norm)
        sum_area_md = sum_area_md + mesh%sub_face(id_sub_face)%area
        if (re > 0) min_apf_md = min(min_apf_md, mesh%sub_face(id_sub_face)%area)
      end do
      grad_md = grad_md / mesh%vert(id_vert)%volume
      h_node = mesh%vert(id_vert)%volume / sum_area_md

      if (adv_mode == 5) then
        ! Haenel-consistent energy row: grad(rho*E) -> h_p * grad(rho)
        sol_w = conserv_to_primit(sol_p_md)
        h_bar = (sol_p_md(5) + sol_w(5))/sol_p_md(1)
        grad_md(5, :) = h_bar*grad_md(1, :)
      end if

      smax_md = 0.0_DOUBLE
      smax_md(1, 1) = abs(v_md(1))
      smax_md(2, 2) = abs(v_md(2))
      smax_md(3, 3) = abs(v_md(3))
      fp_md = tensor_product(sol_p_md, v_md) &
        - 0.5_DOUBLE*h_node*matmul(grad_md, smax_md)
    end if

    ! --- flux assembly ---
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if (rse > 0) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if
      norm = mesh%sub_face(id_sub_face)%norm

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      sol_l = primit_to_conserv(sol_w_l)
      rhol = sol_w_l(1)
      vnl = dot_product(sol_w_l(2:4), norm)
      pl = sol_w_l(5)
      al = sound_speed_w(sol_w_l)

      sol_r = primit_to_conserv(sol_w_r)
      rhor = sol_w_r(1)
      vnr = dot_product(sol_w_r(2:4), norm)
      pr = sol_w_r(5)
      ar = sound_speed_w(sol_w_r)

      lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(vnr - vnl))
      lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(vnr - vnl))

      ! symmetric, impedance-weighted interface velocity (conservative: fplus=fminus)
      vstar = (lambda_l*vnl + lambda_r*vnr - (pr - pl))/(lambda_l + lambda_r)

      sol_m = 0.5_DOUBLE*(sol_r + sol_l)

      ff_lag(1) = 0.0_DOUBLE
      if (tp_fix) then
        ! one-sided EUCCLHYD pressure at the shock end, nodal pressure at the
        ! low-Mach end; energy velocity blends the same way so theta_n=0 is
        ! exactly the WIP2_NOLM flux.
        vn_node = dot_product(u_node, norm)
        pflux_l = theta_n*(pl - lambda_l*(vn_node - vnl)) + (1.0_DOUBLE - theta_n)*pp
        pflux_r = theta_n*(pr + lambda_r*(vn_node - vnr)) + (1.0_DOUBLE - theta_n)*pp
        ff_lag(5) = 0.0_DOUBLE
      else
        ff_lag(2:4) = pp * norm
        ff_lag(5) = pp * vstar
      end if

      ! --- dissipation tangentielle (voir zb_tang) -------------------------
      ft_t = 0.0_DOUBLE ; ust_t = 0.0_DOUBLE
      if (zb_tang > 0 .and. re > 0) then
        utl_t = sol_w_l(2:4) - dot_product(sol_w_l(2:4), norm)*norm
        utr_t = sol_w_r(2:4) - dot_product(sol_w_r(2:4), norm)*norm
        dut_t = sqrt(dot_product(utr_t - utl_t, utr_t - utl_t))
        if (zb_tang == 1) then
          nul_t = lambda_l ; nur_t = lambda_r
        else
          nul_t = rhol*dut_t ; nur_t = rhor*dut_t
        end if
        if (nul_t + nur_t > 0.0_DOUBLE) then
          ust_t = (nul_t*utl_t + nur_t*utr_t)/(nul_t + nur_t)
          ft_t  = nul_t*(ust_t - utl_t)   ! = -nur_t*(ust_t - utr_t) : UNIVALUE
          ff_lag(2:4) = ff_lag(2:4) - ft_t
          ff_lag(5)   = ff_lag(5) - dot_product(ft_t, ust_t)
        end if
      end if

      ! jump used by the advection dissipation
      djump = sol_r - sol_l
      if (enth_fix) then
        h_l = (sol_l(5) + pl)/rhol
        h_r = (sol_r(5) + pr)/rhor
        h_bar = 0.5_DOUBLE*(h_l + h_r)
        djump(5) = h_bar*(rhor - rhol)
      end if

      ! Rusanov-at-vstar advection (WIP form)
      if (need_ducros) then
        if (eps_mode == 4 .or. eps_mode == 5) then
          if (eps_mode == 5) then
            ! Mode 5 fixes mode 4's scaling. sqrt(|dp|/rho) has the dimension of
            ! a velocity but the WRONG Mach scaling: at low Mach dp across a face
            ! is O(rho v^2) = O(h), independent of Ma, so sqrt(dp/rho) ~ sqrt(h)
            ! never vanishes and puts a floor on L2(rho) (measured: 1.06e-7 at
            ! Ma=1e-4 against 6.0e-8 without it). The impedance conversion
            ! |dp|/(rho*a) is O(h*Ma) instead, and is the same conversion already
            ! used everywhere else here, e.g. in vstar's -(p_r-p_l)/(lam_l+lam_r).
            eps_face = max(0.0_DOUBLE, -(vnr - vnl), &
              abs(pr - pl)/(0.5_DOUBLE*(rhol + rhor)*0.5_DOUBLE*(al + ar)))
          else
            eps_face = max(0.0_DOUBLE, -(vnr - vnl), &
              sqrt(abs(pr - pl)/(0.5_DOUBLE*(rhol + rhor))))
          end if
          ff_adv = vstar*sol_m - 0.5_DOUBLE*(abs(vstar) + eps_face)*djump
        else
          ff_adv = vstar*sol_m - 0.5_DOUBLE*(abs(vstar) + corr)*djump
        end if
      end if

      ! AMISO advection (asymmetric by construction: the nodal upwind state
      ! sol_p is approached from each side, so fadv_minus /= fadv_plus)
      if (need_amiso) then
        rhom = 0.5_DOUBLE*(rhol + rhor)
        am = 0.5_DOUBLE*(al + ar)
        vm = 0.5_DOUBLE*(vnr + vnl) - 0.5_DOUBLE/(rhom*am)*(pr - pl)
        lambda_a = max(1e-8_DOUBLE, -vnl, vnr) + corr_a
        ff_a = vm*sol_m - 0.5_DOUBLE*(abs(vm) + corr_a)*djump

        if (re > 0) then
          wpcf = min_apf/mesh%sub_face(id_sub_face)%area * (1.0_DOUBLE/3.0_DOUBLE)
          fadv_minus = wpcf*(sol_l*vnl - lambda_a*(sol_p - sol_l)) &
            + (1.0_DOUBLE - wpcf)*ff_a
          fadv_plus = wpcf*(sol_r*vnr + lambda_a*(sol_p - sol_r)) &
            + (1.0_DOUBLE - wpcf)*ff_a
        else
          fadv_minus = ff_a
          fadv_plus = ff_a
        end if
      end if

      if (tp_fix) then
        ! per-side Lagrange flux (the EUCCLHYD part is one-sided by construction)
        ff_lag_l(1) = 0.0_DOUBLE
        ff_lag_l(2:4) = pflux_l * norm
        ff_lag_l(5) = pflux_l * (theta_n*vn_node + (1.0_DOUBLE - theta_n)*vstar)
        ff_lag_r(1) = 0.0_DOUBLE
        ff_lag_r(2:4) = pflux_r * norm
        ff_lag_r(5) = pflux_r * (theta_n*vn_node + (1.0_DOUBLE - theta_n)*vstar)
      else
        ff_lag_l = ff_lag
        ff_lag_r = ff_lag
      end if

      ! --- partie LAGRANGE du flux d energie, voir lagrange_energy_flux ------
      ! Exclusif de ep_mode : si une energie nodale est active, c est elle qui
      ! fixe le flux d energie et lag_mode ne doit pas la recouvrir.
      if (lag_mode > 0 .and. .not. tp_fix .and. ep_mode == 0) then
        call lagrange_energy_flux(lag_mode, pp, pl, pr, vnl, vnr, &
          lambda_l, lambda_r, dot_product(u_node, norm), ff_lag_l(5), ff_lag_r(5))
      end if

      if (ep_active .and. re > 0 .and. &
          .not. (boundary_2d .and. abs(norm(3)) > 1e-3)) then
        ! one-sided energy flux on the nodal energy; same lambda_l/lambda_r and the
        ! same scalar X as in the e_p pass, so the node telescopes exactly
        e_l = sol_l(5)/rhol
        e_r = sol_r(5)/rhor
        if (ep_mode == 2 .or. ep_mode == 4) then
          e_l = e_l - 0.5_DOUBLE*dot_product(sol_w_l(2:4), sol_w_l(2:4))
          e_r = e_r - 0.5_DOUBLE*dot_product(sol_w_r(2:4), sol_w_r(2:4))
        end if
        if (ep_mode >= 3) then
          mu_l = rhol*abs(vnl - vstar)
          mu_r = rhor*abs(vnr - vstar)
        else
          mu_l = lambda_l
          mu_r = lambda_r
        end if
        if (zb_ppus) then
          ! (pv)*_s = p_p (v_s + ds_s) -/+ mu_s (X_p - X_s) : la pression du
          ! travail est celle de la QdM et la vitesse est l etoile de son cote.
          ff_lag_l(5) = pp*(vnl + (pl - pp)/lambda_l) - mu_l*(e_p - e_l)
          ff_lag_r(5) = pp*(vnr + (pp - pr)/lambda_r) + mu_r*(e_p - e_r)
        else if (zb_pivp) then
          vn_node = dot_product(u_node, norm)
          ff_lag_l(5) = (pl - lambda_l*(vn_node - vnl))*vn_node - mu_l*(e_p - e_l)
          ff_lag_r(5) = (pr + lambda_r*(vn_node - vnr))*vn_node + mu_r*(e_p - e_r)
        else
          ff_lag_l(5) = pl*vnl - mu_l*(e_p - e_l)
          ff_lag_r(5) = pr*vnr + mu_r*(e_p - e_r)
        end if
      end if

      ! --- SELECTEUR COMPLET de la partie lagrangienne (voir zb_lag) --------
      ! Surcharge TOUT ce qui precede : les trois composantes sont posees ici
      ! ensemble, donc aucun melange avec lag_mode, ep_mode ou tp_fix.
      if (zb_lag > 0) then
        vn_node = dot_product(u_node, norm)
        select case (zb_lag)
        case (1)   ! vp : pressions de sous-face UNILATERALES
          ff_lag_l(1) = 0.0_DOUBLE
          ff_lag_l(2:4) = (pl - lambda_l*(vn_node - vnl))*norm
          ff_lag_l(5) = (pl - lambda_l*(vn_node - vnl))*vn_node
          ff_lag_r(1) = 0.0_DOUBLE
          ff_lag_r(2:4) = (pr + lambda_r*(vn_node - vnr))*norm
          ff_lag_r(5) = (pr + lambda_r*(vn_node - vnr))*vn_node
        case (2)   ! vp_pp : pression NODALE, travail a la vitesse nodale
          ff_lag_l(1) = 0.0_DOUBLE
          ff_lag_l(2:4) = pp*norm
          ff_lag_l(5) = pp*vn_node
          ff_lag_r = ff_lag_l
        case (3)   ! epv : pression nodale, travail a la vitesse etoile du cote
          ff_lag_l(1) = 0.0_DOUBLE
          ff_lag_l(2:4) = pp*norm
          ff_lag_l(5) = pp*(vnl + (pl - pp)/lambda_l)
          ff_lag_r(1) = 0.0_DOUBLE
          ff_lag_r(2:4) = pp*norm
          ff_lag_r(5) = pp*(vnr + (pp - pr)/lambda_r)
        end select
      end if

      if (adv_mode == 3 .or. adv_mode == 4 .or. adv_mode == 5) then
        if ((adv_mode == 3 .or. adv_mode == 5) .and. re > 0) then
          w_md = min_apf_md/mesh%sub_face(id_sub_face)%area*(1.0_DOUBLE/3.0_DOUBLE)
        else if (adv_mode == 4) then
          w_md = 0.5_DOUBLE
        else
          w_md = 0.0_DOUBLE
        end if
        ff_adv = w_md*matmul(fp_md, norm) + (1.0_DOUBLE - w_md)*ff_adv
      end if

      select case (adv_mode)
      case (1)
        fminus = fadv_minus + ff_lag_l
        fplus = fadv_plus + ff_lag_r
      case (2)
        if (need_amiso) then
          fminus = (1.0_DOUBLE - w_shock)*ff_adv + w_shock*fadv_minus + ff_lag_l
          fplus = (1.0_DOUBLE - w_shock)*ff_adv + w_shock*fadv_plus + ff_lag_r
        else
          fminus = ff_adv + ff_lag_l
          fplus = ff_adv + ff_lag_r
        end if
      case default
        fminus = ff_adv + ff_lag_l
        fplus = ff_adv + ff_lag_r
      end select

      if (boundary_2d &
        .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1e-3) then
        fminus = 0.0_DOUBLE
        fplus = 0.0_DOUBLE
      end if

      lambda_lts = max(abs(vnl), abs(vnr)) + max(al, ar)

      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        sum_lambda_vert(lse_loc) = sum_lambda_vert(lse_loc) &
          + mesh%sub_face(id_sub_face)%area*lambda_lts
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) &
          + mesh%sub_face(id_sub_face)%area*fminus

        if (rse > 0) then
        if (mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda_lts
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) &
            - mesh%sub_face(id_sub_face)%area*fplus
        end if
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_WIP2

  ! Partie LAGRANGE du flux d ENERGIE, isolee pour pouvoir etre remplacee sans
  ! toucher au reste de WIP2 : la partie ADVECTION et tous les selecteurs
  ! existants de la famille ZB (adv_mode, eps_mode, enth_fix, tp_fix) sont
  ! reutilises tels quels. La partie QUANTITE DE MOUVEMENT reste pp*norm dans
  ! tous les cas ; seule l energie change ici.
  !
  !   lag_mode 1 : pp * (v_p.n)   tout nodal, la forme lagrangienne telle quelle
  !   lag_mode 2 : pp * u*_s      par cote, u*_l = vnl + (pl-pp)/lambda_l
  !   lag_mode 3 : Pi_s * (v_p.n) par cote, pression de sous-face CLASSIQUE
  !
  ! Les energies internes intermediaires correspondantes, en variables
  ! specifiques (forme lagrangienne de Gallice), verifiees a 2e-16 :
  !   2 -> eps* = eps + (pp^2 - p^2)/(2 lambda^2) : AUCUNE VITESSE n y entre,
  !        donc positivite automatique par Einfeldt-Menikoff-Plohr. C est la
  !        seule forme de la famille qui ait cette propriete.
  !   1 -> eps* = eps + pp (u - v_p.n)/lambda - (p-pp)^2/(2 lambda^2) : le terme
  !        de vitesse est invariant galileen mais n a ni signe ni borne.
  !   3 -> eps* = eps + u (d_Pi - d) - ... : terme de vitesse ET non galileen,
  !        puisque la pression du flux d energie (Pi_s) differe de celle du flux
  !        de quantite de mouvement (pp). Ce terme est exactement l ecart a
  !        l invariance, et il s annule quand Pi_s = pp.
  ! Mesure en lagrangien : 1 et 3 tiennent Gresho et cassent Sedov (pas 14 et
  ! pas 1023), 2 tient Sedov et Gresho sur quadrangles.
  subroutine lagrange_energy_flux(lag_mode, pp, pl, pr, vnl, vnr, &
      lambda_l, lambda_r, vn_node, fl5, fr5)
    implicit none
    integer(kind=ENTIER), intent(in) :: lag_mode
    real(kind=DOUBLE), intent(in) :: pp, pl, pr, vnl, vnr
    real(kind=DOUBLE), intent(in) :: lambda_l, lambda_r, vn_node
    real(kind=DOUBLE), intent(out) :: fl5, fr5

    select case (lag_mode)
    case (1)
      fl5 = pp*vn_node
      fr5 = pp*vn_node
    case (2)
      fl5 = pp*(vnl + (pl - pp)/lambda_l)
      fr5 = pp*(vnr + (pp - pr)/lambda_r)
    case (3)
      fl5 = (pl - lambda_l*(vn_node - vnl))*vn_node
      fr5 = (pr + lambda_r*(vn_node - vnr))*vn_node
    case default
      ! lag_mode = 0 n appelle jamais cette routine ; si on arrive ici c est une
      ! valeur non prevue, et rendre le flux inchange masquerait l erreur.
      write(*, *) "lagrange_energy_flux : lag_mode non prevu =", lag_mode
      stop 1
    end select
  end subroutine lagrange_energy_flux

  subroutine compute_rhs_around_vert_usi3d(mesh, sol, grad, &
      nsen, flux_sum_vert, sum_lambda_vert, &
      id_vert, vp, h_p, second_order)
    use ns_global_data_module, only: bc_style, scheme, &
      exclude_bound_vert, boundary_2d
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: h_p
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    integer(kind=ENTIER), intent(in) :: nsen
    real(kind=DOUBLE), dimension(nsen), intent(inout) :: sum_lambda_vert
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: flux_sum_vert
    integer(kind=ENTIER), intent(in) :: id_vert
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, d, id_sub_face, n_neigh
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc, id_sub_elem, id_elem
    real(kind=DOUBLE), dimension(3) :: norm, vp_node, vface, gradp_h
    real(kind=DOUBLE) :: pp_node, corrp, cp_node, lambda_lts, hp_vert, rho_p
    real(kind=DOUBLE) :: pp_avg, div_v_h, area_f, vol
    real(kind=DOUBLE) :: rhol, rhor, vnl, vnr, pl, pr, cL, cR
    real(kind=DOUBLE) :: rhom, am, up, pp_face, smax
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, sol_l, sol_r, sol_m, sol_w, Qface
    real(kind=DOUBLE), dimension(5) :: ff_1d, fminus, fplus, Qp, Fp_n
    real(kind=DOUBLE), dimension(5, 3) :: Fp, gradQ
    real(kind=DOUBLE) :: cfweight, sum_area

    cfweight = 1.0_DOUBLE / 3.0_DOUBLE
    !cfweight = 0.0_DOUBLE
    !cfweight = 1e-3_DOUBLE
    rse_loc  = 0
    hp_vert  = h_p(id_vert)
    vol      = mesh%vert(id_vert)%volume
    n_neigh  = mesh%vert(id_vert)%n_sub_elems_neigh

    ! --- volume-weighted averages over neighbouring sub-elements ---
    Qp      = 0.0_DOUBLE
    pp_avg  = 0.0_DOUBLE
    vp_node = 0.0_DOUBLE
    do j = 1, n_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem     = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_w       = conserv_to_primit(sol(:, id_elem))
      Qp          = Qp      + mesh%sub_elem(id_sub_elem)%volume * sol(:, id_elem)
      pp_avg      = pp_avg  + mesh%sub_elem(id_sub_elem)%volume * sol_w(5)
      vp_node     = vp_node + mesh%sub_elem(id_sub_elem)%volume * sol_w(2:4)
    end do
    Qp      = Qp      / vol
    pp_avg  = pp_avg  / vol
    vp_node = vp_node / vol
    sol_w   = conserv_to_primit(Qp)
    cp_node = sound_speed_w(sol_w)
    rho_p   = sol_w(1)

    ! --- Gauss-based div(v), grad(p), grad(Q) using face-average values ---
    ! All faces included (boundary faces handled via reconstruct_lr_w) so that
    ! sum(area*norm) = 0 over the closed dual cell and gradQ = 0 for uniform flow.
    div_v_h  = 0.0_DOUBLE
    gradp_h  = 0.0_DOUBLE
    gradQ    = 0.0_DOUBLE
    sum_area = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le     = mesh%sub_face(id_sub_face)%left_elem_neigh
      re     = mesh%sub_face(id_sub_face)%right_elem_neigh
      norm   = mesh%sub_face(id_sub_face)%norm
      area_f = mesh%sub_face(id_sub_face)%area
      if (boundary_2d .and. abs(norm(3)) > 1e-3_DOUBLE) cycle

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      div_v_h  = div_v_h + area_f * dot_product(sol_w_r(2:4)-sol_w_l(2:4), norm)
      gradp_h  = gradp_h + area_f * (sol_w_r(5) - sol_w_l(5)) * norm
      gradQ    = gradQ   + area_f * tensor_product(primit_to_conserv(sol_w_r)- primit_to_conserv(sol_w_l), norm)
      sum_area = sum_area + area_f
    end do

    ! Green-Gauss: divide by dual cell volume; multiply by h_p to get dimensionless quantities
    div_v_h = div_v_h / vol * hp_vert
    gradp_h = gradp_h / vol * hp_vert
    gradQ   = gradQ   / sum_area

    ! --- acoustic nodal pressure and velocity (original Sidilkover formulas) ---
    pp_node = pp_avg  - 0.5_DOUBLE * rho_p * cp_node * div_v_h


    ! --- dimensionless shock sensor (original corrp, in [0,1]) ---
    corrp = max(0.0_DOUBLE, min(1.0_DOUBLE, &
      abs(div_v_h) / cp_node + norm2(gradp_h) / (rho_p * cp_node**2)))
    !call compute_corr2(mesh, id_vert, sol, grad, corrp, second_order)

    vp(:, id_vert) = vp_node
    if( mesh%vert(id_vert)%is_bound ) then
      vface = wall_normal(mesh, id_vert)
      if( norm2(vface) > 1e-12_DOUBLE ) then
        vface   = vface / norm2(vface)
        vp_node = vp_node - dot_product(vp_node, vface)*vface
        vp(:, id_vert) = vp_node
      end if
    end if

    ! --- nodal flux tensor: coefficient 1/(nDim+1)=1/3, per direction d ---
    do d = 1, 3
      Fp(:, d) = vp_node(d) * Qp(:) &
        - 0.5_DOUBLE * (abs(vp_node(d)) + corrp * cp_node) * gradQ(:, d)
    end do
    Fp(2:4, :) = Fp(2:4, :) + pp_node * eye3
    Fp(5, :)   = Fp(5, :)   + vp_node * pp_node

    ! --- flux loop ---
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le   = mesh%sub_face(id_sub_face)%left_elem_neigh
      re   = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse  = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse  = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if (rse > 0) rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      norm = mesh%sub_face(id_sub_face)%norm

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      sol_l = primit_to_conserv(sol_w_l)
      sol_r = primit_to_conserv(sol_w_r)
      sol_m = 0.5_DOUBLE * (sol_l + sol_r)

      rhol = sol_w_l(1);  rhor = sol_w_r(1)
      vnl  = dot_product(sol_w_l(2:4), norm)
      vnr  = dot_product(sol_w_r(2:4), norm)
      pl   = sol_w_l(5);  pr   = sol_w_r(5)
      cL   = sound_speed_w(sol_w_l)
      cR   = sound_speed_w(sol_w_r)

      rhom    = 0.5_DOUBLE * (rhol + rhor)
      am      = max(cL, cR)
      up      = 0.5_DOUBLE * (vnl + vnr) - 0.5_DOUBLE / (rhom * am) * (pr - pl)
      pp_face = 0.5_DOUBLE * (pl  + pr ) - 0.5_DOUBLE * rhom * am * (vnr - vnl)
      ! corrp is dimensionless, am is [m/s]: consistent with original smax = |up| + cmax*am
      smax    = abs(up) + corrp * am

      ! 1D Sidilkover flux
      ff_1d      = up * sol_m - 0.5_DOUBLE * smax * (sol_r - sol_l)
      ff_1d(2:4) = ff_1d(2:4) + pp_face * norm
      ff_1d(5)   = ff_1d(5)   + up * pp_face

      ! composite: (1/3)*nodal + (2/3)*1D
      Fp_n   = matmul(Fp, norm)
      fminus = cfweight * Fp_n + (1.0_DOUBLE - cfweight) * ff_1d
      fplus  = fminus

      if (boundary_2d &
        .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1e-3_DOUBLE) then
        fminus = 0.0_DOUBLE
        fplus  = 0.0_DOUBLE
      end if

      lambda_lts = max(norm2(vp_node), abs(vnl), abs(vnr)) + max(am, cp_node)

      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        sum_lambda_vert(lse_loc) = sum_lambda_vert(lse_loc) &
          + mesh%sub_face(id_sub_face)%area * lambda_lts
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) &
          + mesh%sub_face(id_sub_face)%area * fminus

        if (rse > 0) then
        if (mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area * lambda_lts
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) &
            - mesh%sub_face(id_sub_face)%area * fplus
        end if
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_usi3d
end module ns_euler_zb_module