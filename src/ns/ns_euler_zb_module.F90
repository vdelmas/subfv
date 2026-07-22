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

  subroutine compute_nodal_velocity_LSM(mesh, id_vert, sol, grad, vp, mat_h_p, second_order)
    use ns_global_data_module, only : boundary_2d, scheme
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(3, 3, mesh%n_vert), intent(in) :: mat_h_p
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
        * mesh%sub_face(id_sub_face)%area&
        * matmul(mat_h_p(:, :, id_vert), mesh%sub_face(id_sub_face)%norm)
    end do
    grad_p = grad_p / mesh%vert(id_vert)%volume

    if( mesh%vert(id_vert)%is_bound ) then
      Bp = wall_normal(mesh, id_vert)
      if( norm2(Bp) > 1e-12_DOUBLE ) then
        Bp = Bp/norm2(Bp)
        grad_p = grad_p - dot_product(grad_p, Bp)*Bp
        vp(:, id_vert) = vp(:, id_vert) - dot_product(vp(:, id_vert), Bp)*Bp
      end if
    end if

    vp(:, id_vert) = vp(:, id_vert) - 0.5_DOUBLE/(rho_p*a_p)*grad_p
  end subroutine compute_nodal_velocity_LSM

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

  subroutine compute_nodal_pressure_LSM(mesh, id_vert, sol, grad, pp, mat_h_p, second_order)
    use ns_global_data_module, only : scheme
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, 3, mesh%n_vert), intent(in) :: mat_h_p
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
        * matmul(mat_h_p(:, :, id_vert), mesh%sub_face(id_sub_face)%norm))
    end do
    div_v = div_v / mesh%vert(id_vert)%volume

    pp = pp - 0.5_DOUBLE*rho_p*a_p*div_v
  end subroutine compute_nodal_pressure_LSM

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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_ARMDUMAT

  subroutine compute_rhs_around_vert_ARMDMMAT(mesh, sol, grad, &
      nsen, sum_lambda_vert, flux_sum_vert, &
      second_order, id_vert, vp, mat_h_p)
    use ns_global_data_module, only: bc_style, scheme, &
      exclude_bound_vert, boundary_2d
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
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

    real(kind=DOUBLE) :: vnl, vnr, al, ar, lambda, wpcf
    real(kind=DOUBLE), dimension(5) :: sol_p, sol_l, sol_r, ff
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, fminus, fplus
    real(kind=DOUBLE), dimension(5,3) :: fp, grad_sol
    real(kind=DOUBLE), dimension(3, 3) :: SMAX

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
        mesh%sub_face(id_sub_face)%area&
        *matmul(mat_h_p(:, :, id_vert), mesh%sub_face(id_sub_face)%norm))
    end do
    grad_sol = grad_sol / mesh%vert(id_vert)%volume

    call compute_nodal_velocity_LSM(mesh, id_vert, sol, grad, vp, mat_h_p, second_order)

    SMAX = 0.0_DOUBLE
    SMAX(1, 1) = abs(vp(1, id_vert))
    SMAX(2, 2) = abs(vp(2, id_vert))
    SMAX(3, 3) = abs(vp(3, id_vert))
    fp = tensor_product(sol_p, vp(:, id_vert)) &
      - 0.5_DOUBLE*matmul(grad_sol, SMAX)

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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_ARMDMMAT

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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_ARMDU

  subroutine compute_rhs_around_vert_ARMDM(mesh, sol, grad, &
      nsen, sum_lambda_vert, flux_sum_vert, &
      second_order, id_vert, vp, mat_h_p)
    use ns_global_data_module, only: bc_style, scheme, &
      exclude_bound_vert, boundary_2d
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
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
        mesh%sub_face(id_sub_face)%area&
        *matmul(mat_h_p(:, :, id_vert), mesh%sub_face(id_sub_face)%norm))
    end do
    grad_sol = grad_sol / mesh%vert(id_vert)%volume

    call compute_nodal_velocity_LSM(mesh, id_vert, sol, grad, vp, mat_h_p, second_order)

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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_ARMDM

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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_LS

  subroutine compute_rhs_around_vert_LSM(mesh, sol, grad, &
      nsen, flux_sum_vert, &
      id_vert, vp, mat_h_p, second_order)
    use ns_global_data_module, only: bc_style, scheme, &
      exclude_bound_vert, boundary_2d
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
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

    rse_loc = 0
    call compute_nodal_velocity_LSM(mesh, id_vert, sol, grad, vp, mat_h_p, second_order)
    call compute_nodal_pressure_LSM(mesh, id_vert, sol, grad, pp, mat_h_p, second_order)

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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_LSM

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
    real(kind=DOUBLE), dimension(3) :: norm, Bp, vl, vr

    real(kind=DOUBLE) :: sum_area, vnl, vnr, al, ar, lambda, corr
    real(kind=DOUBLE), dimension(5) :: sol_p, sol_l, sol_r, u_bar
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, fminus, fplus, ff

    logical :: vert_is_wall 

    rse_loc = 0
    Bp = wall_normal(mesh, id_vert)
    if( norm2(Bp) > 1e12_DOUBLE ) then
      Bp = Bp / norm2(Bp)
      vert_is_wall = .true.
    else
      vert_is_wall = .false.
    end if

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

      if( vert_is_wall ) then
        vl = sol_w_l(2:4)
        vl = vl - dot_product(vl, Bp)/dot_product(Bp, Bp)*Bp
        vnl = dot_product(sol_w_l(2:4), norm)
        vr = sol_w_r(2:4)
        vr = vr - dot_product(vr, Bp)/dot_product(Bp, Bp)*Bp
        vnr = dot_product(sol_w_r(2:4), norm)
      else
        vnl = dot_product(sol_w_l(2:4), norm)
        vnr = dot_product(sol_w_r(2:4), norm)
      end if

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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
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
    real(kind=DOUBLE), dimension(3) :: norm, Bp, vl, vr

    real(kind=DOUBLE) :: sum_area, vnl, vnr, al, ar, lambda, corr
    real(kind=DOUBLE) :: wpcf, min_apf, pr, pl, rhol, rhor
    real(kind=DOUBLE), dimension(5) :: sol_p, sol_l, sol_r, u_bar, ff, sol_m
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, fminus, fplus
    real(kind=DOUBLE) :: vm, rhom, am

    logical :: vert_is_wall 

    rse_loc = 0
    Bp = wall_normal(mesh, id_vert)
    if( norm2(Bp) > 1e12_DOUBLE ) then
      Bp = Bp / norm2(Bp)
      vert_is_wall = .true.
    else
      vert_is_wall = .false.
    end if

    min_apf = huge(1.0_DOUBLE)
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      if( re > 0 ) then
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

      if( vert_is_wall ) then
        vl = sol_w_l(2:4)
        vl = vl - dot_product(vl, Bp)/dot_product(Bp, Bp)*Bp
        vnl = dot_product(sol_w_l(2:4), norm)
        vr = sol_w_r(2:4)
        vr = vr - dot_product(vr, Bp)/dot_product(Bp, Bp)*Bp
        vnr = dot_product(sol_w_r(2:4), norm)
      else
        vnl = dot_product(sol_w_l(2:4), norm)
        vnr = dot_product(sol_w_r(2:4), norm)
      end if

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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
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

  function compute_ellip(mesh, id_vert) result(mat_h_p)
    use ns_global_data_module, only: boundary_2d, error_2d_h
    use linear_solver_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(3, 3) :: mat_h_p

    integer(kind=ENTIER) :: iter, maxiter
    integer(kind=ENTIER) :: nc
    integer(kind=ENTIER) :: j, id_sub_face, id_face
    integer(kind=ENTIER) :: id_elem, id_sub_elem
    integer(kind=ENTIER) :: idvp, idvm, k, idv, kp, km
    real(kind=DOUBLE), dimension(3) :: ce1, ce2, cp, Ltx, cf, ce
    real(kind=DOUBLE), dimension(6) :: L, delta_L
    real(kind=DOUBLE), dimension(6, 6) :: LLt
    real(kind=DOUBLE), dimension(3, 3) :: Q
    real(kind=DOUBLE), dimension(3) :: Lambda
    real(kind=DOUBLE), dimension(:), allocatable :: fL
    real(kind=DOUBLE), dimension(:, :), allocatable :: xc
    real(kind=DOUBLE), dimension(:, :), allocatable :: JL

    cp = mesh%vert(id_vert)%coord

    nc = 3*mesh%vert(id_vert)%n_sub_faces_neigh + mesh%vert(id_vert)%n_sub_elems_neigh
    !nc = 3*mesh%vert(id_vert)%n_sub_faces_neigh
    allocate(xc(3, nc))
    allocate(JL(nc, 6))
    allocate(fL(nc))

    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      idvp = 0
      idvm = 0
      do k = 1, mesh%face(id_face)%n_vert
        idv = mesh%face(id_face)%vert(k)
        if (idv == id_vert) then
          kp = 1 + mod((k+1)-1, mesh%face(id_face)%n_vert)
          km = 1 + mod((k-1)-1 + mesh%face(id_face)%n_vert, mesh%face(id_face)%n_vert)
          idvm = mesh%face(id_face)%vert(km)
          idvp = mesh%face(id_face)%vert(kp)
          exit
        end if
      end do
      ce1 = 0.5_DOUBLE*(mesh%vert(idvp)%coord + cp)
      ce2 = 0.5_DOUBLE*(mesh%vert(idvm)%coord + cp)
      cf = mesh%face(id_face)%coord
      xc(:, 3*(j-1) + 1) = 2*(ce1 - cp)
      xc(:, 3*(j-1) + 2) = 2*(ce2 - cp)
      xc(:, 3*(j-1) + 3) = 2*(cf - cp)
    end do

    do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      ce = mesh%elem(id_elem)%coord
      xc(:, 3*mesh%vert(id_vert)%n_sub_faces_neigh + j) = 2*(ce - cp)
    end do

    L = 0.0_DOUBLE
    if( boundary_2d ) then
      L(1) = (mesh%vert(id_vert)%volume/error_2d_h)**(-1.0_DOUBLE/2.0_DOUBLE)
      L(4) = (mesh%vert(id_vert)%volume/error_2d_h)**(-1.0_DOUBLE/2.0_DOUBLE)
      L(6) = 1.0_DOUBLE/error_2d_h
    else
      L(1) = mesh%vert(id_vert)%volume**(-1.0_DOUBLE/3.0_DOUBLE)
      L(4) = mesh%vert(id_vert)%volume**(-1.0_DOUBLE/3.0_DOUBLE)
      L(6) = mesh%vert(id_vert)%volume**(-1.0_DOUBLE/3.0_DOUBLE)
    end if

    maxiter = 100
    iter = 0
    delta_L = 1.0_DOUBLE
    do while (maxval(abs(delta_L)) > 1e-4_DOUBLE .and. iter <= maxiter)
      iter = iter + 1
      do j = 1, nc
        Ltx(1) = L(1)*xc(1, j) + L(2)*xc(2, j) + L(3)*xc(3, j)
        Ltx(2) =                  L(4)*xc(2, j) + L(5)*xc(3, j)
        Ltx(3) =                                   L(6)*xc(3, j)
        JL(j, 1) = 2*xc(1, j)*Ltx(1)
        JL(j, 2) = 2*xc(2, j)*Ltx(1)
        JL(j, 3) = 2*xc(3, j)*Ltx(1)
        JL(j, 4) = 2*xc(2, j)*Ltx(2)
        JL(j, 5) = 2*xc(3, j)*Ltx(2)
        JL(j, 6) = 2*xc(3, j)*Ltx(3)
        fL(j) = dot_product(Ltx, Ltx) - 1.0_DOUBLE
      end do
      LLt = matmul(transpose(JL), JL)
      call pseudo_inverse_inplace_lapack(6, LLt)
      delta_L = -matmul(LLt, matmul(transpose(JL), fL))
      L = L + delta_L
    end do

    if( iter >= maxiter ) then
      print*, id_vert, "DID NOT CONVERGE ELLIP", maxval(abs(delta_L))
    end if

    mat_h_p = 0.0_DOUBLE
    mat_h_p(1, 1) = L(1)
    mat_h_p(1, 2) = L(2)
    mat_h_p(1, 3) = L(3)
    mat_h_p(2, 2) = L(4)
    mat_h_p(2, 3) = L(5)
    mat_h_p(3, 3) = L(6)

    mat_h_p = matmul(transpose(mat_h_p), mat_h_p)

    call spectral_decomposition(mat_h_p, Q, Lambda)
    Lambda = 1.0_DOUBLE / sqrt(Lambda)
    mat_h_p = matmul(Q, matmul(diag(Lambda), transpose(Q)))
  end function compute_ellip

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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
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
      vbar = (- (pr -pl))/(lambda_l+lambda_r)
      rhs = rhs + mesh%sub_face(id_sub_face)%area&
      *(lambda_l + lambda_r)*vbar*norm

      !Div(v)
      invlamb = (1.0_DOUBLE/lambda_l + 1.0_DOUBLE/lambda_r)
      pbar = (- (vnr - vnl))/invlamb
      div_v = div_v + mesh%sub_face(id_sub_face)%area*invlamb*pbar
      denomsum = denomsum + mesh%sub_face(id_sub_face)%area*invlamb
    end do

    call inv_lapack(3, mat)
    grad_p = matmul(mat, rhs) * rho_p * a_p

    div_v = div_v/(denomsum*a_p*rho_p)

    corr_div_v = max(0.0_DOUBLE, min(1.0_DOUBLE, abs(div_v)/a_p))*a_p
    corr_grad_p = max(0.0_DOUBLE, min(1.0_DOUBLE, norm2(grad_p)/a_p**2))*a_p
    corr = max(corr_div_v, corr_grad_p)
  end subroutine compute_corr2

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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) &
            - mesh%sub_face(id_sub_face)%area*fplus
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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_LPF

  subroutine compute_rhs_around_vert_WIP(mesh, sol, grad, &
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
    real(kind=DOUBLE), dimension(nsen), intent(inout) :: &
      sum_lambda_vert
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: &
      flux_sum_vert
    integer(kind=ENTIER), intent(in) :: id_vert
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, id_sub_face, id_face
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE), dimension(3) :: norm

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
      sol_r = primit_to_conserv(sol_w_r)
      grad_sol_p = grad_sol_p + mesh%sub_face(id_sub_face)%area&
        *tensor_product(sol_r - sol_l, mesh%sub_face(id_sub_face)%norm)
      sum_area = sum_area + mesh%sub_face(id_sub_face)%area
    end do
    grad_sol_p = grad_sol_p / (0.5_DOUBLE*sum_area)

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

    call compute_corr2(mesh, id_vert, sol, grad, corr, second_order)

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

      lambda_l = rhol*al
      lambda_r = rhor*ar
      vbar = (lambda_l*vnl + lambda_r*vnr - (pr-pl))/(lambda_l+lambda_r)
      pbar = (lambda_r * pl + lambda_l * pr - lambda_l*lambda_r*(vnr-vnl))/(lambda_l+lambda_r)

      !vbar = 0.5_DOUBLE*(vnl+vnr) - 0.5_DOUBLE/(rhom*am) * (pr - pl)
      !pbar = 0.5_DOUBLE*(pl+pr) - 0.5_DOUBLE*rhom*am*(vnr-vnl)
      !ff_adv = vbar*sol_m - 0.5_DOUBLE*max(abs(vnl),abs(vnr))*(sol_r-sol_l)
      !ff_adv = 0.5_DOUBLE*(vnl*sol_l+vnr*sol_r) &
      !  - 0.5_DOUBLE*max(abs(vnl),abs(vnr))*(sol_r-sol_l)

      !ff_adv = vbar*sol_m - 0.5_DOUBLE*(abs(vbar)+min(1.0_DOUBLE, machm)*am)*(sol_r-sol_l)
      ff_adv = vbar*sol_m - 0.5_DOUBLE*(abs(vbar)+corr)*(sol_r-sol_l)

      ff_lag(1) = 0.0_DOUBLE
      !ff_lag(2:4) = pbar * norm
      !ff_lag(5) = pbar * vbar
      ff_lag(2:4) = pp * norm
      ff_lag(5) = pp * vbar

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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*lambda_lts
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) - mesh%sub_face(id_sub_face)%area*fplus
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_WIP

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

    integer(kind=ENTIER) :: j, d, id_sub_face
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc, id_sub_elem, id_elem
    real(kind=DOUBLE), dimension(3) :: norm
    real(kind=DOUBLE) :: pp_node, corrp, cp_node, lambda_lts
    real(kind=DOUBLE) :: rhol, rhor, vnl, vnr, pl, pr, cL, cR
    real(kind=DOUBLE) :: rhom, am, up, pp_face, smax
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, sol_l, sol_r, sol_m
    real(kind=DOUBLE), dimension(5) :: ff_1d, fminus, fplus, Qp, Fp_n
    real(kind=DOUBLE), dimension(3) :: vp_avg, vp_node
    real(kind=DOUBLE), dimension(5, 3) :: Fp, gradQ
    real(kind=DOUBLE) :: cfweight, sum_area

    cfweight = 1.0_DOUBLE / 3.0_DOUBLE
    !cfweight = 0.0_DOUBLE
    rse_loc  = 0

    call compute_nodal_velocity_LVP(mesh, id_vert, sol, grad, vp_node, second_order)
    call compute_nodal_pressure_LPP(mesh, id_vert, sol, grad, pp_node, second_order)
    call compute_corr2(mesh, id_vert, sol, grad, corrp, second_order)

    vp(:, id_vert) = vp_node

    ! --- volume-weighted averages over sub-elements ---
    Qp      = 0.0_DOUBLE
    vp_avg  = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      Qp = Qp + mesh%sub_elem(id_sub_elem)%volume * sol(:, id_elem)
      vp_avg = vp_avg &
        + mesh%sub_elem(id_sub_elem)%volume * sol(2:4, id_elem) / sol(1, id_elem)
    end do
    Qp = Qp / mesh%vert(id_vert)%volume
    vp_avg = vp_avg / mesh%vert(id_vert)%volume
    cp_node = sound_speed_w(conserv_to_primit(Qp))

    ! --- gradient of Q via jump across sub-faces ---
    gradQ    = 0.0_DOUBLE
    sum_area = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      le   = mesh%sub_face(id_sub_face)%left_elem_neigh
      re   = mesh%sub_face(id_sub_face)%right_elem_neigh
      norm = mesh%sub_face(id_sub_face)%norm

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_l, sol_w_r)

      sol_l = primit_to_conserv(sol_w_l)
      sol_r = primit_to_conserv(sol_w_r)

      if( .not. (boundary_2d .and. abs(norm(3)) > 1e-8_DOUBLE) ) then
        gradQ    = gradQ &
          + mesh%sub_face(id_sub_face)%area * tensor_product(sol_r - sol_l, norm)
        sum_area = sum_area + mesh%sub_face(id_sub_face)%area
      end if
    end do
    gradQ = gradQ / (sum_area + 1e-12_DOUBLE)

    ! --- build nodal flux tensor Fp(5,3) ---
    ! advection: vp_node(d)*Qp - 1/(nDim+1)*(|vp_avg(d)| + corrp)*gradQ(:,d)
    ! Lagrangian: pp_node on momentum diagonal and energy column
    !Fp = tensor_product(Qp, vp_node) &
    !    - 0.5_DOUBLE * matmul(diag(abs(vp_node)) + corrp*eye3, gradQ)
    Fp = tensor_product(Qp, vp_node) &
        - 0.5_DOUBLE * (norm2(vp_node) + corrp)*gradQ
    Fp(2:4, :) = Fp(2:4, :) + pp_node * eye3
    Fp(5, :) = Fp(5, :) + vp_node * pp_node

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

      rhol = sol_w_l(1);  rhoR = sol_w_r(1)
      vnl  = dot_product(sol_w_l(2:4), norm)
      vnr  = dot_product(sol_w_r(2:4), norm)
      pl   = sol_w_l(5);  pr  = sol_w_r(5)
      cL   = sound_speed_w(sol_w_l)
      cR   = sound_speed_w(sol_w_r)

      rhom   = 0.5_DOUBLE * (rhol + rhoR)
      am     = max(cL, cR)
      up     = 0.5_DOUBLE * (vnl + vnr) - 0.5_DOUBLE / (rhom * am) * (pr - pl)
      pp_face= 0.5_DOUBLE * (pl  + pr ) - 0.5_DOUBLE * rhom * am * (vnr - vnl)
      smax   = abs(up) + corrp * am

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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area * lambda_lts
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) &
            - mesh%sub_face(id_sub_face)%area * fplus
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_usi3d

  subroutine compute_rhs_around_vert_LPC(mesh, sol, grad, &
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

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) &
            - mesh%sub_face(id_sub_face)%area*fplus
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert_LPC
end module ns_euler_zb_module