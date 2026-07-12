module ns_euler_recon_module
  use precision_module
  use mesh_module
  use ns_euler_primitives_module
  implicit none

contains
  subroutine compute_cell_grad_from_nodal_grad(mesh, sol, grad, all_nodal_grad, method)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, 5, mesh%n_elems), intent(inout) :: grad
    real(kind=DOUBLE), dimension(3, 5, mesh%n_vert), intent(inout) :: all_nodal_grad
    integer(kind=ENTIER), intent(in) :: method

    integer(kind=ENTIER) :: i

    do i=1, mesh%n_vert
      call compute_nodal_grad(mesh, i, sol, all_nodal_grad(:, :, i))
    end do

    if( method == 0 .or. method == 1 ) then
      grad = 0.0_DOUBLE
    else if( method == 2 ) then
      do i=1, mesh%n_elems
        call compute_cell_grad_method_2(mesh, all_nodal_grad, i, grad(:, :, i))
        call physical_slope_limiter(mesh, sol, i, grad)
        call relaxed_maximum_principle(mesh, sol, i, grad)
      end do
    else if( method == 3 ) then
      do i=1, mesh%n_elems
        call compute_cell_grad_method_3(mesh, all_nodal_grad, i, grad(:, :, i))
        call physical_slope_limiter(mesh, sol, i, grad)
        call relaxed_maximum_principle(mesh, sol, i, grad)
      end do
    else if( method == 4 ) then
      do i=1, mesh%n_elems
        call compute_cell_grad_method_4(mesh, all_nodal_grad, i, grad(:, :, i))
        call physical_slope_limiter(mesh, sol, i, grad)
        call relaxed_maximum_principle(mesh, sol, i, grad)
      end do
    else if( method == 5 ) then
      do i=1, mesh%n_elems
        call compute_cell_grad_method_5(mesh, all_nodal_grad, i, grad(:, :, i))
        call physical_slope_limiter(mesh, sol, i, grad)
        call relaxed_maximum_principle(mesh, sol, i, grad)
      end do
    else if( method == 6 ) then
      do i=1, mesh%n_elems
        call compute_cell_grad_method_6(mesh, all_nodal_grad, i, grad(:, :, i))
        call physical_slope_limiter(mesh, sol, i, grad)
        call relaxed_maximum_principle(mesh, sol, i, grad)
      end do
    else
      print*, "Method for second order unknown !"
      error stop
    end if
  end subroutine compute_cell_grad_from_nodal_grad

  subroutine relaxed_maximum_principle(mesh, sol, i, grad)
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE), dimension(3, 5, mesh%n_elems), intent(inout) :: grad

    integer(kind=ENTIER) :: j, id_vert, k, id_elem
    real(kind=DOUBLE) :: omega
    real(kind=DOUBLE), dimension(5) :: sol_w, delta_sol_w

    real(kind=DOUBLE) :: rho_min, rho_max, p_min, p_max

    real(kind=DOUBLE), parameter :: tol = 1e-8_DOUBLE

    integer(kind=ENTIER) :: me, mpi_ierr

    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    !Compute min max
    rho_min = 1e100_DOUBLE
    rho_max = 0.0_DOUBLE
    p_min = 1e100_DOUBLE
    p_max = 0.0_DOUBLE
    do j=1, mesh%elem(i)%n_vert
      id_vert = mesh%elem(i)%vert(j)
      do k=1, mesh%vert(id_vert)%n_elems_neigh
        id_elem = mesh%vert(id_vert)%elem_neigh(k)
        sol_w = conserv_to_primit(sol(:, id_elem))
        p_min = min(p_min, sol_w(5))
        p_max = max(p_max, sol_w(5))
        rho_min = min(rho_min, sol_w(1))
        rho_max = max(rho_max, sol_w(1))
      end do
    end do

    !Relax
    p_min = (1.0_DOUBLE-1e-2_DOUBLE)*p_min
    rho_min = (1.0_DOUBLE-1e-2_DOUBLE)*rho_min
    p_max = (1.0_DOUBLE+1e-2_DOUBLE)*p_max
    rho_max = (1.0_DOUBLE+1e-2_DOUBLE)*rho_max

    sol_w = conserv_to_primit(sol(:, i))
    omega = 1.0_DOUBLE
    do j=1, mesh%elem(i)%n_vert
      id_vert = mesh%elem(i)%vert(j)
      delta_sol_w = matmul(transpose(grad(:, :, i)), &
        mesh%vert(id_vert)%coord - mesh%elem(i)%coord)
      if( delta_sol_w(1) < 0.0_DOUBLE ) then
        omega = min(omega, &
          max(0.0_DOUBLE, (sol_w(1) - rho_min)/(tol - delta_sol_w(1))))
      else 
        omega = min(omega, &
          max(0.0_DOUBLE, -(sol_w(1) - rho_max)/(tol + delta_sol_w(1))))
      end if
      if( delta_sol_w(5) < 0.0_DOUBLE ) then
        omega = min(omega, &
          max(0.0_DOUBLE, (sol_w(5) - p_min)/(tol - delta_sol_w(5))))
      else
        omega = min(omega, &
          max(0.0_DOUBLE, -(sol_w(5) - p_max)/(tol + delta_sol_w(5))))
      end if
    end do
    grad(:, :, i) = omega * grad(:, :, i)
  end subroutine relaxed_maximum_principle

  subroutine physical_slope_limiter(mesh, sol, i, grad)
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE), dimension(3, 5, mesh%n_elems), intent(inout) :: grad

    integer(kind=ENTIER) :: j, id_vert
    real(kind=DOUBLE) :: omega
    real(kind=DOUBLE), dimension(5) :: sol_w, delta_sol_w

    real(kind=DOUBLE), parameter :: tol = 1e-8_DOUBLE

    integer(kind=ENTIER) :: me, mpi_ierr

    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    sol_w = conserv_to_primit(sol(:, i))
    omega = 1.0_DOUBLE
    do j=1, mesh%elem(i)%n_vert
      id_vert = mesh%elem(i)%vert(j)
      delta_sol_w = matmul(transpose(grad(:, :, i)), &
        mesh%vert(id_vert)%coord - mesh%elem(i)%coord)
      if( delta_sol_w(1) < 0.0_DOUBLE ) then
        omega = min(omega, &
          max(0.0_DOUBLE, sol_w(1)/(tol - delta_sol_w(1))))
      end if
      if( delta_sol_w(5) < 0.0_DOUBLE ) then
        omega = min(omega, &
          max(0.0_DOUBLE, sol_w(5)/(tol - delta_sol_w(5))))
      end if
    end do
    grad(:, :, i) = omega * grad(:, :, i)
  end subroutine physical_slope_limiter

  subroutine compute_cell_grad_method_2(mesh, all_nodal_grad, id_elem, grad)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(3, 5, mesh%n_vert), intent(in) :: all_nodal_grad
    integer(kind=ENTIER), intent(in) :: id_elem
    real(kind=DOUBLE), dimension(3, 5), intent(inout) :: grad

    integer(kind=ENTIER) :: j, id_vert
    real(kind=DOUBLE) :: dual_vol_sum

    grad(:, :) = 0.0_DOUBLE
    dual_vol_sum = 0.0_DOUBLE
    do j=1, mesh%elem(id_elem)%n_vert
      id_vert = mesh%elem(id_elem)%vert(j)
      grad(:, :) = grad(:, :) &
        + mesh%vert(id_vert)%volume*all_nodal_grad(:, :, id_vert)
      dual_vol_sum = dual_vol_sum &
        + mesh%vert(id_vert)%volume
    end do
    grad = grad/dual_vol_sum
  end subroutine compute_cell_grad_method_2

  subroutine compute_cell_grad_method_3(mesh, all_nodal_grad, id_elem, grad)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(3, 5, mesh%n_vert), intent(in) :: all_nodal_grad
    integer(kind=ENTIER), intent(in) :: id_elem
    real(kind=DOUBLE), dimension(3, 5), intent(inout) :: grad

    integer(kind=ENTIER) :: j, id_vert, var

    id_vert = mesh%elem(id_elem)%vert(1)
    grad(:, :) = all_nodal_grad(:, :, id_vert)
    do j=2, mesh%elem(id_elem)%n_vert
      id_vert = mesh%elem(id_elem)%vert(j)
      do var=1, 5
        if(OI(all_nodal_grad(:, var, id_vert)) < OI(grad(:, var))) then
          grad(:, var) = all_nodal_grad(:, var, id_vert)
        end if
      end do
    end do
  end subroutine compute_cell_grad_method_3

  subroutine compute_cell_grad_method_4(mesh, all_nodal_grad, id_elem, grad)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(3, 5, mesh%n_vert), intent(in) :: all_nodal_grad
    integer(kind=ENTIER), intent(in) :: id_elem
    real(kind=DOUBLE), dimension(3, 5), intent(inout) :: grad

    real(kind=DOUBLE), parameter :: eps=1e-3_DOUBLE

    integer(kind=ENTIER) :: j, id_vert, var
    real(kind=DOUBLE), dimension(5) :: weight_sum

    !Non-linear combination with alpha=1
    grad(:, :) = 0.0_DOUBLE
    weight_sum(:) = 0.0_DOUBLE
    do j=1, mesh%elem(id_elem)%n_vert
      id_vert = mesh%elem(id_elem)%vert(j)
      do var=1, 5
        grad(:, var) = grad(:, var) &
          + mesh%vert(id_vert)%volume/(eps+OI(all_nodal_grad(:, var, id_vert))) &
          * all_nodal_grad(:, var, id_vert)
        weight_sum(var) = weight_sum(var) &
          + mesh%vert(id_vert)%volume/(eps+OI(all_nodal_grad(:, var, id_vert)))
      end do
    end do

    do var=1, 5
      grad(:, var) = grad(:, var)/weight_sum(var)
    end do
  end subroutine compute_cell_grad_method_4

  subroutine compute_cell_grad_method_5(mesh, all_nodal_grad, id_elem, grad)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(3, 5, mesh%n_vert), intent(in) :: all_nodal_grad
    integer(kind=ENTIER), intent(in) :: id_elem
    real(kind=DOUBLE), dimension(3, 5), intent(inout) :: grad

    real(kind=DOUBLE), parameter :: eps=1e-3_DOUBLE

    integer(kind=ENTIER) :: j, id_vert, var
    real(kind=DOUBLE), dimension(5) :: weight_sum

    !Non-linear combination with alpha=1
    grad(:, :) = 0.0_DOUBLE
    weight_sum(:) = 0.0_DOUBLE
    do j=1, mesh%elem(id_elem)%n_vert
      id_vert = mesh%elem(id_elem)%vert(j)
      do var=1, 5
        grad(:, var) = grad(:, var) &
          + mesh%vert(id_vert)%volume/(eps+OI(all_nodal_grad(:, var, id_vert)))**2 &
          * all_nodal_grad(:, var, id_vert)
        weight_sum(var) = weight_sum(var) &
          + mesh%vert(id_vert)%volume/(eps+OI(all_nodal_grad(:, var, id_vert)))**2
      end do
    end do

    do var=1, 5
      grad(:, var) = grad(:, var)/weight_sum(var)
    end do
  end subroutine compute_cell_grad_method_5

  subroutine compute_cell_grad_method_6(mesh, all_nodal_grad, id_elem, grad)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(3, 5, mesh%n_vert), intent(in) :: all_nodal_grad
    integer(kind=ENTIER), intent(in) :: id_elem
    real(kind=DOUBLE), dimension(3, 5), intent(inout) :: grad

    real(kind=DOUBLE), parameter :: eps=1e-3_DOUBLE

    integer(kind=ENTIER) :: j, id_vert, var
    real(kind=DOUBLE), dimension(5) :: weight_sum
    real(kind=DOUBLE) :: central_weight_sum
    real(kind=DOUBLE), dimension(3, 5) :: central_grad

    !Non-linear combination with alpha=1
    grad(:, :) = 0.0_DOUBLE
    weight_sum(:) = 0.0_DOUBLE
    central_grad(:, :) = 0.0_DOUBLE
    central_weight_sum = 0.0_DOUBLE
    do j=1, mesh%elem(id_elem)%n_vert
      id_vert = mesh%elem(id_elem)%vert(j)
      do var=1, 5
        grad(:, var) = grad(:, var) &
          + mesh%vert(id_vert)%volume/(eps+OI(all_nodal_grad(:, var, id_vert)))**2 &
          * all_nodal_grad(:, var, id_vert)
        weight_sum(var) = weight_sum(var) &
          + mesh%vert(id_vert)%volume/(eps+OI(all_nodal_grad(:, var, id_vert)))**2
      end do
      central_grad(:, :) = central_grad(:, :) &
        + mesh%vert(id_vert)%volume*all_nodal_grad(:, :, id_vert)
      central_weight_sum = central_weight_sum &
        + mesh%vert(id_vert)%volume
    end do

    central_grad = central_grad/central_weight_sum
    grad(:, :) = grad(:, :) &
      + central_weight_sum/(eps+OI(central_grad(:, :)))**2 &
      * central_grad(:, :)
    weight_sum(:) = weight_sum(:) &
      + central_weight_sum/(eps+OI(central_grad(:, :)))**2

    do var=1, 5
      grad(:, var) = grad(:, var)/weight_sum(var)
    end do
  end subroutine compute_cell_grad_method_6

  pure function OI(grad)
    implicit none

    real(kind=DOUBLE), dimension(3), intent(in) :: grad
    real(kind=DOUBLE) :: OI

    OI = norm2(grad(:))
  end function OI

  subroutine compute_right_state(mesh, id_sub_face, re, sol_l, sol_r)
    use ns_global_data_module, only: bc_val, t, &
      bc_euler_id, BC_EULER_WALL, BC_EULER_ADHERENCE_WALL, BC_EULER_FREESTREAM, &
      BC_EULER_OUTFLOWSUPERSONIC, BC_EULER_INFLOW_POND, &
      BC_EULER_INOUT_DOUBLE_MACH, BC_EULER_DOUBLE_MACH_BOTTOM, &
      BC_EULER_POTENTIAL_FLOW_2D
    use ns_euler_exact_sol_module, only: sol_potential_flow_2d
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_sub_face, re
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_l
    real(kind=DOUBLE), dimension(5), intent(inout) :: sol_r

    integer(kind=ENTIER) :: id_face
    real(kind=DOUBLE), dimension(3) :: face_coord
    real(kind=DOUBLE), dimension(5) :: w_pf

    if (re > 0) then
      print*, "re > 0, bad use of compute right sate !"
      error stop
    else if (re == 0) then !Default option wall
      sol_r(:) = sol_l
      sol_r(2:4) = sol_r(2:4) &
        - 2.0_DOUBLE*dot_product(sol_r(2:4), mesh%sub_face(id_sub_face)%norm)*mesh%sub_face(id_sub_face)%norm
    else !!Boundary
      select case (bc_euler_id(-re))
      case (BC_EULER_OUTFLOWSUPERSONIC)
        sol_r = sol_l
      case (BC_EULER_FREESTREAM)
        sol_r = primit_to_conserv(bc_val(:, -re))
      case (BC_EULER_INFLOW_POND)
        sol_r = primit_to_conserv(bc_val(:, -re))
        if (dot_product(sol_r(2:4), mesh%sub_face(id_sub_face)%norm) > 0.0_DOUBLE) then
          sol_r = sol_l
        end if
      case (BC_EULER_WALL)
        sol_r(:) = sol_l
        sol_r(2:4) = sol_r(2:4) &
          - 2.0_DOUBLE*dot_product(sol_r(2:4), mesh%sub_face(id_sub_face)%norm)*mesh%sub_face(id_sub_face)%norm
      case (BC_EULER_ADHERENCE_WALL)
        sol_r(:) = sol_l
        sol_r(2:4) = 2.0_DOUBLE*sol_l(1)*bc_val(2:4, -re) - sol_l(2:4)
        sol_r(5) = sol_l(5) - 0.5_DOUBLE*sol_l(1)*norm2(sol_l(2:4)/sol_l(1))**2 &
          + 0.5_DOUBLE*sol_l(1)*norm2(sol_r(2:4)/sol_r(1))**2
      case (BC_EULER_INOUT_DOUBLE_MACH)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        face_coord = mesh%face(id_face)%coord
        if (face_coord(2) > 1.732_DOUBLE*(face_coord(1) - 0.1667_DOUBLE - 10.0_DOUBLE*t)) then
          sol_r = primit_to_conserv((/8.0_DOUBLE, 7.145_DOUBLE, -4.125_DOUBLE, 0.0_DOUBLE, 116.5_DOUBLE/))
        else
          sol_r = primit_to_conserv((/1.4_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE/))
        end if
      case (BC_EULER_DOUBLE_MACH_BOTTOM)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        face_coord = mesh%face(id_face)%coord
        if (face_coord(1) < 0.1667_DOUBLE) then
          sol_r = primit_to_conserv((/8.0_DOUBLE, 7.145_DOUBLE, -4.125_DOUBLE, 0.0_DOUBLE, 116.5_DOUBLE/))
        else
          sol_r(:) = sol_l
          sol_r(2:4) = sol_r(2:4) &
            - 2.0_DOUBLE*dot_product(sol_r(2:4), mesh%sub_face(id_sub_face)%norm)*mesh%sub_face(id_sub_face)%norm
        end if
      case (BC_EULER_POTENTIAL_FLOW_2D)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        face_coord = mesh%face(id_face)%coord
        call sol_potential_flow_2d(face_coord, w_pf)
        sol_r = primit_to_conserv(w_pf)
      case default
        print *, "BC TYPE NOT RECOGNIZED !"
        error stop
      end select
    end if
  end subroutine compute_right_state

  pure function is_wall(re)
    use ns_global_data_module, only: bc_euler_id, BC_EULER_WALL
    implicit none

    integer(kind=ENTIER), intent(in) :: re
    logical :: is_wall

    if (re > 0) then
      is_wall = .FALSE.
    else if (re == 0) then
      is_wall = .TRUE.
    else
      is_wall = bc_euler_id(-re) == BC_EULER_WALL
    end if
  end function is_wall

  subroutine reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
      second_order, sol_w_l, sol_w_r)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    integer(kind=ENTIER), intent(in) :: id_vert, id_sub_face, le, re
    logical, intent(in) :: second_order
    real(kind=DOUBLE), dimension(5), intent(out) :: sol_w_l, sol_w_r

    real(kind=DOUBLE), dimension(5) :: sol_ghost

    sol_w_l = conserv_to_primit(sol(:, le))
    if (second_order) sol_w_l = sol_w_l + &
      matmul(transpose(grad(:, :, le)), mesh%vert(id_vert)%coord - mesh%elem(le)%coord)
    if (re > 0) then
      sol_w_r = conserv_to_primit(sol(:, re))
      if (second_order) sol_w_r = sol_w_r + &
        matmul(transpose(grad(:, :, re)), mesh%vert(id_vert)%coord - mesh%elem(re)%coord)
    else
      call compute_right_state(mesh, id_sub_face, re, primit_to_conserv(sol_w_l), sol_ghost)
      sol_w_r = conserv_to_primit(sol_ghost)
    end if
  end subroutine reconstruct_lr_w

  subroutine compute_nodal_grad(mesh, id_vert, sol, nodal_grad)
    use linear_solver_module, only: tensor_product, &
      pseudo_inverse_inplace_lapack, print_mat
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, 5), intent(inout) :: nodal_grad

    integer(kind=ENTIER) :: i, id_sub_face, id_face, le, re
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(3, 3) :: mat

    nodal_grad = 0.0_DOUBLE
    mat = 0.0_DOUBLE
    do i=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(i)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      le = mesh%face(id_face)%left_neigh
      re = mesh%face(id_face)%right_neigh
      if( re > 0 ) then
        sol_w_l = conserv_to_primit(sol(:, le))
        sol_w_r = conserv_to_primit(sol(:, re))
        nodal_grad = nodal_grad &
          + mesh%sub_face(id_sub_face)%area&
          * tensor_product(mesh%sub_face(id_sub_face)%norm, &
          sol_w_r - sol_w_l)
        mat = mat &
          + mesh%sub_face(id_sub_face)%area&
          * tensor_product(mesh%sub_face(id_sub_face)%norm, &
          mesh%elem(re)%coord - mesh%elem(le)%coord)
      end if
    end do
    call pseudo_inverse_inplace_lapack(3, mat)
    nodal_grad = matmul(mat, nodal_grad)
    if( mesh%vert(id_vert)%is_bound ) nodal_grad = 0.0_DOUBLE
  end subroutine compute_nodal_grad

  subroutine compute_vert_div_rot(mesh, id_vert, sol, div, rot)
    use linear_solver_module, only: cross_product
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3), intent(inout) :: rot
    real(kind=DOUBLE), intent(inout) :: div

    integer(kind=ENTIER) :: i, id_sub_face, id_face, le, re
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r

    div = 0.0_DOUBLE
    rot = 0.0_DOUBLE
    do i=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(i)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      le = mesh%face(id_face)%left_neigh
      re = mesh%face(id_face)%right_neigh
      if( re > 0 ) then
        sol_w_l = conserv_to_primit(sol(:, le))
        sol_w_r = conserv_to_primit(sol(:, re))
        div = div &
          + mesh%sub_face(id_sub_face)%area&
          * dot_product(mesh%sub_face(id_sub_face)%norm, &
          sol_w_r(2:4) - sol_w_l(2:4))
        rot = rot &
          + mesh%sub_face(id_sub_face)%area&
          * cross_product(mesh%sub_face(id_sub_face)%norm, &
          sol_w_r(2:4) - sol_w_l(2:4))
      end if
    end do

    div = div / mesh%vert(i)%volume
    rot = rot / mesh%vert(i)%volume

    if( mesh%vert(id_vert)%is_bound ) then
      div = 0.0_DOUBLE
      rot = 0.0_DOUBLE
    end if
  end subroutine compute_vert_div_rot

  function omega_ducros(mesh, sol, id_elem) result(omega)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_elem
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE) :: omega

    integer(kind=ENTIER) :: j, id_sub_elem, id_vert
    real(kind=DOUBLE) :: div, a, b
    real(kind=DOUBLE), dimension(3) :: rot

    a = 0.0_DOUBLE
    b = 0.0_DOUBLE
    do j=1, mesh%elem(id_elem)%n_sub_elems
      id_sub_elem = mesh%elem(id_elem)%sub_elem(j)
      id_vert = mesh%sub_elem(id_sub_elem)%mesh_vert
      call compute_vert_div_rot(mesh, id_vert, sol, div, rot)
      a = a + mesh%sub_elem(id_sub_elem)%volume * abs(div)**2
      b = b + mesh%sub_elem(id_sub_elem)%volume &
        * (abs(div)**2 + 5*norm2(rot)**2)
    end do
    omega = a/(b+1e-4_DOUBLE)
  end function omega_ducros

  subroutine flux_limiting(mesh, sol, rhs, rhs_o2, sum_lambda, cfl)
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: rhs_o2
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: sum_lambda
    real(kind=DOUBLE), intent(in) :: cfl

    integer(kind=ENTIER) :: i
    integer :: mpi_ierr
    real(kind=DOUBLE) :: omega, dt

    dt = 1e10
    do i = 1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        dt = min(dt, mesh%elem(i)%volume/sum_lambda(i))
      end if
    end do
    dt = cfl * dt

    call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_DOUBLE, &
      MPI_MIN, MPI_COMM_WORLD, mpi_ierr)

    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        omega = max_omega_for_positivity(sol(:, i) &
          + rhs(:, i)/(mesh%elem(i)%volume/dt + sum_lambda(i)), &
          (rhs_o2(:, i) - rhs(:, i))/(mesh%elem(i)%volume/dt + sum_lambda(i)))
        rhs(:, i) = rhs(:, i) + min(1.0_DOUBLE, omega) * (rhs_o2(:, i) - rhs(:, i))
      end if
    end do

    !do i=1, mesh%n_elems
    !  if( .not. mesh%elem(i)%is_ghost ) then
    !    rhs(:, i) = rhs(:, i) &
    !      + (1.0_DOUBLE - omega_ducros(mesh, sol, i)) * (rhs_o2(:, i) - rhs(:, i))
    !  end if
    !end do
  end subroutine flux_limiting
end module ns_euler_recon_module
