module ns_mesh_metric_module
  use precision_module
  use mesh_module
  implicit none

contains
  subroutine allocate_neigh(mesh, id_vert, n_min, n_neigh, neigh)
    use sort_module, only: add_sort_unique_int
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert, n_min
    integer(kind=ENTIER), intent(inout) :: n_neigh
    integer(kind=ENTIER), dimension(:), allocatable, intent(inout) :: neigh

    integer(kind=ENTIER) :: i
    integer(kind=ENTIER) :: n_pot_vert
    integer(kind=ENTIER), dimension(:), allocatable :: pot_vert

    call add_sort_unique_int(neigh, n_neigh, mesh%vert(id_vert)%elem_neigh)

    do while (n_neigh < n_min)
      do i = 1, n_neigh
        call add_sort_unique_int(pot_vert, n_pot_vert, mesh%elem(neigh(i))%vert)
      end do
      do i = 1, n_pot_vert
        call add_sort_unique_int(neigh, n_neigh, mesh%vert(pot_vert(i))%elem_neigh)
      end do
    end do
  end subroutine allocate_neigh

  subroutine compute_hessian(mesh, id_vert, phi, hessian)
    use linear_solver_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(3, 3), intent(inout) :: hessian

    integer(kind=ENTIER) :: i
    real(kind=DOUBLE), dimension(10) :: v
    real(kind=DOUBLE), dimension(10, 10) :: mat
    real(kind=DOUBLE), dimension(10) :: rhs

    integer(kind=ENTIER) :: n_neigh
    integer(kind=ENTIER), dimension(:), allocatable :: neigh
    real(kind=DOUBLE), dimension(3) :: d
    real(kind=DOUBLE) :: weight

    !call allocate_neigh(mesh, id_vert, 20, n_neigh, neigh)

    !mat = 0.0_DOUBLE
    !rhs = 0.0_DOUBLE
    !do i = 1, n_neigh
    !  d = mesh%elem(neigh(i))%coord - mesh%vert(id_vert)%coord
    !  weight = norm2(d)
    !  v(1) = 1.0_DOUBLE
    !  v(2) = d(1)
    !  v(3) = d(2)
    !  v(4) = d(3)
    !  v(5) = d(1)*d(1)
    !  v(6) = 2*d(1)*d(2)
    !  v(7) = 2*d(1)*d(3)
    !  v(8) = d(2)*d(2)
    !  v(9) = 2*d(2)*d(3)
    !  v(10) = d(3)*d(3)

    !  mat = mat + (1.0_DOUBLE/weight)**2*tensor_product(v, v)
    !  rhs = rhs + (1.0_DOUBLE/weight)**2*phi(neigh(i))*v
    !end do

    !call inv_lapack(10, mat)
    !v = matmul(mat, rhs)

    v = 0.0_DOUBLE

    hessian(1, 1) = v(5)
    hessian(1, 2) = v(6)
    hessian(1, 3) = v(7)
    hessian(2, 2) = v(8)
    hessian(2, 3) = v(9)
    hessian(3, 3) = v(10)

    hessian(2, 1) = hessian(1, 2)
    hessian(3, 1) = hessian(1, 3)
    hessian(3, 2) = hessian(2, 3)
  end subroutine compute_hessian

  subroutine compute_cell_metric(mesh, sol, cell_metric)
    use ns_global_data_module, only: coeffs_surf
    use ns_euler_primitives_module, only: conserv_to_primit, mach_w, temp_w
    use linear_solver_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems) :: sol
    real(kind=DOUBLE), dimension(3, 3, mesh%n_vert) :: cell_metric

    integer(kind=ENTIER) :: i
    real(kind=DOUBLE), dimension(:), allocatable :: sol_mach, sol_temp
    real(kind=DOUBLE), dimension(:, :, :), allocatable :: hessian
    real(kind=DOUBLE), dimension(5) :: w
    real(kind=DOUBLE), dimension(3, 3) :: Q
    real(kind=DOUBLE), dimension(3) :: Lambda

    allocate(sol_mach(mesh%n_elems))
    allocate(sol_temp(mesh%n_elems))
    do i = 1, mesh%n_elems
      w = conserv_to_primit(sol(:, i))
      sol_mach(i) = mach_w(w)
      sol_temp(i) = temp_w(w)
    end do

    allocate(hessian(3, 3, mesh%n_vert))
    do i = 1, mesh%n_vert
      call compute_hessian(mesh, i, sol_temp, hessian(:, :, i))
    end do

    do i = 1, mesh%n_vert
      call spectral_decomposition(hessian(:, :, i), Q, Lambda)
      cell_metric(:, :, i) = matmul(Q, matmul(diag(abs(Lambda)), transpose(Q)))
    end do
  end subroutine compute_cell_metric
end module ns_mesh_metric_module