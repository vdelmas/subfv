module sparse_bcsr_linear_module
  use linear_solver_module
  use precision_module
  implicit none

  type bcsr_matrix_type
    integer(kind=ENTIER) :: m, n, nnz !m rows, n columns, nnz non zero elems
    integer(kind=ENTIER) :: blc_size
    integer(kind=ENTIER), dimension(:), allocatable :: col, row
    real(kind=DOUBLE), dimension(:, :, :), allocatable :: v
  end type bcsr_matrix_type
contains
  function bcsr_matmul(A, x) result(y)
    implicit none

    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(A%blc_size*A%n), intent(in) :: x
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: y

    integer(kind=ENTIER) :: i, j

    y = 0.0_DOUBLE
    do i = 1, A%m
      do j = A%row(i), A%row(i + 1) - 1
        y(A%blc_size*(i - 1) + 1:A%blc_size*i) = y(A%blc_size*(i - 1) + 1:A%blc_size*i) &
          + matmul(A%v(:, :, j), x(A%blc_size*(A%col(j) - 1) + 1:A%blc_size*A%col(j)))
      end do
    end do
  end function bcsr_matmul

  function mpi_bcsr_matmul(A, x, mpi_send_recv) result(y)
    use mpi
    use mpi_module
    implicit none

    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(A%blc_size*A%n), intent(inout) :: x
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: y

    integer(kind=ENTIER) :: i, j

    call mpi_memory_exchange(mpi_send_recv, A%m, A%blc_size, x)

    y = 0.0_DOUBLE
    do i = 1, A%m
      if( .not. mpi_send_recv%is_ghost(i) ) then
        do j = A%row(i), A%row(i + 1) - 1
          y(A%blc_size*(i - 1) + 1:A%blc_size*i) = y(A%blc_size*(i - 1) + 1:A%blc_size*i) &
            + matmul(A%v(:, :, j), x(A%blc_size*(A%col(j) - 1) + 1:A%blc_size*A%col(j)))
        end do
      end if
    end do
  end function mpi_bcsr_matmul

  subroutine add_block_to_bcsr(mat, i, j, M)
    implicit none

    type(bcsr_matrix_type), intent(inout) :: mat
    integer(kind=ENTIER), intent(in) :: i, j
    real(kind=DOUBLE), dimension(mat%blc_size, mat%blc_size) :: M

    integer(kind=ENTIER) :: k

    do k = mat%row(i), mat%row(i + 1) - 1
      if (mat%col(k) == j) then
        mat%v(:, :, k) = mat%v(:, :, k) + M(:, :)
        exit
      end if
    end do
  end subroutine add_block_to_bcsr

  subroutine add_block_to_bcsr_thread_safe(mat, i, j, M)
    implicit none

    type(bcsr_matrix_type), intent(inout) :: mat
    integer(kind=ENTIER), intent(in) :: i, j
    real(kind=DOUBLE), dimension(mat%blc_size, mat%blc_size) :: M

    integer(kind=ENTIER) :: k, ii, jj

    do k = mat%row(i), mat%row(i + 1) - 1
      if (mat%col(k) == j) then
        do ii = 1, mat%blc_size
          do jj = 1, mat%blc_size
            mat%v(ii, jj, k) = mat%v(ii, jj, k) + M(ii, jj)
          end do
        end do
        exit
      end if
    end do
  end subroutine add_block_to_bcsr_thread_safe

  subroutine mpi_bcsr_gmres(A, x, b, tol, m, conv, H, V, r, y, c, s, z, mpi_send_recv)
    use mpi
    use mpi_module
    use sparse_csr_linear_module
    implicit none

    integer(kind=ENTIER), intent(in) :: m
    real(kind=DOUBLE), intent(in) :: tol
    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(in) :: b
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(inout) :: x
    logical, intent(inout) :: conv
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    integer(kind=ENTIER), parameter :: maxit = 20
    integer(kind=ENTIER) :: i, j, k, l, nr, n

    real(kind=DOUBLE), dimension(m + 1, m) :: H
    real(kind=DOUBLE), dimension(A%blc_size*A%m, m + 1) :: V
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: r
    real(kind=DOUBLE), dimension(m + 1) :: y, c, s
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: z

    real(kind=DOUBLE) :: delta, rho, tmp

    conv = .false.

    n = A%blc_size*A%m
    x = 0.0_DOUBLE
    r = b - mpi_bcsr_matmul(A, x, mpi_send_recv)
    rho = mpi_norm2(A%m, A%blc_size, r, mpi_send_recv)
    if (rho < tol) then
      conv = .true.
      return
    end if

    do j = 1, maxit
      y = 0.0_DOUBLE
      y(1) = rho
      r = r/y(1)
      nr = m

      call mpi_modifiedArnoldi(n, m, r, A, H, V, mpi_send_recv)

      !Givens rotations
      do i = 1, m
        do k = 2, i
          tmp = H(k - 1, i)
          H(k - 1, i) = c(k - 1)*H(k - 1, i) + s(k - 1)*H(k, i)
          H(k, i) = -s(k - 1)*tmp + c(k - 1)*H(k, i)
        end do

        delta = sqrt(H(i, i)**2 + H(i + 1, i)*H(i + 1, i))
        c(i) = H(i, i)/delta
        s(i) = H(i + 1, i)/delta

        H(i, i) = c(i)*H(i, i) + s(i)*H(i + 1, i)

        do k = i + 1, m + 1
          H(k, i) = 0.0_DOUBLE
        end do

        y(i + 1) = -s(i)*y(i)
        y(i) = c(i)*y(i)
        rho = abs(y(i + 1))
        if (rho < tol) then
          nr = i + 1
          exit
        end if
      end do

      !Backsolve
      z = 0.0_DOUBLE
      if (i >= m) then
        nr = m
        z(nr) = y(nr)/H(nr, nr)
      end if

      do k = nr - 1, 1, -1
        z(k) = y(k)
        do l = k + 1, nr
          z(k) = z(k) - H(k, l)*z(l)
        end do
        z(k) = z(k)/H(k, k)
      end do

      !Linear combination of basis vector
      do i = 1, nr
        x = x + z(i)*V(:, i)
      end do

      r = b - mpi_bcsr_matmul(A, x, mpi_send_recv)
      rho = mpi_norm2(A%m, A%blc_size, r, mpi_send_recv)
      if (rho < tol) then
        conv = .true.
        exit
      end if
    end do
  end subroutine mpi_bcsr_gmres

  subroutine bcsr_gmres(A, x, b, tol, m, conv, H, V, r, y, c, s, z)
    implicit none

    integer(kind=ENTIER), intent(in) :: m
    real(kind=DOUBLE), intent(in) :: tol
    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(in) :: b
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(inout) :: x
    logical, intent(inout) :: conv

    integer(kind=ENTIER), parameter :: maxit = 20
    integer(kind=ENTIER) :: i, j, k, l, nr, n

    real(kind=DOUBLE), dimension(m + 1, m) :: H
    real(kind=DOUBLE), dimension(A%blc_size*A%m, m + 1) :: V
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: r
    real(kind=DOUBLE), dimension(m + 1) :: y, c, s
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: z

    real(kind=DOUBLE) :: delta, rho, tmp

    n = A%blc_size*A%m
    x = 0.0_DOUBLE
    r = b - bcsr_matmul(A, x)
    rho = norm2(r)
    if (rho < tol) then
      conv = .true.
      return
    end if

    do j = 1, maxit
      y = 0.0_DOUBLE
      y(1) = rho
      r = r/y(1)
      nr = m

      call modifiedArnoldi(n, m, r, A, H, V)

      !Givens rotations
      do i = 1, m
        do k = 2, i
          tmp = H(k - 1, i)
          H(k - 1, i) = c(k - 1)*H(k - 1, i) + s(k - 1)*H(k, i)
          H(k, i) = -s(k - 1)*tmp + c(k - 1)*H(k, i)
        end do

        delta = sqrt(H(i, i)**2 + H(i + 1, i)*H(i + 1, i))
        c(i) = H(i, i)/delta
        s(i) = H(i + 1, i)/delta

        H(i, i) = c(i)*H(i, i) + s(i)*H(i + 1, i)

        do k = i + 1, m + 1
          H(k, i) = 0.0_DOUBLE
        end do

        y(i + 1) = -s(i)*y(i)
        y(i) = c(i)*y(i)
        rho = abs(y(i + 1))
        if (rho < tol) then
          nr = i + 1
          exit
        end if
      end do

      !Backsolve
      z = 0.0_DOUBLE
      if (i >= m) then
        nr = m
        z(nr) = y(nr)/H(nr, nr)
      end if

      do k = nr - 1, 1, -1
        z(k) = y(k)
        do l = k + 1, nr
          z(k) = z(k) - H(k, l)*z(l)
        end do
        z(k) = z(k)/H(k, k)
      end do

      !Linear combination of basis vector
      do i = 1, nr
        x = x + z(i)*V(:, i)
      end do

      r = b - bcsr_matmul(A, x)
      rho = norm2(r)
      if (rho < tol) exit
    end do

    if (norm2(b - bcsr_matmul(A, x)) > tol) then
      ! print*, "GMRES DID NOT FULLY CONVERGE", norm2(b - bcsr_matmul(A, x))
      conv = .false.
    else
      conv = .true.
    end if
  end subroutine bcsr_gmres

  subroutine mpi_modifiedArnoldi(n, m, r, A, H, V, mpi_send_recv)
    use mpi
    use mpi_module
    use sparse_csr_linear_module, only: mpi_dot_product, mpi_norm2
    implicit none

    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), intent(in), dimension(n) :: r
    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(m + 1, m), intent(inout) :: H
    real(kind=DOUBLE), dimension(n, m + 1), intent(inout) :: V
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    integer(kind=ENTIER) :: i, j

    V(:, 1) = r
    do j = 2, m + 1
      V(:, j) = mpi_bcsr_matmul(A, V(:, j - 1), mpi_send_recv)
      do i = 1, j - 1
        H(i, j - 1) = mpi_dot_product(A%m, A%blc_size, V(:, i), V(:, j), mpi_send_recv)
        V(:, j) = V(:, j) - H(i, j - 1)*V(:, i)
      end do
      H(j, j - 1) = mpi_norm2(A%m, A%blc_size, V(:, j), mpi_send_recv)
      V(:, j) = V(:, j)/H(j, j - 1)
    end do
  end subroutine mpi_modifiedArnoldi

  subroutine modifiedArnoldi(n, m, r, A, H, V)
    implicit none

    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), intent(in), dimension(n) :: r
    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(m + 1, m), intent(inout) :: H
    real(kind=DOUBLE), dimension(n, m + 1), intent(inout) :: V

    integer(kind=ENTIER) :: i, j

    V(:, 1) = r
    do j = 2, m + 1
      V(:, j) = bcsr_matmul(A, V(:, j - 1))
      do i = 1, j - 1
        H(i, j - 1) = dot_product(V(:, i), V(:, j))
        V(:, j) = V(:, j) - H(i, j - 1)*V(:, i)
      end do
      H(j, j - 1) = norm2(V(:, j))
      V(:, j) = V(:, j)/H(j, j - 1)
    end do
  end subroutine modifiedArnoldi

  subroutine mpi_bcsr_gmres_lusgs(A, x, b, tol, m, conv, H, V, r, y, c, s, z, zp, mpi_send_recv)
    use mpi
    use mpi_module
    use sparse_csr_linear_module, only: mpi_dot_product, mpi_norm2
    implicit none

    integer(kind=ENTIER), intent(in) :: m
    real(kind=DOUBLE), intent(in) :: tol
    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(in) :: b
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(inout) :: x
    logical, intent(inout) :: conv
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    integer(kind=ENTIER), parameter :: maxit = 20
    integer(kind=ENTIER) :: i, j, k, l, nr, n

    real(kind=DOUBLE), dimension(m + 1, m) :: H
    real(kind=DOUBLE), dimension(A%blc_size*A%m, m + 1) :: V
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: r
    real(kind=DOUBLE), dimension(m + 1) :: y, c, s
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: z, zp

    real(kind=DOUBLE) :: delta, rho, tmp

    conv = .false.

    n = A%blc_size*A%m
    x = 0.0_DOUBLE
    r = b - mpi_bcsr_matmul(A, x, mpi_send_recv)
    rho = mpi_norm2(A%m, A%blc_size, r, mpi_send_recv)
    if (rho < tol) then
      conv = .true.
      return
    end if

    do j = 1, maxit
      y = 0.0_DOUBLE
      y(1) = rho
      r = r/y(1)
      nr = m

      call mpi_modifiedArnoldi_lusgs(n, m, r, A, H, V, mpi_send_recv)

      !Givens rotations
      do i = 1, m
        do k = 2, i
          tmp = H(k - 1, i)
          H(k - 1, i) = c(k - 1)*H(k - 1, i) + s(k - 1)*H(k, i)
          H(k, i) = -s(k - 1)*tmp + c(k - 1)*H(k, i)
        end do

        delta = sqrt(H(i, i)**2 + H(i + 1, i)*H(i + 1, i))
        c(i) = H(i, i)/delta
        s(i) = H(i + 1, i)/delta

        H(i, i) = c(i)*H(i, i) + s(i)*H(i + 1, i)

        do k = i + 1, m + 1
          H(k, i) = 0.0_DOUBLE
        end do

        y(i + 1) = -s(i)*y(i)
        y(i) = c(i)*y(i)
        rho = abs(y(i + 1))
        if (rho < tol) then
          nr = i + 1
          exit
        end if
      end do

      !Backsolve
      z = 0.0_DOUBLE
      if (i >= m) then
        nr = m
        z(nr) = y(nr)/H(nr, nr)
      end if

      do k = nr - 1, 1, -1
        z(k) = y(k)
        do l = k + 1, nr
          z(k) = z(k) - H(k, l)*z(l)
        end do
        z(k) = z(k)/H(k, k)
      end do

      !Linear combination of basis vector
      r = 0.0_DOUBLE !reuse r to save memory
      do i = 1, nr
        !x = x + z(i)*V(:, i)
        r = r + z(i)*V(:, i)
      end do

      call bcsr_lusgs(A, zp, r)
      x = x + zp

      r = b - mpi_bcsr_matmul(A, x, mpi_send_recv)
      rho = mpi_norm2(A%m, A%blc_size, r, mpi_send_recv)
      if (rho < tol) then
        conv = .true.
        exit
      end if
    end do
  end subroutine mpi_bcsr_gmres_lusgs

  subroutine bcsr_gmres_lusgs(A, x, b, tol, m, conv, H, V, r, y, c, s, z, zp)
    implicit none

    integer(kind=ENTIER), intent(in) :: m
    real(kind=DOUBLE), intent(in) :: tol
    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(in) :: b
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(inout) :: x
    logical, intent(inout) :: conv

    integer(kind=ENTIER), parameter :: maxit = 20
    integer(kind=ENTIER) :: i, j, k, l, nr, n

    real(kind=DOUBLE), dimension(m + 1, m) :: H
    real(kind=DOUBLE), dimension(A%blc_size*A%m, m + 1) :: V
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: r
    real(kind=DOUBLE), dimension(m + 1) :: y, c, s
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: z, zp

    real(kind=DOUBLE) :: delta, rho, tmp

    n = A%blc_size*A%m
    x = 0.0_DOUBLE
    r = b - bcsr_matmul(A, x)
    rho = norm2(r)
    if (rho < tol) then
      conv = .true.
      return
    end if

    do j = 1, maxit
      y = 0.0_DOUBLE
      y(1) = rho
      r = r/y(1)
      nr = m

      call modifiedArnoldi_lusgs(n, m, r, A, H, V)

      !Givens rotations
      do i = 1, m
        do k = 2, i
          tmp = H(k - 1, i)
          H(k - 1, i) = c(k - 1)*H(k - 1, i) + s(k - 1)*H(k, i)
          H(k, i) = -s(k - 1)*tmp + c(k - 1)*H(k, i)
        end do

        delta = sqrt(H(i, i)**2 + H(i + 1, i)*H(i + 1, i))
        c(i) = H(i, i)/delta
        s(i) = H(i + 1, i)/delta

        H(i, i) = c(i)*H(i, i) + s(i)*H(i + 1, i)

        do k = i + 1, m + 1
          H(k, i) = 0.0_DOUBLE
        end do

        y(i + 1) = -s(i)*y(i)
        y(i) = c(i)*y(i)
        rho = abs(y(i + 1))
        if (rho < tol) then
          nr = i + 1
          exit
        end if
      end do

      !Backsolve
      z = 0.0_DOUBLE
      if (i >= m) then
        nr = m
        z(nr) = y(nr)/H(nr, nr)
      end if

      do k = nr - 1, 1, -1
        z(k) = y(k)
        do l = k + 1, nr
          z(k) = z(k) - H(k, l)*z(l)
        end do
        z(k) = z(k)/H(k, k)
      end do

      !Linear combination of basis vector
      r = 0.0_DOUBLE !reuse r to save memory
      do i = 1, nr
        !x = x + z(i)*V(:, i)
        r = r + z(i)*V(:, i)
      end do

      call bcsr_lusgs(A, zp, r)
      x = x + zp

      r = b - bcsr_matmul(A, x)
      rho = norm2(r)
      if (rho < tol) exit
    end do

    if (norm2(b - bcsr_matmul(A, x)) > tol) then
      ! print*, "GMRES DID NOT FULLY CONVERGE", norm2(b - bcsr_matmul(A, x))
      conv = .false.
    else
      conv = .true.
    end if
  end subroutine bcsr_gmres_lusgs

  subroutine mpi_modifiedArnoldi_parallel_lusgs(n, m, r, A, H, V, &
      n_color_sets, color_set, color, mpi_send_recv)
    use mpi
    use mpi_module
    use sparse_csr_linear_module, only: color_set_type, mpi_dot_product, mpi_norm2
    implicit none

    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), intent(in), dimension(n) :: r
    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(m + 1, m), intent(inout) :: H
    real(kind=DOUBLE), dimension(n, m + 1), intent(inout) :: V
    integer(kind=ENTIER), intent(in) :: n_color_sets
    type(color_set_type), dimension(n_color_sets), intent(in) :: color_set
    integer(kind=ENTIER), dimension(n), intent(in) :: color
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    integer(kind=ENTIER) :: i, j

    V(:, 1) = r
    do j = 2, m + 1
      call mpi_bcsr_parallel_lusgs(A, V(:, j), V(:, j - 1), &
        n_color_sets, color_set, color, mpi_send_recv)
      V(:, j) = mpi_bcsr_matmul(A, V(:, j), mpi_send_recv)
      do i = 1, j - 1
        H(i, j - 1) = mpi_dot_product(A%m, A%blc_size, V(:, i), V(:, j), mpi_send_recv)
        V(:, j) = V(:, j) - H(i, j - 1)*V(:, i)
      end do
      H(j, j - 1) = mpi_norm2(A%m, A%blc_size, V(:, j), mpi_send_recv)
      V(:, j) = V(:, j)/H(j, j - 1)
    end do
  end subroutine mpi_modifiedArnoldi_parallel_lusgs

  subroutine modifiedArnoldi_parallel_lusgs(n, m, r, A, H, V, &
      n_color_sets, color_set, color)
    use sparse_csr_linear_module, only: color_set_type
    implicit none

    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), intent(in), dimension(n) :: r
    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(m + 1, m), intent(inout) :: H
    real(kind=DOUBLE), dimension(n, m + 1), intent(inout) :: V
    integer(kind=ENTIER), intent(in) :: n_color_sets
    type(color_set_type), dimension(n_color_sets), intent(in) :: color_set
    integer(kind=ENTIER), dimension(n), intent(in) :: color

    integer(kind=ENTIER) :: i, j

    V(:, 1) = r
    do j = 2, m + 1
      call bcsr_parallel_lusgs(A, V(:, j), V(:, j - 1), &
        n_color_sets, color_set, color)
      V(:, j) = bcsr_matmul(A, V(:, j))
      do i = 1, j - 1
        H(i, j - 1) = dot_product(V(:, i), V(:, j))
        V(:, j) = V(:, j) - H(i, j - 1)*V(:, i)
      end do
      H(j, j - 1) = norm2(V(:, j))
      V(:, j) = V(:, j)/H(j, j - 1)
    end do
  end subroutine modifiedArnoldi_parallel_lusgs

  subroutine mpi_modifiedArnoldi_lusgs(n, m, r, A, H, V, mpi_send_recv)
    use mpi
    use mpi_module
    use sparse_csr_linear_module, only: mpi_dot_product, mpi_norm2
    implicit none

    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), intent(in), dimension(n) :: r
    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(m + 1, m), intent(inout) :: H
    real(kind=DOUBLE), dimension(n, m + 1), intent(inout) :: V
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    integer(kind=ENTIER) :: i, j

    V(:, 1) = r
    do j = 2, m + 1
      !V(:, j) = bcsr_matmul(A, V(:, j - 1))
      call bcsr_lusgs(A, V(:, j), V(:, j - 1))
      V(:, j) = mpi_bcsr_matmul(A, V(:, j), mpi_send_recv)
      do i = 1, j - 1
        H(i, j - 1) = mpi_dot_product(A%m, A%blc_size, V(:, i), V(:, j), mpi_send_recv)
        V(:, j) = V(:, j) - H(i, j - 1)*V(:, i)
      end do
      H(j, j - 1) = mpi_norm2(A%m, A%blc_size, V(:, j), mpi_send_recv)
      V(:, j) = V(:, j)/H(j, j - 1)
    end do
  end subroutine mpi_modifiedArnoldi_lusgs

  subroutine modifiedArnoldi_lusgs(n, m, r, A, H, V)
    implicit none

    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), intent(in), dimension(n) :: r
    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(m + 1, m), intent(inout) :: H
    real(kind=DOUBLE), dimension(n, m + 1), intent(inout) :: V

    integer(kind=ENTIER) :: i, j

    V(:, 1) = r
    do j = 2, m + 1
      !V(:, j) = bcsr_matmul(A, V(:, j - 1))
      call bcsr_lusgs(A, V(:, j), V(:, j - 1))
      V(:, j) = bcsr_matmul(A, V(:, j))
      do i = 1, j - 1
        H(i, j - 1) = dot_product(V(:, i), V(:, j))
        V(:, j) = V(:, j) - H(i, j - 1)*V(:, i)
      end do
      H(j, j - 1) = norm2(V(:, j))
      V(:, j) = V(:, j)/H(j, j - 1)
    end do
  end subroutine modifiedArnoldi_lusgs

  subroutine mpi_bcsr_gmres_parallel_lusgs(A, x, b, tol, m, &
      n_color_sets, color_set, color, conv, &
      H, V, r, y, c, s, z, zp, mpi_send_recv)
    use mpi
    use mpi_module
    use sparse_csr_linear_module, only: mpi_dot_product, color_set_type, mpi_norm2
    implicit none

    integer(kind=ENTIER), intent(in) :: m
    real(kind=DOUBLE), intent(in) :: tol
    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(in) :: b
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(inout) :: x
    logical, intent(inout) :: conv
    integer(kind=ENTIER), intent(in) :: n_color_sets
    type(color_set_type), dimension(n_color_sets), intent(in) :: color_set
    integer(kind=ENTIER), dimension(A%blc_size*A%m), intent(in) :: color
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    integer(kind=ENTIER), parameter :: maxit = 40
    integer(kind=ENTIER) :: i, j, k, l, nr, n

    real(kind=DOUBLE), dimension(m + 1, m) :: H
    real(kind=DOUBLE), dimension(A%blc_size*A%m, m + 1) :: V
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: r
    real(kind=DOUBLE), dimension(m + 1) :: y, c, s
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: z, zp

    real(kind=DOUBLE) :: delta, rho, tmp

    conv = .false.

    n = A%blc_size*A%m
    x = 0.0_DOUBLE
    r = b - mpi_bcsr_matmul(A, x, mpi_send_recv)
    rho = mpi_norm2(A%m, A%blc_size, r, mpi_send_recv)
    if (rho < tol) then
      conv = .true.
      return
    end if

    do j = 1, maxit
      y = 0.0_DOUBLE
      y(1) = rho
      r = r/y(1)
      nr = m

      call mpi_modifiedArnoldi_parallel_lusgs(n, m, r, A, H, V, &
        n_color_sets, color_set, color, mpi_send_recv)

      !Givens rotations
      do i = 1, m
        do k = 2, i
          tmp = H(k - 1, i)
          H(k - 1, i) = c(k - 1)*H(k - 1, i) + s(k - 1)*H(k, i)
          H(k, i) = -s(k - 1)*tmp + c(k - 1)*H(k, i)
        end do

        delta = sqrt(H(i, i)**2 + H(i + 1, i)*H(i + 1, i))
        c(i) = H(i, i)/delta
        s(i) = H(i + 1, i)/delta

        H(i, i) = c(i)*H(i, i) + s(i)*H(i + 1, i)

        do k = i + 1, m + 1
          H(k, i) = 0.0_DOUBLE
        end do

        y(i + 1) = -s(i)*y(i)
        y(i) = c(i)*y(i)
        rho = abs(y(i + 1))
        if (rho < tol) then
          nr = i + 1
          exit
        end if
      end do

      !Backsolve
      z = 0.0_DOUBLE
      if (i >= m) then
        nr = m
        z(nr) = y(nr)/H(nr, nr)
      end if

      do k = nr - 1, 1, -1
        z(k) = y(k)
        do l = k + 1, nr
          z(k) = z(k) - H(k, l)*z(l)
        end do
        z(k) = z(k)/H(k, k)
      end do

      !Linear combination of basis vector
      r = 0.0_DOUBLE !reuse r to save memory
      do i = 1, nr
        r = r + z(i)*V(:, i)
      end do

      call mpi_bcsr_parallel_lusgs(A, zp, r, &
        n_color_sets, color_set, color, mpi_send_recv)
      x = x + zp

      r = b - mpi_bcsr_matmul(A, x, mpi_send_recv)
      rho = mpi_norm2(A%m, A%blc_size, r, mpi_send_recv)
      if (rho < tol) then
        conv = .true.
        exit
      end if
    end do
  end subroutine mpi_bcsr_gmres_parallel_lusgs

  subroutine mpi_bcsr_gmres_diag(A, x, b, tol, m, &
      conv, H, V, r, y, c, s, z, zp, mpi_send_recv)
    use mpi
    use mpi_module
    use sparse_csr_linear_module, only: mpi_dot_product, color_set_type, mpi_norm2
    implicit none

    integer(kind=ENTIER), intent(in) :: m
    real(kind=DOUBLE), intent(in) :: tol
    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(in) :: b
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(inout) :: x
    logical, intent(inout) :: conv
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    integer(kind=ENTIER), parameter :: maxit = 20
    integer(kind=ENTIER) :: i, j, k, l, nr, n, me, mpi_ierr

    real(kind=DOUBLE), dimension(m + 1, m) :: H
    real(kind=DOUBLE), dimension(A%blc_size*A%m, m + 1) :: V
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: r
    real(kind=DOUBLE), dimension(m + 1) :: y, c, s
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: z, zp

    real(kind=DOUBLE) :: delta, rho, tmp

    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)
    conv = .false.

    n = A%blc_size*A%m
    x = 0.0_DOUBLE
    r = b - mpi_bcsr_matmul(A, x, mpi_send_recv)
    rho = mpi_norm2(A%m, A%blc_size, r, mpi_send_recv)
    if (rho < tol) then
      conv = .true.
      return
    end if

    do j = 1, maxit
      y = 0.0_DOUBLE
      y(1) = rho
      r = r/y(1)
      nr = m

      call mpi_modifiedArnoldi_diag(n, m, r, A, H, V, mpi_send_recv)

      !Givens rotations
      do i = 1, m
        do k = 2, i
          tmp = H(k - 1, i)
          H(k - 1, i) = c(k - 1)*H(k - 1, i) + s(k - 1)*H(k, i)
          H(k, i) = -s(k - 1)*tmp + c(k - 1)*H(k, i)
        end do

        delta = sqrt(H(i, i)**2 + H(i + 1, i)*H(i + 1, i))
        c(i) = H(i, i)/delta
        s(i) = H(i + 1, i)/delta

        H(i, i) = c(i)*H(i, i) + s(i)*H(i + 1, i)

        do k = i + 1, m + 1
          H(k, i) = 0.0_DOUBLE
        end do

        y(i + 1) = -s(i)*y(i)
        y(i) = c(i)*y(i)
        rho = abs(y(i + 1))
        if (rho < tol) then
          nr = i + 1
          exit
        end if
      end do

      !Backsolve
      z = 0.0_DOUBLE
      if (i >= m) then
        nr = m
        z(nr) = y(nr)/H(nr, nr)
      end if

      do k = nr - 1, 1, -1
        z(k) = y(k)
        do l = k + 1, nr
          z(k) = z(k) - H(k, l)*z(l)
        end do
        z(k) = z(k)/H(k, k)
      end do

      !Linear combination of basis vector
      r = 0.0_DOUBLE !reuse r to save memory
      do i = 1, nr
        r = r + z(i)*V(:, i)
      end do

      call bcsr_diag(A, zp, r)
      x = x + zp

      r = b - mpi_bcsr_matmul(A, x, mpi_send_recv)
      rho = mpi_norm2(A%m, A%blc_size, r, mpi_send_recv)
      !if ( me == 0 ) print*, j, rho
      if (rho < tol) then
        conv = .true.
        exit
      end if
    end do
  end subroutine mpi_bcsr_gmres_diag

  subroutine bcsr_gmres_parallel_lusgs(A, x, b, tol, m, &
      n_color_sets, color_set, color, conv, &
      H, V, r, y, c, s, z, zp)
    use sparse_csr_linear_module, only: color_set_type
    implicit none

    integer(kind=ENTIER), intent(in) :: m
    real(kind=DOUBLE), intent(in) :: tol
    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(in) :: b
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(inout) :: x
    logical, intent(inout) :: conv
    integer(kind=ENTIER), intent(in) :: n_color_sets
    type(color_set_type), dimension(n_color_sets), intent(in) :: color_set
    integer(kind=ENTIER), dimension(A%blc_size*A%m), intent(in) :: color

    integer(kind=ENTIER), parameter :: maxit = 20
    integer(kind=ENTIER) :: i, j, k, l, nr, n

    real(kind=DOUBLE), dimension(m + 1, m) :: H
    real(kind=DOUBLE), dimension(A%blc_size*A%m, m + 1) :: V
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: r
    real(kind=DOUBLE), dimension(m + 1) :: y, c, s
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: z, zp

    real(kind=DOUBLE) :: delta, rho, tmp

    n = A%blc_size*A%m
    x = 0.0_DOUBLE
    r = b - bcsr_matmul(A, x)
    rho = norm2(r)
    if (rho < tol) then
      conv = .true.
      return
    end if

    do j = 1, maxit
      y = 0.0_DOUBLE
      y(1) = rho
      r = r/y(1)
      nr = m

      call modifiedArnoldi_parallel_lusgs(n, m, r, A, H, V, &
        n_color_sets, color_set, color)

      !Givens rotations
      do i = 1, m
        do k = 2, i
          tmp = H(k - 1, i)
          H(k - 1, i) = c(k - 1)*H(k - 1, i) + s(k - 1)*H(k, i)
          H(k, i) = -s(k - 1)*tmp + c(k - 1)*H(k, i)
        end do

        delta = sqrt(H(i, i)**2 + H(i + 1, i)*H(i + 1, i))
        c(i) = H(i, i)/delta
        s(i) = H(i + 1, i)/delta

        H(i, i) = c(i)*H(i, i) + s(i)*H(i + 1, i)

        do k = i + 1, m + 1
          H(k, i) = 0.0_DOUBLE
        end do

        y(i + 1) = -s(i)*y(i)
        y(i) = c(i)*y(i)
        rho = abs(y(i + 1))
        if (rho < tol) then
          nr = i + 1
          exit
        end if
      end do

      !Backsolve
      z = 0.0_DOUBLE
      if (i >= m) then
        nr = m
        z(nr) = y(nr)/H(nr, nr)
      end if

      do k = nr - 1, 1, -1
        z(k) = y(k)
        do l = k + 1, nr
          z(k) = z(k) - H(k, l)*z(l)
        end do
        z(k) = z(k)/H(k, k)
      end do

      !Linear combination of basis vector
      r = 0.0_DOUBLE !reuse r to save memory
      do i = 1, nr
        r = r + z(i)*V(:, i)
      end do

      call bcsr_parallel_lusgs(A, zp, r, &
        n_color_sets, color_set, color)

      x = x + zp

      r = b - bcsr_matmul(A, x)
      rho = norm2(r)
      if (rho < tol) exit
    end do

    if (norm2(b - bcsr_matmul(A, x)) > tol) then
      conv = .false.
    else
      conv = .true.
    end if
  end subroutine bcsr_gmres_parallel_lusgs

  !subroutine bcsr_pcgm_cholesky(n, A, L, x, b, tol)
  !  implicit none
  !
  !  integer(kind=ENTIER), intent(in) :: n
  !  real(kind=DOUBLE), intent(in) :: tol
  !  type(bcsr_matrix_type), intent(in) :: A, L
  !  real(kind=DOUBLE), dimension(n), intent(in) :: b
  !  real(kind=DOUBLE), dimension(n), intent(inout) :: x
  !
  !  integer(kind=ENTIER) :: k
  !  real(kind=DOUBLE) :: alpha, rdr, rpdrp
  !  real(kind=DOUBLE), dimension(n) :: r, p, Ap, mr
  !
  !  r = b - bcsr_matmul(A, x)
  !
  !  call bcsr_solve_cholesky(L, mr, r)
  !
  !  rdr = dot_product(mr, r)
  !  if (sqrt(rdr) < tol) then
  !    return
  !  end if
  !  p = mr
  !
  !  do k = 1, n
  !    Ap = bcsr_matmul(A, p)
  !    alpha = rdr/dot_product(p, Ap)
  !    x = x + alpha*p
  !    r = r - alpha*Ap
  !    call bcsr_solve_cholesky(L, mr, r)
  !    rpdrp = dot_product(mr, r)
  !    if (sqrt(rpdrp) <= tol) then
  !      if (norm2(r) <= tol) exit
  !    end if
  !    p = mr + (rpdrp/rdr)*p
  !    rdr = rpdrp
  !  end do
  !
  !  if (sqrt(rpdrp) > tol) then
  !    print *, "PCGM did not converge fully ! Redisu = ", sqrt(rpdrp)
  !    error stop
  !  end if
  !end subroutine bcsr_pcgm_cholesky
  !

  subroutine mpi_bcsr_pcgm(A, x, b, tol, mpi_send_recv)
    use mpi
    use mpi_module
    use sparse_csr_linear_module, only: mpi_dot_product, mpi_norm2
    implicit none

    real(kind=DOUBLE), intent(in) :: tol
    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(in) :: b
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(inout) :: x
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    integer(kind=ENTIER), parameter :: maxit = 1000
    integer(kind=ENTIER) :: k, i, n
    real(kind=DOUBLE) :: alpha, rdr, rpdrp
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: r, p, Ap, mr
    integer(kind=ENTIER), dimension(A%m) :: d

    n = A%blc_size*A%m
    mr = 0.0_DOUBLE

    !Get inverse diagonal of A
    do i = 1, A%m
      if( .not. mpi_send_recv%is_ghost(i) ) then
        do k = A%row(i), A%row(i + 1) - 1
          if (A%col(k) == i) then
            d(i) = k
            exit
          end if
        end do
      end if
    end do

    r = b - mpi_bcsr_matmul(A, x, mpi_send_recv)
    do i = 1, A%m
      if( .not. mpi_send_recv%is_ghost(i) ) then
        call lu_solve(3, A%v(:, :, d(i)), mr(A%blc_size*(i - 1) + 1:A%blc_size*i), &
          r(A%blc_size*(i - 1) + 1:A%blc_size*i))
      end if
    end do
    rdr = mpi_dot_product(A%m, A%blc_size, mr, r, mpi_send_recv)
    if (sqrt(rdr) < tol) return

    p = mr
    do k = 1, maxit
      Ap = mpi_bcsr_matmul(A, p, mpi_send_recv)
      alpha = rdr/mpi_dot_product(A%m, A%blc_size, p, Ap, mpi_send_recv)
      x = x + alpha*p
      r = r - alpha*Ap
      do i = 1, A%m
        if( .not. mpi_send_recv%is_ghost(i) ) then
          call lu_solve(3, A%v(:, :, d(i)), mr(A%blc_size*(i - 1) + 1:A%blc_size*i), &
            r(A%blc_size*(i - 1) + 1:A%blc_size*i))
        end if
      end do
      rpdrp = mpi_dot_product(A%m, A%blc_size, mr, r, mpi_send_recv)
      ! print*, "Residu ", n, k, sqrt(rpdrp)
      if (sqrt(rpdrp) <= tol) then
        if (mpi_dot_product(A%m, A%blc_size, mr, r, mpi_send_recv) <= tol) exit
      end if
      p = mr + (rpdrp/rdr)*p
      rdr = rpdrp
    end do

    !print*, norm2(b - bcsr_matmul(A, x))
    if (sqrt(rpdrp) > tol) then
      print *, "PCGM did not converge fully ! Redisu = ", sqrt(rpdrp)
      error stop
    end if
  end subroutine mpi_bcsr_pcgm

  subroutine bcsr_pcgm(A, x, b, tol)
    implicit none

    real(kind=DOUBLE), intent(in) :: tol
    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(in) :: b
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(inout) :: x

    integer(kind=ENTIER) :: k, i, n
    real(kind=DOUBLE) :: alpha, rdr, rpdrp
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: r, p, Ap, mr
    integer(kind=ENTIER), dimension(A%m) :: d

    n = A%blc_size*A%m
    rpdrp = 0.0_DOUBLE

    !Get inverse diagonal of A
    do i = 1, A%m
      do k = A%row(i), A%row(i + 1) - 1
        if (A%col(k) == i) then
          d(i) = k
          exit
        end if
      end do
    end do

    r = b - bcsr_matmul(A, x)
    mr = 0.0_DOUBLE
    do i = 1, A%m
      call lu_solve(3, A%v(:, :, d(i)), mr(A%blc_size*(i - 1) + 1:A%blc_size*i), &
        r(A%blc_size*(i - 1) + 1:A%blc_size*i))
    end do
    rdr = dot_product(mr, r)
    if (sqrt(rdr) < tol) return
    p = mr

    do k = 1, n
      Ap = bcsr_matmul(A, p)
      alpha = rdr/dot_product(p, Ap)
      x = x + alpha*p
      r = r - alpha*Ap
      do i = 1, A%m
        call lu_solve(3, A%v(:, :, d(i)), mr(A%blc_size*(i - 1) + 1:A%blc_size*i), &
          r(A%blc_size*(i - 1) + 1:A%blc_size*i))
      end do
      rpdrp = dot_product(mr, r)
      ! print*, "Residu ", n, k, sqrt(rpdrp)
      if (sqrt(rpdrp) <= tol) then
        if (norm2(r) <= tol) exit
      end if
      p = mr + (rpdrp/rdr)*p
      rdr = rpdrp
    end do

    !print*, norm2(b - bcsr_matmul(A, x))
    if (sqrt(rpdrp) > tol) then
      print *, "PCGM did not converge fully ! Redisu = ", sqrt(rpdrp)
      error stop
    end if
  end subroutine bcsr_pcgm

  !subroutine bcsr_solve_cholesky(L, x, b)
  !  implicit none
  !
  !  type(bcsr_matrix_type), intent(in) :: L
  !  real(kind=DOUBLE), dimension(L%n) :: x
  !  real(kind=DOUBLE), dimension(L%m) :: b
  !
  !  integer(kind=ENTIER) :: i, j, kdiag, k
  !  real(kind=DOUBLE), dimension(L%n) :: y
  !  real(kind=DOUBLE) :: s
  !
  !  if (L%n /= L%m) then
  !    print *, 'Bad sizes bcsr_solve_cholesky'
  !    error stop
  !  end if
  !
  !  y = 0.0_DOUBLE
  !  !L y = b
  !  do i = 1, L%n
  !    s = 0.0_DOUBLE
  !    do k = L%row(i), L%row(i + 1) - 1
  !      j = L%col(k)
  !      if (j < i) then
  !        s = s + L%v(k)*y(j)
  !      else if (j == i) then
  !        kdiag = k
  !        exit
  !      end if
  !    end do
  !    y(i) = (b(i) - s)/L%v(kdiag)
  !  end do
  !
  !  !L^t x = y
  !  x = y
  !  do i = L%n, 1, -1
  !    do k = L%row(i + 1) - 1, L%row(i), -1
  !      j = L%col(k)
  !      if (j < i) then
  !        x(j) = x(j) - L%v(k)*x(i)
  !      else if (j == i) then
  !        x(i) = x(i)/L%v(k)
  !      end if
  !    end do
  !  end do
  !end subroutine bcsr_solve_cholesky

  !subroutine bcsr_incomplete_cholesky(A, L)
  !  implicit none
  !
  !  type(bcsr_matrix_type), intent(in) :: A
  !  type(bcsr_matrix_type), intent(inout) :: L
  !
  !  integer(kind=ENTIER) :: i, j, k, i1, i2
  !  real(kind=DOUBLE) :: s
  !  integer(kind=ENTIER), dimension(A%m) :: diag
  !
  !  if (A%m /= A%n) then
  !    print *, "BAD SIZES CSR INCOMP CHOL"
  !    error stop
  !  end if
  !
  !  !Bsed on ILUPP python library
  !  L%v(:) = 0.0_DOUBLE
  !  do i = 1, A%m
  !    do k = L%row(i), L%row(i + 1) - 1
  !      j = L%col(k)
  !      if (j <= i) then
  !
  !        i1 = L%row(i)
  !        i2 = L%row(j)
  !        s = 0.0_DOUBLE
  !        do while (L%col(i1) < i .and. &
  !            L%col(i2) < i)
  !          if (L%col(i1) == L%col(i2)) then
  !            s = s + L%v(i1)*L%v(i2)
  !            i1 = i1 + 1
  !            i2 = i2 + 1
  !          else if (L%col(i1) < L%col(i2)) then
  !            i1 = i1 + 1
  !          else
  !            i2 = i2 + 1
  !          end if
  !        end do
  !
  !        if (j == i) then
  !          L%v(k) = sqrt(A%v(k) - s)
  !          diag(i) = k
  !        else if (j < i) then
  !          L%v(k) = (A%v(k) - s)/L%v(diag(j))
  !        end if
  !      else
  !        exit
  !      end if
  !    end do
  !  end do
  !end subroutine bcsr_incomplete_cholesky

  subroutine bcsr_lusgs(A, x, b)
    use linear_solver_module, only: lu_solve_noswp
    implicit none

    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(in) :: b
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(inout) :: x

    integer(kind=ENTIER) :: i, j
    integer(kind=ENTIER), dimension(A%blc_size*A%m) :: d
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: y
    real(kind=DOUBLE), dimension(A%blc_size, A%blc_size) :: tmp

    !Get diagonal of A
    do i = 1, A%m
      do j = A%row(i), A%row(i + 1) - 1
        if (A%col(j) == i) then
          d(i) = j
          exit
        end if
      end do
    end do

    !Diag
    !do i = 1, A%m
    !call lu_solve_noswp(A%blc_size, A%v(:, :, d(i)), &
    !x(A%blc_size*(i - 1) + 1:A%blc_size*i), b(A%blc_size*(i - 1) + 1:A%blc_size*i), tmp)
    !end do

    !(L+D)D-1(U+D) x = b
    !(L+D) y = b
    !(U+D) x = D y
    do i = 1, A%m
      y(A%blc_size*(i - 1) + 1:A%blc_size*i) = b(A%blc_size*(i - 1) + 1:A%blc_size*i)
      j = A%row(i)
      do while (A%col(j) < i)
        y(A%blc_size*(i - 1) + 1:A%blc_size*i) = y(A%blc_size*(i - 1) + 1:A%blc_size*i) &
          - matmul(A%v(:, :, j), y(A%blc_size*(A%col(j) - 1) + 1:A%blc_size*A%col(j)))
        j = j + 1
      end do
    end do

    do i = 1, A%m
      x(A%blc_size*(i - 1) + 1:A%blc_size*i) = y(A%blc_size*(i - 1) + 1:A%blc_size*i)
      j = A%row(i)
      do while (A%col(j) < i)
        x(A%blc_size*(i - 1) + 1:A%blc_size*i) = x(A%blc_size*(i - 1) + 1:A%blc_size*i) &
          - matmul(A%v(:, :, j), x(A%blc_size*(A%col(j) - 1) + 1:A%blc_size*A%col(j)))
        j = j + 1
      end do
      y(A%blc_size*(i - 1) + 1:A%blc_size*i) = x(A%blc_size*(i - 1) + 1:A%blc_size*i)
      call lu_solve_noswp(A%blc_size, A%v(:, :, d(i)), &
        x(A%blc_size*(i - 1) + 1:A%blc_size*i), y(A%blc_size*(i - 1) + 1:A%blc_size*i), tmp)
    end do

    !do i = A%m, 1, -1
    !  x(A%blc_size*(i - 1) + 1:A%blc_size*i) = y(A%blc_size*(i - 1) + 1:A%blc_size*i)
    !  j = A%row(i+1)-1
    !  do while (j > d(i))
    !    x(A%blc_size*(i - 1) + 1:A%blc_size*i) = x(A%blc_size*(i - 1) + 1:A%blc_size*i) &
    !      - matmul(A%v(:, :, j), x(A%blc_size*(A%col(j) - 1) + 1:A%blc_size*A%col(j)))
    !    j = j - 1
    !  end do
    !  y(A%blc_size*(i - 1) + 1:A%blc_size*i) = x(A%blc_size*(i - 1) + 1:A%blc_size*i)
    !  call lu_solve_noswp(A%blc_size, A%v(:, :, d(i)), &
    !    x(A%blc_size*(i - 1) + 1:A%blc_size*i), y(A%blc_size*(i - 1) + 1:A%blc_size*i), tmp)
    !end do
  end subroutine bcsr_lusgs

  subroutine mpi_bcsr_parallel_lusgs(A, x, b, &
      n_color_sets, color_set, color, mpi_send_recv)
    use sparse_csr_linear_module, only: color_set_type
    use linear_solver_module, only: lu_solve_noswp
    use mpi
    use mpi_module
    implicit none

    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(in) :: b
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(inout) :: x
    integer(kind=ENTIER), intent(in) :: n_color_sets
    type(color_set_type), dimension(n_color_sets), intent(in) :: color_set
    integer(kind=ENTIER), dimension(A%m), intent(in) :: color
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    integer(kind=ENTIER) :: i, j, cs, ics
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: y
    integer(kind=ENTIER), dimension(A%blc_size*A%m) :: d
    real(kind=DOUBLE), dimension(A%blc_size, A%blc_size) :: tmp

    !Get diagonal of A
    do i = 1, A%m
      if( .not. mpi_send_recv%is_ghost(i) ) then
        do j = A%row(i), A%row(i + 1) - 1
          if (A%col(j) == i) then
            d(i) = j
            exit
          end if
        end do
      end if
    end do

    !(L+D)D-1(U+D) x = b
    !(L+D) y = b
    !(U+D) x = D y

    do cs = 1, n_color_sets
      do ics = 1, color_set(cs)%n_elems
        i = color_set(cs)%color(ics)
        if( .not. mpi_send_recv%is_ghost(i) ) then
          y(A%blc_size*(i - 1) + 1:A%blc_size*i) = b(A%blc_size*(i - 1) + 1:A%blc_size*i)
          do j = A%row(i), A%row(i + 1) - 1
            if (color(A%col(j)) < cs) then
              !y(i) = y(i) - A%v(j)*y(A%col(j))
              y(A%blc_size*(i - 1) + 1:A%blc_size*i) = y(A%blc_size*(i - 1) + 1:A%blc_size*i) &
                - matmul(A%v(:, :, j), y(A%blc_size*(A%col(j) - 1) + 1:A%blc_size*A%col(j)))
            end if
          end do
        end if
      end do
      call mpi_memory_exchange(mpi_send_recv, A%m, A%blc_size, y)
    end do

    do cs = n_color_sets, 1, -1
      do ics = 1, color_set(cs)%n_elems
        i = color_set(cs)%color(ics)
        if( .not. mpi_send_recv%is_ghost(i) ) then
          !x(i) = y(i)
          x(A%blc_size*(i - 1) + 1:A%blc_size*i) = y(A%blc_size*(i - 1) + 1:A%blc_size*i)
          do j = A%row(i), A%row(i + 1) - 1
            if (color(A%col(j)) > cs) then
              ! x(i) = x(i) - A%v(j)*x(A%col(j))
              x(A%blc_size*(i - 1) + 1:A%blc_size*i) = x(A%blc_size*(i - 1) + 1:A%blc_size*i) &
                - matmul(A%v(:, :, j), x(A%blc_size*(A%col(j) - 1) + 1:A%blc_size*A%col(j)))
            end if
          end do
          !x(i) = x(i)/d(i)
          y(A%blc_size*(i - 1) + 1:A%blc_size*i) = x(A%blc_size*(i - 1) + 1:A%blc_size*i)
          call lu_solve_noswp(A%blc_size, A%v(:, :, d(i)), &
            x(A%blc_size*(i - 1) + 1:A%blc_size*i), y(A%blc_size*(i - 1) + 1:A%blc_size*i), tmp)
        end if
      end do
      call mpi_memory_exchange(mpi_send_recv, A%m, A%blc_size, x)
    end do
  end subroutine mpi_bcsr_parallel_lusgs

  subroutine bcsr_diag(A, x, b)
    use linear_solver_module, only: lu_solve_noswp
    use mpi
    use mpi_module
    implicit none

    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(in) :: b
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(inout) :: x

    integer(kind=ENTIER) :: i, j
    real(kind=DOUBLE), dimension(A%blc_size,A%blc_size) :: tmp

    do i = 1, A%m
      do j = A%row(i), A%row(i + 1) - 1
        if (A%col(j) == i) then
          call lu_solve_noswp(A%blc_size, A%v(:, :, j), &
            x(A%blc_size*(i - 1) + 1:A%blc_size*i), b(A%blc_size*(i - 1) + 1:A%blc_size*i), tmp)
          exit
        end if
      end do
    end do
  end subroutine bcsr_diag

  subroutine bcsr_parallel_lusgs(A, x, b, &
      n_color_sets, color_set, color)
    use sparse_csr_linear_module, only: color_set_type
    use linear_solver_module, only: lu_solve_noswp
    implicit none

    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(in) :: b
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(inout) :: x
    integer(kind=ENTIER), intent(in) :: n_color_sets
    type(color_set_type), dimension(n_color_sets), intent(in) :: color_set
    integer(kind=ENTIER), dimension(A%m), intent(in) :: color

    integer(kind=ENTIER) :: i, j, cs, ics
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: y
    integer(kind=ENTIER), dimension(A%blc_size*A%m) :: d
    real(kind=DOUBLE), dimension(A%blc_size, A%blc_size) :: tmp

    !Get diagonal of A
    do i = 1, A%m
      do j = A%row(i), A%row(i + 1) - 1
        if (A%col(j) == i) then
          d(i) = j
          exit
        end if
      end do
    end do

    !(L+D)D-1(U+D) x = b
    !(L+D) y = b
    !(U+D) x = D y

    do cs = 1, n_color_sets
      do ics = 1, color_set(cs)%n_elems
        i = color_set(cs)%color(ics)
        y(A%blc_size*(i - 1) + 1:A%blc_size*i) = b(A%blc_size*(i - 1) + 1:A%blc_size*i)
        do j = A%row(i), A%row(i + 1) - 1
          if (color(A%col(j)) < cs) then
            !y(i) = y(i) - A%v(j)*y(A%col(j))
            y(A%blc_size*(i - 1) + 1:A%blc_size*i) = y(A%blc_size*(i - 1) + 1:A%blc_size*i) &
              - matmul(A%v(:, :, j), y(A%blc_size*(A%col(j) - 1) + 1:A%blc_size*A%col(j)))
          end if
        end do
      end do
    end do

    do cs = n_color_sets, 1, -1
      do ics = 1, color_set(cs)%n_elems
        i = color_set(cs)%color(ics)
        !x(i) = y(i)
        x(A%blc_size*(i - 1) + 1:A%blc_size*i) = y(A%blc_size*(i - 1) + 1:A%blc_size*i)
        do j = A%row(i), A%row(i + 1) - 1
          if (color(A%col(j)) > cs) then
            ! x(i) = x(i) - A%v(j)*x(A%col(j))
            x(A%blc_size*(i - 1) + 1:A%blc_size*i) = x(A%blc_size*(i - 1) + 1:A%blc_size*i) &
              - matmul(A%v(:, :, j), x(A%blc_size*(A%col(j) - 1) + 1:A%blc_size*A%col(j)))
          end if
        end do
        !x(i) = x(i)/d(i)
        y(A%blc_size*(i - 1) + 1:A%blc_size*i) = x(A%blc_size*(i - 1) + 1:A%blc_size*i)
        call lu_solve_noswp(A%blc_size, A%v(:, :, d(i)), &
          x(A%blc_size*(i - 1) + 1:A%blc_size*i), y(A%blc_size*(i - 1) + 1:A%blc_size*i), tmp)
      end do
    end do
  end subroutine bcsr_parallel_lusgs

  subroutine bcsr_restore_diag(A, diag)
    implicit none

    type(bcsr_matrix_type), intent(inout) :: A
    real(kind=DOUBLE), dimension(A%blc_size, A%blc_size, A%m), intent(in) :: diag

    integer(kind=ENTIER) :: i, j

    do i = 1, A%m
      do j = A%row(i), A%row(i + 1) - 1
        if (A%col(j) == i) then
          A%v(:, :, j) = diag(:, :, i)
          exit
        end if
      end do
    end do
  end subroutine bcsr_restore_diag

  subroutine bcsr_save_diag(A, diag)
    implicit none

    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(A%blc_size, A%blc_size, A%m), intent(inout) :: diag

    integer(kind=ENTIER) :: i, j

    do i = 1, A%m
      do j = A%row(i), A%row(i + 1) - 1
        if (A%col(j) == i) then
          diag(:, :, i) = A%v(:, :, j)
          exit
        end if
      end do
    end do
  end subroutine bcsr_save_diag

  subroutine mpi_modifiedArnoldi_diag(n, m, r, A, H, V, mpi_send_recv)
    use mpi
    use mpi_module
    use sparse_csr_linear_module, only: mpi_dot_product, mpi_norm2
    implicit none

    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), intent(in), dimension(n) :: r
    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(m + 1, m), intent(inout) :: H
    real(kind=DOUBLE), dimension(n, m + 1), intent(inout) :: V
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    integer(kind=ENTIER) :: i, j

    V(:, 1) = r
    do j = 2, m + 1
      call bcsr_diag(A, V(:, j), V(:, j - 1))
      V(:, j) = mpi_bcsr_matmul(A, V(:, j), mpi_send_recv)
      do i = 1, j - 1
        H(i, j - 1) = mpi_dot_product(A%m, A%blc_size, V(:, i), V(:, j), mpi_send_recv)
        V(:, j) = V(:, j) - H(i, j - 1)*V(:, i)
      end do
      H(j, j - 1) = mpi_norm2(A%m, A%blc_size, V(:, j), mpi_send_recv)
      V(:, j) = V(:, j)/H(j, j - 1)
    end do
  end subroutine mpi_modifiedArnoldi_diag

  subroutine modifiedArnoldi_diag(n, m, r, A, H, V)
    implicit none

    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), intent(in), dimension(n) :: r
    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(m + 1, m), intent(inout) :: H
    real(kind=DOUBLE), dimension(n, m + 1), intent(inout) :: V

    integer(kind=ENTIER) :: i, j

    V(:, 1) = r
    do j = 2, m + 1
      call bcsr_diag(A, V(:, j), V(:, j - 1))
      V(:, j) = bcsr_matmul(A, V(:, j))
      do i = 1, j - 1
        H(i, j - 1) = dot_product(V(:, i), V(:, j))
        V(:, j) = V(:, j) - H(i, j - 1)*V(:, i)
      end do
      H(j, j - 1) = norm2(V(:, j))
      V(:, j) = V(:, j)/H(j, j - 1)
    end do
  end subroutine modifiedArnoldi_diag

  subroutine bcsr_gmres_diag(A, x, b, tol, m, &
      conv, H, V, r, y, c, s, z, zp)
    implicit none

    integer(kind=ENTIER), intent(in) :: m
    real(kind=DOUBLE), intent(in) :: tol
    type(bcsr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(in) :: b
    real(kind=DOUBLE), dimension(A%blc_size*A%m), intent(inout) :: x
    logical, intent(inout) :: conv

    integer(kind=ENTIER), parameter :: maxit = 40
    integer(kind=ENTIER) :: i, j, k, l, nr, n

    real(kind=DOUBLE), dimension(m + 1, m) :: H
    real(kind=DOUBLE), dimension(A%blc_size*A%m, m + 1) :: V
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: r
    real(kind=DOUBLE), dimension(m + 1) :: y, c, s
    real(kind=DOUBLE), dimension(A%blc_size*A%m) :: z, zp

    real(kind=DOUBLE) :: delta, rho, tmp

    conv = .false.

    n = A%blc_size*A%m
    x = 0.0_DOUBLE
    r = b - bcsr_matmul(A, x)
    rho = norm2(r)
    if (rho < tol) then
      conv = .true.
      return
    end if

    do j = 1, maxit
      y = 0.0_DOUBLE
      y(1) = rho
      r = r/y(1)
      nr = m

      call modifiedArnoldi_diag(n, m, r, A, H, V)

      !Givens rotations
      do i = 1, m
        do k = 2, i
          tmp = H(k - 1, i)
          H(k - 1, i) = c(k - 1)*H(k - 1, i) + s(k - 1)*H(k, i)
          H(k, i) = -s(k - 1)*tmp + c(k - 1)*H(k, i)
        end do

        delta = sqrt(H(i, i)**2 + H(i + 1, i)*H(i + 1, i))
        c(i) = H(i, i)/delta
        s(i) = H(i + 1, i)/delta

        H(i, i) = c(i)*H(i, i) + s(i)*H(i + 1, i)

        do k = i + 1, m + 1
          H(k, i) = 0.0_DOUBLE
        end do

        y(i + 1) = -s(i)*y(i)
        y(i) = c(i)*y(i)
        rho = abs(y(i + 1))
        if (rho < tol) then
          nr = i + 1
          exit
        end if
      end do

      !Backsolve
      z = 0.0_DOUBLE
      if (i >= m) then
        nr = m
        z(nr) = y(nr)/H(nr, nr)
      end if

      do k = nr - 1, 1, -1
        z(k) = y(k)
        do l = k + 1, nr
          z(k) = z(k) - H(k, l)*z(l)
        end do
        z(k) = z(k)/H(k, k)
      end do

      !Linear combination of basis vector
      r = 0.0_DOUBLE !reuse r to save memory
      do i = 1, nr
        r = r + z(i)*V(:, i)
      end do

      call bcsr_diag(A, zp, r)
      x = x + zp

      r = b - bcsr_matmul(A, x)
      rho = norm2(r)
      if (rho < tol) then
        conv = .true.
        exit
      end if
    end do
  end subroutine bcsr_gmres_diag
end module sparse_bcsr_linear_module
