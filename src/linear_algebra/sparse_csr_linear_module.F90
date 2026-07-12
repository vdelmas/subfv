module sparse_csr_linear_module
  use precision_module
  use mpi_module
  implicit none

  type csr_matrix_type
    integer(kind=ENTIER) :: m, n, nnz !m rows, n columns, nnz non zero elems
    integer(kind=ENTIER), dimension(:), allocatable :: col, row
    real(kind=DOUBLE), dimension(:), allocatable :: v
    type(mpi_send_recv_type) :: mpi_send_recv
  end type csr_matrix_type

  type color_set_type
    integer(kind=ENTIER) :: n_elems
    integer(kind=ENTIER), dimension(:), allocatable :: color
  end type color_set_type
contains
  subroutine modifiedArnoldi(n, m, r, A, H, V)
    implicit none

    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), intent(in), dimension(n) :: r
    type(csr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(m + 1, m), intent(inout) :: H
    real(kind=DOUBLE), dimension(n, m + 1), intent(inout) :: V

    integer(kind=ENTIER) :: i, j

    H = 0.0_DOUBLE
    V = 0.0_DOUBLE

    V(:, 1) = r
    do j = 2, m + 1
      V(:, j) = csr_matmul(A, V(:, j - 1))
      do i = 1, j - 1
        H(i, j - 1) = dot_product(V(:, i), V(:, j))
        V(:, j) = V(:, j) - H(i, j - 1)*V(:, i)
      end do
      H(j, j - 1) = norm2(V(:, j))
      V(:, j) = V(:, j)/H(j, j - 1)
    end do
  end subroutine modifiedArnoldi

  subroutine modifiedArnoldi_parallel_lusgs(n, m, r, A, H, V, &
      n_color_sets, color_set, color)
    implicit none

    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), intent(in), dimension(n) :: r
    type(csr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(m + 1, m), intent(inout) :: H
    real(kind=DOUBLE), dimension(n, m + 1), intent(inout) :: V
    integer(kind=ENTIER), intent(in) :: n_color_sets
    type(color_set_type), dimension(n_color_sets), intent(in) :: color_set
    integer(kind=ENTIER), dimension(n), intent(in) :: color

    integer(kind=ENTIER) :: i, j

    H = 0.0_DOUBLE
    V = 0.0_DOUBLE

    V(:, 1) = r
    do j = 2, m + 1
      !V(:, j) = csr_matmul(A, V(:, j - 1))
      call csr_parallel_lusgs(n, A, V(:, j), V(:, j - 1), &
        n_color_sets, color_set, color)
      V(:, j) = csr_matmul(A, V(:, j))
      do i = 1, j - 1
        H(i, j - 1) = dot_product(V(:, i), V(:, j))
        V(:, j) = V(:, j) - H(i, j - 1)*V(:, i)
      end do
      H(j, j - 1) = norm2(V(:, j))
      V(:, j) = V(:, j)/H(j, j - 1)
    end do
  end subroutine modifiedArnoldi_parallel_lusgs

  subroutine modifiedArnoldi_lusgs(n, m, r, A, H, V)
    implicit none

    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), intent(in), dimension(n) :: r
    type(csr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(m + 1, m), intent(inout) :: H
    real(kind=DOUBLE), dimension(n, m + 1), intent(inout) :: V

    integer(kind=ENTIER) :: i, j

    H = 0.0_DOUBLE
    V = 0.0_DOUBLE

    V(:, 1) = r
    do j = 2, m + 1
      !V(:, j) = csr_matmul(A, V(:, j - 1))
      call csr_lusgs(n, A, V(:, j), V(:, j - 1))
      V(:, j) = csr_matmul(A, V(:, j))
      do i = 1, j - 1
        H(i, j - 1) = dot_product(V(:, i), V(:, j))
        V(:, j) = V(:, j) - H(i, j - 1)*V(:, i)
      end do
      H(j, j - 1) = norm2(V(:, j))
      V(:, j) = V(:, j)/H(j, j - 1)
    end do
  end subroutine modifiedArnoldi_lusgs

  subroutine csr_gmres_lusgs(n, A, x, b, tol, m, conv, H, V, r, y, c, s, z, zp)
    implicit none

    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), intent(in) :: tol
    type(csr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(n), intent(in) :: b
    real(kind=DOUBLE), dimension(n), intent(inout) :: x
    logical, intent(inout) :: conv

    integer(kind=ENTIER), parameter :: maxit = 20
    integer(kind=ENTIER) :: i, j, k, l, nr

    real(kind=DOUBLE), dimension(m + 1, m) :: H
    real(kind=DOUBLE), dimension(n, m + 1) :: V
    real(kind=DOUBLE), dimension(n) :: r
    real(kind=DOUBLE), dimension(m + 1) :: y, c, s
    real(kind=DOUBLE), dimension(n) :: z, zp

    real(kind=DOUBLE) :: delta, rho, tmp

    x = 0.0_DOUBLE
    r = b - csr_matmul(A, x)
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

      call csr_lusgs(n, A, zp, r)
      x = x + zp

      r = b - csr_matmul(A, x)
      rho = norm2(r)
      if (rho < tol) exit
    end do

    if (norm2(b - csr_matmul(A, x)) > tol) then
      ! print*, "GMRES DID NOT FULLY CONVERGE", norm2(b - csr_matmul(A, x))
      conv = .false.
    else
      conv = .true.
    end if
  end subroutine csr_gmres_lusgs

  subroutine csr_gmres_parallel_lusgs(n, A, x, b, tol, m, &
      n_color_sets, color_set, color, conv, &
      H, V, r, y, c, s, z, zp)
    implicit none

    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), intent(in) :: tol
    type(csr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(n), intent(in) :: b
    real(kind=DOUBLE), dimension(n), intent(inout) :: x
    logical, intent(inout) :: conv
    integer(kind=ENTIER), intent(in) :: n_color_sets
    type(color_set_type), dimension(n_color_sets), intent(in) :: color_set
    integer(kind=ENTIER), dimension(n), intent(in) :: color

    integer(kind=ENTIER), parameter :: maxit = 40
    integer(kind=ENTIER) :: i, j, k, l, nr

    real(kind=DOUBLE), dimension(m + 1, m) :: H
    real(kind=DOUBLE), dimension(n, m + 1) :: V
    real(kind=DOUBLE), dimension(n) :: r
    real(kind=DOUBLE), dimension(m + 1) :: y, c, s
    real(kind=DOUBLE), dimension(n) :: z, zp

    real(kind=DOUBLE) :: delta, rho, tmp

    x = 0.0_DOUBLE
    r = b - csr_matmul(A, x)
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

      call csr_parallel_lusgs(n, A, zp, r, &
        n_color_sets, color_set, color)
      x = x + zp

      r = b - csr_matmul(A, x)
      rho = norm2(r)
      if (rho < tol) exit
    end do

    if (norm2(b - csr_matmul(A, x)) > tol) then
      conv = .false.
    else
      conv = .true.
    end if
  end subroutine csr_gmres_parallel_lusgs

  subroutine csr_gmres(n, A, x, b, tol, m, conv, H, V, r, y, c, s, z)
    implicit none

    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), intent(in) :: tol
    type(csr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(n), intent(in) :: b
    real(kind=DOUBLE), dimension(n), intent(inout) :: x
    logical, intent(inout) :: conv

    integer(kind=ENTIER), parameter :: maxit = 40
    integer(kind=ENTIER) :: i, j, k, l, nr

    real(kind=DOUBLE), dimension(m + 1, m) :: H
    real(kind=DOUBLE), dimension(n, m + 1) :: V
    real(kind=DOUBLE), dimension(n) :: r
    real(kind=DOUBLE), dimension(m + 1) :: y, c, s
    real(kind=DOUBLE), dimension(n) :: z

    real(kind=DOUBLE) :: delta, rho, tmp

    x = 0.0_DOUBLE
    r = b - csr_matmul(A, x)
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

      ! if (rho < tol) exit

      r = b - csr_matmul(A, x)

      rho = norm2(r)
      if (rho < tol) exit
    end do

    if (norm2(b - csr_matmul(A, x)) > tol) then
      ! print*, "GMRES DID NOT FULLY CONVERGE", norm2(b - csr_matmul(A, x))
      conv = .false.
    else
      conv = .true.
    end if
  end subroutine csr_gmres

  subroutine csr_pcgm_cholesky(n, A, L, x, b, tol)
    implicit none

    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), intent(in) :: tol
    type(csr_matrix_type), intent(in) :: A, L
    real(kind=DOUBLE), dimension(n), intent(in) :: b
    real(kind=DOUBLE), dimension(n), intent(inout) :: x

    integer(kind=ENTIER) :: k
    real(kind=DOUBLE) :: alpha, rdr, rpdrp
    real(kind=DOUBLE), dimension(n) :: r, p, Ap, mr

    r = b - csr_matmul(A, x)
    rpdrp = 0.0_DOUBLE

    call csr_solve_cholesky(L, mr, r)

    rdr = dot_product(mr, r)
    if (sqrt(rdr) < tol) then
      return
    end if
    p = mr

    do k = 1, n
      Ap = csr_matmul(A, p)
      alpha = rdr/dot_product(p, Ap)
      x = x + alpha*p
      r = r - alpha*Ap
      call csr_solve_cholesky(L, mr, r)
      rpdrp = dot_product(mr, r)
      if (sqrt(rpdrp) <= tol) then
        if (norm2(r) <= tol) exit
      end if
      p = mr + (rpdrp/rdr)*p
      rdr = rpdrp
    end do

    if (sqrt(rpdrp) > tol) then
      print *, "PCGM did not converge fully ! Redisu = ", sqrt(rpdrp)
      error stop
    end if
  end subroutine csr_pcgm_cholesky

  subroutine csr_pcgm(n, A, x, b, tol)
    implicit none

    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), intent(in) :: tol
    type(csr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(n), intent(in) :: b
    real(kind=DOUBLE), dimension(n), intent(inout) :: x

    integer(kind=ENTIER) :: k, i
    real(kind=DOUBLE) :: alpha, rdr, rpdrp
    real(kind=DOUBLE), dimension(n) :: r, p, Ap, mr, d

    alpha = 0.0_DOUBLE
    rdr = 0.0_DOUBLE
    rpdrp = 0.0_DOUBLE

    !Get diagonal of A
    do i = 1, A%m
      do k = A%row(i), A%row(i + 1) - 1
        if (A%col(k) == i) then
          d(i) = A%v(k)
          exit
        end if
      end do
    end do

    r = b - csr_matmul(A, x)
    do i = 1, n
      mr(i) = r(i)/d(i)
    end do
    rdr = dot_product(mr, r)
    if (sqrt(rdr) < tol) then
      ! print *, "ALREADY CONVERGED AT START", sqrt(rdr)
      return
    end if
    p = mr

    do k = 1, n
      Ap = csr_matmul(A, p)
      alpha = rdr/dot_product(p, Ap)
      x = x + alpha*p
      r = r - alpha*Ap
      do i = 1, n
        mr(i) = r(i)/d(i)
      end do
      rpdrp = dot_product(mr, r)
      ! print*, "Residu ", n, k, sqrt(rpdrp)
      if (sqrt(rpdrp) <= tol) then
        if (norm2(r) <= tol) exit
      end if
      p = mr + (rpdrp/rdr)*p
      rdr = rpdrp
    end do

    if (sqrt(rpdrp) > tol) then
      print *, "PCGM did not converge fully ! Redisu = ", sqrt(rpdrp)
      error stop
    end if
  end subroutine csr_pcgm

  function csr_matmul(A, x) result(y)
    implicit none

    type(csr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(A%n), intent(in) :: x
    real(kind=DOUBLE), dimension(A%m) :: y

    integer(kind=ENTIER) :: i, j

    do i = 1, A%m
      y(i) = 0.0_DOUBLE
      do j = A%row(i), A%row(i + 1) - 1
        y(i) = y(i) + A%v(j)*x(A%col(j))
      end do
    end do
  end function csr_matmul

  subroutine mpi_csr_pcgm(n, A, x, b, tol, mpi_send_recv)
    use mpi
    use mpi_module
    implicit none

    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), intent(in) :: tol
    type(csr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(n), intent(in) :: b
    real(kind=DOUBLE), dimension(n), intent(inout) :: x
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    integer(kind=ENTIER), parameter :: maxit = 1000

    integer(kind=ENTIER) :: k, i
    real(kind=DOUBLE) :: alpha, rdr, rpdrp
    real(kind=DOUBLE), dimension(n) :: r, p, Ap, mr, d

    alpha = 0.0_DOUBLE
    rdr = 0.0_DOUBLE
    rpdrp = 0.0_DOUBLE

    !Get diagonal of A
    do i = 1, A%m
      if( .not. mpi_send_recv%is_ghost(i) ) then
        do k = A%row(i), A%row(i + 1) - 1
          if (A%col(k) == i) then
            d(i) = A%v(k)
            exit
          end if
        end do
      end if
    end do

    r = b - mpi_csr_matmul(A, x, mpi_send_recv)

    mr = 0.0_DOUBLE
    do i = 1, n
      if( .not. mpi_send_recv%is_ghost(i) ) mr(i) = r(i)/d(i)
    end do

    rdr = mpi_dot_product(A%m, 1, mr, r, mpi_send_recv)

    if (sqrt(rdr) < tol) return

    p = mr
    do k = 1, maxit
      Ap = mpi_csr_matmul(A, p, mpi_send_recv)
      alpha = rdr/mpi_dot_product(A%m, 1, p, Ap, mpi_send_recv)
      x = x + alpha*p
      r = r - alpha*Ap
      do i = 1, n
        if( .not. mpi_send_recv%is_ghost(i) ) mr(i) = r(i)/d(i)
      end do
      rpdrp = mpi_dot_product(A%m, 1, mr, r, mpi_send_recv)
      !print*, "Residu ", n, k, sqrt(rpdrp), rpdrp/rdr
      if (sqrt(rpdrp) <= tol) then
        if (mpi_norm2(A%m, 1, r, mpi_send_recv) <= tol) exit
      end if
      p = mr + (rpdrp/rdr)*p
      rdr = rpdrp
    end do

    if (sqrt(rpdrp) > tol) then
      print *, "PCGM did not converge fully ! Redisu = ", sqrt(rpdrp)
      error stop
    end if
  end subroutine mpi_csr_pcgm

  function mpi_norm2(n, blc_size, a, mpi_send_recv)
    use mpi
    use mpi_module
    implicit none

    integer(kind=ENTIER), intent(in) :: n, blc_size
    real(kind=DOUBLE), dimension(n*blc_size), intent(in) :: a
    type(mpi_send_recv_type), intent(in) :: mpi_send_recv
    real(kind=DOUBLE) :: mpi_norm2

    integer(kind=ENTIER) :: i, mpi_ierr

    mpi_norm2 = 0.0_DOUBLE
    do i = 1, n
      if( .not. mpi_send_recv%is_ghost(i) ) then
        mpi_norm2 = mpi_norm2 &
          + dot_product(a((i-1)*blc_size+1:i*blc_size),a((i-1)*blc_size+1:i*blc_size))
      end if
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE, mpi_norm2, 1, &
      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
    mpi_norm2 = sqrt(mpi_norm2)
  end function mpi_norm2

  function mpi_dot_product(n, blc_size, a, b, mpi_send_recv)
    use mpi
    use mpi_module
    implicit none

    integer(kind=ENTIER), intent(in) :: n, blc_size
    real(kind=DOUBLE), dimension(n*blc_size), intent(in) :: a, b
    type(mpi_send_recv_type), intent(in) :: mpi_send_recv
    real(kind=DOUBLE) :: mpi_dot_product

    integer(kind=ENTIER) :: i, mpi_ierr

    mpi_dot_product = 0.0_DOUBLE
    do i = 1, n
      if( .not. mpi_send_recv%is_ghost(i) ) then
        mpi_dot_product = mpi_dot_product &
          + dot_product(a((i-1)*blc_size+1:i*blc_size),b((i-1)*blc_size+1:i*blc_size))
      end if
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE, mpi_dot_product, 1, &
      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
  end function mpi_dot_product

  function mpi_csr_matmul(A, x, mpi_send_recv) result(y)
    implicit none

    type(csr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(A%n), intent(inout) :: x
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv
    real(kind=DOUBLE), dimension(A%m) :: y

    integer(kind=ENTIER) :: i, j

    call mpi_memory_exchange(mpi_send_recv, A%m, 1, x)

    y = 0.0_DOUBLE
    do i = 1, A%m
      if( .not. mpi_send_recv%is_ghost(i) ) then
        do j = A%row(i), A%row(i + 1) - 1
          y(i) = y(i) + A%v(j)*x(A%col(j))
        end do
      end if
    end do
  end function mpi_csr_matmul

  subroutine csr_solve_cholesky(L, x, b)
    implicit none

    type(csr_matrix_type), intent(in) :: L
    real(kind=DOUBLE), dimension(L%n) :: x
    real(kind=DOUBLE), dimension(L%m) :: b

    integer(kind=ENTIER) :: i, j, k
    real(kind=DOUBLE), dimension(L%n) :: y
    real(kind=DOUBLE) :: s, diag

    if (L%n /= L%m) then
      print *, 'Bad sizes csr_solve_cholesky'
      error stop
    end if

    diag = 0.0_DOUBLE
    y = 0.0_DOUBLE
    !L y = b
    do i = 1, L%n
      s = 0.0_DOUBLE
      do k = L%row(i), L%row(i + 1) - 1
        j = L%col(k)
        if (j < i) then
          s = s + L%v(k)*y(j)
        else if (j == i) then
          diag = L%v(k)
          exit
        end if
      end do
      y(i) = (b(i) - s)/diag
    end do

    !L^t x = y
    x = y
    do i = L%n, 1, -1
      do k = L%row(i + 1) - 1, L%row(i), -1
        j = L%col(k)
        if (j < i) then
          x(j) = x(j) - L%v(k)*x(i)
        else if (j == i) then
          x(i) = x(i)/L%v(k)
        end if
      end do
    end do
  end subroutine csr_solve_cholesky

  subroutine csr_incomplete_cholesky(A, L)
    implicit none

    type(csr_matrix_type), intent(in) :: A
    type(csr_matrix_type), intent(inout) :: L

    integer(kind=ENTIER) :: i, j, k, i1, i2
    real(kind=DOUBLE) :: s
    integer(kind=ENTIER), dimension(A%m) :: diag

    if (A%m /= A%n) then
      print *, "BAD SIZES CSR INCOMP CHOL"
      error stop
    end if
    diag = 0.0_DOUBLE

    !Bsed on ILUPP python library
    L%v(:) = 0.0_DOUBLE
    do i = 1, A%m
      do k = L%row(i), L%row(i + 1) - 1
        j = L%col(k)
        if (j <= i) then

          i1 = L%row(i)
          i2 = L%row(j)
          s = 0.0_DOUBLE
          do while (L%col(i1) < i .and. &
              L%col(i2) < i)
            if (L%col(i1) == L%col(i2)) then
              s = s + L%v(i1)*L%v(i2)
              i1 = i1 + 1
              i2 = i2 + 1
            else if (L%col(i1) < L%col(i2)) then
              i1 = i1 + 1
            else
              i2 = i2 + 1
            end if
          end do

          if (j == i) then
            L%v(k) = sqrt(A%v(k) - s)
            diag(i) = k
          else if (j < i) then
            L%v(k) = (A%v(k) - s)/L%v(diag(j))
          end if
        else
          exit
        end if
      end do
    end do
  end subroutine csr_incomplete_cholesky

  subroutine print_csr_matrix(mat)
    implicit none

    type(csr_matrix_type), intent(in) :: mat

    integer(kind=ENTIER) :: i, j, k, nnz, fu

    open (newunit=fu, file="mat.dat", status="unknown")

    print *, "----------"
    nnz = 1
    do i = 1, mat%n
      do j = mat%row(i), mat%row(i + 1) - 1
        if (j == mat%row(i)) then
          do k = 1, mat%col(nnz) - 1
            write (fu, "(f18.2)", advance="no") 0.0_DOUBLE
            write (*, "(f18.2)", advance="no") 0.0_DOUBLE
          end do
        else
          do k = mat%col(nnz - 1), mat%col(nnz) - 2
            write (fu, "(f18.2)", advance="no") 0.0_DOUBLE
            write (*, "(f18.2)", advance="no") 0.0_DOUBLE
          end do
        end if
        write (fu, "(f18.6)", advance="no") mat%v(nnz)
        write (*, "(f18.6)", advance="no") mat%v(nnz)
        nnz = nnz + 1
      end do
      do k = mat%col(nnz - 1), mat%m - 1
        write (fu, "(f18.2)", advance="no") 0.0_DOUBLE
        write (*, "(f18.2)", advance="no") 0.0_DOUBLE
      end do
      write (fu, *) ""
      write (*, *) ""
    end do
    print *, "----------"

    close (fu)
  end subroutine print_csr_matrix

  subroutine csr_add_elem(mat, i, j, v)
    implicit none

    type(csr_matrix_type), intent(inout) :: mat
    integer(kind=ENTIER), intent(in) :: i, j
    real(kind=DOUBLE) :: v

    integer(kind=ENTIER) :: m

    do m = mat%row(i), mat%row(i + 1) - 1
      if (mat%col(m) == j) then
        mat%v(m) = mat%v(m) + v
      end if
    end do
  end subroutine csr_add_elem

  subroutine csr_lusgs(n, A, x, b)
    implicit none

    integer(kind=ENTIER), intent(in) :: n
    type(csr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(n), intent(in) :: b
    real(kind=DOUBLE), dimension(n), intent(inout) :: x

    integer(kind=ENTIER) :: i, j
    real(kind=DOUBLE), dimension(n) :: y, d

    !Get diagonal of A
    do i = 1, n
      do j = A%row(i), A%row(i + 1) - 1
        if (A%col(j) == i) then
          d(i) = A%v(j)
          exit
        end if
      end do
    end do

    !(L+D)D-1(U+D) x = b
    !(L+D) D-1 y = b
    !(U+D) x = y
    do i = 1, n
      y(i) = b(i)
      j = A%row(i)
      do while (A%col(j) < i)
        y(i) = y(i) - A%v(j)*y(A%col(j))/d(i)
        j = j + 1
      end do
    end do

    do i = 1, n
      x(i) = y(i)
      j = A%row(i)
      do while (A%col(j) < i)
        x(i) = x(i) - A%v(j)*x(A%col(j))
        j = j + 1
      end do
      x(i) = x(i)/A%v(j)
    end do
  end subroutine csr_lusgs

  subroutine csr_parallel_lusgs(n, A, x, b, &
      n_color_sets, color_set, color)
    implicit none

    integer(kind=ENTIER), intent(in) :: n
    type(csr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(n), intent(in) :: b
    real(kind=DOUBLE), dimension(n), intent(inout) :: x
    integer(kind=ENTIER), intent(in) :: n_color_sets
    type(color_set_type), dimension(n_color_sets), intent(in) :: color_set
    integer(kind=ENTIER), dimension(n), intent(in) :: color

    integer(kind=ENTIER), parameter :: blc_size = 5

    integer(kind=ENTIER) :: i, j, cs, ics, k
    real(kind=DOUBLE), dimension(n) :: y, d

    !Get diagonal of A
    do i = 1, n
      do j = A%row(i), A%row(i + 1) - 1
        if (A%col(j) == i) then
          d(i) = A%v(j)
          exit
        end if
      end do
    end do

    !(L+D)D-1(U+D) x = b
    !(L+D) D-1 y = b
    !(U+D) x = y

    do cs = 1, n_color_sets
      do ics = 1, color_set(cs)%n_elems
        i = color_set(cs)%color(ics)
        do k = 1 + blc_size*(i - 1), blc_size*i
          y(k) = b(k)
          do j = A%row(k), A%row(k + 1) - 1
            if (color(ceiling(real(A%col(j))/real(blc_size))) < cs) then
              y(k) = y(k) - A%v(j)*y(A%col(j))/d(k)
            end if
          end do
        end do
      end do
    end do

    do cs = n_color_sets, 1, -1
      do ics = 1, color_set(cs)%n_elems
        i = color_set(cs)%color(ics)
        do k = 1 + blc_size*(i - 1), blc_size*i
          x(k) = y(k)
          do j = A%row(k), A%row(k + 1) - 1
            if (color(ceiling(real(A%col(j))/real(blc_size))) > cs) then
              x(k) = x(k) - A%v(j)*x(A%col(j))
            end if
          end do
          x(k) = x(k)/d(k)
        end do
      end do
    end do
  end subroutine csr_parallel_lusgs

  subroutine csr_parallel_gauss_seidel(n, A, x, b, tol, &
      n_color_sets, color_set, color, conv)
    implicit none

    integer(kind=ENTIER), intent(in) :: n
    type(csr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(n), intent(in) :: b
    real(kind=DOUBLE), dimension(n), intent(inout) :: x
    real(kind=DOUBLE), intent(in) :: tol
    integer(kind=ENTIER), intent(in) :: n_color_sets
    type(color_set_type), dimension(n_color_sets), intent(in) :: color_set
    integer(kind=ENTIER), dimension(n), intent(in) :: color
    logical, intent(inout) :: conv

    integer(kind=ENTIER), parameter :: maxit = 20

    integer(kind=ENTIER) :: i, j, p, cs, ics
    real(kind=DOUBLE) :: res, diag
    real(kind=DOUBLE), dimension(n) :: b2

    res = norm2(b - csr_matmul(A, x))
    if (res < tol) return

    diag = 0.0_DOUBLE

    do p = 1, maxit
      !(L+D) x = -U x + b
      !(L+D) x = b2
      !b2 = b - csr_matmul_upper(A, x)

      do cs = 1, n_color_sets
        do ics = 1, color_set(cs)%n_elems
          i = color_set(cs)%color(ics)
          b2(i) = b(i)
          do j = A%row(i), A%row(i + 1) - 1
            if (color(A%col(j)) > cs) then
              b2(i) = b2(i) - A%v(j)*x(A%col(j))
            end if
          end do
        end do
      end do

      do cs = 1, n_color_sets
        do ics = 1, color_set(cs)%n_elems
          i = color_set(cs)%color(ics)
          x(i) = b2(i)
          do j = A%row(i), A%row(i + 1) - 1
            if (color(A%col(j)) < cs) then
              x(i) = x(i) - A%v(j)*x(A%col(j))
            else if (i == A%col(j)) then
              diag = A%v(j)
            end if
          end do
          x(i) = x(i)/diag
        end do
      end do

      res = norm2(b - csr_matmul(A, x))
      if (res < tol) then
        conv = .true.
        return
      end if
    end do

    res = norm2(b - csr_matmul(A, x))
    if (res > tol) then
      print *, "GS DID NOT CONV"
      conv = .false.
    else
      conv = .true.
    end if
  end subroutine csr_parallel_gauss_seidel

  subroutine csr_restore_diag(A, diag)
    implicit none

    type(csr_matrix_type), intent(inout) :: A
    real(kind=DOUBLE), dimension(A%m), intent(in) :: diag

    integer(kind=ENTIER) :: i, j

    do i = 1, A%m
      do j = A%row(i), A%row(i + 1) - 1
        if (A%col(j) == i) then
          A%v(j) = diag(i)
          exit
        end if
      end do
    end do
  end subroutine csr_restore_diag

  subroutine csr_save_diag(A, diag)
    implicit none

    type(csr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(A%m), intent(inout) :: diag

    integer(kind=ENTIER) :: i, j

    do i = 1, A%m
      do j = A%row(i), A%row(i + 1) - 1
        if (A%col(j) == i) then
          diag(i) = A%v(j)
          exit
        end if
      end do
    end do
  end subroutine csr_save_diag
end module sparse_csr_linear_module
