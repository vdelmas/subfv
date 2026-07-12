module sparse_linear_module
  use precision_module
  implicit none

  type csr_matrix_type
    integer(kind=ENTIER) :: m, n
    integer(kind=ENTIER_D) :: nnz !m rows, n columns, nnz non zero elems
    integer(kind=ENTIER), dimension(:), allocatable :: col
    integer(kind=ENTIER_D), dimension(:), allocatable :: row
    real(kind=DOUBLE), dimension(:), allocatable :: v
  end type csr_matrix_type

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

  subroutine csr_gmres(n, A, x, b, tol, m)
    implicit none

    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), intent(in) :: tol
    type(csr_matrix_type), intent(in) :: A
    real(kind=DOUBLE), dimension(n), intent(in) :: b
    real(kind=DOUBLE), dimension(n), intent(inout) :: x

    integer(kind=ENTIER), parameter :: maxit = 100
    integer(kind=ENTIER) :: i, j, k, l, nr

    real(kind=DOUBLE), dimension(m + 1, m) :: H
    real(kind=DOUBLE), dimension(n, m + 1) :: V
    real(kind=DOUBLE), dimension(n) :: r
    real(kind=DOUBLE), dimension(m + 1) :: y, c, s
    real(kind=DOUBLE), dimension(n) :: z

    real(kind=DOUBLE) :: delta, rho, tmp

    x = 0.0_DOUBLE

    r = b - csr_matmul(A, x)

    do j = 1, maxit
      y = 0.0_DOUBLE
      y(1) = norm2(r)
      r = r/y(1)

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

      ! print*, j, rho
    end do

    if (norm2(b - csr_matmul(A, x)) > tol) then
      print *, "GMRES DID NOT FULLY CONVERGE", norm2(b - csr_matmul(A, x))

      ! do i = 1, A%m
      !   do j = A%row(i), A%row(i+1) - 1
      !     write(1, *) i-1, A%col(j)-1, A%v(j)
      !   end do
      ! end do

      ! do i = 1, A%m
      !   write(2, *) b(i)
      ! end do

      ! error stop
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

  subroutine csr_solve_cholesky(L, x, b)
    implicit none

    type(csr_matrix_type), intent(in) :: L
    real(kind=DOUBLE), dimension(L%n) :: x
    real(kind=DOUBLE), dimension(L%m) :: b

    integer(kind=ENTIER) :: i, j, kdiag, k
    real(kind=DOUBLE), dimension(L%n) :: y
    real(kind=DOUBLE) :: s

    if (L%n /= L%m) then
      print *, 'Bad sizes csr_solve_cholesky'
      error stop
    end if

    y = 0.0_DOUBLE
    !L y = b
    do i = 1, L%n
      s = 0.0_DOUBLE
      do k = L%row(i), L%row(i + 1) - 1
        j = L%col(k)
        if (j < i) then
          s = s + L%v(k)*y(j)
        else if (j == i) then
          kdiag = k
          exit
        end if
      end do
      y(i) = (b(i) - s)/L%v(kdiag)
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

    integer(kind=ENTIER) :: i, k, fu
    integer(kind=ENTIER_D) :: nnz, j

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

  subroutine add_block_to_csr(mat, i, j, di, dj, M)
    implicit none

    type(csr_matrix_type), intent(inout) :: mat
    integer(kind=ENTIER), intent(in) :: i, di, j, dj
    real(kind=DOUBLE), dimension(di, dj) :: M

    integer(kind=ENTIER) :: k, p

    do p = 1, dj
      do k = 1, di
        ! N(i+k-1,j+p-1) = N(i+k-1,j+p-1) + M(k,p)
        call csr_add_elem(mat, i + k - 1, j + p - 1, M(k, p))
      end do
    end do
  end subroutine add_block_to_csr

  subroutine csr_add_elem(mat, i, j, v)
    implicit none

    type(csr_matrix_type), intent(inout) :: mat
    integer(kind=ENTIER), intent(in) :: i, j
    real(kind=DOUBLE) :: v

    integer(kind=ENTIER_D) :: m
    logical :: found

    found = .false.
    do m = mat%row(i), mat%row(i + 1) - 1
      if (mat%col(m) == j) then
        mat%v(m) = mat%v(m) + v
        found = .true.
      end if
    end do

    if (.not. found) then
      print *, "BAD CSR FORMAT ADD NOT FOUND!!!!", i, j
    end if
  end subroutine csr_add_elem
end module sparse_linear_module
