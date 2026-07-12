module linear_solver_module
  use precision_module
  implicit none

  real(kind=DOUBLE), dimension(3, 3), parameter :: eye3 = reshape( &
    [1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, &
     0.0_DOUBLE, 1.0_DOUBLE, 0.0_DOUBLE, &
     0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE], [3, 3])

  real(kind=DOUBLE), dimension(5, 5), parameter :: eye5 = reshape( &
    [1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, &
     0.0_DOUBLE, 1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, &
     0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, &
     0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE, 0.0_DOUBLE, &
     0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE], [5, 5])

contains
  subroutine inverse(mat, mat_inv, n)
    implicit none
    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), dimension(n, n), intent(in) :: mat
    real(kind=DOUBLE), dimension(n, n), intent(inout) :: mat_inv

    integer(kind=ENTIER) :: i, j, k
    integer(kind=ENTIER) :: iidx, idx_tmp
    integer(kind=ENTIER), dimension(n) :: idx
    real(kind=DOUBLE) :: coeff
    real(kind=DOUBLE), dimension(n) :: b, d, x
    real(kind=DOUBLE), dimension(n, n) :: a, L, U

    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    L = 0.0_DOUBLE
    U = 0.0_DOUBLE
    b = 0.0_DOUBLE
    a = mat

    do k = 1, n
      idx(k) = k
    end do

    ! step 1: forward elimination
    do k = 1, n - 1

      !!Pivot
      iidx = 1
      do while (abs(a(idx(k), k)) < 1e-12_DOUBLE)
        idx_tmp = idx(k)
        idx(k) = idx(k + iidx)
        idx(k + iidx) = idx_tmp
        do j = 1, k - 1
          coeff = L(idx(k), j)
          L(idx(k), j) = L(idx(k + iidx), j)
          L(idx(k + iidx), j) = coeff
        end do
        iidx = iidx + 1
        if (k + iidx > n) then
          print *, "LU PIVOT NOT WORKING"
          error stop
        end if
      end do

      do i = k + 1, n
        coeff = a(idx(i), k)/a(idx(k), k)
        L(i, k) = coeff
        do j = k + 1, n
          a(idx(i), j) = a(idx(i), j) - coeff*a(idx(k), j)
        end do
      end do
    end do

    ! Step 2: prepare L and U matrices
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i = 1, n
      L(i, i) = 1.0_DOUBLE
    end do
    ! U matrix is the upper triangular part of A
    do j = 1, n
      if (abs(a(j, j)) > 1e-12) then
        do i = 1, j
          U(i, j) = a(idx(i), j)
        end do
      else
        U(i, j) = 1.0_DOUBLE
      end if
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k = 1, n

      do iidx = 1, n
        if (idx(iidx) == k) then
          b(iidx) = 1.0_DOUBLE
        else
          b(iidx) = 0.0_DOUBLE
        end if
      end do

      d(1) = b(1)
      ! Step 3a: Solve Ld=b using the forward substitution
      do i = 2, n
        d(i) = b(i)
        do j = 1, i - 1
          d(i) = d(i) - L(i, j)*d(j)
        end do
      end do
      ! Step 3b: Solve Ux=d using the back substitution
      x(n) = d(n)/U(n, n)
      do i = n - 1, 1, -1
        x(i) = d(i)
        do j = n, i + 1, -1
          x(i) = x(i) - U(i, j)*x(j)
        end do
        x(i) = x(i)/u(i, i)
      end do
      ! Step 3c: fill the solutions x(n) into column k of C
      do i = 1, n
        mat_inv(i, k) = x(i)
      end do
      b(k) = 0.0
    end do
  end subroutine inverse

  pure function det(a, b, c, d)
    implicit none

    real(kind=DOUBLE), intent(in) :: a, b, c, d
    real(kind=DOUBLE) :: det

    det = a*d - c*b
  end function det

  pure subroutine inverse_3_by_3(mat, mat_inv)
    implicit none

    real(kind=DOUBLE), dimension(3, 3), intent(in) :: mat
    real(kind=DOUBLE), dimension(3, 3), intent(inout) :: mat_inv

    real(kind=DOUBLE) :: det_mat
    real(kind=DOUBLE), dimension(3, 3) :: mat_tmp, mat_adj

    mat_tmp = transpose(mat)

    mat_adj(1, 1) = det(mat_tmp(2, 2), mat_tmp(2, 3), &
      mat_tmp(3, 2), mat_tmp(3, 3))
    mat_adj(1, 2) = -det(mat_tmp(2, 1), mat_tmp(2, 3), &
      mat_tmp(3, 1), mat_tmp(3, 3))
    mat_adj(1, 3) = det(mat_tmp(2, 1), mat_tmp(2, 2), &
      mat_tmp(3, 1), mat_tmp(3, 2))

    mat_adj(2, 1) = -det(mat_tmp(1, 2), mat_tmp(1, 3), &
      mat_tmp(3, 2), mat_tmp(3, 3))
    mat_adj(2, 2) = det(mat_tmp(1, 1), mat_tmp(1, 3), &
      mat_tmp(3, 1), mat_tmp(3, 3))
    mat_adj(2, 3) = -det(mat_tmp(1, 1), mat_tmp(1, 2), &
      mat_tmp(3, 1), mat_tmp(3, 2))

    mat_adj(3, 1) = det(mat_tmp(1, 2), mat_tmp(1, 3), &
      mat_tmp(2, 2), mat_tmp(2, 3))
    mat_adj(3, 2) = -det(mat_tmp(1, 1), mat_tmp(1, 3), &
      mat_tmp(2, 1), mat_tmp(2, 3))
    mat_adj(3, 3) = det(mat_tmp(1, 1), mat_tmp(1, 2), &
      mat_tmp(2, 1), mat_tmp(2, 2))

    det_mat = &
      mat(1, 1)*mat_adj(1, 1) &
      + mat(1, 2)*mat_adj(2, 1) &
      + mat(1, 3)*mat_adj(3, 1)

    mat_inv = 1.0_DOUBLE/det_mat*mat_adj
  end subroutine inverse_3_by_3

  pure function tensor_product_3(a, b) result(tp)
    implicit none

    real(kind=DOUBLE), intent(in) :: a(3), b(3)
    real(kind=DOUBLE), dimension(3, 3) :: tp

    tp(:, 1) = a(:)*b(1)
    tp(:, 2) = a(:)*b(2)
    tp(:, 3) = a(:)*b(3)
  end function tensor_product_3

  pure function tensor_product(a, b) result(tp)
    implicit none

    real(kind=DOUBLE), intent(in) :: a(:), b(:)
    real(kind=DOUBLE), dimension(size(a), size(b)) :: tp

    integer(kind=ENTIER) :: j

    do j = 1, size(b)
      tp(:, j) = a(:)*b(j)
    end do
  end function tensor_product

  pure function cross_product(a, b)
    implicit none

    real(kind=DOUBLE), intent(in) :: a(3), b(3)
    real(kind=DOUBLE), dimension(3) :: cross_product

    cross_product(1) = a(2)*b(3) - a(3)*b(2)
    cross_product(2) = -(a(1)*b(3) - a(3)*b(1))
    cross_product(3) = a(1)*b(2) - a(2)*b(1)
  end function cross_product

  pure function det_3_by_3(mat)
    implicit none

    real(kind=DOUBLE), intent(in) :: mat(3, 3)
    real(kind=DOUBLE) :: det_3_by_3

    det_3_by_3 = mat(1, 1)*(mat(2, 2)*mat(3, 3) - mat(3, 2)*mat(2, 3)) &
      - mat(1, 2)*(mat(2, 1)*mat(3, 3) - mat(3, 1)*mat(2, 3)) &
      + mat(1, 3)*(mat(2, 1)*mat(3, 2) - mat(3, 1)*mat(2, 2))
  end function det_3_by_3

  pure subroutine qr_decomp(n, m, A, q, r, p)
    implicit none

    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), dimension(n, m), intent(in) :: A
    real(kind=DOUBLE), dimension(n, n), intent(inout) :: q
    real(kind=DOUBLE), dimension(n, m), intent(inout) :: r
    integer(kind=ENTIER), dimension(m), intent(inout) :: p

    integer(kind=ENTIER) :: i, j
    real(kind=DOUBLE), dimension(m) :: norms
    integer(kind=ENTIER), dimension(1) :: dummy
    real(kind=DOUBLE), dimension(n) :: y

    r = 0.0_DOUBLE
    q = 0.0_DOUBLE

    do i = 1, m
      norms(i) = norm2(A(:, i))
    end do
    do i = 1, m
      dummy = maxloc(norms)
      p(i) = dummy(1)
      norms(p(i)) = -1.0_DOUBLE
    end do

    r(1, 1) = norm2(A(:, p(1)))
    if (r(1, 1) < 1e-14) then
      r(:, 1) = 0.0_DOUBLE
      r(1, 1) = 1.0_DOUBLE
      q(:, 1) = 0.0_DOUBLE
    else
      q(:, 1) = A(:, p(1))/r(1, 1)
    end if

    do i = 2, m
      y = A(:, p(i))
      do j = 1, i - 1
        r(j, i) = dot_product(y, q(:, j))
        y = y - r(j, i)*q(:, j)
      end do

      r(i, i) = norm2(y)
      if (r(i, i) < 1e-14) then
        r(:, i) = 0.0_DOUBLE
        r(i, i) = 1.0_DOUBLE
        q(:, i) = 0.0_DOUBLE
      else
        q(:, i) = y/r(i, i)
      end if
    end do
  end subroutine qr_decomp

  pure subroutine qr_solve(n, m, A, x, b)
    implicit none

    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), dimension(n, m), intent(in) :: A
    real(kind=DOUBLE), dimension(n), intent(in) :: b
    real(kind=DOUBLE), dimension(m), intent(inout) :: x

    integer(kind=ENTIER) :: i, j
    integer(kind=ENTIER), dimension(m) :: p
    real(kind=DOUBLE), dimension(n) :: y
    real(kind=DOUBLE), dimension(n, n) :: q
    real(kind=DOUBLE), dimension(n, m) :: r

    call qr_decomp(n, m, A, q, r, p)
    y = matmul(transpose(q), b)

    do i = m, 1, -1
      x(p(i)) = y(i)
      if (r(i, i) > 1e-12_DOUBLE) then
        do j = i + 1, m
          x(p(i)) = x(p(i)) - r(i, j)*x(p(j))
        end do
        x(p(i)) = x(p(i))/r(i, i)
      else
        x(p(i)) = 0.0_DOUBLE
      end if
    end do
  end subroutine qr_solve

  pure subroutine qrdcmp(a, n, m, c, d, sing)
    implicit none
    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), intent(inout) :: a(n, m), c(m), d(m)
    logical, intent(inout) :: sing

    integer(kind=ENTIER) :: i, j, k
    real(kind=DOUBLE) :: sscale, sigma, ssum, tau
    sing = .FALSE.

    do k = 1, m
      sscale = 0.0_DOUBLE
      do i = k, n
        sscale = max(sscale, abs(a(i, k)))
      end do
      if (sscale .eq. 0.0_DOUBLE) then
        sing = .TRUE.
        c(k) = 0.0_DOUBLE
        d(k) = 0.0_DOUBLE
      else
        do i = k, n
          a(i, k) = a(i, k)/sscale
        end do
        ssum = 0.0_DOUBLE
        do i = k, n
          ssum = ssum + a(i, k)**2
        end do
        sigma = sign(sqrt(ssum), a(k, k))
        a(k, k) = a(k, k) + sigma
        c(k) = sigma*a(k, k)
        d(k) = -sscale*sigma
        do j = k + 1, m
          ssum = 0.0_DOUBLE
          do i = k, n
            ssum = ssum + a(i, k)*a(i, j)
          end do
          tau = ssum/c(k)
          do i = k, n
            a(i, j) = a(i, j) - tau*a(i, k)
          end do
        end do
      end if
    end do
    d(m) = a(m, m)
    if (d(m) .eq. 0.0_DOUBLE) sing = .TRUE.
  end subroutine qrdcmp

  subroutine qrsolv(a, n, m, c, d, b)
    implicit none
    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), intent(in) :: a(n, m), c(m), d(m)
    real(kind=DOUBLE), intent(inout) :: b(n)

    integer(kind=ENTIER) :: i, j
    real(kind=DOUBLE) :: ssum, tau

    do j = 1, m - 1
      ssum = 0.0_DOUBLE
      do i = 1, n
        ssum = ssum + a(i, j)*b(i)
      end do
      tau = ssum/c(j)
      do i = j, n
        b(i) = b(i) - tau*a(i, j)
      end do
    end do
    call rsolv(a, n, m, d, b)
  end subroutine qrsolv

  subroutine rsolv(a, n, m, d, b)
    implicit none
    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), intent(in) :: a(n, m), d(m)
    real(kind=DOUBLE), intent(inout) :: b(n)

    integer(kind=ENTIER) :: i, j
    real(kind=DOUBLE) :: ssum

    b(m) = b(m)/d(m)
    do i = m - 1, 1, -1
      ssum = 0.0_DOUBLE
      do j = i + 1, m
        ssum = ssum + a(i, j)*b(j)
      end do
      b(i) = (b(i) - ssum)/d(i)
    end do
  end subroutine rsolv

  subroutine print_mat(mat)
    implicit none

    integer(kind=ENTIER) :: i, j
    real(kind=DOUBLE), intent(in) :: mat(:, :)

    print *, ""
    do i = 1, size(mat(:, 1))
      do j = 1, size(mat(i, :))
        if (mat(i, j) > 1e-12_DOUBLE) then
          write (*, '(a,f18.8,a)', advance="no") ""//achar(27)//"[36m", mat(i, j), ""//achar(27)//"[0m"
        else if (mat(i, j) < -1e-12_DOUBLE) then
          write (*, '(a,f18.8,a)', advance="no") ""//achar(27)//"[34m", mat(i, j), ""//achar(27)//"[0m"
        else
          write (*, '(f18.8)', advance="no") mat(i, j)
        end if
      end do
      write (*, *) " "
    end do
    print *, ""
  end subroutine print_mat

  pure subroutine lu_solve_noswp(n, A, x, b, C)
    implicit none

    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), dimension(n, n), intent(in) :: A
    real(kind=DOUBLE), dimension(n, n), intent(inout) :: C
    real(kind=DOUBLE), dimension(n), intent(inout) :: x
    real(kind=DOUBLE), dimension(n), intent(in) :: b

    integer(kind=ENTIER) :: k, i, j
    real(kind=DOUBLE) :: invdiag, lik

    C(:, :) = A(:, :)
    x(:) = b(:)
    do k = 1, n
      invdiag = 1.0_DOUBLE/C(k, k)
      do i = k + 1, n
        lik = C(i, k)*invdiag
        do j = k, n
          C(i, j) = C(i, j) - lik*C(k, j)
        end do
        x(i) = x(i) - lik*x(k)
      end do
    end do

    x(n) = x(n)/C(n, n)
    do k = n - 1, 1, -1
      do j = k + 1, n
        x(k) = x(k) - C(k, j)*x(j)
      end do
      x(k) = x(k)/C(k, k)
    end do
  end subroutine lu_solve_noswp

  subroutine lu_solve(n, A, x, b)
    implicit none

    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), dimension(n, n), intent(in) :: A
    real(kind=DOUBLE), dimension(n), intent(inout) :: x
    real(kind=DOUBLE), dimension(n), intent(in) :: b

    integer(kind=ENTIER) :: k, pivot, i, j
    real(kind=DOUBLE) :: invdiag, lik
    real(kind=DOUBLE), dimension(n, n) :: C

    C = A
    x = b

    do k = 1, n
      pivot = maxloc(abs(C(k:n, k)), dim=1)
      call swap(C(k, :), C(pivot + k - 1, :))
      call swap(x(k), x(pivot + k - 1))
      invdiag = 1.0_DOUBLE/C(k, k)
      do i = k + 1, n
        lik = C(i, k)*invdiag
        do j = k, n
          C(i, j) = C(i, j) - lik*C(k, j)
        end do
        x(i) = x(i) - lik*x(k)
      end do
    end do

    if (abs(C(n, n)) < 1e-12_DOUBLE) error stop
    x(n) = x(n)/C(n, n)
    do k = n - 1, 1, -1
      do j = k + 1, n
        x(k) = x(k) - C(k, j)*x(j)
      end do
      x(k) = x(k)/C(k, k)
    end do
  end subroutine lu_solve

  pure subroutine lu_solve_inplace(n, A, b)
    implicit none

    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), dimension(n, n), intent(inout) :: A
    real(kind=DOUBLE), dimension(n), intent(inout) :: b

    integer(kind=ENTIER) :: k, pivot, i, j
    real(kind=DOUBLE) :: invdiag, lik

    do k = 1, n
      pivot = maxloc(abs(A(k:n, k)), dim=1)
      call swap(A(k, :), A(pivot + k - 1, :))
      call swap(b(k), b(pivot + k - 1))
      invdiag = 1.0_DOUBLE/A(k, k)
      do i = k + 1, n
        lik = A(i, k)*invdiag
        do j = k, n
          A(i, j) = A(i, j) - lik*A(k, j)
        end do
        b(i) = b(i) - lik*b(k)
      end do
    end do

    b(n) = b(n)/A(n, n)
    do k = n - 1, 1, -1
      do j = k + 1, n
        b(k) = b(k) - A(k, j)*b(j)
      end do
      b(k) = b(k)/A(k, k)
    end do
  end subroutine lu_solve_inplace

  elemental subroutine swap(a, b)
    real(kind=DOUBLE), intent(inout) :: a, b
    real(kind=DOUBLE) :: save_val
    save_val = a
    a = b
    b = save_val
  end subroutine swap

  pure function eye(n)
    implicit none
    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), dimension(n, n) :: eye

    integer(kind=ENTIER) :: i

    eye = 0.0_DOUBLE
    do i = 1, n
      eye(i, i) = 1.0_DOUBLE
    end do
  end function eye

  pure function eye_vector(i, n)
    implicit none
    integer(kind=ENTIER), intent(in) :: i, n
    real(kind=DOUBLE), dimension(n) :: eye_vector

    eye_vector = 0.0_DOUBLE
    eye_vector(i) = 1.0_DOUBLE
  end function eye_vector

  subroutine lu_inverse(n, A, Ainv)
    implicit none

    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), dimension(n, n), intent(in) :: A
    real(kind=DOUBLE), dimension(n, n), intent(inout) :: Ainv

    integer(kind=ENTIER) :: k, pivot, i, j, p
    real(kind=DOUBLE) :: invdiag, lik
    real(kind=DOUBLE), dimension(n, n) :: C

    C(:, :) = A(:, :)
    Ainv(:, :) = eye(n)

    do k = 1, n
      pivot = maxloc(abs(C(k:n, k)), dim=1)
      call swap(C(k, :), C(pivot + k - 1, :))
      call swap(Ainv(k, :), Ainv(pivot + k - 1, :))
      invdiag = 1.0_DOUBLE/C(k, k)
      do i = k + 1, n
        lik = C(i, k)*invdiag
        do j = k, n
          C(i, j) = C(i, j) - lik*C(k, j)
        end do
        do p = 1, n
          Ainv(i, p) = Ainv(i, p) - lik*Ainv(k, p)
        end do
      end do
    end do

    do p = 1, n
      Ainv(n, p) = Ainv(n, p)/C(n, n)
      do k = n - 1, 1, -1
        do j = k + 1, n
          Ainv(k, p) = Ainv(k, p) - C(k, j)*Ainv(j, p)
        end do
        Ainv(k, p) = Ainv(k, p)/C(k, k)
      end do
    end do
  end subroutine lu_inverse

  subroutine lu_inverse_modif_input(n, A, Ainv)
    implicit none

    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), dimension(n, n), intent(inout) :: A
    real(kind=DOUBLE), dimension(n, n), intent(inout) :: Ainv

    integer(kind=ENTIER) :: k, pivot, i, j, p
    real(kind=DOUBLE) :: invdiag, lik

    Ainv(:, :) = eye(n)

    do k = 1, n
      pivot = maxloc(abs(A(k:n, k)), dim=1)
      call swap(A(k, :), A(pivot + k - 1, :))
      call swap(Ainv(k, :), Ainv(pivot + k - 1, :))
      invdiag = 1.0_DOUBLE/A(k, k)
      do i = k + 1, n
        lik = A(i, k)*invdiag
        do j = k, n
          A(i, j) = A(i, j) - lik*A(k, j)
        end do
        do p = 1, n
          Ainv(i, p) = Ainv(i, p) - lik*Ainv(k, p)
        end do
      end do
    end do

    do p = 1, n
      Ainv(n, p) = Ainv(n, p)/A(n, n)
      do k = n - 1, 1, -1
        do j = k + 1, n
          Ainv(k, p) = Ainv(k, p) - A(k, j)*Ainv(j, p)
        end do
        Ainv(k, p) = Ainv(k, p)/A(k, k)
      end do
    end do
  end subroutine lu_inverse_modif_input

  subroutine inv_lapack(n, A)
    implicit none

    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), dimension(n, n), intent(inout) :: A

    integer(kind=ENTIER) :: info
    integer(kind=ENTIER), dimension(n) :: ipiv
    real(kind=DOUBLE), dimension(n) :: work

    call dgetrf(n, n, A, n, ipiv, info)
    call dgetri(n, A, n, ipiv, work, n, info)
  end subroutine inv_lapack

  subroutine lu_solve_inplace_lapack(n, A, x, b)
    implicit none
    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), dimension(n,n), intent(inout):: A
    real(kind=DOUBLE), dimension(n), intent(in)   :: b
    real(kind=DOUBLE), dimension(n), intent(out)  :: x
    integer(kind=ENTIER), dimension(n) :: ipiv
    integer(kind=ENTIER) :: info
    real(kind=DOUBLE), dimension(n) :: bcopy

    bcopy = b
    call dgetrf(n, n, A, n, ipiv, info)
    call dgetrs('N', n, 1, A, n, ipiv, bcopy, n, info)
    x = bcopy
  end subroutine lu_solve_inplace_lapack

  subroutine svd_solve_lapack(A, x, b, n)
    implicit none
    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), dimension(n,n), intent(inout) :: A
    real(kind=DOUBLE), dimension(n), intent(in) :: b
    real(kind=DOUBLE), dimension(n), intent(out) :: x

    real(kind=DOUBLE), dimension(n) :: bcopy
    real(kind=DOUBLE), dimension(n) :: S
    real(kind=DOUBLE), dimension(1) :: work_query
    real(kind=DOUBLE), allocatable, dimension(:) :: work
    integer(kind=ENTIER) :: info, lwork, rank
    real(kind=DOUBLE) :: rcond

    bcopy = b
    lwork = -1
    call dgelss(n, n, 1, A, n, bcopy, n, S, -1.0_DOUBLE, rank, work_query, lwork, info)
    lwork = int(work_query(1), kind=ENTIER)
    allocate(work(lwork))
    rcond = -1.0_DOUBLE
    call dgelss(n, n, 1, A, n, bcopy, n, S, rcond, rank, work, lwork, info)
    x = bcopy(1:n)
  end subroutine svd_solve_lapack

  subroutine pseudo_inverse_inplace_lapack(n, A)
    implicit none
    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), intent(inout) :: A(n,n)
    real(kind=DOUBLE), dimension(n, n) :: U, Sigma_plus
    real(kind=DOUBLE), dimension(n) :: Sigma
    real(kind=DOUBLE), dimension(:), allocatable :: WORK
    real(kind=DOUBLE) :: TOL
    integer(kind=ENTIER) :: INFO, i, LWORK

    LWORK = -1
    allocate(WORK(1))
    call dgesvd('A','O', n, n, A, n, Sigma, U, n, A, n, WORK, LWORK, INFO)

    LWORK = int(WORK(1))
    deallocate(WORK)
    allocate(WORK(LWORK))

    call dgesvd('A','O', n, n, A, n, Sigma, U, n, A, n, WORK, 5*n, INFO)

    TOL = n * 1e-12_DOUBLE * maxval(Sigma)
    Sigma_plus = 0.0d0
    do i = 1, n
      if(Sigma(i) > TOL) then
        Sigma_plus(i,i) = 1.0d0 / Sigma(i)
      else
        Sigma_plus(i,i) = 0.0d0
      end if
    end do

    A = matmul(transpose(A), matmul(Sigma_plus, transpose(U)))
    deallocate(WORK)
  end subroutine pseudo_inverse_inplace_lapack

  pure subroutine lu_inverse_modif_input_no_swp_inplace(n, A)
    implicit none

    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), dimension(n, n), intent(inout) :: A
    real(kind=DOUBLE), dimension(n, n) :: Ainv

    integer(kind=ENTIER) :: k, i, j, p
    real(kind=DOUBLE) :: invdiag, lik

    Ainv(:, :) = eye(n)

    do k = 1, n
      invdiag = 1.0_DOUBLE/A(k, k)
      do i = k + 1, n
        lik = A(i, k)*invdiag
        do j = k, n
          A(i, j) = A(i, j) - lik*A(k, j)
        end do
        do p = 1, n
          Ainv(i, p) = Ainv(i, p) - lik*Ainv(k, p)
        end do
      end do
    end do

    do p = 1, n
      Ainv(n, p) = Ainv(n, p)/A(n, n)
      do k = n - 1, 1, -1
        do j = k + 1, n
          Ainv(k, p) = Ainv(k, p) - A(k, j)*Ainv(j, p)
        end do
        Ainv(k, p) = Ainv(k, p)/A(k, k)
      end do
    end do

    A = Ainv
  end subroutine lu_inverse_modif_input_no_swp_inplace

  pure subroutine lu_inverse_modif_input_no_swp(n, A, Ainv)
    implicit none

    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), dimension(n, n), intent(inout) :: A
    real(kind=DOUBLE), dimension(n, n), intent(inout) :: Ainv

    integer(kind=ENTIER) :: k, i, j, p
    real(kind=DOUBLE) :: invdiag, lik

    Ainv(:, :) = eye(n)

    do k = 1, n
      invdiag = 1.0_DOUBLE/A(k, k)
      do i = k + 1, n
        lik = A(i, k)*invdiag
        do j = k, n
          A(i, j) = A(i, j) - lik*A(k, j)
        end do
        do p = 1, n
          Ainv(i, p) = Ainv(i, p) - lik*Ainv(k, p)
        end do
      end do
    end do

    do p = 1, n
      Ainv(n, p) = Ainv(n, p)/A(n, n)
      do k = n - 1, 1, -1
        do j = k + 1, n
          Ainv(k, p) = Ainv(k, p) - A(k, j)*Ainv(j, p)
        end do
        Ainv(k, p) = Ainv(k, p)/A(k, k)
      end do
    end do
  end subroutine lu_inverse_modif_input_no_swp

  pure subroutine qr_solve_with_q_r(n, m, q, r, p, x, b)
    implicit none

    integer(kind=ENTIER), intent(in) :: n, m
    real(kind=DOUBLE), dimension(n), intent(in) :: b
    real(kind=DOUBLE), dimension(m), intent(inout) :: x

    integer(kind=ENTIER) :: i, j
    integer(kind=ENTIER), dimension(m), intent(in) :: p
    real(kind=DOUBLE), dimension(n) :: y
    real(kind=DOUBLE), dimension(n, n), intent(in) :: q
    real(kind=DOUBLE), dimension(n, m), intent(in) :: r

    y = matmul(transpose(q), b)

    do i = m, 1, -1
      x(p(i)) = y(i)
      if (r(i, i) > 1e-12_DOUBLE) then
        do j = i + 1, m
          x(p(i)) = x(p(i)) - r(i, j)*x(p(j))
        end do
        x(p(i)) = x(p(i))/r(i, i)
      else
        x(p(i)) = 0.0_DOUBLE
      end if
    end do
  end subroutine qr_solve_with_q_r

  subroutine matmul_blas(A, B, C, m, n, k)
    implicit none
    integer(kind=ENTIER), intent(in) :: m, n, k
    real(kind=DOUBLE), dimension(m,k), intent(in) :: A
    real(kind=DOUBLE), dimension(k,n), intent(in) :: B
    real(kind=DOUBLE), dimension(m,n), intent(out) :: C

    real(kind=DOUBLE), parameter :: alpha = 1.0_DOUBLE
    real(kind=DOUBLE), parameter :: beta = 0.0_DOUBLE

    call dgemm('N','N', m, n, k, alpha, A, m, B, k, beta, C, m)
  end subroutine matmul_blas

  subroutine matmul3_blas(A, B, C, D, m, k, n, p)
    implicit none
    integer(kind=ENTIER), intent(in) :: m, k, n, p
    real(kind=DOUBLE), dimension(m,k), intent(in) :: A
    real(kind=DOUBLE), dimension(k,n), intent(in) :: B
    real(kind=DOUBLE), dimension(n,p), intent(in) :: C
    real(kind=DOUBLE), dimension(m,p), intent(out) :: D

    real(kind=DOUBLE), allocatable, dimension(:,:) :: T1, T2
    integer(kind=ENTIER) :: cost1, cost2

    cost1 = m*n*k + m*p*n
    cost2 = k*n*p + m*p*k
    if (cost1 <= cost2) then
      allocate(T1(m,n))
      call dgemm('N','N', m, n, k, 1.0_DOUBLE, A, m, B, k, 0.0_DOUBLE, T1, m)
      call dgemm('N','N', m, p, n, 1.0_DOUBLE, T1, m, C, n, 0.0_DOUBLE, D, m)
      deallocate(T1)
    else
      allocate(T2(k,p))
      call dgemm('N','N', k, p, n, 1.0_DOUBLE, B, k, C, n, 0.0_DOUBLE, T2, k)
      call dgemm('N','N', m, p, k, 1.0_DOUBLE, A, m, T2, k, 0.0_DOUBLE, D, m)
      deallocate(T2)
    end if
  end subroutine matmul3_blas

  subroutine matvec_blas(A, B, C, m, n)
    implicit none
    integer(kind=ENTIER), intent(in) :: m, n
    real(kind=DOUBLE), dimension(m,n), intent(in) :: A
    real(kind=DOUBLE), dimension(n), intent(in)   :: B
    real(kind=DOUBLE), dimension(m), intent(out)  :: C

    call dgemv('N', m, n, 1.0_DOUBLE, A, m, B, 1, 0.0_DOUBLE, C, 1)
  end subroutine matvec_blas

  subroutine matvec2_blas(A, B, C, D, m, n, p)
    implicit none
    integer(kind=ENTIER), intent(in) :: m, n, p
    real(kind=DOUBLE), dimension(m,n), intent(in) :: A
    real(kind=DOUBLE), dimension(n,p), intent(in) :: B
    real(kind=DOUBLE), dimension(p), intent(in)   :: C
    real(kind=DOUBLE), dimension(m), intent(out)  :: D

    real(kind=DOUBLE), allocatable, dimension(:) :: tmp

    allocate(tmp(n))
    call dgemv('N', n, p, 1.0_DOUBLE, B, n, C, 1, 0.0_DOUBLE, tmp, 1)
    call dgemv('N', m, n, 1.0_DOUBLE, A, m, tmp, 1, 0.0_DOUBLE, D, 1)
    deallocate(tmp)
  end subroutine matvec2_blas

  subroutine matmul2_blas(A, B, C, m, k, n)
    implicit none
    integer(kind=ENTIER), intent(in) :: m, k, n
    real(kind=DOUBLE), dimension(m, k), intent(in)  :: A
    real(kind=DOUBLE), dimension(k, n), intent(in)  :: B
    real(kind=DOUBLE), dimension(m, n), intent(out) :: C
    call dgemm('N', 'N', m, n, k, 1.0_DOUBLE, A, m, B, k, 0.0_DOUBLE, C, m)
  end subroutine matmul2_blas

  subroutine add_block(N, i0, j0, di, dj, M)
    implicit none
    integer(kind=ENTIER), intent(in) :: i0, j0, di, dj
    real(kind=DOUBLE), dimension(:,:), intent(inout) :: N
    real(kind=DOUBLE), dimension(:,:), intent(in)    :: M
    integer(kind=ENTIER) :: i, j

    if (i0 < 1 .or. j0 < 1 .or. i0+di-1 > size(N,1) .or. j0+dj-1 > size(N,2)) then
      error stop "add_block: indices ou taille du bloc invalides"
    end if
    if (size(M,1) /= di .or. size(M,2) /= dj) then
      error stop "add_block: dimensions de M incorrectes"
    end if

    do j = 1, dj
      do i = 1, di
        N(i0+i-1, j0+j-1) = N(i0+i-1, j0+j-1) + M(i,j)
      end do
    end do
  end subroutine add_block

  subroutine add_block_3(N, i0, j0, M)
    implicit none
    integer(kind=ENTIER), intent(in) :: i0, j0
    real(kind=DOUBLE), dimension(:, :), intent(inout) :: N
    real(kind=DOUBLE), dimension(3, 3), intent(in)    :: M

    N(i0:i0+2, j0  ) = N(i0:i0+2, j0  ) + M(:,1)
    N(i0:i0+2, j0+1) = N(i0:i0+2, j0+1) + M(:,2)
    N(i0:i0+2, j0+2) = N(i0:i0+2, j0+2) + M(:,3)
  end subroutine add_block_3

  subroutine add_block_3_transpose(N, i0, j0, M)
    implicit none
    integer(kind=ENTIER), intent(in) :: i0, j0
    real(kind=DOUBLE), dimension(:, :), intent(inout) :: N
    real(kind=DOUBLE), dimension(3, 3), intent(in)    :: M

    N(i0:i0+2, j0  ) = N(i0:i0+2, j0  ) + M(1,:)
    N(i0:i0+2, j0+1) = N(i0:i0+2, j0+1) + M(2,:)
    N(i0:i0+2, j0+2) = N(i0:i0+2, j0+2) + M(3,:)
  end subroutine add_block_3_transpose

  subroutine lu_factor_lapack(n, A, ipiv)
    implicit none
    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), dimension(n, n), intent(inout) :: A
    integer(kind=ENTIER), dimension(n), intent(out) :: ipiv
    integer(kind=ENTIER) :: info
    call dgetrf(n, n, A, n, ipiv, info)
  end subroutine lu_factor_lapack

  subroutine lu_solve_mat_lapack(n, A, ipiv, nrhs, B)
    implicit none
    integer(kind=ENTIER), intent(in) :: n, nrhs
    real(kind=DOUBLE), dimension(n, n), intent(in) :: A
    integer(kind=ENTIER), dimension(n), intent(in) :: ipiv
    real(kind=DOUBLE), dimension(n, nrhs), intent(inout) :: B
    integer(kind=ENTIER) :: info
    call dgetrs('N', n, nrhs, A, n, ipiv, B, n, info)
  end subroutine lu_solve_mat_lapack

  subroutine lu_solve_vec_lapack(n, A, ipiv, b)
    implicit none
    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), dimension(n, n), intent(in) :: A
    integer(kind=ENTIER), dimension(n), intent(in) :: ipiv
    real(kind=DOUBLE), dimension(n), intent(inout) :: b
    integer(kind=ENTIER) :: info
    call dgetrs('N', n, 1, A, n, ipiv, b, n, info)
  end subroutine lu_solve_vec_lapack

  function matrix_abs_Newton_Schulz(A) result(B)
    implicit none
    real(kind=DOUBLE), dimension(:, :), intent(in) :: A
    real(kind=DOUBLE), dimension(size(A, 1), size(A, 2)) :: B

    integer(kind=ENTIER), parameter :: max_iter = 100
    real(kind=DOUBLE), parameter :: tol = 1e-8_DOUBLE

    integer(kind=ENTIER) :: iter
    real(kind=DOUBLE) :: res
    real(kind=DOUBLE), dimension(size(A, 1), size(A, 2)) :: C

    res = 1e6_DOUBLE
    iter = 0
    B = A / (1.01_DOUBLE * norm2(A))
    do while (res > tol .and. iter < max_iter)
      C = 0.5_DOUBLE * matmul(B, 3.0_DOUBLE * eye(size(B, 1)) - matmul(B, B))
      res = norm2(B - C)
      B = C
      iter = iter + 1
    end do
    B = matmul(B, A)
  end function matrix_abs_Newton_Schulz

  pure function diag(a) result(M)
    implicit none
    real(kind=DOUBLE), dimension(:), intent(in) :: a
    real(kind=DOUBLE), dimension(size(a), size(a)) :: M

    integer(kind=ENTIER) :: i

    M = 0.0_DOUBLE
    do i = 1, size(a)
      M(i, i) = a(i)
    end do
  end function diag

  subroutine spectral_decomposition(M, Q, Lambda)
    implicit none

    real(kind=DOUBLE), intent(in) :: M(:, :)
    real(kind=DOUBLE), intent(out) :: Q(:, :), Lambda(:)

    integer(kind=ENTIER) :: n, lda, lwork, info
    real(kind=DOUBLE) :: work_query(1)
    real(kind=DOUBLE), allocatable :: work(:)

    n = size(M, 1)
    lda = n
    Q = M
    lwork = -1
    call dsyev('V', 'U', n, Q, lda, Lambda, work_query, lwork, info)
    lwork = int(work_query(1), ENTIER)
    allocate(work(lwork))
    call dsyev('V', 'U', n, Q, lda, Lambda, work, lwork, info)
    deallocate(work)
  end subroutine spectral_decomposition

end module linear_solver_module
