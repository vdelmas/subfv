module rbf_module
  use precision_module
  implicit none
contains
  pure function dof_offset(i, nimp) result(off)
    !flat-storage offset (0-based) before point i's block in a (3*nimp+nslide)
    !packed dof vector: hard points 1:nimp each own 3 slots, sliding points
    !nimp+1:n each own 1 (their scalar displacement along no)
    implicit none
    integer(kind=4), intent(in) :: i, nimp
    integer(kind=4) :: off

    if( i <= nimp ) then
      off = 3*(i-1)
    else
      off = 3*nimp + (i-nimp-1)
    end if
  end function dof_offset

  subroutine compute_rbf_field(nimp, nslide, x, no, delta_imp, rad, delta)
    !nimp      = number of points where the displacement is imposed strongly (hard)
    !nslide    = number of points that are free to slide on a surface (weak)
    !x         = coordinates of all n=nimp+nslide points, 1:nimp hard then nimp+1:n sliding
    !no        = local (unit) normal at the sliding points, (3,n)-shaped for
    !            uniformity, meaningless/unused on the first nimp columns
    !delta_imp = imposed displacement, packed (3*nimp+nslide): 3 components
    !            per hard point (1:nimp), then 1 component per sliding point
    !            (displacement along no, typically 0 -- the point stays on
    !            its surface, free to move however it wants tangentially
    !            since that's not constrained/solved for at all)
    !rad       = radius of the compact support basis function
    !delta     = resulting coefficients to evaluate the field, same packing
    !            as delta_imp -- a sliding point only ever has 1 real dof
    !            (its position along no), matching the 1 equation actually
    !            enforced there, so the system is full rank (unlike an
    !            otimes(no) projector on a full 3-vector, which would leave 2
    !            unconstrained gauge directions per slip point)

    implicit none

    integer(kind=4), intent(in) :: nimp, nslide
    real(kind=8), dimension(3, nimp+nslide), intent(in) :: x, no
    real(kind=8), dimension(3*nimp+nslide), intent(in) :: delta_imp
    real(kind=8), intent(in) :: rad
    real(kind=8), dimension(3*nimp+nslide), intent(inout) :: delta

    real(kind=8), parameter :: tol=1e-8
    integer(kind=4), parameter :: maxiter=5000
    integer(kind=4) :: n, iter
    real(kind=8), dimension(3*nimp+nslide) :: r, r0hat, p, v, h, res, tt
    real(kind=8) :: rho, rho_old, alpha, omega, beta

    n = nimp + nslide

    r = delta_imp - matmul_rbf(n, nimp, x, no, delta, rad)
    r0hat = r
    rho_old = 1.d0; alpha = 1.d0; omega = 1.d0
    p = 0.d0; v = 0.d0

    iter = 0
    do while ((sqrt(sum(r*r)) > tol) .and. iter < maxiter)
      rho = sum(r0hat*r)
      beta = (rho/(1e-14+rho_old))*(alpha/(1e-14+omega))
      p = r + beta*(p - omega*v)
      v = matmul_rbf(n, nimp, x, no, p, rad)
      alpha = rho/(1e-14+sum(r0hat*v))
      h = delta + alpha*p
      res = r - alpha*v
      if( sqrt(sum(res*res)) < tol ) then
        delta = h
        r = res
        exit
      end if
      tt = matmul_rbf(n, nimp, x, no, res, rad)
      omega = sum(tt*res)/(1e-14+sum(tt*tt))
      delta = h + omega*res
      r = res - omega*tt
      rho_old = rho
      iter = iter + 1
    end do
    if( iter >= maxiter ) print*,"Failed to converge in maxiter", maxiter
  end subroutine compute_rbf_field

  function matmul_rbf(n, nimp, x, no, a, rad) result(b)
    !n     = nimp+nslide, total number of rbf points
    !nimp  = points 1:nimp are hard (3 packed dof), nimp+1:n are sliding (1
    !        packed dof, along the local normal no)
    !a, b  = packed (3*nimp+nslide) dof vectors, see dof_offset
    !
    !delta(x) = sum_q delta_q*phi(x-x_q) + sum_p no_p*delta_p*phi(x-x_p)
    !where q runs over the hard points (delta_q full 3-vector), p over the
    !sliding points (delta_p a scalar along no_p) -- this function applies
    !that same block structure to build the matrix-vector product, per pair
    !(i,j):
    !  hard-hard   : plain 3x3 identity/kernel
    !  hard-slide  : the sliding column's scalar dof is expanded to 3d via
    !                no_j before being added to the hard row
    !  slide-hard  : the hard column's 3d contribution is projected onto
    !                the sliding row's no_i to get its scalar dof
    !  slide-slide : both sides go through their own no_i/no_j
    implicit none

    integer(kind=4), intent(in) :: n, nimp
    real(kind=8), dimension(3, n), intent(in) :: x, no
    real(kind=8), dimension(3*nimp+(n-nimp)), intent(in) :: a
    real(kind=8), dimension(3*nimp+(n-nimp)) :: b
    real(kind=8) :: rad

    integer(kind=4) :: i, j, oi, oj
    real(kind=8) :: kij
    real(kind=8), dimension(3) :: ai3, aj3

    b = 0.d0
    do i=1, n
      oi = dof_offset(i, nimp)
      if( i <= nimp ) then
        b(oi+1:oi+3) = b(oi+1:oi+3) + a(oi+1:oi+3)
      else
        b(oi+1) = b(oi+1) + a(oi+1)
      end if

      do j=i+1, n
        oj = dof_offset(j, nimp)
        kij = rbf_func(norm2(x(:, j)-x(:, i))/rad)

        if( j <= nimp ) then
          aj3 = a(oj+1:oj+3)
        else
          aj3 = a(oj+1)*no(:, j)
        end if
        if( i <= nimp ) then
          b(oi+1:oi+3) = b(oi+1:oi+3) + kij*aj3
        else
          b(oi+1) = b(oi+1) + kij*dot_product(no(:, i), aj3)
        end if

        if( i <= nimp ) then
          ai3 = a(oi+1:oi+3)
        else
          ai3 = a(oi+1)*no(:, i)
        end if
        if( j <= nimp ) then
          b(oj+1:oj+3) = b(oj+1:oj+3) + kij*ai3
        else
          b(oj+1) = b(oj+1) + kij*dot_product(no(:, j), ai3)
        end if
      end do
    end do
  end function matmul_rbf

  function rbf_func(x) result(y)
    implicit none

    real(kind=8), intent(in) :: x
    real(kind=8) :: y

    y = max(0.d0, 1.d0-x)
  end function rbf_func

  subroutine eval_rbf_field(n, nimp, x, no, rad, delta, xi, di)
    !n, nimp, x, no, delta : see matmul_rbf/compute_rbf_field
    !xi = point at which the field is evaluated
    !di = delta(xi) = sum_q delta_q*phi(xi-x_q) + sum_p no_p*delta_p*phi(xi-x_p)
    implicit none

    integer(kind=4), intent(in) :: n, nimp
    real(kind=8), dimension(3, n), intent(in) :: x, no
    real(kind=8), intent(in) :: rad
    real(kind=8), dimension(3*nimp+(n-nimp)), intent(in) :: delta
    real(kind=8), dimension(3), intent(in) :: xi
    real(kind=8), dimension(3), intent(inout) :: di

    integer(kind=4) :: i, oi
    real(kind=8) :: kij

    di = 0.d0
    do i=1, n
      oi = dof_offset(i, nimp)
      kij = rbf_func(norm2(xi-x(:, i))/rad)
      if( i <= nimp ) then
        di = di + kij*delta(oi+1:oi+3)
      else
        di = di + kij*delta(oi+1)*no(:, i)
      end if
    end do
  end subroutine eval_rbf_field

  function imposed_vel(x, t)
    implicit none

    real(kind=8), intent(in) :: t
    real(kind=8), dimension(3), intent(in) :: x

    real(kind=8), dimension(3) :: imposed_vel

    imposed_vel = (/-1.d0/300.D0, 0.d0, 0.d0/)
  end function imposed_vel
end module rbf_module
