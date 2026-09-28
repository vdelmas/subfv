! Local copy of src/rbf/rbf_module.F90's hard/sliding RBF field solve,
! adapted so ale/ never depends on the rbf/ library (see ale_module.F90's
! header comment for why: ale/ only links core). Used to interpolate a
! smooth ALE grid velocity from a set of "hard" nodes (velocity imposed
! exactly, e.g. a tracked material interface or a piston) and "sliding"
! nodes (free to move tangentially along a prescribed surface normal, e.g.
! a domain wall) out to every other mesh vertex.
module ale_rbf_module
  use precision_module
  implicit none
contains
  pure function dof_offset(i, nimp) result(off)
    !flat-storage offset (0-based) before point i's block in a (3*nimp+nslide)
    !packed dof vector: hard points 1:nimp each own 3 slots, sliding points
    !nimp+1:n each own 1 (their scalar value along no)
    implicit none
    integer(kind=ENTIER), intent(in) :: i, nimp
    integer(kind=ENTIER) :: off

    if (i <= nimp) then
      off = 3*(i - 1)
    else
      off = 3*nimp + (i - nimp - 1)
    end if
  end function dof_offset

  subroutine compute_rbf_field(nimp, nslide, x, no, delta_imp, rad, delta)
    !nimp      = number of points where the value is imposed strongly (hard)
    !nslide    = number of points free to slide on a surface (weak)
    !x         = coordinates of all n=nimp+nslide points, 1:nimp hard then
    !            nimp+1:n sliding
    !no        = local (unit) normal at the sliding points, (3,n)-shaped for
    !            uniformity, meaningless/unused on the first nimp columns
    !delta_imp = imposed value, packed (3*nimp+nslide): 3 components per
    !            hard point (1:nimp), then 1 component per sliding point
    !            (value along no, typically 0 -- the point stays on its
    !            surface, free tangentially since that's not constrained)
    !rad       = radius of the compact support basis function
    !delta     = resulting coefficients to evaluate the field, same packing
    implicit none

    integer(kind=ENTIER), intent(in) :: nimp, nslide
    real(kind=DOUBLE), dimension(3, nimp + nslide), intent(in) :: x, no
    real(kind=DOUBLE), dimension(3*nimp + nslide), intent(in) :: delta_imp
    real(kind=DOUBLE), intent(in) :: rad
    real(kind=DOUBLE), dimension(3*nimp + nslide), intent(inout) :: delta

    ! CG (not BiCGSTAB): the Wendland-C2 kernel used by rbf_func below is
    ! strictly positive definite in R^3 (unlike the earlier compact "tent"
    ! kernel max(0,1-x), which is only guaranteed PD in 1D), so the packed
    ! hard/sliding system assembled by matmul_rbf is SPD. Plain CG has a
    ! single matvec/iteration, a monotonically decreasing residual, and none
    ! of BiCGSTAB's rho/omega breakdown modes, so it reaches machine-precision
    ! residuals reliably instead of stalling.
    real(kind=DOUBLE), parameter :: tol = 1e-16_DOUBLE
    integer(kind=ENTIER), parameter :: maxiter = 20000
    integer(kind=ENTIER) :: n, iter
    real(kind=DOUBLE), dimension(3*nimp + nslide) :: r, p, ap
    real(kind=DOUBLE) :: rs_old, rs_new, alpha, beta, r0norm, target_res

    n = nimp + nslide

    r = delta_imp - matmul_rbf(n, nimp, x, no, delta, rad)
    p = r
    rs_old = sum(r*r)
    r0norm = sqrt(rs_old)
    target_res = tol*max(1.0_DOUBLE, r0norm)

    iter = 0
    do while (sqrt(rs_old) > target_res .and. iter < maxiter)
      ap = matmul_rbf(n, nimp, x, no, p, rad)
      alpha = rs_old/(1e-300_DOUBLE + sum(p*ap))
      delta = delta + alpha*p
      r = r - alpha*ap
      rs_new = sum(r*r)
      if (sqrt(rs_new) <= target_res) then
        rs_old = rs_new
        exit
      end if
      beta = rs_new/rs_old
      p = r + beta*p
      rs_old = rs_new
      iter = iter + 1
    end do
    if (iter >= maxiter) print*, "ale_rbf: Failed to converge in maxiter", maxiter, &
      "residual =", sqrt(rs_old)
  end subroutine compute_rbf_field

  function matmul_rbf(n, nimp, x, no, a, rad) result(b)
    !n     = nimp+nslide, total number of rbf points
    !nimp  = points 1:nimp are hard (3 packed dof), nimp+1:n are sliding (1
    !        packed dof, along the local normal no)
    !a, b  = packed (3*nimp+nslide) dof vectors, see dof_offset
    implicit none

    integer(kind=ENTIER), intent(in) :: n, nimp
    real(kind=DOUBLE), dimension(3, n), intent(in) :: x, no
    real(kind=DOUBLE), dimension(3*nimp + (n - nimp)), intent(in) :: a
    real(kind=DOUBLE), dimension(3*nimp + (n - nimp)) :: b
    real(kind=DOUBLE) :: rad

    integer(kind=ENTIER) :: i, j, oi, oj
    real(kind=DOUBLE) :: kij
    real(kind=DOUBLE), dimension(3) :: ai3, aj3

    b = 0.0_DOUBLE
    do i = 1, n
      oi = dof_offset(i, nimp)
      if (i <= nimp) then
        b(oi + 1:oi + 3) = b(oi + 1:oi + 3) + a(oi + 1:oi + 3)
      else
        b(oi + 1) = b(oi + 1) + a(oi + 1)
      end if

      do j = i + 1, n
        oj = dof_offset(j, nimp)
        kij = rbf_func(norm2(x(:, j) - x(:, i))/rad)

        if (j <= nimp) then
          aj3 = a(oj + 1:oj + 3)
        else
          aj3 = a(oj + 1)*no(:, j)
        end if
        if (i <= nimp) then
          b(oi + 1:oi + 3) = b(oi + 1:oi + 3) + kij*aj3
        else
          b(oi + 1) = b(oi + 1) + kij*dot_product(no(:, i), aj3)
        end if

        if (i <= nimp) then
          ai3 = a(oi + 1:oi + 3)
        else
          ai3 = a(oi + 1)*no(:, i)
        end if
        if (j <= nimp) then
          b(oj + 1:oj + 3) = b(oj + 1:oj + 3) + kij*ai3
        else
          b(oj + 1) = b(oj + 1) + kij*dot_product(no(:, j), ai3)
        end if
      end do
    end do
  end function matmul_rbf

  function rbf_func(x) result(y)
    !Wendland C2 compactly-supported RBF, strictly positive definite in R^3
    !(unlike the linear "tent" kernel max(0,1-x), which is PD only in 1D and
    !left the assembled system only weakly/indefinitely conditioned).
    implicit none

    real(kind=DOUBLE), intent(in) :: x
    real(kind=DOUBLE) :: y, xm

    if (x >= 1.0_DOUBLE) then
      y = 0.0_DOUBLE
    else
      xm = 1.0_DOUBLE - x
      y = (xm**4)*(4.0_DOUBLE*x + 1.0_DOUBLE)
    end if
  end function rbf_func

  subroutine eval_rbf_field(n, nimp, x, no, rad, delta, xi, di)
    !n, nimp, x, no, delta : see matmul_rbf/compute_rbf_field
    !xi = point at which the field is evaluated
    implicit none

    integer(kind=ENTIER), intent(in) :: n, nimp
    real(kind=DOUBLE), dimension(3, n), intent(in) :: x, no
    real(kind=DOUBLE), intent(in) :: rad
    real(kind=DOUBLE), dimension(3*nimp + (n - nimp)), intent(in) :: delta
    real(kind=DOUBLE), dimension(3), intent(in) :: xi
    real(kind=DOUBLE), dimension(3), intent(inout) :: di

    integer(kind=ENTIER) :: i, oi
    real(kind=DOUBLE) :: kij

    di = 0.0_DOUBLE
    do i = 1, n
      oi = dof_offset(i, nimp)
      kij = rbf_func(norm2(xi - x(:, i))/rad)
      if (i <= nimp) then
        di = di + kij*delta(oi + 1:oi + 3)
      else
        di = di + kij*delta(oi + 1)*no(:, i)
      end if
    end do
  end subroutine eval_rbf_field
end module ale_rbf_module
