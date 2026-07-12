module sort_module
  use precision_module
  implicit none

  private
  public :: qsort_key, qsort_int, qsort_real, qsort_mat
  public :: add_sort_unique_int

contains

  subroutine qsort_key(n, ra)
    implicit none

    integer(kind=ENTIER), intent(in) :: n
    integer(kind=ENTIER), dimension(n), intent(inout) :: ra

    if (n < 2) return
    call qsort_key_r(ra, 1, n)
  end subroutine qsort_key

  recursive subroutine qsort_key_r(ra, lo, hi)
    implicit none

    integer(kind=ENTIER), dimension(:), intent(inout) :: ra
    integer(kind=ENTIER), intent(in) :: lo, hi

    integer(kind=ENTIER) :: i, j, pivot, tmp

    if (lo >= hi) return
    pivot = ra((lo + hi) / 2)
    i = lo
    j = hi
    do
      do while (ra(i) < pivot); i = i + 1; end do
      do while (ra(j) > pivot); j = j - 1; end do
      if (i >= j) exit
      tmp = ra(i); ra(i) = ra(j); ra(j) = tmp
      i = i + 1; j = j - 1
    end do
    call qsort_key_r(ra, lo, j)
    call qsort_key_r(ra, j + 1, hi)
  end subroutine qsort_key_r

  subroutine qsort_int(n, ra, rb)
    implicit none

    integer(kind=ENTIER), intent(in) :: n
    integer(kind=ENTIER), dimension(n), intent(inout) :: ra
    integer(kind=ENTIER), dimension(n), intent(inout) :: rb

    if (n < 2) return
    call qsort_int_r(ra, rb, 1, n)
  end subroutine qsort_int

  recursive subroutine qsort_int_r(ra, rb, lo, hi)
    implicit none

    integer(kind=ENTIER), dimension(:), intent(inout) :: ra
    integer(kind=ENTIER), dimension(:), intent(inout) :: rb
    integer(kind=ENTIER), intent(in) :: lo, hi

    integer(kind=ENTIER) :: i, j, pivot, tmp

    if (lo >= hi) return
    pivot = ra((lo + hi) / 2)
    i = lo
    j = hi
    do
      do while (ra(i) < pivot); i = i + 1; end do
      do while (ra(j) > pivot); j = j - 1; end do
      if (i >= j) exit
      tmp = ra(i); ra(i) = ra(j); ra(j) = tmp
      tmp = rb(i); rb(i) = rb(j); rb(j) = tmp
      i = i + 1; j = j - 1
    end do
    call qsort_int_r(ra, rb, lo, j)
    call qsort_int_r(ra, rb, j + 1, hi)
  end subroutine qsort_int_r

  subroutine qsort_real(n, ra, rb)
    implicit none

    integer(kind=ENTIER), intent(in) :: n
    integer(kind=ENTIER), dimension(n), intent(inout) :: ra
    real(kind=DOUBLE), dimension(n), intent(inout) :: rb

    if (n < 2) return
    call qsort_real_r(ra, rb, 1, n)
  end subroutine qsort_real

  recursive subroutine qsort_real_r(ra, rb, lo, hi)
    implicit none

    integer(kind=ENTIER), dimension(:), intent(inout) :: ra
    real(kind=DOUBLE), dimension(:), intent(inout) :: rb
    integer(kind=ENTIER), intent(in) :: lo, hi

    integer(kind=ENTIER) :: i, j, pivot, tmp
    real(kind=DOUBLE) :: tmp_rb

    if (lo >= hi) return
    pivot = ra((lo + hi) / 2)
    i = lo
    j = hi
    do
      do while (ra(i) < pivot); i = i + 1; end do
      do while (ra(j) > pivot); j = j - 1; end do
      if (i >= j) exit
      tmp = ra(i); ra(i) = ra(j); ra(j) = tmp
      tmp_rb = rb(i); rb(i) = rb(j); rb(j) = tmp_rb
      i = i + 1; j = j - 1
    end do
    call qsort_real_r(ra, rb, lo, j)
    call qsort_real_r(ra, rb, j + 1, hi)
  end subroutine qsort_real_r

  subroutine qsort_mat(n, ra, rb)
    implicit none

    integer(kind=ENTIER), intent(in) :: n
    integer(kind=ENTIER), dimension(n), intent(inout) :: ra
    real(kind=DOUBLE), dimension(3, 3, n), intent(inout) :: rb

    if (n < 2) return
    call qsort_mat_r(ra, rb, 1, n)
  end subroutine qsort_mat

  recursive subroutine qsort_mat_r(ra, rb, lo, hi)
    implicit none

    integer(kind=ENTIER), dimension(:), intent(inout) :: ra
    real(kind=DOUBLE), dimension(:, :, :), intent(inout) :: rb
    integer(kind=ENTIER), intent(in) :: lo, hi

    integer(kind=ENTIER) :: i, j, pivot, tmp
    real(kind=DOUBLE), dimension(3, 3) :: tmp_rb

    if (lo >= hi) return
    pivot = ra((lo + hi) / 2)
    i = lo
    j = hi
    do
      do while (ra(i) < pivot); i = i + 1; end do
      do while (ra(j) > pivot); j = j - 1; end do
      if (i >= j) exit
      tmp = ra(i); ra(i) = ra(j); ra(j) = tmp
      tmp_rb = rb(:, :, i); rb(:, :, i) = rb(:, :, j); rb(:, :, j) = tmp_rb
      i = i + 1; j = j - 1
    end do
    call qsort_mat_r(ra, rb, lo, j)
    call qsort_mat_r(ra, rb, j + 1, hi)
  end subroutine qsort_mat_r

  subroutine add_sort_unique_int(a, na, b)
    implicit none

    integer(kind=ENTIER), dimension(:), allocatable, intent(inout) :: a
    integer(kind=ENTIER), intent(inout) :: na
    integer(kind=ENTIER), dimension(:), intent(in) :: b

    integer(kind=ENTIER) :: nb, i, iloc
    integer(kind=ENTIER), dimension(:), allocatable :: tmp

    nb = size(b)

    if (allocated(a)) then
      na = size(a)
      allocate(tmp(na))
      tmp = a
      deallocate(a)
      allocate(a(na + nb))
      a(:na) = tmp(:)
      a(na+1:) = b
      na = na + nb
      deallocate(tmp)
    else
      allocate(a(nb))
      a = b
      na = nb
    end if

    allocate(tmp(na))
    tmp = a

    call qsort_key(size(tmp), tmp)

    iloc = 1
    do i = 2, size(tmp)
      if (tmp(i) /= tmp(i-1)) iloc = iloc + 1
    end do

    deallocate(a)
    na = iloc
    allocate(a(na))
    iloc = 1
    a(1) = tmp(1)
    do i = 2, size(tmp)
      if (tmp(i) /= tmp(i-1)) then
        iloc = iloc + 1
        a(iloc) = tmp(i)
      end if
    end do

    deallocate(tmp)
  end subroutine add_sort_unique_int

end module sort_module
