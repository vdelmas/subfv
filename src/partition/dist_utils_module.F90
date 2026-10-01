! Distributed helpers: block distribution of a global index range, Alltoallv
! exchanges and stable O(n log n) merge sorts on permutations.
module dist_utils_module
  use mpi
  use precision_module
  implicit none

  private
  public :: block_starts, block_owner
  public :: exchange_int, exchange_dbl
  public :: sort_perm_int, sort_perm_records, same_record

contains

  ! Global 1-based index range 1..n split into P contiguous blocks:
  ! rank r owns starts(r)+1 .. starts(r+1).
  subroutine block_starts(n, num_procs, starts)
    integer(kind=ENTIER_D), intent(in) :: n
    integer(kind=ENTIER), intent(in) :: num_procs
    integer(kind=ENTIER_D), dimension(0:num_procs), intent(out) :: starts

    integer(kind=ENTIER) :: r

    do r = 0, num_procs
      starts(r) = (n*int(r, ENTIER_D))/int(num_procs, ENTIER_D)
    end do
  end subroutine block_starts

  pure function block_owner(g, num_procs, starts) result(owner)
    integer(kind=ENTIER_D), intent(in) :: g
    integer(kind=ENTIER), intent(in) :: num_procs
    integer(kind=ENTIER_D), dimension(0:num_procs), intent(in) :: starts
    integer(kind=ENTIER) :: owner

    integer(kind=ENTIER) :: lo, hi, mid

    ! Find r with starts(r) < g <= starts(r+1); empty blocks are skipped
    ! naturally since starts is non-decreasing.
    lo = 0
    hi = num_procs - 1
    do while (lo < hi)
      mid = (lo + hi + 1)/2
      if (starts(mid) < g) then
        lo = mid
      else
        hi = mid - 1
      end if
    end do
    owner = lo
  end function block_owner

  ! sendbuf is ordered by destination rank (sendcounts(0:P-1) entries each);
  ! recvbuf comes back ordered by source rank.
  subroutine exchange_int(num_procs, sendcounts, sendbuf, recvcounts, recvbuf)
    integer(kind=ENTIER), intent(in) :: num_procs
    integer(kind=ENTIER), dimension(0:num_procs - 1), intent(in) :: sendcounts
    integer(kind=ENTIER), dimension(:), intent(in) :: sendbuf
    integer(kind=ENTIER), dimension(0:num_procs - 1), intent(out) :: recvcounts
    integer(kind=ENTIER), dimension(:), allocatable, intent(out) :: recvbuf

    integer(kind=ENTIER), dimension(0:num_procs - 1) :: sdispl, rdispl
    integer(kind=ENTIER) :: r, mpi_ierr
    integer(kind=ENTIER), dimension(1) :: dummy

    call MPI_ALLTOALL(sendcounts, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, mpi_ierr)

    sdispl(0) = 0
    rdispl(0) = 0
    do r = 1, num_procs - 1
      sdispl(r) = sdispl(r - 1) + sendcounts(r - 1)
      rdispl(r) = rdispl(r - 1) + recvcounts(r - 1)
    end do

    allocate (recvbuf(sum(int(recvcounts, ENTIER_D))))
    if (size(sendbuf) > 0 .and. size(recvbuf) > 0) then
      call MPI_ALLTOALLV(sendbuf, sendcounts, sdispl, MPI_INTEGER, &
        recvbuf, recvcounts, rdispl, MPI_INTEGER, MPI_COMM_WORLD, mpi_ierr)
    else if (size(sendbuf) > 0) then
      call MPI_ALLTOALLV(sendbuf, sendcounts, sdispl, MPI_INTEGER, &
        dummy, recvcounts, rdispl, MPI_INTEGER, MPI_COMM_WORLD, mpi_ierr)
    else if (size(recvbuf) > 0) then
      call MPI_ALLTOALLV(dummy, sendcounts, sdispl, MPI_INTEGER, &
        recvbuf, recvcounts, rdispl, MPI_INTEGER, MPI_COMM_WORLD, mpi_ierr)
    else
      call MPI_ALLTOALLV(dummy, sendcounts, sdispl, MPI_INTEGER, &
        dummy, recvcounts, rdispl, MPI_INTEGER, MPI_COMM_WORLD, mpi_ierr)
    end if
  end subroutine exchange_int

  subroutine exchange_dbl(num_procs, sendcounts, sendbuf, recvcounts, recvbuf)
    integer(kind=ENTIER), intent(in) :: num_procs
    integer(kind=ENTIER), dimension(0:num_procs - 1), intent(in) :: sendcounts
    real(kind=DOUBLE), dimension(:), intent(in) :: sendbuf
    integer(kind=ENTIER), dimension(0:num_procs - 1), intent(out) :: recvcounts
    real(kind=DOUBLE), dimension(:), allocatable, intent(out) :: recvbuf

    integer(kind=ENTIER), dimension(0:num_procs - 1) :: sdispl, rdispl
    integer(kind=ENTIER) :: r, mpi_ierr
    real(kind=DOUBLE), dimension(1) :: dummy

    call MPI_ALLTOALL(sendcounts, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, &
      MPI_COMM_WORLD, mpi_ierr)

    sdispl(0) = 0
    rdispl(0) = 0
    do r = 1, num_procs - 1
      sdispl(r) = sdispl(r - 1) + sendcounts(r - 1)
      rdispl(r) = rdispl(r - 1) + recvcounts(r - 1)
    end do

    allocate (recvbuf(sum(int(recvcounts, ENTIER_D))))
    if (size(sendbuf) > 0 .and. size(recvbuf) > 0) then
      call MPI_ALLTOALLV(sendbuf, sendcounts, sdispl, MPI_DOUBLE_PRECISION, &
        recvbuf, recvcounts, rdispl, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpi_ierr)
    else if (size(sendbuf) > 0) then
      call MPI_ALLTOALLV(sendbuf, sendcounts, sdispl, MPI_DOUBLE_PRECISION, &
        dummy, recvcounts, rdispl, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpi_ierr)
    else if (size(recvbuf) > 0) then
      call MPI_ALLTOALLV(dummy, sendcounts, sdispl, MPI_DOUBLE_PRECISION, &
        recvbuf, recvcounts, rdispl, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpi_ierr)
    else
      call MPI_ALLTOALLV(dummy, sendcounts, sdispl, MPI_DOUBLE_PRECISION, &
        dummy, recvcounts, rdispl, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpi_ierr)
    end if
  end subroutine exchange_dbl

  ! Stable merge sort: perm(1:n) such that key(perm) is non-decreasing.
  subroutine sort_perm_int(n, key, perm)
    integer(kind=ENTIER), intent(in) :: n
    integer(kind=ENTIER), dimension(n), intent(in) :: key
    integer(kind=ENTIER), dimension(n), intent(out) :: perm

    integer(kind=ENTIER), dimension(:), allocatable :: tmp
    integer(kind=ENTIER) :: i, width, lo, mid, hi, a, b, k

    do i = 1, n
      perm(i) = i
    end do
    if (n < 2) return
    allocate (tmp(n))

    width = 1
    do while (width < n)
      lo = 1
      do while (lo <= n)
        mid = min(lo + width - 1, n)
        hi = min(lo + 2*width - 1, n)
        a = lo
        b = mid + 1
        k = lo
        do while (a <= mid .and. b <= hi)
          if (key(perm(b)) < key(perm(a))) then
            tmp(k) = perm(b)
            b = b + 1
          else
            tmp(k) = perm(a)
            a = a + 1
          end if
          k = k + 1
        end do
        do while (a <= mid)
          tmp(k) = perm(a)
          a = a + 1
          k = k + 1
        end do
        do while (b <= hi)
          tmp(k) = perm(b)
          b = b + 1
          k = k + 1
        end do
        lo = lo + 2*width
      end do
      perm(1:n) = tmp(1:n)
      width = 2*width
    end do
  end subroutine sort_perm_int

  ! Stable merge sort of variable-length integer records
  ! dat(ptr(i):ptr(i+1)-1), ordered by length first, then lexicographically.
  subroutine sort_perm_records(n, ptr, dat, perm)
    integer(kind=ENTIER), intent(in) :: n
    integer(kind=ENTIER), dimension(n + 1), intent(in) :: ptr
    integer(kind=ENTIER), dimension(:), intent(in) :: dat
    integer(kind=ENTIER), dimension(n), intent(out) :: perm

    integer(kind=ENTIER), dimension(:), allocatable :: tmp
    integer(kind=ENTIER) :: i, width, lo, mid, hi, a, b, k

    do i = 1, n
      perm(i) = i
    end do
    if (n < 2) return
    allocate (tmp(n))

    width = 1
    do while (width < n)
      lo = 1
      do while (lo <= n)
        mid = min(lo + width - 1, n)
        hi = min(lo + 2*width - 1, n)
        a = lo
        b = mid + 1
        k = lo
        do while (a <= mid .and. b <= hi)
          if (record_less(perm(b), perm(a))) then
            tmp(k) = perm(b)
            b = b + 1
          else
            tmp(k) = perm(a)
            a = a + 1
          end if
          k = k + 1
        end do
        do while (a <= mid)
          tmp(k) = perm(a)
          a = a + 1
          k = k + 1
        end do
        do while (b <= hi)
          tmp(k) = perm(b)
          b = b + 1
          k = k + 1
        end do
        lo = lo + 2*width
      end do
      perm(1:n) = tmp(1:n)
      width = 2*width
    end do

  contains

    logical function record_less(i1, i2)
      integer(kind=ENTIER), intent(in) :: i1, i2
      integer(kind=ENTIER) :: l1, l2, j

      l1 = ptr(i1 + 1) - ptr(i1)
      l2 = ptr(i2 + 1) - ptr(i2)
      if (l1 /= l2) then
        record_less = l1 < l2
        return
      end if
      do j = 0, l1 - 1
        if (dat(ptr(i1) + j) /= dat(ptr(i2) + j)) then
          record_less = dat(ptr(i1) + j) < dat(ptr(i2) + j)
          return
        end if
      end do
      record_less = .false.
    end function record_less
  end subroutine sort_perm_records

  pure logical function same_record(ptr, dat, i1, i2)
    integer(kind=ENTIER), dimension(:), intent(in) :: ptr, dat
    integer(kind=ENTIER), intent(in) :: i1, i2

    integer(kind=ENTIER) :: l1, j

    l1 = ptr(i1 + 1) - ptr(i1)
    same_record = .false.
    if (l1 /= ptr(i2 + 1) - ptr(i2)) return
    do j = 0, l1 - 1
      if (dat(ptr(i1) + j) /= dat(ptr(i2) + j)) return
    end do
    same_record = .true.
  end function same_record
end module dist_utils_module
