module mpi_module
  use precision_module
  implicit none

  type mpi_elems_type
    integer(kind=ENTIER) :: partition_id = -1, n_elems = 0
    integer(kind=ENTIER), dimension(:), allocatable :: elem_id
    real(kind=DOUBLE), dimension(:, :), allocatable :: sol
    INTEGER(kind=ENTIER), dimension(:, :), allocatable :: sol_int
  end type mpi_elems_type

  type mpi_send_recv_type
    integer(kind=ENTIER) :: n_mpi_send_neigh, n_mpi_recv_neigh
    type(mpi_elems_type), dimension(:), allocatable :: mpi_send_neigh, mpi_recv_neigh
    integer, dimension(:), allocatable :: mpi_reqsend, mpi_reqrecv
    integer, dimension(:, :), allocatable :: mpi_sendstat, mpi_recvstat
    logical, dimension(:), allocatable :: is_ghost
  end type mpi_send_recv_type

contains
  subroutine mpi_memory_exchange(mpi_send_recv, n_elems, n_var, sol)
    use mpi
    implicit none

    integer(kind=ENTIER) :: n_elems, n_var
    real(kind=DOUBLE), dimension(n_var, n_elems), intent(inout) :: sol
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    integer :: mpi_ierr
    integer(kind=ENTIER) :: i, k
    integer(kind=ENTIER) :: id_elem

    do i = 1, mpi_send_recv%n_mpi_send_neigh
      allocate(mpi_send_recv%mpi_send_neigh(i)%sol(n_var, &
        mpi_send_recv%mpi_send_neigh(i)%n_elems))
      do k = 1, mpi_send_recv%mpi_send_neigh(i)%n_elems
        id_elem = mpi_send_recv%mpi_send_neigh(i)%elem_id(k)
        mpi_send_recv%mpi_send_neigh(i)%sol(:, k) = sol(:, id_elem)
      end do

      call mpi_isend(mpi_send_recv%mpi_send_neigh(i)%sol(1, 1), &
        n_var*mpi_send_recv%mpi_send_neigh(i)%n_elems, MPI_DOUBLE, &
        mpi_send_recv%mpi_send_neigh(i)%partition_id, &
        0, MPI_COMM_WORLD, mpi_send_recv%mpi_reqsend(i), mpi_ierr)
    end do

    do i = 1, mpi_send_recv%n_mpi_recv_neigh
      allocate(mpi_send_recv%mpi_recv_neigh(i)%sol(n_var, &
        mpi_send_recv%mpi_recv_neigh(i)%n_elems))
      call mpi_irecv(mpi_send_recv%mpi_recv_neigh(i)%sol(1, 1), &
        n_var*mpi_send_recv%mpi_recv_neigh(i)%n_elems, MPI_DOUBLE, &
        mpi_send_recv%mpi_recv_neigh(i)%partition_id, &
        MPI_ANY_TAG, MPI_COMM_WORLD, mpi_send_recv%mpi_reqrecv(i), mpi_ierr)
    end do

    call mpi_waitall(mpi_send_recv%n_mpi_send_neigh, &
      mpi_send_recv%mpi_reqsend, mpi_send_recv%mpi_sendstat, mpi_ierr)
    call mpi_waitall(mpi_send_recv%n_mpi_recv_neigh, &
      mpi_send_recv%mpi_reqrecv, mpi_send_recv%mpi_recvstat, mpi_ierr)

    do i = 1, mpi_send_recv%n_mpi_recv_neigh
      do k = 1, mpi_send_recv%mpi_recv_neigh(i)%n_elems
        id_elem = mpi_send_recv%mpi_recv_neigh(i)%elem_id(k)
        sol(:, id_elem) = mpi_send_recv%mpi_recv_neigh(i)%sol(:, k)
      end do
    end do

    do i = 1, mpi_send_recv%n_mpi_send_neigh
      deallocate(mpi_send_recv%mpi_send_neigh(i)%sol)
    end do
    do i = 1, mpi_send_recv%n_mpi_recv_neigh
      deallocate(mpi_send_recv%mpi_recv_neigh(i)%sol)
    end do
  end subroutine mpi_memory_exchange

  subroutine mpi_memory_exchange_partial_int(mpi_send_recv, n_elems, n_var, sol, id)
    use mpi
    implicit none

    integer(kind=ENTIER) :: n_elems, n_var, id
    integer(kind=ENTIER), dimension(n_var, n_elems), intent(inout) :: sol
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    integer :: mpi_ierr
    integer(kind=ENTIER) :: i, k, me
    integer(kind=ENTIER) :: id_elem

    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    if( me == id ) then
      do i = 1, mpi_send_recv%n_mpi_send_neigh
        allocate(mpi_send_recv%mpi_send_neigh(i)%sol_int(n_var, &
          mpi_send_recv%mpi_send_neigh(i)%n_elems))

        do k = 1, mpi_send_recv%mpi_send_neigh(i)%n_elems
          id_elem = mpi_send_recv%mpi_send_neigh(i)%elem_id(k)
          mpi_send_recv%mpi_send_neigh(i)%sol_int(:, k) = sol(:, id_elem)
        end do

        call mpi_isend(mpi_send_recv%mpi_send_neigh(i)%sol_int(1, 1), &
          n_var*mpi_send_recv%mpi_send_neigh(i)%n_elems, MPI_INT, &
          mpi_send_recv%mpi_send_neigh(i)%partition_id, &
          0, MPI_COMM_WORLD, mpi_send_recv%mpi_reqsend(i), mpi_ierr)
      end do

      call mpi_waitall(mpi_send_recv%n_mpi_send_neigh, &
        mpi_send_recv%mpi_reqsend, mpi_send_recv%mpi_sendstat, mpi_ierr)

      do i = 1, mpi_send_recv%n_mpi_send_neigh
        deallocate(mpi_send_recv%mpi_send_neigh(i)%sol_int)
      end do
    else

      do i = 1, mpi_send_recv%n_mpi_recv_neigh
        if ( mpi_send_recv%mpi_recv_neigh(i)%partition_id == id ) then
          allocate(mpi_send_recv%mpi_recv_neigh(i)%sol_int(n_var, &
            mpi_send_recv%mpi_recv_neigh(i)%n_elems))
          call mpi_irecv(mpi_send_recv%mpi_recv_neigh(i)%sol_int(1, 1), &
            n_var*mpi_send_recv%mpi_recv_neigh(i)%n_elems, MPI_INT, &
            mpi_send_recv%mpi_recv_neigh(i)%partition_id, &
            MPI_ANY_TAG, MPI_COMM_WORLD, mpi_send_recv%mpi_reqrecv(i), mpi_ierr)
        end if
      end do

      do i = 1, mpi_send_recv%n_mpi_recv_neigh
        if ( mpi_send_recv%mpi_recv_neigh(i)%partition_id == id ) then
          call mpi_wait(mpi_send_recv%mpi_reqrecv(i), mpi_send_recv%mpi_recvstat(:, i), mpi_ierr)

          do k = 1, mpi_send_recv%mpi_recv_neigh(i)%n_elems
            id_elem = mpi_send_recv%mpi_recv_neigh(i)%elem_id(k)
            sol(:, id_elem) = mpi_send_recv%mpi_recv_neigh(i)%sol_int(:, k)
          end do

        end if
      end do

      do i = 1, mpi_send_recv%n_mpi_recv_neigh
        if ( mpi_send_recv%mpi_recv_neigh(i)%partition_id == id ) then
          deallocate(mpi_send_recv%mpi_recv_neigh(i)%sol_int)
        end if
      end do
    end if

    call mpi_barrier(mpi_comm_world, mpi_ierr)
  end subroutine mpi_memory_exchange_partial_int
end module mpi_module
