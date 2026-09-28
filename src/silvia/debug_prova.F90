program main
    implicit none

    integer :: i
    real :: x

    x = 1.0

    do i = 1, 5
        x = x * 2.0
        print *, "i =", i, "x =", x
    end do
end program main
