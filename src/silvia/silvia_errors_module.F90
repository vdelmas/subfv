module silvia_errors_module
	use precision_module
  use mesh_module
	use silvia_base_module
	implicit none

	public :: compute_error_test
	public :: compute_error_isentropic
	public :: compute_error_gresho

contains

	subroutine compute_error_test(mesh, sol, t, mach, h_err, l2_err)
		type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), intent(in) :: t
		real(kind=DOUBLE), intent(in) :: mach
		real(kind=DOUBLE), intent(inout) :: h_err, l2_err

		if (init_uniform) then
      ! still nothign to do
    else if (init_1drp) then
      ! still nothign to do
    else if (init_isentropic_vortex) then
			call compute_error_isentropic(mesh,sol,t,h_err,l2_err)
    else if (init_gresho) then
      call compute_error_gresho(mesh,sol,mach,h_err,l2_err)
    end if

	end subroutine compute_error_test

  subroutine compute_error_isentropic(mesh, sol, t, volume, error)
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), intent(in) :: t
		real(kind=DOUBLE), intent(inout) :: error, volume
    
    integer(kind=ENTIER) :: i
    real(kind=DOUBLE), dimension(5) :: wexact, wsol

    error = 0.0_DOUBLE
    volume = 0.0_DOUBLE
    do i = 1, mesh%n_elems
      if (mesh%elem(i)%is_ghost) then
        cycle
      else if (abs(mesh%elem(i)%coord(1)) < 3.0_DOUBLE &
        .and. abs(mesh%elem(i)%coord(2)) < 3.0_DOUBLE &
        .and. abs(mesh%elem(i)%coord(3)) < 3.0_DOUBLE) then
        call sol_isentropic_vortex(mesh%elem(i)%coord, wexact, t)
        wsol = conserv_to_primit(sol(:, i))
        error = error + mesh%elem(i)%volume*(wsol(1) - wexact(1))**2
        volume = volume + mesh%elem(i)%volume
      end if
    end do

    ! if (error_2d) then
    !   print *, "Error Vortex: ", sqrt((volume/error_2d_h)/mesh%n_elems), sqrt(error)
    ! else
    !   print *, "Error Vortex: ", (volume/mesh%n_elems)**(1.0_DOUBLE/3.0_DOUBLE), sqrt(error)
    ! end if
  end subroutine compute_error_isentropic

	subroutine compute_error_gresho(mesh, sol, mach, volume, error)
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
		real(kind=DOUBLE), intent(in) :: mach
		real(kind=DOUBLE), intent(inout) :: error, volume

    integer(kind=ENTIER) :: i
    real(kind=DOUBLE), dimension(5) :: wexact, wsol

    error = 0.0_DOUBLE
    volume = 0.0_DOUBLE
    do i = 1, mesh%n_elems
      if (mesh%elem(i)%is_ghost) then
        cycle
      ! else if (abs(mesh%elem(i)%coord(1)) < 0.35_DOUBLE &
      !   .and. abs(mesh%elem(i)%coord(2)) < 0.35_DOUBLE &
      !   .and. abs(mesh%elem(i)%coord(3)) < 0.35_DOUBLE) then
      else
        call sol_gresho_mach_C2(mesh%elem(i)%coord, wexact, mach)
        wsol = conserv_to_primit(sol(:, i))
        error = error + mesh%elem(i)%volume*(wsol(1) - wexact(1))**2
        volume = volume + mesh%elem(i)%volume
      end if
    end do

    ! if (error_2d) then
    !   print *, "Error Gresho: ", sqrt((volume/error_2d_h)/mesh%n_elems), sqrt(error)
    ! else
    !   print *, "Error Gresho: ", (volume/mesh%n_elems)**(1.0_DOUBLE/3.0_DOUBLE), sqrt(error)
    ! end if
  end subroutine compute_error_gresho

end module silvia_errors_module