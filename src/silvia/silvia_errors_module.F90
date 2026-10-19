module silvia_errors_module
	use precision_module
  use mesh_module
	use silvia_base_module
	implicit none

	public :: compute_error_test
	public :: compute_error_isentropic
	public :: compute_error_gresho

contains

	subroutine compute_error_test(mesh, qpw, sol, t, mach, volume, h_err, l2_err)
		type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(4,6,mesh%n_elems), intent(in) :: qpw
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), intent(in) :: t
		real(kind=DOUBLE), intent(in) :: mach
		real(kind=DOUBLE), intent(inout) :: h_err, volume
    real(kind=DOUBLE), dimension(5), intent(inout) :: l2_err

		if (init_uniform) then
      ! still nothign to do
    else if (init_1drp) then
      ! still nothign to do
    else if (init_isentropic_vortex) then
			call compute_error_isentropic(mesh,qpw,sol,t,volume,h_err,l2_err)
    else if (init_gresho) then
      call compute_error_gresho(mesh,qpw,sol,mach,h_err,l2_err)
    end if

	end subroutine compute_error_test

  subroutine compute_error_isentropic(mesh, qpw, sol, t, volume, diam, error)
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(4,6,mesh%n_elems), intent(in) :: qpw
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), intent(in) :: t
    real(kind=DOUBLE), dimension(5), intent(inout) :: error
		real(kind=DOUBLE), intent(inout) :: volume, diam
    
    integer(kind=ENTIER) :: i,q
    real(kind=DOUBLE), dimension(5) :: wexact, wsol, wq, uexact_avg

    error = 0.0_DOUBLE
    volume = 0.0_DOUBLE
    diam = 1.0e5
    do i = 1, mesh%n_elems
      if (mesh%elem(i)%is_ghost) then
        cycle
      else 
        if (abs(mesh%elem(i)%coord(1)) < 3.0_DOUBLE .and. abs(mesh%elem(i)%coord(2)) < 3.0_DOUBLE .and. abs(mesh%elem(i)%coord(3)) < 3.0_DOUBLE) then
          if (order<=2) then
            call sol_isentropic_vortex(mesh%elem(i)%coord, wexact, t)
          else
            ! sol stores conservative cell averages at order 3. Average the
            ! exact conservative state first, then convert that average to
            ! primitive variables for a like-for-like comparison.
            uexact_avg = 0.0_DOUBLE
            do q = 1, size(qpw,DIM=2)
              call sol_isentropic_vortex(qpw(1:3,q,i), wq, t)
              uexact_avg = uexact_avg + qpw(4,q,i) * primit_to_conserv(wq)
            end do
            uexact_avg = uexact_avg / mesh%elem(i)%volume
            wexact = conserv_to_primit(uexact_avg)
          end if
          wsol = conserv_to_primit(sol(:, i))
          error = error + mesh%elem(i)%volume*abs(wsol - wexact)
          volume = volume + mesh%elem(i)%volume
          if (error_2d) then
            diam = min( diam , sqrt(mesh%elem(i)%volume) )
          else 
            diam = min( diam , (mesh%elem(i)%volume)**(1.0_DOUBLE/3.0_DOUBLE) )
          end if
        end if
      end if
    end do
  end subroutine compute_error_isentropic

	subroutine compute_error_gresho(mesh, qpw, sol, mach, volume, error)
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(4,6,mesh%n_elems), intent(in) :: qpw
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
		real(kind=DOUBLE), intent(in) :: mach
		real(kind=DOUBLE), dimension(5), intent(inout) :: error
    real(kind=DOUBLE), intent(inout) :: volume

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
        call sol_gresho_mach(mesh%elem(i)%coord, wexact, mach)
        wsol = conserv_to_primit(sol(:, i))
        error = error + mesh%elem(i)%volume*abs((wsol - wexact)**1)
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
