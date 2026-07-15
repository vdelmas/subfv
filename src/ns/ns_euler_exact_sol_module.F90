module ns_euler_exact_sol_module
  use precision_module
  use mesh_module
  implicit none

contains
  pure subroutine sol_gresho(x, w)
    use ns_global_data_module, only: gamma
    implicit none

    real(kind=DOUBLE), dimension(3), intent(in) :: x
    real(kind=DOUBLE), dimension(5), intent(inout) :: w

    real(kind=DOUBLE) :: w_gresho, m_gresho, p_gresho, r
    real(kind=DOUBLE), dimension(3) :: coord2

    w_gresho = 0.2_DOUBLE
    m_gresho = 1e-5
    p_gresho = 1.0_DOUBLE/(gamma*m_gresho**2)

    coord2(:) = x - (/0._DOUBLE, 0._DOUBLE, 0._DOUBLE/)
    r = norm2(coord2(:2))

    if (r < w_gresho) then
      w(1) = 1.0_DOUBLE
      w(2) = -5.0_DOUBLE*coord2(2)
      w(3) = 5.0_DOUBLE*coord2(1)
      w(4) = 0.0_DOUBLE
      w(5) = p_gresho + 12.5*r**2
    else if (r > 1.0_DOUBLE*w_gresho .and. r < 2.0_DOUBLE*w_gresho) then
      w(1) = 1.0_DOUBLE
      w(2) = (5.0_DOUBLE - 2.0_DOUBLE/r)*coord2(2)
      w(3) = (-5.0_DOUBLE + 2.0_DOUBLE/r)*coord2(1)
      w(4) = 0.0_DOUBLE
      w(5) = p_gresho + 12.5_DOUBLE*r**2 + &
        4 - 20.0_DOUBLE*r + 4.0_DOUBLE*log(5.0_DOUBLE*r)
    else
      w(1) = 1.0_DOUBLE
      w(2) = 0.0_DOUBLE
      w(3) = 0.0_DOUBLE
      w(4) = 0.0_DOUBLE
      w(5) = p_gresho - 2.0_DOUBLE + 4.0_DOUBLE*log(2.0_DOUBLE)
    end if
    w(2) = w(2) + 100000.0_DOUBLE
  end subroutine sol_gresho

  pure subroutine sol_isentropic_vortex(coord, w, t)
    use ns_global_data_module, only:gamma, pi
    implicit none

    real(kind=DOUBLE), intent(in) :: t
    real(kind=DOUBLE), dimension(3), intent(in) :: coord
    real(kind=DOUBLE), dimension(5), intent(inout) :: w

    real(kind=DOUBLE), dimension(3) :: center_coord, vel
    real(kind=DOUBLE) :: beta, r

    beta = 5.0_DOUBLE
    vel(:) = (/0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE/)
    center_coord(:) = (/0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE/) + vel*t
    r = (coord(1) - center_coord(1))**2 + (coord(2) - center_coord(2))**2
    w(1) = 1.0_DOUBLE*(1.0_DOUBLE - ((gamma - 1)*beta**2)/(8.0_DOUBLE*gamma*pi**2)* &
      exp(1.0_DOUBLE - r))**(1.0_DOUBLE/(gamma - 1.0_DOUBLE))
    w(2) = vel(1) &
      - (coord(2) - center_coord(2))*beta/(2.0_DOUBLE*pi)*exp(0.5_DOUBLE*(1.0_DOUBLE - r))
    w(3) = vel(2) &
      + (coord(1) - center_coord(1))*beta/(2.0_DOUBLE*pi)*exp(0.5_DOUBLE*(1.0_DOUBLE - r))
    w(4) = 0.0_DOUBLE
    w(5) = w(1)**gamma
  end subroutine sol_isentropic_vortex

  subroutine compute_error_isentropic(mesh, sol, t)
    use ns_global_data_module, only: error_2d, error_2d_h
    use ns_euler_primitives_module, only: conserv_to_primit
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), intent(in) :: t

    integer(kind=ENTIER) :: i
    real(kind=DOUBLE) :: error, volume
    real(kind=DOUBLE), dimension(5) :: wexact, wsol

    error = 0.0_DOUBLE
    volume = 0.0_DOUBLE
    do i = 1, mesh%n_elems
      if (abs(mesh%elem(i)%coord(1)) < 3.0_DOUBLE &
        .and. abs(mesh%elem(i)%coord(2)) < 3.0_DOUBLE &
        .and. abs(mesh%elem(i)%coord(3)) < 3.0_DOUBLE) then
        call sol_isentropic_vortex(mesh%elem(i)%coord, wexact, t)
        wsol = conserv_to_primit(sol(:, i))
        error = error + mesh%elem(i)%volume*(wsol(1) - wexact(1))**2
        volume = volume + mesh%elem(i)%volume
      end if
    end do

    if (error_2d) then
      print *, "Error Vortex: ", sqrt((volume/error_2d_h)/mesh%n_elems), sqrt(error)
    else
      print *, "Error Vortex: ", (volume/mesh%n_elems)**(1.0_DOUBLE/3.0_DOUBLE), sqrt(error)
    end if
  end subroutine compute_error_isentropic

  pure subroutine sol_potential_flow_2d(x, w)
    use ns_global_data_module, only: gamma
    implicit none

    real(kind=DOUBLE), dimension(3), intent(in) :: x
    real(kind=DOUBLE), dimension(5), intent(inout) :: w

    real(kind=DOUBLE), parameter :: r0 = 0.5_DOUBLE, r1 = 5.5_DOUBLE
    real(kind=DOUBLE) :: r, theta, vr, vtheta, ux, uy, vmag2

    real(kind=DOUBLE), parameter :: rhoinf = 1.4_DOUBLE
    real(kind=DOUBLE), parameter :: mach = 1e-4
    real(kind=DOUBLE), parameter :: pinf = 1.0_DOUBLE
    real(kind=DOUBLE) :: vinf

    vinf = mach*sqrt(gamma*pinf/rhoinf)

    r     = sqrt(x(1)**2 + x(2)**2)
    theta = atan2(x(2), x(1))

    vr     = (r1**2/(r1**2 - r0**2)) * (1.0_DOUBLE - (r0/r)**2) * cos(theta)
    vtheta = -(r1**2/(r1**2 - r0**2)) * (1.0_DOUBLE + (r0/r)**2) * sin(theta)

    ux = vinf * (cos(theta)*vr - sin(theta)*vtheta)
    uy = vinf * (sin(theta)*vr + cos(theta)*vtheta)

    vmag2 = ux**2 + uy**2

    w(1) = rhoinf
    w(2) = ux
    w(3) = uy
    w(4) = 0.0_DOUBLE
    w(5) = pinf + 0.5_DOUBLE*rhoinf*(vinf**2 - vmag2)
  end subroutine sol_potential_flow_2d

  pure subroutine sol_potential_flow_3d(x, w)
    use ns_global_data_module, only: gamma
    implicit none

    real(kind=DOUBLE), dimension(3), intent(in) :: x
    real(kind=DOUBLE), dimension(5), intent(inout) :: w

    real(kind=DOUBLE), parameter :: r0 = 0.5_DOUBLE
    real(kind=DOUBLE) :: r, f, ux, uy, uz, vmag2

    real(kind=DOUBLE), parameter :: rhoinf = 1.4_DOUBLE
    real(kind=DOUBLE), parameter :: mach = 1e-4_DOUBLE
    real(kind=DOUBLE), parameter :: pinf = 1.0_DOUBLE
    real(kind=DOUBLE) :: vinf

    vinf = mach*sqrt(gamma*pinf/rhoinf)

    r = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
    f = r0**3 / (2.0_DOUBLE * r**5)

    ux = vinf * (1.0_DOUBLE + f*(r**2 - 3.0_DOUBLE*x(1)**2))
    uy = -3.0_DOUBLE * vinf * f * x(1) * x(2)
    uz = -3.0_DOUBLE * vinf * f * x(1) * x(3)
    vmag2 = ux**2 + uy**2 + uz**2

    w(1) = rhoinf
    w(2) = ux
    w(3) = uy
    w(4) = uz
    w(5) = pinf + 0.5_DOUBLE*rhoinf*(vinf**2 - vmag2)
  end subroutine sol_potential_flow_3d

  function kelvin_helmholtz(y) result(v)
    use ns_global_data_module, only: pi
    implicit none

    real(kind=DOUBLE), intent(in) :: y
    real(kind=DOUBLE) :: v

    real(kind=DOUBLE) :: w_kh = 1.0_DOUBLE / 16.0_DOUBLE

    if (y > -0.25_DOUBLE - w_kh/2.0_DOUBLE .and. &
        y < -0.25_DOUBLE + w_kh/2.0_DOUBLE) then
      v = -sin(pi * (y + 0.25_DOUBLE) / w_kh)
    else if (y > -0.25_DOUBLE + w_kh/2.0_DOUBLE .and. &
             y <  0.25_DOUBLE - w_kh/2.0_DOUBLE) then
      v = -1.0_DOUBLE
    else if (y >  0.25_DOUBLE - w_kh/2.0_DOUBLE .and. &
             y <  0.25_DOUBLE + w_kh/2.0_DOUBLE) then
      v =  sin(pi * (y - 0.25_DOUBLE) / w_kh)
    else
      v = 1.0_DOUBLE
    end if
  end function kelvin_helmholtz

end module ns_euler_exact_sol_module