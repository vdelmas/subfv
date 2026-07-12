module ns_euler_primitives_module
  use precision_module
  use mesh_module
  implicit none

contains
  pure function primit_to_conserv(w) result(u)
    use ns_global_data_module, only: gamma
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE), dimension(5) :: u

    u(1) = w(1)
    u(2:4) = w(2:4)*w(1)
    u(5) = w(5)/(gamma - 1) &
      + 0.5_DOUBLE*w(1)*(w(2)**2 + w(3)**2 + w(4)**2)
  end function primit_to_conserv

  pure function conserv_to_primit(u) result(w)
    use ns_global_data_module, only: gamma
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: u
    real(kind=DOUBLE), dimension(5) :: w

    w(1) = u(1)
    w(2:4) = u(2:4)/u(1)
    w(5) = (gamma - 1)*(u(5) &
      - 0.5_DOUBLE*w(1)*(w(2)**2 + w(3)**2 + w(4)**2))
  end function conserv_to_primit

  pure function sound_speed_w(w) result(a)
    use ns_global_data_module, only: gamma
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE) :: a

    a = sqrt(gamma*w(5)/w(1))
  end function sound_speed_w

  pure function sound_speed_u(u) result(a)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: u
    real(kind=DOUBLE) :: a

    a = sound_speed_w(conserv_to_primit(u))
  end function sound_speed_u

  pure function temp_u(u) result(theta)
    use ns_global_data_module, only: gamma, Cv_p
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: u
    real(kind=DOUBLE) :: theta

    real(kind=DOUBLE) :: tau

    tau = 1.0_DOUBLE/u(1)
    theta = (u(5)*tau - 0.5_DOUBLE*dot_product(u(2:4),u(2:4))*tau**2)/Cv_p
  end function temp_u

  pure function temp_w(w) result(theta)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE) :: theta

    theta = temp_u(primit_to_conserv(w))
  end function temp_w

  pure function mach_u(u) result(mach)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: u
    real(kind=DOUBLE) :: mach

    mach = mach_w(conserv_to_primit(u))
  end function mach_u

  pure function mach_w(w) result(mach)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE) :: mach

    mach = norm2(w(2:4))/sound_speed_w(w)
  end function mach_w

  pure function max_omega_for_positivity(U, deltaU) result(omega)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: U, deltaU
    real(kind=DOUBLE) :: omega

    real(kind=DOUBLE), dimension(5) :: WS, US
    logical :: pos

    omega = 1.0_DOUBLE
    US = U + omega * deltaU
    WS = conserv_to_primit(US)
    pos = (WS(1) > 1e-17_DOUBLE .and. WS(5) > 1e-17_DOUBLE)
    do while (.not. pos .and. omega > 1e-6_DOUBLE)
      omega = omega * 0.9_DOUBLE
      US = U + omega * deltaU
      WS = conserv_to_primit(US)
      pos = (WS(1) > 1e-17_DOUBLE .and. WS(5) > 1e-17_DOUBLE)
    end do
  end function max_omega_for_positivity

  pure subroutine base_change(v1, v2, n, forward)
    implicit none

    logical, intent(in) :: forward
    real(kind=DOUBLE), dimension(3), intent(inout) :: v1, v2
    real(kind=DOUBLE), dimension(3), intent(in) :: n

    real(kind=DOUBLE) :: e1(3), e2(3), e3(3), tmp(3), norm
    real(kind=DOUBLE) :: t1(3), t2(3)

    norm = sqrt(n(1)*n(1) + n(2)*n(2) + n(3)*n(3))
    e1 = n/norm

    if (abs(e1(1)) <= abs(e1(2)) .and. abs(e1(1)) <= abs(e1(3))) then
      tmp = [1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE]
    else if (abs(e1(2)) <= abs(e1(3))) then
      tmp = [0.0_DOUBLE, 1.0_DOUBLE, 0.0_DOUBLE]
    else
      tmp = [0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE]
    end if

    e2 = [tmp(2)*e1(3) - tmp(3)*e1(2), &
      tmp(3)*e1(1) - tmp(1)*e1(3), &
      tmp(1)*e1(2) - tmp(2)*e1(1)]
    e2 = e2/sqrt(e2(1)*e2(1) + e2(2)*e2(2) + e2(3)*e2(3))
    e3 = [e1(2)*e2(3) - e1(3)*e2(2), &
      e1(3)*e2(1) - e1(1)*e2(3), &
      e1(1)*e2(2) - e1(2)*e2(1)]

    if (forward) then
      t1 = [v1(1)*e1(1) + v1(2)*e1(2) + v1(3)*e1(3), &
        v1(1)*e2(1) + v1(2)*e2(2) + v1(3)*e2(3), &
        v1(1)*e3(1) + v1(2)*e3(2) + v1(3)*e3(3)]
      t2 = [v2(1)*e1(1) + v2(2)*e1(2) + v2(3)*e1(3), &
        v2(1)*e2(1) + v2(2)*e2(2) + v2(3)*e2(3), &
        v2(1)*e3(1) + v2(2)*e3(2) + v2(3)*e3(3)]
    else
      t1 = [e1(1)*v1(1) + e2(1)*v1(2) + e3(1)*v1(3), &
        e1(2)*v1(1) + e2(2)*v1(2) + e3(2)*v1(3), &
        e1(3)*v1(1) + e2(3)*v1(2) + e3(3)*v1(3)]
      t2 = [e1(1)*v2(1) + e2(1)*v2(2) + e3(1)*v2(3), &
        e1(2)*v2(1) + e2(2)*v2(2) + e3(2)*v2(3), &
        e1(3)*v2(1) + e2(3)*v2(2) + e3(3)*v2(3)]
    end if

    v1 = t1
    v2 = t2
  end subroutine base_change

  function wall_normal(mesh, id_vert) result(Bp)
    use ns_global_data_module, only: bc_euler_id, BC_EULER_WALL, boundary_2d
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert

    real(kind=DOUBLE), dimension(3) :: Bp
    integer(kind=ENTIER) :: j, re
    integer(kind=ENTIER) :: id_sub_face, id_face
    real(kind=DOUBLE), dimension(3) :: norm

    Bp = 0.0_DOUBLE
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      re = mesh%face(id_face)%right_neigh
      if (re <= 0) then
        norm = mesh%sub_face(id_sub_face)%norm
        if (.not. (boundary_2d .and. abs(norm(3)) > 1e-12)) then
          if (re == 0) then
            Bp = Bp + mesh%sub_face(id_sub_face)%area*mesh%sub_face(id_sub_face)%norm
          else if (bc_euler_id(-re) == BC_EULER_WALL) then
            Bp = Bp + mesh%sub_face(id_sub_face)%area*mesh%sub_face(id_sub_face)%norm
          end if
        end if
      end if
    end do
  end function wall_normal
end module ns_euler_primitives_module