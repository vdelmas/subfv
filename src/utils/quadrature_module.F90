! Quadrature rules for faces and volumes.
!
! Conventions:
!   - coords(3, n_nodes): physical vertex coordinates of the element
!   - weights include the Jacobian (sum of wts = area or volume)
!   - order: 1 = exact for degree 1, 2 = exact for degree 2, 3 = exact for degree 3
!     (order is clamped to the highest available rule if exceeded)
!
! Supported faces (face_quad_pts, n_nodes):
!   3 → triangle, 4 → bilinear quad
!
! Supported volumes (volume_quad_pts, n_nodes):
!   4 → tetrahedron, 5 → pyramid (split into 2 tets)
!   6 → prism/wedge,  8 → hexahedron
!
! Sub-face quadrature: a sub-face is always a quad — pass its 4 vertices
!   [cpm, cp, cpp, cf] to face_quad_pts with n_nodes=4.
!
! Sub-element quadrature: decompose into tetrahedra whose apex is the element
!   centroid ce and whose base is each triangular half of a sub-face; call
!   volume_quad_pts with n_nodes=4 on each tetrahedron and accumulate.
module quadrature_module
  use precision_module
  implicit none
  private

  public :: n_face_quad_pts, face_quad_pts
  public :: n_volume_quad_pts, volume_quad_pts
  public :: n_face_gl_quad_pts, face_gl_quad_pts
  public :: n_face_gl_iso_quad_pts, face_gl_iso_quad_pts

contains

  ! ----------------------------------------------------------------
  ! Query: number of quadrature points before allocating
  ! ----------------------------------------------------------------

  pure function n_face_quad_pts(n_nodes, order) result(n)
    integer(kind=ENTIER), intent(in) :: n_nodes, order
    integer(kind=ENTIER) :: n
    select case (n_nodes)
    case (3)
      select case (order)
      case (1);     n = 1
      case (2);     n = 3
      case default; n = 4   ! order >= 3
      end select
    case (4)
      select case (order)
      case (1);     n = 1
      case (2);     n = 4
      case default; n = 9   ! order >= 3
      end select
    case default; n = 0
    end select
  end function n_face_quad_pts

  pure function n_volume_quad_pts(n_nodes, order) result(n)
    integer(kind=ENTIER), intent(in) :: n_nodes, order
    integer(kind=ENTIER) :: n
    integer(kind=ENTIER) :: n_tet
    select case (n_nodes)
    case (4)   ! tetrahedron
      select case (order)
      case (1);     n = 1
      case (2);     n = 4
      case default; n = 5   ! order >= 3
      end select
    case (5)   ! pyramid = 2 tets
      select case (order)
      case (1);     n = 2 * 1
      case (2);     n = 2 * 4
      case default; n = 2 * 5
      end select
    case (6)   ! prism = triangle x line
      select case (order)
      case (1);     n = 1   ! 1 tri pt x 1 line pt
      case (2);     n = 6   ! 3 tri pts x 2 line pts
      case default; n = 12  ! 4 tri pts x 3 line pts
      end select
    case (8)   ! hexahedron = line x line x line
      select case (order)
      case (1);     n = 1
      case (2);     n = 8
      case default; n = 27   ! order >= 3
      end select
    case default; n = 0
    end select
  end function n_volume_quad_pts

  ! ----------------------------------------------------------------
  ! face_quad_pts: physical quadrature on a face polygon
  ! ----------------------------------------------------------------

  subroutine face_quad_pts(n_nodes, coords, order, pts, wts)
    integer(kind=ENTIER), intent(in) :: n_nodes, order
    real(kind=DOUBLE), dimension(3, n_nodes), intent(in) :: coords
    real(kind=DOUBLE), dimension(:, :), allocatable, intent(out) :: pts
    real(kind=DOUBLE), dimension(:),    allocatable, intent(out) :: wts
    integer(kind=ENTIER) :: n

    n = n_face_quad_pts(n_nodes, order)
    allocate(pts(3, n), wts(n))
    select case (n_nodes)
    case (3); call tri_face_rule(coords, order, pts, wts)
    case (4); call quad_face_rule(coords, order, pts, wts)
    end select
  end subroutine face_quad_pts

  ! ----------------------------------------------------------------
  ! volume_quad_pts: physical quadrature on a volume element
  ! ----------------------------------------------------------------

  subroutine volume_quad_pts(n_nodes, coords, order, pts, wts)
    integer(kind=ENTIER), intent(in) :: n_nodes, order
    real(kind=DOUBLE), dimension(3, n_nodes), intent(in) :: coords
    real(kind=DOUBLE), dimension(:, :), allocatable, intent(out) :: pts
    real(kind=DOUBLE), dimension(:),    allocatable, intent(out) :: wts
    integer(kind=ENTIER) :: n

    n = n_volume_quad_pts(n_nodes, order)
    allocate(pts(3, n), wts(n))
    select case (n_nodes)
    case (4); call tet_rule(coords, order, pts, wts)
    case (5); call pyr_rule(coords, order, pts, wts)
    case (6); call prism_rule(coords, order, pts, wts)
    case (8); call hex_rule(coords, order, pts, wts)
    end select
  end subroutine volume_quad_pts

  ! ================================================================
  ! Triangle face rule (Dunavant)
  !   vertices: coords(:,1..3)
  !   x_gauss = L1*v1 + L2*v2 + L3*v3
  !   w_gauss = w_norm * area
  ! ================================================================
  subroutine tri_face_rule(coords, order, pts, wts)
    real(kind=DOUBLE), dimension(3, 3), intent(in)  :: coords
    integer(kind=ENTIER),               intent(in)  :: order
    real(kind=DOUBLE), dimension(:, :), intent(out) :: pts
    real(kind=DOUBLE), dimension(:),    intent(out) :: wts

    real(kind=DOUBLE), dimension(3) :: v1, v2, v3
    real(kind=DOUBLE) :: area

    v1 = coords(:, 1); v2 = coords(:, 2); v3 = coords(:, 3)
    area = 0.5_DOUBLE * norm2(cross3(v2-v1, v3-v1))

    select case (order)

    case (1)   ! 1-point rule, exact for degree 1
      pts(:, 1) = (v1 + v2 + v3) / 3.0_DOUBLE
      wts(1) = area

    case (2)   ! 3-point Dunavant, exact for degree 2
      !   (L1, L2, L3) = (2/3, 1/6, 1/6) and cyclic permutations, w = 1/3
      pts(:, 1) = (4*v1 +   v2 +   v3) / 6.0_DOUBLE
      pts(:, 2) = (  v1 + 4*v2 +   v3) / 6.0_DOUBLE
      pts(:, 3) = (  v1 +   v2 + 4*v3) / 6.0_DOUBLE
      wts = area / 3.0_DOUBLE

    case default   ! 4-point Dunavant, exact for degree 3
      ! (1/3, 1/3, 1/3), w = -27/48
      pts(:, 1) = (v1 + v2 + v3) / 3.0_DOUBLE
      wts(1) = -27.0_DOUBLE / 48.0_DOUBLE * area
      ! (3/5, 1/5, 1/5) and cyclic, w = 25/48
      pts(:, 2) = (3*v1 +   v2 +   v3) / 5.0_DOUBLE
      pts(:, 3) = (  v1 + 3*v2 +   v3) / 5.0_DOUBLE
      pts(:, 4) = (  v1 +   v2 + 3*v3) / 5.0_DOUBLE
      wts(2:4) = 25.0_DOUBLE / 48.0_DOUBLE * area

    end select
  end subroutine tri_face_rule

  ! ================================================================
  ! Bilinear quad face rule (Gauss-Legendre)
  !   vertices in order: v1(--), v2(+-), v3(++), v4(-+) in (xi,eta)
  !   x(xi,eta) = sum_k N_k(xi,eta) * v_k,  xi,eta in [-1,1]
  !   w_gauss = w_xi * w_eta * |dx/dxi x dx/deta|
  ! ================================================================
  subroutine quad_face_rule(coords, order, pts, wts)
    real(kind=DOUBLE), dimension(3, 4), intent(in)  :: coords
    integer(kind=ENTIER),               intent(in)  :: order
    real(kind=DOUBLE), dimension(:, :), intent(out) :: pts
    real(kind=DOUBLE), dimension(:),    intent(out) :: wts

    ! 1D Gauss nodes and weights on [-1,1]
    real(kind=DOUBLE), dimension(3) :: xi_g, w_g
    integer(kind=ENTIER) :: ng, i, j, k
    real(kind=DOUBLE) :: xi, eta, wx, wy
    real(kind=DOUBLE), dimension(3) :: dxdxi, dxdeta
    real(kind=DOUBLE), dimension(4) :: N

    call gauss1d(order, ng, xi_g, w_g)

    k = 0
    do i = 1, ng
      do j = 1, ng
        xi  = xi_g(i); wx = w_g(i)
        eta = xi_g(j); wy = w_g(j)
        k = k + 1
        call quad_shape(xi, eta, N, dxdxi, dxdeta, coords)
        pts(:, k) = matmul(coords, N)
        wts(k) = wx * wy * norm2(cross3(dxdxi, dxdeta))
      end do
    end do
  end subroutine quad_face_rule

  ! ================================================================
  ! Tetrahedron rule
  !   vertices: coords(:,1..4)
  !   x = L1*v1+L2*v2+L3*v3+L4*v4, L1=1-L2-L3-L4
  !   Jacobian: J = [v2-v1 | v3-v1 | v4-v1], constant
  !   vol = |det J| / 6
  !   w_gauss = w_norm * vol
  ! ================================================================
  subroutine tet_rule(coords, order, pts, wts)
    real(kind=DOUBLE), dimension(3, 4), intent(in)  :: coords
    integer(kind=ENTIER),               intent(in)  :: order
    real(kind=DOUBLE), dimension(:, :), intent(out) :: pts
    real(kind=DOUBLE), dimension(:),    intent(out) :: wts

    real(kind=DOUBLE), dimension(3) :: v1, v2, v3, v4
    real(kind=DOUBLE) :: vol, a, b

    v1 = coords(:, 1); v2 = coords(:, 2)
    v3 = coords(:, 3); v4 = coords(:, 4)
    vol = abs(det3(reshape([v2-v1, v3-v1, v4-v1], [3, 3]))) / 6.0_DOUBLE

    select case (order)

    case (1)   ! 1 pt, centroid
      pts(:, 1) = (v1+v2+v3+v4) / 4.0_DOUBLE
      wts(1) = vol

    case (2)   ! 4 pts, Dunavant–Keast
      ! a = (5+3√5)/20,  b = (5-√5)/20
      a = (5.0_DOUBLE + 3.0_DOUBLE*sqrt(5.0_DOUBLE)) / 20.0_DOUBLE
      b = (5.0_DOUBLE -            sqrt(5.0_DOUBLE)) / 20.0_DOUBLE
      pts(:, 1) = a*v1 + b*v2 + b*v3 + b*v4
      pts(:, 2) = b*v1 + a*v2 + b*v3 + b*v4
      pts(:, 3) = b*v1 + b*v2 + a*v3 + b*v4
      pts(:, 4) = b*v1 + b*v2 + b*v3 + a*v4
      wts = vol / 4.0_DOUBLE

    case default   ! 5 pts, Keast order 3  (Felippa Table 9.4)
      ! centroid with negative weight, then 4 face-biased pts
      pts(:, 1) = (v1+v2+v3+v4) / 4.0_DOUBLE
      wts(1) = -4.0_DOUBLE / 5.0_DOUBLE * vol
      ! L = (1/2, 1/6, 1/6, 1/6) and permutations
      pts(:, 2) = 0.5_DOUBLE*v1 + (v2+v3+v4)/6.0_DOUBLE
      pts(:, 3) = 0.5_DOUBLE*v2 + (v1+v3+v4)/6.0_DOUBLE
      pts(:, 4) = 0.5_DOUBLE*v3 + (v1+v2+v4)/6.0_DOUBLE
      pts(:, 5) = 0.5_DOUBLE*v4 + (v1+v2+v3)/6.0_DOUBLE
      wts(2:5) = 9.0_DOUBLE / 20.0_DOUBLE * vol

    end select
  end subroutine tet_rule

  ! ================================================================
  ! Pyramid rule: split into 2 tetrahedra
  !   GMSH/VTK node order: v1..v4 base (quad), v5 apex
  !   Tet 1: v1, v2, v3, v5
  !   Tet 2: v1, v3, v4, v5
  ! ================================================================
  subroutine pyr_rule(coords, order, pts, wts)
    real(kind=DOUBLE), dimension(3, 5), intent(in)  :: coords
    integer(kind=ENTIER),               intent(in)  :: order
    real(kind=DOUBLE), dimension(:, :), intent(out) :: pts
    real(kind=DOUBLE), dimension(:),    intent(out) :: wts

    real(kind=DOUBLE), dimension(3, 4) :: tet_coords
    integer(kind=ENTIER) :: n_tet
    real(kind=DOUBLE), dimension(:, :), allocatable :: pts_t
    real(kind=DOUBLE), dimension(:),    allocatable :: wts_t

    n_tet = n_volume_quad_pts(4_ENTIER, order)

    tet_coords(:, 1:3) = coords(:, 1:3); tet_coords(:, 4) = coords(:, 5)
    call tet_rule(tet_coords, order, pts_t, wts_t)
    pts(:, 1:n_tet)           = pts_t; wts(1:n_tet)           = wts_t
    deallocate(pts_t, wts_t)

    tet_coords(:, 1) = coords(:, 1); tet_coords(:, 2) = coords(:, 3)
    tet_coords(:, 3) = coords(:, 4); tet_coords(:, 4) = coords(:, 5)
    call tet_rule(tet_coords, order, pts_t, wts_t)
    pts(:, n_tet+1:2*n_tet) = pts_t; wts(n_tet+1:2*n_tet) = wts_t
    deallocate(pts_t, wts_t)
  end subroutine pyr_rule

  ! ================================================================
  ! Prism (wedge) rule: tensor product triangle x line
  !   GMSH node order: v1,v2,v3 bottom triangle, v4,v5,v6 top
  !   x(L,t) = (1-t)/2*(L1*v1+L2*v2+L3*v3) + (1+t)/2*(L1*v4+L2*v5+L3*v6)
  !   t in [-1,1], L barycentric
  !   Physical weight = w_tri_norm/2 * w_line * |det J(L,t)|
  ! ================================================================
  subroutine prism_rule(coords, order, pts, wts)
    real(kind=DOUBLE), dimension(3, 6), intent(in)  :: coords
    integer(kind=ENTIER),               intent(in)  :: order
    real(kind=DOUBLE), dimension(:, :), intent(out) :: pts
    real(kind=DOUBLE), dimension(:),    intent(out) :: wts

    ! Triangle Gauss reference data (barycentric, weights sum to 1)
    integer(kind=ENTIER), parameter :: max_tri = 4, max_line = 3
    integer(kind=ENTIER) :: n_tri, n_line, i, j, k
    real(kind=DOUBLE), dimension(3, max_tri) :: L_tri
    real(kind=DOUBLE), dimension(max_tri)    :: w_tri
    real(kind=DOUBLE), dimension(max_line)   :: t_line, w_line
    real(kind=DOUBLE) :: L1, L2, L3, t, wl
    real(kind=DOUBLE), dimension(3) :: xbot, xtop
    real(kind=DOUBLE), dimension(3) :: dxdL1, dxdL2, dxdt
    real(kind=DOUBLE) :: jac

    call tri_ref_rule(order, n_tri, L_tri, w_tri)
    call gauss1d(order, n_line, t_line, w_line)

    k = 0
    do i = 1, n_tri
      L1 = L_tri(1, i); L2 = L_tri(2, i); L3 = L_tri(3, i)
      do j = 1, n_line
        t  = t_line(j); wl = w_line(j)
        k  = k + 1
        xbot = L1*coords(:,1) + L2*coords(:,2) + L3*coords(:,3)
        xtop = L1*coords(:,4) + L2*coords(:,5) + L3*coords(:,6)
        pts(:, k) = 0.5_DOUBLE*((1-t)*xbot + (1+t)*xtop)
        ! Jacobian columns
        dxdL1 = 0.5_DOUBLE*((1-t)*(coords(:,1)-coords(:,3)) &
                            +(1+t)*(coords(:,4)-coords(:,6)))
        dxdL2 = 0.5_DOUBLE*((1-t)*(coords(:,2)-coords(:,3)) &
                            +(1+t)*(coords(:,5)-coords(:,6)))
        dxdt  = 0.5_DOUBLE*(-xbot + xtop)
        jac   = abs(det3(reshape([dxdL1, dxdL2, dxdt], [3,3])))
        ! Factor 1/2: area of reference triangle in (L1,L2) coords
        wts(k) = 0.5_DOUBLE * w_tri(i) * wl * jac
      end do
    end do
  end subroutine prism_rule

  ! ================================================================
  ! Hexahedron rule: tensor product Gauss-Legendre
  !   GMSH node order (xi,eta,zeta):
  !   v1(-,-,-), v2(+,-,-), v3(+,+,-), v4(-,+,-)
  !   v5(-,-,+), v6(+,-,+), v7(+,+,+), v8(-,+,+)
  ! ================================================================
  subroutine hex_rule(coords, order, pts, wts)
    real(kind=DOUBLE), dimension(3, 8), intent(in)  :: coords
    integer(kind=ENTIER),               intent(in)  :: order
    real(kind=DOUBLE), dimension(:, :), intent(out) :: pts
    real(kind=DOUBLE), dimension(:),    intent(out) :: wts

    integer(kind=ENTIER) :: ng, i, j, k, m
    real(kind=DOUBLE), dimension(3) :: xi_g, w_g
    real(kind=DOUBLE) :: xi, eta, zeta, wx, wy, wz, jac
    real(kind=DOUBLE), dimension(3) :: dxdxi, dxdeta, dxdzeta
    real(kind=DOUBLE), dimension(8) :: N
    integer(kind=ENTIER) :: idx

    call gauss1d(order, ng, xi_g, w_g)

    idx = 0
    do i = 1, ng
      do j = 1, ng
        do k = 1, ng
          xi   = xi_g(i); wx = w_g(i)
          eta  = xi_g(j); wy = w_g(j)
          zeta = xi_g(k); wz = w_g(k)
          idx  = idx + 1
          call hex_shape(xi, eta, zeta, N, dxdxi, dxdeta, dxdzeta, coords)
          pts(:, idx) = matmul(coords, N)
          jac = abs(det3(reshape([dxdxi, dxdeta, dxdzeta], [3,3])))
          wts(idx) = wx * wy * wz * jac
        end do
      end do
    end do
  end subroutine hex_rule

  ! ================================================================
  ! Internal helpers
  ! ================================================================

  ! 1D Gauss-Legendre nodes and weights on [-1,1], n=1,2,3 points
  subroutine gauss1d(order, ng, xi, w)
    integer(kind=ENTIER), intent(in)  :: order
    integer(kind=ENTIER), intent(out) :: ng
    real(kind=DOUBLE), dimension(3), intent(out) :: xi, w

    real(kind=DOUBLE), parameter :: s3 = 1.0_DOUBLE / sqrt(3.0_DOUBLE)
    real(kind=DOUBLE), parameter :: s35 = sqrt(3.0_DOUBLE / 5.0_DOUBLE)

    select case (order)
    case (1)
      ng = 1; xi(1) = 0.0_DOUBLE; w(1) = 2.0_DOUBLE
    case (2)
      ng = 2
      xi(1) = -s3; xi(2) = s3
      w(1)  = 1.0_DOUBLE; w(2) = 1.0_DOUBLE
    case default   ! order >= 3 → 3-point rule
      ng = 3
      xi(1) = -s35; xi(2) = 0.0_DOUBLE; xi(3) = s35
      w(1)  = 5.0_DOUBLE/9.0_DOUBLE
      w(2)  = 8.0_DOUBLE/9.0_DOUBLE
      w(3)  = 5.0_DOUBLE/9.0_DOUBLE
    end select
  end subroutine gauss1d

  ! Triangle reference Gauss rule in barycentric coords (L1,L2,L3)
  ! w_tri are normalized so sum(w_tri) = 1
  subroutine tri_ref_rule(order, n, L, w)
    integer(kind=ENTIER), intent(in)  :: order
    integer(kind=ENTIER), intent(out) :: n
    real(kind=DOUBLE), dimension(3, 4), intent(out) :: L
    real(kind=DOUBLE), dimension(4),    intent(out) :: w

    select case (order)
    case (1)
      n = 1
      L(:, 1) = [1.0_DOUBLE/3, 1.0_DOUBLE/3, 1.0_DOUBLE/3]
      w(1) = 1.0_DOUBLE
    case (2)
      n = 3
      L(:, 1) = [2.0_DOUBLE/3, 1.0_DOUBLE/6, 1.0_DOUBLE/6]
      L(:, 2) = [1.0_DOUBLE/6, 2.0_DOUBLE/3, 1.0_DOUBLE/6]
      L(:, 3) = [1.0_DOUBLE/6, 1.0_DOUBLE/6, 2.0_DOUBLE/3]
      w(1:3) = 1.0_DOUBLE/3
    case default
      n = 4
      L(:, 1) = [1.0_DOUBLE/3, 1.0_DOUBLE/3, 1.0_DOUBLE/3]
      L(:, 2) = [3.0_DOUBLE/5, 1.0_DOUBLE/5, 1.0_DOUBLE/5]
      L(:, 3) = [1.0_DOUBLE/5, 3.0_DOUBLE/5, 1.0_DOUBLE/5]
      L(:, 4) = [1.0_DOUBLE/5, 1.0_DOUBLE/5, 3.0_DOUBLE/5]
      w(1)   = -27.0_DOUBLE/48
      w(2:4) =  25.0_DOUBLE/48
    end select
  end subroutine tri_ref_rule

  ! Bilinear quad shape functions and physical Jacobian columns
  ! Node order: v1(-,-), v2(+,-), v3(+,+), v4(-,+)
  subroutine quad_shape(xi, eta, N, dxdxi, dxdeta, coords)
    real(kind=DOUBLE), intent(in)  :: xi, eta
    real(kind=DOUBLE), dimension(4), intent(out) :: N
    real(kind=DOUBLE), dimension(3), intent(out) :: dxdxi, dxdeta
    real(kind=DOUBLE), dimension(3, 4), intent(in) :: coords

    real(kind=DOUBLE), dimension(4) :: dNdxi, dNdeta

    N(1) = (1-xi)*(1-eta)/4.0_DOUBLE
    N(2) = (1+xi)*(1-eta)/4.0_DOUBLE
    N(3) = (1+xi)*(1+eta)/4.0_DOUBLE
    N(4) = (1-xi)*(1+eta)/4.0_DOUBLE

    dNdxi(1) = -(1-eta)/4.0_DOUBLE; dNdxi(2) =  (1-eta)/4.0_DOUBLE
    dNdxi(3) =  (1+eta)/4.0_DOUBLE; dNdxi(4) = -(1+eta)/4.0_DOUBLE

    dNdeta(1) = -(1-xi)/4.0_DOUBLE; dNdeta(2) = -(1+xi)/4.0_DOUBLE
    dNdeta(3) =  (1+xi)/4.0_DOUBLE; dNdeta(4) =  (1-xi)/4.0_DOUBLE

    dxdxi  = matmul(coords, dNdxi)
    dxdeta = matmul(coords, dNdeta)
  end subroutine quad_shape

  ! Trilinear hex shape functions and Jacobian columns
  ! Node order: v1(---)...v8(-++) see GMSH convention
  subroutine hex_shape(xi, eta, zeta, N, dxdxi, dxdeta, dxdzeta, coords)
    real(kind=DOUBLE), intent(in) :: xi, eta, zeta
    real(kind=DOUBLE), dimension(8), intent(out) :: N
    real(kind=DOUBLE), dimension(3), intent(out) :: dxdxi, dxdeta, dxdzeta
    real(kind=DOUBLE), dimension(3, 8), intent(in) :: coords

    real(kind=DOUBLE), dimension(8) :: dNdxi, dNdeta, dNdzeta
    real(kind=DOUBLE) :: xm, xp, ym, yp, zm, zp

    xm=1-xi; xp=1+xi; ym=1-eta; yp=1+eta; zm=1-zeta; zp=1+zeta

    N(1)=xm*ym*zm/8; N(2)=xp*ym*zm/8; N(3)=xp*yp*zm/8; N(4)=xm*yp*zm/8
    N(5)=xm*ym*zp/8; N(6)=xp*ym*zp/8; N(7)=xp*yp*zp/8; N(8)=xm*yp*zp/8

    dNdxi(1)=-ym*zm/8; dNdxi(2)= ym*zm/8; dNdxi(3)= yp*zm/8; dNdxi(4)=-yp*zm/8
    dNdxi(5)=-ym*zp/8; dNdxi(6)= ym*zp/8; dNdxi(7)= yp*zp/8; dNdxi(8)=-yp*zp/8

    dNdeta(1)=-xm*zm/8; dNdeta(2)=-xp*zm/8; dNdeta(3)= xp*zm/8; dNdeta(4)= xm*zm/8
    dNdeta(5)=-xm*zp/8; dNdeta(6)=-xp*zp/8; dNdeta(7)= xp*zp/8; dNdeta(8)= xm*zp/8

    dNdzeta(1)=-xm*ym/8; dNdzeta(2)=-xp*ym/8; dNdzeta(3)=-xp*yp/8; dNdzeta(4)=-xm*yp/8
    dNdzeta(5)= xm*ym/8; dNdzeta(6)= xp*ym/8; dNdzeta(7)= xp*yp/8; dNdzeta(8)= xm*yp/8

    dxdxi   = matmul(coords, dNdxi)
    dxdeta  = matmul(coords, dNdeta)
    dxdzeta = matmul(coords, dNdzeta)
  end subroutine hex_shape

  ! Cross product of two 3-vectors
  pure function cross3(a, b) result(c)
    real(kind=DOUBLE), dimension(3), intent(in) :: a, b
    real(kind=DOUBLE), dimension(3) :: c
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross3

  ! Determinant of a 3x3 matrix
  pure function det3(A) result(d)
    real(kind=DOUBLE), dimension(3, 3), intent(in) :: A
    real(kind=DOUBLE) :: d
    d = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) &
      - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) &
      + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
  end function det3


  ! ================================================================
  ! Gauss-Lobatto vertex rules: quadrature nodes at face vertices
  !
  ! face_gl_quad_pts (n_nodes=3, order=1):
  !   3 points at triangle vertices, w = A/3 each. Exact for P_1.
  !
  ! face_gl_quad_pts (n_nodes=4, order=1):
  !   4 points at quad vertices, weights from triangle split
  !   T1=[v1,v2,v3] + T2=[v1,v3,v4]. Exact for P_1 on any quad.
  !   w1 = (A1+A2)/3,  w2 = A1/3,  w3 = (A1+A2)/3,  w4 = A2/3
  !
  ! face_gl_iso_quad_pts: iso variant — vertex i_iso has weight iso_wt
  !   (iso_wt <= natural GL weight at that vertex, enforced externally).
  !   Other weights adjusted to maintain P_1 exactness.
  !
  !   n_nodes=3: adds 1 extra interior point at the face centroid;
  !     the 3 remaining unknowns (2 free vertices + centroid weight)
  !     are solved from the 3 P_1 conditions (∫1, ∫c1, ∫c2).
  !
  !   n_nodes=4: no extra point needed; fixing iso_wt uses the
  !     one free DOF of the underdetermined GL system, giving a unique
  !     solution for the other 3 vertex weights from 3 P_1 conditions.
  !
  ! The 2 in-plane coordinates (c1, c2) are chosen automatically from
  ! the face normal (largest-component rule) for 3D robustness.
  ! ================================================================

  pure function n_face_gl_quad_pts(n_nodes, order) result(n)
    integer(kind=ENTIER), intent(in) :: n_nodes, order
    integer(kind=ENTIER) :: n
    select case (n_nodes)
    case (3); n = 3
    case (4); n = 4
    case default; n = 0
    end select
  end function n_face_gl_quad_pts

  pure function n_face_gl_iso_quad_pts(n_nodes, order) result(n)
    integer(kind=ENTIER), intent(in) :: n_nodes, order
    integer(kind=ENTIER) :: n
    select case (n_nodes)
    case (3); n = 4   ! 3 vertices + 1 centroid interior point
    case (4); n = 4   ! 4 vertices, iso constraint fills the free DOF
    case default; n = 0
    end select
  end function n_face_gl_iso_quad_pts

  subroutine face_gl_quad_pts(n_nodes, coords, order, pts, wts)
    integer(kind=ENTIER), intent(in) :: n_nodes, order
    real(kind=DOUBLE), dimension(3, n_nodes), intent(in) :: coords
    real(kind=DOUBLE), dimension(:, :), allocatable, intent(out) :: pts
    real(kind=DOUBLE), dimension(:),    allocatable, intent(out) :: wts
    integer(kind=ENTIER) :: n

    n = n_face_gl_quad_pts(n_nodes, order)
    allocate(pts(3, n), wts(n))
    select case (n_nodes)
    case (3); call tri_gl_rule(coords, pts, wts)
    case (4); call quad_gl_rule(coords, pts, wts)
    end select
  end subroutine face_gl_quad_pts

  subroutine face_gl_iso_quad_pts(n_nodes, coords, order, i_iso, iso_wt, pts, wts)
    integer(kind=ENTIER), intent(in) :: n_nodes, order, i_iso
    real(kind=DOUBLE),    intent(in) :: iso_wt
    real(kind=DOUBLE), dimension(3, n_nodes), intent(in) :: coords
    real(kind=DOUBLE), dimension(:, :), allocatable, intent(out) :: pts
    real(kind=DOUBLE), dimension(:),    allocatable, intent(out) :: wts
    integer(kind=ENTIER) :: n

    n = n_face_gl_iso_quad_pts(n_nodes, order)
    allocate(pts(3, n), wts(n))
    select case (n_nodes)
    case (3); call tri_gl_iso_rule(coords, i_iso, iso_wt, pts, wts)
    case (4); call quad_gl_iso_rule(coords, i_iso, iso_wt, pts, wts)
    end select
  end subroutine face_gl_iso_quad_pts

  ! ----------------------------------------------------------------
  ! Triangle GL vertex rule: w = A/3 at each vertex. Exact for P_1.
  ! ----------------------------------------------------------------
  subroutine tri_gl_rule(coords, pts, wts)
    real(kind=DOUBLE), dimension(3, 3), intent(in)  :: coords
    real(kind=DOUBLE), dimension(3, 3), intent(out) :: pts
    real(kind=DOUBLE), dimension(3),    intent(out) :: wts
    real(kind=DOUBLE) :: area

    area = 0.5_DOUBLE * norm2(cross3(coords(:,2)-coords(:,1), coords(:,3)-coords(:,1)))
    pts  = coords
    wts  = area / 3.0_DOUBLE
  end subroutine tri_gl_rule

  ! ----------------------------------------------------------------
  ! Quad GL vertex rule via triangle split T1=[v1,v2,v3]+T2=[v1,v3,v4].
  ! w1=(A1+A2)/3, w2=A1/3, w3=(A1+A2)/3, w4=A2/3. Exact for P_1.
  ! ----------------------------------------------------------------
  subroutine quad_gl_rule(coords, pts, wts)
    real(kind=DOUBLE), dimension(3, 4), intent(in)  :: coords
    real(kind=DOUBLE), dimension(3, 4), intent(out) :: pts
    real(kind=DOUBLE), dimension(4),    intent(out) :: wts
    real(kind=DOUBLE) :: a1, a2

    a1 = 0.5_DOUBLE * norm2(cross3(coords(:,2)-coords(:,1), coords(:,3)-coords(:,1)))
    a2 = 0.5_DOUBLE * norm2(cross3(coords(:,3)-coords(:,1), coords(:,4)-coords(:,1)))
    pts    = coords
    wts(1) = (a1 + a2) / 3.0_DOUBLE
    wts(2) = a1 / 3.0_DOUBLE
    wts(3) = (a1 + a2) / 3.0_DOUBLE
    wts(4) = a2 / 3.0_DOUBLE
  end subroutine quad_gl_rule

  ! ----------------------------------------------------------------
  ! Triangle GL iso rule: vertex i_iso has weight iso_wt.
  ! Adds 1 extra interior point (face centroid). Solves 3x3 system
  ! for the 2 free vertex weights + centroid weight.
  ! ----------------------------------------------------------------
  subroutine tri_gl_iso_rule(coords, i_iso, iso_wt, pts, wts)
    real(kind=DOUBLE), dimension(3, 3), intent(in)  :: coords
    integer(kind=ENTIER),               intent(in)  :: i_iso
    real(kind=DOUBLE),                  intent(in)  :: iso_wt
    real(kind=DOUBLE), dimension(3, 4), intent(out) :: pts
    real(kind=DOUBLE), dimension(4),    intent(out) :: wts

    real(kind=DOUBLE) :: area
    real(kind=DOUBLE), dimension(3) :: centroid, int_xyz, sol
    integer(kind=ENTIER) :: ic1, ic2, j, k, free(2)
    real(kind=DOUBLE), dimension(3, 3) :: mat
    real(kind=DOUBLE), dimension(3) :: rhs

    area     = 0.5_DOUBLE * norm2(cross3(coords(:,2)-coords(:,1), coords(:,3)-coords(:,1)))
    centroid = (coords(:,1) + coords(:,2) + coords(:,3)) / 3.0_DOUBLE
    int_xyz  = area * centroid

    call face_coord_pair(cross3(coords(:,2)-coords(:,1), coords(:,3)-coords(:,1)), ic1, ic2)

    ! Identify the 2 free vertex indices
    k = 0
    do j = 1, 3
      if (j /= i_iso) then
        k = k + 1
        free(k) = j
      end if
    end do

    ! Place points: 3 vertices then centroid (index 4)
    pts(:, 1:3) = coords
    pts(:, 4)   = centroid
    wts(i_iso)  = iso_wt

    ! 3x3 system for w_{free(1)}, w_{free(2)}, w_centroid
    mat(1, 1) = 1.0_DOUBLE
    mat(1, 2) = 1.0_DOUBLE
    mat(1, 3) = 1.0_DOUBLE
    mat(2, 1) = coords(ic1, free(1))
    mat(2, 2) = coords(ic1, free(2))
    mat(2, 3) = centroid(ic1)
    mat(3, 1) = coords(ic2, free(1))
    mat(3, 2) = coords(ic2, free(2))
    mat(3, 3) = centroid(ic2)

    rhs(1) = area         - iso_wt
    rhs(2) = int_xyz(ic1) - iso_wt * coords(ic1, i_iso)
    rhs(3) = int_xyz(ic2) - iso_wt * coords(ic2, i_iso)

    call solve3x3(mat, rhs, sol)

    wts(free(1)) = sol(1)
    wts(free(2)) = sol(2)
    wts(4)       = sol(3)
  end subroutine tri_gl_iso_rule

  ! ----------------------------------------------------------------
  ! Quad GL iso rule: vertex i_iso has weight iso_wt.
  ! Solves 3x3 system for the 3 remaining vertex weights.
  ! No extra point needed (iso constraint fills the one free DOF).
  ! ----------------------------------------------------------------
  subroutine quad_gl_iso_rule(coords, i_iso, iso_wt, pts, wts)
    real(kind=DOUBLE), dimension(3, 4), intent(in)  :: coords
    integer(kind=ENTIER),               intent(in)  :: i_iso
    real(kind=DOUBLE),                  intent(in)  :: iso_wt
    real(kind=DOUBLE), dimension(3, 4), intent(out) :: pts
    real(kind=DOUBLE), dimension(4),    intent(out) :: wts

    real(kind=DOUBLE) :: a1, a2
    real(kind=DOUBLE), dimension(3) :: int_xyz, sol
    integer(kind=ENTIER) :: ic1, ic2, j, k, free(3)
    real(kind=DOUBLE), dimension(3, 3) :: mat
    real(kind=DOUBLE), dimension(3) :: rhs

    ! Integrals via triangle split T1=[v1,v2,v3] + T2=[v1,v3,v4]
    a1 = 0.5_DOUBLE * norm2(cross3(coords(:,2)-coords(:,1), coords(:,3)-coords(:,1)))
    a2 = 0.5_DOUBLE * norm2(cross3(coords(:,3)-coords(:,1), coords(:,4)-coords(:,1)))
    int_xyz = a1 * (coords(:,1) + coords(:,2) + coords(:,3)) / 3.0_DOUBLE &
            + a2 * (coords(:,1) + coords(:,3) + coords(:,4)) / 3.0_DOUBLE

    call face_coord_pair(cross3(coords(:,2)-coords(:,1), coords(:,3)-coords(:,1)), ic1, ic2)

    pts       = coords
    wts(i_iso) = iso_wt

    k = 0
    do j = 1, 4
      if (j /= i_iso) then
        k = k + 1
        free(k) = j
      end if
    end do

    mat(1, 1) = 1.0_DOUBLE
    mat(1, 2) = 1.0_DOUBLE
    mat(1, 3) = 1.0_DOUBLE
    mat(2, 1) = coords(ic1, free(1))
    mat(2, 2) = coords(ic1, free(2))
    mat(2, 3) = coords(ic1, free(3))
    mat(3, 1) = coords(ic2, free(1))
    mat(3, 2) = coords(ic2, free(2))
    mat(3, 3) = coords(ic2, free(3))

    rhs(1) = (a1 + a2)    - iso_wt
    rhs(2) = int_xyz(ic1) - iso_wt * coords(ic1, i_iso)
    rhs(3) = int_xyz(ic2) - iso_wt * coords(ic2, i_iso)

    call solve3x3(mat, rhs, sol)

    do k = 1, 3
      wts(free(k)) = sol(k)
    end do
  end subroutine quad_gl_iso_rule

  ! ----------------------------------------------------------------
  ! Pick the 2 in-plane coordinate indices from a face normal vector.
  ! Avoids the direction where the normal is largest (most degenerate).
  ! ----------------------------------------------------------------
  subroutine face_coord_pair(normal, ic1, ic2)
    real(kind=DOUBLE), dimension(3), intent(in)  :: normal
    integer(kind=ENTIER),            intent(out) :: ic1, ic2
    real(kind=DOUBLE), dimension(3) :: an

    an = abs(normal)
    if (an(1) >= an(2) .and. an(1) >= an(3)) then
      ic1 = 2; ic2 = 3
    else if (an(2) >= an(1) .and. an(2) >= an(3)) then
      ic1 = 1; ic2 = 3
    else
      ic1 = 1; ic2 = 2
    end if
  end subroutine face_coord_pair

  ! ----------------------------------------------------------------
  ! 3x3 linear system via Cramer's rule: A*x = b
  ! ----------------------------------------------------------------
  subroutine solve3x3(A, b, x)
    real(kind=DOUBLE), dimension(3, 3), intent(in)  :: A
    real(kind=DOUBLE), dimension(3),    intent(in)  :: b
    real(kind=DOUBLE), dimension(3),    intent(out) :: x
    real(kind=DOUBLE) :: d
    real(kind=DOUBLE), dimension(3, 3) :: Ak

    d = det3(A)
    Ak = A; Ak(:, 1) = b; x(1) = det3(Ak) / d
    Ak = A; Ak(:, 2) = b; x(2) = det3(Ak) / d
    Ak = A; Ak(:, 3) = b; x(3) = det3(Ak) / d
  end subroutine solve3x3

end module quadrature_module
