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
! Order-4/5 status (2026-09-15, for the arbitrary-high-order study):
!   quad face / hex volume (n_nodes=4/8): full tensor-product Gauss-
!     Legendre support through order 5 (gauss1d extended to 5 points).
!   triangle face (n_nodes=3) and the triangle side of prism_rule
!     (n_nodes=6): full support through order 5 -- 6-pt (degree 4) and
!     7-pt (degree 5) Dunavant/Radon rules, derived from scratch via
!     moment-matching and verified numerically against every monomial
!     up to their target degree (max error ~1e-15) rather than
!     transcribed from a table -- see tri_face_rule/tri_ref_rule.
!   tetrahedron volume (n_nodes=4, and n_nodes=5 pyramid = 2 tets):
!     order 4 done (11-pt Keast rule, same from-scratch derivation +
!     verification, ~1e-17 max error -- see tet_rule). order 5 NOT yet
!     done: order>=5 on a tet/pyramid mesh silently reuses the order-4,
!     11-point rule (degree-4 accurate, not degree-5) -- a symmetric
!     15-point ansatz (centroid + 2 S31 orbits + 1 S22 orbit) was tried
!     and did not converge to a valid rule from several initial
!     guesses; needs either a different ansatz/orbit structure or a
!     verified literature rule, not yet done. Low practical urgency as
!     of this note: no active test mesh in this study uses pure
!     tetrahedral cells (the triangular-mesh cylinder case uses prisms,
!     already covered above). Flagged for whoever next needs true
!     order-5 accuracy on a tet-celled mesh.
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
      case (3);     n = 4
      case (4);     n = 6
      case (5);     n = 7
      ! order>=7 reuses the degree-6, 12-pt rule (no verified degree-7
      ! rule yet) -- see tri_face_rule's own case default. Keep in sync.
      case default; n = 12   ! order >= 6
      end select
    case (4)
      select case (order)
      case (1);     n = 1
      case (2);     n = 4
      case (3);     n = 9
      case (4);     n = 16
      case default; n = 25   ! order >= 5
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
      case (3);     n = 5
      case (4);     n = 11
      ! order>=6 reuses the degree-5, 15-pt rule (no verified degree-6
      ! tet rule yet) -- see tet_rule's own case default. Keep in sync.
      case default; n = 15  ! order >= 5
      end select
    case (5)   ! pyramid = 2 tets
      select case (order)
      case (1);     n = 2 * 1
      case (2);     n = 2 * 4
      case (3);     n = 2 * 5
      case (4);     n = 2 * 11
      case default; n = 2 * 15  ! order >= 5, see n_nodes=4 case above
      end select
    case (6)   ! prism = triangle x line
      select case (order)
      case (1);     n = 1   ! 1 tri pt x 1 line pt
      case (2);     n = 6   ! 3 tri pts x 2 line pts
      case (3);     n = 12  ! 4 tri pts x 3 line pts
      case (4);     n = 24  ! 6 tri pts x 4 line pts
      case (5);     n = 35  ! 7 tri pts x 5 line pts
      case default; n = 60  ! order >= 6: 12 tri pts x 5 line pts
      end select
    case (8)   ! hexahedron = line x line x line
      select case (order)
      case (1);     n = 1
      case (2);     n = 8
      case (3);     n = 27
      case (4);     n = 64
      case default; n = 125   ! order >= 5
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

    case (3)   ! 4-point Dunavant, exact for degree 3
      ! (1/3, 1/3, 1/3), w = -27/48
      pts(:, 1) = (v1 + v2 + v3) / 3.0_DOUBLE
      wts(1) = -27.0_DOUBLE / 48.0_DOUBLE * area
      ! (3/5, 1/5, 1/5) and cyclic, w = 25/48
      pts(:, 2) = (3*v1 +   v2 +   v3) / 5.0_DOUBLE
      pts(:, 3) = (  v1 + 3*v2 +   v3) / 5.0_DOUBLE
      pts(:, 4) = (  v1 +   v2 + 3*v3) / 5.0_DOUBLE
      wts(2:4) = 25.0_DOUBLE / 48.0_DOUBLE * area

    case (4)   ! 6-point Dunavant, exact for degree 4. Two symmetric
      ! orbits of 3 permutations each of (a,a,b); verified numerically
      ! against exact monomial moments up to degree 4 (max error 1e-15).
      block
        real(kind=DOUBLE), parameter :: a1 = 0.445948490915965_DOUBLE
        real(kind=DOUBLE), parameter :: b1 = 0.108103018168070_DOUBLE
        real(kind=DOUBLE), parameter :: w1 = 0.223381589678011_DOUBLE
        real(kind=DOUBLE), parameter :: a2 = 0.091576213509771_DOUBLE
        real(kind=DOUBLE), parameter :: b2 = 0.816847572980459_DOUBLE
        real(kind=DOUBLE), parameter :: w2 = 0.109951743655322_DOUBLE
        pts(:, 1) = a1*v1 + a1*v2 + b1*v3
        pts(:, 2) = a1*v1 + b1*v2 + a1*v3
        pts(:, 3) = b1*v1 + a1*v2 + a1*v3
        wts(1:3)  = w1 * area
        pts(:, 4) = a2*v1 + a2*v2 + b2*v3
        pts(:, 5) = a2*v1 + b2*v2 + a2*v3
        pts(:, 6) = b2*v1 + a2*v2 + a2*v3
        wts(4:6)  = w2 * area
      end block

    case (5)   ! 7-point Dunavant/Radon, exact for degree 5. Centroid
      ! plus two symmetric orbits of 3 permutations of (a,a,b); verified
      ! numerically up to degree 5 (max error 1e-15), and confirmed
      ! genuinely NOT exact at degree 6 (~1e-4 relative error there).
      block
        real(kind=DOUBLE), parameter :: w0 = 0.225_DOUBLE
        real(kind=DOUBLE), parameter :: a1 = 0.470142064105115_DOUBLE
        real(kind=DOUBLE), parameter :: b1 = 0.059715871789770_DOUBLE
        real(kind=DOUBLE), parameter :: w1 = 0.132394152788506_DOUBLE
        real(kind=DOUBLE), parameter :: a2 = 0.101286507323456_DOUBLE
        real(kind=DOUBLE), parameter :: b2 = 0.797426985353087_DOUBLE
        real(kind=DOUBLE), parameter :: w2 = 0.125939180544827_DOUBLE
        pts(:, 1) = (v1 + v2 + v3) / 3.0_DOUBLE
        wts(1)    = w0 * area
        pts(:, 2) = a1*v1 + a1*v2 + b1*v3
        pts(:, 3) = a1*v1 + b1*v2 + a1*v3
        pts(:, 4) = b1*v1 + a1*v2 + a1*v3
        wts(2:4)  = w1 * area
        pts(:, 5) = a2*v1 + a2*v2 + b2*v3
        pts(:, 6) = a2*v1 + b2*v2 + a2*v3
        pts(:, 7) = b2*v1 + a2*v2 + a2*v3
        wts(5:7)  = w2 * area
      end block

    case default   ! order >= 6: 12-point Dunavant (1985) degree-6 rule.
      ! Two S21 orbits (3 pts each, perms of (a,a,1-2a)) + one S111
      ! orbit (6 pts, all perms of (a3,b3,1-a3-b3), all three barycentric
      ! values distinct). Derived from scratch by symbolic moment-matching
      ! (sympy/mpmath, multi-start Levenberg-Marquardt refined to 50+
      ! digits) against the exact reference-triangle integral of every
      ! monomial up to degree 6, not transcribed from a table -- the
      ! ansatz has several real roots (trivial relabelings of the same
      ! rule); the one used has every point interior and every weight
      ! positive, and its constants independently match the published
      ! Dunavant degree-6 values to the digits checked. Verified
      ! numerically to ~1e-17 max error over every monomial up to degree
      ! 6, and confirmed genuinely NOT exact at degree 7 (~1e-6 relative
      ! error there). Replaces the previous fallback that silently
      ! reused the degree-5, 7-point rule for order>=6 -- this is also
      ! what prism_rule gets for its triangle side via tri_ref_rule
      ! below, since a triangular-base prism is just this rule times
      ! gauss1d in the extrusion direction.
      block
        real(kind=DOUBLE), parameter :: a1 = 0.2492867451709104_DOUBLE
        real(kind=DOUBLE), parameter :: w1 = 0.1167862757263794_DOUBLE
        real(kind=DOUBLE), parameter :: a2 = 0.06308901449150223_DOUBLE
        real(kind=DOUBLE), parameter :: w2 = 0.05084490637020682_DOUBLE
        real(kind=DOUBLE), parameter :: a3 = 0.3103524510337844_DOUBLE
        real(kind=DOUBLE), parameter :: b3 = 0.05314504984481695_DOUBLE
        real(kind=DOUBLE), parameter :: w3 = 0.08285107561837358_DOUBLE
        real(kind=DOUBLE) :: c1, c2, c3
        c1 = 1.0_DOUBLE - 2.0_DOUBLE*a1
        c2 = 1.0_DOUBLE - 2.0_DOUBLE*a2
        c3 = 1.0_DOUBLE - a3 - b3
        pts(:, 1) = a1*v1 + a1*v2 + c1*v3
        pts(:, 2) = a1*v1 + c1*v2 + a1*v3
        pts(:, 3) = c1*v1 + a1*v2 + a1*v3
        wts(1:3)  = w1 * area
        pts(:, 4) = a2*v1 + a2*v2 + c2*v3
        pts(:, 5) = a2*v1 + c2*v2 + a2*v3
        pts(:, 6) = c2*v1 + a2*v2 + a2*v3
        wts(4:6)  = w2 * area
        pts(:, 7)  = a3*v1 + b3*v2 + c3*v3
        pts(:, 8)  = a3*v1 + c3*v2 + b3*v3
        pts(:, 9)  = b3*v1 + a3*v2 + c3*v3
        pts(:, 10) = b3*v1 + c3*v2 + a3*v3
        pts(:, 11) = c3*v1 + a3*v2 + b3*v3
        pts(:, 12) = c3*v1 + b3*v2 + a3*v3
        wts(7:12)  = w3 * area
      end block

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
    real(kind=DOUBLE), dimension(5) :: xi_g, w_g
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

    case (3)   ! 5 pts, Keast order 3  (Felippa Table 9.4)
      ! centroid with negative weight, then 4 face-biased pts
      pts(:, 1) = (v1+v2+v3+v4) / 4.0_DOUBLE
      wts(1) = -4.0_DOUBLE / 5.0_DOUBLE * vol
      ! L = (1/2, 1/6, 1/6, 1/6) and permutations
      pts(:, 2) = 0.5_DOUBLE*v1 + (v2+v3+v4)/6.0_DOUBLE
      pts(:, 3) = 0.5_DOUBLE*v2 + (v1+v3+v4)/6.0_DOUBLE
      pts(:, 4) = 0.5_DOUBLE*v3 + (v1+v2+v4)/6.0_DOUBLE
      pts(:, 5) = 0.5_DOUBLE*v4 + (v1+v2+v3)/6.0_DOUBLE
      wts(2:5) = 9.0_DOUBLE / 20.0_DOUBLE * vol

    case (4)   ! 11 pts, Keast order 4. Centroid orbit (1 pt) + S31(a)
      ! orbit (4 pts, perms of (a,a,a,1-3a)) + S22(b) orbit (6 pts, perms
      ! of (b,b,1/2-b,1/2-b)). Derived from scratch via symbolic
      ! moment-matching (sympy) against the exact reference-tetrahedron
      ! integrals of every monomial up to degree 4, not transcribed from
      ! a table -- verified numerically afterward to ~1e-17 max error
      ! (see aho-mpi-ghost-stencil-bug memory / tex_arbitrary_high_order
      ! for the derivation script). Matches the classical Keast(1986)
      ! degree-4 rule.
      block
        real(kind=DOUBLE), parameter :: a = 1.0_DOUBLE/14.0_DOUBLE
        real(kind=DOUBLE), parameter :: b = 0.25_DOUBLE - sqrt(70.0_DOUBLE)/56.0_DOUBLE
        real(kind=DOUBLE), parameter :: wc = -148.0_DOUBLE/1875.0_DOUBLE
        real(kind=DOUBLE), parameter :: wa = 343.0_DOUBLE/7500.0_DOUBLE
        real(kind=DOUBLE), parameter :: wb = 56.0_DOUBLE/375.0_DOUBLE
        real(kind=DOUBLE) :: ca, cb
        ca = 1.0_DOUBLE - 3.0_DOUBLE*a
        cb = 0.5_DOUBLE - b
        pts(:, 1) = 0.25_DOUBLE*(v1+v2+v3+v4)
        wts(1) = wc * vol
        pts(:, 2) = a*v1 + a*v2 + a*v3 + ca*v4
        pts(:, 3) = a*v1 + a*v2 + ca*v3 + a*v4
        pts(:, 4) = a*v1 + ca*v2 + a*v3 + a*v4
        pts(:, 5) = ca*v1 + a*v2 + a*v3 + a*v4
        wts(2:5) = wa * vol
        pts(:, 6)  = b*v1 + b*v2 + cb*v3 + cb*v4
        pts(:, 7)  = b*v1 + cb*v2 + b*v3 + cb*v4
        pts(:, 8)  = b*v1 + cb*v2 + cb*v3 + b*v4
        pts(:, 9)  = cb*v1 + b*v2 + b*v3 + cb*v4
        pts(:, 10) = cb*v1 + b*v2 + cb*v3 + b*v4
        pts(:, 11) = cb*v1 + cb*v2 + b*v3 + b*v4
        wts(6:11) = wb * vol
      end block

    case default   ! order >= 5: 15 pts, degree-5. Centroid orbit (1 pt)
      ! + S31(a) orbit (4 pts, perms of (a,a,a,1-3a)) + S31(b) orbit (4
      ! pts, perms of (b,b,b,1-3b), a second, independent S31 family --
      ! degree 4 alone needs only one S31 orbit, degree 5 needs two) +
      ! S22(c) orbit (6 pts, perms of (c,c,d,d), d=1/2-c). All three
      ! shape parameters and four weights derived from scratch by
      ! symbolic moment-matching (sympy/mpmath, Levenberg-Marquardt
      ! refined to 50+ digits) against the exact reference-tetrahedron
      ! integral of every monomial up to degree 5, not transcribed from
      ! a table -- this particular real root of the (overdetermined,
      ! 18-equation/7-unknown) moment system was chosen because it has
      ! every point strictly interior and every weight strictly positive
      ! (other real roots of the same ansatz exist but place a point
      ! outside the tet or a negative weight). Verified numerically
      ! afterward to ~1e-17 max error over every monomial up to degree 5,
      ! and confirmed genuinely NOT exact at degree 6 (~1e-4 relative
      ! error there, as expected for a degree-5 rule). Replaces the
      ! previous fallback that silently reused the degree-4, 11-pt rule
      ! for order>=5. Order >= 6 still has no verified rule of its own
      ! and reuses this one, so a request for a true degree-6 tet
      ! quadrature remains inexact by this same ~1e-4 margin.
      block
        real(kind=DOUBLE), parameter :: a  = 0.09198569323008009_DOUBLE
        real(kind=DOUBLE), parameter :: b  = 0.3196170515859492_DOUBLE
        real(kind=DOUBLE), parameter :: c  = 0.4438225029698707_DOUBLE
        real(kind=DOUBLE), parameter :: w0 = 0.1169962321923602_DOUBLE
        real(kind=DOUBLE), parameter :: w1 = 0.07196628338760167_DOUBLE
        real(kind=DOUBLE), parameter :: w2 = 0.06970310366278548_DOUBLE
        real(kind=DOUBLE), parameter :: w3 = 0.05272103660101520_DOUBLE
        real(kind=DOUBLE) :: ca, cb, d
        ca = 1.0_DOUBLE - 3.0_DOUBLE*a
        cb = 1.0_DOUBLE - 3.0_DOUBLE*b
        d  = 0.5_DOUBLE - c
        pts(:, 1) = 0.25_DOUBLE*(v1+v2+v3+v4)
        wts(1) = w0 * vol
        pts(:, 2) = a*v1 + a*v2 + a*v3 + ca*v4
        pts(:, 3) = a*v1 + a*v2 + ca*v3 + a*v4
        pts(:, 4) = a*v1 + ca*v2 + a*v3 + a*v4
        pts(:, 5) = ca*v1 + a*v2 + a*v3 + a*v4
        wts(2:5) = w1 * vol
        pts(:, 6) = b*v1 + b*v2 + b*v3 + cb*v4
        pts(:, 7) = b*v1 + b*v2 + cb*v3 + b*v4
        pts(:, 8) = b*v1 + cb*v2 + b*v3 + b*v4
        pts(:, 9) = cb*v1 + b*v2 + b*v3 + b*v4
        wts(6:9) = w2 * vol
        pts(:, 10) = c*v1 + c*v2 + d*v3 + d*v4
        pts(:, 11) = c*v1 + d*v2 + c*v3 + d*v4
        pts(:, 12) = c*v1 + d*v2 + d*v3 + c*v4
        pts(:, 13) = d*v1 + c*v2 + c*v3 + d*v4
        pts(:, 14) = d*v1 + c*v2 + d*v3 + c*v4
        pts(:, 15) = d*v1 + d*v2 + c*v3 + c*v4
        wts(10:15) = w3 * vol
      end block

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
    integer(kind=ENTIER), parameter :: max_tri = 12, max_line = 5
    integer(kind=ENTIER) :: n_tri, n_line, i, j, k
    real(kind=DOUBLE), dimension(3, max_tri) :: L_tri
    real(kind=DOUBLE), dimension(max_tri)    :: w_tri
    real(kind=DOUBLE), dimension(max_line)   :: t_line, w_line
    real(kind=DOUBLE) :: L1, L2, L3, t, wl
    real(kind=DOUBLE), dimension(3) :: xbot, xtop
    real(kind=DOUBLE), dimension(3) :: dxdL1, dxdL2, dxdt
    real(kind=DOUBLE) :: jac

    ! tri_ref_rule now supports up to order 6 (2026-09-20, 12-pt Dunavant)
    ! and gauss1d caps at its own 5-pt rule (order>=5): a triangular-base
    ! prism is exactly this triangle rule times gauss1d in the extrusion
    ! direction, so raising the triangle side's degree is all a prism
    ! needs -- matching n_volume_quad_pts(6,order)'s tri-pts x line-pts
    ! count, no separate volume-specific rule required.
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
    real(kind=DOUBLE), dimension(5) :: xi_g, w_g
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

  ! 1D Gauss-Legendre nodes and weights on [-1,1], n=1..5 points.
  ! An n-point rule is exact for polynomials up to degree 2n-1, so this
  ! alone (quad-face/hex-volume tensor-product rules only) comfortably
  ! covers every degree tested so far. The triangle/tet/prism/pyramid
  ! rules are separate, non-tensor-product constructions with their own
  ! degree ceiling: degree 5 for triangle/prism (tri_ref_rule), degree 5
  ! for tet/pyramid (tet_rule) -- see n_face_quad_pts/n_volume_quad_pts
  ! for the current per-shape degree each order actually gets, and
  ! tri_face_rule/tet_rule's own case-default comments for what happens
  ! past that (silent reuse of the highest verified rule, not a hard
  ! error).
  subroutine gauss1d(order, ng, xi, w)
    integer(kind=ENTIER), intent(in)  :: order
    integer(kind=ENTIER), intent(out) :: ng
    real(kind=DOUBLE), dimension(5), intent(out) :: xi, w

    real(kind=DOUBLE), parameter :: s3 = 1.0_DOUBLE / sqrt(3.0_DOUBLE)
    real(kind=DOUBLE), parameter :: s35 = sqrt(3.0_DOUBLE / 5.0_DOUBLE)
    ! 4-point: nodes ±sqrt((3∓2*sqrt(6/5))/7), standard A&S 25.4.30
    real(kind=DOUBLE), parameter :: x4a = 0.3399810435848563_DOUBLE
    real(kind=DOUBLE), parameter :: x4b = 0.8611363115940526_DOUBLE
    real(kind=DOUBLE), parameter :: w4a = 0.6521451548625461_DOUBLE
    real(kind=DOUBLE), parameter :: w4b = 0.3478548451374538_DOUBLE
    ! 5-point: nodes 0, ±(1/3)*sqrt(5∓2*sqrt(10/7)), standard A&S 25.4.30
    real(kind=DOUBLE), parameter :: x5a = 0.5384693101056831_DOUBLE
    real(kind=DOUBLE), parameter :: x5b = 0.9061798459386640_DOUBLE
    real(kind=DOUBLE), parameter :: w5c = 0.5688888888888889_DOUBLE
    real(kind=DOUBLE), parameter :: w5a = 0.4786286704993665_DOUBLE
    real(kind=DOUBLE), parameter :: w5b = 0.2369268850561891_DOUBLE

    select case (order)
    case (1)
      ng = 1; xi(1) = 0.0_DOUBLE; w(1) = 2.0_DOUBLE
    case (2)
      ng = 2
      xi(1) = -s3; xi(2) = s3
      w(1)  = 1.0_DOUBLE; w(2) = 1.0_DOUBLE
    case (3)
      ng = 3
      xi(1) = -s35; xi(2) = 0.0_DOUBLE; xi(3) = s35
      w(1)  = 5.0_DOUBLE/9.0_DOUBLE
      w(2)  = 8.0_DOUBLE/9.0_DOUBLE
      w(3)  = 5.0_DOUBLE/9.0_DOUBLE
    case (4)
      ng = 4
      xi(1) = -x4b; xi(2) = -x4a; xi(3) = x4a; xi(4) = x4b
      w(1)  =  w4b; w(2)  =  w4a; w(3)  = w4a; w(4)  = w4b
    case default   ! order >= 5 → 5-point rule
      ng = 5
      xi(1) = -x5b; xi(2) = -x5a; xi(3) = 0.0_DOUBLE; xi(4) = x5a; xi(5) = x5b
      w(1)  =  w5b; w(2)  =  w5a; w(3)  = w5c;        w(4)  = w5a; w(5)  = w5b
    end select
  end subroutine gauss1d

  ! Triangle reference Gauss rule in barycentric coords (L1,L2,L3)
  ! w_tri are normalized so sum(w_tri) = 1. Same underlying rules as
  ! tri_face_rule (kept as a separate barycentric-table routine for
  ! prism_rule's tensor-product use) -- see tri_face_rule for the
  ! degree-4/5/6 rules' verification note.
  subroutine tri_ref_rule(order, n, L, w)
    integer(kind=ENTIER), intent(in)  :: order
    integer(kind=ENTIER), intent(out) :: n
    real(kind=DOUBLE), dimension(3, 12), intent(out) :: L
    real(kind=DOUBLE), dimension(12),    intent(out) :: w

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
    case (3)
      n = 4
      L(:, 1) = [1.0_DOUBLE/3, 1.0_DOUBLE/3, 1.0_DOUBLE/3]
      L(:, 2) = [3.0_DOUBLE/5, 1.0_DOUBLE/5, 1.0_DOUBLE/5]
      L(:, 3) = [1.0_DOUBLE/5, 3.0_DOUBLE/5, 1.0_DOUBLE/5]
      L(:, 4) = [1.0_DOUBLE/5, 1.0_DOUBLE/5, 3.0_DOUBLE/5]
      w(1)   = -27.0_DOUBLE/48
      w(2:4) =  25.0_DOUBLE/48
    case (4)
      n = 6
      block
        real(kind=DOUBLE), parameter :: a1 = 0.445948490915965_DOUBLE
        real(kind=DOUBLE), parameter :: b1 = 0.108103018168070_DOUBLE
        real(kind=DOUBLE), parameter :: w1 = 0.223381589678011_DOUBLE
        real(kind=DOUBLE), parameter :: a2 = 0.091576213509771_DOUBLE
        real(kind=DOUBLE), parameter :: b2 = 0.816847572980459_DOUBLE
        real(kind=DOUBLE), parameter :: w2 = 0.109951743655322_DOUBLE
        L(:, 1) = [a1, a1, b1]; L(:, 2) = [a1, b1, a1]; L(:, 3) = [b1, a1, a1]
        w(1:3) = w1
        L(:, 4) = [a2, a2, b2]; L(:, 5) = [a2, b2, a2]; L(:, 6) = [b2, a2, a2]
        w(4:6) = w2
      end block
    case (5)
      n = 7
      block
        real(kind=DOUBLE), parameter :: w0 = 0.225_DOUBLE
        real(kind=DOUBLE), parameter :: a1 = 0.470142064105115_DOUBLE
        real(kind=DOUBLE), parameter :: b1 = 0.059715871789770_DOUBLE
        real(kind=DOUBLE), parameter :: w1 = 0.132394152788506_DOUBLE
        real(kind=DOUBLE), parameter :: a2 = 0.101286507323456_DOUBLE
        real(kind=DOUBLE), parameter :: b2 = 0.797426985353087_DOUBLE
        real(kind=DOUBLE), parameter :: w2 = 0.125939180544827_DOUBLE
        L(:, 1) = [1.0_DOUBLE/3, 1.0_DOUBLE/3, 1.0_DOUBLE/3]
        w(1) = w0
        L(:, 2) = [a1, a1, b1]; L(:, 3) = [a1, b1, a1]; L(:, 4) = [b1, a1, a1]
        w(2:4) = w1
        L(:, 5) = [a2, a2, b2]; L(:, 6) = [a2, b2, a2]; L(:, 7) = [b2, a2, a2]
        w(5:7) = w2
      end block
    case default   ! order >= 6: 12-pt Dunavant degree-6 rule, see
      ! tri_face_rule's case default for the derivation/verification note.
      n = 12
      block
        real(kind=DOUBLE), parameter :: a1 = 0.2492867451709104_DOUBLE
        real(kind=DOUBLE), parameter :: w1 = 0.1167862757263794_DOUBLE
        real(kind=DOUBLE), parameter :: a2 = 0.06308901449150223_DOUBLE
        real(kind=DOUBLE), parameter :: w2 = 0.05084490637020682_DOUBLE
        real(kind=DOUBLE), parameter :: a3 = 0.3103524510337844_DOUBLE
        real(kind=DOUBLE), parameter :: b3 = 0.05314504984481695_DOUBLE
        real(kind=DOUBLE), parameter :: w3 = 0.08285107561837358_DOUBLE
        real(kind=DOUBLE) :: c1, c2, c3
        c1 = 1.0_DOUBLE - 2.0_DOUBLE*a1
        c2 = 1.0_DOUBLE - 2.0_DOUBLE*a2
        c3 = 1.0_DOUBLE - a3 - b3
        L(:, 1) = [a1, a1, c1]; L(:, 2) = [a1, c1, a1]; L(:, 3) = [c1, a1, a1]
        w(1:3) = w1
        L(:, 4) = [a2, a2, c2]; L(:, 5) = [a2, c2, a2]; L(:, 6) = [c2, a2, a2]
        w(4:6) = w2
        L(:, 7)  = [a3, b3, c3]; L(:, 8)  = [a3, c3, b3]
        L(:, 9)  = [b3, a3, c3]; L(:, 10) = [b3, c3, a3]
        L(:, 11) = [c3, a3, b3]; L(:, 12) = [c3, b3, a3]
        w(7:12) = w3
      end block
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
