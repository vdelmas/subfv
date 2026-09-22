! Mesh-node movement that aligns a shock with mesh lines, factored out of
! shock_adapt_module so a second solver (euler_ho) can drive the same
! algorithm. Everything here is solver-agnostic: the sensor takes a cell
! density array rather than a conservative-state layout, and the thresholds
! are arguments rather than namelist globals, because the two callers keep
! their state in different shapes (5 components here, 6 with gamma there).
module shock_adapt_move_module
  use precision_module
  use mesh_module
  implicit none

  private
  public :: compute_shock_sensor_grad_rho
  public :: build_vert_adjacency
  public :: compute_node_displacement_curvature
  public :: compute_local_scale
  public :: move_mesh
  public :: min_elem_volume
  public :: mpi_memory_exchange_vert

contains

  subroutine compute_shock_sensor_grad_rho(mesh, rho, threshold, node_sensor, node_flagged)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: rho
    real(kind=DOUBLE), intent(in) :: threshold
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(out) :: node_sensor
    logical, dimension(mesh%n_vert), intent(out) :: node_flagged

    integer(kind=ENTIER) :: iv, j, ide
    real(kind=DOUBLE), dimension(2) :: dx, rhs2
    real(kind=DOUBLE), dimension(2, 2) :: mat2
    real(kind=DOUBLE) :: wj, sum_w, rho_bar, gx, gy, local_scale, grad_mag

    node_sensor = 0.0_DOUBLE
    node_flagged = .FALSE.

    do iv = 1, mesh%n_vert
      if (mesh%vert(iv)%is_bound) cycle ! never flag/move boundary-geometry nodes
      if (mesh%vert(iv)%n_elems_neigh < 3) cycle ! need >=3 for a well-posed 2x2 WLS fit

      sum_w = 0.0_DOUBLE; rho_bar = 0.0_DOUBLE
      do j = 1, mesh%vert(iv)%n_elems_neigh
        ide = mesh%vert(iv)%elem_neigh(j)
        dx = mesh%elem(ide)%coord(1:2) - mesh%vert(iv)%coord(1:2)
        wj = 1.0_DOUBLE/max(dot_product(dx, dx), 1e-30_DOUBLE)
        sum_w = sum_w + wj
        rho_bar = rho_bar + wj*rho(ide)
      end do
      rho_bar = rho_bar/sum_w

      mat2 = 0.0_DOUBLE; rhs2 = 0.0_DOUBLE
      do j = 1, mesh%vert(iv)%n_elems_neigh
        ide = mesh%vert(iv)%elem_neigh(j)
        dx = mesh%elem(ide)%coord(1:2) - mesh%vert(iv)%coord(1:2)
        wj = 1.0_DOUBLE/max(dot_product(dx, dx), 1e-30_DOUBLE)
        mat2(1, 1) = mat2(1, 1) + wj*dx(1)*dx(1)
        mat2(1, 2) = mat2(1, 2) + wj*dx(1)*dx(2)
        mat2(2, 2) = mat2(2, 2) + wj*dx(2)*dx(2)
        rhs2(1) = rhs2(1) + wj*dx(1)*(rho(ide) - rho_bar)
        rhs2(2) = rhs2(2) + wj*dx(2)*(rho(ide) - rho_bar)
      end do
      mat2(2, 1) = mat2(1, 2)

      call solve2(mat2, rhs2, gx, gy)
      grad_mag = sqrt(gx*gx + gy*gy)

      local_scale = huge(1.0_DOUBLE)
      do j = 1, mesh%vert(iv)%n_elems_neigh
        ide = mesh%vert(iv)%elem_neigh(j)
        local_scale = min(local_scale, norm2(mesh%vert(iv)%coord - mesh%elem(ide)%coord))
      end do

      node_sensor(iv) = grad_mag*local_scale/rho_bar
      node_flagged(iv) = node_sensor(iv) > threshold
    end do
  end subroutine compute_shock_sensor_grad_rho

  ! Topological vertex-vertex adjacency (in-plane mesh edges only). Built
  ! from the z-normal "cap" faces (abs(norm(3))>1-1e-8) -- under the
  ! boundary_2d 1-cell extrusion convention, BOTH the top and bottom caps
  ! are themselves domain-boundary faces (there is only one layer), and
  ! each cap face's vertex list is exactly the true 2D triangle/quad
  ! connectivity: consecutive vertices in a cap face's vert() list are
  ! genuine in-plane mesh edges, with no cross-z-layer contamination (top
  ! and bottom layers get their own, separate adjacency, as they should
  ! since each layer moves independently -- though identically, by mesh
  ! symmetry). Fixed max degree per vertex (16, generous for a 2D tri/quad
  ! mesh); errors loudly if exceeded rather than silently truncating.
  subroutine build_vert_adjacency(mesh, n_neigh, vneigh)
    implicit none
    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), dimension(mesh%n_vert), intent(out) :: n_neigh
    integer(kind=ENTIER), dimension(16, mesh%n_vert), intent(out) :: vneigh

    integer(kind=ENTIER) :: iface, k, n, a, b, kk
    logical :: already

    n_neigh = 0

    do iface = 1, mesh%n_faces
      if (abs(mesh%face(iface)%norm(3)) < 1.0_DOUBLE - 1e-8_DOUBLE) cycle ! not a cap face
      n = mesh%face(iface)%n_vert
      do k = 1, n
        a = mesh%face(iface)%vert(k)
        b = mesh%face(iface)%vert(mod(k, n) + 1)

        already = .FALSE.
        do kk = 1, n_neigh(a)
          if (vneigh(kk, a) == b) already = .TRUE.
        end do
        if (.not. already) then
          n_neigh(a) = n_neigh(a) + 1
          if (n_neigh(a) > 16) then
            print*, "[-] build_vert_adjacency: vertex degree exceeds 16 at vert", a
            error stop
          end if
          vneigh(n_neigh(a), a) = b
        end if

        already = .FALSE.
        do kk = 1, n_neigh(b)
          if (vneigh(kk, b) == a) already = .TRUE.
        end do
        if (.not. already) then
          n_neigh(b) = n_neigh(b) + 1
          if (n_neigh(b) > 16) then
            print*, "[-] build_vert_adjacency: vertex degree exceeds 16 at vert", b
            error stop
          end if
          vneigh(n_neigh(b), b) = a
        end if
      end do
    end do
  end subroutine build_vert_adjacency

  ! Move each flagged node along the local normal of a quadratic curve
  ! fit through its flagged topological neighbors, to reduce the front's
  ! local (discrete) curvature -- a node is predicted from its neighbors'
  ! own quadratic trend, and moved toward that prediction (never using its
  ! own position in the fit, so this is a genuine "flatten toward
  ! neighbors" step, not self-reinforcing).
  !
  ! For node i with flagged neighbors {j}: local 2D frame at the neighbor
  ! centroid, tangent = principal axis of the neighbor point cloud (2x2
  ! eigenproblem, closed-form), normal = perpendicular. Project neighbors
  ! (and i) into (s,h) local coords. Fit h(s)=a*s^2+b*s+c through the
  ! neighbors only (least squares if >3 points, exact if 3, linear-only
  ! i.e. a=0 if only 2). Predicted height at i's own s: h_pred(s_i). Move
  ! i along the normal by (h_pred(s_i) - h_i), capped as usual.
  subroutine compute_node_displacement_curvature(mesh, node_flagged, n_neigh, vneigh, &
      max_move_frac, curvature_relax, disp, n_moved, max_disp, &
      orig_local_scale, cum_disp)
    implicit none
    type(mesh_type), intent(in) :: mesh
    logical, dimension(mesh%n_vert), intent(in) :: node_flagged
    integer(kind=ENTIER), dimension(mesh%n_vert), intent(in) :: n_neigh
    integer(kind=ENTIER), dimension(16, mesh%n_vert), intent(in) :: vneigh
    real(kind=DOUBLE), intent(in) :: max_move_frac, curvature_relax
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(out) :: disp
    integer(kind=ENTIER), intent(out) :: n_moved
    real(kind=DOUBLE), intent(out) :: max_disp
    ! Supply both to cap the cumulative movement across repeated cycles
    ! against the pre-movement geometry; omit both for a single-shot move.
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in), optional :: orig_local_scale
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout), optional :: cum_disp

    integer(kind=ENTIER) :: iv, j, jn, jn2, k, kk
    integer(kind=ENTIER), dimension(16) :: nbr_id
    real(kind=DOUBLE), dimension(2, 16) :: nbr_xy
    real(kind=DOUBLE), dimension(2) :: centroid, t_hat, n_hat, d, xy_i
    real(kind=DOUBLE), dimension(16) :: s, h
    real(kind=DOUBLE) :: cxx, cxy, cyy, tr, det_, lam1, ex, ey, enorm
    real(kind=DOUBLE) :: s_i, h_i, h_pred, a_q, b_q, c_q, dn, local_scale, cap, mag
    real(kind=DOUBLE), dimension(3, 3) :: mat
    real(kind=DOUBLE), dimension(3) :: rhs_v, new_cum
    real(kind=DOUBLE) :: cum_mag_before
    logical :: cumulative

    disp = 0.0_DOUBLE
    n_moved = 0
    max_disp = 0.0_DOUBLE
    cumulative = present(orig_local_scale) .and. present(cum_disp)

    do iv = 1, mesh%n_vert
      if (.not. node_flagged(iv)) cycle

      ! Collect flagged neighbors up to 2 hops along the flagged-node
      ! subgraph (not just direct topological neighbors): a clean shock
      ! front is topologically a thin CHAIN, so an interior node typically
      ! has exactly 2 direct flagged neighbors (one on each side) -- never
      ! 3 -- so requiring k>=3 from direct neighbors alone would leave
      ! every chain-interior node permanently unmovable (confirmed: an
      ! early version flagged 40 nodes but moved 0 of them). Extending to
      ! each direct flagged neighbor's OWN flagged neighbors too (still
      ! excluding iv itself) gives an interior chain node up to 4 points
      ! (2 on each side), enough for a meaningful quadratic/curvature fit.
      k = 0
      nbr_id = 0
      do j = 1, n_neigh(iv)
        jn = vneigh(j, iv)
        if (.not. node_flagged(jn) .or. k >= 16) cycle
        if (any(nbr_id(1:k) == jn)) cycle
        k = k + 1
        nbr_xy(:, k) = mesh%vert(jn)%coord(1:2)
        nbr_id(k) = jn

        do kk = 1, n_neigh(jn)
          jn2 = vneigh(kk, jn)
          if (jn2 == iv .or. .not. node_flagged(jn2) .or. k >= 16) cycle
          if (any(nbr_id(1:k) == jn2)) cycle
          k = k + 1
          nbr_xy(:, k) = mesh%vert(jn2)%coord(1:2)
          nbr_id(k) = jn2
        end do
      end do
      ! Need >=3 neighbors for a curvature-bearing (quadratic) fit -- a
      ! 2-point line has zero curvature everywhere and gives no reliable
      ! direction to move toward (a plane/line is exactly the degenerate
      ! case that motivated requiring a curvature term in the first place).
      if (k < 3) cycle

      centroid = 0.0_DOUBLE
      do kk = 1, k
        centroid = centroid + nbr_xy(:, kk)
      end do
      centroid = centroid/real(k, DOUBLE)

      ! Principal direction of the neighbor point cloud: 2x2 covariance,
      ! closed-form dominant eigenvector.
      cxx = 0.0_DOUBLE; cxy = 0.0_DOUBLE; cyy = 0.0_DOUBLE
      do kk = 1, k
        d = nbr_xy(:, kk) - centroid
        cxx = cxx + d(1)*d(1); cxy = cxy + d(1)*d(2); cyy = cyy + d(2)*d(2)
      end do
      tr = cxx + cyy
      det_ = cxx*cyy - cxy*cxy
      lam1 = 0.5_DOUBLE*tr + sqrt(max(0.0_DOUBLE, 0.25_DOUBLE*tr*tr - det_))
      if (abs(cxy) > 1e-14_DOUBLE) then
        ex = lam1 - cyy; ey = cxy
      else if (cxx >= cyy) then
        ex = 1.0_DOUBLE; ey = 0.0_DOUBLE
      else
        ex = 0.0_DOUBLE; ey = 1.0_DOUBLE
      end if
      enorm = sqrt(ex*ex + ey*ey)
      if (enorm < 1e-14_DOUBLE) cycle ! degenerate (coincident neighbors)
      t_hat = (/ex, ey/)/enorm
      n_hat = (/-t_hat(2), t_hat(1)/)

      do kk = 1, k
        d = nbr_xy(:, kk) - centroid
        s(kk) = dot_product(d, t_hat)
        h(kk) = dot_product(d, n_hat)
      end do
      xy_i = mesh%vert(iv)%coord(1:2) - centroid
      s_i = dot_product(xy_i, t_hat)
      h_i = dot_product(xy_i, n_hat)

      ! Least-squares (or exact if k==3) quadratic fit through the
      ! neighbors: normal equations for [a,b,c].
      mat = 0.0_DOUBLE; rhs_v = 0.0_DOUBLE
      do kk = 1, k
        mat(1, 1) = mat(1, 1) + s(kk)**4
        mat(1, 2) = mat(1, 2) + s(kk)**3
        mat(1, 3) = mat(1, 3) + s(kk)**2
        mat(2, 3) = mat(2, 3) + s(kk)
        mat(3, 3) = mat(3, 3) + 1.0_DOUBLE
        rhs_v(1) = rhs_v(1) + s(kk)**2*h(kk)
        rhs_v(2) = rhs_v(2) + s(kk)*h(kk)
        rhs_v(3) = rhs_v(3) + h(kk)
      end do
      mat(2, 1) = mat(1, 2); mat(2, 2) = mat(1, 3)
      mat(3, 1) = mat(1, 3); mat(3, 2) = mat(2, 3)
      call solve3(mat, rhs_v, a_q, b_q, c_q)

      h_pred = a_q*s_i*s_i + b_q*s_i + c_q
      dn = curvature_relax*(h_pred - h_i)

      if (cumulative) then
        ! Cap the CUMULATIVE total against the ORIGINAL (pre-movement) local
        ! scale. Capping each cycle's step against the CURRENT geometry (the
        ! branch below) bounds one step but not the drift: the scale shrinks
        ! as the cells compress, so repeated cycles keep squeezing the same
        ! cells. That is what collapsed dt by 300x over 4 cycles here.
        cap = max_move_frac*orig_local_scale(iv)
        if (abs(dn) > cap) dn = sign(cap, dn)
        disp(1:2, iv) = dn*n_hat
        disp(3, iv) = 0.0_DOUBLE

        cum_mag_before = norm2(cum_disp(:, iv))
        if (cum_mag_before >= cap) then
          disp(:, iv) = 0.0_DOUBLE
        else
          new_cum = cum_disp(:, iv) + disp(:, iv)
          if (norm2(new_cum) > cap) then
            ! Scale the increment so the new cumulative total lands exactly on
            ! the cap, rather than scaling the increment on its own.
            disp(:, iv) = disp(:, iv)*max(0.0_DOUBLE, &
              (cap - cum_mag_before)/max(norm2(disp(:, iv)), 1e-300_DOUBLE))
          end if
        end if
        cum_disp(:, iv) = cum_disp(:, iv) + disp(:, iv)
        dn = norm2(disp(:, iv))
      else
        ! Cap by max_move_frac of the local neighbor-centroid distance
        ! (same safety spirit as the snap-based mover).
        local_scale = huge(1.0_DOUBLE)
        do j = 1, mesh%vert(iv)%n_elems_neigh
          local_scale = min(local_scale, &
            norm2(mesh%vert(iv)%coord - mesh%elem(mesh%vert(iv)%elem_neigh(j))%coord))
        end do
        cap = max_move_frac*local_scale
        if (abs(dn) > cap) dn = sign(cap, dn)

        disp(1:2, iv) = dn*n_hat
        disp(3, iv) = 0.0_DOUBLE
      end if

      mag = abs(dn)
      if (mag > 1e-14_DOUBLE) then
        n_moved = n_moved + 1
        max_disp = max(max_disp, mag)
      end if
    end do
  end subroutine compute_node_displacement_curvature

  ! Tiny local linear solves (3x3/2x2), no LAPACK dependency needed for
  ! per-vertex systems this small.
  subroutine solve3(mat, rhs, x1, x2, x3)
    implicit none
    real(kind=DOUBLE), dimension(3, 3), intent(in) :: mat
    real(kind=DOUBLE), dimension(3), intent(in) :: rhs
    real(kind=DOUBLE), intent(out) :: x1, x2, x3
    real(kind=DOUBLE) :: det_
    real(kind=DOUBLE), dimension(3, 3) :: m1, m2, m3

    det_ = mat(1, 1)*(mat(2, 2)*mat(3, 3) - mat(2, 3)*mat(3, 2)) &
      - mat(1, 2)*(mat(2, 1)*mat(3, 3) - mat(2, 3)*mat(3, 1)) &
      + mat(1, 3)*(mat(2, 1)*mat(3, 2) - mat(2, 2)*mat(3, 1))

    if (abs(det_) < 1e-30_DOUBLE) then
      x1 = 0.0_DOUBLE; x2 = 0.0_DOUBLE; x3 = 0.0_DOUBLE
      return
    end if

    m1 = mat; m1(:, 1) = rhs
    m2 = mat; m2(:, 2) = rhs
    m3 = mat; m3(:, 3) = rhs

    x1 = (m1(1, 1)*(m1(2, 2)*m1(3, 3) - m1(2, 3)*m1(3, 2)) &
      - m1(1, 2)*(m1(2, 1)*m1(3, 3) - m1(2, 3)*m1(3, 1)) &
      + m1(1, 3)*(m1(2, 1)*m1(3, 2) - m1(2, 2)*m1(3, 1)))/det_
    x2 = (m2(1, 1)*(m2(2, 2)*m2(3, 3) - m2(2, 3)*m2(3, 2)) &
      - m2(1, 2)*(m2(2, 1)*m2(3, 3) - m2(2, 3)*m2(3, 1)) &
      + m2(1, 3)*(m2(2, 1)*m2(3, 2) - m2(2, 2)*m2(3, 1)))/det_
    x3 = (m3(1, 1)*(m3(2, 2)*m3(3, 3) - m3(2, 3)*m3(3, 2)) &
      - m3(1, 2)*(m3(2, 1)*m3(3, 3) - m3(2, 3)*m3(3, 1)) &
      + m3(1, 3)*(m3(2, 1)*m3(3, 2) - m3(2, 2)*m3(3, 1)))/det_
  end subroutine solve3

  subroutine solve2(mat, rhs, x1, x2)
    implicit none
    real(kind=DOUBLE), dimension(2, 2), intent(in) :: mat
    real(kind=DOUBLE), dimension(2), intent(in) :: rhs
    real(kind=DOUBLE), intent(out) :: x1, x2
    real(kind=DOUBLE) :: det_

    det_ = mat(1, 1)*mat(2, 2) - mat(1, 2)*mat(2, 1)
    if (abs(det_) < 1e-30_DOUBLE) then
      x1 = 0.0_DOUBLE; x2 = 0.0_DOUBLE
      return
    end if
    x1 = (rhs(1)*mat(2, 2) - mat(1, 2)*rhs(2))/det_
    x2 = (mat(1, 1)*rhs(2) - rhs(1)*mat(2, 1))/det_
  end subroutine solve2

  ! centroid), computed once on the given mesh state -- call this right
  ! after Stage 1 (before any movement) to get `orig_local_scale` for the
  ! iterative mode's cumulative displacement cap.
  subroutine compute_local_scale(mesh, local_scale)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(out) :: local_scale

    integer(kind=ENTIER) :: iv, j

    do iv = 1, mesh%n_vert
      local_scale(iv) = huge(1.0_DOUBLE)
      do j = 1, mesh%vert(iv)%n_elems_neigh
        local_scale(iv) = min(local_scale(iv), &
          norm2(mesh%vert(iv)%coord - mesh%elem(mesh%vert(iv)%elem_neigh(j))%coord))
      end do
    end do
  end subroutine compute_local_scale

  subroutine move_mesh(mesh, disp)
    implicit none
    type(mesh_type), intent(inout) :: mesh
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: disp

    integer(kind=ENTIER) :: i

    do i = 1, mesh%n_vert
      mesh%vert(i)%coord = mesh%vert(i)%coord + disp(:, i)
    end do
  end subroutine move_mesh

  function min_elem_volume(mesh) result(vmin)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE) :: vmin

    integer(kind=ENTIER) :: i

    vmin = huge(1.0_DOUBLE)
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) vmin = min(vmin, mesh%elem(i)%volume)
    end do
  end function min_elem_volume


  ! Ported verbatim from lagrange_module.F90 / ale_module.F90 (identical
  ! routine in both): send/recv vert%coord for every vertex belonging to a
  ! shared/ghost element, so a moved mesh stays consistent across MPI
  ! partition boundaries. Call after every move_mesh, before
  ! compute_geometry_mesh, whenever num_procs>1.
  subroutine mpi_memory_exchange_vert(mesh, mpi_send_recv)
    use mpi
    use mpi_module
    implicit none

    type(mesh_type), intent(inout) :: mesh
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    integer :: mpi_ierr
    integer(kind=ENTIER) :: i, k, j
    integer(kind=ENTIER) :: id_elem, n_vert_tot, id_vert
    logical, dimension(:), allocatable :: vert_averaged

    do i = 1, mpi_send_recv%n_mpi_send_neigh
      n_vert_tot = 0
      do k = 1, mpi_send_recv%mpi_send_neigh(i)%n_elems
        id_elem = mpi_send_recv%mpi_send_neigh(i)%elem_id(k)
        n_vert_tot = n_vert_tot + mesh%elem(id_elem)%n_vert
      end do

      allocate(mpi_send_recv%mpi_send_neigh(i)%sol(3, n_vert_tot))

      n_vert_tot = 1
      do k = 1, mpi_send_recv%mpi_send_neigh(i)%n_elems
        id_elem = mpi_send_recv%mpi_send_neigh(i)%elem_id(k)
        do j = 1, mesh%elem(id_elem)%n_vert
          id_vert = mesh%elem(id_elem)%vert(j)
          mpi_send_recv%mpi_send_neigh(i)%sol(:, n_vert_tot) = mesh%vert(id_vert)%coord
          n_vert_tot = n_vert_tot + 1
        end do
      end do

      call mpi_isend(mpi_send_recv%mpi_send_neigh(i)%sol(1, 1), &
        3*(n_vert_tot - 1), MPI_DOUBLE, &
        mpi_send_recv%mpi_send_neigh(i)%partition_id, &
        0, MPI_COMM_WORLD, mpi_send_recv%mpi_reqsend(i), mpi_ierr)
    end do

    do i = 1, mpi_send_recv%n_mpi_recv_neigh
      n_vert_tot = 0
      do k = 1, mpi_send_recv%mpi_recv_neigh(i)%n_elems
        id_elem = mpi_send_recv%mpi_recv_neigh(i)%elem_id(k)
        n_vert_tot = n_vert_tot + mesh%elem(id_elem)%n_vert
      end do
      allocate(mpi_send_recv%mpi_recv_neigh(i)%sol(3, n_vert_tot))
      call mpi_irecv(mpi_send_recv%mpi_recv_neigh(i)%sol(1, 1), &
        3*n_vert_tot, MPI_DOUBLE, &
        mpi_send_recv%mpi_recv_neigh(i)%partition_id, &
        MPI_ANY_TAG, MPI_COMM_WORLD, mpi_send_recv%mpi_reqrecv(i), mpi_ierr)
    end do

    call mpi_waitall(mpi_send_recv%n_mpi_send_neigh, &
      mpi_send_recv%mpi_reqsend, mpi_send_recv%mpi_sendstat, mpi_ierr)
    call mpi_waitall(mpi_send_recv%n_mpi_recv_neigh, &
      mpi_send_recv%mpi_reqrecv, mpi_send_recv%mpi_recvstat, mpi_ierr)

    ! Unpack: a vertex genuinely owned only as a ghost here (is_ghost=.true.
    ! -- ALL its neighboring elements are ghost elements, see
    ! mesh_connectivity_module.F90's compute_ghost_vert) has no legitimate
    ! local computation to protect, so take the sender's value outright.
    ! A vertex on the actual partition CUT, though, touches at least one
    ! local (owned) element on BOTH sides of the boundary, so
    ! is_ghost=.FALSE. there on every rank that touches it -- each such
    ! rank independently computed its own (generally different, since each
    ! only sees its own side's local neighbor subgraph for the curvature
    ! fit) displacement for the SAME physical vertex. Confirmed by the
    ! user directly inspecting a cluster run: shared boundary nodes ended
    ! up at different positions on each side. Average with the received
    ! value there instead of overwriting, so both ranks converge on the
    ! same final position (order-independent, symmetric).
    !
    ! A vertex on the cut is typically referenced by SEVERAL ghost elements
    ! from the same neighbor rank (and possibly, at a multi-rank corner, by
    ! ghost elements from more than one neighbor), so it appears more than
    ! once across these nested loops. The averaging update is NOT
    ! idempotent (applying it twice in place biases the result toward the
    ! received value: 0.5, then 0.75, then 0.875, ...), and the occurrence
    ! count generally differs between the two ranks sharing a vertex, so
    ! repeating it silently broke the intended 50/50 symmetry -- this was
    ! the actual cause of a visible mesh crack at partition cuts after a
    ! second movement cycle. Guard so each vertex is only averaged once per
    ! call, regardless of how many times it is encountered.
    allocate(vert_averaged(size(mesh%vert)))
    vert_averaged = .false.
    do i = 1, mpi_send_recv%n_mpi_recv_neigh
      n_vert_tot = 1
      do k = 1, mpi_send_recv%mpi_recv_neigh(i)%n_elems
        id_elem = mpi_send_recv%mpi_recv_neigh(i)%elem_id(k)
        do j = 1, mesh%elem(id_elem)%n_vert
          id_vert = mesh%elem(id_elem)%vert(j)
          if (.not. vert_averaged(id_vert)) then
            if (mesh%vert(id_vert)%is_ghost) then
              mesh%vert(id_vert)%coord = mpi_send_recv%mpi_recv_neigh(i)%sol(:, n_vert_tot)
            else
              mesh%vert(id_vert)%coord = 0.5_DOUBLE*(mesh%vert(id_vert)%coord &
                + mpi_send_recv%mpi_recv_neigh(i)%sol(:, n_vert_tot))
            end if
            vert_averaged(id_vert) = .true.
          end if
          n_vert_tot = n_vert_tot + 1
        end do
      end do
    end do
    deallocate(vert_averaged)

    do i = 1, mpi_send_recv%n_mpi_send_neigh
      deallocate(mpi_send_recv%mpi_send_neigh(i)%sol)
    end do
    do i = 1, mpi_send_recv%n_mpi_recv_neigh
      deallocate(mpi_send_recv%mpi_recv_neigh(i)%sol)
    end do
  end subroutine mpi_memory_exchange_vert

end module shock_adapt_move_module
