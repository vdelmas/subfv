module silvia_reconstruction_module
	! use silvia_base_module
	use precision_module
  use mesh_module
  
	implicit none

	public :: reconstruct, ls_reconstruction

	contains

  ! Polynomial reconstruction of primitive variable w at point xq from cell i.
  ! Hierarchical fallback: if an order-p reconstruction gives unphysical
  ! rho or p, it is replaced by the order-(p-1) result.
  function reconstruct(prim, grad, hess, i, xq, xc) result(w)
    real(kind=DOUBLE), dimension(:, :),          intent(in) :: prim    ! (5, n_elems)
    real(kind=DOUBLE), dimension(:, :, :),       intent(in) :: grad    ! (5, 3, n_elems)
    real(kind=DOUBLE), dimension(:, :, :, :), allocatable, intent(in) :: hess
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE), dimension(3), intent(in) :: xq, xc
    real(kind=DOUBLE), dimension(5) :: w

    real(kind=DOUBLE), dimension(3) :: dx
    real(kind=DOUBLE), dimension(5) :: w_try
    integer(kind=ENTIER) :: j, k

    dx = xq - xc
    w  = prim(:, i)

    if (order >= 2) then
      w_try = w + matmul(grad(:, :, i), dx)
      if (physical_state(w_try, w)) w = w_try
    end if

    if (order >= 3 .and. allocated(hess)) then
      if (.not. allocated(cell_moment)) then
        print *, 'FATAL: reconstruct() called at order>=3 but compute_cell_moments ', &
          'was never called -- the order-3 Taylor polynomial would silently carry the ', &
          'wrong cell average (see header comment above compute_cell_moments). Call ', &
          'compute_cell_moments(mesh) once, right after compute_geometry_mesh.'
        error stop 1
      end if
      ! Always build from scratch: cell value + gradient + Hessian, MINUS
      ! the cell's own curvature-average bias (see compute_cell_moments)
      ! so that this polynomial's cell average is prim(:,i), not its
      ! value at xc plus a spurious O(h^2) offset.
      w_try = prim(:, i) + matmul(grad(:, :, i), dx)
      do j = 1, 3
        do k = 1, 3
          w_try = w_try + 0.5_DOUBLE * hess(:, j, k, i) &
            * (dx(j) * dx(k) - cell_moment(j, k, i))
        end do
      end do
      if (physical_state(w_try, prim(:, i))) w = w_try
    end if

    ! Enforce positivity (safety net)
    w(1) = max(w(1), 1.0e-12_DOUBLE)
    w(5) = max(w(5), 1.0e-12_DOUBLE)
  end function reconstruct

	! Returns .true. if w_cand is a physically acceptable reconstruction
  ! relative to the reference state w_ref.
  ! Rejects: negative rho or p; velocity more than 20x the reference speed + c.
  pure function physical_state(w_cand, w_ref) result(ok)
    real(kind=DOUBLE), dimension(5), intent(in) :: w_cand, w_ref
    logical :: ok
    real(kind=DOUBLE) :: spd_ref, spd_cand, c_ref

    ok = .false.
    if (w_cand(1) <= 0.0_DOUBLE) return
    if (w_cand(5) <= 0.0_DOUBLE) return

    ! Velocity magnitude check
    spd_ref  = w_ref(2)**2  + w_ref(3)**2  + w_ref(4)**2
    spd_cand = w_cand(2)**2 + w_cand(3)**2 + w_cand(4)**2
    c_ref    = gamma_gas * w_ref(5) / max(w_ref(1), 1.0e-16_DOUBLE)
    if (spd_cand > 400.0_DOUBLE * (spd_ref + c_ref)) return   ! 20x speed limit

    ok = .true.
  end function physical_state

	! ================================================================
  ! Internal subroutines
  ! ================================================================

  ! Cell-centred least-squares polynomial reconstruction.
  !
  ! For each non-ghost cell i, collects all vertex-neighbours (cells sharing
  ! at least one vertex), then fits a degree-(order-1) polynomial to their
  ! mean values by solving the normal equations (A^T A) x = A^T b.
  !
  ! Polynomial basis (boundary_2d=T, i.e. 2D):
  !   order 2: [dx, dy]                                    (2 coefficients)
  !   order 3: [dx, dy, dx²/2, dx·dy, dy²/2]              (5 coefficients)
  ! Full 3D:
  !   order 2: [dx, dy, dz]                                (3 coefficients)
  !   order 3: [dx, dy, dz, dx²/2, dx·dy, dx·dz,
  !              dy²/2, dy·dz, dz²/2]                      (9 coefficients)
  !
  ! Coefficients are stored in the existing grad/hess arrays so that the
  ! reconstruct() function needs no changes.
  subroutine ls_reconstruction(mesh, prim, grad, hess)
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: prim
    real(kind=DOUBLE), dimension(5, 3, mesh%n_elems), intent(inout) :: grad
    real(kind=DOUBLE), dimension(:, :, :, :), allocatable, intent(inout) :: hess

    integer(kind=ENTIER) :: i, iv, mf, iface, il, ir, jn, kv, nv
    integer(kind=ENTIER) :: n_coeff, n_neigh, n_needed
    integer(kind=ENTIER) :: ring_start, ring_end, ic, ivv
    integer(kind=ENTIER), dimension(500) :: neigh_list
    real(kind=DOUBLE), allocatable :: Amat(:, :)
    real(kind=DOUBLE), dimension(500) :: wt_arr   ! sqrt(1/dist) per neighbour
    real(kind=DOUBLE), dimension(9, 9) :: ATA
    real(kind=DOUBLE), dimension(9)    :: ATb, x_coeff
    real(kind=DOUBLE) :: ddx, ddy, ddz, dist, wt, Lref, sdx, sdy, sdz
    real(kind=DOUBLE), dimension(3) :: xci, xcj

    if (boundary_2d) then
      n_coeff = merge(5, 2, order >= 3)
    else
      n_coeff = merge(9, 3, order >= 3)
    end if

    do i = 1, mesh%n_elems
      if (mesh%elem(i)%is_ghost) cycle

      xci    = mesh%elem(i)%coord
      n_neigh = 0

      ! Collect unique real vertex-neighbours of cell i (ring 1)
      do nv = 1, mesh%elem(i)%n_vert
        iv = mesh%elem(i)%vert(nv)
        do mf = 1, mesh%vert(iv)%n_faces_neigh
          iface = mesh%vert(iv)%face_neigh(mf)
          call add_unique_neighbor(mesh%face(iface)%left_neigh)
          call add_unique_neighbor(mesh%face(iface)%right_neigh)
        end do
      end do

      ! Ring-expand (neighbours of neighbours, etc.) when the immediate
      ! vertex-neighbour ring doesn't have enough cells to determine the
      ! unknowns actually alive on THIS mesh -- e.g. a nominal order-3 fit
      ! has 5 unknowns (boundary_2d: dx,dy,dxx,dxy,dyy), but on a
      ! single-row-in-y mesh (subfv's Sod tube convention) every neighbour
      ! shares the same y, so dy/dxy/dyy are identically zero columns and
      ! only 2 unknowns (dx,dxx) are really determinable -- n_needed
      ! (below) tracks that reduced count, recomputed as the ring grows,
      ! so the stencil stays as narrow (and thus as locally accurate) as
      ! the data actually require instead of always demanding the full,
      ! nominal n_coeff neighbours even when most of those unknowns are
      ! degenerate. Previously this used a fixed n_coeff target and simply
      ! gave up (`cycle`, leaving grad=hess at their initialized 0) when
      ! the immediate ring fell short of it, silently degrading every
      ! order-3 ls_reconstruction cell to order 1 on the whole Sod-tube
      ! mesh family. Mirrors arbitrary_high_order_module's
      ! gather_ls_neighbors ring growth; the degeneracy detection mirrors
      ! its compute_nodal_derivative_at_vertex direction-spread check.
      ring_start = 1
      do
        n_needed = count_needed_coeffs()
        if (n_neigh >= n_needed .or. n_neigh >= 500) exit
        ring_end = n_neigh
        do ic = ring_start, ring_end
          do ivv = 1, mesh%elem(neigh_list(ic))%n_vert
            iv = mesh%elem(neigh_list(ic))%vert(ivv)
            do mf = 1, mesh%vert(iv)%n_faces_neigh
              iface = mesh%vert(iv)%face_neigh(mf)
              call add_unique_neighbor(mesh%face(iface)%left_neigh)
              call add_unique_neighbor(mesh%face(iface)%right_neigh)
            end do
          end do
        end do
        if (ring_end == n_neigh) exit   ! stencil can't grow further (tiny/disconnected mesh)
        ring_start = ring_end + 1
      end do

      if (n_neigh < n_needed) cycle   ! still under-determined after ring expansion: leave grad=0

      allocate(Amat(n_neigh, n_coeff))

      ! Characteristic length for THIS cell, used to nondimensionalize the
      ! design matrix columns before solving (see Lref note below) --
      ! cube root of the cell volume is a reasonable, cheap proxy for its
      ! own size regardless of cell shape.
      Lref = max(mesh%elem(i)%volume, 1.0e-300_DOUBLE)**(1.0_DOUBLE/3.0_DOUBLE)

      ! Build the design matrix A weighted by w_j = 1/dist_j (inverse-distance).
      ! Row j of Amat is scaled by sqrt(w_j) = 1/sqrt(dist_j) so that
      ! ATA = A^T W A  and  ATb = A^T W b  (with W = diag(w_j)). Columns are
      ! built from ddx/Lref etc (order-1 in the SCALED offset), not raw
      ! ddx -- found necessary 2026-09-14: with raw offsets, the linear-term
      ! columns (~dx*wt) and quadratic-term columns (~dx^2*wt) differ in
      ! scale by a factor of ~1/dx, which at a fine mesh (dx=0.01) inflates
      ! ATA's condition number enough to genuinely corrupt the solve (order-3
      ! ls_reconstruction was found NOT exact -- Linf 5e-4 instead of machine
      ! precision -- on a degree-2 test polynomial on a genuinely-3D
      ! pseudo-1D mesh; a plain absolute-vs-relative pivot-tolerance fix,
      ! tried first, made no difference, pointing at conditioning rather
      ! than a threshold). x_coeff comes back in these SCALED units and is
      ! divided back down by the matching power of Lref where grad/hess are
      ! assigned below.
      do jn = 1, n_neigh
        xcj = mesh%elem(neigh_list(jn))%coord
        ddx = xcj(1) - xci(1)
        ddy = xcj(2) - xci(2)
        ddz = xcj(3) - xci(3)
        if (boundary_2d) then
          dist = sqrt(ddx**2 + ddy**2)
        else
          dist = sqrt(ddx**2 + ddy**2 + ddz**2)
        end if
        wt = 1.0_DOUBLE / sqrt(max(dist, 1.0e-14_DOUBLE))   ! sqrt(1/dist)
        wt_arr(jn) = wt
        sdx = ddx / Lref
        sdy = ddy / Lref
        sdz = ddz / Lref

        if (boundary_2d) then
          Amat(jn, 1) = sdx * wt
          Amat(jn, 2) = sdy * wt
          if (order >= 3) then
            Amat(jn, 3) = sdx * sdx * 0.5_DOUBLE * wt
            Amat(jn, 4) = sdx * sdy             * wt
            Amat(jn, 5) = sdy * sdy * 0.5_DOUBLE * wt
          end if
        else
          Amat(jn, 1) = sdx * wt
          Amat(jn, 2) = sdy * wt
          Amat(jn, 3) = sdz * wt
          if (order >= 3) then
            Amat(jn, 4) = sdx * sdx * 0.5_DOUBLE * wt
            Amat(jn, 5) = sdx * sdy             * wt
            Amat(jn, 6) = sdx * sdz             * wt
            Amat(jn, 7) = sdy * sdy * 0.5_DOUBLE * wt
            Amat(jn, 8) = sdy * sdz             * wt
            Amat(jn, 9) = sdz * sdz * 0.5_DOUBLE * wt
          end if
        end if
      end do

      ! Weighted normal matrix A^T W A  (n_coeff x n_coeff)
      ATA(1:n_coeff, 1:n_coeff) = matmul(transpose(Amat), Amat)

      ! Solve once per primitive variable (ATb = A^T W b)
      do kv = 1, 5
        ATb(1:n_coeff) = 0.0_DOUBLE
        do jn = 1, n_neigh
          ATb(1:n_coeff) = ATb(1:n_coeff) + Amat(jn, :) * wt_arr(jn) * &
            (prim(kv, neigh_list(jn)) - prim(kv, i))
        end do

        call solve_ls(ATA, ATb, x_coeff, n_coeff)

        ! x_coeff came back in Lref-scaled units (see Amat construction
        ! above): d/dx = (1/Lref) d/d(sdx), d^2/dx^2 = (1/Lref^2) d^2/d(sdx)^2.
        grad(kv, 1, i) = x_coeff(1) / Lref
        grad(kv, 2, i) = x_coeff(2) / Lref
        if (boundary_2d) then
          grad(kv, 3, i) = 0.0_DOUBLE
        else
          grad(kv, 3, i) = x_coeff(3) / Lref
        end if

        if (order >= 3 .and. allocated(hess)) then
          hess(kv, :, :, i) = 0.0_DOUBLE
          if (boundary_2d) then
            hess(kv, 1, 1, i) = x_coeff(3) / Lref**2
            hess(kv, 1, 2, i) = x_coeff(4) / Lref**2
            hess(kv, 2, 1, i) = x_coeff(4) / Lref**2
            hess(kv, 2, 2, i) = x_coeff(5) / Lref**2
          else
            hess(kv, 1, 1, i) = x_coeff(4) / Lref**2
            hess(kv, 1, 2, i) = x_coeff(5) / Lref**2
            hess(kv, 2, 1, i) = x_coeff(5) / Lref**2
            hess(kv, 1, 3, i) = x_coeff(6) / Lref**2
            hess(kv, 3, 1, i) = x_coeff(6) / Lref**2
            hess(kv, 2, 2, i) = x_coeff(7) / Lref**2
            hess(kv, 2, 3, i) = x_coeff(8) / Lref**2
            hess(kv, 3, 2, i) = x_coeff(8) / Lref**2
            hess(kv, 3, 3, i) = x_coeff(9) / Lref**2
          end if
        end if

      end do

      deallocate(Amat)
    end do

  contains

    ! Add cell `cand` to neigh_list(1:n_neigh) (host-associated, declared
    ! above in ls_reconstruction) if it is a real, non-ghost cell other
    ! than i itself and not already present.
    subroutine add_unique_neighbor(cand)
      integer(kind=ENTIER), intent(in) :: cand
      integer(kind=ENTIER) :: jn
      logical :: found

      if (cand <= 0 .or. cand == i .or. mesh%elem(cand)%is_ghost) return
      found = .false.
      do jn = 1, n_neigh
        if (neigh_list(jn) == cand) then
          found = .true.
          exit
        end if
      end do
      if (.not. found .and. n_neigh < 500) then
        n_neigh = n_neigh + 1
        neigh_list(n_neigh) = cand
      end if
    end subroutine add_unique_neighbor

    ! Number of unknowns actually determinable from the CURRENT
    ! neigh_list(1:n_neigh) (host-associated), given which raw coordinate
    ! directions have any spread among these neighbours relative to i.
    ! x is assumed always alive (a mesh degenerate even in x is not
    ! handled specially). A degenerate y drops y itself and every
    ! quadratic term involving it (xy, yy); same for z. Reduces exactly
    ! to the nominal n_coeff (5 or 9 for order>=3, 2 or 3 otherwise) on
    ! any ordinary 2D/3D mesh where all directions are alive.
    function count_needed_coeffs() result(n_needed)
      integer(kind=ENTIER) :: n_needed
      real(kind=DOUBLE) :: spread_x, spread_y, spread_z
      real(kind=DOUBLE) :: ddxj, ddyj, ddzj
      logical :: y_alive, z_alive
      integer(kind=ENTIER) :: jn2, n_lin

      spread_x = 0.0_DOUBLE
      spread_y = 0.0_DOUBLE
      spread_z = 0.0_DOUBLE
      do jn2 = 1, n_neigh
        ddxj = mesh%elem(neigh_list(jn2))%coord(1) - xci(1)
        ddyj = mesh%elem(neigh_list(jn2))%coord(2) - xci(2)
        ddzj = mesh%elem(neigh_list(jn2))%coord(3) - xci(3)
        spread_x = max(spread_x, abs(ddxj))
        spread_y = max(spread_y, abs(ddyj))
        spread_z = max(spread_z, abs(ddzj))
      end do

      y_alive = spread_y >= 1.0e-8_DOUBLE * max(spread_x, 1.0e-300_DOUBLE)
      z_alive = (.not. boundary_2d) .and. &
        (spread_z >= 1.0e-8_DOUBLE * max(spread_x, 1.0e-300_DOUBLE))

      n_lin = 1
      if (y_alive) n_lin = n_lin + 1
      if (z_alive) n_lin = n_lin + 1

      n_needed = n_lin
      if (order >= 3) n_needed = n_needed + n_lin * (n_lin + 1) / 2
    end function count_needed_coeffs

  end subroutine ls_reconstruction

end module silvia_reconstruction_module