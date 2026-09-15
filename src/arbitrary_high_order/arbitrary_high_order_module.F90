module arbitrary_high_order_module
  ! Arbitrary-order cell derivatives on unstructured meshes, built as a
  ! recursive multi-D generalization of 1D divided differences.
  !
  ! Given the order-(k-1) derivative tensor in every cell (order 0 being the
  ! field itself), one order-k step (compute_next_order_derivative) is:
  !   1) for every vertex, a weighted-least-squares affine fit of the
  !      order-(k-1) cell tensor over the elements touching that vertex
  !      gives its gradient, i.e. the order-k tensor at the vertex. Being a
  !      genuine least-squares fit (solved via LAPACK LU, one factorization
  !      per vertex reused for all nc_in right-hand sides), it is exact
  !      (zero residual) whenever the order-(k-1) field is itself locally
  !      affine, regardless of mesh irregularity -- unlike the Green-Gauss
  !      jump formula this replaced (2026-09-12), whose consistency on a
  !      general mesh had not actually been verified and produced
  !      unreliable convergence rates. See compute_nodal_derivative_at_vertex.
  !   2) that vertex tensor is immediately scattered into every cell around
  !      the vertex as one nonlinear-WENO-weighted contribution, weight
  !      1/(eps+|D^k phi_v|^2) (strongly de-weighting a vertex whose
  !      estimate is large/oscillatory), with a running weight sum, so no
  !      (nc_out, n_vert) array of nodal values is ever stored -- one vertex
  !      is visited once, its accumulator is updated in every adjacent cell,
  !      and only the final per-cell sums are kept;
  !   3) once every vertex has been swept, each cell normalizes its
  !      accumulator by its weight sum. (2026-09-13: this replaced an
  !      earlier two-candidate blend against a separate central/unweighted
  !      average -- see the JCP paper for why that extra step was dropped:
  !      it didn't change the achieved order, only added complexity and a
  !      now-fixed indicator bug.)
  !
  ! Only the order-k tensor needs a ghost exchange before it can serve as the
  ! order-(k-1) input of the next step (one flat exchange per order, one
  ! ghost layer deep since the dual stencil around a vertex only reaches
  ! elements sharing that vertex). compute_next_order_derivative performs a
  ! single such step and leaves the exchange to the caller, so that
  ! communication of order k can be overlapped with unrelated work while
  ! order k+1 is prepared -- see compute_derivative_hierarchy for the
  ! straightforward (blocking) way to chain all orders.
  use precision_module
  use mesh_module
  use mpi_module
  implicit none

  private

  public :: derivative_field_type
  public :: n_derivative_components
  public :: compute_next_order_derivative
  public :: compute_next_order_derivative_cweno
  public :: compute_derivative_hierarchy
  public :: compute_next_order_derivative_overlap
  public :: compute_derivative_hierarchy_overlap
  public :: timeline_event_type
  public :: compute_derivative_hierarchy_timed
  public :: compute_derivative_hierarchy_overlap_timed

  ! WENO-indicator regularization. eps_weno (used elsewhere, unchanged)
  ! is just enough to avoid a literal division by zero when a candidate's
  ! tensor is exactly zero -- per the user's original request, replacing
  ! ad hoc 1e-6/1e-8 floors, and still correct for that use.
  !
  ! eps_weight_num is different (2026-09-15): it regularizes oi_v, the
  ! per-vertex WENO oscillation indicator in scatter_weno_weighted -- see
  ! compute_nodal_derivative_at_vertex's header. oi_v is now a DIMENSIONLESS
  ! relative residual (the fit's own mean-square misfit divided by the
  ! neighbors' own mean-square phi), so a fixed, scale-independent
  ! constant is the right kind of floor here, unlike a raw physical
  ! quantity -- tiny(1.0) provides no regularization at all against
  ! floating-point-level residual noise once oi_v itself is this small for
  ! smooth data (confirmed: it left the Hessian error on a smooth Gaussian
  ! test field completely flat under mesh refinement, not converging at
  ! all -- 1e-6 restored clean O(h^2) convergence, verified against the
  ! uniform-weight reference). 1e-6 matches the standard WENO/Jiang-Shu
  ! convention for a dimensionless epsilon.
  real(kind=DOUBLE), parameter :: eps_weno = tiny(1.0_DOUBLE)
  real(kind=DOUBLE), parameter :: eps_weight_num = 1.0e-6_DOUBLE
  ! Exponent on the oscillation indicator in the WENO weight,
  ! weight = 1/(eps+|dphi_v|^weno_power) -- 2026-09-13 experiment, per the
  ! user's request, to see whether a more aggressive de-centering (in the
  ! absence of any real positivity/slope limiter) reduces the small
  ! overshoots seen e.g. just behind a Sod rarefaction tail. Standard WENO
  ! practice is power=2; try higher powers here directly rather than
  ! building an actual limiter.
  integer(kind=ENTIER), parameter :: weno_power = 2

  ! Linear-blend mode (2026-09-15, per the user's request for a baseline
  ! comparison against the pure-ls reconstruction source): when .false.,
  ! the nonlinear WENO weight below is replaced by a plain volume weight
  ! (weight=1, so vertex_weno_weight cancels out of both the numerator and
  ! denominator, leaving a straight sub_elem_volume-weighted average of
  ! nodal derivatives -- no oscillation-adaptive de-centering at all,
  ! analogous to what ls_reconstruction's raw, unlimited fit represents,
  ! but built from aho's nodal-derivative machinery instead of a
  ! cell-centered polynomial fit). Set by euler_ho_module's
  ! use_weno_blend namelist flag before calling aho_reconstruction;
  ! defaults to .true. (the normal, published WENO behavior) everywhere
  ! else, including arbitrary_high_order_main's own standalone tests.
  logical, public :: use_weno_blend = .true.

  ! CWENO-style central-candidate experiment (2026-09-15, per the user's
  ! request): classical CWENO (Semplice & Visconti 2020, the same
  ! reference already cited above for the OI fix) doesn't just blend
  ! per-vertex candidates against each other -- it ALSO includes one
  ! "optimal"/central candidate (here: the plain linear/volume-weighted
  ! blend of the same per-vertex estimates, i.e. exactly what
  ! use_weno_blend=.false. already computes) in the SAME nonlinear
  ! blend, given a large ideal/linear weight so it dominates whenever
  ! the neighborhood is smooth, but still yields to the more robust
  ! per-vertex candidates when its own oscillation indicator blows up
  ! near a real discontinuity. Off by default (use
  ! compute_next_order_derivative_cweno explicitly to opt in); not yet
  ! wired into aho_reconstruction/euler_ho_module.
  logical, public :: use_cweno_center = .false.
  real(kind=DOUBLE), public :: cweno_center_weight = 1000.0_DOUBLE

  ! Per-vertex neighbor-list cache (2026-09-14): gather_ls_neighbors'
  ! ring-expansion + sort-based uniqueness only depends on mesh topology
  ! (id_vert, mesh) -- never on phi, nc_in, or which RK stage/hierarchy
  ! order is calling -- yet compute_nodal_derivative_at_vertex called it
  ! fresh every single time. Profiling (perf, order-2 vortex, aho) found
  ! the sort_module calls it makes (qsort_key_r + add_sort_unique_int)
  ! consuming ~58% of TOTAL solver runtime, versus ~7% for the actual
  ! LAPACK LU solve -- almost the entire per-vertex cost was re-deriving
  ! the same static answer over and over. Cached here in CSR form,
  ! (re)built once per distinct mesh (by n_elems as a cheap fingerprint)
  ! and reused for the life of the run; the per-call cost is now just an
  ! O(1) index lookup into it.
  integer(kind=ENTIER), save :: neigh_cache_n_vert = -1
  integer(kind=ENTIER), dimension(:), allocatable, save :: neigh_cache_start
  integer(kind=ENTIER), dimension(:), allocatable, save :: neigh_cache_list

  type :: derivative_field_type
    real(kind=DOUBLE), dimension(:, :), allocatable :: val ! (d**order, n_elems)
  end type derivative_field_type

  ! One phase of one order's work, timestamped with MPI_WTIME() by the
  ! *_timed drivers below, purely for the blocking-vs-overlap timeline
  ! figure in tex_arbitrary_high_order -- no effect on any other routine.
  type :: timeline_event_type
    integer(kind=ENTIER) :: order
    character(len=20)    :: phase
    real(kind=DOUBLE)    :: t0, t1
  end type timeline_event_type

contains

  pure function n_derivative_components(d, order) result(nc)
    implicit none

    integer(kind=ENTIER), intent(in) :: d, order
    integer(kind=ENTIER) :: nc

    nc = d**order
  end function n_derivative_components

  ! One step of the hierarchy: nc_in = d**(order-1) components/cell in, and
  ! nc_out = nc_in*d components/cell out. phi must already be valid on the
  ! single ghost layer (see compute_derivative_hierarchy). boundary_2d must
  ! match the flag the mesh was built with (mesh_geometry_module/build_mesh):
  ! on a mesh extruded as a single, thin layer in z (subfv's convention for
  ! 2D cases), every vertex's element neighbors sit at the same z offset
  ! from it, which makes the full 3D least-squares fit of Step~1 singular
  ! -- boundary_2d=.true. drops z from the fit basis instead.
  subroutine compute_next_order_derivative(mesh, d, nc_in, boundary_2d, phi, dphi)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(nc_in, mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(nc_in*d, mesh%n_elems), intent(out) :: dphi

    integer(kind=ENTIER) :: nc_out, id_vert, id_elem
    real(kind=DOUBLE), dimension(:, :), allocatable :: weno_num
    real(kind=DOUBLE), dimension(:), allocatable :: weno_den
    real(kind=DOUBLE), dimension(:, :), allocatable :: dphi_v_cache
    logical, dimension(:), allocatable :: valid_cache
    real(kind=DOUBLE), dimension(:), allocatable :: oi_cache

    nc_out = nc_in*d

    allocate(weno_num(nc_out, mesh%n_elems))
    allocate(weno_den(mesh%n_elems))
    allocate(dphi_v_cache(nc_out, mesh%n_vert))
    allocate(valid_cache(mesh%n_vert))
    allocate(oi_cache(mesh%n_vert))
    weno_num = 0.0_DOUBLE
    weno_den = 0.0_DOUBLE

    ! Single pass: each vertex's nodal derivative is computed once and
    ! immediately scattered, WENO-weighted, into every cell touching it --
    ! see accumulate_weno_contribution. A vertex ON THE PHYSICAL DOMAIN
    ! BOUNDARY (mesh%vert%is_bound) is skipped here (2026-09-15, per the
    ! user): its dual-cell neighbor gather only sees interior elements (no
    ! mirror/ghost cell across the wall), so the weighted-least-squares fit
    ! at that vertex is a one-sided sample -- correct along the wall-normal
    ! direction (real physical variation to fit), but along directions
    ! tangent to the wall it turns a field with genuinely zero tangential
    ! variation into a spurious nonzero derivative purely from the
    ! neighbor set's asymmetry (confirmed on Sod's 100x10x10 mesh: the
    ! order-3 aho Hessian picked up O(0.1-0.8) spurious y/z cross-terms at
    ! a y=const/z=const wall-adjacent cell, where ls_reconstruction's
    ! direct per-cell quadratic fit correctly gives ~1e-9). Every cell has
    ! at least one strictly-interior vertex (not touching any domain
    ! boundary face) unless the mesh is only 1 cell deep in some direction
    ! that boundary_2d was not told to drop -- boundary_2d already handles
    ! the single-thin-layer-in-z case by excluding z from the fit basis
    ! entirely (see header), so is_bound is only skipped when .not.
    ! boundary_2d; skipping it there too would starve every vertex at once
    ! (both z faces touch every vertex) and divide by zero below.
    do id_vert = 1, mesh%n_vert
      if (mesh%vert(id_vert)%is_bound) cycle
      call accumulate_weno_contribution(mesh, d, nc_in, boundary_2d, id_vert, &
        phi, dphi_v_cache, valid_cache, oi_cache, weno_num, weno_den)
    end do

    ! Fallback for a cell left with zero total weight (e.g. every one of
    ! its vertices sits on the domain boundary -- possible in a corner
    ! region of a very coarse mesh, or every touching vertex was invalid,
    ! see compute_nodal_derivative_at_vertex): re-admit boundary vertices
    ! for that cell only, rather than dividing by zero.
    do id_elem = 1, mesh%n_elems
      if (weno_den(id_elem) == 0.0_DOUBLE) then
        call rescue_zero_weight_cell(mesh, d, nc_in, boundary_2d, id_elem, &
          phi, weno_num, weno_den)
      end if
    end do

    do id_elem = 1, mesh%n_elems
      ! A fully isolated cell (every touching vertex invalid, even after
      ! rescue) has no resolvable direction anywhere around it -- report a
      ! zero derivative rather than divide by zero (see rescue's own note).
      if (weno_den(id_elem) == 0.0_DOUBLE) then
        dphi(:, id_elem) = 0.0_DOUBLE
      else
        dphi(:, id_elem) = weno_num(:, id_elem) / weno_den(id_elem)
      end if
    end do

    deallocate(weno_num, weno_den, dphi_v_cache, valid_cache, oi_cache)
  end subroutine compute_next_order_derivative

  ! EXPERIMENTAL (2026-09-15): CWENO variant of compute_next_order_derivative
  ! above -- see use_cweno_center's header for the idea. Accumulates, in one
  ! pass (one LS solve per vertex, same as the plain routine), THREE running
  ! sums per cell: the usual nonlinear WENO num/den, a LINEAR (unweighted by
  ! OI) num/den built from the exact same per-vertex candidates, and a
  ! volume-weighted running average of the vertices' own OI. The final
  ! per-cell derivative folds the linear candidate back into the nonlinear
  ! blend with weight cweno_center_weight/(eps+OI_c**power), OI_c being that
  ! cell's own volume-weighted-average OI -- i.e. the central/linear
  ! candidate acts as one more entry in the same WENO sum, just with a much
  ! larger ideal weight so it wins whenever OI_c is small (smooth data).
  subroutine compute_next_order_derivative_cweno(mesh, d, nc_in, boundary_2d, phi, dphi)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(nc_in, mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(nc_in*d, mesh%n_elems), intent(out) :: dphi

    integer(kind=ENTIER) :: nc_out, id_vert, id_elem
    real(kind=DOUBLE), dimension(:, :), allocatable :: weno_num, lin_num
    real(kind=DOUBLE), dimension(:), allocatable :: weno_den, lin_den, oi_num, oi_den
    real(kind=DOUBLE), dimension(:, :), allocatable :: dphi_v_cache
    logical, dimension(:), allocatable :: valid_cache
    real(kind=DOUBLE), dimension(:), allocatable :: oi_cache
    real(kind=DOUBLE) :: oi_c, w_c
    real(kind=DOUBLE), dimension(:), allocatable :: dphi_lin_elem

    nc_out = nc_in*d

    allocate(weno_num(nc_out, mesh%n_elems), lin_num(nc_out, mesh%n_elems))
    allocate(weno_den(mesh%n_elems), lin_den(mesh%n_elems))
    allocate(oi_num(mesh%n_elems), oi_den(mesh%n_elems))
    allocate(dphi_v_cache(nc_out, mesh%n_vert))
    allocate(valid_cache(mesh%n_vert))
    allocate(oi_cache(mesh%n_vert))
    allocate(dphi_lin_elem(nc_out))
    weno_num = 0.0_DOUBLE; weno_den = 0.0_DOUBLE
    lin_num = 0.0_DOUBLE; lin_den = 0.0_DOUBLE
    oi_num = 0.0_DOUBLE; oi_den = 0.0_DOUBLE

    do id_vert = 1, mesh%n_vert
      if (mesh%vert(id_vert)%is_bound) cycle
      call accumulate_cweno_contribution(mesh, d, nc_in, boundary_2d, id_vert, &
        phi, dphi_v_cache, valid_cache, oi_cache, weno_num, weno_den, lin_num, lin_den, &
        oi_num, oi_den)
    end do

    do id_elem = 1, mesh%n_elems
      if (weno_den(id_elem) == 0.0_DOUBLE) then
        call rescue_zero_weight_cell(mesh, d, nc_in, boundary_2d, id_elem, &
          phi, weno_num, weno_den)
        ! Rescue's one-sided fallback has no linear/OI counterpart -- just
        ! mirror it into the linear accumulators too so this cell's final
        ! blend below still has a well-defined (if degenerate) center.
        lin_num(:, id_elem) = weno_num(:, id_elem)
        lin_den(id_elem) = weno_den(id_elem)
      end if
    end do

    do id_elem = 1, mesh%n_elems
      if (weno_den(id_elem) == 0.0_DOUBLE) then
        dphi(:, id_elem) = 0.0_DOUBLE
        cycle
      end if
      if (lin_den(id_elem) > 0.0_DOUBLE .and. oi_den(id_elem) > 0.0_DOUBLE) then
        oi_c = oi_num(id_elem) / oi_den(id_elem)
        w_c = cweno_center_weight / (eps_weight_num + oi_c**weno_power)
        dphi_lin_elem = lin_num(:, id_elem) / lin_den(id_elem)
        dphi(:, id_elem) = (weno_num(:, id_elem) + w_c * dphi_lin_elem) / (weno_den(id_elem) + w_c)
      else
        dphi(:, id_elem) = weno_num(:, id_elem) / weno_den(id_elem)
      end if
    end do

    deallocate(weno_num, weno_den, lin_num, lin_den, oi_num, oi_den)
    deallocate(dphi_v_cache, valid_cache, oi_cache, dphi_lin_elem)
  end subroutine compute_next_order_derivative_cweno

  ! Same per-vertex LS solve as accumulate_weno_contribution, but scatters
  ! into three running sums at once (nonlinear WENO, linear/center, and a
  ! volume-weighted running OI average) -- see
  ! compute_next_order_derivative_cweno's header.
  subroutine accumulate_cweno_contribution(mesh, d, nc_in, boundary_2d, id_vert, &
      phi, dphi_v_cache, valid_cache, oi_cache, weno_num, weno_den, lin_num, lin_den, &
      oi_num, oi_den)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in, id_vert
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(nc_in, mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(:, :), intent(inout) :: dphi_v_cache
    logical, dimension(:), intent(inout) :: valid_cache
    real(kind=DOUBLE), dimension(:), intent(inout) :: oi_cache
    real(kind=DOUBLE), dimension(:, :), intent(inout) :: weno_num, lin_num
    real(kind=DOUBLE), dimension(:), intent(inout) :: weno_den, lin_den, oi_num, oi_den

    real(kind=DOUBLE), dimension(nc_in*d) :: dphi_v
    logical :: valid
    real(kind=DOUBLE) :: oi_v
    integer(kind=ENTIER) :: j, id_elem, id_sub_elem
    real(kind=DOUBLE) :: sub_elem_volume, vertex_weno_weight

    call compute_nodal_derivative_at_vertex(mesh, d, nc_in, boundary_2d, &
      id_vert, phi, dphi_v, valid, oi_v)
    dphi_v_cache(:, id_vert) = dphi_v
    valid_cache(id_vert) = valid
    oi_cache(id_vert) = oi_v
    if (.not. valid) return

    vertex_weno_weight = 1.0_DOUBLE / (eps_weight_num + oi_v**weno_power)

    do j = 1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sub_elem_volume = mesh%sub_elem(id_sub_elem)%volume

      weno_num(:, id_elem) = weno_num(:, id_elem) + (sub_elem_volume * vertex_weno_weight) * dphi_v
      weno_den(id_elem) = weno_den(id_elem) + sub_elem_volume * vertex_weno_weight

      lin_num(:, id_elem) = lin_num(:, id_elem) + sub_elem_volume * dphi_v
      lin_den(id_elem) = lin_den(id_elem) + sub_elem_volume

      oi_num(id_elem) = oi_num(id_elem) + sub_elem_volume * oi_v
      oi_den(id_elem) = oi_den(id_elem) + sub_elem_volume
    end do
  end subroutine accumulate_cweno_contribution

  ! Rescue path for a cell that ended up with zero total WENO weight after
  ! is_bound vertices were skipped above -- only possible if every one of
  ! the cell's own vertices sits on the domain boundary. Re-admits just
  ! this cell's own boundary vertices (their nodal derivative is still a
  ! valid, if one-sided, estimate) so the caller never divides by zero.
  subroutine rescue_zero_weight_cell(mesh, d, nc_in, boundary_2d, id_elem, &
      phi, weno_num, weno_den)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in, id_elem
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(nc_in, mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(:, :), intent(inout) :: weno_num
    real(kind=DOUBLE), dimension(:), intent(inout) :: weno_den

    integer(kind=ENTIER) :: j, id_vert, id_sub_elem
    real(kind=DOUBLE), dimension(nc_in*d) :: dphi_v
    real(kind=DOUBLE) :: sub_elem_volume, vertex_weno_weight, oi_v
    logical :: valid

    do j = 1, mesh%elem(id_elem)%n_vert
      id_vert = mesh%elem(id_elem)%vert(j)
      id_sub_elem = mesh%elem(id_elem)%sub_elem(j)
      call compute_nodal_derivative_at_vertex(mesh, d, nc_in, boundary_2d, &
        id_vert, phi, dphi_v, valid, oi_v)
      if (.not. valid) cycle ! see compute_nodal_derivative_at_vertex's header
      sub_elem_volume = mesh%sub_elem(id_sub_elem)%volume
      if (use_weno_blend) then
        vertex_weno_weight = 1.0_DOUBLE / (eps_weight_num + oi_v**weno_power)
      else
        vertex_weno_weight = 1.0_DOUBLE
      end if
      weno_num(:, id_elem) = weno_num(:, id_elem) &
        + (sub_elem_volume * vertex_weno_weight) * dphi_v
      weno_den(id_elem) = weno_den(id_elem) + sub_elem_volume * vertex_weno_weight
    end do
    ! Every one of this cell's vertices was invalid too (a fully isolated
    ! cell/vertex cluster with no resolvable direction anywhere) -- leave
    ! weno_den at 0; the caller must not divide by it (see its own guard).
  end subroutine rescue_zero_weight_cell

  ! Same order-k step as compute_next_order_derivative, but overlapping the
  ! ghost exchange of its own output with local work, as follows:
  !   1) sweep only the vertices touching a cell that must be sent to
  !      another MPI rank (is_send_cell(v)'s union), and finalize just
  !      those cells into dphi;
  !   2) post that exchange non-blockingly (mpi_memory_exchange_post);
  !   3) sweep every remaining vertex and finalize every remaining cell
  !      while the exchange is in flight;
  !   4) wait for the exchange to complete, filling in this rank's ghost
  !      cells of dphi.
  ! With num_procs=1 (no send/recv neighbors) this degenerates to exactly
  ! compute_next_order_derivative (step 1 touches nothing, everything runs
  ! in step 3, and post/wait are no-ops). See compute_derivative_hierarchy_overlap
  ! for the driver and Section on 3D MPI overlap performance in
  ! tex_arbitrary_high_order for the measured benefit.
  subroutine compute_next_order_derivative_overlap(mesh, mpi_send_recv, num_procs, &
      d, nc_in, boundary_2d, phi, dphi)
    use mpi_module, only: mpi_memory_exchange_post, mpi_memory_exchange_wait
    implicit none

    type(mesh_type), intent(in) :: mesh
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv
    integer(kind=ENTIER), intent(in) :: num_procs, d, nc_in
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(nc_in, mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(nc_in*d, mesh%n_elems), intent(out) :: dphi

    integer(kind=ENTIER) :: nc_out, id_vert, id_elem, i, j
    real(kind=DOUBLE), dimension(:, :), allocatable :: weno_num
    real(kind=DOUBLE), dimension(:), allocatable :: weno_den
    real(kind=DOUBLE), dimension(:, :), allocatable :: dphi_v_cache
    logical, dimension(:), allocatable :: valid_cache
    real(kind=DOUBLE), dimension(:), allocatable :: oi_cache
    logical, dimension(:), allocatable :: is_send_cell, is_boundary_vertex

    nc_out = nc_in*d

    allocate(weno_num(nc_out, mesh%n_elems))
    allocate(weno_den(mesh%n_elems))
    allocate(dphi_v_cache(nc_out, mesh%n_vert))
    allocate(valid_cache(mesh%n_vert))
    allocate(oi_cache(mesh%n_vert))
    weno_num = 0.0_DOUBLE
    weno_den = 0.0_DOUBLE

    ! With a single rank there is nothing to send/receive and mpi_send_recv
    ! is not meaningfully populated (same convention as
    ! compute_derivative_hierarchy's num_procs>1 guard on
    ! mpi_memory_exchange): treat every cell as "not sent", so everything
    ! runs in Step 3 below and this degenerates to compute_next_order_derivative.
    allocate(is_send_cell(mesh%n_elems))
    is_send_cell = .false.
    if (num_procs > 1) then
      do i = 1, mpi_send_recv%n_mpi_send_neigh
        do j = 1, mpi_send_recv%mpi_send_neigh(i)%n_elems
          is_send_cell(mpi_send_recv%mpi_send_neigh(i)%elem_id(j)) = .true.
        end do
      end do
    end if

    allocate(is_boundary_vertex(mesh%n_vert))
    is_boundary_vertex = .false.
    do id_elem = 1, mesh%n_elems
      if (is_send_cell(id_elem)) then
        do j = 1, mesh%elem(id_elem)%n_vert
          is_boundary_vertex(mesh%elem(id_elem)%vert(j)) = .true.
        end do
      end if
    end do

    ! Step 1: vertices touching a to-be-sent cell -- compute+cache their
    ! nodal derivative and accumulate it, WENO-weighted, into the send
    ! cells only (skip_cell filters at the cell level). Physical-domain-
    ! boundary vertices are skipped throughout (Steps 1 and 3 alike, see
    ! compute_next_order_derivative for why), so dphi_v_cache is never
    ! read for one in Step 3's cached-reuse pass below.
    do id_vert = 1, mesh%n_vert
      if (.not. is_boundary_vertex(id_vert)) cycle
      if (mesh%vert(id_vert)%is_bound) cycle
      call accumulate_weno_contribution(mesh, d, nc_in, boundary_2d, id_vert, &
        phi, dphi_v_cache, valid_cache, oi_cache, weno_num, weno_den, skip_cell=(.not. is_send_cell))
    end do
    do id_elem = 1, mesh%n_elems
      if (.not. is_send_cell(id_elem)) cycle
      if (weno_den(id_elem) == 0.0_DOUBLE) then
        call rescue_zero_weight_cell(mesh, d, nc_in, boundary_2d, id_elem, &
          phi, weno_num, weno_den)
      end if
      if (weno_den(id_elem) == 0.0_DOUBLE) then
        dphi(:, id_elem) = 0.0_DOUBLE
      else
        dphi(:, id_elem) = weno_num(:, id_elem) / weno_den(id_elem)
      end if
    end do

    ! Step 2: hand the just-finalized boundary cells to MPI and move on
    ! without waiting.
    if (num_procs > 1) call mpi_memory_exchange_post(mpi_send_recv, mesh%n_elems, nc_out, dphi)

    ! Step 3: every other vertex/cell, computed while the exchange above is
    ! in flight. Non-boundary vertices compute+cache+accumulate fresh;
    ! boundary vertices reuse their Step 1 cached dphi_v (accumulate_cached_
    ! weno_contribution) since a non-send cell can share a boundary vertex
    ! with a send cell and still needs that vertex's contribution.
    do id_vert = 1, mesh%n_vert
      if (is_boundary_vertex(id_vert)) cycle
      if (mesh%vert(id_vert)%is_bound) cycle
      call accumulate_weno_contribution(mesh, d, nc_in, boundary_2d, id_vert, &
        phi, dphi_v_cache, valid_cache, oi_cache, weno_num, weno_den, skip_cell=is_send_cell)
    end do
    do id_vert = 1, mesh%n_vert
      if (.not. is_boundary_vertex(id_vert)) cycle
      if (mesh%vert(id_vert)%is_bound) cycle
      call accumulate_cached_weno_contribution(mesh, id_vert, dphi_v_cache, &
        valid_cache, oi_cache, weno_num, weno_den, skip_cell=is_send_cell)
    end do
    do id_elem = 1, mesh%n_elems
      if (is_send_cell(id_elem)) cycle
      if (weno_den(id_elem) == 0.0_DOUBLE) then
        call rescue_zero_weight_cell(mesh, d, nc_in, boundary_2d, id_elem, &
          phi, weno_num, weno_den)
      end if
      if (weno_den(id_elem) == 0.0_DOUBLE) then
        dphi(:, id_elem) = 0.0_DOUBLE
      else
        dphi(:, id_elem) = weno_num(:, id_elem) / weno_den(id_elem)
      end if
    end do

    ! Step 4: this rank's ghost cells of dphi are only valid past this point.
    if (num_procs > 1) call mpi_memory_exchange_wait(mpi_send_recv, mesh%n_elems, nc_out, dphi)

    deallocate(weno_num, weno_den, dphi_v_cache, valid_cache, oi_cache)
    deallocate(is_send_cell, is_boundary_vertex)
  end subroutine compute_next_order_derivative_overlap

  ! Reconstructs every order 1..max_order from the scalar field phi0.
  ! Each order's cell-centered input is exchanged across the single ghost
  ! layer (mpi_memory_exchange, reused as-is from mpi_module) before it is
  ! used to build the next, node-centered-then-primal, order -- so only a
  ! d**(order-1)-component field ever needs to cross an MPI boundary to
  ! advance one order. ns/ can call compute_next_order_derivative directly
  ! instead of this driver, and pipeline its own non-blocking exchange of
  ! one order with other work while the next order is prepared.
  subroutine compute_derivative_hierarchy(mesh, mpi_send_recv, num_procs, &
      d, boundary_2d, max_order, phi0, dfield)
    implicit none

    type(mesh_type), intent(in) :: mesh
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv
    integer(kind=ENTIER), intent(in) :: num_procs, d, max_order
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: phi0
    type(derivative_field_type), dimension(max_order), intent(out) :: dfield

    integer(kind=ENTIER) :: order, nc_in
    real(kind=DOUBLE), dimension(:, :), allocatable :: phi_prev

    allocate(phi_prev(1, mesh%n_elems))
    phi_prev(1, :) = phi0

    do order = 1, max_order
      nc_in = d**(order-1)

      if (num_procs > 1) then
        call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, nc_in, phi_prev)
      end if

      allocate(dfield(order)%val(nc_in*d, mesh%n_elems))
      call compute_next_order_derivative(mesh, d, nc_in, boundary_2d, &
        phi_prev, dfield(order)%val)

      deallocate(phi_prev)
      allocate(phi_prev(nc_in*d, mesh%n_elems))
      phi_prev = dfield(order)%val
    end do

    deallocate(phi_prev)
  end subroutine compute_derivative_hierarchy

  ! Same as compute_derivative_hierarchy, but using
  ! compute_next_order_derivative_overlap for each order: since that
  ! routine already waits for its own output's exchange before returning,
  ! phi_prev is guaranteed ghost-complete for the next order without a
  ! separate top-of-loop exchange -- except for phi0 itself (order 0),
  ! which nothing "produces" and so is exchanged once, plainly, up front.
  subroutine compute_derivative_hierarchy_overlap(mesh, mpi_send_recv, num_procs, &
      d, boundary_2d, max_order, phi0, dfield)
    implicit none

    type(mesh_type), intent(in) :: mesh
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv
    integer(kind=ENTIER), intent(in) :: num_procs, d, max_order
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: phi0
    type(derivative_field_type), dimension(max_order), intent(out) :: dfield

    integer(kind=ENTIER) :: order, nc_in
    real(kind=DOUBLE), dimension(:, :), allocatable :: phi_prev

    allocate(phi_prev(1, mesh%n_elems))
    phi_prev(1, :) = phi0
    if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 1, phi_prev)

    do order = 1, max_order
      nc_in = d**(order-1)

      allocate(dfield(order)%val(nc_in*d, mesh%n_elems))
      call compute_next_order_derivative_overlap(mesh, mpi_send_recv, num_procs, &
        d, nc_in, boundary_2d, phi_prev, dfield(order)%val)

      deallocate(phi_prev)
      allocate(phi_prev(nc_in*d, mesh%n_elems))
      phi_prev = dfield(order)%val
    end do

    deallocate(phi_prev)
  end subroutine compute_derivative_hierarchy_overlap

  ! Same as compute_derivative_hierarchy, but timestamping the two phases
  ! of every order (blocking exchange, then compute) into events(1:n_events)
  ! via MPI_WTIME() -- for the blocking-vs-overlap timeline figure only.
  ! events must be preallocated by the caller to at least 2*max_order.
  subroutine compute_derivative_hierarchy_timed(mesh, mpi_send_recv, num_procs, &
      d, boundary_2d, max_order, phi0, dfield, events, n_events)
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv
    integer(kind=ENTIER), intent(in) :: num_procs, d, max_order
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: phi0
    type(derivative_field_type), dimension(max_order), intent(out) :: dfield
    type(timeline_event_type), dimension(:), intent(inout) :: events
    integer(kind=ENTIER), intent(inout) :: n_events

    integer(kind=ENTIER) :: order, nc_in
    real(kind=DOUBLE), dimension(:, :), allocatable :: phi_prev
    real(kind=DOUBLE) :: t0, t1

    allocate(phi_prev(1, mesh%n_elems))
    phi_prev(1, :) = phi0

    do order = 1, max_order
      nc_in = d**(order-1)

      t0 = MPI_WTIME()
      if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, nc_in, phi_prev)
      t1 = MPI_WTIME()
      n_events = n_events + 1
      events(n_events) = timeline_event_type(order, 'exchange', t0, t1)

      allocate(dfield(order)%val(nc_in*d, mesh%n_elems))
      t0 = MPI_WTIME()
      call compute_next_order_derivative(mesh, d, nc_in, boundary_2d, &
        phi_prev, dfield(order)%val)
      t1 = MPI_WTIME()
      n_events = n_events + 1
      events(n_events) = timeline_event_type(order, 'compute', t0, t1)

      deallocate(phi_prev)
      allocate(phi_prev(nc_in*d, mesh%n_elems))
      phi_prev = dfield(order)%val
    end do

    deallocate(phi_prev)
  end subroutine compute_derivative_hierarchy_timed

  ! Same as compute_derivative_hierarchy_overlap, but timestamping the four
  ! phases of every order (compute boundary cells, post, compute interior,
  ! wait) into events(1:n_events) -- for the timeline figure only. events
  ! must be preallocated by the caller to at least 4*max_order.
  subroutine compute_derivative_hierarchy_overlap_timed(mesh, mpi_send_recv, num_procs, &
      d, boundary_2d, max_order, phi0, dfield, events, n_events)
    use mpi
    use mpi_module, only: mpi_memory_exchange_post, mpi_memory_exchange_wait
    implicit none

    type(mesh_type), intent(in) :: mesh
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv
    integer(kind=ENTIER), intent(in) :: num_procs, d, max_order
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: phi0
    type(derivative_field_type), dimension(max_order), intent(out) :: dfield
    type(timeline_event_type), dimension(:), intent(inout) :: events
    integer(kind=ENTIER), intent(inout) :: n_events

    integer(kind=ENTIER) :: order, nc_in, nc_out, id_vert, id_elem, i, j
    real(kind=DOUBLE), dimension(:, :), allocatable :: phi_prev
    real(kind=DOUBLE), dimension(:, :), allocatable :: weno_num
    real(kind=DOUBLE), dimension(:), allocatable :: weno_den
    real(kind=DOUBLE), dimension(:, :), allocatable :: dphi_v_cache
    logical, dimension(:), allocatable :: valid_cache
    real(kind=DOUBLE), dimension(:), allocatable :: oi_cache
    logical, dimension(:), allocatable :: is_send_cell, is_boundary_vertex
    real(kind=DOUBLE) :: t0, t1

    allocate(phi_prev(1, mesh%n_elems))
    phi_prev(1, :) = phi0
    if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 1, phi_prev)

    do order = 1, max_order
      nc_in  = d**(order-1)
      nc_out = nc_in*d
      allocate(dfield(order)%val(nc_in*d, mesh%n_elems))

      allocate(weno_num(nc_out, mesh%n_elems))
      allocate(weno_den(mesh%n_elems))
      allocate(dphi_v_cache(nc_out, mesh%n_vert))
      allocate(valid_cache(mesh%n_vert))
      allocate(oi_cache(mesh%n_vert))
      weno_num = 0.0_DOUBLE; weno_den = 0.0_DOUBLE

      allocate(is_send_cell(mesh%n_elems))
      is_send_cell = .false.
      if (num_procs > 1) then
        do i = 1, mpi_send_recv%n_mpi_send_neigh
          do j = 1, mpi_send_recv%mpi_send_neigh(i)%n_elems
            is_send_cell(mpi_send_recv%mpi_send_neigh(i)%elem_id(j)) = .true.
          end do
        end do
      end if
      allocate(is_boundary_vertex(mesh%n_vert))
      is_boundary_vertex = .false.
      do id_elem = 1, mesh%n_elems
        if (is_send_cell(id_elem)) then
          do j = 1, mesh%elem(id_elem)%n_vert
            is_boundary_vertex(mesh%elem(id_elem)%vert(j)) = .true.
          end do
        end if
      end do

      t0 = MPI_WTIME()
      do id_vert = 1, mesh%n_vert
        if (.not. is_boundary_vertex(id_vert)) cycle
        if (mesh%vert(id_vert)%is_bound) cycle
        call accumulate_weno_contribution(mesh, d, nc_in, boundary_2d, id_vert, &
          phi_prev, dphi_v_cache, valid_cache, oi_cache, weno_num, weno_den, skip_cell=(.not. is_send_cell))
      end do
      do id_elem = 1, mesh%n_elems
        if (.not. is_send_cell(id_elem)) cycle
        if (weno_den(id_elem) == 0.0_DOUBLE) then
          call rescue_zero_weight_cell(mesh, d, nc_in, boundary_2d, id_elem, &
            phi_prev, weno_num, weno_den)
        end if
        if (weno_den(id_elem) == 0.0_DOUBLE) then
          dfield(order)%val(:, id_elem) = 0.0_DOUBLE
        else
          dfield(order)%val(:, id_elem) = weno_num(:, id_elem) / weno_den(id_elem)
        end if
      end do
      t1 = MPI_WTIME()
      n_events = n_events + 1
      events(n_events) = timeline_event_type(order, 'compute_bnd', t0, t1)

      t0 = MPI_WTIME()
      if (num_procs > 1) call mpi_memory_exchange_post(mpi_send_recv, mesh%n_elems, nc_out, dfield(order)%val)
      t1 = MPI_WTIME()
      n_events = n_events + 1
      events(n_events) = timeline_event_type(order, 'post', t0, t1)

      t0 = MPI_WTIME()
      do id_vert = 1, mesh%n_vert
        if (is_boundary_vertex(id_vert)) cycle
        if (mesh%vert(id_vert)%is_bound) cycle
        call accumulate_weno_contribution(mesh, d, nc_in, boundary_2d, id_vert, &
          phi_prev, dphi_v_cache, valid_cache, oi_cache, weno_num, weno_den, skip_cell=is_send_cell)
      end do
      do id_vert = 1, mesh%n_vert
        if (.not. is_boundary_vertex(id_vert)) cycle
        if (mesh%vert(id_vert)%is_bound) cycle
        call accumulate_cached_weno_contribution(mesh, id_vert, dphi_v_cache, &
          valid_cache, oi_cache, weno_num, weno_den, skip_cell=is_send_cell)
      end do
      do id_elem = 1, mesh%n_elems
        if (is_send_cell(id_elem)) cycle
        if (weno_den(id_elem) == 0.0_DOUBLE) then
          call rescue_zero_weight_cell(mesh, d, nc_in, boundary_2d, id_elem, &
            phi_prev, weno_num, weno_den)
        end if
        if (weno_den(id_elem) == 0.0_DOUBLE) then
          dfield(order)%val(:, id_elem) = 0.0_DOUBLE
        else
          dfield(order)%val(:, id_elem) = weno_num(:, id_elem) / weno_den(id_elem)
        end if
      end do
      t1 = MPI_WTIME()
      n_events = n_events + 1
      events(n_events) = timeline_event_type(order, 'compute_int', t0, t1)

      t0 = MPI_WTIME()
      if (num_procs > 1) call mpi_memory_exchange_wait(mpi_send_recv, mesh%n_elems, nc_out, dfield(order)%val)
      t1 = MPI_WTIME()
      n_events = n_events + 1
      events(n_events) = timeline_event_type(order, 'wait', t0, t1)

      deallocate(weno_num, weno_den, dphi_v_cache, valid_cache, oi_cache)
      deallocate(is_send_cell, is_boundary_vertex)
      deallocate(phi_prev)
      allocate(phi_prev(nc_in*d, mesh%n_elems))
      phi_prev = dfield(order)%val
    end do

    deallocate(phi_prev)
  end subroutine compute_derivative_hierarchy_overlap_timed

  ! Computes (D^k phi)_v at vertex id_vert, caches it in dphi_v_cache(:,
  ! id_vert) (for a later accumulate_cached_weno_contribution call on the
  ! same vertex, used by the MPI-overlap driver -- see below), and
  ! immediately scatters it into every cell touching id_vert as a
  ! nonlinear-WENO-weighted contribution: weight = 1/(eps+OI^2), OI a
  ! mesh-size-scaled oscillation indicator (see scatter_weno_weighted), the
  ! same de-weight-the-large/oscillatory-candidate idea as a classical WENO
  ! reconstruction, applied directly to each vertex's own nodal derivative
  ! with no separate central/unweighted candidate to blend against.
  !
  ! (2026-09-13: this replaced a two-candidate scheme that also computed a
  ! plain volume-weighted "central" average per cell and blended it in,
  ! weighting THIS candidate by deviation from that central one rather than
  ! by its own magnitude -- found necessary at the time because raw
  ! magnitude-weighting mis-read smooth-but-large curvature as a shock
  ! signal. That risk is unchanged here; simplicity was preferred once it
  ! was confirmed, on the vortex/Sod/cylinder suite, that the two-candidate
  ! blend did not actually improve the achieved order over this simpler
  ! form. If a future test resurfaces the old symptom -- order stuck below
  ! design value on a smooth, strongly-curved field -- that's the first
  ! place to look. 2026-09-15: that symptom resurfaced -- see
  ! scatter_weno_weighted for the actual fix, which keeps this single-
  ! candidate structure but corrects the weight formula itself instead of
  ! reverting to the two-candidate blend.)
  !
  ! skip_cell, if present, excludes cells where skip_cell(id_elem) is
  ! .true. from accumulation -- used by compute_next_order_derivative_overlap
  ! to accumulate into only the cells that should receive this vertex's
  ! contribution at the point this is called (Step 1: only send cells;
  ! Step 3: only non-send cells).
  subroutine accumulate_weno_contribution(mesh, d, nc_in, boundary_2d, id_vert, &
      phi, dphi_v_cache, valid_cache, oi_cache, weno_num, weno_den, skip_cell)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in, id_vert
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(nc_in, mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(:, :), intent(inout) :: dphi_v_cache
    logical, dimension(:), intent(inout) :: valid_cache
    real(kind=DOUBLE), dimension(:), intent(inout) :: oi_cache
    real(kind=DOUBLE), dimension(:, :), intent(inout) :: weno_num
    real(kind=DOUBLE), dimension(:), intent(inout) :: weno_den
    logical, dimension(:), intent(in), optional :: skip_cell

    real(kind=DOUBLE), dimension(nc_in*d) :: dphi_v
    logical :: valid
    real(kind=DOUBLE) :: oi_v

    call compute_nodal_derivative_at_vertex(mesh, d, nc_in, boundary_2d, &
      id_vert, phi, dphi_v, valid, oi_v)
    dphi_v_cache(:, id_vert) = dphi_v
    valid_cache(id_vert) = valid
    oi_cache(id_vert) = oi_v
    if (.not. valid) return ! see compute_nodal_derivative_at_vertex's header

    call scatter_weno_weighted(mesh, id_vert, dphi_v, oi_v, weno_num, weno_den, skip_cell)
  end subroutine accumulate_weno_contribution

  ! Same scatter as accumulate_weno_contribution, but reusing an
  ! already-cached dphi_v instead of recomputing it -- used by
  ! compute_next_order_derivative_overlap's Step 3 to fold a boundary
  ! vertex's (Step-1-computed) contribution into the non-send cells it also
  ! touches, without a second LAPACK solve.
  subroutine accumulate_cached_weno_contribution(mesh, id_vert, dphi_v_cache, &
      valid_cache, oi_cache, weno_num, weno_den, skip_cell)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(:, :), intent(in) :: dphi_v_cache
    logical, dimension(:), intent(in) :: valid_cache
    real(kind=DOUBLE), dimension(:), intent(in) :: oi_cache
    real(kind=DOUBLE), dimension(:, :), intent(inout) :: weno_num
    real(kind=DOUBLE), dimension(:), intent(inout) :: weno_den
    logical, dimension(:), intent(in), optional :: skip_cell

    if (.not. valid_cache(id_vert)) return

    call scatter_weno_weighted(mesh, id_vert, dphi_v_cache(:, id_vert), &
      oi_cache(id_vert), weno_num, weno_den, skip_cell)
  end subroutine accumulate_cached_weno_contribution

  ! Shared scatter: weight = omega_p/(eps+OI^2), omega_p=sub_elem_volume (the
  ! geometric/"linear weight" role, as in classical CWENO -- see Semplice &
  ! Visconti 2020, Def. 3) and OI=oi_v, the vertex's own LS-fit weighted
  ! residual computed in compute_nodal_derivative_at_vertex.
  !
  ! FIX (2026-09-15, per the user -- pointed at Semplice & Visconti 2020's
  ! CWENO smoothness indicator as the model to follow): a WENO oscillation
  ! indicator must measure DISAGREEMENT with a smooth local model, not raw
  ! derivative magnitude, and it must vanish under mesh refinement for
  ! smooth data (Semplice & Visconti's Definition 2/Eq.(1): "I[P] -> 0
  ! under grid refinement if P is associated to smooth data", via an
  ! explicit Delta-x^(2i-1) factor). Two formulas were tried and rejected
  ! before this one:
  !   1) OI=|dphi_v| (no correction at all): converges to the true, finite,
  !      generally-nonzero local derivative as h->0 for any field with
  !      genuine spatial variation -- never vanishes, so it keeps
  !      penalizing a vertex for having a legitimately larger derivative
  !      than its neighbors, not for being actually oscillatory.
  !   2) OI=(neighbor coordinate spread)*|dphi_v|: multiplying by a
  !      per-vertex length scale that is roughly CONSTANT across one cell's
  !      touching vertices on a uniform mesh does not change the *relative*
  !      weighting between those vertices at all (a common factor cancels
  !      in the weighted average) -- confirmed by testing it directly:
  !      identical errors to formula 1, no improvement whatsoever.
  ! The fix that actually works: use the LS fit's own weighted residual
  ! (how much the actual neighbor data disagrees with the fitted affine
  ! model, computed alongside dphi_v) as OI. This is exactly the "deviation
  ! from a smooth reference" that classical WENO/CWENO indicators measure,
  ! rather than a candidate's own absolute size: it is identically zero for
  ! genuinely affine data on any mesh at any resolution, shrinks with h for
  ! smooth-but-curved data (the affine model's own residual against a
  ! curved function shrinks as the neighborhood shrinks), and stays O(1)
  ! at a genuine discontinuity in the stencil. Confirmed by the standalone
  ! analytic quadratic-reconstruction test (see aho-mpi-ghost-stencil-bug
  ! memory): forcing a uniform (non-WENO) weight as a diagnostic dropped
  ! the interior Hessian error from O(1) (not shrinking with h) to ~1e-12
  ! -- i.e. the per-vertex fit itself was already exact, only the
  ! averaging weight was wrong; this residual-based OI achieves the same
  ! effect without disabling WENO's actual shock-detection role.
  subroutine scatter_weno_weighted(mesh, id_vert, dphi_v, oi_v, weno_num, weno_den, skip_cell)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(:), intent(in) :: dphi_v
    real(kind=DOUBLE), intent(in) :: oi_v
    real(kind=DOUBLE), dimension(:, :), intent(inout) :: weno_num
    real(kind=DOUBLE), dimension(:), intent(inout) :: weno_den
    logical, dimension(:), intent(in), optional :: skip_cell

    integer(kind=ENTIER) :: j, id_elem, id_sub_elem
    real(kind=DOUBLE) :: sub_elem_volume, vertex_weno_weight

    if (use_weno_blend) then
      vertex_weno_weight = 1.0_DOUBLE / (eps_weight_num + oi_v**weno_power)
    else
      ! Linear-blend mode: weight=1 cancels out of numerator/denominator,
      ! leaving a plain sub_elem_volume-weighted average (see use_weno_blend's
      ! header comment).
      vertex_weno_weight = 1.0_DOUBLE
    end if

    do j = 1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      if (present(skip_cell)) then
        if (skip_cell(id_elem)) cycle
      end if
      sub_elem_volume = mesh%sub_elem(id_sub_elem)%volume

      weno_num(:, id_elem) = weno_num(:, id_elem) &
        + (sub_elem_volume * vertex_weno_weight) * dphi_v
      weno_den(id_elem) = weno_den(id_elem) + sub_elem_volume * vertex_weno_weight
    end do
  end subroutine scatter_weno_weighted

  ! Weighted-least-squares gradient of phi at id_vert, over the elements
  ! touching that vertex (mesh%vert(id_vert)%elem_neigh): fits
  ! phi_c(x) ~= a0_c + G_c . (x - x_vert) for each of the nc_in scalar
  ! components c independently, all sharing the same (n_basis, n_basis)
  ! geometry-only normal matrix (one LU factorization per vertex, reused as
  ! nc_in right-hand sides).
  !
  ! The candidate directions are x,y (boundary_2d) or x,y,z (else); among
  ! those, only the ones with a genuinely resolvable spread among the
  ! gathered neighbors (relative to the largest candidate spread) enter the
  ! basis -- e.g. z is dropped by boundary_2d for a mesh extruded as a
  ! single thin layer in z (every neighbor then sits at the same z offset
  ! from id_vert), but the SAME degeneracy can hit any other direction too:
  ! a quasi-1D mesh (a single row of cells across y, as in a shock-tube
  ! test modeled as a thin 2D strip) has every cell centroid at the same y
  ! regardless of x, which silently made the fixed-2-unknown (dx,dy) fit
  ! singular here before this dynamic basis was added (2026-09-12) --
  ! dropped directions get an exact 0 derivative, same convention as
  ! boundary_2d's z. This is exact (zero residual) whenever phi is itself
  ! locally affine around id_vert in the active directions, regardless of
  ! mesh irregularity -- replaces (2026-09-12) a Green-Gauss jump formula
  ! whose consistency on a general mesh had not been verified and which
  ! produced unreliable convergence rates.
  ! The stencil is always exactly mesh%vert(id_vert)%elem_neigh -- the cells
  ! directly touching the vertex, never ring-expanded (2026-09-15, per the
  ! user): a hex vertex always has 8 of these, a quad vertex 4, a tet/tri
  ! vertex typically more -- always enough for a well-posed fit at every
  ! INTERIOR vertex. The only vertices short of that are on the physical
  ! domain boundary, and those never contribute to the reconstruction
  ! anyway (mesh%vert%is_bound is skipped by the caller, see
  ! compute_next_order_derivative) except through rescue_zero_weight_cell's
  ! deliberately one-sided fallback. Ring-expanding past the immediate
  ! elem_neigh (an earlier version of this code did, up to n_min=2*max_basis)
  ! was the actual bug behind an MPI-rank-count-dependent vortex accuracy
  ! regression: a single node-based ghost layer (subfv-gmsh's partitioning)
  ! only guarantees a *complete* elem_neigh for vertices directly on the
  ! partition seam, not for the second ring reached by expansion, so seam
  ! vertices silently lost neighbors under MPI that a serial run still had.
  ! See mesh-is-bound-timing-bug memory / tex_arbitrary_high_order for the
  ! diagnosis; confirmed via subfvns (single-ring stencil, unaffected by
  ! rank count) showing identical error at 1 and 2 ranks on the same case.
  subroutine compute_nodal_derivative_at_vertex(mesh, d, nc_in, boundary_2d, &
      id_vert, phi, dphi_v, valid, oi_v)
    use linear_solver_module, only: lu_factor_lapack, lu_solve_mat_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in, id_vert
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(nc_in, mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(nc_in*d), intent(out) :: dphi_v
    ! .false. when not even one direction could be resolved (e.g. a mesh
    ! corner vertex with too few immediate elem_neigh to fit anything) --
    ! see the caller: such a vertex must be excluded from the WENO scatter
    ! entirely, not scattered as dphi_v=0, which scatter_weno_weighted's
    ! 1/(eps+|dphi_v|^2) weight would otherwise read as a perfectly smooth
    ! (near-infinite-confidence) sample and let it swamp every genuinely
    ! resolved neighbor's contribution to the cells touching it.
    logical, intent(out) :: valid
    ! WENO oscillation indicator for this vertex's fit: the fit's own
    ! weighted-mean-squared residual (how much the actual neighbor data
    ! disagrees with the fitted affine model), maximized over the nc_in
    ! components and made DIMENSIONLESS by dividing by the neighbors' own
    ! weighted mean-square phi (so a fixed, scale-independent eps works in
    ! scatter_weno_weighted regardless of the field's physical units/
    ! magnitude -- e.g. density~O(1) vs. pressure~O(1e5)) -- see
    ! scatter_weno_weighted's header for why a residual, not |dphi_v|
    ! itself, is the right quantity, and why it must be relative.
    real(kind=DOUBLE), intent(out) :: oi_v

    integer(kind=ENTIER), parameter :: max_basis = 4
    real(kind=DOUBLE), parameter :: rel_spread_tol = 1.0e-8_DOUBLE
    integer(kind=ENTIER) :: n_basis, n_cand, n_active, j, a, b, i1, id_elem, n_neigh
    integer(kind=ENTIER), dimension(3) :: active_dim
    integer(kind=ENTIER), dimension(:), allocatable :: neigh
    integer(kind=ENTIER), dimension(max_basis) :: ipiv
    real(kind=DOUBLE), dimension(3) :: dx, dmin, dmax, spread
    real(kind=DOUBLE) :: weight, max_spread, weight_sum, phi_scale2
    real(kind=DOUBLE), dimension(max_basis) :: basis
    real(kind=DOUBLE), dimension(max_basis, max_basis) :: mat
    real(kind=DOUBLE), dimension(max_basis, nc_in) :: rhs
    real(kind=DOUBLE), dimension(nc_in) :: predicted, resid_sq, phi_sq_sum

    call ensure_neighbor_cache(mesh)
    n_neigh = neigh_cache_start(id_vert+1) - neigh_cache_start(id_vert)
    allocate(neigh(n_neigh))
    neigh = neigh_cache_list(neigh_cache_start(id_vert):neigh_cache_start(id_vert+1)-1)

    ! Which of the candidate directions (x,y[,z]) actually vary among the
    ! gathered neighbors -- see header note.
    n_cand = merge(2_ENTIER, 3_ENTIER, boundary_2d)
    dmin(1:n_cand) = huge(1.0_DOUBLE)
    dmax(1:n_cand) = -huge(1.0_DOUBLE)
    do j = 1, n_neigh
      dx = mesh%elem(neigh(j))%coord - mesh%vert(id_vert)%coord
      do a = 1, n_cand
        dmin(a) = min(dmin(a), dx(a))
        dmax(a) = max(dmax(a), dx(a))
      end do
    end do
    spread(1:n_cand) = dmax(1:n_cand) - dmin(1:n_cand)
    max_spread = maxval(spread(1:n_cand))

    n_active = 0
    do a = 1, n_cand
      if (spread(a) > rel_spread_tol * max(max_spread, 1.0e-300_DOUBLE)) then
        n_active = n_active + 1
        active_dim(n_active) = a
      end if
    end do
    n_basis = 1 + n_active

    mat = 0.0_DOUBLE
    rhs = 0.0_DOUBLE
    weight_sum = 0.0_DOUBLE
    phi_sq_sum = 0.0_DOUBLE

    do j = 1, n_neigh
      id_elem = neigh(j)
      dx = mesh%elem(id_elem)%coord - mesh%vert(id_vert)%coord
      weight = 1.0_DOUBLE / max(dot_product(dx, dx), 1.0e-24_DOUBLE)
      weight_sum = weight_sum + weight
      phi_sq_sum = phi_sq_sum + weight * phi(:, id_elem)**2

      basis(1) = 1.0_DOUBLE
      do a = 1, n_active
        basis(1+a) = dx(active_dim(a))
      end do

      do a = 1, n_basis
        do b = 1, n_basis
          mat(a, b) = mat(a, b) + weight * basis(a) * basis(b)
        end do
        rhs(a, :) = rhs(a, :) + (weight * basis(a)) * phi(:, id_elem)
      end do
    end do

    dphi_v = 0.0_DOUBLE
    oi_v = 0.0_DOUBLE
    valid = (n_active > 0)
    if (.not. valid) return ! no resolvable direction at all: report zero gradient

    call lu_factor_lapack(n_basis, mat(1:n_basis, 1:n_basis), ipiv(1:n_basis))
    call lu_solve_mat_lapack(n_basis, mat(1:n_basis, 1:n_basis), &
      ipiv(1:n_basis), nc_in, rhs(1:n_basis, :))

    ! WENO oscillation indicator: the fit's own weighted-mean-squared
    ! residual against the actual neighbor data, using the SAME weight as
    ! the fit itself -- see scatter_weno_weighted's header. Unlike |dphi_v|,
    ! this vanishes for genuinely affine-in-the-active-directions data
    ! (any mesh, any resolution) and shrinks with h for smooth-but-curved
    ! data, while staying O(1) at a real discontinuity in the stencil --
    ! exactly the property a WENO smoothness indicator needs.
    resid_sq = 0.0_DOUBLE
    do j = 1, n_neigh
      id_elem = neigh(j)
      dx = mesh%elem(id_elem)%coord - mesh%vert(id_vert)%coord
      weight = 1.0_DOUBLE / max(dot_product(dx, dx), 1.0e-24_DOUBLE)
      predicted = rhs(1, :)
      do a = 1, n_active
        predicted = predicted + rhs(1+a, :) * dx(active_dim(a))
      end do
      resid_sq = resid_sq + weight * (phi(:, id_elem) - predicted)**2
    end do
    ! phi_scale2 floored at a tiny absolute value only to avoid literal 0/0
    ! when phi is exactly zero everywhere near this vertex (e.g. a void) --
    ! in that case resid_sq is 0 too, so oi_v correctly comes out as 0
    ! regardless of the floor's exact value.
    phi_scale2 = maxval(phi_sq_sum) / max(weight_sum, 1.0e-300_DOUBLE)
    oi_v = sqrt((maxval(resid_sq) / max(weight_sum, 1.0e-300_DOUBLE)) &
      / max(phi_scale2, 1.0e-300_DOUBLE))

    ! rhs(1+a, i1) now holds d(phi_i1)/dx_{active_dim(a)}; flatten with i1
    ! (the "carried" direction, from phi's own components) fast-varying and
    ! the direction slow-varying, matching the rest of the module's tensor
    ! layout. Inactive directions stay 0.
    do a = 1, n_active
      do i1 = 1, nc_in
        dphi_v((active_dim(a)-1)*nc_in + i1) = rhs(1+a, i1)
      end do
    end do
  end subroutine compute_nodal_derivative_at_vertex

  ! Builds neigh_cache_start/neigh_cache_list (CSR) by calling
  ! gather_ls_neighbors once per vertex -- the one-time cost that used to
  ! be paid on every one of compute_nodal_derivative_at_vertex's many
  ! calls per vertex (see the cache arrays' declaration). A no-op if the
  ! cache already matches this mesh's vertex count.
  subroutine ensure_neighbor_cache(mesh)
    implicit none

    type(mesh_type), intent(in) :: mesh

    integer(kind=ENTIER) :: v, n_neigh, total
    integer(kind=ENTIER), dimension(:), allocatable :: neigh

    if (neigh_cache_n_vert == mesh%n_vert) return

    if (allocated(neigh_cache_start)) deallocate(neigh_cache_start)
    if (allocated(neigh_cache_list))  deallocate(neigh_cache_list)
    allocate(neigh_cache_start(mesh%n_vert + 1))

    neigh_cache_start(1) = 1
    do v = 1, mesh%n_vert
      call gather_ls_neighbors(mesh, v, n_neigh, neigh)
      neigh_cache_start(v+1) = neigh_cache_start(v) + n_neigh
      deallocate(neigh)
    end do

    total = neigh_cache_start(mesh%n_vert + 1) - 1
    allocate(neigh_cache_list(total))
    do v = 1, mesh%n_vert
      call gather_ls_neighbors(mesh, v, n_neigh, neigh)
      neigh_cache_list(neigh_cache_start(v):neigh_cache_start(v+1)-1) = neigh
      deallocate(neigh)
    end do

    neigh_cache_n_vert = mesh%n_vert
  end subroutine ensure_neighbor_cache

  ! Element neighbors of id_vert for the least-squares fit: exactly
  ! mesh%vert(id_vert)%elem_neigh, the cells directly touching the vertex --
  ! see the header comment on compute_nodal_derivative_at_vertex for why
  ! this is never ring-expanded.
  subroutine gather_ls_neighbors(mesh, id_vert, n_neigh, neigh)
    use sort_module, only: add_sort_unique_int
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    integer(kind=ENTIER), intent(out) :: n_neigh
    integer(kind=ENTIER), dimension(:), allocatable, intent(out) :: neigh

    n_neigh = 0
    call add_sort_unique_int(neigh, n_neigh, mesh%vert(id_vert)%elem_neigh)
  end subroutine gather_ls_neighbors
end module arbitrary_high_order_module
