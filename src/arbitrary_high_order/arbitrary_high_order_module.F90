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
  public :: apply_grad_bias_correction

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
  ! all).
  real(kind=DOUBLE), parameter :: eps_weno = tiny(1.0_DOUBLE)
  ! Made runtime-settable (2026-09-17, was `parameter`), per the user's
  ! request to sweep OI/epsilon choices without a rebuild per value.
  !
  ! DEFAULT RAISED 1e-6 -> 1e-2 (2026-09-17), per a precise per-level
  ! diagnostic (test_tensor_levels.F90, scratchpad -- compares grad/hess/
  ! third against the ISENTROPIC VORTEX's own analytic derivatives, not
  ! the synthetic cubic field 1e-6 was originally calibrated against).
  ! Root cause: on genuinely smooth (non-polynomial) data, oi_v's
  ! residual term shrinks like O(h^2) as expected, but its actual
  ! magnitude at PRACTICAL mesh resolutions (h down to 0.0125, N=80) sits
  ! around 1e-2 to 1e-4 for this field's real curvature -- nowhere near
  ! small enough for a 1e-6 floor to ever dominate the weight, so the
  ! "WENO degenerates to the linear blend on smooth data" property never
  ! actually engaged at any resolution that's practical to run. Measured
  ! directly: at 1e-6, grad's own order caps at ~1.5 (not 2) and
  ! compounds through hess (~0.7-1.1, degrading with refinement) down to
  ! third (order goes NEGATIVE, N=40->80) -- confirmed NOT caused by
  ! eps_weight_num_deep (swept 1 to 1e4, negligible effect) nor by the
  ! gradient-norm term (grad_norm_derate=1e12, i.e. effectively disabled,
  ! made no difference either). Raising eps_weight_num itself to 1e-2
  ! recovers each level to within a few percent of the pure-linear-blend
  ! reference at the same resolutions (grad order ~2.0-2.04, hess/third
  ! visually overlapping the linear curve), and the full solver's
  ! order-4 vortex convergence rate rises from ~2.1-2.7 to ~3.4-3.9
  ! (N=20->40->80, matched cfl/tmax to Table~\ref{tab:vortex-solver-conv}
  ! in the paper). Cost: Sod's order-4 overshoot grows from
  ! rho_max=1.00048 to 1.00758 (p_max 1.00067->1.01063) -- roughly an
  ! order of magnitude worse in relative terms, but still comparable to
  ! or better than the already-accepted pure-linear-blend overshoot
  ! (rho_max=1.0103, p_max=1.0145) reported elsewhere in this paper, and
  ! nowhere near unstable (no NaN, no qualitative blowup at any order
  ! tested). Net: a real, better-motivated calibration than 1e-6 was,
  ! not merely a different arbitrary constant -- 1e-6 was calibrated
  ! against a field (exactly cubic, so oi_v's residual term is
  ! identically zero away from any true discontinuity) that could not
  ! have exposed this problem in the first place.
  real(kind=DOUBLE), public :: eps_weight_num = 1.0e-2_DOUBLE
  ! Separate, LARGER eps for the hess/third recursion levels only
  ! (deriv_order>=2) -- 2026-09-16 experiment, per the user's request to
  ! keep pushing on the order-4 vs order-3 gap while staying within the
  ! existing single-affine-fit-per-vertex architecture (no quadratic fit,
  ! no wider stencil). Motivation: order 2 (grad only, no hess/third) was
  ! found to regress on Sod as soon as eps_weight_num was raised past
  ! 1e-6, proving the Sod cost comes specifically from the GRAD level's
  ! own de-centering being weakened -- so a single shared eps forces an
  ! all-or-nothing tradeoff between Sod robustness and smooth-field
  ! hess/third accuracy. Splitting them lets hess/third (only exercised
  ! at order>=3, and whose own smooth-data OI values this session found
  ! are already tiny/negligible except right at a genuine discontinuity)
  ! use a more forgiving floor without touching grad's own protection.
  real(kind=DOUBLE), public :: eps_weight_num_deep = 1.0_DOUBLE
  ! Exponent on the oscillation indicator in the WENO weight,
  ! weight = omega_p/(eps+OI_v^weno_power). Standard WENO/Jiang-Shu
  ! practice is power=2; RECALIBRATED to power=1 on 2026-09-16, per the
  ! user's request for more aggressive de-centering on Sod. Since OI_v
  ! (Eq. in the paper, combined residual+gradient-norm indicator) is
  ! dimensionless and typically < 1 even at a discontinuity (raw
  ! oi_gradnorm divided by grad_norm_derate=1e4), RAISING the power
  ! SHRINKS oi_v^power faster than it shrinks eps -- i.e. a higher power
  ! makes the indicator LESS discriminating here, not more, once oi_v^power
  ! drops below eps_weight_num. Confirmed directly: sweeping
  ! power=1,2,4,6,8 on Sod (order 2/3/4, local run) found power=1 gives
  ! by far the smallest pressure overshoot/undershoot (order 3:
  ! p_min=0.099994/p_max=1.000494, vs power=2's 0.099904/1.006215 and
  ! power=4's 0.093734/1.013196 -- degrading monotonically as power
  ! increases from there). Also verified power=1 does NOT regress (in
  ! fact modestly IMPROVES) the smooth stationary-vortex benchmark
  ! (order 3 L2: 1.0511e-2 vs power=2's 1.1233e-2; order 4: 1.2081e-2 vs
  ! 1.2426e-2) and leaves the synthetic-cubic grad/hess convergence
  ! order unchanged (hess still ~order 1, third still flat but smaller
  ! in absolute magnitude). A real, dual-validated improvement, not a
  ! narrow overfit to one case.
  integer(kind=ENTIER), public :: weno_power = 1

  ! Calibration divisor on the gradient-norm term added to oi_v (see
  ! compute_nodal_derivative_at_vertex) -- calibrated 2026-09-16 by an A/B
  ! sweep on the real order-4 isentropic vortex solver (N=80, cfl=0.1,
  ! two_wave), not just the synthetic cubic test. At derate=1 the raw term
  ! regressed the smooth vortex by 26-33% (L2 1.24e-2->1.56e-2 at order 4)
  ! because a REAL flow's gradient magnitude is far from the tiny synthetic
  ! test's -- the term was firing on ordinary smooth variation everywhere,
  ! not just at a discontinuity. Swept 1/1e2/1e4/1e6/1e8: the regression is
  ! fully gone (matches the pre-fix L2 to 4 significant figures) by 1e4 and
  ! stays flat at 1e6, so 1e4 is picked as the smallest value that already
  ! plateaus (no reason to derate further and needlessly shrink the shock
  ! margin below). Confirmed still comfortably detects a Sod-like
  ! axis-aligned jump at this value: that field's own oi_gradnorm/derate is
  ! ~1.5e-4 (N=20) growing to ~1.2e-3 (N=80) -- 150x-1200x above
  ! eps_weight_num=1e-6, a margin that only IMPROVES with mesh refinement
  ! -- while the smooth vortex/cubic field's own oi_gradnorm/derate (~1e-9
  ! at N=40) sits safely BELOW eps, i.e. genuinely negligible there, unlike
  ! at derate=1.
  real(kind=DOUBLE), public :: grad_norm_derate = 1.0e4_DOUBLE

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

  ! Per-cell second geometric moment (M2, about each cell's own centroid),
  ! used by apply_grad_bias_correction -- purely a function of mesh
  ! geometry, so cached once per distinct mesh (same n_elems fingerprint
  ! convention as neigh_cache above) rather than recomputed via quadrature
  ! on every call (every RK stage, every iteration): the un-cached version
  ! measurably slowed down an order-4 run (found 2026-09-16).
  integer(kind=ENTIER), save :: m2_cache_n_elems = -1
  real(kind=DOUBLE), dimension(:), allocatable, save :: m2_cache_xx, m2_cache_xy, m2_cache_yy

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
  subroutine compute_next_order_derivative(mesh, d, nc_in, boundary_2d, phi, dphi, deriv_order, &
      dphi_v_out, valid_v_out, oi_v_out)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(nc_in, mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(nc_in*d, mesh%n_elems), intent(out) :: dphi
    ! Which recursion level this call computes: 1=grad, 2=hess, 3=third --
    ! selects eps_weight_num (level 1) vs eps_weight_num_deep (level>=2)
    ! in accumulate_weno_contribution/scatter_weno_weighted. Optional,
    ! defaults to 1 (the more protective/conservative choice) for any
    ! caller that doesn't specify it.
    integer(kind=ENTIER), intent(in), optional :: deriv_order
    ! Optional: exposes the PER-VERTEX nodal estimate (before WENO
    ! scatter/blend into cells) that this call computes internally anyway
    ! -- used by the order-4 grad bias correction (apply_grad_bias_
    ! correction) to get an accurate per-vertex Hessian (from the hess-
    ! level call) and third-derivative tensor (from the third-level call)
    ! without a second, redundant LS solve.
    real(kind=DOUBLE), dimension(:, :), allocatable, intent(out), optional :: dphi_v_out
    logical, dimension(:), allocatable, intent(out), optional :: valid_v_out
    ! Optional: exposes the PER-VERTEX oscillation indicator oi_v this
    ! call computes internally anyway. Used by apply_grad_bias_correction
    ! (2026-09-17) to recompute the SAME per-vertex WENO weight
    ! (1/(eps_weight_num+oi_v**weno_power), or 1 if use_weno_blend is
    ! .false.) that this call itself used to scatter the gradient into
    ! cells -- see that subroutine's header for why this consistency
    ! matters (the correction's own vertex-to-cell scatter must match the
    ! weighting grad_cell was actually built with, or the two don't
    ! cancel correctly under a non-uniform WENO blend).
    real(kind=DOUBLE), dimension(:), allocatable, intent(out), optional :: oi_v_out

    integer(kind=ENTIER) :: nc_out, id_vert, id_elem, deriv_order_eff
    real(kind=DOUBLE), dimension(:, :), allocatable :: weno_num
    real(kind=DOUBLE), dimension(:), allocatable :: weno_den
    real(kind=DOUBLE), dimension(:, :), allocatable :: dphi_v_cache
    logical, dimension(:), allocatable :: valid_cache
    real(kind=DOUBLE), dimension(:), allocatable :: oi_cache

    nc_out = nc_in*d
    deriv_order_eff = 1
    if (present(deriv_order)) deriv_order_eff = deriv_order

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
        phi, dphi_v_cache, valid_cache, oi_cache, weno_num, weno_den, &
        deriv_order=deriv_order_eff)
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

    if (present(dphi_v_out)) then
      allocate(dphi_v_out(nc_out, mesh%n_vert))
      dphi_v_out = dphi_v_cache
    end if
    if (present(valid_v_out)) then
      allocate(valid_v_out(mesh%n_vert))
      valid_v_out = valid_cache
    end if
    if (present(oi_v_out)) then
      allocate(oi_v_out(mesh%n_vert))
      oi_v_out = oi_cache
    end if

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

    ! BUG FIX (2026-09-15): this was unconditionally nonlinear, ignoring
    ! use_weno_blend entirely -- a linear-mode run (use_weno_blend=.false.)
    ! silently still got the full nonlinear per-vertex weighting here, so
    ! the "linear blend" ablation was never actually linear once CWENO
    ! shipped. Mirrors scatter_weno_weighted's own use_weno_blend gate.
    ! When use_weno_blend=.false., weno_num/weno_den reduce to exactly
    ! lin_num/lin_den (same formula, weight=1), so the center-candidate
    ! blend below is provably a no-op in that mode (algebraically
    ! (X + w_c*(X/Y))/(Y + w_c) = X/Y for any w_c) -- CWENO stays
    ! correct, it just has nothing left to do when there is no nonlinear
    ! candidate to protect against.
    if (use_weno_blend) then
      vertex_weno_weight = 1.0_DOUBLE / (eps_weight_num + oi_v**weno_power)
    else
      vertex_weno_weight = 1.0_DOUBLE
    end if

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
        phi_prev, dfield(order)%val, deriv_order=order)

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
        phi_prev, dfield(order)%val, deriv_order=order)
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
      phi, dphi_v_cache, valid_cache, oi_cache, weno_num, weno_den, skip_cell, deriv_order)
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
    integer(kind=ENTIER), intent(in), optional :: deriv_order

    real(kind=DOUBLE), dimension(nc_in*d) :: dphi_v
    logical :: valid
    real(kind=DOUBLE) :: oi_v, eps_here
    integer(kind=ENTIER) :: deriv_order_eff

    call compute_nodal_derivative_at_vertex(mesh, d, nc_in, boundary_2d, &
      id_vert, phi, dphi_v, valid, oi_v)
    dphi_v_cache(:, id_vert) = dphi_v
    valid_cache(id_vert) = valid
    oi_cache(id_vert) = oi_v
    if (.not. valid) return ! see compute_nodal_derivative_at_vertex's header

    deriv_order_eff = 1
    if (present(deriv_order)) deriv_order_eff = deriv_order
    eps_here = merge(eps_weight_num, eps_weight_num_deep, deriv_order_eff <= 1)

    call scatter_weno_weighted(mesh, id_vert, dphi_v, oi_v, weno_num, weno_den, skip_cell, eps_here)
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
  subroutine scatter_weno_weighted(mesh, id_vert, dphi_v, oi_v, weno_num, weno_den, skip_cell, eps_in)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(:), intent(in) :: dphi_v
    real(kind=DOUBLE), intent(in) :: oi_v
    real(kind=DOUBLE), dimension(:, :), intent(inout) :: weno_num
    real(kind=DOUBLE), dimension(:), intent(inout) :: weno_den
    logical, dimension(:), intent(in), optional :: skip_cell
    ! Overrides eps_weight_num when present -- lets the caller pick a
    ! level-dependent eps (grad vs hess/third), see eps_weight_num_deep's
    ! header and accumulate_weno_contribution.
    real(kind=DOUBLE), intent(in), optional :: eps_in

    integer(kind=ENTIER) :: j, id_elem, id_sub_elem
    real(kind=DOUBLE) :: sub_elem_volume, vertex_weno_weight, eps_use

    eps_use = eps_weight_num
    if (present(eps_in)) eps_use = eps_in

    if (use_weno_blend) then
      vertex_weno_weight = 1.0_DOUBLE / (eps_use + oi_v**weno_power)
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
    !
    ! FIX (2026-09-16, per the user): this residual is BLIND to a
    ! discontinuity that happens to be locally aligned with the mesh (e.g.
    ! a 1D shock like Sod's, or any shock front locally parallel to a grid
    ! line): an interior quad/hex vertex's dual stencil samples only 2
    ! distinct coordinate values per active direction, so a jump that is
    ! constant along the OTHER active direction(s) is fit *exactly* by the
    ! affine model regardless of its size -- confirmed numerically
    ! (test_oi_stats in the aho-oi-redesign investigation): resid-based
    ! OI = O(1e-16) machine noise at a genuine Sod-like 8x jump straddling
    ! a vertex, indistinguishable from smooth data. Combined here with a
    ! second, classical (Jiang-Shu-style) indicator that does NOT share
    ! this blind spot: beta = h_local*|grad|^2/phi_scale^2, using this
    ! same vertex's own just-solved gradient (rhs(2:,:)) and the local
    ! neighbor coordinate spread (max_spread) already computed above. A
    ! smooth, converged gradient gives beta -> 0 as h_local -> 0 (finite
    ! |grad|, shrinking prefactor); a vertex whose stencil straddles a
    ! real jump has an affine-fit slope that blows up like 1/h_local
    ! (fitting a fixed jump over a shrinking baseline), so beta ~
    ! h_local*(1/h_local)^2 = 1/h_local DIVERGES under refinement --
    ! confirmed to catch exactly the Sod-like case the residual alone
    ! missed (beta ~1.2-12, growing with resolution, vs ~1e-5 max on the
    ! smooth cubic verification field), while adding negligible cost
    ! (reuses this vertex's own fit, no extra neighbor pass). Taking the
    ! max of the two keeps whichever indicator is already doing its job:
    ! the residual for genuine multi-directional curvature/oscillation,
    ! the gradient-norm term for an axis-aligned jump the residual can't
    ! see. eps_weight_num/weno_power (both unchanged, classical WENO
    ! values) apply to this combined oi_v exactly as before.
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
    oi_v = max(oi_v, (max_spread * sum(rhs(2:n_basis, :)**2) / max(phi_scale2, 1.0e-300_DOUBLE)) &
      / grad_norm_derate)

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

  ! Builds m2_cache_xx/xy/yy (each cell's own second geometric moment
  ! about its own centroid) once per distinct mesh, via the same
  ! degree-5 quadrature apply_grad_bias_correction used to recompute this
  ! on every call before this cache existed (2026-09-16 perf fix -- see
  ! that cache's own header). A no-op if already built for this mesh.
  subroutine ensure_m2_cache(mesh)
    use quadrature_module, only: volume_quad_pts
    implicit none

    type(mesh_type), intent(in) :: mesh

    integer(kind=ENTIER) :: i, kv, n_v
    real(kind=DOUBLE), dimension(:, :), allocatable :: vcoords, pq
    real(kind=DOUBLE), dimension(:), allocatable :: wq
    real(kind=DOUBLE) :: xc, yc

    if (m2_cache_n_elems == mesh%n_elems) return

    if (allocated(m2_cache_xx)) deallocate(m2_cache_xx)
    if (allocated(m2_cache_xy)) deallocate(m2_cache_xy)
    if (allocated(m2_cache_yy)) deallocate(m2_cache_yy)
    allocate(m2_cache_xx(mesh%n_elems), m2_cache_xy(mesh%n_elems), m2_cache_yy(mesh%n_elems))

    do i = 1, mesh%n_elems
      n_v = mesh%elem(i)%n_vert
      allocate(vcoords(3, n_v))
      do kv = 1, n_v
        vcoords(:, kv) = mesh%vert(mesh%elem(i)%vert(kv))%coord
      end do
      call volume_quad_pts(n_v, vcoords, 5_ENTIER, pq, wq)
      xc = mesh%elem(i)%coord(1); yc = mesh%elem(i)%coord(2)
      m2_cache_xx(i) = sum(wq*(pq(1,:)-xc)**2) / sum(wq)
      m2_cache_xy(i) = sum(wq*(pq(1,:)-xc)*(pq(2,:)-yc)) / sum(wq)
      m2_cache_yy(i) = sum(wq*(pq(2,:)-yc)**2) / sum(wq)
      deallocate(vcoords, pq, wq)
    end do

    m2_cache_n_elems = mesh%n_elems
  end subroutine ensure_m2_cache

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

  ! Corrects the grad step's own O(h^2) bias against a genuinely cubic
  ! field, IN PLACE, using the per-vertex nodal Hessian (Hv, from the
  ! hess-level compute_next_order_derivative call's dphi_v_out) and
  ! per-vertex nodal third-derivative tensor (Tv, from the third-level
  ! call's dphi_v_out) -- both already computed anyway as intermediate
  ! quantities, no extra neighbors or LS solves needed beyond a handful
  ! of quadrature calls for each cell's own second geometric moment.
  !
  ! Derivation (2026-09-16, validated on the synthetic cubic field to
  ! machine precision under a plain/linear blend, and to ~order 4 under
  ! WENO with a sufficiently large eps_weight_num_deep): the vertex's own
  ! weighted-LS affine fit (Step 1, same 4-neighbor stencil) is biased by
  ! THREE distinct, additive contributions, all linear in the (assumed
  ! locally cubic) field and its own H,T:
  !   (1) the vertex's OWN Taylor terms beyond affine (quadratic+cubic),
  !       evaluated exactly at the vertex using its own Hv,Tv;
  !   (2) each NEIGHBOR cell's own cell-average-vs-point-value gap
  !       (phi_j fed to the fit is a genuine cell AVERAGE, not the point
  !       value at that neighbor's centroid the Taylor expansion in (1)
  !       assumes) -- 0.5*H(centroid_j):M2_j, H(centroid_j)=Hv+Tv.dx_j
  !       (exact, H linear for a cubic field), M2_j that neighbor cell's
  !       own second geometric moment (same quantity compute_cell_moments
  !       computes elsewhere in this codebase, just for every INPUT
  !       neighbor here, not only the cell being reconstructed);
  !   (3) a further gap once the (now Taylor-corrected) per-vertex
  !       estimates are themselves blended (volume-weighted average, same
  !       weights as the ordinary grad blend) into a CELL value: since
  !       the true gradient is itself curved across the cell (its own
  !       "Hessian" is exactly Tv), averaging several corner (vertex)
  !       samples does not equal the value at the cell's own centroid --
  !       corrected the same way, using the DISCRETE second moment of the
  !       touching vertices' own positions (not a continuous quadrature
  !       moment) weighted by the same sub_elem_volume blend weights.
  ! All three are exactly zero for genuinely affine or quadratic data (no
  ! bias introduced where none exists); for a cubic field they combine to
  ! reproduce the FULL observed grad bias exactly (confirmed by direct
  ! comparison against the analytic bias at several individual vertices).
  !
  ! Only implemented for boundary_2d=.true. (2D active directions x,y);
  ! a .false. call is a no-op (grad returned unchanged) -- extending to
  ! genuine 3D would need the analogous but algebraically larger 3D
  ! Taylor/moment formulas, not attempted here.
  subroutine apply_grad_bias_correction(mesh, d, nc_in, boundary_2d, grad_cell, &
      hess_v, third_v, valid_hess_v, valid_third_v, grad_oi_v)
    use linear_solver_module, only: lu_factor_lapack, lu_solve_mat_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(nc_in*d, mesh%n_elems), intent(inout) :: grad_cell
    real(kind=DOUBLE), dimension(nc_in*d*d, mesh%n_vert), intent(in) :: hess_v
    real(kind=DOUBLE), dimension(nc_in*d*d*d, mesh%n_vert), intent(in) :: third_v
    logical, dimension(mesh%n_vert), intent(in) :: valid_hess_v, valid_third_v
    ! Per-vertex oscillation indicator from the GRADIENT-level call of
    ! compute_next_order_derivative (its oi_v_out) -- used to recompute
    ! the exact per-vertex weight (sub_elem_volume*vertex_weno_weight)
    ! that grad_cell was actually scattered with, so the two
    ! vertex-to-cell scatters below (bias and cell-blend-curvature) match
    ! how grad_cell was built instead of silently assuming a plain
    ! sub_elem_volume (uniform) blend. Under use_weno_blend=.false. this
    ! reduces to weight=1 exactly as before (bit-identical to the
    ! pre-2026-09-17 code); under WENO it is what closes the gap between
    ! this correction reaching machine precision (linear blend) and only
    ! ~order 2.7 (WENO blend) on the same smooth field -- see
    ! sec:results-order4's "Fix, attempt 3" discussion.
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: grad_oi_v

    integer(kind=ENTIER) :: iv, j, id_elem, id_sub_elem, kv, n_v, ic, i
    integer(kind=ENTIER) :: n_neigh, n_basis
    integer(kind=ENTIER), dimension(:), allocatable :: neigh
    real(kind=DOUBLE), dimension(3, 3) :: mat
    real(kind=DOUBLE), dimension(3, nc_in) :: rhs
    real(kind=DOUBLE), dimension(3) :: basis
    integer(kind=ENTIER), dimension(3) :: ipiv
    real(kind=DOUBLE), dimension(3) :: dx_v
    real(kind=DOUBLE) :: weight, dxx, dyy, sub_elem_volume, dvx, dvy, vweight
    real(kind=DOUBLE), dimension(:), allocatable :: Hxx, Hxy, Hyy, Txxx, Txxy, Txyy, Tyyy
    real(kind=DOUBLE), dimension(:), allocatable :: HxxJ, HxyJ, HyyJ, moment_j
    real(kind=DOUBLE), dimension(:, :), allocatable :: bias_num_x, bias_num_y
    real(kind=DOUBLE), dimension(:), allocatable :: bias_den

    if (.not. boundary_2d) return ! 3D not implemented, see header

    call ensure_m2_cache(mesh)
    call ensure_neighbor_cache(mesh)

    allocate(Hxx(nc_in), Hxy(nc_in), Hyy(nc_in), Txxx(nc_in), Txxy(nc_in), Txyy(nc_in), Tyyy(nc_in))
    allocate(HxxJ(nc_in), HxyJ(nc_in), HyyJ(nc_in), moment_j(nc_in))
    allocate(bias_num_x(nc_in, mesh%n_elems), bias_num_y(nc_in, mesh%n_elems), bias_den(mesh%n_elems))
    bias_num_x = 0.0_DOUBLE; bias_num_y = 0.0_DOUBLE; bias_den = 0.0_DOUBLE

    n_basis = 3 ! 1, dx, dy -- SAME affine basis as the Step-1 fit
    do iv = 1, mesh%n_vert
      if (mesh%vert(iv)%is_bound) cycle
      if (.not. valid_hess_v(iv) .or. .not. valid_third_v(iv)) cycle

      do ic = 1, nc_in
        Hxx(ic) = hess_v(ic, iv)                                    ! (dir1,dir2)=(1,1)
        Hxy(ic) = 0.5_DOUBLE*(hess_v(2*nc_in+ic, iv) + hess_v(nc_in+ic, iv)) ! (1,2)&(2,1) avg
        Hyy(ic) = hess_v(3*nc_in+ic, iv)                             ! (2,2)  [nc_in*d=3*nc_in]
        Txxx(ic) = third_v(ic, iv)                                   ! (1,1,1)
        Txxy(ic) = (third_v(3*nc_in+ic, iv) + third_v(9*nc_in+ic, iv) &
          + third_v(nc_in+ic, iv)) / 3.0_DOUBLE                      ! (1,1,2)&(1,2,1)&(2,1,1)
        Txyy(ic) = (third_v(4*nc_in+ic, iv) + third_v(10*nc_in+ic, iv) &
          + third_v(12*nc_in+ic, iv)) / 3.0_DOUBLE                   ! (2,2,1)&(2,1,2)&(1,2,2)
        Tyyy(ic) = third_v(13*nc_in+ic, iv)                          ! (2,2,2) [nc_in*d*d=9*nc_in]
      end do

      n_neigh = neigh_cache_start(iv+1) - neigh_cache_start(iv)
      allocate(neigh(n_neigh))
      neigh = neigh_cache_list(neigh_cache_start(iv):neigh_cache_start(iv+1)-1)
      mat = 0.0_DOUBLE; rhs = 0.0_DOUBLE
      do j = 1, n_neigh
        id_elem = neigh(j)
        dx_v = mesh%elem(id_elem)%coord - mesh%vert(iv)%coord
        weight = 1.0_DOUBLE / max(dot_product(dx_v, dx_v), 1.0e-24_DOUBLE)
        dxx = dx_v(1); dyy = dx_v(2)
        basis(1) = 1.0_DOUBLE; basis(2) = dxx; basis(3) = dyy
        mat = mat + weight * spread(basis, 2, 3) * spread(basis, 1, 3)

        ! (1) vertex's own Taylor terms beyond affine
        moment_j = 0.5_DOUBLE*(Hxx*dxx**2 + 2.0_DOUBLE*Hxy*dxx*dyy + Hyy*dyy**2) &
          + (Txxx*dxx**3 + 3.0_DOUBLE*Txxy*dxx**2*dyy + 3.0_DOUBLE*Txyy*dxx*dyy**2 + Tyyy*dyy**3)/6.0_DOUBLE

        ! (2) this neighbor's own cell-average-vs-point-value gap
        HxxJ = Hxx + Txxx*dxx + Txxy*dyy
        HxyJ = Hxy + Txxy*dxx + Txyy*dyy
        HyyJ = Hyy + Txyy*dxx + Tyyy*dyy
        moment_j = moment_j + 0.5_DOUBLE*(HxxJ*m2_cache_xx(id_elem) &
          + 2.0_DOUBLE*HxyJ*m2_cache_xy(id_elem) + HyyJ*m2_cache_yy(id_elem))

        do ic = 1, nc_in
          rhs(:, ic) = rhs(:, ic) + weight*basis*moment_j(ic)
        end do
      end do
      deallocate(neigh)

      call lu_factor_lapack(n_basis, mat, ipiv)
      call lu_solve_mat_lapack(n_basis, mat, ipiv, nc_in, rhs)

      ! Match grad_cell's OWN vertex-to-cell weighting exactly (see this
      ! subroutine's header on grad_oi_v) instead of assuming a plain
      ! sub_elem_volume blend.
      if (use_weno_blend) then
        vweight = 1.0_DOUBLE / (eps_weight_num + grad_oi_v(iv)**weno_power)
      else
        vweight = 1.0_DOUBLE
      end if

      do j = 1, mesh%vert(iv)%n_sub_elems_neigh
        id_sub_elem = mesh%vert(iv)%sub_elem_neigh(j)
        id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
        sub_elem_volume = mesh%sub_elem(id_sub_elem)%volume * vweight
        bias_num_x(:, id_elem) = bias_num_x(:, id_elem) + sub_elem_volume*rhs(2, :)
        bias_num_y(:, id_elem) = bias_num_y(:, id_elem) + sub_elem_volume*rhs(3, :)
        bias_den(id_elem) = bias_den(id_elem) + sub_elem_volume
      end do
    end do

    do i = 1, mesh%n_elems
      if (bias_den(i) <= 0.0_DOUBLE) cycle
      grad_cell(1:nc_in, i) = grad_cell(1:nc_in, i) - bias_num_x(:, i)/bias_den(i)
      grad_cell(nc_in+1:2*nc_in, i) = grad_cell(nc_in+1:2*nc_in, i) - bias_num_y(:, i)/bias_den(i)
    end do

    ! (3) cell-blend curvature: even with each vertex's own grad exactly
    ! unbiased, blending several corner samples of a CURVED (quadratic)
    ! gradient field does not equal its value at the cell centroid --
    ! correct using the touching vertices' own Tv, blended the same way,
    ! and the DISCRETE second moment of their positions about the cell
    ! centroid (not a quadrature moment).
    block
      real(kind=DOUBLE), dimension(:, :), allocatable :: m2d_num_xx, m2d_num_xy, m2d_num_yy
      real(kind=DOUBLE), dimension(:), allocatable :: m2d_den
      real(kind=DOUBLE), dimension(:, :), allocatable :: t_num_xxx, t_num_xxy, t_num_xyy, t_num_yyy
      allocate(m2d_num_xx(1, mesh%n_elems), m2d_num_xy(1, mesh%n_elems), m2d_num_yy(1, mesh%n_elems))
      allocate(m2d_den(mesh%n_elems))
      allocate(t_num_xxx(nc_in, mesh%n_elems), t_num_xxy(nc_in, mesh%n_elems))
      allocate(t_num_xyy(nc_in, mesh%n_elems), t_num_yyy(nc_in, mesh%n_elems))
      m2d_num_xx = 0.0_DOUBLE; m2d_num_xy = 0.0_DOUBLE; m2d_num_yy = 0.0_DOUBLE; m2d_den = 0.0_DOUBLE
      t_num_xxx = 0.0_DOUBLE; t_num_xxy = 0.0_DOUBLE; t_num_xyy = 0.0_DOUBLE; t_num_yyy = 0.0_DOUBLE
      do iv = 1, mesh%n_vert
        if (mesh%vert(iv)%is_bound) cycle
        if (.not. valid_third_v(iv)) cycle
        do ic = 1, nc_in
          Txxx(ic) = third_v(ic, iv)
          Txxy(ic) = (third_v(3*nc_in+ic, iv) + third_v(9*nc_in+ic, iv) + third_v(nc_in+ic, iv)) / 3.0_DOUBLE
          Txyy(ic) = (third_v(4*nc_in+ic, iv) + third_v(10*nc_in+ic, iv) + third_v(12*nc_in+ic, iv)) / 3.0_DOUBLE
          Tyyy(ic) = third_v(13*nc_in+ic, iv)
        end do
        ! Same grad_cell-consistent weight as the bias scatter above --
        ! this loop redistributes grad's OWN per-vertex samples, so it
        ! must match grad's own blend too, not third_v's.
        if (use_weno_blend) then
          vweight = 1.0_DOUBLE / (eps_weight_num + grad_oi_v(iv)**weno_power)
        else
          vweight = 1.0_DOUBLE
        end if

        do j = 1, mesh%vert(iv)%n_sub_elems_neigh
          id_sub_elem = mesh%vert(iv)%sub_elem_neigh(j)
          id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
          sub_elem_volume = mesh%sub_elem(id_sub_elem)%volume * vweight
          dvx = mesh%vert(iv)%coord(1) - mesh%elem(id_elem)%coord(1)
          dvy = mesh%vert(iv)%coord(2) - mesh%elem(id_elem)%coord(2)
          m2d_num_xx(1,id_elem) = m2d_num_xx(1,id_elem) + sub_elem_volume*dvx*dvx
          m2d_num_xy(1,id_elem) = m2d_num_xy(1,id_elem) + sub_elem_volume*dvx*dvy
          m2d_num_yy(1,id_elem) = m2d_num_yy(1,id_elem) + sub_elem_volume*dvy*dvy
          m2d_den(id_elem) = m2d_den(id_elem) + sub_elem_volume
          t_num_xxx(:,id_elem) = t_num_xxx(:,id_elem) + sub_elem_volume*Txxx
          t_num_xxy(:,id_elem) = t_num_xxy(:,id_elem) + sub_elem_volume*Txxy
          t_num_xyy(:,id_elem) = t_num_xyy(:,id_elem) + sub_elem_volume*Txyy
          t_num_yyy(:,id_elem) = t_num_yyy(:,id_elem) + sub_elem_volume*Tyyy
        end do
      end do
      do i = 1, mesh%n_elems
        if (m2d_den(i) <= 0.0_DOUBLE) cycle
        block
          real(kind=DOUBLE) :: M2xx, M2xy, M2yy
          real(kind=DOUBLE), dimension(nc_in) :: Tx1, Tx2, Tx3, Tx4, extra_x, extra_y
          M2xx = m2d_num_xx(1,i)/m2d_den(i); M2xy = m2d_num_xy(1,i)/m2d_den(i); M2yy = m2d_num_yy(1,i)/m2d_den(i)
          Tx1 = t_num_xxx(:,i)/m2d_den(i); Tx2 = t_num_xxy(:,i)/m2d_den(i)
          Tx3 = t_num_xyy(:,i)/m2d_den(i); Tx4 = t_num_yyy(:,i)/m2d_den(i)
          extra_x = 0.5_DOUBLE*(Tx1*M2xx + 2.0_DOUBLE*Tx2*M2xy + Tx3*M2yy)
          extra_y = 0.5_DOUBLE*(Tx2*M2xx + 2.0_DOUBLE*Tx3*M2xy + Tx4*M2yy)
          grad_cell(1:nc_in, i) = grad_cell(1:nc_in, i) - extra_x
          grad_cell(nc_in+1:2*nc_in, i) = grad_cell(nc_in+1:2*nc_in, i) - extra_y
        end block
      end do
    end block

    deallocate(Hxx, Hxy, Hyy, Txxx, Txxy, Txyy, Tyyy, HxxJ, HxyJ, HyyJ, moment_j)
    deallocate(bias_num_x, bias_num_y, bias_den)
  end subroutine apply_grad_bias_correction
end module arbitrary_high_order_module
