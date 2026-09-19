! Arbitrary-order cell derivatives on unstructured meshes: a recursive multi-D
! generalization of 1D divided differences. Given the order-(k-1) derivative
! tensor in every cell (order 0 = the field itself), one order-k step
! (compute_next_order_derivative) is:
!   1) at every vertex, a weighted-least-squares (or Green-Gauss) affine fit
!      of the order-(k-1) tensor over the elements touching that vertex gives
!      the order-k tensor there (compute_nodal_derivative_at_vertex[_green_gauss]);
!   2) that vertex tensor is scattered into every cell around the vertex as a
!      nonlinear-WENO-weighted contribution, weight 1/(eps+OI_v^p);
!   3) each cell normalizes its accumulator by its weight sum.
! Only the order-k tensor needs a ghost exchange before serving as the
! order-(k-1) input of the next step (one layer deep) -- compute_next_order_
! derivative performs a single step and leaves the exchange to the caller, so
! communication of order k can overlap with unrelated work while order k+1 is
! prepared; see compute_derivative_hierarchy for the blocking chain of all orders.
module arbitrary_high_order_module
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
  public :: test_green_gauss_nodal
  public :: test_green_gauss_vortex

  ! eps_weno floors a WENO indicator ratio against literal division by zero.
  real(kind=DOUBLE), parameter :: eps_weno = tiny(1.0_DOUBLE)
  ! Regularizes oi_v (the dimensionless per-vertex WENO oscillation indicator)
  ! in scatter_weno_weighted; eps_weight_num_deep is a separate, more forgiving
  ! floor for the hess/third recursion levels. LS-path values; GG has its own
  ! decoupled eps_weight_num_gg/eps_weight_num_deep_gg below.
  real(kind=DOUBLE), public :: eps_weight_num = 1.0e-2_DOUBLE
  real(kind=DOUBLE), public :: eps_weight_num_deep = 1.0_DOUBLE
  real(kind=DOUBLE), public :: eps_weight_num_gg = 1.0e-2_DOUBLE
  real(kind=DOUBLE), public :: eps_weight_num_deep_gg = 1.0_DOUBLE
  ! Exponent on the oscillation indicator in the WENO weight, weight =
  ! omega_p/(eps+OI_v^weno_power).
  integer(kind=ENTIER), public :: weno_power = 1
  ! Calibration divisor on the gradient-norm term added to oi_v (LS) / oi_v's
  ! only term (GG) -- see compute_nodal_derivative_at_vertex[_green_gauss].
  real(kind=DOUBLE), public :: grad_norm_derate = 1.0e4_DOUBLE
  real(kind=DOUBLE), public :: grad_norm_derate_gg = 1.0e4_DOUBLE
  ! When .true. AND use_green_gauss, scatter_weno_weighted/rescue_zero_weight_
  ! cell use omega_p/(eps+(1+OI)^p) instead of omega_p/(eps+OI^p).
  logical, public :: use_alt_gg_weight = .false.
  ! .false.: nonlinear WENO weight is replaced by a plain volume weight
  ! (weight=1), i.e. a straight sub_elem_volume-weighted average of nodal
  ! derivatives with no oscillation-adaptive de-centering.
  logical, public :: use_weno_blend = .true.
  ! Selects the per-vertex nodal fit compute_next_order_derivative uses at
  ! every recursion level: .false. (default) is the weighted least-squares
  ! fit ("aho ls"); .true. is the Green-Gauss divergence-theorem fit
  ! ("aho gg"). Both produce the same (dphi_v, oi_v) shape and are consumed
  ! identically downstream (WENO scatter, apply_grad_bias_correction).
  logical, public :: use_green_gauss = .false.
  ! CWENO-style central-candidate experiment (Semplice & Visconti 2020):
  ! blends one extra "optimal" candidate (the plain linear/volume-weighted
  ! average) into the same nonlinear WENO sum, weight cweno_center_weight/
  ! (eps+OI_c**cweno_center_power). Off by default; opt in via
  ! compute_next_order_derivative_cweno explicitly.
  logical, public :: use_cweno_center = .false.
  real(kind=DOUBLE), public :: cweno_center_weight = 1000.0_DOUBLE
  integer(kind=ENTIER), public :: cweno_center_power = 4

  ! Per-vertex neighbor-list cache (CSR), built once per distinct mesh:
  ! gather_ls_neighbors' ring-expansion is pure topology, independent of phi.
  integer(kind=ENTIER), save :: neigh_cache_n_vert = -1
  integer(kind=ENTIER), dimension(:), allocatable, save :: neigh_cache_start
  integer(kind=ENTIER), dimension(:), allocatable, save :: neigh_cache_list

  ! Per-cell second geometric moment (M2, about each cell's own centroid),
  ! used by apply_grad_bias_correction: pure mesh geometry, cached once.
  integer(kind=ENTIER), save :: m2_cache_n_elems = -1
  real(kind=DOUBLE), dimension(:), allocatable, save :: m2_cache_xx, m2_cache_xy, m2_cache_yy

  ! Per-vertex Green-Gauss fit matrix, ALREADY INVERTED: pure mesh geometry,
  ! identical at every recursion level and every call for a given vertex, so
  ! cached once per distinct mesh instead of rebuilt+inverted from scratch
  ! every call. boundary_2d is baked in (row/column 3 zeroed).
  integer(kind=ENTIER), save :: gg_mat_cache_n_vert = -1
  logical, save :: gg_mat_cache_boundary_2d = .false.
  real(kind=DOUBLE), dimension(:, :, :), allocatable, save :: gg_mat_inv_cache
  logical, dimension(:), allocatable, save :: gg_mat_valid_cache

  ! Flattened (CSR-style) per-vertex list of the (area*normal, source-cell)
  ! pairs walked by compute_nodal_derivative_at_vertex_green_gauss's flux
  ! sum: also pure geometry, built alongside gg_mat_inv_cache above.
  ! gg_flux_offset_cache(v):gg_flux_offset_cache(v+1)-1 indexes into the two
  ! flat arrays for vertex v. gg_oi_hlocal_cache is the OI indicator's
  ! geometric h_local term, cached the same way.
  integer(kind=ENTIER), dimension(:), allocatable, save :: gg_flux_offset_cache
  real(kind=DOUBLE), dimension(:, :), allocatable, save :: gg_flux_w_cache
  integer(kind=ENTIER), dimension(:), allocatable, save :: gg_flux_elem_cache
  real(kind=DOUBLE), dimension(:), allocatable, save :: gg_oi_hlocal_cache

  type :: derivative_field_type
    real(kind=DOUBLE), dimension(:, :), allocatable :: val ! (d**order, n_elems)
  end type derivative_field_type

  ! One phase of one order's work, timestamped by the *_timed drivers below,
  ! for the blocking-vs-overlap timeline figure in tex_arbitrary_high_order.
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
  ! single ghost layer. boundary_2d must match the mesh's own build flag: on
  ! a mesh extruded as one thin z-layer, every vertex's neighbors share the
  ! same z, which makes the full 3D fit singular -- boundary_2d=.true. drops
  ! z from the fit basis instead.
  subroutine compute_next_order_derivative(mesh, d, nc_in, boundary_2d, phi, dphi, deriv_order, &
      dphi_v_out, valid_v_out, oi_v_out)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(nc_in, mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(nc_in*d, mesh%n_elems), intent(out) :: dphi
    ! Recursion level this call computes: 1=grad, 2=hess, 3=third -- selects
    ! eps_weight_num (level 1) vs eps_weight_num_deep (level>=2). Optional,
    ! defaults to 1.
    integer(kind=ENTIER), intent(in), optional :: deriv_order
    ! Optional: exposes the per-vertex nodal estimate this call computes
    ! internally anyway -- used by apply_grad_bias_correction to get an
    ! accurate per-vertex Hessian/third tensor without a redundant LS solve.
    real(kind=DOUBLE), dimension(:, :), allocatable, intent(out), optional :: dphi_v_out
    logical, dimension(:), allocatable, intent(out), optional :: valid_v_out
    ! Optional: exposes the per-vertex oscillation indicator oi_v, used by
    ! apply_grad_bias_correction to recompute the same per-vertex WENO
    ! weight grad_cell was actually scattered with.
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

    ! Each vertex's nodal derivative is computed once and immediately
    ! scattered, WENO-weighted, into every cell touching it. A vertex on the
    ! physical domain boundary is skipped (its dual-cell neighbor gather is
    ! one-sided, spurious tangential derivatives) unless .not. boundary_2d,
    ! where skipping it too would starve every vertex at once.
    do id_vert = 1, mesh%n_vert
      if (mesh%vert(id_vert)%is_bound) cycle
      call accumulate_weno_contribution(mesh, d, nc_in, boundary_2d, id_vert, &
        phi, dphi_v_cache, valid_cache, oi_cache, weno_num, weno_den, &
        deriv_order=deriv_order_eff)
    end do

    ! Fallback for a cell left with zero total weight (e.g. every touching
    ! vertex sits on the domain boundary): re-admit boundary vertices for
    ! that cell only.
    do id_elem = 1, mesh%n_elems
      if (weno_den(id_elem) == 0.0_DOUBLE) then
        call rescue_zero_weight_cell(mesh, d, nc_in, boundary_2d, id_elem, &
          phi, weno_num, weno_den)
      end if
    end do

    do id_elem = 1, mesh%n_elems
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

  ! CWENO variant of compute_next_order_derivative: accumulates, in one pass,
  ! the usual nonlinear WENO num/den, a LINEAR (unweighted) num/den from the
  ! same per-vertex candidates, and a volume-weighted running average of the
  ! vertices' own OI. The final per-cell derivative folds the linear
  ! candidate back into the nonlinear blend with weight cweno_center_weight/
  ! (eps+OI_c**power) -- see use_cweno_center's header.
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
        ! Rescue's fallback has no linear/OI counterpart -- mirror it into
        ! the linear accumulators too.
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
        w_c = cweno_center_weight / (eps_weight_num + oi_c**cweno_center_power)
        dphi_lin_elem = lin_num(:, id_elem) / lin_den(id_elem)
        dphi(:, id_elem) = (weno_num(:, id_elem) + w_c * dphi_lin_elem) / (weno_den(id_elem) + w_c)
      else
        dphi(:, id_elem) = weno_num(:, id_elem) / weno_den(id_elem)
      end if
    end do

    deallocate(weno_num, weno_den, lin_num, lin_den, oi_num, oi_den)
    deallocate(dphi_v_cache, valid_cache, oi_cache, dphi_lin_elem)
  end subroutine compute_next_order_derivative_cweno

  ! Same per-vertex LS/GG solve as accumulate_weno_contribution, but scatters
  ! into three running sums at once (nonlinear WENO, linear/center, and a
  ! volume-weighted running OI average).
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

    if (use_green_gauss) then
      call compute_nodal_derivative_at_vertex_green_gauss(mesh, d, nc_in, boundary_2d, &
        id_vert, phi, dphi_v, valid, oi_v)
    else
      call compute_nodal_derivative_at_vertex(mesh, d, nc_in, boundary_2d, &
        id_vert, phi, dphi_v, valid, oi_v)
    end if
    dphi_v_cache(:, id_vert) = dphi_v
    valid_cache(id_vert) = valid
    oi_cache(id_vert) = oi_v
    if (.not. valid) return

    ! use_weno_blend=.false. must reduce weno_num/weno_den to exactly
    ! lin_num/lin_den (weight=1), so the center-candidate blend below is a
    ! provable no-op in that mode -- mirrors scatter_weno_weighted's own gate.
    if (use_weno_blend) then
      if (use_green_gauss) then
        if (use_alt_gg_weight) then
          vertex_weno_weight = 1.0_DOUBLE / (eps_weight_num_gg + (1.0_DOUBLE + oi_v)**weno_power)
        else
          vertex_weno_weight = 1.0_DOUBLE / (eps_weight_num_gg + oi_v**weno_power)
        end if
      else
        vertex_weno_weight = 1.0_DOUBLE / (eps_weight_num + oi_v**weno_power)
      end if
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

  ! Rescue path for a cell left with zero total WENO weight after is_bound
  ! vertices were skipped -- only possible if every one of the cell's own
  ! vertices sits on the domain boundary. Re-admits just this cell's own
  ! boundary vertices so the caller never divides by zero.
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
      ! Under use_green_gauss this is currently a no-op for a fully-boundary
      ! cell: the GG fit always returns valid=.false. at a boundary vertex
      ! (no one-sided estimate implemented for it), so such a cell is
      ! correctly left at weno_den=0 rather than silently borrowing an LS estimate.
      if (use_green_gauss) then
        call compute_nodal_derivative_at_vertex_green_gauss(mesh, d, nc_in, boundary_2d, &
          id_vert, phi, dphi_v, valid, oi_v)
      else
        call compute_nodal_derivative_at_vertex(mesh, d, nc_in, boundary_2d, &
          id_vert, phi, dphi_v, valid, oi_v)
      end if
      if (.not. valid) cycle
      sub_elem_volume = mesh%sub_elem(id_sub_elem)%volume
      if (use_weno_blend) then
        if (use_green_gauss) then
          if (use_alt_gg_weight) then
            vertex_weno_weight = 1.0_DOUBLE / (eps_weight_num_gg + (1.0_DOUBLE + oi_v)**weno_power)
          else
            vertex_weno_weight = 1.0_DOUBLE / (eps_weight_num_gg + oi_v**weno_power)
          end if
        else
          vertex_weno_weight = 1.0_DOUBLE / (eps_weight_num + oi_v**weno_power)
        end if
      else
        vertex_weno_weight = 1.0_DOUBLE
      end if
      weno_num(:, id_elem) = weno_num(:, id_elem) &
        + (sub_elem_volume * vertex_weno_weight) * dphi_v
      weno_den(id_elem) = weno_den(id_elem) + sub_elem_volume * vertex_weno_weight
    end do
    ! If every vertex was invalid too (isolated cell cluster), leave
    ! weno_den at 0; the caller must not divide by it.
  end subroutine rescue_zero_weight_cell

  ! Same order-k step as compute_next_order_derivative, but overlapping the
  ! ghost exchange of its own output with local work: (1) sweep vertices
  ! touching a to-be-sent cell and finalize just those cells; (2) post that
  ! exchange non-blockingly; (3) sweep every remaining vertex/cell while the
  ! exchange is in flight; (4) wait for it, filling in this rank's ghost
  ! cells. With num_procs=1 this degenerates to compute_next_order_derivative.
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

    ! Step 1: vertices touching a to-be-sent cell.
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

    ! Step 2: hand the just-finalized boundary cells to MPI and move on.
    if (num_procs > 1) call mpi_memory_exchange_post(mpi_send_recv, mesh%n_elems, nc_out, dphi)

    ! Step 3: everything else, computed while the exchange is in flight.
    ! Boundary vertices reuse their Step-1 cached dphi_v.
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

  ! Reconstructs every order 1..max_order from the scalar field phi0. Each
  ! order's cell-centered input is exchanged across the single ghost layer
  ! before it is used to build the next order.
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

  ! Same as compute_derivative_hierarchy, using compute_next_order_derivative_
  ! overlap for each order: that routine already waits for its own output's
  ! exchange, so phi_prev is ghost-complete for the next order without a
  ! separate top-of-loop exchange -- except phi0 itself, exchanged once up front.
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

  ! Same as compute_derivative_hierarchy, timestamping the two phases of
  ! every order (blocking exchange, then compute) into events(1:n_events) --
  ! for the blocking-vs-overlap timeline figure only. events must be
  ! preallocated by the caller to at least 2*max_order.
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

  ! Same as compute_derivative_hierarchy_overlap, timestamping the four
  ! phases of every order (compute boundary cells, post, compute interior,
  ! wait) -- for the timeline figure only. events must be preallocated to
  ! at least 4*max_order.
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
  ! same vertex, used by the MPI-overlap driver), and immediately scatters
  ! it into every cell touching id_vert as a nonlinear-WENO-weighted
  ! contribution: weight = 1/(eps+OI^p), OI a mesh-size-scaled oscillation
  ! indicator (see scatter_weno_weighted).
  ! skip_cell, if present, excludes cells where skip_cell(id_elem) is .true.
  ! -- used by compute_next_order_derivative_overlap's two-pass split.
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

    if (use_green_gauss) then
      call compute_nodal_derivative_at_vertex_green_gauss(mesh, d, nc_in, boundary_2d, &
        id_vert, phi, dphi_v, valid, oi_v)
    else
      call compute_nodal_derivative_at_vertex(mesh, d, nc_in, boundary_2d, &
        id_vert, phi, dphi_v, valid, oi_v)
    end if
    dphi_v_cache(:, id_vert) = dphi_v
    valid_cache(id_vert) = valid
    oi_cache(id_vert) = oi_v
    if (.not. valid) return

    deriv_order_eff = 1
    if (present(deriv_order)) deriv_order_eff = deriv_order
    if (use_green_gauss) then
      eps_here = merge(eps_weight_num_gg, eps_weight_num_deep_gg, deriv_order_eff <= 1)
    else
      eps_here = merge(eps_weight_num, eps_weight_num_deep, deriv_order_eff <= 1)
    end if

    call scatter_weno_weighted(mesh, id_vert, dphi_v, oi_v, weno_num, weno_den, skip_cell, eps_here)
  end subroutine accumulate_weno_contribution

  ! Same scatter as accumulate_weno_contribution, but reusing an
  ! already-cached dphi_v instead of recomputing it -- used by the MPI
  ! overlap driver's Step 3 to fold a boundary vertex's Step-1 contribution
  ! into non-send cells it also touches, without a second LAPACK solve.
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

  ! Shared scatter: weight = omega_p/(eps+OI^p), omega_p=sub_elem_volume (the
  ! geometric/"linear weight" role, as in classical CWENO -- Semplice &
  ! Visconti 2020, Def. 3) and OI=oi_v, the vertex's own fit residual/
  ! gradient-norm indicator (see compute_nodal_derivative_at_vertex). OI
  ! must measure disagreement with a smooth local model (vanishing under
  ! refinement for smooth data), not raw derivative magnitude -- see that
  ! routine's own header for the indicator's derivation.
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
    ! level-dependent eps (grad vs hess/third).
    real(kind=DOUBLE), intent(in), optional :: eps_in

    integer(kind=ENTIER) :: j, id_elem, id_sub_elem
    real(kind=DOUBLE) :: sub_elem_volume, vertex_weno_weight, eps_use

    eps_use = eps_weight_num
    if (present(eps_in)) eps_use = eps_in

    if (use_weno_blend) then
      ! GG-only alternative: omega_p/(eps+(1+OI)^p) bounds the smooth-data
      ! weight near 1/(eps+1) instead of letting it blow up to 1/eps as OI->0.
      if (use_green_gauss .and. use_alt_gg_weight) then
        vertex_weno_weight = 1.0_DOUBLE / (eps_use + (1.0_DOUBLE + oi_v)**weno_power)
      else
        vertex_weno_weight = 1.0_DOUBLE / (eps_use + oi_v**weno_power)
      end if
    else
      ! weight=1 cancels out of numerator/denominator, leaving a plain
      ! sub_elem_volume-weighted average of nodal derivatives.
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

  ! Weighted-least-squares gradient of phi at id_vert, over mesh%vert(id_vert)%
  ! elem_neigh (never ring-expanded, see below): fits phi_c(x) ~= a0_c + G_c.
  ! (x-x_vert) for each of the nc_in components independently, sharing one
  ! (n_basis,n_basis) geometry-only normal matrix (one LU factorization per
  ! vertex, reused for nc_in right-hand sides).
  ! Candidate directions are x,y (boundary_2d) or x,y,z; only ones with a
  ! genuinely resolvable spread among the gathered neighbors enter the basis
  ! (dropped directions get an exact 0 derivative) -- handles both
  ! boundary_2d's z and a quasi-1D mesh's degenerate y/z. Exact (zero
  ! residual) whenever phi is locally affine in the active directions,
  ! regardless of mesh irregularity.
  ! The stencil is always exactly elem_neigh, never ring-expanded: a wider
  ! ring is not guaranteed complete for a vertex on an MPI partition seam
  ! (only a single node-based ghost layer), which previously caused a
  ! rank-count-dependent accuracy regression.
  subroutine compute_nodal_derivative_at_vertex(mesh, d, nc_in, boundary_2d, &
      id_vert, phi, dphi_v, valid, oi_v)
    use linear_solver_module, only: lu_factor_lapack, lu_solve_mat_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in, id_vert
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(nc_in, mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(nc_in*d), intent(out) :: dphi_v
    ! .false. when not even one direction could be resolved: the caller must
    ! exclude such a vertex from the WENO scatter entirely, not scatter
    ! dphi_v=0 (which the 1/(eps+OI^p) weight would read as perfect confidence).
    logical, intent(out) :: valid
    ! WENO oscillation indicator: the fit's own weighted-mean-squared
    ! residual against the neighbor data, maximized over components and made
    ! dimensionless by the neighbors' own weighted mean-square phi, combined
    ! (via max) with a gradient-norm term beta=h_local*|grad|^2/phi_scale^2
    ! that catches an axis-aligned jump the residual alone is blind to (a
    ! jump exactly fit by the affine model when it doesn't vary along the
    ! stencil's other active directions).
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
    if (.not. valid) return

    call lu_factor_lapack(n_basis, mat(1:n_basis, 1:n_basis), ipiv(1:n_basis))
    call lu_solve_mat_lapack(n_basis, mat(1:n_basis, 1:n_basis), &
      ipiv(1:n_basis), nc_in, rhs(1:n_basis, :))

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
    phi_scale2 = maxval(phi_sq_sum) / max(weight_sum, 1.0e-300_DOUBLE)
    oi_v = sqrt((maxval(resid_sq) / max(weight_sum, 1.0e-300_DOUBLE)) &
      / max(phi_scale2, 1.0e-300_DOUBLE))
    oi_v = max(oi_v, (max_spread * sum(rhs(2:n_basis, :)**2) / max(phi_scale2, 1.0e-300_DOUBLE)) &
      / grad_norm_derate)

    ! rhs(1+a,i1) holds d(phi_i1)/dx_{active_dim(a)}; flatten with i1 (the
    ! carried component) fast-varying and direction slow-varying, matching
    ! the rest of the module's tensor layout. Inactive directions stay 0.
    do a = 1, n_active
      do i1 = 1, nc_in
        dphi_v((active_dim(a)-1)*nc_in + i1) = rhs(1+a, i1)
      end do
    end do
  end subroutine compute_nodal_derivative_at_vertex

  ! Builds gg_mat_inv_cache/gg_flux_*_cache/gg_oi_hlocal_cache (all pure mesh
  ! geometry) once per distinct mesh/boundary_2d combination -- a no-op if
  ! already built. Uses pseudo_inverse_inplace_lapack (SVD) unconditionally:
  ! this is now a one-time-per-vertex cost, so the inversion method no longer
  ! matters for performance, and the SVD's robustness (handles the
  ! boundary_2d case where mat's z-column is identically zero, and any
  ! near-degenerate 3D element) is used instead of a closed-form inverse's
  ! own singularity edge cases.
  subroutine ensure_green_gauss_mat_cache(mesh, boundary_2d)
    use linear_solver_module, only: pseudo_inverse_inplace_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d

    integer(kind=ENTIER) :: v, i, j, a, b, k, n_cand
    integer(kind=ENTIER) :: id_sub_elem, id_elem, id_sub_face, total_pairs
    real(kind=DOUBLE), dimension(3) :: dx, norm, dminn, dmaxn
    real(kind=DOUBLE), dimension(3, 3) :: mat, mat_inv

    ! NB: .eqv. binds LOOSER than .and., so this must be parenthesized --
    ! `A .and. B .eqv. C` parses as `(A .and. B) .eqv. C`. Without the
    ! parens, a differently-sized mesh with boundary_2d=.false. (matching
    ! the cache's own default) can false-positive as "already built".
    if (gg_mat_cache_n_vert == mesh%n_vert .and. (gg_mat_cache_boundary_2d .eqv. boundary_2d)) return

    if (allocated(gg_mat_inv_cache))  deallocate(gg_mat_inv_cache)
    if (allocated(gg_mat_valid_cache)) deallocate(gg_mat_valid_cache)
    if (allocated(gg_flux_offset_cache)) deallocate(gg_flux_offset_cache)
    if (allocated(gg_flux_w_cache)) deallocate(gg_flux_w_cache)
    if (allocated(gg_flux_elem_cache)) deallocate(gg_flux_elem_cache)
    if (allocated(gg_oi_hlocal_cache)) deallocate(gg_oi_hlocal_cache)
    allocate(gg_mat_inv_cache(3, 3, mesh%n_vert))
    allocate(gg_mat_valid_cache(mesh%n_vert))
    allocate(gg_oi_hlocal_cache(mesh%n_vert))
    allocate(gg_flux_offset_cache(mesh%n_vert + 1))

    total_pairs = 0
    do v = 1, mesh%n_vert
      if (mesh%vert(v)%is_bound) cycle
      do i = 1, mesh%vert(v)%n_sub_elems_neigh
        total_pairs = total_pairs + mesh%sub_elem(mesh%vert(v)%sub_elem_neigh(i))%n_sub_faces
      end do
    end do
    allocate(gg_flux_w_cache(3, total_pairs))
    allocate(gg_flux_elem_cache(total_pairs))

    n_cand = merge(2_ENTIER, 3_ENTIER, boundary_2d)
    k = 0
    do v = 1, mesh%n_vert
      gg_flux_offset_cache(v) = k + 1
      if (mesh%vert(v)%is_bound) then
        gg_mat_valid_cache(v) = .false.
        gg_mat_inv_cache(:, :, v) = 0.0_DOUBLE
        gg_oi_hlocal_cache(v) = 0.0_DOUBLE
        cycle
      end if
      gg_mat_valid_cache(v) = .true.

      mat = 0.0_DOUBLE
      do i = 1, mesh%vert(v)%n_sub_elems_neigh
        id_sub_elem = mesh%vert(v)%sub_elem_neigh(i)
        id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
        dx = mesh%elem(id_elem)%coord - mesh%vert(v)%coord

        do j = 1, mesh%sub_elem(id_sub_elem)%n_sub_faces
          id_sub_face = mesh%sub_elem(id_sub_elem)%sub_face(j)
          ! sub_face%norm points from left_elem_neigh to right_elem_neigh;
          ! flip it to point OUTWARD from this sub_elem's own cell.
          if (mesh%sub_face(id_sub_face)%left_elem_neigh == id_elem) then
            norm = mesh%sub_face(id_sub_face)%norm
          else
            norm = -mesh%sub_face(id_sub_face)%norm
          end if
          do a = 1, 3
            do b = 1, 3
              mat(a, b) = mat(a, b) &
                + mesh%sub_face(id_sub_face)%area * norm(a) * dx(b)
            end do
          end do
          k = k + 1
          gg_flux_w_cache(:, k) = mesh%sub_face(id_sub_face)%area * norm
          gg_flux_elem_cache(k) = id_elem
        end do
      end do

      dminn(1:n_cand) = huge(1.0_DOUBLE)
      dmaxn(1:n_cand) = -huge(1.0_DOUBLE)
      do i = 1, size(mesh%vert(v)%elem_neigh)
        dx = mesh%elem(mesh%vert(v)%elem_neigh(i))%coord - mesh%vert(v)%coord
        do a = 1, n_cand
          dminn(a) = min(dminn(a), dx(a))
          dmaxn(a) = max(dmaxn(a), dx(a))
        end do
      end do
      gg_oi_hlocal_cache(v) = maxval(dmaxn(1:n_cand) - dminn(1:n_cand))

      ! mat's third column is identically zero for boundary_2d (every
      ! neighbor shares the vertex's own z on a thin single-layer mesh),
      ! singular by construction; the SVD handles that and any
      ! near-degenerate 3D element the same way.
      mat_inv = mat
      call pseudo_inverse_inplace_lapack(3_ENTIER, mat_inv)
      if (boundary_2d) mat_inv(3, :) = 0.0_DOUBLE
      gg_mat_inv_cache(:, :, v) = mat_inv
    end do
    gg_flux_offset_cache(mesh%n_vert + 1) = total_pairs + 1

    gg_mat_cache_n_vert = mesh%n_vert
    gg_mat_cache_boundary_2d = boundary_2d
  end subroutine ensure_green_gauss_mat_cache

  ! Alternative to compute_nodal_derivative_at_vertex ("aho ls"): a discrete
  ! Green-Gauss/divergence-theorem sum over the vertex's own dual control
  ! volume (median-dual sub_elem/sub_face decomposition), corrected by the
  ! inverse of the local geometric moment matrix so the result is exact for
  ! a linear field on any (including skewed) stencil -- "aho green gauss".
  !
  ! Derivation: for phi linear, phi(x)=phi(x_p)+A.(x-x_p), Green-Gauss over
  ! the closed interior dual cell gives grad_raw = phi(x_p)*sum_f(A_f*n_f) +
  ! sum_f A_f*n_f*(A.(x_elem_f-x_p)). The first term vanishes for an
  ! interior vertex (closed surface); boundary vertices are excluded
  ! (valid=.false.) rather than corrected for the term not vanishing there.
  ! The second term is mat.A with mat = sum_f A_f*(n_f (x) (x_elem_f-x_p)),
  ! so A = mat^-1.grad_raw. mat depends only on geometry, so the same
  ! mat^-1 (gg_mat_inv_cache) is reused unchanged at every recursion order
  ! (grad->hess->third): the extra carried indices of a higher-order phi
  ! are just extra right-hand-side columns.
  subroutine compute_nodal_derivative_at_vertex_green_gauss(mesh, d, nc_in, &
      boundary_2d, id_vert, phi, dphi_v, valid, oi_v)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in, id_vert
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(nc_in, mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(nc_in*d), intent(out) :: dphi_v
    logical, intent(out) :: valid
    real(kind=DOUBLE), intent(out) :: oi_v

    integer(kind=ENTIER) :: i1, a, k, k0, k1, id_elem, n_cand
    real(kind=DOUBLE), dimension(3, nc_in) :: grad_raw, grad_true
    real(kind=DOUBLE) :: phi_scale2, weight_sum
    real(kind=DOUBLE), dimension(nc_in) :: phi_sq_sum

    oi_v = 0.0_DOUBLE
    dphi_v = 0.0_DOUBLE

    call ensure_green_gauss_mat_cache(mesh, boundary_2d)

    valid = gg_mat_valid_cache(id_vert)
    if (.not. valid) return

    ! Pure geometry (weighted normal, source cell) is precomputed once per
    ! vertex in ensure_green_gauss_mat_cache; only this phi-dependent gather
    ! is redone at every call (grad/hess/third x every RK stage x every step).
    grad_raw = 0.0_DOUBLE
    k0 = gg_flux_offset_cache(id_vert)
    k1 = gg_flux_offset_cache(id_vert + 1) - 1
    do k = k0, k1
      id_elem = gg_flux_elem_cache(k)
      do a = 1, 3
        do i1 = 1, nc_in
          grad_raw(a, i1) = grad_raw(a, i1) + gg_flux_w_cache(a, k) * phi(i1, id_elem)
        end do
      end do
    end do

    grad_true = matmul(gg_mat_inv_cache(:, :, id_vert), grad_raw)

    ! Gradient-norm oscillation indicator, same formula/calibration idea as
    ! compute_nodal_derivative_at_vertex's own gradient-norm term (no
    ! residual counterpart here: GG doesn't explicitly fit an affine model,
    ! so there is no natural "residual against the fit"). h_local is pure
    ! geometry (cached above); phi_scale2 depends on phi and is gathered here.
    n_cand = merge(2_ENTIER, 3_ENTIER, boundary_2d)
    weight_sum = 0.0_DOUBLE
    phi_sq_sum = 0.0_DOUBLE
    do i1 = 1, size(mesh%vert(id_vert)%elem_neigh)
      weight_sum = weight_sum + 1.0_DOUBLE
      phi_sq_sum = phi_sq_sum + phi(:, mesh%vert(id_vert)%elem_neigh(i1))**2
    end do
    phi_scale2 = maxval(phi_sq_sum) / max(weight_sum, 1.0e-300_DOUBLE)
    oi_v = (gg_oi_hlocal_cache(id_vert) * sum(grad_true(1:n_cand, :)**2) &
      / max(phi_scale2, 1.0e-300_DOUBLE)) / grad_norm_derate_gg

    ! Same "direction slow-varying, carried-component fast-varying" flat
    ! layout as compute_nodal_derivative_at_vertex's own dphi_v.
    do a = 1, d
      do i1 = 1, nc_in
        dphi_v((a-1)*nc_in + i1) = grad_true(a, i1)
      end do
    end do
  end subroutine compute_nodal_derivative_at_vertex_green_gauss

  ! Exactness/accuracy test for compute_nodal_derivative_at_vertex_green_
  ! gauss, degree by degree: builds a polynomial phi of degree k whose exact
  ! order-k derivative is a known constant tensor, calls the routine at
  ! every non-boundary vertex, and prints max/RMS error for degree 1 (grad),
  ! 2 (hess, fed the exact linear gradient field), and 3 (third, fed the
  ! exact linear Hessian field), plus a degree-mismatch check (grad of a
  ! quadratic field, not exact by construction, to track convergence order).
  subroutine test_green_gauss_nodal(mesh, boundary_2d)
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d

    real(kind=DOUBLE), dimension(3) :: v1, v3, x0, dx
    real(kind=DOUBLE), dimension(3, 3) :: h2
    real(kind=DOUBLE), dimension(:, :), allocatable :: phi1, phi2, phi3
    real(kind=DOUBLE), dimension(:), allocatable :: dphi_v
    real(kind=DOUBLE) :: err, max_err, rms_err
    integer(kind=ENTIER) :: i, id_vert, n_valid, a, i1
    logical :: valid
    real(kind=DOUBLE) :: oi_v, s

    x0 = 0.0_DOUBLE
    v1 = [1.3_DOUBLE, -0.7_DOUBLE, 0.5_DOUBLE]
    v3 = [1.3_DOUBLE, -0.7_DOUBLE, 0.5_DOUBLE]
    h2 = reshape([2.0_DOUBLE, 0.3_DOUBLE, -0.1_DOUBLE, &
                  0.3_DOUBLE, -1.5_DOUBLE, 0.2_DOUBLE, &
                  -0.1_DOUBLE, 0.2_DOUBLE, 0.8_DOUBLE], [3, 3])
    if (boundary_2d) then
      v1(3) = 0.0_DOUBLE
      v3(3) = 0.0_DOUBLE
      h2(3, :) = 0.0_DOUBLE
      h2(:, 3) = 0.0_DOUBLE
    end if

    ! Degree 1: phi1 = v1.(x-x0), exact grad = v1 (constant).
    allocate(phi1(1, mesh%n_elems), dphi_v(3))
    do i = 1, mesh%n_elems
      phi1(1, i) = dot_product(v1, mesh%elem(i)%coord - x0)
    end do

    max_err = 0.0_DOUBLE
    rms_err = 0.0_DOUBLE
    n_valid = 0
    do id_vert = 1, mesh%n_vert
      if (mesh%vert(id_vert)%is_ghost) cycle
      call compute_nodal_derivative_at_vertex_green_gauss(mesh, 3_ENTIER, &
        1_ENTIER, boundary_2d, id_vert, phi1, dphi_v, valid, oi_v)
      if (.not. valid) cycle
      n_valid = n_valid + 1
      err = norm2(dphi_v - v1)
      max_err = max(max_err, err)
      rms_err = rms_err + err**2
    end do
    rms_err = sqrt(rms_err / max(n_valid, 1))
    print *, 'green_gauss degree=1 (grad of linear): n_valid=', n_valid, &
      ' max_err=', max_err, ' rms_err=', rms_err
    deallocate(phi1, dphi_v)

    ! Degree 2: phi2 = 1/2 (x-x0).H2.(x-x0), exact grad = H2.(x-x0), exact
    ! hess = H2 (constant).
    allocate(phi2(3, mesh%n_elems), dphi_v(9))
    do i = 1, mesh%n_elems
      dx = mesh%elem(i)%coord - x0
      phi2(:, i) = matmul(h2, dx)
    end do

    max_err = 0.0_DOUBLE
    rms_err = 0.0_DOUBLE
    n_valid = 0
    do id_vert = 1, mesh%n_vert
      if (mesh%vert(id_vert)%is_ghost) cycle
      call compute_nodal_derivative_at_vertex_green_gauss(mesh, 3_ENTIER, &
        3_ENTIER, boundary_2d, id_vert, phi2, dphi_v, valid, oi_v)
      if (.not. valid) cycle
      n_valid = n_valid + 1
      err = 0.0_DOUBLE
      do a = 1, 3
        do i1 = 1, 3
          err = err + (dphi_v((a-1)*3+i1) - h2(a, i1))**2
        end do
      end do
      err = sqrt(err)
      max_err = max(max_err, err)
      rms_err = rms_err + err**2
    end do
    rms_err = sqrt(rms_err / max(n_valid, 1))
    print *, 'green_gauss degree=2 (hess, from exact grad field): n_valid=', &
      n_valid, ' max_err=', max_err, ' rms_err=', rms_err
    deallocate(phi2, dphi_v)

    ! Degree 3: phi3 = 1/6 (v3.(x-x0))^3, exact hess = (v3.(x-x0))*v3(x)v3,
    ! exact third = v3(x)v3(x)v3 (constant, fully symmetric).
    allocate(phi3(9, mesh%n_elems), dphi_v(27))
    do i = 1, mesh%n_elems
      dx = mesh%elem(i)%coord - x0
      s = dot_product(v3, dx)
      do a = 1, 3
        do i1 = 1, 3
          phi3((a-1)*3+i1, i) = s * v3(a) * v3(i1)
        end do
      end do
    end do

    max_err = 0.0_DOUBLE
    rms_err = 0.0_DOUBLE
    n_valid = 0
    do id_vert = 1, mesh%n_vert
      if (mesh%vert(id_vert)%is_ghost) cycle
      call compute_nodal_derivative_at_vertex_green_gauss(mesh, 3_ENTIER, &
        9_ENTIER, boundary_2d, id_vert, phi3, dphi_v, valid, oi_v)
      if (.not. valid) cycle
      n_valid = n_valid + 1
      err = 0.0_DOUBLE
      do a = 1, 3
        do i1 = 1, 9
          err = err + (dphi_v((a-1)*9+i1) &
            - v3(a)*v3(1+(i1-1)/3)*v3(1+mod(i1-1, 3)))**2
        end do
      end do
      err = sqrt(err)
      max_err = max(max_err, err)
      rms_err = rms_err + err**2
    end do
    rms_err = sqrt(rms_err / max(n_valid, 1))
    print *, 'green_gauss degree=3 (third, from exact hess field): n_valid=', &
      n_valid, ' max_err=', max_err, ' rms_err=', rms_err
    deallocate(phi3, dphi_v)

    ! Degree-mismatch: grad (order=1 call) of the SAME quadratic phi2q used
    ! above -- not exact (scheme is only linear-exact); report the error so
    ! the caller can confirm the expected shrink rate under refinement.
    allocate(phi1(1, mesh%n_elems), dphi_v(3))
    do i = 1, mesh%n_elems
      dx = mesh%elem(i)%coord - x0
      phi1(1, i) = 0.5_DOUBLE * dot_product(dx, matmul(h2, dx))
    end do

    max_err = 0.0_DOUBLE
    rms_err = 0.0_DOUBLE
    n_valid = 0
    do id_vert = 1, mesh%n_vert
      if (mesh%vert(id_vert)%is_ghost) cycle
      call compute_nodal_derivative_at_vertex_green_gauss(mesh, 3_ENTIER, &
        1_ENTIER, boundary_2d, id_vert, phi1, dphi_v, valid, oi_v)
      if (.not. valid) cycle
      n_valid = n_valid + 1
      err = norm2(dphi_v - matmul(h2, mesh%vert(id_vert)%coord - x0))
      max_err = max(max_err, err)
      rms_err = rms_err + err**2
    end do
    rms_err = sqrt(rms_err / max(n_valid, 1))
    print *, 'green_gauss degree-mismatch (grad of quadratic, order=1 call): &
      &n_valid=', n_valid, ' max_err=', max_err, ' rms_err=', rms_err
    deallocate(phi1, dphi_v)
  end subroutine test_green_gauss_nodal

  ! Same idea as test_green_gauss_nodal's degree-mismatch check, but on the
  ! real (non-polynomial) stationary isentropic vortex density profile
  ! (domain [-10,10]^2, beta=5, gamma=1.4, t=0), compared against its exact
  ! analytic grad/hess/third. Static test only (no RK3, no flux, no WENO):
  ! isolates whether the Green-Gauss nodal operator alone gets the right
  ! order on smooth, non-polynomial data -- not the same as the paper's own
  ! vortex L2(rho) convergence table, which measures the full solver.
  subroutine test_green_gauss_vortex(mesh, boundary_2d)
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d

    real(kind=DOUBLE), dimension(:, :), allocatable :: phi1, phi2, phi3
    real(kind=DOUBLE), dimension(:), allocatable :: dphi_v
    real(kind=DOUBLE) :: err, max_err, rms_err
    real(kind=DOUBLE) :: rho_v, gx, gy, hxx, hxy, hyy, txxx, txxy, txyy, tyyy
    real(kind=DOUBLE) :: xv, yv
    integer(kind=ENTIER) :: i, id_vert, n_valid, a, i1
    logical :: valid
    real(kind=DOUBLE) :: oi_v

    ! Degree 1: grad(rho) at every cell centroid, exact vs analytic gx,gy.
    allocate(phi1(1, mesh%n_elems), dphi_v(3))
    do i = 1, mesh%n_elems
      call vortex_ref(mesh%elem(i)%coord(1), mesh%elem(i)%coord(2), &
        rho_v, gx, gy, hxx, hxy, hyy, txxx, txxy, txyy, tyyy)
      phi1(1, i) = rho_v
    end do

    max_err = 0.0_DOUBLE
    rms_err = 0.0_DOUBLE
    n_valid = 0
    do id_vert = 1, mesh%n_vert
      if (mesh%vert(id_vert)%is_ghost) cycle
      call compute_nodal_derivative_at_vertex_green_gauss(mesh, 3_ENTIER, &
        1_ENTIER, boundary_2d, id_vert, phi1, dphi_v, valid, oi_v)
      if (.not. valid) cycle
      xv = mesh%vert(id_vert)%coord(1); yv = mesh%vert(id_vert)%coord(2)
      call vortex_ref(xv, yv, rho_v, gx, gy, hxx, hxy, hyy, txxx, txxy, txyy, tyyy)
      n_valid = n_valid + 1
      err = norm2(dphi_v - [gx, gy, 0.0_DOUBLE])
      max_err = max(max_err, err)
      rms_err = rms_err + err**2
    end do
    rms_err = sqrt(rms_err / max(n_valid, 1))
    print *, 'green_gauss vortex order=1 (grad of rho): n_valid=', n_valid, &
      ' max_err=', max_err, ' rms_err=', rms_err
    deallocate(phi1, dphi_v)

    ! Degree 2: hess(rho), fed the exact analytic grad field.
    allocate(phi2(3, mesh%n_elems), dphi_v(9))
    do i = 1, mesh%n_elems
      call vortex_ref(mesh%elem(i)%coord(1), mesh%elem(i)%coord(2), &
        rho_v, gx, gy, hxx, hxy, hyy, txxx, txxy, txyy, tyyy)
      phi2(:, i) = [gx, gy, 0.0_DOUBLE]
    end do

    max_err = 0.0_DOUBLE
    rms_err = 0.0_DOUBLE
    n_valid = 0
    do id_vert = 1, mesh%n_vert
      if (mesh%vert(id_vert)%is_ghost) cycle
      call compute_nodal_derivative_at_vertex_green_gauss(mesh, 3_ENTIER, &
        3_ENTIER, boundary_2d, id_vert, phi2, dphi_v, valid, oi_v)
      if (.not. valid) cycle
      xv = mesh%vert(id_vert)%coord(1); yv = mesh%vert(id_vert)%coord(2)
      call vortex_ref(xv, yv, rho_v, gx, gy, hxx, hxy, hyy, txxx, txxy, txyy, tyyy)
      n_valid = n_valid + 1
      ! true Hess = [[hxx,hxy,0],[hxy,hyy,0],[0,0,0]], layout (a-1)*3+i1
      err = (dphi_v(1)-hxx)**2 + (dphi_v(2)-hxy)**2 + dphi_v(3)**2 &
          + (dphi_v(4)-hxy)**2 + (dphi_v(5)-hyy)**2 + dphi_v(6)**2 &
          + dphi_v(7)**2 + dphi_v(8)**2 + dphi_v(9)**2
      err = sqrt(err)
      max_err = max(max_err, err)
      rms_err = rms_err + err**2
    end do
    rms_err = sqrt(rms_err / max(n_valid, 1))
    print *, 'green_gauss vortex order=2 (hess, from exact grad field): &
      &n_valid=', n_valid, ' max_err=', max_err, ' rms_err=', rms_err
    deallocate(phi2, dphi_v)

    ! Degree 3: third(rho), fed the exact analytic hess field.
    allocate(phi3(9, mesh%n_elems), dphi_v(27))
    do i = 1, mesh%n_elems
      call vortex_ref(mesh%elem(i)%coord(1), mesh%elem(i)%coord(2), &
        rho_v, gx, gy, hxx, hxy, hyy, txxx, txxy, txyy, tyyy)
      phi3(:, i) = [hxx, hxy, 0.0_DOUBLE, hxy, hyy, 0.0_DOUBLE, &
                    0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE]
    end do

    max_err = 0.0_DOUBLE
    rms_err = 0.0_DOUBLE
    n_valid = 0
    do id_vert = 1, mesh%n_vert
      if (mesh%vert(id_vert)%is_ghost) cycle
      call compute_nodal_derivative_at_vertex_green_gauss(mesh, 3_ENTIER, &
        9_ENTIER, boundary_2d, id_vert, phi3, dphi_v, valid, oi_v)
      if (.not. valid) cycle
      xv = mesh%vert(id_vert)%coord(1); yv = mesh%vert(id_vert)%coord(2)
      call vortex_ref(xv, yv, rho_v, gx, gy, hxx, hxy, hyy, txxx, txxy, txyy, tyyy)
      n_valid = n_valid + 1
      ! true third tensor: fully symmetric in (x,y), zero if any index is z.
      ! dphi_v((a-1)*9+(a2-1)*3+i1): a=new (3rd) derivative direction,
      ! a2/i1=Hess indices carried from the previous step.
      err = 0.0_DOUBLE
      err = err + (dphi_v(1)-txxx)**2 + (dphi_v(2)-txxy)**2 + dphi_v(3)**2
      err = err + (dphi_v(4)-txxy)**2 + (dphi_v(5)-txyy)**2 + dphi_v(6)**2
      err = err + dphi_v(7)**2 + dphi_v(8)**2 + dphi_v(9)**2
      err = err + (dphi_v(10)-txxy)**2 + (dphi_v(11)-txyy)**2 + dphi_v(12)**2
      err = err + (dphi_v(13)-txyy)**2 + (dphi_v(14)-tyyy)**2 + dphi_v(15)**2
      err = err + dphi_v(16)**2 + dphi_v(17)**2 + dphi_v(18)**2
      err = err + dphi_v(19)**2 + dphi_v(20)**2 + dphi_v(21)**2
      err = err + dphi_v(22)**2 + dphi_v(23)**2 + dphi_v(24)**2
      err = err + dphi_v(25)**2 + dphi_v(26)**2 + dphi_v(27)**2
      err = sqrt(err)
      max_err = max(max_err, err)
      rms_err = rms_err + err**2
    end do
    rms_err = sqrt(rms_err / max(n_valid, 1))
    print *, 'green_gauss vortex order=3 (third, from exact hess field): &
      &n_valid=', n_valid, ' max_err=', max_err, ' rms_err=', rms_err
    deallocate(phi3, dphi_v)
  end subroutine test_green_gauss_vortex

  ! Exact analytic rho and its grad/hess/third (x,y components only -- the
  ! vortex is z-invariant) for the stationary (t=0) isentropic vortex used
  ! throughout this paper (vortex_prim in euler_ho_module: beta=5, gamma=1.4,
  ! domain-centred, no background advection at t=0). Generated via sympy
  ! (diff+cse) from the same closed-form rho, K=(gamma-1)*beta^2/(8*gamma*pi^2).
  pure subroutine vortex_ref(x, y, rho_v, gx, gy, hxx, hxy, hyy, &
      txxx, txxy, txyy, tyyy)
    implicit none

    real(kind=DOUBLE), intent(in) :: x, y
    real(kind=DOUBLE), intent(out) :: rho_v, gx, gy, hxx, hxy, hyy
    real(kind=DOUBLE), intent(out) :: txxx, txxy, txyy, tyyy

    real(kind=DOUBLE) :: cse0, cse1, cse2, cse3, cse4, cse5, cse6, cse7, &
      cse8, cse9, cse10, cse11, cse12, cse13, cse14, cse15, cse16, cse17, &
      cse18, cse19, cse20, cse21

    include "vortex_ref_cse.inc"
  end subroutine vortex_ref

  ! Builds neigh_cache_start/neigh_cache_list (CSR) by calling
  ! gather_ls_neighbors once per vertex -- the one-time cost that used to be
  ! paid on every one of compute_nodal_derivative_at_vertex's many calls per
  ! vertex. A no-op if the cache already matches this mesh's vertex count.
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

  ! Builds m2_cache_xx/xy/yy (each cell's own second geometric moment about
  ! its own centroid) once per distinct mesh, via degree-5 quadrature. A
  ! no-op if already built for this mesh.
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
  ! mesh%vert(id_vert)%elem_neigh, never ring-expanded (see
  ! compute_nodal_derivative_at_vertex's header).
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
  ! field, in place, using the per-vertex nodal Hessian (hess_v) and
  ! third-derivative tensor (third_v) -- both already computed as
  ! intermediates elsewhere, no extra neighbors or LS solves needed beyond
  ! each cell's own second geometric moment.
  !
  ! The vertex's own weighted-LS affine fit is biased by three additive
  ! contributions, all linear in the (assumed locally cubic) field and its
  ! own H,T: (1) the vertex's own Taylor terms beyond affine, evaluated
  ! exactly using its own Hv,Tv; (2) each neighbor cell's own
  ! cell-average-vs-point-value gap, 0.5*H(centroid_j):M2_j; (3) the gap
  ! from blending several corner (vertex) samples of a curved gradient
  ! field into one cell value, corrected via the discrete second moment of
  ! the touching vertices' own positions. All three vanish for affine or
  ! quadratic data; for cubic data they reproduce the full observed grad bias.
  !
  ! Only implemented for boundary_2d=.true.; a .false. call is a no-op.
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
    ! Per-vertex oscillation indicator from the gradient-level call of
    ! compute_next_order_derivative -- used to recompute the exact
    ! per-vertex weight (sub_elem_volume*vertex_weno_weight) that grad_cell
    ! was actually scattered with, matching its own vertex-to-cell blend
    ! instead of assuming a plain sub_elem_volume average.
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
    real(kind=DOUBLE), dimension(:, :), allocatable :: m2d_num_xx, m2d_num_xy, m2d_num_yy
    real(kind=DOUBLE), dimension(:), allocatable :: m2d_den
    real(kind=DOUBLE), dimension(:, :), allocatable :: t_num_xxx, t_num_xxy, t_num_xyy, t_num_yyy
    real(kind=DOUBLE) :: M2xx, M2xy, M2yy
    real(kind=DOUBLE), dimension(:), allocatable :: Tx1, Tx2, Tx3, Tx4, extra_x, extra_y

    if (.not. boundary_2d) return

    call ensure_m2_cache(mesh)
    call ensure_neighbor_cache(mesh)

    allocate(Hxx(nc_in), Hxy(nc_in), Hyy(nc_in), Txxx(nc_in), Txxy(nc_in), Txyy(nc_in), Tyyy(nc_in))
    allocate(HxxJ(nc_in), HxyJ(nc_in), HyyJ(nc_in), moment_j(nc_in))
    allocate(bias_num_x(nc_in, mesh%n_elems), bias_num_y(nc_in, mesh%n_elems), bias_den(mesh%n_elems))
    allocate(Tx1(nc_in), Tx2(nc_in), Tx3(nc_in), Tx4(nc_in), extra_x(nc_in), extra_y(nc_in))
    bias_num_x = 0.0_DOUBLE; bias_num_y = 0.0_DOUBLE; bias_den = 0.0_DOUBLE

    n_basis = 3 ! 1, dx, dy -- same affine basis as the Step-1 fit
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

      ! Match grad_cell's own vertex-to-cell weighting exactly.
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
    ! unbiased, blending several corner samples of a curved gradient field
    ! does not equal its value at the cell centroid -- correct using the
    ! touching vertices' own Tv and the discrete second moment of their
    ! positions about the cell centroid.
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
      ! Same grad_cell-consistent weight as the bias scatter above -- this
      ! redistributes grad's own per-vertex samples, so it must match
      ! grad's own blend, not third_v's.
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
      M2xx = m2d_num_xx(1,i)/m2d_den(i); M2xy = m2d_num_xy(1,i)/m2d_den(i); M2yy = m2d_num_yy(1,i)/m2d_den(i)
      Tx1 = t_num_xxx(:,i)/m2d_den(i); Tx2 = t_num_xxy(:,i)/m2d_den(i)
      Tx3 = t_num_xyy(:,i)/m2d_den(i); Tx4 = t_num_yyy(:,i)/m2d_den(i)
      extra_x = 0.5_DOUBLE*(Tx1*M2xx + 2.0_DOUBLE*Tx2*M2xy + Tx3*M2yy)
      extra_y = 0.5_DOUBLE*(Tx2*M2xx + 2.0_DOUBLE*Tx3*M2xy + Tx4*M2yy)
      grad_cell(1:nc_in, i) = grad_cell(1:nc_in, i) - extra_x
      grad_cell(nc_in+1:2*nc_in, i) = grad_cell(nc_in+1:2*nc_in, i) - extra_y
    end do

    deallocate(Hxx, Hxy, Hyy, Txxx, Txxy, Txyy, Tyyy, HxxJ, HxyJ, HyyJ, moment_j)
    deallocate(bias_num_x, bias_num_y, bias_den)
    deallocate(Tx1, Tx2, Tx3, Tx4, extra_x, extra_y)
    deallocate(m2d_num_xx, m2d_num_xy, m2d_num_yy, m2d_den)
    deallocate(t_num_xxx, t_num_xxy, t_num_xyy, t_num_yyy)
  end subroutine apply_grad_bias_correction

end module arbitrary_high_order_module
