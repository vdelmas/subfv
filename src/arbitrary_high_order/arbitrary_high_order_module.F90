! Arbitrary-order cell derivatives on unstructured meshes: recursive multi-D divided differences, one nodal LS/Green-Gauss fit + WENO scatter per recursion order (grad->hess->third), overlappable MPI exchange between orders.
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
  public :: apply_local_taylor_correction
  public :: apply_gradient_node_correction
  public :: compute_gradient_node_bias
  public :: compute_node_derivative_bias
  public :: recombine_derivative_regression
  public :: compute_2exact_hessian_aho_cls
  public :: apply_discrete_vertex_moment_correction

  ! eps_weno floors a WENO indicator ratio against literal division by zero.
  real(kind=DOUBLE), parameter :: eps_weno = tiny(1.0_DOUBLE)
  ! Regularizes oi_v in scatter_weno_weighted; eps_weight_num_deep is a looser floor for hess/third; GG has its own decoupled eps_weight_num_gg/_deep_gg.
  real(kind=DOUBLE), public :: eps_weight_num = 1.0e-2_DOUBLE
  real(kind=DOUBLE), public :: eps_weight_num_deep = 1.0_DOUBLE
  real(kind=DOUBLE), public :: eps_weight_num_gg = 1.0e-2_DOUBLE
  real(kind=DOUBLE), public :: eps_weight_num_deep_gg = 1.0_DOUBLE
  ! Exponent on the oscillation indicator: weight = omega_p/(eps+OI_v^weno_power).
  integer(kind=ENTIER), public :: weno_power = 1
  ! Calibration divisor on the gradient-norm term added to oi_v (LS) / oi_v's only term (GG).
  real(kind=DOUBLE), public :: grad_norm_derate = 1.0e4_DOUBLE
  real(kind=DOUBLE), public :: grad_norm_derate_gg = 1.0e4_DOUBLE
  ! GG-only: when .true., use omega_p/(eps+(1+OI)^p) instead of omega_p/(eps+OI^p).
  logical, public :: use_alt_gg_weight = .false.
  ! .false.: plain sub_elem_volume-weighted average of nodal derivatives, no oscillation-adaptive de-centering.
  logical, public :: use_weno_blend = .true.
  ! .false.=weighted least-squares nodal fit ("aho ls"); .true.=Green-Gauss divergence-theorem fit ("aho gg"); both produce the same (dphi_v,oi_v) shape.
  logical, public :: use_green_gauss = .false.
  ! CWENO-style central candidate (Semplice & Visconti 2020): blends the plain linear average into the WENO sum, weight cweno_center_weight/(eps+OI_c** cweno_center_power); off by default, opt in via compute_next_order_derivative_cweno.
  logical, public :: use_cweno_center = .false.
  real(kind=DOUBLE), public :: cweno_center_weight = 1000.0_DOUBLE
  integer(kind=ENTIER), public :: cweno_center_power = 4

  ! Per-vertex neighbor-list cache (CSR), built once per mesh: gather_ls_neighbors is pure topology.
  integer(kind=ENTIER), save :: neigh_cache_n_vert = -1
  integer(kind=ENTIER), dimension(:), allocatable, save :: neigh_cache_start
  integer(kind=ENTIER), dimension(:), allocatable, save :: neigh_cache_list

  ! Vertex-to-vertex adjacency (CSR), built once per mesh: v' is a neighbor of v if they share at
  ! least one cell. Pure topology (no such array exists on vert_type/mesh_type itself -- every
  ! _neigh field there points to elements/faces, never other vertices). Needed for aho_cls
  ! (Haider, Croisille & Courbet 2011): raising a vertex's k-exact k-th derivative to (k+1)-exact
  ! uses differences against NEIGHBORING VERTICES' own k-exact derivative, staying at the vertex
  ! level between orders instead of compute_next_order_derivative's node->cell->node chain.
  integer(kind=ENTIER), save :: vv_neigh_cache_n_vert = -1
  integer(kind=ENTIER), dimension(:), allocatable, save :: vv_neigh_cache_start
  integer(kind=ENTIER), dimension(:), allocatable, save :: vv_neigh_cache_list

  ! Per-cell second geometric moment, used by apply_grad_bias_correction: pure mesh geometry, cached once.
  integer(kind=ENTIER), save :: m2_cache_n_elems = -1
  real(kind=DOUBLE), dimension(:), allocatable, save :: m2_cache_xx, m2_cache_yy, m2_cache_zz
  real(kind=DOUBLE), dimension(:), allocatable, save :: m2_cache_xy, m2_cache_xz, m2_cache_yz

  ! Per-vertex Green-Gauss fit matrix, already inverted: pure geometry, cached once per mesh (boundary_2d baked in).
  integer(kind=ENTIER), save :: gg_mat_cache_n_vert = -1
  logical, save :: gg_mat_cache_boundary_2d = .false.
  real(kind=DOUBLE), dimension(:, :, :), allocatable, save :: gg_mat_inv_cache
  logical, dimension(:), allocatable, save :: gg_mat_valid_cache

  ! Flattened per-vertex (area*normal, source-cell) list for the GG flux sum, plus the OI indicator's h_local term: also pure geometry, built alongside gg_mat_inv_cache.
  integer(kind=ENTIER), dimension(:), allocatable, save :: gg_flux_offset_cache
  real(kind=DOUBLE), dimension(:, :), allocatable, save :: gg_flux_w_cache
  integer(kind=ENTIER), dimension(:), allocatable, save :: gg_flux_elem_cache
  real(kind=DOUBLE), dimension(:), allocatable, save :: gg_oi_hlocal_cache

  type :: derivative_field_type
    real(kind=DOUBLE), dimension(:, :), allocatable :: val ! (d**order, n_elems)
  end type derivative_field_type

  ! Per-cell moments about the cell's own centroid, M_c^{(m)} = (1/V_c) int_c (x-x_c)^{tensor m} dV,
  ! full flat tensor of size 3**m (row-major, same convention as hess_flat/third_flat). Cached once
  ! per mesh, m=2..max_order; feeds apply_local_taylor_correction.
  type :: cell_moment_ptr_type
    real(kind=DOUBLE), dimension(:, :), allocatable :: m ! (3**order, n_elems)
  end type cell_moment_ptr_type

  integer(kind=ENTIER), save :: cell_moment_cache_n_elems = -1
  integer(kind=ENTIER), save :: cell_moment_cache_max_order = 0
  type(cell_moment_ptr_type), dimension(:), allocatable, save :: cell_moment_cache ! indexed 2:max_order

  ! Discrete analogue of cell_moment_cache: mu_m^discrete(c) = average over the cell's own corner
  ! vertices of (x_v-x_c)^{tensor m}, for averaging a per-vertex SAMPLE (not a continuous field)
  ! over a cell's corners. See apply_discrete_vertex_moment_correction.
  integer(kind=ENTIER), save :: discrete_vmom_cache_n_elems = -1
  integer(kind=ENTIER), save :: discrete_vmom_cache_max_order = 0
  type(cell_moment_ptr_type), dimension(:), allocatable, save :: discrete_vmom_cache ! indexed 2:max_order

  ! Per-vertex 1-exact-gradient of the geometric field (x_elem-x_vert)^{tensor m} (Pont et al. 2017,
  ! JCP 350, eq. 56-61, generalized to arbitrary m): same GG flux-sum + gg_mat_inv_cache machinery
  ! used for phi, fed the geometric field instead. Flat tensor (3*3**m, n_vert), component t*3+i
  ! (t=geometric multi-index, i=gradient direction). Pure 1-ring geometry, cached once per mesh;
  ! corrects the gradient from any higher true derivative m=2..max_order (apply_gradient_node_correction).
  integer(kind=ENTIER), save :: gg_grad_h1_cache_n_vert = -1
  integer(kind=ENTIER), save :: gg_grad_h1_cache_max_order = 0
  logical, save :: gg_grad_h1_cache_boundary_2d = .false.
  type(cell_moment_ptr_type), dimension(:), allocatable, save :: gg_grad_h1_cache ! indexed 2:max_order

  ! One timestamped phase of one order's work, for the blocking-vs-overlap timeline figure.
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

  ! One step of the hierarchy: nc_in=d**(order-1) in, nc_out=nc_in*d out; boundary_2d must match the mesh's own build flag (drops z from the fit basis on a thin-extruded mesh).
  subroutine compute_next_order_derivative(mesh, d, nc_in, boundary_2d, phi, dphi, deriv_order, &
      dphi_v_out, valid_v_out, oi_v_out)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(nc_in, mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(nc_in*d, mesh%n_elems), intent(out) :: dphi
    ! Recursion level: 1=grad, 2=hess, 3=third -- selects eps_weight_num vs eps_weight_num_deep. Defaults to 1.
    integer(kind=ENTIER), intent(in), optional :: deriv_order
    ! Optional: exposes the per-vertex nodal estimate, used by apply_grad_bias_correction to avoid a redundant LS solve.
    real(kind=DOUBLE), dimension(:, :), allocatable, intent(out), optional :: dphi_v_out
    logical, dimension(:), allocatable, intent(out), optional :: valid_v_out
    ! Optional: exposes the per-vertex oscillation indicator, used by apply_grad_bias_correction to match grad_cell's own WENO weight.
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

    ! A boundary vertex is skipped (one-sided neighbor gather) unless .not. boundary_2d, where skipping it would starve every vertex.
    do id_vert = 1, mesh%n_vert
      if (mesh%vert(id_vert)%is_bound) cycle
      call accumulate_weno_contribution(mesh, d, nc_in, boundary_2d, id_vert, &
        phi, dphi_v_cache, valid_cache, oi_cache, weno_num, weno_den, &
        deriv_order=deriv_order_eff)
    end do

    ! Fallback for a cell left with zero weight (e.g. every touching vertex on the boundary): re-admit its own boundary vertices.
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

  ! CWENO variant: also accumulates a linear (unweighted) num/den and a volume-weighted OI average, folded back with weight cweno_center_weight/(eps+OI_c**power) -- see use_cweno_center.
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
        ! Rescue's fallback has no linear/OI counterpart -- mirror it into the linear accumulators too.
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

  ! Same per-vertex LS/GG solve as accumulate_weno_contribution, scattering into three running sums at once (WENO, linear/center, OI average).
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

    ! use_weno_blend=.false. must reduce weno_num/den to exactly lin_num/den (weight=1) so the center blend below is a provable no-op.
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

  ! Rescue for a cell left with zero WENO weight (every vertex on the boundary): re-admits just this cell's own boundary vertices.
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
      ! Under use_green_gauss this is a no-op for a fully-boundary cell: the GG fit always returns valid=.false. there.
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
    ! If every vertex is invalid too, weno_den is left at 0; the caller must not divide by it.
  end subroutine rescue_zero_weight_cell

  ! Same step as compute_next_order_derivative, overlapping its own ghost exchange with local work (send-cells first, post, remaining work, wait); degenerates to it when num_procs=1.
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

    ! Step 3: everything else, computed while the exchange is in flight (boundary vertices reuse their Step-1 cached dphi_v).
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

  ! Reconstructs every order 1..max_order from phi0, exchanging each order's ghost layer before building the next.
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

  ! Same as compute_derivative_hierarchy using the overlap step: phi_prev is already ghost-complete except phi0, exchanged once up front.
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

  ! Same as compute_derivative_hierarchy, timestamping exchange/compute into events (preallocate to >=2*max_order) for the timeline figure.
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

  ! Same as compute_derivative_hierarchy_overlap, timestamping all four phases into events (preallocate to >=4*max_order).
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

  ! Computes+caches (D^k phi)_v at id_vert and scatters it WENO-weighted into every touching cell; skip_cell excludes cells (used by the MPI-overlap two-pass split).
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

  ! Same scatter as accumulate_weno_contribution, reusing an already-cached dphi_v (MPI-overlap Step 3) instead of a second LAPACK solve.
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

  ! Shared scatter: weight=omega_p/(eps+OI^p), omega_p=sub_elem_volume, OI=oi_v (the vertex's own fit residual/gradient-norm indicator, see compute_nodal_derivative_at_vertex).
  subroutine scatter_weno_weighted(mesh, id_vert, dphi_v, oi_v, weno_num, weno_den, skip_cell, eps_in)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(:), intent(in) :: dphi_v
    real(kind=DOUBLE), intent(in) :: oi_v
    real(kind=DOUBLE), dimension(:, :), intent(inout) :: weno_num
    real(kind=DOUBLE), dimension(:), intent(inout) :: weno_den
    logical, dimension(:), intent(in), optional :: skip_cell
    ! Overrides eps_weight_num when present (level-dependent eps, grad vs hess/third).
    real(kind=DOUBLE), intent(in), optional :: eps_in

    integer(kind=ENTIER) :: j, id_elem, id_sub_elem
    real(kind=DOUBLE) :: sub_elem_volume, vertex_weno_weight, eps_use

    eps_use = eps_weight_num
    if (present(eps_in)) eps_use = eps_in

    if (use_weno_blend) then
      ! GG-only alternative bounds the smooth-data weight near 1/(eps+1) instead of blowing up to 1/eps as OI->0.
      if (use_green_gauss .and. use_alt_gg_weight) then
        vertex_weno_weight = 1.0_DOUBLE / (eps_use + (1.0_DOUBLE + oi_v)**weno_power)
      else
        vertex_weno_weight = 1.0_DOUBLE / (eps_use + oi_v**weno_power)
      end if
    else
      ! weight=1 cancels out of num/den, leaving a plain sub_elem_volume-weighted average.
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

  ! Weighted-LS gradient of phi at id_vert over elem_neigh (never ring-expanded: a wider ring isn't guaranteed complete at an MPI partition seam), dynamic basis dropping any direction without a resolvable spread; exact for a locally affine phi.
  subroutine compute_nodal_derivative_at_vertex(mesh, d, nc_in, boundary_2d, &
      id_vert, phi, dphi_v, valid, oi_v)
    use linear_solver_module, only: lu_factor_lapack, lu_solve_mat_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in, id_vert
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(nc_in, mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(nc_in*d), intent(out) :: dphi_v
    ! .false. when no direction could be resolved -- caller must exclude the vertex from the WENO scatter, not scatter dphi_v=0.
    logical, intent(out) :: valid
    ! WENO oscillation indicator: fit residual (dimensionless) combined via max with a gradient-norm term that catches an axis-aligned jump the residual is blind to.
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

    ! rhs(1+a,i1)=d(phi_i1)/dx_active_dim(a); flattened component-fast/direction-slow, inactive directions left at 0.
    do a = 1, n_active
      do i1 = 1, nc_in
        dphi_v((active_dim(a)-1)*nc_in + i1) = rhs(1+a, i1)
      end do
    end do
  end subroutine compute_nodal_derivative_at_vertex

  ! Builds gg_mat_inv_cache/gg_flux_*_cache/gg_oi_hlocal_cache (pure geometry) once per mesh; uses the LAPACK SVD pseudo-inverse unconditionally since inversion is now a one-time cost.
  subroutine ensure_green_gauss_mat_cache(mesh, boundary_2d)
    use linear_solver_module, only: pseudo_inverse_inplace_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d

    integer(kind=ENTIER) :: v, i, j, a, b, k, n_cand
    integer(kind=ENTIER) :: id_sub_elem, id_elem, id_sub_face, total_pairs
    real(kind=DOUBLE), dimension(3) :: dx, norm, dminn, dmaxn
    real(kind=DOUBLE), dimension(3, 3) :: mat, mat_inv

    ! NB: .eqv. binds looser than .and. -- must be parenthesized or a differently-sized mesh can false-positive as already cached.
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
          ! sub_face%norm points left_elem_neigh->right_elem_neigh; flip it outward from this sub_elem's own cell.
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

      ! mat's third column is identically zero for boundary_2d (singular by construction); the SVD handles that and any near-degenerate 3D element.
      mat_inv = mat
      call pseudo_inverse_inplace_lapack(3_ENTIER, mat_inv)
      if (boundary_2d) mat_inv(3, :) = 0.0_DOUBLE
      gg_mat_inv_cache(:, :, v) = mat_inv
    end do
    gg_flux_offset_cache(mesh%n_vert + 1) = total_pairs + 1

    gg_mat_cache_n_vert = mesh%n_vert
    gg_mat_cache_boundary_2d = boundary_2d
  end subroutine ensure_green_gauss_mat_cache

  ! Builds cell_moment_cache(m), m=2..max_order: each cell's own moment tensor M_c^{(m)} =
  ! (1/V_c) int_c (x-x_c)^{tensor m} dV, full flat tensor of size 3**m (row-major, same
  ! convention as hess_flat/third_flat). Pure geometry, cached once per mesh.
  subroutine ensure_cell_moment_cache(mesh, max_order)
    use quadrature_module, only: volume_quad_pts
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: max_order

    integer(kind=ENTIER) :: i, kv, n_v, m, t, idx, r, dloc, n_quad, quad_deg
    real(kind=DOUBLE), dimension(:, :), allocatable :: vcoords, pq
    real(kind=DOUBLE), dimension(:), allocatable :: wq, prod
    real(kind=DOUBLE), dimension(:, :), allocatable :: dxq
    real(kind=DOUBLE) :: xc, yc, zc

    if (cell_moment_cache_n_elems == mesh%n_elems .and. cell_moment_cache_max_order >= max_order) return

    if (allocated(cell_moment_cache)) deallocate(cell_moment_cache)
    allocate(cell_moment_cache(2:max_order))
    do m = 2, max_order
      allocate(cell_moment_cache(m)%m(3**m, mesh%n_elems))
      cell_moment_cache(m)%m = 0.0_DOUBLE
    end do

    quad_deg = max(5_ENTIER, max_order)

    do i = 1, mesh%n_elems
      n_v = mesh%elem(i)%n_vert
      allocate(vcoords(3, n_v))
      do kv = 1, n_v
        vcoords(:, kv) = mesh%vert(mesh%elem(i)%vert(kv))%coord
      end do
      call volume_quad_pts(n_v, vcoords, quad_deg, pq, wq)
      xc = mesh%elem(i)%coord(1); yc = mesh%elem(i)%coord(2); zc = mesh%elem(i)%coord(3)

      n_quad = size(wq)
      allocate(dxq(3, n_quad), prod(n_quad))
      dxq(1, :) = pq(1, :) - xc
      dxq(2, :) = pq(2, :) - yc
      dxq(3, :) = pq(3, :) - zc

      do m = 2, max_order
        do t = 0, 3**m - 1
          ! unravel t into its m base-3 digits (order irrelevant: product is commutative)
          idx = t
          prod = 1.0_DOUBLE
          do r = 1, m
            dloc = mod(idx, 3) + 1
            idx = idx / 3
            prod = prod * dxq(dloc, :)
          end do
          cell_moment_cache(m)%m(t+1, i) = sum(wq * prod) / mesh%elem(i)%volume
        end do
      end do
      deallocate(dxq, prod, vcoords, pq, wq)
    end do

    cell_moment_cache_n_elems = mesh%n_elems
    cell_moment_cache_max_order = max_order
  end subroutine ensure_cell_moment_cache

  ! z_vK^(m) = (1/|T_K|) int_{T_K} (x-x_v)^{tensor m} dV, for ONE neighbor cell id_elem and shift
  ! s = x_elem - x_v, full flat tensor of size 3**m. Binomial/moment-shift expansion (Haider,
  ! Croisille & Courbet 2011, eq. 1/13's z_{alpha,beta}): writing x-x_v = s + (x-x_K),
  !   z^(m)[i_1..i_m] = sum over subsets S of {1..m} ( prod_{l not in S} s_{i_l} ) * M_K^{(|S|)}[i_l, l in S]
  ! with M_K^(0)=1, M_K^(1)=0, M_K^(j)=cell_moment_cache(j) for j>=2 -- folds the node-stencil
  ! geometry bias (S={} term) and K's own cell-average-vs-point-value gap (S!={} terms) into one
  ! moment, reusing cell_moment_cache algebraically instead of a second GG pass per correction.
  subroutine shifted_cell_moment_full(id_elem, s, m, z)
    implicit none

    integer(kind=ENTIER), intent(in) :: id_elem, m
    real(kind=DOUBLE), dimension(3), intent(in) :: s
    real(kind=DOUBLE), dimension(3**m), intent(out) :: z

    integer(kind=ENTIER) :: t, idx, l, mask, j, sub_idx
    integer(kind=ENTIER), dimension(m) :: digits
    real(kind=DOUBLE) :: prod_s, mval

    z = 0.0_DOUBLE
    do t = 0, 3**m - 1
      idx = t
      do l = m, 1, -1
        digits(l) = mod(idx, 3)
        idx = idx / 3
      end do
      do mask = 0, 2**m - 1
        j = popcnt(mask)
        prod_s = 1.0_DOUBLE
        sub_idx = 0
        do l = 1, m
          if (btest(mask, l-1)) then
            sub_idx = sub_idx*3 + digits(l)
          else
            prod_s = prod_s * s(digits(l)+1)
          end if
        end do
        if (j == 0) then
          mval = 1.0_DOUBLE
        else if (j == 1) then
          mval = 0.0_DOUBLE
        else
          mval = cell_moment_cache(j)%m(sub_idx+1, id_elem)
        end if
        z(t+1) = z(t+1) + prod_s * mval
      end do
    end do
  end subroutine shifted_cell_moment_full

  ! Builds gg_grad_h1_cache(m), m=2..max_order: for each vertex v, the SAME 1-exact GG gradient
  ! operator (flux-sum + gg_mat_inv_cache) applied to the geometric field {z_vK^(m)}_K
  ! (shifted_cell_moment_full) instead of phi. Pure geometry (v's own 1-ring plus each touching
  ! cell's own cached moments), cached once per mesh.
  subroutine ensure_gg_gradient_h1_cache(mesh, boundary_2d, max_order)
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d
    integer(kind=ENTIER), intent(in) :: max_order

    integer(kind=ENTIER) :: v, k, id_elem, m, tt
    real(kind=DOUBLE), dimension(3) :: dx_v, tmp
    real(kind=DOUBLE), dimension(:, :), allocatable :: raw
    real(kind=DOUBLE), dimension(:), allocatable :: z

    call ensure_green_gauss_mat_cache(mesh, boundary_2d)
    call ensure_cell_moment_cache(mesh, max_order)

    if (gg_grad_h1_cache_n_vert == mesh%n_vert .and. gg_grad_h1_cache_max_order >= max_order &
        .and. (gg_grad_h1_cache_boundary_2d .eqv. boundary_2d)) return

    if (allocated(gg_grad_h1_cache)) deallocate(gg_grad_h1_cache)
    allocate(gg_grad_h1_cache(2:max_order))
    do m = 2, max_order
      allocate(gg_grad_h1_cache(m)%m(3*3**m, mesh%n_vert))
      gg_grad_h1_cache(m)%m = 0.0_DOUBLE
    end do

    do v = 1, mesh%n_vert
      if (.not. gg_mat_valid_cache(v)) cycle
      do m = 2, max_order
        allocate(raw(3, 3**m), z(3**m))
        raw = 0.0_DOUBLE
        do k = gg_flux_offset_cache(v), gg_flux_offset_cache(v+1) - 1
          id_elem = gg_flux_elem_cache(k)
          dx_v = mesh%elem(id_elem)%coord - mesh%vert(v)%coord
          call shifted_cell_moment_full(id_elem, dx_v, m, z)
          raw = raw + spread(gg_flux_w_cache(:, k), 2, 3**m) * spread(z, 1, 3)
        end do
        do tt = 1, 3**m
          tmp = matmul(gg_mat_inv_cache(:, :, v), raw(:, tt))
          gg_grad_h1_cache(m)%m((tt-1)*3+1 : (tt-1)*3+3, v) = tmp
        end do
        deallocate(raw, z)
      end do
    end do

    gg_grad_h1_cache_n_vert = mesh%n_vert
    gg_grad_h1_cache_max_order = max_order
    gg_grad_h1_cache_boundary_2d = boundary_2d
  end subroutine ensure_gg_gradient_h1_cache

  ! Alternative to the LS fit ("aho gg"): Green-Gauss sum over the vertex's dual control volume, corrected by mat^-1=gg_mat_inv_cache (see that cache's own header for the derivation) so the result is exact for a linear field on any stencil.
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

    ! Geometry is precomputed in ensure_green_gauss_mat_cache; only this phi gather is redone every call.
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

    ! Gradient-norm oscillation indicator (no residual counterpart for GG); h_local is cached, phi_scale2 depends on phi and is gathered here.
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

    ! Same direction-slow/component-fast flat layout as compute_nodal_derivative_at_vertex's dphi_v.
    do a = 1, d
      do i1 = 1, nc_in
        dphi_v((a-1)*nc_in + i1) = grad_true(a, i1)
      end do
    end do
  end subroutine compute_nodal_derivative_at_vertex_green_gauss

  ! Builds the CSR neighbor cache once per mesh so compute_nodal_derivative_at_vertex doesn't re-derive it on every call.
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

  ! Builds vv_neigh_cache: for each vertex v, the set of OTHER vertices sharing at least one cell
  ! with v (deduplicated), via a two-pass CSR build using a reusable "last touched by v" mark
  ! array (avoids an O(n_vert) reset per vertex). Pure topology, cached once per mesh.
  subroutine ensure_vv_neighbor_cache(mesh)
    implicit none

    type(mesh_type), intent(in) :: mesh

    integer(kind=ENTIER) :: v, i, kv, id_elem, n_v, w, n_uniq, p
    integer(kind=ENTIER), dimension(:), allocatable :: mark

    if (vv_neigh_cache_n_vert == mesh%n_vert) return

    if (allocated(vv_neigh_cache_start)) deallocate(vv_neigh_cache_start)
    if (allocated(vv_neigh_cache_list))  deallocate(vv_neigh_cache_list)
    allocate(vv_neigh_cache_start(mesh%n_vert + 1))
    allocate(mark(mesh%n_vert))
    mark = 0

    vv_neigh_cache_start(1) = 1
    do v = 1, mesh%n_vert
      n_uniq = 0
      do i = 1, mesh%vert(v)%n_elems_neigh
        id_elem = mesh%vert(v)%elem_neigh(i)
        n_v = mesh%elem(id_elem)%n_vert
        do kv = 1, n_v
          w = mesh%elem(id_elem)%vert(kv)
          if (w == v) cycle
          if (mark(w) /= v) then
            mark(w) = v
            n_uniq = n_uniq + 1
          end if
        end do
      end do
      vv_neigh_cache_start(v+1) = vv_neigh_cache_start(v) + n_uniq
    end do

    allocate(vv_neigh_cache_list(vv_neigh_cache_start(mesh%n_vert+1) - 1))
    mark = 0
    do v = 1, mesh%n_vert
      p = vv_neigh_cache_start(v)
      do i = 1, mesh%vert(v)%n_elems_neigh
        id_elem = mesh%vert(v)%elem_neigh(i)
        n_v = mesh%elem(id_elem)%n_vert
        do kv = 1, n_v
          w = mesh%elem(id_elem)%vert(kv)
          if (w == v) cycle
          if (mark(w) /= v) then
            mark(w) = v
            vv_neigh_cache_list(p) = w
            p = p + 1
          end if
        end do
      end do
    end do
    deallocate(mark)

    vv_neigh_cache_n_vert = mesh%n_vert
  end subroutine ensure_vv_neighbor_cache

  ! Builds each cell's own second geometric moment once per mesh via degree-5 quadrature.
  subroutine ensure_m2_cache(mesh)
    use quadrature_module, only: volume_quad_pts
    implicit none

    type(mesh_type), intent(in) :: mesh

    integer(kind=ENTIER) :: i, kv, n_v
    real(kind=DOUBLE), dimension(:, :), allocatable :: vcoords, pq
    real(kind=DOUBLE), dimension(:), allocatable :: wq
    real(kind=DOUBLE) :: xc, yc, zc

    if (m2_cache_n_elems == mesh%n_elems) return

    if (allocated(m2_cache_xx)) deallocate(m2_cache_xx)
    if (allocated(m2_cache_yy)) deallocate(m2_cache_yy)
    if (allocated(m2_cache_zz)) deallocate(m2_cache_zz)
    if (allocated(m2_cache_xy)) deallocate(m2_cache_xy)
    if (allocated(m2_cache_xz)) deallocate(m2_cache_xz)
    if (allocated(m2_cache_yz)) deallocate(m2_cache_yz)
    allocate(m2_cache_xx(mesh%n_elems), m2_cache_yy(mesh%n_elems), m2_cache_zz(mesh%n_elems))
    allocate(m2_cache_xy(mesh%n_elems), m2_cache_xz(mesh%n_elems), m2_cache_yz(mesh%n_elems))

    do i = 1, mesh%n_elems
      n_v = mesh%elem(i)%n_vert
      allocate(vcoords(3, n_v))
      do kv = 1, n_v
        vcoords(:, kv) = mesh%vert(mesh%elem(i)%vert(kv))%coord
      end do
      call volume_quad_pts(n_v, vcoords, 5_ENTIER, pq, wq)
      xc = mesh%elem(i)%coord(1); yc = mesh%elem(i)%coord(2); zc = mesh%elem(i)%coord(3)
      m2_cache_xx(i) = sum(wq*(pq(1,:)-xc)**2) / sum(wq)
      m2_cache_yy(i) = sum(wq*(pq(2,:)-yc)**2) / sum(wq)
      m2_cache_zz(i) = sum(wq*(pq(3,:)-zc)**2) / sum(wq)
      m2_cache_xy(i) = sum(wq*(pq(1,:)-xc)*(pq(2,:)-yc)) / sum(wq)
      m2_cache_xz(i) = sum(wq*(pq(1,:)-xc)*(pq(3,:)-zc)) / sum(wq)
      m2_cache_yz(i) = sum(wq*(pq(2,:)-yc)*(pq(3,:)-zc)) / sum(wq)
      deallocate(vcoords, pq, wq)
    end do

    m2_cache_n_elems = mesh%n_elems
  end subroutine ensure_m2_cache

  ! n!/(n-kk)! for small non-negative integers (kk<=n): converts Taylor coefficients to/from raw
  ! partial derivatives without a separate factorial table.
  pure function fact_ratio(n, kk) result(r)
    implicit none

    integer(kind=ENTIER), intent(in) :: n, kk
    real(kind=DOUBLE) :: r
    integer(kind=ENTIER) :: i

    r = 1.0_DOUBLE
    do i = n-kk+1, n
      r = r * real(i, kind=DOUBLE)
    end do
  end function fact_ratio

  ! Element neighbors of id_vert: exactly elem_neigh, never ring-expanded (see compute_nodal_derivative_at_vertex).
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

  ! Corrects the grad step's own O(h^2) bias against a cubic field in place, using each
  ! vertex's nodal Hessian/third tensor and per-cell second moment. Handles both
  ! boundary_2d=.true. (2 active directions x,y) and genuinely-3D meshes (3 active
  ! directions x,y,z) through one unified formula: the basis size (3 or 4) is the only
  ! branch, everything else is the full 3D tensor contraction, which degenerates exactly
  ! to the 2D case when z-related H/T components are zero (as they are for boundary_2d).
  ! hess_v/third_v flat layout: index (0-based, in units of nc_in) for a component whose
  ! derivative directions are taken in order (dir1,dir2[,dir3]) is
  ! (dir3-1)*9+(dir2-1)*3+(dir1-1) (third_v) or (dir2-1)*3+(dir1-1) (hess_v), 1=x,2=y,3=z;
  ! symmetric components are averaged over every valid ordering (verified numerically
  ! against manufactured cubic fields with fully distinct coefficients).
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
    ! Per-vertex OI from the gradient-level call, used to match grad_cell's own vertex-to-cell WENO weight exactly.
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: grad_oi_v

    integer(kind=ENTIER) :: iv, j, id_elem, id_sub_elem, ic, i
    integer(kind=ENTIER) :: n_neigh, n_basis, n_cand
    integer(kind=ENTIER), dimension(:), allocatable :: neigh
    real(kind=DOUBLE), dimension(4, 4) :: mat
    real(kind=DOUBLE), dimension(4, nc_in) :: rhs
    real(kind=DOUBLE), dimension(4) :: basis
    integer(kind=ENTIER), dimension(4) :: ipiv
    real(kind=DOUBLE), dimension(3) :: dx_v
    real(kind=DOUBLE) :: weight, dxx, dyy, dzz, sub_elem_volume, dvx, dvy, dvz, vweight
    real(kind=DOUBLE), dimension(:), allocatable :: Hxx, Hyy, Hzz, Hxy, Hxz, Hyz
    real(kind=DOUBLE), dimension(:), allocatable :: Txxx, Tyyy, Tzzz, Txxy, Txxz, Txyy, Tyyz, Txzz, Tyzz, Txyz
    real(kind=DOUBLE), dimension(:), allocatable :: HxxJ, HyyJ, HzzJ, HxyJ, HxzJ, HyzJ, moment_j
    real(kind=DOUBLE), dimension(:, :), allocatable :: bias_num_x, bias_num_y, bias_num_z
    real(kind=DOUBLE), dimension(:), allocatable :: bias_den
    real(kind=DOUBLE), dimension(:), allocatable :: m2d_num_xx, m2d_num_yy, m2d_num_zz
    real(kind=DOUBLE), dimension(:), allocatable :: m2d_num_xy, m2d_num_xz, m2d_num_yz, m2d_den
    real(kind=DOUBLE), dimension(:, :), allocatable :: t_num_xxx, t_num_yyy, t_num_zzz, t_num_xxy, t_num_xxz
    real(kind=DOUBLE), dimension(:, :), allocatable :: t_num_xyy, t_num_yyz, t_num_xzz, t_num_yzz, t_num_xyz
    real(kind=DOUBLE) :: M2xx, M2yy, M2zz, M2xy, M2xz, M2yz
    real(kind=DOUBLE), dimension(:), allocatable :: Tx1, Tx2, Tx3, Tx4, Tx5, Tx6, Tx7, Tx8, Tx9, Tx10
    real(kind=DOUBLE), dimension(:), allocatable :: extra_x, extra_y, extra_z

    call ensure_m2_cache(mesh)
    call ensure_neighbor_cache(mesh)

    n_cand = merge(2_ENTIER, 3_ENTIER, boundary_2d)
    n_basis = 1 + n_cand

    allocate(Hxx(nc_in), Hyy(nc_in), Hzz(nc_in), Hxy(nc_in), Hxz(nc_in), Hyz(nc_in))
    allocate(Txxx(nc_in), Tyyy(nc_in), Tzzz(nc_in), Txxy(nc_in), Txxz(nc_in))
    allocate(Txyy(nc_in), Tyyz(nc_in), Txzz(nc_in), Tyzz(nc_in), Txyz(nc_in))
    allocate(HxxJ(nc_in), HyyJ(nc_in), HzzJ(nc_in), HxyJ(nc_in), HxzJ(nc_in), HyzJ(nc_in), moment_j(nc_in))
    allocate(bias_num_x(nc_in, mesh%n_elems), bias_num_y(nc_in, mesh%n_elems), bias_num_z(nc_in, mesh%n_elems))
    allocate(bias_den(mesh%n_elems))
    bias_num_x = 0.0_DOUBLE; bias_num_y = 0.0_DOUBLE; bias_num_z = 0.0_DOUBLE; bias_den = 0.0_DOUBLE

    do iv = 1, mesh%n_vert
      if (mesh%vert(iv)%is_bound) cycle
      if (.not. valid_hess_v(iv) .or. .not. valid_third_v(iv)) cycle

      do ic = 1, nc_in
        Hxx(ic) = hess_v(0*nc_in+ic, iv)
        Hyy(ic) = hess_v(4*nc_in+ic, iv)
        Hzz(ic) = hess_v(8*nc_in+ic, iv)
        Hxy(ic) = 0.5_DOUBLE*(hess_v(1*nc_in+ic, iv) + hess_v(3*nc_in+ic, iv))
        Hxz(ic) = 0.5_DOUBLE*(hess_v(2*nc_in+ic, iv) + hess_v(6*nc_in+ic, iv))
        Hyz(ic) = 0.5_DOUBLE*(hess_v(5*nc_in+ic, iv) + hess_v(7*nc_in+ic, iv))
        Txxx(ic) = third_v(0*nc_in+ic, iv)
        Tyyy(ic) = third_v(13*nc_in+ic, iv)
        Tzzz(ic) = third_v(26*nc_in+ic, iv)
        Txxy(ic) = (third_v(1*nc_in+ic, iv) + third_v(3*nc_in+ic, iv) + third_v(9*nc_in+ic, iv)) / 3.0_DOUBLE
        Txxz(ic) = (third_v(2*nc_in+ic, iv) + third_v(6*nc_in+ic, iv) + third_v(18*nc_in+ic, iv)) / 3.0_DOUBLE
        Txyy(ic) = (third_v(4*nc_in+ic, iv) + third_v(10*nc_in+ic, iv) + third_v(12*nc_in+ic, iv)) / 3.0_DOUBLE
        Tyyz(ic) = (third_v(14*nc_in+ic, iv) + third_v(16*nc_in+ic, iv) + third_v(22*nc_in+ic, iv)) / 3.0_DOUBLE
        Txzz(ic) = (third_v(8*nc_in+ic, iv) + third_v(20*nc_in+ic, iv) + third_v(24*nc_in+ic, iv)) / 3.0_DOUBLE
        Tyzz(ic) = (third_v(17*nc_in+ic, iv) + third_v(23*nc_in+ic, iv) + third_v(25*nc_in+ic, iv)) / 3.0_DOUBLE
        Txyz(ic) = (third_v(5*nc_in+ic, iv) + third_v(7*nc_in+ic, iv) + third_v(11*nc_in+ic, iv) &
          + third_v(15*nc_in+ic, iv) + third_v(19*nc_in+ic, iv) + third_v(21*nc_in+ic, iv)) / 6.0_DOUBLE
      end do

      n_neigh = neigh_cache_start(iv+1) - neigh_cache_start(iv)
      allocate(neigh(n_neigh))
      neigh = neigh_cache_list(neigh_cache_start(iv):neigh_cache_start(iv+1)-1)
      mat(1:n_basis, 1:n_basis) = 0.0_DOUBLE
      rhs(1:n_basis, :) = 0.0_DOUBLE
      do j = 1, n_neigh
        id_elem = neigh(j)
        dx_v = mesh%elem(id_elem)%coord - mesh%vert(iv)%coord
        weight = 1.0_DOUBLE / max(dot_product(dx_v, dx_v), 1.0e-24_DOUBLE)
        dxx = dx_v(1); dyy = dx_v(2); dzz = dx_v(3)
        basis(1) = 1.0_DOUBLE; basis(2) = dxx; basis(3) = dyy
        if (n_basis == 4) basis(4) = dzz
        mat(1:n_basis, 1:n_basis) = mat(1:n_basis, 1:n_basis) &
          + weight * spread(basis(1:n_basis), 2, n_basis) * spread(basis(1:n_basis), 1, n_basis)

        ! (1) vertex's own Taylor terms beyond affine (full 3D; z-terms vanish
        ! identically when boundary_2d, since dzz=0 for every neighbor there).
        moment_j = 0.5_DOUBLE*(Hxx*dxx**2 + Hyy*dyy**2 + Hzz*dzz**2 &
          + 2.0_DOUBLE*Hxy*dxx*dyy + 2.0_DOUBLE*Hxz*dxx*dzz + 2.0_DOUBLE*Hyz*dyy*dzz) &
          + (Txxx*dxx**3 + Tyyy*dyy**3 + Tzzz*dzz**3 &
             + 3.0_DOUBLE*Txxy*dxx**2*dyy + 3.0_DOUBLE*Txxz*dxx**2*dzz &
             + 3.0_DOUBLE*Txyy*dxx*dyy**2 + 3.0_DOUBLE*Tyyz*dyy**2*dzz &
             + 3.0_DOUBLE*Txzz*dxx*dzz**2 + 3.0_DOUBLE*Tyzz*dyy*dzz**2 &
             + 6.0_DOUBLE*Txyz*dxx*dyy*dzz) / 6.0_DOUBLE

        ! (2) this neighbor's own cell-average-vs-point-value gap
        HxxJ = Hxx + Txxx*dxx + Txxy*dyy + Txxz*dzz
        HyyJ = Hyy + Txyy*dxx + Tyyy*dyy + Tyyz*dzz
        HzzJ = Hzz + Txzz*dxx + Tyzz*dyy + Tzzz*dzz
        HxyJ = Hxy + Txxy*dxx + Txyy*dyy + Txyz*dzz
        HxzJ = Hxz + Txxz*dxx + Txyz*dyy + Txzz*dzz
        HyzJ = Hyz + Txyz*dxx + Tyyz*dyy + Tyzz*dzz
        moment_j = moment_j + 0.5_DOUBLE*(HxxJ*m2_cache_xx(id_elem) + HyyJ*m2_cache_yy(id_elem) &
          + HzzJ*m2_cache_zz(id_elem) + 2.0_DOUBLE*HxyJ*m2_cache_xy(id_elem) &
          + 2.0_DOUBLE*HxzJ*m2_cache_xz(id_elem) + 2.0_DOUBLE*HyzJ*m2_cache_yz(id_elem))

        do ic = 1, nc_in
          rhs(1:n_basis, ic) = rhs(1:n_basis, ic) + weight*basis(1:n_basis)*moment_j(ic)
        end do
      end do
      deallocate(neigh)

      call lu_factor_lapack(n_basis, mat(1:n_basis, 1:n_basis), ipiv(1:n_basis))
      call lu_solve_mat_lapack(n_basis, mat(1:n_basis, 1:n_basis), ipiv(1:n_basis), nc_in, rhs(1:n_basis, :))

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
        if (n_basis == 4) bias_num_z(:, id_elem) = bias_num_z(:, id_elem) + sub_elem_volume*rhs(4, :)
        bias_den(id_elem) = bias_den(id_elem) + sub_elem_volume
      end do
    end do

    do i = 1, mesh%n_elems
      if (bias_den(i) <= 0.0_DOUBLE) cycle
      grad_cell(1:nc_in, i) = grad_cell(1:nc_in, i) - bias_num_x(:, i)/bias_den(i)
      grad_cell(nc_in+1:2*nc_in, i) = grad_cell(nc_in+1:2*nc_in, i) - bias_num_y(:, i)/bias_den(i)
      if (.not. boundary_2d) grad_cell(2*nc_in+1:3*nc_in, i) = grad_cell(2*nc_in+1:3*nc_in, i) - bias_num_z(:, i)/bias_den(i)
    end do

    ! (3) cell-blend curvature: correct the gap from blending several corner samples of a curved gradient field via the touching vertices' own Tv and the discrete second moment of their positions about the cell centroid.
    allocate(m2d_num_xx(mesh%n_elems), m2d_num_yy(mesh%n_elems), m2d_num_zz(mesh%n_elems))
    allocate(m2d_num_xy(mesh%n_elems), m2d_num_xz(mesh%n_elems), m2d_num_yz(mesh%n_elems), m2d_den(mesh%n_elems))
    allocate(t_num_xxx(nc_in, mesh%n_elems), t_num_yyy(nc_in, mesh%n_elems), t_num_zzz(nc_in, mesh%n_elems))
    allocate(t_num_xxy(nc_in, mesh%n_elems), t_num_xxz(nc_in, mesh%n_elems), t_num_xyy(nc_in, mesh%n_elems))
    allocate(t_num_yyz(nc_in, mesh%n_elems), t_num_xzz(nc_in, mesh%n_elems), t_num_yzz(nc_in, mesh%n_elems))
    allocate(t_num_xyz(nc_in, mesh%n_elems))
    allocate(Tx1(nc_in), Tx2(nc_in), Tx3(nc_in), Tx4(nc_in), Tx5(nc_in))
    allocate(Tx6(nc_in), Tx7(nc_in), Tx8(nc_in), Tx9(nc_in), Tx10(nc_in))
    allocate(extra_x(nc_in), extra_y(nc_in), extra_z(nc_in))
    m2d_num_xx = 0.0_DOUBLE; m2d_num_yy = 0.0_DOUBLE; m2d_num_zz = 0.0_DOUBLE
    m2d_num_xy = 0.0_DOUBLE; m2d_num_xz = 0.0_DOUBLE; m2d_num_yz = 0.0_DOUBLE; m2d_den = 0.0_DOUBLE
    t_num_xxx = 0.0_DOUBLE; t_num_yyy = 0.0_DOUBLE; t_num_zzz = 0.0_DOUBLE
    t_num_xxy = 0.0_DOUBLE; t_num_xxz = 0.0_DOUBLE; t_num_xyy = 0.0_DOUBLE
    t_num_yyz = 0.0_DOUBLE; t_num_xzz = 0.0_DOUBLE; t_num_yzz = 0.0_DOUBLE; t_num_xyz = 0.0_DOUBLE
    do iv = 1, mesh%n_vert
      if (mesh%vert(iv)%is_bound) cycle
      if (.not. valid_third_v(iv)) cycle
      do ic = 1, nc_in
        Txxx(ic) = third_v(0*nc_in+ic, iv)
        Tyyy(ic) = third_v(13*nc_in+ic, iv)
        Tzzz(ic) = third_v(26*nc_in+ic, iv)
        Txxy(ic) = (third_v(1*nc_in+ic, iv) + third_v(3*nc_in+ic, iv) + third_v(9*nc_in+ic, iv)) / 3.0_DOUBLE
        Txxz(ic) = (third_v(2*nc_in+ic, iv) + third_v(6*nc_in+ic, iv) + third_v(18*nc_in+ic, iv)) / 3.0_DOUBLE
        Txyy(ic) = (third_v(4*nc_in+ic, iv) + third_v(10*nc_in+ic, iv) + third_v(12*nc_in+ic, iv)) / 3.0_DOUBLE
        Tyyz(ic) = (third_v(14*nc_in+ic, iv) + third_v(16*nc_in+ic, iv) + third_v(22*nc_in+ic, iv)) / 3.0_DOUBLE
        Txzz(ic) = (third_v(8*nc_in+ic, iv) + third_v(20*nc_in+ic, iv) + third_v(24*nc_in+ic, iv)) / 3.0_DOUBLE
        Tyzz(ic) = (third_v(17*nc_in+ic, iv) + third_v(23*nc_in+ic, iv) + third_v(25*nc_in+ic, iv)) / 3.0_DOUBLE
        Txyz(ic) = (third_v(5*nc_in+ic, iv) + third_v(7*nc_in+ic, iv) + third_v(11*nc_in+ic, iv) &
          + third_v(15*nc_in+ic, iv) + third_v(19*nc_in+ic, iv) + third_v(21*nc_in+ic, iv)) / 6.0_DOUBLE
      end do
      ! Same grad_cell-consistent weight as the bias scatter above (redistributes grad's own samples).
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
        dvz = mesh%vert(iv)%coord(3) - mesh%elem(id_elem)%coord(3)
        m2d_num_xx(id_elem) = m2d_num_xx(id_elem) + sub_elem_volume*dvx*dvx
        m2d_num_yy(id_elem) = m2d_num_yy(id_elem) + sub_elem_volume*dvy*dvy
        m2d_num_zz(id_elem) = m2d_num_zz(id_elem) + sub_elem_volume*dvz*dvz
        m2d_num_xy(id_elem) = m2d_num_xy(id_elem) + sub_elem_volume*dvx*dvy
        m2d_num_xz(id_elem) = m2d_num_xz(id_elem) + sub_elem_volume*dvx*dvz
        m2d_num_yz(id_elem) = m2d_num_yz(id_elem) + sub_elem_volume*dvy*dvz
        m2d_den(id_elem) = m2d_den(id_elem) + sub_elem_volume
        t_num_xxx(:,id_elem) = t_num_xxx(:,id_elem) + sub_elem_volume*Txxx
        t_num_yyy(:,id_elem) = t_num_yyy(:,id_elem) + sub_elem_volume*Tyyy
        t_num_zzz(:,id_elem) = t_num_zzz(:,id_elem) + sub_elem_volume*Tzzz
        t_num_xxy(:,id_elem) = t_num_xxy(:,id_elem) + sub_elem_volume*Txxy
        t_num_xxz(:,id_elem) = t_num_xxz(:,id_elem) + sub_elem_volume*Txxz
        t_num_xyy(:,id_elem) = t_num_xyy(:,id_elem) + sub_elem_volume*Txyy
        t_num_yyz(:,id_elem) = t_num_yyz(:,id_elem) + sub_elem_volume*Tyyz
        t_num_xzz(:,id_elem) = t_num_xzz(:,id_elem) + sub_elem_volume*Txzz
        t_num_yzz(:,id_elem) = t_num_yzz(:,id_elem) + sub_elem_volume*Tyzz
        t_num_xyz(:,id_elem) = t_num_xyz(:,id_elem) + sub_elem_volume*Txyz
      end do
    end do
    do i = 1, mesh%n_elems
      if (m2d_den(i) <= 0.0_DOUBLE) cycle
      M2xx = m2d_num_xx(i)/m2d_den(i); M2yy = m2d_num_yy(i)/m2d_den(i); M2zz = m2d_num_zz(i)/m2d_den(i)
      M2xy = m2d_num_xy(i)/m2d_den(i); M2xz = m2d_num_xz(i)/m2d_den(i); M2yz = m2d_num_yz(i)/m2d_den(i)
      Tx1 = t_num_xxx(:,i)/m2d_den(i); Tx2 = t_num_yyy(:,i)/m2d_den(i); Tx3 = t_num_zzz(:,i)/m2d_den(i)
      Tx4 = t_num_xxy(:,i)/m2d_den(i); Tx5 = t_num_xxz(:,i)/m2d_den(i); Tx6 = t_num_xyy(:,i)/m2d_den(i)
      Tx7 = t_num_yyz(:,i)/m2d_den(i); Tx8 = t_num_xzz(:,i)/m2d_den(i); Tx9 = t_num_yzz(:,i)/m2d_den(i)
      Tx10 = t_num_xyz(:,i)/m2d_den(i)
      extra_x = 0.5_DOUBLE*(Tx1*M2xx + Tx6*M2yy + Tx8*M2zz + 2.0_DOUBLE*Tx4*M2xy + 2.0_DOUBLE*Tx5*M2xz + 2.0_DOUBLE*Tx10*M2yz)
      extra_y = 0.5_DOUBLE*(Tx4*M2xx + Tx2*M2yy + Tx9*M2zz + 2.0_DOUBLE*Tx6*M2xy + 2.0_DOUBLE*Tx10*M2xz + 2.0_DOUBLE*Tx7*M2yz)
      extra_z = 0.5_DOUBLE*(Tx5*M2xx + Tx7*M2yy + Tx3*M2zz + 2.0_DOUBLE*Tx10*M2xy + 2.0_DOUBLE*Tx8*M2xz + 2.0_DOUBLE*Tx9*M2yz)
      grad_cell(1:nc_in, i) = grad_cell(1:nc_in, i) - extra_x
      grad_cell(nc_in+1:2*nc_in, i) = grad_cell(nc_in+1:2*nc_in, i) - extra_y
      if (.not. boundary_2d) grad_cell(2*nc_in+1:3*nc_in, i) = grad_cell(2*nc_in+1:3*nc_in, i) - extra_z
    end do

    deallocate(Hxx, Hyy, Hzz, Hxy, Hxz, Hyz)
    deallocate(Txxx, Tyyy, Tzzz, Txxy, Txxz, Txyy, Tyyz, Txzz, Tyzz, Txyz)
    deallocate(HxxJ, HyyJ, HzzJ, HxyJ, HxzJ, HyzJ, moment_j)
    deallocate(bias_num_x, bias_num_y, bias_num_z, bias_den)
    deallocate(Tx1, Tx2, Tx3, Tx4, Tx5, Tx6, Tx7, Tx8, Tx9, Tx10, extra_x, extra_y, extra_z)
    deallocate(m2d_num_xx, m2d_num_yy, m2d_num_zz, m2d_num_xy, m2d_num_xz, m2d_num_yz, m2d_den)
    deallocate(t_num_xxx, t_num_yyy, t_num_zzz, t_num_xxy, t_num_xxz)
    deallocate(t_num_xyy, t_num_yyz, t_num_xzz, t_num_yzz, t_num_xyz)
  end subroutine apply_grad_bias_correction

  ! Contracts the LAST m=(k-q) indices of the rank-k tensor D_k (flat, row-major, size d**k*nc_in
  ! per cell) against the rank-m geometric tensor M_m (flat, row-major, size d**m per cell),
  ! leaving a rank-q tensor (flat, size d**q*nc_in per cell). Since both operands are stored in the
  ! FULL (redundant) flat convention already used for hess_flat/third_flat, this plain index-matched
  ! sum reproduces the correct symmetric-contraction multinomial weights automatically -- no
  ! separate combinatorial bookkeeping needed, unlike an independent-component (Hxx/Hxy/...)
  ! representation.
  subroutine contract_last_indices(d, nc_in, n_elems, k, m, D_k, M_m, corr_q)
    implicit none

    integer(kind=ENTIER), intent(in) :: d, nc_in, n_elems, k, m
    real(kind=DOUBLE), dimension(d**k*nc_in, n_elems), intent(in) :: D_k
    real(kind=DOUBLE), dimension(d**m, n_elems), intent(in) :: M_m
    real(kind=DOUBLE), dimension(d**(k-m)*nc_in, n_elems), intent(out) :: corr_q

    integer(kind=ENTIER) :: tq, tm, ic, n_tq, n_tm, off

    n_tq = d**(k-m)
    n_tm = d**m

    corr_q = 0.0_DOUBLE
    do tm = 0, n_tm - 1
      do tq = 0, n_tq - 1
        off = (tq*n_tm + tm) * nc_in
        do ic = 1, nc_in
          corr_q(tq*nc_in+ic, :) = corr_q(tq*nc_in+ic, :) + D_k(off+ic, :) * M_m(tm+1, :)
        end do
      end do
    end do
  end subroutine contract_last_indices

  ! Corrects, in place, the cell-level Taylor coefficients D^(0)=phi,...,D^(k_max) held in
  ! dfield(0:k_max)%val: each is a RAW cell-AVERAGE estimate, not the point value at x_c.
  ! Taylor-expanding about x_c and averaging over the cell (M_c^(1)=0 identically):
  !   D^(q)(x_c) = avg(D^(q))_c - sum_{m=2}^{k_max-q} (1/m!) D^(q+m)(x_c) : M_c^(m)
  ! applied incrementally k=2,...,k_max using the RAW dfield(k)%val (Pont et al. 2017, JCP 350,
  ! sec. 3.4's successive-correction idea, but entirely local to the cell: only cell_moment_cache).
  subroutine apply_local_taylor_correction(mesh, d, nc_in, k_max, dfield)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in, k_max
    type(derivative_field_type), dimension(0:k_max), intent(inout) :: dfield

    integer(kind=ENTIER) :: kk, q, m
    real(kind=DOUBLE), dimension(:, :), allocatable :: corr

    if (k_max < 2) return

    call ensure_cell_moment_cache(mesh, k_max)

    do kk = 2, k_max
      do q = 0, kk - 2
        m = kk - q
        allocate(corr(nc_in*d**q, mesh%n_elems))
        call contract_last_indices(d, nc_in, mesh%n_elems, kk, m, &
          dfield(kk)%val, cell_moment_cache(m)%m, corr)
        dfield(q)%val = dfield(q)%val - corr / fact_ratio(m, m)
        deallocate(corr)
      end do
    end do
  end subroutine apply_local_taylor_correction

  ! Contracts the m geometric indices of H_m^(1)(v) (flat, layout t*3+i, t=0..3**m-1, i=1..3)
  ! against the rank-m derivative tensor D_m_v (flat, layout t*nc_in+ic, same row-major
  ! convention as hess_v/third_v), leaving the rank-1 (gradient-direction i, component ic) bias
  ! bias_v(i,ic,vertex). Same "full flat tensor, plain index-matched sum" trick as
  ! contract_last_indices: the symmetric-contraction multinomial weights fall out automatically.
  subroutine contract_grad_node_bias(nc_in, m, n_vert, D_m_v, H_m1_v, bias_v)
    implicit none

    integer(kind=ENTIER), intent(in) :: nc_in, m, n_vert
    real(kind=DOUBLE), dimension(3**m*nc_in, n_vert), intent(in) :: D_m_v
    real(kind=DOUBLE), dimension(3**m*3, n_vert), intent(in) :: H_m1_v
    real(kind=DOUBLE), dimension(3*nc_in, n_vert), intent(out) :: bias_v

    integer(kind=ENTIER) :: t, i, ic, n_t

    n_t = 3**m
    bias_v = 0.0_DOUBLE
    do t = 0, n_t - 1
      do i = 1, 3
        do ic = 1, nc_in
          bias_v((i-1)*nc_in+ic, :) = bias_v((i-1)*nc_in+ic, :) &
            + D_m_v(t*nc_in+ic, :) * H_m1_v(t*3+i, :)
        end do
      end do
    end do
  end subroutine contract_grad_node_bias

  ! Corrects a RAW order-k vertex derivative (built by one application of the GG operator L_v to
  ! the order-(k-1) cell field) using every available higher vertex derivative
  ! dfield_v(k+1)%val,...,dfield_v(k_max)%val:
  !   (D^k phi)_v^corrected = (D^k phi)_v^raw - sum_{m=2}^{k_max-k+1} (1/m!) D^(k-1+m)(x_v) : H_m^(1)(v)
  ! (Pont et al. 2017, JCP 350, sec. 3.4 generalized to arbitrary order, folded into one moment
  ! per Haider, Croisille & Courbet 2011 eq. 13). H_m^(1)(v)=gg_grad_h1_cache(m) is the SAME
  ! operator L_v applied to every level of the aho_gg recursion, so this ONE cache corrects any
  ! order k by treating D^(k-1)'s d**(k-1) components as independent scalar fields (nc_in
  ! replaced by nc_in*d**(k-1) in eq:grad-bias-gg). k=1 is the gradient case.
  subroutine compute_node_derivative_bias(mesh, d, nc_in, boundary_2d, k, k_max, dfield_v, bias_v_total)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in, k, k_max
    logical, intent(in) :: boundary_2d
    type(derivative_field_type), dimension(k+1:k_max), intent(in) :: dfield_v
    real(kind=DOUBLE), dimension(nc_in*d**k, mesh%n_vert), intent(out) :: bias_v_total

    integer(kind=ENTIER) :: m, nc_eff
    real(kind=DOUBLE), dimension(:, :), allocatable :: bias_m

    bias_v_total = 0.0_DOUBLE
    if (k_max < k+1) return

    call ensure_gg_gradient_h1_cache(mesh, boundary_2d, k_max-k+1)

    nc_eff = nc_in * d**(k-1)
    do m = 2, k_max-k+1
      allocate(bias_m(nc_eff*d, mesh%n_vert))
      call contract_grad_node_bias(nc_eff, m, mesh%n_vert, dfield_v(k-1+m)%val, gg_grad_h1_cache(m)%m, bias_m)
      bias_v_total = bias_v_total + bias_m / fact_ratio(m, m)
      deallocate(bias_m)
    end do
  end subroutine compute_node_derivative_bias

  ! k=1 special case of compute_node_derivative_bias, kept for backward compatibility with
  ! existing callers (apply_gradient_node_correction).
  subroutine compute_gradient_node_bias(mesh, d, nc_in, boundary_2d, k_max, dfield_v, bias_v_total)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in, k_max
    logical, intent(in) :: boundary_2d
    type(derivative_field_type), dimension(2:k_max), intent(in) :: dfield_v
    real(kind=DOUBLE), dimension(nc_in*d, mesh%n_vert), intent(out) :: bias_v_total

    call compute_node_derivative_bias(mesh, d, nc_in, boundary_2d, 1_ENTIER, k_max, dfield_v, bias_v_total)
  end subroutine compute_gradient_node_bias

  ! Recombines a per-vertex quantity into cells via a weighted least-squares affine regression
  ! (val_v(x) ~= val_cell + G.(x-x_c), inverse-square-distance weight) over each cell's own
  ! touching vertices, instead of a naive average (only exact when the sampled quantity is
  ! constant across the cell). Falls back to a plain average when a cell has too few valid
  ! vertices for the fit (e.g. near a boundary).
  subroutine recombine_derivative_regression(mesh, boundary_2d, n_comp, val_v, valid_v, val_cell_out)
    use linear_solver_module, only: lu_factor_lapack, lu_solve_mat_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d
    integer(kind=ENTIER), intent(in) :: n_comp
    real(kind=DOUBLE), dimension(n_comp, mesh%n_vert), intent(in) :: val_v
    logical, dimension(mesh%n_vert), intent(in) :: valid_v
    real(kind=DOUBLE), dimension(n_comp, mesh%n_elems), intent(out) :: val_cell_out

    integer(kind=ENTIER) :: i, kv, n_v, id_v, n_cand, n_basis, n_valid, cnt
    real(kind=DOUBLE), dimension(4, 4) :: mat
    real(kind=DOUBLE), dimension(4, n_comp) :: rhs
    real(kind=DOUBLE), dimension(4) :: basis
    integer(kind=ENTIER), dimension(4) :: ipiv
    real(kind=DOUBLE), dimension(3) :: dx_v
    real(kind=DOUBLE) :: weight
    real(kind=DOUBLE), dimension(n_comp) :: s

    n_cand = merge(2_ENTIER, 3_ENTIER, boundary_2d)
    n_basis = 1 + n_cand
    val_cell_out = 0.0_DOUBLE

    do i = 1, mesh%n_elems
      n_v = mesh%elem(i)%n_vert
      n_valid = 0
      do kv = 1, n_v
        if (valid_v(mesh%elem(i)%vert(kv))) n_valid = n_valid + 1
      end do

      if (n_valid < n_basis) then
        s = 0.0_DOUBLE; cnt = 0
        do kv = 1, n_v
          id_v = mesh%elem(i)%vert(kv)
          if (.not. valid_v(id_v)) cycle
          s = s + val_v(:, id_v)
          cnt = cnt + 1
        end do
        if (cnt > 0) val_cell_out(:, i) = s / real(cnt, kind=DOUBLE)
        cycle
      end if

      mat(1:n_basis, 1:n_basis) = 0.0_DOUBLE
      rhs(1:n_basis, :) = 0.0_DOUBLE
      do kv = 1, n_v
        id_v = mesh%elem(i)%vert(kv)
        if (.not. valid_v(id_v)) cycle
        dx_v = mesh%vert(id_v)%coord - mesh%elem(i)%coord
        weight = 1.0_DOUBLE / max(dot_product(dx_v, dx_v), 1.0e-24_DOUBLE)
        basis(1) = 1.0_DOUBLE; basis(2) = dx_v(1); basis(3) = dx_v(2)
        if (n_basis == 4) basis(4) = dx_v(3)
        mat(1:n_basis, 1:n_basis) = mat(1:n_basis, 1:n_basis) &
          + weight * spread(basis(1:n_basis), 2, n_basis) * spread(basis(1:n_basis), 1, n_basis)
        rhs(1:n_basis, :) = rhs(1:n_basis, :) &
          + weight * spread(basis(1:n_basis), 2, n_comp) * spread(val_v(:, id_v), 1, n_basis)
      end do

      call lu_factor_lapack(n_basis, mat(1:n_basis, 1:n_basis), ipiv(1:n_basis))
      call lu_solve_mat_lapack(n_basis, mat(1:n_basis, 1:n_basis), ipiv(1:n_basis), n_comp, rhs(1:n_basis, :))

      val_cell_out(:, i) = rhs(1, :)
    end do
  end subroutine recombine_derivative_regression

  subroutine apply_gradient_node_correction(mesh, d, nc_in, boundary_2d, grad_cell, k_max, &
      dfield_v, valid_v)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in, k_max
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(nc_in*d, mesh%n_elems), intent(inout) :: grad_cell
    type(derivative_field_type), dimension(2:k_max), intent(in) :: dfield_v
    logical, dimension(mesh%n_vert), intent(in) :: valid_v

    integer(kind=ENTIER) :: iv, j, id_elem, id_sub_elem
    real(kind=DOUBLE), dimension(:, :), allocatable :: bias_v_total
    real(kind=DOUBLE), dimension(:, :), allocatable :: bias_num
    real(kind=DOUBLE), dimension(:), allocatable :: bias_den
    real(kind=DOUBLE) :: sub_elem_volume

    if (k_max < 2) return

    allocate(bias_v_total(nc_in*d, mesh%n_vert))
    call compute_gradient_node_bias(mesh, d, nc_in, boundary_2d, k_max, dfield_v, bias_v_total)

    allocate(bias_num(nc_in*d, mesh%n_elems), bias_den(mesh%n_elems))
    bias_num = 0.0_DOUBLE; bias_den = 0.0_DOUBLE

    do iv = 1, mesh%n_vert
      if (mesh%vert(iv)%is_bound) cycle
      if (.not. valid_v(iv)) cycle
      do j = 1, mesh%vert(iv)%n_sub_elems_neigh
        id_sub_elem = mesh%vert(iv)%sub_elem_neigh(j)
        id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
        sub_elem_volume = mesh%sub_elem(id_sub_elem)%volume
        bias_num(:, id_elem) = bias_num(:, id_elem) + sub_elem_volume * bias_v_total(:, iv)
        bias_den(id_elem) = bias_den(id_elem) + sub_elem_volume
      end do
    end do

    do j = 1, mesh%n_elems
      if (bias_den(j) <= 0.0_DOUBLE) cycle
      grad_cell(:, j) = grad_cell(:, j) - bias_num(:, j) / bias_den(j)
    end do

    deallocate(bias_v_total, bias_num, bias_den)
  end subroutine apply_gradient_node_correction

  ! aho_cls (Haider, Croisille & Courbet 2011, "Efficient Implementation of High Order
  ! Reconstruction in Finite Volume Methods"): contracts H_2^(1) (gg_grad_h1_cache(2), full flat,
  ! layout t*3+i) with a rank-2 tensor E9 (full flat, layout (i-1)*3+j), leaving the rank-1
  ! (direction i) result -- used to build the functional-identity matrix below.
  pure function contract_geom_h2(Hcol, E9) result(vec3)
    implicit none

    real(kind=DOUBLE), dimension(27), intent(in) :: Hcol
    real(kind=DOUBLE), dimension(9), intent(in) :: E9
    real(kind=DOUBLE), dimension(3) :: vec3

    integer(kind=ENTIER) :: t, i

    vec3 = 0.0_DOUBLE
    do t = 0, 8
      do i = 1, 3
        vec3(i) = vec3(i) + Hcol(t*3+i) * E9(t+1)
      end do
    end do
  end function contract_geom_h2

  ! aho_cls, step k=1->2: builds a GENUINELY 2-exact Hessian at every vertex directly from
  ! ONE-RING vertex data, following Haider, Croisille & Courbet (2011) eq. 15-18 ("functional
  ! identity"), adapted from their cell-based setting to our vertex-based dual stencil. For a
  ! genuinely quadratic field u with true Hessian H, eq:grad-bias-gg with k_max=2 gives at any
  ! vertex w: w_w^(1|1)[u] = grad(u)(x_w) + (1/2) H:H_2^(1)(w) exactly. Subtracting this relation
  ! at a neighbor v' from the one at v, and using grad(u)(x_v')-grad(u)(x_v)=H.(x_v'-x_v) exactly:
  !   w_v'^(1|1)[u] - w_v^(1|1)[u] = H.(x_v'-x_v) + (1/2) H : [H_2^(1)(v') - H_2^(1)(v)]
  ! a known linear map J_v(b) in a candidate tensor b -- solving J_v(b)={...}_v' by least squares
  ! recovers b=H exactly when u is quadratic. Needs only vv_neigh_cache and gg_grad_h1_cache(2).
  subroutine compute_2exact_hessian_aho_cls(mesh, boundary_2d, grad_v, valid_grad_v, hess_v_out, valid_out)
    use linear_solver_module, only: lu_factor_lapack, lu_solve_mat_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: grad_v
    logical, dimension(mesh%n_vert), intent(in) :: valid_grad_v
    real(kind=DOUBLE), dimension(6, mesh%n_vert), intent(out) :: hess_v_out ! Hxx,Hxy,Hxz,Hyy,Hyz,Hzz
    logical, dimension(mesh%n_vert), intent(out) :: valid_out

    integer(kind=ENTIER) :: v, i, vp, ell, row, n_valid
    real(kind=DOUBLE), dimension(:, :), allocatable :: Jmat
    real(kind=DOUBLE), dimension(:), allocatable :: rhsvec
    real(kind=DOUBLE), dimension(3) :: h, gcorr_v, gcorr_vp
    real(kind=DOUBLE), dimension(9) :: E9
    real(kind=DOUBLE), dimension(6, 6) :: normal_mat
    real(kind=DOUBLE), dimension(6, 1) :: normal_rhs
    integer(kind=ENTIER), dimension(6) :: ipiv
    integer(kind=ENTIER), parameter :: basis_i(6) = (/1,1,1,2,2,3/)
    integer(kind=ENTIER), parameter :: basis_j(6) = (/1,2,3,2,3,3/)

    call ensure_vv_neighbor_cache(mesh)
    call ensure_gg_gradient_h1_cache(mesh, boundary_2d, 2_ENTIER)

    hess_v_out = 0.0_DOUBLE
    valid_out = .false.

    do v = 1, mesh%n_vert
      if (mesh%vert(v)%is_bound) cycle
      if (.not. valid_grad_v(v)) cycle

      n_valid = 0
      do i = vv_neigh_cache_start(v), vv_neigh_cache_start(v+1) - 1
        if (valid_grad_v(vv_neigh_cache_list(i))) n_valid = n_valid + 1
      end do
      if (n_valid < 2) cycle

      allocate(Jmat(3*n_valid, 6), rhsvec(3*n_valid))
      row = 0
      do i = vv_neigh_cache_start(v), vv_neigh_cache_start(v+1) - 1
        vp = vv_neigh_cache_list(i)
        if (.not. valid_grad_v(vp)) cycle
        h = mesh%vert(vp)%coord - mesh%vert(v)%coord
        do ell = 1, 6
          E9 = 0.0_DOUBLE
          E9((basis_i(ell)-1)*3 + basis_j(ell)) = 1.0_DOUBLE
          E9((basis_j(ell)-1)*3 + basis_i(ell)) = 1.0_DOUBLE
          gcorr_v  = contract_geom_h2(gg_grad_h1_cache(2)%m(:, v),  E9)
          gcorr_vp = contract_geom_h2(gg_grad_h1_cache(2)%m(:, vp), E9)
          Jmat(row*3+1, ell) = E9(1)*h(1) + E9(2)*h(2) + E9(3)*h(3) + 0.5_DOUBLE*(gcorr_vp(1)-gcorr_v(1))
          Jmat(row*3+2, ell) = E9(4)*h(1) + E9(5)*h(2) + E9(6)*h(3) + 0.5_DOUBLE*(gcorr_vp(2)-gcorr_v(2))
          Jmat(row*3+3, ell) = E9(7)*h(1) + E9(8)*h(2) + E9(9)*h(3) + 0.5_DOUBLE*(gcorr_vp(3)-gcorr_v(3))
        end do
        rhsvec(row*3+1:row*3+3) = grad_v(:, vp) - grad_v(:, v)
        row = row + 1
      end do

      normal_mat = matmul(transpose(Jmat), Jmat)
      normal_rhs(:, 1) = matmul(transpose(Jmat), rhsvec)
      call lu_factor_lapack(6_ENTIER, normal_mat, ipiv)
      call lu_solve_mat_lapack(6_ENTIER, normal_mat, ipiv, 1_ENTIER, normal_rhs)

      hess_v_out(:, v) = normal_rhs(:, 1)
      valid_out(v) = .true.
      deallocate(Jmat, rhsvec)
    end do
  end subroutine compute_2exact_hessian_aho_cls

  ! Builds discrete_vmom_cache(m), m=2..max_order: for each cell, mu_m^discrete(c) = the average
  ! over its OWN corner vertices v of (x_v-x_c)^{tensor m}, as a full flat tensor of size 3**m
  ! (same row-major convention as cell_moment_cache). Pure topology+geometry, cached once per
  ! mesh -- unlike cell_moment_cache's quadrature, this is a plain average over mesh%elem(:)%vert.
  subroutine ensure_discrete_vertex_moment_cache(mesh, max_order)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: max_order

    integer(kind=ENTIER) :: i, kv, n_v, m, t, idx, r, dloc, id_v
    real(kind=DOUBLE), dimension(3) :: dx_v
    real(kind=DOUBLE) :: prodval

    if (discrete_vmom_cache_n_elems == mesh%n_elems .and. discrete_vmom_cache_max_order >= max_order) return

    if (allocated(discrete_vmom_cache)) deallocate(discrete_vmom_cache)
    allocate(discrete_vmom_cache(2:max_order))
    do m = 2, max_order
      allocate(discrete_vmom_cache(m)%m(3**m, mesh%n_elems))
      discrete_vmom_cache(m)%m = 0.0_DOUBLE
    end do

    do i = 1, mesh%n_elems
      n_v = mesh%elem(i)%n_vert
      do m = 2, max_order
        do kv = 1, n_v
          id_v = mesh%elem(i)%vert(kv)
          dx_v = mesh%vert(id_v)%coord - mesh%elem(i)%coord
          do t = 0, 3**m - 1
            idx = t
            prodval = 1.0_DOUBLE
            do r = 1, m
              dloc = mod(idx, 3) + 1
              idx = idx / 3
              prodval = prodval * dx_v(dloc)
            end do
            discrete_vmom_cache(m)%m(t+1, i) = discrete_vmom_cache(m)%m(t+1, i) + prodval
          end do
        end do
        discrete_vmom_cache(m)%m(:, i) = discrete_vmom_cache(m)%m(:, i) / real(n_v, kind=DOUBLE)
      end do
    end do

    discrete_vmom_cache_n_elems = mesh%n_elems
    discrete_vmom_cache_max_order = max_order
  end subroutine ensure_discrete_vertex_moment_cache

  ! Corrects, in place, a cell-level quantity D^(q)_cell (built by recombine_derivative_regression,
  ! which already removes the m=1 term via its affine fit) using every available higher blended
  ! derivative D^(q+2)_cell,...,D^(k_max)_cell -- same Taylor argument and loop structure as
  ! apply_local_taylor_correction, but for the DISCRETE corner-vertex average discrete_vmom_cache
  ! instead of a continuous volume integral, and starting at m=2 since the regression already
  ! handles m=1. Verified to close the vertex-to-cell recombination gap to ~1e-10.
  subroutine apply_discrete_vertex_moment_correction(mesh, d, nc_in, q, k_max, q_val, dfield_cell)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in, q, k_max
    real(kind=DOUBLE), dimension(nc_in*d**q, mesh%n_elems), intent(inout) :: q_val
    type(derivative_field_type), dimension(q+2:k_max), intent(in) :: dfield_cell

    integer(kind=ENTIER) :: m
    real(kind=DOUBLE), dimension(:, :), allocatable :: corr

    if (k_max < q+2) return
    call ensure_discrete_vertex_moment_cache(mesh, k_max-q)

    do m = 2, k_max - q
      allocate(corr(nc_in*d**q, mesh%n_elems))
      call contract_last_indices(d, nc_in, mesh%n_elems, q+m, m, dfield_cell(q+m)%val, &
        discrete_vmom_cache(m)%m, corr)
      q_val = q_val - corr / fact_ratio(m, m)
      deallocate(corr)
    end do
  end subroutine apply_discrete_vertex_moment_correction

end module arbitrary_high_order_module
