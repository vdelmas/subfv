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
  public :: apply_derivative_node_correction
  public :: scatter_bias_to_cells
  public :: compute_gradient_node_bias
  public :: compute_node_derivative_bias
  public :: ensure_ls_grad_h1_cache
  public :: ensure_gg_gradient_h1_cache
  public :: recombine_derivative_regression
  public :: compute_2exact_hessian_aho_cls
  public :: apply_discrete_vertex_moment_correction
  public :: invalidate_geometry_caches
  public :: compute_derivatives_aho_cls
  public :: compute_next_order_node_aho_cls
  public :: blend_node_to_cell_aho_cls
  public :: compute_derivatives_cls_classic
  public :: compute_derivatives_cls_classic_order3
  public :: ensure_cell_ls_grad_mat_cache
  public :: compute_cell_ls_grad
  public :: compute_cell_cls_2exact_hessian
  public :: apply_cell_cls_hess_operator
  public :: compute_cell_cls_3exact_third
  public :: apply_cell_cls_eq13_correction
  public :: third_red_to_full27
  public :: grad_of_geom_field_at
  public :: hess_of_geom_field_at
  public :: third_of_geom_field_at
  public :: compute_cell_cls_4exact_fourth
  public :: canon10_index
  public :: canon15_index
  public :: ensure_cell_moment_cache
  public :: set_wall_mirror_data
  public :: mirror_wall_vec_start

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
  ! Debug-only: when > 0, compute_node_derivative_bias prints intermediate norms for this vertex.
  integer(kind=ENTIER), public :: debug_bias_vertex = 0
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

  ! Wall-tangent fit for boundary-vertex reconstruction (2026-09-20; supersedes an earlier
  ! mirror-ghost-cell attempt that turned out numerically unstable -- reflecting sub-face geometry
  ! across the wall could drive the LS/GG fit matrix arbitrarily close to singular, producing huge
  ! finite gradients that only showed up as NaN several RK stages downstream, well after any
  ! dt/t-based health check would catch it): the caller (euler_ho_module) computes, once per mesh,
  ! a per-vertex outward wall normal from the vertex's touching WALL-type boundary faces and
  ! registers it here via set_wall_mirror_data. mirror_wall_vec_start=0 disables the whole
  ! mechanism (default); any nonzero value enables it (the specific value is no longer used to pick
  ! which nc_in components are velocity -- this approach never touches phi's components at all).
  ! Enabled, a wall vertex's fit is restricted to the plane tangent to its wall normal (1 tangent
  ! direction for boundary_2d, 2 otherwise) instead of the global x[,y[,z]] directions: a small,
  ! well-posed problem using only the well-resolved along-wall neighbor spread, leaving the
  ! (poorly-resolved, one-sided) wall-normal gradient component at exactly 0 rather than
  ! extrapolating it. See compute_nodal_derivative_at_vertex (LS) and
  ! ensure_green_gauss_mat_cache/compute_nodal_derivative_at_vertex_green_gauss (GG). Only
  ! consulted at deriv_order=1 (nc_in=5, primitives) -- compute_next_order_derivative's boundary
  ! loop keeps the earlier phantom-zero-gradient hack for deriv_order=2/3 (hess/third) and for any
  ! boundary vertex without a clean single-wall normal (mixed wall+inflow/outflow vertices, or
  ! wall_mirror_valid=.false.).
  integer(kind=ENTIER), save :: mirror_wall_vec_start = 0
  integer(kind=ENTIER), save :: wall_mirror_n_vert = -1
  real(kind=DOUBLE), dimension(:, :), allocatable, save :: wall_mirror_norm
  logical, dimension(:), allocatable, save :: wall_mirror_valid

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

  ! Persistent scratch buffers for compute_next_order_derivative's WENO-blend accumulation --
  ! previously allocated and zeroed fresh on EVERY call (3x per RK stage: grad/hess/third, nc_out
  ! up to 135), a real allocator + first-touch-page-fault cost on large meshes (profiled: memset
  ! alone was ~9% of total runtime on a 6000-vertex tet mesh). Sized once to the largest nc_out
  ! actually used (135, from third's nc_in=45, d=3) and reused via a (1:nc_out,:) slice --
  ! assumed-shape dummies in accumulate_weno_contribution/scatter_weno_weighted/
  ! rescue_zero_weight_cell accept that slice by descriptor, no copy.
  integer(kind=ENTIER), parameter :: WENO_BUF_MAX_NC_OUT = 135
  real(kind=DOUBLE), dimension(:, :), allocatable, save :: weno_num_buf
  real(kind=DOUBLE), dimension(:), allocatable, save :: weno_den_buf
  real(kind=DOUBLE), dimension(:, :), allocatable, save :: dphi_v_cache_buf
  logical, dimension(:), allocatable, save :: valid_cache_buf
  real(kind=DOUBLE), dimension(:), allocatable, save :: oi_cache_buf
  integer(kind=ENTIER), save :: weno_buf_n_elems = -1
  integer(kind=ENTIER), save :: weno_buf_n_vert = -1

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

  ! LS analogue of gg_mat_inv_cache: per-vertex weighted-LS fit matrix, reproducing
  ! compute_nodal_derivative_at_vertex's own operator exactly (weight=1/|dx|^2, basis={1,x[,y[,z]]},
  ! same dynamic active-dimension dropping) so ls_grad_h1_cache below applies the EXACT SAME
  ! operator aho_ls itself uses -- required for Haider's eq. 13 correction to cancel the bias
  ! exactly. Skips boundary/wall vertices (matches gg_mat_inv_cache; eq. 13 only ever corrects
  ! interior-vertex-fed derivatives). Pure geometry, cached once per mesh.
  integer(kind=ENTIER), save :: ls_mat_cache_n_vert = -1
  logical, save :: ls_mat_cache_boundary_2d = .false.
  real(kind=DOUBLE), dimension(:, :, :), allocatable, save :: ls_mat_inv_cache ! (4,4,n_vert)
  integer(kind=ENTIER), dimension(:, :), allocatable, save :: ls_active_dim_cache ! (3,n_vert)
  integer(kind=ENTIER), dimension(:), allocatable, save :: ls_n_active_cache
  logical, dimension(:), allocatable, save :: ls_mat_valid_cache

  ! LS analogue of gg_grad_h1_cache: H_m^(1)(v), m=2..max_order, built by applying the SAME
  ! ls_mat_inv_cache operator to the geometric field {z_vK^(m)}_K (shifted_cell_moment_full)
  ! instead of phi. Pure geometry, cached once per mesh.
  integer(kind=ENTIER), save :: ls_grad_h1_cache_n_vert = -1
  integer(kind=ENTIER), save :: ls_grad_h1_cache_max_order = 0
  logical, save :: ls_grad_h1_cache_boundary_2d = .false.
  type(cell_moment_ptr_type), dimension(:), allocatable, save :: ls_grad_h1_cache ! indexed 2:max_order

  ! Cell-to-cell 1-exact gradient LS fit matrix (Haider, Croisille & Courbet 2011, Definition 1,
  ! step 1): per cell alpha, unweighted normal-equations inverse over the cell's own
  ! mesh%elem(alpha)%neigh_by_vert stencil. Pure geometry, cached once per mesh (boundary_2d baked
  ! in as for gg_mat_inv_cache). This is the genuinely cell-centered CLS reference -- no vertex
  ! fit, no cell blend anywhere -- distinct from aho_gg/aho_ls's vertex-fit+WENO-scatter chain.
  integer(kind=ENTIER), save :: cell_ls_mat_cache_n_elems = -1
  logical, save :: cell_ls_mat_cache_boundary_2d = .false.
  real(kind=DOUBLE), dimension(:, :, :), allocatable, save :: cell_ls_mat_inv_cache

  ! Cell-to-cell analogue of gg_grad_h1_cache: H_m^(1)(alpha) = cell_ls_mat_inv_cache(alpha)
  ! applied to the geometric field {z_{alpha,beta}^(m)-z_{alpha,alpha}^(m)}_beta
  ! (shifted_cell_moment_full) instead of phi, over the SAME neigh_by_vert stencil as the grad
  ! fit. Reproduces Haider's eq. 15-16 J-operator correction term faithfully (their w_beta^(k|k)
  ! applied to z_beta^(k+1), for k=1). Pure geometry, cached once per mesh.
  integer(kind=ENTIER), save :: cell_ls_h1_cache_n_elems = -1
  integer(kind=ENTIER), save :: cell_ls_h1_cache_max_order = 0
  logical, save :: cell_ls_h1_cache_boundary_2d = .false.
  type(cell_moment_ptr_type), dimension(:), allocatable, save :: cell_ls_h1_cache ! indexed 2:max_order

  ! Per-cell X's OWN hess operator applied to the geometric field z_X^(3) instead of phi
  ! (Haider's w_beta^(2|2)[z_beta^(3)], eq. 15-16 m=2 case) -- see
  ! ensure_cell_cls_hess_of_z3_cache's own header comment for the derivation. Reduced 6-component
  ! hess basis x 27-component z basis (162), cached once per mesh.
  integer(kind=ENTIER), save :: cell_cls_hess_of_z3_n_elems = -1
  logical, save :: cell_cls_hess_of_z3_boundary_2d = .false.
  real(kind=DOUBLE), dimension(:, :), allocatable, save :: cell_cls_hess_of_z3_cache

  ! One timestamped phase of one order's work, for the blocking-vs-overlap timeline figure.
  type :: timeline_event_type
    integer(kind=ENTIER) :: order
    character(len=20)    :: phase
    real(kind=DOUBLE)    :: t0, t1
  end type timeline_event_type

contains

  ! Drops every cached quantity that was computed from vertex COORDINATES, so
  ! the next reconstruction call rebuilds them from the current geometry.
  !
  ! Each cache above is guarded only by a size check (cache_n_vert ==
  ! mesh%n_vert, cache_n_elems == mesh%n_elems). That is sound for a fixed
  ! mesh, but a solver that MOVES nodes changes neither count, so without this
  ! call every cache silently survives the move and the reconstruction is then
  ! built on the pre-move geometry -- wrong at order >= 2, and silent.
  !
  ! Call it after move_mesh/compute_geometry_mesh, before the next
  ! reconstruction. Topology caches (neigh/vv_neigh) are reset too: they are
  ! connectivity-only and a pure move leaves them valid, but resetting is cheap
  ! next to a full rebuild and keeps this routine a single honest statement
  ! ("geometry changed") rather than a list the caller must keep in sync.
  subroutine invalidate_geometry_caches()
    implicit none

    neigh_cache_n_vert = -1
    vv_neigh_cache_n_vert = -1
    m2_cache_n_elems = -1
    gg_mat_cache_n_vert = -1
    wall_mirror_n_vert = -1
    cell_moment_cache_n_elems = -1
    cell_moment_cache_max_order = 0
    discrete_vmom_cache_n_elems = -1
    discrete_vmom_cache_max_order = 0
    gg_grad_h1_cache_n_vert = -1
    gg_grad_h1_cache_max_order = 0
    ls_mat_cache_n_vert = -1
    ls_grad_h1_cache_n_vert = -1
    ls_grad_h1_cache_max_order = 0
    cell_ls_mat_cache_n_elems = -1
    cell_ls_h1_cache_n_elems = -1
    cell_ls_h1_cache_max_order = 0
    cell_cls_hess_of_z3_n_elems = -1
  end subroutine invalidate_geometry_caches

  ! Registers the per-vertex wall normal used by the mirror-ghost-cell fix (see the cache block's
  ! header comment above). Call once per mesh, before the first compute_next_order_derivative/
  ! ensure_green_gauss_mat_cache call that should see it (mesh%n_vert-sized re-registration is
  ! cheap; the GG cache only rebuilds when mesh%n_vert or boundary_2d changes, so calling this
  ! AFTER ensure_green_gauss_mat_cache has already run for this mesh size would not retroactively
  ! rebuild it -- register before the first reconstruction call of a run).
  subroutine set_wall_mirror_data(mesh, norm_v, valid_v)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(:, :), intent(in) :: norm_v
    logical, dimension(:), intent(in) :: valid_v

    if (allocated(wall_mirror_norm))  deallocate(wall_mirror_norm)
    if (allocated(wall_mirror_valid)) deallocate(wall_mirror_valid)
    allocate(wall_mirror_norm(3, mesh%n_vert))
    allocate(wall_mirror_valid(mesh%n_vert))
    wall_mirror_norm  = norm_v
    wall_mirror_valid = valid_v
    wall_mirror_n_vert = mesh%n_vert
  end subroutine set_wall_mirror_data

  ! .true. iff v is a registered wall-mirror vertex (a valid unit normal is available and mirroring is enabled).
  pure function is_wall_mirror_vertex(v, n_vert_mesh) result(r)
    implicit none

    integer(kind=ENTIER), intent(in) :: v, n_vert_mesh
    logical :: r

    r = mirror_wall_vec_start > 0 .and. wall_mirror_n_vert == n_vert_mesh
    if (r) r = allocated(wall_mirror_valid)
    if (r) r = wall_mirror_valid(v)
  end function is_wall_mirror_vertex

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

    nc_out = nc_in*d
    deriv_order_eff = 1
    if (present(deriv_order)) deriv_order_eff = deriv_order

    ! nc_out cycles through a small fixed set of values (15/45/135 for grad/hess/third) every RK
    ! stage; resized once to the largest ever requested (135 in production; a test-only caller in
    ! arbitrary_high_order_main.F90 could ask for more, hence the max() rather than a hardcoded
    ! bound) and reused via a (1:nc_out,:) ASSOCIATE alias -- unlike POINTER association this adds
    ! no indirection (the compiler substitutes the slice expression at compile time), so it avoids
    ! both the allocate/deallocate/first-touch cost this used to pay every single call AND any
    ! aliasing-analysis penalty a real pointer would add.
    if (weno_buf_n_elems /= mesh%n_elems .or. weno_buf_n_vert /= mesh%n_vert &
        .or. (allocated(weno_num_buf) .and. size(weno_num_buf, 1) < nc_out)) then
      if (allocated(weno_num_buf)) then
        deallocate(weno_num_buf, weno_den_buf, dphi_v_cache_buf, valid_cache_buf, oi_cache_buf)
      end if
      allocate(weno_num_buf(max(WENO_BUF_MAX_NC_OUT, nc_out), mesh%n_elems))
      allocate(weno_den_buf(mesh%n_elems))
      allocate(dphi_v_cache_buf(max(WENO_BUF_MAX_NC_OUT, nc_out), mesh%n_vert))
      allocate(valid_cache_buf(mesh%n_vert))
      allocate(oi_cache_buf(mesh%n_vert))
      weno_buf_n_elems = mesh%n_elems
      weno_buf_n_vert = mesh%n_vert
    end if

    associate (weno_num => weno_num_buf(1:nc_out, :), weno_den => weno_den_buf, &
        dphi_v_cache => dphi_v_cache_buf(1:nc_out, :), valid_cache => valid_cache_buf, &
        oi_cache => oi_cache_buf)
    weno_num = 0.0_DOUBLE
    weno_den = 0.0_DOUBLE

    ! A boundary vertex is skipped (one-sided neighbor gather) unless .not. boundary_2d, where skipping it would starve every vertex.
    ! 2026-09-20: at deriv_order=1 (the primitive gradient), a registered wall-mirror vertex
    ! (is_wall_mirror_vertex) instead gets a symmetrized fit via a mirrored ghost-neighbor
    ! contribution (LS: wall_norm passed into compute_nodal_derivative_at_vertex; GG: baked into
    ! ensure_green_gauss_mat_cache's own per-vertex matrix) -- see the wall-mirror cache block's
    ! header comment. Every other boundary vertex (deriv_order=2/3's hess/third, or any boundary
    ! vertex without a clean single-wall normal) keeps the earlier phantom-zero-gradient hack:
    ! scatter a phantom near-zero gradient (dphi_v~1e-16, oi_v=0) with the standard
    ! weight=omega_p/(eps+OI^p) formula -- OI=0 makes this phantom sample dominate the blend for
    ! any touching cell, pulling the reconstructed derivative toward flat/first-order near walls.
    do id_vert = 1, mesh%n_vert
      if (mesh%vert(id_vert)%is_bound) then
        if (deriv_order_eff == 1 .and. is_wall_mirror_vertex(id_vert, mesh%n_vert)) then
          if (use_green_gauss) then
            call accumulate_weno_contribution(mesh, d, nc_in, boundary_2d, id_vert, &
              phi, dphi_v_cache, valid_cache, oi_cache, weno_num, weno_den, &
              deriv_order=deriv_order_eff)
          else
            call accumulate_weno_contribution(mesh, d, nc_in, boundary_2d, id_vert, &
              phi, dphi_v_cache, valid_cache, oi_cache, weno_num, weno_den, &
              deriv_order=deriv_order_eff, wall_norm=wall_mirror_norm(:, id_vert))
          end if
          cycle
        end if
        block
          real(kind=DOUBLE), dimension(nc_out) :: dphi_v_zero
          dphi_v_zero = 1.0e-16_DOUBLE
          dphi_v_cache(:, id_vert) = dphi_v_zero
          valid_cache(id_vert) = .true.
          oi_cache(id_vert) = 0.0_DOUBLE
          call scatter_weno_weighted(mesh, id_vert, dphi_v_zero, 0.0_DOUBLE, weno_num, weno_den, &
            eps_in=1.0e-16_DOUBLE)
        end block
        cycle
      end if
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
    end associate
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
      phi, dphi_v_cache, valid_cache, oi_cache, weno_num, weno_den, skip_cell, deriv_order, wall_norm)
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
    ! LS-only wall-mirror normal for this vertex (GG's mirror handling is entirely baked into
    ! ensure_green_gauss_mat_cache/compute_nodal_derivative_at_vertex_green_gauss, keyed off the
    ! module-level wall_mirror_norm/valid, so it needs no argument here).
    real(kind=DOUBLE), dimension(3), intent(in), optional :: wall_norm

    real(kind=DOUBLE), dimension(nc_in*d) :: dphi_v
    logical :: valid
    real(kind=DOUBLE) :: oi_v, eps_here
    integer(kind=ENTIER) :: deriv_order_eff

    if (use_green_gauss) then
      call compute_nodal_derivative_at_vertex_green_gauss(mesh, d, nc_in, boundary_2d, &
        id_vert, phi, dphi_v, valid, oi_v)
    else if (present(wall_norm)) then
      call compute_nodal_derivative_at_vertex(mesh, d, nc_in, boundary_2d, &
        id_vert, phi, dphi_v, valid, oi_v, wall_norm)
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
      id_vert, phi, dphi_v, valid, oi_v, wall_norm)
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
    ! Wall-tangent fit (2026-09-20, replaces an earlier mirror-ghost-cell attempt that turned out
    ! numerically unstable): when present (and mirror_wall_vec_start>0), the fit basis is restricted
    ! to the plane tangent to this unit wall normal (1 tangent direction for boundary_2d, 2
    ! otherwise) instead of the global x[,y[,z]] directions -- a small, well-posed problem using
    ! only the well-resolved along-wall neighbor spread. The wall-normal gradient component is never
    ! solved for and is left at exactly 0 rather than extrapolated from a one-sided stencil.
    real(kind=DOUBLE), dimension(3), intent(in), optional :: wall_norm

    integer(kind=ENTIER), parameter :: max_basis = 4
    real(kind=DOUBLE), parameter :: rel_spread_tol = 1.0e-8_DOUBLE
    integer(kind=ENTIER) :: n_basis, n_cand, n_active, j, a, b, i1, id_elem, n_neigh, bdir, n_tan
    integer(kind=ENTIER), dimension(3) :: active_dim
    integer(kind=ENTIER), dimension(:), allocatable :: neigh
    integer(kind=ENTIER), dimension(max_basis) :: ipiv
    real(kind=DOUBLE), dimension(3) :: dx, dmin, dmax, spread, bp, refv
    real(kind=DOUBLE), dimension(3, 2) :: tvec
    real(kind=DOUBLE) :: weight, max_spread, weight_sum, phi_scale2, dcoord
    real(kind=DOUBLE), dimension(max_basis) :: basis
    real(kind=DOUBLE), dimension(max_basis, max_basis) :: mat
    real(kind=DOUBLE), dimension(max_basis, nc_in) :: rhs
    real(kind=DOUBLE), dimension(nc_in) :: predicted, resid_sq, phi_sq_sum
    logical :: is_wall

    is_wall = present(wall_norm) .and. mirror_wall_vec_start > 0
    bp = 0.0_DOUBLE
    n_tan = 0
    if (is_wall) then
      if (dot_product(wall_norm, wall_norm) < 1.0e-24_DOUBLE) then
        is_wall = .false.
      else
        bp = wall_norm / sqrt(dot_product(wall_norm, wall_norm))
      end if
    end if

    if (is_wall) then
      ! Orthonormal tangent basis {t1[,t2]} spanning the plane perp to bp (Gram-Schmidt off a
      ! reference axis not near-parallel to bp; boundary_2d keeps t1 in-plane, z=0, since bp itself
      ! has bp(3)=0 for any side wall of a z-extruded 2D mesh).
      n_tan = merge(1_ENTIER, 2_ENTIER, boundary_2d)
      if (abs(bp(1)) < 0.9_DOUBLE) then
        refv = (/1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE/)
      else
        refv = (/0.0_DOUBLE, 1.0_DOUBLE, 0.0_DOUBLE/)
      end if
      tvec(:, 1) = refv - dot_product(refv, bp) * bp
      tvec(:, 1) = tvec(:, 1) / sqrt(dot_product(tvec(:, 1), tvec(:, 1)))
      if (n_tan == 2) then
        tvec(1, 2) = bp(2)*tvec(3,1) - bp(3)*tvec(2,1)
        tvec(2, 2) = bp(3)*tvec(1,1) - bp(1)*tvec(3,1)
        tvec(3, 2) = bp(1)*tvec(2,1) - bp(2)*tvec(1,1)
      end if
    end if

    call ensure_neighbor_cache(mesh)
    n_neigh = neigh_cache_start(id_vert+1) - neigh_cache_start(id_vert)
    allocate(neigh(n_neigh))
    neigh = neigh_cache_list(neigh_cache_start(id_vert):neigh_cache_start(id_vert+1)-1)

    if (is_wall) then
      n_cand = n_tan
    else
      n_cand = merge(2_ENTIER, 3_ENTIER, boundary_2d)
    end if
    dmin(1:n_cand) = huge(1.0_DOUBLE)
    dmax(1:n_cand) = -huge(1.0_DOUBLE)
    do j = 1, n_neigh
      dx = mesh%elem(neigh(j))%coord - mesh%vert(id_vert)%coord
      do a = 1, n_cand
        if (is_wall) then
          dcoord = dot_product(dx, tvec(:, a))
        else
          dcoord = dx(a)
        end if
        dmin(a) = min(dmin(a), dcoord)
        dmax(a) = max(dmax(a), dcoord)
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
        if (is_wall) then
          basis(1+a) = dot_product(dx, tvec(:, active_dim(a)))
        else
          basis(1+a) = dx(active_dim(a))
        end if
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
        if (is_wall) then
          predicted = predicted + rhs(1+a, :) * dot_product(dx, tvec(:, active_dim(a)))
        else
          predicted = predicted + rhs(1+a, :) * dx(active_dim(a))
        end if
      end do
      resid_sq = resid_sq + weight * (phi(:, id_elem) - predicted)**2
    end do
    phi_scale2 = maxval(phi_sq_sum) / max(weight_sum, 1.0e-300_DOUBLE)
    oi_v = sqrt((maxval(resid_sq) / max(weight_sum, 1.0e-300_DOUBLE)) &
      / max(phi_scale2, 1.0e-300_DOUBLE))
    oi_v = max(oi_v, (max_spread * sum(rhs(2:n_basis, :)**2) / max(phi_scale2, 1.0e-300_DOUBLE)) &
      / grad_norm_derate)

    ! rhs(1+a,i1)=d(phi_i1)/dx_active_dim(a); flattened component-fast/direction-slow, inactive
    ! directions left at 0. is_wall: active_dim(a) instead indexes a SURVIVING TANGENT direction
    ! (tvec(:,active_dim(a))), expanded back into global x[,y[,z]] components -- the wall-normal
    ! component is never touched and stays 0 (dphi_v was zeroed above).
    if (is_wall) then
      do a = 1, n_active
        do i1 = 1, nc_in
          do bdir = 1, d
            dphi_v((bdir-1)*nc_in + i1) = dphi_v((bdir-1)*nc_in + i1) &
              + rhs(1+a, i1) * tvec(bdir, active_dim(a))
          end do
        end do
      end do
    else
      do a = 1, n_active
        do i1 = 1, nc_in
          dphi_v((active_dim(a)-1)*nc_in + i1) = rhs(1+a, i1)
        end do
      end do
    end if
  end subroutine compute_nodal_derivative_at_vertex

  ! Builds gg_mat_inv_cache/gg_flux_*_cache/gg_oi_hlocal_cache (pure geometry) once per mesh; uses the LAPACK SVD pseudo-inverse unconditionally since inversion is now a one-time cost.
  subroutine ensure_green_gauss_mat_cache(mesh, boundary_2d)
    use linear_solver_module, only: pseudo_inverse_inplace_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d

    integer(kind=ENTIER) :: v, i, j, a, b, k, n_cand
    integer(kind=ENTIER) :: id_sub_elem, id_elem, id_sub_face, total_pairs
    real(kind=DOUBLE), dimension(3) :: dx, norm, dminn, dmaxn, bp, refv
    real(kind=DOUBLE), dimension(3, 3) :: mat, mat_inv
    real(kind=DOUBLE), dimension(3, 2) :: tvec, mt
    real(kind=DOUBLE), dimension(2, 2) :: tam, tam_inv
    integer(kind=ENTIER) :: n_tan
    logical :: is_wm

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

    ! A registered wall vertex (see the wall-mirror cache block's header comment -- name kept for
    ! the shared is_wall_mirror_vertex helper, though this GG path no longer mirrors anything) gets
    ! a valid matrix too, built from real geometry only, same as any interior vertex.
    total_pairs = 0
    do v = 1, mesh%n_vert
      is_wm = is_wall_mirror_vertex(v, mesh%n_vert)
      if (mesh%vert(v)%is_bound .and. .not. is_wm) cycle
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
      is_wm = is_wall_mirror_vertex(v, mesh%n_vert)
      if (mesh%vert(v)%is_bound .and. .not. is_wm) then
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

      if (is_wm) then
        ! Wall-tangent fit (see compute_nodal_derivative_at_vertex's own header comment for the
        ! LS twin of this same idea): restrict the flux-matching unknown to grad=T@grad_t (T's
        ! columns an orthonormal tangent basis), solved by least squares --
        ! grad_t=(T'mat'mat T)^-1 T'mat' grad_raw -- then folded into one effective 3x3 operator
        ! mat_inv_eff=T@(T'mat'mat T)^-1@T'mat' so compute_nodal_derivative_at_vertex_green_gauss's
        ! plain grad_true=mat_inv_eff@grad_raw needs no changes. The wall-normal component is never
        ! solved for, so it comes out exactly 0.
        bp = wall_mirror_norm(:, v)
        n_tan = merge(1_ENTIER, 2_ENTIER, boundary_2d)
        if (abs(bp(1)) < 0.9_DOUBLE) then
          refv = (/1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE/)
        else
          refv = (/0.0_DOUBLE, 1.0_DOUBLE, 0.0_DOUBLE/)
        end if
        tvec(:, 1) = refv - dot_product(refv, bp) * bp
        tvec(:, 1) = tvec(:, 1) / sqrt(dot_product(tvec(:, 1), tvec(:, 1)))
        if (n_tan == 2) then
          tvec(1, 2) = bp(2)*tvec(3,1) - bp(3)*tvec(2,1)
          tvec(2, 2) = bp(3)*tvec(1,1) - bp(1)*tvec(3,1)
          tvec(3, 2) = bp(1)*tvec(2,1) - bp(2)*tvec(1,1)
        end if

        mt(:, 1:n_tan) = matmul(mat, tvec(:, 1:n_tan))
        tam(1:n_tan, 1:n_tan) = matmul(transpose(mt(:, 1:n_tan)), mt(:, 1:n_tan))
        ! Tikhonov ridge, relative to mat's own scale: a genuinely near-singular tam (t happens to
        ! sit close to mat's null space too, seen on a handful of mesh-quality-degenerate wall
        ! vertices) then damps smoothly toward a near-0 tangential gradient there instead of the
        ! plain SVD pseudo-inverse's huge-but-finite blowup (1/tiny-sigma) -- confirmed by direct
        ! measurement to reach ~1e6-1e7 on this mesh's worst vertices, enough to detonate the whole
        ! solve within ~100 iterations despite looking like a merely "ill-conditioned", not exactly
        ! singular, 1x1/2x2 system.
        do a = 1, n_tan
          tam(a, a) = tam(a, a) + 1.0e-6_DOUBLE * sum(mat**2)
        end do
        tam_inv(1:n_tan, 1:n_tan) = tam(1:n_tan, 1:n_tan)
        call pseudo_inverse_inplace_lapack(n_tan, tam_inv(1:n_tan, 1:n_tan))
        gg_mat_inv_cache(:, :, v) = matmul(tvec(:, 1:n_tan), &
          matmul(tam_inv(1:n_tan, 1:n_tan), transpose(mt(:, 1:n_tan))))
      else
        ! mat's third column is identically zero for boundary_2d (singular by construction); the SVD handles that and any near-degenerate 3D element.
        mat_inv = mat
        call pseudo_inverse_inplace_lapack(3_ENTIER, mat_inv)
        if (boundary_2d) mat_inv(3, :) = 0.0_DOUBLE
        gg_mat_inv_cache(:, :, v) = mat_inv
      end if
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

  ! Builds ls_mat_inv_cache/ls_active_dim_cache/ls_n_active_cache/ls_mat_valid_cache: per-vertex
  ! weighted-LS fit matrix reproducing compute_nodal_derivative_at_vertex's own operator exactly
  ! (weight=1/|dx|^2, basis={1,x[,y[,z]]}, same dynamic active-dimension dropping by spread), minus
  ! the wall-tangent special case (wall/boundary vertices are simply left invalid, matching
  ! gg_mat_valid_cache's own convention -- eq. 13's correction only ever applies at interior
  ! vertices).
  subroutine ensure_ls_mat_cache(mesh, boundary_2d)
    use linear_solver_module, only: pseudo_inverse_inplace_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d

    real(kind=DOUBLE), parameter :: rel_spread_tol = 1.0e-8_DOUBLE
    integer(kind=ENTIER) :: v, j, id_elem, n_neigh, n_cand, n_active, n_basis, a, b
    integer(kind=ENTIER), dimension(:), allocatable :: neigh
    integer(kind=ENTIER), dimension(3) :: active_dim
    real(kind=DOUBLE), dimension(3) :: dx, dmin, dmax, spread_v
    real(kind=DOUBLE) :: weight, max_spread
    real(kind=DOUBLE), dimension(4) :: basis
    real(kind=DOUBLE), dimension(4, 4) :: mat, mat_inv

    if (ls_mat_cache_n_vert == mesh%n_vert .and. (ls_mat_cache_boundary_2d .eqv. boundary_2d)) return

    call ensure_neighbor_cache(mesh)

    if (allocated(ls_mat_inv_cache)) deallocate(ls_mat_inv_cache)
    if (allocated(ls_active_dim_cache)) deallocate(ls_active_dim_cache)
    if (allocated(ls_n_active_cache)) deallocate(ls_n_active_cache)
    if (allocated(ls_mat_valid_cache)) deallocate(ls_mat_valid_cache)
    allocate(ls_mat_inv_cache(4, 4, mesh%n_vert))
    allocate(ls_active_dim_cache(3, mesh%n_vert))
    allocate(ls_n_active_cache(mesh%n_vert))
    allocate(ls_mat_valid_cache(mesh%n_vert))
    ls_mat_inv_cache = 0.0_DOUBLE
    ls_active_dim_cache = 0
    ls_n_active_cache = 0
    ls_mat_valid_cache = .false.

    n_cand = merge(2_ENTIER, 3_ENTIER, boundary_2d)

    do v = 1, mesh%n_vert
      if (mesh%vert(v)%is_bound) cycle

      n_neigh = neigh_cache_start(v+1) - neigh_cache_start(v)
      allocate(neigh(n_neigh))
      neigh = neigh_cache_list(neigh_cache_start(v):neigh_cache_start(v+1)-1)

      dmin(1:n_cand) = huge(1.0_DOUBLE)
      dmax(1:n_cand) = -huge(1.0_DOUBLE)
      do j = 1, n_neigh
        dx = mesh%elem(neigh(j))%coord - mesh%vert(v)%coord
        do a = 1, n_cand
          dmin(a) = min(dmin(a), dx(a))
          dmax(a) = max(dmax(a), dx(a))
        end do
      end do
      spread_v(1:n_cand) = dmax(1:n_cand) - dmin(1:n_cand)
      max_spread = maxval(spread_v(1:n_cand))

      n_active = 0
      do a = 1, n_cand
        if (spread_v(a) > rel_spread_tol * max(max_spread, 1.0e-300_DOUBLE)) then
          n_active = n_active + 1
          active_dim(n_active) = a
        end if
      end do
      n_basis = 1 + n_active

      if (n_active == 0) then
        deallocate(neigh)
        cycle
      end if

      mat = 0.0_DOUBLE
      do j = 1, n_neigh
        id_elem = neigh(j)
        dx = mesh%elem(id_elem)%coord - mesh%vert(v)%coord
        weight = 1.0_DOUBLE / max(dot_product(dx, dx), 1.0e-24_DOUBLE)
        basis(1) = 1.0_DOUBLE
        do a = 1, n_active
          basis(1+a) = dx(active_dim(a))
        end do
        do a = 1, n_basis
          do b = 1, n_basis
            mat(a, b) = mat(a, b) + weight * basis(a) * basis(b)
          end do
        end do
      end do

      mat_inv(1:n_basis, 1:n_basis) = mat(1:n_basis, 1:n_basis)
      call pseudo_inverse_inplace_lapack(n_basis, mat_inv(1:n_basis, 1:n_basis))

      ls_mat_inv_cache(1:n_basis, 1:n_basis, v) = mat_inv(1:n_basis, 1:n_basis)
      ls_active_dim_cache(1:n_active, v) = active_dim(1:n_active)
      ls_n_active_cache(v) = n_active
      ls_mat_valid_cache(v) = .true.

      deallocate(neigh)
    end do

    ls_mat_cache_n_vert = mesh%n_vert
    ls_mat_cache_boundary_2d = boundary_2d
  end subroutine ensure_ls_mat_cache

  ! LS analogue of ensure_gg_gradient_h1_cache: H_m^(1)(v), m=2..max_order, built by applying the
  ! SAME ls_mat_inv_cache weighted-LS operator to the geometric field {z_vK^(m)}_K
  ! (shifted_cell_moment_full) instead of phi -- only the gradient rows of the fit (dropping the
  ! constant/intercept row) are kept, matching compute_nodal_derivative_at_vertex's own dphi_v.
  ! Pure geometry, cached once per mesh.
  subroutine ensure_ls_grad_h1_cache(mesh, boundary_2d, max_order)
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d
    integer(kind=ENTIER), intent(in) :: max_order

    integer(kind=ENTIER) :: v, j, id_elem, n_neigh, n_active, n_basis, m, tt, n_tt, a
    real(kind=DOUBLE), dimension(3) :: dx
    real(kind=DOUBLE) :: weight
    real(kind=DOUBLE), dimension(:), allocatable :: z_k
    real(kind=DOUBLE), dimension(:, :), allocatable :: rhs, sol

    call ensure_ls_mat_cache(mesh, boundary_2d)
    call ensure_cell_moment_cache(mesh, max_order)
    call ensure_neighbor_cache(mesh)

    if (ls_grad_h1_cache_n_vert == mesh%n_vert .and. ls_grad_h1_cache_max_order >= max_order &
        .and. (ls_grad_h1_cache_boundary_2d .eqv. boundary_2d)) return

    if (allocated(ls_grad_h1_cache)) deallocate(ls_grad_h1_cache)
    allocate(ls_grad_h1_cache(2:max_order))
    do m = 2, max_order
      allocate(ls_grad_h1_cache(m)%m(3*3**m, mesh%n_vert))
      ls_grad_h1_cache(m)%m = 0.0_DOUBLE
    end do

    do m = 2, max_order
      n_tt = 3**m
      allocate(z_k(n_tt))
      do v = 1, mesh%n_vert
        if (.not. ls_mat_valid_cache(v)) cycle
        n_active = ls_n_active_cache(v)
        n_basis = 1 + n_active

        n_neigh = neigh_cache_start(v+1) - neigh_cache_start(v)
        allocate(rhs(n_basis, n_tt), sol(n_basis, n_tt))
        rhs = 0.0_DOUBLE
        do j = 1, n_neigh
          id_elem = neigh_cache_list(neigh_cache_start(v)+j-1)
          dx = mesh%elem(id_elem)%coord - mesh%vert(v)%coord
          weight = 1.0_DOUBLE / max(dot_product(dx, dx), 1.0e-24_DOUBLE)
          call shifted_cell_moment_full(id_elem, dx, m, z_k)
          rhs(1, :) = rhs(1, :) + weight * z_k
          do a = 1, n_active
            rhs(1+a, :) = rhs(1+a, :) + (weight * dx(ls_active_dim_cache(a, v))) * z_k
          end do
        end do

        sol = matmul(ls_mat_inv_cache(1:n_basis, 1:n_basis, v), rhs)

        do tt = 1, n_tt
          do a = 1, n_active
            ls_grad_h1_cache(m)%m((tt-1)*3+ls_active_dim_cache(a, v), v) = sol(1+a, tt)
          end do
        end do

        deallocate(rhs, sol)
      end do
      deallocate(z_k)
    end do

    ls_grad_h1_cache_n_vert = mesh%n_vert
    ls_grad_h1_cache_max_order = max_order
    ls_grad_h1_cache_boundary_2d = boundary_2d
  end subroutine ensure_ls_grad_h1_cache

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

    ! Geometry is precomputed in ensure_green_gauss_mat_cache (gg_mat_inv_cache is the plain
    ! pseudo-inverse for an interior vertex, or the wall-tangent-projected effective operator for a
    ! wall vertex -- see that cache's own header comment); only this phi gather is redone every call.
    ! i1 outer / a inner (a=1:3 is grad_raw's and gg_flux_w_cache's own contiguous leading
    ! dimension): each inner step is a length-3 axpy on contiguous memory instead of a
    ! stride-3 write, same total sum either way.
    grad_raw = 0.0_DOUBLE
    k0 = gg_flux_offset_cache(id_vert)
    k1 = gg_flux_offset_cache(id_vert + 1) - 1
    do k = k0, k1
      id_elem = gg_flux_elem_cache(k)
      do i1 = 1, nc_in
        grad_raw(:, i1) = grad_raw(:, i1) + phi(i1, id_elem) * gg_flux_w_cache(:, k)
      end do
    end do

    grad_true = matmul(gg_mat_inv_cache(:, :, id_vert), grad_raw)

    ! Gradient-norm oscillation indicator (no residual counterpart for GG); h_local is cached, phi_scale2 depends on phi and is gathered here.
    n_cand = merge(2_ENTIER, 3_ENTIER, boundary_2d)
    weight_sum = real(size(mesh%vert(id_vert)%elem_neigh), kind=DOUBLE)
    phi_sq_sum = 0.0_DOUBLE
    do i1 = 1, size(mesh%vert(id_vert)%elem_neigh)
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

    ! Single merged pass over vertices: the original code ran two separate do iv=1,mesh%n_vert
    ! loops back to back -- one for the bias correction (needs hess AND third valid), one for the
    ! cell-blend curvature correction (needs only third valid) -- each independently re-extracting
    ! the SAME 10-component Txxx..Txyz reduction from third_v and recomputing the SAME vweight.
    ! Merged into one pass (extract once, guard the hess-dependent part on valid_hess_v(iv) exactly
    ! as the original outer cycle did, then feed both accumulations from a single sub_elem_neigh
    ! loop since sub_elem_volume is identical in both) -- purely eliminates duplicate work, no
    ! change to any formula or to which vertices contribute to which correction.
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

      ! Match grad_cell's own vertex-to-cell weighting exactly.
      if (use_weno_blend) then
        vweight = 1.0_DOUBLE / (eps_weight_num + grad_oi_v(iv)**weno_power)
      else
        vweight = 1.0_DOUBLE
      end if

      if (valid_hess_v(iv)) then
        do ic = 1, nc_in
          Hxx(ic) = hess_v(0*nc_in+ic, iv)
          Hyy(ic) = hess_v(4*nc_in+ic, iv)
          Hzz(ic) = hess_v(8*nc_in+ic, iv)
          Hxy(ic) = 0.5_DOUBLE*(hess_v(1*nc_in+ic, iv) + hess_v(3*nc_in+ic, iv))
          Hxz(ic) = 0.5_DOUBLE*(hess_v(2*nc_in+ic, iv) + hess_v(6*nc_in+ic, iv))
          Hyz(ic) = 0.5_DOUBLE*(hess_v(5*nc_in+ic, iv) + hess_v(7*nc_in+ic, iv))
        end do

        ! Index neigh_cache_list directly instead of copying this vertex's neighbor slice into a
        ! freshly allocate()'d local array every single vertex (n_vert allocate/deallocate cycles
        ! per call, a real cost on a tet mesh where a vertex's neigh_by_vert stencil is large).
        n_neigh = neigh_cache_start(iv+1) - neigh_cache_start(iv)
        mat(1:n_basis, 1:n_basis) = 0.0_DOUBLE
        rhs(1:n_basis, :) = 0.0_DOUBLE
        do j = 1, n_neigh
          id_elem = neigh_cache_list(neigh_cache_start(iv)+j-1)
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

        call lu_factor_lapack(n_basis, mat(1:n_basis, 1:n_basis), ipiv(1:n_basis))
        call lu_solve_mat_lapack(n_basis, mat(1:n_basis, 1:n_basis), ipiv(1:n_basis), nc_in, rhs(1:n_basis, :))
      end if

      do j = 1, mesh%vert(iv)%n_sub_elems_neigh
        id_sub_elem = mesh%vert(iv)%sub_elem_neigh(j)
        id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
        sub_elem_volume = mesh%sub_elem(id_sub_elem)%volume * vweight

        if (valid_hess_v(iv)) then
          bias_num_x(:, id_elem) = bias_num_x(:, id_elem) + sub_elem_volume*rhs(2, :)
          bias_num_y(:, id_elem) = bias_num_y(:, id_elem) + sub_elem_volume*rhs(3, :)
          if (n_basis == 4) bias_num_z(:, id_elem) = bias_num_z(:, id_elem) + sub_elem_volume*rhs(4, :)
          bias_den(id_elem) = bias_den(id_elem) + sub_elem_volume
        end if

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
      if (bias_den(i) <= 0.0_DOUBLE) cycle
      grad_cell(1:nc_in, i) = grad_cell(1:nc_in, i) - bias_num_x(:, i)/bias_den(i)
      grad_cell(nc_in+1:2*nc_in, i) = grad_cell(nc_in+1:2*nc_in, i) - bias_num_y(:, i)/bias_den(i)
      if (.not. boundary_2d) grad_cell(2*nc_in+1:3*nc_in, i) = grad_cell(2*nc_in+1:3*nc_in, i) - bias_num_z(:, i)/bias_den(i)
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
  ! Loop nest is vertex-outermost so each inner pass over t/i/ic stays within one contiguous
  ! column of D_m_v/H_m1_v/bias_v (arrays are (component, vertex), vertex trailing/column-major) --
  ! the previous component-outermost order made every access a full-n_vert stride, one useful
  ! element per cache line, which dominated compute_node_derivative_bias's cost (~0.23s/call at
  ! N=160, vs a ~10ms flop-count estimate) and was the main driver of order 4's slowdown.
  subroutine contract_grad_node_bias(nc_in, m, n_vert, D_m_v, H_m1_v, bias_v)
    implicit none

    integer(kind=ENTIER), intent(in) :: nc_in, m, n_vert
    real(kind=DOUBLE), dimension(3**m*nc_in, n_vert), intent(in) :: D_m_v
    real(kind=DOUBLE), dimension(3**m*3, n_vert), intent(in) :: H_m1_v
    real(kind=DOUBLE), dimension(3*nc_in, n_vert), intent(out) :: bias_v

    integer(kind=ENTIER) :: t, i, ic, n_t, v

    n_t = 3**m
    do v = 1, n_vert
      bias_v(:, v) = 0.0_DOUBLE
      do t = 0, n_t - 1
        do i = 1, 3
          do ic = 1, nc_in
            bias_v((i-1)*nc_in+ic, v) = bias_v((i-1)*nc_in+ic, v) &
              + D_m_v(t*nc_in+ic, v) * H_m1_v(t*3+i, v)
          end do
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

    ! Dispatches on the SAME module-level flag aho_reconstruction uses to pick aho_gg vs aho_ls:
    ! the h1-cache must be built from whichever operator actually produced dfield_v's own m-exact
    ! derivatives, or eq. 13's correction does not cancel the right bias.
    if (use_green_gauss) then
      call ensure_gg_gradient_h1_cache(mesh, boundary_2d, k_max-k+1)
    else
      call ensure_ls_grad_h1_cache(mesh, boundary_2d, k_max-k+1)
    end if

    nc_eff = nc_in * d**(k-1)
    do m = 2, k_max-k+1
      allocate(bias_m(nc_eff*d, mesh%n_vert))
      if (use_green_gauss) then
        call contract_grad_node_bias(nc_eff, m, mesh%n_vert, dfield_v(k-1+m)%val, gg_grad_h1_cache(m)%m, bias_m)
        if (debug_bias_vertex > 0) print '(A,I2,A,I2,A,ES14.4,A,ES14.4)', '    [dbg k=', k, ' m=', m, &
          '] |H_m1|=', sqrt(sum(gg_grad_h1_cache(m)%m(:,debug_bias_vertex)**2)), &
          ' |Dm|=', sqrt(sum(dfield_v(k-1+m)%val(:,debug_bias_vertex)**2))
      else
        call contract_grad_node_bias(nc_eff, m, mesh%n_vert, dfield_v(k-1+m)%val, ls_grad_h1_cache(m)%m, bias_m)
        if (debug_bias_vertex > 0) print '(A,I2,A,I2,A,ES14.4,A,ES14.4)', '    [dbg k=', k, ' m=', m, &
          '] |H_m1|=', sqrt(sum(ls_grad_h1_cache(m)%m(:,debug_bias_vertex)**2)), &
          ' |Dm|=', sqrt(sum(dfield_v(k-1+m)%val(:,debug_bias_vertex)**2))
      end if
      if (debug_bias_vertex > 0) print '(A,ES14.4)', '    [dbg] |bias_m|=', sqrt(sum(bias_m(:,debug_bias_vertex)**2))
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

    call apply_derivative_node_correction(mesh, d, nc_in, boundary_2d, 1_ENTIER, k_max, &
      grad_cell, dfield_v, valid_v)
  end subroutine apply_gradient_node_correction

  ! Generalizes apply_gradient_node_correction to any order k (k=1 is the gradient case): corrects
  ! dfield_cell (a k-exact k-th derivative in cell-average form) in place using every available
  ! higher vertex derivative dfield_v(k+1)%val,...,dfield_v(k_max)%val, per Haider, Croisille &
  ! Courbet 2011 eq. 13 (compute_node_derivative_bias, which itself dispatches on use_green_gauss
  ! to reproduce whichever operator -- aho_gg or aho_ls -- actually built dfield_v), then blends
  ! the correction into cells via a plain sub_elem_volume average (unweighted, matching the
  ! original grad-only routine -- the correction is a geometry+derivative bias, not a re-blend of
  ! samples, so it does not need grad_cell's own WENO weight).
  subroutine apply_derivative_node_correction(mesh, d, nc_in, boundary_2d, k, k_max, dfield_cell, &
      dfield_v, valid_v)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in, k, k_max
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(nc_in*d**k, mesh%n_elems), intent(inout) :: dfield_cell
    type(derivative_field_type), dimension(k+1:k_max), intent(in) :: dfield_v
    logical, dimension(mesh%n_vert), intent(in) :: valid_v

    integer(kind=ENTIER) :: iv, j, id_elem, id_sub_elem
    real(kind=DOUBLE), dimension(:, :), allocatable :: bias_v_total
    real(kind=DOUBLE), dimension(:, :), allocatable :: bias_num
    real(kind=DOUBLE), dimension(:), allocatable :: bias_den
    real(kind=DOUBLE) :: sub_elem_volume

    if (k_max < k+1) return

    allocate(bias_v_total(nc_in*d**k, mesh%n_vert))
    call compute_node_derivative_bias(mesh, d, nc_in, boundary_2d, k, k_max, dfield_v, bias_v_total)

    allocate(bias_num(nc_in*d**k, mesh%n_elems), bias_den(mesh%n_elems))
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
      dfield_cell(:, j) = dfield_cell(:, j) - bias_num(:, j) / bias_den(j)
    end do

    deallocate(bias_v_total, bias_num, bias_den)
  end subroutine apply_derivative_node_correction

  ! Same cell-blend as apply_derivative_node_correction's second half, but takes an already-computed
  ! per-vertex bias instead of calling compute_node_derivative_bias itself -- for a caller (the
  ! eq.13 cascade) that also needs the per-vertex bias on its own (to correct a lower derivative's
  ! own vertex values before using them to correct the next one down), so the bias is computed once
  ! and reused rather than recomputed from scratch for the cell-level pass.
  subroutine scatter_bias_to_cells(mesh, boundary_2d, n_comp, bias_v_total, valid_v, dfield_cell)
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d
    integer(kind=ENTIER), intent(in) :: n_comp
    real(kind=DOUBLE), dimension(n_comp, mesh%n_vert), intent(in) :: bias_v_total
    logical, dimension(mesh%n_vert), intent(in) :: valid_v
    real(kind=DOUBLE), dimension(n_comp, mesh%n_elems), intent(inout) :: dfield_cell

    integer(kind=ENTIER) :: iv, j, id_elem, id_sub_elem
    real(kind=DOUBLE), dimension(:, :), allocatable :: bias_num
    real(kind=DOUBLE), dimension(:), allocatable :: bias_den
    real(kind=DOUBLE) :: sub_elem_volume

    allocate(bias_num(n_comp, mesh%n_elems), bias_den(mesh%n_elems))
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
      dfield_cell(:, j) = dfield_cell(:, j) - bias_num(:, j) / bias_den(j)
    end do

    deallocate(bias_num, bias_den)
  end subroutine scatter_bias_to_cells

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

  ! aho_cls, fully standalone (no dependence on aho_gg/aho_ls's own fit),
  ! recursive to arbitrary order. Builds p(1),p(2),...,p(k_max) -- each
  ! GENUINELY exact to its own order -- alternating a node step and a
  ! cell-blend step:
  !   p(0)_cell = phi (given)
  !     -> [order-1 node step, below]      -> p(1)_node   (grad_v)
  !     -> [linear/WENO blend]             -> p(1)_cell   (grad_out)
  !     -> [functional-identity node step] -> p(2)_node   (hess_v)
  !     -> [linear/WENO blend]             -> p(2)_cell   (hess_out)
  !     -> [functional-identity node step] -> p(3)_node   (third_v)
  !     -> [linear/WENO blend]             -> p(3)_cell   (third_out)
  ! compute_order1_node_aho_cls (k=1) and compute_next_order_node_aho_cls
  ! (k>=2, SAME routine for k=2,3,...,40) both rest on the identical idea --
  ! a quantity built by a SIMPLE, EXPLICIT interpolation/blend of known data
  ! has an EXACTLY COMPUTABLE bias against the true field, a plain discrete
  ! geometric moment (no quadrature, no gg_grad_h1_cache) -- but they differ
  ! because the interpolation runs in OPPOSITE directions at k=1 (node built
  ! FROM cells, order 0) vs k>=2 (cell built FROM nodes, order k-1, from the
  ! previous blend step), so the discrete moment (nu_v vs mu_C(v)) and the
  ! delta each corrects are not the same formula -- seeded from an ordinary
  ! GG/LS fit instead, grad_v carries an UNCOMPUTABLE-here bias (it needs
  ! gg_grad_h1_cache to characterize), which is why that was tried and
  ! abandoned in favor of this fully self-contained construction.
  subroutine compute_derivatives_aho_cls(mesh, nc_in, phi, grad_out, hess_out, third_out)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: nc_in
    real(kind=DOUBLE), dimension(nc_in, mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(nc_in*3, mesh%n_elems), intent(out) :: grad_out
    real(kind=DOUBLE), dimension(nc_in*9, mesh%n_elems), intent(out) :: hess_out
    real(kind=DOUBLE), dimension(nc_in*27, mesh%n_elems), intent(out) :: third_out

    real(kind=DOUBLE), dimension(:, :), allocatable :: grad_v, hess_v, third_v
    real(kind=DOUBLE), dimension(:), allocatable :: grad_oi_v, hess_oi_v, third_oi_v_unused

    allocate(grad_v(nc_in*3, mesh%n_vert), grad_oi_v(mesh%n_vert))
    call compute_order1_node_aho_cls(mesh, nc_in, phi, grad_v, grad_oi_v)
    call blend_node_to_cell_aho_cls(mesh, nc_in*3_ENTIER, grad_v, grad_oi_v, grad_out)

    allocate(hess_v(nc_in*9, mesh%n_vert), hess_oi_v(mesh%n_vert))
    call compute_next_order_node_aho_cls(mesh, nc_in, 2_ENTIER, grad_v, grad_out, grad_oi_v, hess_v, hess_oi_v)
    call blend_node_to_cell_aho_cls(mesh, nc_in*9_ENTIER, hess_v, hess_oi_v, hess_out)

    allocate(third_v(nc_in*27, mesh%n_vert), third_oi_v_unused(mesh%n_vert))
    call compute_next_order_node_aho_cls(mesh, nc_in, 3_ENTIER, hess_v, hess_out, hess_oi_v, &
      third_v, third_oi_v_unused)
    call blend_node_to_cell_aho_cls(mesh, nc_in*27_ENTIER, third_v, third_oi_v_unused, third_out)

    deallocate(grad_v, hess_v, third_v, grad_oi_v, hess_oi_v, third_oi_v_unused)
  end subroutine compute_derivatives_aho_cls

  ! Special k=1 node step (no p(0)_node exists to difference against, unlike
  ! k>=2): phi_v(v), a SIMPLE inverse-square-distance-weighted interpolation
  ! of phi over v's own neighbor cells, has an exactly computable bias
  ! against phi(x_v) for a genuinely affine field: phi_v(v) = phi(x_v) +
  ! grad(x_v).nu_v, nu_v the SAME-weighted average of (x_C-x_v). So for each
  ! neighbor cell C, phi(C) - phi_v(v) = grad(x_v).[(x_C-x_v) - nu_v] exactly
  ! -- solved by weighted LS over v's 1-ring for grad_v, without ever
  ! forming phi_v as a separate output (it cancels out of the fit). Verified
  ! to reproduce grad_v exactly for a synthetic affine field even though
  ! phi_v itself is measurably biased.
  subroutine compute_order1_node_aho_cls(mesh, nc_in, phi, grad_v, oi_v)
    use linear_solver_module, only: lu_factor_lapack, lu_solve_mat_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: nc_in
    real(kind=DOUBLE), dimension(nc_in, mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(nc_in*3, mesh%n_vert), intent(out) :: grad_v
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(out) :: oi_v

    integer(kind=ENTIER) :: id_vert, n_neigh, gi, id_elem, r
    real(kind=DOUBLE), dimension(:, :), allocatable :: Mmat, rhs, rhs_normal
    real(kind=DOUBLE), dimension(:), allocatable :: beta
    real(kind=DOUBLE), dimension(3, 3) :: normal_mat
    integer(kind=ENTIER), dimension(3) :: ipiv
    real(kind=DOUBLE), dimension(3) :: nu, dxloc
    real(kind=DOUBLE) :: bsum, resid_e, scale2
    real(kind=DOUBLE), dimension(nc_in) :: phi_v_vec

    grad_v = 0.0_DOUBLE
    oi_v = 0.0_DOUBLE

    do id_vert = 1, mesh%n_vert
      if (mesh%vert(id_vert)%is_bound) cycle
      n_neigh = mesh%vert(id_vert)%n_elems_neigh
      if (n_neigh < 3) cycle

      allocate(Mmat(n_neigh, 3), rhs(n_neigh, nc_in), beta(n_neigh))

      nu = 0.0_DOUBLE
      bsum = 0.0_DOUBLE
      do gi = 1, n_neigh
        id_elem = mesh%vert(id_vert)%elem_neigh(gi)
        dxloc = mesh%elem(id_elem)%coord - mesh%vert(id_vert)%coord
        beta(gi) = 1.0_DOUBLE / max(dot_product(dxloc, dxloc), 1.0e-24_DOUBLE)
        nu = nu + beta(gi) * dxloc
        bsum = bsum + beta(gi)
      end do
      if (bsum > 0.0_DOUBLE) nu = nu / bsum

      phi_v_vec = 0.0_DOUBLE
      do gi = 1, n_neigh
        id_elem = mesh%vert(id_vert)%elem_neigh(gi)
        phi_v_vec = phi_v_vec + beta(gi) * phi(:, id_elem)
      end do
      if (bsum > 0.0_DOUBLE) phi_v_vec = phi_v_vec / bsum

      do gi = 1, n_neigh
        id_elem = mesh%vert(id_vert)%elem_neigh(gi)
        dxloc = mesh%elem(id_elem)%coord - mesh%vert(id_vert)%coord
        Mmat(gi, :) = dxloc - nu
        rhs(gi, :) = phi(:, id_elem) - phi_v_vec
      end do

      normal_mat = 0.0_DOUBLE
      do gi = 1, n_neigh
        normal_mat = normal_mat + beta(gi) * spread(Mmat(gi,:),2,3) * spread(Mmat(gi,:),1,3)
      end do

      allocate(rhs_normal(3, nc_in))
      rhs_normal = 0.0_DOUBLE
      do gi = 1, n_neigh
        do r = 1, 3
          rhs_normal(r, :) = rhs_normal(r, :) + beta(gi) * Mmat(gi, r) * rhs(gi, :)
        end do
      end do

      call lu_factor_lapack(3_ENTIER, normal_mat, ipiv)
      call lu_solve_mat_lapack(3_ENTIER, normal_mat, ipiv, nc_in, rhs_normal)

      do r = 1, 3
        grad_v((r-1)*nc_in+1 : r*nc_in, id_vert) = rhs_normal(r, :)
      end do

      resid_e = 0.0_DOUBLE
      scale2 = 0.0_DOUBLE
      do gi = 1, n_neigh
        block
          real(kind=DOUBLE), dimension(nc_in) :: pred
          pred = 0.0_DOUBLE
          do r = 1, 3
            pred = pred + Mmat(gi, r) * grad_v((r-1)*nc_in+1 : r*nc_in, id_vert)
          end do
          resid_e = resid_e + sum((pred - rhs(gi, :))**2)
        end block
        id_elem = mesh%vert(id_vert)%elem_neigh(gi)
        scale2 = max(scale2, maxval(phi(:, id_elem)**2))
      end do
      oi_v(id_vert) = (resid_e / real(n_neigh*nc_in, kind=DOUBLE)) / max(scale2, 1.0e-300_DOUBLE)

      deallocate(Mmat, rhs, beta, rhs_normal)
    end do
  end subroutine compute_order1_node_aho_cls

  ! Blends a per-vertex tensor (ncomp components, any order -- k=1..40 alike)
  ! into cells via the SAME sub_elem_volume/WENO weight used everywhere else
  ! in this module: weight = sub_elem_volume (use_weno_blend=.false.) or
  ! sub_elem_volume/(eps_weight_num+oi_node) (use_weno_blend=.true.).
  subroutine blend_node_to_cell_aho_cls(mesh, ncomp, p_node, oi_node, p_cell)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: ncomp
    real(kind=DOUBLE), dimension(ncomp, mesh%n_vert), intent(in) :: p_node
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: oi_node
    real(kind=DOUBLE), dimension(ncomp, mesh%n_elems), intent(out) :: p_cell

    real(kind=DOUBLE), dimension(:, :), allocatable :: num
    real(kind=DOUBLE), dimension(:), allocatable :: den
    integer(kind=ENTIER) :: id_vert, j, id_elem, id_sub_elem
    real(kind=DOUBLE) :: sub_elem_volume, w

    allocate(num(ncomp, mesh%n_elems), den(mesh%n_elems))
    num = 0.0_DOUBLE; den = 0.0_DOUBLE

    do id_vert = 1, mesh%n_vert
      if (mesh%vert(id_vert)%is_bound) cycle
      do j = 1, mesh%vert(id_vert)%n_sub_elems_neigh
        id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
        id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
        sub_elem_volume = mesh%sub_elem(id_sub_elem)%volume
        if (use_weno_blend) then
          w = sub_elem_volume / (eps_weight_num + oi_node(id_vert))
        else
          w = sub_elem_volume
        end if
        num(:, id_elem) = num(:, id_elem) + w * p_node(:, id_vert)
        den(id_elem) = den(id_elem) + w
      end do
    end do

    do id_elem = 1, mesh%n_elems
      if (den(id_elem) > 0.0_DOUBLE) then
        p_cell(:, id_elem) = num(:, id_elem) / den(id_elem)
      else
        p_cell(:, id_elem) = 0.0_DOUBLE
      end if
    end do

    deallocate(num, den)
  end subroutine blend_node_to_cell_aho_cls

  ! GENERAL k-exact node step, identical for every k>=2 (k=2,3,...,40): builds
  ! p_next_node = D^k(x_v), genuinely k-exact, at every non-boundary vertex,
  ! from p_prev_node = D^(k-1)(x_v) (already known, built at the previous
  ! step) and p_prev_cell = D^(k-1) blended into cells (blend_node_to_cell_
  ! aho_cls's own output, using oi_prev_node's weights). For a genuine
  ! degree-k field, D^(k-1) is exactly affine, so each neighbor cell C's own
  ! blend gap delta_C := p_prev_cell(C) - p_prev_node(v) equals exactly
  ! D^k(x_v) contracted with mu_C(v), the SAME sub_elem_volume/WENO-weighted
  ! average of (x_v'-x_v) over C's own corner vertices v' used to build that
  ! blend (a DISCRETE moment, no continuous quadrature cache needed here --
  ! this is Haider's functional identity, eq. 15-18, generalized from
  ! vertex-vertex differences (compute_2exact_hessian_aho_cls) to
  ! vertex-cell differences, so it composes without a growing stencil or a
  ! growing per-step linear system). Stacking over v's own 1-ring of cells
  ! (n_neigh equations) and contracting only the LAST tensor index of D^k
  ! against mu_C(v) decouples into nc_in*3**(k-1) INDEPENDENT 3-unknown
  ! weighted least-squares fits sharing one (n_neigh x 3) geometry matrix --
  ! exactly the size of an ordinary gradient fit, at any k. oi_next_node is
  ! this fit's own normalized residual, used as the NEXT blend's WENO weight.
  subroutine compute_next_order_node_aho_cls(mesh, nc_in, k, p_prev_node, p_prev_cell, oi_prev_node, &
      p_next_node, oi_next_node)
    use linear_solver_module, only: lu_factor_lapack, lu_solve_mat_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: nc_in, k
    real(kind=DOUBLE), dimension(nc_in*3_ENTIER**(k-1), mesh%n_vert), intent(in) :: p_prev_node
    real(kind=DOUBLE), dimension(nc_in*3_ENTIER**(k-1), mesh%n_elems), intent(in) :: p_prev_cell
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: oi_prev_node
    real(kind=DOUBLE), dimension(nc_in*3_ENTIER**k, mesh%n_vert), intent(out) :: p_next_node
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(out) :: oi_next_node

    integer(kind=ENTIER) :: ncomp, id_vert, n_neigh, gi, id_elem, jc, id_v2, n_v2, r
    real(kind=DOUBLE), dimension(:, :), allocatable :: Mmat, rhs, rhs_normal
    real(kind=DOUBLE), dimension(:), allocatable :: wvec
    real(kind=DOUBLE), dimension(3, 3) :: normal_mat
    integer(kind=ENTIER), dimension(3) :: ipiv
    real(kind=DOUBLE), dimension(3) :: mu, dxloc
    real(kind=DOUBLE) :: wsum, wloc, sub_elem_volume, resid_e, scale2

    ncomp = nc_in * 3_ENTIER**(k-1)
    p_next_node = 0.0_DOUBLE
    oi_next_node = 0.0_DOUBLE

    do id_vert = 1, mesh%n_vert
      if (mesh%vert(id_vert)%is_bound) cycle
      n_neigh = mesh%vert(id_vert)%n_elems_neigh
      if (n_neigh < 3) cycle

      allocate(Mmat(n_neigh, 3), rhs(n_neigh, ncomp), wvec(n_neigh))

      do gi = 1, n_neigh
        id_elem = mesh%vert(id_vert)%elem_neigh(gi)
        mu = 0.0_DOUBLE
        wsum = 0.0_DOUBLE
        n_v2 = mesh%elem(id_elem)%n_vert
        do jc = 1, n_v2
          id_v2 = mesh%elem(id_elem)%vert(jc)
          sub_elem_volume = mesh%sub_elem(mesh%elem(id_elem)%sub_elem(jc))%volume
          if (use_weno_blend) then
            wloc = sub_elem_volume / (eps_weight_num + oi_prev_node(id_v2))
          else
            wloc = sub_elem_volume
          end if
          dxloc = mesh%vert(id_v2)%coord - mesh%vert(id_vert)%coord
          mu = mu + wloc * dxloc
          wsum = wsum + wloc
        end do
        if (wsum > 0.0_DOUBLE) mu = mu / wsum

        Mmat(gi, :) = mu
        wvec(gi) = 1.0_DOUBLE / max(dot_product(mu, mu), 1.0e-24_DOUBLE)
        rhs(gi, :) = p_prev_cell(:, id_elem) - p_prev_node(:, id_vert)
      end do

      normal_mat = 0.0_DOUBLE
      do gi = 1, n_neigh
        normal_mat = normal_mat + wvec(gi) * spread(Mmat(gi,:),2,3) * spread(Mmat(gi,:),1,3)
      end do

      allocate(rhs_normal(3, ncomp))
      rhs_normal = 0.0_DOUBLE
      do gi = 1, n_neigh
        do r = 1, 3
          rhs_normal(r, :) = rhs_normal(r, :) + wvec(gi) * Mmat(gi, r) * rhs(gi, :)
        end do
      end do

      call lu_factor_lapack(3_ENTIER, normal_mat, ipiv)
      call lu_solve_mat_lapack(3_ENTIER, normal_mat, ipiv, ncomp, rhs_normal)

      ! rhs_normal(r,:) = D^k(x_v)[..., r] -- r is the newly-appended, most
      ! significant tensor index (same convention as hess_flat/third_flat).
      do r = 1, 3
        p_next_node((r-1)*ncomp+1 : r*ncomp, id_vert) = rhs_normal(r, :)
      end do

      resid_e = 0.0_DOUBLE
      scale2 = 0.0_DOUBLE
      do gi = 1, n_neigh
        block
          real(kind=DOUBLE), dimension(ncomp) :: pred
          pred = 0.0_DOUBLE
          do r = 1, 3
            pred = pred + Mmat(gi, r) * p_next_node((r-1)*ncomp+1 : r*ncomp, id_vert)
          end do
          resid_e = resid_e + sum((pred - rhs(gi, :))**2)
        end block
        scale2 = max(scale2, maxval(p_prev_cell(:, mesh%vert(id_vert)%elem_neigh(gi))**2))
      end do
      oi_next_node(id_vert) = (resid_e / real(n_neigh*ncomp, kind=DOUBLE)) / max(scale2, 1.0e-300_DOUBLE)

      deallocate(Mmat, rhs, wvec, rhs_normal)
    end do
  end subroutine compute_next_order_node_aho_cls

  ! Builds cell_ls_mat_inv_cache (Haider, Croisille & Courbet 2011, Definition 1, step 1: a
  ! 1-exact 1st derivative computed directly from cell averages on a small stencil). Unweighted LS
  ! normal-equations inverse over mesh%elem(alpha)%neigh_by_vert (the cell's own vertex-stencil
  ! neighbor cells) -- purely cell-centered, no vertex fit or blend involved anywhere.
  subroutine ensure_cell_ls_grad_mat_cache(mesh, boundary_2d)
    use linear_solver_module, only: pseudo_inverse_inplace_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d

    integer(kind=ENTIER) :: alpha, i, beta, a, b
    real(kind=DOUBLE), dimension(3) :: h
    real(kind=DOUBLE), dimension(3, 3) :: mat, mat_inv

    if (cell_ls_mat_cache_n_elems == mesh%n_elems .and. (cell_ls_mat_cache_boundary_2d .eqv. boundary_2d)) return

    if (allocated(cell_ls_mat_inv_cache)) deallocate(cell_ls_mat_inv_cache)
    allocate(cell_ls_mat_inv_cache(3, 3, mesh%n_elems))

    do alpha = 1, mesh%n_elems
      mat = 0.0_DOUBLE
      do i = 1, mesh%elem(alpha)%n_neigh_by_vert
        beta = mesh%elem(alpha)%neigh_by_vert(i)
        h = mesh%elem(beta)%coord - mesh%elem(alpha)%coord
        do a = 1, 3
          do b = 1, 3
            mat(a, b) = mat(a, b) + h(a)*h(b)
          end do
        end do
      end do
      mat_inv = mat
      call pseudo_inverse_inplace_lapack(3_ENTIER, mat_inv)
      if (boundary_2d) mat_inv(3, :) = 0.0_DOUBLE
      cell_ls_mat_inv_cache(:, :, alpha) = mat_inv
    end do

    cell_ls_mat_cache_n_elems = mesh%n_elems
    cell_ls_mat_cache_boundary_2d = boundary_2d
  end subroutine ensure_cell_ls_grad_mat_cache

  ! Cell-to-cell 1-exact gradient (Definition 1, step 1): grad_cell(:,alpha) =
  ! cell_ls_mat_inv_cache(alpha) @ sum_beta h_{alpha,beta}*(phi(beta)-phi(alpha)), beta ranging
  ! over neigh_by_vert(alpha). Caller must have called ensure_cell_ls_grad_mat_cache first.
  subroutine compute_cell_ls_grad(mesh, nc_in, phi, grad_cell)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: nc_in
    real(kind=DOUBLE), dimension(nc_in, mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(nc_in*3, mesh%n_elems), intent(out) :: grad_cell

    integer(kind=ENTIER) :: alpha, i, beta, ic, p
    real(kind=DOUBLE), dimension(3) :: h
    real(kind=DOUBLE), dimension(3, nc_in) :: rhs, res

    do alpha = 1, mesh%n_elems
      rhs = 0.0_DOUBLE
      do i = 1, mesh%elem(alpha)%n_neigh_by_vert
        beta = mesh%elem(alpha)%neigh_by_vert(i)
        h = mesh%elem(beta)%coord - mesh%elem(alpha)%coord
        do ic = 1, nc_in
          rhs(:, ic) = rhs(:, ic) + h * (phi(ic, beta) - phi(ic, alpha))
        end do
      end do
      res = matmul(cell_ls_mat_inv_cache(:, :, alpha), rhs)
      do ic = 1, nc_in
        do p = 1, 3
          grad_cell((p-1)*nc_in+ic, alpha) = res(p, ic)
        end do
      end do
    end do
  end subroutine compute_cell_ls_grad

  ! Cell-to-cell analogue of ensure_gg_gradient_h1_cache: H_m^(1)(alpha), m=2..max_order, built by
  ! applying the SAME cell_ls_mat_inv_cache grad operator to the geometric field
  ! {z_{alpha,beta}^(m)-z_{alpha,alpha}^(m)}_beta (shifted_cell_moment_full) instead of phi, over
  ! alpha's own neigh_by_vert stencil. This is Haider's w_beta^(k|k) applied to z_beta^(k+1) (eq.
  ! 15-16) for the k=1 case. Pure geometry, cached once per mesh.
  subroutine ensure_cell_ls_h1_cache(mesh, boundary_2d, max_order)
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d
    integer(kind=ENTIER), intent(in) :: max_order

    integer(kind=ENTIER) :: alpha, i, beta, m, tt, n_tt
    real(kind=DOUBLE), dimension(3) :: h, tmp
    real(kind=DOUBLE), dimension(:, :), allocatable :: rhs
    real(kind=DOUBLE), dimension(:), allocatable :: z_beta, z_alpha

    call ensure_cell_ls_grad_mat_cache(mesh, boundary_2d)
    call ensure_cell_moment_cache(mesh, max_order)

    if (cell_ls_h1_cache_n_elems == mesh%n_elems .and. cell_ls_h1_cache_max_order >= max_order &
        .and. (cell_ls_h1_cache_boundary_2d .eqv. boundary_2d)) return

    if (allocated(cell_ls_h1_cache)) deallocate(cell_ls_h1_cache)
    allocate(cell_ls_h1_cache(2:max_order))
    do m = 2, max_order
      allocate(cell_ls_h1_cache(m)%m(3*3**m, mesh%n_elems))
      cell_ls_h1_cache(m)%m = 0.0_DOUBLE
    end do

    do m = 2, max_order
      n_tt = 3**m
      allocate(z_beta(n_tt), z_alpha(n_tt), rhs(3, n_tt))
      do alpha = 1, mesh%n_elems
        z_alpha = cell_moment_cache(m)%m(:, alpha)
        rhs = 0.0_DOUBLE
        do i = 1, mesh%elem(alpha)%n_neigh_by_vert
          beta = mesh%elem(alpha)%neigh_by_vert(i)
          h = mesh%elem(beta)%coord - mesh%elem(alpha)%coord
          call shifted_cell_moment_full(beta, h, m, z_beta)
          do tt = 1, n_tt
            rhs(:, tt) = rhs(:, tt) + h * (z_beta(tt) - z_alpha(tt))
          end do
        end do
        do tt = 1, n_tt
          tmp = matmul(cell_ls_mat_inv_cache(:, :, alpha), rhs(:, tt))
          cell_ls_h1_cache(m)%m((tt-1)*3+1:(tt-1)*3+3, alpha) = tmp
        end do
      end do
      deallocate(z_beta, z_alpha, rhs)
    end do

    cell_ls_h1_cache_n_elems = mesh%n_elems
    cell_ls_h1_cache_max_order = max_order
    cell_ls_h1_cache_boundary_2d = boundary_2d
  end subroutine ensure_cell_ls_h1_cache

  ! Cell-to-cell CLS Hessian (Haider, Croisille & Courbet 2011, Definition 1, step 2, k=1->2):
  ! builds a genuinely 2-exact Hessian directly at cell centers from a 1-exact gradient
  ! (grad_cell, compute_cell_ls_grad) and the SAME operator's own bias on the geometric moment
  ! field (cell_ls_h1_cache(2)), following the functional identity (eq. 15-18) with i indexing
  ! neigh_by_vert(alpha) instead of vv_neigh_cache. Direct cell-to-cell port of the already-proven
  ! vertex-vertex compute_2exact_hessian_aho_cls -- this is Haider's OWN algorithm (which never
  ! involves vertices/nodes at all), as opposed to Setzwein's vertex-centered adaptation of it.
  subroutine compute_cell_cls_2exact_hessian(mesh, boundary_2d, grad_cell, hess_cell)
    use linear_solver_module, only: lu_factor_lapack, lu_solve_mat_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(3, mesh%n_elems), intent(in) :: grad_cell
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(out) :: hess_cell ! Hxx,Hxy,Hxz,Hyy,Hyz,Hzz

    integer(kind=ENTIER) :: alpha, i, beta, ell, row, n_neigh
    real(kind=DOUBLE), dimension(:, :), allocatable :: Jmat
    real(kind=DOUBLE), dimension(:), allocatable :: rhsvec
    real(kind=DOUBLE), dimension(3) :: h, gcorr_alpha, gcorr_beta
    real(kind=DOUBLE), dimension(9) :: E9
    real(kind=DOUBLE), dimension(6, 6) :: normal_mat
    real(kind=DOUBLE), dimension(6, 1) :: normal_rhs
    integer(kind=ENTIER), dimension(6) :: ipiv
    integer(kind=ENTIER), parameter :: basis_i(6) = (/1,1,1,2,2,3/)
    integer(kind=ENTIER), parameter :: basis_j(6) = (/1,2,3,2,3,3/)

    call ensure_cell_ls_h1_cache(mesh, boundary_2d, 2_ENTIER)

    hess_cell = 0.0_DOUBLE

    do alpha = 1, mesh%n_elems
      n_neigh = mesh%elem(alpha)%n_neigh_by_vert
      if (n_neigh < 2) cycle
      allocate(Jmat(3*n_neigh, 6), rhsvec(3*n_neigh))
      row = 0
      do i = 1, n_neigh
        beta = mesh%elem(alpha)%neigh_by_vert(i)
        h = mesh%elem(beta)%coord - mesh%elem(alpha)%coord
        do ell = 1, 6
          E9 = 0.0_DOUBLE
          E9((basis_i(ell)-1)*3 + basis_j(ell)) = 1.0_DOUBLE
          E9((basis_j(ell)-1)*3 + basis_i(ell)) = 1.0_DOUBLE
          gcorr_alpha = contract_geom_h2(cell_ls_h1_cache(2)%m(:, alpha), E9)
          gcorr_beta  = contract_geom_h2(cell_ls_h1_cache(2)%m(:, beta), E9)
          Jmat(row*3+1, ell) = E9(1)*h(1) + E9(2)*h(2) + E9(3)*h(3) + 0.5_DOUBLE*(gcorr_beta(1)-gcorr_alpha(1))
          Jmat(row*3+2, ell) = E9(4)*h(1) + E9(5)*h(2) + E9(6)*h(3) + 0.5_DOUBLE*(gcorr_beta(2)-gcorr_alpha(2))
          Jmat(row*3+3, ell) = E9(7)*h(1) + E9(8)*h(2) + E9(9)*h(3) + 0.5_DOUBLE*(gcorr_beta(3)-gcorr_alpha(3))
        end do
        rhsvec(row*3+1:row*3+3) = grad_cell(:, beta) - grad_cell(:, alpha)
        row = row + 1
      end do
      normal_mat = matmul(transpose(Jmat), Jmat)
      normal_rhs(:, 1) = matmul(transpose(Jmat), rhsvec)
      call lu_factor_lapack(6_ENTIER, normal_mat, ipiv)
      call lu_solve_mat_lapack(6_ENTIER, normal_mat, ipiv, 1_ENTIER, normal_rhs)
      hess_cell(:, alpha) = normal_rhs(:, 1)
      deallocate(Jmat, rhsvec)
    end do
  end subroutine compute_cell_cls_2exact_hessian

  ! Driver: Haider, Croisille & Courbet (2011) Definition 1 CLS algorithm, faithfully cell-to-cell
  ! (reference implementation; k=1->2 only for now, eq. 13's final k-exactness correction and the
  ! k=2->3 step are not yet implemented). nc_in=1 only for now. grad_out/hess_out use the full flat
  ! (redundant) tensor convention, (dir2-1)*3*nc_in+(dir1-1)*nc_in+ic, matching hess_flat elsewhere.
  subroutine compute_derivatives_cls_classic(mesh, boundary_2d, phi, grad_out, hess_out)
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(1, mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(3, mesh%n_elems), intent(out) :: grad_out
    real(kind=DOUBLE), dimension(9, mesh%n_elems), intent(out) :: hess_out

    real(kind=DOUBLE), dimension(6, mesh%n_elems) :: hess_red
    integer(kind=ENTIER) :: alpha, ell, i_, j_, t1, t2
    integer(kind=ENTIER), parameter :: basis_i(6) = (/1,1,1,2,2,3/)
    integer(kind=ENTIER), parameter :: basis_j(6) = (/1,2,3,2,3,3/)

    call ensure_cell_ls_grad_mat_cache(mesh, boundary_2d)
    call compute_cell_ls_grad(mesh, 1_ENTIER, phi, grad_out)
    call compute_cell_cls_2exact_hessian(mesh, boundary_2d, grad_out, hess_red)

    hess_out = 0.0_DOUBLE
    do ell = 1, 6
      i_ = basis_i(ell); j_ = basis_j(ell)
      t1 = (j_-1)*3 + i_
      t2 = (i_-1)*3 + j_
      do alpha = 1, mesh%n_elems
        hess_out(t1, alpha) = hess_red(ell, alpha)
        hess_out(t2, alpha) = hess_red(ell, alpha)
      end do
    end do
  end subroutine compute_derivatives_cls_classic

  ! GENERAL reusable "apply the SAME grad-to-hess functional-identity OPERATOR (Haider eq. 16,
  ! m=1 case) to an ARBITRARY nc_in-component grad-like field", instead of the real 1-exact
  ! gradient -- since the operator's own geometry (Jmat/normal_mat, built from
  ! cell_ls_h1_cache(2) and cell centroid displacements alone) does NOT depend on which field is
  ! being processed, this is the SAME construction as compute_cell_cls_2exact_hessian, just
  ! batched over nc_in independent "grad-like" fields at once. This is the key building block
  ! that makes the k=2->3 step (needing beta's OWN hess operator applied to a geometric moment
  ! field, not phi) tractable: call this with grad_like_cell = grad-of-the-geometric-field.
  ! Output is the reduced 6-component (Hxx,Hxy,Hxz,Hyy,Hyz,Hzz) basis, times nc_in.
  subroutine apply_cell_cls_hess_operator(mesh, boundary_2d, nc_in, grad_like_cell, hess_like_out)
    use linear_solver_module, only: lu_factor_lapack, lu_solve_mat_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d
    integer(kind=ENTIER), intent(in) :: nc_in
    real(kind=DOUBLE), dimension(3*nc_in, mesh%n_elems), intent(in) :: grad_like_cell
    real(kind=DOUBLE), dimension(6*nc_in, mesh%n_elems), intent(out) :: hess_like_out

    integer(kind=ENTIER) :: alpha, i, beta, ell, row, n_neigh, ic, p
    real(kind=DOUBLE), dimension(:, :), allocatable :: Jmat, rhsvec, normal_rhs
    real(kind=DOUBLE), dimension(3) :: h, gcorr_alpha, gcorr_beta
    real(kind=DOUBLE), dimension(9) :: E9
    real(kind=DOUBLE), dimension(6, 6) :: normal_mat
    integer(kind=ENTIER), dimension(6) :: ipiv
    integer(kind=ENTIER), parameter :: basis_i(6) = (/1,1,1,2,2,3/)
    integer(kind=ENTIER), parameter :: basis_j(6) = (/1,2,3,2,3,3/)

    call ensure_cell_ls_h1_cache(mesh, boundary_2d, 2_ENTIER)

    hess_like_out = 0.0_DOUBLE

    do alpha = 1, mesh%n_elems
      n_neigh = mesh%elem(alpha)%n_neigh_by_vert
      if (n_neigh < 2) cycle
      allocate(Jmat(3*n_neigh, 6), rhsvec(3*n_neigh, nc_in))
      row = 0
      do i = 1, n_neigh
        beta = mesh%elem(alpha)%neigh_by_vert(i)
        h = mesh%elem(beta)%coord - mesh%elem(alpha)%coord
        do ell = 1, 6
          E9 = 0.0_DOUBLE
          E9((basis_i(ell)-1)*3 + basis_j(ell)) = 1.0_DOUBLE
          E9((basis_j(ell)-1)*3 + basis_i(ell)) = 1.0_DOUBLE
          gcorr_alpha = contract_geom_h2(cell_ls_h1_cache(2)%m(:, alpha), E9)
          gcorr_beta  = contract_geom_h2(cell_ls_h1_cache(2)%m(:, beta), E9)
          Jmat(row*3+1, ell) = E9(1)*h(1) + E9(2)*h(2) + E9(3)*h(3) + 0.5_DOUBLE*(gcorr_beta(1)-gcorr_alpha(1))
          Jmat(row*3+2, ell) = E9(4)*h(1) + E9(5)*h(2) + E9(6)*h(3) + 0.5_DOUBLE*(gcorr_beta(2)-gcorr_alpha(2))
          Jmat(row*3+3, ell) = E9(7)*h(1) + E9(8)*h(2) + E9(9)*h(3) + 0.5_DOUBLE*(gcorr_beta(3)-gcorr_alpha(3))
        end do
        do ic = 1, nc_in
          do p = 1, 3
            rhsvec(row*3+p, ic) = grad_like_cell((p-1)*nc_in+ic, beta) - grad_like_cell((p-1)*nc_in+ic, alpha)
          end do
        end do
        row = row + 1
      end do
      normal_mat = matmul(transpose(Jmat), Jmat)
      allocate(normal_rhs(6, nc_in))
      normal_rhs = matmul(transpose(Jmat), rhsvec)
      call lu_factor_lapack(6_ENTIER, normal_mat, ipiv)
      call lu_solve_mat_lapack(6_ENTIER, normal_mat, ipiv, nc_in, normal_rhs)
      do ell = 1, 6
        do ic = 1, nc_in
          hess_like_out((ell-1)*nc_in+ic, alpha) = normal_rhs(ell, ic)
        end do
      end do
      deallocate(Jmat, rhsvec, normal_rhs)
    end do
  end subroutine apply_cell_cls_hess_operator

  ! Builds cell_cls_hess_of_z3_cache(X) := X's OWN hess operator (apply_cell_cls_hess_operator)
  ! applied to the geometric field {z_{X,gamma}^(3)}_gamma (shifted_cell_moment_full) instead of
  ! phi -- Haider's w_beta^(2|2) applied to z_beta^(3) (eq. 15-16, m=2 case), the piece needed to
  ! build a genuinely 3-exact third derivative from a 2-exact Hessian (Definition 1, step 2,
  ! k=2->3) that a simple reapplication of grad's OWN operator (as used for k=1->2) cannot give.
  ! For each X, first fits the geometric field's OWN "grad" at X and each of X's neigh_by_vert
  ! neighbors (a local, bounded 2-ring computation via cell_ls_mat_inv_cache), then applies X's
  ! hess operator to those local grad values. Pure geometry, expensive (O(n_elems * n_neigh^2)),
  ! cached once per mesh. Output: reduced 6-component hess basis x 27-component z basis (162
  ! total), layout (ell_hess-1)*27+ell_z.
  subroutine ensure_cell_cls_hess_of_z3_cache(mesh, boundary_2d)
    use linear_solver_module, only: lu_factor_lapack, lu_solve_mat_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d

    integer(kind=ENTIER) :: X, Y, iy, n1, ny, j, gamma
    real(kind=DOUBLE), dimension(3) :: h
    real(kind=DOUBLE), dimension(27) :: z_gamma_at_X, z_Y_at_X
    real(kind=DOUBLE), dimension(:, :), allocatable :: rhs, grad_of_z3_x, hess_of_z3_x

    call ensure_cell_ls_grad_mat_cache(mesh, boundary_2d)
    call ensure_cell_ls_h1_cache(mesh, boundary_2d, 2_ENTIER)
    call ensure_cell_moment_cache(mesh, 3_ENTIER)

    if (cell_cls_hess_of_z3_n_elems == mesh%n_elems .and. (cell_cls_hess_of_z3_boundary_2d .eqv. boundary_2d)) return

    if (allocated(cell_cls_hess_of_z3_cache)) deallocate(cell_cls_hess_of_z3_cache)
    allocate(cell_cls_hess_of_z3_cache(162, mesh%n_elems))
    cell_cls_hess_of_z3_cache = 0.0_DOUBLE

    do X = 1, mesh%n_elems
      n1 = mesh%elem(X)%n_neigh_by_vert
      if (n1 < 2) cycle

      ! grad-of-z3_X at X itself (iy=0) and at each of X's own neigh_by_vert neighbors (iy=1..n1).
      allocate(grad_of_z3_x(3*27, 0:n1))
      do iy = 0, n1
        Y = merge(X, mesh%elem(X)%neigh_by_vert(iy), iy == 0)
        ny = mesh%elem(Y)%n_neigh_by_vert
        allocate(rhs(3, 27))
        rhs = 0.0_DOUBLE
        call shifted_cell_moment_full(Y, mesh%elem(Y)%coord - mesh%elem(X)%coord, 3_ENTIER, z_Y_at_X)
        do j = 1, ny
          gamma = mesh%elem(Y)%neigh_by_vert(j)
          h = mesh%elem(gamma)%coord - mesh%elem(Y)%coord
          call shifted_cell_moment_full(gamma, mesh%elem(gamma)%coord - mesh%elem(X)%coord, 3_ENTIER, z_gamma_at_X)
          rhs(1, :) = rhs(1, :) + h(1) * (z_gamma_at_X - z_Y_at_X)
          rhs(2, :) = rhs(2, :) + h(2) * (z_gamma_at_X - z_Y_at_X)
          rhs(3, :) = rhs(3, :) + h(3) * (z_gamma_at_X - z_Y_at_X)
        end do
        block
          real(kind=DOUBLE), dimension(3, 27) :: res
          res = matmul(cell_ls_mat_inv_cache(:, :, Y), rhs)
          grad_of_z3_x((0)*27+1:(0)*27+27, iy) = res(1, :)
          grad_of_z3_x((1)*27+1:(1)*27+27, iy) = res(2, :)
          grad_of_z3_x((2)*27+1:(2)*27+27, iy) = res(3, :)
        end block
        deallocate(rhs)
      end do

      ! Apply X's OWN hess operator to grad_of_z3_x (nc_in=27), batched over the 27 z-components.
      block
        real(kind=DOUBLE), dimension(:, :), allocatable :: grad_like_local
        real(kind=DOUBLE), dimension(:, :), allocatable :: Jmat, rhsvec, normal_rhs
        real(kind=DOUBLE), dimension(3) :: hh, gcorr_alpha, gcorr_beta
        real(kind=DOUBLE), dimension(9) :: E9
        real(kind=DOUBLE), dimension(6, 6) :: normal_mat
        integer(kind=ENTIER), dimension(6) :: ipiv
        integer(kind=ENTIER), parameter :: basis_i(6) = (/1,1,1,2,2,3/)
        integer(kind=ENTIER), parameter :: basis_j(6) = (/1,2,3,2,3,3/)
        integer(kind=ENTIER) :: i2, beta2, ell2, row2, ic2, p2

        allocate(Jmat(3*n1, 6), rhsvec(3*n1, 27))
        row2 = 0
        do i2 = 1, n1
          beta2 = mesh%elem(X)%neigh_by_vert(i2)
          hh = mesh%elem(beta2)%coord - mesh%elem(X)%coord
          do ell2 = 1, 6
            E9 = 0.0_DOUBLE
            E9((basis_i(ell2)-1)*3 + basis_j(ell2)) = 1.0_DOUBLE
            E9((basis_j(ell2)-1)*3 + basis_i(ell2)) = 1.0_DOUBLE
            gcorr_alpha = contract_geom_h2(cell_ls_h1_cache(2)%m(:, X), E9)
            gcorr_beta  = contract_geom_h2(cell_ls_h1_cache(2)%m(:, beta2), E9)
            Jmat(row2*3+1, ell2) = E9(1)*hh(1)+E9(2)*hh(2)+E9(3)*hh(3) + 0.5_DOUBLE*(gcorr_beta(1)-gcorr_alpha(1))
            Jmat(row2*3+2, ell2) = E9(4)*hh(1)+E9(5)*hh(2)+E9(6)*hh(3) + 0.5_DOUBLE*(gcorr_beta(2)-gcorr_alpha(2))
            Jmat(row2*3+3, ell2) = E9(7)*hh(1)+E9(8)*hh(2)+E9(9)*hh(3) + 0.5_DOUBLE*(gcorr_beta(3)-gcorr_alpha(3))
          end do
          do ic2 = 1, 27
            do p2 = 1, 3
              rhsvec(row2*3+p2, ic2) = grad_of_z3_x((p2-1)*27+ic2, i2) - grad_of_z3_x((p2-1)*27+ic2, 0)
            end do
          end do
          row2 = row2 + 1
        end do
        normal_mat = matmul(transpose(Jmat), Jmat)
        allocate(normal_rhs(6, 27))
        normal_rhs = matmul(transpose(Jmat), rhsvec)
        call lu_factor_lapack(6_ENTIER, normal_mat, ipiv)
        call lu_solve_mat_lapack(6_ENTIER, normal_mat, ipiv, 27_ENTIER, normal_rhs)
        do ell2 = 1, 6
          do ic2 = 1, 27
            cell_cls_hess_of_z3_cache((ell2-1)*27+ic2, X) = normal_rhs(ell2, ic2)
          end do
        end do
        deallocate(Jmat, rhsvec, normal_rhs)
      end block

      deallocate(grad_of_z3_x)
    end do

    cell_cls_hess_of_z3_n_elems = mesh%n_elems
    cell_cls_hess_of_z3_boundary_2d = boundary_2d
  end subroutine ensure_cell_cls_hess_of_z3_cache

  ! Cell-to-cell CLS third derivative (Haider, Croisille & Courbet 2011, Definition 1, step 2,
  ! k=2->3): builds a genuinely 3-exact third derivative directly at cell centers from a 2-exact
  ! Hessian (hess_red, reduced 6-component Hxx,Hxy,Hxz,Hyy,Hyz,Hzz, from
  ! compute_cell_cls_2exact_hessian) and cell_cls_hess_of_z3_cache (beta's OWN hess operator
  ! applied to z_beta^(3), NOT grad's operator -- the genuinely faithful eq. 16 construction,
  ! unlike compute_node_derivative_bias's approximation of always reusing the order-1 operator).
  ! third_out is the full flat (redundant) 27-component tensor; NOT explicitly symmetrized (each
  ! of the 27 unknowns is fit independently against whichever equations reference it, so
  ! permutation-equivalent entries may differ by a small residual amount rather than being
  ! forced exactly equal).
  ! third_red is the GENUINELY REDUCED 10-component basis of S^3(R^3) (independent components of a
  ! symmetric rank-3 3D tensor: xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz), NOT the full-flat
  ! 27-component redundant array used elsewhere in this module. Using the redundant basis for an
  ! UNKNOWN being solved for (as opposed to a KNOWN geometric field, where redundancy is harmless)
  ! leaves the normal-equations system rank-deficient: its null space lets the Moore-Penrose
  ! pseudo-inverse return the MINIMUM-NORM solution among many consistent ones, not the unique
  ! true (symmetric) value -- a real, load-bearing bug, not a cosmetic asymmetry, found by
  ! rereading Haider, Croisille & Courbet 2011 in full: their own S^m(R^d) notation IS the reduced
  ! symmetric tensor space, never a redundant flat array.
  subroutine compute_cell_cls_3exact_third(mesh, boundary_2d, hess_red, third_red)
    use linear_solver_module, only: pseudo_inverse_inplace_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(in) :: hess_red
    real(kind=DOUBLE), dimension(10, mesh%n_elems), intent(out) :: third_red

    integer(kind=ENTIER) :: alpha, i, beta, n_neigh, row, ell_h, ell_b, c, p, q
    real(kind=DOUBLE), dimension(:, :), allocatable :: Jmat, normal_rhs
    real(kind=DOUBLE), dimension(:), allocatable :: rhsvec
    real(kind=DOUBLE), dimension(3) :: h
    real(kind=DOUBLE), dimension(10, 10) :: normal_mat
    integer(kind=ENTIER), parameter :: basis_i(6) = (/1,1,1,2,2,3/)
    integer(kind=ENTIER), parameter :: basis_j(6) = (/1,2,3,2,3,3/)
    ! Sorted (a,b,c) triple for each of the 10 canonical third-derivative components, and their
    ! full-flat 27-index (c-1)*9+(b-1)*3+a (matching cell_cls_hess_of_z3_cache's own z-input
    ! layout) -- z is a KNOWN geometric field there, so picking one canonical representative per
    ! permutation-equivalent group is exact, no averaging needed.
    integer(kind=ENTIER), parameter :: red10_a(10) = (/1,1,1,1,1,1,2,2,2,3/)
    integer(kind=ENTIER), parameter :: red10_b(10) = (/1,1,1,2,2,3,2,2,3,3/)
    integer(kind=ENTIER), parameter :: red10_c(10) = (/1,2,3,2,3,3,2,3,3,3/)
    integer(kind=ENTIER), parameter :: red10_full27(10) = (/1,10,19,13,22,25,14,23,26,27/)

    call ensure_cell_cls_hess_of_z3_cache(mesh, boundary_2d)

    third_red = 0.0_DOUBLE

    do alpha = 1, mesh%n_elems
      n_neigh = mesh%elem(alpha)%n_neigh_by_vert
      ! 10 unknowns, 6 equations/neighbor: n_neigh>=3 already over-determines the system.
      if (n_neigh < 5) cycle
      allocate(Jmat(6*n_neigh, 10), rhsvec(6*n_neigh))
      Jmat = 0.0_DOUBLE
      row = 0
      do i = 1, n_neigh
        beta = mesh%elem(alpha)%neigh_by_vert(i)
        h = mesh%elem(beta)%coord - mesh%elem(alpha)%coord
        do ell_h = 1, 6
          p = basis_i(ell_h); q = basis_j(ell_h)
          ! "b.h" (eq. 9): (b.h)_{p,q} = sum_c b_{p,q,c} h(c) -- b_{p,q,c} maps to whichever of
          ! the 10 canonical unknowns matches sorted(p,q,c).
          do c = 1, 3
            ell_b = canon10_index(p, q, c)
            Jmat(row*6+ell_h, ell_b) = Jmat(row*6+ell_h, ell_b) + h(c)
          end do
          ! Geometric correction term (eq. 16): (1/3!) [HESS_OF_Z3(beta)-HESS_OF_Z3(alpha)],
          ! read off at the SAME canonical z-representative used for each of the 10 unknowns.
          do ell_b = 1, 10
            Jmat(row*6+ell_h, ell_b) = Jmat(row*6+ell_h, ell_b) + (1.0_DOUBLE/6.0_DOUBLE) * &
              (cell_cls_hess_of_z3_cache((ell_h-1)*27+red10_full27(ell_b), beta) &
               - cell_cls_hess_of_z3_cache((ell_h-1)*27+red10_full27(ell_b), alpha))
          end do
        end do
        rhsvec(row*6+1:row*6+6) = hess_red(:, beta) - hess_red(:, alpha)
        row = row + 1
      end do
      normal_mat = matmul(transpose(Jmat), Jmat)
      allocate(normal_rhs(10, 1))
      normal_rhs(:, 1) = matmul(transpose(Jmat), rhsvec)
      call pseudo_inverse_inplace_lapack(10_ENTIER, normal_mat)
      third_red(:, alpha) = matmul(normal_mat, normal_rhs(:, 1))
      deallocate(Jmat, rhsvec, normal_rhs)
    end do
  end subroutine compute_cell_cls_3exact_third

  ! Sorts (a,b,c) and returns the matching index (1..10) into the canonical S^3(R^3) basis
  ! (xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz). Module-level so every routine building or applying
  ! a third-derivative operator uses the SAME canonical convention.
  pure function canon10_index(a_in, b_in, c_in) result(ell)
    implicit none
    integer(kind=ENTIER), intent(in) :: a_in, b_in, c_in
    integer(kind=ENTIER) :: ell
    integer(kind=ENTIER), dimension(3) :: s
    integer(kind=ENTIER) :: tmp
    integer(kind=ENTIER), parameter :: red10_a(10) = (/1,1,1,1,1,1,2,2,2,3/)
    integer(kind=ENTIER), parameter :: red10_b(10) = (/1,1,1,2,2,3,2,2,3,3/)
    integer(kind=ENTIER), parameter :: red10_c(10) = (/1,2,3,2,3,3,2,3,3,3/)
    s = (/a_in, b_in, c_in/)
    if (s(1) > s(2)) then; tmp=s(1); s(1)=s(2); s(2)=tmp; end if
    if (s(2) > s(3)) then; tmp=s(2); s(2)=s(3); s(3)=tmp; end if
    if (s(1) > s(2)) then; tmp=s(1); s(1)=s(2); s(2)=tmp; end if
    do ell = 1, 10
      if (red10_a(ell) == s(1) .and. red10_b(ell) == s(2) .and. red10_c(ell) == s(3)) return
    end do
    ell = -1
  end function canon10_index

  ! Sorts (a,b,c,d) and returns the matching index (1..15) into the canonical S^4(R^3) basis
  ! (xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,xyyz,xyzz,xzzz,yyyy,yyyz,yyzz,yzzz,zzzz).
  pure function canon15_index(a_in, b_in, c_in, d_in) result(ell)
    implicit none
    integer(kind=ENTIER), intent(in) :: a_in, b_in, c_in, d_in
    integer(kind=ENTIER) :: ell
    integer(kind=ENTIER), dimension(4) :: s
    integer(kind=ENTIER) :: tmp, i, j
    integer(kind=ENTIER), parameter :: red15_a(15) = (/1,1,1,1,1,1,1,1,1,1,2,2,2,2,3/)
    integer(kind=ENTIER), parameter :: red15_b(15) = (/1,1,1,1,1,1,2,2,2,3,2,2,2,3,3/)
    integer(kind=ENTIER), parameter :: red15_c(15) = (/1,1,1,2,2,3,2,2,3,3,2,2,3,3,3/)
    integer(kind=ENTIER), parameter :: red15_d(15) = (/1,2,3,2,3,3,2,3,3,3,2,3,3,3,3/)
    s = (/a_in, b_in, c_in, d_in/)
    do i = 1, 3
      do j = 1, 4-i
        if (s(j) > s(j+1)) then; tmp=s(j); s(j)=s(j+1); s(j+1)=tmp; end if
      end do
    end do
    do ell = 1, 15
      if (red15_a(ell) == s(1) .and. red15_b(ell) == s(2) .and. red15_c(ell) == s(3) &
          .and. red15_d(ell) == s(4)) return
    end do
    ell = -1
  end function canon15_index

  ! Expands a reduced 10-component (S^3(R^3): xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz) tensor
  ! into the full flat 27-component convention, for feeding routines (contract_grad_node_bias)
  ! that expect the general full-flat layout used elsewhere in this module. Since third_red's 10
  ! components are already the UNIQUE, exactly-symmetric solution (no redundant unknowns solved
  ! for), this expansion is exact by construction -- no averaging needed.
  pure function third_red_to_full27(third_red, n_elems) result(third_full)
    implicit none
    integer(kind=ENTIER), intent(in) :: n_elems
    real(kind=DOUBLE), dimension(10, n_elems), intent(in) :: third_red
    real(kind=DOUBLE), dimension(27, n_elems) :: third_full
    integer(kind=ENTIER) :: ell, a, b, c, mask
    integer(kind=ENTIER), parameter :: red10_a(10) = (/1,1,1,1,1,1,2,2,2,3/)
    integer(kind=ENTIER), parameter :: red10_b(10) = (/1,1,1,2,2,3,2,2,3,3/)
    integer(kind=ENTIER), parameter :: red10_c(10) = (/1,2,3,2,3,3,2,3,3,3/)

    third_full = 0.0_DOUBLE
    do ell = 1, 10
      a = red10_a(ell); b = red10_b(ell); c = red10_c(ell)
      ! Scatter into every raw permutation of (a,b,c) (1, 3, or 6 slots depending on repeats).
      do mask = 0, 5
        block
          integer(kind=ENTIER), dimension(3) :: perm
          select case (mask)
          case (0); perm = (/a,b,c/)
          case (1); perm = (/a,c,b/)
          case (2); perm = (/b,a,c/)
          case (3); perm = (/b,c,a/)
          case (4); perm = (/c,a,b/)
          case (5); perm = (/c,b,a/)
          end select
          third_full((perm(3)-1)*9+(perm(2)-1)*3+perm(1), :) = third_red(ell, :)
        end block
      end do
    end do
  end function third_red_to_full27

  ! Averages every permutation-equivalent group of components of a full-flat (27, n_elems)
  ! rank-3 tensor (flat index (c-1)*9+(b-1)*3+a for (a,b,c) in 1..3) so the result is EXACTLY
  ! symmetric under any permutation of its 3 indices, instead of merely approximately so.
  subroutine symmetrize_third_flat(third_out, n_elems)
    implicit none

    integer(kind=ENTIER), intent(in) :: n_elems
    real(kind=DOUBLE), dimension(27, n_elems), intent(inout) :: third_out

    integer(kind=ENTIER) :: a, b, c, idx0, k, m, np, idxk
    integer(kind=ENTIER), dimension(6, 3) :: perms
    integer(kind=ENTIER), dimension(6) :: perm_idx
    logical :: dup
    logical, dimension(27) :: done
    real(kind=DOUBLE), dimension(n_elems) :: avg

    done = .false.
    do a = 1, 3
      do b = 1, 3
        do c = 1, 3
          idx0 = (c-1)*9 + (b-1)*3 + a
          if (done(idx0)) cycle
          perms(1, :) = (/a, b, c/)
          perms(2, :) = (/a, c, b/)
          perms(3, :) = (/b, a, c/)
          perms(4, :) = (/b, c, a/)
          perms(5, :) = (/c, a, b/)
          perms(6, :) = (/c, b, a/)

          np = 0
          do k = 1, 6
            idxk = (perms(k,3)-1)*9 + (perms(k,2)-1)*3 + perms(k,1)
            dup = .false.
            do m = 1, np
              if (perm_idx(m) == idxk) dup = .true.
            end do
            if (.not. dup) then
              np = np + 1
              perm_idx(np) = idxk
            end if
          end do

          avg = 0.0_DOUBLE
          do k = 1, np
            avg = avg + third_out(perm_idx(k), :)
          end do
          avg = avg / real(np, kind=DOUBLE)
          do k = 1, np
            third_out(perm_idx(k), :) = avg
            done(perm_idx(k)) = .true.
          end do
        end do
      end do
    end do
  end subroutine symmetrize_third_flat

  ! Driver: full Haider, Croisille & Courbet (2011) Definition 1 CLS algorithm, faithfully
  ! cell-to-cell, up to the third derivative (grad 1-exact, hess 2-exact via eq. 16's functional
  ! identity, third 3-exact via the SAME identity reusing hess's OWN operator, not grad's).
  ! nc_in=1 only. grad_out/hess_out use the full flat (redundant) tensor convention; third_out is
  ! full flat too and not explicitly symmetrized (see compute_cell_cls_3exact_third).
  subroutine compute_derivatives_cls_classic_order3(mesh, boundary_2d, phi, grad_out, hess_out, third_out)
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(1, mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(3, mesh%n_elems), intent(out) :: grad_out
    real(kind=DOUBLE), dimension(9, mesh%n_elems), intent(out) :: hess_out
    real(kind=DOUBLE), dimension(27, mesh%n_elems), intent(out) :: third_out

    real(kind=DOUBLE), dimension(6, mesh%n_elems) :: hess_red
    real(kind=DOUBLE), dimension(10, mesh%n_elems) :: third_red
    integer(kind=ENTIER) :: alpha, ell, i_, j_, t1, t2
    integer(kind=ENTIER), parameter :: basis_i(6) = (/1,1,1,2,2,3/)
    integer(kind=ENTIER), parameter :: basis_j(6) = (/1,2,3,2,3,3/)

    call ensure_cell_ls_grad_mat_cache(mesh, boundary_2d)
    call compute_cell_ls_grad(mesh, 1_ENTIER, phi, grad_out)
    call compute_cell_cls_2exact_hessian(mesh, boundary_2d, grad_out, hess_red)

    hess_out = 0.0_DOUBLE
    do ell = 1, 6
      i_ = basis_i(ell); j_ = basis_j(ell)
      t1 = (j_-1)*3 + i_
      t2 = (i_-1)*3 + j_
      do alpha = 1, mesh%n_elems
        hess_out(t1, alpha) = hess_red(ell, alpha)
        hess_out(t2, alpha) = hess_red(ell, alpha)
      end do
    end do

    call compute_cell_cls_3exact_third(mesh, boundary_2d, hess_red, third_red)
    third_out = third_red_to_full27(third_red, mesh%n_elems)

    ! Definition 1, step 3 (eq. 13): correct grad (k=1) and hess (k=2) using third (k+1=3), now
    ! that it is available, so grad/hess ALSO reach 3-exactness instead of stopping at their own
    ! 1-/2-exact construction. For a genuinely cubic field this should drive grad and hess to
    ! (near) machine precision.
    call apply_cell_cls_eq13_correction(mesh, boundary_2d, grad_out, hess_red, third_out)

    hess_out = 0.0_DOUBLE
    do ell = 1, 6
      i_ = basis_i(ell); j_ = basis_j(ell)
      t1 = (j_-1)*3 + i_
      t2 = (i_-1)*3 + j_
      do alpha = 1, mesh%n_elems
        hess_out(t1, alpha) = hess_red(ell, alpha)
        hess_out(t2, alpha) = hess_red(ell, alpha)
      end do
    end do
  end subroutine compute_derivatives_cls_classic_order3

  ! Haider, Croisille & Courbet (2011) Definition 1, step 3 (eq. 13): corrects grad (k=1, in
  ! place) and hess_red (k=2, reduced 6-component, in place) using the now-available third
  ! derivative, reusing the SAME geometric operators already built for the k->k+1 recursion
  ! (cell_ls_h1_cache for grad's own correction, cell_cls_hess_of_z3_cache for hess's own
  ! correction -- NOT always grad's operator, unlike compute_node_derivative_bias's
  ! approximation for the vertex-based aho_gg/aho_ls chain).
  subroutine apply_cell_cls_eq13_correction(mesh, boundary_2d, grad_out, hess_red, third_out)
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(3, mesh%n_elems), intent(inout) :: grad_out
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(inout) :: hess_red
    real(kind=DOUBLE), dimension(27, mesh%n_elems), intent(in) :: third_out

    real(kind=DOUBLE), dimension(3, mesh%n_elems) :: bias_grad_m2, bias_grad_m3
    real(kind=DOUBLE), dimension(6, mesh%n_elems) :: bias_hess
    integer(kind=ENTIER) :: alpha, ell_h, ell_z

    call ensure_cell_ls_h1_cache(mesh, boundary_2d, 3_ENTIER)
    call ensure_cell_cls_hess_of_z3_cache(mesh, boundary_2d)

    ! Hess correction FIRST: (1/3!) [hess's own operator on z^3] : third (Haider eq. 13, k=2) --
    ! so grad's OWN correction below reuses the NOW-corrected hess, not the raw 2-exact one.
    bias_hess = 0.0_DOUBLE
    do alpha = 1, mesh%n_elems
      do ell_h = 1, 6
        do ell_z = 1, 27
          bias_hess(ell_h, alpha) = bias_hess(ell_h, alpha) &
            + cell_cls_hess_of_z3_cache((ell_h-1)*27+ell_z, alpha) * third_out(ell_z, alpha)
        end do
      end do
    end do
    hess_red = hess_red - bias_hess/6.0_DOUBLE

    ! Grad correction: (1/2!) H_2^(1):hess + (1/3!) H_3^(1):third (Haider eq. 13, k=1), using the
    ! now-corrected hess_red.
    call contract_grad_node_bias(1_ENTIER, 2_ENTIER, mesh%n_elems, hess_red_to_full9(hess_red, mesh%n_elems), &
      cell_ls_h1_cache(2)%m, bias_grad_m2)
    call contract_grad_node_bias(1_ENTIER, 3_ENTIER, mesh%n_elems, third_out, cell_ls_h1_cache(3)%m, bias_grad_m3)
    grad_out = grad_out - bias_grad_m2/2.0_DOUBLE - bias_grad_m3/6.0_DOUBLE
  end subroutine apply_cell_cls_eq13_correction

  ! Expands a reduced 6-component (Hxx,Hxy,Hxz,Hyy,Hyz,Hzz) hess field into the full flat
  ! 9-component convention, for feeding contract_grad_node_bias (which expects the general
  ! full-flat layout used everywhere else in the module).
  pure function hess_red_to_full9(hess_red, n_elems) result(hess_full)
    implicit none
    integer(kind=ENTIER), intent(in) :: n_elems
    real(kind=DOUBLE), dimension(6, n_elems), intent(in) :: hess_red
    real(kind=DOUBLE), dimension(9, n_elems) :: hess_full
    integer(kind=ENTIER) :: ell, i_, j_, t1, t2
    integer(kind=ENTIER), parameter :: basis_i(6) = (/1,1,1,2,2,3/)
    integer(kind=ENTIER), parameter :: basis_j(6) = (/1,2,3,2,3,3/)

    hess_full = 0.0_DOUBLE
    do ell = 1, 6
      i_ = basis_i(ell); j_ = basis_j(ell)
      t1 = (j_-1)*3 + i_
      t2 = (i_-1)*3 + j_
      hess_full(t1, :) = hess_red(ell, :)
      hess_full(t2, :) = hess_red(ell, :)
    end do
  end function hess_red_to_full9

  ! ============================================================================================
  ! Order-5 extension (Haider, Croisille & Courbet 2011, Definition 1, k=3->4 step): builds a
  ! genuinely 4-exact 4th derivative from the 3-exact third derivative, needing "third's OWN
  ! operator applied to z^(4)" (eq. 16, m=3 case) -- ONE MORE nesting level than the k=2->3 step
  ! (which needed "hess's own operator applied to z^(3)"). Implemented ON DEMAND (no persistent
  ! module-level cache): every quantity is recomputed fresh from mesh geometry + cell_ls_mat_inv_
  ! cache + cell_ls_h1_cache(2) + cell_cls_hess_of_z3_cache each call, mirroring exactly how a
  ! human would hand-evaluate eq. 15-18 recursively. This is expensive (each call to
  ! third_of_geom_field_at costs O(n_neigh^3): a grad-fit at every point of a 2-ring, a hess-fit
  ! at every point of a 1-ring, then one more solve) -- acceptable for validating correctness on
  ! a small mesh, not yet a production-ready order-5 implementation.
  ! ============================================================================================

  ! id_eval's own 1-exact grad operator (Definition 1, step 1) applied to the KNOWN geometric
  ! field {z_{id_home,gamma}^(m_order)}_gamma (shifted_cell_moment_full) instead of phi, using
  ! id_eval's own neigh_by_vert stencil. Full-flat output (3, 3**m_order).
  subroutine grad_of_geom_field_at(mesh, id_home, m_order, id_eval, g)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_home, m_order, id_eval
    real(kind=DOUBLE), dimension(3, 3_ENTIER**m_order), intent(out) :: g

    integer(kind=ENTIER) :: n_eval, j, gamma, nc
    real(kind=DOUBLE), dimension(3) :: h
    real(kind=DOUBLE), dimension(3_ENTIER**m_order) :: z_gamma, z_eval
    real(kind=DOUBLE), dimension(3, 3_ENTIER**m_order) :: rhs

    nc = 3_ENTIER**m_order
    n_eval = mesh%elem(id_eval)%n_neigh_by_vert
    rhs = 0.0_DOUBLE
    call shifted_cell_moment_full(id_eval, mesh%elem(id_eval)%coord - mesh%elem(id_home)%coord, m_order, z_eval)
    do j = 1, n_eval
      gamma = mesh%elem(id_eval)%neigh_by_vert(j)
      h = mesh%elem(gamma)%coord - mesh%elem(id_eval)%coord
      call shifted_cell_moment_full(gamma, mesh%elem(gamma)%coord - mesh%elem(id_home)%coord, m_order, z_gamma)
      rhs(1, :) = rhs(1, :) + h(1) * (z_gamma - z_eval)
      rhs(2, :) = rhs(2, :) + h(2) * (z_gamma - z_eval)
      rhs(3, :) = rhs(3, :) + h(3) * (z_gamma - z_eval)
    end do
    g = matmul(cell_ls_mat_inv_cache(:, :, id_eval), rhs)
  end subroutine grad_of_geom_field_at

  ! id_eval's own 2-exact hess operator (eq. 16, m=1 case) applied to the SAME geometric field,
  ! via grad_of_geom_field_at evaluated over id_eval's own neigh_by_vert stencil. Reduced
  ! 6-component output (6, 3**m_order).
  subroutine hess_of_geom_field_at(mesh, boundary_2d, id_home, m_order, id_eval, h6)
    use linear_solver_module, only: pseudo_inverse_inplace_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d
    integer(kind=ENTIER), intent(in) :: id_home, m_order, id_eval
    real(kind=DOUBLE), dimension(6, 3_ENTIER**m_order), intent(out) :: h6

    integer(kind=ENTIER) :: nc, n_neigh, i, beta, ell, row
    real(kind=DOUBLE), dimension(3) :: hh, gcorr_a, gcorr_b
    real(kind=DOUBLE), dimension(9) :: E9
    real(kind=DOUBLE), dimension(:, :), allocatable :: Jmat, rhsvec, normal_rhs, g_eval, g_beta
    real(kind=DOUBLE), dimension(6, 6) :: normal_mat
    integer(kind=ENTIER), parameter :: basis_i(6) = (/1,1,1,2,2,3/)
    integer(kind=ENTIER), parameter :: basis_j(6) = (/1,2,3,2,3,3/)

    nc = 3_ENTIER**m_order
    n_neigh = mesh%elem(id_eval)%n_neigh_by_vert
    allocate(g_eval(3, nc), g_beta(3, nc))
    call grad_of_geom_field_at(mesh, id_home, m_order, id_eval, g_eval)
    allocate(Jmat(3*n_neigh, 6), rhsvec(3*n_neigh, nc))
    row = 0
    do i = 1, n_neigh
      beta = mesh%elem(id_eval)%neigh_by_vert(i)
      hh = mesh%elem(beta)%coord - mesh%elem(id_eval)%coord
      do ell = 1, 6
        E9 = 0.0_DOUBLE
        E9((basis_i(ell)-1)*3 + basis_j(ell)) = 1.0_DOUBLE
        E9((basis_j(ell)-1)*3 + basis_i(ell)) = 1.0_DOUBLE
        gcorr_a = contract_geom_h2(cell_ls_h1_cache(2)%m(:, id_eval), E9)
        gcorr_b = contract_geom_h2(cell_ls_h1_cache(2)%m(:, beta), E9)
        Jmat(row*3+1, ell) = E9(1)*hh(1)+E9(2)*hh(2)+E9(3)*hh(3) + 0.5_DOUBLE*(gcorr_b(1)-gcorr_a(1))
        Jmat(row*3+2, ell) = E9(4)*hh(1)+E9(5)*hh(2)+E9(6)*hh(3) + 0.5_DOUBLE*(gcorr_b(2)-gcorr_a(2))
        Jmat(row*3+3, ell) = E9(7)*hh(1)+E9(8)*hh(2)+E9(9)*hh(3) + 0.5_DOUBLE*(gcorr_b(3)-gcorr_a(3))
      end do
      call grad_of_geom_field_at(mesh, id_home, m_order, beta, g_beta)
      rhsvec(row*3+1, :) = g_beta(1, :) - g_eval(1, :)
      rhsvec(row*3+2, :) = g_beta(2, :) - g_eval(2, :)
      rhsvec(row*3+3, :) = g_beta(3, :) - g_eval(3, :)
      row = row + 1
    end do
    normal_mat = matmul(transpose(Jmat), Jmat)
    allocate(normal_rhs(6, nc))
    normal_rhs = matmul(transpose(Jmat), rhsvec)
    call pseudo_inverse_inplace_lapack(6_ENTIER, normal_mat)
    h6 = matmul(normal_mat, normal_rhs)
    deallocate(Jmat, rhsvec, normal_rhs, g_eval, g_beta)
  end subroutine hess_of_geom_field_at

  ! id_eval's own 3-exact third operator (eq. 16, m=2 case, reduced-10 basis) applied to the SAME
  ! geometric field, via hess_of_geom_field_at evaluated over id_eval's own neigh_by_vert
  ! stencil, using cell_cls_hess_of_z3_cache for the SAME geometric correction term used to build
  ! the real third derivative. Reduced 10-component output (10, 3**m_order).
  subroutine third_of_geom_field_at(mesh, boundary_2d, id_home, m_order, id_eval, t10)
    use linear_solver_module, only: pseudo_inverse_inplace_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d
    integer(kind=ENTIER), intent(in) :: id_home, m_order, id_eval
    real(kind=DOUBLE), dimension(10, 3_ENTIER**m_order), intent(out) :: t10

    integer(kind=ENTIER) :: nc, n_neigh, i, beta, ell_h, ell_b, c, p, q, row
    real(kind=DOUBLE), dimension(3) :: hh
    real(kind=DOUBLE), dimension(:, :), allocatable :: Jmat, rhsvec, normal_rhs, hess_eval, hess_beta
    real(kind=DOUBLE), dimension(10, 10) :: normal_mat
    integer(kind=ENTIER), parameter :: basis_i(6) = (/1,1,1,2,2,3/)
    integer(kind=ENTIER), parameter :: basis_j(6) = (/1,2,3,2,3,3/)
    integer(kind=ENTIER), parameter :: red10_full27(10) = (/1,10,19,13,22,25,14,23,26,27/)

    nc = 3_ENTIER**m_order
    n_neigh = mesh%elem(id_eval)%n_neigh_by_vert
    allocate(hess_eval(6, nc), hess_beta(6, nc))
    call hess_of_geom_field_at(mesh, boundary_2d, id_home, m_order, id_eval, hess_eval)
    allocate(Jmat(6*n_neigh, 10), rhsvec(6*n_neigh, nc))
    Jmat = 0.0_DOUBLE
    row = 0
    do i = 1, n_neigh
      beta = mesh%elem(id_eval)%neigh_by_vert(i)
      hh = mesh%elem(beta)%coord - mesh%elem(id_eval)%coord
      do ell_h = 1, 6
        p = basis_i(ell_h); q = basis_j(ell_h)
        do c = 1, 3
          ell_b = canon10_index(p, q, c)
          Jmat(row*6+ell_h, ell_b) = Jmat(row*6+ell_h, ell_b) + hh(c)
        end do
        do ell_b = 1, 10
          Jmat(row*6+ell_h, ell_b) = Jmat(row*6+ell_h, ell_b) + (1.0_DOUBLE/6.0_DOUBLE) * &
            (cell_cls_hess_of_z3_cache((ell_h-1)*27+red10_full27(ell_b), beta) &
             - cell_cls_hess_of_z3_cache((ell_h-1)*27+red10_full27(ell_b), id_eval))
        end do
      end do
      call hess_of_geom_field_at(mesh, boundary_2d, id_home, m_order, beta, hess_beta)
      rhsvec(row*6+1:row*6+6, :) = hess_beta - hess_eval
      row = row + 1
    end do
    normal_mat = matmul(transpose(Jmat), Jmat)
    allocate(normal_rhs(10, nc))
    normal_rhs = matmul(transpose(Jmat), rhsvec)
    call pseudo_inverse_inplace_lapack(10_ENTIER, normal_mat)
    t10 = matmul(normal_mat, normal_rhs)
    deallocate(Jmat, rhsvec, normal_rhs, hess_eval, hess_beta)
  end subroutine third_of_geom_field_at

  ! Cell-to-cell CLS fourth derivative (Definition 1, step 2, k=3->4): builds a genuinely
  ! 4-exact fourth derivative directly at cell centers from a 3-exact third derivative
  ! (third_red, reduced 10-component) and third_of_geom_field_at(alpha,4,alpha) (the SAME
  ! geometric correction mechanism, one nesting level deeper). fourth_red is the reduced
  ! 15-component S^4(R^3) basis (xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,xyyz,xyzz,xzzz,yyyy,yyyz,
  ! yyzz,yzzz,zzzz). No caching: third_of_geom_field_at is recomputed from scratch for every
  ! (alpha,beta) pair -- expensive, see this section's header comment.
  subroutine compute_cell_cls_4exact_fourth(mesh, boundary_2d, third_red, fourth_red)
    use linear_solver_module, only: pseudo_inverse_inplace_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(10, mesh%n_elems), intent(in) :: third_red
    real(kind=DOUBLE), dimension(15, mesh%n_elems), intent(out) :: fourth_red

    integer(kind=ENTIER) :: alpha, i, beta, n_neigh, row, ell_t, ell_b, s, p, q, r
    real(kind=DOUBLE), dimension(:, :), allocatable :: Jmat, normal_rhs
    real(kind=DOUBLE), dimension(:), allocatable :: rhsvec
    real(kind=DOUBLE), dimension(:, :), allocatable :: third_of_z4_alpha, third_of_z4_beta
    real(kind=DOUBLE), dimension(3) :: h
    real(kind=DOUBLE), dimension(15, 15) :: normal_mat
    integer(kind=ENTIER), parameter :: red10_a(10) = (/1,1,1,1,1,1,2,2,2,3/)
    integer(kind=ENTIER), parameter :: red10_b(10) = (/1,1,1,2,2,3,2,2,3,3/)
    integer(kind=ENTIER), parameter :: red10_c(10) = (/1,2,3,2,3,3,2,3,3,3/)
    integer(kind=ENTIER), parameter :: red15_full81(15) = (/1,28,55,37,64,73,40,67,76,79,41,68,77,80,81/)

    call ensure_cell_moment_cache(mesh, 4_ENTIER)
    fourth_red = 0.0_DOUBLE

    do alpha = 1, mesh%n_elems
      n_neigh = mesh%elem(alpha)%n_neigh_by_vert
      if (n_neigh < 9) cycle
      allocate(third_of_z4_alpha(10, 81))
      call third_of_geom_field_at(mesh, boundary_2d, alpha, 4_ENTIER, alpha, third_of_z4_alpha)
      allocate(Jmat(10*n_neigh, 15), rhsvec(10*n_neigh), third_of_z4_beta(10, 81))
      Jmat = 0.0_DOUBLE
      row = 0
      do i = 1, n_neigh
        beta = mesh%elem(alpha)%neigh_by_vert(i)
        h = mesh%elem(beta)%coord - mesh%elem(alpha)%coord
        call third_of_geom_field_at(mesh, boundary_2d, beta, 4_ENTIER, beta, third_of_z4_beta)
        do ell_t = 1, 10
          p = red10_a(ell_t); q = red10_b(ell_t); r = red10_c(ell_t)
          do s = 1, 3
            ell_b = canon15_index(p, q, r, s)
            Jmat(row*10+ell_t, ell_b) = Jmat(row*10+ell_t, ell_b) + h(s)
          end do
          do ell_b = 1, 15
            Jmat(row*10+ell_t, ell_b) = Jmat(row*10+ell_t, ell_b) + (1.0_DOUBLE/24.0_DOUBLE) * &
              (third_of_z4_beta(ell_t, red15_full81(ell_b)) - third_of_z4_alpha(ell_t, red15_full81(ell_b)))
          end do
        end do
        rhsvec(row*10+1:row*10+10) = third_red(:, beta) - third_red(:, alpha)
        row = row + 1
      end do
      normal_mat = matmul(transpose(Jmat), Jmat)
      allocate(normal_rhs(15, 1))
      normal_rhs(:, 1) = matmul(transpose(Jmat), rhsvec)
      call pseudo_inverse_inplace_lapack(15_ENTIER, normal_mat)
      fourth_red(:, alpha) = matmul(normal_mat, normal_rhs(:, 1))
      deallocate(Jmat, rhsvec, normal_rhs, third_of_z4_beta, third_of_z4_alpha)
    end do
  end subroutine compute_cell_cls_4exact_fourth

end module arbitrary_high_order_module
