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
  public :: compute_derivative_hierarchy
  public :: compute_next_order_derivative_overlap
  public :: compute_derivative_hierarchy_overlap
  public :: timeline_event_type
  public :: compute_derivative_hierarchy_timed
  public :: compute_derivative_hierarchy_overlap_timed

  ! WENO-indicator regularization: just enough to avoid a literal division
  ! by zero when a candidate's tensor is exactly zero, no more -- per the
  ! user's request, replacing the earlier ad hoc 1e-6/1e-8 floors.
  real(kind=DOUBLE), parameter :: eps_weno = tiny(1.0_DOUBLE)
  real(kind=DOUBLE), parameter :: eps_weight_num = eps_weno
  ! Exponent on the oscillation indicator in the WENO weight,
  ! weight = 1/(eps+|dphi_v|^weno_power) -- 2026-09-13 experiment, per the
  ! user's request, to see whether a more aggressive de-centering (in the
  ! absence of any real positivity/slope limiter) reduces the small
  ! overshoots seen e.g. just behind a Sod rarefaction tail. Standard WENO
  ! practice is power=2; try higher powers here directly rather than
  ! building an actual limiter.
  integer(kind=ENTIER), parameter :: weno_power = 2

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

    nc_out = nc_in*d

    allocate(weno_num(nc_out, mesh%n_elems))
    allocate(weno_den(mesh%n_elems))
    allocate(dphi_v_cache(nc_out, mesh%n_vert))
    weno_num = 0.0_DOUBLE
    weno_den = 0.0_DOUBLE

    ! Single pass: each vertex's nodal derivative is computed once and
    ! immediately scattered, WENO-weighted, into every cell touching it --
    ! see accumulate_weno_contribution.
    do id_vert = 1, mesh%n_vert
      call accumulate_weno_contribution(mesh, d, nc_in, boundary_2d, id_vert, &
        phi, dphi_v_cache, weno_num, weno_den)
    end do

    do id_elem = 1, mesh%n_elems
      dphi(:, id_elem) = weno_num(:, id_elem) / weno_den(id_elem)
    end do

    deallocate(weno_num, weno_den, dphi_v_cache)
  end subroutine compute_next_order_derivative

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
    logical, dimension(:), allocatable :: is_send_cell, is_boundary_vertex

    nc_out = nc_in*d

    allocate(weno_num(nc_out, mesh%n_elems))
    allocate(weno_den(mesh%n_elems))
    allocate(dphi_v_cache(nc_out, mesh%n_vert))
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
    ! cells only (skip_cell filters at the cell level).
    do id_vert = 1, mesh%n_vert
      if (.not. is_boundary_vertex(id_vert)) cycle
      call accumulate_weno_contribution(mesh, d, nc_in, boundary_2d, id_vert, &
        phi, dphi_v_cache, weno_num, weno_den, skip_cell=(.not. is_send_cell))
    end do
    do id_elem = 1, mesh%n_elems
      if (.not. is_send_cell(id_elem)) cycle
      dphi(:, id_elem) = weno_num(:, id_elem) / weno_den(id_elem)
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
      call accumulate_weno_contribution(mesh, d, nc_in, boundary_2d, id_vert, &
        phi, dphi_v_cache, weno_num, weno_den, skip_cell=is_send_cell)
    end do
    do id_vert = 1, mesh%n_vert
      if (.not. is_boundary_vertex(id_vert)) cycle
      call accumulate_cached_weno_contribution(mesh, id_vert, dphi_v_cache, &
        weno_num, weno_den, skip_cell=is_send_cell)
    end do
    do id_elem = 1, mesh%n_elems
      if (is_send_cell(id_elem)) cycle
      dphi(:, id_elem) = weno_num(:, id_elem) / weno_den(id_elem)
    end do

    ! Step 4: this rank's ghost cells of dphi are only valid past this point.
    if (num_procs > 1) call mpi_memory_exchange_wait(mpi_send_recv, mesh%n_elems, nc_out, dphi)

    deallocate(weno_num, weno_den, dphi_v_cache)
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
        call accumulate_weno_contribution(mesh, d, nc_in, boundary_2d, id_vert, &
          phi_prev, dphi_v_cache, weno_num, weno_den, skip_cell=(.not. is_send_cell))
      end do
      do id_elem = 1, mesh%n_elems
        if (.not. is_send_cell(id_elem)) cycle
        dfield(order)%val(:, id_elem) = weno_num(:, id_elem) / weno_den(id_elem)
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
        call accumulate_weno_contribution(mesh, d, nc_in, boundary_2d, id_vert, &
          phi_prev, dphi_v_cache, weno_num, weno_den, skip_cell=is_send_cell)
      end do
      do id_vert = 1, mesh%n_vert
        if (.not. is_boundary_vertex(id_vert)) cycle
        call accumulate_cached_weno_contribution(mesh, id_vert, dphi_v_cache, &
          weno_num, weno_den, skip_cell=is_send_cell)
      end do
      do id_elem = 1, mesh%n_elems
        if (is_send_cell(id_elem)) cycle
        dfield(order)%val(:, id_elem) = weno_num(:, id_elem) / weno_den(id_elem)
      end do
      t1 = MPI_WTIME()
      n_events = n_events + 1
      events(n_events) = timeline_event_type(order, 'compute_int', t0, t1)

      t0 = MPI_WTIME()
      if (num_procs > 1) call mpi_memory_exchange_wait(mpi_send_recv, mesh%n_elems, nc_out, dfield(order)%val)
      t1 = MPI_WTIME()
      n_events = n_events + 1
      events(n_events) = timeline_event_type(order, 'wait', t0, t1)

      deallocate(weno_num, weno_den, dphi_v_cache)
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
  ! nonlinear-WENO-weighted contribution: weight = 1/(eps+|dphi_v|^2), the
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
  ! place to look.)
  !
  ! skip_cell, if present, excludes cells where skip_cell(id_elem) is
  ! .true. from accumulation -- used by compute_next_order_derivative_overlap
  ! to accumulate into only the cells that should receive this vertex's
  ! contribution at the point this is called (Step 1: only send cells;
  ! Step 3: only non-send cells).
  subroutine accumulate_weno_contribution(mesh, d, nc_in, boundary_2d, id_vert, &
      phi, dphi_v_cache, weno_num, weno_den, skip_cell)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in, id_vert
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(nc_in, mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(:, :), intent(inout) :: dphi_v_cache
    real(kind=DOUBLE), dimension(:, :), intent(inout) :: weno_num
    real(kind=DOUBLE), dimension(:), intent(inout) :: weno_den
    logical, dimension(:), intent(in), optional :: skip_cell

    real(kind=DOUBLE), dimension(nc_in*d) :: dphi_v

    call compute_nodal_derivative_at_vertex(mesh, d, nc_in, boundary_2d, &
      id_vert, phi, dphi_v)
    dphi_v_cache(:, id_vert) = dphi_v

    call scatter_weno_weighted(mesh, id_vert, dphi_v, weno_num, weno_den, skip_cell)
  end subroutine accumulate_weno_contribution

  ! Same scatter as accumulate_weno_contribution, but reusing an
  ! already-cached dphi_v instead of recomputing it -- used by
  ! compute_next_order_derivative_overlap's Step 3 to fold a boundary
  ! vertex's (Step-1-computed) contribution into the non-send cells it also
  ! touches, without a second LAPACK solve.
  subroutine accumulate_cached_weno_contribution(mesh, id_vert, dphi_v_cache, &
      weno_num, weno_den, skip_cell)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(:, :), intent(in) :: dphi_v_cache
    real(kind=DOUBLE), dimension(:, :), intent(inout) :: weno_num
    real(kind=DOUBLE), dimension(:), intent(inout) :: weno_den
    logical, dimension(:), intent(in), optional :: skip_cell

    call scatter_weno_weighted(mesh, id_vert, dphi_v_cache(:, id_vert), &
      weno_num, weno_den, skip_cell)
  end subroutine accumulate_cached_weno_contribution

  ! Shared scatter: weight = 1/(eps+max_c|dphi_v|^2), volume-weighted into
  ! every cell touching id_vert not excluded by skip_cell.
  subroutine scatter_weno_weighted(mesh, id_vert, dphi_v, weno_num, weno_den, skip_cell)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(:), intent(in) :: dphi_v
    real(kind=DOUBLE), dimension(:, :), intent(inout) :: weno_num
    real(kind=DOUBLE), dimension(:), intent(inout) :: weno_den
    logical, dimension(:), intent(in), optional :: skip_cell

    integer(kind=ENTIER) :: j, id_elem, id_sub_elem
    real(kind=DOUBLE) :: sub_elem_volume, vertex_weno_weight

    vertex_weno_weight = 1.0_DOUBLE / (eps_weight_num + maxval(abs(dphi_v))**weno_power)

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
  ! A vertex with fewer than n_min_neighbors immediate elem_neigh (mesh
  ! corners/edges) would otherwise give a singular/ill-conditioned system;
  ! gather_ls_neighbors expands the stencil by rings until there are enough,
  ! the same technique as ns_mesh_metric_module's (disabled) allocate_neigh.
  subroutine compute_nodal_derivative_at_vertex(mesh, d, nc_in, boundary_2d, &
      id_vert, phi, dphi_v)
    use linear_solver_module, only: lu_factor_lapack, lu_solve_mat_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: d, nc_in, id_vert
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(nc_in, mesh%n_elems), intent(in) :: phi
    real(kind=DOUBLE), dimension(nc_in*d), intent(out) :: dphi_v

    integer(kind=ENTIER), parameter :: max_basis = 4
    integer(kind=ENTIER), parameter :: n_min_neighbors = 2*max_basis
    real(kind=DOUBLE), parameter :: rel_spread_tol = 1.0e-8_DOUBLE
    integer(kind=ENTIER) :: n_basis, n_cand, n_active, j, a, b, i1, id_elem, n_neigh
    integer(kind=ENTIER), dimension(3) :: active_dim
    integer(kind=ENTIER), dimension(:), allocatable :: neigh
    integer(kind=ENTIER), dimension(max_basis) :: ipiv
    real(kind=DOUBLE), dimension(3) :: dx, dmin, dmax, spread
    real(kind=DOUBLE) :: weight, max_spread
    real(kind=DOUBLE), dimension(max_basis) :: basis
    real(kind=DOUBLE), dimension(max_basis, max_basis) :: mat
    real(kind=DOUBLE), dimension(max_basis, nc_in) :: rhs

    call ensure_neighbor_cache(mesh, n_min_neighbors)
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

    do j = 1, n_neigh
      id_elem = neigh(j)
      dx = mesh%elem(id_elem)%coord - mesh%vert(id_vert)%coord
      weight = 1.0_DOUBLE / max(dot_product(dx, dx), 1.0e-24_DOUBLE)

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
    if (n_active == 0) return ! no resolvable direction at all: report zero gradient

    call lu_factor_lapack(n_basis, mat(1:n_basis, 1:n_basis), ipiv(1:n_basis))
    call lu_solve_mat_lapack(n_basis, mat(1:n_basis, 1:n_basis), &
      ipiv(1:n_basis), nc_in, rhs(1:n_basis, :))

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

  ! Element neighbors of id_vert for the least-squares fit, ring-expanded
  ! (via each ring's own vertices) until there are at least n_min -- a
  ! near-boundary or corner vertex can otherwise have too few immediate
  ! elem_neigh to determine the fit, which silently produces a singular or
  ! near-singular system (and was observed to blow up the LU solve before
  ! this was added). Same technique as ns_mesh_metric_module's (disabled)
  ! allocate_neigh, reimplemented here to avoid a dependency from this
  ! core-level module onto the ns library.
  ! Builds neigh_cache_start/neigh_cache_list (CSR) by calling
  ! gather_ls_neighbors once per vertex -- the one-time cost that used to
  ! be paid on every one of compute_nodal_derivative_at_vertex's many
  ! calls per vertex (see the cache arrays' declaration). A no-op if the
  ! cache already matches this mesh's vertex count.
  subroutine ensure_neighbor_cache(mesh, n_min)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: n_min

    integer(kind=ENTIER) :: v, n_neigh, total
    integer(kind=ENTIER), dimension(:), allocatable :: neigh

    if (neigh_cache_n_vert == mesh%n_vert) return

    if (allocated(neigh_cache_start)) deallocate(neigh_cache_start)
    if (allocated(neigh_cache_list))  deallocate(neigh_cache_list)
    allocate(neigh_cache_start(mesh%n_vert + 1))

    neigh_cache_start(1) = 1
    do v = 1, mesh%n_vert
      call gather_ls_neighbors(mesh, v, n_min, n_neigh, neigh)
      neigh_cache_start(v+1) = neigh_cache_start(v) + n_neigh
      deallocate(neigh)
    end do

    total = neigh_cache_start(mesh%n_vert + 1) - 1
    allocate(neigh_cache_list(total))
    do v = 1, mesh%n_vert
      call gather_ls_neighbors(mesh, v, n_min, n_neigh, neigh)
      neigh_cache_list(neigh_cache_start(v):neigh_cache_start(v+1)-1) = neigh
      deallocate(neigh)
    end do

    neigh_cache_n_vert = mesh%n_vert
  end subroutine ensure_neighbor_cache

  subroutine gather_ls_neighbors(mesh, id_vert, n_min, n_neigh, neigh)
    use sort_module, only: add_sort_unique_int
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert, n_min
    integer(kind=ENTIER), intent(out) :: n_neigh
    integer(kind=ENTIER), dimension(:), allocatable, intent(out) :: neigh

    integer(kind=ENTIER) :: i, n_prev, n_pot_vert
    integer(kind=ENTIER), dimension(:), allocatable :: pot_vert

    n_neigh = 0
    call add_sort_unique_int(neigh, n_neigh, mesh%vert(id_vert)%elem_neigh)

    do while (n_neigh < min(n_min, mesh%n_elems))
      n_prev = n_neigh
      n_pot_vert = 0
      if (allocated(pot_vert)) deallocate(pot_vert)
      do i = 1, n_neigh
        call add_sort_unique_int(pot_vert, n_pot_vert, mesh%elem(neigh(i))%vert)
      end do
      do i = 1, n_pot_vert
        call add_sort_unique_int(neigh, n_neigh, mesh%vert(pot_vert(i))%elem_neigh)
      end do
      if (n_neigh == n_prev) exit ! stencil can't grow any further (tiny/disconnected mesh)
    end do
  end subroutine gather_ls_neighbors
end module arbitrary_high_order_module
