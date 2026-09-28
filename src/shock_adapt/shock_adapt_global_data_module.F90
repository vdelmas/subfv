! Standalone shock-adaptive node-movement solver: minimal, local mesh-node
! displacement near a detected shock front, X-Mesh-style (node-movement +
! fixed-topology mechanism from the "Mesh3D" crack-tracking paper, applied
! here to a first-order Euler shock instead of a crack surface), so that
! the shock ends up sitting on a mesh line/face instead of cutting through
! cell interiors. Fully self-contained: only depends on `core` (mesh,
! precision, mpi), no dependency on ns/lagrange/euler_ho/ale, so it can
! never affect (or be affected by) those solvers.
module shock_adapt_global_data_module
  use precision_module
  implicit none

  integer, parameter :: mnbc = 10
  real(kind=DOUBLE), parameter :: gamma = 1.4_DOUBLE

  !Mesh
  character(len=255) :: meshfile_path = "", meshfile
  logical :: boundary_2d = .FALSE.

  !Boundary conditions: 'wall' (default, mirror normal velocity), 'freestream'
  !(fixed bc_val state), 'outflowsupersonic'/'outflow' (zero-gradient ghost =
  !interior state, standard supersonic-outflow treatment -- same convention
  !as euler_ho_module.F90's BC_OUTFLOW / ns's BC_EULER_OUTFLOWSUPERSONIC).
  integer(kind=ENTIER), parameter :: BC_WALL = 0, BC_FREESTREAM = 1, BC_OUTFLOW = 2
  integer(kind=ENTIER) :: n_bc = 0
  character(len=255), dimension(mnbc) :: bc_name = ""
  character(len=255), dimension(mnbc) :: bc_type = ""
  real(kind=DOUBLE), dimension(5, mnbc) :: bc_val = 0.0_DOUBLE
  integer(kind=ENTIER), dimension(mnbc) :: bc_kind = BC_WALL

  !Flux scheme for the 2-state face Riemann solver (all three are plain
  !two-state fluxes, ported from ns_euler_rs_module.F90 -- no nodal
  !lambda-iteration solve needed for a first-order, fixed-mesh-between-moves
  !scheme).
  character(len=32) :: flux_scheme = 'modified_three_wave'
  integer(kind=ENTIER), parameter :: FLUX_THREE_WAVE = 0, FLUX_MODIFIED_THREE_WAVE = 1, FLUX_TWO_WAVE = 2
  integer(kind=ENTIER) :: flux_scheme_id = FLUX_MODIFIED_THREE_WAVE

  !Time
  real(kind=DOUBLE) :: cfl = 0.5_DOUBLE

  ! Reconstruction order. 1 = plain cell-average face states (the original
  ! scheme). 2 = linear (MUSCL-type) extrapolation to each face using a
  ! Green-Gauss gradient from src/arbitrary_high_order/ (the same "aho_gg"
  ! reconstruction euler_ho_module.F90 uses, reused here as a library --
  ! only links core, so this stays consistent with shock_adapt's own
  ! independence). See shock_adapt_module.F90's compute_rhs.
  integer(kind=ENTIER) :: order = 1
  ! Run the first order_ramp_iter Stage-1 iterations at order 1 regardless
  ! of `order`, to get past the sharp startup transient from a uniform
  ! initial condition before switching to the requested order -- see
  ! shock_adapt_module.F90's compute_rhs header comment.
  integer(kind=ENTIER) :: order_ramp_iter = 0

  !Init: uniform state only (this prototype targets the cylinder-tunnel case,
  !whose initial condition is a uniform freestream).
  logical :: init_uniform = .FALSE.
  real(kind=DOUBLE), dimension(5) :: sol_uniform = 0.0_DOUBLE

  !--- Shock-adaptive node movement (Stage 2) ---
  ! Face is flagged "shocked" if the normalized pressure jump across it
  ! exceeds this threshold AND the face is compressive (converging normal
  ! velocity) -- a standard shock sensor, gated to avoid flagging contacts/
  ! expansions.
  real(kind=DOUBLE) :: shock_sensor_threshold = 0.15_DOUBLE
  ! Optional extra AND-condition on the normalized density jump
  ! (|rhoR-rhoL|/(rhoR+rhoL)) -- pressure+compression already excludes
  ! contact discontinuities (no pressure jump there), so this only makes
  ! the detector stricter, not a replacement criterion.
  logical :: use_density_sensor = .FALSE.
  real(kind=DOUBLE) :: density_sensor_threshold = 0.15_DOUBLE
  ! Density-GRADIENT-magnitude nodal detector (what a numerical schlieren
  ! actually shows), normalized dimensionless by local cell size/rho_bar --
  ! see compute_shock_sensor_grad in shock_adapt_module.F90. This is the
  ! detector actually used by use_curvature_move/detect_only (supersedes
  ! the earlier pressure+velocity-jump compute_shock_sensor_nodal, kept in
  ! the code but unused by default).
  real(kind=DOUBLE) :: grad_sensor_threshold = 0.5_DOUBLE
  ! Cap on a vertex's displacement, as a fraction of its local minimum
  ! distance to a neighboring cell centroid -- keeps the movement "minimal"
  ! (bounded, local) as opposed to a full mesh relocation.
  real(kind=DOUBLE) :: max_move_frac = 0.3_DOUBLE
  ! Controllable pull toward the downstream (compressed) side, as a
  ! fraction of the flagged face's own inter-centroid distance -- this is
  ! the actual lever for how visible the node movement is (max_move_frac
  ! above is only a safety cap, see shock_adapt_module.F90's
  ! compute_node_displacement for why a pressure-interpolated target
  ! self-defeatingly canceled out regardless of the cap).
  real(kind=DOUBLE) :: shock_snap_frac = 0.4_DOUBLE

  !--- Stage iteration budgets ---
  ! Stage 1: march to a (quasi-)steady baseline on the unmodified mesh.
  integer(kind=ENTIER) :: n_iter_steady = 20000
  ! Single-shot mode (use_iterative_moves=.FALSE.): reconvergence budget
  ! after the one Stage-2 node movement.
  integer(kind=ENTIER) :: n_iter_post_move = 5000

  !--- Iterative mode (use_iterative_moves=.TRUE.) ---
  ! Repeat move->quick-restabilize n_cycles times instead of moving once:
  ! after the Stage-1 baseline, each cycle re-detects the shock, applies a
  ! bounded move, then runs n_iter_restab iterations (much shorter than
  ! n_iter_post_move -- only enough to let the *local* perturbation from
  ! that one move settle, not a full reconvergence) before the next cycle's
  ! detection. A VTU snapshot is written at cycle 1, every 30th cycle, and
  ! the final cycle, so the cumulative mesh drift can be inspected.
  logical :: use_iterative_moves = .FALSE.
  integer(kind=ENTIER) :: n_cycles = 100
  integer(kind=ENTIER) :: n_iter_restab = 300

  !--- Curvature-based move (new, per-node) ---
  ! Detect shock nodes directly (nodal min/max pressure(+density) jump
  ! among a vertex's neighboring cells, reusing shock_sensor_threshold/
  ! density_sensor_threshold), then move each flagged node along the local
  ! normal of a quadratic curve fit through its flagged topological
  ! neighbors, to reduce the front's local curvature -- see
  ! shock_adapt_module.F90's compute_node_displacement_curvature. Single
  ! pass only for now (no repeated geometric smoothing, no CFD-in-the-loop).
  logical :: use_curvature_move = .FALSE.
  ! Relaxation on the curvature-based correction: dn = curvature_relax *
  ! (h_pred(s_i) - h_i), instead of moving the node all the way to fully
  ! flatten the local fit in one step. 1.0 = full correction (the original
  ! behavior); user-confirmed 0.5 is a good default -- under-relaxing
  ! avoids overshoot, standard practice for this kind of curvature-flow
  ! smoothing.
  real(kind=DOUBLE) :: curvature_relax = 0.5_DOUBLE
  ! If .TRUE., stop right after Stage 1 and write the nodal sensor/flagged
  ! diagnostic fields, without moving anything -- for inspecting which
  ! nodes get flagged before trusting the move step.
  logical :: detect_only = .FALSE.

  !Output
  integer(kind=ENTIER) :: n_iter_print = 200

contains
  subroutine read_input_parameters(filename)
    implicit none

    character(len=*), intent(in) :: filename

    integer(kind=ENTIER) :: funit

    namelist /INPUT_PARAM/ &
      meshfile_path, meshfile, boundary_2d, &
      n_bc, bc_name, bc_type, bc_val, &
      flux_scheme, cfl, order, order_ramp_iter, &
      init_uniform, sol_uniform, &
      shock_sensor_threshold, max_move_frac, shock_snap_frac, &
      use_density_sensor, density_sensor_threshold, grad_sensor_threshold, &
      n_iter_steady, n_iter_post_move, n_iter_print, &
      use_iterative_moves, n_cycles, n_iter_restab, &
      use_curvature_move, curvature_relax, detect_only

    open(newunit=funit, file=trim(adjustl(filename)))
    read(unit=funit, nml=INPUT_PARAM)
    close(funit)

    select case (trim(adjustl(flux_scheme)))
    case ('three_wave')
      flux_scheme_id = FLUX_THREE_WAVE
    case ('modified_three_wave')
      flux_scheme_id = FLUX_MODIFIED_THREE_WAVE
    case ('two_wave')
      flux_scheme_id = FLUX_TWO_WAVE
    case default
      print*, "ERROR: unknown flux_scheme '", trim(flux_scheme), &
        "' (only 'three_wave'/'modified_three_wave'/'two_wave' supported)"
      error stop
    end select
  end subroutine read_input_parameters

  subroutine init_bc_kind()
    implicit none

    integer(kind=ENTIER) :: i
    character(len=255) :: t

    do i = 1, n_bc
      t = trim(adjustl(bc_type(i)))
      if (t == "wall" .or. t == "") then
        bc_kind(i) = BC_WALL
      else if (t == "freestream") then
        bc_kind(i) = BC_FREESTREAM
      else if (t == "outflowsupersonic" .or. t == "outflow") then
        bc_kind(i) = BC_OUTFLOW
      else
        print*, "ERROR: bc_type unrecognized for BC ", i, ": '", trim(t), &
          "' (only 'wall'/'freestream'/'outflowsupersonic' supported)"
        error stop
      end if
    end do
  end subroutine init_bc_kind
end module shock_adapt_global_data_module
