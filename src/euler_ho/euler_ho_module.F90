! High-order explicit Euler solver (finite volume, SSP-RK3 in time).
!
! Conservative variables: sol(1:5,i) = [rho, rho*u, rho*v, rho*w, rho*E]
! Primitive  variables:   prim(1:5,i) = [rho, u, v, w, p]
!
! Spatial reconstruction up to order 4 (aho reconstruction source only --
! ls_reconstruction still tops out at order 3, see its own header) using
! cell-centred least-squares polynomial fit over vertex-neighbours, with
! face quadrature.
! boundary_2d=T → 2D polynomial (x,y only); boundary_2d=F → full 3D.
! Numerical flux: Rusanov (local Lax-Friedrichs).
module euler_ho_module
  use precision_module
  use mesh_module
  use quadrature_module
  use arbitrary_high_order_module, only: compute_next_order_derivative, &
    compute_next_order_derivative_cweno, apply_grad_bias_correction, &
    aho_module_use_weno_blend => use_weno_blend, &
    aho_module_use_green_gauss => use_green_gauss, &
    aho_module_use_cweno_center => use_cweno_center, &
    aho_module_cweno_center_weight => cweno_center_weight, &
    aho_module_cweno_center_power => cweno_center_power, &
    aho_module_eps_weight_num => eps_weight_num, &
    aho_module_eps_weight_num_deep => eps_weight_num_deep, &
    aho_module_eps_weight_num_gg => eps_weight_num_gg, &
    aho_module_eps_weight_num_deep_gg => eps_weight_num_deep_gg, &
    aho_module_weno_power => weno_power, &
    aho_module_grad_norm_derate => grad_norm_derate, &
    aho_module_grad_norm_derate_gg => grad_norm_derate_gg, &
    aho_module_use_alt_gg_weight => use_alt_gg_weight
  use mpi_module, only: mpi_send_recv_type, mpi_memory_exchange
  implicit none
  private

  real(kind=DOUBLE), parameter :: PI = 4.0_DOUBLE * datan(1.0_DOUBLE)

  ! ----------------------------------------------------------------
  ! Global parameters (set by read_params)
  ! ----------------------------------------------------------------
  integer(kind=ENTIER), parameter :: MAX_BC = 10
  character(len=255), public :: meshfile_path = ''
  character(len=255), public :: meshfile      = ''
  real(kind=DOUBLE),  public :: gamma_gas     = 1.4_DOUBLE
  ! Per-cell ratio of specific heats, for multi-material cases (e.g. a
  ! helium bubble in air). Allocated and default-filled with gamma_gas by
  ! init_sol; specific init cases (e.g. shock_bubble, init=6) overwrite it
  ! per-cell from a geometric material criterion. Every physics function
  ! that used to read the module-global scalar gamma_gas directly now
  ! takes the relevant cell(s)' gamma_arr value(s) as an explicit argument
  ! instead, mirroring src/lagrange/lagrange_module.F90's gamma_arr/gl/gr
  ! pattern (gamma_gas itself is kept only as the uniform-fill default and
  ! for the single-material analytic vortex state, vortex_prim).
  real(kind=DOUBLE), dimension(:), allocatable, public :: gamma_arr
  ! Multi-material gamma TRANSPORT (Abgrall 1996 quasi-conservative
  ! approach): rho*Gamma, Gamma = 1/(gamma-1), is carried as sol's 6th
  ! component (see sol's declaration in init_sol for the full rationale),
  ! evolved by euler_ho_main's ordinary RK3/MPI-exchange machinery
  ! together with the 5 physical equations -- no separate bookkeeping.
  ! Advecting Gamma itself (rather than gamma directly) avoids spurious
  ! pressure oscillations at a material interface, a well-known pitfall of
  ! the naive approach. gamma_arr is re-derived from sol(6,:)/sol(1,:) via
  ! sync_gamma_arr, called once per RK stage with that stage's own sol.
  ! For a single-material run, Gamma is uniform everywhere so this
  ! reduces to a no-op and gamma_arr stays exactly gamma_gas, as before --
  ! zero behaviour change.
  real(kind=DOUBLE),  public :: cfl           = 0.4_DOUBLE
  real(kind=DOUBLE),  public :: tmax          = 1.0_DOUBLE
  integer(kind=ENTIER), public :: order       = 1     ! 1, 2, 3, or 4 (aho only; see ls_reconstruction)
  logical, public :: boundary_2d              = .false.
  integer(kind=ENTIER), public :: n_sol_vtu  = 10
  logical, public :: compute_error            = .false.
  logical, public :: error_2d                 = .true.  ! vortex in xy plane
  ! At order>=2, use arbitrary_high_order_module's least-squares dual/primal
  ! hierarchy (compute_next_order_derivative) for grad/hess instead of this
  ! module's own ls_reconstruction -- same face_flux_loop/reconstruct/RK3
  ! downstream, so this isolates the reconstruction source for a like-for-like
  ! comparison (e.g. on a shock case).
  logical, public :: use_aho_reconstruction   = .false.
  ! At use_aho_reconstruction=.true. only: .false. disables the nonlinear
  ! WENO blend (arbitrary_high_order_module's own use_weno_blend flag,
  ! set from this one in aho_reconstruction), leaving a plain
  ! volume-weighted linear combination of nodal derivatives -- a baseline
  ! comparison point analogous to ls_reconstruction's raw, unlimited fit,
  ! but built from aho's nodal-derivative machinery instead of a
  ! cell-centered polynomial fit.
  logical, public :: use_weno_blend           = .true.
  ! Selects the per-vertex nodal fit compute_next_order_derivative uses
  ! (arbitrary_high_order_module's own use_green_gauss flag, set from this
  ! one in aho_reconstruction, exactly like use_weno_blend): .false.
  ! (default) is the weighted least-squares fit ("aho ls", the original,
  ! fully-validated path); .true. is the Green-Gauss divergence-theorem
  ! fit ("aho gg", 2026-09-17), roughly 4x cheaper per call once its own
  ! geometry-only fit matrix is cached (see that module's ensure_green_
  ! gauss_mat_cache), reaches the same order 2/3 as aho ls on the vortex,
  ! and reuses apply_grad_bias_correction unchanged for order 4 -- but its
  ! own oscillation indicator is a gradient-norm term only (no residual
  ! counterpart), not yet validated on non-axis-aligned discontinuities
  ! the way aho ls's combined indicator has been.
  logical, public :: use_green_gauss          = .false.
  ! CWENO-style central candidate (2026-09-15): at use_weno_blend=.true.
  ! only, folds a large-weight plain/linear candidate into the SAME
  ! nonlinear blend at EVERY recursion level (grad, hess, third). Verified
  ! on a synthetic cubic field to let hess/third converge with a real
  ! order instead of being stuck at O(h^0) under plain WENO -- but
  ! reverted to .false. as the production default on 2026-09-16: on real
  ! shock benchmarks (DMR, cylinder, FFS) it made positivity/overshoot
  ! measurably WORSE than plain WENO at orders 2 and 3 too (not just 4),
  ! because the large always-on linear-candidate weight dilutes classical
  ! WENO's de-centering near genuine discontinuities everywhere, not just
  ! where CWENO's own OI_c indicator judges it safe to do so -- and it
  ! still never reached true order-4 accuracy (capped at order 3's rate
  ! by the grad step's own bias, untouched by this fix). Kept available
  ! behind this flag for reference/experimentation; the production path
  ! is plain compute_next_order_derivative at every level. See
  ! arbitrary-high-order memory for the full writeup.
  logical, public :: use_cweno_center          = .false.
  real(kind=DOUBLE), public :: cweno_center_weight = 1000.0_DOUBLE
  integer(kind=ENTIER), public :: cweno_center_power = 4
  ! At order=4 only: corrects the grad step's own O(h^2) affine-fit bias
  ! in place, using the per-vertex nodal Hessian/third-derivative tensors
  ! already computed as intermediates while building hess/third (see
  ! arbitrary_high_order_module's apply_grad_bias_correction for the full
  ! derivation) -- no extra neighbors, no higher-degree fit, same 4-point
  ! stencil throughout. Verified on the synthetic cubic field: reduces
  ! grad's error by ~4-6x under the current WENO scheme (would reach
  ! machine precision under a plain linear blend, where hess/third are
  ! themselves already exact); the residual gap under WENO tracks how far
  ! hess/third themselves still are from exact (see eps_weight_num_deep).
  ! Only implemented for boundary_2d=.true.; a .false. mesh silently gets
  ! no correction (apply_grad_bias_correction is a documented no-op
  ! there), not yet extended to genuine 3D.
  logical, public :: use_grad_bias_correction  = .true.
  ! Classical (non-SSP) 4-stage RK4 in place of SSP-RK3 (2026-09-18, per
  ! the user, to rule out RK3's own O(dt^3) temporal error as the reason
  ! order 4's spatial accuracy doesn't clearly separate from order 3 on
  ! some benchmarks -- a direct, unambiguous check to run alongside the
  ! cheaper dt-shrinking diagnostic already tried for the same question).
  ! .false. (default) keeps the existing SSP-RK3 exactly as before; see
  ! euler_ho_main.F90's time loop for both branches, selected once at
  ! the top of the timestep, everything else (compute_rhs, MPI exchange,
  ! sync_gamma_arr, output/error cadence) shared identically either way.
  logical, public :: use_rk4 = .false.
  ! CFL time-step criterion: .false. (default) keeps the existing dt =
  ! cfl*volume/sum_lambda(i), summing each face's A_pcf*lambda_pcf over
  ! every face of the cell (a domain-of-dependence bound on the cell's
  ! full boundary). .true. instead uses dt = cfl*volume/max_lambda(i),
  ! the MAX over the cell's faces of that same per-face quantity, a much
  ! less restrictive criterion for cells with many faces (2026-09-18, per
  ! the user, to test alongside cfl=0.95 locally whether it still holds
  ! stable at every order on a coarse shock-bubble mesh). face_flux_loop
  ! accumulates a face-local total across that face's quadrature points
  ! first, then folds it into sum_lambda(il)/sum_lambda(ir) via + (sum
  ! mode) or max() (max mode) -- same array either way, so compute_dt and
  ! every call site are untouched; sum mode is bit-for-bit identical to
  ! the pre-existing behaviour.
  logical, public :: use_max_lambda_dt = .false.
  ! Runtime-settable mirrors of arbitrary_high_order_module's own
  ! eps_weight_num/eps_weight_num_deep/weno_power/grad_norm_derate
  ! (2026-09-17, made settable there too -- see that module's own
  ! comments for what each controls, including why eps_weight_num's
  ! default was raised 1e-6 -> 1e-2 the same day). Defaults match the
  ! module's defaults exactly; exposed here, and propagated into the aho
  ! module inside aho_reconstruction exactly like use_weno_blend/
  ! use_cweno_center already are, so an input_data.f can sweep OI/epsilon
  ! calibration without a rebuild.
  real(kind=DOUBLE), public :: eps_weight_num      = 1.0e-2_DOUBLE
  real(kind=DOUBLE), public :: eps_weight_num_deep = 1.0_DOUBLE
  real(kind=DOUBLE), public :: eps_weight_num_gg = 1.0e-2_DOUBLE
  real(kind=DOUBLE), public :: eps_weight_num_deep_gg = 1.0_DOUBLE
  integer(kind=ENTIER), public :: weno_power       = 1
  real(kind=DOUBLE), public :: grad_norm_derate    = 1.0e4_DOUBLE
  real(kind=DOUBLE), public :: grad_norm_derate_gg = 1.0e4_DOUBLE
  logical, public :: use_alt_gg_weight = .false.
  ! Numerical flux at each face quadrature point: 'rusanov' (local
  ! Lax-Friedrichs, the original default) or 'three_wave' (an HLLC-family,
  ! 3-wave approximate Riemann solver -- same algorithm as
  ! ns_euler_rs_module's three_wave, ported here rather than linked, to
  ! keep this solver's dependency footprint self-contained). three_wave
  ! resolves the contact wave much more sharply than Rusanov, which is
  ! the wave any reconstruction artifact near a contact discontinuity
  ! would otherwise be partly masked by.
  character(len=32), public :: flux_scheme    = 'rusanov'
  ! flux_scheme resolved to an integer once in read_params, instead of a
  ! trim+adjustl+string-compare per quadrature point inside face_flux_loop
  ! (profiling found this string comparison alone costing ~9.5% of total
  ! solver runtime at order 2 -- flux_scheme never changes after startup,
  ! so comparing it as a string on every one of the (many quad points) x
  ! (faces) x (RK stages) x (iterations) calls was pure waste).
  integer(kind=ENTIER), parameter :: FLUX_RUSANOV = 0, FLUX_TWO_WAVE = 1, FLUX_THREE_WAVE = 2
  integer(kind=ENTIER) :: flux_scheme_id = FLUX_RUSANOV

  ! Boundary conditions (n_bc <= MAX_BC)
  integer(kind=ENTIER), public :: n_bc = 0
  character(len=255), dimension(MAX_BC), public :: bc_name = ''
  character(len=255), dimension(MAX_BC), public :: bc_type = ''
  real(kind=DOUBLE), dimension(5, MAX_BC), public :: bc_val  = 0.0_DOUBLE
  ! bc_type resolved to integers once in read_params, for the same reason
  ! as flux_scheme_id: ghost_prim's string select-case re-parsed
  ! bc_type(id_bc) on every boundary quadrature-point call.
  integer(kind=ENTIER), parameter :: BC_WALL = 0, BC_FREESTREAM = 1, &
    BC_OUTFLOW = 2, BC_DMR_TOP = 3
  integer(kind=ENTIER), dimension(MAX_BC) :: bc_type_id = BC_WALL

  ! Init
  !   0 = uniform (sol_uniform: primitive w)
  !   1 = Sod 1D (x1drp, sol_w_1drp_l, sol_w_1drp_r: primitive)
  !   2 = isentropic vortex (beta=5, centre at origin at t=0, moves with (u_bg,v_bg))
  !   6 = shock-bubble interaction (Haas & Sturtevant 1987 JFM 181; see
  !       init_sol's case(6) for the full geometry/state -- multi-material,
  !       sets gamma_arr per-cell)
  real(kind=DOUBLE), public :: u_bg_vortex  = 0.0_DOUBLE  ! background x-velocity
  real(kind=DOUBLE), public :: v_bg_vortex  = 0.0_DOUBLE  ! background y-velocity
  integer(kind=ENTIER), public :: init         = 0
  real(kind=DOUBLE), dimension(5), public :: sol_uniform    = [1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE]
  real(kind=DOUBLE), public :: x1drp          = 0.5_DOUBLE
  real(kind=DOUBLE), dimension(5), public :: sol_w_1drp_l   = [1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE]
  real(kind=DOUBLE), dimension(5), public :: sol_w_1drp_r   = [0.125_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.1_DOUBLE]

  public :: read_params
  public :: init_sol
  public :: compute_prim
  public :: compute_rhs
  public :: compute_dt
  public :: sync_gamma_arr
  public :: compute_error_vortex
  public :: test_reconstruction_exactness
  public :: compute_cell_moments
  ! Exposed for the standalone quadrature+reconstruction consistency test
  ! (scratchpad/quad_recon_test) -- not otherwise called outside this module.
  public :: reconstruct, ls_reconstruction, aho_reconstruction

  ! Cache of each cell's own second geometric moment about its centroid,
  ! M_jk(i) = (1/V_i) * int_cell (x_j - xc_j)(x_k - xc_k) dV -- a pure
  ! mesh-geometry quantity, independent of order/field, computed once by
  ! compute_cell_moments (called from euler_ho_main right after
  ! compute_geometry_mesh) and read by reconstruct() at order 3. See the
  ! header comment above reconstruct() for why it is needed.
  real(kind=DOUBLE), dimension(:, :, :), allocatable :: cell_moment
  ! Same idea, one degree up: the third geometric moment
  ! M_jkl(i) = (1/V_i) * int_cell (x_j-xc_j)(x_k-xc_k)(x_l-xc_l) dV, needed
  ! by reconstruct()'s order-4 (cubic) term for the same reason -- only
  ! allocated/filled when order>=4 (compute_cell_moments needs degree-3
  ! quadrature to get this one exactly, vs degree-2 for cell_moment alone).
  real(kind=DOUBLE), dimension(:, :, :, :), allocatable :: cell_moment3

contains

  ! ----------------------------------------------------------------
  ! Namelist I/O
  ! ----------------------------------------------------------------
  subroutine read_params(filename)
    character(len=*), intent(in) :: filename
    integer(kind=ENTIER) :: funit, i_bc
    namelist /INPUT_PARAM/ &
      meshfile_path, meshfile, &
      n_bc, bc_name, bc_type, bc_val, &
      boundary_2d, &
      init, sol_uniform, x1drp, sol_w_1drp_l, sol_w_1drp_r, &
      u_bg_vortex, v_bg_vortex, &
      gamma_gas, order, cfl, tmax, n_sol_vtu, &
      compute_error, error_2d, use_aho_reconstruction, use_weno_blend, use_green_gauss, &
      use_cweno_center, cweno_center_weight, cweno_center_power, use_grad_bias_correction, use_rk4, &
      use_max_lambda_dt, flux_scheme, &
      eps_weight_num, eps_weight_num_deep, weno_power, grad_norm_derate, &
      grad_norm_derate_gg, use_alt_gg_weight, eps_weight_num_gg, eps_weight_num_deep_gg
    open(newunit=funit, file=trim(adjustl(filename)))
    read(nml=INPUT_PARAM, unit=funit)
    close(funit)

    ! Resolve flux_scheme to an integer once here, matching subfvns's own
    ! scheme_id convention (ns_global_data_module) -- see flux_scheme_id's
    ! declaration for why. Same fallback as the old string select case:
    ! anything not recognized defaults to Rusanov.
    select case (trim(adjustl(flux_scheme)))
    case ('three_wave')
      flux_scheme_id = FLUX_THREE_WAVE
    case ('two_wave')
      flux_scheme_id = FLUX_TWO_WAVE
    case default
      flux_scheme_id = FLUX_RUSANOV
    end select

    ! Same treatment for bc_type -- see bc_type_id's declaration. Same
    ! fallback as ghost_prim's old string select case: anything
    ! unrecognized (including blank) defaults to a slip wall.
    do i_bc = 1, n_bc
      select case (trim(adjustl(bc_type(i_bc))))
      case ('freestream')
        bc_type_id(i_bc) = BC_FREESTREAM
      case ('outflowsupersonic', 'outflow')
        bc_type_id(i_bc) = BC_OUTFLOW
      case ('dmr_top')
        bc_type_id(i_bc) = BC_DMR_TOP
      case default
        bc_type_id(i_bc) = BC_WALL
      end select
    end do
  end subroutine read_params

  ! ----------------------------------------------------------------
  ! Solution initialisation
  ! ----------------------------------------------------------------
  subroutine init_sol(mesh, sol)
    type(mesh_type), intent(in)    :: mesh
    ! sol(1:5,:) = [rho, rho*u, rho*v, rho*w, rho*E] as before; sol(6,:) =
    ! rho*Gamma, Gamma = 1/(gamma-1), the multi-material gamma-transport
    ! variable (Abgrall 1996 quasi-conservative approach -- advecting
    ! Gamma itself, rather than gamma directly, avoids spurious pressure
    ! oscillations at a material interface). Folded into sol (rather than
    ! kept as a separate parallel array) so it rides through euler_ho_main's
    ! existing RK3/MPI-exchange machinery for free, with no separate
    ! bookkeeping to keep in sync -- a real bug from an earlier, separate-
    ! array version (module-level state going stale mid-RK-step, found by
    ! subfv-c1, 2026-09-16) is structurally impossible this way. Every
    ! function that only knows the physical 5-equation system (conserv_to_
    ! primit, reconstruct, three_wave/two_wave/rusanov, cs, euler_flux)
    ! takes an explicit sol(1:5,i)/wL(1:5)/wR(1:5) slice; gamma_arr is
    ! re-derived from sol(6,i) via sync_gamma_arr after every RK stage
    ! (see euler_ho_main.F90). For a single-material run, Gamma is uniform
    ! everywhere so this is a no-op and gamma_arr stays exactly gamma_gas,
    ! as before -- zero behaviour change.
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(out) :: sol

    integer(kind=ENTIER) :: i
    real(kind=DOUBLE), dimension(5) :: w
    real(kind=DOUBLE), dimension(3) :: xc
    real(kind=DOUBLE) :: rb

    if (.not. allocated(gamma_arr)) allocate(gamma_arr(mesh%n_elems))
    gamma_arr = gamma_gas

    do i = 1, mesh%n_elems
      xc = mesh%elem(i)%coord
      select case (init)
      case (0)   ! uniform
        w = sol_uniform
      case (1)   ! Sod 1D (x-direction)
        if (xc(1) < x1drp) then
          w = sol_w_1drp_l
        else
          w = sol_w_1drp_r
        end if
      case (3)   ! double Mach reflection: shock at t=0, see dmr_state
        w = dmr_state(xc(1), xc(2), 0.0_DOUBLE)
      case (4)   ! Shu-Osher shock/entropy-wave interaction, see shu_osher_state
        w = shu_osher_state(xc(1))
      case (5)   ! Woodward-Colella interacting blast waves, see woodward_colella_state
        w = woodward_colella_state(xc(1))
      case (2)   ! isentropic vortex: high-order cell average via volume quadrature
        call vortex_cell_average(mesh, i, 0.0_DOUBLE, w)
      case (6)   ! Shock-bubble interaction (Haas & Sturtevant 1987, JFM 181;
                 ! also Quirk & Karni 1996 JFM 318, Razmi et al. 2019 JAFM
                 ! 12(2) 631-645). Ms=1.22 planar shock in air hits a
                 ! cylindrical He+28%air bubble. Full domain 500x89mm (not
                 ! symmetry-reduced): shock initially at x=400mm, bubble
                 ! centred at (350mm,44.5mm), radius 25mm. The post-shock
                 ! block (x>400mm) set here is ALSO imposed continuously at
                 ! the right boundary via a 'freestream' BC (see
                 ! input_data.f / bc_val) since the shock propagates left,
                 ! away from that edge, over the run -- the IC and the BC
                 ! must both carry the same post-shock state for
                 ! consistency at t=0+.
        block
          real(kind=DOUBLE), parameter :: xc_bub    = 0.350_DOUBLE
          real(kind=DOUBLE), parameter :: yc_bub    = 0.0445_DOUBLE
          real(kind=DOUBLE), parameter :: r_bub     = 0.025_DOUBLE
          real(kind=DOUBLE), parameter :: x_shock   = 0.400_DOUBLE
          real(kind=DOUBLE), parameter :: gamma_air = 1.4_DOUBLE
          ! He contaminated 28% by mass with air (Haas & Sturtevant):
          ! effective gamma and density ratio vs. ambient air.
          real(kind=DOUBLE), parameter :: gamma_bub = 1.648_DOUBLE
          real(kind=DOUBLE), parameter :: rho_bub   = 0.287_DOUBLE / 1.578_DOUBLE
          ! Post-shock air state from the normal-shock relations for
          ! gamma=1.4, Ms=1.22 (rho2/rho1, p2/p1, and u2 via the standard
          ! shock-jump formulas); shock moves in -x, so u2 < 0.
          real(kind=DOUBLE), parameter :: rho_post  = 1.376364_DOUBLE
          real(kind=DOUBLE), parameter :: p_post    = 1.569800_DOUBLE
          real(kind=DOUBLE), parameter :: vx_post   = -0.394729_DOUBLE
          rb = sqrt((xc(1) - xc_bub)**2 + (xc(2) - yc_bub)**2)
          if (xc(1) > x_shock) then
            gamma_arr(i) = gamma_air
            w = [rho_post, vx_post, 0.0_DOUBLE, 0.0_DOUBLE, p_post]
          else if (rb <= r_bub) then
            gamma_arr(i) = gamma_bub
            w = [rho_bub, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE]
          else
            gamma_arr(i) = gamma_air
            w = [1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE]
          end if
        end block
      case default
        w = sol_uniform
      end select
      sol(1:5, i) = primit_to_conserv(w, gamma_arr(i))
      sol(6, i)   = sol(1, i) / (gamma_arr(i) - 1.0_DOUBLE)
    end do
  end subroutine init_sol

  ! ----------------------------------------------------------------
  ! Convert all cells from conservative to primitive
  ! ----------------------------------------------------------------
  subroutine compute_prim(mesh, sol, prim)
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(in)  :: sol
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(out) :: prim
    integer(kind=ENTIER) :: i
    do i = 1, mesh%n_elems
      prim(:, i) = conserv_to_primit(sol(1:5, i), gamma_arr(i))
    end do
  end subroutine compute_prim

  ! ----------------------------------------------------------------
  ! Recompute gamma_arr(i) = 1 + sol(1,i)/sol(6,i) (rho/(rho*Gamma), i.e.
  ! 1+1/Gamma) from the given RK-stage's own sol -- both the density and
  ! the just-transported Gamma field come from the SAME stage array (e.g.
  ! sol1), which is what keeps this correctly synchronised stage-by-stage
  ! (see sol's declaration in init_sol for why this used to be a separate,
  ! easy-to-desync array). Call after every RK-stage state update, before
  ! the next compute_prim/compute_rhs call reads gamma_arr for its flux
  ! gL/gR -- mirrors how sol's own ghost cells get refreshed via
  ! mpi_memory_exchange after each stage (see euler_ho_main.F90).
  ! ----------------------------------------------------------------
  subroutine sync_gamma_arr(mesh, sol)
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(in) :: sol
    integer(kind=ENTIER) :: i
    do i = 1, mesh%n_elems
      gamma_arr(i) = 1.0_DOUBLE + sol(1, i) / sol(6, i)
    end do
  end subroutine sync_gamma_arr

  ! ----------------------------------------------------------------
  ! RHS: -1/V * sum_faces(int_face F·n dA), using face quadrature
  ! ----------------------------------------------------------------
  subroutine compute_rhs(mesh, sol, prim, rhs, sum_lambda, t, num_procs, mpi_send_recv)
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(in)  :: sol
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in)  :: prim
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(out) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems),    intent(out) :: sum_lambda
    ! Current stage's time, for a time-dependent BC (e.g. 'dmr_top'). Not
    ! tracked per-RK-substage (all 3 SSP-RK3 stages of one step reuse the
    ! step's start time) -- a minor, deliberate simplification for a
    ! qualitative test, not exact substage timing.
    real(kind=DOUBLE), intent(in) :: t
    ! Optional: only needed by aho_reconstruction under MPI (see its header
    ! comment) -- absent, compute_rhs behaves exactly as before (serial, or
    ! the ls_reconstruction path which needs no cross-rank exchange here).
    integer(kind=ENTIER), intent(in), optional :: num_procs
    type(mpi_send_recv_type), intent(inout), optional :: mpi_send_recv

    ! Per-cell gradients (5 primitives x 3 spatial dims)
    real(kind=DOUBLE), allocatable :: grad(:, :, :)  ! (5, 3, n_elems)
    real(kind=DOUBLE), allocatable :: hess(:, :, :, :) ! (5, 3, 3, n_elems)
    real(kind=DOUBLE), allocatable :: third(:, :, :, :, :) ! (5, 3, 3, 3, n_elems), order 4 only

    allocate(grad(5, 3, mesh%n_elems))
    grad = 0.0_DOUBLE

    if (order >= 3) then
      allocate(hess(5, 3, 3, mesh%n_elems))
      hess = 0.0_DOUBLE
    end if

    if (order >= 4) then
      ! ls_reconstruction only ever fits up to a quadratic (order-3)
      ! polynomial (see its own header) -- no silent order-3-in-disguise
      ! fallback here if someone asks for order 4 with it.
      if (.not. use_aho_reconstruction) then
        print *, 'FATAL: order>=4 requires use_aho_reconstruction=.true. -- ', &
          'ls_reconstruction does not fit a cubic term (see its header comment).'
        error stop 1
      end if
      allocate(third(5, 3, 3, 3, mesh%n_elems))
      third = 0.0_DOUBLE
    end if

    if (order >= 2) then
      if (use_aho_reconstruction) then
        call aho_reconstruction(mesh, prim, grad, hess, third, num_procs, mpi_send_recv)
      else
        call ls_reconstruction(mesh, prim, grad, hess)
      end if
    end if

    rhs        = 0.0_DOUBLE
    sum_lambda = 0.0_DOUBLE
    call face_flux_loop(mesh, sol, prim, grad, hess, rhs, sum_lambda, t, third)

    ! Divide by cell volume (all 6 components, including the gamma-
    ! transport one)
    block
      integer(kind=ENTIER) :: i
      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) then
          rhs(:, i) = rhs(:, i) / mesh%elem(i)%volume
        end if
      end do
    end block


    if (allocated(third)) deallocate(third)
    if (allocated(hess)) deallocate(hess)
    deallocate(grad)
  end subroutine compute_rhs

  ! ----------------------------------------------------------------
  ! CFL time step
  ! ----------------------------------------------------------------
  function compute_dt(mesh, sum_lambda) result(dt)
    use mpi
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: sum_lambda
    real(kind=DOUBLE) :: dt, dt_i
    integer(kind=ENTIER) :: i
    integer :: mpi_ierr

    dt = huge(1.0_DOUBLE)
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost .and. sum_lambda(i) > 0.0_DOUBLE) then
        dt_i = cfl * mesh%elem(i)%volume / sum_lambda(i)
        dt = min(dt, dt_i)
      end if
    end do

    ! Found 2026-09-14: this was a purely local (per-rank) minimum, with no
    ! cross-rank reduction -- unlike ns_euler_module's own compute_dt
    ! (ns_euler_module.F90, MPI_ALLREDUCE/MPI_MIN), which is why subfvns
    ! never showed this problem on the same meshes. Under MPI, each rank's
    ! own dt differs (most sharply near a strong feature like a bow shock,
    ! where sum_lambda is far larger than in quiescent far-field cells), so
    ! t = t + dt drifts apart rank-to-rank; a far-field rank can reach
    ! t>=tmax and hit MPI_FINALIZE while a shock-region rank is still deep
    ! in its own loop expecting to exchange with it -- an eternal
    ! PMPI_Waitall inside mpi_memory_exchange on the ranks left behind.
    ! Reproduced deterministically (same physical t, not the same iteration
    ! count, across a cfl change) on cylinder-tri at 6 ranks locally, no
    ! cluster involved; confirmed via gdb backtraces on multiple stuck
    ! cluster runs (WC, DMR, Shu-Osher, cylinder-quad, cylinder-tri, FFS)
    ! all landing in the same mpi_memory_exchange/PMPI_Waitall frame.
    call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_DOUBLE, &
      MPI_MIN, MPI_COMM_WORLD, mpi_ierr)
  end function compute_dt

  ! ----------------------------------------------------------------
  ! Diagnostic: k-exactness of reconstruct() on a synthetic polynomial
  ! field of known degree, for both ls_reconstruction and
  ! aho_reconstruction. degree=1 (affine field) exercises the gradient
  ! path (euler_ho order>=2); degree=2 (quadratic field) exercises the
  ! grad+Hessian path (order 3). mesh must be a uniform structured
  ! quad/hex grid of spacing (dx,dy) in x,y (boundary_2d convention, z
  ! dropped) so that a quadratic field's exact cell average has the
  ! closed form used below (avg(x^2) = xc^2 + dx^2/12, etc.) -- no
  ! quadrature needed, so this is independent of quadrature_module.
  ! Prints the Linf error, over interior cells only (|xc|,|yc| < 0.6 *
  ! half-domain-width, away from the one-sided-fit boundary layer) and
  ! 4 sample points per cell, for each reconstruction source.
  ! This does NOT touch module-global `order` except to set it for the
  ! duration of the call (both `ls_reconstruction` and
  ! `aho_reconstruction` branch on it), restoring it before returning.
  subroutine test_reconstruction_exactness(mesh, degree, dx, dy, half_width, order_override)
    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: degree  ! 1 or 2
    real(kind=DOUBLE), intent(in) :: dx, dy, half_width
    ! Normally the method's order is degree+1 (the "matched", k-exact case
    ! this test exists to check). Passing order_override lets a caller
    ! deliberately MISMATCH them instead -- e.g. order_override=2 with
    ! degree=2 checks that a degree-2 field is NOT reproduced exactly by
    ! the gradient-only (order-2) path, which has no Hessian term to
    ! represent curvature with. A non-zero, O(h^2)-scale error there is
    ! the expected, correct outcome -- it is a negative control confirming
    ! this test would actually catch a real exactness failure, rather than
    ! trivially reporting ~0 for any input.
    integer(kind=ENTIER), intent(in), optional :: order_override

    real(kind=DOUBLE), dimension(5, mesh%n_elems) :: prim
    real(kind=DOUBLE), dimension(5, 3, mesh%n_elems) :: grad_ls, grad_aho
    real(kind=DOUBLE), dimension(:, :, :, :), allocatable :: hess_ls, hess_aho
    real(kind=DOUBLE), dimension(5) :: A, B, C, D, E, F
    real(kind=DOUBLE), dimension(3) :: xq, xc_vec
    real(kind=DOUBLE), dimension(5) :: wexact, w_ls, w_aho
    real(kind=DOUBLE) :: xc, yc, err_ls, err_aho
    integer(kind=ENTIER) :: i, v, s, saved_order
    real(kind=DOUBLE), dimension(4, 2) :: offsets

    saved_order = order
    if (present(order_override)) then
      order = order_override
    else
      order = degree + 1
    end if

    ! Large positive base for rho (v=1) and p (v=5) so physical_state's
    ! positivity check never rejects the reconstructed candidate over the
    ! sampled domain (|x|,|y| up to half_width) -- small slopes/curvatures
    ! keep the field close to its base value everywhere it is sampled.
    A = [100.0_DOUBLE, -0.7_DOUBLE, 2.3_DOUBLE, 0.4_DOUBLE, 100.0_DOUBLE]
    B = [0.07_DOUBLE, 0.13_DOUBLE, -0.05_DOUBLE, 0.09_DOUBLE, 0.02_DOUBLE]
    C = [-0.03_DOUBLE, 0.06_DOUBLE, 0.11_DOUBLE, -0.08_DOUBLE, 0.04_DOUBLE]
    D = [0.013_DOUBLE, -0.009_DOUBLE, 0.021_DOUBLE, 0.005_DOUBLE, -0.017_DOUBLE]
    E = [-0.021_DOUBLE, 0.011_DOUBLE, -0.008_DOUBLE, 0.019_DOUBLE, 0.006_DOUBLE]
    F = [0.017_DOUBLE, 0.008_DOUBLE, -0.014_DOUBLE, -0.006_DOUBLE, 0.012_DOUBLE]

    do i = 1, mesh%n_elems
      xc = mesh%elem(i)%coord(1)
      yc = mesh%elem(i)%coord(2)
      do v = 1, 5
        prim(v, i) = A(v) + B(v)*xc + C(v)*yc
        if (degree >= 2) then
          prim(v, i) = prim(v, i) + D(v)*(xc**2 + dx**2/12.0_DOUBLE) &
                                   + E(v)*(xc*yc) &
                                   + F(v)*(yc**2 + dy**2/12.0_DOUBLE)
        end if
      end do
    end do

    if (order >= 3) then
      allocate(hess_ls(5, 3, 3, mesh%n_elems), hess_aho(5, 3, 3, mesh%n_elems))
      hess_ls = 0.0_DOUBLE; hess_aho = 0.0_DOUBLE
    end if
    grad_ls = 0.0_DOUBLE; grad_aho = 0.0_DOUBLE

    call ls_reconstruction(mesh, prim, grad_ls, hess_ls)
    call aho_reconstruction(mesh, prim, grad_aho, hess_aho)

    offsets(1, :) = [ 0.30_DOUBLE,  0.20_DOUBLE]
    offsets(2, :) = [-0.40_DOUBLE,  0.10_DOUBLE]
    offsets(3, :) = [ 0.15_DOUBLE, -0.35_DOUBLE]
    offsets(4, :) = [-0.20_DOUBLE, -0.25_DOUBLE]

    err_ls = 0.0_DOUBLE; err_aho = 0.0_DOUBLE

    do i = 1, mesh%n_elems
      if (mesh%elem(i)%is_ghost) cycle
      xc = mesh%elem(i)%coord(1)
      yc = mesh%elem(i)%coord(2)
      if (abs(xc) > 0.6_DOUBLE*half_width .or. abs(yc) > 0.6_DOUBLE*half_width) cycle
      xc_vec = mesh%elem(i)%coord

      do s = 1, 4
        xq(1) = xc + offsets(s, 1)*dx
        xq(2) = yc + offsets(s, 2)*dy
        xq(3) = xc_vec(3)

        do v = 1, 5
          wexact(v) = A(v) + B(v)*xq(1) + C(v)*xq(2)
          if (degree >= 2) then
            wexact(v) = wexact(v) + D(v)*xq(1)**2 + E(v)*xq(1)*xq(2) + F(v)*xq(2)**2
          end if
        end do

        w_ls  = reconstruct(prim, grad_ls,  hess_ls,  i, xq, xc_vec, gamma_gas)
        w_aho = reconstruct(prim, grad_aho, hess_aho, i, xq, xc_vec, gamma_gas)

        err_ls  = max(err_ls,  maxval(abs(w_ls  - wexact)))
        err_aho = max(err_aho, maxval(abs(w_aho - wexact)))
      end do
    end do

    print *, 'RECON_EXACTNESS degree=', degree, ' order=', order, &
      ' Linf(ls)=', err_ls, ' Linf(aho)=', err_aho

    if (allocated(hess_ls))  deallocate(hess_ls)
    if (allocated(hess_aho)) deallocate(hess_aho)
    order = saved_order
  end subroutine test_reconstruction_exactness

  ! ----------------------------------------------------------------
  ! Error on isentropic vortex (L2 on rho)
  ! ----------------------------------------------------------------
  subroutine compute_error_vortex(mesh, sol, t, h, l2err)
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), intent(in) :: t
    real(kind=DOUBLE), intent(out) :: h, l2err

    integer(kind=ENTIER) :: i, kf, iface_loc, n_fv, k, n_inner
    real(kind=DOUBLE) :: err2, area_tot, area_q, area_2d
    real(kind=DOUBLE), dimension(5) :: wexact, wsol
    real(kind=DOUBLE), dimension(3) :: xc
    real(kind=DOUBLE), dimension(:, :), allocatable :: qpts_v, fc
    real(kind=DOUBLE), dimension(:),    allocatable :: qwts_v
    logical :: found_zface

    err2     = 0.0_DOUBLE
    area_tot = 0.0_DOUBLE
    n_inner  = 0

    do i = 1, mesh%n_elems
      if (mesh%elem(i)%is_ghost) cycle
      xc = mesh%elem(i)%coord
      if (abs(xc(1)) > 5.0_DOUBLE .or. abs(xc(2)) > 5.0_DOUBLE) cycle

      ! Cell average of the exact solution, via genuine volume quadrature
      ! (works for any cell type, hex/tet/prism/pyramid alike -- see
      ! vortex_cell_average's header for why this replaced the previous
      ! z-face-projection approach, which relied on a guaranteed z-normal
      ! face that a tet cell essentially never has). Mathematically
      ! identical to the old approach for a z-invariant field on a
      ! z-extruded boundary_2d mesh (the two are the same integral divided
      ! by the same volume, just computed via a full 3D rule instead of a
      ! face rule x cancelled z-thickness), so this does not change any
      ! already-published boundary_2d result.
      call vortex_cell_average(mesh, i, t, wexact)

      ! h/area_2d bookkeeping is untouched (still needs a real z-face
      ! for the boundary_2d convention's h = sqrt(area) definition): find
      ! one purely to measure area_q below, independent of wexact now.
      found_zface = .false.
      do kf = 1, mesh%elem(i)%n_faces
        iface_loc = mesh%elem(i)%face(kf)
        if (abs(abs(mesh%face(iface_loc)%norm(3)) - 1.0_DOUBLE) < 1.0e-6_DOUBLE) then
          n_fv = mesh%face(iface_loc)%n_vert
          if (n_fv >= 3 .and. allocated(mesh%face(iface_loc)%vert)) then
            allocate(fc(3, n_fv))
            do k = 1, n_fv
              fc(:, k) = mesh%vert(mesh%face(iface_loc)%vert(k))%coord
            end do
            call face_quad_pts(n_fv, fc, max(order, 3), qpts_v, qwts_v)
            deallocate(fc)
            area_q = sum(qwts_v)
            deallocate(qpts_v, qwts_v)
            found_zface = .true.
            exit
          end if
        end if
      end do

      wsol = conserv_to_primit(sol(1:5, i), gamma_arr(i))
      ! Use z-face area as 2D cell area so h = sqrt(avg_area) = h_2D
      ! independent of dz. Fallback to vol^(2/3) if no z-face was found.
      if (found_zface .and. area_q > 0.0_DOUBLE) then
        area_2d = area_q
      else
        area_2d = mesh%elem(i)%volume**(2.0_DOUBLE/3.0_DOUBLE)
      end if
      err2     = err2     + area_2d * (wsol(1) - wexact(1))**2
      area_tot = area_tot + area_2d
      n_inner  = n_inner  + 1
    end do

    ! Under MPI, each rank only owns a slice of the domain (is_ghost cells
    ! are excluded above, so no cell is double-counted across ranks) --
    ! err2/area_tot/n_inner must be summed over all ranks before taking the
    ! sqrt, or a multi-rank run silently reports only rank 0's own local
    ! partition's error (found 2026-09-13 while setting up 6-rank vortex
    ! convergence runs: L2(rho) came back orders of magnitude too small).
    call mpi_sum_error_vortex(err2, area_tot, n_inner)

    if (n_inner > 0) then
      h = sqrt(area_tot / n_inner)   ! = h_2D for any dz
    else
      h = 0.0_DOUBLE
    end if
    l2err = sqrt(err2)
  end subroutine compute_error_vortex

  ! In-place MPI_SUM of the three running accumulators compute_error_vortex
  ! builds from this rank's own (non-ghost) cells only. A no-op under a
  ! single rank (MPI_ALLREDUCE with a 1-rank communicator returns its input
  ! unchanged), so this is safe to call unconditionally.
  subroutine mpi_sum_error_vortex(err2, area_tot, n_inner)
    use mpi
    implicit none

    real(kind=DOUBLE), intent(inout) :: err2, area_tot
    integer(kind=ENTIER), intent(inout) :: n_inner

    real(kind=DOUBLE) :: buf_real(2)
    integer(kind=ENTIER) :: n_inner_sum
    integer :: mpi_ierr

    call MPI_ALLREDUCE((/err2, area_tot/), buf_real, 2, MPI_DOUBLE_PRECISION, &
      MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
    err2     = buf_real(1)
    area_tot = buf_real(2)

    call MPI_ALLREDUCE(n_inner, n_inner_sum, 1, MPI_INTEGER, MPI_SUM, &
      MPI_COMM_WORLD, mpi_ierr)
    n_inner = n_inner_sum
  end subroutine mpi_sum_error_vortex

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
  !
  ! Tops out at order 3 (quadratic, n_coeff above): no cubic/order-4 basis
  ! is fitted here. compute_rhs enforces this explicitly (error stop if
  ! order>=4 and use_aho_reconstruction is false) rather than silently
  ! reusing the order-3 fit under an order-4 label -- growing this to a
  ! cubic fit would also need a wider neighbor stencil (ls's ring widens
  ! with order), not just more basis functions.
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

      ! Same non-guaranteed-short-circuit hazard as face_flux_loop's ir
      ! check (see 2026-09-15 comment there): nest the sign guard so
      ! mesh%elem(cand) is never indexed for a non-positive (boundary
      ! marker) cand.
      if (cand <= 0) return
      if (cand == i .or. mesh%elem(cand)%is_ghost) return
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

  ! Drop-in alternative to ls_reconstruction: grad/hess of the 5 primitive
  ! variables from arbitrary_high_order_module's least-squares dual/primal
  ! hierarchy. The 5 variables are differentiated together as one nc_in=5
  ! "vector field" (the module's least-squares fit already treats each
  ! component as an independent right-hand side against a shared,
  ! geometry-only normal matrix, so this is exactly as one call per
  ! variable would be, just batched).
  !
  ! compute_next_order_derivative flattens its output with the *input*
  ! component fast-varying and the new derivative direction slow-varying:
  ! grad_flat((dir-1)*5 + v, e) = d(prim_v)/dx_dir. Differentiating that
  ! again gives hess_flat((dir2-1)*15 + (dir1-1)*5 + v, e) = d2(prim_v)/(dx_dir1 dx_dir2).
  !
  ! MPI (2026-09-15, per the user's fresh 1-proc-vs-2-proc vortex test):
  ! compute_next_order_derivative computes grad_flat/hess_flat for EVERY
  ! cell it's given, ghost cells included, using this rank's OWN local
  ! vertex-neighbor topology (mesh%vert%elem_neigh) for that cell's
  ! vertices. A ghost cell's own vertices sit one hop further from the
  ! interface than the interface itself -- exactly where subfv-gmsh's
  ! single node-based ghost layer no longer guarantees a *complete*
  ! elem_neigh (it only guarantees that for vertices directly touching an
  ! owned cell). So a ghost cell's locally-computed grad/hess can be wrong
  ! (built from a genuinely incomplete local stencil), yet it's exactly
  ! what face_flux_loop uses for that ghost cell's own Taylor extrapolation
  ! at the shared face with a real, owned cell -- corrupting that owned
  ! cell's flux, then its evolution, then (over many RK3 steps) the whole
  ! rank-adjacent region. compute_derivative_hierarchy (this module's own
  ! sibling driver, used by arbitrary_high_order_main, not by this solver)
  ! already exchanges every order's tensor before using it as the next
  ! order's input for exactly this reason -- this routine just wasn't
  ! doing the analogous thing, for either stage. Confirmed by two direct
  ! tests: (1) subfvns, a classical single-ring-stencil solver sharing the
  ! same partitioned mesh, gives IDENTICAL error at 1 and 2 ranks; (2) this
  ! routine's order=2 (grad only, no hess) already showed the same
  ! MPI-dependent error inflation as order=3, ruling out the hess stage
  ! specifically and pointing at grad_flat itself being ghost-incomplete.
  ! Fix: exchange grad_flat right after it's computed (both because
  ! face_flux_loop reads it directly at order 2, and because it's hess's
  ! own input at order 3), and exchange hess_flat right after IT is
  ! computed, before face_flux_loop reads it. num_procs/mpi_send_recv are
  ! optional so this remains a plain serial call wherever nothing MPI is
  ! in play (verify_reconstruction_order's synthetic single-rank tests).
  ! third (order-4 cubic term) is optional and follows the exact same
  ! recursive pattern one step further: hess_flat (45 components/cell,
  ! nc_in for this step) goes into compute_next_order_derivative to give
  ! third_flat (135 components/cell), exchanged across ranks the same way,
  ! then unflattened with the third (newest) direction slow-varying over
  ! hess_flat's own index -- same convention as hess_flat's own layout
  ! relative to grad_flat, one level up.
  subroutine aho_reconstruction(mesh, prim, grad, hess, third, num_procs, mpi_send_recv)
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: prim
    real(kind=DOUBLE), dimension(5, 3, mesh%n_elems), intent(inout) :: grad
    real(kind=DOUBLE), dimension(:, :, :, :), allocatable, intent(inout) :: hess
    real(kind=DOUBLE), dimension(:, :, :, :, :), allocatable, intent(inout), optional :: third
    integer(kind=ENTIER), intent(in), optional :: num_procs
    type(mpi_send_recv_type), intent(inout), optional :: mpi_send_recv

    real(kind=DOUBLE), dimension(:, :), allocatable :: grad_flat, hess_flat, third_flat
    integer(kind=ENTIER) :: e, v, dir1, dir2, dir3
    logical :: do_exchange
    ! Per-vertex nodal Hessian/third-derivative estimates (before WENO
    ! scatter into cells), needed by apply_grad_bias_correction -- see
    ! use_grad_bias_correction's header.
    real(kind=DOUBLE), dimension(:, :), allocatable :: hess_v, third_v
    logical, dimension(:), allocatable :: valid_hess_v, valid_third_v
    ! Per-vertex oscillation indicator from the GRADIENT-level call,
    ! needed by apply_grad_bias_correction so its own vertex-to-cell
    ! scatter matches grad_flat's actual WENO weighting -- see that
    ! subroutine's header. Zero-filled (-> uniform weight) when
    ! use_cweno_center=.true., since compute_next_order_derivative_cweno
    ! does not expose an oi_v_out (that path is retired/off by default;
    ! not fixed here).
    real(kind=DOUBLE), dimension(:), allocatable :: grad_oi_v

    do_exchange = present(num_procs)
    if (do_exchange) do_exchange = num_procs > 1

    ! CWENO-style central candidate (2026-09-15, replaces an earlier,
    ! cruder "WENO only on the grad step" interim fix): rather than
    ! skipping WENO's nonlinear oscillation-indicator weight at the
    ! deeper recursion levels (hess from grad, third from hess), fold a
    ! large-weight linear/center candidate into the SAME nonlinear blend
    ! at EVERY level via compute_next_order_derivative_cweno -- see that
    ! subroutine's header and arbitrary_high_order_module's
    ! use_cweno_center/cweno_center_weight. This keeps the per-vertex
    ! nonlinear WENO candidates genuinely present (for real shock
    ! detection) at every level, unlike the interim fix, while still
    ! letting hess/third converge with a real order on smooth data
    ! (verified via test_order4_cweno.F90, scratchpad, on a synthetic
    ! cubic field: hess/third residual ~1/cweno_center_weight, real
    ! O(h^2-3) convergence, vs. O(h^0)/non-convergent with plain WENO at
    ! every level). Falls back to the original compute_next_order_derivative
    ! (WENO or linear per use_weno_blend, uniformly at every level) when
    ! use_cweno_center=.false.
    aho_module_use_weno_blend = use_weno_blend
    aho_module_use_green_gauss = use_green_gauss
    aho_module_use_cweno_center = use_cweno_center
    aho_module_cweno_center_weight = cweno_center_weight
    aho_module_cweno_center_power = cweno_center_power
    aho_module_eps_weight_num = eps_weight_num
    aho_module_eps_weight_num_deep = eps_weight_num_deep
    aho_module_weno_power = weno_power
    aho_module_grad_norm_derate = grad_norm_derate
    aho_module_grad_norm_derate_gg = grad_norm_derate_gg
    aho_module_use_alt_gg_weight = use_alt_gg_weight
    aho_module_eps_weight_num_gg = eps_weight_num_gg
    aho_module_eps_weight_num_deep_gg = eps_weight_num_deep_gg

    allocate(grad_flat(15, mesh%n_elems))
    allocate(grad_oi_v(mesh%n_vert))
    grad_oi_v = 0.0_DOUBLE
    if (use_cweno_center) then
      call compute_next_order_derivative_cweno(mesh, 3_ENTIER, 5_ENTIER, boundary_2d, &
        prim, grad_flat)
    else
      call compute_next_order_derivative(mesh, 3_ENTIER, 5_ENTIER, boundary_2d, &
        prim, grad_flat, deriv_order=1_ENTIER, oi_v_out=grad_oi_v)
    end if

    if (do_exchange) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 15_ENTIER, grad_flat)

    ! grad_flat is unpacked into `grad` further below, AFTER the order-4
    ! bias correction (if any) has had a chance to correct it in place --
    ! see use_grad_bias_correction's header.
    if (order >= 3 .and. allocated(hess)) then
      allocate(hess_flat(45, mesh%n_elems))
      if (use_cweno_center) then
        call compute_next_order_derivative_cweno(mesh, 3_ENTIER, 15_ENTIER, boundary_2d, &
          grad_flat, hess_flat)
      else
        block
          real(kind=DOUBLE), dimension(:, :), allocatable :: hess_v_local
          logical, dimension(:), allocatable :: valid_hess_v_local
          logical :: need_hess_v
          need_hess_v = (order >= 4 .and. present(third) .and. use_grad_bias_correction)
          if (need_hess_v) then
            call compute_next_order_derivative(mesh, 3_ENTIER, 15_ENTIER, boundary_2d, &
              grad_flat, hess_flat, deriv_order=2_ENTIER, &
              dphi_v_out=hess_v_local, valid_v_out=valid_hess_v_local)
            if (allocated(hess_v)) deallocate(hess_v)
            if (allocated(valid_hess_v)) deallocate(valid_hess_v)
            allocate(hess_v(size(hess_v_local,1), size(hess_v_local,2)))
            allocate(valid_hess_v(size(valid_hess_v_local)))
            hess_v = hess_v_local
            valid_hess_v = valid_hess_v_local
          else
            call compute_next_order_derivative(mesh, 3_ENTIER, 15_ENTIER, boundary_2d, &
              grad_flat, hess_flat, deriv_order=2_ENTIER)
          end if
        end block
      end if
      if (do_exchange) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 45_ENTIER, hess_flat)

      do e = 1, mesh%n_elems
        do dir2 = 1, 3
          do dir1 = 1, 3
            do v = 1, 5
              hess(v, dir1, dir2, e) = hess_flat((dir2-1)*15 + (dir1-1)*5 + v, e)
            end do
          end do
        end do
      end do

      if (order >= 4 .and. present(third)) then
        if (allocated(third)) then
          allocate(third_flat(135, mesh%n_elems))
          if (use_cweno_center) then
            call compute_next_order_derivative_cweno(mesh, 3_ENTIER, 45_ENTIER, boundary_2d, &
              hess_flat, third_flat)
          else
            if (use_grad_bias_correction) then
              block
                real(kind=DOUBLE), dimension(:, :), allocatable :: third_v_local
                logical, dimension(:), allocatable :: valid_third_v_local
                call compute_next_order_derivative(mesh, 3_ENTIER, 45_ENTIER, boundary_2d, &
                  hess_flat, third_flat, deriv_order=3_ENTIER, &
                  dphi_v_out=third_v_local, valid_v_out=valid_third_v_local)
                if (allocated(third_v)) deallocate(third_v)
                if (allocated(valid_third_v)) deallocate(valid_third_v)
                allocate(third_v(size(third_v_local,1), size(third_v_local,2)))
                allocate(valid_third_v(size(valid_third_v_local)))
                third_v = third_v_local
                valid_third_v = valid_third_v_local
              end block
              ! Correct grad_flat IN PLACE using the per-vertex nodal
              ! hess_v/third_v just computed above -- see
              ! apply_grad_bias_correction's header for the full
              ! derivation. Only implemented/safe for boundary_2d; a
              ! .false. call is a documented no-op.
              call apply_grad_bias_correction(mesh, 3_ENTIER, 5_ENTIER, boundary_2d, &
                grad_flat, hess_v, third_v, valid_hess_v, valid_third_v, grad_oi_v)
              ! grad_flat was corrected IN PLACE above, for this rank's own
              ! owned cells only (apply_grad_bias_correction has no notion
              ! of ownership, but a cell's bias is only ever accumulated
              ! from ITS OWN touching vertices, computed identically on
              ! whichever rank owns each vertex). Found 2026-09-17 (MPI
              ! order-4 vs serial comparison, isentropic vortex): without
              ! re-exchanging here, a rank's GHOST copy of a neighbor
              ! rank's cell still holds the PRE-correction grad_flat, so a
              ! face straddling a partition boundary reconstructs wL
              ! (locally owned, corrected) against wR (ghost, uncorrected)
              ! -- a real, partition-count-dependent ~4-6% L2 discrepancy
              ! on the same vortex test at 6 MPI ranks vs 1 (confirmed
              ! absent with use_grad_bias_correction=.false.: no
              ! correction, no discrepancy, at any rank count -- isolating
              ! this call as the cause before this fix).
              if (do_exchange) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 15_ENTIER, grad_flat)
            else
              call compute_next_order_derivative(mesh, 3_ENTIER, 45_ENTIER, boundary_2d, &
                hess_flat, third_flat, deriv_order=3_ENTIER)
            end if
          end if
          if (do_exchange) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 135_ENTIER, third_flat)

          do e = 1, mesh%n_elems
            do dir3 = 1, 3
              do dir2 = 1, 3
                do dir1 = 1, 3
                  do v = 1, 5
                    third(v, dir1, dir2, dir3, e) = &
                      third_flat((dir3-1)*45 + (dir2-1)*15 + (dir1-1)*5 + v, e)
                  end do
                end do
              end do
            end do
          end do
          deallocate(third_flat)
        end if
      end if
      deallocate(hess_flat)
    end if

    ! grad's own unpacking, deferred until here (regardless of order) so
    ! the order-4 bias correction above (if it ran) has already updated
    ! grad_flat in place.
    do e = 1, mesh%n_elems
      do dir1 = 1, 3
        do v = 1, 5
          grad(v, dir1, e) = grad_flat((dir1-1)*5 + v, e)
        end do
      end do
    end do

    deallocate(grad_flat)
    if (allocated(grad_oi_v)) deallocate(grad_oi_v)
    if (allocated(hess_v)) deallocate(hess_v)
    if (allocated(valid_hess_v)) deallocate(valid_hess_v)
    if (allocated(third_v)) deallocate(third_v)
    if (allocated(valid_third_v)) deallocate(valid_third_v)
  end subroutine aho_reconstruction

  ! Solves ATA_in * x = rhs_in (ATA_in = A^T W A, a Gram/normal matrix,
  ! always symmetric PSD). Drops directions with no information in the
  ! data -- e.g. the y-gradient unknown on a mesh that is a single row of
  ! cells in y (ls_reconstruction's Sod-tube mesh: boundary_2d, every
  ! neighbour at the same y, so the whole y-row/column of ATA is exactly
  ! zero) -- and solves only the reduced, well-posed subsystem for the
  ! rest, rather than giving up on the whole vector.
  !
  ! BUG FIXED (previously): a plain Gaussian elimination with partial
  ! pivoting returned x=0 for EVERY unknown, including well-determined
  ! ones, the instant it hit any near-zero pivot -- so on the Sod mesh
  ! above, the perfectly-determined x-gradient was ALSO silently zeroed
  ! alongside the genuinely-degenerate y one, making ls_reconstruction
  ! produce bit-identical output to plain order 1 at every order (2 and
  ! 3) on that whole mesh family. Because ATA is a Gram matrix, a
  ! genuinely-degenerate direction j has its ENTIRE row/column j equal to
  ! zero (not just a small pivot reached after row operations), so this
  ! is detected directly from the untouched input, up front, rather than
  ! discovered mid-elimination -- and only actually-degenerate directions
  ! are dropped.
  subroutine solve_ls(ATA_in, rhs_in, x, n)
    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), dimension(9, 9), intent(in) :: ATA_in
    real(kind=DOUBLE), dimension(9),    intent(in) :: rhs_in
    real(kind=DOUBLE), dimension(9),    intent(out) :: x

    integer(kind=ENTIER), dimension(9) :: idx
    integer(kind=ENTIER) :: n_keep, ii, jj, fail_at
    real(kind=DOUBLE), dimension(9, 9) :: ATA_red
    real(kind=DOUBLE), dimension(9) :: rhs_red, x_red
    real(kind=DOUBLE) :: col_tol

    x = 0.0_DOUBLE

    ! Relative to the matrix's own scale, same reasoning as gauss_solve's
    ! piv_tol: an absolute floor here misclassifies a genuinely-alive but
    ! naturally-small quadratic-term column (scale ~dx^2*wt, vs ~dx*wt for
    ! the linear-term columns) as degenerate purely from that scale gap.
    col_tol = 1.0e-12_DOUBLE * maxval(abs(ATA_in(1:n, 1:n)))

    n_keep = 0
    do jj = 1, n
      if (any(abs(ATA_in(1:n, jj)) > col_tol)) then
        n_keep = n_keep + 1
        idx(n_keep) = jj
      end if
    end do

    ! Iteratively solve the currently-active subset (idx(1:n_keep)); if
    ! gauss_solve hits a genuinely-degenerate pivot MID-elimination (a
    ! direction that looked alive column-wise above but turns out to be a
    ! linear combination of others once the others are eliminated against
    ! it -- e.g. a symmetric structured-mesh stencil where some quadratic
    ! cross-term direction is exactly indeterminate), drop just that one
    ! direction and retry with the rest, rather than discarding the whole
    ! vector. Found necessary 2026-09-14: on a genuinely-3D pseudo-1D mesh
    ! (regular hex cells), the old single-shot gauss_solve threw away
    ! EVERY coefficient, including perfectly well-determined ones like the
    ! plain x-gradient, the moment ANY one of the 9 order-3 directions hit
    ! a degenerate pivot -- caught because order-3 ls_reconstruction was
    ! found not exact (Linf 5e-4, not machine precision) on a degree-2
    ! test polynomial on that mesh.
    do
      if (n_keep == 0) return

      do ii = 1, n_keep
        do jj = 1, n_keep
          ATA_red(ii, jj) = ATA_in(idx(ii), idx(jj))
        end do
        rhs_red(ii) = rhs_in(idx(ii))
      end do

      call gauss_solve(ATA_red, rhs_red, x_red, n_keep, fail_at)

      if (fail_at == 0) then
        do ii = 1, n_keep
          x(idx(ii)) = x_red(ii)
        end do
        return
      end if

      ! Drop the direction that failed and retry with the rest.
      do ii = fail_at, n_keep - 1
        idx(ii) = idx(ii + 1)
      end do
      n_keep = n_keep - 1
    end do
  end subroutine solve_ls

  ! Gaussian elimination with partial pivoting for small dense systems.
  ! Solves ATA_in * x = rhs_in. On a degenerate pivot at step jj, returns
  ! immediately with fail_at=jj (x left at 0) instead of guessing --
  ! solve_ls's caller loop drops direction jj and retries with the rest,
  ! rather than discarding every other, well-determined direction along
  ! with the one genuinely-degenerate one (see solve_ls's own comment).
  ! fail_at=0 signals a complete, successful solve.
  subroutine gauss_solve(ATA_in, rhs_in, x, n, fail_at)
    integer(kind=ENTIER), intent(in) :: n
    real(kind=DOUBLE), dimension(9, 9), intent(in) :: ATA_in
    real(kind=DOUBLE), dimension(9),    intent(in) :: rhs_in
    real(kind=DOUBLE), dimension(9),    intent(out) :: x
    integer(kind=ENTIER), intent(out) :: fail_at

    real(kind=DOUBLE), dimension(9, 10) :: aug
    integer(kind=ENTIER) :: ii, jj, kk, piv_row
    real(kind=DOUBLE) :: pivot_val, fac, tmp, piv_tol

    x = 0.0_DOUBLE
    fail_at = 0
    aug(:, 1:n)  = ATA_in(1:9, 1:n)
    aug(:, n+1)  = rhs_in

    ! Degeneracy threshold RELATIVE to the matrix's own scale, not an
    ! absolute constant -- an absolute floor misclassifies a legitimately
    ! nonzero but naturally-small pivot (e.g. from a quadratic-term
    ! column, scale ~dx^2*wt, vs ~dx*wt for linear-term columns) as
    ! degenerate purely from that scale gap, not from true rank-deficiency.
    piv_tol = 1.0e-12_DOUBLE * maxval(abs(ATA_in(1:n, 1:n)))

    do jj = 1, n
      piv_row   = jj
      pivot_val = abs(aug(jj, jj))
      do ii = jj + 1, n
        if (abs(aug(ii, jj)) > pivot_val) then
          pivot_val = abs(aug(ii, jj))
          piv_row   = ii
        end if
      end do

      if (piv_row /= jj) then
        do kk = jj, n + 1
          tmp            = aug(jj,      kk)
          aug(jj,      kk) = aug(piv_row, kk)
          aug(piv_row, kk) = tmp
        end do
      end if

      if (abs(aug(jj, jj)) < piv_tol) then
        fail_at = jj
        return
      end if

      do ii = jj + 1, n
        fac = aug(ii, jj) / aug(jj, jj)
        aug(ii, jj:n+1) = aug(ii, jj:n+1) - fac * aug(jj, jj:n+1)
      end do
    end do

    do ii = n, 1, -1
      x(ii) = aug(ii, n + 1)
      do kk = ii + 1, n
        x(ii) = x(ii) - aug(ii, kk) * x(kk)
      end do
      x(ii) = x(ii) / aug(ii, ii)
    end do
  end subroutine gauss_solve

  ! Main face flux loop
  subroutine face_flux_loop(mesh, sol, prim, grad, hess, rhs, sum_lambda, t, third)
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(in)  :: sol
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in)  :: prim
    real(kind=DOUBLE), dimension(5, 3, mesh%n_elems), intent(in) :: grad
    real(kind=DOUBLE), dimension(:, :, :, :), allocatable, intent(in) :: hess
    real(kind=DOUBLE), dimension(6, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems),    intent(inout) :: sum_lambda
    real(kind=DOUBLE), intent(in) :: t
    real(kind=DOUBLE), dimension(:, :, :, :, :), allocatable, intent(in), optional :: third

    integer(kind=ENTIER) :: iface, il, ir, iv, k, n_fvert, n_qpts, q
    real(kind=DOUBLE), dimension(3) :: norm, xface
    real(kind=DOUBLE), dimension(:, :), allocatable :: face_coords, qpts
    real(kind=DOUBLE), dimension(:),    allocatable :: qwts
    real(kind=DOUBLE), dimension(3, 1) :: qpts_single
    real(kind=DOUBLE), dimension(1)    :: qwts_single
    real(kind=DOUBLE) :: qwt
    real(kind=DOUBLE), dimension(5) :: wL, wR, flux
    real(kind=DOUBLE) :: lambda, gL, gR, face_lambda_accum
    ! Passive-scalar flux for the gamma-transport variable sol(6,:) =
    ! rho*Gamma (see gamma_arr's declaration): upwinded by the sign of the
    ! mass flux (flux(1)) already computed below -- same convention as a
    ! species mass fraction in a Godunov scheme, the material identity
    ! moves with the local mass flux.
    real(kind=DOUBLE) :: flux_rgm1, Gl_rgm1, Gr_rgm1
    logical :: is_zface, use_face_quad

    do iface = 1, mesh%n_faces
      il   = mesh%face(iface)%left_neigh
      ir   = mesh%face(iface)%right_neigh
      norm = mesh%face(iface)%norm

      ! Skip ghost-owned faces
      if (mesh%elem(il)%is_ghost) cycle

      ! Ghost/boundary cells (ir <= 0) have no gamma_arr entry of their
      ! own -- assume the same material as the interior cell il. Correct
      ! for a wall/outflow/freestream BC where the boundary does not
      ! coincide with a material interface (the case for every BC used so
      ! far, including shock_bubble's right-edge post-shock inflow, which
      ! is pure air on both sides); would need a per-BC gamma if a future
      ! case put a material interface directly on a domain boundary.
      gL = gamma_arr(il)
      if (ir > 0) then
        gR = gamma_arr(ir)
      else
        gR = gL
      end if

      n_fvert = mesh%face(iface)%n_vert

      ! For z-faces on 2D-extruded meshes (norm_z ≈ ±1, boundary_2d=T), the
      ! solution has no z-variation. A single centroid point is exact AND
      ! avoids spurious net z-flux from GMSH's vertex ordering (top/bottom
      ! faces get different (x,y) Gauss points due to opposite windings).
      is_zface = boundary_2d .and. abs(abs(norm(3)) - 1.0_DOUBLE) < 1.0e-6_DOUBLE

      ! For z-faces in 2D mode: zero net physics, skip entirely -- before
      ! touching qpts/qwts at all, so the single-point fast path below
      ! never allocates/deallocates anything for the (order>=2) branch
      ! either.
      if (is_zface) cycle

      ! Order 2 (linear reconstruction) only needs a single centroid flux
      ! evaluation -- the classical, textbook 2nd-order FV convention.
      ! Found 2026-09-15: this used to take the multi-point face_quad_pts
      ! branch at order 2 as well (n_face_quad_pts gives 4 points on a
      ! quad face there, matching a much higher exactness degree than a
      ! linear reconstruction needs). Even though wL/wR are each exactly
      ! linear along the face, the *numerical flux* F(wL,wR) is a
      ! nonlinear function of them, so the extra points were integrating
      ! that nonlinearity more accurately than a 2nd-order scheme is
      ! supposed to -- silently giving order 2 a quadrature-driven
      ! accuracy boost past its design order, which is why its observed
      ! rate/error sat suspiciously close to order 3's. Only order>=3
      ! (a genuinely quadratic reconstruction) needs the multi-point rule.
      use_face_quad = order >= 3 .and. n_fvert > 0 .and. allocated(mesh%face(iface)%vert)

      if (use_face_quad) then
        ! Quadrature points on face. face_quad_pts allocates qpts/qwts
        ! itself (their size varies with n_fvert/order) -- unavoidable
        ! per-face allocate/deallocate here, but this branch only runs at
        ! order>=3, and only on genuinely curved/multi-point faces.
        !
        ! Found 2026-09-15: convergence order is always polynomial degree
        ! + 1 (order 1 = degree-0/constant, order 2 = degree-1/affine,
        ! order 3 = degree-2/quadratic via the Hessian term), and
        ! quadrature_module's own "order" argument is defined to match
        ! the polynomial DEGREE being integrated ("order N = exact for
        ! degree N"), not the solver's convergence order. Passing `order`
        ! itself here (3 at order=3) asks for degree-3 exactness -- a full
        ! convergence-order too many -- when the reconstructed state is
        ! only ever degree `order-1`. Pass `order-1` instead: a quad face
        ! goes from 9 points (3x3 Gauss, exact to degree 5) down to 4
        ! (2x2, exact to degree 3), still exact for the degree-2
        ! reconstructed state, with no accuracy loss.
        allocate(face_coords(3, n_fvert))
        do k = 1, n_fvert
          iv = mesh%face(iface)%vert(k)
          face_coords(:, k) = mesh%vert(iv)%coord
        end do
        call face_quad_pts(int(n_fvert, ENTIER), face_coords, &
          int(order - 1, ENTIER), qpts, qwts)
        deallocate(face_coords)
        n_qpts = size(qwts)
      else
        ! Order 1 (or a non-curved face at any order): single point at the
        ! face centroid. Fixed-size local arrays -- no allocate/deallocate
        ! at all, since this is the hot path (every face, every RK stage,
        ! every iteration): millions of alloc/dealloc cycles here were
        ! found (2026-09-14) to eventually stall OpenMPI's shared-memory
        ! transport after tens of thousands of iterations (`gdb` on a
        ! stuck rank always landed in `mca_btl_sm_component_progress` /
        ! `opal_progress`, inside a `PMPI_Waitall` that never returns).
        n_qpts = 1
        qpts_single(:, 1) = mesh%face(iface)%coord
        qwts_single(1)    = mesh%face(iface)%area
      end if

      face_lambda_accum = 0.0_DOUBLE
      do q = 1, n_qpts
        if (use_face_quad) then
          xface = qpts(:, q)
          qwt   = qwts(q)
        else
          xface = qpts_single(:, q)
          qwt   = qwts_single(q)
        end if

        wL = reconstruct(prim, grad, hess, il, xface, mesh%elem(il)%coord, gL, third)
        if (ir > 0) then
          wR = reconstruct(prim, grad, hess, ir, xface, mesh%elem(ir)%coord, gR, third)
        else
          ! Boundary: reconstruct left then apply BC
          wR = ghost_prim(xface, norm, wL, -ir, t)
        end if

        ! Numerical flux (physical: per-area * area_weight_at_qpt).
        ! flux_scheme_id resolved once in read_params -- see its
        ! declaration for why this is an integer compare, not a string one.
        select case (flux_scheme_id)
        case (FLUX_THREE_WAVE)
          flux = three_wave_flux(wL, wR, norm, gL, gR) * qwt
        case (FLUX_TWO_WAVE)
          flux = two_wave_flux(wL, wR, norm, gL, gR) * qwt
        case default
          flux = rusanov(wL, wR, norm, gL, gR) * qwt
        end select
        lambda = max(abs(dot_product(wL(2:4), norm)) + cs(wL, gL), &
                     abs(dot_product(wR(2:4), norm)) + cs(wR, gR)) * qwt
        face_lambda_accum = face_lambda_accum + lambda

        ! Gamma-transport passive-scalar flux: upwind Gamma = sol(6,:)/
        ! sol(1,:) (rho*Gamma/rho) by the sign of the mass flux flux(1)
        ! just computed above (the same mass flux that already carries
        ! rhs(1,:)'s continuity equation). Ghost/boundary (ir<=0) reuses
        ! Gl_rgm1, matching the existing gR=gL "same material as interior
        ! neighbour" convention used for gL/gR just above.
        Gl_rgm1 = sol(6, il) / sol(1, il)
        if (ir > 0) then
          Gr_rgm1 = sol(6, ir) / sol(1, ir)
        else
          Gr_rgm1 = Gl_rgm1
        end if
        if (flux(1) >= 0.0_DOUBLE) then
          flux_rgm1 = flux(1) * Gl_rgm1
        else
          flux_rgm1 = flux(1) * Gr_rgm1
        end if

        rhs(1:5, il) = rhs(1:5, il) - flux
        rhs(6, il)   = rhs(6, il)   - flux_rgm1
        if (ir > 0) then
          ! Fortran's .and. does not guarantee short-circuit evaluation
          ! (unlike C) -- a combined "ir > 0 .and. .not. mesh%elem(ir)%is_ghost"
          ! let an unoptimized (-O0) build evaluate mesh%elem(ir) even when
          ! ir was a negative boundary-condition marker, indexing e.g.
          ! mesh%elem(-3) (found 2026-09-15, cluster debug build with
          ! -fcheck=bounds; silently read out-of-bounds memory instead of
          ! crashing in every optimized build up to now). Nest the checks
          ! instead so ir's sign always gates the array access.
          if (.not. mesh%elem(ir)%is_ghost) then
            rhs(1:5, ir) = rhs(1:5, ir) + flux
            rhs(6, ir)   = rhs(6, ir)   + flux_rgm1
          end if
        end if
      end do

      ! Fold this face's total (summed over its own quad points just
      ! above) into sum_lambda(il)/sum_lambda(ir): += for the original
      ! sum-over-faces CFL bound (bit-for-bit identical to summing every
      ! quad point directly, since + is associative), or max() for the
      ! looser max-over-faces bound (use_max_lambda_dt) -- see that
      ! flag's declaration.
      if (use_max_lambda_dt) then
        sum_lambda(il) = max(sum_lambda(il), face_lambda_accum)
      else
        sum_lambda(il) = sum_lambda(il) + face_lambda_accum
      end if
      if (ir > 0) then
        if (.not. mesh%elem(ir)%is_ghost) then
          if (use_max_lambda_dt) then
            sum_lambda(ir) = max(sum_lambda(ir), face_lambda_accum)
          else
            sum_lambda(ir) = sum_lambda(ir) + face_lambda_accum
          end if
        end if
      end if

      if (use_face_quad) deallocate(qpts, qwts)
    end do
  end subroutine face_flux_loop

  ! Returns .true. if w_cand is a physically acceptable reconstruction
  ! relative to the reference state w_ref.
  ! Rejects: negative rho or p; velocity more than 20x the reference speed + c.
  pure function physical_state(w_cand, w_ref, g) result(ok)
    real(kind=DOUBLE), dimension(5), intent(in) :: w_cand, w_ref
    real(kind=DOUBLE), intent(in) :: g
    logical :: ok
    real(kind=DOUBLE) :: spd_ref, spd_cand, c_ref

    ok = .false.
    if (w_cand(1) <= 0.0_DOUBLE) return
    if (w_cand(5) <= 0.0_DOUBLE) return

    ! Velocity magnitude check
    spd_ref  = w_ref(2)**2  + w_ref(3)**2  + w_ref(4)**2
    spd_cand = w_cand(2)**2 + w_cand(3)**2 + w_cand(4)**2
    c_ref    = g * w_ref(5) / max(w_ref(1), 1.0e-16_DOUBLE)
    if (spd_cand > 400.0_DOUBLE * (spd_ref + c_ref)) return   ! 20x speed limit

    ok = .true.
  end function physical_state

  ! Second (and, at order>=4, third) geometric moment of every cell about
  ! its own centroid, M_jk(i) = (1/V_i) * int_cell (x_j-xc_j)(x_k-xc_k) dV
  ! [and M_jkl(i), the same one degree up], needed so that reconstruct()'s
  ! order-3 (quadratic) and order-4 (cubic) Taylor terms have the correct
  ! CELL AVERAGE (equal to the given prim(i)), not merely the correct
  ! VALUE at the centroid -- for an affine field the two coincide (the
  ! moment of a linear term vanishes by definition of the centroid), but
  ! for a field with real curvature the cell average of a degree-p monomial
  ! differs from its centroid value by exactly this moment, an O(h^p) bias
  ! that a plain Taylor expansion would otherwise silently carry into every
  ! reconstructed value. Each Taylor term's own moment correction is
  ! independent of the others (integration is linear, so the mean of the
  ! full polynomial is just the sum of each term's own mean) -- adding the
  ! cubic term never changes the existing quadratic correction. A pure
  ! mesh-geometry quantity, independent of field/reconstruction source (ls
  ! or aho) -- computed once (call from euler_ho_main right after
  ! compute_geometry_mesh) and cached in the module-level cell_moment /
  ! cell_moment3 arrays read by reconstruct(). Quadrature order 2 is exact
  ! for M_jk (a degree-2 integrand); M_jkl needs order 3 (degree-3), so
  ! cell_moment3 is only allocated/filled at order>=4, using that pricier
  ! rule for both moments together rather than allocating two separate
  ! quadratures per cell.
  subroutine compute_cell_moments(mesh)
    type(mesh_type), intent(in) :: mesh

    integer(kind=ENTIER) :: i, k, n_v, quad_order
    real(kind=DOUBLE), dimension(:, :), allocatable :: vcoords, pts
    real(kind=DOUBLE), dimension(:), allocatable :: wts
    real(kind=DOUBLE), dimension(3) :: xc
    integer(kind=ENTIER) :: jj, kk, ll
    logical :: need_third

    need_third = (order >= 4)
    quad_order = merge(3_ENTIER, 2_ENTIER, need_third)

    if (allocated(cell_moment)) deallocate(cell_moment)
    allocate(cell_moment(3, 3, mesh%n_elems))
    cell_moment = 0.0_DOUBLE

    if (allocated(cell_moment3)) deallocate(cell_moment3)
    if (need_third) then
      allocate(cell_moment3(3, 3, 3, mesh%n_elems))
      cell_moment3 = 0.0_DOUBLE
    end if

    do i = 1, mesh%n_elems
      n_v = mesh%elem(i)%n_vert
      if (n_v /= 4 .and. n_v /= 5 .and. n_v /= 6 .and. n_v /= 8) cycle  ! unsupported cell kind: leave moment=0
      allocate(vcoords(3, n_v))
      do k = 1, n_v
        vcoords(:, k) = mesh%vert(mesh%elem(i)%vert(k))%coord
      end do
      call volume_quad_pts(n_v, vcoords, quad_order, pts, wts)
      deallocate(vcoords)

      xc = mesh%elem(i)%coord
      do jj = 1, 3
        do kk = 1, 3
          cell_moment(jj, kk, i) = sum(wts * (pts(jj, :) - xc(jj)) * (pts(kk, :) - xc(kk))) &
            / max(mesh%elem(i)%volume, 1.0e-300_DOUBLE)
        end do
      end do
      if (need_third) then
        do jj = 1, 3
          do kk = 1, 3
            do ll = 1, 3
              cell_moment3(jj, kk, ll, i) = sum(wts * (pts(jj, :) - xc(jj)) &
                * (pts(kk, :) - xc(kk)) * (pts(ll, :) - xc(ll))) &
                / max(mesh%elem(i)%volume, 1.0e-300_DOUBLE)
            end do
          end do
        end do
      end if
      deallocate(pts, wts)
    end do
  end subroutine compute_cell_moments

  ! Polynomial reconstruction of primitive variable w at point xq from cell i.
  ! Hierarchical fallback: if an order-p reconstruction gives unphysical
  ! rho or p, it is replaced by the order-(p-1) result. `third` (the
  ! order-4 cubic term) is optional, allocatable-but-absent-when-unused
  ! exactly like `hess` already is, so every existing call site that only
  ! ever ran at order<=3 compiles and behaves unchanged without passing it.
  function reconstruct(prim, grad, hess, i, xq, xc, g, third) result(w)
    real(kind=DOUBLE), dimension(:, :),          intent(in) :: prim    ! (5, n_elems)
    real(kind=DOUBLE), dimension(:, :, :),       intent(in) :: grad    ! (5, 3, n_elems)
    real(kind=DOUBLE), dimension(:, :, :, :), allocatable, intent(in) :: hess
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE), dimension(3), intent(in) :: xq, xc
    ! Owning cell i's ratio of specific heats -- only used to bound
    ! physical_state's reference sound speed, see gamma_arr's declaration.
    real(kind=DOUBLE), intent(in) :: g
    real(kind=DOUBLE), dimension(:, :, :, :, :), allocatable, intent(in), optional :: third
    real(kind=DOUBLE), dimension(5) :: w

    real(kind=DOUBLE), dimension(3) :: dx
    real(kind=DOUBLE), dimension(5) :: w_try, w_prev
    integer(kind=ENTIER) :: j, k, l

    dx = xq - xc
    w  = prim(:, i)

    if (order >= 2) then
      w_try = w + matmul(grad(:, :, i), dx)
      if (physical_state(w_try, w, g)) w = w_try
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
      if (physical_state(w_try, prim(:, i), g)) w = w_try
    end if

    if (order >= 4 .and. present(third)) then
      if (allocated(third)) then
        if (.not. allocated(cell_moment3)) then
          print *, 'FATAL: reconstruct() called at order>=4 but compute_cell_moments ', &
            'never filled cell_moment3 (order was not yet 4 when it ran, or it was ', &
            'never called) -- the order-4 cubic term would silently carry the wrong ', &
            'cell average. Call compute_cell_moments(mesh) once, at order=4, right ', &
            'after compute_geometry_mesh.'
          error stop 1
        end if
        ! Same pattern one degree up: cell value + gradient + Hessian
        ! (with its own moment correction, unchanged by adding this term --
        ! see compute_cell_moments' header) + cubic term, minus the
        ! cell's own third-moment bias, so this polynomial's cell average
        ! is still exactly prim(:,i). w_prev is the order-3 candidate
        ! from just above (already validated against physical_state); the
        ! order-4 candidate falls back to it, not all the way to order 2,
        ! if the cubic term alone pushes the state unphysical.
        w_prev = w
        w_try = w_prev
        do j = 1, 3
          do k = 1, 3
            do l = 1, 3
              w_try = w_try + (1.0_DOUBLE / 6.0_DOUBLE) * third(:, j, k, l, i) &
                * (dx(j) * dx(k) * dx(l) - cell_moment3(j, k, l, i))
            end do
          end do
        end do
        if (physical_state(w_try, w_prev, g)) w = w_try
      end if
    end if

    ! Enforce positivity (safety net)
    w(1) = max(w(1), 1.0e-12_DOUBLE)
    w(5) = max(w(5), 1.0e-12_DOUBLE)
  end function reconstruct

  ! Ghost cell primitive state for boundary condition
  function ghost_prim(xf, norm, wL, id_bc, t) result(wR)
    real(kind=DOUBLE), dimension(3), intent(in) :: xf, norm
    real(kind=DOUBLE), dimension(5), intent(in) :: wL
    integer(kind=ENTIER), intent(in) :: id_bc
    real(kind=DOUBLE), intent(in) :: t
    real(kind=DOUBLE), dimension(5) :: wR

    real(kind=DOUBLE) :: vn

    if (id_bc < 1 .or. id_bc > n_bc) then
      ! Default: slip wall (mirror normal velocity)
      wR    = wL
      vn    = dot_product(wL(2:4), norm)
      wR(2) = wL(2) - 2.0_DOUBLE * vn * norm(1)
      wR(3) = wL(3) - 2.0_DOUBLE * vn * norm(2)
      wR(4) = wL(4) - 2.0_DOUBLE * vn * norm(3)
      return
    end if

    ! bc_type_id resolved once in read_params -- see its declaration.
    select case (bc_type_id(id_bc))
    case (BC_FREESTREAM)
      wR = bc_val(:, id_bc)
    case (BC_OUTFLOW)
      ! Zero-gradient extrapolation: valid where every characteristic
      ! leaves the domain (locally supersonic outflow).
      wR = wL
    case (BC_DMR_TOP)
      ! Double Mach reflection: exact post/pre-shock state at (xf(1), t),
      ! following the shock's known straight-line motion. See dmr_state.
      wR = dmr_state(xf(1), xf(2), t)
    case default
      ! BC_WALL (also the fallback for blank/unrecognized bc_type,
      ! matching the old string select case's behavior)
      wR    = wL
      vn    = dot_product(wL(2:4), norm)
      wR(2) = wL(2) - 2.0_DOUBLE * vn * norm(1)
      wR(3) = wL(3) - 2.0_DOUBLE * vn * norm(2)
      wR(4) = wL(4) - 2.0_DOUBLE * vn * norm(3)
    end select
  end function ghost_prim

  ! ================================================================
  ! Euler physics
  ! ================================================================

  pure function conserv_to_primit(u, g) result(w)
    real(kind=DOUBLE), dimension(5), intent(in) :: u
    real(kind=DOUBLE), intent(in) :: g
    real(kind=DOUBLE), dimension(5) :: w
    real(kind=DOUBLE) :: rho_inv, ke
    rho_inv = 1.0_DOUBLE / max(u(1), 1.0e-16_DOUBLE)
    w(1) = u(1)
    w(2) = u(2) * rho_inv
    w(3) = u(3) * rho_inv
    w(4) = u(4) * rho_inv
    ke   = 0.5_DOUBLE * (w(2)**2 + w(3)**2 + w(4)**2)
    w(5) = (g - 1.0_DOUBLE) * (u(5) - u(1) * ke)
  end function conserv_to_primit

  pure function primit_to_conserv(w, g) result(u)
    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE), intent(in) :: g
    real(kind=DOUBLE), dimension(5) :: u
    real(kind=DOUBLE) :: ke
    ke   = 0.5_DOUBLE * (w(2)**2 + w(3)**2 + w(4)**2)
    u(1) = w(1)
    u(2) = w(1) * w(2)
    u(3) = w(1) * w(3)
    u(4) = w(1) * w(4)
    u(5) = w(1) * (ke + w(5) / ((g - 1.0_DOUBLE) * w(1)))
  end function primit_to_conserv

  ! Sound speed from primitive state
  pure function cs(w, g) result(c)
    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE), intent(in) :: g
    real(kind=DOUBLE) :: c
    c = sqrt(max(g * w(5) / max(w(1), 1.0e-16_DOUBLE), 0.0_DOUBLE))
  end function cs

  ! Euler flux F(w)·n in conservative variables
  pure function euler_flux(w, n, g) result(F)
    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), intent(in) :: g
    real(kind=DOUBLE), dimension(5) :: F
    real(kind=DOUBLE) :: vn, rho, p, e

    rho = w(1); p = w(5)
    vn  = w(2)*n(1) + w(3)*n(2) + w(4)*n(3)
    e   = p / ((g - 1.0_DOUBLE) * rho) + 0.5_DOUBLE*(w(2)**2+w(3)**2+w(4)**2)

    F(1) = rho * vn
    F(2) = rho * w(2) * vn + p * n(1)
    F(3) = rho * w(3) * vn + p * n(2)
    F(4) = rho * w(4) * vn + p * n(3)
    F(5) = rho * e * vn + p * vn
  end function euler_flux

  ! Rusanov (local Lax-Friedrichs) numerical flux per unit area
  pure function rusanov(wL, wR, n, gL, gR) result(F)
    real(kind=DOUBLE), dimension(5), intent(in) :: wL, wR
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), intent(in) :: gL, gR
    real(kind=DOUBLE), dimension(5) :: F
    real(kind=DOUBLE) :: lam, vnL, vnR

    vnL = wL(2)*n(1) + wL(3)*n(2) + wL(4)*n(3)
    vnR = wR(2)*n(1) + wR(3)*n(2) + wR(4)*n(3)
    lam = max(abs(vnL) + cs(wL, gL), abs(vnR) + cs(wR, gR))

    F = 0.5_DOUBLE * (euler_flux(wL, n, gL) + euler_flux(wR, n, gR)) &
      - 0.5_DOUBLE * lam * (primit_to_conserv(wR, gR) - primit_to_conserv(wL, gL))
  end function rusanov

  ! Three-wave (HLLC-family) approximate Riemann solver: resolves left,
  ! contact and right waves separately instead of Rusanov's single-speed
  ! dissipation, so it captures a contact discontinuity with much less
  ! smearing. Same algorithm as ns_euler_rs_module's three_wave (a
  ! flux-difference form of HLLC -- Toro, "Riemann Solvers and Numerical
  ! Methods for Fluid Dynamics", ch. 10), kept as its own local copy here
  ! (rather than a shared routine with optional gl/gr) so this solver's
  ! gamma_arr multi-material support never touches ns_euler_rs_module or
  ! any other consumer of it.
  pure function three_wave_flux(wL, wR, n, gL, gR) result(F)
    real(kind=DOUBLE), dimension(5), intent(in) :: wL, wR
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), intent(in) :: gL, gR
    real(kind=DOUBLE), dimension(5) :: F

    real(kind=DOUBLE), dimension(5) :: uL, uR, fL, fR, uL_star, uR_star
    real(kind=DOUBLE) :: rhoL, rhoR, vnL, vnR, pL, pR, eL, eR, aL, aR
    real(kind=DOUBLE) :: lambdaL, lambdaR, v_star
    real(kind=DOUBLE) :: rhoL_star, rhoR_star, pL_star, pR_star, sL_wave, sR_wave

    rhoL = wL(1); vnL = dot_product(wL(2:4), n); pL = wL(5)
    uL = primit_to_conserv(wL, gL); eL = uL(5)/rhoL; aL = cs(wL, gL)

    rhoR = wR(1); vnR = dot_product(wR(2:4), n); pR = wR(5)
    uR = primit_to_conserv(wR, gR); eR = uR(5)/rhoR; aR = cs(wR, gR)

    fL(1)   = vnL*uL(1)
    fL(2:4) = vnL*uL(2:4) + pL*n
    fL(5)   = (uL(5) + pL)*vnL

    fR(1)   = vnR*uR(1)
    fR(2:4) = vnR*uR(2:4) + pR*n
    fR(5)   = (uR(5) + pR)*vnR

    lambdaL = max(aL*rhoL, sqrt(rhoL*max(0.0_DOUBLE, pR - pL)), -rhoL*(vnR - vnL))
    lambdaR = max(aR*rhoR, sqrt(rhoR*max(0.0_DOUBLE, pL - pR)), -rhoR*(vnR - vnL))
    v_star  = (lambdaL*vnL + lambdaR*vnR - (pR - pL)) / (lambdaR + lambdaL)

    rhoL_star = 1.0_DOUBLE / (1.0_DOUBLE/rhoL + (v_star - vnL)/lambdaL)
    pL_star   = pL - lambdaL*(v_star - vnL)
    uL_star(1)   = rhoL_star
    uL_star(2:4) = rhoL_star*(wL(2:4) + (v_star - vnL)*n)
    uL_star(5)   = rhoL_star*(eL + (pL*vnL - pL_star*v_star)/lambdaL)

    rhoR_star = 1.0_DOUBLE / (1.0_DOUBLE/rhoR + (vnR - v_star)/lambdaR)
    pR_star   = pR + lambdaR*(v_star - vnR)
    uR_star(1)   = rhoR_star
    uR_star(2:4) = rhoR_star*(wR(2:4) + (v_star - vnR)*n)
    uR_star(5)   = rhoR_star*(eR + (pR_star*v_star - pR*vnR)/lambdaR)

    sL_wave = vnL - lambdaL/rhoL
    sR_wave = vnR + lambdaR/rhoR

    F = 0.5_DOUBLE*(fL + fR) - 0.5_DOUBLE * ( &
      abs(sL_wave)*(uL_star - uL) + &
      abs(v_star)*(uR_star - uL_star) + &
      abs(sR_wave)*(uR - uR_star))
  end function three_wave_flux

  ! Two-wave (HLL-family) approximate Riemann solver: a single intermediate
  ! star state between the left and right acoustic waves (no separate
  ! contact wave), so it is more dissipative than three_wave at a contact
  ! but cheaper. Same algorithm as ns_euler_rs_module's two_wave, kept as
  ! its own local copy for the same reason as three_wave_flux above.
  pure function two_wave_flux(wL, wR, n, gL, gR) result(F)
    real(kind=DOUBLE), dimension(5), intent(in) :: wL, wR
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), intent(in) :: gL, gR
    real(kind=DOUBLE), dimension(5) :: F

    real(kind=DOUBLE), dimension(5) :: uL, uR, fL, fR, u_star
    real(kind=DOUBLE) :: rhoL, rhoR, vnL, vnR, pL, pR, aL, aR
    real(kind=DOUBLE) :: lambdaL, lambdaR, sL_wave, sR_wave

    rhoL = wL(1); vnL = dot_product(wL(2:4), n); pL = wL(5)
    aL = cs(wL, gL); uL = primit_to_conserv(wL, gL)

    rhoR = wR(1); vnR = dot_product(wR(2:4), n); pR = wR(5)
    aR = cs(wR, gR); uR = primit_to_conserv(wR, gR)

    fL(1)   = vnL*uL(1)
    fL(2:4) = vnL*uL(2:4) + pL*n
    fL(5)   = (uL(5) + pL)*vnL

    fR(1)   = vnR*uR(1)
    fR(2:4) = vnR*uR(2:4) + pR*n
    fR(5)   = (uR(5) + pR)*vnR

    lambdaL = max(aL*rhoL, sqrt(rhoL*max(0.0_DOUBLE, pR - pL)), -rhoL*(vnR - vnL))
    lambdaR = max(aR*rhoR, sqrt(rhoR*max(0.0_DOUBLE, pL - pR)), -rhoR*(vnR - vnL))

    sL_wave = vnL - lambdaL/rhoL
    sR_wave = vnR + lambdaR/rhoR

    u_star = (sR_wave*uR - sL_wave*uL - (fR - fL)) / (sR_wave - sL_wave)

    F = 0.5_DOUBLE*(fL + fR) - 0.5_DOUBLE*( &
      abs(sL_wave)*(u_star - uL) + abs(sR_wave)*(uR - u_star))
  end function two_wave_flux

  ! Isentropic vortex primitive state at (x, t)
  pure subroutine vortex_prim(x, t, w)
    real(kind=DOUBLE), dimension(3), intent(in)  :: x
    real(kind=DOUBLE), intent(in) :: t
    real(kind=DOUBLE), dimension(5), intent(out) :: w

    real(kind=DOUBLE) :: beta, r2, tmp, dx, dy

    beta = 5.0_DOUBLE
    dx   = x(1) - u_bg_vortex * t   ! position relative to moving centre
    dy   = x(2) - v_bg_vortex * t
    r2   = dx**2 + dy**2

    tmp  = (gamma_gas - 1.0_DOUBLE) * beta**2 / (8.0_DOUBLE * gamma_gas * PI**2) * exp(1.0_DOUBLE - r2)
    w(1) = (1.0_DOUBLE - tmp)**(1.0_DOUBLE / (gamma_gas - 1.0_DOUBLE))
    w(2) = u_bg_vortex - dy * beta / (2.0_DOUBLE * PI) * exp(0.5_DOUBLE * (1.0_DOUBLE - r2))
    w(3) = v_bg_vortex + dx * beta / (2.0_DOUBLE * PI) * exp(0.5_DOUBLE * (1.0_DOUBLE - r2))
    w(4) = 0.0_DOUBLE
    w(5) = w(1)**gamma_gas
  end subroutine vortex_prim

  ! High-order cell average of the exact vortex field, via genuine volume
  ! quadrature -- works for any cell type volume_quad_pts supports (hex,
  ! tet, prism, pyramid), unlike the z-face-projection trick this
  ! replaced (init_sol's case(2) and compute_error_vortex both used to
  ! find a face with norm ~= +-z and integrate over IT instead, exact
  ! only because a hex/quad cell on the boundary_2d single-layer-extrusion
  ! convention is guaranteed to have one; a tet cell essentially never
  ! does, so that approach silently fell back to a point value at the
  ! centroid there -- a real, if asymptotically vanishing, bias on a
  ! genuinely 3D tet mesh). Degree-3 quadrature is far more than this
  ! smooth (Gaussian-profile) field's own accuracy needs at any mesh size
  ! tested, so it is not the bottleneck at any order; unsupported cell
  ! kinds fall back to the point value at the centroid, same as
  ! compute_cell_moments does for its own moments.
  subroutine vortex_cell_average(mesh, i, t, w)
    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE), intent(in) :: t
    real(kind=DOUBLE), dimension(5), intent(out) :: w

    integer(kind=ENTIER) :: n_v, k, q
    real(kind=DOUBLE), dimension(:, :), allocatable :: vcoords, pts
    real(kind=DOUBLE), dimension(:), allocatable :: wts
    real(kind=DOUBLE), dimension(5) :: w_sum, wq
    real(kind=DOUBLE) :: vol_q

    n_v = mesh%elem(i)%n_vert
    if (n_v /= 4 .and. n_v /= 5 .and. n_v /= 6 .and. n_v /= 8) then
      call vortex_prim(mesh%elem(i)%coord, t, w)
      return
    end if

    allocate(vcoords(3, n_v))
    do k = 1, n_v
      vcoords(:, k) = mesh%vert(mesh%elem(i)%vert(k))%coord
    end do
    call volume_quad_pts(n_v, vcoords, 3_ENTIER, pts, wts)
    deallocate(vcoords)

    w_sum = 0.0_DOUBLE; vol_q = 0.0_DOUBLE
    do q = 1, size(wts)
      call vortex_prim(pts(:, q), t, wq)
      w_sum = w_sum + wts(q) * wq
      vol_q = vol_q + wts(q)
    end do
    if (vol_q > 0.0_DOUBLE) then
      w = w_sum / vol_q
    else
      call vortex_prim(mesh%elem(i)%coord, t, w)
    end if
    deallocate(pts, wts)
  end subroutine vortex_cell_average

  ! Classical double Mach reflection (Woodward & Colella 1984): a Mach-10
  ! shock, inclined 60 deg to the x-axis, initially touching the x-axis at
  ! x=1/6 and moving in +x. Ahead of the shock: gamma=1.4, rho=1.4, u=v=0,
  ! p=1 (a1=1, so the shock speed normal to itself is exactly Ms*a1=10).
  ! Behind the shock (exact 2-state Rankine-Hugoniot values, as used
  ! throughout the DMR literature): rho=8, u=8.25*cos(30deg),
  ! v=-8.25*sin(30deg), p=116.5. The shock's horizontal speed is
  ! Ms*a1/sin(60deg); at time t and height y its x-position is
  ! 1/6 + y/tan(60deg) + that speed * t. Only the domain's top boundary
  ! needs this as a (time-dependent) BC ('dmr_top'); init=3 uses the same
  ! function at t=0 for the initial condition; the bottom boundary is a
  ! fixed post-shock inflow for x<1/6 and a reflecting wall for x>=1/6
  ! (two ordinary BC groups, no time dependence needed there).
  pure function dmr_state(x, y, t) result(w)
    real(kind=DOUBLE), intent(in) :: x, y, t
    real(kind=DOUBLE), dimension(5) :: w

    real(kind=DOUBLE), parameter :: shock_angle = PI/3.0_DOUBLE ! 60 deg
    real(kind=DOUBLE), parameter :: mach_shock = 10.0_DOUBLE
    real(kind=DOUBLE), parameter :: x0 = 1.0_DOUBLE/6.0_DOUBLE
    real(kind=DOUBLE), parameter :: a1 = 1.0_DOUBLE ! pre-shock sound speed (rho=1.4,p=1,gamma=1.4)
    real(kind=DOUBLE) :: shock_speed_x, x_shock

    shock_speed_x = mach_shock*a1/sin(shock_angle)
    x_shock = x0 + y/tan(shock_angle) + shock_speed_x*t

    if (x < x_shock) then
      ! Post-shock
      w = (/ 8.0_DOUBLE, 8.25_DOUBLE*cos(PI/6.0_DOUBLE), -8.25_DOUBLE*sin(PI/6.0_DOUBLE), &
             0.0_DOUBLE, 116.5_DOUBLE /)
    else
      ! Pre-shock (undisturbed)
      w = (/ 1.4_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE /)
    end if
  end function dmr_state

  ! Shu-Osher shock/entropy-wave interaction (Shu & Osher, JCP 1989,
  ! "Efficient implementation of essentially non-oscillatory
  ! shock-capturing schemes, II"): a Mach-3 shock at x=-4 moving into a
  ! sinusoidal density field at rest. Domain [-5,5], gamma=1.4, run to
  ! t=1.8. No closed-form solution exists once the shock crosses the
  ! sine wave (a genuine nonlinear shock/entropy-wave interaction, not a
  ! simple Riemann problem); the standard convention is a very
  ! fine-mesh numerical solution as the "exact" reference.
  pure function shu_osher_state(x) result(w)
    real(kind=DOUBLE), intent(in) :: x
    real(kind=DOUBLE), dimension(5) :: w

    if (x < -4.0_DOUBLE) then
      w = (/ 3.857143_DOUBLE, 2.629369_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 10.33333_DOUBLE /)
    else
      w = (/ 1.0_DOUBLE + 0.2_DOUBLE*sin(5.0_DOUBLE*x), 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE /)
    end if
  end function shu_osher_state

  ! Woodward-Colella interacting blast waves (Woodward & Colella, JCP 1984,
  ! "The numerical simulation of two-dimensional fluid flow with strong
  ! shocks"): domain [0,1], gamma=1.4, density and velocity uniform
  ! (rho=1, u=0), three uniform pressure zones with a 1e5 pressure ratio
  ! at the extremes, reflecting (wall) boundaries at both ends, run to
  ! t=0.038. Two strong shocks form at the pressure jumps, propagate
  ! inward, collide near the domain center, reflect off each other and
  ! off the walls, producing a dense, transient multi-wave structure --
  ! one of the most demanding standard 1D shock-interaction benchmarks,
  ! considerably more severe than Shu-Osher's single Mach-3 shock.
  pure function woodward_colella_state(x) result(w)
    real(kind=DOUBLE), intent(in) :: x
    real(kind=DOUBLE), dimension(5) :: w

    if (x < 0.1_DOUBLE) then
      w = (/ 1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1000.0_DOUBLE /)
    else if (x < 0.9_DOUBLE) then
      w = (/ 1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.01_DOUBLE /)
    else
      w = (/ 1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 100.0_DOUBLE /)
    end if
  end function woodward_colella_state

end module euler_ho_module
