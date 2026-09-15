module silvia_module
  use precision_module
  use mesh_module
  use quadrature_module
  use arbitrary_high_order_module, only: compute_next_order_derivative
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
  real(kind=DOUBLE),  public :: cfl           = 0.4_DOUBLE
  real(kind=DOUBLE),  public :: tmax          = 1.0_DOUBLE
  integer(kind=ENTIER), public :: order       = 1     ! 1, 2, or 3
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
  real(kind=DOUBLE), public :: u_bg_vortex  = 0.0_DOUBLE  ! background x-velocity
  real(kind=DOUBLE), public :: v_bg_vortex  = 0.0_DOUBLE  ! background y-velocity
  integer(kind=ENTIER), public :: init         = 0
  real(kind=DOUBLE), dimension(5), public :: sol_uniform    = [1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE]
  real(kind=DOUBLE), public :: x1drp          = 0.5_DOUBLE
  real(kind=DOUBLE), dimension(5), public :: sol_w_1drp_l   = [1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE]
  real(kind=DOUBLE), dimension(5), public :: sol_w_1drp_r   = [0.125_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.1_DOUBLE]

  ! public :: read_params
  ! public :: init_sol
  ! public :: compute_prim
  ! public :: compute_rhs
  ! public :: compute_dt
  ! public :: compute_error_vortex
  ! public :: test_reconstruction_exactness
  ! public :: compute_cell_moments
  ! Exposed for the standalone quadrature+reconstruction consistency test
  ! (scratchpad/quad_recon_test) -- not otherwise called outside this module.
  !public :: reconstruct, ls_reconstruction, aho_reconstruction
  public :: primit_to_conserv
  public :: conserv_to_primit

  ! Cache of each cell's own second geometric moment about its centroid,
  ! M_jk(i) = (1/V_i) * int_cell (x_j - xc_j)(x_k - xc_k) dV -- a pure
  ! mesh-geometry quantity, independent of order/field, computed once by
  ! compute_cell_moments (called from euler_ho_main right after
  ! compute_geometry_mesh) and read by reconstruct() at order 3. See the
  ! header comment above reconstruct() for why it is needed.
  real(kind=DOUBLE), dimension(:, :, :), allocatable :: cell_moment

contains
  pure function primit_to_conserv(w) result(u)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE), dimension(5) :: u

    u(1) = w(1)
    u(2:4) = w(2:4)*w(1)
    u(5) = w(5)/(gamma_gas - 1) &
      + 0.5_DOUBLE*w(1)*(w(2)**2 + w(3)**2 + w(4)**2)
  end function primit_to_conserv

  pure function conserv_to_primit(u) result(w)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: u
    real(kind=DOUBLE), dimension(5) :: w

    w(1) = u(1)
    w(2:4) = u(2:4)/u(1)
    w(5) = (gamma_gas - 1)*(u(5) &
      - 0.5_DOUBLE*w(1)*(w(2)**2 + w(3)**2 + w(4)**2))
  end function conserv_to_primit
end module silvia_module