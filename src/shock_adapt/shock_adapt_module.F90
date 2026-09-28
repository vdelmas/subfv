module shock_adapt_module
  use precision_module
  use mesh_module
  use shock_adapt_global_data_module, only: gamma, n_bc, bc_kind, bc_val, &
    BC_WALL, BC_FREESTREAM, BC_OUTFLOW, boundary_2d, &
    flux_scheme_id, FLUX_THREE_WAVE, FLUX_MODIFIED_THREE_WAVE, FLUX_TWO_WAVE
  implicit none

  ! Cached by setup_wall_mirror, consumed by compute_rhs at order>=2: cells
  ! touching a wall-only boundary vertex get grad forced to zero (flattened
  ! to first-order there), on top of/regardless of the wall-mirror fit
  ! itself -- user-requested extra safety margin near walls specifically
  ! (not all boundaries), since the wall-mirror fit alone (a real, nonzero
  ! gradient) was not enough to survive the sharp uniform-IC startup
  ! transient on this Mach-3 cylinder-tunnel case.
  logical, dimension(:), allocatable, save :: is_wall_vert_cached

contains

  pure function primit_to_conserv(w) result(u)
    implicit none
    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE), dimension(5) :: u

    u(1) = w(1)
    u(2:4) = w(2:4)*w(1)
    u(5) = w(5)/(gamma - 1.0_DOUBLE) + 0.5_DOUBLE*w(1)*(w(2)**2 + w(3)**2 + w(4)**2)
  end function primit_to_conserv

  pure function conserv_to_primit(u) result(w)
    implicit none
    real(kind=DOUBLE), dimension(5), intent(in) :: u
    real(kind=DOUBLE), dimension(5) :: w

    w(1) = u(1)
    w(2:4) = u(2:4)/u(1)
    w(5) = (gamma - 1.0_DOUBLE)*(u(5) - 0.5_DOUBLE*w(1)*(w(2)**2 + w(3)**2 + w(4)**2))
  end function conserv_to_primit

  pure function sound_speed_w(w) result(a)
    implicit none
    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE) :: a

    a = sqrt(gamma*w(5)/w(1))
  end function sound_speed_w

  ! Ghost primitive state for a boundary face, id_bc = -right_neigh (>=1).
  ! re==0 (no Physical Surface tag) falls back to 'wall' too, same
  ! convention as ale_module.F90's is_wall.
  pure function ghost_state(wL, n, re) result(wR)
    implicit none
    real(kind=DOUBLE), dimension(5), intent(in) :: wL
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    integer(kind=ENTIER), intent(in) :: re
    real(kind=DOUBLE), dimension(5) :: wR

    real(kind=DOUBLE) :: vn
    integer(kind=ENTIER) :: id_bc, kind_

    id_bc = -re
    kind_ = BC_WALL
    if (id_bc >= 1 .and. id_bc <= n_bc) kind_ = bc_kind(id_bc)

    select case (kind_)
    case (BC_FREESTREAM)
      wR = bc_val(:, id_bc)
    case (BC_OUTFLOW)
      wR = wL
    case default ! BC_WALL
      wR = wL
      vn = dot_product(wL(2:4), n)
      wR(2:4) = wL(2:4) - 2.0_DOUBLE*vn*n
    end select
  end function ghost_state

  ! --- 2-state face Riemann solvers (ported from ns_euler_rs_module.F90,
  ! self-contained: use this module's own primit_to_conserv/sound_speed_w,
  ! single constant gamma, no dependency on ns). ---

  subroutine three_wave(sol_w_l, sol_w_r, n, lr_flux, sl, sr)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), dimension(5, 2), intent(inout) :: lr_flux
    real(kind=DOUBLE), intent(inout) :: sl, sr

    real(kind=DOUBLE), dimension(5) :: sol_l, sol_r
    real(kind=DOUBLE) :: rhol, rhor, vn_l, vn_r
    real(kind=DOUBLE) :: pl, pr, rhol_et, rhor_et
    real(kind=DOUBLE) :: v_et, el, er, al, ar
    real(kind=DOUBLE) :: v_bar, pl_bar, pr_bar
    real(kind=DOUBLE) :: lambda_l, lambda_r
    real(kind=DOUBLE), dimension(5) :: fl, fr, sol_l_et, sol_r_et

    rhol = sol_w_l(1)
    vn_l = dot_product(sol_w_l(2:4), n)
    pl = sol_w_l(5)
    sol_l = primit_to_conserv(sol_w_l)
    el = sol_l(5)/rhol
    al = sound_speed_w(sol_w_l)

    rhor = sol_w_r(1)
    vn_r = dot_product(sol_w_r(2:4), n)
    pr = sol_w_r(5)
    sol_r = primit_to_conserv(sol_w_r)
    er = sol_r(5)/rhor
    ar = sound_speed_w(sol_w_r)

    fl(1)   = vn_l*sol_l(1)
    fl(2:4) = vn_l*sol_l(2:4) + pl*n
    fl(5)   = (sol_l(5) + pl)*vn_l

    fr(1)   = vn_r*sol_r(1)
    fr(2:4) = vn_r*sol_r(2:4) + pr*n
    fr(5)   = (sol_r(5) + pr)*vn_r

    lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(vn_r - vn_l))
    lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(vn_r - vn_l))
    v_bar = (lambda_l*vn_l + lambda_r*vn_r - (pr - pl))/(lambda_r + lambda_l)
    v_et = v_bar

    rhol_et = 1.0_DOUBLE/(1.0_DOUBLE/rhol + (v_et - vn_l)/lambda_l)
    pl_bar = pl - lambda_l*(v_et - vn_l)

    sol_l_et(1)   = rhol_et
    sol_l_et(2:4) = rhol_et*(sol_w_l(2:4) + (v_et - vn_l)*n)
    sol_l_et(5)   = rhol_et*(el + (pl*vn_l - pl_bar*v_et)/lambda_l)

    rhor_et = 1.0_DOUBLE/(1.0_DOUBLE/rhor + (vn_r - v_et)/lambda_r)
    pr_bar = pr + lambda_r*(v_et - vn_r)

    sol_r_et(1)   = rhor_et
    sol_r_et(2:4) = rhor_et*(sol_w_r(2:4) + (v_et - vn_r)*n)
    sol_r_et(5)   = rhor_et*(er + (pr_bar*v_et - pr*vn_r)/lambda_r)

    sl = vn_l - lambda_l/rhol
    sr = vn_r + lambda_r/rhor

    lr_flux(:, 1) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(v_et)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et))

    lr_flux(:, 2) = -lr_flux(:, 1)
  end subroutine three_wave

  subroutine modified_three_wave(sol_w_l, sol_w_r, n, lr_flux, sl, sr)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), dimension(5, 2), intent(inout) :: lr_flux
    real(kind=DOUBLE), intent(inout) :: sl, sr

    real(kind=DOUBLE), dimension(5) :: sol_l, sol_r
    real(kind=DOUBLE) :: rhol, rhor, vn_l, vn_r
    real(kind=DOUBLE) :: pl, pr, rhol_et, rhor_et
    real(kind=DOUBLE) :: v_et, el, er, al, ar
    real(kind=DOUBLE) :: v_bar, pl_bar, pr_bar
    real(kind=DOUBLE) :: lambda_l, lambda_r, tang_coeff, tang_energy
    real(kind=DOUBLE), dimension(3) :: vt_l, vt_r, vt_et, dvt
    real(kind=DOUBLE), dimension(5) :: fl, fr, sol_l_et, sol_r_et

    rhol = sol_w_l(1)
    vn_l = dot_product(sol_w_l(2:4), n)
    pl = sol_w_l(5)
    sol_l = primit_to_conserv(sol_w_l)
    el = sol_l(5)/rhol
    al = sound_speed_w(sol_w_l)

    rhor = sol_w_r(1)
    vn_r = dot_product(sol_w_r(2:4), n)
    pr = sol_w_r(5)
    sol_r = primit_to_conserv(sol_w_r)
    er = sol_r(5)/rhor
    ar = sound_speed_w(sol_w_r)

    fl(1)   = vn_l*sol_l(1)
    fl(2:4) = vn_l*sol_l(2:4) + pl*n
    fl(5)   = (sol_l(5) + pl)*vn_l

    fr(1)   = vn_r*sol_r(1)
    fr(2:4) = vn_r*sol_r(2:4) + pr*n
    fr(5)   = (sol_r(5) + pr)*vn_r

    lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(vn_r - vn_l))
    lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(vn_r - vn_l))
    v_bar = (lambda_l*vn_l + lambda_r*vn_r - (pr - pl))/(lambda_r + lambda_l)
    v_et = v_bar

    vt_l = sol_w_l(2:4) - vn_l*n
    vt_r = sol_w_r(2:4) - vn_r*n
    vt_et = (lambda_l*vt_l + lambda_r*vt_r)/(lambda_l + lambda_r)
    dvt = vt_r - vt_l

    tang_coeff = -lambda_l*lambda_r/(lambda_l + lambda_r)
    tang_energy = tang_coeff*dot_product(dvt, vt_et)

    rhol_et = 1.0_DOUBLE/(1.0_DOUBLE/rhol + (v_et - vn_l)/lambda_l)
    pl_bar = pl - lambda_l*(v_et - vn_l)

    sol_l_et(1)   = rhol_et
    sol_l_et(2:4) = rhol_et*(v_et*n + vt_et)
    sol_l_et(5)   = rhol_et*(el + (pl*vn_l - pl_bar*v_et - tang_energy)/lambda_l)

    rhor_et = 1.0_DOUBLE/(1.0_DOUBLE/rhor + (vn_r - v_et)/lambda_r)
    pr_bar = pr + lambda_r*(v_et - vn_r)

    sol_r_et(1)   = rhor_et
    sol_r_et(2:4) = rhor_et*(v_et*n + vt_et)
    sol_r_et(5)   = rhor_et*(er + (pr_bar*v_et - pr*vn_r + tang_energy)/lambda_r)

    sl = vn_l - lambda_l/rhol
    sr = vn_r + lambda_r/rhor

    lr_flux(:, 1) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(v_et)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et))

    lr_flux(:, 2) = -lr_flux(:, 1)
  end subroutine modified_three_wave

  subroutine two_wave(sol_w_l, sol_w_r, n, lr_flux, sl, sr)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), dimension(5, 2), intent(inout) :: lr_flux
    real(kind=DOUBLE), intent(inout) :: sl, sr

    real(kind=DOUBLE), dimension(5) :: sol_l, sol_r
    real(kind=DOUBLE) :: rhol, rhor, vn_l, vn_r, pl, pr, al, ar, lambda_l, lambda_r
    real(kind=DOUBLE), dimension(5) :: fl, fr, sol_et

    rhol = sol_w_l(1)
    vn_l = dot_product(sol_w_l(2:4), n)
    pl = sol_w_l(5)
    al = sound_speed_w(sol_w_l)
    sol_l = primit_to_conserv(sol_w_l)

    rhor = sol_w_r(1)
    vn_r = dot_product(sol_w_r(2:4), n)
    pr = sol_w_r(5)
    ar = sound_speed_w(sol_w_r)
    sol_r = primit_to_conserv(sol_w_r)

    fl(1)   = vn_l*sol_l(1)
    fl(2:4) = vn_l*sol_l(2:4) + pl*n
    fl(5)   = (sol_l(5) + pl)*vn_l

    fr(1)   = vn_r*sol_r(1)
    fr(2:4) = vn_r*sol_r(2:4) + pr*n
    fr(5)   = (sol_r(5) + pr)*vn_r

    lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(vn_r - vn_l))
    lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(vn_r - vn_l))

    sl = vn_l - lambda_l/rhol
    sr = vn_r + lambda_r/rhor

    sol_et = (sr*sol_r - sl*sol_l - (fr - fl))/(sr - sl)

    lr_flux(:, 1) = 0.5_DOUBLE*(fl + fr) &
      - 0.5_DOUBLE*(abs(sl)*(sol_et - sol_l) + abs(sr)*(sol_r - sol_et))
    lr_flux(:, 2) = -lr_flux(:, 1)
  end subroutine two_wave

  ! Cell-centered face-loop RHS: rhs(:,i) accumulates -sum_faces(area*flux),
  ! sum_lambda(:,i) accumulates sum_faces(area*max(|sl|,|sr|)).
  !
  ! order=1: plain cell-average face states (unchanged from the original
  ! scheme). order=2: linear (MUSCL-type) extrapolation to each face using
  ! a Green-Gauss gradient from arbitrary_high_order_module (the same
  ! "aho_gg" reconstruction euler_ho_module.F90 uses at order 2, reused
  ! here as a library -- see the header comment on `order` in
  ! shock_adapt_global_data_module.F90). A reconstructed state with
  ! negative density/pressure falls back to the cell-average there (same
  ! positivity safeguard convention as euler_ho_module.F90's
  ! `physical_state`) -- not a silent fix of an anomaly, standard MUSCL
  ! practice given no limiter is applied on top of the WENO-blended
  ! Green-Gauss gradient.
  ! `order_override`, when present, replaces the module's `order` for this
  ! call only (used by shock_adapt_main.F90's `advance` to ramp: run the
  ! first `order_ramp_iter` Stage-1 iterations at order 1 regardless of the
  ! requested order, since the sharp startup transient from a uniform
  ! initial condition -- before any real shock has formed -- is too steep
  ! for the unlimited order-2 reconstruction: confirmed to blow up to NaN
  ! around iter 165 (t~0.047, well before steady state) without ramping,
  ! even with setup_wall_mirror registered. Standard practice in high-order
  ! CFD startup, not a workaround for a bug in the reconstruction itself.
  subroutine compute_rhs(mesh, sol, rhs, sum_lambda, order_override)
    use shock_adapt_global_data_module, only: order
    use arbitrary_high_order_module, only: compute_next_order_derivative, use_green_gauss
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(out) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(out) :: sum_lambda
    integer(kind=ENTIER), intent(in), optional :: order_override

    integer(kind=ENTIER) :: iface, il, ir, ie, v, dir, order_eff, iv_loc
    real(kind=DOUBLE), dimension(3) :: n, dxL, dxR
    real(kind=DOUBLE) :: area, sl, sr, wave_speed
    real(kind=DOUBLE), dimension(5) :: wL, wR, flux, wL0, wR0
    real(kind=DOUBLE), dimension(5, 2) :: lr_flux
    real(kind=DOUBLE), dimension(:, :), allocatable :: prim, grad_flat
    real(kind=DOUBLE), dimension(:, :, :), allocatable :: grad

    order_eff = order
    if (present(order_override)) order_eff = order_override

    rhs = 0.0_DOUBLE
    sum_lambda = 1e-12_DOUBLE

    if (order_eff >= 2) then
      use_green_gauss = .TRUE.
      allocate(prim(5, mesh%n_elems), grad_flat(15, mesh%n_elems))
      allocate(grad(5, 3, mesh%n_elems))
      do ie = 1, mesh%n_elems
        prim(:, ie) = conserv_to_primit(sol(:, ie))
      end do
      call compute_next_order_derivative(mesh, 3_ENTIER, 5_ENTIER, boundary_2d, &
        prim, grad_flat, deriv_order=1_ENTIER)
      do ie = 1, mesh%n_elems
        do dir = 1, 3
          do v = 1, 5
            grad(v, dir, ie) = grad_flat((dir - 1)*5 + v, ie)
          end do
        end do
      end do

      ! Force grad to zero (flatten to first-order) for any cell touching
      ! a wall vertex -- see the module-level is_wall_vert_cached comment.
      if (allocated(is_wall_vert_cached)) then
        do ie = 1, mesh%n_elems
          do iv_loc = 1, mesh%elem(ie)%n_vert
            if (is_wall_vert_cached(mesh%elem(ie)%vert(iv_loc))) then
              grad(:, :, ie) = 0.0_DOUBLE
              exit
            end if
          end do
        end do
      end if

      ! Positivity kill (ported from euler_ho_module.F90's
      ! apply_positivity_kill, per subfv-de's finding): if the linear
      ! (grad-only) reconstruction predicts negative rho or p at any of a
      ! cell's own VERTICES (not just its face centroids -- a stricter,
      ! more comprehensive check than the per-face fallback below, since it
      ! bounds the extrapolation over the cell's whole vertex set and kills
      ! the gradient consistently for every face of that cell, not just the
      ! one that happened to look bad), zero that cell's grad entirely.
      ! Always on here (no kill_recons opt-in flag): this is the missing
      ! stability ingredient identified via the cross-session debugging
      ! thread on this exact test case, not an optional safety net.
      do ie = 1, mesh%n_elems
        if (mesh%elem(ie)%is_ghost) cycle
        do iv_loc = 1, mesh%elem(ie)%n_vert
          dxL = mesh%vert(mesh%elem(ie)%vert(iv_loc))%coord - mesh%elem(ie)%coord
          if (prim(1, ie) + dot_product(grad(1, :, ie), dxL) <= 0.0_DOUBLE .or. &
            prim(5, ie) + dot_product(grad(5, :, ie), dxL) <= 0.0_DOUBLE) then
            grad(:, :, ie) = 0.0_DOUBLE
            exit
          end if
        end do
      end do
    end if

    do iface = 1, mesh%n_faces
      il = mesh%face(iface)%left_neigh
      ir = mesh%face(iface)%right_neigh
      n = mesh%face(iface)%norm
      area = mesh%face(iface)%area

      wL0 = conserv_to_primit(sol(:, il))
      if (ir > 0) then
        wR0 = conserv_to_primit(sol(:, ir))
      else
        wR0 = ghost_state(wL0, n, ir)
      end if

      wL = wL0; wR = wR0
      if (order_eff >= 2) then
        dxL = mesh%face(iface)%coord - mesh%elem(il)%coord
        wL = wL0 + matmul(grad(:, :, il), dxL)
        ! `x /= x` is the standard finite-value test (false for a normal
        ! number, true for NaN) -- needed because a NaN reconstructed state
        ! (e.g. from a degenerate/boundary-vertex GG fit) silently passes
        ! a "<=0" positivity check (any comparison against NaN is .FALSE.
        ! in Fortran), so a plain positivity guard alone does not catch it.
        ! Confirmed: without this, a NaN state appeared around iter
        ! 200-250 on the very-coarse order-2 run and silently corrupted
        ! the whole solve (dt pinned at its 1e10*cfl sentinel from then on,
        ! since sum_lambda itself went NaN and every subsequent min()
        ! comparison against it was also .FALSE.).
        if (wL(1) <= 0.0_DOUBLE .or. wL(5) <= 0.0_DOUBLE &
          .or. any(wL /= wL)) wL = wL0

        if (ir > 0) then
          dxR = mesh%face(iface)%coord - mesh%elem(ir)%coord
          wR = wR0 + matmul(grad(:, :, ir), dxR)
          if (wR(1) <= 0.0_DOUBLE .or. wR(5) <= 0.0_DOUBLE &
            .or. any(wR /= wR)) wR = wR0
        else
          wR = ghost_state(wL, n, ir)
        end if
      end if

      select case (flux_scheme_id)
      case (FLUX_THREE_WAVE)
        call three_wave(wL, wR, n, lr_flux, sl, sr)
      case (FLUX_TWO_WAVE)
        call two_wave(wL, wR, n, lr_flux, sl, sr)
      case default ! FLUX_MODIFIED_THREE_WAVE
        call modified_three_wave(wL, wR, n, lr_flux, sl, sr)
      end select

      flux = area*lr_flux(:, 1)
      wave_speed = area*max(abs(sl), abs(sr))

      if (.not. mesh%elem(il)%is_ghost) then
        rhs(:, il) = rhs(:, il) - flux
        sum_lambda(il) = sum_lambda(il) + wave_speed
      end if
      if (ir > 0) then
        if (.not. mesh%elem(ir)%is_ghost) then
          rhs(:, ir) = rhs(:, ir) + flux
          sum_lambda(ir) = sum_lambda(ir) + wave_speed
        end if
      end if
    end do
  end subroutine compute_rhs

  ! Registers the aho module's wall-mirror ghost-cell fix (ported from
  ! euler_ho_module.F90's setup_wall_mirror, adapted to shock_adapt's own
  ! BC naming) -- required for order>=2: without it, `wall`-boundary
  ! vertices only get the generic phantom-zero-gradient fallback, which is
  ! not enough to keep this Mach-3 cylinder-tunnel case's cell averages
  ! positivity-preserving (confirmed: order=2 diverged to NaN around
  ! iter 200-250 without this, on the very-coarse mesh, matching the
  ! already-documented `wall-node-reconstruction-todo` gap this exact case
  ! hit before). For every boundary vertex whose touching boundary faces
  ! are ALL wall-type (never mixed with freestream/outflow), average their
  ! area-weighted outward normal into a unit wall normal and register it;
  ! `mirror_wall_vec_start` must also be set >0 to activate the mechanism.
  ! Call once, after compute_geometry_mesh, before the first compute_rhs
  ! call (only matters when order>=2).
  subroutine setup_wall_mirror(mesh)
    use arbitrary_high_order_module, only: set_wall_mirror_data, mirror_wall_vec_start
    implicit none
    type(mesh_type), intent(in) :: mesh

    integer(kind=ENTIER) :: iface, ir, id_bc, iv, v
    real(kind=DOUBLE), dimension(:, :), allocatable :: norm_acc, wall_norm_v
    logical, dimension(:), allocatable :: touches_wall, touches_nonwall, wall_valid_v
    real(kind=DOUBLE), dimension(3) :: n

    allocate(norm_acc(3, mesh%n_vert))
    allocate(touches_wall(mesh%n_vert), touches_nonwall(mesh%n_vert))
    norm_acc = 0.0_DOUBLE
    touches_wall = .FALSE.
    touches_nonwall = .FALSE.

    do iface = 1, mesh%n_faces
      ir = mesh%face(iface)%right_neigh
      if (ir >= 0) cycle
      id_bc = -ir
      if (id_bc < 1 .or. id_bc > n_bc) cycle
      if (.not. allocated(mesh%face(iface)%vert)) cycle
      if (boundary_2d .and. abs(abs(mesh%face(iface)%norm(3)) - 1.0_DOUBLE) < 1.0e-6_DOUBLE) cycle
      do iv = 1, mesh%face(iface)%n_vert
        v = mesh%face(iface)%vert(iv)
        if (bc_kind(id_bc) == BC_WALL) then
          norm_acc(:, v) = norm_acc(:, v) + mesh%face(iface)%area*mesh%face(iface)%norm
          touches_wall(v) = .TRUE.
        else
          touches_nonwall(v) = .TRUE.
        end if
      end do
    end do

    allocate(wall_norm_v(3, mesh%n_vert), wall_valid_v(mesh%n_vert))
    wall_norm_v = 0.0_DOUBLE
    wall_valid_v = .FALSE.
    do v = 1, mesh%n_vert
      if (.not. mesh%vert(v)%is_bound) cycle
      if (.not. touches_wall(v) .or. touches_nonwall(v)) cycle
      n = norm_acc(:, v)
      if (dot_product(n, n) < 1.0e-24_DOUBLE) cycle
      wall_norm_v(:, v) = n/sqrt(dot_product(n, n))
      wall_valid_v(v) = .TRUE.
    end do

    mirror_wall_vec_start = 2_ENTIER
    call set_wall_mirror_data(mesh, wall_norm_v, wall_valid_v)

    if (allocated(is_wall_vert_cached)) deallocate(is_wall_vert_cached)
    allocate(is_wall_vert_cached(mesh%n_vert))
    is_wall_vert_cached = touches_wall

    deallocate(norm_acc, touches_wall, touches_nonwall, wall_norm_v, wall_valid_v)
  end subroutine setup_wall_mirror

  subroutine compute_dt(mesh, sum_lambda, cfl, dt)
    use mpi
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: sum_lambda
    real(kind=DOUBLE), intent(in) :: cfl
    real(kind=DOUBLE), intent(out) :: dt

    integer(kind=ENTIER) :: i, mpi_ierr

    dt = 1e10_DOUBLE
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) dt = min(dt, mesh%elem(i)%volume/sum_lambda(i))
    end do
    dt = cfl*dt

    call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, mpi_ierr)
  end subroutine compute_dt

  subroutine init_sol(mesh, sol)
    use shock_adapt_global_data_module, only: init_uniform, sol_uniform
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: sol

    integer(kind=ENTIER) :: i

    if (init_uniform) then
      do i = 1, mesh%n_elems
        sol(:, i) = primit_to_conserv(sol_uniform)
      end do
    else
      print*, "[-] No init chosen! (only init_uniform is supported by this standalone module)"
      error stop
    end if
  end subroutine init_sol

  ! --- Shock sensor + minimal node movement (Stage 2) ---
  ! Per-face normalized pressure jump, gated by compression, for every
  ! interior in-plane face (abs(norm(3))<1e-8 under boundary_2d, same filter
  ! already used for wall classification in ale_main.F90 -- skips the
  ! z-normal cap faces of the 1-cell extrusion). Optionally ANDs in a
  ! normalized density-jump threshold too (`density_sensor_threshold`):
  ! pressure alone already excludes contact discontinuities/shear layers
  ! via the compression gate (density can jump there with no pressure
  ! jump), so density is an extra, stricter requirement on top of
  ! pressure+compression, not a replacement for it -- makes the detector
  ! more conservative (less sensitive to numerical noise in weak jumps),
  ! at the cost of potentially missing very weak shocks. `sensor` (used for
  ! the vertex diagnostic field) stays the pressure jump either way.
  subroutine compute_shock_sensor(mesh, sol, sensor, flagged)
    use shock_adapt_global_data_module, only: shock_sensor_threshold, &
      use_density_sensor, density_sensor_threshold
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(mesh%n_faces), intent(out) :: sensor
    logical, dimension(mesh%n_faces), intent(out) :: flagged

    integer(kind=ENTIER) :: iface, il, ir
    real(kind=DOUBLE), dimension(5) :: wL, wR
    real(kind=DOUBLE) :: pL, pR, vnL, vnR, rho_jump
    logical :: in_plane, density_ok

    sensor = 0.0_DOUBLE
    flagged = .FALSE.

    do iface = 1, mesh%n_faces
      il = mesh%face(iface)%left_neigh
      ir = mesh%face(iface)%right_neigh
      if (ir <= 0) cycle ! only interior faces are eligible

      in_plane = .TRUE.
      if (boundary_2d) in_plane = abs(mesh%face(iface)%norm(3)) < 1e-8_DOUBLE
      if (.not. in_plane) cycle

      wL = conserv_to_primit(sol(:, il))
      wR = conserv_to_primit(sol(:, ir))
      pL = wL(5); pR = wR(5)
      vnL = dot_product(wL(2:4), mesh%face(iface)%norm)
      vnR = dot_product(wR(2:4), mesh%face(iface)%norm)

      sensor(iface) = abs(pR - pL)/(pR + pL)

      density_ok = .TRUE.
      if (use_density_sensor) then
        rho_jump = abs(wR(1) - wL(1))/(wR(1) + wL(1))
        density_ok = rho_jump > density_sensor_threshold
      end if

      flagged(iface) = (vnR - vnL < 0.0_DOUBLE) .and. (sensor(iface) > shock_sensor_threshold) &
        .and. density_ok
    end do
  end subroutine compute_shock_sensor

  ! Nodal shock detector: for each vertex, the normalized pressure jump is
  ! the (max-min) range of pressure among its neighboring cells (same
  ! thresholds as the face-based compute_shock_sensor, applied directly at
  ! the vertex instead of aggregated from faces).
  !
  ! Compression gate, generalized from the face-based one: among the
  ! vertex's neighbor cells, take the one with p_min ("upstream") and the
  ! one with p_max ("downstream"), n = normalize(x_downstream-x_upstream);
  ! flagged only if (v_downstream-v_upstream).n < 0 (converging along that
  ! direction) -- same physical signature as the 2-cell face gate
  ! (vnR-vnL<0), just picking the two extremal-pressure cells instead of a
  ! face's fixed left/right pair. Without this, a first version flagged a
  ! broad band wrapping the whole cylinder (including the shoulder
  ! expansion fans, which also have a large pressure range but are
  ! divergent, not convergent) -- confirmed visually on a very coarse mesh
  ! before adding this gate.
  subroutine compute_shock_sensor_nodal(mesh, sol, node_sensor, node_flagged)
    use shock_adapt_global_data_module, only: shock_sensor_threshold, &
      use_density_sensor, density_sensor_threshold
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(out) :: node_sensor
    logical, dimension(mesh%n_vert), intent(out) :: node_flagged

    integer(kind=ENTIER) :: iv, j, ide, ide_min, ide_max
    real(kind=DOUBLE), dimension(5) :: w, w_min, w_max
    real(kind=DOUBLE) :: p_min, p_max, rho_min, rho_max, rho_jump
    real(kind=DOUBLE), dimension(3) :: n_dir
    real(kind=DOUBLE) :: n_norm
    logical :: density_ok, compression_ok

    node_sensor = 0.0_DOUBLE
    node_flagged = .FALSE.

    do iv = 1, mesh%n_vert
      p_min = huge(1.0_DOUBLE); p_max = -huge(1.0_DOUBLE)
      rho_min = huge(1.0_DOUBLE); rho_max = -huge(1.0_DOUBLE)
      ide_min = 0; ide_max = 0
      do j = 1, mesh%vert(iv)%n_elems_neigh
        ide = mesh%vert(iv)%elem_neigh(j)
        w = conserv_to_primit(sol(:, ide))
        if (w(5) < p_min) then; p_min = w(5); ide_min = ide; end if
        if (w(5) > p_max) then; p_max = w(5); ide_max = ide; end if
        rho_min = min(rho_min, w(1)); rho_max = max(rho_max, w(1))
      end do
      if (mesh%vert(iv)%n_elems_neigh < 2 .or. ide_min == ide_max) cycle ! nothing to compare

      node_sensor(iv) = (p_max - p_min)/(p_max + p_min)

      n_dir = mesh%elem(ide_max)%coord - mesh%elem(ide_min)%coord
      n_norm = norm2(n_dir)
      compression_ok = .FALSE.
      if (n_norm > 1e-14_DOUBLE) then
        n_dir = n_dir/n_norm
        w_min = conserv_to_primit(sol(:, ide_min))
        w_max = conserv_to_primit(sol(:, ide_max))
        compression_ok = dot_product(w_max(2:4) - w_min(2:4), n_dir) < 0.0_DOUBLE
      end if

      density_ok = .TRUE.
      if (use_density_sensor) then
        rho_jump = (rho_max - rho_min)/(rho_max + rho_min)
        density_ok = rho_jump > density_sensor_threshold
      end if

      node_flagged(iv) = (node_sensor(iv) > shock_sensor_threshold) .and. density_ok &
        .and. compression_ok
    end do
  end subroutine compute_shock_sensor_nodal

  ! Density-gradient-based nodal detector: what a numerical schlieren
  ! actually visualizes (|grad(rho)|), computed directly at the vertex via
  ! a weighted least-squares ("Green-Gauss-like") fit through the
  ! neighboring cell densities -- for cell centroids x_j (relative to the
  ! vertex) and densities rho_j, weight w_j=1/|dx_j|^2 (inverse-distance-
  ! squared, standard WLS gradient weighting), rho_bar = the w-weighted
  ! mean (acting as the reference value at the vertex since rho is only
  ! cell-centered), then solve the 2x2 normal equations for grad=(gx,gy)
  ! minimizing sum_j w_j*[(rho_j-rho_bar) - grad.dx_j]^2. Sensor is
  ! |grad|*local_scale/rho_bar (dimensionless, same normalization spirit as
  ! the other sensors here) so `grad_sensor_threshold` is mesh-independent.
  subroutine compute_shock_sensor_grad(mesh, sol, node_sensor, node_flagged)
    use shock_adapt_global_data_module, only: grad_sensor_threshold
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(out) :: node_sensor
    logical, dimension(mesh%n_vert), intent(out) :: node_flagged

    integer(kind=ENTIER) :: iv, j, ide
    real(kind=DOUBLE), dimension(5) :: w
    real(kind=DOUBLE), dimension(2) :: dx, rhs2
    real(kind=DOUBLE), dimension(2, 2) :: mat2
    real(kind=DOUBLE) :: wj, sum_w, rho_bar, gx, gy, local_scale, grad_mag

    node_sensor = 0.0_DOUBLE
    node_flagged = .FALSE.

    do iv = 1, mesh%n_vert
      if (mesh%vert(iv)%is_bound) cycle ! never flag/move boundary-geometry nodes
      if (mesh%vert(iv)%n_elems_neigh < 3) cycle ! need >=3 for a well-posed 2x2 WLS fit

      sum_w = 0.0_DOUBLE; rho_bar = 0.0_DOUBLE
      do j = 1, mesh%vert(iv)%n_elems_neigh
        ide = mesh%vert(iv)%elem_neigh(j)
        w = conserv_to_primit(sol(:, ide))
        dx = mesh%elem(ide)%coord(1:2) - mesh%vert(iv)%coord(1:2)
        wj = 1.0_DOUBLE/max(dot_product(dx, dx), 1e-30_DOUBLE)
        sum_w = sum_w + wj
        rho_bar = rho_bar + wj*w(1)
      end do
      rho_bar = rho_bar/sum_w

      mat2 = 0.0_DOUBLE; rhs2 = 0.0_DOUBLE
      do j = 1, mesh%vert(iv)%n_elems_neigh
        ide = mesh%vert(iv)%elem_neigh(j)
        w = conserv_to_primit(sol(:, ide))
        dx = mesh%elem(ide)%coord(1:2) - mesh%vert(iv)%coord(1:2)
        wj = 1.0_DOUBLE/max(dot_product(dx, dx), 1e-30_DOUBLE)
        mat2(1, 1) = mat2(1, 1) + wj*dx(1)*dx(1)
        mat2(1, 2) = mat2(1, 2) + wj*dx(1)*dx(2)
        mat2(2, 2) = mat2(2, 2) + wj*dx(2)*dx(2)
        rhs2(1) = rhs2(1) + wj*dx(1)*(w(1) - rho_bar)
        rhs2(2) = rhs2(2) + wj*dx(2)*(w(1) - rho_bar)
      end do
      mat2(2, 1) = mat2(1, 2)

      call solve2(mat2, rhs2, gx, gy)
      grad_mag = sqrt(gx*gx + gy*gy)

      local_scale = huge(1.0_DOUBLE)
      do j = 1, mesh%vert(iv)%n_elems_neigh
        ide = mesh%vert(iv)%elem_neigh(j)
        local_scale = min(local_scale, norm2(mesh%vert(iv)%coord - mesh%elem(ide)%coord))
      end do

      node_sensor(iv) = grad_mag*local_scale/rho_bar
      node_flagged(iv) = node_sensor(iv) > grad_sensor_threshold
    end do
  end subroutine compute_shock_sensor_grad

  ! Topological vertex-vertex adjacency (in-plane mesh edges only). Built
  ! from the z-normal "cap" faces (abs(norm(3))>1-1e-8) -- under the
  ! boundary_2d 1-cell extrusion convention, BOTH the top and bottom caps
  ! are themselves domain-boundary faces (there is only one layer), and
  ! each cap face's vertex list is exactly the true 2D triangle/quad
  ! connectivity: consecutive vertices in a cap face's vert() list are
  ! genuine in-plane mesh edges, with no cross-z-layer contamination (top
  ! and bottom layers get their own, separate adjacency, as they should
  ! since each layer moves independently -- though identically, by mesh
  ! symmetry). Fixed max degree per vertex (16, generous for a 2D tri/quad
  ! mesh); errors loudly if exceeded rather than silently truncating.
  subroutine build_vert_adjacency(mesh, n_neigh, vneigh)
    implicit none
    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), dimension(mesh%n_vert), intent(out) :: n_neigh
    integer(kind=ENTIER), dimension(16, mesh%n_vert), intent(out) :: vneigh

    integer(kind=ENTIER) :: iface, k, n, a, b, kk
    logical :: already

    n_neigh = 0

    do iface = 1, mesh%n_faces
      if (abs(mesh%face(iface)%norm(3)) < 1.0_DOUBLE - 1e-8_DOUBLE) cycle ! not a cap face
      n = mesh%face(iface)%n_vert
      do k = 1, n
        a = mesh%face(iface)%vert(k)
        b = mesh%face(iface)%vert(mod(k, n) + 1)

        already = .FALSE.
        do kk = 1, n_neigh(a)
          if (vneigh(kk, a) == b) already = .TRUE.
        end do
        if (.not. already) then
          n_neigh(a) = n_neigh(a) + 1
          if (n_neigh(a) > 16) then
            print*, "[-] build_vert_adjacency: vertex degree exceeds 16 at vert", a
            error stop
          end if
          vneigh(n_neigh(a), a) = b
        end if

        already = .FALSE.
        do kk = 1, n_neigh(b)
          if (vneigh(kk, b) == a) already = .TRUE.
        end do
        if (.not. already) then
          n_neigh(b) = n_neigh(b) + 1
          if (n_neigh(b) > 16) then
            print*, "[-] build_vert_adjacency: vertex degree exceeds 16 at vert", b
            error stop
          end if
          vneigh(n_neigh(b), b) = a
        end if
      end do
    end do
  end subroutine build_vert_adjacency

  ! Move each flagged node along the local normal of a quadratic curve
  ! fit through its flagged topological neighbors, to reduce the front's
  ! local (discrete) curvature -- a node is predicted from its neighbors'
  ! own quadratic trend, and moved toward that prediction (never using its
  ! own position in the fit, so this is a genuine "flatten toward
  ! neighbors" step, not self-reinforcing).
  !
  ! For node i with flagged neighbors {j}: local 2D frame at the neighbor
  ! centroid, tangent = principal axis of the neighbor point cloud (2x2
  ! eigenproblem, closed-form), normal = perpendicular. Project neighbors
  ! (and i) into (s,h) local coords. Fit h(s)=a*s^2+b*s+c through the
  ! neighbors only (least squares if >3 points, exact if 3, linear-only
  ! i.e. a=0 if only 2). Predicted height at i's own s: h_pred(s_i). Move
  ! i along the normal by (h_pred(s_i) - h_i), capped as usual.
  subroutine compute_node_displacement_curvature(mesh, node_flagged, n_neigh, vneigh, &
      disp, n_moved, max_disp)
    use shock_adapt_global_data_module, only: max_move_frac, curvature_relax
    implicit none
    type(mesh_type), intent(in) :: mesh
    logical, dimension(mesh%n_vert), intent(in) :: node_flagged
    integer(kind=ENTIER), dimension(mesh%n_vert), intent(in) :: n_neigh
    integer(kind=ENTIER), dimension(16, mesh%n_vert), intent(in) :: vneigh
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(out) :: disp
    integer(kind=ENTIER), intent(out) :: n_moved
    real(kind=DOUBLE), intent(out) :: max_disp

    integer(kind=ENTIER) :: iv, j, jn, jn2, k, kk
    integer(kind=ENTIER), dimension(16) :: nbr_id
    real(kind=DOUBLE), dimension(2, 16) :: nbr_xy
    real(kind=DOUBLE), dimension(2) :: centroid, t_hat, n_hat, d, xy_i
    real(kind=DOUBLE), dimension(16) :: s, h
    real(kind=DOUBLE) :: cxx, cxy, cyy, tr, det_, lam1, ex, ey, enorm
    real(kind=DOUBLE) :: s_i, h_i, h_pred, a_q, b_q, c_q, dn, local_scale, cap, mag
    real(kind=DOUBLE), dimension(3, 3) :: mat
    real(kind=DOUBLE), dimension(3) :: rhs_v

    disp = 0.0_DOUBLE
    n_moved = 0
    max_disp = 0.0_DOUBLE

    do iv = 1, mesh%n_vert
      if (.not. node_flagged(iv)) cycle

      ! Collect flagged neighbors up to 2 hops along the flagged-node
      ! subgraph (not just direct topological neighbors): a clean shock
      ! front is topologically a thin CHAIN, so an interior node typically
      ! has exactly 2 direct flagged neighbors (one on each side) -- never
      ! 3 -- so requiring k>=3 from direct neighbors alone would leave
      ! every chain-interior node permanently unmovable (confirmed: an
      ! early version flagged 40 nodes but moved 0 of them). Extending to
      ! each direct flagged neighbor's OWN flagged neighbors too (still
      ! excluding iv itself) gives an interior chain node up to 4 points
      ! (2 on each side), enough for a meaningful quadratic/curvature fit.
      k = 0
      nbr_id = 0
      do j = 1, n_neigh(iv)
        jn = vneigh(j, iv)
        if (.not. node_flagged(jn) .or. k >= 16) cycle
        if (any(nbr_id(1:k) == jn)) cycle
        k = k + 1
        nbr_xy(:, k) = mesh%vert(jn)%coord(1:2)
        nbr_id(k) = jn

        do kk = 1, n_neigh(jn)
          jn2 = vneigh(kk, jn)
          if (jn2 == iv .or. .not. node_flagged(jn2) .or. k >= 16) cycle
          if (any(nbr_id(1:k) == jn2)) cycle
          k = k + 1
          nbr_xy(:, k) = mesh%vert(jn2)%coord(1:2)
          nbr_id(k) = jn2
        end do
      end do
      ! Need >=3 neighbors for a curvature-bearing (quadratic) fit -- a
      ! 2-point line has zero curvature everywhere and gives no reliable
      ! direction to move toward (a plane/line is exactly the degenerate
      ! case that motivated requiring a curvature term in the first place).
      if (k < 3) cycle

      centroid = 0.0_DOUBLE
      do kk = 1, k
        centroid = centroid + nbr_xy(:, kk)
      end do
      centroid = centroid/real(k, DOUBLE)

      ! Principal direction of the neighbor point cloud: 2x2 covariance,
      ! closed-form dominant eigenvector.
      cxx = 0.0_DOUBLE; cxy = 0.0_DOUBLE; cyy = 0.0_DOUBLE
      do kk = 1, k
        d = nbr_xy(:, kk) - centroid
        cxx = cxx + d(1)*d(1); cxy = cxy + d(1)*d(2); cyy = cyy + d(2)*d(2)
      end do
      tr = cxx + cyy
      det_ = cxx*cyy - cxy*cxy
      lam1 = 0.5_DOUBLE*tr + sqrt(max(0.0_DOUBLE, 0.25_DOUBLE*tr*tr - det_))
      if (abs(cxy) > 1e-14_DOUBLE) then
        ex = lam1 - cyy; ey = cxy
      else if (cxx >= cyy) then
        ex = 1.0_DOUBLE; ey = 0.0_DOUBLE
      else
        ex = 0.0_DOUBLE; ey = 1.0_DOUBLE
      end if
      enorm = sqrt(ex*ex + ey*ey)
      if (enorm < 1e-14_DOUBLE) cycle ! degenerate (coincident neighbors)
      t_hat = (/ex, ey/)/enorm
      n_hat = (/-t_hat(2), t_hat(1)/)

      do kk = 1, k
        d = nbr_xy(:, kk) - centroid
        s(kk) = dot_product(d, t_hat)
        h(kk) = dot_product(d, n_hat)
      end do
      xy_i = mesh%vert(iv)%coord(1:2) - centroid
      s_i = dot_product(xy_i, t_hat)
      h_i = dot_product(xy_i, n_hat)

      ! Least-squares (or exact if k==3) quadratic fit through the
      ! neighbors: normal equations for [a,b,c].
      mat = 0.0_DOUBLE; rhs_v = 0.0_DOUBLE
      do kk = 1, k
        mat(1, 1) = mat(1, 1) + s(kk)**4
        mat(1, 2) = mat(1, 2) + s(kk)**3
        mat(1, 3) = mat(1, 3) + s(kk)**2
        mat(2, 3) = mat(2, 3) + s(kk)
        mat(3, 3) = mat(3, 3) + 1.0_DOUBLE
        rhs_v(1) = rhs_v(1) + s(kk)**2*h(kk)
        rhs_v(2) = rhs_v(2) + s(kk)*h(kk)
        rhs_v(3) = rhs_v(3) + h(kk)
      end do
      mat(2, 1) = mat(1, 2); mat(2, 2) = mat(1, 3)
      mat(3, 1) = mat(1, 3); mat(3, 2) = mat(2, 3)
      call solve3(mat, rhs_v, a_q, b_q, c_q)

      h_pred = a_q*s_i*s_i + b_q*s_i + c_q
      dn = curvature_relax*(h_pred - h_i)

      ! Cap by max_move_frac of the local neighbor-centroid distance
      ! (same safety spirit as the snap-based mover).
      local_scale = huge(1.0_DOUBLE)
      do j = 1, mesh%vert(iv)%n_elems_neigh
        local_scale = min(local_scale, &
          norm2(mesh%vert(iv)%coord - mesh%elem(mesh%vert(iv)%elem_neigh(j))%coord))
      end do
      cap = max_move_frac*local_scale
      if (abs(dn) > cap) dn = sign(cap, dn)

      disp(1:2, iv) = dn*n_hat
      disp(3, iv) = 0.0_DOUBLE

      mag = abs(dn)
      if (mag > 1e-14_DOUBLE) then
        n_moved = n_moved + 1
        max_disp = max(max_disp, mag)
      end if
    end do
  end subroutine compute_node_displacement_curvature

  ! Tiny local linear solves (3x3/2x2), no LAPACK dependency needed for
  ! per-vertex systems this small.
  subroutine solve3(mat, rhs, x1, x2, x3)
    implicit none
    real(kind=DOUBLE), dimension(3, 3), intent(in) :: mat
    real(kind=DOUBLE), dimension(3), intent(in) :: rhs
    real(kind=DOUBLE), intent(out) :: x1, x2, x3
    real(kind=DOUBLE) :: det_
    real(kind=DOUBLE), dimension(3, 3) :: m1, m2, m3

    det_ = mat(1, 1)*(mat(2, 2)*mat(3, 3) - mat(2, 3)*mat(3, 2)) &
      - mat(1, 2)*(mat(2, 1)*mat(3, 3) - mat(2, 3)*mat(3, 1)) &
      + mat(1, 3)*(mat(2, 1)*mat(3, 2) - mat(2, 2)*mat(3, 1))

    if (abs(det_) < 1e-30_DOUBLE) then
      x1 = 0.0_DOUBLE; x2 = 0.0_DOUBLE; x3 = 0.0_DOUBLE
      return
    end if

    m1 = mat; m1(:, 1) = rhs
    m2 = mat; m2(:, 2) = rhs
    m3 = mat; m3(:, 3) = rhs

    x1 = (m1(1, 1)*(m1(2, 2)*m1(3, 3) - m1(2, 3)*m1(3, 2)) &
      - m1(1, 2)*(m1(2, 1)*m1(3, 3) - m1(2, 3)*m1(3, 1)) &
      + m1(1, 3)*(m1(2, 1)*m1(3, 2) - m1(2, 2)*m1(3, 1)))/det_
    x2 = (m2(1, 1)*(m2(2, 2)*m2(3, 3) - m2(2, 3)*m2(3, 2)) &
      - m2(1, 2)*(m2(2, 1)*m2(3, 3) - m2(2, 3)*m2(3, 1)) &
      + m2(1, 3)*(m2(2, 1)*m2(3, 2) - m2(2, 2)*m2(3, 1)))/det_
    x3 = (m3(1, 1)*(m3(2, 2)*m3(3, 3) - m3(2, 3)*m3(3, 2)) &
      - m3(1, 2)*(m3(2, 1)*m3(3, 3) - m3(2, 3)*m3(3, 1)) &
      + m3(1, 3)*(m3(2, 1)*m3(3, 2) - m3(2, 2)*m3(3, 1)))/det_
  end subroutine solve3

  subroutine solve2(mat, rhs, x1, x2)
    implicit none
    real(kind=DOUBLE), dimension(2, 2), intent(in) :: mat
    real(kind=DOUBLE), dimension(2), intent(in) :: rhs
    real(kind=DOUBLE), intent(out) :: x1, x2
    real(kind=DOUBLE) :: det_

    det_ = mat(1, 1)*mat(2, 2) - mat(1, 2)*mat(2, 1)
    if (abs(det_) < 1e-30_DOUBLE) then
      x1 = 0.0_DOUBLE; x2 = 0.0_DOUBLE
      return
    end if
    x1 = (rhs(1)*mat(2, 2) - mat(1, 2)*rhs(2))/det_
    x2 = (mat(1, 1)*rhs(2) - rhs(1)*mat(2, 1))/det_
  end subroutine solve2

  ! Bounded per-vertex displacement pulling flagged faces toward the
  ! downstream (compressed) side of the shock, by a controllable fraction
  ! of the inter-centroid distance (`shock_snap_frac`).
  !
  ! An earlier version targeted the point where a *linear interpolation of
  ! pressure* between the two neighbor centroids crossed the arithmetic
  ! mean (pL+pR)/2 -- the paper's Fig. 4b "shortest-distance projection"
  ! idea. That target turned out to be self-defeating in practice: because
  ! the mean-pressure crossing point of a roughly-linear interpolant sits
  ! very close to frac=0.5 by construction, and the face itself already
  ! sits close to the geometric midpoint of its two neighbor centroids on a
  ! fairly regular unstructured mesh, `x_star - face%coord` came out
  ! numerically tiny regardless of `max_move_frac` (confirmed: an mf=0.3
  ! and an mf=0.7 run produced bit-identical output -- the cap was never
  ! the binding constraint, the target itself barely moved). Replaced with
  ! a direct, user-controllable pull of `shock_snap_frac` of the whole
  ! inter-centroid segment, still safety-capped by `max_move_frac` below.
  ! `orig_local_scale`/`cum_disp` are optional and only meant for the
  ! iterative mode (many repeated calls on the same, progressively-moved
  ! mesh): without them, the safety cap is `max_move_frac` of the CURRENT
  ! local neighbor-centroid distance, recomputed fresh each call -- fine
  ! for a single move, but unsafe under repetition. The cap is a MIN over
  ! *all* neighboring cell centroids, including ones off the compression
  ! direction that never shrink, so it does not actually shrink at the same
  ! rate as the specific upstream-downstream gap being compressed call
  ! after call -- confirmed by direct test: repeating the single-shot
  ! formula for a handful of cycles collapsed a cell to negative volume
  ! (by cycle 4 of a 10-cycle smoke test) even with the cap nominally
  ! active. With `orig_local_scale` (each vertex's local scale computed
  ! ONCE before any movement) and `cum_disp` (running total displacement
  ! per vertex, updated in place), the cap instead bounds the CUMULATIVE
  ! displacement against the ORIGINAL geometry -- so a vertex that has
  ! already moved close to `max_move_frac` of its original neighborhood
  ! size simply stops moving, regardless of how many more cycles run.
  subroutine compute_node_displacement(mesh, sol, flagged, disp, n_moved, max_disp, &
      orig_local_scale, cum_disp)
    use shock_adapt_global_data_module, only: max_move_frac, shock_snap_frac
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    logical, dimension(mesh%n_faces), intent(in) :: flagged
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(out) :: disp
    integer(kind=ENTIER), intent(out) :: n_moved
    real(kind=DOUBLE), intent(out) :: max_disp
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in), optional :: orig_local_scale
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout), optional :: cum_disp

    integer(kind=ENTIER) :: iface, il, ir, k, iv, j
    real(kind=DOUBLE) :: local_scale, cap, mag, cum_mag_before, room
    real(kind=DOUBLE), dimension(3) :: dx, new_cum
    integer(kind=ENTIER), dimension(mesh%n_vert) :: wcount
    logical :: cumulative

    cumulative = present(orig_local_scale) .and. present(cum_disp)

    disp = 0.0_DOUBLE
    wcount = 0

    do iface = 1, mesh%n_faces
      if (.not. flagged(iface)) cycle
      il = mesh%face(iface)%left_neigh
      ir = mesh%face(iface)%right_neigh

      dx = shock_snap_frac*(mesh%elem(ir)%coord - mesh%elem(il)%coord)
      dx(3) = 0.0_DOUBLE ! in-plane only (boundary_2d convention)

      do k = 1, mesh%face(iface)%n_vert
        iv = mesh%face(iface)%vert(k)
        disp(:, iv) = disp(:, iv) + dx
        wcount(iv) = wcount(iv) + 1
      end do
    end do

    n_moved = 0
    max_disp = 0.0_DOUBLE
    do iv = 1, mesh%n_vert
      if (wcount(iv) == 0) cycle
      disp(:, iv) = disp(:, iv)/real(wcount(iv), DOUBLE)

      if (cumulative) then
        ! Cap the CUMULATIVE total against the ORIGINAL (pre-movement)
        ! local scale, not the current (possibly already-shrunk) geometry.
        cap = max_move_frac*orig_local_scale(iv)
        cum_mag_before = norm2(cum_disp(:, iv))
        room = cap - cum_mag_before
        if (room <= 0.0_DOUBLE) then
          disp(:, iv) = 0.0_DOUBLE
        else
          new_cum = cum_disp(:, iv) + disp(:, iv)
          mag = norm2(new_cum)
          if (mag > cap .and. mag > 0.0_DOUBLE) then
            ! Scale the *increment* so the new cumulative total lands
            ! exactly on the cap, not the increment itself.
            disp(:, iv) = disp(:, iv)*max(0.0_DOUBLE, (cap - cum_mag_before)/norm2(disp(:, iv)))
          end if
        end if
        cum_disp(:, iv) = cum_disp(:, iv) + disp(:, iv)
      else
        ! Local length scale: min distance from this vertex to any of its
        ! neighboring cell centroids -- always available, a robust proxy
        ! for "local mesh size" without needing an explicit vertex-edge
        ! list. Recomputed fresh each call (safe for a single move only).
        local_scale = huge(1.0_DOUBLE)
        do j = 1, mesh%vert(iv)%n_elems_neigh
          local_scale = min(local_scale, &
            norm2(mesh%vert(iv)%coord - mesh%elem(mesh%vert(iv)%elem_neigh(j))%coord))
        end do
        cap = max_move_frac*local_scale

        mag = norm2(disp(:, iv))
        if (mag > cap .and. mag > 0.0_DOUBLE) disp(:, iv) = disp(:, iv)*(cap/mag)
      end if

      mag = norm2(disp(:, iv))
      if (mag > 1e-14_DOUBLE) then
        n_moved = n_moved + 1
        max_disp = max(max_disp, mag)
      end if
    end do
  end subroutine compute_node_displacement

  ! Per-vertex local length scale (min distance to a neighboring cell
  ! centroid), computed once on the given mesh state -- call this right
  ! after Stage 1 (before any movement) to get `orig_local_scale` for the
  ! iterative mode's cumulative displacement cap.
  subroutine compute_local_scale(mesh, local_scale)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(out) :: local_scale

    integer(kind=ENTIER) :: iv, j

    do iv = 1, mesh%n_vert
      local_scale(iv) = huge(1.0_DOUBLE)
      do j = 1, mesh%vert(iv)%n_elems_neigh
        local_scale(iv) = min(local_scale(iv), &
          norm2(mesh%vert(iv)%coord - mesh%elem(mesh%vert(iv)%elem_neigh(j))%coord))
      end do
    end do
  end subroutine compute_local_scale

  subroutine move_mesh(mesh, disp)
    implicit none
    type(mesh_type), intent(inout) :: mesh
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: disp

    integer(kind=ENTIER) :: i

    do i = 1, mesh%n_vert
      mesh%vert(i)%coord = mesh%vert(i)%coord + disp(:, i)
    end do
  end subroutine move_mesh

  function min_elem_volume(mesh) result(vmin)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE) :: vmin

    integer(kind=ENTIER) :: i

    vmin = huge(1.0_DOUBLE)
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) vmin = min(vmin, mesh%elem(i)%volume)
    end do
  end function min_elem_volume

  ! Vertex-scalar diagnostic field for VTU output: max sensor value over the
  ! faces touching each vertex (0 where no flagged/near-shock face touches
  ! it).
  subroutine compute_vert_sensor(mesh, sensor, vert_sensor)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(mesh%n_faces), intent(in) :: sensor
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(out) :: vert_sensor

    integer(kind=ENTIER) :: iv, j, ifc

    vert_sensor = 0.0_DOUBLE
    do iv = 1, mesh%n_vert
      do j = 1, mesh%vert(iv)%n_faces_neigh
        ifc = mesh%vert(iv)%face_neigh(j)
        vert_sensor(iv) = max(vert_sensor(iv), sensor(ifc))
      end do
    end do
  end subroutine compute_vert_sensor

  ! Ported verbatim from lagrange_module.F90 / ale_module.F90 (identical
  ! routine in both): send/recv vert%coord for every vertex belonging to a
  ! shared/ghost element, so a moved mesh stays consistent across MPI
  ! partition boundaries. Call after every move_mesh, before
  ! compute_geometry_mesh, whenever num_procs>1.
  subroutine mpi_memory_exchange_vert(mesh, mpi_send_recv)
    use mpi
    use mpi_module
    implicit none

    type(mesh_type), intent(inout) :: mesh
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    integer :: mpi_ierr
    integer(kind=ENTIER) :: i, k, j
    integer(kind=ENTIER) :: id_elem, n_vert_tot, id_vert
    logical, dimension(:), allocatable :: vert_averaged

    do i = 1, mpi_send_recv%n_mpi_send_neigh
      n_vert_tot = 0
      do k = 1, mpi_send_recv%mpi_send_neigh(i)%n_elems
        id_elem = mpi_send_recv%mpi_send_neigh(i)%elem_id(k)
        n_vert_tot = n_vert_tot + mesh%elem(id_elem)%n_vert
      end do

      allocate(mpi_send_recv%mpi_send_neigh(i)%sol(3, n_vert_tot))

      n_vert_tot = 1
      do k = 1, mpi_send_recv%mpi_send_neigh(i)%n_elems
        id_elem = mpi_send_recv%mpi_send_neigh(i)%elem_id(k)
        do j = 1, mesh%elem(id_elem)%n_vert
          id_vert = mesh%elem(id_elem)%vert(j)
          mpi_send_recv%mpi_send_neigh(i)%sol(:, n_vert_tot) = mesh%vert(id_vert)%coord
          n_vert_tot = n_vert_tot + 1
        end do
      end do

      call mpi_isend(mpi_send_recv%mpi_send_neigh(i)%sol(1, 1), &
        3*(n_vert_tot - 1), MPI_DOUBLE, &
        mpi_send_recv%mpi_send_neigh(i)%partition_id, &
        0, MPI_COMM_WORLD, mpi_send_recv%mpi_reqsend(i), mpi_ierr)
    end do

    do i = 1, mpi_send_recv%n_mpi_recv_neigh
      n_vert_tot = 0
      do k = 1, mpi_send_recv%mpi_recv_neigh(i)%n_elems
        id_elem = mpi_send_recv%mpi_recv_neigh(i)%elem_id(k)
        n_vert_tot = n_vert_tot + mesh%elem(id_elem)%n_vert
      end do
      allocate(mpi_send_recv%mpi_recv_neigh(i)%sol(3, n_vert_tot))
      call mpi_irecv(mpi_send_recv%mpi_recv_neigh(i)%sol(1, 1), &
        3*n_vert_tot, MPI_DOUBLE, &
        mpi_send_recv%mpi_recv_neigh(i)%partition_id, &
        MPI_ANY_TAG, MPI_COMM_WORLD, mpi_send_recv%mpi_reqrecv(i), mpi_ierr)
    end do

    call mpi_waitall(mpi_send_recv%n_mpi_send_neigh, &
      mpi_send_recv%mpi_reqsend, mpi_send_recv%mpi_sendstat, mpi_ierr)
    call mpi_waitall(mpi_send_recv%n_mpi_recv_neigh, &
      mpi_send_recv%mpi_reqrecv, mpi_send_recv%mpi_recvstat, mpi_ierr)

    ! Unpack: a vertex genuinely owned only as a ghost here (is_ghost=.true.
    ! -- ALL its neighboring elements are ghost elements, see
    ! mesh_connectivity_module.F90's compute_ghost_vert) has no legitimate
    ! local computation to protect, so take the sender's value outright.
    ! A vertex on the actual partition CUT, though, touches at least one
    ! local (owned) element on BOTH sides of the boundary, so
    ! is_ghost=.FALSE. there on every rank that touches it -- each such
    ! rank independently computed its own (generally different, since each
    ! only sees its own side's local neighbor subgraph for the curvature
    ! fit) displacement for the SAME physical vertex. Confirmed by the
    ! user directly inspecting a cluster run: shared boundary nodes ended
    ! up at different positions on each side. Average with the received
    ! value there instead of overwriting, so both ranks converge on the
    ! same final position (order-independent, symmetric).
    !
    ! A vertex on the cut is typically referenced by SEVERAL ghost elements
    ! from the same neighbor rank (and possibly, at a multi-rank corner, by
    ! ghost elements from more than one neighbor), so it appears more than
    ! once across these nested loops. The averaging update is NOT
    ! idempotent (applying it twice in place biases the result toward the
    ! received value: 0.5, then 0.75, then 0.875, ...), and the occurrence
    ! count generally differs between the two ranks sharing a vertex, so
    ! repeating it silently broke the intended 50/50 symmetry -- this was
    ! the actual cause of a visible mesh crack at partition cuts after a
    ! second movement cycle. Guard so each vertex is only averaged once per
    ! call, regardless of how many times it is encountered.
    allocate(vert_averaged(size(mesh%vert)))
    vert_averaged = .false.
    do i = 1, mpi_send_recv%n_mpi_recv_neigh
      n_vert_tot = 1
      do k = 1, mpi_send_recv%mpi_recv_neigh(i)%n_elems
        id_elem = mpi_send_recv%mpi_recv_neigh(i)%elem_id(k)
        do j = 1, mesh%elem(id_elem)%n_vert
          id_vert = mesh%elem(id_elem)%vert(j)
          if (.not. vert_averaged(id_vert)) then
            if (mesh%vert(id_vert)%is_ghost) then
              mesh%vert(id_vert)%coord = mpi_send_recv%mpi_recv_neigh(i)%sol(:, n_vert_tot)
            else
              mesh%vert(id_vert)%coord = 0.5_DOUBLE*(mesh%vert(id_vert)%coord &
                + mpi_send_recv%mpi_recv_neigh(i)%sol(:, n_vert_tot))
            end if
            vert_averaged(id_vert) = .true.
          end if
          n_vert_tot = n_vert_tot + 1
        end do
      end do
    end do
    deallocate(vert_averaged)

    do i = 1, mpi_send_recv%n_mpi_send_neigh
      deallocate(mpi_send_recv%mpi_send_neigh(i)%sol)
    end do
    do i = 1, mpi_send_recv%n_mpi_recv_neigh
      deallocate(mpi_send_recv%mpi_recv_neigh(i)%sol)
    end do
  end subroutine mpi_memory_exchange_vert

end module shock_adapt_module
