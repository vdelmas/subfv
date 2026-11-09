module ns_euler_module
  use precision_module
  use mesh_module
  use ns_euler_primitives_module
  use ns_euler_recon_module
  use ns_euler_zb_module
  use ns_euler_exact_sol_module
  use ns_euler_rs_module
  implicit none

  ! multi_point_eps diag, per cell, rebuilt at every flux call: number of its nodes with
  ! eps_p <= 0, and number of its boundary_2d z-face sides with eps* < 0 (printed on the
  ! [omega] culprit line)
  integer, allocatable, dimension(:) :: eps_cell_neg, eps_cell_zf
  ! ener_dump_sides: unit of the per-rank side dump, open during the first flux call only
  ! (newunit returns NEGATIVE unit numbers in gfortran: track the open state explicitly)
  integer :: ener_dump_unit = -1
  logical :: ener_dumping = .false.
  integer(kind=8) :: ener_calls = 0_8
  real(kind=DOUBLE), dimension(10), parameter :: ener_edges = (/ -100.0_DOUBLE, -10.0_DOUBLE, &
    -1.0_DOUBLE, 0.0_DOUBLE, 0.1_DOUBLE, 0.5_DOUBLE, 0.9_DOUBLE, 1.1_DOUBLE, 2.0_DOUBLE, 10.0_DOUBLE /)

contains
  subroutine init_sol(mesh, sol, me, num_procs, mpi_send_recv)
    use mpi
    use ns_global_data_module
    use mpi_module
    use ns_io_module, only: restart_from_vtu_file
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), &
      intent(inout) :: sol
    integer(kind=ENTIER), intent(in) :: me, num_procs
    type(mpi_send_recv_type) :: mpi_send_recv

    integer(kind=ENTIER) :: i, mpi_ierr
    real(kind=DOUBLE), dimension(5) :: w
    real(kind=DOUBLE) :: M_helm, r_helm, delta_helm, H_helm
    real(kind=DOUBLE), dimension(5) :: sol_left, sol_right
    real(kind=DOUBLE) :: rmin, rmin_share
    integer(kind=ENTIER) :: imin

    if (init_uniform) then
      do i = 1, mesh%n_elems
        sol(:, i) = primit_to_conserv(sol_uniform)
      end do
    else if (init_1drp) then
      do i = 1, mesh%n_elems
        if( mesh%elem(i)%coord(1) < x1drp ) then
          sol(:, i) = primit_to_conserv(sol_w_1drp_l)
        else 
          sol(:, i) = primit_to_conserv(sol_w_1drp_r)
        end if
      end do
    else if (init_isentropic_vortex) then
      do i = 1, mesh%n_elems
        call sol_isentropic_vortex(mesh%elem(i)%coord, w, 0.0_DOUBLE)
        sol(:, i) = primit_to_conserv(w)
      end do
    else if (init_gresho) then
      do i = 1, mesh%n_elems
        call sol_gresho_mach(mesh%elem(i)%coord, w, 1.0e-5_DOUBLE)
        sol(:, i) = primit_to_conserv(w)
      end do
    else if (init_potential_flow_2d) then
      do i = 1, mesh%n_elems
        call sol_potential_flow_2d(mesh%elem(i)%coord, w)
        sol(:, i) = primit_to_conserv(w)
      end do
    else if (init_potential_flow_3d) then
      do i = 1, mesh%n_elems
        call sol_potential_flow_3d(mesh%elem(i)%coord, w)
        sol(:, i) = primit_to_conserv(w)
      end do
    else if (init_sedov) then
      rmin = huge(1.0_DOUBLE)
      imin = -1
      do i=1, mesh%n_elems
        if( boundary_2d ) then
          if( norm2(mesh%elem(i)%coord(:2)) < rmin ) then
            rmin = norm2(mesh%elem(i)%coord(:2))
            imin = i
          end if
        else
          if( norm2(mesh%elem(i)%coord) < rmin ) then
            rmin = norm2(mesh%elem(i)%coord)
            imin = i
          end if
        end if
      end do

      if (num_procs > 1) then
        rmin_share = rmin
        call MPI_ALLREDUCE(MPI_IN_PLACE, rmin_share, 1, MPI_DOUBLE, &
          MPI_MIN, MPI_COMM_WORLD, mpi_ierr)
        if( abs(rmin - rmin_share) > 1e-14_DOUBLE ) then
          imin = -1
        end if
      end if

      do i = 1, mesh%n_elems
        sol(1, i) = 1.0_DOUBLE
        sol(2, i) = 0.0_DOUBLE
        sol(3, i) = 0.0_DOUBLE
        sol(4, i) = 0.0_DOUBLE
        sol(5, i) = 1e-12_DOUBLE/(gamma - 1.0_DOUBLE)
      end do

      if( imin > 0 ) then
        if( boundary_2d ) then
          sol(1, imin) = 1.0_DOUBLE
          sol(2, imin) = 0.0_DOUBLE
          sol(3, imin) = 0.0_DOUBLE
          sol(4, imin) = 0.0_DOUBLE
          sol(5, imin) = 0.984042_DOUBLE/mesh%elem(imin)%volume
        else
          sol(1, imin) = 1.0_DOUBLE
          sol(2, imin) = 0.0_DOUBLE
          sol(3, imin) = 0.0_DOUBLE
          sol(4, imin) = 0.0_DOUBLE
          sol(5, imin) = 0.851072_DOUBLE/mesh%elem(imin)%volume
        end if
      end if
    else if (init_kelvin) then
      M_helm   = 1e-2_DOUBLE
      r_helm   = M_helm
      delta_helm = 0.1_DOUBLE
      tmax = 0.8_DOUBLE / M_helm
      print *, "T_max =", tmax
      do i = 1, mesh%n_elems
        H_helm = kelvin_helmholtz(mesh%elem(i)%coord(2))
        w(1) = gamma + r_helm * H_helm
        w(2) = M_helm * H_helm
        w(3) = delta_helm * M_helm * sin(2.0_DOUBLE*pi*mesh%elem(i)%coord(1))
        w(4) = 0.0_DOUBLE
        w(5) = 1.0_DOUBLE
        sol(:, i) = primit_to_conserv(w)
      end do
    else if (init_double_mach) then
      sol_left  = (/1.4_DOUBLE, 0.0_DOUBLE,    0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE/)
      sol_right = (/8.0_DOUBLE, 7.145_DOUBLE, -4.125_DOUBLE, 0.0_DOUBLE, 116.5_DOUBLE/)
      tmax = 0.2_DOUBLE
      print *, "T_max =", tmax
      do i = 1, mesh%n_elems
        if (mesh%elem(i)%coord(2) < 1.732_DOUBLE*(mesh%elem(i)%coord(1) - 0.1667_DOUBLE)) then
          sol(:, i) = primit_to_conserv(sol_left)
        else
          sol(:, i) = primit_to_conserv(sol_right)
        end if
      end do
    else if (init_triple_point) then
      ! Triple-point shock interaction (single material, gamma = 1.4).
      ! Domain [0,7] x [0,3], fluid at rest, reflecting walls everywhere.
      !   x <= 1            : rho = 1,     p = 1
      !   x >  1, y <= 1.5  : rho = 1,     p = 0.1
      !   x >  1, y >  1.5  : rho = 0.125, p = 0.1
      w(2:4) = 0.0_DOUBLE
      do i = 1, mesh%n_elems
        if (mesh%elem(i)%coord(1) <= 1.0_DOUBLE) then
          w(1) = 1.0_DOUBLE
          w(5) = 1.0_DOUBLE
        else if (mesh%elem(i)%coord(2) <= 1.5_DOUBLE) then
          w(1) = 1.0_DOUBLE
          w(5) = 0.1_DOUBLE
        else
          w(1) = 0.125_DOUBLE
          w(5) = 0.1_DOUBLE
        end if
        sol(:, i) = primit_to_conserv(w)
      end do
    else if (init_restart) then
      call restart_from_vtu_file(mesh, sol, restart_file, &
        me, num_procs)
    else
      print*, "[-] No init chosen!"
    end if

    if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, &
      mesh%n_elems, 5, sol)
  end subroutine init_sol

  subroutine compute_euler_flux(mesh, sol, rhs, &
      grad, sum_lambda, second_order, vp, h_p, mat_h_p)
    use mpi
    use ns_global_data_module, only: zb_lag_mode, zb_adv_mode, zb_eps_mode, &
      bc_style, scheme, exclude_bound_vert, ener_dump_sides, &
      scheme_id, scheme_adv_id, scheme_lag_id, &
      SCHEME_ZB, SCHEME_WIP, SCHEME_WIP2, SCHEME_WIP2_NOLM, SCHEME_WIP2_AADV, SCHEME_WIP2_HYB, &
      SCHEME_WIP2_EJUMP, SCHEME_WIP2_EDIV, SCHEME_WIP2_EMAX, &
      SCHEME_WIP2_ENTH, SCHEME_WIP2_ENTH_EJUMP, &
      SCHEME_WIP2_TP, SCHEME_WIP2_TP_ENTH, SCHEME_WIP2_EPOS, SCHEME_WIP2_EPOS2, &
      SCHEME_WIP2_ED1, SCHEME_WIP2_ED2, &
      SCHEME_WIP2_MD, SCHEME_WIP2_MD5, SCHEME_WIP2_MDH, SCHEME_USI3D, &
      SCHEME_WIP2_PPEP, SCHEME_WIP2_PPEPI, SCHEME_WIP2_PPEPV, SCHEME_WIP2_PPEPVI, &
      SCHEME_MULTI_POINT_ENER, SCHEME_MULTI_POINT_ENER_ISO, &
      SCHEME_MULTI_POINT_ENER_V, SCHEME_MULTI_POINT_ENER_V_ISO, &
      SCHEME_MULTI_POINT_EPS, SCHEME_MULTI_POINT_EPS_ISO, &
      SCHEME_MULTI_POINT_EPS2, SCHEME_MULTI_POINT_EPS2M, SCHEME_MULTI_POINT_PS, SCHEME_MULTI_POINT_PSW, SCHEME_MULTI_POINT_PV, &
      SCHEME_ADV_AR1D, SCHEME_ADV_AM, SCHEME_ADV_AMISO, &
      SCHEME_ADV_ARMD, SCHEME_ADV_ARMDU, &
      SCHEME_ADV_ARMDMAT, SCHEME_ADV_ARMDUMAT, &
      SCHEME_LAG_LS, SCHEME_LAG_LSU, SCHEME_ADV_ARMDWIP, &
      SCHEME_LAG_LSWIP, SCHEME_LAG_LPP, SCHEME_LAG_LVPPP, SCHEME_LAG_LS1D, SCHEME_LAG_LPF
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(3, mesh%n_elems), intent(inout) :: vp
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(inout) :: sum_lambda
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: h_p
    real(kind=DOUBLE), dimension(3, 3, mesh%n_vert), intent(in) :: mat_h_p
    logical, intent(in) :: second_order

    integer(kind=ENTIER) :: j, idse, ide, id_vert, nsfn, nsen, max_nsen
    integer(kind=ENTIER) :: adv_mode, eps_mode, ep_mode
    integer(kind=ENTIER) :: lag_mode
    real(kind=DOUBLE), dimension(3) :: v_vert
    real(kind=DOUBLE), dimension(:), allocatable :: sum_lambda_vert
    real(kind=DOUBLE), dimension(:, :), allocatable :: flux_sum_vert

    rhs = 0.0_DOUBLE
    v_vert = 0.0_DOUBLE
    sum_lambda = 1e-12_DOUBLE

    max_nsen = 0
    do id_vert = 1, mesh%n_vert
      max_nsen = max(max_nsen, mesh%vert(id_vert)%n_sub_elems_neigh)
    end do
    allocate(sum_lambda_vert(max_nsen))
    allocate(flux_sum_vert(5, max_nsen))

    if (scheme_id == SCHEME_ZB) then
      do id_vert = 1, mesh%n_vert
        nsfn = mesh%vert(id_vert)%n_sub_faces_neigh
        nsen = mesh%vert(id_vert)%n_sub_elems_neigh
        sum_lambda_vert(1:nsen) = 0.0_DOUBLE
        flux_sum_vert(:, 1:nsen) = 0.0_DOUBLE

        select case (scheme_adv_id)
        case (SCHEME_ADV_AR1D)
          call compute_rhs_around_vert_AR1D(mesh, sol, grad, &
            nsen, sum_lambda_vert, &
            flux_sum_vert, second_order, id_vert)
        case (SCHEME_ADV_AM)
          call compute_rhs_around_vert_AM(mesh, sol, grad, &
            nsen, v_vert, sum_lambda_vert, &
            flux_sum_vert, second_order, id_vert, h_p)
        case (SCHEME_ADV_AMISO)
          call compute_rhs_around_vert_AMISO(mesh, sol, grad, &
            nsen, v_vert, sum_lambda_vert, &
            flux_sum_vert, second_order, id_vert, h_p)
        case (SCHEME_ADV_ARMD)
          call compute_rhs_around_vert_ARMD(mesh, sol, grad, &
            nsen, sum_lambda_vert, &
            flux_sum_vert, second_order, id_vert, vp, h_p)
        case (SCHEME_ADV_ARMDWIP)
          call compute_rhs_around_vert_ARMDWIP(mesh, sol, grad, &
            nsen, sum_lambda_vert, &
            flux_sum_vert, second_order, id_vert, vp, h_p, mat_h_p)
        case (SCHEME_ADV_ARMDU)
          call compute_rhs_around_vert_ARMDU(mesh, sol, grad, &
            nsen, sum_lambda_vert, &
            flux_sum_vert, second_order, id_vert, vp)
        case (SCHEME_ADV_ARMDMAT)
          call compute_rhs_around_vert_ARMDMAT(mesh, sol, grad, &
            nsen, sum_lambda_vert, &
            flux_sum_vert, second_order, id_vert, vp, h_p)
        case (SCHEME_ADV_ARMDUMAT)
          call compute_rhs_around_vert_ARMDUMAT(mesh, sol, grad, &
            nsen, sum_lambda_vert, &
            flux_sum_vert, second_order, id_vert, vp)
        case default
          print*,"Unknown ZB advection"
          error stop
        end select

        select case (scheme_lag_id)
        case (SCHEME_LAG_LS)
          call compute_rhs_around_vert_LS(mesh, sol, grad, &
            nsen, flux_sum_vert, &
            id_vert, vp, h_p, second_order)
        case (SCHEME_LAG_LSWIP)
          call compute_rhs_around_vert_LSWIP(mesh, sol, grad, &
            nsen, flux_sum_vert, &
            id_vert, vp, h_p, mat_h_p, second_order)
        case (SCHEME_LAG_LSU)
          call compute_rhs_around_vert_LSU(mesh, sol, grad, &
            nsen, flux_sum_vert, &
            id_vert, vp, second_order)
        case (SCHEME_LAG_LPP)
          call compute_rhs_around_vert_LPP(mesh, sol, grad, &
            nsen, flux_sum_vert, id_vert, second_order)
        case (SCHEME_LAG_LVPPP)
          call compute_rhs_around_vert_LVPPP(mesh, sol, grad, &
            nsen, flux_sum_vert, id_vert, vp, h_p, second_order)
        case (SCHEME_LAG_LS1D)
          call compute_rhs_around_vert_LS1D(mesh, sol, grad, &
            nsen, flux_sum_vert, id_vert, second_order)
        case (SCHEME_LAG_LPF)
          call compute_rhs_around_vert_LPF(mesh, sol, grad, &
            nsen, flux_sum_vert, id_vert, second_order)
        case default
          print*,"Unknown ZB Lagrange"
          error stop
        end select

        do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
          idse = mesh%vert(id_vert)%sub_elem_neigh(j)
          ide = mesh%sub_elem(idse)%mesh_elem
          rhs(:, ide) = rhs(:, ide) - flux_sum_vert(:, j)
          sum_lambda(ide) = sum_lambda(ide) + sum_lambda_vert(j)
        end do
      end do
    else if (scheme_id == SCHEME_WIP) then
      do id_vert = 1, mesh%n_vert
        nsfn = mesh%vert(id_vert)%n_sub_faces_neigh
        nsen = mesh%vert(id_vert)%n_sub_elems_neigh
        sum_lambda_vert(1:nsen) = 0.0_DOUBLE
        flux_sum_vert(:, 1:nsen) = 0.0_DOUBLE
        call compute_rhs_around_vert_WIP(mesh, sol, grad, &
          nsen, flux_sum_vert, sum_lambda_vert, &
          id_vert, vp, h_p, second_order)

        do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
          idse = mesh%vert(id_vert)%sub_elem_neigh(j)
          ide = mesh%sub_elem(idse)%mesh_elem
          rhs(:, ide) = rhs(:, ide) - flux_sum_vert(:, j)
          sum_lambda(ide) = sum_lambda(ide) + sum_lambda_vert(j)
        end do
      end do
    else if (scheme_id == SCHEME_WIP2 .or. scheme_id == SCHEME_WIP2_NOLM &
        .or. scheme_id == SCHEME_WIP2_AADV .or. scheme_id == SCHEME_WIP2_HYB &
        .or. scheme_id == SCHEME_WIP2_EJUMP .or. scheme_id == SCHEME_WIP2_EDIV &
        .or. scheme_id == SCHEME_WIP2_EMAX .or. scheme_id == SCHEME_WIP2_ENTH &
        .or. scheme_id == SCHEME_WIP2_ENTH_EJUMP .or. scheme_id == SCHEME_WIP2_TP &
        .or. scheme_id == SCHEME_WIP2_TP_ENTH .or. scheme_id == SCHEME_WIP2_EPOS .or. scheme_id == SCHEME_WIP2_EPOS2 &
        .or. scheme_id == SCHEME_WIP2_ED1 .or. scheme_id == SCHEME_WIP2_ED2 &
        .or. scheme_id == SCHEME_WIP2_MD .or. scheme_id == SCHEME_WIP2_MD5 &
        .or. scheme_id == SCHEME_WIP2_MDH .or. scheme_id == SCHEME_WIP2_PPEP &
        .or. scheme_id == SCHEME_WIP2_PPEPI &
        .or. scheme_id == SCHEME_WIP2_PPEPV .or. scheme_id == SCHEME_WIP2_PPEPVI) then
      do id_vert = 1, mesh%n_vert
        nsfn = mesh%vert(id_vert)%n_sub_faces_neigh
        nsen = mesh%vert(id_vert)%n_sub_elems_neigh
        sum_lambda_vert(1:nsen) = 0.0_DOUBLE
        flux_sum_vert(:, 1:nsen) = 0.0_DOUBLE
        if (scheme_id == SCHEME_WIP2_AADV) then
          adv_mode = 1
        else if (scheme_id == SCHEME_WIP2_HYB) then
          adv_mode = 2
        else if (scheme_id == SCHEME_WIP2_MD) then
          adv_mode = 3
        else if (scheme_id == SCHEME_WIP2_MD5) then
          adv_mode = 4
        else if (scheme_id == SCHEME_WIP2_MDH) then
          adv_mode = 5
        else
          adv_mode = 0
        end if
        if (scheme_id == SCHEME_WIP2_EJUMP .or. scheme_id == SCHEME_WIP2_ENTH_EJUMP) then
          eps_mode = 1
        else if (scheme_id == SCHEME_WIP2_EDIV) then
          eps_mode = 2
        else if (scheme_id == SCHEME_WIP2_EMAX) then
          eps_mode = 3
        else if (scheme_id == SCHEME_WIP2_EPOS) then
          eps_mode = 4
        else if (scheme_id == SCHEME_WIP2_EPOS2) then
          eps_mode = 5
        else if (scheme_id == SCHEME_WIP2_ED1) then
          eps_mode = 6
        else if (scheme_id == SCHEME_WIP2_ED2) then
          eps_mode = 7
        else
          eps_mode = 0
        end if
        ! La namelist l emporte sur le scheme_id, pour balayer les capteurs sans
        ! creer une constante SCHEME_ par combinaison.
        if (zb_eps_mode >= 0) eps_mode = zb_eps_mode
        ! lag_mode : forme du flux d energie de la partie Lagrange. Aucun
        ! scheme_id ne le fixe pour l instant ; il se pilote par la namelist
        ! zb_lag_mode, ce qui evite de creer une constante SCHEME_ par
        ! combinaison lag_mode x adv_mode (il y en aurait 24).
        lag_mode = 0
        if (zb_lag_mode >= 0) lag_mode = zb_lag_mode
        ! Meme chose pour adv_mode : la namelist l emporte sur le scheme_id,
        ! pour pouvoir balayer les flux d advection sans toucher au code.
        if (zb_adv_mode >= 0) adv_mode = zb_adv_mode

        if (scheme_id == SCHEME_WIP2_PPEP) then
          ep_mode = 1
        else if (scheme_id == SCHEME_WIP2_PPEPI) then
          ep_mode = 2
        else if (scheme_id == SCHEME_WIP2_PPEPV) then
          ep_mode = 3
        else if (scheme_id == SCHEME_WIP2_PPEPVI) then
          ep_mode = 4
        else
          ep_mode = 0
        end if
        call compute_rhs_around_vert_WIP2(mesh, sol, grad, &
          nsen, flux_sum_vert, sum_lambda_vert, &
          id_vert, second_order, scheme_id == SCHEME_WIP2, adv_mode, eps_mode, &
          scheme_id == SCHEME_WIP2_ENTH .or. scheme_id == SCHEME_WIP2_ENTH_EJUMP &
          .or. scheme_id == SCHEME_WIP2_TP_ENTH, &
          scheme_id == SCHEME_WIP2_TP .or. scheme_id == SCHEME_WIP2_TP_ENTH, &
          ep_mode, lag_mode)

        do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
          idse = mesh%vert(id_vert)%sub_elem_neigh(j)
          ide = mesh%sub_elem(idse)%mesh_elem
          rhs(:, ide) = rhs(:, ide) - flux_sum_vert(:, j)
          sum_lambda(ide) = sum_lambda(ide) + sum_lambda_vert(j)
        end do
      end do
    else if (scheme_id == SCHEME_USI3D) then
      do id_vert = 1, mesh%n_vert
        nsfn = mesh%vert(id_vert)%n_sub_faces_neigh
        nsen = mesh%vert(id_vert)%n_sub_elems_neigh
        sum_lambda_vert(1:nsen) = 0.0_DOUBLE
        flux_sum_vert(:, 1:nsen) = 0.0_DOUBLE
        call compute_rhs_around_vert_usi3d(mesh, sol, grad, &
          nsen, flux_sum_vert, sum_lambda_vert, &
          id_vert, vp, h_p, second_order)

        do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
          idse = mesh%vert(id_vert)%sub_elem_neigh(j)
          ide = mesh%sub_elem(idse)%mesh_elem
          rhs(:, ide) = rhs(:, ide) - flux_sum_vert(:, j)
          sum_lambda(ide) = sum_lambda(ide) + sum_lambda_vert(j)
        end do
      end do
    else
      if (ener_dump_sides .and. ener_calls == 0) then
        block
          integer :: rk, ie
          character(len=64) :: fn
          call MPI_COMM_RANK(MPI_COMM_WORLD, rk, ie)
          write(fn, '(a,i0,a)') 'ener_sides_', rk, '.txt'
          open(newunit=ener_dump_unit, file=trim(fn), status='replace')
          ener_dumping = .true.
          write(ener_dump_unit, '(a)') '# vert j side lambda eps Delta u_n ut2 E_p E_s ratio is_ghost_vert'
        end block
      end if
      if (scheme_id == SCHEME_MULTI_POINT_EPS .or. scheme_id == SCHEME_MULTI_POINT_EPS_ISO .or. &
          scheme_id == SCHEME_MULTI_POINT_EPS2 .or. scheme_id == SCHEME_MULTI_POINT_EPS2M .or. &
          scheme_id == SCHEME_MULTI_POINT_PS .or. scheme_id == SCHEME_MULTI_POINT_PSW .or. &
          scheme_id == SCHEME_MULTI_POINT_PV) then
        if (.not. allocated(eps_cell_neg)) allocate(eps_cell_neg(mesh%n_elems), eps_cell_zf(mesh%n_elems))
        eps_cell_neg = 0
        eps_cell_zf = 0
      end if
      do id_vert = 1, mesh%n_vert
        nsfn = mesh%vert(id_vert)%n_sub_faces_neigh
        nsen = mesh%vert(id_vert)%n_sub_elems_neigh
        sum_lambda_vert(1:nsen) = 0.0_DOUBLE
        flux_sum_vert(:, 1:nsen) = 0.0_DOUBLE
        call compute_rhs_around_vert(mesh, sol, grad, &
          nsfn, nsen, v_vert, sum_lambda_vert, &
          flux_sum_vert, second_order, id_vert, vp)

        do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
          idse = mesh%vert(id_vert)%sub_elem_neigh(j)
          ide = mesh%sub_elem(idse)%mesh_elem
          rhs(:, ide) = rhs(:, ide) - flux_sum_vert(:, j)
          sum_lambda(ide) = sum_lambda(ide) + sum_lambda_vert(j)
        end do
      end do
      if (scheme_id == SCHEME_MULTI_POINT_ENER .or. scheme_id == SCHEME_MULTI_POINT_ENER_ISO .or. &
          scheme_id == SCHEME_MULTI_POINT_ENER_V .or. scheme_id == SCHEME_MULTI_POINT_ENER_V_ISO .or. &
          scheme_id == SCHEME_MULTI_POINT_EPS .or. scheme_id == SCHEME_MULTI_POINT_EPS_ISO .or. &
          scheme_id == SCHEME_MULTI_POINT_EPS2 .or. scheme_id == SCHEME_MULTI_POINT_EPS2M .or. &
          scheme_id == SCHEME_MULTI_POINT_PS .or. scheme_id == SCHEME_MULTI_POINT_PSW .or. &
          scheme_id == SCHEME_MULTI_POINT_PV) then
        ener_calls = ener_calls + 1
        if (ener_dumping) then
          close(ener_dump_unit)
          ener_dumping = .false.
        end if
      end if
    end if
    deallocate(sum_lambda_vert, flux_sum_vert)
  end subroutine compute_euler_flux

  subroutine compute_rhs_around_vert(mesh, sol, grad, &
      nsfn, nsen, v_vert, sum_lambda_vert, flux_sum_vert, &
      second_order, id_vert, vp)
    use ns_global_data_module, only: bc_style, scheme, exclude_bound_vert, &
      scheme_id, ener_lambda_pos, ener_internal, ener_pstar, &
      SCHEME_MULTI_POINT, SCHEME_MULTI_POINT_ISO, SCHEME_MULTI_POINT_PRESSURE, &
      SCHEME_MULTI_POINT_PRESSURE_PH, &
      SCHEME_THREE_WAVE, SCHEME_TWO_WAVE, SCHEME_MODIFIED_THREE_WAVE, &
      SCHEME_THREE_WAVE_ENTHALPY, SCHEME_MULTI_POINT_ENTHALPY, SCHEME_THREE_WAVE_ENTHALPY2, &
      SCHEME_MULTI_POINT_VILAR, &
      SCHEME_MULTI_POINT_ENER, SCHEME_MULTI_POINT_ENER_ISO, &
      SCHEME_MULTI_POINT_ENER_V, SCHEME_MULTI_POINT_ENER_V_ISO, &
      SCHEME_MULTI_POINT_EPS, SCHEME_MULTI_POINT_EPS_ISO, &
      SCHEME_MULTI_POINT_EPS2, SCHEME_MULTI_POINT_EPS2M, SCHEME_MULTI_POINT_PS, SCHEME_MULTI_POINT_PSW, &
      SCHEME_MULTI_POINT_PV, boundary_2d
    use linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    integer(kind=ENTIER), intent(in) :: nsfn, nsen
    real(kind=DOUBLE), dimension(3), intent(inout) :: v_vert
    real(kind=DOUBLE), dimension(nsen), intent(inout) :: &
      sum_lambda_vert
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: &
      flux_sum_vert
    logical, intent(in) :: second_order
    integer(kind=ENTIER), intent(in) :: id_vert

    integer(kind=ENTIER) :: j, id_sub_face, id_face
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE) :: p_bound, p_nodal, sl, sr
    real(kind=DOUBLE), dimension(5, 2, nsfn) :: sol_w_lr
    real(kind=DOUBLE), dimension(5, 2, nsfn) :: lr_flux
    real(kind=DOUBLE), dimension(2, nsfn) :: lambda
    real(kind=DOUBLE), dimension(nsfn) :: v_bars, p_bars, warea
    real(kind=DOUBLE), dimension(5, 2, nsfn) :: lr_flux_3w
    real(kind=DOUBLE), dimension(5) :: sol_p_vilar
    ! multi_point_ener* nodal passes
    real(kind=DOUBLE) :: p_ener, pp_den, e_ener, ener_num, ener_den, invlam
    real(kind=DOUBLE) :: rho_l, rho_r, vn_l, vn_r, p_l, p_r, e_l, e_r, mu_l, mu_r
    real(kind=DOUBLE), dimension(nsfn) :: pv_l_ener, pv_r_ener, vs_ener
    logical :: ener_active, ener_v, ener_iso, is_zf, ener_raised, ener_eps
    real(kind=DOUBLE) :: eps_l, eps_r, dsl_l, dsl_r
    integer(kind=ENTIER) :: eps_it, eps_nr, eps_ns, eps_kind
    real(kind=DOUBLE) :: eps_min, eps_tgt
    logical, dimension(2, nsfn) :: ener_raised_side
    integer(kind=ENTIER) :: ener_pass, k_side
    real(kind=DOUBLE) :: eps_s, d_s, b_s, ratio, lam_req, lam_ratio
    real(kind=DOUBLE), dimension(2, nsfn) :: lam_cl, lam_wave
    real(kind=DOUBLE) :: w_ener, ustar_l, ustar_r
    real(kind=DOUBLE), dimension(5) :: u_tmp
    rse_loc = 0
    ener_eps = .false.
    ener_iso = .false.
    eps_kind = 0
    flux_sum_vert = 0.0_DOUBLE
    v_vert = 0.0_DOUBLE
    sum_lambda_vert = 0.0_DOUBLE
    sol_w_lr = 0.0_DOUBLE
    lr_flux = 0.0_DOUBLE
    lr_flux_3w= 0.0_DOUBLE
    v_bars = 0.0_DOUBLE
    warea = 0.0_DOUBLE

    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face

      le = mesh%face(id_face)%left_neigh
      re = mesh%face(id_face)%right_neigh

      warea(j) = mesh%sub_face(id_sub_face)%area

      call reconstruct_lr_w(mesh, sol, grad, id_vert, id_sub_face, le, re, &
        second_order, sol_w_lr(:, 1, j), sol_w_lr(:, 2, j))

    end do

    warea = minval(warea)/warea

    select case (scheme_id)
    case (SCHEME_MULTI_POINT, SCHEME_MULTI_POINT_ENTHALPY, &
        SCHEME_MULTI_POINT_ENER, SCHEME_MULTI_POINT_ENER_V, SCHEME_MULTI_POINT_EPS, &
        SCHEME_MULTI_POINT_EPS2, SCHEME_MULTI_POINT_EPS2M, SCHEME_MULTI_POINT_PS, &
        SCHEME_MULTI_POINT_PSW, SCHEME_MULTI_POINT_PV)
      ! The nodal velocity comes from the paper's linear system (16), which involves only the
      ! slopes and the face velocities -- no energy -- so the enthalpy-preserving variant
      ! shares this step with multi_point unchanged.
      lambda(:, :) = 0.0_DOUBLE
      if( mesh%vert(id_vert)%is_bound .and. exclude_bound_vert) then
        v_vert = 0.0_DOUBLE
        vp(:, id_vert) = v_vert
      else
        call compute_lambdas_and_solve_nodal_velocity(mesh, id_vert, &
          sol_w_lr, lambda, v_bars, v_vert, p_bound)
        vp(:, id_vert) = v_vert
      end if
    case (SCHEME_MULTI_POINT_ISO, SCHEME_MULTI_POINT_ENER_ISO, SCHEME_MULTI_POINT_ENER_V_ISO, &
        SCHEME_MULTI_POINT_EPS_ISO)
      lambda(:, :) = 0.0_DOUBLE
      if( mesh%vert(id_vert)%is_bound .and. exclude_bound_vert) then
        v_vert = 0.0_DOUBLE
        vp(:, id_vert) = v_vert
      else
        call compute_lambdas_and_solve_nodal_velocity_iso(mesh, id_vert, &
          sol_w_lr, lambda, v_bars, v_vert, p_bound)
        vp(:, id_vert) = v_vert
      end if
    case (SCHEME_MULTI_POINT_VILAR)
      sol_p_vilar = 0.0_DOUBLE
      call compute_nodal_state_VILAR(mesh, id_vert, sol_w_lr, sol_p_vilar)
    case (SCHEME_MULTI_POINT_PRESSURE)
      lambda(:, :) = 0.0_DOUBLE
      p_nodal = 0.0_DOUBLE
      call compute_lambdas_and_solve_nodal_pressure(mesh, id_vert, &
        sol_w_lr, lambda, p_bars, p_nodal)
    case (SCHEME_MULTI_POINT_PRESSURE_PH)
      lambda(:, :) = 0.0_DOUBLE
      p_nodal = 0.0_DOUBLE
      call compute_lambdas_and_solve_nodal_pressure(mesh, id_vert, &
        sol_w_lr, lambda, p_bars, p_nodal)
    end select

    ! --- multi_point_ener*: nodal pressure p_p and nodal energy e_p ---
    ! Same lambdas as the star states (the iterated ones from the nodal velocity solve).
    ! p_p : 1/lambda-weighted mean of the acoustic interface pressures, all sub-faces but
    !       the boundary_2d z-faces (the boundary ones through their ghost state), exactly
    !       as the Lagrangian vp_pp pass.
    ! e_p : telescoping definition, INTERIOR sub-faces only. A boundary sub-face has no
    !       partner to telescope with and lambda (e_p - e) carries no velocity factor: kept
    !       in, it leaks energy through the boundary and is a fake inviscid wall heat flux.
    !       Boundary sub-faces keep the energy flux p_p v* (zero at a wall, v* = 0).
    ! mu  : lambda (multi_point_ener*) or rho |v_s - v*| (multi_point_ener_v*), v* = v_p.n.
    !       A degenerate node (sum A (mu_l + mu_r) = 0) uses p_p v* on every sub-face, which
    !       is conservative since both sides then carry the same energy flux.
    ener_active = .false.
    if (scheme_id == SCHEME_MULTI_POINT_ENER .or. scheme_id == SCHEME_MULTI_POINT_ENER_ISO .or. &
        scheme_id == SCHEME_MULTI_POINT_ENER_V .or. scheme_id == SCHEME_MULTI_POINT_ENER_V_ISO .or. &
          scheme_id == SCHEME_MULTI_POINT_EPS .or. scheme_id == SCHEME_MULTI_POINT_EPS_ISO .or. &
          scheme_id == SCHEME_MULTI_POINT_EPS2 .or. scheme_id == SCHEME_MULTI_POINT_EPS2M .or. &
          scheme_id == SCHEME_MULTI_POINT_PS .or. scheme_id == SCHEME_MULTI_POINT_PSW .or. &
          scheme_id == SCHEME_MULTI_POINT_PV) then
      ener_v = (scheme_id == SCHEME_MULTI_POINT_ENER_V .or. scheme_id == SCHEME_MULTI_POINT_ENER_V_ISO)
      ener_iso = (scheme_id == SCHEME_MULTI_POINT_ENER_ISO .or. scheme_id == SCHEME_MULTI_POINT_ENER_V_ISO .or. &
        scheme_id == SCHEME_MULTI_POINT_EPS_ISO)
      ! eps_kind: 1 eps (eps*_l = eps*_r = eps_p), 2 eps2 (eps*_s = eps_s + kappa),
      ! 3 eps2m (eps*_s = eps_s (1 + kappa)); e_ener holds eps_p (1) or kappa (2, 3)
      eps_kind = 0
      if (scheme_id == SCHEME_MULTI_POINT_EPS .or. scheme_id == SCHEME_MULTI_POINT_EPS_ISO) eps_kind = 1
      if (scheme_id == SCHEME_MULTI_POINT_EPS2) eps_kind = 2
      if (scheme_id == SCHEME_MULTI_POINT_EPS2M) eps_kind = 3
      ! 4: multi_point_ps, (pv)*_s = p_p u*_s: eps*_s = eps_s + (p_p^2 - p_s^2)/(2 lam_s^2) >= 0 with the
      ! classic lambda (lam >= rho a); no bound, no nodal energy (e_ener unused, kept 0)
      if (scheme_id == SCHEME_MULTI_POINT_PS) eps_kind = 4
      ! 5: multi_point_psw, ps with p_p over interior sub-faces and classic multi_point on boundary/z
      if (scheme_id == SCHEME_MULTI_POINT_PSW) eps_kind = 5
      ! 6: multi_point_pv, psw with the classic energy flux Pi_s v* on interior sub-faces
      if (scheme_id == SCHEME_MULTI_POINT_PV) eps_kind = 6
      ener_eps = (eps_kind > 0)
      ener_raised_side = .false.
      ! lambda positivity bound (ener_lambda_pos). With option A the star internal energy is
      !   eps*_l = eps_l + B_l/lambda_l - D_l^2/(2 lambda_l^2),  B_l = p_p u_l - pv_l,  D_l = p_l - p_p
      !   eps*_r = eps_r + B_r/lambda_r - D_r^2/(2 lambda_r^2),  B_r = pv_r - p_p u_r,  D_r = p_r - p_p
      ! (u_s the normal velocity of the cell state), checked against the star states the solver
      ! builds. eps*_s >= 0 is a convex parabola in lambda with one negative and one positive
      ! root, so lambda_req,s = (-B_s + sqrt(B_s^2 + 2 eps_s D_s^2))/(2 eps_s) always exists.
      ! p_p, v* and e_p depend on every lambda of the node, hence a fixed point: at most 5
      ! passes, lambdas only increase; the last pass's violators are counted, never clipped.
      ener_pass = 0
      do while (.true.)
        ener_pass = ener_pass + 1
        if (ener_pass > 1) then
          if (ener_iso) then
            call compute_lambdas_and_solve_nodal_velocity_iso(mesh, id_vert, &
              sol_w_lr, lambda, v_bars, v_vert, p_bound)
          else
            call compute_lambdas_and_solve_nodal_velocity(mesh, id_vert, &
              sol_w_lr, lambda, v_bars, v_vert, p_bound)
          end if
          vp(:, id_vert) = v_vert
        end if
        p_ener = 0.0_DOUBLE
        pp_den = 0.0_DOUBLE
        ener_num = 0.0_DOUBLE
        ener_den = 0.0_DOUBLE
        do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
          id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
          id_face = mesh%sub_face(id_sub_face)%mesh_face
          re = mesh%face(id_face)%right_neigh
          if (is_wall(re)) then
            vs_ener(j) = 0.0_DOUBLE
          else
            vs_ener(j) = dot_product(v_vert, mesh%sub_face(id_sub_face)%norm)
          end if
          is_zf = boundary_2d .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE
          if (is_zf) cycle
          rho_l = sol_w_lr(1, 1, j)
          vn_l = dot_product(sol_w_lr(2:4, 1, j), mesh%sub_face(id_sub_face)%norm)
          p_l = sol_w_lr(5, 1, j)
          rho_r = sol_w_lr(1, 2, j)
          vn_r = dot_product(sol_w_lr(2:4, 2, j), mesh%sub_face(id_sub_face)%norm)
          p_r = sol_w_lr(5, 2, j)
          invlam = 1.0_DOUBLE/lambda(1, j) + 1.0_DOUBLE/lambda(2, j)
          if (eps_kind >= 5 .and. re <= 0) cycle
          p_ener = p_ener + mesh%sub_face(id_sub_face)%area &
            *(p_l/lambda(1, j) + p_r/lambda(2, j) - (vn_r - vn_l))
          pp_den = pp_den + mesh%sub_face(id_sub_face)%area*invlam
          if (re <= 0) cycle
          u_tmp = primit_to_conserv(sol_w_lr(:, 1, j))
          e_l = u_tmp(5)/rho_l
          u_tmp = primit_to_conserv(sol_w_lr(:, 2, j))
          e_r = u_tmp(5)/rho_r
          if (ener_internal .and. .not. ener_eps) then
            e_l = e_l - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 1, j), sol_w_lr(2:4, 1, j))
            e_r = e_r - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 2, j), sol_w_lr(2:4, 2, j))
          end if
          if (ener_v) then
            mu_l = rho_l*abs(vn_l - vs_ener(j))
            mu_r = rho_r*abs(vn_r - vs_ener(j))
          else
            mu_l = lambda(1, j)
            mu_r = lambda(2, j)
          end if
          ! iso energy solver (the _iso variants): no areas, the counterpart of the iso
          ! velocity solver. Their ener flux enters weighted by A_j warea_j = min_j A_j, a
          ! constant over the node, so it telescopes exactly when e_p carries no areas;
          ! with the areas it leaked energy (6e-10 on Sedov).
          if (ener_iso) then
            w_ener = 1.0_DOUBLE
          else
            w_ener = mesh%sub_face(id_sub_face)%area
          end if
          ener_num = ener_num + w_ener*(mu_l*e_l + mu_r*e_r + (p_l*vn_l - p_r*vn_r))
          ener_den = ener_den + w_ener*(mu_l + mu_r)
        end do
        if (pp_den > 0.0_DOUBLE) p_ener = p_ener/pp_den
        ! multi_point_eps: lambda lower bound per SIDE (Vincent, 2026-10-29). With
        ! D_l = p_l - p_p, D_r = p_p - p_r, the bracket of each sub-face in the eps_p numerator splits as
        !   [lam_l eps_l + p_p u_l - D_l^2/(2 lam_l)] + [lam_r eps_r - p_p u_r - D_r^2/(2 lam_r)],
        ! so requiring each half >= 0, i.e. the convex parabolas
        !   eps_l lam_l^2 + p_p u_l lam_l - D_l^2/2 >= 0,   eps_r lam_r^2 - p_p u_r lam_r - D_r^2/2 >= 0,
        ! gives eps_p >= 0, hence eps*_l = eps*_r = eps_p >= 0 on every interior sub-face. A boundary or
        ! z-face sub-face (flux L_l, pv_l = p_p v*) has eps*_l lam^2 = eps_l lam^2 + B lam - D^2/2 with
        ! B = p_p (u_l - v*): the same parabola (B = p_p u_l at a wall or z-face, v* = 0).
        ! lam <- max(lam, root): lambdas only increase, each side keeps its own. p_p depends on the
        ! lambdas, so it is recomputed twice; the third raise is made against the FROZEN p_p, which then
        ! guarantees the condition exactly (p_p need not be the 1/lambda mean: p_p n and the eps_p
        ! telescoping are conservative for any p_p). v* is the classic one (eps*_s = eps_p does not
        ! depend on v*).
        if (ener_eps) then
          lam_cl = lambda
          ! eps2: the per-side parabolas use the node's smallest interior internal energy eps_min
          ! (kappa >= -eps_min <=> every eps*_s = eps_s + kappa >= 0); eps / eps2m use eps_s
          eps_min = huge(1.0_DOUBLE)
          do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
            id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
            id_face = mesh%sub_face(id_sub_face)%mesh_face
            re = mesh%face(id_face)%right_neigh
            is_zf = boundary_2d .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE
            if (re <= 0 .or. is_zf) cycle
            do k_side = 1, 2
              u_tmp = primit_to_conserv(sol_w_lr(:, k_side, j))
              eps_min = min(eps_min, u_tmp(5)/sol_w_lr(1, k_side, j) &
                - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, k_side, j), sol_w_lr(2:4, k_side, j)))
            end do
          end do
          do eps_it = 1, 3
            if (eps_kind >= 4) exit
            eps_nr = 0
            do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
              id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
              id_face = mesh%sub_face(id_sub_face)%mesh_face
              re = mesh%face(id_face)%right_neigh
              is_zf = boundary_2d .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE
              do k_side = 1, 2
                if (k_side == 2 .and. (re <= 0 .or. is_zf)) cycle
                rho_l = sol_w_lr(1, k_side, j)
                vn_l = dot_product(sol_w_lr(2:4, k_side, j), mesh%sub_face(id_sub_face)%norm)
                u_tmp = primit_to_conserv(sol_w_lr(:, k_side, j))
                eps_s = u_tmp(5)/rho_l - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, k_side, j), sol_w_lr(2:4, k_side, j))
                if (eps_kind == 2 .and. re > 0 .and. .not. is_zf) eps_s = eps_min
                d_s = sol_w_lr(5, k_side, j) - p_ener
                if (k_side == 1) then
                  if (re <= 0 .or. is_zf) then
                    b_s = p_ener*(vn_l - vs_ener(j))
                  else
                    b_s = p_ener*vn_l
                  end if
                else
                  b_s = -p_ener*vn_l
                end if
                ! positive root, in the form without cancellation (b >= 0: D^2/(b + sq))
                lam_req = sqrt(b_s**2 + 2.0_DOUBLE*eps_s*d_s**2)
                if (b_s >= 0.0_DOUBLE) then
                  if (lam_req + b_s > 0.0_DOUBLE) then
                    lam_req = d_s**2/(b_s + lam_req)
                  else
                    lam_req = 0.0_DOUBLE
                  end if
                else
                  lam_req = (lam_req - b_s)/(2.0_DOUBLE*eps_s)
                end if
                if (lam_req > lambda(k_side, j)) then
                  lambda(k_side, j) = lam_req
                  eps_nr = eps_nr + 1
                end if
              end do
            end do
            if (eps_nr == 0 .or. eps_it == 3) exit
            ! p_p with the raised lambdas (same sub-faces as above: all but the z-faces)
            p_ener = 0.0_DOUBLE
            pp_den = 0.0_DOUBLE
            do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
              id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
              is_zf = boundary_2d .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE
              if (is_zf) cycle
              vn_l = dot_product(sol_w_lr(2:4, 1, j), mesh%sub_face(id_sub_face)%norm)
              vn_r = dot_product(sol_w_lr(2:4, 2, j), mesh%sub_face(id_sub_face)%norm)
              p_ener = p_ener + mesh%sub_face(id_sub_face)%area &
                *(sol_w_lr(5, 1, j)/lambda(1, j) + sol_w_lr(5, 2, j)/lambda(2, j) - (vn_r - vn_l))
              pp_den = pp_den + mesh%sub_face(id_sub_face)%area*(1.0_DOUBLE/lambda(1, j) + 1.0_DOUBLE/lambda(2, j))
            end do
            if (pp_den > 0.0_DOUBLE) p_ener = p_ener/pp_den
          end do
          ! nodes whose EVERY eligible side ended up raised: only there is eps_p driven to 0
          if (.not. mesh%vert(id_vert)%is_ghost) then
            eps_nr = 0
            eps_ns = 0
            do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
              id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
              id_face = mesh%sub_face(id_sub_face)%mesh_face
              re = mesh%face(id_face)%right_neigh
              is_zf = boundary_2d .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE
              do k_side = 1, 2
                if (k_side == 2 .and. (re <= 0 .or. is_zf)) cycle
                eps_ns = eps_ns + 1
                if (lambda(k_side, j) > lam_cl(k_side, j)) eps_nr = eps_nr + 1
              end do
            end do
          end if
        end if
        ! ener_pstar: X_p from sum_f w [(pv)*_r - (pv)*_l] = 0 with the star-velocity pressure work,
        ! X_p sum w (mu_l + mu_r) = sum w [mu_l X_l + mu_r X_r - p_p (u*_r - u*_l)]; it needs p_p,
        ! hence this second pass over the sub-faces.
        if (ener_pstar .or. ener_eps) then
          ener_num = 0.0_DOUBLE
          ener_den = 0.0_DOUBLE
          do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
            id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
            id_face = mesh%sub_face(id_sub_face)%mesh_face
            re = mesh%face(id_face)%right_neigh
            is_zf = boundary_2d .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE
            if (is_zf .or. re <= 0) cycle
            rho_l = sol_w_lr(1, 1, j)
            vn_l = dot_product(sol_w_lr(2:4, 1, j), mesh%sub_face(id_sub_face)%norm)
            p_l = sol_w_lr(5, 1, j)
            rho_r = sol_w_lr(1, 2, j)
            vn_r = dot_product(sol_w_lr(2:4, 2, j), mesh%sub_face(id_sub_face)%norm)
            p_r = sol_w_lr(5, 2, j)
            u_tmp = primit_to_conserv(sol_w_lr(:, 1, j))
            e_l = u_tmp(5)/rho_l
            u_tmp = primit_to_conserv(sol_w_lr(:, 2, j))
            e_r = u_tmp(5)/rho_r
            if (ener_internal .and. .not. ener_eps) then
              e_l = e_l - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 1, j), sol_w_lr(2:4, 1, j))
              e_r = e_r - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 2, j), sol_w_lr(2:4, 2, j))
            end if
            if (ener_v) then
              mu_l = rho_l*abs(vn_l - vs_ener(j))
              mu_r = rho_r*abs(vn_r - vs_ener(j))
            else
              mu_l = lambda(1, j)
              mu_r = lambda(2, j)
            end if
            if (ener_iso) then
              w_ener = 1.0_DOUBLE
            else
              w_ener = mesh%sub_face(id_sub_face)%area
            end if
            if (eps_kind >= 4) then
              ener_den = ener_den + w_ener*(lambda(1, j) + lambda(2, j))
            else if (eps_kind >= 2) then
              ! eps2 / eps2m: with eps*_s imposed, the jump relation gives
              !   (pv)*_l = p_p u_l - D_l^2/(2 lam_l) - lam_l c_l kappa
              !   (pv)*_r = p_p u_r + D_r^2/(2 lam_r) + lam_r c_r kappa
              ! (c_s = 1 for eps2, eps_s for eps2m), D_l = p_l - p_p, D_r = p_p - p_r: no relaxation term
              ! left, hence no conduction. sum_f w [(pv)*_r - (pv)*_l] = 0 fixes
              ! kappa = -N / sum w (lam_l c_l + lam_r c_r), N = sum w [p_p (u_r - u_l) + D_r^2/(2 lam_r)
              ! + D_l^2/(2 lam_l)], accumulated here as num = -N.
              eps_l = e_l - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 1, j), sol_w_lr(2:4, 1, j))
              eps_r = e_r - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 2, j), sol_w_lr(2:4, 2, j))
              ener_num = ener_num - w_ener*(p_ener*(vn_r - vn_l) + 0.5_DOUBLE*(p_ener - p_r)**2/lambda(2, j) &
                + 0.5_DOUBLE*(p_l - p_ener)**2/lambda(1, j))
              if (eps_kind == 2) then
                ener_den = ener_den + w_ener*(lambda(1, j) + lambda(2, j))
              else
                ener_den = ener_den + w_ener*(lambda(1, j)*eps_l + lambda(2, j)*eps_r)
              end if
            else if (ener_eps) then
              ! multi_point_eps: e*_s = eps_p + |u*_s|^2/2 in the Lagrangian energy line, with the
              ! acoustic star velocity u*_s = u_s + d_s n (tangential part unchanged),
              ! d_l = (p_l - p_p)/lambda_l, d_r = (p_p - p_r)/lambda_r, so that
              ! e_s - |u*_s|^2/2 = eps_s - u_s d_s - d_s^2/2 (same sign on both sides), and
              ! eps_p sum w (lam_l + lam_r) = sum w [lam_l (eps_l - u_l d_l - d_l^2/2)
              !                              + lam_r (eps_r - u_r d_r - d_r^2/2) + p_l u_l - p_r u_r]
              ! from sum_f w [(pv)*_r - (pv)*_l] = 0. Then eps*_l = eps*_r = eps_p identically
              ! (option A star states): positivity is the single condition eps_p > 0.
              eps_l = e_l - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 1, j), sol_w_lr(2:4, 1, j))
              eps_r = e_r - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 2, j), sol_w_lr(2:4, 2, j))
              dsl_l = (p_l - p_ener)/lambda(1, j)
              dsl_r = (p_ener - p_r)/lambda(2, j)
              ener_num = ener_num + w_ener*(lambda(1, j)*(eps_l - vn_l*dsl_l - 0.5_DOUBLE*dsl_l**2) &
                + lambda(2, j)*(eps_r - vn_r*dsl_r - 0.5_DOUBLE*dsl_r**2) + (p_l*vn_l - p_r*vn_r))
              ener_den = ener_den + w_ener*(lambda(1, j) + lambda(2, j))
            else
              ustar_l = vn_l + (p_l - p_ener)/lambda(1, j)
              ustar_r = vn_r + (p_ener - p_r)/lambda(2, j)
              ener_num = ener_num + w_ener*(mu_l*e_l + mu_r*e_r - p_ener*(ustar_r - ustar_l))
              ener_den = ener_den + w_ener*(mu_l + mu_r)
            end if
          end do
        end if
        e_ener = 0.0_DOUBLE
        ener_active = (ener_den > 0.0_DOUBLE)
        if (ener_active .and. eps_kind < 4) e_ener = ener_num/ener_den
        ! multi_point_eps: e_ener is eps_p, the only positivity condition of the node; eps2: some
        ! eps*_s = eps_s + kappa <= 0 iff kappa <= -eps_min; eps2m: iff kappa <= -1
        eps_tgt = 0.0_DOUBLE
        if (eps_kind == 2) eps_tgt = -eps_min
        if (eps_kind == 3) eps_tgt = -1.0_DOUBLE
        if (eps_kind >= 4) eps_tgt = -huge(1.0_DOUBLE)
        if (ener_eps .and. ener_pass == 1 .and. ener_active .and. .not. mesh%vert(id_vert)%is_ghost) then
          if (e_ener <= eps_tgt) then
            do j = 1, mesh%vert(id_vert)%n_sub_elems_neigh
              le = mesh%sub_elem(mesh%vert(id_vert)%sub_elem_neigh(j))%mesh_elem
              eps_cell_neg(le) = eps_cell_neg(le) + 1
            end do
          end if
        end if

        ! one-sided energy fluxes per sub-face
        do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
          id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
          id_face = mesh%sub_face(id_sub_face)%mesh_face
          re = mesh%face(id_face)%right_neigh
          is_zf = boundary_2d .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE
          pv_l_ener(j) = p_ener*vs_ener(j)
          pv_r_ener(j) = p_ener*vs_ener(j)
          if (re <= 0 .or. is_zf .or. .not. ener_active) cycle
          rho_l = sol_w_lr(1, 1, j)
          vn_l = dot_product(sol_w_lr(2:4, 1, j), mesh%sub_face(id_sub_face)%norm)
          p_l = sol_w_lr(5, 1, j)
          rho_r = sol_w_lr(1, 2, j)
          vn_r = dot_product(sol_w_lr(2:4, 2, j), mesh%sub_face(id_sub_face)%norm)
          p_r = sol_w_lr(5, 2, j)
          u_tmp = primit_to_conserv(sol_w_lr(:, 1, j))
          e_l = u_tmp(5)/rho_l
          u_tmp = primit_to_conserv(sol_w_lr(:, 2, j))
          e_r = u_tmp(5)/rho_r
          if (ener_internal .and. .not. ener_eps) then
            e_l = e_l - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 1, j), sol_w_lr(2:4, 1, j))
            e_r = e_r - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 2, j), sol_w_lr(2:4, 2, j))
          end if
          if (ener_v) then
            mu_l = rho_l*abs(vn_l - vs_ener(j))
            mu_r = rho_r*abs(vn_r - vs_ener(j))
          else
            mu_l = lambda(1, j)
            mu_r = lambda(2, j)
          end if
          if (eps_kind == 6) then
            ! classic energy flux: Pi_l = p_l + lam_l (u_l - v*), Pi_r = p_r - lam_r (u_r - v*), times v*
            pv_l_ener(j) = (p_l + lambda(1, j)*(vn_l - vs_ener(j)))*vs_ener(j)
            pv_r_ener(j) = (p_r - lambda(2, j)*(vn_r - vs_ener(j)))*vs_ener(j)
          else if (eps_kind >= 4) then
            pv_l_ener(j) = p_ener*(vn_l + (p_l - p_ener)/lambda(1, j))
            pv_r_ener(j) = p_ener*(vn_r + (p_ener - p_r)/lambda(2, j))
          else if (eps_kind == 2) then
            pv_l_ener(j) = p_ener*vn_l - 0.5_DOUBLE*(p_l - p_ener)**2/lambda(1, j) - lambda(1, j)*e_ener
            pv_r_ener(j) = p_ener*vn_r + 0.5_DOUBLE*(p_ener - p_r)**2/lambda(2, j) + lambda(2, j)*e_ener
          else if (eps_kind == 3) then
            eps_l = e_l - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 1, j), sol_w_lr(2:4, 1, j))
            eps_r = e_r - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 2, j), sol_w_lr(2:4, 2, j))
            pv_l_ener(j) = p_ener*vn_l - 0.5_DOUBLE*(p_l - p_ener)**2/lambda(1, j) - lambda(1, j)*eps_l*e_ener
            pv_r_ener(j) = p_ener*vn_r + 0.5_DOUBLE*(p_ener - p_r)**2/lambda(2, j) + lambda(2, j)*eps_r*e_ener
          else if (ener_eps) then
            ! (pv)*_l = p_l u_l - lam_l (eps_p + |u*_l|^2/2 - e_l),  (pv)*_r = p_r u_r + lam_r (...)
            ! (e_l, e_r are the TOTAL energies here: ener_internal plays no role)
            eps_l = e_l - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 1, j), sol_w_lr(2:4, 1, j))
            eps_r = e_r - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 2, j), sol_w_lr(2:4, 2, j))
            dsl_l = (p_l - p_ener)/lambda(1, j)
            dsl_r = (p_ener - p_r)/lambda(2, j)
            pv_l_ener(j) = p_l*vn_l - lambda(1, j)*(e_ener - eps_l + vn_l*dsl_l + 0.5_DOUBLE*dsl_l**2)
            pv_r_ener(j) = p_r*vn_r + lambda(2, j)*(e_ener - eps_r + vn_r*dsl_r + 0.5_DOUBLE*dsl_r**2)
          else if (ener_pstar) then
            pv_l_ener(j) = p_ener*(vn_l + (p_l - p_ener)/lambda(1, j)) - mu_l*(e_ener - e_l)
            pv_r_ener(j) = p_ener*(vn_r + (p_ener - p_r)/lambda(2, j)) + mu_r*(e_ener - e_r)
          else
            pv_l_ener(j) = p_l*vn_l - mu_l*(e_ener - e_l)
            pv_r_ener(j) = p_r*vn_r + mu_r*(e_ener - e_r)
          end if
        end do

        ! nodal energy balance: energy created at the node = sum over interior sub-faces of
        ! A ((pv)*_r - (pv)*_l) (a boundary sub-face's pv is a genuine outflow). Relative to
        ! sum A p_p (lam_l/rho_l + lam_r/rho_r)/2, a p*a scale that does not vanish at rest (sum |pv|
        ! does: 0/0 in a fluid at rest); interior and boundary nodes kept apart. Not for _iso (A warea).
        if (ener_pass == 1 .and. .not. ener_iso .and. .not. mesh%vert(id_vert)%is_ghost) then
          ener_num = 0.0_DOUBLE
          ener_den = 0.0_DOUBLE
          do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
            id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
            id_face = mesh%sub_face(id_sub_face)%mesh_face
            re = mesh%face(id_face)%right_neigh
            is_zf = boundary_2d .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE
            if (re <= 0 .or. is_zf) cycle
            ener_num = ener_num + mesh%sub_face(id_sub_face)%area*(pv_r_ener(j) - pv_l_ener(j))
            p_l = sol_w_lr(5, 1, j)
            p_r = sol_w_lr(5, 2, j)
            ratio = 0.5_DOUBLE*((1.0_DOUBLE/lambda(1, j) + 1.0_DOUBLE/lambda(2, j))*p_ener**2 &
              - p_l**2/lambda(1, j) - p_r**2/lambda(2, j)) &
              /(max(p_l, p_r)*0.5_DOUBLE*(lambda(1, j)/sol_w_lr(1, 1, j) + lambda(2, j)/sol_w_lr(1, 2, j)))
            ener_den = ener_den + mesh%sub_face(id_sub_face)%area*abs(p_ener) &
              *0.5_DOUBLE*(lambda(1, j)/sol_w_lr(1, 1, j) + lambda(2, j)/sol_w_lr(1, 2, j))
          end do
          if (ener_den > 0.0_DOUBLE) then
            k_side = 1
            if (mesh%vert(id_vert)%is_bound) k_side = 2
          end if
        end if

        ! star internal energies with the current lambdas, and the bound
        ener_raised = .false.
        lam_wave = lambda
        do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
          id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
          id_face = mesh%sub_face(id_sub_face)%mesh_face
          re = mesh%face(id_face)%right_neigh
          is_zf = boundary_2d .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE
          if (re <= 0 .or. is_zf) then
            ! not bounded: boundary and boundary_2d z-face sub-faces. Their LEFT star state
            ! still enters the cell's convex combination (s_l < 0), so count it.
            ! (psw, pv: these sub-faces carry the classic multi_point flux, not p_p: not counted)
            if (ener_pass == 1 .and. .not. mesh%vert(id_vert)%is_ghost .and. eps_kind < 5) then
              rho_l = sol_w_lr(1, 1, j)
              vn_l = dot_product(sol_w_lr(2:4, 1, j), mesh%sub_face(id_sub_face)%norm)
              u_tmp = primit_to_conserv(sol_w_lr(:, 1, j))
              eps_s = u_tmp(5)/rho_l - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 1, j), sol_w_lr(2:4, 1, j))
              d_s = sol_w_lr(5, 1, j) - p_ener
              b_s = p_ener*vn_l - pv_l_ener(j)
              ratio = (eps_s + b_s/lambda(1, j) - 0.5_DOUBLE*d_s**2/lambda(1, j)**2)/eps_s
              if (ratio < -1e-12_DOUBLE) then
                if (is_zf) then
                  if (ener_eps) then
                    le = mesh%face(id_face)%left_neigh
                    eps_cell_zf(le) = eps_cell_zf(le) + 1
                  end if
                else
                end if
              end if
              ! bound the left wave of these sub-faces too. With v* = 0 (wall, z-face) the
              ! flux is L_l whatever the wave lambda, so this only lowers the cell's local
              ! step through s_l (sum_lambda): exactly what the convex-combination argument
              ! needs, an admissible U*_l with a compatible step. Same closed root, B frozen.
              if (ener_lambda_pos .and. ratio < 0.0_DOUBLE) then
                lam_req = (-b_s + sqrt(b_s**2 + 2.0_DOUBLE*eps_s*d_s**2))/(2.0_DOUBLE*eps_s)
                if (lam_req > lam_wave(1, j)) then
                  lam_wave(1, j) = lam_req
                end if
              end if
            end if
            cycle
          end if
          do k_side = 1, 2
            rho_l = sol_w_lr(1, k_side, j)
            vn_l = dot_product(sol_w_lr(2:4, k_side, j), mesh%sub_face(id_sub_face)%norm)
            u_tmp = primit_to_conserv(sol_w_lr(:, k_side, j))
            eps_s = u_tmp(5)/rho_l - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, k_side, j), sol_w_lr(2:4, k_side, j))
            d_s = sol_w_lr(5, k_side, j) - p_ener
            if (k_side == 1) then
              b_s = p_ener*vn_l - pv_l_ener(j)
            else
              b_s = pv_r_ener(j) - p_ener*vn_l
            end if
            ratio = (eps_s + b_s/lambda(k_side, j) - 0.5_DOUBLE*d_s**2/lambda(k_side, j)**2)/eps_s
            if (ener_pass == 1 .and. ener_dumping) then
              write(ener_dump_unit, '(i8,i4,i2,8es16.8,l2)') id_vert, j, k_side, lambda(k_side, j), eps_s, &
                d_s, vn_l, dot_product(sol_w_lr(2:4, k_side, j), sol_w_lr(2:4, k_side, j)) - vn_l**2, &
                e_ener, u_tmp(5)/rho_l, ratio, mesh%vert(id_vert)%is_ghost
            end if
            ! ghost vertices: their fluxes only feed ghost cells, refreshed by the MPI exchange;
            ! neither bounded nor counted (they inflated the first counts ~5x)
            if (mesh%vert(id_vert)%is_ghost) cycle
            ! multi_point_eps: the closed-form eps*_s must equal eps_p on every interior side
            if (ener_eps .and. ener_pass == 1 .and. ener_active) then
              eps_tgt = e_ener
              if (eps_kind == 2) eps_tgt = eps_s + e_ener
              if (eps_kind == 3) eps_tgt = eps_s*(1.0_DOUBLE + e_ener)
              if (eps_kind == 6) eps_tgt = eps_s*ratio
              if (eps_kind == 4 .or. eps_kind == 5) &
                eps_tgt = eps_s + 0.5_DOUBLE*(p_ener**2 - (d_s + p_ener)**2)/lambda(k_side, j)**2
            end if
            if (ener_pass == 1) then
              if (.not. ener_eps) lam_cl(k_side, j) = lambda(k_side, j)
            end if
            ! DECOUPLED bound (Vincent): the nodal stage (p_p, v*, e_p, mu) keeps the classic
            ! lambda; only the WAVE lambda of the star states (s_l, s_r) is raised. With the nodal
            ! stage frozen, pv_s is a constant, so B_s and D_s do not depend on the wave lambda:
            ! eps*(lam) = eps + B/lam - D^2/(2 lam^2) is convex in 1/lam for ALL four variants
            ! and its positive root always exists (eps > 0). One pass, no fixed point. The jump
            ! F^+ - F^- = L_r - L_l depends on the nodal stage only, so conservation and the
            ! telescoping are untouched.
            if (ener_lambda_pos .and. ratio < 0.0_DOUBLE) then
              lam_req = (-b_s + sqrt(b_s**2 + 2.0_DOUBLE*eps_s*d_s**2))/(2.0_DOUBLE*eps_s)
              if (lam_req > lam_wave(k_side, j)) then
                lam_wave(k_side, j) = lam_req
              end if
            end if
          end do
        end do
        if (.not. ener_lambda_pos .or. .not. ener_raised .or. ener_pass == 5) exit
      end do
      if (.not. mesh%vert(id_vert)%is_ghost) then
      end if
      ! NB: when the 5th pass still raised a lambda, p_p/e_p/pv are those of the lambdas BEFORE
      ! that raise: recompute them once more so that fluxes and lambdas are consistent.
      if (ener_lambda_pos .and. ener_raised) then
        if (ener_iso) then
          call compute_lambdas_and_solve_nodal_velocity_iso(mesh, id_vert, &
            sol_w_lr, lambda, v_bars, v_vert, p_bound)
        else
          call compute_lambdas_and_solve_nodal_velocity(mesh, id_vert, &
            sol_w_lr, lambda, v_bars, v_vert, p_bound)
        end if
        vp(:, id_vert) = v_vert
          p_ener = 0.0_DOUBLE
          pp_den = 0.0_DOUBLE
          ener_num = 0.0_DOUBLE
          ener_den = 0.0_DOUBLE
          do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
            id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
            id_face = mesh%sub_face(id_sub_face)%mesh_face
            re = mesh%face(id_face)%right_neigh
            if (is_wall(re)) then
              vs_ener(j) = 0.0_DOUBLE
            else
              vs_ener(j) = dot_product(v_vert, mesh%sub_face(id_sub_face)%norm)
            end if
            is_zf = boundary_2d .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE
            if (is_zf) cycle
            rho_l = sol_w_lr(1, 1, j)
            vn_l = dot_product(sol_w_lr(2:4, 1, j), mesh%sub_face(id_sub_face)%norm)
            p_l = sol_w_lr(5, 1, j)
            rho_r = sol_w_lr(1, 2, j)
            vn_r = dot_product(sol_w_lr(2:4, 2, j), mesh%sub_face(id_sub_face)%norm)
            p_r = sol_w_lr(5, 2, j)
            invlam = 1.0_DOUBLE/lambda(1, j) + 1.0_DOUBLE/lambda(2, j)
            p_ener = p_ener + mesh%sub_face(id_sub_face)%area &
              *(p_l/lambda(1, j) + p_r/lambda(2, j) - (vn_r - vn_l))
            pp_den = pp_den + mesh%sub_face(id_sub_face)%area*invlam
            if (re <= 0) cycle
            u_tmp = primit_to_conserv(sol_w_lr(:, 1, j))
            e_l = u_tmp(5)/rho_l
            u_tmp = primit_to_conserv(sol_w_lr(:, 2, j))
            e_r = u_tmp(5)/rho_r
            if (ener_internal) then
              e_l = e_l - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 1, j), sol_w_lr(2:4, 1, j))
              e_r = e_r - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 2, j), sol_w_lr(2:4, 2, j))
            end if
            if (ener_v) then
              mu_l = rho_l*abs(vn_l - vs_ener(j))
              mu_r = rho_r*abs(vn_r - vs_ener(j))
            else
              mu_l = lambda(1, j)
              mu_r = lambda(2, j)
            end if
            ! iso energy solver (the _iso variants): no areas, the counterpart of the iso
            ! velocity solver. Their ener flux enters weighted by A_j warea_j = min_j A_j, a
            ! constant over the node, so it telescopes exactly when e_p carries no areas;
            ! with the areas it leaked energy (6e-10 on Sedov).
            if (ener_iso) then
              w_ener = 1.0_DOUBLE
            else
              w_ener = mesh%sub_face(id_sub_face)%area
            end if
            ener_num = ener_num + w_ener*(mu_l*e_l + mu_r*e_r + (p_l*vn_l - p_r*vn_r))
            ener_den = ener_den + w_ener*(mu_l + mu_r)
          end do
          if (pp_den > 0.0_DOUBLE) p_ener = p_ener/pp_den
          ! ener_pstar: X_p from sum_f w [(pv)*_r - (pv)*_l] = 0 with the star-velocity pressure work,
          ! X_p sum w (mu_l + mu_r) = sum w [mu_l X_l + mu_r X_r - p_p (u*_r - u*_l)]; it needs p_p,
          ! hence this second pass over the sub-faces.
          if (ener_pstar) then
            ener_num = 0.0_DOUBLE
            ener_den = 0.0_DOUBLE
            do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
              id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
              id_face = mesh%sub_face(id_sub_face)%mesh_face
              re = mesh%face(id_face)%right_neigh
              is_zf = boundary_2d .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE
              if (is_zf .or. re <= 0) cycle
              rho_l = sol_w_lr(1, 1, j)
              vn_l = dot_product(sol_w_lr(2:4, 1, j), mesh%sub_face(id_sub_face)%norm)
              p_l = sol_w_lr(5, 1, j)
              rho_r = sol_w_lr(1, 2, j)
              vn_r = dot_product(sol_w_lr(2:4, 2, j), mesh%sub_face(id_sub_face)%norm)
              p_r = sol_w_lr(5, 2, j)
              u_tmp = primit_to_conserv(sol_w_lr(:, 1, j))
              e_l = u_tmp(5)/rho_l
              u_tmp = primit_to_conserv(sol_w_lr(:, 2, j))
              e_r = u_tmp(5)/rho_r
              if (ener_internal) then
                e_l = e_l - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 1, j), sol_w_lr(2:4, 1, j))
                e_r = e_r - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 2, j), sol_w_lr(2:4, 2, j))
              end if
              if (ener_v) then
                mu_l = rho_l*abs(vn_l - vs_ener(j))
                mu_r = rho_r*abs(vn_r - vs_ener(j))
              else
                mu_l = lambda(1, j)
                mu_r = lambda(2, j)
              end if
              if (ener_iso) then
                w_ener = 1.0_DOUBLE
              else
                w_ener = mesh%sub_face(id_sub_face)%area
              end if
              ustar_l = vn_l + (p_l - p_ener)/lambda(1, j)
              ustar_r = vn_r + (p_ener - p_r)/lambda(2, j)
              ener_num = ener_num + w_ener*(mu_l*e_l + mu_r*e_r - p_ener*(ustar_r - ustar_l))
              ener_den = ener_den + w_ener*(mu_l + mu_r)
            end do
          end if
          e_ener = 0.0_DOUBLE
          ener_active = (ener_den > 0.0_DOUBLE)
          if (ener_active) e_ener = ener_num/ener_den

          ! one-sided energy fluxes per sub-face
          do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
            id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
            id_face = mesh%sub_face(id_sub_face)%mesh_face
            re = mesh%face(id_face)%right_neigh
            is_zf = boundary_2d .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE
            pv_l_ener(j) = p_ener*vs_ener(j)
            pv_r_ener(j) = p_ener*vs_ener(j)
            if (re <= 0 .or. is_zf .or. .not. ener_active) cycle
            rho_l = sol_w_lr(1, 1, j)
            vn_l = dot_product(sol_w_lr(2:4, 1, j), mesh%sub_face(id_sub_face)%norm)
            p_l = sol_w_lr(5, 1, j)
            rho_r = sol_w_lr(1, 2, j)
            vn_r = dot_product(sol_w_lr(2:4, 2, j), mesh%sub_face(id_sub_face)%norm)
            p_r = sol_w_lr(5, 2, j)
            u_tmp = primit_to_conserv(sol_w_lr(:, 1, j))
            e_l = u_tmp(5)/rho_l
            u_tmp = primit_to_conserv(sol_w_lr(:, 2, j))
            e_r = u_tmp(5)/rho_r
            if (ener_internal) then
              e_l = e_l - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 1, j), sol_w_lr(2:4, 1, j))
              e_r = e_r - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 2, j), sol_w_lr(2:4, 2, j))
            end if
            if (ener_v) then
              mu_l = rho_l*abs(vn_l - vs_ener(j))
              mu_r = rho_r*abs(vn_r - vs_ener(j))
            else
              mu_l = lambda(1, j)
              mu_r = lambda(2, j)
            end if
            if (ener_pstar) then
              pv_l_ener(j) = p_ener*(vn_l + (p_l - p_ener)/lambda(1, j)) - mu_l*(e_ener - e_l)
              pv_r_ener(j) = p_ener*(vn_r + (p_ener - p_r)/lambda(2, j)) + mu_r*(e_ener - e_r)
            else
              pv_l_ener(j) = p_l*vn_l - mu_l*(e_ener - e_l)
              pv_r_ener(j) = p_r*vn_r + mu_r*(e_ener - e_r)
            end if
          end do
      end if
      ! post-bound histogram with the final lambdas
      do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
        id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        re = mesh%face(id_face)%right_neigh
        is_zf = boundary_2d .and. abs(mesh%sub_face(id_sub_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE
        if (re <= 0 .or. is_zf) then
          ! boundary / z-face left star state after the bound
          if ((ener_lambda_pos .or. ener_eps) .and. .not. mesh%vert(id_vert)%is_ghost) then
            rho_l = sol_w_lr(1, 1, j)
            vn_l = dot_product(sol_w_lr(2:4, 1, j), mesh%sub_face(id_sub_face)%norm)
            u_tmp = primit_to_conserv(sol_w_lr(:, 1, j))
            eps_s = u_tmp(5)/rho_l - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, 1, j), sol_w_lr(2:4, 1, j))
            d_s = sol_w_lr(5, 1, j) - p_ener
            b_s = p_ener*vn_l - pv_l_ener(j)
            ratio = (eps_s + b_s/lam_wave(1, j) - 0.5_DOUBLE*d_s**2/lam_wave(1, j)**2)/eps_s
          end if
          cycle
        end if
        do k_side = 1, 2
          rho_l = sol_w_lr(1, k_side, j)
          vn_l = dot_product(sol_w_lr(2:4, k_side, j), mesh%sub_face(id_sub_face)%norm)
          u_tmp = primit_to_conserv(sol_w_lr(:, k_side, j))
          eps_s = u_tmp(5)/rho_l - 0.5_DOUBLE*dot_product(sol_w_lr(2:4, k_side, j), sol_w_lr(2:4, k_side, j))
          d_s = sol_w_lr(5, k_side, j) - p_ener
          if (k_side == 1) then
            b_s = p_ener*vn_l - pv_l_ener(j)
          else
            b_s = pv_r_ener(j) - p_ener*vn_l
          end if
          ratio = (eps_s + b_s/lam_wave(k_side, j) - 0.5_DOUBLE*d_s**2/lam_wave(k_side, j)**2)/eps_s
          if (mesh%vert(id_vert)%is_ghost) cycle
          ! lambda_wave/lambda_nodal, bins 1 | (1,1.1] | (1.1,2] | (2,10] | (10,100] | >100
          lam_ratio = lam_wave(k_side, j)/lam_cl(k_side, j)
          if (lam_ratio <= 1.0_DOUBLE) then
          else
          end if
          ! sides still negative after the bound, by cause
          ! must stay 0 with the decoupled bound: a negative here is a wiring bug
        end do
      end do
    end if

    !!Compute flux across each sub_face
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face

      le = mesh%face(id_face)%left_neigh
      re = mesh%face(id_face)%right_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if

      select case (scheme_id)
      case (SCHEME_MULTI_POINT)
        if( mesh%vert(id_vert)%is_bound .and. exclude_bound_vert) then
          call two_wave(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
            mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), sl, sr)
        else
          if (is_wall(re) ) then
            call multi_point(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
              mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), 0.0_DOUBLE, &
              lambda(1, j), lambda(2, j), sl, sr)
          else
            call multi_point(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
              mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), &
              dot_product(v_vert, mesh%sub_face(id_sub_face)%norm), &
              lambda(1, j), lambda(2, j), sl, sr)
          end if
        end if
      case (SCHEME_MULTI_POINT_ISO)
        if( mesh%vert(id_vert)%is_bound .and. exclude_bound_vert) then
          call two_wave(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
            mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), sl, sr)
        else
          if (is_wall(re) ) then
            call multi_point(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
              mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), 0.0_DOUBLE, &
              lambda(1, j), lambda(2, j), sl, sr)
          else
            call multi_point(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
              mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), &
              dot_product(v_vert, mesh%sub_face(id_sub_face)%norm), &
              lambda(1, j), lambda(2, j), sl, sr)
          end if
        end if
        call three_wave(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
          mesh%sub_face(id_sub_face)%norm, lr_flux_3w(:, :, j), sl, sr)
        lr_flux(:, :, j) = warea(j) * lr_flux(:, :, j) &
          + (1.0_DOUBLE - warea(j)) * lr_flux_3w(:, :, j)
      case (SCHEME_THREE_WAVE)
        call three_wave(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
          mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), sl, sr)
      case (SCHEME_TWO_WAVE)
        call two_wave(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
          mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), sl, sr)
      case (SCHEME_MODIFIED_THREE_WAVE)
        call modified_three_wave(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
          mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), sl, sr)
      case (SCHEME_THREE_WAVE_ENTHALPY)
        call three_wave_enthalpy(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
          mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), sl, sr)
      case (SCHEME_THREE_WAVE_ENTHALPY2)
        call three_wave_enthalpy2(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
          mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), sl, sr)
      case (SCHEME_MULTI_POINT_ENTHALPY)
        ! Mirrors the SCHEME_MULTI_POINT branch, except that the boundary-vertex fallback is
        ! the enthalpy-preserving 1D solver rather than two_wave: the paper applies the
        ! one-dimensional scheme at boundary nodes ("ensuring conservation, which is critical
        ! for accuracy"), and falling back to two_wave there would reintroduce an enthalpy
        ! defect on exactly the wall cells the property is wanted on.
        if( mesh%vert(id_vert)%is_bound .and. exclude_bound_vert) then
          call three_wave_enthalpy(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
            mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), sl, sr)
        else
          if (is_wall(re) ) then
            call multi_point_enthalpy(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
              mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), 0.0_DOUBLE, &
              lambda(1, j), lambda(2, j), sl, sr)
          else
            call multi_point_enthalpy(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
              mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), &
              dot_product(v_vert, mesh%sub_face(id_sub_face)%norm), &
              lambda(1, j), lambda(2, j), sl, sr)
          end if
        end if
      case (SCHEME_MULTI_POINT_ENER, SCHEME_MULTI_POINT_ENER_V, &
          SCHEME_MULTI_POINT_ENER_ISO, SCHEME_MULTI_POINT_ENER_V_ISO, &
          SCHEME_MULTI_POINT_EPS, SCHEME_MULTI_POINT_EPS_ISO, &
          SCHEME_MULTI_POINT_EPS2, SCHEME_MULTI_POINT_EPS2M, SCHEME_MULTI_POINT_PS, &
        SCHEME_MULTI_POINT_PSW, SCHEME_MULTI_POINT_PV)
        ! multi_point_ener, with the same boundary-vertex fallback and wall v* = 0 as
        ! multi_point; the _iso variants blend with three_wave exactly as multi_point_iso
        if( mesh%vert(id_vert)%is_bound .and. exclude_bound_vert) then
          call two_wave(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
            mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), sl, sr)
        else
          ener_ent_interior = (re > 0) .and. .not. (boundary_2d .and. &
            abs(mesh%sub_face(id_sub_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE)
          if (eps_kind >= 5 .and. .not. ener_ent_interior) then
            ! psw: boundary / z-face sub-face, classic multi_point (one-sided pbar, v* = 0 at a wall)
            call multi_point(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
              mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), vs_ener(j), &
              lam_wave(1, j), lam_wave(2, j), sl, sr)
          else
          call multi_point_ener(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
            mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), vs_ener(j), &
            lam_wave(1, j), lam_wave(2, j), sl, sr, p_ener, pv_l_ener(j), pv_r_ener(j))
          end if
        end if
        if (scheme_id == SCHEME_MULTI_POINT_ENER_ISO .or. scheme_id == SCHEME_MULTI_POINT_ENER_V_ISO .or. &
            scheme_id == SCHEME_MULTI_POINT_EPS_ISO) then
          call three_wave(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
            mesh%sub_face(id_sub_face)%norm, lr_flux_3w(:, :, j), sl, sr)
          lr_flux(:, :, j) = warea(j) * lr_flux(:, :, j) &
            + (1.0_DOUBLE - warea(j)) * lr_flux_3w(:, :, j)
        end if
      case (SCHEME_MULTI_POINT_VILAR)
        call multi_point_vilar(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
          mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), sol_p_vilar, sl, sr)
      case (SCHEME_MULTI_POINT_PRESSURE)
        call multi_point_pressure(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
          mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), p_nodal, &
          lambda(1, j), lambda(2, j), sl, sr)
      case (SCHEME_MULTI_POINT_PRESSURE_PH)
        call multi_point_pressure_ph(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
          mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), p_nodal, &
          lambda(1, j), lambda(2, j), sl, sr)
      case default
        print*,"Unknown scheme"
        error stop
      end select

      ! multi_point_eps diag: the same wave sums with the classic lambdas (no effect on the solution).
      ! The _iso variants take sl, sr from three_wave, which never sees the bound: S_cl = S there.
      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        sum_lambda_vert(lse_loc) = sum_lambda_vert(lse_loc) &
          + mesh%sub_face(id_sub_face)%area*max(0.0_DOUBLE, -sl)
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) &
          + mesh%sub_face(id_sub_face)%area*lr_flux(:, 1, j)

        if (rse > 0) then
        if (mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*max(0.0_DOUBLE, sr)
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) + mesh%sub_face(id_sub_face)%area*lr_flux(:, 2, j)
        end if
        end if
      end if
    end do
  end subroutine compute_rhs_around_vert

  subroutine add_timedisc_and_solve(mesh, sol, rhs, diag, sum_lambda, &
      delta_sol, res, res_vect, num_procs, me, a1, a2, cfl, t, dt)
    use mpi
    use mpi_module
    use linear_solver_module, only: lu_solve_inplace_lapack, eye, eye5
    use sparse_csr_linear_module, only: mpi_norm2
    use ns_global_data_module, only: local_time_step, tmax, ener_theta
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: sol, rhs
    real(kind=DOUBLE), dimension(5, 5, mesh%n_elems), intent(inout) :: diag
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(inout) :: sum_lambda
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: delta_sol
    real(kind=DOUBLE), intent(inout) :: res
    real(kind=DOUBLE), dimension(5), intent(inout) :: res_vect
    integer(kind=ENTIER), intent(in) :: num_procs, me
    real(kind=DOUBLE), intent(in) :: a1, a2 !For BDF2
    real(kind=DOUBLE), intent(in) :: cfl, t
    real(kind=DOUBLE), intent(inout) :: dt

    integer(kind=ENTIER) :: i, mpi_ierr, k
    real(kind=DOUBLE) :: omega, omega_i
    integer, save :: n_culprit = 0
    real(kind=DOUBLE), dimension(5) :: new_sol

    dt = 1e10
    do i = 1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        dt = min(dt, mesh%elem(i)%volume/sum_lambda(i))
      end if
    end do
    dt = cfl * dt

    call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_DOUBLE, &
      MPI_MIN, MPI_COMM_WORLD, mpi_ierr)

    if( .not. local_time_step ) then
      dt = min(tmax - t + 1e-12_DOUBLE, dt)
    end if


    if( .not. local_time_step ) then
      sum_lambda = 0.0_DOUBLE
    end if

    res = 0.0_DOUBLE
    res_vect = 0.0_DOUBLE
    omega = 1.0_DOUBLE
    do i = 1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        do k=1, 5
          res_vect(k) = res_vect(k) + rhs(k, i)**2
        end do
        diag(:, :, i) = diag(:, :, i) &
          + a1*(mesh%elem(i)%volume/dt + sum_lambda(i)/ener_theta)*eye5
        rhs(:, i) = rhs(:, i) &
          + a1*(mesh%elem(i)%volume/dt + sum_lambda(i)/ener_theta)*sol(:, i)
        new_sol = sol(:, i)
        call lu_solve_inplace_lapack(5, diag(:, :, i), new_sol, rhs(:, i))
        delta_sol(:, i) = new_sol - sol(:, i)
        omega_i = max_omega_for_positivity(sol(:, i), a2*delta_sol(:, i))
        if (omega_i < 1.0_DOUBLE) then
          ! culprit of an omega reduction: where it is and what the full update would give
          block
            real(kind=DOUBLE), dimension(5) :: w0, w1
            real(kind=DOUBLE) :: ac
            w0 = conserv_to_primit(sol(:, i))
            w1 = conserv_to_primit(sol(:, i) + a2*delta_sol(:, i))
            if (local_time_step) then
              ac = 1.0_DOUBLE - sum_lambda(i)/(mesh%elem(i)%volume/dt + sum_lambda(i)/ener_theta)
            else
              ac = -1.0_DOUBLE
            end if
            if (n_culprit < 20) print '(a,i3,a,i8,a,3f9.4,a,2es11.3,a,2es11.3,a,f7.4,a,es10.3)', &
              ' [omega] rank', me, ' cell', i, ' x,y,z=', mesh%elem(i)%coord, &
              ' rho,p=', w0(1), w0(5), ' -> ', w1(1), w1(5), ' alpha_c=', ac, ' omega_i=', omega_i
            if (n_culprit < 20 .and. allocated(eps_cell_neg)) print '(a,i8,a,i3,a,i3)', &
              ' [omega] cell', i, '  nodes with eps_p<=0:', eps_cell_neg(i), '  z-face sides eps*<0:', eps_cell_zf(i)
            n_culprit = n_culprit + 1
          end block
        end if
        omega = min(omega, omega_i)
        res = res + norm2(delta_sol(:, i))
      end if
    end do

    if (num_procs > 1) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, res, 1, MPI_DOUBLE, &
        MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, res_vect, 5, MPI_DOUBLE, &
        MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, omega, 1, MPI_DOUBLE, &
        MPI_MIN, MPI_COMM_WORLD, mpi_ierr)
    end if

    do k=1, 5
      res_vect(k) = sqrt(res_vect(k))
    end do

    if( omega <= 1e-6 ) then
      print*, "BAD OMEGA", omega
      error stop
    end if

    if( omega < 1.0_DOUBLE ) then
      if( me == 0 ) print*, "OMEGA", omega
    end if

    if( .not. local_time_step .and. omega < 1.0_DOUBLE ) then
      dt = (1.0_DOUBLE-omega)*dt
    end if

    sol = sol + a2*omega*delta_sol
  end subroutine add_timedisc_and_solve

end module ns_euler_module