module ns_euler_module
  use precision_module
  use mesh_module
  use ns_euler_primitives_module
  use ns_euler_recon_module
  use ns_euler_zb_module
  use ns_euler_exact_sol_module
  use ns_euler_rs_module
  implicit none

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
        sol(5, i) = 1e-16_DOUBLE/(gamma - 1.0_DOUBLE)
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
    use ns_global_data_module, only: bc_style, scheme, exclude_bound_vert, &
      scheme_id, scheme_adv_id, scheme_lag_id, &
      SCHEME_ZB, &
      SCHEME_ADV_AR1D, SCHEME_ADV_AM, SCHEME_ADV_AMISO, &
      SCHEME_ADV_ARMD, SCHEME_ADV_ARMDU, SCHEME_ADV_ARMDM, &
      SCHEME_ADV_ARMDMAT, SCHEME_ADV_ARMDUMAT, SCHEME_ADV_ARMDMMAT, &
      SCHEME_LAG_LS, SCHEME_LAG_LSU, SCHEME_LAG_LSM, SCHEME_ADV_ARMDWIP, &
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
        case (SCHEME_ADV_ARMDM)
          call compute_rhs_around_vert_ARMDM(mesh, sol, grad, &
            nsen, sum_lambda_vert, &
            flux_sum_vert, second_order, id_vert, vp, mat_h_p)
        case (SCHEME_ADV_ARMDMAT)
          call compute_rhs_around_vert_ARMDMAT(mesh, sol, grad, &
            nsen, sum_lambda_vert, &
            flux_sum_vert, second_order, id_vert, vp, h_p)
        case (SCHEME_ADV_ARMDUMAT)
          call compute_rhs_around_vert_ARMDUMAT(mesh, sol, grad, &
            nsen, sum_lambda_vert, &
            flux_sum_vert, second_order, id_vert, vp)
        case (SCHEME_ADV_ARMDMMAT)
          call compute_rhs_around_vert_ARMDMMAT(mesh, sol, grad, &
            nsen, sum_lambda_vert, &
            flux_sum_vert, second_order, id_vert, vp, mat_h_p)
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
        case (SCHEME_LAG_LSM)
          call compute_rhs_around_vert_LSM(mesh, sol, grad, &
            nsen, flux_sum_vert, &
            id_vert, vp, mat_h_p, second_order)
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
    else
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
    end if
    deallocate(sum_lambda_vert, flux_sum_vert)
  end subroutine compute_euler_flux

  subroutine compute_rhs_around_vert(mesh, sol, grad, &
      nsfn, nsen, v_vert, sum_lambda_vert, flux_sum_vert, &
      second_order, id_vert, vp)
    use ns_global_data_module, only: bc_style, scheme, exclude_bound_vert, &
      scheme_id, &
      SCHEME_MULTI_POINT, SCHEME_MULTI_POINT_ISO, SCHEME_MULTI_POINT_PRESSURE, &
      SCHEME_THREE_WAVE, SCHEME_TWO_WAVE, SCHEME_MODIFIED_THREE_WAVE
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
    rse_loc = 0
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
    case (SCHEME_MULTI_POINT)
      lambda(:, :) = 0.0_DOUBLE
      if( mesh%vert(id_vert)%is_bound .and. exclude_bound_vert) then
        v_vert = 0.0_DOUBLE
        vp(:, id_vert) = v_vert
      else
        call compute_lambdas_and_solve_nodal_velocity(mesh, id_vert, &
          sol_w_lr, lambda, v_bars, v_vert, p_bound)
        vp(:, id_vert) = v_vert
      end if
    case (SCHEME_MULTI_POINT_ISO)
      lambda(:, :) = 0.0_DOUBLE
      if( mesh%vert(id_vert)%is_bound .and. exclude_bound_vert) then
        v_vert = 0.0_DOUBLE
        vp(:, id_vert) = v_vert
      else
        call compute_lambdas_and_solve_nodal_velocity_iso(mesh, id_vert, &
          sol_w_lr, lambda, v_bars, v_vert, p_bound)
        vp(:, id_vert) = v_vert
      end if
    case (SCHEME_MULTI_POINT_PRESSURE)
      lambda(:, :) = 0.0_DOUBLE
      p_nodal = 0.0_DOUBLE
      call compute_lambdas_and_solve_nodal_pressure(mesh, id_vert, &
        sol_w_lr, lambda, p_bars, p_nodal)
    end select

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
      case (SCHEME_MULTI_POINT_PRESSURE)
        call multi_point_pressure(sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
          mesh%sub_face(id_sub_face)%norm, lr_flux(:, :, j), p_nodal, &
          lambda(1, j), lambda(2, j), sl, sr)
      case default
        print*,"Unknown scheme"
        error stop
      end select

      if (mesh%sub_elem(lse)%mesh_vert == id_vert) then
        sum_lambda_vert(lse_loc) = sum_lambda_vert(lse_loc) &
          + mesh%sub_face(id_sub_face)%area*max(0.0_DOUBLE, -sl)
        flux_sum_vert(:, lse_loc) = flux_sum_vert(:, lse_loc) &
          + mesh%sub_face(id_sub_face)%area*lr_flux(:, 1, j)

        if (rse > 0 .and. mesh%sub_elem(rse)%mesh_vert == id_vert) then
          sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) &
            + mesh%sub_face(id_sub_face)%area*max(0.0_DOUBLE, sr)
          flux_sum_vert(:, rse_loc) = flux_sum_vert(:, rse_loc) + mesh%sub_face(id_sub_face)%area*lr_flux(:, 2, j)
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
    use ns_global_data_module, only: local_time_step, tmax
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
    real(kind=DOUBLE) :: omega
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
          + a1*(mesh%elem(i)%volume/dt + sum_lambda(i))*eye5
        rhs(:, i) = rhs(:, i) &
          + a1*(mesh%elem(i)%volume/dt + sum_lambda(i))*sol(:, i)
        new_sol = sol(:, i)
        call lu_solve_inplace_lapack(5, diag(:, :, i), new_sol, rhs(:, i))
        delta_sol(:, i) = new_sol - sol(:, i)
        omega = min(omega, max_omega_for_positivity(sol(:, i), a2*delta_sol(:, i)))
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