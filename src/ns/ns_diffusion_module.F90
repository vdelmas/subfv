module ns_diffusion_module
  use precision_module
  use mesh_module
  implicit none

contains
  subroutine build_diff_nodal_system(mesh, id_vert, N, S, S_tilde_T, B, mu_arr)
    use ns_global_data_module, &
      only: bc_name, bc_type_T, bc_val_T, bc_T_id, BC_T_DIRICHLET
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_faces_neigh, mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(inout) :: N
    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_faces_neigh, mesh%vert(id_vert)%n_sub_elems_neigh), &
      intent(inout) :: S
    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_elems_neigh, mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(inout) :: S_tilde_T
    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_faces_neigh), intent(inout) :: B
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mu_arr

    integer(kind=ENTIER) :: l, p
    integer(kind=ENTIER) :: id_sub_elem, id_elem, id_sub_elem_loc
    integer(kind=ENTIER) :: nsfn, nsen
    integer(kind=ENTIER) :: isfl, isfp, isfl_loc, isfp_loc, ifl, ifp
    real(kind=DOUBLE) :: kap, blc_1
    real(kind=DOUBLE), dimension(3) :: norml, normp

    nsfn = mesh%vert(id_vert)%n_sub_faces_neigh
    nsen = mesh%vert(id_vert)%n_sub_elems_neigh

    N = 0.0_DOUBLE
    S = 0.0_DOUBLE
    S_tilde_T = 0.0_DOUBLE
    B = 0.0_DOUBLE

    !Build NT = ST + B
    do l = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      isfl = mesh%vert(id_vert)%sub_face_neigh(l)
      isfl_loc = mesh%sub_face(isfl)%id_loc_around_node
      ifl = mesh%sub_face(isfl)%mesh_face

      !Left sub elem
      id_sub_elem = mesh%sub_face(isfl)%left_sub_elem_neigh
      id_elem = mesh%face(ifl)%left_neigh
      id_sub_elem_loc = mesh%sub_elem(id_sub_elem)%id_loc_around_node
      norml = mesh%sub_face(isfl)%norm

      kap = kappa(mu_arr(id_elem))

      !Add terms in N, S, and S_tilde
      do p = 1, mesh%sub_elem(id_sub_elem)%n_sub_faces
        isfp = mesh%sub_elem(id_sub_elem)%sub_face(p)
        isfp_loc = mesh%sub_face(isfp)%id_loc_around_node
        ifp = mesh%sub_face(isfp)%mesh_face
        if (mesh%sub_face(isfp)%left_sub_elem_neigh == id_sub_elem) then
          normp = mesh%sub_face(isfp)%norm
        else
          normp = -mesh%sub_face(isfp)%norm
        end if
        blc_1 = kap*mesh%sub_face(isfl)%area*mesh%sub_face(isfp)%area &
          *(1.0_DOUBLE/mesh%sub_elem(id_sub_elem)%volume) &
          *dot_product(normp, norml)
        N(isfl_loc, isfp_loc) = N(isfl_loc, isfp_loc) + blc_1
        S(isfl_loc, id_sub_elem_loc) = S(isfl_loc, id_sub_elem_loc) + blc_1
        S_tilde_T(id_sub_elem_loc, isfl_loc) = &
          S_tilde_T(id_sub_elem_loc, isfl_loc) + blc_1
      end do

      !Right sub elem
      norml = -mesh%sub_face(isfl)%norm
      id_sub_elem = mesh%sub_face(isfl)%right_sub_elem_neigh
      id_elem = mesh%face(ifl)%right_neigh
      if (id_sub_elem > 0) then
        id_sub_elem_loc = mesh%sub_elem(id_sub_elem)%id_loc_around_node

        kap = kappa(mu_arr(id_elem))

        !Add terms in N, S, and S_tilde
        do p = 1, mesh%sub_elem(id_sub_elem)%n_sub_faces
          isfp = mesh%sub_elem(id_sub_elem)%sub_face(p)
          isfp_loc = mesh%sub_face(isfp)%id_loc_around_node
          ifp = mesh%sub_face(isfp)%mesh_face
          if (mesh%sub_face(isfp)%left_sub_elem_neigh == id_sub_elem) then
            normp = mesh%sub_face(isfp)%norm
          else
            normp = -mesh%sub_face(isfp)%norm
          end if
          blc_1 = kap*mesh%sub_face(isfl)%area*mesh%sub_face(isfp)%area &
            *(1.0_DOUBLE/mesh%sub_elem(id_sub_elem)%volume) &
            *dot_product(normp, norml)
          N(isfl_loc, isfp_loc) = N(isfl_loc, isfp_loc) + blc_1
          S(isfl_loc, id_sub_elem_loc) = S(isfl_loc, id_sub_elem_loc) + blc_1
          S_tilde_T(id_sub_elem_loc, isfl_loc) = &
            S_tilde_T(id_sub_elem_loc, isfl_loc) + blc_1
        end do
      else
        id_elem = mesh%face(ifl)%right_neigh
        if( id_elem < 0 ) then
          if (bc_T_id(-id_elem) == BC_T_DIRICHLET) then
            B(isfl_loc) = mesh%sub_face(isfl)%area*bc_val_T(-id_elem)
            S(isfl_loc, :) = 0.0_DOUBLE
            N(isfl_loc, :) = 0.0_DOUBLE
            N(isfl_loc, isfl_loc) = mesh%sub_face(isfl)%area
          end if
        end if
      end if
    end do
  end subroutine build_diff_nodal_system

  subroutine build_diff_mat_into_global_mat_diag(mesh, diag, rhs, sol, mu_arr)
    use ns_global_data_module, only: bc_name, bc_type_T, bc_val_T
    use ns_euler_primitives_module, only: temp_u
    use linear_solver_module, only: inv_lapack, &
    matmul3_blas, matvec2_blas, matvec_blas
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, 5, mesh%n_elems), intent(inout) :: diag
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mu_arr

    integer(kind=ENTIER) :: l, p
    integer(kind=ENTIER) :: id_vert
    integer(kind=ENTIER) :: nsfn, nsen
    integer(kind=ENTIER) :: iel, iep, isel, isep
    real(kind=DOUBLE), dimension(:), allocatable :: B, SB, theta, blc_1
    real(kind=DOUBLE), dimension(:, :), allocatable :: N, S, S_tilde_T, SNS
    real(kind=DOUBLE), dimension(:, :), allocatable :: mat_P

    do id_vert = 1, mesh%n_vert
      if ( .not. mesh%vert(id_vert)%is_ghost ) then
        nsfn = mesh%vert(id_vert)%n_sub_faces_neigh
        nsen = mesh%vert(id_vert)%n_sub_elems_neigh

        allocate (N(nsfn, nsfn))
        allocate (S(nsfn, nsen))
        allocate (S_tilde_T(nsen, nsfn))
        allocate (B(nsfn))
        allocate (SNS(nsen, nsen))
        allocate (SB(nsen))
        allocate (theta(nsen))
        allocate (mat_P(5, nsen))
        allocate (blc_1(nsen))

        call build_diff_nodal_system(mesh, id_vert, N, S, S_tilde_T, B, mu_arr)

        call inv_lapack(nsfn, N)
        call matmul3_blas(S_tilde_T, N, S, SNS, nsen, nsfn, nsfn, nsen)
        call matvec2_blas(S_tilde_T, N, B, SB, nsen, nsfn, nsfn)

        !Precompute
        do l = 1, nsen
          isel = mesh%vert(id_vert)%sub_elem_neigh(l)
          iel = mesh%sub_elem(isel)%mesh_elem
          theta(l) = temp_u(sol(:, iel))
          mat_P(:, l) = dtheta_dU(sol(:, iel))
          blc_1(l) = sum(S_tilde_T(l, :))
        end do

        !Scatter into matrix
        do l = 1, nsen
          isel = mesh%vert(id_vert)%sub_elem_neigh(l)
          iel = mesh%sub_elem(isel)%mesh_elem
          do p = 1, nsen
            isep = mesh%vert(id_vert)%sub_elem_neigh(p)
            iep = mesh%sub_elem(isep)%mesh_elem
            !call add_block_to_bcsr_thread_safe(mat, iel, iep, -SNS(l,p)*mat_P)
            rhs(5, iel) = rhs(5, iel) + SNS(l, p)*theta(p)
          end do
          !call add_block_to_bcsr_thread_safe(mat, iel, iel, mat_Q*mat_P*blc_1)
          diag(5, :, iel) = diag(5, :, iel) + (blc_1(l)-SNS(l,l))*mat_P(:, l)
          rhs(5, iel) = rhs(5, iel) &
            - SNS(l, l)*dot_product(mat_P(:, l),sol(:, iel)) &
            + (blc_1(l)-SNS(l,l))*dot_product(mat_P(:, l),sol(:, iel)) &
            - blc_1(l)*theta(l) + SB(l)
        end do

        deallocate (N, S, S_tilde_T, B, SNS, SB, theta, mat_P, blc_1)
      end if
    end do
  end subroutine build_diff_mat_into_global_mat_diag

  pure function kappa(mu)
    use ns_global_data_module, only: Cv_p, Prandtl, gamma
    implicit none

    real(kind=DOUBLE), intent(in) :: mu
    real(kind=DOUBLE) :: kappa

    kappa = mu*Cv_p*gamma/Prandtl
  end function kappa

  pure function dtheta_dU(u) result(mat_P)
    use ns_global_data_module, only: Cv_p
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: u
    real(kind=DOUBLE), dimension(5) :: mat_P

    mat_P(1)   = (norm2(u(2:4)/u(1))**2 - u(5)/u(1))/(u(1)*Cv_p)
    mat_P(2:4) = -u(2:4)/(u(1)**2*Cv_p)
    mat_P(5)   = 1.0_DOUBLE/(u(1)*Cv_p)
  end function dtheta_dU
end module ns_diffusion_module