module ns_vectorial_diffusion_module
  use precision_module
  use mesh_module
  implicit none

contains
  subroutine build_momentum_diff_nodal_system(mesh, mu_arr, id_vert, &
      N, S, S_tilde_T, B, mu0corr)
    use ns_global_data_module, only: bc_name, bc_type_V, bc_val_V, bc_V_id, BC_V_DIRICHLET, BC_V_NEUMANN
    use linear_solver_module, only: eye3, tensor_product_3, add_block_3, add_block_3_transpose
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mu_arr
    integer(kind=ENTIER), intent(in) :: id_vert

    real(kind=DOUBLE), &
      dimension(3*mesh%vert(id_vert)%n_sub_faces_neigh, &
      3*mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(inout) :: N
    real(kind=DOUBLE), &
      dimension(3*mesh%vert(id_vert)%n_sub_faces_neigh, &
      3*mesh%vert(id_vert)%n_sub_elems_neigh), &
      intent(inout) :: S
    real(kind=DOUBLE), &
      dimension(3*mesh%vert(id_vert)%n_sub_elems_neigh, &
      3*mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(inout) :: S_tilde_T
    real(kind=DOUBLE), dimension(3*mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(inout) :: B
    real(kind=DOUBLE), intent(in) :: mu0corr

    integer(kind=ENTIER) :: l, p
    integer(kind=ENTIER) :: id_sub_elem, id_elem, id_sub_elem_loc
    integer(kind=ENTIER) :: nsfn, nsen
    integer(kind=ENTIER) :: isfl, isfp, isfl_loc, isfp_loc
    integer(kind=ENTIER) :: ifl0, ifp0, ise0
    real(kind=DOUBLE), dimension(3, 3) :: tp_np_nl, blc_3
    real(kind=DOUBLE), dimension(3) :: norml, normp
    real(kind=DOUBLE) :: mu_elem
    logical :: is_dirichlet, is_neumann

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
      ifl0 = 3*(isfl_loc - 1)

      ! Determine face type before processing
      id_sub_elem = mesh%sub_face(isfl)%right_sub_elem_neigh
      id_elem = mesh%sub_face(isfl)%right_elem_neigh
      if (id_elem > 0) then
        is_neumann   = .false.
        is_dirichlet = .false.
      else if (id_elem < 0) then
        is_dirichlet = bc_V_id(-id_elem) == BC_V_DIRICHLET
        is_neumann   = bc_V_id(-id_elem) == BC_V_NEUMANN
      else
        is_neumann   = .true.
        is_dirichlet = .false.
      end if

      if (.not. is_neumann .and. .not. is_dirichlet) then
        !Left sub elem
        id_sub_elem = mesh%sub_face(isfl)%left_sub_elem_neigh
        id_elem = mesh%sub_face(isfl)%left_elem_neigh
        id_sub_elem_loc = mesh%sub_elem(id_sub_elem)%id_loc_around_node
        ise0 = 3*(id_sub_elem_loc - 1)
        norml = mesh%sub_face(isfl)%norm

        mu_elem = mu_arr(id_elem)

        do p = 1, mesh%sub_elem(id_sub_elem)%n_sub_faces
          isfp = mesh%sub_elem(id_sub_elem)%sub_face(p)
          isfp_loc = mesh%sub_face(isfp)%id_loc_around_node
          ifp0 = 3*(isfp_loc - 1)
          if (mesh%sub_face(isfp)%left_sub_elem_neigh == id_sub_elem) then
            normp = mesh%sub_face(isfp)%norm
          else
            normp = -mesh%sub_face(isfp)%norm
          end if
          tp_np_nl = tensor_product_3(normp, norml)
          blc_3 = mesh%sub_face(isfl)%area*mesh%sub_face(isfp)%area &
            *(1.0_DOUBLE/mesh%sub_elem(id_sub_elem)%volume) &
            * (mu_elem*dot_product(norml, normp)*eye3 &
            + (mu_elem - mu0corr)*tp_np_nl &
            + (mu0corr - (2._DOUBLE/3.0_DOUBLE)*mu_elem)*transpose(tp_np_nl))
          call add_block_3_transpose(S_tilde_T, ise0+1, ifl0+1, blc_3)
          call add_block_3(N, ifl0+1, ifp0+1, blc_3)
          call add_block_3(S, ifl0+1, ise0+1, blc_3)
        end do

        !!Right sub elem
        norml = -mesh%sub_face(isfl)%norm
        id_sub_elem = mesh%sub_face(isfl)%right_sub_elem_neigh
        id_elem = mesh%sub_face(isfl)%right_elem_neigh
        id_sub_elem_loc = mesh%sub_elem(id_sub_elem)%id_loc_around_node
        ise0 = 3*(id_sub_elem_loc - 1)

        mu_elem = mu_arr(id_elem)

        do p = 1, mesh%sub_elem(id_sub_elem)%n_sub_faces
          isfp = mesh%sub_elem(id_sub_elem)%sub_face(p)
          isfp_loc = mesh%sub_face(isfp)%id_loc_around_node
          ifp0 = 3*(isfp_loc - 1)
          if (mesh%sub_face(isfp)%left_sub_elem_neigh == id_sub_elem) then
            normp = mesh%sub_face(isfp)%norm
          else
            normp = -mesh%sub_face(isfp)%norm
          end if
          tp_np_nl = tensor_product_3(normp, norml)
          blc_3 = mesh%sub_face(isfl)%area*mesh%sub_face(isfp)%area &
            *(1.0_DOUBLE/mesh%sub_elem(id_sub_elem)%volume) &
            * (mu_elem*dot_product(norml, normp)*eye3 &
            + (mu_elem - mu0corr)*tp_np_nl &
            + (mu0corr - (2._DOUBLE/3.0_DOUBLE)*mu_elem)*transpose(tp_np_nl))
          call add_block_3(N, ifl0+1, ifp0+1, blc_3)
          call add_block_3(S, ifl0+1, ise0+1, blc_3)
          call add_block_3_transpose(S_tilde_T, ise0+1, ifl0+1, blc_3)
        end do
      else if (is_neumann) then
        !Left sub elem
        id_sub_elem = mesh%sub_face(isfl)%left_sub_elem_neigh
        id_elem = mesh%sub_face(isfl)%left_elem_neigh
        id_sub_elem_loc = mesh%sub_elem(id_sub_elem)%id_loc_around_node
        ise0 = 3*(id_sub_elem_loc - 1)
        norml = mesh%sub_face(isfl)%norm

        mu_elem = mu_arr(id_elem)

        do p = 1, mesh%sub_elem(id_sub_elem)%n_sub_faces
          isfp = mesh%sub_elem(id_sub_elem)%sub_face(p)
          isfp_loc = mesh%sub_face(isfp)%id_loc_around_node
          ifp0 = 3*(isfp_loc - 1)
          if (mesh%sub_face(isfp)%left_sub_elem_neigh == id_sub_elem) then
            normp = mesh%sub_face(isfp)%norm
          else
            normp = -mesh%sub_face(isfp)%norm
          end if
          tp_np_nl = tensor_product_3(normp, norml)
          blc_3 = mesh%sub_face(isfl)%area*mesh%sub_face(isfp)%area &
            *(1.0_DOUBLE/mesh%sub_elem(id_sub_elem)%volume) &
            *(dot_product(norml, normp)*eye3 &
            + transpose(tp_np_nl)/3.0_DOUBLE)
          call add_block_3(N, ifl0+1, ifp0+1, blc_3)
          call add_block_3(S, ifl0+1, ise0+1, blc_3)
          blc_3 = mesh%sub_face(isfl)%area*mesh%sub_face(isfp)%area &
            *(1.0_DOUBLE/mesh%sub_elem(id_sub_elem)%volume) &
            * (mu_elem*dot_product(norml, normp)*eye3 &
            + (mu_elem - mu0corr)*tp_np_nl &
            + (mu0corr - (2._DOUBLE/3.0_DOUBLE)*mu_elem)*transpose(tp_np_nl))
          call add_block_3_transpose(S_tilde_T, ise0+1, ifl0+1, blc_3)
        end do
      else if (is_dirichlet) then
        !Left contribution
        id_sub_elem = mesh%sub_face(isfl)%left_sub_elem_neigh
        id_elem = mesh%sub_face(isfl)%left_elem_neigh
        id_sub_elem_loc = mesh%sub_elem(id_sub_elem)%id_loc_around_node
        ise0 = 3*(id_sub_elem_loc - 1)
        norml = mesh%sub_face(isfl)%norm

        mu_elem = mu_arr(id_elem)

        do p = 1, mesh%sub_elem(id_sub_elem)%n_sub_faces
          isfp = mesh%sub_elem(id_sub_elem)%sub_face(p)
          isfp_loc = mesh%sub_face(isfp)%id_loc_around_node
          ifp0 = 3*(isfp_loc - 1)
          if (mesh%sub_face(isfp)%left_sub_elem_neigh == id_sub_elem) then
            normp = mesh%sub_face(isfp)%norm
          else
            normp = -mesh%sub_face(isfp)%norm
          end if
          tp_np_nl = tensor_product_3(normp, norml)
          blc_3 = mesh%sub_face(isfl)%area*mesh%sub_face(isfp)%area &
            *(1.0_DOUBLE/mesh%sub_elem(id_sub_elem)%volume) &
            * (mu_elem*dot_product(norml, normp)*eye3 &
            + (mu_elem - mu0corr)*tp_np_nl &
            + (mu0corr - (2._DOUBLE/3.0_DOUBLE)*mu_elem)*transpose(tp_np_nl))
          call add_block_3_transpose(S_tilde_T, ise0+1, ifl0+1, blc_3)
        end do

        B(3*(isfl_loc - 1) + 1:3*isfl_loc) = &
          mesh%sub_face(isfl)%area*bc_val_V(:, -mesh%sub_face(isfl)%right_elem_neigh)
        S(3*(isfl_loc - 1) + 1:3*isfl_loc, :) = 0.0_DOUBLE
        N(3*(isfl_loc - 1) + 1:3*isfl_loc, :) = 0.0_DOUBLE
        N(3*(isfl_loc - 1) + 1:3*isfl_loc, 3*(isfl_loc - 1) + 1:3*isfl_loc) = &
          mesh%sub_face(isfl)%area*eye3
      end if
    end do
  end subroutine build_momentum_diff_nodal_system

  subroutine build_momentum_diff_mat_into_global_mat_diag(mesh, diag, rhs, sol, mu_arr, mu0corr)
    use linear_solver_module, &
      only: eye3, tensor_product_3, inv_lapack, matmul3_blas, &
      matmul2_blas, matvec2_blas, matvec_blas
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, 5, mesh%n_elems) :: diag
    real(kind=DOUBLE), dimension(5, mesh%n_elems) :: rhs
    real(kind=DOUBLE), dimension(5, mesh%n_elems) :: sol
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mu_arr
    real(kind=DOUBLE), intent(in) :: mu0corr

    integer(kind=ENTIER) :: id_vert, l, p, nsfn, nsen
    integer(kind=ENTIER) :: iel, iep, isel, isep
    integer(kind=ENTIER) :: isfl, isfl_loc, ifl
    real(kind=DOUBLE), dimension(:), allocatable :: B, SB, SNS_Ve, Ve, Vf
    real(kind=DOUBLE), dimension(:, :), allocatable :: N, S, S_tilde_T, SNS
    real(kind=DOUBLE), dimension(:, :, :), allocatable :: blc_3
    real(kind=DOUBLE), dimension(3, 3) :: Spc, tp
    real(kind=DOUBLE), dimension(3) :: norml, SpcVF, VfmVe

    do id_vert = 1, mesh%n_vert
      if( .not. mesh%vert(id_vert)%is_ghost ) then
        nsfn = mesh%vert(id_vert)%n_sub_faces_neigh
        nsen = mesh%vert(id_vert)%n_sub_elems_neigh

        allocate (N(3*nsfn, 3*nsfn))
        allocate (S(3*nsfn, 3*nsen))
        allocate (S_tilde_T(3*nsen, 3*nsfn))
        allocate (B(3*nsfn))
        allocate (SB(3*nsen))
        allocate (SNS(3*nsen, 3*nsen))
        allocate (Ve(3*nsen))
        allocate (Vf(3*nsfn))
        allocate (blc_3(3,3,nsen))
        allocate (SNS_Ve(3*nsen))

        call build_momentum_diff_nodal_system(mesh, mu_arr, id_vert, &
          N, S, S_tilde_T, B, mu0corr)

        ! Precompute Ve before overwriting S and B
        do l = 1, nsen
          isel = mesh%vert(id_vert)%sub_elem_neigh(l)
          iel = mesh%sub_elem(isel)%mesh_elem
          Ve(3*(l-1)+1:3*l) = sol(2:4, iel)/sol(1, iel)
        end do

        ! Compute blc_3
        do p = 1, nsen
          isep = mesh%vert(id_vert)%sub_elem_neigh(p)
          iep = mesh%sub_elem(isep)%mesh_elem
          Spc = 0.0_DOUBLE
          do l = 1, nsfn
            Spc = Spc + S_tilde_T(3*(p-1)+1:3*p, 3*(l-1)+1:3*l)
          end do
          blc_3(:, :, p) = transpose(Spc) / sol(1, iep)
        end do

        ! Compute Vf = S*Ve + B, then invert N in-place, compute all products with dgemm
        call matvec_blas(S, Ve, Vf, 3*nsfn, 3*nsen)
        Vf = Vf + B

        call inv_lapack(3*nsfn, N)
        call matmul3_blas(S_tilde_T, N, S, SNS, 3*nsen, 3*nsfn, 3*nsfn, 3*nsen)
        call matvec2_blas(S_tilde_T, N, B, SB, 3*nsen, 3*nsfn, 3*nsfn)
        ! B is no longer needed; reuse as temp for N^{-1}*(S*Ve+B)
        call matvec_blas(N, Vf, B, 3*nsfn, 3*nsfn)
        Vf = B
        call matvec_blas(SNS, Ve, SNS_Ve, 3*nsen, 3*nsen)

        do l = 1, nsen
          isel = mesh%vert(id_vert)%sub_elem_neigh(l)
          iel = mesh%sub_elem(isel)%mesh_elem
          diag(2:4, 2:4, iel) = diag(2:4, 2:4, iel) &
            - SNS(3*(l-1)+1:3*l, 3*(l-1)+1:3*l)/sol(1, iel) + blc_3(:, :, l)
          rhs(2:4, iel) = rhs(2:4, iel) &
            + SNS_Ve(3*(l-1)+1:3*l) &
            - matmul(SNS(3*(l-1)+1:3*l, 3*(l-1)+1:3*l), Ve(3*(l-1)+1:3*l)) &
            + SB(3*(l-1)+1:3*l)
        end do

        !Add viscous work (Vf = N^{-1}*(S*Ve+B) already computed)
        do p = 1, nsen
          isep = mesh%vert(id_vert)%sub_elem_neigh(p)
          iep = mesh%sub_elem(isep)%mesh_elem
          Spc = 0.0_DOUBLE
          do l = 1, mesh%sub_elem(isep)%n_sub_faces
            isfl = mesh%sub_elem(isep)%sub_face(l)
            ifl = mesh%sub_face(isfl)%mesh_face
            isfl_loc = mesh%sub_face(isfl)%id_loc_around_node
            if (mesh%sub_face(isfl)%left_sub_elem_neigh == isep) then
              norml = mesh%sub_face(isfl)%norm
            else
              norml = -mesh%sub_face(isfl)%norm
            end if
            VfmVe = Vf(3*(isfl_loc-1)+1:3*isfl_loc)- Ve(3*(p-1)+1:3*p)
            tp = tensor_product_3(VfmVe, norml)
            Spc = Spc + mesh%sub_face(isfl)%area*(tp + transpose(tp) &
              - (2.0_DOUBLE/3.0_DOUBLE)*dot_product(VfmVe, norml)*eye3)
          end do
          Spc = Spc * mu_arr(iep)/mesh%sub_elem(isep)%volume
          do l = 1, mesh%sub_elem(isep)%n_sub_faces
            isfl = mesh%sub_elem(isep)%sub_face(l)
            ifl = mesh%sub_face(isfl)%mesh_face
            isfl_loc = mesh%sub_face(isfl)%id_loc_around_node
            if (mesh%sub_face(isfl)%left_sub_elem_neigh == isep) then
              norml = mesh%sub_face(isfl)%norm
            else
              norml = -mesh%sub_face(isfl)%norm
            end if
            SpcVF = matmul(Spc, Vf(3*(isfl_loc-1)+1:3*isfl_loc))
            rhs(5, iep) = rhs(5, iep) &
              + mesh%sub_face(isfl)%area*dot_product(SpcVF, norml)
          end do
        end do

        deallocate (N, S, SNS, B, S_tilde_T, Ve, Vf, SB, SNS_Ve, blc_3)
      end if
    end do
  end subroutine build_momentum_diff_mat_into_global_mat_diag

  subroutine compute_mu(mesh, sol, num_procs, mu_arr, mu0corr)
    use mpi
    use ns_global_data_module, only: use_sutherland, mu_p
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    integer(kind=ENTIER), intent(in) :: num_procs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(out) :: mu_arr
    real(kind=DOUBLE), intent(out) :: mu0corr

    integer(kind=ENTIER) :: id_elem, mpi_ierr

    if (use_sutherland) then
      mu0corr = 1e12_DOUBLE
      do id_elem = 1, mesh%n_elems
        mu_arr(id_elem) = mu_sutherland(sol(:, id_elem))
        mu0corr = min(mu0corr, mu_arr(id_elem))
      end do
      if (num_procs > 1) call MPI_ALLREDUCE(MPI_IN_PLACE, mu0corr, 1, &
        MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, mpi_ierr)
    else
      mu_arr  = mu_p
      mu0corr = mu_p
    end if
  end subroutine compute_mu

  pure function mu_sutherland(sol)
    use ns_global_data_module, only: mu0, T0, C0
    use ns_euler_primitives_module, only: temp_u
    implicit none
    real(kind=DOUBLE), dimension(5), intent(in) :: sol

    real(kind=DOUBLE) :: mu_sutherland, T

    T = temp_u(sol)
    mu_sutherland = mu0*(T/T0)**(1.5_DOUBLE)*(T0 + C0)/(T + C0)
  end function mu_sutherland
end module ns_vectorial_diffusion_module