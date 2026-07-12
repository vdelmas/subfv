module ns_io_module
  use mpi_module
  use precision_module
  use mesh_module
  implicit none

contains
  subroutine write_sol(mesh, sol, delta_sol, grad, vp, h_p, mat_h_p, me, iaff)
    use ns_global_data_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol, delta_sol
    integer(kind=ENTIER), intent(in) :: me, iaff

    real(kind=DOUBLE), dimension(:, :, :), intent(inout) :: grad
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(inout) :: h_p
    real(kind=DOUBLE), dimension(3, 3, mesh%n_vert), intent(in) :: mat_h_p

    character(len=255) :: iaff_str, me_str

    write (me_str, *) me
    write (iaff_str, *) iaff
    call write_sol_vtu(mesh, &
      trim(adjustl(me_str))//"_output_"//trim(adjustl(iaff_str))//".vtu", &
      sol, delta_sol, grad, vp, h_p, mat_h_p)
    call write_sol_meta_pvtu("output_"//trim(adjustl(iaff_str))//".pvtu", iaff_str)

    if (plot_solution_dat) call write_sol_dat(mesh, &
      "output_"//trim(adjustl(iaff_str))//".dat", sol)

    if (compute_coeffs) then
      call write_sol_meta_pvtu_coeffs("coeffs_"//trim(adjustl(iaff_str))//".pvtu", &
        iaff_str)
      call write_coeffs_vtu(mesh, &
        trim(adjustl(me_str))//"_coeffs_"//trim(adjustl(iaff_str))//".vtu", sol)
    end if

    if( write_cell_size ) then
      call write_cell_size_dat("cell_sizes_"//trim(adjustl(iaff_str))//".dat", &
        mesh, sol, grad)
      !call write_cell_size_mmg("cell_sizes_"//trim(adjustl(iaff_str))//".sol", &
      !  mesh, sol, grad, all_nodal_grad)
    end if

    if( me == 0 ) print *, "[+] Sol written ", iaff
  end subroutine write_sol

  subroutine write_sol_dat(mesh, filename, sol)
    use mpi
    use ns_global_data_module
    use ns_euler_exact_sol_module, only: sol_isentropic_vortex
    use ns_euler_primitives_module, only: conserv_to_primit
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    character(len=*) :: filename

    integer(kind=ENTIER) :: i, vtkout
    integer(kind=ENTIER) :: me, num_procs, mpi_ierr, idproc
    real(kind=DOUBLE), dimension(5) :: w


    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    if( me == 0 ) then
      open (newunit=vtkout, file=trim(adjustl(filename)), status="unknown")
      write(vtkout, *) ""
      close(vtkout)
    end if

    do idproc=1, num_procs
      if( me == idproc - 1 ) then
        open (newunit=vtkout, file=trim(adjustl(filename)), position="append")
        do i = 1, mesh%n_elems
          if (.not. mesh%elem(i)%is_ghost) then
            if ((mesh%elem(i)%coord(1) < xmax_dat .and. &
              mesh%elem(i)%coord(1) > xmin_dat) .and. &
              (mesh%elem(i)%coord(2) < ymax_dat .and. &
              mesh%elem(i)%coord(2) > ymin_dat) .and. &
              (mesh%elem(i)%coord(3) < zmax_dat .and. &
              mesh%elem(i)%coord(3) > zmin_dat)) then

              w = conserv_to_primit(sol(:, i))
              write (vtkout, *) mesh%elem(i)%coord, sol(:, i), w(:)
            end if
          end if
        end do
        close (vtkout)
      end if
      call mpi_barrier(mpi_comm_world, mpi_ierr)
    end do
  end subroutine write_sol_dat

  subroutine write_coeffs_vtu(mesh, filename, sol)
    use ns_global_data_module, only: compute_coeffs, coeffs_surf, &
      pinf, rhoinf, vinf, T0, bc_val_T, use_sutherland, mu_p
    use mpi
    use mpi_module
    use ns_euler_primitives_module, only: conserv_to_primit, temp_u
    use ns_diffusion_module, only: build_diff_nodal_system, kappa
    use ns_vectorial_diffusion_module, only: mu_sutherland
    use linear_solver_module, only: inv_lapack, tensor_product
    implicit none

    type(mesh_type), intent(in) :: mesh
    character(len=*), intent(in) :: filename
    real(kind=DOUBLE), dimension(5, mesh%n_elems) :: sol

    integer(kind=ENTIER) :: fn
    integer(kind=ENTIER) :: i, j, k
    integer(kind=ENTIER) :: id_elem, id_face, id_sub_face, id_vert
    integer(kind=ENTIER) :: id_sub_elem, id_sub_elem_loc, id_sub_face_k_loc
    integer(kind=ENTIER) :: id_elem_k, id_face_k, id_sub_face_k, id_sub_elem_k
    integer(kind=ENTIER) :: n_bdy_faces
    real(kind=DOUBLE), dimension(5) :: w
    integer(kind=ENTIER), dimension(:), allocatable :: face_ids
    real(kind=DOUBLE), dimension(:), allocatable :: cp, q

    integer(kind=ENTIER) :: nsfn, nsen
    real(kind=DOUBLE), dimension(:, :), allocatable :: N, S, S_tilde
    real(kind=DOUBLE), dimension(:), allocatable :: B, Tf, Te, mu_arr
    real(kind=DOUBLE), dimension(3) :: qpc
    real(kind=DOUBLE), dimension(3, 3) :: mat
    real(kind=DOUBLE) :: kap

    integer(kind=ENTIER), dimension(:), allocatable :: id_vert_no_ghost
    real(kind=DOUBLE) :: dmin, dmax
    integer(kind=ENTIER) :: me, num_procs, mpi_ierr, sm
    integer(kind=ENTIER) :: n_interior_elems, n_interior_vert

    !Get the face ids from the boundary coeffs_surf
    n_bdy_faces = 0
    do i = 1, mesh%n_faces
      if (-mesh%face(i)%right_neigh == coeffs_surf) n_bdy_faces = n_bdy_faces + 1
    end do
    allocate (face_ids(n_bdy_faces))
    n_bdy_faces = 0
    do i = 1, mesh%n_faces
      if (-mesh%face(i)%right_neigh == coeffs_surf) then
        n_bdy_faces = n_bdy_faces + 1
        face_ids(n_bdy_faces) = i
      end if
    end do

    allocate (cp(n_bdy_faces))
    !Compute Cp for the face
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      id_elem = mesh%face(id_face)%left_neigh
      w = conserv_to_primit(sol(:, id_elem))
      cp(i) = (w(5) - pinf)/(0.5_DOUBLE*rhoinf*vinf**2)
    end do

    allocate (mu_arr(mesh%n_elems))
    if (use_sutherland) then
      do i = 1, mesh%n_elems
        mu_arr(i) = mu_sutherland(sol(:, i))
      end do
    else
      mu_arr = mu_p
    end if

    !Compute Q for the subfaces
    allocate (q(n_bdy_faces))
    q = 0.0_DOUBLE
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      id_elem = mesh%face(id_face)%left_neigh
      do j = 1, mesh%face(id_face)%n_sub_faces
        id_sub_face = mesh%face(id_face)%sub_face(j)
        id_vert = mesh%sub_face(id_sub_face)%mesh_vert
        id_sub_elem = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
        id_sub_elem_loc = mesh%sub_elem(id_sub_elem)%id_loc_around_node

        nsfn = mesh%vert(id_vert)%n_sub_faces_neigh
        nsen = mesh%vert(id_vert)%n_sub_elems_neigh

        allocate (N(nsfn, nsfn))
        N = 0.0_DOUBLE
        allocate (S(nsfn, nsen))
        S = 0.0_DOUBLE
        allocate (S_tilde(nsfn, nsen))
        S_tilde = 0.0_DOUBLE
        allocate (B(nsfn))
        B = 0.0_DOUBLE

        allocate (Te(nsen))
        Te = 0.0_DOUBLE
        allocate (Tf(nsfn))
        Tf = 0.0_DOUBLE

        call build_diff_nodal_system(mesh, id_vert, N, S, S_tilde, B, mu_arr)
        call inv_lapack(nsfn, N)

        do k = 1, mesh%vert(id_vert)%n_sub_elems_neigh
          id_sub_elem_k = mesh%vert(id_vert)%sub_elem_neigh(k)
          id_elem_k = mesh%sub_elem(id_sub_elem_k)%mesh_elem
          Te(k) = temp_u(sol(:, id_elem_k))
        end do

        Tf = matmul(N, matmul(S, Te) + B)

        mat = 0.0_DOUBLE
        qpc = 0.0_DOUBLE
        do k = 1, mesh%sub_elem(id_sub_elem)%n_sub_faces
          id_sub_face_k = mesh%sub_elem(id_sub_elem)%sub_face(k)
          id_face_k = mesh%sub_face(id_sub_face_k)%mesh_face
          id_sub_face_k_loc = mesh%sub_face(id_sub_face_k)%id_loc_around_node

          kap = kappa(mu_arr(id_elem))

          qpc = qpc - kap*(Tf(id_sub_face_k_loc) - Te(id_sub_elem_loc))*mesh%face(id_face_k)%norm &
            *mesh%sub_face(id_sub_face_k)%area/mesh%sub_elem(id_sub_elem)%volume
          mat = mat + tensor_product(mesh%face(id_face_k)%norm, mesh%face(id_face_k)%norm)
        end do

        call inv_lapack(3, mat)
        qpc = matmul(mat, qpc)
        q(i) = q(i) + mesh%sub_face(id_sub_face)%area*dot_product(qpc, mesh%face(id_face)%norm)
        deallocate (N, S, S_tilde, B, Te, Tf)
      end do
    end do

    !Write VTU Data
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    allocate(id_vert_no_ghost(mesh%n_vert))
    id_vert_no_ghost = 0
    n_interior_vert = 0
    do i=1, mesh%n_vert
      if ( .not. mesh%vert(i)%is_ghost ) then
        n_interior_vert = n_interior_vert + 1
        id_vert_no_ghost(i) = n_interior_vert
      end if
    end do

    n_interior_elems = 0
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) n_interior_elems = n_interior_elems + 1
    end do

    open(newunit=fn, file=trim(adjustl(filename)))

    write(fn, *) "<VTKFile type='UnstructuredGrid' version='1.0' &
      &byte_order='LittleEndian' header_type='UInt64'>"
    write(fn, *) "<UnstructuredGrid>\n<Piece NumberOfPoints='", n_interior_vert, &
      &"' NumberOfCells='", n_interior_elems, "'>"

    dmin = 1e100_DOUBLE
    dmax = -1e100_DOUBLE
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        do j=1, 3
          if( mesh%vert(i)%coord(j) > dmax ) dmax = mesh%vert(i)%coord(j)
          if( mesh%vert(i)%coord(j) < dmin ) dmin = mesh%vert(i)%coord(j)
        end do
      end if
    end do

    write(fn, *) "<Points>"
    write(fn, *) " <DataArray type='Float64' Name='Points' NumberOfComponents='3' &
      & format='ascii' RangeMin='", dmin, "' RangeMax='", dmax, "'>"

    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        write(fn, *) mesh%vert(i)%coord(:)
      end if
    end do

    write(fn, *) "</DataArray>"
    write(fn, *) "</Points>"

    write(fn, *) "<Cells>"
    write(fn, *) "<DataArray type='Int64' Name='connectivity' format='ascii' &
      &RangeMin='0' RangeMax='", mesh%n_vert-1, "'>"
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        write(fn, *) id_vert_no_ghost(mesh%face(id_face)%vert) - 1
      end if
    end do
    write(fn, *) "</DataArray>"

    sm = 0
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        sm = sm + mesh%face(id_face)%n_vert
      end if
    end do

    write(fn, *) "<DataArray type='Int64' Name='offsets' format='ascii' &
      &RangeMin='", 0, "' RangeMax='", sm, "'>"

    sm = 0
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        sm = sm + mesh%face(id_face)%n_vert
        write(fn, *) sm
      end if
    end do

    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='UInt8' Name='types' format='ascii' &
      &RangeMin='7' RangeMax='7'>"

    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        write(fn, *) 7 !Type for polygon
      end if
    end do

    write(fn, *) "</DataArray>"
    write(fn, *) "</Cells>"

    write(fn, *) "<CellData>"

    write(fn, *) "<DataArray type='Float64' Name='Centroid' format='ascii' NumberOfComponents='3'>"
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        write (fn, *) mesh%face(id_face)%coord
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Cp' format='ascii' NumberOfComponents='1'>"
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        write (fn, *) cp(i)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='q' format='ascii' NumberOfComponents='1'>"
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        write (fn, *) q(i)/mesh%face(id_face)%area
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Density' format='ascii' NumberOfComponents='1'>"
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        id_elem = mesh%face(id_face)%left_neigh
        write (fn, *) sol(1, id_elem)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Temperature' format='ascii' NumberOfComponents='1'>"
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        id_elem = mesh%face(id_face)%left_neigh
        write (fn, *) temp_u(sol(:, id_elem))
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Pressure' format='ascii' NumberOfComponents='1'>"
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        id_elem = mesh%face(id_face)%left_neigh
        w = conserv_to_primit(sol(:, id_elem))
        write (fn, *) w(5)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Theta3' format='ascii' NumberOfComponents='1'>"
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        id_elem = mesh%face(id_face)%left_neigh
        write (fn, *) atan2(-mesh%elem(id_elem)%coord(2), -mesh%elem(id_elem)%coord(1))
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "</CellData>"

    write(fn, *) "</Piece>"
    write(fn, *) "</UnstructuredGrid>"
    write(fn, *) "</VTKFile>"
    close(fn)
  end subroutine write_coeffs_vtu

  subroutine write_sol_vtu(mesh, filename, sol, delta_sol, grad, vp, h_p, mat_h_p)
    use mpi
    use, intrinsic :: ieee_arithmetic
    use ns_global_data_module
    use ns_euler_exact_sol_module, only: sol_isentropic_vortex
    use ns_euler_primitives_module, only: conserv_to_primit, sound_speed_w, &
      mach_u, temp_u
    use ns_euler_recon_module, only: compute_nodal_grad
    use ns_euler_zb_module, only: compute_nodal_pressure_LS, compute_nodal_velocity_LS
    use ns_euler_recon_module, only: omega_ducros
    use ns_vectorial_diffusion_module, only: mu_sutherland
    use linear_solver_module
    use ns_mesh_metric_module, only: compute_hessian
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol, delta_sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(inout) :: h_p
    real(kind=DOUBLE), dimension(3, 3, mesh%n_vert), intent(in) :: mat_h_p

    character(len=*), intent(in) :: filename


    integer(kind=ENTIER) :: me, num_procs, mpi_ierr
    integer(kind=ENTIER) :: n_interior_elems, n_interior_vert
    integer(kind=ENTIER) :: fn, i, j, s, id_face, k
    real(kind=DOUBLE) :: dmin, dmax

    real(kind=DOUBLE) :: pp
    real(kind=DOUBLE), dimension(5) :: w
    real(kind=DOUBLE), dimension(3, 5) :: nodal_grad
    integer(kind=ENTIER), dimension(:), allocatable :: id_vert_no_ghost
    integer(kind=ENTIER), dimension(:), allocatable :: local_id_vert_no_ghost
    real(kind=DOUBLE), dimension(3, 3) :: hessian
    real(kind=DOUBLE), dimension(mesh%n_elems) :: sol_phi

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    allocate(id_vert_no_ghost(mesh%n_vert))
    id_vert_no_ghost = 0
    n_interior_vert = 0
    do i=1, mesh%n_vert
      if ( .not. mesh%vert(i)%is_ghost ) then
        n_interior_vert = n_interior_vert + 1
        id_vert_no_ghost(i) = n_interior_vert
      end if
    end do

    n_interior_elems = 0
    do i=1, mesh%n_elems
      if ( .not. mesh%elem(i)%is_ghost ) n_interior_elems = n_interior_elems + 1
    end do

    open(newunit=fn, file=trim(adjustl(filename)))

    write(fn, *) "<VTKFile type='UnstructuredGrid' version='1.0' &
      &byte_order='LittleEndian' header_type='UInt64'>"
    write(fn, *) "<UnstructuredGrid>\n<Piece NumberOfPoints='", n_interior_vert, &
      &"' NumberOfCells='", n_interior_elems, "'>"

    dmin = 1e100_DOUBLE
    dmax = -1e100_DOUBLE
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        do j=1, 3
          if( mesh%vert(i)%coord(j) > dmax ) dmax = mesh%vert(i)%coord(j)
          if( mesh%vert(i)%coord(j) < dmin ) dmin = mesh%vert(i)%coord(j)
        end do
      end if
    end do

    write(fn, *) "<Points>"
    write(fn, *) " <DataArray type='Float64' Name='Points' NumberOfComponents='3' &
      & format='ascii' RangeMin='", dmin, "' RangeMax='", dmax, "'>"

    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        write(fn, *) mesh%vert(i)%coord(:)
      end if
    end do

    write(fn, *) "</DataArray>"
    write(fn, *) "</Points>"
    write(fn, *) "<Cells>"
    write(fn, *) "<DataArray type='Int64' Name='connectivity' format='ascii' &
      &RangeMin='0' RangeMax='", mesh%n_vert-1, "'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) id_vert_no_ghost(mesh%elem(i)%vert) - 1
      end if
    end do
    write(fn, *) "</DataArray>"

    s = 0
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        s = s + mesh%elem(i)%n_vert
      end if
    end do

    write(fn, *) "<DataArray type='Int64' Name='offsets' format='ascii' &
      &RangeMin='", 0, "' RangeMax='", s, "'>"

    s = 0
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        s = s + mesh%elem(i)%n_vert
        write(fn, *) s
      end if
    end do

    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='UInt8' Name='types' format='ascii' &
      &RangeMin='42' RangeMax='42'>"

    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) 42 !Type for polyhedra
      end if
    end do

    write(fn, *) "</DataArray>"
    write(fn, *) "<DataArray type='Int64' Name='faces' format='ascii' &
      &RangeMin='0' RangeMax='", mesh%n_vert-1, "'>"

    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) mesh%elem(i)%n_faces
        do j=1, mesh%elem(i)%n_faces
          id_face = mesh%elem(i)%face(j)
          allocate(local_id_vert_no_ghost(mesh%face(id_face)%n_vert))
          do k=1, mesh%face(id_face)%n_vert
            local_id_vert_no_ghost(k) = id_vert_no_ghost(mesh%face(id_face)%vert(k))
          end do
          !write(fn, *) mesh%face(id_face)%n_vert, mesh%face(id_face)%vert(:) - 1
          write(fn, *) mesh%face(id_face)%n_vert, local_id_vert_no_ghost - 1
          deallocate(local_id_vert_no_ghost)
        end do
      end if
    end do

    write(fn, *) "</DataArray>"

    s = 0
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        do j=1, mesh%elem(i)%n_faces
          id_face = mesh%elem(i)%face(j)
          s = s + mesh%face(id_face)%n_vert + 1
        end do
        s = s + 1
      end if
    end do

    write(fn, *) "<DataArray type='Int64' Name='faceoffsets' format='ascii' &
      &RangeMin='", 0, "' RangeMax='", s, "'>"

    s = 0
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        do j=1, mesh%elem(i)%n_faces
          id_face = mesh%elem(i)%face(j)
          s = s + mesh%face(id_face)%n_vert + 1
        end do
        s = s + 1
        write(fn, *) s
      end if
    end do

    write(fn, *) "</DataArray>"
    write(fn, *) "</Cells>"

    !Point data
    write(fn, *) "<PointData>"

    write(fn, *) "<DataArray type='Float64' Name='Nodal_Velocity' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        call compute_nodal_velocity_LS(mesh, i, sol, grad, vp, h_p, second_order)
        if (.not. all(ieee_is_finite(vp(:, i)))) vp(:, i) = 0.0_DOUBLE
        write (fn, *) vp(:, i)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Nodal_Pressure' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        call compute_nodal_pressure_LS(mesh, i, sol, grad, pp, h_p, second_order)
        if (.not. ieee_is_finite(pp)) pp = 0.0_DOUBLE
        write (fn, *) pp
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Nodal_Length' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        write (fn, *) mat_h_p(1,1,i), mat_h_p(2,2,i), mat_h_p(3,3,i)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Nodal_Grad_Density' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        call compute_nodal_grad(mesh, i, sol, nodal_grad)
        write (fn, *) nodal_grad(:, 1)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Nodal_Grad_Velocity_X' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        call compute_nodal_grad(mesh, i, sol, nodal_grad)
        write (fn, *) nodal_grad(:, 2)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Nodal_Grad_Velocity_Y' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        call compute_nodal_grad(mesh, i, sol, nodal_grad)
        write (fn, *) nodal_grad(:, 3)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Nodal_Grad_Velocity_Z' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        call compute_nodal_grad(mesh, i, sol, nodal_grad)
        write (fn, *) nodal_grad(:, 4)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Nodal_Grad_Pressure' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        call compute_nodal_grad(mesh, i, sol, nodal_grad)
        write (fn, *) nodal_grad(:, 5)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Hessian Mach' format='ascii' NumberOfComponents='9'>"
    do i=1, mesh%n_elems
      sol_phi(i) = mach_u(sol(:, i))
    end do
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        call compute_hessian(mesh, i, sol_phi, hessian)
        write (fn, *) hessian
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Hessian Temp' format='ascii' NumberOfComponents='9'>"
    do i=1, mesh%n_elems
      sol_phi(i) = temp_u(sol(:, i))
    end do
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        call compute_hessian(mesh, i, sol_phi, hessian)
        write (fn, *) hessian
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "</PointData>"

    !Cell data
    write(fn, *) "<CellData>"

    write(fn, *) "<DataArray type='Float64' Name='Centroid' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) mesh%elem(i)%coord
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Int32' Name='Tag' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write (fn, *) mesh%elem(i)%tag
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Int32' Name='MPI_Color' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) me
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Density' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) sol(1, i)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Momentum' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) sol(2:4, i)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Energy' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) sol(5, i)/sol(1, i)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Velocity' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) sol(2:4, i)/sol(1, i)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Pressure' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        w = conserv_to_primit(sol(:, i))
        write(fn, *) w(5)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Internal_energy' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write (fn, *) sol(5, i)/sol(1, i) - 0.5_DOUBLE*norm2(sol(2:4, i)/sol(1, i))**2
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='H' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        w = conserv_to_primit(sol(:, i))
        write (fn, *) sol(5, i)/sol(1, i) + w(5)/sol(1, i)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Mach' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        w = conserv_to_primit(sol(:, i))
        write (fn, *) norm2(sol(2:4, i)/sol(1, i))/sound_speed_w(w)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Temperature' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write (fn, *) temp_u(sol(:, i))
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='mu' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        if( use_sutherland ) then
          write (fn, *) mu_sutherland(sol(:, i))
        else
          write (fn, *) mu_p
        end if
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Ducros' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write (fn, *) omega_ducros(mesh, sol, i)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='r' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write (fn, *) norm2(mesh%elem(i)%coord)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='rxy' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write (fn, *) norm2(mesh%elem(i)%coord(1:2))
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Theta3' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write (fn, *) atan2(-mesh%elem(i)%coord(2), -mesh%elem(i)%coord(1))
      end if
    end do
    write(fn, *) "</DataArray>"

    if( second_order ) then
      write(fn, *) "<DataArray type='Float64' Name='Density_Grad' format='ascii' NumberOfComponents='3'>"
      do i=1, mesh%n_elems
        if( .not. mesh%elem(i)%is_ghost ) then
          write (fn, *) grad(:, 1, i)
        end if
      end do
      write(fn, *) "</DataArray>"

      write(fn, *) "<DataArray type='Float64' Name='Velocity_X_Grad' format='ascii' NumberOfComponents='3'>"
      do i=1, mesh%n_elems
        if( .not. mesh%elem(i)%is_ghost ) then
          write (fn, *) grad(:, 2, i)
        end if
      end do
      write(fn, *) "</DataArray>"

      write(fn, *) "<DataArray type='Float64' Name='Velocity_Y_Grad' format='ascii' NumberOfComponents='3'>"
      do i=1, mesh%n_elems
        if( .not. mesh%elem(i)%is_ghost ) then
          write (fn, *) grad(:, 3, i)
        end if
      end do
      write(fn, *) "</DataArray>"

      write(fn, *) "<DataArray type='Float64' Name='Velocity_Z_Grad' format='ascii' NumberOfComponents='3'>"
      do i=1, mesh%n_elems
        if( .not. mesh%elem(i)%is_ghost ) then
          write (fn, *) grad(:, 4, i)
        end if
      end do
      write(fn, *) "</DataArray>"

      write(fn, *) "<DataArray type='Float64' Name='Pressure_Grad' format='ascii' NumberOfComponents='3'>"
      do i=1, mesh%n_elems
        if( .not. mesh%elem(i)%is_ghost ) then
          write (fn, *) grad(:, 5, i)
        end if
      end do
      write(fn, *) "</DataArray>"
    end if

    write(fn, *) "<DataArray type='Float64' Name='Residual' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write (fn, *) norm2(delta_sol(:, i))
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "</CellData>"

    write(fn, *) "</Piece>"
    write(fn, *) "</UnstructuredGrid>"
    write(fn, *) "</VTKFile>"

    close(fn)
  end subroutine write_sol_vtu

  subroutine write_sol_meta_pvtu(filename, iaff_char)
    use mpi
    use mpi_module
    use ns_global_data_module, only: second_order, scheme
    implicit none

    character(len=*), intent(in) :: filename, iaff_char

    integer(kind=ENTIER) :: fn, i, num_procs, me, mpi_ierr
    character(len=255) :: i_char

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    if( me == 0 ) then
      open(newunit=fn, file=trim(adjustl(filename)))

      write(fn, *) "<VTKFile type='PUnstructuredGrid' version='0.1' byte_order='LittleEndian' header_type='UInt64'>"
      write(fn, *) "<PUnstructuredGrid GhostLevel='0'>"
      write(fn, *) "<PPoints>"
      write(fn, *) "<PDataArray type='Float64' Name='Points' NumberOfComponents='3'/>"
      write(fn, *) "</PPoints>"
      write(fn, *) "<PCells>"
      write(fn, *) "<PDataArray type='Int64' Name='connectivity' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Int64' Name='offsets'      NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='UInt8' Name='types'        NumberOfComponents='1'/>"
      write(fn, *) "</PCells>"
      write(fn, *) "<PPointData>"
      write(fn, *) "<PDataArray type='Float64' Name='Nodal_Velocity' format='ascii' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Nodal_Pressure' format='ascii' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Nodal_Length' format='ascii' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Nodal_Grad_Density' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Nodal_Grad_Velocity_X' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Nodal_Grad_Velocity_Y' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Nodal_Grad_Velocity_Z' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Nodal_Grad_Pressure' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Hessian Mach' NumberOfComponents='9'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Hessian Temp' NumberOfComponents='9'/>"
      write(fn, *) "</PPointData>"
      write(fn, *) "<PCellData>"
      write(fn, *) "<PDataArray type='Float64' Name='Centroid' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Int32' Name='Tag' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Int32' Name='MPI_Color' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Density' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Momentum' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Energy' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Velocity' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Pressure' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Internal_energy' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='H' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Mach' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Temperature' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='mu' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Ducros' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='r' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='rxy' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Theta3' NumberOfComponents='1'/>"
      if( second_order ) then
        write(fn, *) "<PDataArray type='Float64' Name='Density_Grad' NumberOfComponents='3'/>"
        write(fn, *) "<PDataArray type='Float64' Name='Velocity_X_Grad' NumberOfComponents='3'/>"
        write(fn, *) "<PDataArray type='Float64' Name='Velocity_Y_Grad' NumberOfComponents='3'/>"
        write(fn, *) "<PDataArray type='Float64' Name='Velocity_Z_Grad' NumberOfComponents='3'/>"
        write(fn, *) "<PDataArray type='Float64' Name='Pressure_Grad' NumberOfComponents='3'/>"
      end if
      write(fn, *) "<PDataArray type='Float64' Name='Residual' NumberOfComponents='1'/>"
      write(fn, *) "</PCellData>"

      do i=0, num_procs-1
        write(i_char, *) i
        write(fn, *) "<Piece Source='"//trim(adjustl(i_char))//"_output_"//&
          trim(adjustl(iaff_char))//".vtu'/>"
      end do

      write(fn, *) "</PUnstructuredGrid>"
      write(fn, *) "</VTKFile>"

      close(fn)
    end if

    call mpi_barrier(mpi_comm_world, mpi_ierr)
  end subroutine write_sol_meta_pvtu

  subroutine write_sol_meta_pvtu_coeffs(filename, iaff_char)
    use mpi
    use mpi_module
    use ns_global_data_module, only: second_order
    implicit none

    character(len=*), intent(in) :: filename, iaff_char

    integer(kind=ENTIER) :: fn, i, num_procs, me, mpi_ierr
    character(len=255) :: i_char

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    if( me == 0 ) then
      open(newunit=fn, file=trim(adjustl(filename)))

      write(fn, *) "<VTKFile type='PUnstructuredGrid' version='0.1' byte_order='LittleEndian' header_type='UInt64'>"
      write(fn, *) "<PUnstructuredGrid GhostLevel='0'>"
      write(fn, *) "<PPoints>"
      write(fn, *) "<PDataArray type='Float64' Name='Points' NumberOfComponents='3'/>"
      write(fn, *) "</PPoints>"
      write(fn, *) "<PCells>"
      write(fn, *) "<PDataArray type='Int64' Name='connectivity' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Int64' Name='offsets'      NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='UInt8' Name='types'        NumberOfComponents='1'/>"
      write(fn, *) "</PCells>"
      write(fn, *) "<PCellData>"
      write(fn, *) "<PDataArray type='Float64' Name='Centroid' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Density' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Pressure' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Temperature' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Theta3' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Cp' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='q' NumberOfComponents='1'/>"
      write(fn, *) "</PCellData>"

      do i=0, num_procs-1
        write(i_char, *) i
        write(fn, *) "<Piece Source='"//trim(adjustl(i_char))//"_coeffs_"//&
          trim(adjustl(iaff_char))//".vtu'/>"
      end do

      write(fn, *) "</PUnstructuredGrid>"
      write(fn, *) "</VTKFile>"

      close(fn)
    end if

    call mpi_barrier(mpi_comm_world, mpi_ierr)
  end subroutine write_sol_meta_pvtu_coeffs

  subroutine write_cell_size_dat(filename, mesh, sol, grad)
    use mpi
    use mpi_module
    use linear_solver_module, only: spectral_decomposition
    use ns_mesh_metric_module, only: compute_cell_metric
    implicit none

    character(len=*), intent(in) :: filename

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, 5, mesh%n_elems), intent(inout) :: grad
    integer(kind=ENTIER) :: fn, i, num_procs, me, mpi_ierr, p, n_cells

    real(kind=DOUBLE), dimension(mesh%n_vert) :: cell_size
    real(kind=DOUBLE), dimension(3, 3, mesh%n_vert) :: cell_metric

    real(kind=DOUBLE), dimension(3, 3) :: Q
    real(kind=DOUBLE), dimension(3) :: Lambda


    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    call compute_cell_metric(mesh, sol, cell_metric)
    do i=1, mesh%n_vert
      call spectral_decomposition(cell_metric(:, :, i),Q,Lambda)
      !cell_size(i) = max(2e-2, sqrt(1.0_DOUBLE/(1e-12+2*maxval(Lambda))))!Mach
      cell_size(i) = max(2e-2, 4e2*sqrt(1.0_DOUBLE/(1e-16+maxval(Lambda))))!Temp
      if( cell_size(i) > 1.0_DOUBLE ) cell_size(i) = -1
    end do

    do p=0, num_procs-1

      n_cells = 0
      do i=1, mesh%n_vert
        if( cell_size(i) > 0.0_DOUBLE) then
          n_cells = n_cells + 1
        end if
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE, n_cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)

      if( me == p ) then
        if( p==0 ) then
          open(newunit=fn, file=trim(adjustl(filename)))
          write(fn, *) n_cells
        else
          open(newunit=fn, file=trim(adjustl(filename)), position='append', status='old')
        end if

        n_cells = 0
        do i=1, mesh%n_vert
          if( cell_size(i) > 0.0_DOUBLE) then
            n_cells = n_cells + 1
          end if
        end do

        do i=1, mesh%n_vert
          if( cell_size(i) > 0.0_DOUBLE) then
            write(fn, *) mesh%vert(i)%coord, cell_size(i)
          end if
        end do

        if( p==num_procs-1 ) write(fn, *) 0
        close(fn)
      end if
      call mpi_barrier(mpi_comm_world, mpi_ierr)
    end do
  end subroutine write_cell_size_dat

  subroutine write_cell_size_mmg(filename, mesh, sol, grad, all_nodal_grad)
    use mpi
    use mpi_module
    use ns_global_data_module, only: second_order
    use ns_euler_recon_module, only: compute_cell_grad_from_nodal_grad
    use sort_module, only: qsort_mat
    use linear_solver_module, only: spectral_decomposition
    use ns_mesh_metric_module, only: compute_cell_metric
    implicit none

    character(len=*), intent(in) :: filename

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol 
    real(kind=DOUBLE), dimension(3, 5, mesh%n_elems), intent(inout) :: grad
    real(kind=DOUBLE), dimension(3, 5, mesh%n_vert), intent(inout) :: all_nodal_grad
    integer(kind=ENTIER) :: fn, i, num_procs, me, mpi_ierr, n_uniq
    integer(kind=ENTIER) :: j, last_glob, iloc

    integer(kind=ENTIER), dimension(:), allocatable :: n_vert_vect
    real(kind=DOUBLE), dimension(mesh%n_vert) :: cell_size
    real(kind=DOUBLE), dimension(3, 3, mesh%n_vert) :: cell_metric
    integer(kind=ENTIER), dimension(mesh%n_vert) :: id_glob

    integer, dimension(:, :), allocatable :: status

    real(kind=DOUBLE), dimension(:), allocatable :: all_cell_size
    real(kind=DOUBLE), dimension(:,:,:), allocatable :: all_cell_metric
    integer(kind=ENTIER), dimension(:), allocatable :: all_id_glob

    type neigh_type
      real(kind=DOUBLE), dimension(:), allocatable :: cell_size
      real(kind=DOUBLE), dimension(:,:,:), allocatable :: cell_metric
      integer(kind=ENTIER), dimension(:), allocatable :: id_glob
    end type neigh_type

    type(neigh_type), dimension(:), allocatable :: neigh

    real(kind=DOUBLE), dimension(3, 3) :: Q
    real(kind=DOUBLE), dimension(3) :: Lambda

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    call compute_cell_grad_from_nodal_grad(mesh, sol, grad, all_nodal_grad, 2)
    call compute_cell_metric(mesh, sol, cell_metric)

    do i=1, mesh%n_vert
      call spectral_decomposition(cell_metric(:, :, i),Q,Lambda)
      cell_size(i) = maxval(Lambda)
    end do

    do i=1, mesh%n_vert
      id_glob(i) = mesh%vert(i)%id_glob
    end do

    if( num_procs > 1 ) then
      allocate(n_vert_vect(num_procs))
      n_vert_vect = 0
      call MPI_Gather(mesh%n_vert, 1, MPI_INT, &
        n_vert_vect, 1, MPI_INT, &
        0, MPI_COMM_WORLD, mpi_ierr)

      if( me == 0 ) then
        allocate(status(MPI_STATUS_SIZE, num_procs))
        allocate(neigh(num_procs))
        do i=1, num_procs
          allocate(neigh(i)%cell_size(n_vert_vect(i)))
          allocate(neigh(i)%cell_metric(3, 3, n_vert_vect(i)))
          allocate(neigh(i)%id_glob(n_vert_vect(i)))
        end do
        do i=1, mesh%n_vert
          neigh(1)%cell_size(i) = cell_size(i)
          neigh(1)%cell_metric(:,:,i) = cell_metric(:,:,i)
          neigh(1)%id_glob(i) = mesh%vert(i)%id_glob
        end do
        do i=2, num_procs
          call MPI_RECV(neigh(i)%cell_size, n_vert_vect(i), MPI_DOUBLE, &
            i-1, i-1, MPI_COMM_WORLD, status(:, i), mpi_ierr)
          call MPI_RECV(neigh(i)%cell_metric(1, 1, 1), 3*3*n_vert_vect(i), MPI_DOUBLE, &
            i-1, i-1, MPI_COMM_WORLD, status(:, i), mpi_ierr)
          call MPI_RECV(neigh(i)%id_glob, n_vert_vect(i), MPI_INT, &
            i-1, i-1, MPI_COMM_WORLD, status(:, i), mpi_ierr)
        end do
      else
        call MPI_SEND(cell_size, mesh%n_vert, MPI_DOUBLE, &
          0, me, MPI_COMM_WORLD, mpi_ierr)
        call MPI_SEND(cell_metric, 3*3*mesh%n_vert, MPI_DOUBLE, &
          0, me, MPI_COMM_WORLD, mpi_ierr)
        call MPI_SEND(id_glob, mesh%n_vert, MPI_INT, &
          0, me, MPI_COMM_WORLD, mpi_ierr)
      end if

      call mpi_barrier(mpi_comm_world, mpi_ierr)

      if( me == 0 ) then
        allocate( all_cell_size(sum(n_vert_vect)) )
        allocate( all_cell_metric(3,3,sum(n_vert_vect)) )
        allocate( all_id_glob(sum(n_vert_vect)) )
        iloc = 1
        do i=1, num_procs
          do j=1, n_vert_vect(i)
            all_cell_size(iloc) = neigh(i)%cell_size(j)
            all_cell_metric(:,:,iloc) = neigh(i)%cell_metric(:,:,j)
            all_id_glob(iloc) = neigh(i)%id_glob(j)
            iloc = iloc + 1
          end do
        end do

        !call hpsort_real(sum(n_vert_vect), all_id_glob, all_cell_size)
        call qsort_mat(sum(n_vert_vect), all_id_glob, all_cell_metric)

        n_uniq = 1
        last_glob = all_id_glob(1)
        do i=2, sum(n_vert_vect)
          if(all_id_glob(i) /= last_glob) then
            last_glob = all_id_glob(i)
            n_uniq = n_uniq + 1
          end if
        end do

        open(newunit=fn, file=trim(adjustl(filename)))
        write(fn, *) "MeshVersionFormatted 2"
        write(fn, *) ""
        write(fn, *) "Dimension 3"
        write(fn, *) ""
        write(fn, *) "SolAtVertices"
        write(fn, *) ""
        write(fn, *) n_uniq
        write(fn, *) "1 6"
        write(fn, *) ""

        n_uniq = 1
        last_glob = all_id_glob(1)
        write(fn, *) all_cell_metric(1, 1, 1), &
          all_cell_metric(1, 2, 1), &
          all_cell_metric(1, 3, 1), &
          all_cell_metric(2, 2, 1), &
          all_cell_metric(2, 3, 1), &
          all_cell_metric(3, 3, 1)
        do i=2, sum(n_vert_vect)
          if(all_id_glob(i) /= last_glob) then
            last_glob = all_id_glob(i)
            n_uniq = n_uniq + 1
            write(fn, *) all_cell_metric(1, 1, i), &
              all_cell_metric(1, 2, i), &
              all_cell_metric(1, 3, i), &
              all_cell_metric(2, 2, i), &
              all_cell_metric(2, 3, i), &
              all_cell_metric(3, 3, i)
          end if
        end do

        close(fn)
      end if

      call mpi_barrier(mpi_comm_world, mpi_ierr)
    else 

      !call hpsort_real(mesh%n_vert, id_glob, cell_size)
      call qsort_mat(mesh%n_vert, id_glob, cell_metric)

      open(newunit=fn, file=trim(adjustl(filename)))
      write(fn, *) "MeshVersionFormatted 2"
      write(fn, *) ""
      write(fn, *) "Dimension 3"
      write(fn, *) ""
      write(fn, *) "SolAtVertices"
      write(fn, *) ""
      write(fn, *) mesh%n_vert
      write(fn, *) "1 6"
      write(fn, *) ""

      do i=1, mesh%n_vert
        write(fn, *) cell_metric(1, 1, i), &
          cell_metric(1, 2, i), &
          cell_metric(1, 3, i), &
          cell_metric(2, 2, i), &
          cell_metric(2, 3, i), &
          cell_metric(3, 3, i)
      end do
      close(fn)
    end if

    call mpi_barrier(mpi_comm_world, mpi_ierr)
  end subroutine write_cell_size_mmg

  subroutine init_residual_file(me, fn_residual)
    use ns_global_data_module, only: write_residual, init_restart, &
      restart_iter, restart_time, restart_cpu_time
    implicit none
    integer(kind=ENTIER), intent(in) :: me
    integer(kind=ENTIER), intent(out) :: fn_residual
    integer(kind=ENTIER) :: fn_residual_error
    character(len=255) :: dummy_read
    logical :: res_exists
    if (.not. (write_residual .and. me == 0)) return
    inquire(file="residual.dat", exist=res_exists)
    open(newunit=fn_residual, file="residual.dat")
    if (init_restart .and. res_exists) then
      do
        read(fn_residual, *, iostat=fn_residual_error) &
          restart_iter, restart_time, dummy_read(1:3), restart_cpu_time
        if (fn_residual_error == iostat_end) exit
      end do
      close(fn_residual)
      open(newunit=fn_residual, file="residual.dat", position="append")
    end if
  end subroutine init_residual_file

  subroutine write_iter_info(iter, res, res_vect, me, fn_residual, start_count, count_rate, &
      mesh, sol, delta_sol, cell_grad, nodal_grad, vp, h_p, mat_h_p)
    use ns_global_data_module, only: n_iter_print, local_time_step, t, &
      write_residual, n_iter_residual, restart_iter, &
      restart_cpu_time, n_iter_write_sol
    implicit none
    integer(kind=ENTIER), intent(in) :: iter, me, fn_residual, start_count, count_rate
    real(kind=DOUBLE), intent(in) :: res
    real(kind=DOUBLE), dimension(5), intent(in) :: res_vect
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol, delta_sol
    real(kind=DOUBLE), dimension(3, 5, mesh%n_elems), intent(inout) :: cell_grad
    real(kind=DOUBLE), dimension(3, 5, mesh%n_vert), intent(inout) :: nodal_grad
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(inout) :: h_p
    real(kind=DOUBLE), dimension(3, 3, mesh%n_vert), intent(in) :: mat_h_p
    integer(kind=ENTIER) :: end_count
    if (me == 0 .and. mod(iter, n_iter_print) == 0) then
      if (local_time_step) then
        print *, iter, res
      else
        print *, iter, t
      end if
    end if
    if (me == 0 .and. write_residual .and. mod(iter, n_iter_residual) == 0) then
      call system_clock(end_count)
      write(fn_residual, *) iter + restart_iter, &
        restart_cpu_time + real(end_count-start_count, DOUBLE)/real(count_rate, DOUBLE), res, res_vect
    end if
    if (n_iter_write_sol > 0 .and. mod(iter+1, n_iter_write_sol) == 0) &
      call write_sol(mesh, sol, delta_sol, cell_grad, vp, h_p, mat_h_p, me, (iter+1)/n_iter_write_sol)
  end subroutine write_iter_info

  subroutine restart_from_vtu_file(mesh, sol, restart_file, me, num_procs)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: sol
    integer(kind=ENTIER), intent(in) :: me, num_procs
    character(len=*), intent(in) :: restart_file

    integer(kind=ENTIER) :: i, vtkin
    character(len=255) :: text, me_str, mpi_restart_file

    if (num_procs > 1) then
      write (me_str, *) me
      mpi_restart_file = trim(adjustl(me_str))//"_"//trim(adjustl(restart_file))
      open (newunit=vtkin, file=trim(mpi_restart_file), status="old")
    else
      open (newunit=vtkin, file=trim(restart_file), status="old")
    end if

    read (vtkin, '(a)') text
    do while (" <DataArray type='Float64' Name='Density' format='ascii' NumberOfComponents='1'>" /= text(:80))
      read (vtkin, '(a)') text
    end do
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        read (vtkin, *) sol(1, i)
      end if
    end do

    do while (" <DataArray type='Float64' Name='Momentum' format='ascii' NumberOfComponents='3'>" /= text(:81))
      read (vtkin, '(a)') text
    end do
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        read (vtkin, *) sol(2, i), sol(3, i), sol(4, i)
      end if
    end do

    do while (" <DataArray type='Float64' Name='Energy' format='ascii' NumberOfComponents='1'>" /= text(:79))
      read (vtkin, '(a)') text
    end do
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        read (vtkin, *) sol(5, i)
        sol(5, i) = sol(5, i) * sol(1, i)
      end if
    end do

    close (vtkin)
  end subroutine restart_from_vtu_file
end module ns_io_module