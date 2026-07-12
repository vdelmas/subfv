module lagrange_io_module
  use mpi
  use precision_module
  use mpi_module
  use mesh_module
  use mesh_reading_module
  use mesh_geometry_module
  use mesh_connectivity_module
  use lagrange_module
  implicit none
contains
  subroutine write_sol_lag(mesh, filename, sol, vp, pp)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: vp
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: pp
    character(len=*), intent(in) :: filename

    call write_sol_meta_pvtu_lag(filename)
    call write_sol_vtu_lag(mesh, filename, sol, vp, pp)

  end subroutine write_sol_lag

  subroutine write_sol_vtu_lag(mesh, basename, sol, vp, pp)
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: vp
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: pp
    character(len=*), intent(in) :: basename

    integer(kind=ENTIER) :: me, num_procs, mpi_ierr
    integer(kind=ENTIER) :: n_interior_elems, n_interior_vert
    integer(kind=ENTIER) :: fn, i, j, s, id_face, size_tot, k
    real(kind=DOUBLE) :: dmin, dmax
    character(len=255) :: filename, me_char

    real(kind=DOUBLE), dimension(5) :: w
    integer(kind=ENTIER), dimension(:), allocatable :: id_vert_no_ghost
    integer(kind=ENTIER), dimension(:), allocatable :: local_id_vert_no_ghost

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

    write(me_char, *) me
    write(filename, *) trim(adjustl(me_char))//"_"//trim(adjustl(basename))//".vtu"
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

    write(fn, *) "<DataArray type='Float64' Name='Vp' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        write (fn, *) vp(:, i)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Pp' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        write (fn, *) pp(i)
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

    write(fn, *) "<DataArray type='Float64' Name='Rxy' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) norm2(mesh%elem(i)%coord(:2))
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='R' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) norm2(mesh%elem(i)%coord)
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

    write(fn, *) "<DataArray type='Int32' Name='MPI_Color' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write (fn, *) me
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Density' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) 1.0_DOUBLE/sol(1, i)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Velocity' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) sol(2:4, i)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Pressure' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) pressure(sol(:, i))
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Internal_energy' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) pressure(sol(:, i))*sol(1, i)/(gamma-1.0_DOUBLE)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "</CellData>"

    write(fn, *) "</Piece>"
    write(fn, *) "</UnstructuredGrid>"
    write(fn, *) "</VTKFile>"

    close(fn)
  end subroutine write_sol_vtu_lag

  subroutine write_sol_meta_pvtu_lag(basename)
    use mpi
    use mpi_module
    use lagrange_global_data_module, only: second_order
    implicit none

    character(len=*), intent(in) :: basename

    integer(kind=ENTIER) :: fn, i, num_procs, me, mpi_ierr
    character(len=255) :: i_char, filename

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    if( me == 0 ) then
      write(filename, *) trim(adjustl(basename))//".pvtu"
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
      write(fn, *) "<DataArray type='Float64' Name='Vp' format='ascii' NumberOfComponents='3'/>"
      write(fn, *) "<DataArray type='Float64' Name='Pp' format='ascii' NumberOfComponents='1'/>"
      write(fn, *) "</PPointData>"
      write(fn, *) "<PCellData>"
      write(fn, *) "<PDataArray type='Float64' Name='Centroid' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Rxy' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='R' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Int32' Name='Tag' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Int32' Name='MPI_Color' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Density' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Velocity' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Pressure' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Internal_energy' NumberOfComponents='1'/>"
      write(fn, *) "</PCellData>"

      do i=0, num_procs-1
        write(i_char, *) i
        write(fn, *) "<Piece Source='"//trim(adjustl(i_char))//"_"//trim(adjustl(basename))//".vtu'/>"
      end do

      write(fn, *) "</PUnstructuredGrid>"
      write(fn, *) "</VTKFile>"

      close(fn)
    end if

    call mpi_barrier(mpi_comm_world, mpi_ierr)
  end subroutine write_sol_meta_pvtu_lag
end module lagrange_io_module
