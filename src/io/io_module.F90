module io_module
  implicit none
contains
  subroutine write_file_vtu_start_cell_data(mesh, filename, fn_vtu, fn_pvtu)
    use mpi
    use mpi_module
    use mesh_module
    implicit none

    type(mesh_type) :: mesh
    character(len=*), intent(in) :: filename
    integer(kind=ENTIER), intent(inout) :: fn_vtu, fn_pvtu

    integer(kind=ENTIER) :: num_procs, me, mpi_ierr
    character(len=255) :: filename_vtu, filename_pvtu, me_str

    integer(kind=ENTIER) :: i

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    if( me == 0 ) then
      !PVTU header
      write(fn_pvtu, *) "<PCellData>"
    end if

    !VTU
    write(fn_vtu, *) "<CellData>"
  end subroutine write_file_vtu_start_cell_data

  subroutine write_file_vtu_end_cell_data(mesh, filename, fn_vtu, fn_pvtu)
    use mpi
    use mpi_module
    use mesh_module
    implicit none

    type(mesh_type) :: mesh
    character(len=*), intent(in) :: filename
    integer(kind=ENTIER), intent(inout) :: fn_vtu, fn_pvtu

    integer(kind=ENTIER) :: num_procs, me, mpi_ierr
    character(len=255) :: filename_vtu, filename_pvtu, me_str

    integer(kind=ENTIER) :: i

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    if( me == 0 ) then
      !PVTU header
      write(fn_pvtu, *) "</PCellData>"
    end if

    !VTU
    write(fn_vtu, *) "</CellData>"
  end subroutine write_file_vtu_end_cell_data

  subroutine write_file_vtu_cell_scalar(mesh, filename, fn_vtu, fn_pvtu, data, name)
    use mpi
    use mpi_module
    use mesh_module
    implicit none

    type(mesh_type) :: mesh
    character(len=*), intent(in) :: filename, name
    integer(kind=ENTIER), intent(inout) :: fn_vtu, fn_pvtu
    real(kind=DOUBLE), dimension(mesh%n_elems) :: data

    integer(kind=ENTIER) :: num_procs, me, mpi_ierr
    character(len=255) :: filename_vtu, filename_pvtu, me_str

    integer(kind=ENTIER) :: i

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    if( me == 0 ) then
      !PVTU header
      write(fn_pvtu, *) "<DataArray type='Float64' Name='"&
        //trim(adjustl(name))//"' format='ascii' NumberOfComponents='1'/>"
    end if

    !VTU
    write(fn_vtu, *) "<DataArray type='Float64' Name='"&
      //trim(adjustl(name))//"' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write (fn_vtu, *) data(i)
      end if
    end do
    write(fn_vtu, *) "</DataArray>"
  end subroutine write_file_vtu_cell_scalar

  subroutine write_file_vtu_cell_vector(mesh, filename, fn_vtu, fn_pvtu, data, name)
    use mpi
    use mpi_module
    use mesh_module
    implicit none

    type(mesh_type) :: mesh
    character(len=*), intent(in) :: filename, name
    integer(kind=ENTIER), intent(inout) :: fn_vtu, fn_pvtu
    real(kind=DOUBLE), dimension(3, mesh%n_elems) :: data

    integer(kind=ENTIER) :: num_procs, me, mpi_ierr
    character(len=255) :: filename_vtu, filename_pvtu, me_str

    integer(kind=ENTIER) :: i

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    if( me == 0 ) then
      !PVTU header
      write(fn_pvtu, *) "<DataArray type='Float64' Name='"&
        //trim(adjustl(name))//"' format='ascii' NumberOfComponents='3'/>"
    end if

    !VTU
    write(fn_vtu, *) "<DataArray type='Float64' Name='"&
      //trim(adjustl(name))//"' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write (fn_vtu, *) data(:, i)
      end if
    end do
    write(fn_vtu, *) "</DataArray>"
  end subroutine write_file_vtu_cell_vector

  subroutine write_file_vtu_start_vert_data(mesh, filename, fn_vtu, fn_pvtu)
    use mpi
    use mpi_module
    use mesh_module
    implicit none

    type(mesh_type) :: mesh
    character(len=*), intent(in) :: filename
    integer(kind=ENTIER), intent(inout) :: fn_vtu, fn_pvtu

    integer(kind=ENTIER) :: num_procs, me, mpi_ierr
    character(len=255) :: filename_vtu, filename_pvtu, me_str

    integer(kind=ENTIER) :: i

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    if( me == 0 ) then
      !PVTU header
      write(fn_pvtu, *) "<PPointData>"
    end if

    !VTU
    write(fn_vtu, *) "<PointData>"
  end subroutine write_file_vtu_start_vert_data

  subroutine write_file_vtu_end_vert_data(mesh, filename, fn_vtu, fn_pvtu)
    use mpi
    use mpi_module
    use mesh_module
    implicit none

    type(mesh_type) :: mesh
    character(len=*), intent(in) :: filename
    integer(kind=ENTIER), intent(inout) :: fn_vtu, fn_pvtu

    integer(kind=ENTIER) :: num_procs, me, mpi_ierr
    character(len=255) :: filename_vtu, filename_pvtu, me_str

    integer(kind=ENTIER) :: i

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    if( me == 0 ) then
      !PVTU header
      write(fn_pvtu, *) "</PPointData>"
    end if

    !VTU
    write(fn_vtu, *) "</PointData>"
  end subroutine write_file_vtu_end_vert_data

  subroutine write_file_vtu_vert_scalar(mesh, filename, fn_vtu, fn_pvtu, data, name)
    use mpi
    use mpi_module
    use mesh_module
    implicit none

    type(mesh_type) :: mesh
    character(len=*), intent(in) :: filename, name
    integer(kind=ENTIER), intent(inout) :: fn_vtu, fn_pvtu
    real(kind=DOUBLE), dimension(mesh%n_vert) :: data

    integer(kind=ENTIER) :: num_procs, me, mpi_ierr
    character(len=255) :: filename_vtu, filename_pvtu, me_str

    integer(kind=ENTIER) :: i

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    if( me == 0 ) then
      !PVTU header
      write(fn_pvtu, *) "<DataArray type='Float64' Name='"&
        //trim(adjustl(name))//"' format='ascii' NumberOfComponents='1'/>"
    end if

    !VTU
    write(fn_vtu, *) "<DataArray type='Float64' Name='"&
      //trim(adjustl(name))//"' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        write (fn_vtu, *) data(i)
      end if
    end do
    write(fn_vtu, *) "</DataArray>"
  end subroutine write_file_vtu_vert_scalar

  subroutine write_file_vtu_vert_vector(mesh, filename, fn_vtu, fn_pvtu, data, name)
    use mpi
    use mpi_module
    use mesh_module
    implicit none

    type(mesh_type) :: mesh
    character(len=*), intent(in) :: filename, name
    integer(kind=ENTIER), intent(inout) :: fn_vtu, fn_pvtu
    real(kind=DOUBLE), dimension(3, mesh%n_vert) :: data

    integer(kind=ENTIER) :: num_procs, me, mpi_ierr
    character(len=255) :: filename_vtu, filename_pvtu, me_str

    integer(kind=ENTIER) :: i

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    if( me == 0 ) then
      !PVTU header
      write(fn_pvtu, *) "<DataArray type='Float64' Name='"&
        //trim(adjustl(name))//"' format='ascii' NumberOfComponents='3'/>"
    end if

    !VTU
    write(fn_vtu, *) "<DataArray type='Float64' Name='"&
      //trim(adjustl(name))//"' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        write (fn_vtu, *) data(:, i)
      end if
    end do
    write(fn_vtu, *) "</DataArray>"
  end subroutine write_file_vtu_vert_vector

  subroutine open_file_vtu(mesh, filename, fn_vtu, fn_pvtu)
    use mpi
    use mpi_module
    use mesh_module
    implicit none

    type(mesh_type) :: mesh
    character(len=*), intent(in) :: filename
    integer(kind=ENTIER), intent(inout) :: fn_vtu, fn_pvtu

    integer(kind=ENTIER) :: num_procs, me, mpi_ierr
    character(len=255) :: filename_vtu, filename_pvtu, me_str

    integer(kind=ENTIER) :: i, j, k, s, id_face
    integer(kind=ENTIER) :: n_interior_elems, n_interior_vert
    real(kind=DOUBLE) :: dmin, dmax

    integer(kind=ENTIER), dimension(:), allocatable :: id_vert_no_ghost
    integer(kind=ENTIER), dimension(:), allocatable :: local_id_vert_no_ghost

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    if( me == 0 ) then
      !PVTU header
      filename_pvtu = trim(adjustl(filename))//".pvtu"
      open(newunit=fn_pvtu, file=trim(adjustl(filename_pvtu)))

      write(fn_pvtu, *) "<VTKFile type='PUnstructuredGrid' version='0.1' byte_order='LittleEndian' header_type='UInt64'>"
      write(fn_pvtu, *) "<PUnstructuredGrid GhostLevel='0'>"
      write(fn_pvtu, *) "<PPoints>"
      write(fn_pvtu, *) "<PDataArray type='Float64' Name='Points' NumberOfComponents='3'/>"
      write(fn_pvtu, *) "</PPoints>"
      write(fn_pvtu, *) "<PCells>"
      write(fn_pvtu, *) "<PDataArray type='Int64' Name='connectivity' NumberOfComponents='1'/>"
      write(fn_pvtu, *) "<PDataArray type='Int64' Name='offsets'      NumberOfComponents='1'/>"
      write(fn_pvtu, *) "<PDataArray type='UInt8' Name='types'        NumberOfComponents='1'/>"
      write(fn_pvtu, *) "</PCells>"
    end if

    !VTU
    write(me_str, *) me
    filename_vtu = trim(adjustl(me_str))//"_"//trim(adjustl(filename))//".vtu"
    open(newunit=fn_vtu, file=trim(adjustl(filename_vtu)))

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

    write(fn_vtu, *) "<VTKFile type='UnstructuredGrid' version='1.0' &
      &byte_order='LittleEndian' header_type='UInt64'>"
    write(fn_vtu, *) "<UnstructuredGrid>\n<Piece NumberOfPoints='", n_interior_vert, &
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

    write(fn_vtu, *) "<Points>"
    write(fn_vtu, *) " <DataArray type='Float64' Name='Points' NumberOfComponents='3' &
      & format='ascii' RangeMin='", dmin, "' RangeMax='", dmax, "'>"

    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        write(fn_vtu, *) mesh%vert(i)%coord(:)
      end if
    end do

    write(fn_vtu, *) "</DataArray>"
    write(fn_vtu, *) "</Points>"
    write(fn_vtu, *) "<Cells>"
    write(fn_vtu, *) "<DataArray type='Int64' Name='connectivity' format='ascii' &
      &RangeMin='0' RangeMax='", mesh%n_vert-1, "'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn_vtu, *) id_vert_no_ghost(mesh%elem(i)%vert) - 1
      end if
    end do
    write(fn_vtu, *) "</DataArray>"

    s = 0
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        s = s + mesh%elem(i)%n_vert
      end if
    end do

    write(fn_vtu, *) "<DataArray type='Int64' Name='offsets' format='ascii' &
      &RangeMin='", 0, "' RangeMax='", s, "'>"

    s = 0
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        s = s + mesh%elem(i)%n_vert
        write(fn_vtu, *) s
      end if
    end do

    write(fn_vtu, *) "</DataArray>"

    write(fn_vtu, *) "<DataArray type='UInt8' Name='types' format='ascii' &
      &RangeMin='42' RangeMax='42'>"

    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn_vtu, *) 42 !Type for polyhedra
      end if
    end do

    write(fn_vtu, *) "</DataArray>"
    write(fn_vtu, *) "<DataArray type='Int64' Name='faces' format='ascii' &
      &RangeMin='0' RangeMax='", mesh%n_vert-1, "'>"

    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn_vtu, *) mesh%elem(i)%n_faces
        do j=1, mesh%elem(i)%n_faces
          id_face = mesh%elem(i)%face(j)
          allocate(local_id_vert_no_ghost(mesh%face(id_face)%n_vert))
          do k=1, mesh%face(id_face)%n_vert
            local_id_vert_no_ghost(k) = id_vert_no_ghost(mesh%face(id_face)%vert(k))
          end do
          !write(fn_vtu, *) mesh%face(id_face)%n_vert, mesh%face(id_face)%vert(:) - 1
          write(fn_vtu, *) mesh%face(id_face)%n_vert, local_id_vert_no_ghost - 1
          deallocate(local_id_vert_no_ghost)
        end do
      end if
    end do

    write(fn_vtu, *) "</DataArray>"

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

    write(fn_vtu, *) "<DataArray type='Int64' Name='faceoffsets' format='ascii' &
      &RangeMin='", 0, "' RangeMax='", s, "'>"

    s = 0
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        do j=1, mesh%elem(i)%n_faces
          id_face = mesh%elem(i)%face(j)
          s = s + mesh%face(id_face)%n_vert + 1
        end do
        s = s + 1
        write(fn_vtu, *) s
      end if
    end do

    write(fn_vtu, *) "</DataArray>"
    write(fn_vtu, *) "</Cells>"

  end subroutine open_file_vtu

  subroutine close_file_vtu(mesh, filename, fn_vtu, fn_pvtu)
    use mpi
    use mpi_module
    use mesh_module
    implicit none

    type(mesh_type) :: mesh
    character(len=*), intent(in) :: filename
    integer(kind=ENTIER), intent(in) :: fn_vtu, fn_pvtu

    integer(kind=ENTIER) :: me, num_procs, mpi_ierr

    integer(kind=ENTIER) :: i
    character(len=255) :: i_char

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    if( me == 0 ) then
      !PVTU header
      do i=0, num_procs-1
        write(i_char, *) i
        write(fn_pvtu, *) "<Piece Source='"//trim(adjustl(i_char))//&
          "_"//trim(adjustl(filename))//".vtu'/>"
      end do
      write(fn_pvtu, *) "</PUnstructuredGrid>"
      write(fn_pvtu, *) "</VTKFile>"
      close(fn_pvtu)
    end if

    !VTU
    write(fn_vtu, *) "</Piece>"
    write(fn_vtu, *) "</UnstructuredGrid>"
    write(fn_vtu, *) "</VTKFile>"
    close(fn_vtu)
  end subroutine close_file_vtu
end module io_module
