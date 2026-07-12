module mesh_reading_module
  use mpi_module
  use precision_module
  use mesh_module
  use sort_module
  implicit none

  private

  public :: read_mesh_msh
  public :: read_mesh_wasilij
  public :: rescale_mesh

contains
  subroutine read_mesh_wasilij(mesh, meshfile_path, meshfile)
    implicit none

    type(mesh_type), intent(inout) :: mesh
    character(len=255), intent(in) :: meshfile, meshfile_path

    integer(kind=ENTIER) :: funit
    integer(kind=ENTIER) :: i, j
    integer(kind=ENTIER) :: n_faces, dummy
    integer(kind=ENTIER), dimension(10) :: dummy_vert
    character(len=255) :: text

    open (newunit=funit, file=trim(adjustl(meshfile_path))//trim(adjustl(meshfile)), &
      status="old")

    read (funit, '(a)') text
    do while ("#" /= text(1:1))
      read (funit, '(a)') text
    end do

    read (funit, '(a)') text
    do while ("#" /= text(1:1))
      read (funit, '(a)') text
    end do

    read (funit, *) mesh%n_vert
    mesh%n_vert = 2*mesh%n_vert
    read (funit, '(a)') text
    do i = 1, mesh%n_vert/2
      read (funit, *) text
    end do

    read (funit, '(a)') text ! #total number of cells
    read (funit, *) mesh%n_elems
    read (funit, '(a)') text !junk
    mesh%n_interior_elems = mesh%n_elems

    mesh%n_faces = 0
    do i = 1, mesh%n_elems
      read (funit, *) dummy, dummy_vert(:dummy)
      mesh%n_faces = mesh%n_faces + dummy + 2
    end do

    if (mesh%n_elems < 1 .or. mesh%n_faces < 1) then
      print *, achar(27)//"[31m[-] No 3D elements in mesh !"//achar(27)//"[0m", mesh%n_elems, mesh%n_faces
      error stop
    end if

    close (funit)

    allocate (mesh%vert(mesh%n_vert))
    allocate (mesh%elem(mesh%n_elems))

    !!Each face is counted twice this is fixed in build face function
    allocate (mesh%minimal_face(mesh%n_faces))
    mesh%n_boundary_faces = 0

    open (newunit=funit, file=trim(adjustl(meshfile_path))//trim(adjustl(meshfile)), &
      status="old")

    read (funit, '(a)') text
    do while ("#" /= text(1:1))
      read (funit, '(a)') text
    end do

    read (funit, '(a)') text
    do while ("#" /= text(1:1))
      read (funit, '(a)') text
    end do

    read (funit, '(a)') text !n_vert
    read (funit, '(a)') text
    do i = 1, mesh%n_vert/2
      read (funit, *) mesh%vert(i)%coord(1), mesh%vert(i)%coord(2)
      mesh%vert(i)%coord(3) = 0.0_DOUBLE
      mesh%vert(i + mesh%n_vert/2)%coord(:) = mesh%vert(i)%coord(:)
      mesh%vert(i + mesh%n_vert/2)%coord(3) = 0.1_DOUBLE
    end do

    read (funit, '(a)') text ! #total number of cells
    read (funit, '(a)') text ! n_elems
    read (funit, '(a)') text ! #junk

    n_faces = 1
    do i = 1, mesh%n_elems
      read (funit, *) dummy, dummy_vert(:dummy)
      do j = 1, dummy
        dummy_vert(j) = dummy_vert(j) + 1
      end do

      mesh%elem(i)%n_vert = 2*dummy
      allocate (mesh%elem(i)%vert(2*dummy))
      do j = 1, dummy
        mesh%elem(i)%vert(j) = dummy_vert(j)
        mesh%elem(i)%vert(j + dummy) = mesh%n_vert/2 + dummy_vert(j)
      end do

      mesh%elem(i)%n_faces = dummy + 2
      allocate (mesh%elem(i)%face(dummy + 2))

      !Bottom face
      mesh%elem(i)%face(1) = n_faces
      allocate (mesh%minimal_face(n_faces)%vert(dummy))
      mesh%minimal_face(n_faces)%left_neigh = i
      mesh%minimal_face(n_faces)%left_neigh_face = 1
      mesh%minimal_face(n_faces)%n_vert = dummy
      do j = 1, dummy
        mesh%minimal_face(n_faces)%vert(j) = dummy_vert(j)
      end do
      n_faces = n_faces + 1

      mesh%elem(i)%face(2) = n_faces
      allocate (mesh%minimal_face(n_faces)%vert(dummy))
      mesh%minimal_face(n_faces)%left_neigh = i
      mesh%minimal_face(n_faces)%left_neigh_face = 2
      mesh%minimal_face(n_faces)%n_vert = dummy
      do j = 1, dummy
        mesh%minimal_face(n_faces)%vert(j) = mesh%n_vert/2 + dummy_vert(j)
      end do
      n_faces = n_faces + 1

      !Side faces
      do j = 1, dummy
        mesh%elem(i)%face(2 + j) = n_faces
        mesh%minimal_face(n_faces)%left_neigh = i
        mesh%minimal_face(n_faces)%n_vert = 4
        mesh%minimal_face(n_faces)%left_neigh_face = 2 + j
        allocate (mesh%minimal_face(n_faces)%vert(4))
        mesh%minimal_face(n_faces)%vert(1) = dummy_vert(1 + mod(j - 1, dummy))
        mesh%minimal_face(n_faces)%vert(2) = dummy_vert(1 + mod(j, dummy))
        mesh%minimal_face(n_faces)%vert(3) = mesh%n_vert/2 + dummy_vert(1 + mod(j, dummy))
        mesh%minimal_face(n_faces)%vert(4) = mesh%n_vert/2 + dummy_vert(1 + mod(j - 1, dummy))
        n_faces = n_faces + 1
      end do
    end do

    close (funit)
  end subroutine read_mesh_wasilij

  subroutine read_mesh_msh(mesh, meshfile_path, meshfile, &
      n_bc, bc_name, me, num_procs, mpi_send_recv)
    implicit none

    type(mesh_type), intent(inout) :: mesh
    character(len=255), intent(in) :: meshfile, meshfile_path
    integer(kind=ENTIER), intent(in) :: num_procs, me
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv
    integer(kind=ENTIER), intent(in) :: n_bc
    character(len=255), dimension(:), intent(in) :: bc_name

    integer(kind=ENTIER) :: funit, ppos
    character(len=255) :: text, me_str, mpi_meshfile
    real(kind=DOUBLE) :: version

    if (num_procs > 1) then
      ppos = scan(trim(meshfile), ".", BACK=.TRUE.)
      write (me_str, *) me + 1

      mpi_meshfile = trim(adjustl(meshfile(1:ppos - 1)))//"_"// &
        trim(adjustl(me_str))//trim(adjustl(meshfile(ppos:)))
      open (newunit=funit, file=trim(adjustl(meshfile_path))//trim(adjustl(mpi_meshfile)), &
        status="old")
      read (funit, *) text
      read (funit, *) version, text
      close (funit)
    else
      open (newunit=funit, file=trim(adjustl(meshfile_path))//trim(adjustl(meshfile)), &
        status="old")
      read (funit, *) text
      read (funit, *) version, text
      close (funit)
    end if

    if (floor(version) == 4) then
      if (num_procs > 1) then
        call read_mesh_msh_4(mesh, meshfile_path, mpi_meshfile, &
          n_bc, bc_name, me, num_procs, mpi_send_recv)
      else
        call read_mesh_msh_4(mesh, meshfile_path, meshfile, &
          n_bc, bc_name, me, num_procs, mpi_send_recv)
      end if
    else if (floor(version) == 2) then
      if (num_procs > 1) then
        print *, "More than one MPI process with msh file version 2 not implemented !"
        print *, "Use MSH4 format with one file per partition."
        error stop
      else
        call read_mesh_msh_2(mesh, meshfile_path, meshfile, n_bc, bc_name)
      end if
    end if
  end subroutine read_mesh_msh

  subroutine read_mesh_msh_4(mesh, meshfile_path, meshfile, &
      n_bc, bc_name, me, num_procs, mpi_send_recv)
    use mpi
    implicit none

    type(mesh_type), intent(inout) :: mesh
    character(len=255), intent(in) :: meshfile, meshfile_path
    integer(kind=ENTIER), intent(in) :: num_procs, me
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv
    integer(kind=ENTIER), intent(in) :: n_bc
    character(len=255), dimension(:), intent(in) :: bc_name
    integer(kind=ENTIER), dimension(n_bc) :: bc_tag

    integer(kind=ENTIER) :: funit
    integer(kind=ENTIER) :: i, j, k, iloc, l, iloc2
    integer(kind=ENTIER) :: n_nodes_in_entity, n_entity_blocks
    integer(kind=ENTIER) :: n_partitions, n_ghost_entities
    integer(kind=ENTIER) :: n_entity_points, n_entity_curves
    integer(kind=ENTIER) :: n_entity_surfs, n_entity_vols
    integer(kind=ENTIER) :: entity_dim, entity_elem_kind, n_elems_in_entity
    integer(kind=ENTIER) :: n_elems, n_faces, n_boundary_faces, n_bc_tmp
    integer(kind=ENTIER) :: dummy, tag, dummy_n_tag2, dummy_n_tag3
    real(kind=DOUBLE) :: dummy_real
    integer(kind=ENTIER) :: entity_tag, entity_phys_tag
    integer(kind=ENTIER), dimension(10) :: dummy_vert
    integer(kind=ENTIER), dimension(100) :: dummy_tag, dummy_tag2, dummy_tag3
    integer(kind=ENTIER), dimension(:), allocatable :: entity_to_phys_tag
    integer(kind=ENTIER), dimension(:), allocatable :: entity_loc_tag
    integer(kind=ENTIER), dimension(:), allocatable :: node_list_to_read
    character(len=255) :: text
    logical, dimension(:), allocatable :: do_read_node

    integer(kind=ENTIER), dimension(:), allocatable :: send_elem_sizes, recv_elem_sizes
    integer(kind=ENTIER) :: owner, n_receiver, n_ghost_elems, id_glob_elem, id_loc_elem
    integer(kind=ENTIER), dimension(:), allocatable :: receiver
    integer(kind=ENTIER), dimension(:, :), allocatable :: id_vert, id_elems
    integer(kind=ENTIER) :: kloc, iloc_send, iloc_recv, iloc_vert

    real(kind=DOUBLE), dimension(3) :: entity_coord_min, entity_coord_max

    !!https://gmsh.info/doc/texinfo/gmsh.html#Gmsh-file-formats

    !! http://www.manpagez.com/info/gmsh/gmsh-2.2.6/gmsh_63.php
    integer(kind=ENTIER), dimension(31), parameter :: n_vert_kind = (/ &
      2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, &
      27, 18, 14, 1, 8, 20, 15, 13, 9, 10, 12, &
      15, 15, 21, 4, 5, 6, 20, 35, 56/)
    integer(kind=ENTIER), dimension(31), parameter :: n_face_kind = (/ &
      0, 0, 0, 4, 6, 5, 5, 0, 0, 0, 0, &
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      0, 0, 0, 0, 0, 0, 0, 0, 0/)
    integer(kind=ENTIER), dimension(31), parameter :: n_face_vert_kind = (/ &
      0, 0, 0, 4*3, 6*4, 3*4 + 2*3, 1*4 + 4*3, 0, 0, 0, 0, &
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      0, 0, 0, 0, 0, 0, 0, 0, 0/)

    open (newunit=funit, file=trim(adjustl(meshfile_path))//trim(adjustl(meshfile)), &
      status="old")

    if (n_bc /= 0) then
      read (funit, '(a)') text
      do while ("$PhysicalNames" /= text(1:14))
        read (funit, '(a)') text
        if ("$Nodes" == text(1:6)) then
          print *, "No Physical names inside mesh!!!"
          error stop
        end if
      end do

      !!Read PhysicalNames of BC
      read (funit, *) n_bc_tmp
      do i = 1, n_bc_tmp
        read (funit, *) dummy, tag, text
        do j = 1, n_bc
          if (trim(adjustl(bc_name(j))) == trim(adjustl(text))) then
            bc_tag(j) = tag
            exit
          end if
        end do
      end do
    end if

    if (num_procs > 1) then
      read (funit, '(a)') text
      do while ("$PartitionedEntities" /= text(1:20))
        read (funit, '(a)') text
      end do

      read (funit, *) n_partitions
      read (funit, *) n_ghost_entities
      do j = 1, n_ghost_entities
        read (funit, '(a)') text
      end do

      read (funit, *) n_entity_points, n_entity_curves, n_entity_surfs, n_entity_vols

      do j = 1, n_entity_points
        read (funit, '(a)') text
      end do

      do j = 1, n_entity_curves
        read (funit, '(a)') text
      end do

      allocate (entity_loc_tag(n_entity_surfs))
      allocate (entity_to_phys_tag(n_entity_surfs))
      do j = 1, n_entity_surfs
        read (funit, *) entity_loc_tag(j), &
          dummy, dummy, & !Parent dim, Parent tag
          dummy, dummy_tag(:dummy), & !Numpartitions, partitiontags
          entity_coord_min(1), entity_coord_min(2), entity_coord_min(3), &
          entity_coord_max(1), entity_coord_max(2), entity_coord_max(3), &
          dummy_n_tag2, dummy_tag2(:dummy_n_tag2), dummy_n_tag3, dummy_tag3(:dummy_n_tag3)

        if (dummy_n_tag2 > 0) then
          ! Takes the first phys tag that is found in bc_tag
          entity_to_phys_tag(j) = 0
          do k=1, dummy_n_tag2
            if( entity_to_phys_tag(j) /= 0 ) exit
            do l=1, n_bc
              if( dummy_tag2(k) == bc_tag(l) ) then
                entity_to_phys_tag(j) = dummy_tag2(k)
                exit
              end if
            end do
          end do
        else
          entity_to_phys_tag(j) = 0
        end if
      end do
    else
      !Reading Entities
      read (funit, '(a)') text
      do while ("$Entities" /= text(1:9))
        read (funit, '(a)') text
      end do

      read (funit, *) n_entity_points, n_entity_curves, n_entity_surfs, n_entity_vols

      do j = 1, n_entity_points
        read (funit, '(a)') text
      end do

      do j = 1, n_entity_curves
        read (funit, '(a)') text
      end do

      allocate (entity_loc_tag(n_entity_surfs))
      allocate (entity_to_phys_tag(n_entity_surfs))
      do j = 1, n_entity_surfs
        read (funit, *) entity_loc_tag(j), &
          entity_coord_min(1), entity_coord_min(2), entity_coord_min(3), &
          entity_coord_max(1), entity_coord_max(2), entity_coord_max(3), &
          dummy_n_tag2, dummy_tag2(:dummy_n_tag2)

        if (dummy_n_tag2 > 0) then
          ! Takes the first phys tag that is found in bc_tag
          entity_to_phys_tag(j) = 0
          do k=1, dummy_n_tag2
            if( entity_to_phys_tag(j) /= 0 ) exit
            do l=1, n_bc
              if( dummy_tag2(k) == bc_tag(l) ) then
                entity_to_phys_tag(j) = dummy_tag2(k)
                exit
              end if
            end do
          end do
        else
          entity_to_phys_tag(j) = 0
        end if
      end do
    end if

    do j = 1, n_entity_vols
      read (funit, '(a)') text
    end do

    !Reading Elements
    do while ("$Elements" /= text(1:9))
      read (funit, '(a)') text
    end do

    read (funit, *) n_entity_blocks, n_elems, dummy, dummy

    mesh%n_elems = 0
    mesh%n_faces = 0
    mesh%n_boundary_faces = 0
    do i = 1, n_entity_blocks
      read (funit, *) entity_dim, entity_tag, entity_elem_kind, n_elems_in_entity

      if (entity_dim == 0) then
        do j = 1, n_elems_in_entity
          read (funit, *) text
        end do
      else if (entity_dim == 1) then
        do j = 1, n_elems_in_entity
          read (funit, *) text
        end do
      else if (entity_dim == 2) then
        if (entity_elem_kind == 2 .or. entity_elem_kind == 3) then
          do j = 1, n_elems_in_entity
            read (funit, *) dummy, dummy_vert(:n_vert_kind(entity_elem_kind))
            mesh%n_boundary_faces = mesh%n_boundary_faces + 1
          end do
        else
          do j = 1, n_elems_in_entity
            read (funit, *) text
          end do
        end if
      else if (entity_dim == 3) then
        if ((entity_elem_kind == 4 .or. entity_elem_kind == 5) .or. &
          (entity_elem_kind == 6 .or. entity_elem_kind == 7)) then
          do j = 1, n_elems_in_entity
            read (funit, *) text
            mesh%n_faces = mesh%n_faces + n_face_kind(entity_elem_kind)
            mesh%n_elems = mesh%n_elems + 1
          end do
        else
          do j = 1, n_elems_in_entity
            read (funit, *) text
          end do
        end if
      end if
    end do

    !Read ghost elems to allocate
    if (num_procs > 1) then
      ! send_elem_sizes
      allocate (send_elem_sizes(num_procs))
      send_elem_sizes = 0
      allocate (recv_elem_sizes(num_procs))
      recv_elem_sizes = 0

      allocate (receiver(num_procs))

      read (funit, '(a)') text
      do while ("$GhostElements" /= text(1:14))
        read (funit, '(a)') text
      end do

      read (funit, *) n_ghost_elems

      do j = 1, n_ghost_elems
        read (funit, *) dummy, owner, n_receiver, receiver(:n_receiver)
        if (me + 1 == owner) then
          do i = 1, n_receiver
            send_elem_sizes(receiver(i)) = send_elem_sizes(receiver(i)) + 1
          end do
        else
          recv_elem_sizes(owner) = recv_elem_sizes(owner) + 1
        end if
      end do

      ! allocate proper sizes for send/recv operations
      mpi_send_recv%n_mpi_send_neigh = 0
      mpi_send_recv%n_mpi_recv_neigh = 0
      do i = 1, num_procs
        if (send_elem_sizes(i) > 0) mpi_send_recv%n_mpi_send_neigh = mpi_send_recv%n_mpi_send_neigh + 1
        if (recv_elem_sizes(i) > 0) mpi_send_recv%n_mpi_recv_neigh = mpi_send_recv%n_mpi_recv_neigh + 1
      end do

      allocate (mpi_send_recv%mpi_send_neigh(mpi_send_recv%n_mpi_send_neigh))
      allocate (mpi_send_recv%mpi_recv_neigh(mpi_send_recv%n_mpi_recv_neigh))
      allocate (mpi_send_recv%mpi_reqsend(mpi_send_recv%n_mpi_send_neigh))
      allocate (mpi_send_recv%mpi_reqrecv(mpi_send_recv%n_mpi_recv_neigh))

      allocate (mpi_send_recv%mpi_sendstat(MPI_STATUS_SIZE, mpi_send_recv%n_mpi_send_neigh))
      allocate (mpi_send_recv%mpi_recvstat(MPI_STATUS_SIZE, mpi_send_recv%n_mpi_recv_neigh))

      iloc_send = 1
      iloc_recv = 1
      do i = 1, num_procs
        if (send_elem_sizes(i) > 0) then
          mpi_send_recv%mpi_send_neigh(iloc_send)%partition_id = i - 1
          mpi_send_recv%mpi_send_neigh(iloc_send)%n_elems = send_elem_sizes(i)
          allocate (mpi_send_recv%mpi_send_neigh(iloc_send)%elem_id(send_elem_sizes(i)))
          mpi_send_recv%mpi_send_neigh(iloc_send)%elem_id(:) = -1
          iloc_send = iloc_send + 1
        end if
        if (recv_elem_sizes(i) > 0) then
          mpi_send_recv%mpi_recv_neigh(iloc_recv)%partition_id = i - 1
          mpi_send_recv%mpi_recv_neigh(iloc_recv)%n_elems = recv_elem_sizes(i)
          allocate (mpi_send_recv%mpi_recv_neigh(iloc_recv)%elem_id(recv_elem_sizes(i)))
          mpi_send_recv%mpi_recv_neigh(iloc_recv)%elem_id(:) = -1
          iloc_recv = iloc_recv + 1
        end if
      end do
    end if

    close (funit)

    if (mesh%n_elems < 1 .or. mesh%n_faces < 1) then
      print *, achar(27)//"[31m[-] No 3D elements in mesh !"//achar(27)//"[0m", mesh%n_elems, mesh%n_faces
      error stop
    end if

    mesh%n_interior_elems = mesh%n_elems
    allocate (mesh%elem(mesh%n_elems))
    allocate (id_elems(mesh%n_elems, 2))

    !!Each face is counted twice this is fixed in build face function
    allocate (mesh%minimal_face(mesh%n_faces))
    allocate (mesh%boundary_face(mesh%n_boundary_faces))

    open (newunit=funit, file=trim(adjustl(meshfile_path))//trim(adjustl(meshfile)), &
      status="old")

    !Reading Elements
    do while ("$Elements" /= text(1:9))
      read (funit, '(a)') text
    end do

    read (funit, *) n_entity_blocks, dummy, dummy, dummy

    iloc = 1
    n_faces = 1
    n_boundary_faces = 1
    do i = 1, n_entity_blocks

      read (funit, *) entity_dim, entity_tag, entity_elem_kind, n_elems_in_entity

      if (entity_dim == 0) then
        do j = 1, n_elems_in_entity
          read (funit, *) text
        end do
      else if (entity_dim == 1) then
        !!Skip
        do j = 1, n_elems_in_entity
          read (funit, *) text
        end do
      else if (entity_dim == 2) then
        if (entity_elem_kind == 2 .or. entity_elem_kind == 3) then

          !Retreive physical tag from entity tag
          entity_phys_tag = -1
          do j = 1, n_entity_surfs
            if (entity_loc_tag(j) == entity_tag) then
              entity_phys_tag = entity_to_phys_tag(j)
            end if
          end do

          if (entity_phys_tag < 0) then
            print *, "Entity phys tag not found !"
            error stop
          end if

          do j = 1, n_elems_in_entity
            read (funit, *) dummy, dummy_vert(:n_vert_kind(entity_elem_kind))

            mesh%boundary_face(n_boundary_faces)%tag = 0
            do k = 1, n_bc
              if (bc_tag(k) == entity_phys_tag) then
                mesh%boundary_face(n_boundary_faces)%tag = k
              end if
            end do

            mesh%boundary_face(n_boundary_faces)%n_vert = n_vert_kind(entity_elem_kind)
            allocate (mesh%boundary_face(n_boundary_faces)%vert(n_vert_kind(entity_elem_kind)))
            do k = 1, n_vert_kind(entity_elem_kind)
              mesh%boundary_face(n_boundary_faces)%vert(k) = dummy_vert(k)
            end do
            n_boundary_faces = n_boundary_faces + 1
          end do

        else
          !Skip
          do j = 1, n_elems_in_entity
            read (funit, *) text
          end do
        end if

      else if (entity_dim == 3) then
        if ((entity_elem_kind == 4 .or. entity_elem_kind == 5) .or. &
          (entity_elem_kind == 6 .or. entity_elem_kind == 7)) then
          do k = 1, n_elems_in_entity
            read (funit, *) id_elems(iloc, 1), dummy_vert(:n_vert_kind(entity_elem_kind))
            id_elems(iloc, 2) = iloc

            mesh%elem(iloc)%elem_kind = entity_elem_kind

            mesh%elem(iloc)%n_vert = n_vert_kind(entity_elem_kind)
            mesh%elem(iloc)%n_faces = n_face_kind(entity_elem_kind)
            allocate (mesh%elem(iloc)%vert(n_vert_kind(entity_elem_kind)))

            do j = 1, n_vert_kind(entity_elem_kind)
              mesh%elem(iloc)%vert(j) = dummy_vert(j)
            end do

            do j = 1, n_face_kind(entity_elem_kind)
              mesh%minimal_face(n_faces + j - 1)%left_neigh = iloc
              mesh%minimal_face(n_faces + j - 1)%left_neigh_face = j
            end do

            allocate (mesh%elem(iloc)%face(n_face_kind(entity_elem_kind)))
            if (entity_elem_kind == 4) then !!tet
              mesh%minimal_face(n_faces)%n_vert = 3
              allocate (mesh%minimal_face(n_faces)%vert(3))
              mesh%elem(iloc)%face(1) = n_faces
              mesh%minimal_face(n_faces)%vert(1) = mesh%elem(iloc)%vert(1)
              mesh%minimal_face(n_faces)%vert(2) = mesh%elem(iloc)%vert(3)
              mesh%minimal_face(n_faces)%vert(3) = mesh%elem(iloc)%vert(2)

              mesh%minimal_face(n_faces + 1)%n_vert = 3
              allocate (mesh%minimal_face(n_faces + 1)%vert(3))
              mesh%elem(iloc)%face(2) = n_faces + 1
              mesh%minimal_face(n_faces + 1)%vert(1) = mesh%elem(iloc)%vert(1)
              mesh%minimal_face(n_faces + 1)%vert(2) = mesh%elem(iloc)%vert(2)
              mesh%minimal_face(n_faces + 1)%vert(3) = mesh%elem(iloc)%vert(4)

              mesh%minimal_face(n_faces + 2)%n_vert = 3
              allocate (mesh%minimal_face(n_faces + 2)%vert(3))
              mesh%elem(iloc)%face(3) = n_faces + 2
              mesh%minimal_face(n_faces + 2)%vert(1) = mesh%elem(iloc)%vert(1)
              mesh%minimal_face(n_faces + 2)%vert(2) = mesh%elem(iloc)%vert(4)
              mesh%minimal_face(n_faces + 2)%vert(3) = mesh%elem(iloc)%vert(3)

              mesh%minimal_face(n_faces + 3)%n_vert = 3
              allocate (mesh%minimal_face(n_faces + 3)%vert(3))
              mesh%elem(iloc)%face(4) = n_faces + 3
              mesh%minimal_face(n_faces + 3)%vert(1) = mesh%elem(iloc)%vert(2)
              mesh%minimal_face(n_faces + 3)%vert(2) = mesh%elem(iloc)%vert(3)
              mesh%minimal_face(n_faces + 3)%vert(3) = mesh%elem(iloc)%vert(4)

            else if (entity_elem_kind == 5) then !! Hexa
              mesh%minimal_face(n_faces)%n_vert = 4
              allocate (mesh%minimal_face(n_faces)%vert(4))
              mesh%elem(iloc)%face(1) = n_faces
              mesh%minimal_face(n_faces)%vert(1) = mesh%elem(iloc)%vert(1)
              mesh%minimal_face(n_faces)%vert(2) = mesh%elem(iloc)%vert(4)
              mesh%minimal_face(n_faces)%vert(3) = mesh%elem(iloc)%vert(3)
              mesh%minimal_face(n_faces)%vert(4) = mesh%elem(iloc)%vert(2)

              mesh%minimal_face(n_faces + 1)%n_vert = 4
              allocate (mesh%minimal_face(n_faces + 1)%vert(4))
              mesh%elem(iloc)%face(2) = n_faces + 1
              mesh%minimal_face(n_faces + 1)%vert(1) = mesh%elem(iloc)%vert(3)
              mesh%minimal_face(n_faces + 1)%vert(2) = mesh%elem(iloc)%vert(4)
              mesh%minimal_face(n_faces + 1)%vert(3) = mesh%elem(iloc)%vert(8)
              mesh%minimal_face(n_faces + 1)%vert(4) = mesh%elem(iloc)%vert(7)

              mesh%minimal_face(n_faces + 2)%n_vert = 4
              allocate (mesh%minimal_face(n_faces + 2)%vert(4))
              mesh%elem(iloc)%face(3) = n_faces + 2
              mesh%minimal_face(n_faces + 2)%vert(1) = mesh%elem(iloc)%vert(7)
              mesh%minimal_face(n_faces + 2)%vert(2) = mesh%elem(iloc)%vert(8)
              mesh%minimal_face(n_faces + 2)%vert(3) = mesh%elem(iloc)%vert(5)
              mesh%minimal_face(n_faces + 2)%vert(4) = mesh%elem(iloc)%vert(6)

              mesh%minimal_face(n_faces + 3)%n_vert = 4
              allocate (mesh%minimal_face(n_faces + 3)%vert(4))
              mesh%elem(iloc)%face(4) = n_faces + 3
              mesh%minimal_face(n_faces + 3)%vert(1) = mesh%elem(iloc)%vert(1)
              mesh%minimal_face(n_faces + 3)%vert(2) = mesh%elem(iloc)%vert(2)
              mesh%minimal_face(n_faces + 3)%vert(3) = mesh%elem(iloc)%vert(6)
              mesh%minimal_face(n_faces + 3)%vert(4) = mesh%elem(iloc)%vert(5)

              mesh%minimal_face(n_faces + 4)%n_vert = 4
              allocate (mesh%minimal_face(n_faces + 4)%vert(4))
              mesh%elem(iloc)%face(5) = n_faces + 4
              mesh%minimal_face(n_faces + 4)%vert(1) = mesh%elem(iloc)%vert(1)
              mesh%minimal_face(n_faces + 4)%vert(2) = mesh%elem(iloc)%vert(5)
              mesh%minimal_face(n_faces + 4)%vert(3) = mesh%elem(iloc)%vert(8)
              mesh%minimal_face(n_faces + 4)%vert(4) = mesh%elem(iloc)%vert(4)

              mesh%minimal_face(n_faces + 5)%n_vert = 4
              allocate (mesh%minimal_face(n_faces + 5)%vert(4))
              mesh%elem(iloc)%face(6) = n_faces + 5
              mesh%minimal_face(n_faces + 5)%vert(1) = mesh%elem(iloc)%vert(2)
              mesh%minimal_face(n_faces + 5)%vert(2) = mesh%elem(iloc)%vert(3)
              mesh%minimal_face(n_faces + 5)%vert(3) = mesh%elem(iloc)%vert(7)
              mesh%minimal_face(n_faces + 5)%vert(4) = mesh%elem(iloc)%vert(6)

            else if (entity_elem_kind == 6) then !! Prism
              mesh%minimal_face(n_faces)%n_vert = 3
              allocate (mesh%minimal_face(n_faces)%vert(3))
              mesh%elem(iloc)%face(1) = n_faces
              mesh%minimal_face(n_faces)%vert(1) = mesh%elem(iloc)%vert(1)
              mesh%minimal_face(n_faces)%vert(2) = mesh%elem(iloc)%vert(3)
              mesh%minimal_face(n_faces)%vert(3) = mesh%elem(iloc)%vert(2)

              mesh%minimal_face(n_faces + 1)%n_vert = 3
              allocate (mesh%minimal_face(n_faces + 1)%vert(3))
              mesh%elem(iloc)%face(2) = n_faces + 1
              mesh%minimal_face(n_faces + 1)%vert(1) = mesh%elem(iloc)%vert(4)
              mesh%minimal_face(n_faces + 1)%vert(2) = mesh%elem(iloc)%vert(5)
              mesh%minimal_face(n_faces + 1)%vert(3) = mesh%elem(iloc)%vert(6)

              mesh%minimal_face(n_faces + 2)%n_vert = 4
              allocate (mesh%minimal_face(n_faces + 2)%vert(4))
              mesh%elem(iloc)%face(3) = n_faces + 2
              mesh%minimal_face(n_faces + 2)%vert(1) = mesh%elem(iloc)%vert(1)
              mesh%minimal_face(n_faces + 2)%vert(2) = mesh%elem(iloc)%vert(2)
              mesh%minimal_face(n_faces + 2)%vert(3) = mesh%elem(iloc)%vert(5)
              mesh%minimal_face(n_faces + 2)%vert(4) = mesh%elem(iloc)%vert(4)

              mesh%minimal_face(n_faces + 3)%n_vert = 4
              allocate (mesh%minimal_face(n_faces + 3)%vert(4))
              mesh%elem(iloc)%face(4) = n_faces + 3
              mesh%minimal_face(n_faces + 3)%vert(1) = mesh%elem(iloc)%vert(2)
              mesh%minimal_face(n_faces + 3)%vert(2) = mesh%elem(iloc)%vert(3)
              mesh%minimal_face(n_faces + 3)%vert(3) = mesh%elem(iloc)%vert(6)
              mesh%minimal_face(n_faces + 3)%vert(4) = mesh%elem(iloc)%vert(5)

              mesh%minimal_face(n_faces + 4)%n_vert = 4
              allocate (mesh%minimal_face(n_faces + 4)%vert(4))
              mesh%elem(iloc)%face(5) = n_faces + 4
              mesh%minimal_face(n_faces + 4)%vert(1) = mesh%elem(iloc)%vert(1)
              mesh%minimal_face(n_faces + 4)%vert(2) = mesh%elem(iloc)%vert(4)
              mesh%minimal_face(n_faces + 4)%vert(3) = mesh%elem(iloc)%vert(6)
              mesh%minimal_face(n_faces + 4)%vert(4) = mesh%elem(iloc)%vert(3)

            else if (entity_elem_kind == 7) then !! Pyramid
              mesh%minimal_face(n_faces)%n_vert = 4
              allocate (mesh%minimal_face(n_faces)%vert(4))
              mesh%elem(iloc)%face(1) = n_faces
              mesh%minimal_face(n_faces)%vert(1) = mesh%elem(iloc)%vert(1)
              mesh%minimal_face(n_faces)%vert(2) = mesh%elem(iloc)%vert(4)
              mesh%minimal_face(n_faces)%vert(3) = mesh%elem(iloc)%vert(3)
              mesh%minimal_face(n_faces)%vert(4) = mesh%elem(iloc)%vert(2)

              mesh%minimal_face(n_faces + 1)%n_vert = 3
              allocate (mesh%minimal_face(n_faces + 1)%vert(3))
              mesh%elem(iloc)%face(2) = n_faces + 1
              mesh%minimal_face(n_faces + 1)%vert(1) = mesh%elem(iloc)%vert(1)
              mesh%minimal_face(n_faces + 1)%vert(2) = mesh%elem(iloc)%vert(2)
              mesh%minimal_face(n_faces + 1)%vert(3) = mesh%elem(iloc)%vert(5)

              mesh%minimal_face(n_faces + 2)%n_vert = 3
              allocate (mesh%minimal_face(n_faces + 2)%vert(3))
              mesh%elem(iloc)%face(3) = n_faces + 2
              mesh%minimal_face(n_faces + 2)%vert(1) = mesh%elem(iloc)%vert(2)
              mesh%minimal_face(n_faces + 2)%vert(2) = mesh%elem(iloc)%vert(3)
              mesh%minimal_face(n_faces + 2)%vert(3) = mesh%elem(iloc)%vert(5)

              mesh%minimal_face(n_faces + 3)%n_vert = 3
              allocate (mesh%minimal_face(n_faces + 3)%vert(3))
              mesh%elem(iloc)%face(4) = n_faces + 3
              mesh%minimal_face(n_faces + 3)%vert(1) = mesh%elem(iloc)%vert(3)
              mesh%minimal_face(n_faces + 3)%vert(2) = mesh%elem(iloc)%vert(4)
              mesh%minimal_face(n_faces + 3)%vert(3) = mesh%elem(iloc)%vert(5)

              mesh%minimal_face(n_faces + 4)%n_vert = 3
              allocate (mesh%minimal_face(n_faces + 4)%vert(3))
              mesh%elem(iloc)%face(5) = n_faces + 4
              mesh%minimal_face(n_faces + 4)%vert(1) = mesh%elem(iloc)%vert(1)
              mesh%minimal_face(n_faces + 4)%vert(2) = mesh%elem(iloc)%vert(5)
              mesh%minimal_face(n_faces + 4)%vert(3) = mesh%elem(iloc)%vert(4)
            end if

            iloc = iloc + 1
            n_faces = n_faces + n_face_kind(entity_elem_kind)
          end do
        else
          do j = 1, n_elems_in_entity
            read (funit, *) text
          end do
        end if
      end if
    end do

    if (num_procs > 1) then

      read (funit, '(a)') text
      do while ("$GhostElements" /= text(1:14))
        read (funit, '(a)') text
      end do

      read (funit, *) n_ghost_elems

      do j = 1, n_ghost_elems
        read (funit, *) id_glob_elem, owner, n_receiver, receiver(:n_receiver)
        if (me + 1 == owner) then
          do i = 1, n_receiver
            do k = 1, mpi_send_recv%n_mpi_send_neigh
              !Matching receiver
              if (mpi_send_recv%mpi_send_neigh(k)%partition_id + 1 == receiver(i)) then
                !Add elem to vector
                do l = 1, mpi_send_recv%mpi_send_neigh(k)%n_elems
                  if (mpi_send_recv%mpi_send_neigh(k)%elem_id(l) < 0) then
                    mpi_send_recv%mpi_send_neigh(k)%elem_id(l) = id_glob_elem
                    exit
                  end if
                end do
                exit
              end if
            end do
          end do
        else
          do k = 1, mpi_send_recv%n_mpi_recv_neigh
            !Matching sender
            if (mpi_send_recv%mpi_recv_neigh(k)%partition_id + 1 == owner) then
              do l = 1, mpi_send_recv%mpi_recv_neigh(k)%n_elems
                if (mpi_send_recv%mpi_recv_neigh(k)%elem_id(l) < 0) then
                  mpi_send_recv%mpi_recv_neigh(k)%elem_id(l) = id_glob_elem
                  exit
                end if
              end do
              exit
            end if
          end do
        end if
      end do

    end if

    close (funit)

    !Make a list of global node numbers to read
    !Should use C++ set to optimize
    iloc = 0
    do i = 1, mesh%n_elems
      iloc = iloc + mesh%elem(i)%n_vert
    end do

    allocate (node_list_to_read(iloc)) !!Too large, should be reduced

    iloc = 1
    do i = 1, mesh%n_elems
      do j = 1, mesh%elem(i)%n_vert
        node_list_to_read(iloc) = mesh%elem(i)%vert(j)
        iloc = iloc + 1
      end do
    end do

    call qsort_key(iloc - 1, node_list_to_read)

    mesh%n_vert = 1
    do i = 2, iloc - 1
      if (node_list_to_read(i) /= node_list_to_read(i - 1)) then
        mesh%n_vert = mesh%n_vert + 1
      end if
    end do

    open (newunit=funit, file=trim(adjustl(meshfile_path))//trim(adjustl(meshfile)), &
      status="old")

    !Reading Nodes
    read (funit, '(a)') text
    do while ("$Nodes" /= text(1:6))
      read (funit, '(a)') text
    end do

    read (funit, *) n_entity_blocks, dummy, dummy, dummy

    allocate (mesh%vert(mesh%n_vert))
    allocate (id_vert(mesh%n_vert, 2))
    id_vert(:, :) = -1
    iloc = 1
    do i = 1, n_entity_blocks
      read (funit, *) entity_dim, dummy, dummy, n_nodes_in_entity !Entity stuff

      if (n_nodes_in_entity > 0) then
        allocate (do_read_node(n_nodes_in_entity))
        iloc2 = iloc
        do j = 1, n_nodes_in_entity
          read (funit, *) dummy !Global node tag
          call get_loc_id(size(node_list_to_read), node_list_to_read, dummy, iloc_vert)
          if (iloc_vert > 0) then
            id_vert(iloc, 1) = dummy !Global node tag
            id_vert(iloc, 2) = iloc !Local node tag
            iloc = iloc + 1
            do_read_node(j) = .TRUE.
          else
            do_read_node(j) = .FALSE.
          end if
        end do

        iloc = iloc2 !Reset iloc back to start

        do j = 1, n_nodes_in_entity
          if (do_read_node(j)) then
            read (funit, *) mesh%vert(iloc)%coord(1), mesh%vert(iloc)%coord(2), mesh%vert(iloc)%coord(3)
            iloc = iloc + 1
          else
            read (funit, *) dummy_real, dummy_real, dummy_real
          end if
        end do
        deallocate (do_read_node)
      end if
    end do

    close (funit)

    if( num_procs == 1 ) then
      !!Renumber nodes
      call qsort_int(mesh%n_vert, id_vert(:, 1), id_vert(:, 2))

      do i=1, mesh%n_vert
        mesh%vert(id_vert(i, 2))%id_glob = id_vert(i, 1)
      end do

      do i = 1, mesh%n_elems
        do j = 1, mesh%elem(i)%n_vert
          call get_loc_id(mesh%n_vert, id_vert, mesh%elem(i)%vert(j), kloc)
          mesh%elem(i)%vert(j) = id_vert(kloc, 2)
        end do
      end do

      do i = 1, mesh%n_faces
        do j = 1, mesh%minimal_face(i)%n_vert
          call get_loc_id(mesh%n_vert, id_vert, mesh%minimal_face(i)%vert(j), kloc)
          mesh%minimal_face(i)%vert(j) = id_vert(kloc, 2)
        end do
      end do

      do i = 1, mesh%n_boundary_faces
        do j = 1, mesh%boundary_face(i)%n_vert
          call get_loc_id(mesh%n_vert, id_vert, mesh%boundary_face(i)%vert(j), kloc)
          mesh%boundary_face(i)%vert(j) = id_vert(kloc, 2)
        end do
      end do
    else 
      !!Renumber nodes
      call qsort_int(mesh%n_vert, id_vert(:, 1), id_vert(:, 2))

      do i=1, mesh%n_vert
        mesh%vert(id_vert(i, 2))%id_glob = id_vert(i, 1)
      end do

      do i = 1, mesh%n_elems
        do j = 1, mesh%elem(i)%n_vert
          call get_loc_id(mesh%n_vert, id_vert, mesh%elem(i)%vert(j), kloc)
          mesh%elem(i)%vert(j) = id_vert(kloc, 2)
        end do
      end do

      do i = 1, mesh%n_faces
        do j = 1, mesh%minimal_face(i)%n_vert
          call get_loc_id(mesh%n_vert, id_vert, mesh%minimal_face(i)%vert(j), kloc)
          mesh%minimal_face(i)%vert(j) = id_vert(kloc, 2)
        end do
      end do

      do i = 1, mesh%n_boundary_faces
        do j = 1, mesh%boundary_face(i)%n_vert
          call get_loc_id(mesh%n_vert, id_vert, mesh%boundary_face(i)%vert(j), kloc)
          mesh%boundary_face(i)%vert(j) = id_vert(kloc, 2)
        end do
      end do

      allocate(mpi_send_recv%is_ghost(mesh%n_elems))
      mpi_send_recv%is_ghost = .FALSE.

      !Get local number of elems to send and receive
      call qsort_int(mesh%n_elems, id_elems(:, 1), id_elems(:, 2))

      mesh%n_interior_elems = mesh%n_elems
      do k = 1, mpi_send_recv%n_mpi_recv_neigh
        !Sort elems by global number before renumbering
        call qsort_key(mpi_send_recv%mpi_recv_neigh(k)%n_elems, &
          mpi_send_recv%mpi_recv_neigh(k)%elem_id)
        do l = 1, mpi_send_recv%mpi_recv_neigh(k)%n_elems
          ! print*, mpi_send_recv%mpi_recv_neigh(k)%elem_id(l)
          if (mpi_send_recv%mpi_recv_neigh(k)%elem_id(l) <= 0) then
            print *, "ERROR !", mpi_send_recv%mpi_recv_neigh(k)%elem_id(l)
            error stop
          end if
          call get_loc_id(mesh%n_elems, id_elems, &
            mpi_send_recv%mpi_recv_neigh(k)%elem_id(l), id_loc_elem)
          mpi_send_recv%mpi_recv_neigh(k)%elem_id(l) = id_elems(id_loc_elem, 2)

          if (.not. mesh%elem(id_elems(id_loc_elem, 2))%is_ghost) then
            mesh%elem(id_elems(id_loc_elem, 2))%is_ghost = .TRUE.
            mpi_send_recv%is_ghost(id_elems(id_loc_elem, 2)) = .TRUE.
            mesh%n_interior_elems = mesh%n_interior_elems - 1
          end if
        end do
      end do

      do k = 1, mpi_send_recv%n_mpi_send_neigh
        !Sort elems by global number before renumbering
        call qsort_key(mpi_send_recv%mpi_send_neigh(k)%n_elems, &
          mpi_send_recv%mpi_send_neigh(k)%elem_id)
        do l = 1, mpi_send_recv%mpi_send_neigh(k)%n_elems
          call get_loc_id(mesh%n_elems, id_elems, &
            mpi_send_recv%mpi_send_neigh(k)%elem_id(l), id_loc_elem)
          mpi_send_recv%mpi_send_neigh(k)%elem_id(l) = id_elems(id_loc_elem, 2)
        end do
      end do
    end if
  end subroutine read_mesh_msh_4

  subroutine read_mesh_msh_2(mesh, meshfile_path, meshfile, n_bc, bc_name)
    implicit none

    type(mesh_type), intent(inout) :: mesh
    character(len=255), intent(in) :: meshfile, meshfile_path

    integer(kind=ENTIER), intent(in) :: n_bc
    character(len=255), dimension(:), intent(in) :: bc_name
    integer(kind=ENTIER), dimension(n_bc) :: bc_tag

    integer(kind=ENTIER) :: funit
    integer(kind=ENTIER) :: i, j, iloc
    integer(kind=ENTIER) :: n_elems, n_faces, n_boundary_faces, n_bc_tmp
    integer(kind=ENTIER) :: dummy, dummy_n_tag, elem_kind, tag
    integer(kind=ENTIER), dimension(10) :: dummy_vert, dummy_tag
    character(len=255) :: text

    !! http://www.manpagez.com/info/gmsh/gmsh-2.2.6/gmsh_63.php
    integer(kind=ENTIER), dimension(31), parameter :: n_vert_kind = (/ &
      2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, &
      27, 18, 14, 1, 8, 20, 15, 13, 9, 10, 12, &
      15, 15, 21, 4, 5, 6, 20, 35, 56/)
    integer(kind=ENTIER), dimension(31), parameter :: n_face_kind = (/ &
      0, 0, 0, 4, 6, 5, 5, 0, 0, 0, 0, &
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      0, 0, 0, 0, 0, 0, 0, 0, 0/)
    integer(kind=ENTIER), dimension(31), parameter :: n_face_vert_kind = (/ &
      0, 0, 0, 4*3, 6*4, 3*4 + 2*3, 1*4 + 4*3, 0, 0, 0, 0, &
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      0, 0, 0, 0, 0, 0, 0, 0, 0/)

    open (newunit=funit, file=trim(adjustl(meshfile_path))//trim(adjustl(meshfile)), &
      status="old")

    if (n_bc /= 0) then
      read (funit, '(a)') text
      do while ("$PhysicalNames" /= text(1:14))
        read (funit, '(a)') text
        if ("$Nodes" == text(1:6)) then
          print *, "No Physical names inside mesh!!!"
          error stop
        end if
      end do

      !!Read PhysicalNames of BC
      read (funit, *) n_bc_tmp
      do i = 1, n_bc_tmp
        read (funit, *) dummy, tag, text
        do j = 1, n_bc
          if (trim(adjustl(bc_name(j))) == trim(adjustl(text))) then
            bc_tag(j) = tag
          end if
        end do
      end do
    end if

    read (funit, '(a)') text
    do while ("$Nodes" /= text(1:6))
      read (funit, '(a)') text
    end do

    read (funit, *) mesh%n_vert
    do i = 1, mesh%n_vert
      read (funit, *) text
    end do

    read (funit, '(a)') text
    do while ("$Elements" /= text(1:9))
      read (funit, '(a)') text
    end do

    read (funit, *) n_elems
    iloc = 1
    n_faces = 0
    n_boundary_faces = 0
    do i = 1, n_elems
      read (funit, *) dummy, elem_kind, &
        dummy_n_tag, dummy_tag(:dummy_n_tag), &
        dummy_vert(:n_vert_kind(elem_kind))

      if ((elem_kind == 2 .or. elem_kind == 3)) then
        n_boundary_faces = n_boundary_faces + 1
      end if

      if ((elem_kind == 4 .or. elem_kind == 5) .or. &
        (elem_kind == 6 .or. elem_kind == 7)) then
        iloc = iloc + 1
        n_faces = n_faces + n_face_kind(elem_kind)
      end if
    end do

    mesh%n_elems = iloc - 1
    mesh%n_interior_elems = mesh%n_elems
    mesh%n_faces = n_faces
    mesh%n_boundary_faces = n_boundary_faces

    if (mesh%n_elems < 1 .or. mesh%n_faces < 1) then
      print *, achar(27)//"[31m[-] No 3D elements in mesh !"//achar(27)//"[0m", mesh%n_elems, mesh%n_faces
      error stop
    end if

    close (funit)

    allocate (mesh%vert(mesh%n_vert))
    allocate (mesh%elem(mesh%n_elems))

    !!Each face is counted twice this is fixed in build face function
    allocate (mesh%minimal_face(mesh%n_faces))
    allocate (mesh%boundary_face(mesh%n_boundary_faces))

    open (newunit=funit, file=trim(adjustl(meshfile_path))//trim(adjustl(meshfile)), &
      status="old")

    read (funit, '(a)') text
    do while ("$Nodes" /= text(1:6))
      read (funit, '(a)') text
    end do

    read (funit, *) mesh%n_vert
    do i = 1, mesh%n_vert
      read (funit, *) dummy, mesh%vert(i)%coord(1), mesh%vert(i)%coord(2), mesh%vert(i)%coord(3)
    end do

    read (funit, '(a)') text
    do while ("$Elements" /= text(1:9))
      read (funit, '(a)') text
    end do

    read (funit, *) n_elems
    iloc = 1
    n_faces = 1
    n_boundary_faces = 1

    do i = 1, n_elems
      read (funit, *) dummy, elem_kind, &
        dummy_n_tag, dummy_tag(:dummy_n_tag), &
        dummy_vert(:n_vert_kind(elem_kind))

      if (elem_kind == 2 .or. elem_kind == 3) then

        mesh%boundary_face(n_boundary_faces)%tag = 0
        do j = 1, n_bc
          if (bc_tag(j) == dummy_tag(1)) then
            mesh%boundary_face(n_boundary_faces)%tag = j
          end if
        end do

        mesh%boundary_face(n_boundary_faces)%n_vert = n_vert_kind(elem_kind)
        allocate (mesh%boundary_face(n_boundary_faces)%vert(n_vert_kind(elem_kind)))
        do j = 1, n_vert_kind(elem_kind)
          mesh%boundary_face(n_boundary_faces)%vert(j) = dummy_vert(j)
        end do
        n_boundary_faces = n_boundary_faces + 1
      end if

      if ((elem_kind == 4 .or. elem_kind == 5) .or. &
        (elem_kind == 6 .or. elem_kind == 7)) then

        mesh%elem(iloc)%elem_kind = elem_kind

        mesh%elem(iloc)%n_vert = n_vert_kind(elem_kind)
        mesh%elem(iloc)%n_faces = n_face_kind(elem_kind)
        allocate (mesh%elem(iloc)%vert(n_vert_kind(elem_kind)))

        do j = 1, n_vert_kind(elem_kind)
          mesh%elem(iloc)%vert(j) = dummy_vert(j)
        end do

        do j = 1, n_face_kind(elem_kind)
          mesh%minimal_face(n_faces + j - 1)%left_neigh = iloc
          mesh%minimal_face(n_faces + j - 1)%left_neigh_face = j
        end do

        allocate (mesh%elem(iloc)%face(n_face_kind(elem_kind)))
        if (elem_kind == 4) then !!tet
          mesh%minimal_face(n_faces)%n_vert = 3
          allocate (mesh%minimal_face(n_faces)%vert(3))
          mesh%elem(iloc)%face(1) = n_faces
          mesh%minimal_face(n_faces)%vert(1) = mesh%elem(iloc)%vert(1)
          mesh%minimal_face(n_faces)%vert(2) = mesh%elem(iloc)%vert(3)
          mesh%minimal_face(n_faces)%vert(3) = mesh%elem(iloc)%vert(2)

          mesh%minimal_face(n_faces + 1)%n_vert = 3
          allocate (mesh%minimal_face(n_faces + 1)%vert(3))
          mesh%elem(iloc)%face(2) = n_faces + 1
          mesh%minimal_face(n_faces + 1)%vert(1) = mesh%elem(iloc)%vert(1)
          mesh%minimal_face(n_faces + 1)%vert(2) = mesh%elem(iloc)%vert(2)
          mesh%minimal_face(n_faces + 1)%vert(3) = mesh%elem(iloc)%vert(4)

          mesh%minimal_face(n_faces + 2)%n_vert = 3
          allocate (mesh%minimal_face(n_faces + 2)%vert(3))
          mesh%elem(iloc)%face(3) = n_faces + 2
          mesh%minimal_face(n_faces + 2)%vert(1) = mesh%elem(iloc)%vert(1)
          mesh%minimal_face(n_faces + 2)%vert(2) = mesh%elem(iloc)%vert(4)
          mesh%minimal_face(n_faces + 2)%vert(3) = mesh%elem(iloc)%vert(3)

          mesh%minimal_face(n_faces + 3)%n_vert = 3
          allocate (mesh%minimal_face(n_faces + 3)%vert(3))
          mesh%elem(iloc)%face(4) = n_faces + 3
          mesh%minimal_face(n_faces + 3)%vert(1) = mesh%elem(iloc)%vert(2)
          mesh%minimal_face(n_faces + 3)%vert(2) = mesh%elem(iloc)%vert(3)
          mesh%minimal_face(n_faces + 3)%vert(3) = mesh%elem(iloc)%vert(4)

        else if (elem_kind == 5) then !! Hexa
          mesh%minimal_face(n_faces)%n_vert = 4
          allocate (mesh%minimal_face(n_faces)%vert(4))
          mesh%elem(iloc)%face(1) = n_faces
          mesh%minimal_face(n_faces)%vert(1) = mesh%elem(iloc)%vert(1)
          mesh%minimal_face(n_faces)%vert(2) = mesh%elem(iloc)%vert(4)
          mesh%minimal_face(n_faces)%vert(3) = mesh%elem(iloc)%vert(3)
          mesh%minimal_face(n_faces)%vert(4) = mesh%elem(iloc)%vert(2)

          mesh%minimal_face(n_faces + 1)%n_vert = 4
          allocate (mesh%minimal_face(n_faces + 1)%vert(4))
          mesh%elem(iloc)%face(2) = n_faces + 1
          mesh%minimal_face(n_faces + 1)%vert(1) = mesh%elem(iloc)%vert(3)
          mesh%minimal_face(n_faces + 1)%vert(2) = mesh%elem(iloc)%vert(4)
          mesh%minimal_face(n_faces + 1)%vert(3) = mesh%elem(iloc)%vert(8)
          mesh%minimal_face(n_faces + 1)%vert(4) = mesh%elem(iloc)%vert(7)

          mesh%minimal_face(n_faces + 2)%n_vert = 4
          allocate (mesh%minimal_face(n_faces + 2)%vert(4))
          mesh%elem(iloc)%face(3) = n_faces + 2
          mesh%minimal_face(n_faces + 2)%vert(1) = mesh%elem(iloc)%vert(7)
          mesh%minimal_face(n_faces + 2)%vert(2) = mesh%elem(iloc)%vert(8)
          mesh%minimal_face(n_faces + 2)%vert(3) = mesh%elem(iloc)%vert(5)
          mesh%minimal_face(n_faces + 2)%vert(4) = mesh%elem(iloc)%vert(6)

          mesh%minimal_face(n_faces + 3)%n_vert = 4
          allocate (mesh%minimal_face(n_faces + 3)%vert(4))
          mesh%elem(iloc)%face(4) = n_faces + 3
          mesh%minimal_face(n_faces + 3)%vert(1) = mesh%elem(iloc)%vert(1)
          mesh%minimal_face(n_faces + 3)%vert(2) = mesh%elem(iloc)%vert(2)
          mesh%minimal_face(n_faces + 3)%vert(3) = mesh%elem(iloc)%vert(6)
          mesh%minimal_face(n_faces + 3)%vert(4) = mesh%elem(iloc)%vert(5)

          mesh%minimal_face(n_faces + 4)%n_vert = 4
          allocate (mesh%minimal_face(n_faces + 4)%vert(4))
          mesh%elem(iloc)%face(5) = n_faces + 4
          mesh%minimal_face(n_faces + 4)%vert(1) = mesh%elem(iloc)%vert(1)
          mesh%minimal_face(n_faces + 4)%vert(2) = mesh%elem(iloc)%vert(5)
          mesh%minimal_face(n_faces + 4)%vert(3) = mesh%elem(iloc)%vert(8)
          mesh%minimal_face(n_faces + 4)%vert(4) = mesh%elem(iloc)%vert(4)

          mesh%minimal_face(n_faces + 5)%n_vert = 4
          allocate (mesh%minimal_face(n_faces + 5)%vert(4))
          mesh%elem(iloc)%face(6) = n_faces + 5
          mesh%minimal_face(n_faces + 5)%vert(1) = mesh%elem(iloc)%vert(2)
          mesh%minimal_face(n_faces + 5)%vert(2) = mesh%elem(iloc)%vert(3)
          mesh%minimal_face(n_faces + 5)%vert(3) = mesh%elem(iloc)%vert(7)
          mesh%minimal_face(n_faces + 5)%vert(4) = mesh%elem(iloc)%vert(6)

        else if (elem_kind == 6) then !! Prism
          mesh%minimal_face(n_faces)%n_vert = 3
          allocate (mesh%minimal_face(n_faces)%vert(3))
          mesh%elem(iloc)%face(1) = n_faces
          mesh%minimal_face(n_faces)%vert(1) = mesh%elem(iloc)%vert(1)
          mesh%minimal_face(n_faces)%vert(2) = mesh%elem(iloc)%vert(3)
          mesh%minimal_face(n_faces)%vert(3) = mesh%elem(iloc)%vert(2)

          mesh%minimal_face(n_faces + 1)%n_vert = 3
          allocate (mesh%minimal_face(n_faces + 1)%vert(3))
          mesh%elem(iloc)%face(2) = n_faces + 1
          mesh%minimal_face(n_faces + 1)%vert(1) = mesh%elem(iloc)%vert(4)
          mesh%minimal_face(n_faces + 1)%vert(2) = mesh%elem(iloc)%vert(5)
          mesh%minimal_face(n_faces + 1)%vert(3) = mesh%elem(iloc)%vert(6)

          mesh%minimal_face(n_faces + 2)%n_vert = 4
          allocate (mesh%minimal_face(n_faces + 2)%vert(4))
          mesh%elem(iloc)%face(3) = n_faces + 2
          mesh%minimal_face(n_faces + 2)%vert(1) = mesh%elem(iloc)%vert(1)
          mesh%minimal_face(n_faces + 2)%vert(2) = mesh%elem(iloc)%vert(2)
          mesh%minimal_face(n_faces + 2)%vert(3) = mesh%elem(iloc)%vert(5)
          mesh%minimal_face(n_faces + 2)%vert(4) = mesh%elem(iloc)%vert(4)

          mesh%minimal_face(n_faces + 3)%n_vert = 4
          allocate (mesh%minimal_face(n_faces + 3)%vert(4))
          mesh%elem(iloc)%face(4) = n_faces + 3
          mesh%minimal_face(n_faces + 3)%vert(1) = mesh%elem(iloc)%vert(2)
          mesh%minimal_face(n_faces + 3)%vert(2) = mesh%elem(iloc)%vert(3)
          mesh%minimal_face(n_faces + 3)%vert(3) = mesh%elem(iloc)%vert(6)
          mesh%minimal_face(n_faces + 3)%vert(4) = mesh%elem(iloc)%vert(5)

          mesh%minimal_face(n_faces + 4)%n_vert = 4
          allocate (mesh%minimal_face(n_faces + 4)%vert(4))
          mesh%elem(iloc)%face(5) = n_faces + 4
          mesh%minimal_face(n_faces + 4)%vert(1) = mesh%elem(iloc)%vert(1)
          mesh%minimal_face(n_faces + 4)%vert(2) = mesh%elem(iloc)%vert(4)
          mesh%minimal_face(n_faces + 4)%vert(3) = mesh%elem(iloc)%vert(6)
          mesh%minimal_face(n_faces + 4)%vert(4) = mesh%elem(iloc)%vert(3)

        else if (elem_kind == 7) then !! Pyramid
          mesh%minimal_face(n_faces)%n_vert = 4
          allocate (mesh%minimal_face(n_faces)%vert(4))
          mesh%elem(iloc)%face(1) = n_faces
          mesh%minimal_face(n_faces)%vert(1) = mesh%elem(iloc)%vert(1)
          mesh%minimal_face(n_faces)%vert(2) = mesh%elem(iloc)%vert(4)
          mesh%minimal_face(n_faces)%vert(3) = mesh%elem(iloc)%vert(3)
          mesh%minimal_face(n_faces)%vert(4) = mesh%elem(iloc)%vert(2)

          mesh%minimal_face(n_faces + 1)%n_vert = 3
          allocate (mesh%minimal_face(n_faces + 1)%vert(3))
          mesh%elem(iloc)%face(2) = n_faces + 1
          mesh%minimal_face(n_faces + 1)%vert(1) = mesh%elem(iloc)%vert(1)
          mesh%minimal_face(n_faces + 1)%vert(2) = mesh%elem(iloc)%vert(2)
          mesh%minimal_face(n_faces + 1)%vert(3) = mesh%elem(iloc)%vert(5)

          mesh%minimal_face(n_faces + 2)%n_vert = 3
          allocate (mesh%minimal_face(n_faces + 2)%vert(3))
          mesh%elem(iloc)%face(3) = n_faces + 2
          mesh%minimal_face(n_faces + 2)%vert(1) = mesh%elem(iloc)%vert(2)
          mesh%minimal_face(n_faces + 2)%vert(2) = mesh%elem(iloc)%vert(3)
          mesh%minimal_face(n_faces + 2)%vert(3) = mesh%elem(iloc)%vert(5)

          mesh%minimal_face(n_faces + 3)%n_vert = 3
          allocate (mesh%minimal_face(n_faces + 3)%vert(3))
          mesh%elem(iloc)%face(4) = n_faces + 3
          mesh%minimal_face(n_faces + 3)%vert(1) = mesh%elem(iloc)%vert(3)
          mesh%minimal_face(n_faces + 3)%vert(2) = mesh%elem(iloc)%vert(4)
          mesh%minimal_face(n_faces + 3)%vert(3) = mesh%elem(iloc)%vert(5)

          mesh%minimal_face(n_faces + 4)%n_vert = 3
          allocate (mesh%minimal_face(n_faces + 4)%vert(3))
          mesh%elem(iloc)%face(5) = n_faces + 4
          mesh%minimal_face(n_faces + 4)%vert(1) = mesh%elem(iloc)%vert(1)
          mesh%minimal_face(n_faces + 4)%vert(2) = mesh%elem(iloc)%vert(5)
          mesh%minimal_face(n_faces + 4)%vert(3) = mesh%elem(iloc)%vert(4)
        end if

        iloc = iloc + 1
        n_faces = n_faces + n_face_kind(elem_kind)
      end if
    end do

    close (funit)
  end subroutine read_mesh_msh_2

  subroutine get_loc_id(n, vect, kglob, kloc)
    implicit none

    !vect should be sorted ascending
    integer(kind=ENTIER), intent(in) :: n, kglob
    integer(kind=ENTIER), intent(inout) :: kloc
    integer(kind=ENTIER), dimension(n), intent(in) :: vect

    integer(kind=ENTIER) :: lmin, lmax, l

    lmin = 1
    lmax = n

    do while (lmax > lmin + 1)
      l = (lmax + lmin - 1)/2 + 1
      if (vect(l) == kglob) then
        kloc = l
        return
      else if (vect(lmin) == kglob) then
        kloc = lmin
        return
      else if (vect(lmax) == kglob) then
        kloc = lmax
        return
      else if (vect(l) < kglob) then
        lmin = l
      else if (vect(l) > kglob) then
        lmax = l
      end if
    end do

    if (vect(lmin) == kglob) then
      kloc = lmin
    else if (vect(lmax) == kglob) then
      kloc = lmax
    else
      kloc = -1
    end if
  end subroutine get_loc_id

  subroutine rescale_mesh(mesh, sf)
    implicit none

    type(mesh_type), intent(inout) :: mesh
    real(kind=DOUBLE), intent(in) :: sf

    integer(kind=ENTIER) :: i

    do i=1, mesh%n_vert
      mesh%vert(i)%coord = sf*mesh%vert(i)%coord
    end do
  end subroutine rescale_mesh

end module mesh_reading_module
