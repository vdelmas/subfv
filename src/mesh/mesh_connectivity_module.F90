module mesh_connectivity_module
  use precision_module
  use mesh_module
  use mpi_module
  implicit none

  private

  public :: build_mesh
contains

  subroutine build_mesh(mesh, num_procs, mpi_send_recv, use_sub_entities, b2d)
    use mpi_module
    implicit none

    type(mesh_type), intent(inout) :: mesh
    integer, intent(in) :: num_procs
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv
    logical, intent(in) :: use_sub_entities, b2d

    call build_node_elem_neigh(mesh)
    call renumber_vert_mesh(mesh)
    call renumber_elem_mesh(mesh, num_procs, mpi_send_recv)
    call build_faces(mesh)
    if (num_procs > 1) call mpi_exchange_fix_boundary(mesh, mpi_send_recv)
    call build_node_face_connectivity(mesh)
    call renumber_faces_mesh(mesh)
    if (use_sub_entities) call build_sub_faces(mesh)
    if (use_sub_entities) call build_sub_elems(mesh)
    if (use_sub_entities) call find_sub_elems_around_sub_face(mesh)
    if (use_sub_entities) call build_node_sub_face_connectivity(mesh)
    call compute_elem_tag(mesh)
    if (use_sub_entities) call compute_n_max_sub_faces_around_node(mesh)
    call compute_neigh_by_vert(mesh)
    if (num_procs > 1) call compute_ghost_vert(mesh)
    call find_if_vert_is_bound(mesh, b2d)
  end subroutine build_mesh

  subroutine build_faces(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    !!Temporary container used while detecting duplicates in build_faces
    type face_container_type
      type(minimal_face_type), dimension(:), allocatable :: face
    end type face_container_type

    integer(kind=ENTIER) :: i, j, k, iel, ier, ifl, ifr
    integer(kind=ENTIER) :: id_glob_max, id_min, n_faces
    integer(kind=ENTIER), dimension(:), allocatable :: n_collision
    type(face_container_type), dimension(:), allocatable :: face_container
    logical :: already_contained

    n_faces = 0

    id_glob_max = 0
    do j = 1, mesh%n_faces
      id_min = mesh%minimal_face(j)%vert(1)
      do k = 2, mesh%minimal_face(j)%n_vert
        id_min = min(id_min, mesh%minimal_face(j)%vert(k))
      end do
      id_glob_max = max(id_glob_max, id_min)
    end do

    allocate (n_collision(id_glob_max))
    allocate (face_container(id_glob_max))

    n_collision(:) = 0
    do j = 1, mesh%n_faces
      id_min = mesh%minimal_face(j)%vert(1)
      do k = 2, mesh%minimal_face(j)%n_vert
        id_min = min(id_min, mesh%minimal_face(j)%vert(k))
      end do
      n_collision(id_min) = n_collision(id_min) + 1
    end do

    do j = 1, id_glob_max
      allocate (face_container(j)%face(n_collision(j)))
      do k = 1, n_collision(j)
        face_container(j)%face(k)%n_vert = -1
      end do
    end do

    do j = 1, mesh%n_faces
      id_min = mesh%minimal_face(j)%vert(1)
      do k = 2, mesh%minimal_face(j)%n_vert
        id_min = min(id_min, mesh%minimal_face(j)%vert(k))
      end do

      already_contained = .FALSE.
      do k = 1, n_collision(id_min)
        if (face_container(id_min)%face(k)%n_vert > 0) then
          if (same_faces(mesh%minimal_face(j), face_container(id_min)%face(k))) then
            already_contained = .TRUE.

            face_container(id_min)%face(k)%right_neigh = mesh%minimal_face(j)%left_neigh
            face_container(id_min)%face(k)%right_neigh_face = mesh%minimal_face(j)%left_neigh_face
            exit
          end if
        end if
      end do

      if (.not. already_contained) then
        do k = 1, n_collision(id_min)
          if (face_container(id_min)%face(k)%n_vert <= 0) then
            n_faces = n_faces + 1
            face_container(id_min)%face(k)%n_vert = mesh%minimal_face(j)%n_vert
            allocate (face_container(id_min)%face(k)%vert(mesh%minimal_face(j)%n_vert))
            do i = 1, mesh%minimal_face(j)%n_vert
              face_container(id_min)%face(k)%vert(i) = mesh%minimal_face(j)%vert(i)
            end do
            face_container(id_min)%face(k)%left_neigh = mesh%minimal_face(j)%left_neigh
            face_container(id_min)%face(k)%left_neigh_face = mesh%minimal_face(j)%left_neigh_face
            exit
          end if
        end do
      end if
    end do

    !!Check for boundary faces match
    do j = 1, mesh%n_boundary_faces
      id_min = mesh%boundary_face(j)%vert(1)
      do k = 2, mesh%boundary_face(j)%n_vert
        id_min = min(id_min, mesh%boundary_face(j)%vert(k))
      end do

      already_contained = .FALSE.
      do k = 1, n_collision(id_min)
        if (face_container(id_min)%face(k)%n_vert > 0) then
          if (same_faces_boundary(mesh%boundary_face(j), face_container(id_min)%face(k))) then
            already_contained = .TRUE.
            if (face_container(id_min)%face(k)%right_neigh <= 0) then
              face_container(id_min)%face(k)%right_neigh = -mesh%boundary_face(j)%tag
              face_container(id_min)%face(k)%right_neigh_face = 0
            end if
            exit
          end if
        end if
      end do
    end do

    do i = 1, mesh%n_faces
      deallocate (mesh%minimal_face(i)%vert)
    end do
    deallocate (mesh%minimal_face)

    mesh%n_faces = n_faces
    allocate (mesh%face(mesh%n_faces))

    n_faces = 0
    do j = 1, id_glob_max
      do k = 1, n_collision(j)
        if (face_container(j)%face(k)%n_vert > 0) then
          n_faces = n_faces + 1

          mesh%face(n_faces)%n_vert = face_container(j)%face(k)%n_vert
          allocate (mesh%face(n_faces)%vert(mesh%face(n_faces)%n_vert))
          do i = 1, mesh%face(n_faces)%n_vert
            mesh%face(n_faces)%vert(i) = face_container(j)%face(k)%vert(i)
          end do
          mesh%face(n_faces)%left_neigh = face_container(j)%face(k)%left_neigh
          mesh%face(n_faces)%left_neigh_face = face_container(j)%face(k)%left_neigh_face

          mesh%face(n_faces)%right_neigh = face_container(j)%face(k)%right_neigh
          mesh%face(n_faces)%right_neigh_face = face_container(j)%face(k)%right_neigh_face

          mesh%face(n_faces)%area = 0.0_DOUBLE
          mesh%face(n_faces)%norm = 0.0_DOUBLE

          iel = mesh%face(n_faces)%left_neigh
          ifl = mesh%face(n_faces)%left_neigh_face

          ier = mesh%face(n_faces)%right_neigh
          ifr = mesh%face(n_faces)%right_neigh_face

          mesh%elem(iel)%face(ifl) = n_faces

          if (ier > 0) then
            mesh%elem(ier)%face(ifr) = n_faces
          end if

          deallocate (face_container(j)%face(k)%vert)
        end if
      end do
    end do

    deallocate (face_container)
  end subroutine build_faces

  subroutine mpi_exchange_fix_boundary(mesh, info_mpi)
    use mpi
    use mpi_module
    implicit none

    type(mesh_type), intent(inout) :: mesh
    type(mpi_send_recv_type), intent(inout) :: info_mpi

    integer :: mpi_ierr
    integer(kind=ENTIER) :: i, k, j, iloc
    integer(kind=ENTIER) :: id_elem, id_face

    type right_neigh_type
      integer(kind=ENTIER), dimension(:), allocatable :: rn
    end type right_neigh_type

    type(right_neigh_type), dimension(:), allocatable :: send_right_neigh
    type(right_neigh_type), dimension(:), allocatable :: recv_right_neigh

    allocate (send_right_neigh(info_mpi%n_mpi_send_neigh))
    do i = 1, info_mpi%n_mpi_send_neigh
      iloc = 1
      do k = 1, info_mpi%mpi_send_neigh(i)%n_elems
        id_elem = info_mpi%mpi_send_neigh(i)%elem_id(k)
        iloc = iloc + mesh%elem(id_elem)%n_faces
      end do
      allocate (send_right_neigh(i)%rn(iloc - 1))

      iloc = 1
      do k = 1, info_mpi%mpi_send_neigh(i)%n_elems
        id_elem = info_mpi%mpi_send_neigh(i)%elem_id(k)
        do j = 1, mesh%elem(id_elem)%n_faces
          id_face = mesh%elem(id_elem)%face(j)
          if (mesh%face(id_face)%right_neigh < 0) then
            send_right_neigh(i)%rn(iloc) = mesh%face(id_face)%right_neigh
          else
            send_right_neigh(i)%rn(iloc) = 42
          end if
          iloc = iloc + 1
        end do
      end do

      call mpi_isend(send_right_neigh(i)%rn(1), &
        size(send_right_neigh(i)%rn), MPI_INT, &
        info_mpi%mpi_send_neigh(i)%partition_id, &
        0, MPI_COMM_WORLD, info_mpi%mpi_reqsend(i), mpi_ierr)
    end do

    allocate (recv_right_neigh(info_mpi%n_mpi_recv_neigh))
    do i = 1, info_mpi%n_mpi_recv_neigh
      iloc = 1
      do k = 1, info_mpi%mpi_recv_neigh(i)%n_elems
        id_elem = info_mpi%mpi_recv_neigh(i)%elem_id(k)
        iloc = iloc + mesh%elem(id_elem)%n_faces
      end do
      allocate (recv_right_neigh(i)%rn(iloc - 1))
    end do

    do i = 1, info_mpi%n_mpi_recv_neigh
      call mpi_irecv(recv_right_neigh(i)%rn(1), &
        size(recv_right_neigh(i)%rn), MPI_INT, &
        info_mpi%mpi_recv_neigh(i)%partition_id, &
        MPI_ANY_TAG, MPI_COMM_WORLD, info_mpi%mpi_reqrecv(i), mpi_ierr)
    end do

    call mpi_waitall(info_mpi%n_mpi_send_neigh, &
      info_mpi%mpi_reqsend, info_mpi%mpi_sendstat, mpi_ierr)
    call mpi_waitall(info_mpi%n_mpi_recv_neigh, &
      info_mpi%mpi_reqrecv, info_mpi%mpi_recvstat, mpi_ierr)

    do i = 1, info_mpi%n_mpi_recv_neigh
      iloc = 1
      do k = 1, info_mpi%mpi_recv_neigh(i)%n_elems
        id_elem = info_mpi%mpi_recv_neigh(i)%elem_id(k)
        do j = 1, mesh%elem(id_elem)%n_faces
          id_face = mesh%elem(id_elem)%face(j)
          ! print*, id_face, mesh%face(id_face)%right_neigh, recv_right_neigh(i)%rn(iloc)
          if (recv_right_neigh(i)%rn(iloc) < 0) then
            mesh%face(id_face)%right_neigh = recv_right_neigh(i)%rn(iloc)
          end if
          iloc = iloc + 1
        end do
      end do
    end do
  end subroutine mpi_exchange_fix_boundary

  pure subroutine build_node_face_connectivity(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    integer(kind=ENTIER) :: i, id_vert, j, k

    do i = 1, mesh%n_vert
      mesh%vert(i)%n_faces_neigh = 0
    end do

    do i = 1, mesh%n_faces
      do j = 1, mesh%face(i)%n_vert
        id_vert = mesh%face(i)%vert(j)
        mesh%vert(id_vert)%n_faces_neigh = mesh%vert(id_vert)%n_faces_neigh + 1
      end do
    end do

    do i = 1, mesh%n_vert
      allocate (mesh%vert(i)%face_neigh(mesh%vert(i)%n_faces_neigh))
      mesh%vert(i)%face_neigh(:) = -1
    end do

    do i = 1, mesh%n_faces
      do j = 1, mesh%face(i)%n_vert
        id_vert = mesh%face(i)%vert(j)
        do k = 1, mesh%vert(id_vert)%n_faces_neigh
          if (mesh%vert(id_vert)%face_neigh(k) < 0) then
            mesh%vert(id_vert)%face_neigh(k) = i
            exit
          end if
        end do
      end do
    end do
  end subroutine build_node_face_connectivity

  pure subroutine build_sub_faces(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    integer(kind=ENTIER) :: i, j, k
    integer(kind=ENTIER) :: id_face, id_sub_face, id_left, id_right

    mesh%n_sub_faces = 0
    do i = 1, mesh%n_faces
      mesh%face(i)%n_sub_faces = mesh%face(i)%n_vert
      allocate (mesh%face(i)%sub_face(mesh%face(i)%n_sub_faces))
      mesh%face(i)%sub_face(:) = -1
      mesh%n_sub_faces = mesh%n_sub_faces + mesh%face(i)%n_vert
    end do
    allocate (mesh%sub_face(mesh%n_sub_faces))

    do i = 1, mesh%n_elems
      mesh%elem(i)%n_sub_faces = 0
      do j = 1, mesh%elem(i)%n_faces
        id_face = mesh%elem(i)%face(j)
        mesh%elem(i)%n_sub_faces = mesh%elem(i)%n_sub_faces + mesh%face(id_face)%n_vert
      end do
      allocate (mesh%elem(i)%sub_face(mesh%elem(i)%n_sub_faces))
      mesh%elem(i)%sub_face(:) = -1
    end do

    id_sub_face = 1
    do i = 1, mesh%n_vert
      do j = 1, mesh%vert(i)%n_faces_neigh

        id_face = mesh%vert(i)%face_neigh(j)
        ! if (id_face <= 0) error stop

        do k = 1, mesh%face(id_face)%n_sub_faces
          if (mesh%face(id_face)%sub_face(k) < 0) then
            mesh%face(id_face)%sub_face(k) = id_sub_face
            exit
          end if
        end do

        mesh%sub_face(id_sub_face)%mesh_face = id_face
        mesh%sub_face(id_sub_face)%mesh_vert = i
        mesh%sub_face(id_sub_face)%left_elem_neigh = mesh%face(id_face)%left_neigh
        mesh%sub_face(id_sub_face)%right_elem_neigh = mesh%face(id_face)%right_neigh

        id_left = mesh%face(id_face)%left_neigh
        do k = 1, mesh%elem(id_left)%n_sub_faces
          if (mesh%elem(id_left)%sub_face(k) < 0) then
            mesh%elem(id_left)%sub_face(k) = id_sub_face
            exit
          end if
        end do

        id_right = mesh%face(id_face)%right_neigh
        if (id_right > 0) then
          do k = 1, mesh%elem(id_right)%n_sub_faces
            if (mesh%elem(id_right)%sub_face(k) < 0) then
              mesh%elem(id_right)%sub_face(k) = id_sub_face
              exit
            end if
          end do
        end if

        id_sub_face = id_sub_face + 1
      end do
    end do
  end subroutine build_sub_faces

  pure subroutine build_sub_elems(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    integer(kind=ENTIER) :: i, j, k
    integer(kind=ENTIER) :: id_sub_face, id_sub_elem, id_sub_face_loc

    mesh%n_sub_elems = 0
    do i = 1, mesh%n_elems
      mesh%elem(i)%n_sub_elems = mesh%elem(i)%n_vert
      allocate (mesh%elem(i)%sub_elem(mesh%elem(i)%n_sub_elems))
      mesh%n_sub_elems = mesh%n_sub_elems + mesh%elem(i)%n_sub_elems
    end do
    allocate (mesh%sub_elem(mesh%n_sub_elems))

    id_sub_elem = 1
    do i = 1, mesh%n_elems
      do j = 1, mesh%elem(i)%n_vert
        mesh%sub_elem(id_sub_elem)%n_sub_faces = 0
        do k = 1, mesh%elem(i)%n_sub_faces
          id_sub_face = mesh%elem(i)%sub_face(k)
          if (mesh%sub_face(id_sub_face)%mesh_vert == mesh%elem(i)%vert(j)) then
            mesh%sub_elem(id_sub_elem)%n_sub_faces = mesh%sub_elem(id_sub_elem)%n_sub_faces + 1
          end if
        end do
        allocate (mesh%sub_elem(id_sub_elem)%sub_face(mesh%sub_elem(id_sub_elem)%n_sub_faces))
        id_sub_elem = id_sub_elem + 1
      end do
    end do

    id_sub_elem = 1
    do i = 1, mesh%n_elems
      do j = 1, mesh%elem(i)%n_vert
        id_sub_face_loc = 1
        mesh%elem(i)%sub_elem(j) = id_sub_elem
        do k = 1, mesh%elem(i)%n_sub_faces
          id_sub_face = mesh%elem(i)%sub_face(k)
          if (mesh%sub_face(id_sub_face)%mesh_vert == mesh%elem(i)%vert(j)) then
            mesh%sub_elem(id_sub_elem)%mesh_elem = i
            mesh%sub_elem(id_sub_elem)%mesh_vert = mesh%sub_face(id_sub_face)%mesh_vert
            mesh%sub_elem(id_sub_elem)%sub_face(id_sub_face_loc) = id_sub_face
            id_sub_face_loc = id_sub_face_loc + 1
          end if
        end do
        id_sub_elem = id_sub_elem + 1
      end do
    end do
  end subroutine build_sub_elems

  pure subroutine find_sub_elems_around_sub_face(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    integer(kind=ENTIER) :: i, j, id_sub_face, id_elem

    do i = 1, mesh%n_sub_elems
      id_elem = mesh%sub_elem(i)%mesh_elem
      do j = 1, mesh%sub_elem(i)%n_sub_faces
        id_sub_face = mesh%sub_elem(i)%sub_face(j)
        if (mesh%sub_face(id_sub_face)%left_elem_neigh == id_elem) then
          mesh%sub_face(id_sub_face)%left_sub_elem_neigh = i
        end if
        if (mesh%sub_face(id_sub_face)%right_elem_neigh == id_elem) then
          mesh%sub_face(id_sub_face)%right_sub_elem_neigh = i
        end if
      end do
    end do
  end subroutine find_sub_elems_around_sub_face

  pure subroutine build_node_elem_neigh(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    integer(kind=ENTIER) :: i, id_vert, j

    do i = 1, mesh%n_vert
      mesh%vert(i)%n_elems_neigh = 0
    end do

    do i = 1, mesh%n_elems
      do j = 1, mesh%elem(i)%n_vert
        id_vert = mesh%elem(i)%vert(j)
        mesh%vert(id_vert)%n_elems_neigh = &
          mesh%vert(id_vert)%n_elems_neigh + 1
      end do
    end do

    do i = 1, mesh%n_vert
      allocate (mesh%vert(i)%elem_neigh(mesh%vert(i)%n_elems_neigh))
    end do

    do i = 1, mesh%n_vert
      mesh%vert(i)%n_elems_neigh = 0
    end do

    do i = 1, mesh%n_elems
      do j = 1, mesh%elem(i)%n_vert
        id_vert = mesh%elem(i)%vert(j)
        mesh%vert(id_vert)%elem_neigh(mesh%vert(id_vert)%n_elems_neigh + 1) = i
        mesh%vert(id_vert)%n_elems_neigh = mesh%vert(id_vert)%n_elems_neigh + 1
      end do
    end do
  end subroutine build_node_elem_neigh

  pure subroutine build_node_sub_face_connectivity(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    integer(kind=ENTIER) :: i, id_vert

    do i = 1, mesh%n_vert
      mesh%vert(i)%n_sub_faces_neigh = 0
      mesh%vert(i)%n_sub_elems_neigh = 0
    end do

    do i = 1, mesh%n_sub_faces
      id_vert = mesh%sub_face(i)%mesh_vert
      mesh%vert(id_vert)%n_sub_faces_neigh = mesh%vert(id_vert)%n_sub_faces_neigh + 1
    end do

    do i = 1, mesh%n_sub_elems
      id_vert = mesh%sub_elem(i)%mesh_vert
      mesh%vert(id_vert)%n_sub_elems_neigh = mesh%vert(id_vert)%n_sub_elems_neigh + 1
    end do

    do i = 1, mesh%n_vert
      allocate (mesh%vert(i)%sub_face_neigh(mesh%vert(i)%n_sub_faces_neigh))
      allocate (mesh%vert(i)%sub_elem_neigh(mesh%vert(i)%n_sub_elems_neigh))
    end do

    do i = 1, mesh%n_vert
      mesh%vert(i)%n_sub_faces_neigh = 0
      mesh%vert(i)%n_sub_elems_neigh = 0
    end do

    do i = 1, mesh%n_sub_faces
      id_vert = mesh%sub_face(i)%mesh_vert
      mesh%sub_face(i)%id_loc_around_node = mesh%vert(id_vert)%n_sub_faces_neigh + 1
      mesh%vert(id_vert)%sub_face_neigh(mesh%vert(id_vert)%n_sub_faces_neigh + 1) = i
      mesh%vert(id_vert)%n_sub_faces_neigh = mesh%vert(id_vert)%n_sub_faces_neigh + 1
    end do

    do i = 1, mesh%n_sub_elems
      id_vert = mesh%sub_elem(i)%mesh_vert
      mesh%sub_elem(i)%id_loc_around_node = mesh%vert(id_vert)%n_sub_elems_neigh + 1
      mesh%vert(id_vert)%sub_elem_neigh(mesh%vert(id_vert)%n_sub_elems_neigh + 1) = i
      mesh%vert(id_vert)%n_sub_elems_neigh = mesh%vert(id_vert)%n_sub_elems_neigh + 1
    end do
  end subroutine build_node_sub_face_connectivity

  pure subroutine compute_elem_tag(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    integer(kind=ENTIER) :: i, tag, id_elem

    do i = 1, mesh%n_faces
      if (mesh%face(i)%right_neigh < 0) then
        tag = -mesh%face(i)%right_neigh
        id_elem = mesh%face(i)%left_neigh
        if (mesh%elem(id_elem)%tag < 0) then
          mesh%elem(id_elem)%tag = tag
        else
          !Smallest tag number wins
          if (mesh%elem(id_elem)%tag > tag) then
            mesh%elem(id_elem)%tag = tag
          end if
        end if
      end if
    end do
  end subroutine compute_elem_tag

  pure subroutine compute_n_max_sub_faces_around_node(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    integer(kind=ENTIER) :: i

    mesh%n_max_sub_faces_around_node = 0
    do i = 1, mesh%n_vert
      mesh%n_max_sub_faces_around_node = &
        max(mesh%n_max_sub_faces_around_node, mesh%vert(i)%n_sub_faces_neigh)
    end do
  end subroutine compute_n_max_sub_faces_around_node

  subroutine compute_neigh_by_vert(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    integer(kind=ENTIER) :: i, j, k, iloc, id_vert, n_neigh_by_vert, id_face
    integer(kind=ENTIER), dimension(:), allocatable :: temp_neigh_by_vert

    real(kind=DOUBLE) :: rn
    logical :: already_added

    do i = 1, mesh%n_elems
      n_neigh_by_vert = 0
      do j = 1, mesh%elem(i)%n_vert
        id_vert = mesh%elem(i)%vert(j)
        do k = 1, mesh%vert(id_vert)%n_elems_neigh
          if (mesh%vert(id_vert)%elem_neigh(k) /= i) then
            n_neigh_by_vert = n_neigh_by_vert + 1
          end if
        end do
      end do

      allocate (temp_neigh_by_vert(n_neigh_by_vert))

      iloc = 1
      do j = 1, mesh%elem(i)%n_vert
        id_vert = mesh%elem(i)%vert(j)
        do k = 1, mesh%vert(id_vert)%n_elems_neigh
          if (mesh%vert(id_vert)%elem_neigh(k) /= i) then
            temp_neigh_by_vert(iloc) = mesh%vert(id_vert)%elem_neigh(k)
            iloc = iloc + 1
          end if
        end do
      end do

      !Fix duplicate elements
      mesh%elem(i)%n_neigh_by_vert = n_neigh_by_vert
      do j = 1, n_neigh_by_vert
        do k = j + 1, n_neigh_by_vert
          if (temp_neigh_by_vert(k) > 0 .and. &
            temp_neigh_by_vert(j) == temp_neigh_by_vert(k)) then
            temp_neigh_by_vert(k) = -1
            mesh%elem(i)%n_neigh_by_vert = mesh%elem(i)%n_neigh_by_vert - 1
          end if
        end do
      end do

      allocate (mesh%elem(i)%neigh_by_vert(mesh%elem(i)%n_neigh_by_vert))
      iloc = 1
      do j = 1, n_neigh_by_vert
        if (temp_neigh_by_vert(j) > 0) then
          mesh%elem(i)%neigh_by_vert(iloc) = temp_neigh_by_vert(j)
          iloc = iloc + 1
        end if
      end do

      deallocate (temp_neigh_by_vert)
    end do

  end subroutine compute_neigh_by_vert

  subroutine compute_ghost_vert(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    integer(kind=ENTIER) :: i, j, id_elem
    logical :: all_ghost

    do i = 1, mesh%n_vert
      all_ghost = .TRUE.
      do j = 1, mesh%vert(i)%n_elems_neigh
        id_elem = mesh%vert(i)%elem_neigh(j)
        if (.not. mesh%elem(id_elem)%is_ghost) all_ghost = .FALSE.
      end do

      if (all_ghost) then
        mesh%vert(i)%is_ghost = .TRUE.
      else
        mesh%vert(i)%is_ghost = .FALSE.
      end if
    end do
  end subroutine compute_ghost_vert

  pure function same_faces(a, b)
    implicit none

    type(minimal_face_type), intent(in) :: a, b
    logical :: same_faces

    logical :: vert_found
    integer(kind=ENTIER) :: i, j

    if (a%n_vert /= b%n_vert) then
      same_faces = .FALSE.
    else
      same_faces = .TRUE.
      do i = 1, a%n_vert
        vert_found = .FALSE.
        do j = 1, b%n_vert
          if (a%vert(i) == b%vert(j)) vert_found = .TRUE.
        end do
        same_faces = vert_found .AND. same_faces
      end do
    end if
  end function same_faces

  pure function same_faces_boundary(a, b)
    implicit none

    type(boundary_face_type), intent(in) :: a
    type(minimal_face_type), intent(in) :: b
    logical :: same_faces_boundary

    logical :: vert_found
    integer(kind=ENTIER) :: i, j

    if (a%n_vert /= b%n_vert) then
      same_faces_boundary = .FALSE.
    else
      same_faces_boundary = .TRUE.
      do i = 1, a%n_vert
        vert_found = .FALSE.
        do j = 1, b%n_vert
          if (a%vert(i) == b%vert(j)) vert_found = .TRUE.
        end do
        same_faces_boundary = vert_found .AND. same_faces_boundary
      end do
    end if
  end function same_faces_boundary

  subroutine find_if_vert_is_bound(mesh, b2d)
    implicit none

    type(mesh_type), intent(inout) :: mesh
    logical, intent(in) :: b2d

    integer(kind=DOUBLE) :: i, j, id_face

    do i = 1, mesh%n_vert
      mesh%vert(i)%is_bound = .FALSE.
      do j = 1, mesh%vert(i)%n_faces_neigh
        id_face = mesh%vert(i)%face_neigh(j)
        if( b2d ) then
          if (mesh%face(id_face)%right_neigh <= 0 .and. &
            abs(mesh%face(id_face)%norm(3)) < 1e-12_DOUBLE ) then
            mesh%vert(i)%is_bound = .TRUE.
          end if
        else
          if (mesh%face(id_face)%right_neigh <= 0) then
            mesh%vert(i)%is_bound = .TRUE.
          end if
        end if
      end do
    end do
  end subroutine find_if_vert_is_bound

  subroutine renumber_elem_mesh(mesh, num_procs, mpi_send_recv)
    implicit none

    type(mesh_type), intent(inout) :: mesh
    integer(kind=ENTIER), intent(in) :: num_procs
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    integer(kind=ENTIER) :: i, iloc, j, &
      id_face, id_sub_face, le, re, n_ghost_stuff
    integer(kind=ENTIER), dimension(mesh%n_elems) :: otn

    type(elem_type), dimension(:), allocatable :: elem_tmp

    otn = -1 

    !Gen olt to new numbering of elements
    if( num_procs > 1 ) then
      otn = -1
      do i=1, mpi_send_recv%n_mpi_send_neigh
        do j=1, mpi_send_recv%mpi_send_neigh(i)%n_elems
          otn( mpi_send_recv%mpi_send_neigh(i)%elem_id(j) ) = 0
        end do
      end do
      do i=1, mpi_send_recv%n_mpi_recv_neigh
        do j=1, mpi_send_recv%mpi_recv_neigh(i)%n_elems
          otn(mpi_send_recv%mpi_recv_neigh(i)%elem_id(j)) = 0
        end do
      end do

      n_ghost_stuff = 0
      do i=1, mesh%n_elems
        if( otn(i) == 0 ) n_ghost_stuff = n_ghost_stuff + 1
      end do

      otn = -1
      iloc = mesh%n_elems - n_ghost_stuff + 1
      do i=1, mpi_send_recv%n_mpi_send_neigh
        do j=1, mpi_send_recv%mpi_send_neigh(i)%n_elems
          if( otn( mpi_send_recv%mpi_send_neigh(i)%elem_id(j) ) < 0 ) then
            otn( mpi_send_recv%mpi_send_neigh(i)%elem_id(j) ) = iloc
            iloc = iloc + 1
          end if
        end do
      end do
      do i=1, mpi_send_recv%n_mpi_recv_neigh
        do j=1, mpi_send_recv%mpi_recv_neigh(i)%n_elems
          if( otn(mpi_send_recv%mpi_recv_neigh(i)%elem_id(j)) < 0 ) then
            otn(mpi_send_recv%mpi_recv_neigh(i)%elem_id(j)) = iloc
            iloc = iloc + 1
          end if
        end do
      end do
    end if

    iloc = 1
    do i=1, mesh%n_vert
      do j=1, mesh%vert(i)%n_elems_neigh
        le = mesh%vert(i)%elem_neigh(j)
        if( otn(le) <= 0 ) then
          otn(le) = iloc
          iloc = iloc + 1
        end if
      end do
    end do

    allocate(elem_tmp(mesh%n_elems))
    elem_tmp = mesh%elem

    do i=1, mesh%n_elems
      mesh%elem(otn(i)) = elem_tmp(i)
    end do

    !Fix datastructure
    do i=1, mesh%n_faces
      mesh%minimal_face(i)%left_neigh = otn(mesh%minimal_face(i)%left_neigh)
    end do
    do i=1, mesh%n_vert
      do j=1, mesh%vert(i)%n_elems_neigh
        mesh%vert(i)%elem_neigh(j) = otn(mesh%vert(i)%elem_neigh(j))
      end do
    end do

    !Fix mpi stuff
    if( num_procs > 1 ) then
      do i=1, mpi_send_recv%n_mpi_send_neigh
        do j=1, mpi_send_recv%mpi_send_neigh(i)%n_elems
          mpi_send_recv%mpi_send_neigh(i)%elem_id(j) = &
            otn(mpi_send_recv%mpi_send_neigh(i)%elem_id(j))
        end do
      end do
      do i=1, mpi_send_recv%n_mpi_recv_neigh
        do j=1, mpi_send_recv%mpi_recv_neigh(i)%n_elems
          mpi_send_recv%mpi_recv_neigh(i)%elem_id(j) = &
            otn(mpi_send_recv%mpi_recv_neigh(i)%elem_id(j))
        end do
      end do
    end if
  end subroutine renumber_elem_mesh

  subroutine renumber_vert_mesh(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    type(mesh_type) :: mesh_tmp
    integer(kind=ENTIER) :: i, iloc, iloc_before, j, &
      id_face, id_sub_face, le, re, id_vert, &
      id_vert_neigh, id_elem, k, iloc_start, iloc_end
    integer(kind=ENTIER), dimension(mesh%n_vert) :: otn
    integer(kind=ENTIER), dimension(mesh%n_vert) :: nto

    type(vert_type), dimension(:), allocatable :: vert_tmp

    otn = -1
    nto = -1
    otn(1) = 1
    nto(1) = 1
    iloc = 2

    iloc_start = iloc-1
    iloc_end = iloc-1
    do while ( iloc <= mesh%n_vert )
      iloc_before = iloc
      do i=iloc_start, iloc_end
        id_vert = nto(i)
        do j=1, mesh%vert(id_vert)%n_elems_neigh
          id_elem = mesh%vert(id_vert)%elem_neigh(j)
          do k=1, mesh%elem(id_elem)%n_vert
            id_vert_neigh = mesh%elem(id_elem)%vert(k)
            if( otn(id_vert_neigh) < 0 ) then
              otn(id_vert_neigh) = iloc
              nto(iloc) = id_vert_neigh
              iloc = iloc + 1
            end if
          end do
        end do
      end do
      if( iloc == iloc_before ) then
        do i = 1, mesh%n_vert
          if (otn(i) < 0) then
            otn(i) = iloc
            nto(iloc) = i
            iloc = iloc + 1
            exit
          end if
        end do
      end if
      iloc_start = iloc_end+1
      iloc_end = iloc-1
    end do

    allocate(vert_tmp(mesh%n_vert))
    vert_tmp = mesh%vert

    do i=1, mesh%n_vert
      mesh%vert(otn(i)) = vert_tmp(i)
    end do

    !Fix datastructure
    do i=1, mesh%n_elems
      do j=1, mesh%elem(i)%n_vert
        mesh%elem(i)%vert(j) = otn(mesh%elem(i)%vert(j))
      end do
    end do
    do i=1, mesh%n_faces
      do j=1, mesh%minimal_face(i)%n_vert
        mesh%minimal_face(i)%vert(j) = otn(mesh%minimal_face(i)%vert(j))
      end do
    end do
    do i=1, mesh%n_boundary_faces
      do j=1, mesh%boundary_face(i)%n_vert
        mesh%boundary_face(i)%vert(j) = otn(mesh%boundary_face(i)%vert(j))
      end do
    end do
  end subroutine renumber_vert_mesh

  subroutine renumber_faces_mesh(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    type(mesh_type) :: mesh_tmp
    integer(kind=ENTIER) :: i, iloc, j, &
      id_face, id_sub_face, le, re, id_vert, &
      id_vert_neigh, id_elem, k, iloc_start, iloc_end
    integer(kind=ENTIER), dimension(mesh%n_faces) :: otn
    integer(kind=ENTIER), dimension(mesh%n_faces) :: nto

    type(face_type), dimension(:), allocatable :: faces_tmp

    otn = -1
    nto = -1
    iloc = 1
    do i=1, mesh%n_vert
      do j=1, mesh%vert(i)%n_faces_neigh
        id_face = mesh%vert(i)%face_neigh(j)
        if( otn(id_face) < 0 ) then
          otn(id_face) = iloc
          nto(iloc) = id_face
          iloc = iloc + 1
        end if
      end do
    end do

    allocate(faces_tmp(mesh%n_faces))
    faces_tmp = mesh%face

    do i=1, mesh%n_faces
      mesh%face(otn(i)) = faces_tmp(i)
    end do

    !Fix datastructure
    do i=1, mesh%n_elems
      do j=1, mesh%elem(i)%n_faces
        mesh%elem(i)%face(j) = otn(mesh%elem(i)%face(j))
      end do
    end do
    do i=1, mesh%n_vert
      do j=1, mesh%vert(i)%n_faces_neigh
        mesh%vert(i)%face_neigh(j) = otn(mesh%vert(i)%face_neigh(j))
      end do
    end do
  end subroutine renumber_faces_mesh
end module mesh_connectivity_module
