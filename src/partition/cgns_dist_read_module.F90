! Block-wise CGNS reading: each rank opens the file read-only and reads only its
! ranges of cells, faces and nodes. Single base / single unstructured zone only.
module cgns_dist_read_module
  use mpi
  use precision_module
  use cgns, only: CGSIZE_T, cgenum_t, CG_MODE_READ, CG_NODE_NOT_FOUND, &
    cg_open_f, cg_close_f, cg_nbases_f, cg_base_read_f, cg_nzones_f, &
    cg_zone_type_f, cg_zone_read_f, cg_nsections_f, cg_section_read_f, &
    cg_elementpartialsize_f, cg_ncoords_f, cg_coord_info_f, cg_nbocos_f, &
    cg_boco_info_f, cg_boco_gridlocation_read_f, cg_goto_f, cg_famname_read_f, &
    cg_nfamilies_f, cg_family_read_f, cg_nfamily_names_f, cg_family_name_read_f, &
    cg_get_error_f, &
    Unstructured, NODE, BAR_2, TRI_3, QUAD_4, TETRA_4, PYRA_5, PENTA_6, HEXA_8, &
    MIXED, NGON_n, NFACE_n, LongInteger, RealDouble, &
    PointRange, PointList, ElementRange, ElementList, &
    Vertex, CellCenter, FaceCenter
  implicit none

  private
  public :: cgns_mesh_file_type, KIND_POLYHEDRON
  public :: cgns_open_mesh, cgns_close_mesh
  public :: cgns_read_cell_block, cgns_read_face_block, cgns_read_coord_block
  public :: cgns_face_space_index, cgns_read_bc_faces

  ! elem_kind stored for NFACE_n polyhedra (VTK_POLYHEDRON)
  integer(kind=ENTIER), parameter :: KIND_POLYHEDRON = 42

  type :: cgns_section_type
    integer :: id
    integer(cgenum_t) :: etype
    integer(kind=ENTIER_D) :: istart, iend
    logical :: in_vol_space = .false., in_face_space = .false.
    ! offset of this section's first element in the volume/face index spaces
    integer(kind=ENTIER_D) :: vol_off = 0, face_off = 0
  end type cgns_section_type

  type :: cgns_mesh_file_type
    integer :: fn = -1, base = 1, zone = 1
    integer(kind=ENTIER_D) :: n_nodes = 0
    integer(kind=ENTIER_D) :: n_vol_space = 0, n_face_space = 0
    integer :: n_sections = 0
    type(cgns_section_type), dimension(:), allocatable :: sec
    character(len=32), dimension(3) :: coord_name
  end type cgns_mesh_file_type

contains

  subroutine cg_check(ier, what)
    integer, intent(in) :: ier
    character(len=*), intent(in) :: what

    character(len=256) :: msg
    integer :: mpi_ierr

    if (ier /= 0) then
      msg = ""
      call cg_get_error_f(msg)
      print *, achar(27)//"[31m[-] CGNS error in "//trim(what)//": "//trim(msg)//achar(27)//"[0m"
      call MPI_ABORT(MPI_COMM_WORLD, 1, mpi_ierr)
    end if
  end subroutine cg_check

  subroutine fatal(msg)
    character(len=*), intent(in) :: msg
    integer :: mpi_ierr

    print *, achar(27)//"[31m[-] "//trim(msg)//achar(27)//"[0m"
    call MPI_ABORT(MPI_COMM_WORLD, 1, mpi_ierr)
  end subroutine fatal

  ! Number of nodes of a linear CGNS element type, 0 for types we do not
  ! accept as cells/faces (higher-order, BAR, NODE, ...).
  pure integer function npe_of(etype)
    integer(cgenum_t), intent(in) :: etype

    select case (etype)
    case (TRI_3)
      npe_of = 3
    case (QUAD_4, TETRA_4)
      npe_of = 4
    case (PYRA_5)
      npe_of = 5
    case (PENTA_6)
      npe_of = 6
    case (HEXA_8)
      npe_of = 8
    case default
      npe_of = 0
    end select
  end function npe_of

  ! CGNS 3D element type -> gmsh elem_kind used throughout subfv.
  pure integer function gmsh_kind_of(etype)
    integer(cgenum_t), intent(in) :: etype

    select case (etype)
    case (TETRA_4)
      gmsh_kind_of = 4
    case (HEXA_8)
      gmsh_kind_of = 5
    case (PENTA_6)
      gmsh_kind_of = 6
    case (PYRA_5)
      gmsh_kind_of = 7
    case default
      gmsh_kind_of = 0
    end select
  end function gmsh_kind_of

  pure logical function is_3d_type(etype)
    integer(cgenum_t), intent(in) :: etype
    is_3d_type = (etype == TETRA_4 .or. etype == PYRA_5 .or. etype == PENTA_6 &
      .or. etype == HEXA_8)
  end function is_3d_type

  pure logical function is_2d_type(etype)
    integer(cgenum_t), intent(in) :: etype
    is_2d_type = (etype == TRI_3 .or. etype == QUAD_4)
  end function is_2d_type

  ! NODE/BAR elements are ignored; any other unsupported type is an error.
  pure logical function is_ignored_type(etype)
    integer(cgenum_t), intent(in) :: etype
    is_ignored_type = (etype == NODE .or. etype == BAR_2)
  end function is_ignored_type

  subroutine cgns_open_mesh(filename, cf)
    character(len=*), intent(in) :: filename
    type(cgns_mesh_file_type), intent(out) :: cf

    integer :: ier, nbases, nzones, cell_dim, phys_dim, ncoords, s, i, j
    integer :: nbndry, parent_flag
    integer(cgenum_t) :: ztype, etype, dtype
    integer(CGSIZE_T), dimension(9) :: zsize
    integer(CGSIZE_T) :: istart, iend
    character(len=32) :: bname, zname, sname
    type(cgns_section_type) :: tmp

    call cg_open_f(trim(filename), CG_MODE_READ, cf%fn, ier)
    call cg_check(ier, "cg_open "//trim(filename))

    call cg_nbases_f(cf%fn, nbases, ier)
    call cg_check(ier, "cg_nbases")
    if (nbases /= 1) call fatal("CGNS: exactly one base expected, found more (multi-base not supported)")

    call cg_base_read_f(cf%fn, cf%base, bname, cell_dim, phys_dim, ier)
    call cg_check(ier, "cg_base_read")
    if (cell_dim /= 3 .or. phys_dim /= 3) &
      call fatal("CGNS: only 3D meshes (cell_dim = phys_dim = 3) are supported")

    call cg_nzones_f(cf%fn, cf%base, nzones, ier)
    call cg_check(ier, "cg_nzones")
    if (nzones /= 1) call fatal("CGNS: exactly one zone expected (multi-zone meshes not supported)")

    call cg_zone_type_f(cf%fn, cf%base, cf%zone, ztype, ier)
    call cg_check(ier, "cg_zone_type")
    if (ztype /= Unstructured) call fatal("CGNS: only unstructured zones are supported")

    call cg_zone_read_f(cf%fn, cf%base, cf%zone, zname, zsize, ier)
    call cg_check(ier, "cg_zone_read")
    cf%n_nodes = int(zsize(1), ENTIER_D)
    if (cf%n_nodes > int(huge(1_ENTIER), ENTIER_D)) &
      call fatal("CGNS: more nodes than a 32-bit ENTIER can index")

    call cg_ncoords_f(cf%fn, cf%base, cf%zone, ncoords, ier)
    call cg_check(ier, "cg_ncoords")
    if (ncoords /= 3) call fatal("CGNS: 3 coordinate arrays expected")
    do i = 1, 3
      call cg_coord_info_f(cf%fn, cf%base, cf%zone, i, dtype, cf%coord_name(i), ier)
      call cg_check(ier, "cg_coord_info")
    end do

    call cg_nsections_f(cf%fn, cf%base, cf%zone, cf%n_sections, ier)
    call cg_check(ier, "cg_nsections")
    allocate (cf%sec(cf%n_sections))

    do s = 1, cf%n_sections
      call cg_section_read_f(cf%fn, cf%base, cf%zone, s, sname, etype, istart, iend, &
        nbndry, parent_flag, ier)
      call cg_check(ier, "cg_section_read")
      cf%sec(s)%id = s
      cf%sec(s)%etype = etype
      cf%sec(s)%istart = int(istart, ENTIER_D)
      cf%sec(s)%iend = int(iend, ENTIER_D)
      if (cf%sec(s)%iend > int(huge(1_ENTIER), ENTIER_D)) &
        call fatal("CGNS: element ids beyond 32-bit ENTIER range")

      if (is_3d_type(etype) .or. etype == NFACE_n) then
        cf%sec(s)%in_vol_space = .true.
      else if (is_2d_type(etype) .or. etype == NGON_n) then
        cf%sec(s)%in_face_space = .true.
      else if (etype == MIXED) then
        cf%sec(s)%in_vol_space = .true.
        cf%sec(s)%in_face_space = .true.
      else if (.not. is_ignored_type(etype)) then
        call fatal("CGNS: section '"//trim(sname)//"' has an unsupported element type "// &
          "(only linear TRI/QUAD/TETRA/PYRA/PENTA/HEXA, MIXED, NGON_n, NFACE_n)")
      end if
    end do

    ! Sort sections by first element id so that both index spaces follow
    ! ascending element ids (insertion sort, a handful of sections).
    do s = 2, cf%n_sections
      tmp = cf%sec(s)
      j = s - 1
      do while (j >= 1)
        if (cf%sec(j)%istart <= tmp%istart) exit
        cf%sec(j + 1) = cf%sec(j)
        j = j - 1
      end do
      cf%sec(j + 1) = tmp
    end do

    cf%n_vol_space = 0
    cf%n_face_space = 0
    do s = 1, cf%n_sections
      if (cf%sec(s)%in_vol_space) then
        cf%sec(s)%vol_off = cf%n_vol_space
        cf%n_vol_space = cf%n_vol_space + cf%sec(s)%iend - cf%sec(s)%istart + 1
      end if
      if (cf%sec(s)%in_face_space) then
        cf%sec(s)%face_off = cf%n_face_space
        cf%n_face_space = cf%n_face_space + cf%sec(s)%iend - cf%sec(s)%istart + 1
      end if
    end do

    if (cf%n_vol_space == 0) call fatal("CGNS: no 3D element section in mesh")
  end subroutine cgns_open_mesh

  subroutine cgns_close_mesh(cf)
    type(cgns_mesh_file_type), intent(inout) :: cf
    integer :: ier

    call cg_close_f(cf%fn, ier)
    call cg_check(ier, "cg_close")
    cf%fn = -1
  end subroutine cgns_close_mesh

  ! Reads raw poly/mixed connectivity of elements [e1,e2] of section s:
  ! conn(:) and offsets(1:n+1), offsets rebased to start at 0.
  subroutine read_poly_range(cf, s, e1, e2, conn, offsets)
    type(cgns_mesh_file_type), intent(in) :: cf
    integer, intent(in) :: s
    integer(kind=ENTIER_D), intent(in) :: e1, e2
    integer(kind=ENTIER_D), dimension(:), allocatable, intent(out) :: conn, offsets

    integer(CGSIZE_T) :: r1, r2, dsize
    integer :: ier

    r1 = int(e1, CGSIZE_T)
    r2 = int(e2, CGSIZE_T)
    call cg_elementpartialsize_f(cf%fn, cf%base, cf%zone, cf%sec(s)%id, r1, r2, dsize, ier)
    call cg_check(ier, "cg_ElementPartialSize")
    allocate (conn(max(dsize, 1_CGSIZE_T)), offsets(e2 - e1 + 2))
    call cg_poly_elements_general_read_f(cf%fn, cf%base, cf%zone, cf%sec(s)%id, r1, r2, &
      LongInteger, conn, offsets, ier)
    call cg_check(ier, "cg_poly_elements_general_read")
    offsets = offsets - offsets(1)
  end subroutine read_poly_range

  subroutine read_fixed_range(cf, s, e1, e2, npe, conn)
    type(cgns_mesh_file_type), intent(in) :: cf
    integer, intent(in) :: s, npe
    integer(kind=ENTIER_D), intent(in) :: e1, e2
    integer(kind=ENTIER_D), dimension(:), allocatable, intent(out) :: conn

    integer(CGSIZE_T) :: r1, r2
    integer :: ier

    r1 = int(e1, CGSIZE_T)
    r2 = int(e2, CGSIZE_T)
    allocate (conn(npe*(e2 - e1 + 1)))
    call cg_elements_general_read_f(cf%fn, cf%base, cf%zone, cf%sec(s)%id, r1, r2, &
      LongInteger, conn, ier)
    call cg_check(ier, "cg_elements_general_read")
  end subroutine read_fixed_range

  ! Cells with volume-space indices first..last, as CSR: node ids, or signed
  ! NGON ids for KIND_POLYHEDRON. 2D elements inside MIXED sections are skipped.
  subroutine cgns_read_cell_block(cf, first, last, n, gid, kind, ptr, dat)
    type(cgns_mesh_file_type), intent(in) :: cf
    integer(kind=ENTIER_D), intent(in) :: first, last
    integer(kind=ENTIER), intent(out) :: n
    integer(kind=ENTIER), dimension(:), allocatable, intent(out) :: gid, kind, ptr, dat

    integer :: s
    integer(kind=ENTIER_D) :: a, b, e1, e2, j, len_sec, k, np
    integer(kind=ENTIER_D), dimension(:), allocatable :: conn, offsets
    integer(kind=ENTIER) :: npe, nmax, ndat
    integer(cgenum_t) :: et

    ! First pass sizes are bounded by the raw connectivity read; grow simply.
    nmax = int(max(last - first + 1, 0_ENTIER_D), ENTIER)
    allocate (gid(nmax), kind(nmax), ptr(nmax + 1))
    allocate (dat(max(8*nmax, 1)))
    n = 0
    ndat = 0
    ptr(1) = 1

    do s = 1, cf%n_sections
      if (.not. cf%sec(s)%in_vol_space) cycle
      len_sec = cf%sec(s)%iend - cf%sec(s)%istart + 1
      a = max(first, cf%sec(s)%vol_off + 1)
      b = min(last, cf%sec(s)%vol_off + len_sec)
      if (a > b) cycle
      e1 = cf%sec(s)%istart + (a - cf%sec(s)%vol_off - 1)
      e2 = cf%sec(s)%istart + (b - cf%sec(s)%vol_off - 1)
      et = cf%sec(s)%etype

      if (is_3d_type(et)) then
        npe = npe_of(et)
        call read_fixed_range(cf, s, e1, e2, npe, conn)
        do j = 0, e2 - e1
          call push_cell(int(e1 + j, ENTIER), gmsh_kind_of(et), &
            int(conn(j*npe + 1:j*npe + npe), ENTIER))
        end do
      else if (et == NFACE_n) then
        call read_poly_range(cf, s, e1, e2, conn, offsets)
        do j = 0, e2 - e1
          call push_cell(int(e1 + j, ENTIER), KIND_POLYHEDRON, &
            int(conn(offsets(j + 1) + 1:offsets(j + 2)), ENTIER))
        end do
      else if (et == MIXED) then
        call read_poly_range(cf, s, e1, e2, conn, offsets)
        do j = 0, e2 - e1
          k = offsets(j + 1) + 1
          et = int(conn(k), cgenum_t)
          np = offsets(j + 2) - k
          if (is_3d_type(et)) then
            if (np /= npe_of(et)) call fatal("CGNS: inconsistent MIXED element size")
            call push_cell(int(e1 + j, ENTIER), gmsh_kind_of(et), &
              int(conn(k + 1:k + np), ENTIER))
          else if (.not. (is_2d_type(et) .or. is_ignored_type(et))) then
            call fatal("CGNS: unsupported element type inside a MIXED section")
          end if
        end do
      end if
    end do

  contains

    subroutine push_cell(g, kd, list)
      integer(kind=ENTIER), intent(in) :: g, kd
      integer(kind=ENTIER), dimension(:), intent(in) :: list

      integer(kind=ENTIER), dimension(:), allocatable :: tmp

      if (ndat + size(list) > size(dat)) then
        allocate (tmp(2*size(dat) + size(list)))
        tmp(1:ndat) = dat(1:ndat)
        call move_alloc(tmp, dat)
      end if
      n = n + 1
      gid(n) = g
      kind(n) = kd
      dat(ndat + 1:ndat + size(list)) = list
      ndat = ndat + size(list)
      ptr(n + 1) = ndat + 1
    end subroutine push_cell
  end subroutine cgns_read_cell_block

  ! Face-space index of a CGNS element id, 0 if the id is not in a section
  ! that can hold faces.
  pure function cgns_face_space_index(cf, elem_id) result(idx)
    type(cgns_mesh_file_type), intent(in) :: cf
    integer(kind=ENTIER_D), intent(in) :: elem_id
    integer(kind=ENTIER_D) :: idx

    integer :: s

    idx = 0
    do s = 1, cf%n_sections
      if (.not. cf%sec(s)%in_face_space) cycle
      if (elem_id >= cf%sec(s)%istart .and. elem_id <= cf%sec(s)%iend) then
        idx = cf%sec(s)%face_off + (elem_id - cf%sec(s)%istart) + 1
        return
      end if
    end do
  end function cgns_face_space_index

  ! Faces with face-space indices first..last, as CSR node lists (empty for a
  ! 3D element inside a MIXED section).
  subroutine cgns_read_face_block(cf, first, last, ptr, dat)
    type(cgns_mesh_file_type), intent(in) :: cf
    integer(kind=ENTIER_D), intent(in) :: first, last
    integer(kind=ENTIER), dimension(:), allocatable, intent(out) :: ptr, dat

    integer :: s
    integer(kind=ENTIER_D) :: a, b, e1, e2, j, len_sec, k, np
    integer(kind=ENTIER_D), dimension(:), allocatable :: conn, offsets
    integer(kind=ENTIER) :: npe, n, nloc, ndat
    integer(cgenum_t) :: et

    nloc = int(max(last - first + 1, 0_ENTIER_D), ENTIER)
    allocate (ptr(nloc + 1))
    allocate (dat(max(4*nloc, 1)))
    n = 0
    ndat = 0
    ptr(1) = 1

    do s = 1, cf%n_sections
      if (.not. cf%sec(s)%in_face_space) cycle
      len_sec = cf%sec(s)%iend - cf%sec(s)%istart + 1
      a = max(first, cf%sec(s)%face_off + 1)
      b = min(last, cf%sec(s)%face_off + len_sec)
      if (a > b) cycle
      e1 = cf%sec(s)%istart + (a - cf%sec(s)%face_off - 1)
      e2 = cf%sec(s)%istart + (b - cf%sec(s)%face_off - 1)
      et = cf%sec(s)%etype

      if (is_2d_type(et)) then
        npe = npe_of(et)
        call read_fixed_range(cf, s, e1, e2, npe, conn)
        do j = 0, e2 - e1
          call push_face(int(conn(j*npe + 1:j*npe + npe), ENTIER))
        end do
      else if (et == NGON_n) then
        call read_poly_range(cf, s, e1, e2, conn, offsets)
        do j = 0, e2 - e1
          call push_face(int(conn(offsets(j + 1) + 1:offsets(j + 2)), ENTIER))
        end do
      else if (et == MIXED) then
        call read_poly_range(cf, s, e1, e2, conn, offsets)
        do j = 0, e2 - e1
          k = offsets(j + 1) + 1
          et = int(conn(k), cgenum_t)
          np = offsets(j + 2) - k
          if (is_2d_type(et)) then
            call push_face(int(conn(k + 1:k + np), ENTIER))
          else
            call push_face([integer(kind=ENTIER) ::])
          end if
        end do
      end if
    end do

    if (n /= nloc) call fatal("CGNS: face block size mismatch (internal error)")

  contains

    subroutine push_face(list)
      integer(kind=ENTIER), dimension(:), intent(in) :: list

      integer(kind=ENTIER), dimension(:), allocatable :: tmp

      if (ndat + size(list) > size(dat)) then
        allocate (tmp(2*size(dat) + size(list)))
        tmp(1:ndat) = dat(1:ndat)
        call move_alloc(tmp, dat)
      end if
      n = n + 1
      dat(ndat + 1:ndat + size(list)) = list
      ndat = ndat + size(list)
      ptr(n + 1) = ndat + 1
    end subroutine push_face
  end subroutine cgns_read_face_block

  subroutine cgns_read_coord_block(cf, first, last, xyz)
    type(cgns_mesh_file_type), intent(in) :: cf
    integer(kind=ENTIER_D), intent(in) :: first, last
    real(kind=DOUBLE), dimension(:, :), allocatable, intent(out) :: xyz

    real(kind=DOUBLE), dimension(:), allocatable :: buf
    integer(CGSIZE_T) :: r1, r2
    integer :: ier, i, n

    n = int(max(last - first + 1, 0_ENTIER_D))
    allocate (xyz(3, n), buf(max(n, 1)))
    if (n == 0) return
    r1 = int(first, CGSIZE_T)
    r2 = int(last, CGSIZE_T)
    do i = 1, 3
      call cg_coord_read_f(cf%fn, cf%base, cf%zone, trim(cf%coord_name(i)), RealDouble, &
        r1, r2, buf, ier)
      call cg_check(ier, "cg_coord_read "//trim(cf%coord_name(i)))
      xyz(i, :) = buf(1:n)
    end do
  end subroutine cgns_read_coord_block

  ! This rank's share of (face element id, bc index); a BC is matched on its name,
  ! its FamilyName or that family's FamilyName_t children (gmsh). Index 0 = unlisted.
  subroutine cgns_read_bc_faces(cf, n_bc, bc_name, me, num_procs, n_pairs, face_id, bc_idx)
    type(cgns_mesh_file_type), intent(in) :: cf
    integer(kind=ENTIER), intent(in) :: n_bc, me, num_procs
    character(len=255), dimension(:), intent(in) :: bc_name
    integer(kind=ENTIER), intent(out) :: n_pairs
    integer(kind=ENTIER), dimension(:), allocatable, intent(out) :: face_id, bc_idx

    integer :: ier, nbocos, bc, nfam, f, nnames, nb, ng, i, k, kmatch, ndataset
    integer, dimension(3) :: normal_index
    integer(cgenum_t) :: bctype, ptset, loc, ndtype
    integer(CGSIZE_T) :: npnts, nlistsize
    integer(CGSIZE_T), dimension(:), allocatable :: pnts
    real(kind=DOUBLE), dimension(:), allocatable :: normals
    character(len=32) :: boconame, fname, nname, nfamily
    character(len=128) :: famname
    character(len=32), dimension(:), allocatable :: cand
    integer :: ncand
    integer(kind=ENTIER_D) :: n_ids, i0, i1, j, g
    integer(kind=ENTIER), dimension(:), allocatable :: tmp1, tmp2
    logical, dimension(:), allocatable :: bc_used
    integer :: nvol, nchk

    call cg_nbocos_f(cf%fn, cf%base, cf%zone, nbocos, ier)
    call cg_check(ier, "cg_nbocos")
    call cg_nfamilies_f(cf%fn, cf%base, nfam, ier)
    call cg_check(ier, "cg_nfamilies")

    allocate (face_id(16), bc_idx(16))
    allocate (bc_used(max(n_bc, 1)))
    bc_used = .false.
    n_pairs = 0

    do bc = 1, nbocos
      call cg_boco_info_f(cf%fn, cf%base, cf%zone, bc, boconame, bctype, ptset, npnts, &
        normal_index, nlistsize, ndtype, ndataset, ier)
      call cg_check(ier, "cg_boco_info")
      call cg_boco_gridlocation_read_f(cf%fn, cf%base, cf%zone, bc, loc, ier)
      call cg_check(ier, "cg_boco_gridlocation_read")

      ! Candidate names
      allocate (cand(64))
      ncand = 1
      cand(1) = boconame
      famname = ""
      call cg_goto_f(cf%fn, cf%base, ier, 'Zone_t', cf%zone, 'ZoneBC_t', 1, 'BC_t', bc, 'end')
      call cg_check(ier, "cg_goto BC")
      call cg_famname_read_f(famname, ier)
      if (ier == 0) then
        ncand = ncand + 1
        cand(ncand) = famname(1:32)
        do f = 1, nfam
          call cg_family_read_f(cf%fn, cf%base, f, fname, nb, ng, ier)
          call cg_check(ier, "cg_family_read")
          if (trim(fname) /= trim(famname)) cycle
          call cg_nfamily_names_f(cf%fn, cf%base, f, nnames, ier)
          call cg_check(ier, "cg_nfamily_names")
          do i = 1, nnames
            call cg_family_name_read_f(cf%fn, cf%base, f, i, nname, nfamily, ier)
            call cg_check(ier, "cg_family_name_read")
            if (ncand + 2 <= size(cand)) then
              cand(ncand + 1) = nname
              cand(ncand + 2) = nfamily
              ncand = ncand + 2
            end if
          end do
        end do
      else if (ier /= CG_NODE_NOT_FOUND) then
        call cg_check(ier, "cg_famname_read")
      end if

      kmatch = 0
      do i = 1, ncand
        do k = 1, n_bc
          if (trim(adjustl(bc_name(k))) == trim(adjustl(cand(i)))) then
            if (kmatch /= 0 .and. kmatch /= k) then
              call fatal("CGNS: BC '"//trim(boconame)//"' matches two different bc_name entries")
            end if
            kmatch = k
          end if
        end do
      end do
      if (kmatch > 0) bc_used(kmatch) = .true.

      deallocate (cand)

      if (loc == Vertex) then
        call fatal("CGNS: BC '"//trim(boconame)//"' is given at Vertex location; "// &
          "only face-element BCs (FaceCenter/CellCenter on 2D elements) are supported")
      else if (loc /= FaceCenter .and. loc /= CellCenter) then
        call fatal("CGNS: BC '"//trim(boconame)//"' has an unsupported GridLocation")
      end if

      ! Point set. Ranges are split without reading anything; lists are read
      ! whole (boundary-sized) then this rank keeps its slice.
      allocate (normals(max(int(nlistsize), 1)))
      if (ptset == PointRange .or. ptset == ElementRange) then
        allocate (pnts(2))
        call cg_boco_read_f(cf%fn, cf%base, cf%zone, bc, pnts, normals, ier)
        call cg_check(ier, "cg_boco_read")
        n_ids = int(pnts(2) - pnts(1) + 1, ENTIER_D)
      else if (ptset == PointList .or. ptset == ElementList) then
        allocate (pnts(max(npnts, 1_CGSIZE_T)))
        call cg_boco_read_f(cf%fn, cf%base, cf%zone, bc, pnts, normals, ier)
        call cg_check(ier, "cg_boco_read")
        n_ids = int(npnts, ENTIER_D)
      else
        call fatal("CGNS: BC '"//trim(boconame)//"' has an unsupported point set type")
      end if
      deallocate (normals)

      ! gmsh writes a BC_t for volume groups too: skip it, refuse mixed ones.
      if (n_ids > 0) then
        if (ptset == PointRange .or. ptset == ElementRange) then
          nvol = count_vol(int(pnts(1), ENTIER_D)) + count_vol(int(pnts(2), ENTIER_D))
          nchk = 2
        else
          nvol = 0
          do j = 1, n_ids
            nvol = nvol + count_vol(int(pnts(j), ENTIER_D))
          end do
          nchk = int(n_ids)
        end if
        if (nvol == nchk) then
          if (me == 0) print '(a,a,a)', " [cgns] BC '", trim(boconame), &
            "' spans 3D cells (a volume group, not a boundary condition): ignored"
          deallocate (pnts)
          cycle
        else if (nvol > 0) then
          call fatal("CGNS: BC '"//trim(boconame)//"' mixes 3D cells and faces")
        end if
      end if

      if (me == 0) then
        if (kmatch > 0) then
          print '(a,a,a,a,a,i0,a)', " [cgns] BC '", trim(boconame), "' (family '", &
            trim(famname), "') -> bc_name '"//trim(adjustl(bc_name(kmatch)))//"' (tag ", kmatch, ")"
        else
          print '(a,a,a,a,a)', " [cgns] BC '", trim(boconame), "' (family '", trim(famname), &
            "') matches no bc_name entry -> tag 0"
        end if
      end if

      i0 = (n_ids*int(me, ENTIER_D))/int(num_procs, ENTIER_D)
      i1 = (n_ids*int(me + 1, ENTIER_D))/int(num_procs, ENTIER_D)
      do j = i0 + 1, i1
        if (ptset == PointRange .or. ptset == ElementRange) then
          g = int(pnts(1), ENTIER_D) + j - 1
        else
          g = int(pnts(j), ENTIER_D)
        end if
        if (n_pairs >= size(face_id)) then
          allocate (tmp1(2*size(face_id)), tmp2(2*size(face_id)))
          tmp1(1:n_pairs) = face_id(1:n_pairs)
          tmp2(1:n_pairs) = bc_idx(1:n_pairs)
          call move_alloc(tmp1, face_id)
          call move_alloc(tmp2, bc_idx)
        end if
        n_pairs = n_pairs + 1
        face_id(n_pairs) = int(g, ENTIER)
        bc_idx(n_pairs) = kmatch
      end do
      deallocate (pnts)
    end do

    if (me == 0) then
      do k = 1, n_bc
        if (.not. bc_used(k)) print '(a,a,a)', achar(27)//"[33m [cgns] warning: bc_name '", &
          trim(adjustl(bc_name(k))), "' matches no BC in the CGNS file"//achar(27)//"[0m"
      end do
    end if

  contains

    ! 1 if the element id lies in a section holding only 3D cells.
    integer function count_vol(eid)
      integer(kind=ENTIER_D), intent(in) :: eid
      integer :: sec_i
      count_vol = 0
      do sec_i = 1, cf%n_sections
        if (eid >= cf%sec(sec_i)%istart .and. eid <= cf%sec(sec_i)%iend) then
          if (cf%sec(sec_i)%in_vol_space .and. .not. cf%sec(sec_i)%in_face_space) count_vol = 1
          return
        end if
      end do
    end function count_vol
  end subroutine cgns_read_bc_faces
end module cgns_dist_read_module
