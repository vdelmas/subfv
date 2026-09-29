! Distributed mesh setup from CGNS: block read, ParMETIS partition, migration,
! node-directory ghost layers and BC matching; output identical to read_mesh_msh_4.
! Env: SUBFV_PARTITIONER=parmetis|block, SUBFV_GHOST_LAYERS=n (default 1).
module dist_mesh_module
  use mpi
  use precision_module
  use mpi_module
  use mesh_module
  use dist_utils_module
  use cgns_dist_read_module
  use iso_c_binding, only: c_int
  implicit none

  private
  public :: read_mesh_cgns_dist

  ! Cell record: gid, kind, nvert, verts(nvert), nface, (nv, verts(nv)) * nface
  type :: cell_store_type
    integer(kind=ENTIER) :: n = 0
    integer(kind=ENTIER), dimension(:), allocatable :: ptr ! n+1
    integer(kind=ENTIER), dimension(:), allocatable :: dat
  end type cell_store_type

  interface
    ! ParMETIS' own Fortran entry point (frename.c): takes the Fortran communicator.
    ! idx_t/real_t must be 32-bit (IDXTYPEWIDTH/REALTYPEWIDTH, checked in CMake).
    integer(c_int) function parmetis_v3_partmeshkway(elmdist, eptr, eind, elmwgt, wgtflag, &
        numflag, ncon, ncommonnodes, nparts, tpwgts, ubvec, options, edgecut, part, comm) &
        bind(C, name="parmetis_v3_partmeshkway")
      use iso_c_binding, only: c_int, c_float, c_ptr
      integer(c_int), dimension(*), intent(in) :: elmdist, eptr, eind
      type(c_ptr), value :: elmwgt
      integer(c_int), intent(in) :: wgtflag, numflag, ncon, ncommonnodes, nparts
      real(c_float), dimension(*), intent(in) :: tpwgts, ubvec
      integer(c_int), dimension(*), intent(in) :: options
      integer(c_int), intent(out) :: edgecut
      integer(c_int), dimension(*), intent(out) :: part
      integer, intent(in) :: comm
    end function parmetis_v3_partmeshkway
  end interface

contains

  subroutine fatal(msg)
    character(len=*), intent(in) :: msg
    integer :: mpi_ierr

    print *, achar(27)//"[31m[-] "//trim(msg)//achar(27)//"[0m"
    call MPI_ABORT(MPI_COMM_WORLD, 1, mpi_ierr)
  end subroutine fatal

  ! Record accessors
  pure integer(kind=ENTIER) function rec_gid(cs, i)
    type(cell_store_type), intent(in) :: cs
    integer(kind=ENTIER), intent(in) :: i
    rec_gid = cs%dat(cs%ptr(i))
  end function rec_gid

  pure integer(kind=ENTIER) function rec_nvert(cs, i)
    type(cell_store_type), intent(in) :: cs
    integer(kind=ENTIER), intent(in) :: i
    rec_nvert = cs%dat(cs%ptr(i) + 2)
  end function rec_nvert

  ! position of vertex j (1-based) of cell i
  pure integer(kind=ENTIER) function rec_vpos(cs, i, j)
    type(cell_store_type), intent(in) :: cs
    integer(kind=ENTIER), intent(in) :: i, j
    rec_vpos = cs%ptr(i) + 2 + j
  end function rec_vpos

  pure integer(kind=ENTIER) function rec_nface(cs, i)
    type(cell_store_type), intent(in) :: cs
    integer(kind=ENTIER), intent(in) :: i
    rec_nface = cs%dat(cs%ptr(i) + 3 + rec_nvert(cs, i))
  end function rec_nface

  ! position of the first face header (nv) of cell i
  pure integer(kind=ENTIER) function rec_fpos(cs, i)
    type(cell_store_type), intent(in) :: cs
    integer(kind=ENTIER), intent(in) :: i
    rec_fpos = cs%ptr(i) + 4 + rec_nvert(cs, i)
  end function rec_fpos

  subroutine store_append(cs, rec)
    type(cell_store_type), intent(inout) :: cs
    integer(kind=ENTIER), dimension(:), intent(in) :: rec

    integer(kind=ENTIER), dimension(:), allocatable :: tmp
    integer(kind=ENTIER) :: ndat

    if (.not. allocated(cs%ptr)) then
      allocate (cs%ptr(64), cs%dat(1024))
      cs%n = 0
      cs%ptr(1) = 1
    end if
    ndat = cs%ptr(cs%n + 1) - 1
    if (cs%n + 2 > size(cs%ptr)) then
      allocate (tmp(2*size(cs%ptr)))
      tmp(1:cs%n + 1) = cs%ptr(1:cs%n + 1)
      call move_alloc(tmp, cs%ptr)
    end if
    if (ndat + size(rec) > size(cs%dat)) then
      allocate (tmp(2*size(cs%dat) + size(rec)))
      tmp(1:ndat) = cs%dat(1:ndat)
      call move_alloc(tmp, cs%dat)
    end if
    cs%dat(ndat + 1:ndat + size(rec)) = rec
    cs%n = cs%n + 1
    cs%ptr(cs%n + 1) = ndat + size(rec) + 1
  end subroutine store_append

  subroutine store_init(cs)
    type(cell_store_type), intent(inout) :: cs
    if (allocated(cs%ptr)) deallocate (cs%ptr)
    if (allocated(cs%dat)) deallocate (cs%dat)
    allocate (cs%ptr(64), cs%dat(1024))
    cs%n = 0
    cs%ptr(1) = 1
  end subroutine store_init

  ! Reorders a store by ascending global id.
  subroutine store_sort_by_gid(cs)
    type(cell_store_type), intent(inout) :: cs

    type(cell_store_type) :: tmp
    integer(kind=ENTIER), dimension(:), allocatable :: key, perm
    integer(kind=ENTIER) :: i

    allocate (key(cs%n), perm(cs%n))
    do i = 1, cs%n
      key(i) = rec_gid(cs, i)
    end do
    call sort_perm_int(cs%n, key, perm)
    call store_init(tmp)
    do i = 1, cs%n
      call store_append(tmp, cs%dat(cs%ptr(perm(i)):cs%ptr(perm(i) + 1) - 1))
    end do
    call move_alloc(tmp%ptr, cs%ptr)
    call move_alloc(tmp%dat, cs%dat)
  end subroutine store_sort_by_gid

  ! Standard cell faces, same tables and orientation as read_mesh_msh_4.
  subroutine build_std_record(gid, kind, v, rec)
    integer(kind=ENTIER), intent(in) :: gid, kind
    integer(kind=ENTIER), dimension(:), intent(in) :: v
    integer(kind=ENTIER), dimension(:), allocatable, intent(out) :: rec

    integer(kind=ENTIER), dimension(4, 4), parameter :: tet = reshape([ &
      1, 3, 2, 0, 1, 2, 4, 0, 1, 4, 3, 0, 2, 3, 4, 0], [4, 4])
    integer(kind=ENTIER), dimension(4, 6), parameter :: hex = reshape([ &
      1, 4, 3, 2, 3, 4, 8, 7, 7, 8, 5, 6, 1, 2, 6, 5, 1, 5, 8, 4, 2, 3, 7, 6], [4, 6])
    integer(kind=ENTIER), dimension(4, 5), parameter :: pri = reshape([ &
      1, 3, 2, 0, 4, 5, 6, 0, 1, 2, 5, 4, 2, 3, 6, 5, 1, 4, 6, 3], [4, 5])
    integer(kind=ENTIER), dimension(4, 5), parameter :: pyr = reshape([ &
      1, 4, 3, 2, 1, 2, 5, 0, 2, 3, 5, 0, 3, 4, 5, 0, 1, 5, 4, 0], [4, 5])

    select case (kind)
    case (4)
      call fill(tet)
    case (5)
      call fill(hex)
    case (6)
      call fill(pri)
    case (7)
      call fill(pyr)
    case default
      call fatal("dist_mesh: unknown standard cell kind")
    end select

  contains

    subroutine fill(tab)
      integer(kind=ENTIER), dimension(:, :), intent(in) :: tab

      integer(kind=ENTIER) :: nf, f, k, nv, pos, ntot

      nf = size(tab, 2)
      ntot = 4 + size(v)
      do f = 1, nf
        ntot = ntot + 1 + count(tab(:, f) > 0)
      end do
      allocate (rec(ntot))
      rec(1) = gid
      rec(2) = kind
      rec(3) = size(v)
      rec(4:3 + size(v)) = v
      rec(4 + size(v)) = nf
      pos = 5 + size(v)
      do f = 1, nf
        nv = count(tab(:, f) > 0)
        rec(pos) = nv
        do k = 1, nv
          rec(pos + k) = v(tab(k, f))
        end do
        pos = pos + nv + 1
      end do
    end subroutine fill
  end subroutine build_std_record

  ! Fetch variable-length integer data by global index from a block-distributed
  ! store; results come back in request order.
  subroutine fetch_csr(me, num_procs, starts, n_req, req, lptr, ldat, out_ptr, out_dat)
    integer(kind=ENTIER), intent(in) :: me, num_procs, n_req
    integer(kind=ENTIER_D), dimension(0:num_procs), intent(in) :: starts
    integer(kind=ENTIER), dimension(:), intent(in) :: req
    integer(kind=ENTIER), dimension(:), intent(in) :: lptr, ldat
    integer(kind=ENTIER), dimension(:), allocatable, intent(out) :: out_ptr, out_dat

    integer(kind=ENTIER), dimension(:), allocatable :: owner, perm, sendbuf, recvbuf, &
      replybuf, ansbuf, lens, ans_pos
    integer(kind=ENTIER), dimension(0:num_procs - 1) :: scount, rcount, rep_count, acount
    integer(kind=ENTIER) :: i, t, r, l, pos, k, n_in, ln

    allocate (owner(n_req), perm(n_req), sendbuf(n_req))
    do i = 1, n_req
      owner(i) = block_owner(int(req(i), ENTIER_D), num_procs, starts)
    end do
    call sort_perm_int(n_req, owner, perm)
    scount = 0
    do t = 1, n_req
      sendbuf(t) = req(perm(t))
      scount(owner(perm(t))) = scount(owner(perm(t))) + 1
    end do
    call exchange_int(num_procs, scount, sendbuf, rcount, recvbuf)

    ! Answer: len, data... per request, in received order.
    n_in = size(recvbuf)
    rep_count = 0
    pos = 0
    do r = 0, num_procs - 1
      do k = 1, rcount(r)
        pos = pos + 1
        l = int(int(recvbuf(pos), ENTIER_D) - starts(me), ENTIER)
        if (l < 1 .or. l > size(lptr) - 1) call fatal("dist_mesh: fetch out of block (internal error)")
        rep_count(r) = rep_count(r) + 1 + lptr(l + 1) - lptr(l)
      end do
    end do
    allocate (replybuf(sum(rep_count)))
    pos = 0
    k = 0
    do i = 1, n_in
      l = int(int(recvbuf(i), ENTIER_D) - starts(me), ENTIER)
      ln = lptr(l + 1) - lptr(l)
      replybuf(k + 1) = ln
      replybuf(k + 2:k + 1 + ln) = ldat(lptr(l):lptr(l + 1) - 1)
      k = k + 1 + ln
    end do
    call exchange_int(num_procs, rep_count, replybuf, acount, ansbuf)

    ! Answers arrive grouped by owner in the order the requests were sent,
    ! i.e. in perm order.
    allocate (lens(n_req), ans_pos(n_req))
    k = 0
    do t = 1, n_req
      lens(perm(t)) = ansbuf(k + 1)
      ans_pos(perm(t)) = k + 2
      k = k + 1 + ansbuf(k + 1)
    end do
    allocate (out_ptr(n_req + 1))
    out_ptr(1) = 1
    do i = 1, n_req
      out_ptr(i + 1) = out_ptr(i) + lens(i)
    end do
    allocate (out_dat(max(out_ptr(n_req + 1) - 1, 1)))
    do i = 1, n_req
      out_dat(out_ptr(i):out_ptr(i + 1) - 1) = ansbuf(ans_pos(i):ans_pos(i) + lens(i) - 1)
    end do
  end subroutine fetch_csr

  subroutine fetch_coords(me, num_procs, starts, n_req, req, lxyz, xyz)
    integer(kind=ENTIER), intent(in) :: me, num_procs, n_req
    integer(kind=ENTIER_D), dimension(0:num_procs), intent(in) :: starts
    integer(kind=ENTIER), dimension(:), intent(in) :: req
    real(kind=DOUBLE), dimension(:, :), intent(in) :: lxyz
    real(kind=DOUBLE), dimension(:, :), allocatable, intent(out) :: xyz

    integer(kind=ENTIER), dimension(:), allocatable :: owner, perm, sendbuf, recvbuf
    integer(kind=ENTIER), dimension(0:num_procs - 1) :: scount, rcount, rep_count, acount
    real(kind=DOUBLE), dimension(:), allocatable :: replybuf, ansbuf
    integer(kind=ENTIER) :: i, t, l

    allocate (owner(n_req), perm(n_req), sendbuf(n_req))
    do i = 1, n_req
      owner(i) = block_owner(int(req(i), ENTIER_D), num_procs, starts)
    end do
    call sort_perm_int(n_req, owner, perm)
    scount = 0
    do t = 1, n_req
      sendbuf(t) = req(perm(t))
      scount(owner(perm(t))) = scount(owner(perm(t))) + 1
    end do
    call exchange_int(num_procs, scount, sendbuf, rcount, recvbuf)

    allocate (replybuf(3*size(recvbuf)))
    do i = 1, size(recvbuf)
      l = int(int(recvbuf(i), ENTIER_D) - starts(me), ENTIER)
      replybuf(3*i - 2:3*i) = lxyz(:, l)
    end do
    rep_count = 3*rcount
    call exchange_dbl(num_procs, rep_count, replybuf, acount, ansbuf)

    allocate (xyz(3, n_req))
    do t = 1, n_req
      xyz(:, perm(t)) = ansbuf(3*t - 2:3*t)
    end do
  end subroutine fetch_coords

  ! Sends whole records to destination ranks; received records appended to
  ! dst, with their source rank in src (optional).
  subroutine send_records(num_procs, cs, n_items, item_cell, item_dest, dst, src)
    integer(kind=ENTIER), intent(in) :: num_procs, n_items
    type(cell_store_type), intent(in) :: cs
    integer(kind=ENTIER), dimension(:), intent(in) :: item_cell, item_dest
    type(cell_store_type), intent(inout) :: dst
    integer(kind=ENTIER), dimension(:), allocatable, intent(inout), optional :: src

    integer(kind=ENTIER), dimension(:), allocatable :: perm, sendbuf, recvbuf, tmp
    integer(kind=ENTIER), dimension(0:num_procs - 1) :: scount, rcount
    integer(kind=ENTIER) :: t, c, k, len, r, pos, n0, nn

    allocate (perm(n_items))
    call sort_perm_int(n_items, item_dest, perm)
    scount = 0
    do t = 1, n_items
      c = item_cell(perm(t))
      scount(item_dest(perm(t))) = scount(item_dest(perm(t))) + cs%ptr(c + 1) - cs%ptr(c)
    end do
    allocate (sendbuf(sum(scount)))
    k = 0
    do t = 1, n_items
      c = item_cell(perm(t))
      len = cs%ptr(c + 1) - cs%ptr(c)
      sendbuf(k + 1:k + len) = cs%dat(cs%ptr(c):cs%ptr(c + 1) - 1)
      k = k + len
    end do
    call exchange_int(num_procs, scount, sendbuf, rcount, recvbuf)

    n0 = dst%n
    pos = 1
    do r = 0, num_procs - 1
      k = pos + rcount(r)
      do while (pos < k)
        len = record_length(recvbuf, pos)
        call store_append(dst, recvbuf(pos:pos + len - 1))
        pos = pos + len
        if (present(src)) then
          nn = dst%n
          if (.not. allocated(src)) allocate (src(max(64, 2*nn)))
          if (nn > size(src)) then
            allocate (tmp(2*nn))
            tmp(1:size(src)) = src
            call move_alloc(tmp, src)
          end if
          src(nn) = r
        end if
      end do
    end do
  end subroutine send_records

  pure integer(kind=ENTIER) function record_length(buf, pos)
    integer(kind=ENTIER), dimension(:), intent(in) :: buf
    integer(kind=ENTIER), intent(in) :: pos

    integer(kind=ENTIER) :: p, nf, f

    p = pos + 3 + buf(pos + 2)
    nf = buf(p)
    p = p + 1
    do f = 1, nf
      p = p + 1 + buf(p)
    end do
    record_length = p - pos
  end function record_length

  ! Unique sorted list of every vertex of cells 1..n of cs.
  subroutine store_unique_nodes(cs, n, nodes)
    type(cell_store_type), intent(in) :: cs
    integer(kind=ENTIER), intent(in) :: n
    integer(kind=ENTIER), dimension(:), allocatable, intent(out) :: nodes

    integer(kind=ENTIER), dimension(:), allocatable :: all, perm
    integer(kind=ENTIER) :: i, j, k, m

    m = 0
    do i = 1, n
      m = m + rec_nvert(cs, i)
    end do
    allocate (all(m), perm(m))
    k = 0
    do i = 1, n
      do j = 1, rec_nvert(cs, i)
        k = k + 1
        all(k) = cs%dat(rec_vpos(cs, i, j))
      end do
    end do
    call sort_perm_int(m, all, perm)
    allocate (nodes(m))
    k = 0
    do i = 1, m
      if (k == 0) then
        k = 1
        nodes(1) = all(perm(i))
      else if (all(perm(i)) /= nodes(k)) then
        k = k + 1
        nodes(k) = all(perm(i))
      end if
    end do
    nodes = nodes(1:k)
  end subroutine store_unique_nodes

  pure integer(kind=ENTIER) function bsearch(n, sorted, key)
    integer(kind=ENTIER), intent(in) :: n, key
    integer(kind=ENTIER), dimension(:), intent(in) :: sorted

    integer(kind=ENTIER) :: lo, hi, mid

    bsearch = 0
    lo = 1
    hi = n
    do while (lo <= hi)
      mid = (lo + hi)/2
      if (sorted(mid) == key) then
        bsearch = mid
        return
      else if (sorted(mid) < key) then
        lo = mid + 1
      else
        hi = mid - 1
      end if
    end do
  end function bsearch

  ! Insertion sort, for the tiny per-face node lists.
  pure subroutine small_sort(a)
    integer(kind=ENTIER), dimension(:), intent(inout) :: a
    integer(kind=ENTIER) :: i, j, t

    do i = 2, size(a)
      t = a(i)
      j = i - 1
      do while (j >= 1)
        if (a(j) <= t) exit
        a(j + 1) = a(j)
        j = j - 1
      end do
      a(j + 1) = t
    end do
  end subroutine small_sort

  ! =====================================================================
  subroutine read_mesh_cgns_dist(mesh, meshfile_path, meshfile, n_bc, bc_name, me, &
      num_procs, mpi_send_recv)
    type(mesh_type), intent(inout) :: mesh
    character(len=255), intent(in) :: meshfile, meshfile_path
    integer(kind=ENTIER), intent(in) :: n_bc, me, num_procs
    character(len=255), dimension(:), intent(in) :: bc_name
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    type(cgns_mesh_file_type) :: cf
    integer(kind=ENTIER_D), dimension(0:num_procs) :: node_starts, vol_starts, face_starts
    integer(kind=ENTIER) :: n_raw, i, j, k, f, nv, pos, mpi_ierr, n_layers, layer
    integer(kind=ENTIER), dimension(:), allocatable :: raw_gid, raw_kind, raw_ptr, raw_dat
    integer(kind=ENTIER), dimension(:), allocatable :: fs_ptr, fs_dat
    real(kind=DOUBLE), dimension(:, :), allocatable :: node_xyz, held_xyz
    type(cell_store_type) :: loaded, held
    integer(kind=ENTIER), dimension(:), allocatable :: rec, part, ghost_src, held_nodes
    integer(kind=ENTIER) :: n_owned, edgecut
    character(len=64) :: partitioner, env
    integer(kind=ENTIER) :: env_len, env_stat
    ! send pairs accumulated over layers: (dest rank, owned cell index)
    integer(kind=ENTIER), dimension(:), allocatable :: sp_dest, sp_cell
    integer(kind=ENTIER) :: n_sp
    ! global id of every assembled (local) cell
    integer(kind=ENTIER), dimension(:), allocatable :: elem_gid
    real(kind=DOUBLE) :: t0, t_read, t_part, t_ghost, t_bc

    t0 = MPI_WTIME()

    ! ---------------- options ----------------
#ifdef SUBFV_HAVE_PARMETIS
    partitioner = "parmetis"
#else
    partitioner = "block"
#endif
    call get_environment_variable("SUBFV_PARTITIONER", env, env_len, env_stat)
    if (env_stat == 0 .and. env_len > 0) partitioner = trim(env)
    n_layers = 1
    call get_environment_variable("SUBFV_GHOST_LAYERS", env, env_len, env_stat)
    if (env_stat == 0 .and. env_len > 0) read (env, *) n_layers
    if (n_layers < 1) call fatal("SUBFV_GHOST_LAYERS must be >= 1")
#ifndef SUBFV_HAVE_PARMETIS
    if (trim(partitioner) == "parmetis") call fatal("subfv built without ParMETIS")
#endif
    if (trim(partitioner) /= "parmetis" .and. trim(partitioner) /= "block") &
      call fatal("SUBFV_PARTITIONER must be 'parmetis' or 'block'")

    ! ---------------- 1. block reads ----------------
    call cgns_open_mesh(trim(adjustl(meshfile_path))//trim(adjustl(meshfile)), cf)
    call block_starts(cf%n_nodes, num_procs, node_starts)
    call block_starts(cf%n_vol_space, num_procs, vol_starts)
    call block_starts(cf%n_face_space, num_procs, face_starts)

    call cgns_read_cell_block(cf, vol_starts(me) + 1, vol_starts(me + 1), n_raw, &
      raw_gid, raw_kind, raw_ptr, raw_dat)
    call cgns_read_face_block(cf, face_starts(me) + 1, face_starts(me + 1), fs_ptr, fs_dat)
    call cgns_read_coord_block(cf, node_starts(me) + 1, node_starts(me + 1), node_xyz)

    ! ---------------- 2. cell records (polyhedra resolved) ----------------
    call build_records()
    t_read = MPI_WTIME() - t0

    ! ---------------- 3-4. partition + migration ----------------
    call partition_cells()
    call migrate()
    t_part = MPI_WTIME() - t0 - t_read

    ! ---------------- 5. ghost layers ----------------
    n_sp = 0
    allocate (sp_dest(0), sp_cell(0))
    do layer = 1, n_layers
      call ghost_round()
    end do
    t_ghost = MPI_WTIME() - t0 - t_read - t_part

    ! ---------------- 6. coordinates ----------------
    call store_unique_nodes(held, held%n, held_nodes)
    call fetch_coords(me, num_procs, node_starts, size(held_nodes), held_nodes, node_xyz, held_xyz)
    deallocate (node_xyz)

    ! ---------------- 7-8. BC faces + assembly ----------------
    call assemble_mesh()
    call match_boundary_faces()
    t_bc = MPI_WTIME() - t0 - t_read - t_part - t_ghost

    call cgns_close_mesh(cf)
    call check_send_recv_consistency()
    call print_summary()

  contains

    ! Standard cells use the msh face tables; NFACE_n faces are fetched and
    ! reversed when their sign is negative (outward orientation).
    subroutine build_records()
      integer(kind=ENTIER), dimension(:), allocatable :: req, fptr, fdat
      integer(kind=ENTIER) :: c, n_req, q, fid, nf, ntot, nvert
      integer(kind=ENTIER), dimension(:), allocatable :: allv, perm, verts
      integer(kind=ENTIER_D) :: fidx

      ! Every rank takes part in the fetch (collective), even with no polyhedra.
      n_req = 0
      do c = 1, n_raw
        if (raw_kind(c) == KIND_POLYHEDRON) n_req = n_req + raw_ptr(c + 1) - raw_ptr(c)
      end do
      allocate (req(n_req))
      q = 0
      do c = 1, n_raw
        if (raw_kind(c) /= KIND_POLYHEDRON) cycle
        do j = raw_ptr(c), raw_ptr(c + 1) - 1
          q = q + 1
          fidx = cgns_face_space_index(cf, int(abs(raw_dat(j)), ENTIER_D))
          if (fidx == 0) call fatal("CGNS: NFACE_n references an element that is not a face")
          req(q) = int(fidx, ENTIER)
        end do
      end do
      call fetch_csr(me, num_procs, face_starts, n_req, req, fs_ptr, fs_dat, fptr, fdat)

      call store_init(loaded)
      q = 0
      do c = 1, n_raw
        if (raw_kind(c) /= KIND_POLYHEDRON) then
          call build_std_record(raw_gid(c), raw_kind(c), &
            raw_dat(raw_ptr(c):raw_ptr(c + 1) - 1), rec)
        else
          nf = raw_ptr(c + 1) - raw_ptr(c)
          ! unique vertex list = union of face nodes
          ntot = 0
          do f = 1, nf
            ntot = ntot + fptr(q + f + 1) - fptr(q + f)
          end do
          allocate (allv(ntot), perm(ntot), verts(ntot))
          k = 0
          do f = 1, nf
            do j = fptr(q + f), fptr(q + f + 1) - 1
              k = k + 1
              allv(k) = fdat(j)
            end do
          end do
          call sort_perm_int(ntot, allv, perm)
          nvert = 0
          do j = 1, ntot
            if (nvert == 0) then
              nvert = 1
              verts(1) = allv(perm(j))
            else if (allv(perm(j)) /= verts(nvert)) then
              nvert = nvert + 1
              verts(nvert) = allv(perm(j))
            end if
          end do
          allocate (rec(4 + nvert + nf + ntot))
          rec(1) = raw_gid(c)
          rec(2) = KIND_POLYHEDRON
          rec(3) = nvert
          rec(4:3 + nvert) = verts(1:nvert)
          rec(4 + nvert) = nf
          pos = 5 + nvert
          do f = 1, nf
            fid = raw_dat(raw_ptr(c) + f - 1)
            nv = fptr(q + f + 1) - fptr(q + f)
            if (nv < 3) call fatal("CGNS: polyhedron face with fewer than 3 nodes")
            rec(pos) = nv
            if (fid > 0) then
              rec(pos + 1:pos + nv) = fdat(fptr(q + f):fptr(q + f + 1) - 1)
            else
              rec(pos + 1:pos + nv) = fdat(fptr(q + f + 1) - 1:fptr(q + f):-1)
            end if
            pos = pos + nv + 1
          end do
          q = q + nf
          deallocate (allv, perm, verts)
        end if
        call store_append(loaded, rec)
        deallocate (rec)
      end do
      deallocate (raw_gid, raw_kind, raw_ptr, raw_dat)
    end subroutine build_records

    subroutine partition_cells()
      integer(kind=ENTIER), dimension(0:num_procs) :: elmdist
      integer(kind=ENTIER), dimension(:), allocatable :: eptr, eind
      integer(kind=ENTIER) :: c, nloc, stat, m
      integer(kind=ENTIER), dimension(num_procs) :: counts

      nloc = loaded%n
      allocate (part(max(nloc, 1)))
      edgecut = -1
      if (num_procs == 1 .or. trim(partitioner) == "block") then
        part(1:nloc) = me
        return
      end if

      call MPI_ALLGATHER(nloc, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, mpi_ierr)
      if (minval(counts) < 1) call fatal("ParMETIS needs at least one cell per rank")
      elmdist(0) = 0
      do c = 1, num_procs
        elmdist(c) = elmdist(c - 1) + counts(c)
      end do

      m = 0
      do c = 1, nloc
        m = m + rec_nvert(loaded, c)
      end do
      allocate (eptr(nloc + 1), eind(max(m, 1)))
      eptr(1) = 0
      k = 0
      do c = 1, nloc
        do j = 1, rec_nvert(loaded, c)
          k = k + 1
          eind(k) = loaded%dat(rec_vpos(loaded, c, j)) - 1
        end do
        eptr(c + 1) = k
      end do

#ifdef SUBFV_HAVE_PARMETIS
      ! ncommonnodes = 3: two 3D cells are dual-graph neighbours when they
      ! share a face (any face has >= 3 nodes), not merely an edge.
      block
        use iso_c_binding, only: c_null_ptr, c_float
        real(c_float), dimension(num_procs) :: tpwgts
        real(c_float), dimension(1) :: ubvec
        integer(c_int), dimension(3) :: options
        tpwgts = 1.0_c_float/real(num_procs, c_float)
        ubvec = 1.05_c_float
        options = 0
        stat = parmetis_v3_partmeshkway(elmdist, eptr, eind, c_null_ptr, 0, 0, 1, 3, num_procs, &
          tpwgts, ubvec, options, edgecut, part, MPI_COMM_WORLD)
      end block
      if (stat /= 1) call fatal("ParMETIS_V3_PartMeshKway failed")  ! METIS_OK = 1
#else
      stat = 0
#endif
    end subroutine partition_cells

    subroutine migrate()
      integer(kind=ENTIER), dimension(:), allocatable :: items
      integer(kind=ENTIER) :: c

      allocate (items(loaded%n))
      do c = 1, loaded%n
        items(c) = c
      end do
      call store_init(held)
      call send_records(num_procs, loaded, loaded%n, items, part, held)
      deallocate (part)
      call store_init(loaded)
      call store_sort_by_gid(held)
      n_owned = held%n
      allocate (ghost_src(max(n_owned, 64)))
      ghost_src = -1
    end subroutine migrate

    ! One ghost layer: owned cell c goes to every rank touching one of its nodes.
    subroutine ghost_round()
      integer(kind=ENTIER), dimension(:), allocatable :: nodes, owner, perm, sendbuf, recvbuf
      integer(kind=ENTIER), dimension(:), allocatable :: pn, pr, rbuf, dir_node, dir_ptr, dir_rank
      integer(kind=ENTIER), dimension(:), allocatable :: new_dest, new_cell, keyq, keyc, p1, p2
      integer(kind=ENTIER), dimension(:), allocatable :: tmp_d, tmp_c, cand_q, cand_c
      integer(kind=ENTIER), dimension(0:num_procs - 1) :: scount, rcount, rep_count, acount
      integer(kind=ENTIER) :: n, t, r, g0, g1, m, a, b, idx, c, nc, nn, nd, iq, np
      logical :: dup

      ! (a) register held nodes with their directory owner
      call store_unique_nodes(held, held%n, nodes)
      n = size(nodes)
      allocate (owner(n))
      do i = 1, n
        owner(i) = block_owner(int(nodes(i), ENTIER_D), num_procs, node_starts)
      end do
      scount = 0
      do i = 1, n
        scount(owner(i)) = scount(owner(i)) + 1
      end do
      ! nodes sorted => owners non-decreasing => already grouped by owner
      call exchange_int(num_procs, scount, nodes, rcount, recvbuf)

      ! (b) directory: (node, rank) pairs sorted by node (stable => ranks ascending)
      np = size(recvbuf)
      allocate (pn(np), pr(np), perm(np))
      k = 0
      do r = 0, num_procs - 1
        do t = 1, rcount(r)
          k = k + 1
          pn(k) = recvbuf(k)
          pr(k) = r
        end do
      end do
      call sort_perm_int(np, pn, perm)

      ! reply [node, m, ranks(m)] to every rank of each shared node
      rep_count = 0
      g0 = 1
      do while (g0 <= np)
        g1 = g0
        do while (g1 < np)
          if (pn(perm(g1 + 1)) /= pn(perm(g0))) exit
          g1 = g1 + 1
        end do
        m = g1 - g0 + 1
        if (m >= 2) then
          do t = g0, g1
            rep_count(pr(perm(t))) = rep_count(pr(perm(t))) + 2 + m
          end do
        end if
        g0 = g1 + 1
      end do
      allocate (sendbuf(sum(rep_count)))
      block
        integer(kind=ENTIER), dimension(0:num_procs - 1) :: off
        off(0) = 0
        do r = 1, num_procs - 1
          off(r) = off(r - 1) + rep_count(r - 1)
        end do
        g0 = 1
        do while (g0 <= np)
          g1 = g0
          do while (g1 < np)
            if (pn(perm(g1 + 1)) /= pn(perm(g0))) exit
            g1 = g1 + 1
          end do
          m = g1 - g0 + 1
          if (m >= 2) then
            do t = g0, g1
              r = pr(perm(t))
              sendbuf(off(r) + 1) = pn(perm(g0))
              sendbuf(off(r) + 2) = m
              do a = 0, m - 1
                sendbuf(off(r) + 3 + a) = pr(perm(g0 + a))
              end do
              off(r) = off(r) + 2 + m
            end do
          end if
          g0 = g1 + 1
        end do
      end block
      call exchange_int(num_procs, rep_count, sendbuf, acount, rbuf)
      deallocate (sendbuf, pn, pr, perm, recvbuf)

      ! (c) node -> other ranks, sorted by node for binary search
      nd = 0
      k = 1
      do while (k <= size(rbuf))
        nd = nd + 1
        k = k + 2 + rbuf(k + 1)
      end do
      allocate (dir_node(nd), dir_ptr(nd + 1), dir_rank(max(size(rbuf), 1)), perm(nd))
      block
        integer(kind=ENTIER), dimension(:), allocatable :: raw_node, raw_pos
        allocate (raw_node(nd), raw_pos(nd))
        k = 1
        do i = 1, nd
          raw_node(i) = rbuf(k)
          raw_pos(i) = k
          k = k + 2 + rbuf(k + 1)
        end do
        call sort_perm_int(nd, raw_node, perm)
        dir_ptr(1) = 1
        do i = 1, nd
          dir_node(i) = raw_node(perm(i))
          m = rbuf(raw_pos(perm(i)) + 1)
          dir_rank(dir_ptr(i):dir_ptr(i) + m - 1) = &
            rbuf(raw_pos(perm(i)) + 2:raw_pos(perm(i)) + 1 + m)
          dir_ptr(i + 1) = dir_ptr(i) + m
        end do
      end block

      ! (d) candidate (dest, owned cell) pairs
      allocate (cand_q(64), cand_c(64))
      nc = 0
      do c = 1, n_owned
        do j = 1, rec_nvert(held, c)
          idx = bsearch(nd, dir_node, held%dat(rec_vpos(held, c, j)))
          if (idx == 0) cycle
          do t = dir_ptr(idx), dir_ptr(idx + 1) - 1
            if (dir_rank(t) == me) cycle
            if (nc >= size(cand_q)) then
              allocate (tmp_d(2*size(cand_q)), tmp_c(2*size(cand_q)))
              tmp_d(1:nc) = cand_q(1:nc)
              tmp_c(1:nc) = cand_c(1:nc)
              call move_alloc(tmp_d, cand_q)
              call move_alloc(tmp_c, cand_c)
            end if
            nc = nc + 1
            cand_q(nc) = dir_rank(t)
            cand_c(nc) = c
          end do
        end do
      end do

      ! (e) unique, lexicographic (dest, cell); drop pairs already sent
      allocate (keyq(nc), keyc(nc), p1(nc), p2(nc))
      call sort_perm_int(nc, cand_c(1:nc), p1)
      keyq(1:nc) = cand_q(p1)
      call sort_perm_int(nc, keyq, p2)
      keyq = cand_q(p1(p2))
      keyc = cand_c(p1(p2))
      allocate (new_dest(nc), new_cell(nc))
      nn = 0
      iq = 1
      do t = 1, nc
        if (t > 1) then
          if (keyq(t) == keyq(t - 1) .and. keyc(t) == keyc(t - 1)) cycle
        end if
        ! previous pairs are kept sorted by (dest, cell): advance a cursor
        dup = .false.
        do while (iq <= n_sp)
          if (sp_dest(iq) < keyq(t) .or. (sp_dest(iq) == keyq(t) .and. sp_cell(iq) < keyc(t))) then
            iq = iq + 1
          else
            exit
          end if
        end do
        if (iq <= n_sp) then
          if (sp_dest(iq) == keyq(t) .and. sp_cell(iq) == keyc(t)) dup = .true.
        end if
        if (.not. dup) then
          nn = nn + 1
          new_dest(nn) = keyq(t)
          new_cell(nn) = keyc(t)
        end if
      end do

      ! (f) ship new ghosts (copied out first: held also receives them)
      block
        type(cell_store_type) :: outgoing
        integer(kind=ENTIER), dimension(:), allocatable :: items
        call store_init(outgoing)
        allocate (items(nn))
        do t = 1, nn
          call store_append(outgoing, held%dat(held%ptr(new_cell(t)):held%ptr(new_cell(t) + 1) - 1))
          items(t) = t
        end do
        call send_records(num_procs, outgoing, nn, items, new_dest(1:nn), held, ghost_src)
      end block

      ! merge new pairs into the sorted accumulated list
      allocate (tmp_d(n_sp + nn), tmp_c(n_sp + nn))
      a = 1
      b = 1
      k = 0
      do while (a <= n_sp .or. b <= nn)
        k = k + 1
        if (b > nn) then
          tmp_d(k) = sp_dest(a); tmp_c(k) = sp_cell(a); a = a + 1
        else if (a > n_sp) then
          tmp_d(k) = new_dest(b); tmp_c(k) = new_cell(b); b = b + 1
        else if (sp_dest(a) < new_dest(b) .or. &
            (sp_dest(a) == new_dest(b) .and. sp_cell(a) < new_cell(b))) then
          tmp_d(k) = sp_dest(a); tmp_c(k) = sp_cell(a); a = a + 1
        else
          tmp_d(k) = new_dest(b); tmp_c(k) = new_cell(b); b = b + 1
        end if
      end do
      call move_alloc(tmp_d, sp_dest)
      call move_alloc(tmp_c, sp_cell)
      n_sp = n_sp + nn
    end subroutine ghost_round

    ! Fill mesh_type/mpi_send_recv like read_mesh_msh_4: owned then ghost cells by gid.
    subroutine assemble_mesh()
      integer(kind=ENTIER), dimension(:), allocatable :: gkey, perm, order, src_sorted
      integer(kind=ENTIER), dimension(:), allocatable :: loc_of_held
      integer(kind=ENTIER) :: c, ic, nf, fpos, iface, n_ghost, r, nsend, nrecv, is, ir
      integer(kind=ENTIER), dimension(0:num_procs - 1) :: n_to, n_from
      integer(kind=ENTIER) :: ip

      n_ghost = held%n - n_owned

      ! owned first (already by gid), ghosts sorted by gid
      allocate (order(held%n), gkey(n_ghost), perm(n_ghost))
      do c = 1, n_owned
        order(c) = c
      end do
      do c = 1, n_ghost
        gkey(c) = rec_gid(held, n_owned + c)
      end do
      call sort_perm_int(n_ghost, gkey, perm)
      do c = 1, n_ghost
        order(n_owned + c) = n_owned + perm(c)
      end do
      allocate (loc_of_held(held%n), elem_gid(held%n))
      do c = 1, held%n
        loc_of_held(order(c)) = c
        elem_gid(c) = rec_gid(held, order(c))
      end do

      mesh%n_vert = size(held_nodes)
      allocate (mesh%vert(mesh%n_vert))
      do i = 1, mesh%n_vert
        mesh%vert(i)%coord = held_xyz(:, i)
        mesh%vert(i)%id_glob = held_nodes(i)
      end do

      mesh%n_elems = held%n
      mesh%n_interior_elems = n_owned
      allocate (mesh%elem(mesh%n_elems))
      mesh%n_faces = 0
      do c = 1, held%n
        mesh%n_faces = mesh%n_faces + rec_nface(held, c)
      end do
      allocate (mesh%minimal_face(mesh%n_faces))

      iface = 0
      do ic = 1, held%n
        c = order(ic)
        mesh%elem(ic)%elem_kind = held%dat(held%ptr(c) + 1)
        mesh%elem(ic)%n_vert = rec_nvert(held, c)
        allocate (mesh%elem(ic)%vert(mesh%elem(ic)%n_vert))
        do j = 1, mesh%elem(ic)%n_vert
          mesh%elem(ic)%vert(j) = local_node(held%dat(rec_vpos(held, c, j)))
        end do
        nf = rec_nface(held, c)
        mesh%elem(ic)%n_faces = nf
        allocate (mesh%elem(ic)%face(nf))
        fpos = rec_fpos(held, c)
        do f = 1, nf
          iface = iface + 1
          nv = held%dat(fpos)
          mesh%elem(ic)%face(f) = iface
          mesh%minimal_face(iface)%n_vert = nv
          allocate (mesh%minimal_face(iface)%vert(nv))
          do j = 1, nv
            mesh%minimal_face(iface)%vert(j) = local_node(held%dat(fpos + j))
          end do
          mesh%minimal_face(iface)%left_neigh = ic
          mesh%minimal_face(iface)%left_neigh_face = f
          fpos = fpos + nv + 1
        end do
        mesh%elem(ic)%is_ghost = (ic > n_owned)
      end do

      ! ---- mpi_send_recv ----
      n_to = 0
      do ip = 1, n_sp
        n_to(sp_dest(ip)) = n_to(sp_dest(ip)) + 1
      end do
      n_from = 0
      do c = n_owned + 1, held%n
        n_from(ghost_src(c)) = n_from(ghost_src(c)) + 1
      end do
      nsend = count(n_to > 0)
      nrecv = count(n_from > 0)
      mpi_send_recv%n_mpi_send_neigh = nsend
      mpi_send_recv%n_mpi_recv_neigh = nrecv
      allocate (mpi_send_recv%mpi_send_neigh(nsend), mpi_send_recv%mpi_recv_neigh(nrecv))
      allocate (mpi_send_recv%mpi_reqsend(nsend), mpi_send_recv%mpi_reqrecv(nrecv))
      allocate (mpi_send_recv%mpi_sendstat(MPI_STATUS_SIZE, nsend))
      allocate (mpi_send_recv%mpi_recvstat(MPI_STATUS_SIZE, nrecv))

      ! send lists: sp_* sorted by (dest, owned index); owned index order is
      ! gid order, so each list is already sorted by global id.
      is = 0
      k = 1
      do r = 0, num_procs - 1
        if (n_to(r) == 0) cycle
        is = is + 1
        mpi_send_recv%mpi_send_neigh(is)%partition_id = r
        mpi_send_recv%mpi_send_neigh(is)%n_elems = n_to(r)
        allocate (mpi_send_recv%mpi_send_neigh(is)%elem_id(n_to(r)))
        do j = 1, n_to(r)
          mpi_send_recv%mpi_send_neigh(is)%elem_id(j) = loc_of_held(sp_cell(k))
          k = k + 1
        end do
      end do

      ! recv lists: ghosts are laid out by gid, filter per source
      ir = 0
      do r = 0, num_procs - 1
        if (n_from(r) == 0) cycle
        ir = ir + 1
        mpi_send_recv%mpi_recv_neigh(ir)%partition_id = r
        mpi_send_recv%mpi_recv_neigh(ir)%n_elems = n_from(r)
        allocate (mpi_send_recv%mpi_recv_neigh(ir)%elem_id(n_from(r)))
        k = 0
        do ic = n_owned + 1, held%n
          if (ghost_src(order(ic)) == r) then
            k = k + 1
            mpi_send_recv%mpi_recv_neigh(ir)%elem_id(k) = ic
          end if
        end do
      end do

      allocate (mpi_send_recv%is_ghost(mesh%n_elems))
      do ic = 1, mesh%n_elems
        mpi_send_recv%is_ghost(ic) = mesh%elem(ic)%is_ghost
      end do
    end subroutine assemble_mesh

    integer(kind=ENTIER) function local_node(g)
      integer(kind=ENTIER), intent(in) :: g
      local_node = bsearch(size(held_nodes), held_nodes, g)
      if (local_node == 0) call fatal("dist_mesh: node not held (internal error)")
    end function local_node

    ! Face directory keyed on sorted nodes: BC faces are sent to the rank owning
    ! the single cell that has that face.
    subroutine match_boundary_faces()
      integer(kind=ENTIER) :: n_pairs, c, fpos, nf, t, r, g0, g1, ncell, nbcf, kbc, dest
      integer(kind=ENTIER), dimension(:), allocatable :: bc_face, bc_k, req, fptr, fdat
      integer(kind=ENTIER), dimension(:), allocatable :: sendbuf, recvbuf, rbuf, key
      integer(kind=ENTIER), dimension(:), allocatable :: kptr, kdat, ktype, ksrc, kval, perm
      integer(kind=ENTIER), dimension(0:num_procs - 1) :: scount, rcount, rep_count, acount, off
      integer(kind=ENTIER_D) :: fidx
      integer(kind=ENTIER), dimension(8) :: stats, gstats
      integer(kind=ENTIER) :: nk, m, ib

      ! stats: 1 bc faces read, 2 matched, 3 internal (both sides), 4 unmatched,
      !        5 exterior untagged, 6 conflicting tags, 7 non-manifold
      stats = 0

      call cgns_read_bc_faces(cf, n_bc, bc_name, me, num_procs, n_pairs, bc_face, bc_k)
      allocate (req(n_pairs))
      do t = 1, n_pairs
        fidx = cgns_face_space_index(cf, int(bc_face(t), ENTIER_D))
        if (fidx == 0) call fatal("CGNS: a BC references an element that is not in a face section")
        req(t) = int(fidx, ENTIER)
      end do
      call fetch_csr(me, num_procs, face_starts, n_pairs, req, fs_ptr, fs_dat, fptr, fdat)
      deallocate (fs_ptr, fs_dat)
      stats(1) = n_pairs

      ! messages: [type(1=cell face, 2=bc face), value(k), nv, sorted nodes]
      scount = 0
      do c = 1, n_owned
        fpos = rec_fpos(held, c)
        do f = 1, rec_nface(held, c)
          nv = held%dat(fpos)
          dest = block_owner(int(minval(held%dat(fpos + 1:fpos + nv)), ENTIER_D), num_procs, node_starts)
          scount(dest) = scount(dest) + 3 + nv
          fpos = fpos + nv + 1
        end do
      end do
      do t = 1, n_pairs
        nv = fptr(t + 1) - fptr(t)
        if (nv < 3) call fatal("CGNS: a BC element is not a polygon face")
        dest = block_owner(int(minval(fdat(fptr(t):fptr(t + 1) - 1)), ENTIER_D), num_procs, node_starts)
        scount(dest) = scount(dest) + 3 + nv
      end do
      allocate (sendbuf(sum(scount)))
      off(0) = 0
      do r = 1, num_procs - 1
        off(r) = off(r - 1) + scount(r - 1)
      end do
      do c = 1, n_owned
        fpos = rec_fpos(held, c)
        do f = 1, rec_nface(held, c)
          nv = held%dat(fpos)
          allocate (key(nv))
          key = held%dat(fpos + 1:fpos + nv)
          call small_sort(key)
          dest = block_owner(int(key(1), ENTIER_D), num_procs, node_starts)
          sendbuf(off(dest) + 1) = 1
          sendbuf(off(dest) + 2) = 0
          sendbuf(off(dest) + 3) = nv
          sendbuf(off(dest) + 4:off(dest) + 3 + nv) = key
          off(dest) = off(dest) + 3 + nv
          deallocate (key)
          fpos = fpos + nv + 1
        end do
      end do
      do t = 1, n_pairs
        nv = fptr(t + 1) - fptr(t)
        allocate (key(nv))
        key = fdat(fptr(t):fptr(t + 1) - 1)
        call small_sort(key)
        dest = block_owner(int(key(1), ENTIER_D), num_procs, node_starts)
        sendbuf(off(dest) + 1) = 2
        sendbuf(off(dest) + 2) = bc_k(t)
        sendbuf(off(dest) + 3) = nv
        sendbuf(off(dest) + 4:off(dest) + 3 + nv) = key
        off(dest) = off(dest) + 3 + nv
        deallocate (key)
      end do
      call exchange_int(num_procs, scount, sendbuf, rcount, recvbuf)
      deallocate (sendbuf)

      ! parse into keyed entries
      nk = 0
      k = 1
      do while (k <= size(recvbuf))
        nk = nk + 1
        k = k + 3 + recvbuf(k + 2)
      end do
      allocate (kptr(nk + 1), kdat(max(size(recvbuf), 1)), ktype(nk), ksrc(nk), kval(nk), perm(nk))
      kptr(1) = 1
      k = 1
      i = 0
      do r = 0, num_procs - 1
        m = k + rcount(r)
        do while (k < m)
          i = i + 1
          ktype(i) = recvbuf(k)
          kval(i) = recvbuf(k + 1)
          ksrc(i) = r
          nv = recvbuf(k + 2)
          kdat(kptr(i):kptr(i) + nv - 1) = recvbuf(k + 3:k + 2 + nv)
          kptr(i + 1) = kptr(i) + nv
          k = k + 3 + nv
        end do
      end do
      call sort_perm_records(nk, kptr, kdat, perm)

      ! groups of identical faces
      rep_count = 0
      allocate (sendbuf(max(size(recvbuf), 1)))
      deallocate (recvbuf)
      ib = 0
      ! first pass computes the reply for each destination into per-rank
      ! segments; replies are small, so use a two-pass count/fill.
      do ib = 1, 2
        if (ib == 2) then
          deallocate (sendbuf)
          allocate (sendbuf(sum(rep_count)))
          off(0) = 0
          do r = 1, num_procs - 1
            off(r) = off(r - 1) + rep_count(r - 1)
          end do
        end if
        g0 = 1
        do while (g0 <= nk)
          g1 = g0
          do while (g1 < nk)
            if (.not. same_record(kptr, kdat, perm(g1 + 1), perm(g0))) exit
            g1 = g1 + 1
          end do
          ncell = 0
          nbcf = 0
          kbc = -1
          dest = -1
          do t = g0, g1
            if (ktype(perm(t)) == 1) then
              ncell = ncell + 1
              dest = ksrc(perm(t))
            else
              nbcf = nbcf + 1
              if (kbc < 0) then
                kbc = kval(perm(t))
              else if (kval(perm(t)) /= kbc) then
                if (ib == 1) stats(6) = stats(6) + 1
                ! deterministic choice, same as the msh path's first match:
                ! smallest nonzero tag wins
                if (kbc == 0 .or. (kval(perm(t)) /= 0 .and. kval(perm(t)) < kbc)) kbc = kval(perm(t))
              end if
            end if
          end do
          nv = kptr(perm(g0) + 1) - kptr(perm(g0))
          if (ib == 1) then
            if (ncell > 2) stats(7) = stats(7) + 1
            if (nbcf > 0 .and. ncell == 0) stats(4) = stats(4) + nbcf
            if (nbcf > 0 .and. ncell == 2) stats(3) = stats(3) + nbcf
            if (nbcf > 0 .and. ncell == 1) stats(2) = stats(2) + nbcf
            if (nbcf == 0 .and. ncell == 1) stats(5) = stats(5) + 1
            if (nbcf > 0 .and. ncell == 1) rep_count(dest) = rep_count(dest) + 2 + nv
          else if (nbcf > 0 .and. ncell == 1) then
            sendbuf(off(dest) + 1) = kbc
            sendbuf(off(dest) + 2) = nv
            sendbuf(off(dest) + 3:off(dest) + 2 + nv) = kdat(kptr(perm(g0)):kptr(perm(g0) + 1) - 1)
            off(dest) = off(dest) + 2 + nv
          end if
          g0 = g1 + 1
        end do
      end do
      call exchange_int(num_procs, rep_count, sendbuf, acount, rbuf)

      ! my boundary faces
      m = 0
      k = 1
      do while (k <= size(rbuf))
        m = m + 1
        k = k + 2 + rbuf(k + 1)
      end do
      mesh%n_boundary_faces = m
      allocate (mesh%boundary_face(m))
      k = 1
      do t = 1, m
        nv = rbuf(k + 1)
        mesh%boundary_face(t)%tag = rbuf(k)
        mesh%boundary_face(t)%n_vert = nv
        allocate (mesh%boundary_face(t)%vert(nv))
        do j = 1, nv
          mesh%boundary_face(t)%vert(j) = local_node(rbuf(k + 1 + j))
        end do
        k = k + 2 + nv
      end do

      call MPI_REDUCE(stats, gstats, 8, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, mpi_ierr)
      if (me == 0) then
        print '(a,i0,a,i0,a)', " [cgns] boundary faces: ", gstats(1), " in BCs, ", gstats(2), &
          " attached to a cell"
        if (gstats(3) > 0) print '(a,i0,a)', achar(27)//"[33m [cgns] warning: ", gstats(3), &
          " BC faces lie between two cells (internal surface), ignored"//achar(27)//"[0m"
        if (gstats(5) > 0) print '(a,i0,a)', " [cgns] ", gstats(5), &
          " exterior faces without BC -> tag 0 (e.g. z caps of extruded 2D meshes)"
        if (gstats(6) > 0) print '(a,i0,a)', achar(27)//"[33m [cgns] warning: ", gstats(6), &
          " faces appear in several BCs with different tags (smallest nonzero kept)"//achar(27)//"[0m"
      end if
      call MPI_BCAST(gstats, 8, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
      if (gstats(4) > 0) call fatal("CGNS: some BC faces match no face of any cell")
      if (gstats(7) > 0) call fatal("CGNS: non-manifold mesh (a face shared by more than two cells)")
    end subroutine match_boundary_faces

    ! Each send list p->q must equal, element by element (global ids in
    ! order), q's recv list from p.
    subroutine check_send_recv_consistency()
      integer(kind=ENTIER), dimension(0:num_procs - 1) :: scount, rcount
      integer(kind=ENTIER), dimension(:), allocatable :: sendbuf, recvbuf
      integer(kind=ENTIER) :: r, is, nbad, gbad, pos, n

      scount = 0
      do is = 1, mpi_send_recv%n_mpi_send_neigh
        scount(mpi_send_recv%mpi_send_neigh(is)%partition_id) = mpi_send_recv%mpi_send_neigh(is)%n_elems
      end do
      allocate (sendbuf(sum(scount)))
      pos = 0
      do is = 1, mpi_send_recv%n_mpi_send_neigh
        do j = 1, mpi_send_recv%mpi_send_neigh(is)%n_elems
          pos = pos + 1
          sendbuf(pos) = elem_gid(mpi_send_recv%mpi_send_neigh(is)%elem_id(j))
        end do
      end do
      call exchange_int(num_procs, scount, sendbuf, rcount, recvbuf)

      nbad = 0
      pos = 0
      do r = 0, num_procs - 1
        n = 0
        do is = 1, mpi_send_recv%n_mpi_recv_neigh
          if (mpi_send_recv%mpi_recv_neigh(is)%partition_id == r) then
            n = mpi_send_recv%mpi_recv_neigh(is)%n_elems
            if (n /= rcount(r)) then
              nbad = nbad + 1
            else
              do j = 1, n
                if (recvbuf(pos + j) /= elem_gid(mpi_send_recv%mpi_recv_neigh(is)%elem_id(j))) &
                  nbad = nbad + 1
              end do
            end if
          end if
        end do
        if (n == 0 .and. rcount(r) > 0) nbad = nbad + 1
        pos = pos + rcount(r)
      end do
      call MPI_ALLREDUCE(nbad, gbad, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
      if (gbad > 0) call fatal("dist_mesh: send/recv lists disagree between neighbours")
    end subroutine check_send_recv_consistency


    subroutine print_summary()
      integer(kind=ENTIER), dimension(4) :: loc, gmin, gmax
      integer(kind=ENTIER) :: gsum

      loc = [n_owned, held%n - n_owned, mpi_send_recv%n_mpi_send_neigh, mesh%n_vert]
      call MPI_REDUCE(loc, gmin, 4, MPI_INTEGER, MPI_MIN, 0, MPI_COMM_WORLD, mpi_ierr)
      call MPI_REDUCE(loc, gmax, 4, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, mpi_ierr)
      call MPI_REDUCE(n_owned, gsum, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, mpi_ierr)
      if (me == 0) then
        print '(a,i0,a,i0,a,a,a,i0,a)', " [cgns] ", gsum, " cells, ", cf%n_nodes, " nodes, partitioner ", &
          trim(partitioner), ", ", n_layers, " ghost layer(s)"
        if (edgecut >= 0) print '(a,i0)', " [cgns] ParMETIS edgecut: ", edgecut
        print '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)', " [cgns] owned cells min/max ", gmin(1), "/", gmax(1), &
          ", ghosts ", gmin(2), "/", gmax(2), ", neighbours ", gmin(3), "/", gmax(3)
        print '(a,f8.3,a,f8.3,a,f8.3,a,f8.3,a)', " [cgns] time read ", t_read, "s, partition ", t_part, &
          "s, ghosts ", t_ghost, "s, bc+assembly ", t_bc, "s"
      end if
    end subroutine print_summary
  end subroutine read_mesh_cgns_dist
end module dist_mesh_module
