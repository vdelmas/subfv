module mesh_geometry_module
  use precision_module
  use mesh_module
  use mpi_module
  implicit none

  private

  public :: compute_geometry_mesh
  public :: map_cylinder
  public :: odd_even
  public :: project_sol
  public :: project_sol_box
  public :: mpi_project_sol_box
contains
  subroutine compute_geometry_mesh(mesh, use_sub_entities, b2d, me, num_procs)
    implicit none

    type(mesh_type), intent(inout) :: mesh
    logical, intent(in) :: use_sub_entities, b2d
    integer(kind=ENTIER), optional, intent(in) :: me, num_procs

    integer(kind=ENTIER) :: i, id_sub_elem

    do i=1, mesh%n_faces
      mesh%face(i)%coord = compute_face_centroid(mesh, i)
    end do

    do i=1, mesh%n_elems
      mesh%elem(i)%coord = compute_elem_centroid(mesh, i)
    end do

    do i=1, mesh%n_sub_faces
      mesh%sub_face(i)%norm = compute_subface_norm(mesh, i)
      mesh%sub_face(i)%area = norm2(mesh%sub_face(i)%norm)
      mesh%sub_face(i)%norm = mesh%sub_face(i)%norm / mesh%sub_face(i)%area
    end do

    do i=1, mesh%n_faces
      mesh%face(i)%norm = compute_face_norm(mesh, i)
      mesh%face(i)%area = norm2(mesh%face(i)%norm)
      mesh%face(i)%norm = mesh%face(i)%norm / mesh%face(i)%area
    end do

    do i=1, mesh%n_sub_elems
      mesh%sub_elem(i)%volume = compute_subelem_volume(mesh, i)
    end do

    do i=1, mesh%n_elems
      mesh%elem(i)%volume = compute_elem_volume(mesh, i)
    end do

    call check_mesh(mesh, me, num_procs)

    call find_if_vert_is_bound(mesh, b2d)
    call compute_vert_volume(mesh)
  end subroutine compute_geometry_mesh

  pure function compute_elem_volume(mesh, i) result(vol)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE) :: vol

    integer(kind=ENTIER) :: j, id_sub_elem

    vol = 0.0_DOUBLE
    do j=1, mesh%elem(i)%n_sub_elems
      id_sub_elem = mesh%elem(i)%sub_elem(j)
      vol = vol + mesh%sub_elem(id_sub_elem)%volume
    end do
  end function compute_elem_volume

  pure function compute_elem_centroid(mesh, i) result(centroid)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE), dimension(3) :: centroid

    integer(kind=ENTIER) :: j, k, kp
    integer(kind=ENTIER) :: id_face
    integer(kind=ENTIER) :: idv, idvp
    real(kind=DOUBLE), dimension(3) :: cp, cpp, ce, cf
    real(kind=DOUBLE) :: vol, sum_vol

    ce = compute_elem_cgeom(mesh, i)
    centroid = 0.0_DOUBLE
    sum_vol = 0.0_DOUBLE
    do j=1, mesh%elem(i)%n_faces
      id_face = mesh%elem(i)%face(j)
      cf = mesh%face(id_face)%coord
      do k = 1, mesh%face(id_face)%n_vert
        kp = 1+mod((k-1)+1,mesh%face(id_face)%n_vert)
        idv = mesh%face(id_face)%vert(k)
        idvp = mesh%face(id_face)%vert(kp)
        cp = mesh%vert(idv)%coord
        cpp = mesh%vert(idvp)%coord
        if( mesh%face(id_face)%left_neigh == i ) then
          vol = dot_product(cross_product(cp-cf, cpp-cf), cf-ce)/6.0_DOUBLE
        else
          vol = dot_product(cross_product(cpp-cf, cp-cf), cf-ce)/6.0_DOUBLE
        end if
        centroid = centroid + vol*(cp+cpp+cf+ce)/4.0_DOUBLE
        sum_vol = sum_vol + vol
      end do
    end do
    centroid = centroid / sum_vol
  end function compute_elem_centroid

  pure function compute_elem_cgeom(mesh, i) result(cgeom)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE), dimension(3) :: cgeom

    integer(kind=ENTIER) :: j, id_vert

    cgeom = 0.0_DOUBLE
    do j=1, mesh%elem(i)%n_vert
      id_vert = mesh%elem(i)%vert(j)
      cgeom = cgeom + mesh%vert(id_vert)%coord
    end do
    cgeom = cgeom / mesh%elem(i)%n_vert
  end function compute_elem_cgeom

  pure function compute_subelem_volume(mesh, i) result(vol)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE) :: vol

    integer(kind=ENTIER) :: j, k, kp, km, idvp, idv, idvm
    integer(kind=ENTIER) :: id_sub_face, id_vert, id_face, id_elem
    real(kind=DOUBLE), dimension(3) :: cf, cp, cpm, cpp, ce

    id_elem = mesh%sub_elem(i)%mesh_elem
    ce = mesh%elem(id_elem)%coord

    vol = 0.0_DOUBLE
    do j=1, mesh%sub_elem(i)%n_sub_faces
      id_sub_face = mesh%sub_elem(i)%sub_face(j)
      id_vert = mesh%sub_face(id_sub_face)%mesh_vert
      id_face = mesh%sub_face(id_sub_face)%mesh_face

      cf = mesh%face(id_face)%coord
      do k=1, mesh%face(id_face)%n_vert
        idv = mesh%face(id_face)%vert(k)
        if( idv == id_vert ) then
          kp = 1+mod((k-1)+1, mesh%face(id_face)%n_vert)
          km = 1+mod((k-1)-1+mesh%face(id_face)%n_vert,&
            mesh%face(id_face)%n_vert)
          exit
        end if
      end do
      idvp = mesh%face(id_face)%vert(kp)
      idvm = mesh%face(id_face)%vert(km)
      cp = mesh%vert(idv)%coord
      cpm = 0.5_DOUBLE*(mesh%vert(idvm)%coord+cp)
      cpp = 0.5_DOUBLE*(mesh%vert(idvp)%coord+cp)

      if( mesh%face(id_face)%left_neigh == id_elem ) then
        vol = vol &
          + (dot_product(cross_product(cpm-cf, cp-cf), cf-ce))/6.0_DOUBLE&
          + (dot_product(cross_product(cp-cf, cpp-cf), cf-ce))/6.0_DOUBLE
      else
        vol = vol &
          + (dot_product(cross_product(cpp-cf, cp-cf), cf-ce))/6.0_DOUBLE&
          + (dot_product(cross_product(cp-cf, cpm-cf), cf-ce))/6.0_DOUBLE
      end if
    end do
  end function compute_subelem_volume

  pure function compute_face_norm(mesh, i) result(norm)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE), dimension(3) :: norm

    integer(kind=ENTIER) :: j, id_sub_face

    norm = 0.0_DOUBLE
    do j=1, mesh%face(i)%n_sub_faces
      id_sub_face = mesh%face(i)%sub_face(j)
      norm = norm + compute_subface_norm(mesh, id_sub_face)
    end do
  end function compute_face_norm

  pure function compute_face_cgeom(mesh, i) result(cgeom)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE), dimension(3) :: cgeom 

    integer(kind=ENTIER) :: j, id_vert

    cgeom = 0.0_DOUBLE
    do j=1, mesh%face(i)%n_vert
      id_vert = mesh%face(i)%vert(j)
      cgeom = cgeom + mesh%vert(id_vert)%coord
    end do
    cgeom = cgeom / mesh%face(i)%n_vert
  end function compute_face_cgeom

  pure function compute_face_centroid(mesh, i) result(centroid)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE), dimension(3) :: centroid

    integer(kind=ENTIER) :: j, jp, idv, idvp
    real(kind=DOUBLE), dimension(3) :: cgeom, cv, cvp, pdv
    real(kind=DOUBLE) :: area, sum_area

    cgeom = compute_face_cgeom(mesh, i)
    centroid = 0.0_DOUBLE
    sum_area = 0.0_DOUBLE
    do j=1, mesh%face(i)%n_vert
      jp = 1+mod((j-1)+1, mesh%face(i)%n_vert)
      idv = mesh%face(i)%vert(j)
      idvp = mesh%face(i)%vert(jp)
      cv = mesh%vert(idv)%coord
      cvp = mesh%vert(idvp)%coord
      area = 0.5_DOUBLE*norm2(cross_product(cv-cgeom, cvp-cgeom))
      centroid = centroid + area*(cv+cvp+cgeom)/3.0_DOUBLE
      sum_area = sum_area + area
    end do
    centroid = centroid / sum_area
  end function compute_face_centroid

  pure function compute_subface_norm(mesh, i) result(norm)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE), dimension(3) :: norm

    integer(kind=ENTIER) :: id_face, id_vert
    integer(kind=ENTIER) :: j, jp, jm, idvp, idv, idvm

    real(kind=DOUBLE), dimension(3) :: cpm, cp, cpp, cf

    id_vert = mesh%sub_face(i)%mesh_vert
    id_face = mesh%sub_face(i)%mesh_face

    do j=1, mesh%face(id_face)%n_vert
      idv = mesh%face(id_face)%vert(j)
      if( idv == id_vert ) then
        jp = 1+mod((j-1)+1, mesh%face(id_face)%n_vert)
        jm = 1+mod((j-1)-1+mesh%face(id_face)%n_vert,&
          mesh%face(id_face)%n_vert)
        exit
      end if
    end do

    idvp = mesh%face(id_face)%vert(jp)
    idvm = mesh%face(id_face)%vert(jm)
    cf = mesh%face(id_face)%coord
    cp = mesh%vert(idv)%coord
    cpm = 0.5_DOUBLE*(mesh%vert(idvm)%coord+cp)
    cpp = 0.5_DOUBLE*(mesh%vert(idvp)%coord+cp)

    norm = 0.5_DOUBLE*(cross_product(cpm-cf, cp-cf) &
      + cross_product(cp-cf, cpp-cf))
  end function compute_subface_norm

  subroutine compute_vert_volume(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    integer(kind=ENTIER) :: i, j, id_sub_elem

    do i=1, mesh%n_vert
      mesh%vert(i)%volume = 0.0_DOUBLE
      do j=1, mesh%vert(i)%n_sub_elems_neigh
        id_sub_elem = mesh%vert(i)%sub_elem_neigh(j)
        mesh%vert(i)%volume = mesh%vert(i)%volume &
          + mesh%sub_elem(id_sub_elem)%volume
      end do
    end do
  end subroutine compute_vert_volume

  function check_mesh_elem_volume(mesh, i) result(ok)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: i
    logical :: ok

    integer(kind=ENTIER) :: j, id_face
    real(kind=DOUBLE) :: vol
    real(kind=DOUBLE), dimension(3) :: norm

    vol = 0.0_DOUBLE
    do j=1, mesh%elem(i)%n_faces
      id_face = mesh%elem(i)%face(j)
      if( mesh%face(id_face)%left_neigh == i ) then
        norm = mesh%face(id_face)%norm
      else
        norm = -mesh%face(id_face)%norm
      end if
      vol = vol + dot_product(mesh%face(id_face)%coord, &
        mesh%face(id_face)%area*norm)
    end do
    vol = vol / 3.0_DOUBLE

    ok = abs(mesh%elem(i)%volume - vol) < 1e-14_DOUBLE
    if( .not. ok ) print*, "check_mesh_elem_volume", &
      abs(mesh%elem(i)%volume - vol)
  end function check_mesh_elem_volume

  function check_mesh_elem_closed(mesh, i) result(ok)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: i
    logical :: ok

    integer(kind=ENTIER) :: j, id_face
    real(kind=DOUBLE), dimension(3) :: norm, norm_tot

    norm_tot = 0.0_DOUBLE
    do j=1, mesh%elem(i)%n_faces
      id_face = mesh%elem(i)%face(j)
      if( mesh%face(id_face)%left_neigh == i ) then
        norm = mesh%face(id_face)%norm
      else
        norm = -mesh%face(id_face)%norm
      end if
      norm_tot = norm_tot + mesh%face(id_face)%area * norm
    end do

    ok = norm2(norm_tot) < 1e-14_DOUBLE
    if( .not. ok ) print*, "check_mesh_elem_closed", &
      norm2(norm_tot)
  end function check_mesh_elem_closed

  function check_mesh_elem_subface_closed(mesh, i) result(ok)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: i
    logical :: ok

    integer(kind=ENTIER) :: j, id_face, k, id_sub_face
    real(kind=DOUBLE), dimension(3) :: norm, norm_tot

    norm_tot = 0.0_DOUBLE
    do j=1, mesh%elem(i)%n_faces
      id_face = mesh%elem(i)%face(j)
      do k=1, mesh%face(id_face)%n_sub_faces
        id_sub_face = mesh%face(id_face)%sub_face(k)
        if( mesh%face(id_face)%left_neigh == i ) then
          norm = mesh%sub_face(id_sub_face)%norm
        else
          norm = -mesh%sub_face(id_sub_face)%norm
        end if
        norm_tot = norm_tot + mesh%sub_face(id_sub_face)%area * norm
      end do
    end do

    ok = norm2(norm_tot) < 1e-14_DOUBLE
    if( .not. ok ) print*, "check_mesh_elem_subface_closed", &
      norm2(norm_tot)
  end function check_mesh_elem_subface_closed

  subroutine check_mesh(mesh, me, num_procs)
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), optional, intent(in) :: me, num_procs

    integer(kind=ENTIER) :: i, me_, num_procs_, mpi_ierr
    logical :: ok

    me_        = 0; if (present(me))        me_        = me
    num_procs_ = 1; if (present(num_procs)) num_procs_ = num_procs

    ok = .true.
    do i=1, mesh%n_elems
      ok = (ok .and. check_mesh_elem_volume(mesh, i))
      ok = (ok .and. check_mesh_elem_closed(mesh, i))
      ok = (ok .and. check_mesh_elem_subface_closed(mesh, i))
    end do

    do i=1, mesh%n_sub_elems
      ok = (ok .and. (mesh%sub_elem(i)%volume > 0.0_DOUBLE))
    end do

    if (num_procs_ > 1) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, ok, 1, MPI_LOGICAL, MPI_LAND, &
        MPI_COMM_WORLD, mpi_ierr)
    end if

    if (.not. ok) then
      if (me_ == 0) print*, "MESH GEOMETRY IS NOT OK, &
        DO NOT TRY TO CIRCUMVENT THIS WARNING &
        WITH ANY KIND OF HOPE THAT THE CODE WILL WORK !"
      error stop
    else
      if (me_ == 0) print*, "MESH GEOMETRY IS OK"
    end if
  end subroutine check_mesh

  pure subroutine map_cylinder(mesh, dimen)
    implicit none

    type(mesh_type), intent(inout) :: mesh
    integer, intent(in) :: dimen

    integer(kind=ENTIER) :: i

    do i = 1, mesh%n_vert
      mesh%vert(i)%coord = cylinder_mapping(mesh%vert(i)%coord, dimen)
    end do
  end subroutine map_cylinder

  subroutine odd_even(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    integer(kind=ENTIER) :: i
    real(kind=DOUBLE), parameter :: PI = 4.0_DOUBLE*datan(1.0_DOUBLE)
    real(kind=DOUBLE) :: l, theta
    ! l=1e-9
    l=1e-9_DOUBLE

    do i = 1, mesh%n_vert
      if( abs(mesh%vert(i)%coord(2)) < 1e-2_DOUBLE .and. &
        abs(mesh%vert(i)%coord(3)) < 1e-2_DOUBLE) then
        theta = mesh%vert(i)%coord(1)*pi/2.0_DOUBLE
        ! print*,mesh%vert(i)%coord
        mesh%vert(i)%coord(2) = mesh%vert(i)%coord(2) + l*cos(theta)
        mesh%vert(i)%coord(3) = mesh%vert(i)%coord(3) + l*sin(theta)
        ! print*, l*cos(theta), l*sin(theta)
        ! print*,mesh%vert(i)%coord
        ! print*,"-------"
      end if
    end do
  end subroutine odd_even

  pure function cylinder_mapping(coord, dimen) result(new_coord)
    implicit none

    integer, intent(in) :: dimen
    real(kind=DOUBLE), dimension(3), intent(in) :: coord
    real(kind=DOUBLE), dimension(3) :: new_coord

    real(kind=DOUBLE) :: r, new_r, theta, d
    real(kind=DOUBLE) :: r1, new_r1, r2, new_r2

    if (dimen == 2) then
      r = norm2(coord(:2))
      theta = atan2(coord(2), -coord(1))

      r1 = 1.0
      new_r1 = r1

      r2 = 1.3851721918740656_DOUBLE
      new_r2 = r2 &
        + 0.22082080765324227_DOUBLE*theta**2 &
        + 0.07247276314520341_DOUBLE*theta**4 &
        - 0.01622444200319442_DOUBLE*theta**6 &
        + 0.00780545523147566_DOUBLE*theta**8

      new_r = (r - r1)*(new_r2 - new_r1)/(r2 - r1) + new_r1

      new_coord(1) = -new_r*cos(theta)
      new_coord(2) = new_r*sin(theta)
      new_coord(3) = coord(3)
    else if (dimen == 21) then !For cylinder with radius 0.5
      r = norm2(2.0_DOUBLE*coord(:2))
      theta = atan2(2.0_DOUBLE*coord(2), -2.0_DOUBLE*coord(1))

      r1 = 1.0
      new_r1 = r1

      r2 = 1.3851721918740656_DOUBLE
      new_r2 = r2 &
        + 0.22082080765324227_DOUBLE*theta**2 &
        + 0.07247276314520341_DOUBLE*theta**4 &
        - 0.01622444200319442_DOUBLE*theta**6 &
        + 0.00780545523147566_DOUBLE*theta**8

      new_r = (r - r1)*(new_r2 - new_r1)/(r2 - r1) + new_r1

      new_coord(1) = -new_r*cos(theta)/2.0_DOUBLE
      new_coord(2) = new_r*sin(theta)/2.0_DOUBLE
      new_coord(3) = coord(3)
    else if (dimen == 3) then
      r = norm2(coord)
      d = norm2(coord(2:))
      theta = atan2(d, -coord(1))

      r1 = 1.0
      new_r1 = r1

      r2 = 1.3851721918740656_DOUBLE
      new_r2 = r2 &
        + 0.22082080765324227_DOUBLE*theta**2 &
        + 0.07247276314520341_DOUBLE*theta**4 &
        - 0.01622444200319442_DOUBLE*theta**6 &
        + 0.00780545523147566_DOUBLE*theta**8

      new_r = (r - r1)*(new_r2 - new_r1)/(r2 - r1) + new_r1
      new_coord = coord + (new_r - r)*coord
    end if

  end function cylinder_mapping

  pure function tensor_product(a, b)
    implicit none

    real(kind=DOUBLE), intent(in) :: a(:), b(:)
    real(kind=DOUBLE), dimension(size(a), size(a)) :: tensor_product

    integer(kind=ENTIER) :: i, j

    do i = 1, size(a)
      do j = 1, size(a)
        tensor_product(i, j) = a(i)*b(j)
      end do
    end do
  end function tensor_product

  pure function cross_product(a, b)
    implicit none

    real(kind=DOUBLE), intent(in) :: a(3), b(3)
    real(kind=DOUBLE), dimension(3) :: cross_product

    cross_product(1) =   a(2)*b(3) - a(3)*b(2)
    cross_product(2) = -(a(1)*b(3) - a(3)*b(1))
    cross_product(3) =   a(1)*b(2) - a(2)*b(1)
  end function cross_product

  subroutine project_sol_box(sol_size, mesh1, sol1, mesh2, sol2)
    implicit none

    integer(kind=ENTIER) :: sol_size
    type(mesh_type), intent(in) :: mesh1, mesh2
    real(kind=DOUBLE), dimension(sol_size, mesh1%n_elems), intent(in) :: sol1
    real(kind=DOUBLE), dimension(sol_size, mesh2%n_elems), intent(out) :: sol2

    type :: box_type
      integer(kind=ENTIER), dimension(:), allocatable :: elem
    end type box_type

    real(kind=DOUBLE), parameter :: r_per_box=3.0_DOUBLE

    integer(kind=ENTIER) :: i, nx, ny, nz, ix, iy, iz, iclosest, j
    integer(kind=ENTIER) :: kx, ky, kz, idbx, idby, idbz
    real(kind=DOUBLE) :: xmin, xmax, ymin, ymax, zmin, zmax
    real(kind=DOUBLE) :: r, tot_vol, dx, dy, dz, d, d1
    integer(kind=ENTIER), dimension(:,:,:), allocatable :: n_box
    type(box_type), dimension(:,:,:), allocatable :: box

    xmin = mesh1%elem(1)%coord(1)
    xmax = mesh1%elem(1)%coord(1)
    ymin = mesh1%elem(1)%coord(2)
    ymax = mesh1%elem(1)%coord(2)
    zmin = mesh1%elem(1)%coord(3)
    zmax = mesh1%elem(1)%coord(3)

    tot_vol = 0.0_DOUBLE
    do i = 1, mesh1%n_elems
      if(mesh1%elem(i)%coord(1) < xmin) xmin = mesh1%elem(i)%coord(1)
      if(mesh1%elem(i)%coord(1) > xmax) xmax = mesh1%elem(i)%coord(1)
      if(mesh1%elem(i)%coord(2) < ymin) ymin = mesh1%elem(i)%coord(2)
      if(mesh1%elem(i)%coord(2) > ymax) ymax = mesh1%elem(i)%coord(2)
      if(mesh1%elem(i)%coord(3) < zmin) zmin = mesh1%elem(i)%coord(3)
      if(mesh1%elem(i)%coord(3) > zmax) zmax = mesh1%elem(i)%coord(3)
      tot_vol = tot_vol + mesh1%elem(i)%volume
    end do

    r = (tot_vol/mesh1%n_elems)**(1.0_DOUBLE/3.0_DOUBLE)

    xmin = xmin - 1e-1_DOUBLE*r
    xmax = xmax + 1e-1_DOUBLE*r
    ymin = ymin - 1e-1_DOUBLE*r
    ymax = ymax + 1e-1_DOUBLE*r
    zmin = zmin - 1e-1_DOUBLE*r
    zmax = zmax + 1e-1_DOUBLE*r

    if( xmax - xmin < r ) then
      dx = xmax - xmin
    else
      dx = r_per_box*r
    end if
    nx = ceiling( (xmax - xmin)/dx )

    if( ymax - ymin < r ) then
      dy = ymax - ymin
    else
      dy = r_per_box*r
    end if
    ny = ceiling( (ymax - ymin)/dy )

    if( zmax - zmin < r ) then
      dz = zmax - zmin
    else
      dz = r_per_box*r
    end if
    nz = ceiling( (zmax - zmin)/dz )

    allocate(n_box(nx,ny,nz))
    n_box = 0
    do i=1,mesh1%n_elems
      ix = ceiling( (mesh1%elem(i)%coord(1) - xmin)/dx )
      iy = ceiling( (mesh1%elem(i)%coord(2) - ymin)/dy )
      iz = ceiling( (mesh1%elem(i)%coord(3) - zmin)/dz )
      n_box(ix,iy,iz) = n_box(ix,iy,iz) + 1
    end do

    allocate(box(nx,ny,nz))
    do ix=1,nx
      do iy=1,ny
        do iz=1,nz
          if(n_box(ix,iy,iz) > 0) then
            allocate(box(ix,iy,iz)%elem(n_box(ix,iy,iz)))
            box(ix,iy,iz)%elem = 0
          end if
        end do
      end do
    end do

    n_box = 0
    do i=1,mesh1%n_elems
      ix = ceiling( (mesh1%elem(i)%coord(1) - xmin)/dx )
      iy = ceiling( (mesh1%elem(i)%coord(2) - ymin)/dy )
      iz = ceiling( (mesh1%elem(i)%coord(3) - zmin)/dz )
      n_box(ix,iy,iz) = n_box(ix,iy,iz) + 1
      box(ix,iy,iz)%elem(n_box(ix,iy,iz)) = i
    enddo

    do i=1,mesh2%n_elems
      ix = ceiling( (mesh2%elem(i)%coord(1) - xmin)/dx )
      iy = ceiling( (mesh2%elem(i)%coord(2) - ymin)/dy )
      iz = ceiling( (mesh2%elem(i)%coord(3) - zmin)/dz )

      d = 10.0_DOUBLE*r
      iclosest = 0
      do kx=1,3
        do ky=1,3
          do kz=1,3
            idbx = min(max(1,ix+kx-2),nx)
            idby = min(max(1,iy+ky-2),ny)
            idbz = min(max(1,iz+kz-2),nz)
            do j=1,n_box(idbx, idby, idbz)
              d1 = norm2(mesh2%elem(i)%coord&
                -mesh1%elem(box(idbx,idby,idbz)%elem(j))%coord)
              if( d1 < d ) then
                iclosest = box(idbx,idby,idbz)%elem(j)
                d = d1
              end if
            end do
          end do
        end do
      end do

      if(iclosest == 0) then
        print*, "Issue to find closet elems in box"
        error stop
      else
        sol2(:, i) = sol1(:, iclosest)
      end if
    end do
  end subroutine project_sol_box

  subroutine project_sol(sol_size, mesh1, sol1, mesh2, sol2)
    implicit none

    integer(kind=ENTIER) :: sol_size
    type(mesh_type), intent(in) :: mesh1, mesh2
    real(kind=DOUBLE), dimension(sol_size, mesh1%n_elems), intent(in) :: sol1
    real(kind=DOUBLE), dimension(sol_size, mesh2%n_elems), intent(out) :: sol2

    integer(kind=ENTIER) :: i1, i2, iclosest, itest, k

    i1 = 1
    do i2 = 1, mesh2%n_elems
      !Find i1 in mesh1 closest to i2 in mesh2
      do while (.TRUE.)
        iclosest = i1
        do k = 1, mesh1%elem(i1)%n_neigh_by_vert
          itest = mesh1%elem(i1)%neigh_by_vert(k)
          if (norm2(mesh1%elem(itest)%coord - mesh2%elem(i2)%coord) < &
            (norm2(mesh1%elem(iclosest)%coord - mesh2%elem(i2)%coord))) then
            iclosest = itest
          end if
        end do

        if (i1 == iclosest) then
          exit
        else
          i1 = iclosest
        end if
      end do

      !Interpolate solution
      sol2(:, i2) = sol1(:, i1)
    end do
  end subroutine project_sol

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

  subroutine mpi_project_sol_box(sol_size, mesh1, sol1, mesh2, sol2)
    use mpi
    implicit none

    integer(kind=ENTIER) :: sol_size
    type(mesh_type), intent(in) :: mesh1, mesh2
    real(kind=DOUBLE), dimension(sol_size, mesh1%n_elems), intent(in) :: sol1
    real(kind=DOUBLE), dimension(sol_size, mesh2%n_elems), intent(inout) :: sol2

    type :: box_type
      integer(kind=ENTIER), dimension(:), allocatable :: elem
      real(kind=DOUBLE), dimension(:,:), allocatable :: coord
      real(kind=DOUBLE), dimension(:,:), allocatable :: sol
    end type box_type

    real(kind=DOUBLE), parameter :: r_per_box = 5.0_DOUBLE

    integer(kind=ENTIER) :: me, num_procs, p, mpi_ierr
    integer(kind=ENTIER) :: i, nx, ny, nz, ix, iy, iz, iclosest, j
    integer(kind=ENTIER) :: kx, ky, kz, idbx, idby, idbz
    real(kind=DOUBLE) :: xmin, xmax, ymin, ymax, zmin, zmax
    real(kind=DOUBLE) :: r, tot_vol, dx, dy, dz, d1
    real(kind=DOUBLE) :: dx_current, dy_current, dz_current
    integer(kind=ENTIER), dimension(:, :, :), allocatable :: n_box
    type(box_type), dimension(:, :, :), allocatable :: box
    real(kind=double), dimension(mesh2%n_elems) :: d_closest

    real(kind=DOUBLE) :: t1, t2
    real(kind=DOUBLE) :: xmin_current, xmax_current, &
      ymin_current, ymax_current, &
      zmin_current, zmax_current
    integer(kind=ENTIER) :: nx_current, ny_current, nz_current
    integer(kind=ENTIER), dimension(:, :, :), allocatable :: n_box_current
    type(box_type), dimension(:, :, :), allocatable :: box_current

    real(kind=double), dimension(3, 5, mesh2%n_elems) :: grad
    real(kind=double), dimension(mesh2%n_elems) :: residu
    integer(kind=ENTIER), dimension(mesh2%n_elems) :: color
    character(len=255) :: me_str

    logical :: found

    sol2 = -1.
    grad = 0.
    residu = 0.
    color = 0

    call cpu_time(t1)

    xmin = mesh1%elem(1)%coord(1)
    xmax = mesh1%elem(1)%coord(1)
    ymin = mesh1%elem(1)%coord(2)
    ymax = mesh1%elem(1)%coord(2)
    zmin = mesh1%elem(1)%coord(3)
    zmax = mesh1%elem(1)%coord(3)

    tot_vol = 0.0_DOUBLE
    do i = 1, mesh1%n_elems
      if (mesh1%elem(i)%coord(1) < xmin) xmin = mesh1%elem(i)%coord(1)
      if (mesh1%elem(i)%coord(1) > xmax) xmax = mesh1%elem(i)%coord(1)
      if (mesh1%elem(i)%coord(2) < ymin) ymin = mesh1%elem(i)%coord(2)
      if (mesh1%elem(i)%coord(2) > ymax) ymax = mesh1%elem(i)%coord(2)
      if (mesh1%elem(i)%coord(3) < zmin) zmin = mesh1%elem(i)%coord(3)
      if (mesh1%elem(i)%coord(3) > zmax) zmax = mesh1%elem(i)%coord(3)
      tot_vol = tot_vol + mesh1%elem(i)%volume
    end do

    r = (tot_vol/mesh1%n_elems)**(1.0_DOUBLE/3.0_DOUBLE)

    xmin = xmin - r
    xmax = xmax + r
    ymin = ymin - r
    ymax = ymax + r
    zmin = zmin - r
    zmax = zmax + r

    if (xmax - xmin < r) then
      dx = xmax - xmin
    else
      dx = r_per_box*r
    end if
    nx = ceiling((xmax - xmin)/dx)+1

    if (ymax - ymin < r) then
      dy = ymax - ymin
    else
      dy = r_per_box*r
    end if
    ny = ceiling((ymax - ymin)/dy)+1

    if (zmax - zmin < r) then
      dz = zmax - zmin
    else
      dz = r_per_box*r
    end if
    nz = ceiling((zmax - zmin)/dz)+1

    allocate (n_box(nx, ny, nz))
    n_box = 0
    do i = 1, mesh1%n_elems
      ix = ceiling((mesh1%elem(i)%coord(1) - xmin)/dx)
      iy = ceiling((mesh1%elem(i)%coord(2) - ymin)/dy)
      iz = ceiling((mesh1%elem(i)%coord(3) - zmin)/dz)
      n_box(ix, iy, iz) = n_box(ix, iy, iz) + 1
    end do

    allocate (box(nx, ny, nz))
    do ix = 1, nx
      do iy = 1, ny
        do iz = 1, nz
          if (n_box(ix, iy, iz) > 0) then
            allocate (box(ix, iy, iz)%elem(n_box(ix, iy, iz)))
            box(ix, iy, iz)%elem = 0
            allocate (box(ix, iy, iz)%sol(sol_size, n_box(ix, iy, iz)))
            box(ix, iy, iz)%sol = 0
            allocate (box(ix, iy, iz)%coord(3, n_box(ix, iy, iz)))
            box(ix, iy, iz)%coord = 0
          end if
        end do
      end do
    end do

    n_box = 0
    do i = 1, mesh1%n_elems
      ix = ceiling((mesh1%elem(i)%coord(1) - xmin)/dx)
      iy = ceiling((mesh1%elem(i)%coord(2) - ymin)/dy)
      iz = ceiling((mesh1%elem(i)%coord(3) - zmin)/dz)
      n_box(ix, iy, iz) = n_box(ix, iy, iz) + 1
      box(ix, iy, iz)%elem(n_box(ix, iy, iz)) = i
      box(ix, iy, iz)%sol(:, n_box(ix, iy, iz)) = sol1(:, i)
      box(ix, iy, iz)%coord(:, n_box(ix, iy, iz)) = mesh1%elem(i)%coord
    end do

    !MPI
    sol2 = 0.0_DOUBLE
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    d_closest = 1e10

    do p=0, num_procs-1
      if( p == me ) then
        print*, "Proj from", me
        xmin_current = xmin
        xmax_current = xmax
        ymin_current = ymin
        ymax_current = ymax
        zmin_current = zmin
        zmax_current = zmax
        nx_current = nx
        ny_current = ny
        nz_current = nz
        dx_current = dx
        dy_current = dy
        dz_current = dz
      end if

      call mpi_bcast(xmin_current, 1, MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(xmax_current, 1, MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(ymin_current, 1, MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(ymax_current, 1, MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(zmin_current, 1, MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(zmax_current, 1, MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(dx_current, 1, MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(dy_current, 1, MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(dz_current, 1, MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(nx_current, 1, MPI_INT, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(ny_current, 1, MPI_INT, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(nz_current, 1, MPI_INT, p, MPI_COMM_WORLD, mpi_ierr)
      allocate (n_box_current(nx_current, ny_current, nz_current))
      if( me == p ) n_box_current = n_box
      call mpi_bcast(n_box_current, nx_current*ny_current*nz_current, &
        MPI_INT, p, MPI_COMM_WORLD, mpi_ierr)

      allocate (box_current(nx_current, ny_current, nz_current))
      do ix = 1, nx_current
        do iy = 1, ny_current
          do iz = 1, nz_current
            if (n_box_current(ix, iy, iz) > 0) then
              allocate (box_current(ix, iy, iz)%elem(n_box_current(ix, iy, iz)))
              if( me == p ) box_current(ix, iy, iz)%elem = box(ix, iy, iz)%elem
              call mpi_bcast(box_current(ix,iy,iz)%elem, n_box_current(ix,iy,iz), &
                MPI_INT, p, MPI_COMM_WORLD, mpi_ierr)
              allocate (box_current(ix, iy, iz)%sol(sol_size, n_box_current(ix, iy, iz)))
              if( me == p ) box_current(ix, iy, iz)%sol = box(ix, iy, iz)%sol
              call mpi_bcast(box_current(ix,iy,iz)%sol, sol_size*n_box_current(ix,iy,iz), &
                MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
              allocate (box_current(ix, iy, iz)%coord(3, n_box_current(ix, iy, iz)))
              if( me == p ) box_current(ix, iy, iz)%coord = box(ix, iy, iz)%coord
              call mpi_bcast(box_current(ix,iy,iz)%coord, 3*n_box_current(ix,iy,iz), &
                MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
            end if
          end do
        end do
      end do

      do i = 1, mesh2%n_elems
        if( mesh2%elem(i)%coord(1) > xmin_current .and. &
          mesh2%elem(i)%coord(1) < xmax_current .and.&
          mesh2%elem(i)%coord(2) > ymin_current .and. &
          mesh2%elem(i)%coord(2) < ymax_current .and.&
          mesh2%elem(i)%coord(3) > zmin_current .and. &
          mesh2%elem(i)%coord(3) < zmax_current ) then

          ix = ceiling((mesh2%elem(i)%coord(1) - xmin_current)/dx_current)
          iy = ceiling((mesh2%elem(i)%coord(2) - ymin_current)/dy_current)
          iz = ceiling((mesh2%elem(i)%coord(3) - zmin_current)/dz_current)

          found = .false.
          iclosest = 0
          do kx = 1, 5
            do ky = 1, 5
              do kz = 1, 5
                idbx = min(max(1, ix + kx - 3), nx_current)
                idby = min(max(1, iy + ky - 3), ny_current)
                idbz = min(max(1, iz + kz - 3), nz_current)
                do j = 1, n_box_current(idbx, idby, idbz)
                  d1 = norm2(mesh2%elem(i)%coord &
                    - box_current(idbx, idby, idbz)%coord(:, j))
                  if (d1 < d_closest(i)) then
                    iclosest = box_current(idbx, idby, idbz)%elem(j)
                    sol2(:, i) = box_current(idbx, idby, idbz)%sol(:, j)
                    d_closest(i) = d1
                    found = .true.
                  end if
                end do
              end do
            end do
          end do

        end if !bounding box
      end do
      deallocate (n_box_current, box_current)
      call mpi_barrier(mpi_comm_world, mpi_ierr)
    end do

    call mpi_barrier(mpi_comm_world, mpi_ierr)
    call cpu_time(t2)
    if( me == 0 ) then
      print *, ""//achar(27)//"[33m[*] Time for projection onto next mesh :"//achar(27)//"[0m", t2 - t1
    end if
  end subroutine mpi_project_sol_box
end module mesh_geometry_module
