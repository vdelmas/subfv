module acoustic_module
  use precision_module
  use mesh_module
  use mesh_reading_module
  use mesh_connectivity_module
  use mesh_geometry_module
  implicit none

  real(kind=DOUBLE), parameter :: rho = 2.0_DOUBLE, a = 1.0_DOUBLE
contains
  subroutine multi_point_pressure(soll, solr, fluxl, fluxr, smax, p_nodal) 
    implicit none

    real(kind=DOUBLE), dimension(4), intent(in) :: soll, solr
    real(kind=DOUBLE), dimension(4), intent(inout) :: fluxl, fluxr
    real(kind=DOUBLE), intent(inout) :: smax
    real(kind=DOUBLE), intent(in) :: p_nodal

    real(kind=DOUBLE) :: ul_et, ur_et

    smax = a

    ul_et = soll(2) - (p_nodal - soll(1))/(rho*a)
    ur_et = solr(2) + (p_nodal - solr(1))/(rho*a)

    fluxl(1) = rho*a**2*ul_et
    fluxl(2) = p_nodal/rho
    fluxl(3) = 0._DOUBLE
    fluxl(4) = 0._DOUBLE
    fluxl = - fluxl

    fluxr(1) = rho*a**2*ur_et
    fluxr(2) = p_nodal/rho
    fluxr(3) = 0._DOUBLE
    fluxr(4) = 0._DOUBLE
  end subroutine multi_point_pressure

  subroutine rusanov(soll, solr, fluxl, fluxr, smax) 
    implicit none

    real(kind=DOUBLE), dimension(4), intent(in) :: soll, solr
    real(kind=DOUBLE), dimension(4), intent(inout) :: fluxl, fluxr
    real(kind=DOUBLE), intent(inout) :: smax

    smax = a
    fluxr = 0.5_DOUBLE*(acoustic_flux(solr) + acoustic_flux(soll)) &
      - 0.5_DOUBLE*smax*(solr - soll)
    fluxl = - fluxr
  end subroutine rusanov

  function acoustic_flux(sol)
    implicit none

    real(kind=DOUBLE), dimension(4), intent(in) :: sol
    real(kind=DOUBLE), dimension(4) :: acoustic_flux

    acoustic_flux(1) = rho*a**2*sol(2)
    acoustic_flux(2) = sol(1)/rho
    acoustic_flux(3) = 0.0_DOUBLE
    acoustic_flux(4) = 0.0_DOUBLE
  end function acoustic_flux

  subroutine compute_error_sphere(mesh, sol)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(4, mesh%n_elems), intent(in) :: sol

    integer :: i, id_elem
    real(kind=DOUBLE) :: r, theta, phi, error, exact_sol(4), tot_volume

    error = 0.0_DOUBLE
    tot_volume = 0.0_DOUBLE
    do i=1,mesh%n_elems
      r = norm2(mesh%elem(i)%coord)
      theta = atan2(norm2(mesh%elem(i)%coord(2:3)), mesh%elem(i)%coord(1))
      phi = atan2(mesh%elem(id_elem)%coord(3), mesh%elem(id_elem)%coord(2))
      call exact_sphere(r,theta,phi,exact_sol)
      ! error = error + mesh%elem(i)%volume*norm2(sol(2:4,i)-exact_sol(2:4))**2
      error = error + mesh%elem(i)%volume*(norm2(sol(2:4,i))-norm2(exact_sol(2:4)))**2
      tot_volume = tot_volume + mesh%elem(i)%volume
    end do

    error = sqrt(error)
    print*,(tot_volume/mesh%n_elems)**(1.0_DOUBLE/3.0_DOUBLE), error
  end subroutine compute_error_sphere

  pure subroutine exact_sphere(r, theta, phi, sol)
    implicit none

    real(kind=DOUBLE), intent(in) :: r, theta, phi
    real(kind=DOUBLE), dimension(4), intent(inout) :: sol(4)

    real(kind=DOUBLE), parameter :: r0 = 0.5, r1 = 5.5
    real(kind=DOUBLE) :: vr, vtheta, A, B

    B = -1.0_DOUBLE/(2.0_DOUBLE/r0**3 + 1.0_DOUBLE/r1**3)
    A = 2.0_DOUBLE*B/r0**3

    ! vr = -(A-2.0_DOUBLE*B/r**3)*cos(theta)
    ! vtheta = (A+B/r**3)*sin(theta)

    vr = (1.0_DOUBLE-(r0/r)**3)*cos(theta)
    vtheta = -(1.0_DOUBLE+0.5_DOUBLE*(r0/r)**3)*sin(theta)

    sol(1) = 0.0_DOUBLE
    sol(2) = cos(theta)*vr - sin(theta)*vtheta
    sol(3) = (sin(theta)*vr + cos(theta)*vtheta)*cos(phi)
    sol(4) = (sin(theta)*vr + cos(theta)*vtheta)*sin(phi)
  end subroutine exact_sphere

  subroutine compute_error_cyl(mesh, sol)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(4, mesh%n_elems), intent(in) :: sol

    integer :: i
    real(kind=DOUBLE) :: r, theta, error, exact_sol(4), tot_volume

    error = 0.0_DOUBLE
    tot_volume = 0.0_DOUBLE
    do i=1,mesh%n_elems
      r = sqrt(mesh%elem(i)%coord(2)**2 + mesh%elem(i)%coord(1)**2)
      theta = atan2(mesh%elem(i)%coord(2), mesh%elem(i)%coord(1))
      call exact_cyl(r,theta,exact_sol)
      error = error + mesh%elem(i)%volume*norm2(sol(2:4,i)-exact_sol(2:4))**2
      tot_volume = tot_volume + mesh%elem(i)%volume
    end do

    error = sqrt(error)
    print*,(tot_volume/mesh%n_elems)**(1.0_DOUBLE/3.0_DOUBLE), error
  end subroutine compute_error_cyl

  pure subroutine exact_cyl(r, theta, sol)
    implicit none

    real(kind=DOUBLE), intent(in) :: r, theta
    real(kind=DOUBLE), dimension(4), intent(inout) :: sol(4)

    real(kind=DOUBLE), parameter :: r0 = 0.5, r1 = 5.5
    real(kind=DOUBLE) :: vr, vtheta

    vr = (r1**2/(r1**2-r0**2))*(1.0_DOUBLE-(r0/r)**2)*cos(theta)
    vtheta = -(r1**2/(r1**2-r0**2))*(1+(r0/r)**2)*sin(theta)

    sol(1) = 0.0_DOUBLE
    sol(2) = cos(theta)*vr - sin(theta)*vtheta
    sol(3) = sin(theta)*vr + cos(theta)*vtheta
    sol(4) = 0.0_DOUBLE
  end subroutine exact_cyl

  subroutine write_sol_vtk(mesh, sol, filename)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(4, mesh%n_elems), intent(in) :: sol
    character(len=255), intent(in) :: filename

    integer(kind=ENTIER) :: vtkout, i, j, size_tot
    real(kind=DOUBLE), dimension(3, mesh%n_vert) :: vorticity
    real(kind=DOUBLE), dimension(mesh%n_vert) :: divergence
    real(kind=DOUBLE) :: vorticity_norm, divergence_norm

    open(newunit=vtkout, file=filename)
    write (vtkout, '(A)') '# vtk DataFile Version 3.0'
    write (vtkout, '(A)') filename
    write (vtkout, '(A)') 'ASCII'
    write (vtkout, '(A)') 'DATASET UNSTRUCTURED_GRID'

    write (vtkout, *) 'POINTS ', mesh%n_vert, ' double'
    do i = 1, mesh%n_vert
      write (vtkout, *) mesh%vert(i)%coord(1), mesh%vert(i)%coord(2), &
        mesh%vert(i)%coord(3)
    end do

    size_tot = 0
    do i = 1, mesh%n_elems
      size_tot = size_tot + mesh%elem(i)%n_vert + 1
    end do

    write (vtkout, *) 'CELLS ', mesh%n_interior_elems, size_tot
    do i = 1, mesh%n_elems
      write (vtkout, '(i12)', advance='no') mesh%elem(i)%n_vert
      do j = 1, mesh%elem(i)%n_vert - 1
        write (vtkout, '(i12)', advance='no') mesh%elem(i)%vert(j) - 1
      end do
      write (vtkout, '(i12)') mesh%elem(i)%vert(mesh%elem(i)%n_vert) - 1
    end do

    write (vtkout, '(a,i11)') 'CELL_TYPES ', mesh%n_interior_elems
    do i = 1, mesh%n_elems
      if (mesh%elem(i)%elem_kind == 4) then
        write (vtkout, '(i12)') 10
      else if (mesh%elem(i)%elem_kind == 5) then
        write (vtkout, '(i12)') 12
      else if (mesh%elem(i)%elem_kind == 6) then
        write (vtkout, '(i12)') 13
      else if (mesh%elem(i)%elem_kind == 7) then
        write (vtkout, '(i12)') 14
      end if
    end do

    write (vtkout, '(a,i11)') 'CELL_DATA ', mesh%n_elems
    write (vtkout, *) 'SCALARS Pressure double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      write (vtkout, *) sol(1, i)
    end do

    write (vtkout, *) 'VECTORS Velocity double'
    do i = 1, mesh%n_elems
      write (vtkout, *) sol(2, i), sol(3, i), sol(4, i)
    end do

    call compute_operators(mesh, sol, vorticity, vorticity_norm, divergence, divergence_norm)
    write (vtkout, '(a,i11)') 'POINT_DATA ', mesh%n_vert
    write (vtkout, *) 'VECTORS Vorticity double'
    do i=1, mesh%n_vert
      write (vtkout, *) vorticity(:, i)
    end do

    write (vtkout, *) 'SCALARS Divergence double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i=1, mesh%n_vert
      write (vtkout, *) divergence(i)
    end do

    close(vtkout)
  end subroutine write_sol_vtk

  pure subroutine rotate(v, n1, forward)
    implicit none

    real(kind=DOUBLE), dimension(3), intent(inout) :: v
    real(kind=DOUBLE), dimension(3), intent(in) :: n1
    logical, intent(in) :: forward

    real(kind=DOUBLE) :: nn2
    real(kind=DOUBLE), dimension(3) :: n2, n3, vtmp

    n2 = (/n1(3), 0.0_DOUBLE, -n1(1)/)
    nn2 = norm2(n2)
    if( nn2 > 1e-10_DOUBLE ) then
      n2 = n2/nn2
    else
      n2 = (/0.0_DOUBLE, -n1(3), n1(2)/)
      n2 = n2/norm2(n2)
    end if

    n3 = cross_product(n1, n2)
    n3 = n3/norm2(n3)

    if( forward ) then
      vtmp(1) = dot_product(v, n1)
      vtmp(2) = dot_product(v, n2)
      vtmp(3) = dot_product(v, n3)
      v = vtmp
    else
      v = v(1)*n1 + v(2)*n2 + v(3)*n3
    end if
  end subroutine rotate

  pure function cross_product(a, b) result(c)
    implicit none

    real(kind=DOUBLE), dimension(3), intent(in) :: a, b
    real(kind=DOUBLE), dimension(3) :: c

    c(1) =   a(2)*b(3) - a(3)*b(2)
    c(2) = -(a(1)*b(3) - a(3)*b(1))
    c(3) =   a(1)*b(2) - a(2)*b(1)
  end function cross_product

  subroutine smooth_deform_mesh(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    integer(kind=ENTIER) :: i
    real(kind=DOUBLE), parameter :: coeff = 30.0_DOUBLE
    real(kind=DOUBLE), parameter :: xl = -0.5_DOUBLE, xs = 0._DOUBLE, xr = 0.5_DOUBLE
    real(kind=DOUBLE), parameter :: yl = -0.5_DOUBLE, ys = 0._DOUBLE, yr = 0.5_DOUBLE
    real(kind=DOUBLE), parameter :: zl = -0.5_DOUBLE, zs = 0._DOUBLE, zr = 0.5_DOUBLE
    real(kind=DOUBLE) :: x, y, z

    do i=1, mesh%n_vert
      x = mesh%vert(i)%coord(1)
      y = mesh%vert(i)%coord(2)
      z = mesh%vert(i)%coord(3)

      if( y > ys ) then
        mesh%vert(i)%coord(1) = x + (-coeff*(y-ys)*(y-yr))*(x-xl)*(x-xs)*(x-xr)
      else
        mesh%vert(i)%coord(1) = x - (-coeff*(y-yl)*(y-ys))*(x-xl)*(x-xs)*(x-xr)
      end if

      if(z > zs) then
        mesh%vert(i)%coord(2) = y + (-coeff*(z-zs)*(z-zr))*(y-yl)*(y-ys)*(y-yr)
      else 
        mesh%vert(i)%coord(2) = y - (-coeff*(z-zl)*(z-zs))*(y-yl)*(y-ys)*(y-yr)
      end if

      if(x > xs) then
        mesh%vert(i)%coord(3) = z + (-coeff*(x-xs)*(x-xr))*(z-zl)*(z-zs)*(z-zr)
      else
        mesh%vert(i)%coord(3) = z - (-coeff*(x-xl)*(x-xs))*(z-zl)*(z-zs)*(z-zr)
      end if
    end do
  end subroutine smooth_deform_mesh

  subroutine write_sol_vtk_wasilij(mesh, sol, filename)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(4, mesh%n_elems), intent(in) :: sol
    character(len=255), intent(in) :: filename

    integer(kind=ENTIER) :: i, j, size_tot, vtkout
    real(kind=DOUBLE), dimension(3, mesh%n_vert) :: vorticity
    real(kind=DOUBLE), dimension(mesh%n_vert) :: divergence
    real(kind=DOUBLE) :: vorticity_norm, divergence_norm

    open (newunit=vtkout, file=trim(filename), status="unknown")

    write (vtkout, '(A)') '# vtk DataFile Version 3.0'
    write (vtkout, '(A)') filename
    write (vtkout, '(A)') 'ASCII'
    write (vtkout, '(A)') 'DATASET UNSTRUCTURED_GRID'

    write (vtkout, *) 'POINTS ', mesh%n_vert/2, ' double'
    do i = 1, mesh%n_vert/2
      write (vtkout, *) mesh%vert(i)%coord(1), mesh%vert(i)%coord(2), &
        mesh%vert(i)%coord(3)
    end do

    size_tot = 0
    do i = 1, mesh%n_elems
      size_tot = size_tot + mesh%elem(i)%n_vert/2 + 1
    end do

    write (vtkout, *) 'CELLS ', mesh%n_interior_elems, size_tot
    do i = 1, mesh%n_elems
      write (vtkout, '(i12)', advance='no') mesh%elem(i)%n_vert/2
      do j = 1, mesh%elem(i)%n_vert/2 - 1
        write (vtkout, '(i12)', advance='no') mesh%elem(i)%vert(j) - 1
      end do
      write (vtkout, '(i12)') mesh%elem(i)%vert(mesh%elem(i)%n_vert/2) - 1
    end do

    write (vtkout, '(a,i11)') 'CELL_TYPES ', mesh%n_interior_elems
    do i = 1, mesh%n_elems
      write (vtkout, '(i12)') 7
    end do

    write (vtkout, '(a,i11)') 'CELL_DATA ', mesh%n_interior_elems
    write (vtkout, *) 'SCALARS Centroid double', 3
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      write (vtkout, *) mesh%elem(i)%coord
    end do

    write (vtkout, *) 'SCALARS Pressure double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      write (vtkout, *) sol(1, i)
    end do

    write (vtkout, *) 'VECTORS Velocity double'
    do i = 1, mesh%n_elems
      write (vtkout, *) sol(2, i), sol(3, i), sol(4, i)
    end do

    call compute_operators(mesh, sol, vorticity, vorticity_norm, divergence, divergence_norm)
    write (vtkout, '(a,i11)') 'POINT_DATA ', mesh%n_vert
    write (vtkout, *) 'VECTORS Vorticity double'
    do i = 1, mesh%n_vert
      write (vtkout, *) vorticity(:, i)
    end do

    write (vtkout, *) 'SCALARS Divergence double'
    do i = 1, mesh%n_vert
      write (vtkout, *) divergence(i)
    end do

    close (vtkout)
  end subroutine write_sol_vtk_wasilij

  subroutine compute_operators(mesh, sol, vorticity, vorticity_norm, divergence, divergence_norm)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(4, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vorticity
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(inout) :: divergence
    real(kind=DOUBLE), intent(inout) :: vorticity_norm, divergence_norm

    integer(kind=ENTIER) :: i, j, k
    integer(kind=ENTIER) :: id_sub_elem, id_elem, id_face, id_sub_face
    real(kind=DOUBLE) :: vol
    real(kind=DOUBLE), dimension(3) :: Npcf, circ

    vorticity_norm = 0.0_DOUBLE
    divergence_norm = 0.0_DOUBLE
    do i=1, mesh%n_vert
      vorticity(:, i) = 0.0_DOUBLE
      divergence(i) = 0.0_DOUBLE
      if( .not. mesh%vert(i)%is_bound ) then
        vol = 0.0_DOUBLE
        circ = 0.0_DOUBLE
        do j=1, mesh%vert(i)%n_sub_elems_neigh
          id_sub_elem = mesh%vert(i)%sub_elem_neigh(j)
          id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem

          Npcf = 0.0_DOUBLE
          do k=1, mesh%sub_elem(id_sub_elem)%n_sub_faces
            id_sub_face = mesh%sub_elem(id_sub_elem)%sub_face(k)
            id_face = mesh%sub_face(id_sub_face)%mesh_face
            if( mesh%sub_face(id_sub_face)%left_elem_neigh == id_elem ) then
              Npcf = Npcf + mesh%sub_face(id_sub_face)%area*mesh%face(id_face)%norm
            else
              Npcf = Npcf - mesh%sub_face(id_sub_face)%area*mesh%face(id_face)%norm
            end if
          end do

          vorticity(:, i) = vorticity(:, i) - cross_product(Npcf, sol(2:4, id_elem))
          divergence(i) = divergence(i) - dot_product(Npcf, sol(2:4, id_elem))
          vol = vol + mesh%sub_elem(id_sub_elem)%volume
          circ = circ + Npcf
        end do
        vorticity(:, i) = vorticity(:, i)/vol
        vorticity_norm = vorticity_norm + vol*abs(vorticity(1,i))
        vorticity_norm = vorticity_norm + vol*abs(vorticity(2,i))
        vorticity_norm = vorticity_norm + vol*abs(vorticity(3,i))

        divergence(i) = divergence(i)/vol
        divergence_norm = divergence_norm + vol*abs(divergence(i))

        if( norm2(circ) > 1e-14_DOUBLE ) then
          print*,"Bad circ"
          error stop
        end if
      end if
    end do
  end subroutine compute_operators
end module acoustic_module
