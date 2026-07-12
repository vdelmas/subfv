program main
  use mpi
  use precision_module
  use mpi_module
  use mesh_module
  use mesh_reading_module
  use mesh_connectivity_module
  use mesh_geometry_module
  use acoustic_module
  implicit none

  real(kind=DOUBLE), parameter :: PI = 4.0_DOUBLE*datan(1.0_DOUBLE)

  type(mesh_type) :: mesh
  integer(kind=ENTIER) :: n_bc, reduced_neigh, funit
  character(len=255) :: meshfile_path, meshfile, iterchar, output
  character(len=255), dimension(10) :: bc_name, bc_type
  real(kind=DOUBLE), dimension(4, 10) :: bc_val
  logical :: use_sub_entities

  real(kind=DOUBLE), parameter :: gamma = 1.4_DOUBLE
  real(kind=DOUBLE) :: cfl = 0.2_DOUBLE
  real(kind=DOUBLE), dimension(4) :: sol_uniform
  logical :: init_uniform = .false.
  logical :: init_gresho = .false.
  logical :: init_vortex_ring = .false.
  real(kind=DOUBLE) :: distXY, largeR
  logical :: error_cyl = .false., error_sphere = .false.
  logical :: smooth_deform = .false.
  logical :: mesh_was = .false.
  logical :: write_operators = .false.
  integer(kind=ENTIER) :: fno
  real(kind=DOUBLE) :: vorticity_norm, divergence_norm

  integer(kind=ENTIER) :: j, k
  real(kind=DOUBLE) :: max_norm_circ
  real(kind=DOUBLE), dimension(3) :: norm, circ, coord_bound

  integer(kind=ENTIER) :: me, num_procs, mpi_ierr
  type(mpi_send_recv_type) :: mpi_send_recv

  integer(kind=ENTIER) :: iter, il, ir, id_vert, id_sub_face, id_face, id_elem
  real(kind=DOUBLE) :: t, tmax, dt, smax, p_nodal, tot_area, p_tilde, r, theta, phi
  real(kind=DOUBLE) :: velocity, w
  real(kind=DOUBLE), dimension(4) :: soll, solr, fluxl, fluxr
  real(kind=DOUBLE), dimension(:, :), allocatable :: sol, flux, sol_exact, vorticity
  real(kind=DOUBLE), dimension(:), allocatable :: divergence
  real(kind=DOUBLE), dimension(:, :, :), allocatable :: lr_states
  logical :: boundary_2d

  namelist /INPUT_PARAM/ meshfile_path, meshfile, n_bc, bc_name, bc_type, bc_val, &
    tmax, cfl, init_uniform, sol_uniform, error_cyl, error_sphere, init_gresho, &
    smooth_deform, write_operators, init_vortex_ring, boundary_2d

  call MPI_INIT(mpi_ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

  if( num_procs > 1 ) then
    print*,"Not yet fitted for MPI"
    error stop
  end if

  open (newunit=funit, file="input_data.f")
  read (nml=INPUT_PARAM, unit=funit)
  close (unit=funit)

  if( mesh_was ) then
    call read_mesh_wasilij(mesh, meshfile_path, meshfile)
  else
    call read_mesh_msh(mesh, meshfile_path, meshfile, &
      n_bc, bc_name, me, num_procs, mpi_send_recv)
  end if
  call build_mesh(mesh, num_procs, mpi_send_recv, .true., boundary_2d)
  call compute_geometry_mesh(mesh, .true., boundary_2d)
  if(smooth_deform) then
    call smooth_deform_mesh(mesh)
    call compute_geometry_mesh(mesh, .true., boundary_2d)
  end if

  allocate(sol(4, mesh%n_elems))
  allocate(vorticity(3, mesh%n_vert))
  allocate(divergence(mesh%n_vert))
  allocate(flux(4, mesh%n_elems))
  if( error_cyl .or. error_sphere ) allocate(sol_exact(4, mesh%n_elems))

  !Init sol
  if( init_uniform ) then
    do id_elem=1, mesh%n_elems
      sol(:, id_elem) = sol_uniform
    end do
  else if (init_gresho) then
    w = 0.2
    do id_elem=1, mesh%n_elems
      r = norm2(mesh%elem(id_elem)%coord(:2))
      if(mesh%elem(id_elem)%coord(1) > 1e-12_DOUBLE) then
        theta = atan2(mesh%elem(id_elem)%coord(2), mesh%elem(id_elem)%coord(1))
      else
        theta = 0.0_DOUBLE
      end if
      velocity = 5.d0

      if(r < w) then
        velocity = r/w
      else if(r < 2.d0*w) then
        velocity = 1.d0- (r-w)/w
      else 
        velocity = 0.d0
      end if

      if( r < 1e-12_DOUBLE ) then
        sol(1,id_elem) = 1.d0
        sol(2,id_elem) = 0.0_DOUBLE
        sol(3,id_elem) = 0.0_DOUBLE
        sol(4,id_elem) = 0.d0
      else
        sol(1,id_elem) = 1.d0
        sol(2,id_elem) = -mesh%elem(id_elem)%coord(2)/r*velocity
        sol(3,id_elem) = mesh%elem(id_elem)%coord(1)/r*velocity
        sol(4,id_elem) = 0.d0
      end if
    end do
  else if (init_vortex_ring) then
    w = 0.1_DOUBLE
    largeR = 0.25_DOUBLE
    do id_elem=1, mesh%n_elems
      distXY = sqrt(mesh%elem(id_elem)%coord(1)**2 + mesh%elem(id_elem)%coord(2)**2) + 1e-14
      r = sqrt((distXY-largeR)**2 + mesh%elem(id_elem)%coord(3)**2)
      phi = atan2(mesh%elem(id_elem)%coord(2), mesh%elem(id_elem)%coord(1)+1e-14)
      theta = asin(mesh%elem(id_elem)%coord(3)/(r+1e-14_DOUBLE))

      if(distXY < largeR) theta = pi - theta

      if(r < w) then
        velocity = 10.0_DOUBLE*r
      else 
        velocity = max(0.0_DOUBLE, 2.0_DOUBLE - 10.0_DOUBLE*r)
      end if

      sol(1,id_elem) = 1._DOUBLE
      if(distXY < 1e-12) then
        sol(2:4,id_elem) = 0._DOUBLE
      else
        sol(2, id_elem) = - sin(theta)*cos(phi)*velocity/distXY
        sol(3, id_elem) = - sin(theta)*sin(phi)*velocity/distXY
        sol(4, id_elem) =   cos(theta)*velocity/distXY
      end if
    end do
  else
    print*,"Unknown init"
    error stop
  end if

  t = 0.0_DOUBLE
  iter = 0

  write(iterchar, *) iter
  output = "output_"//trim(adjustl(iterchar))//".vtk"
  if( mesh_was ) then
    call write_sol_vtk_wasilij(mesh, sol, output)
  else
    call write_sol_vtk(mesh, sol, output)
  end if

  if ( write_operators ) then
    open(newunit=fno, file="operators.dat")
    call compute_operators(mesh, sol, vorticity, vorticity_norm, divergence, divergence_norm)
    write(fno, *) iter, t, vorticity_norm, divergence_norm
  end if

  do while( t < tmax )

    !Compute flux
    dt = 1e30_DOUBLE
    flux = 0.0_DOUBLE
    do id_vert=1, mesh%n_vert

      !Compute left and right states for the subfaces around the node
      allocate(lr_states(4, 2, mesh%vert(id_vert)%n_sub_faces_neigh))
      do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
        id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face

        il = mesh%sub_face(id_sub_face)%left_elem_neigh
        soll = sol(:, il)

        ir = mesh%sub_face(id_sub_face)%right_elem_neigh
        if( ir > 0 ) then
          solr = sol(:, ir)
        else !Boundary conditions
          if( ir == 0 ) then
            if( mesh_was ) then
              if( dot_product(mesh%face(id_face)%norm,mesh%elem(il)%coord) < 0.0_DOUBLE) then
                solr = bc_val(:,-ir)
              else  if( dot_product(mesh%face(id_face)%norm,mesh%elem(il)%coord) > 0.0_DOUBLE) then
                solr(1) = soll(1)
                solr(2:4) = soll(2:4) &
                  - 2.0_DOUBLE*dot_product(soll(2:4), mesh%face(id_face)%norm)*mesh%face(id_face)%norm
              end if
            else
              !Wall bc with oposite normal velocity
              solr(1) = soll(1)
              solr(2:4) = soll(2:4) &
                - 2.0_DOUBLE*dot_product(soll(2:4), mesh%face(id_face)%norm)*mesh%face(id_face)%norm
            end if
          else
            if( bc_type(-ir) == "freestream" ) then
              solr = bc_val(:,-ir)
            else if (trim(adjustl(bc_type(-ir))) == "inflow_pond") then
              if (dot_product(solr(2:4), mesh%face(id_face)%norm) > 0.0_DOUBLE) then
                solr = soll
              else
                solr = bc_val(:, -ir)
              end if
            else if( bc_type(-ir) == "outflow" ) then
              !Symmetric
              solr = soll
            else if( bc_type(-ir) == "sphere" ) then
              !Symmetric
              coord_bound = mesh%face(id_face)%coord + (mesh%face(id_face)%coord - mesh%elem(il)%coord)
              r = norm2(coord_bound)
              theta = atan2(norm2(coord_bound(2:3)), coord_bound(1))
              phi = atan2(coord_bound(3), coord_bound(2))
              call exact_sphere(r,theta,phi,solr)
            else if( bc_type(-ir) == "outflow" ) then
              !Symmetric
              solr = soll
            else if( bc_type(-ir) == "wall" ) then
              !Wall bc with oposite normal velocity
              solr(1) = soll(1)
              solr(2:4) = soll(2:4) &
                - 2.0_DOUBLE*dot_product(soll(2:4), mesh%face(id_face)%norm)*mesh%face(id_face)%norm
            else
              print*, "Unknown BC type!"
              error stop
            end if
          end if
        end if

        !Rotate sol
        call rotate(soll(2:4), mesh%face(id_face)%norm, .true.)
        call rotate(solr(2:4), mesh%face(id_face)%norm, .true.)

        lr_states(:, 1, j) = soll
        lr_states(:, 2, j) = solr
      end do

      !Compute nodal pressure
      p_nodal = 0.0_DOUBLE
      tot_area = 0.0_DOUBLE
      do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
        id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        p_tilde = 0.5_DOUBLE*(lr_states(1, 1, j) + lr_states(1, 2, j)) &
          -0.5_DOUBLE*rho*a*(lr_states(2, 2, j) - lr_states(2, 1, j))
        p_nodal = p_nodal + mesh%sub_face(id_sub_face)%area*p_tilde
        tot_area = tot_area + mesh%sub_face(id_sub_face)%area
      end do
      p_nodal = p_nodal/tot_area

      !Compute the fluxes through the subfaces
      do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
        id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        il = mesh%sub_face(id_sub_face)%left_elem_neigh
        ir = mesh%sub_face(id_sub_face)%right_elem_neigh
        call multi_point_pressure(lr_states(:, 1, j), lr_states(:, 2, j), fluxl, fluxr, smax, p_nodal)
        ! call rusanov(lr_states(:, 1, j), lr_states(:, 2, j), fluxl, fluxr, smax)

        !Rotate back flux
        call rotate(fluxl(2:4), mesh%face(id_face)%norm, .false.)
        call rotate(fluxr(2:4), mesh%face(id_face)%norm, .false.)

        !Update flux
        do k=1, 4
          flux(k, il) = flux(k, il) + mesh%sub_face(id_sub_face)%area*fluxl(k)
        end do
        if( ir > 0 ) then
          do k=1, 4
            flux(k, ir) = flux(k, ir) + mesh%sub_face(id_sub_face)%area*fluxr(k)
          end do
        end if

        dt = min(dt, mesh%elem(il)%volume/(mesh%face(id_face)%area*smax))
        if( ir > 0 ) then
          dt = min(dt, mesh%elem(ir)%volume/(mesh%face(id_face)%area*smax))
        end if
      end do

      deallocate(lr_states)
    end do

    do id_elem=1, mesh%n_elems
      sol(:, id_elem) = sol(:, id_elem) + (cfl*dt/mesh%elem(id_elem)%volume)*flux(:, id_elem)
    end do

    t = t + dt
    iter = iter + 1

    if ( write_operators .and. mod(iter, 100) == 0) then
      call compute_operators(mesh, sol, vorticity, vorticity_norm, divergence, divergence_norm)
      write(fno, *) iter, t, vorticity_norm, divergence_norm
    end if
  end do

  output = "output_final.vtk"
  if( mesh_was ) then
    call write_sol_vtk_wasilij(mesh, sol, output)
  else
    call write_sol_vtk(mesh, sol, output)
  end if

  if ( write_operators ) then
    close(fno)
  end if

  if(error_cyl) then
    call compute_error_cyl(mesh, sol)

    !Write exact solution into sol
    do id_elem=1, mesh%n_elems
      r = norm2(mesh%elem(id_elem)%coord)
      theta = atan2(mesh%elem(id_elem)%coord(2), mesh%elem(id_elem)%coord(1))
      call exact_cyl(r,theta,sol(:, id_elem))
    end do
    output = "output_exact.vtk"
    if( mesh_was ) then
      call write_sol_vtk_wasilij(mesh, sol, output)
    else
      call write_sol_vtk(mesh, sol, output)
    end if
  else if(error_sphere) then
    call compute_error_sphere(mesh, sol)

    !Write exact solution into sol
    do id_elem=1, mesh%n_elems
      r = norm2(mesh%elem(id_elem)%coord)
      theta = atan2(norm2(mesh%elem(id_elem)%coord(2:3)), mesh%elem(id_elem)%coord(1))
      phi = atan2(mesh%elem(id_elem)%coord(3), mesh%elem(id_elem)%coord(2))
      soll = sol(:, id_elem)
      call exact_sphere(r,theta,phi,solr)
      ! sol(:, id_elem) = mesh%elem(id_elem)%volume*norm2(soll-solr)**2
      sol(:, id_elem) = solr
    end do
    output = "output_exact.vtk"
    if( mesh_was ) then
      call write_sol_vtk_wasilij(mesh, sol, output)
    else
      call write_sol_vtk(mesh, sol, output)
    end if
  end if
end program main
