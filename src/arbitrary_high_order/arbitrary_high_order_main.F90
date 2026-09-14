program arbitrary_high_order_main
  ! Standalone accuracy/cost check for arbitrary_high_order_module: builds
  ! the order 1..max_order cell derivatives of a scalar test field and
  ! dumps everything to a .vtu for inspection. test_case selects the field:
  !   affine       phi = x+2y+3z: Step 1 must recover (1,2,3) to roundoff.
  !   vortex       classical Shu isentropic vortex density (3D mesh; needs
  !                boundary_2d=.false., see arbitrary_high_order_module).
  !   smooth       sin(x)cos(y) on a thin-extruded 2D mesh (boundary_2d=.true.).
  !   discontinuous  a jump in x, for qualitative non-oscillation checks.
  ! Meant as the harness behind the convergence and cost data in
  ! tex_arbitrary_high_order.
  use mpi
  use precision_module
  use mpi_module
  use mesh_module
  use mesh_reading_module
  use mesh_geometry_module
  use mesh_connectivity_module
  use io_module
  use quadrature_module
  use arbitrary_high_order_module

  implicit none

  integer(kind=ENTIER), parameter :: d = 3

  integer(kind=ENTIER) :: n_bc = 0, me, num_procs, fn, mpi_ierr
  character(len=255) :: meshfile_path, meshfile
  character(len=255), dimension(10) :: bc_name
  type(mesh_type) :: mesh
  type(mpi_send_recv_type) :: mpi_send_recv
  logical :: b2d

  integer(kind=ENTIER) :: max_order = 3
  character(len=32) :: test_case = "smooth"
  real(kind=DOUBLE) :: gamma_gas = 1.4_DOUBLE, vortex_beta = 5.0_DOUBLE
  real(kind=DOUBLE) :: interior_radius = 3.0_DOUBLE
  logical :: use_overlap = .FALSE.
  logical :: dump_timeline = .FALSE.

  namelist /params/ meshfile, max_order, test_case, n_bc, bc_name, b2d, &
    gamma_gas, vortex_beta, interior_radius, use_overlap, dump_timeline

  integer(kind=ENTIER) :: fn_vtu, fn_pvtu
  integer(kind=ENTIER) :: i, order, k
  real(kind=DOUBLE) :: x, y
  real(kind=DOUBLE) :: wt1, wt2, wt_local, wt_max
  integer(kind=ENTIER) :: n_interior
  real(kind=DOUBLE), dimension(:), allocatable :: phi
  type(derivative_field_type), dimension(:), allocatable :: dfield
  character(len=16) :: comp_name
  type(timeline_event_type), dimension(:), allocatable :: events
  integer(kind=ENTIER) :: n_events, ev, fn_tl
  character(len=255) :: tl_fname

  call MPI_INIT(mpi_ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

  meshfile_path = ""
  b2d = .TRUE.

  open(newunit=fn, file="input_data.f", status="old")
  read(fn, nml=params)
  close(fn)

  call read_mesh_msh(mesh, meshfile_path, meshfile, &
    n_bc, bc_name, me, num_procs, mpi_send_recv)
  call build_mesh(mesh, num_procs, mpi_send_recv, .TRUE., b2d)
  call compute_geometry_mesh(mesh, .TRUE., b2d)

  allocate(phi(mesh%n_elems))
  do i = 1, mesh%n_elems
    select case (trim(adjustl(test_case)))
    case ("affine")
      ! Unit test for Step 1: phi = x + 2y + 3z is affine everywhere, so
      ! the least-squares gradient must recover (1, 2, 3) exactly (zero
      ! residual) at every vertex and hence every cell, on any mesh.
      phi(i) = mesh%elem(i)%coord(1) + 2.0_DOUBLE*mesh%elem(i)%coord(2) &
        + 3.0_DOUBLE*mesh%elem(i)%coord(3)
    case ("discontinuous")
      x = mesh%elem(i)%coord(1)
      if (x < 0.0_DOUBLE) then
        phi(i) = (x-1.0_DOUBLE)**2 * x**2
      else
        phi(i) = 20.0_DOUBLE + (x+1.0_DOUBLE)**2 * x**2
      end if
    case ("vortex")
      ! Quadrature-accurate volumetric cell average of the classical Shu
      ! isentropic vortex density (same formula and beta/gamma convention
      ! as euler_ho_module's vortex_prim, at rest, t=0), on a genuine 3D
      ! mesh -- reusing quadrature_module's hex volume rule, exactly the
      ! way euler_ho_module's init_sol uses face_quad_pts for its own
      ! (2D, extruded) vortex cell average.
      phi(i) = vortex_density_cell_average(mesh, i, gamma_gas, vortex_beta)
    case default ! "smooth"
      ! Quadrature-accurate cell average of sin(x)cos(y) over the cell's
      ! z-face, rather than a point value at the centroid -- same
      ! quadrature-based cell-averaging convention as the concurrent
      ! high-order Euler work (euler_ho_module's init_sol/quadrature_module),
      ! so phi is a genuine FV cell average, not a point sample.
      phi(i) = smooth_field_cell_average(mesh, i)
    end select
  end do

  n_interior = count(.not. mesh%elem%is_ghost)
  allocate(dfield(max_order))

  ! MPI_WTIME (wall clock), maxed over ranks, is the meaningful "did
  ! overlap help" number -- cpu_time would miss time spent blocked in MPI
  ! waits, which is exactly what use_overlap=.true. is meant to reduce.
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
  wt1 = MPI_WTIME()
  if (dump_timeline) then
    allocate(events(4*max_order))
    n_events = 0
    if (use_overlap) then
      call compute_derivative_hierarchy_overlap_timed(mesh, mpi_send_recv, num_procs, &
        d, b2d, max_order, phi, dfield, events, n_events)
    else
      call compute_derivative_hierarchy_timed(mesh, mpi_send_recv, num_procs, &
        d, b2d, max_order, phi, dfield, events, n_events)
    end if
  else if (use_overlap) then
    call compute_derivative_hierarchy_overlap(mesh, mpi_send_recv, num_procs, &
      d, b2d, max_order, phi, dfield)
  else
    call compute_derivative_hierarchy(mesh, mpi_send_recv, num_procs, &
      d, b2d, max_order, phi, dfield)
  end if
  wt2 = MPI_WTIME()

  if (dump_timeline) then
    write(tl_fname, '(A,I0,A)') 'timeline_', me, '.csv'
    open(newunit=fn_tl, file=trim(tl_fname), status='replace')
    write(fn_tl, '(A)') 'rank,order,phase,t0,t1'
    do ev = 1, n_events
      write(fn_tl, '(I0,A,I0,A,A,A,ES16.9,A,ES16.9)') me, ',', events(ev)%order, ',', &
        trim(events(ev)%phase), ',', events(ev)%t0 - wt1, ',', events(ev)%t1 - wt1
    end do
    close(fn_tl)
    deallocate(events)
  end if
  wt_local = wt2 - wt1
  call MPI_REDUCE(wt_local, wt_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, mpi_ierr)

  if (me == 0) then
    print *, "use_overlap: ", use_overlap
    print *, "n_elems (interior, rank 0): ", n_interior
    print *, "derivative hierarchy wall time, max over ranks (s): ", wt_max
  end if

  ! Report the Linf error against the known exact derivative, away from the
  ! physical boundary where the dual derivative is one-sided by
  ! construction (see the note above compute_nodal_derivative_at_vertex in
  ! arbitrary_high_order_module), not by truncation error.
  select case (trim(adjustl(test_case)))
  case ("affine")
    call report_affine_errors(mesh, max_order, dfield)
  case ("vortex")
    call report_vortex_errors(mesh, max_order, dfield, gamma_gas, vortex_beta, &
      interior_radius)
  case ("smooth")
    call report_smooth_errors(mesh, max_order, dfield, interior_radius)
  end select

  call open_file_vtu(mesh, "arbitrary_high_order", fn_vtu, fn_pvtu)
  call write_file_vtu_start_cell_data(mesh, "arbitrary_high_order", fn_vtu, fn_pvtu)
  call write_file_vtu_cell_scalar(mesh, "arbitrary_high_order", fn_vtu, fn_pvtu, phi, "phi")

  if (max_order >= 1) then
    call write_file_vtu_cell_vector(mesh, "arbitrary_high_order", fn_vtu, fn_pvtu, &
      dfield(1)%val, "dphi")
  end if

  do order = 2, max_order
    do k = 1, d**order
      write(comp_name, '(A,I0,A,I0)') "d", order, "phi_", k
      call write_file_vtu_cell_scalar(mesh, "arbitrary_high_order", fn_vtu, fn_pvtu, &
        dfield(order)%val(k, :), trim(comp_name))
    end do
  end do

  call write_file_vtu_end_cell_data(mesh, "arbitrary_high_order", fn_vtu, fn_pvtu)
  call close_file_vtu(mesh, "arbitrary_high_order", fn_vtu, fn_pvtu)

  call MPI_FINALIZE(mpi_ierr)

contains
  ! Unit test for Step 1 (compute_nodal_derivative_at_vertex): phi = x +
  ! 2y + 3z is affine everywhere, so the least-squares gradient must
  ! recover (1, 2, 3) exactly (zero residual, up to roundoff) at every
  ! vertex and hence every cell, on any mesh -- this is the exactness
  ! property motivating the switch away from the Green-Gauss jump formula
  ! (see arbitrary_high_order_module's header).
  subroutine report_affine_errors(mesh, max_order, dfield)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: max_order
    type(derivative_field_type), dimension(max_order), intent(in) :: dfield

    real(kind=DOUBLE), dimension(3), parameter :: grad_exact = &
      (/ 1.0_DOUBLE, 2.0_DOUBLE, 3.0_DOUBLE /)
    integer(kind=ENTIER) :: i, me, mpi_ierr
    real(kind=DOUBLE) :: linf1

    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    linf1 = 0.0_DOUBLE
    do i = 1, mesh%n_elems
      if (mesh%elem(i)%is_ghost) cycle
      linf1 = max(linf1, maxval(abs(dfield(1)%val(:, i) - grad_exact)))
    end do

    if (me == 0) print *, "Linf error, gradient (order 1, should be ~0): ", linf1
  end subroutine report_affine_errors

  ! Classical Shu isentropic vortex, at rest (u_bg=v_bg=0), t=0, centered
  ! at the origin -- same formula as euler_ho_module's vortex_prim. Returns
  ! the density and its exact gradient/Hessian (flattened with the same
  ! i1-fast/i2-slow convention as arbitrary_high_order_module), analytic
  ! since rho(x,y) = g(r^2)^n with g = 1 - A*exp(1-r^2), n = 1/(gamma-1).
  ! rho does not depend on z, so every z-derivative is exactly 0.
  subroutine vortex_exact(x, y, gamma, beta, rho, grad, hess)
    implicit none

    real(kind=DOUBLE), intent(in) :: x, y, gamma, beta
    real(kind=DOUBLE), intent(out) :: rho
    real(kind=DOUBLE), dimension(3), intent(out) :: grad
    real(kind=DOUBLE), dimension(9), intent(out) :: hess

    real(kind=DOUBLE), parameter :: PI = 4.0_DOUBLE * atan(1.0_DOUBLE)
    real(kind=DOUBLE) :: a_coef, r2, e_val, g_val, h_val, n_exp
    real(kind=DOUBLE) :: rho_x, rho_y, rho_xx, rho_yy, rho_xy

    a_coef = (gamma - 1.0_DOUBLE) * beta**2 / (8.0_DOUBLE * gamma * PI**2)
    n_exp = 1.0_DOUBLE / (gamma - 1.0_DOUBLE)
    r2 = x*x + y*y
    e_val = exp(1.0_DOUBLE - r2)
    g_val = 1.0_DOUBLE - a_coef*e_val
    h_val = 2.0_DOUBLE * a_coef * e_val

    rho = g_val**n_exp

    rho_x = n_exp * g_val**(n_exp-1.0_DOUBLE) * h_val * x
    rho_y = n_exp * g_val**(n_exp-1.0_DOUBLE) * h_val * y

    rho_xx = n_exp*(n_exp-1.0_DOUBLE) * g_val**(n_exp-2.0_DOUBLE) * h_val**2 * x**2 &
      + n_exp * g_val**(n_exp-1.0_DOUBLE) * h_val * (1.0_DOUBLE - 2.0_DOUBLE*x**2)
    rho_yy = n_exp*(n_exp-1.0_DOUBLE) * g_val**(n_exp-2.0_DOUBLE) * h_val**2 * y**2 &
      + n_exp * g_val**(n_exp-1.0_DOUBLE) * h_val * (1.0_DOUBLE - 2.0_DOUBLE*y**2)
    rho_xy = x*y * ( n_exp*(n_exp-1.0_DOUBLE) * g_val**(n_exp-2.0_DOUBLE) * h_val**2 &
      - 2.0_DOUBLE * n_exp * g_val**(n_exp-1.0_DOUBLE) * h_val )

    grad = (/ rho_x, rho_y, 0.0_DOUBLE /)
    hess = (/ rho_xx, rho_xy, 0.0_DOUBLE, &
              rho_xy, rho_yy, 0.0_DOUBLE, &
              0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE /)
  end subroutine vortex_exact

  ! Quadrature-accurate volumetric cell average of the vortex density over
  ! a hexahedral cell, using quadrature_module's hex rule directly on the
  ! cell's 8 vertices (gmsh/VTK node order, as documented in
  ! quadrature_module and assumed preserved by mesh_reading_module).
  function vortex_density_cell_average(mesh, id_elem, gamma, beta) result(phi_avg)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_elem
    real(kind=DOUBLE), intent(in) :: gamma, beta
    real(kind=DOUBLE) :: phi_avg

    real(kind=DOUBLE), dimension(:, :), allocatable :: qpts
    real(kind=DOUBLE), dimension(:), allocatable :: qwts
    real(kind=DOUBLE), dimension(3, 8) :: coords
    real(kind=DOUBLE) :: rho, val_sum, vol_q
    real(kind=DOUBLE), dimension(3) :: grad_dummy
    real(kind=DOUBLE), dimension(9) :: hess_dummy
    integer(kind=ENTIER) :: kv, q

    if (mesh%elem(id_elem)%n_vert /= 8) then
      call vortex_exact(mesh%elem(id_elem)%coord(1), mesh%elem(id_elem)%coord(2), &
        gamma, beta, phi_avg, grad_dummy, hess_dummy)
      return
    end if

    do kv = 1, 8
      coords(:, kv) = mesh%vert(mesh%elem(id_elem)%vert(kv))%coord
    end do
    call volume_quad_pts(8_ENTIER, coords, 3_ENTIER, qpts, qwts)

    val_sum = 0.0_DOUBLE
    vol_q = 0.0_DOUBLE
    do q = 1, size(qwts)
      call vortex_exact(qpts(1, q), qpts(2, q), gamma, beta, rho, grad_dummy, hess_dummy)
      val_sum = val_sum + qwts(q) * rho
      vol_q = vol_q + qwts(q)
    end do
    deallocate(qpts, qwts)

    phi_avg = val_sum / vol_q
  end function vortex_density_cell_average

  ! Linf error of dfield(1)/dfield(2) against the exact vortex
  ! gradient/Hessian, restricted to cells whose centroid is well inside
  ! the box (|x|,|y|,|z| < interior_radius), away from all 6 physical
  ! boundary faces where the dual derivative is one-sided by construction.
  ! Reduced (MPI_MAX) across ranks, so the reported number is the true
  ! global error regardless of how many partitions the mesh was split
  ! into -- necessary for comparing a multi-rank run against a serial
  ! reference, since any single rank only owns a fraction of the domain
  ! and its own local max is not otherwise comparable across partition
  ! counts.
  subroutine report_vortex_errors(mesh, max_order, dfield, gamma, beta, interior_radius)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: max_order
    type(derivative_field_type), dimension(max_order), intent(in) :: dfield
    real(kind=DOUBLE), intent(in) :: gamma, beta, interior_radius

    integer(kind=ENTIER) :: i, me, mpi_ierr
    real(kind=DOUBLE) :: rho_exact, linf1, linf2, linf1_glob, linf2_glob
    real(kind=DOUBLE), dimension(3) :: xc, grad_exact
    real(kind=DOUBLE), dimension(9) :: hess_exact

    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    linf1 = 0.0_DOUBLE
    linf2 = 0.0_DOUBLE
    do i = 1, mesh%n_elems
      if (mesh%elem(i)%is_ghost) cycle
      xc = mesh%elem(i)%coord
      if (maxval(abs(xc)) >= interior_radius) cycle

      call vortex_exact(xc(1), xc(2), gamma, beta, rho_exact, grad_exact, hess_exact)

      linf1 = max(linf1, maxval(abs(dfield(1)%val(:, i) - grad_exact)))
      if (max_order >= 2) then
        linf2 = max(linf2, maxval(abs(dfield(2)%val(:, i) - hess_exact)))
      end if
    end do

    call MPI_REDUCE(linf1, linf1_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_REDUCE(linf2, linf2_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, mpi_ierr)

    if (me == 0) then
      print *, "Linf error, gradient (order 1): ", linf1_glob
      if (max_order >= 2) print *, "Linf error, Hessian (order 2): ", linf2_glob
    end if
  end subroutine report_vortex_errors

  ! Quadrature-accurate cell average of sin(x)cos(y) over cell id_elem's
  ! z-face (these test meshes are 2D slabs extruded in z), using the same
  ! face_quad_pts / quadrature_module machinery as euler_ho_module's
  ! init_sol for its isentropic vortex cell average. Falls back to the
  ! centroid point value if no z-face is found (should not happen on these
  ! meshes).
  function smooth_field_cell_average(mesh, id_elem) result(phi_avg)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_elem
    real(kind=DOUBLE) :: phi_avg

    integer(kind=ENTIER) :: kf, iface_loc, n_fv, kv, q
    real(kind=DOUBLE), dimension(:, :), allocatable :: qpts, fc
    real(kind=DOUBLE), dimension(:), allocatable :: qwts
    real(kind=DOUBLE) :: val_sum, area_q

    do kf = 1, mesh%elem(id_elem)%n_faces
      iface_loc = mesh%elem(id_elem)%face(kf)
      if (abs(abs(mesh%face(iface_loc)%norm(3)) - 1.0_DOUBLE) < 1.0e-6_DOUBLE) then
        n_fv = mesh%face(iface_loc)%n_vert
        if (n_fv >= 3 .and. allocated(mesh%face(iface_loc)%vert)) then
          allocate(fc(3, n_fv))
          do kv = 1, n_fv
            fc(:, kv) = mesh%vert(mesh%face(iface_loc)%vert(kv))%coord
          end do
          call face_quad_pts(n_fv, fc, 3_ENTIER, qpts, qwts)
          deallocate(fc)

          val_sum = 0.0_DOUBLE
          area_q = 0.0_DOUBLE
          do q = 1, size(qwts)
            val_sum = val_sum + qwts(q) * sin(qpts(1, q)) * cos(qpts(2, q))
            area_q = area_q + qwts(q)
          end do
          deallocate(qpts, qwts)

          phi_avg = val_sum / area_q
          return
        end if
      end if
    end do

    phi_avg = sin(mesh%elem(id_elem)%coord(1)) * cos(mesh%elem(id_elem)%coord(2))
  end function smooth_field_cell_average

  ! Linf error of dfield(1) and dfield(2) against the exact gradient/Hessian
  ! of sin(x)*cos(y), restricted to cells within interior_radius of the
  ! origin (the "smooth" test case is run on a mesh with a farfield outer
  ! boundary, whose adjacent cells carry a one-sided-stencil bias by
  ! construction rather than a discretization error -- see the note above
  ! compute_nodal_derivative_at_vertex in arbitrary_high_order_module).
  subroutine report_smooth_errors(mesh, max_order, dfield, interior_radius)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: max_order
    type(derivative_field_type), dimension(max_order), intent(in) :: dfield
    real(kind=DOUBLE), intent(in) :: interior_radius

    integer(kind=ENTIER) :: i, me, mpi_ierr
    real(kind=DOUBLE) :: xc, yc, linf1, linf2
    real(kind=DOUBLE), dimension(3) :: grad_exact
    real(kind=DOUBLE), dimension(9) :: hess_exact

    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    linf1 = 0.0_DOUBLE
    linf2 = 0.0_DOUBLE
    do i = 1, mesh%n_elems
      if (mesh%elem(i)%is_ghost) cycle
      if (norm2(mesh%elem(i)%coord(1:2)) >= interior_radius) cycle

      xc = mesh%elem(i)%coord(1)
      yc = mesh%elem(i)%coord(2)

      grad_exact = (/ cos(xc)*cos(yc), -sin(xc)*sin(yc), 0.0_DOUBLE /)
      linf1 = max(linf1, maxval(abs(dfield(1)%val(:, i) - grad_exact)))

      if (max_order >= 2) then
        hess_exact = (/ -sin(xc)*cos(yc), -cos(xc)*sin(yc), 0.0_DOUBLE, &
                         -cos(xc)*sin(yc), -sin(xc)*cos(yc), 0.0_DOUBLE, &
                          0.0_DOUBLE,       0.0_DOUBLE,       0.0_DOUBLE /)
        linf2 = max(linf2, maxval(abs(dfield(2)%val(:, i) - hess_exact)))
      end if
    end do

    if (me == 0) then
      print *, "Linf error, gradient (order 1): ", linf1
      if (max_order >= 2) print *, "Linf error, Hessian (order 2): ", linf2
    end if
  end subroutine report_smooth_errors
end program arbitrary_high_order_main
