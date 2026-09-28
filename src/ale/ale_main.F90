program main
  use mpi
  use precision_module
  use mpi_module
  use mesh_module
  use mesh_reading_module
  use mesh_geometry_module
  use mesh_connectivity_module
  use io_module

  use ale_global_data_module
  use ale_module
  use ale_rbf_module

  implicit none

  integer(kind=ENTIER) :: me, num_procs, mpi_ierr
  integer(kind=ENTIER) :: iter, i, j, i_sol_vtu
  type(mesh_type) :: mesh
  type(mpi_send_recv_type) :: mpi_send_recv
  integer(kind=ENTIER) :: fn_vtu, fn_pvtu
  character(len=255) :: fln

  real(kind=DOUBLE), dimension(:, :), allocatable :: sol, rhs
  real(kind=DOUBLE), dimension(:), allocatable :: gamma_arr
  real(kind=DOUBLE), dimension(:, :), allocatable :: vp, wp
  real(kind=DOUBLE), dimension(:), allocatable :: sum_lambda, vol_old
  real(kind=DOUBLE), dimension(:), allocatable :: density, pressure, gamma_out
  real(kind=DOUBLE), dimension(:, :), allocatable :: velocity, centroid
  real(kind=DOUBLE), dimension(5) :: w

  ! --- RBF hybrid grid-velocity setup (ale_grid_velocity=='rbf_hybrid'
  ! only): hard points (bubble interface only, w_p=v_p exactly) and sliding
  ! points (every other domain-boundary vertex, constrained to a fixed
  ! normal but otherwise free tangentially -- the piston wall included,
  ! with a nonzero imposed value piston_vel.n instead of the usual 0 for a
  ! static wall) feed a compactly-supported RBF interpolant (ale_rbf_module,
  ! a local copy of src/rbf/rbf_module.F90) whose evaluation at every
  ! vertex gives a smooth w_p elsewhere. Piston nodes are hard (full 3D
  ! vector = piston_vel, not a "sliding" 1-dof point): a corner node
  ! shared with a side wall has an ambiguous single sliding normal (only
  ! one of the two walls' normals can be recorded there), which left
  ! piston/wall corner nodes stuck instead of moving with the piston.
  ! piston_vel is spatially constant, so -- by the same "a smooth/constant
  ! target keeps the RBF system well-behaved regardless of point density"
  ! argument used for the sliding walls -- treating every piston node as
  ! hard with that same constant vector does not reintroduce the
  ! conditioning problem that motivated thinning the bubble interface
  ! (bubble_iface_stride): unlike the bubble, the piston's imposed value
  ! doesn't vary from node to node. Point membership is a mesh-topology
  ! property set once from the initial condition (material placement +
  ! boundary tags); only positions (and, for the bubble's hard points,
  ! values) are refreshed every iteration -- the piston's hard values are
  ! constant and set once.
  integer(kind=ENTIER) :: nimp, nimp_bubble, nslide, n_rbf, k
  integer(kind=ENTIER) :: piston_bc_idx, id_sub_face, id_face, ide, n_bubble_iface_seen
  integer(kind=ENTIER), dimension(4) :: n_wall_seen
  logical, dimension(:), allocatable :: is_hard_vert, is_piston_hard_vert
  integer(kind=ENTIER), dimension(:), allocatable :: imp_vert
  real(kind=DOUBLE), dimension(:, :), allocatable :: x_rbf, no_rbf
  real(kind=DOUBLE), dimension(:), allocatable :: delta_imp_rbf, delta_rbf
  real(kind=DOUBLE), dimension(3) :: piston_vel, no_tmp
  logical :: has_air, has_bub, is_bubble_iface, is_piston_vert, has_lateral_bound_face
  integer(kind=ENTIER) :: wall_dir
  ! Exact axis-aligned unit normals, keyed by wall_dir (1=+x,2=-x,3=+y,4=-y):
  ! the raw sub-face normal computed from mesh geometry is only numerically
  ! close to axis-aligned (mesh generation/geometry round-off), and using it
  ! as-is for a sliding wall's RBF constraint leaves a small residual
  ! velocity component perpendicular to the wall instead of exactly zero.
  ! Snapping to the exact axis vector (once we already know which wall this
  ! is, from the same classification used for the wall_stride counters)
  ! removes that residual.
  real(kind=DOUBLE), dimension(3, 4), parameter :: axis_normal = reshape( &
    [1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, -1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, &
     0.0_DOUBLE, 1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, -1.0_DOUBLE, 0.0_DOUBLE], [3, 4])

  call MPI_INIT(mpi_ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

  call read_input_parameters("input_data.f")
  call init_bc_flags()

  call read_mesh_msh(mesh, meshfile_path, meshfile, &
    n_bc, bc_name, me, num_procs, mpi_send_recv)
  call build_mesh(mesh, num_procs, mpi_send_recv, .true., boundary_2d)
  call compute_geometry_mesh(mesh, .true., boundary_2d, me, num_procs)

  allocate(sol(5, mesh%n_elems))     ; sol = 0.0_DOUBLE
  allocate(gamma_arr(mesh%n_elems))  ; gamma_arr = gamma
  allocate(rhs(5, mesh%n_elems))     ; rhs = 0.0_DOUBLE
  allocate(sum_lambda(mesh%n_elems)) ; sum_lambda = 0.0_DOUBLE
  allocate(vol_old(mesh%n_elems))    ; vol_old = 0.0_DOUBLE
  allocate(vp(3, mesh%n_vert))       ; vp = 0.0_DOUBLE
  allocate(wp(3, mesh%n_vert))       ; wp = 0.0_DOUBLE

  allocate(density(mesh%n_elems))
  allocate(pressure(mesh%n_elems))
  allocate(gamma_out(mesh%n_elems))
  allocate(velocity(3, mesh%n_elems))
  allocate(centroid(3, mesh%n_elems))

  call init_sol(mesh, sol, gamma_arr)
  if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol)

  ! --- RBF hybrid setup: classify vertices once, from the initial gamma
  ! field and the piston boundary tag ---
  if (trim(ale_grid_velocity) == 'rbf_hybrid') then
    piston_bc_idx = 0
    piston_vel = 0.0_DOUBLE
    do i = 1, n_bc
      if (bc_is_piston(i)) then
        piston_bc_idx = i
        piston_vel = bc_val(2:4, i)
        exit
      end if
    end do

    allocate(is_hard_vert(mesh%n_vert)); is_hard_vert = .FALSE.
    allocate(is_piston_hard_vert(mesh%n_vert)); is_piston_hard_vert = .FALSE.
    nimp_bubble = 0
    nimp = 0
    nslide = 0
    n_bubble_iface_seen = 0
    n_wall_seen = 0
    do i = 1, mesh%n_vert
      has_air = .FALSE.; has_bub = .FALSE.
      do j = 1, mesh%vert(i)%n_elems_neigh
        ide = mesh%vert(i)%elem_neigh(j)
        if (abs(gamma_arr(ide) - gamma_air) < 1e-8_DOUBLE) has_air = .TRUE.
        if (abs(gamma_arr(ide) - gamma_bub) < 1e-8_DOUBLE) has_bub = .TRUE.
      end do
      is_bubble_iface = has_air .and. has_bub
      ! Keep only every bubble_iface_stride-th interface vertex as a hard
      ! RBF landmark -- see ale_global_data_module for why (RBF matrix
      ! conditioning at a domain-spanning radius).
      if (is_bubble_iface) then
        n_bubble_iface_seen = n_bubble_iface_seen + 1
        is_bubble_iface = (mod(n_bubble_iface_seen - 1, bubble_iface_stride) == 0)
      end if

      is_piston_vert = .FALSE.
      if (piston_bc_idx > 0) then
        do j = 1, mesh%vert(i)%n_sub_faces_neigh
          id_sub_face = mesh%vert(i)%sub_face_neigh(j)
          id_face = mesh%sub_face(id_sub_face)%mesh_face
          if (mesh%face(id_face)%right_neigh == -piston_bc_idx) then
            is_piston_vert = .TRUE.
            exit
          end if
        end do
      end if

      if (is_bubble_iface) then
        is_hard_vert(i) = .TRUE.
        nimp_bubble = nimp_bubble + 1
      else if (is_piston_vert) then
        is_hard_vert(i) = .TRUE.
        is_piston_hard_vert(i) = .TRUE.
      else if (mesh%vert(i)%is_bound) then
        ! Sliding candidate: exclude the top/bottom (z-normal) faces of the
        ! boundary_2d single-layer extrusion, only the true 2D domain edges
        ! slide.
        has_lateral_bound_face = .FALSE.
        no_tmp = 0.0_DOUBLE
        do j = 1, mesh%vert(i)%n_sub_faces_neigh
          id_sub_face = mesh%vert(i)%sub_face_neigh(j)
          id_face = mesh%sub_face(id_sub_face)%mesh_face
          if (mesh%face(id_face)%right_neigh <= 0) then
            if (abs(mesh%sub_face(id_sub_face)%norm(3)) < 1e-8_DOUBLE) then
              has_lateral_bound_face = .TRUE.
              no_tmp = mesh%sub_face(id_sub_face)%norm
            end if
          end if
        end do
        ! Same thinning idea as the bubble interface, applied to every
        ! (non-piston) wall sliding candidate (see ale_global_data_module).
        ! Counted per wall direction (dominant sign of the local normal),
        ! not with one counter shared across every wall combined: a single
        ! shared counter aliases against the mesh's own vertex-numbering
        ! order (e.g. all of the bottom wall's nodes numbered before the
        ! top's), so "every wall_stride-th" can end up keeping most of one
        ! wall and almost none of another purely by numbering coincidence
        ! -- observed on the refined mesh as 3 kept top-wall points against
        ! 11 bottom-wall ones, leaving the top wall's RBF motion almost
        ! unconstrained and visibly drooping into the domain.
        if (has_lateral_bound_face) then
          if (abs(no_tmp(1)) > abs(no_tmp(2))) then
            wall_dir = merge(1, 2, no_tmp(1) > 0.0_DOUBLE)
          else
            wall_dir = merge(3, 4, no_tmp(2) > 0.0_DOUBLE)
          end if
          n_wall_seen(wall_dir) = n_wall_seen(wall_dir) + 1
          has_lateral_bound_face = (mod(n_wall_seen(wall_dir) - 1, wall_stride) == 0)
        end if
        if (has_lateral_bound_face) nslide = nslide + 1
      end if
    end do

    nimp = count(is_hard_vert)
    n_rbf = nimp + nslide
    if (me == 0) print*, "RBF hybrid: n_vert, nimp_bubble, nimp_piston, nslide (walls) =", &
      mesh%n_vert, nimp_bubble, nimp - nimp_bubble, nslide

    allocate(imp_vert(n_rbf))
    allocate(x_rbf(3, n_rbf))
    allocate(no_rbf(3, n_rbf))
    allocate(delta_imp_rbf(3*nimp + nslide))
    allocate(delta_rbf(3*nimp + nslide))
    no_rbf = 0.0_DOUBLE
    delta_imp_rbf = 0.0_DOUBLE
    delta_rbf = 0.0_DOUBLE

    ! Hard points: bubble interface first (1:nimp_bubble, value refreshed
    ! every iteration from v_p), then piston (nimp_bubble+1:nimp, value
    ! constant = piston_vel, set once here and never touched again).
    j = 0
    do i = 1, mesh%n_vert
      if (is_hard_vert(i) .and. .not. is_piston_hard_vert(i)) then
        j = j + 1
        imp_vert(j) = i
      end if
    end do
    do i = 1, mesh%n_vert
      if (is_piston_hard_vert(i)) then
        j = j + 1
        imp_vert(j) = i
        delta_imp_rbf(dof_offset(j, nimp) + 1:dof_offset(j, nimp) + 3) = piston_vel
      end if
    end do

    k = nimp
    n_wall_seen = 0
    do i = 1, mesh%n_vert
      if (.not. is_hard_vert(i) .and. mesh%vert(i)%is_bound) then
        has_lateral_bound_face = .FALSE.
        no_tmp = 0.0_DOUBLE
        do j = 1, mesh%vert(i)%n_sub_faces_neigh
          id_sub_face = mesh%vert(i)%sub_face_neigh(j)
          id_face = mesh%sub_face(id_sub_face)%mesh_face
          if (mesh%face(id_face)%right_neigh <= 0) then
            if (abs(mesh%sub_face(id_sub_face)%norm(3)) < 1e-8_DOUBLE) then
              has_lateral_bound_face = .TRUE.
              no_tmp = mesh%sub_face(id_sub_face)%norm
            end if
          end if
        end do
        if (has_lateral_bound_face) then
          if (abs(no_tmp(1)) > abs(no_tmp(2))) then
            wall_dir = merge(1, 2, no_tmp(1) > 0.0_DOUBLE)
          else
            wall_dir = merge(3, 4, no_tmp(2) > 0.0_DOUBLE)
          end if
          n_wall_seen(wall_dir) = n_wall_seen(wall_dir) + 1
          has_lateral_bound_face = (mod(n_wall_seen(wall_dir) - 1, wall_stride) == 0)
        end if
        if (has_lateral_bound_face) then
          k = k + 1
          imp_vert(k) = i
          no_rbf(:, k) = axis_normal(:, wall_dir)
          ! delta_imp stays 0 here: a static wall's sliding value along its
          ! normal.
        end if
      end if
    end do
  end if

  i_sol_vtu = 0
  write(fln, *) i_sol_vtu
  fln = "output_"//trim(adjustl(fln))
  call open_file_vtu(mesh, trim(fln), fn_vtu, fn_pvtu)
  call write_file_vtu_start_cell_data(mesh, trim(fln), fn_vtu, fn_pvtu)
  do i = 1, mesh%n_elems
    w = conserv_to_primit(sol(:, i), gamma_arr(i))
    density(i) = w(1); velocity(:, i) = w(2:4); pressure(i) = w(5)
    gamma_out(i) = gamma_arr(i)
    centroid(:, i) = mesh%elem(i)%coord
  end do
  call write_file_vtu_cell_scalar(mesh, trim(fln), fn_vtu, fn_pvtu, density, "Density")
  call write_file_vtu_cell_vector(mesh, trim(fln), fn_vtu, fn_pvtu, velocity, "Velocity")
  call write_file_vtu_cell_scalar(mesh, trim(fln), fn_vtu, fn_pvtu, pressure, "Pressure")
  call write_file_vtu_cell_scalar(mesh, trim(fln), fn_vtu, fn_pvtu, gamma_out, "Gamma")
  call write_file_vtu_cell_vector(mesh, trim(fln), fn_vtu, fn_pvtu, centroid, "Centroid")
  call write_file_vtu_end_cell_data(mesh, trim(fln), fn_vtu, fn_pvtu)
  call close_file_vtu(mesh, trim(fln), fn_vtu, fn_pvtu)
  i_sol_vtu = i_sol_vtu + 1

  t = 0.0_DOUBLE
  iter = 0
  do while (t < tmax)
    if (trim(ale_grid_velocity) == 'lagrangian') then
      call compute_nodal_velocity_field(mesh, sol, gamma_arr, vp)
      wp = vp
    else if (trim(ale_grid_velocity) == 'zero') then
      call compute_grid_velocity(mesh, t, wp)
    else if (trim(ale_grid_velocity) == 'rbf_hybrid') then
      call compute_nodal_velocity_field(mesh, sol, gamma_arr, vp)

      do j = 1, n_rbf
        x_rbf(:, j) = mesh%vert(imp_vert(j))%coord
      end do
      ! Only the bubble's hard points get their target refreshed every
      ! iteration (v_p, the fluid's own nodal velocity there); the piston's
      ! (nimp_bubble+1:nimp) stay at the constant piston_vel set once above.
      do j = 1, nimp_bubble
        delta_imp_rbf(dof_offset(j, nimp) + 1:dof_offset(j, nimp) + 3) = vp(:, imp_vert(j))
      end do

      call compute_rbf_field(nimp, nslide, x_rbf, no_rbf, delta_imp_rbf, rbf_radius, delta_rbf)
      do i = 1, mesh%n_vert
        call eval_rbf_field(n_rbf, nimp, x_rbf, no_rbf, rbf_radius, delta_rbf, mesh%vert(i)%coord, wp(:, i))
      end do
      ! Force hard points exactly (bypass RBF/CG round-off): the bubble
      ! interface moves at v_p, the piston at its own imposed piston_vel
      ! (not v_p, which is just whatever the fluid's nodal solve computed
      ! there and generally differs from the rigid piston's motion).
      do j = 1, nimp_bubble
        wp(:, imp_vert(j)) = vp(:, imp_vert(j))
      end do
      do j = nimp_bubble + 1, nimp
        wp(:, imp_vert(j)) = piston_vel
      end do
      ! Sliding wall points: the RBF fit only *targets* zero velocity along
      ! the local normal (delta_imp=0 there), it doesn't guarantee it
      ! exactly at every evaluation -- CG tolerance and the interpolation
      ! itself leave a small residual normal-component velocity, which
      ! integrates into the node slowly drifting off the wall over many
      ! timesteps. Project it out explicitly so these nodes are mathematically
      ! guaranteed to stay exactly on their wall.
      do j = nimp + 1, n_rbf
        wp(:, imp_vert(j)) = wp(:, imp_vert(j)) &
          - dot_product(wp(:, imp_vert(j)), no_rbf(:, j))*no_rbf(:, j)
      end do
    else
      print*, "Unknown ale_grid_velocity: '", trim(ale_grid_velocity), "'"
      error stop
    end if

    ! In boundary_2d, the mesh is a single flat prism layer: any nonzero
    ! z-component here (the RBF field interpolates a full 3-vector, with
    ! nothing forcing it flat) would slowly warp that layer out of plane,
    ! corrupt the z-normal cap faces, and eventually blow up as a spurious
    ! supersonic z-velocity feeds back into the fluid state. Same reasoning
    ! lagrange_main.F90 documents for its own nodal solver.
    if (boundary_2d) wp(3, :) = 0.0_DOUBLE

    call compute_rhs_ale(mesh, sol, gamma_arr, wp, rhs, sum_lambda, vp)
    call compute_dt_ale(mesh, sum_lambda, cfl, dt)
    ! A crushed (zero/negative-volume) cell eventually turns dt into NaN
    ! (0/0 or a negative sound speed upstream); left unchecked, "t < tmax"
    ! with a NaN t is simply false in IEEE754, so the loop below would exit
    ! silently with a clean "success" exit code and a NaN-filled final vtu,
    ! looking exactly like a normal completion at a much later time. Fail
    ! loudly instead (observed once, coarse mesh, tmax=1.0, t~0.91: dt
    ! collapsed smoothly for ~2000 iterations then finally hit exactly 0).
    if (dt /= dt .or. dt <= 0.0_DOUBLE) then
      print*, "[-] ale_main: invalid dt at iter=", iter, "t=", t, "dt=", dt
      error stop
    end if
    if (t + dt > tmax) dt = tmax - t

    do i = 1, mesh%n_elems
      vol_old(i) = mesh%elem(i)%volume
    end do

    call move_mesh(mesh, wp, dt)
    if (num_procs > 1) call mpi_memory_exchange_vert(mesh, mpi_send_recv)
    call compute_geometry_mesh(mesh, .true., boundary_2d, do_check=.false.)

    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        ! Cell-volume-weighted conservative update: the flux balance rhs
        ! is extensive (area*flux, already signed by compute_rhs_ale), so
        ! forward Euler on vol*sol reads vol^{n+1}*sol^{n+1} = vol^n*sol^n
        ! + dt*rhs. With wp=0 the mesh doesn't move and this reduces to the
        ! plain Eulerian update sol += dt*rhs/vol.
        sol(:, i) = (vol_old(i)*sol(:, i) + dt*rhs(:, i)) / mesh%elem(i)%volume
      end if
    end do
    if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol)

    iter = iter + 1
    t = t + dt
    if (me == 0 .and. mod(iter, n_iter_print) == 0) print*, iter, t, dt

    if (n_iter_write_sol > 0) then
      if (mod(iter, n_iter_write_sol) == 0) then
        write(fln, *) i_sol_vtu
        fln = "output_"//trim(adjustl(fln))
        call open_file_vtu(mesh, trim(fln), fn_vtu, fn_pvtu)
        call write_file_vtu_start_cell_data(mesh, trim(fln), fn_vtu, fn_pvtu)
        do i = 1, mesh%n_elems
          w = conserv_to_primit(sol(:, i), gamma_arr(i))
          density(i) = w(1); velocity(:, i) = w(2:4); pressure(i) = w(5)
          gamma_out(i) = gamma_arr(i)
          centroid(:, i) = mesh%elem(i)%coord
        end do
        call write_file_vtu_cell_scalar(mesh, trim(fln), fn_vtu, fn_pvtu, density, "Density")
        call write_file_vtu_cell_vector(mesh, trim(fln), fn_vtu, fn_pvtu, velocity, "Velocity")
        call write_file_vtu_cell_scalar(mesh, trim(fln), fn_vtu, fn_pvtu, pressure, "Pressure")
        call write_file_vtu_cell_scalar(mesh, trim(fln), fn_vtu, fn_pvtu, gamma_out, "Gamma")
        call write_file_vtu_cell_vector(mesh, trim(fln), fn_vtu, fn_pvtu, centroid, "Centroid")
        call write_file_vtu_end_cell_data(mesh, trim(fln), fn_vtu, fn_pvtu)
        call close_file_vtu(mesh, trim(fln), fn_vtu, fn_pvtu)
        i_sol_vtu = i_sol_vtu + 1
      end if
    end if
  end do

  fln = "output_-1"
  call open_file_vtu(mesh, trim(fln), fn_vtu, fn_pvtu)
  call write_file_vtu_start_cell_data(mesh, trim(fln), fn_vtu, fn_pvtu)
  do i = 1, mesh%n_elems
    w = conserv_to_primit(sol(:, i), gamma_arr(i))
    density(i) = w(1); velocity(:, i) = w(2:4); pressure(i) = w(5)
    gamma_out(i) = gamma_arr(i)
    centroid(:, i) = mesh%elem(i)%coord
  end do
  call write_file_vtu_cell_scalar(mesh, trim(fln), fn_vtu, fn_pvtu, density, "Density")
  call write_file_vtu_cell_vector(mesh, trim(fln), fn_vtu, fn_pvtu, velocity, "Velocity")
  call write_file_vtu_cell_scalar(mesh, trim(fln), fn_vtu, fn_pvtu, pressure, "Pressure")
  call write_file_vtu_cell_scalar(mesh, trim(fln), fn_vtu, fn_pvtu, gamma_out, "Gamma")
  call write_file_vtu_cell_vector(mesh, trim(fln), fn_vtu, fn_pvtu, centroid, "Centroid")
  call write_file_vtu_end_cell_data(mesh, trim(fln), fn_vtu, fn_pvtu)
  call close_file_vtu(mesh, trim(fln), fn_vtu, fn_pvtu)

  call MPI_FINALIZE(mpi_ierr)
end program main
