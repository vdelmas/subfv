! Standalone shock-adaptive node-movement demo/solver.
!
! Single-shot mode (use_iterative_moves=.FALSE., the default): three stages
!   1) march a first-order, fixed-mesh Euler solve to a (quasi-)steady
!      baseline ("before" VTU);
!   2) detect the shock (pressure-jump sensor) and apply a small, bounded,
!      local node displacement so flagged faces relocate toward the shock
!      (X-Mesh-style node movement, minimal/local rather than a full mesh
!      relocation);
!   3) reconverge the same fixed-mesh solve on the now shock-aligned mesh
!      ("after" VTU).
!
! Iterative mode (use_iterative_moves=.TRUE.): after the same Stage-1
! baseline, repeat move->quick-restabilize n_cycles times instead of moving
! once, so the mesh gradually settles onto the shock over many small steps
! rather than one larger one. A VTU snapshot ("output_cycle<N>") is written
! at cycle 1, every 30th cycle, and the final cycle.
!
! Multi-rank MPI is supported: sol is exchanged every iteration (advance),
! and vertex coordinates are exchanged after every node move
! (mpi_memory_exchange_vert, ported from lagrange_module.F90/
! ale_module.F90 -- identical routine in both, send/recv vert%coord for
! shared-element vertices). Needs a pre-partitioned mesh for num_procs>1
! (<meshbase>_<rank>.msh files from subfv-gmsh -part N -part_split
! -part_ghosts), same convention as every other subfv solver.
program main
  use mpi
  use precision_module
  use mpi_module
  use mesh_module
  use mesh_reading_module
  use mesh_geometry_module
  use mesh_connectivity_module
  use io_module

  use shock_adapt_global_data_module
  use shock_adapt_module

  implicit none

  integer(kind=ENTIER) :: me, num_procs, mpi_ierr
  integer(kind=ENTIER) :: icycle
  character(len=255) :: fln_cycle
  type(mesh_type) :: mesh
  type(mpi_send_recv_type) :: mpi_send_recv

  real(kind=DOUBLE), dimension(:, :), allocatable :: sol, rhs
  real(kind=DOUBLE), dimension(:), allocatable :: sum_lambda
  real(kind=DOUBLE), dimension(:), allocatable :: sensor, vert_sensor
  logical, dimension(:), allocatable :: flagged
  real(kind=DOUBLE), dimension(:), allocatable :: node_sensor, node_flagged_out
  logical, dimension(:), allocatable :: node_flagged
  integer(kind=ENTIER), dimension(:), allocatable :: n_neigh
  integer(kind=ENTIER), dimension(:, :), allocatable :: vneigh
  real(kind=DOUBLE), dimension(:, :), allocatable :: disp, cum_disp
  real(kind=DOUBLE), dimension(:), allocatable :: orig_local_scale
  real(kind=DOUBLE) :: t
  integer(kind=ENTIER) :: n_moved
  real(kind=DOUBLE) :: max_disp, vmin
  logical :: do_snapshot
  integer(kind=ENTIER) :: n_moved_g, n_flagged_g
  real(kind=DOUBLE) :: max_disp_g

  call MPI_INIT(mpi_ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

  call read_input_parameters("input_data.f")
  call init_bc_kind()

  call read_mesh_msh(mesh, meshfile_path, meshfile, &
    n_bc, bc_name, me, num_procs, mpi_send_recv)
  call build_mesh(mesh, num_procs, mpi_send_recv, .true., boundary_2d)
  call compute_geometry_mesh(mesh, .true., boundary_2d, me, num_procs)
  if (order >= 2) call setup_wall_mirror(mesh)

  allocate(sol(5, mesh%n_elems))       ; sol = 0.0_DOUBLE
  allocate(rhs(5, mesh%n_elems))       ; rhs = 0.0_DOUBLE
  allocate(sum_lambda(mesh%n_elems))   ; sum_lambda = 0.0_DOUBLE
  allocate(sensor(mesh%n_faces))
  allocate(flagged(mesh%n_faces))
  allocate(vert_sensor(mesh%n_vert))
  allocate(disp(3, mesh%n_vert))
  allocate(node_sensor(mesh%n_vert))
  allocate(node_flagged(mesh%n_vert))
  allocate(node_flagged_out(mesh%n_vert))

  call init_sol(mesh, sol)

  ! --- Stage 1: march to a (quasi-)steady baseline ---
  t = 0.0_DOUBLE
  call advance(n_iter_steady, "stage1")
  call write_snapshot("output_before")

  if (detect_only) then
    if (me == 0) print*, "[detect_only] nodes flagged =", count(node_flagged)
    call MPI_FINALIZE(mpi_ierr)
    stop
  end if

  if (use_curvature_move) then
    ! --- Curvature-based move, repeated n_cycles times (>=2 -- a single
    ! pass only nudges nodes part-way per curvature_relax; several passes
    ! let the front actually settle, each one re-detecting on the
    ! restabilized solution). ---
    allocate(n_neigh(mesh%n_vert))
    allocate(vneigh(16, mesh%n_vert))
    call build_vert_adjacency(mesh, n_neigh, vneigh)

    do icycle = 1, n_cycles
      if (me == 0) print*, "[curvature] cycle", icycle, &
        " min elem volume before move =", min_elem_volume(mesh)

      call compute_shock_sensor_grad(mesh, sol, node_sensor, node_flagged)
      call compute_node_displacement_curvature(mesh, node_flagged, n_neigh, vneigh, &
        disp, n_moved, max_disp)

      call move_mesh(mesh, disp)
      if (num_procs > 1) call mpi_memory_exchange_vert(mesh, mpi_send_recv)
      call compute_geometry_mesh(mesh, .true., boundary_2d, do_check=.false.)
      call check_volume("curvature")
      ! Topology is unchanged by a move, but vertex positions did, and the
      ! adjacency itself only depends on connectivity (not coordinates) --
      ! no need to rebuild n_neigh/vneigh between cycles.

      call advance(n_iter_restab, "stage3")

      call global_stats(count(node_flagged), n_moved, max_disp)
      if (me == 0) print*, "[curvature] cycle", icycle, " nodes flagged =", n_flagged_g, &
        " vertices moved =", n_moved_g, " max displacement =", max_disp_g, " t=", t

      write(fln_cycle, '(A,I0)') "output_cycle", icycle
      call write_snapshot(trim(fln_cycle))
    end do
    call write_snapshot("output_after")
  else if (.not. use_iterative_moves) then
    ! --- Single-shot Stage 2: bounded node movement onto the shock ---
    if (me == 0) print*, "[stage2] min elem volume before move =", min_elem_volume(mesh)

    call compute_shock_sensor(mesh, sol, sensor, flagged)
    call compute_node_displacement(mesh, sol, flagged, disp, n_moved, max_disp)
    call global_stats(0_ENTIER, n_moved, max_disp)
    if (me == 0) print*, "[stage2] vertices moved =", n_moved_g, " max displacement =", max_disp_g

    call move_mesh(mesh, disp)
    if (num_procs > 1) call mpi_memory_exchange_vert(mesh, mpi_send_recv)
    call compute_geometry_mesh(mesh, .true., boundary_2d, do_check=.false.)
    call check_volume("stage2")

    ! --- Stage 3: reconvergence on the shock-aligned mesh ---
    call advance(n_iter_post_move, "stage3")
    call write_snapshot("output_after")
  else
    ! --- Iterative mode: many small move+restabilize cycles ---
    ! Local scale is fixed at its Stage-1 (pre-movement) value, and the
    ! per-vertex CUMULATIVE displacement (not just this cycle's increment)
    ! is capped against it -- see compute_node_displacement's header
    ! comment for why capping against the current, progressively-shrunk
    ! geometry each cycle is NOT safe under repetition (confirmed: it let
    ! a cell collapse to negative volume by cycle 4 of an early test).
    allocate(orig_local_scale(mesh%n_vert))
    allocate(cum_disp(3, mesh%n_vert)); cum_disp = 0.0_DOUBLE
    call compute_local_scale(mesh, orig_local_scale)

    do icycle = 1, n_cycles
      call compute_shock_sensor(mesh, sol, sensor, flagged)
      call compute_node_displacement(mesh, sol, flagged, disp, n_moved, max_disp, &
        orig_local_scale, cum_disp)
      call move_mesh(mesh, disp)
      if (num_procs > 1) call mpi_memory_exchange_vert(mesh, mpi_send_recv)
      call compute_geometry_mesh(mesh, .true., boundary_2d, do_check=.false.)
      call check_volume("cycle")

      call advance(n_iter_restab, "restab")

      call global_stats(0_ENTIER, n_moved, max_disp)
      do_snapshot = (icycle == 1) .or. (mod(icycle, 30) == 0) .or. (icycle == n_cycles)
      if (me == 0) print*, "[cycle]", icycle, "vertices moved =", n_moved_g, &
        " max displacement =", max_disp_g, " t=", t, " snapshot=", do_snapshot
      if (do_snapshot) then
        write(fln_cycle, '(A,I0)') "output_cycle", icycle
        call write_snapshot(trim(fln_cycle))
      end if
    end do
  end if

  call MPI_FINALIZE(mpi_ierr)

contains

  ! Sums n_moved/count(flagged) and max-reduces max_disp across ranks, so
  ! the printed diagnostics are meaningful global totals when num_procs>1
  ! (each rank otherwise only sees its own local partition's vertices --
  ! confirmed with a local mpirun -np 4 test that the un-reduced counts
  ! differed sharply from the serial run's, as expected). Ghost vertices
  ! shared at a partition boundary can be double-counted in the SUM here
  ! (n_moved/n_flagged) -- a known, minor imprecision, same convention as
  ! other diagnostic prints elsewhere in this codebase (e.g. ale_main.F90's
  ! RBF-hybrid vertex counts), not worth an exact dedup for a print line.
  subroutine global_stats(n_flagged_local, n_moved_local, max_disp_local)
    implicit none
    integer(kind=ENTIER), intent(in) :: n_flagged_local, n_moved_local
    real(kind=DOUBLE), intent(in) :: max_disp_local
    integer :: ierr

    n_flagged_g = n_flagged_local
    n_moved_g = n_moved_local
    max_disp_g = max_disp_local
    if (num_procs > 1) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, n_flagged_g, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, n_moved_g, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, max_disp_g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, ierr)
    end if
  end subroutine global_stats

  ! Explicit forward-Euler time march for n_iter iterations on the current
  ! (fixed) mesh, updating the module-level t/sol in place.
  subroutine advance(n_iter, tag)
    implicit none
    integer(kind=ENTIER), intent(in) :: n_iter
    character(len=*), intent(in) :: tag

    integer(kind=ENTIER) :: iter, i
    real(kind=DOUBLE) :: dt

    do iter = 1, n_iter
      if (trim(tag) == "stage1" .and. iter <= order_ramp_iter) then
        call compute_rhs(mesh, sol, rhs, sum_lambda, order_override=1_ENTIER)
      else
        call compute_rhs(mesh, sol, rhs, sum_lambda)
      end if
      call compute_dt(mesh, sum_lambda, cfl, dt)
      if (dt /= dt .or. dt <= 0.0_DOUBLE) then
        print*, "[-] shock_adapt_main (", trim(tag), "): invalid dt at iter=", iter, "dt=", dt
        error stop
      end if

      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) sol(:, i) = sol(:, i) + dt*rhs(:, i)/mesh%elem(i)%volume
      end do
      if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5_ENTIER, sol)

      t = t + dt
      if (me == 0 .and. n_iter_print > 0 .and. mod(iter, n_iter_print) == 0) &
        print*, "["//trim(tag)//"]", iter, t, dt
    end do
  end subroutine advance

  ! Global MIN reduce before deciding to abort: check_volume is called by
  ! every rank (not gated by me==0), so without this a rank whose OWN local
  ! partition happens to have the degenerate cell would error stop alone
  ! while the others carry on into the next collective call (compute_rhs's
  ! MPI exchanges, compute_dt's ALLREDUCE, ...) and hang waiting on a
  ! message that never comes from the now-dead rank -- every rank must
  ! reach the same abort/continue decision together.
  subroutine check_volume(tag)
    implicit none
    character(len=*), intent(in) :: tag
    real(kind=DOUBLE) :: v
    integer :: ierr

    v = min_elem_volume(mesh)
    if (num_procs > 1) call MPI_ALLREDUCE(MPI_IN_PLACE, v, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
    if (v <= 0.0_DOUBLE) then
      print*, "[-] shock_adapt_main (", trim(tag), "): non-positive cell volume after node movement."
      error stop
    end if
  end subroutine check_volume

  subroutine write_snapshot(fln)
    implicit none
    character(len=*), intent(in) :: fln

    integer(kind=ENTIER) :: fn_vtu, fn_pvtu, i
    real(kind=DOUBLE), dimension(5) :: w
    real(kind=DOUBLE), dimension(mesh%n_elems) :: density, pressure
    real(kind=DOUBLE), dimension(3, mesh%n_elems) :: velocity

    call compute_shock_sensor(mesh, sol, sensor, flagged)
    call compute_vert_sensor(mesh, sensor, vert_sensor)
    call compute_shock_sensor_grad(mesh, sol, node_sensor, node_flagged)
    node_flagged_out = merge(1.0_DOUBLE, 0.0_DOUBLE, node_flagged)

    call open_file_vtu(mesh, fln, fn_vtu, fn_pvtu)
    call write_file_vtu_start_cell_data(mesh, fln, fn_vtu, fn_pvtu)
    do i = 1, mesh%n_elems
      w = conserv_to_primit(sol(:, i))
      density(i) = w(1); velocity(:, i) = w(2:4); pressure(i) = w(5)
    end do
    call write_file_vtu_cell_scalar(mesh, fln, fn_vtu, fn_pvtu, density, "Density")
    call write_file_vtu_cell_vector(mesh, fln, fn_vtu, fn_pvtu, velocity, "Velocity")
    call write_file_vtu_cell_scalar(mesh, fln, fn_vtu, fn_pvtu, pressure, "Pressure")
    call write_file_vtu_end_cell_data(mesh, fln, fn_vtu, fn_pvtu)
    call write_file_vtu_start_vert_data(mesh, fln, fn_vtu, fn_pvtu)
    call write_file_vtu_vert_scalar(mesh, fln, fn_vtu, fn_pvtu, vert_sensor, "ShockSensor")
    call write_file_vtu_vert_scalar(mesh, fln, fn_vtu, fn_pvtu, node_sensor, "NodalShockSensor")
    call write_file_vtu_vert_scalar(mesh, fln, fn_vtu, fn_pvtu, node_flagged_out, "NodalFlagged")
    call write_file_vtu_end_vert_data(mesh, fln, fn_vtu, fn_pvtu)
    call close_file_vtu(mesh, fln, fn_vtu, fn_pvtu)
  end subroutine write_snapshot

end program main
