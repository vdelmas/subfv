&INPUT_PARAM

  meshfile_path='/run/media/delmasv/scratch/Codes/subfv/subfv-ns/test/test_error_heat_flux_hybrid/meshes/0.2_2.0000e-05/',
  meshfile='half_cylinder_tri_ns_error_heat.msh'

  use_cylinder_map = .false.
  dim_cylinder_map = 2

  n_iter_print=1
  n_iter_write_sol=999999
  n_max_iter=10
  final_res=1e-12

  second_order = .false.
  method = 1,
  activate_diffusion = .true.

  use_sutherland = .true.
  mu0 = 1.789e-5
  T0 = 288.

  init_uniform=.true.,
  !sol_uniform=1e-3, 5000., 0., 0., 57.615576818289696
  sol_uniform=1e-3, 0., 0., 0., 57.615576818289696

  n_bc=3,
  !bc_name='in_surf','out_surf','cyl_surf',
  bc_name='surf_in','surf_out','surf_cyl',
  bc_type='freestream', 'outflowsupersonic','wall', 
  bc_val=1e-3,5000.,0.,0., 57.615576818289696
  0., 0., 0., 0., 0.
  0., 0., 0., 0., 0.

  bc_type_V='', '','dirichlet', 
  bc_val_V=0., 0., 0.
  0., 0., 0.,
  0., 0., 0.,

  bc_type_T='', '','dirichlet', 
  bc_val_T=0., 0., 500.

  bc_style = 2,
  boundary_2d=.true.
  error_2d_h=1e-3,
  exclude_bound_vert=.false.

  scheme='multi_point_iso',

  write_residual=.true.
  n_iter_residual=10

  compute_coeffs = .true.,
  coeffs_surf = 3,

  rhoinf = 1e-3,
  vinf = 5000.,
  pinf = 57.615576818289696,
  /
