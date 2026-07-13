&INPUT_PARAM

  meshfile_path='',
  meshfile='spherical_sedov_hex.msh'

  use_cylinder_map = .false.
  dim_cylinder_map = 2

  n_iter_print=10
  n_iter_write_sol=100000
  n_max_iter=100000
  final_res=1e-12

  second_order = .false.
  method = 0
  activate_diffusion = .false.

  local_time_step=.false.
  cfl=0.8
  tmax=1.0

  use_sutherland = .true.
  mu0 = 1.789e-5
  T0 = 288.

  init_uniform=.false.,
  sol_uniform=1e-3, 5000., 0., 0., 57.615576818289696

  init_gresho=.false.,

  init_sedov=.true.,
  r_sedov=0.04,

  n_bc=0,
  bc_name='in_surf','out_surf','cyl_surf',
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
  boundary_2d=.false.
  exclude_bound_vert=.false.

  scheme='ZB_ARMD_LS',

  write_residual=.false.
  n_iter_residual=10

  compute_coeffs = .false.,
  coeffs_surf = 3,

  rhoinf = 1e-3,
  vinf = 5000.,
  pinf = 57.615576818289696,
  /
