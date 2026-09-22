&INPUT_PARAM
  meshfile_path='/scratch/vdelmas/subfv/test/euler_ho_cylinder_tunnel/mesh16x_mesh/'
  meshfile='cyl_m3_16x.msh'
  boundary_2d = .true.
  init = 0
  sol_uniform = 1.4, 3.0, 0.0, 0.0, 1.0
  gamma_gas = 1.4
  order = 4
  cfl = 0.95
  tmax = 4.0
  flux_scheme = 'modified_three_wave'
  n_bc = 5
  bc_name = 'in_surf', 'out_surf', 'top_surf', 'bot_surf', 'cyl_surf'
  bc_type = 'freestream', 'outflowsupersonic', 'wall', 'wall', 'wall'
  bc_val =  1.4, 3.0, 0.0, 0.0, 1.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0
  aho_method = 1
  use_weno_blend = .true.
  use_cweno_center = .false.
  use_max_lambda_dt = .true.
  n_sol_vtu = 30
  compute_error = .false.
/
