&INPUT_PARAM
  meshfile_path=''
  meshfile='cyl_m3_fine.msh'
  boundary_2d = .true.
  init = 0
  sol_uniform = 1.4, 3.0, 0.0, 0.0, 1.0
  gamma_gas = 1.4
  order = 2
  cfl = 0.95
  tmax = 4.0
  flux_scheme = 'three_wave'
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
  eps_weight_num_gg = 1.0e-6
  use_max_lambda_dt = .true.
  n_sol_vtu = 10
  compute_error = .false.
/
