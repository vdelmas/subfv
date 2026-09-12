&INPUT_PARAM
  meshfile_path=''
  meshfile='sedov_2d_quad.msh'

  boundary_2d = .true.

  sol_uniform=1.0, 0., 0., 0., 1.0,
  t_max = 1.0
  cfl = 0.5

  scheme = 'classic'

  init = 1

  n_sol_vtu = 51

  n_bc=0
  /
