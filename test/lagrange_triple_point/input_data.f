&INPUT_PARAM
  meshfile_path=''
  meshfile='clean2_triple_point.msh'

  boundary_2d = .true.

  sol_uniform=1.0, 0., 0., 0., 2.0,
  t_max = 0.24
  cfl = 0.2

  scheme = 'sidil'

  init = 7

  n_sol_vtu = 11

  n_bc=0
  /
