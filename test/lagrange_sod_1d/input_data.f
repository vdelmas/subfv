&INPUT_PARAM
  meshfile_path=''
  meshfile='sod_1d.msh'

  boundary_2d = .true.

  sol_uniform=1.0, 0., 0., 0., 1.0,
  t_max = 0.14
  cfl = 0.95

  scheme = 'classic'
  second_order = .true.

  init = 0

  n_sol_vtu = 11

  n_bc=0
  /
