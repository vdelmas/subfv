&INPUT_PARAM
  meshfile_path=''
  meshfile='sod_bimat.msh'

  boundary_2d = .true.

  sol_uniform=1.0, 0., 0., 0., 1.0,
  t_max = 0.14
  cfl = 0.5

  scheme = 'classic'
  second_order = .true.

  init = 6

  n_sol_vtu = 11

  n_bc=0
  /
