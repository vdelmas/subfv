&INPUT_PARAM
  meshfile_path=''
  meshfile='grad_test.msh'

  boundary_2d = .true.

  sol_uniform=1.0, 0., 0., 0., 1.0,
  t_max = 1e-10
  cfl = 0.5

  scheme = 'classic'
  second_order = .true.

  init = 8

  n_sol_vtu = 2

  n_bc=0
  /
