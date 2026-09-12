&INPUT_PARAM
  meshfile_path=''
  meshfile='noh_2d.msh'

  boundary_2d = .true.

  sol_uniform=1.0, 0., 0., 0., 0.5,
  t_max = 0.6
  cfl = 0.5

  scheme = 'classic'
  second_order = .true.

  init = 3

  n_sol_vtu = 51

  n_bc=1
  bc_name='movingbound'
  bc_type='pressure'
  bc_val=1.e-6, 0., 0., 0., 0.
  /
