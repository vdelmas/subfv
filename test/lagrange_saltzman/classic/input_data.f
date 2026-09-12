&INPUT_PARAM
  meshfile_path=''
  meshfile='saltzman.msh'

  boundary_2d = .true.

  sol_uniform=1.0, 0., 0., 0., 1.e-6,
  t_max = 0.9
  cfl = 0.5

  scheme = 'classic'

  init = 5

  n_sol_vtu = 11

  n_bc=2
  bc_name='right_wall','left_wall'
  bc_type='wall','piston'
  bc_val=1.0, 0., 0., 0., 0.,
         1.0, 1., 0., 0., 0.
  /
