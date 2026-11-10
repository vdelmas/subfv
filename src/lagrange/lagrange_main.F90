program main
  use mpi
  use precision_module
  use mpi_module
  use mesh_module
  use mesh_reading_module
  use mesh_geometry_module
  use mesh_connectivity_module
  use lagrange_module
  use lagrange_io_module
  implicit none

  integer(kind=ENTIER) :: fn
  integer :: mpi_ierr, me, num_procs
  integer :: inan_vp
  type(mpi_send_recv_type) :: mpi_send_recv
  character(len=255) :: meshfile, meshfile_path

  character(len=255) :: fln

  integer(kind=ENTIER), parameter :: n_max_bc = 10
  integer(kind=ENTIER) :: n_bc, iter, i, j, method_length = 0
  integer(kind=ENTIER) :: id_sub_face, id_face, piston_bc_idx
  real(kind=DOUBLE) :: cfl, cfl_max = 0.8, b2d_h = 1.0
  real(kind=DOUBLE) :: h_extrude = 1.0
  integer(kind=ENTIER) :: init = 0
  ! Diagnostic: drop boundary vertices from the o2 nodal gradient fit, as the Euler
  ! code does, instead of mirroring. Default .false. = behaviour unchanged.
  logical :: grad_bound_skip = .false.
  ! .true. : gradient exact pour un champ lineaire aux noeuds de mur (voir compute_gradients).
  logical :: grad_bound_interior = .false.
  character(len=255) :: scheme = ""
  integer(kind=ENTIER) :: iso_k
  ! selecteur de base du flux d'energie, extrait du nom de schema et MIS EN
  ! MAJUSCULE : les codes a lettre (WA, WB) sont ecrits en minuscule dans les
  ! noms de cas et pos_base_mode les compare en majuscule.
  character(len=1) :: base_sel
  character(len=255), dimension(n_max_bc) :: bc_name
  character(len=255), dimension(n_max_bc) :: bc_type
  real(kind=DOUBLE), dimension(5, n_max_bc) :: bc_val
  real(kind=DOUBLE), dimension(n_max_bc) :: bc_p
  real(kind=DOUBLE), dimension(5) :: sol_uniform
  real(kind=DOUBLE), dimension(3) :: piston_vel
  type(mesh_type) :: mesh

  logical :: boundary_2d = .true.

  real(kind=DOUBLE) :: pvp_bnd_abs_work = 0.0_DOUBLE
  real(kind=DOUBLE) :: odd_even_amp = 0.0_DOUBLE
  real(kind=DOUBLE) :: t, t_max, dt
  real(kind=DOUBLE), dimension(:,:), allocatable :: vp
  integer(kind=ENTIER) :: fdiag_c1 = -1, fdiag_c2 = -1, fdiag_c3 = -1
  real(kind=DOUBLE), dimension(:), allocatable :: mass, pp, gamma_arr
  real(kind=DOUBLE), dimension(:,:), allocatable :: sol
  real(kind=DOUBLE), dimension(:,:), allocatable :: rhs
  real(kind=DOUBLE), dimension(:,:), allocatable :: new_sol
  logical, dimension(:), allocatable :: vp_is_imposed
  ! Second-order reconstruction: per-cell grad_v (3x3), grad_p (3), div_v (1)
  logical :: second_order = .false.
  ! MUSCL-Hancock coherent : meme dt pour le predicteur et le correcteur.
  ! .false. restaure l'ancien decalage (predicteur sur dt^{n-1}, 0 au pas 1).
  logical :: dt_at_step_start = .true.
  ! Integration en temps. 'muscl_hancock' (DEFAUT) = extrapolation dt/2 dans la
  ! reconstruction. 'rk2' = reconstruction SPATIALE seule + maillage au demi-pas.
  character(len=20) :: time_scheme = 'muscl_hancock'
  logical :: rk2_time
  integer(kind=ENTIER) :: istage, nstage
  real(kind=DOUBLE), dimension(:,:), allocatable :: coord_n, sol_n
  ! Cadence des impressions PAR PAS. 1 = toutes (comportement historique).
  integer(kind=ENTIER) :: trace_every = 1
  logical :: do_trace
  ! Controle a posteriori du retard d'un pas sur la borne de volume.
  real(kind=DOUBLE) :: dt_post, dtviol_max = 0.0_DOUBLE
  integer(kind=ENTIER_D) :: n_dtviol = 0
  ! Iteration du predicteur : 2 = un seul re-predicteur. 1 = pas d'iteration.
  integer(kind=ENTIER) :: dt_n_iter = 2, it_dt
  integer(kind=ENTIER_D) :: n_dtiter2 = 0
  real(kind=DOUBLE) :: dtiter2_max = 0.0_DOUBLE
  real(kind=DOUBLE), dimension(:,:,:), allocatable :: grad_v
  real(kind=DOUBLE), dimension(:,:), allocatable :: grad_p
  real(kind=DOUBLE), dimension(:), allocatable :: div_v
  real(kind=DOUBLE), dimension(:), allocatable :: alpha_p_arr
  integer(kind=ENTIER), dimension(:), allocatable :: p_limited
  real(kind=DOUBLE), dimension(:), allocatable :: h_p_arr

  integer(kind=ENTIER) :: n_sol_vtu=2
  integer(kind=ENTIER) :: i_sol_vtu
  integer(kind=ENTIER) :: iso_weight_mode = 2   ! w_pcf recipe for classic_iso_vol
  ! Isentropic-vortex (init=2) convergence diagnostics
  logical :: compute_error = .false.
  logical :: grad_weno = .true.
  logical :: grad_gg_corrected = .false.
  real(kind=DOUBLE) :: grad_eps_weno = 0.0_DOUBLE   ! <=0 -> tiny(), the historical value
  real(kind=DOUBLE) :: h_err, l1err, l2err, linferr

  namelist /INPUT_PARAM/ &
    meshfile_path, meshfile, &
    t_max, cfl, &
    init, &
    grad_bound_skip, grad_bound_interior, &
    lambda_acoustic_only, diag_entropy, &
    vp_ep_wip_mode, vp_ep_wip_floor, vp_ep_wip_clip, &
    pvp_eps, pvp_lim_iter, pvp_energy_mode, odd_even_amp, pos_base_mode, pos_w_alpha, &
    theta_pi, theta_m_floor, mp_cvisc, diag_divratio, diag_firstprint, &
    md_cvisc, md_first, md_trrs, md_noswitch, md_trace, md_watch, &
    kap_floor, kap_one, kap_check, kap_front_jump, kap_fixed, diag_iso, &
    pif_vp_classic, pif_duk, pif_duk_it, &
    kap_loc, kap_loc_c, kap_loc_mode, diag_dtb, &
    diag_force, fdiag_c1, fdiag_c2, fdiag_c3, &
    diag_area, diag_lamratio, diag_alpha, &
    mdiv_theta, &
    psvp_lam_iter, psvp_lam_marge, &
    fbd_lam_mode, &
    divd_coef, divd_mode, divd_theta, divd_nl, &
    divd_wdbg, divd_wnode, divd_wmax, &
    sedov_nodal_deposit, sedov_energy, &
    tang_mode, pi_bord_psw, convex_limiting, &
    scheme, method_length, b2d_h, &
    n_bc, bc_name, bc_type, bc_val, &
    sol_uniform, boundary_2d, &
    n_sol_vtu, second_order, iso_weight_mode, compute_error, grad_weno, grad_gg_corrected, grad_eps_weno, &
    gresho_mach, gamma_uniform, vac_vortex, vac_wave, dt_at_step_start, dt_n_iter, trace_every, time_scheme, grad_weight_v, &
    wall_mode, grad_eps_vol, vp_gcorr

  call MPI_INIT(mpi_ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

  open(newunit=fn, file="input_data.f")
  read(unit=fn, nml=INPUT_PARAM)
  close(fn)

  ! "classic_rhoa" is the EUCCLHYD-original variant of "classic": the sub-face
  ! swept-mass coefficient uses the pure acoustic impedance lambda = rho*a,
  ! without the pressure-/velocity-jump (Dukowicz-like) terms. Selecting it by
  ! scheme name rather than by the lambda_acoustic_only namelist flag alone
  ! keeps the two variants distinguishable in a campaign's inputs and logs.
  if (scheme == "classic_rhoa") lambda_acoustic_only = .true.
  if (scheme == "vp_pp_pvp_e") pvp_energy_mode = 'bal'
  if (scheme == "vp_pp_pos") pvp_pos = .true.
  ! 'vp_pp_psvp' : quantite de mouvement de vp_pp (p_p nodale), energie de
  ! classic (Pi_c unilaterale). Les trois conservations tiennent, la production
  ! d'entropie n'est plus une somme de carres -- voir pi_psvp dans le module.
  if (scheme == "vp_pp_psvp") pi_psvp = .true.
  ! 'vp_pp_ep' : flux ( -v_p.n , p_p n , (pv)*_s ). Voir pi_pp_ep.
  if (scheme == "vp_pp_ep") pi_pp_ep = .true.
  ! variante a energie INTERNE, voir pi_ep_interne.
  if (scheme == "vp_pp_ep_i") then
    pi_pp_ep = .true.
    pi_ep_interne = .true.
  end if
  if (scheme == "vp_ep_i") pi_ep_interne = .true.
  ! 'vp_pp_epsp' : energie INTERNE au noeud, forme canonique (eps*_c = eps_p).
  ! Voir pi_pp_epsp dans le module : a ne pas confondre avec vp_pp_ep_i, qui
  ! moyenne eps_c sans la correction de glissement.
  if (scheme == "vp_pp_epsp") then
    pi_pp_ep = .true.
    pi_pp_epsp = .true.
  end if
  ! 'vp_pp_eps2' / 'vp_pp_eps2m' : deux eps* differents, kappa nodal unique.
  ! Voir pi_pp_eps2 dans le module. ep_node y porte kappa, pas une energie.
  if (scheme == "vp_pp_eps2") then
    pi_pp_ep = .true.
    pi_pp_eps2 = .true.
  end if
  if (scheme == "vp_pp_eps2m") then
    pi_pp_ep = .true.
    pi_pp_eps2m = .true.
  end if
  ! 'vp_pp_epv' : flux d'energie = pression nodale x vitesse etoile de chaque
  ! cote, sans energie nodale. Conservation et positivite automatiques, voir
  ! pi_pp_epv dans le module.
  if (scheme == "vp_pp_epv") then
    pi_pp_ep = .true.
    pi_pp_epv = .true.
  end if
  ! 'vp_pp_vpps' : flux d energie = pression de sous-face CLASSIQUE Pi_c x
  ! vitesse NODALE, alors que le flux de quantite de mouvement garde p_p.
  ! Combinaison volontairement INCOHERENTE : la pression du flux d energie
  ! differe de celle du flux de QdM, donc eps* depend du REFERENTIEL (au taux
  ! (Pi_c - p_p)/lam_c) et l energie n est PLUS CONSERVEE, puisque l equilibre
  ! nodal porte sur p_p et non sur Pi_c. Le diagnostic VPPSTEL et la derive
  ! d energie totale sont obligatoires pour l interpreter.
  if (scheme == "vp_pp_vpps") then
    pi_pp_ep = .true.
    pi_pp_vpps = .true.
  end if
  ! 'vp_pp_ppvn' : flux d energie = p_p (v_p.n). Pression NODALE des deux cotes
  ! et vitesse NODALE : conservatif trivialement (F_l = F_r), galileen, et sans
  ! terme lineaire en la vitesse de maille dans eps*.
  if (scheme == "vp_pp_ppvn") then
    pi_pp_ep = .true.
    pi_pp_ppvn = .true.
  end if
  ! 'vp_pp_ep_wdiv' : vp_pp_ep dont la pression nodale pondere les sauts de
  ! vitesse en 1/lambda. Voir pi_pp_wdiv dans le module.
  if (scheme == "vp_pp_ep_wdiv") then
    pi_pp_ep = .true.
    pi_pp_wdiv = .true.
  end if
  ! 'vp_pp_ppvn_wdiv' : flux d energie p_p (v_p.n), AUCUNE energie nodale, et
  ! pression nodale fermee par MOINDRES CARRES. Avec vp_pp_ep la relaxation
  ! d energie est bien plus grosse que l effet de la ponderation et le masque ;
  ! ici elle est absente, donc la ponderation est la SEULE source de dissipation.
  if (scheme == "vp_pp_ppvn_wdiv") then
    pi_pp_ep = .true.
    pi_pp_ppvn = .true.
    pi_pp_wdiv = .true.
  end if
  ! 'classic_conv' : `classic` dont l IMPEDANCE DE STABILISATION est convective,
  ! mu_c = rho_c |v_p.n - v_c.n|, au lieu de rho_c a_c. Voir pi_conv_stab.
  if (scheme == "classic_conv") pi_conv_stab = .true.
  ! 'vp_pp_pivp' : QdM en p_p (bas Mach, aucune trainee), energie en
  ! Pi_c (v_p.n) + relaxation nodale (dissipation de classic). Voir pi_pivp.
  if (scheme == "vp_pp_pivp") then
    pi_pp_ep = .true.
    pi_pivp = .true.
  end if
  ! 'vp_pp_ppus' : le flux d energie INVARIANT par changement de referentiel,
  ! (pv)*_c = p_p (u*_c.n) -/+ mu_c (eps_p - eps_c). Voir pi_pp_ppus.
  if (scheme == "vp_pp_ppus") then
    pi_pp_ep = .true.
    pi_ep_interne = .true.
    pi_pp_ppus = .true.
  end if
  ! 'vp_pp_ep_veps' : vp_pp_ep dont la pression nodale est fermee par
  ! somme A (eps*_r - eps*_l) = 0. Voir pi_pp_veps dans le module.
  if (scheme == "vp_pp_ep_veps") then
    pi_pp_ep = .true.
    pi_pp_veps = .true.
  end if
  ! 'vp_hat' : vp dont la pression de sous-face est mixee avec la pression nodale,
  ! LE MEME objet dans les deux flux. theta_pi = 0 doit redonner vp au bit pres.
  if (scheme == "vp_hat") pi_hat = .true.
  ! moyenne de FACE au lieu de la pression nodale : beta = theta_pi.
  ! Pi = p_bar_f - mu_p (div v)_p, le meme dans les deux flux ; v_p vient du
  ! solveur existant et n'est pas determine par ce flux.
  if (scheme == "mean_div") pi_meandiv = .true.
  if (scheme == "vp_kappa") pi_kappa = .true.
  ! 'vp_fbd' : Pi = pbar_f (ponderee par l'impedance, LOCALE a la sous-face)
  ! moins D_p (NODALE, homogene a une pression sans aucune longueur a choisir).
  ! Meme forme que 'vp_mp' mais sans h_p ni coefficient libre, et avec une moyenne
  ! ponderee au lieu de l'arithmetique -- c'est ce qui rend la limite 1D exacte.
  if (scheme == "vp_fbd") pi_fbd = .true.
  if (scheme == "vp_fbw") pi_fbw = .true.
  ! 'vp_mdiv' : la fermeture fbw corrigee par la dilatation du SOUS-ELEMENT.
  ! Elle a besoin de `delta_pc`, donc du chemin `compute_kappa` qui le calcule et
  ! de son echange MPI : on active donc `pi_kappa` avec, mais `kap_one` pour que
  ! le kappa lui-meme reste a 1 et ne modifie pas le flux -- seul `mdv_se` sert.
  if (scheme == "vp_mdiv") then
    pi_fbw = .true.
    pi_mdiv = .true.
    pi_kappa = .true.
    kap_one = .true.
  end if
  if (scheme == "vp_div") pi_divdirect = .true.
  if (scheme == "vp_pif") pi_pif = .true.
  if (scheme == "vp_mp") then
    pi_hat = .true.
    pi_mp = .true.
  end if
  if (scheme == "vp_facebar") then
    pi_hat = .true.
    pi_face = .true.
  end if
  ! 'vp_hat_md' : idem, plus l'amortissement bas Mach theta_M = min(Ma,1) de la
  ! seule part visqueuse. theta_g reste theta_pi.
  if (scheme == "vp_hat_md") then
    pi_hat = .true.
    pi_mach_damp = .true.
  end if
  if (scheme(1:min(11,len_trim(scheme))) == "vp_pp_pos_v") then
    pvp_pos = .true.
    base_sel = scheme(12:12)
    if (base_sel >= 'a' .and. base_sel <= 'z') base_sel = achar(iachar(base_sel)-32)
    pos_base_mode = "V"//base_sel
  end if
  ! meme famille W, mais limiteur ACTIF : c'est la mesure de cavite qui impose
  ! de pouvoir croiser base et limiteur, les deux etant jusqu'ici confondus.
  if (scheme(1:min(11,len_trim(scheme))) == "vp_pp_pos_w") then
    pvp_pos = .true.
    base_sel = scheme(12:12)
    if (base_sel >= 'a' .and. base_sel <= 'z') base_sel = achar(iachar(base_sel)-32)
    pos_base_mode = "W"//base_sel
  end if
  ! mode pur : meme cible, limiteur entierement court-circuite
  if (scheme(1:min(11,len_trim(scheme))) == "vp_pp_pur_v") then
    pvp_pos = .true.
    pos_no_limit = .true.
    base_sel = scheme(12:12)
    if (base_sel >= 'a' .and. base_sel <= 'z') base_sel = achar(iachar(base_sel)-32)
    pos_base_mode = "V"//base_sel
  end if
  if (scheme(1:min(11,len_trim(scheme))) == "vp_pp_pur_w") then
    pvp_pos = .true.
    pos_no_limit = .true.
    base_sel = scheme(12:12)
    if (base_sel >= 'a' .and. base_sel <= 'z') base_sel = achar(iachar(base_sel)-32)
    pos_base_mode = "W"//base_sel
  end if
  ! 'vp_pp_pos_sv' : la variante UNIVALUEE de l'ecretage. Le flux garde une seule
  ! valeur par sous-face, donc la conservation est intacte et il n'y a RIEN a
  ! redistribuer, contrairement aux familles V et W ci-dessus.
  ! 'vp_pp_pos_sva' : la meme, mais toujours le milieu, pour mesurer ce que coute
  ! le milieu systematique la ou le flux naturel etait deja admissible.
  if (scheme == "vp_pp_pos_sv" .or. scheme == "vp_pp_pos_sva") then
    pvp_pos = .true.
    pvp_pos_sv = .true.
    pvp_pos_sv_always = (scheme == "vp_pp_pos_sva")
  end if
  if (scheme == "vp_ep_wip") vp_ep_wip = .true.
  ! Le limiteur convexe (pvp_limit) existe dans lagrange_module mais n'est PAS
  ! atteignable par un nom de schema : il n'a pas ete valide sur la campagne
  ! complete. Pour l'activer il faut reintroduire ici un nom dedie.
  if (pvp_pos_sv .and. me == 0) call pos_sv_selfcheck()
  if (pi_fbd .and. me == 0) call fbd_selfcheck()
  ! Controle du mode sablier : il ne depend d'aucun etat physique, seulement
  ! de la geometrie et de la structure de la fermeture, donc il est lance
  ! pour fbd ET pour fbw (ou D_p est le meme scalaire par noeud).
  if ((pi_fbd .or. pi_fbw) .and. me == 0) call hourglass_selfcheck()
  ! Plafond de positivite du terme de divergence de sous-element : ne depend
  ! d'aucun etat du maillage, donc lance des que la fermeture peut l'utiliser.
  if ((pi_fbd .or. pi_fbw .or. pi_kappa) .and. me == 0) call md_posflux_selfcheck()
  if (me == 0) then
    if (lambda_acoustic_only) then
      print*, "lambda = rho*a (EUCCLHYD original)"
    else
      print*, "lambda = max(rho*a, pressure-jump, velocity-jump)"
    end if
  end if

  call read_mesh_msh(mesh, meshfile_path, meshfile, &
    n_bc, bc_name, me, num_procs, mpi_send_recv)
  call build_mesh(mesh, num_procs, mpi_send_recv, &
    .true., boundary_2d)

  ! Odd-even decoupling : perturbation de la LIGNE CENTRALE du maillage a la
  ! frequence de Nyquist, c'est-a-dire en alternant le signe d'une maille a
  ! l'autre le long de x. C'est le mode que le test doit sonder -- celui ou les
  ! mailles paires et impaires se decouplent, d'ou le nom. La subroutine odd_even
  ! de mesh_geometry_module, jamais appelee, perturbe en helice sur theta =
  ! x*pi/2, soit une periode de 4 mailles : quatre fois trop lente, et un schema
  ! peut l'amortir tout en etant instable sur le vrai mode.
  ! odd_even_amp est une FRACTION de la taille de maille, pour que le reglage ne
  ! depende pas du maillage. Trop petite elle reste dans le bruit d'arrondi et
  ! rien ne pousse ; trop grande elle devient un vrai defaut geometrique et ce
  ! qu'on observe n'est plus une instabilite mais une condition initiale.
  if (odd_even_amp > 0.0_DOUBLE) then
    block
      real(kind=DOUBLE) :: xmin, xmax, ymid, zmid, hx, dx
      integer(kind=ENTIER) :: iv, k, n_moved
      xmin = minval(mesh%vert(1:mesh%n_vert)%coord(1))
      xmax = maxval(mesh%vert(1:mesh%n_vert)%coord(1))
      ymid = 0.5_DOUBLE*(minval(mesh%vert(1:mesh%n_vert)%coord(2)) &
                       + maxval(mesh%vert(1:mesh%n_vert)%coord(2)))
      zmid = 0.5_DOUBLE*(minval(mesh%vert(1:mesh%n_vert)%coord(3)) &
                       + maxval(mesh%vert(1:mesh%n_vert)%coord(3)))
      ! pas en x : plus petit ecart strictement positif entre deux sommets
      hx = xmax - xmin
      do iv = 1, mesh%n_vert
        dx = abs(mesh%vert(iv)%coord(1) - xmin)
        if (dx > 1.0e-12_DOUBLE .and. dx < hx) hx = dx
      end do
      n_moved = 0
      do iv = 1, mesh%n_vert
        if (abs(mesh%vert(iv)%coord(2) - ymid) < 0.25_DOUBLE*hx .and. &
            abs(mesh%vert(iv)%coord(3) - zmid) < 0.25_DOUBLE*hx) then
          k = nint((mesh%vert(iv)%coord(1) - xmin)/hx)
          mesh%vert(iv)%coord(2) = mesh%vert(iv)%coord(2) &
            + odd_even_amp*hx*real(1 - 2*modulo(k, 2), DOUBLE)
          n_moved = n_moved + 1
        end if
      end do
      if (me == 0) print*, "ODD_EVEN: hx =", hx, " amplitude =", odd_even_amp*hx, &
        " sommets deplaces =", n_moved
    end block
  end if

  ! Saltzman (init=5): skewed mesh by mapping X_sk = X + (0.1-Y)*sin(pi*X), Y_sk = Y
  ! (Maire 2007 / Vilar et al.). sin(pi*X)=0 at X=0,1 so left/right boundaries unaffected.
  if (init == 5) then
    block
      real(kind=DOUBLE), parameter :: pi_cst = acos(-1.0_DOUBLE)
      real(kind=DOUBLE) :: xv, yv
      do i = 1, mesh%n_vert
        xv = mesh%vert(i)%coord(1)
        yv = mesh%vert(i)%coord(2)
        mesh%vert(i)%coord(1) = xv + (0.1_DOUBLE - yv) * sin(pi_cst*xv)
      end do
    end block
  end if

  call compute_geometry_mesh(mesh, .true., boundary_2d)

  ! Extrusion height of a boundary_2d mesh (single z-layer): used by the
  ! classic_iso_vol weight to turn a prism sub-element volume back into an
  ! in-plane area scale.  The nodal solver keeps vp_z = 0 in b2d, so this is
  ! constant in time.
  h_extrude = 1.0_DOUBLE
  if (boundary_2d .and. mesh%n_vert > 0) then
    h_extrude = maxval(mesh%vert(1:mesh%n_vert)%coord(3)) &
              - minval(mesh%vert(1:mesh%n_vert)%coord(3))
    if (h_extrude <= 0.0_DOUBLE) h_extrude = 1.0_DOUBLE
  end if

  allocate(vp(3, mesh%n_vert))
  if (diag_alpha) allocate(alpha_sal(mesh%n_elems))
  if (pi_fbw) then
    allocate(fbw_cviol(mesh%n_elems)); fbw_cviol = 0
  end if
  allocate(pp(mesh%n_vert))
  allocate(sol(5, mesh%n_elems))
  allocate(new_sol(5, mesh%n_elems))
  allocate(rhs(5, mesh%n_elems))
  allocate(gamma_arr(mesh%n_elems))
  allocate(vp_is_imposed(mesh%n_vert))

  ! gamma_uniform > 0 overrides the per-init default. Needed for the *classic*
  ! Noh problem, which is posed with gamma = 5/3 (post-shock rho = 16 in 2D
  ! cylindrical, shock at r = t/3) whereas init == 3 defaults to 1.4 here.
  if (gamma_uniform > 0.0_DOUBLE) then
    gamma_arr = gamma_uniform
  else if (init == 5) then
    gamma_arr = 5.0_DOUBLE / 3.0_DOUBLE
  else
    gamma_arr = 1.4_DOUBLE
  end if
  vp = 0.0_DOUBLE
  call init_sol(mesh, sol, sol_uniform, init, me, num_procs, gamma_arr, boundary_2d)
  call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol)

  ! Identify piston nodes from bc_type='piston'
  vp_is_imposed = .false.
  piston_bc_idx = 0
  piston_vel = 0.0_DOUBLE
  do i = 1, n_bc
    if (trim(bc_type(i)) == 'piston') then
      piston_bc_idx = i
      piston_vel = bc_val(2:4, i)
      exit
    end if
  end do
  if (piston_bc_idx > 0) then
    do i = 1, mesh%n_vert
      do j = 1, mesh%vert(i)%n_sub_faces_neigh
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        if (mesh%face(id_face)%right_neigh == -piston_bc_idx) then
          vp_is_imposed(i) = .true.
          exit
        end if
      end do
    end do
  end if

  if (md_watch) then
    allocate(w_prod(mesh%n_elems), source=0.0_DOUBLE)
    allocate(w_dpi_min(mesh%n_elems), source=1.0e300_DOUBLE)
    allocate(w_dpi_max(mesh%n_elems), source=-1.0e300_DOUBLE)
  end if
  if (diag_divratio .or. diag_firstprint) then
    diag_divratio = .true.
    allocate(dr_vdot(mesh%n_elems), source=0.0_DOUBLE)
  end if
  allocate(grad_v(3, 3, mesh%n_elems))
  allocate(grad_p(3, mesh%n_elems))
  allocate(div_v(mesh%n_elems))
  allocate(alpha_p_arr(mesh%n_vert))
  allocate(p_limited(mesh%n_elems))
  allocate(h_p_arr(mesh%n_vert))
  grad_v = 0.0_DOUBLE
  grad_p = 0.0_DOUBLE
  div_v  = 0.0_DOUBLE
  alpha_p_arr = 1.0_DOUBLE
  p_limited   = 0
  h_p_arr     = 0.0_DOUBLE
  if (scheme == "sidil") then
    do i = 1, mesh%n_vert
      h_p_arr(i) = compute_length(mesh, i, method_length, boundary_2d, h_extrude)
    end do
  end if

  i_sol_vtu = 0
  write(fln, *) i_sol_vtu
  write(fln, *) "output_"//trim(adjustl(fln))
  call write_sol_lag(mesh, fln, sol, vp, pp, gamma_arr, grad_v, grad_p, div_v, alpha_p_arr, p_limited, h_p_arr)
  call write_sol_dat_lag(mesh, fln, sol, gamma_arr)
  i_sol_vtu = i_sol_vtu + 1

  allocate(mass(mesh%n_elems))
  mass = 0.0_DOUBLE
  do i=1, mesh%n_elems
    mass(i) = mesh%elem(i)%volume/sol(1, i)
  end do

  ! Pression imposee par condition de bord, pour la borne acoustique de compute_dt.
  bc_p = -1.0_DOUBLE
  do i = 1, n_bc
    if (trim(bc_type(i)) == 'pressure') bc_p(i) = bc_val(1, i)
  end do
  t = 0.0_DOUBLE
  ! Lu par le demi-pas GRP avant d'etre affecte dans l'ancien ordre ; inutilise
  ! quand dt_at_step_start est vrai.
  dt = 0.0_DOUBLE
  rk2_time = (trim(time_scheme) == 'rk2')
  nstage = 1
  if (rk2_time) then
    ! reconstruction spatiale seule : c'est le demi-pas de MAILLAGE qui porte
    ! desormais l'ordre 2 en temps.
    rec_dt_fac = 0.0_DOUBLE
    nstage = 2
    allocate(coord_n(3, mesh%n_vert), sol_n(5, mesh%n_elems))
    if (me == 0) print*, "RK2 : reconstruction spatiale, maillage au demi-pas"
  end if
  iter = 1
  do while ( t < t_max )

    do_trace = (trace_every <= 1) .or. (mod(iter, trace_every) == 0)

    ! RK2 : on repartira de X^n et U^n pour le correcteur.
    if (rk2_time) then
      do i = 1, mesh%n_vert
        coord_n(:, i) = mesh%vert(i)%coord
      end do
      sol_n = sol
    end if

    ! dt en tete de pas : la coupure en aval est idempotente.
    if (dt_at_step_start) then
      ! Phase 1 : les DEUX bornes, la borne de volume servant d'ESTIMATION avec le
      ! v_p precedent -- elle ne decide pas de l'admissibilite, la phase 2 s'en
      ! charge. Avec la seule borne acoustique, le predicteur d'ordre 2 recevait un
      ! dt jusqu'a 9.3x trop grand (mesure sur sedov_tri_rnd).
      call compute_dt(mesh, sol, dt, cfl, vp, me, num_procs, gamma_arr, bc_p=bc_p)
      if (t + dt > t_max) dt = t_max - t
    end if

    ! Set imposed velocities (piston BC) before the RHS computation
    if (piston_bc_idx > 0) then
      do i = 1, mesh%n_vert
        if (vp_is_imposed(i)) vp(:, i) = piston_vel
      end do
    end if

    ! W12 a besoin de grad_v A L'ORDRE 1 : il est initialise a zero ligne 269 et
    ! n'etait rempli que sous second_order. Sans cet appel la viscosite
    ! tensorielle serait identiquement nulle et le schema rendrait du V0 sans
    ! aucun signal. second_order reste .false., donc la RECONSTRUCTION demeure
    ! d'ordre 1 : seul le tenseur visqueux consomme le gradient.
    ! Iteration du predicteur. Le demi-pas GRP tourne avec dt_1 ; si la phase 2
    ! reduit dt on refait le predicteur avec le dt reduit, puis on RE-EVALUE la
    ! borne : sans cette re-evaluation finale la garantie d'admissibilite saute.
    do istage = 1, nstage

    do it_dt = 1, max(1, dt_n_iter)

    if (second_order .or. pos_base_mode == 'WE') then
      call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol)
      call compute_gradients(mesh, sol, gamma_arr, dt, grad_v, grad_p, div_v, p_limited, grad_weno, &
        grad_gg_corrected, grad_eps_weno, bound_skip=grad_bound_skip, &
        bound_interior=grad_bound_interior, b2d=boundary_2d)
      ! Exchange ghost cell gradients so second-order RHS can reconstruct across MPI boundaries
      call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 3, grad_p)
      call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 9, grad_v)
      call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 1, div_v)
    end if
    ! Garde explicite. Elle ne peut PAS tester "grad_v nul" tel quel : a t=0 sur
    ! Sedov le fluide est au repos partout, donc grad_v y est legitimement nul.
    ! Elle teste donc l'implication qui compte : un champ de vitesse non nul avec
    ! un gradient exactement nul signifie que l'appel ci-dessus n'a pas eu lieu.
    if (pos_base_mode == 'WE') then
      if (maxval(abs(sol(2:4, :))) > 0.0_DOUBLE .and. maxval(abs(grad_v)) == 0.0_DOUBLE) then
        if (me == 0) print*, "W12 : grad_v EXACTEMENT nul avec une vitesse non nulle -- ", &
          "compute_gradients n'a pas tourne, le schema serait du V0 deguise. Arret."
        error stop
      end if
    end if

    ! --- 'vp_kappa' : kappa EXPLICITE, calcule sur la vitesse nodale du pas
    ! precedent (vp porte encore ses valeurs d'entree ici) et sur la geometrie
    ! courante. Il traverse l'echange MPI generique par un tableau de maille a
    ! largeur fixe : une maille FANTOME a besoin de tous ses sommets pour sa
    ! divergence de sous-maille, et certains sont hors du halo, donc la calculer
    ! localement pour les fantomes donnerait deux Pi differents de part et
    ! d'autre d'une frontiere de rang.
    ! 'vp_mdiv' passe par le MEME chemin : c'est `compute_kappa` qui remplit
    ! `mdv_c` et `kappa_scatter` qui le redistribue en `mdv_se`. Sans cet appel
    ! `mdv_se` n'est JAMAIS alloue, le garde `allocated(mdv_se)` de la fermeture
    ! est toujours faux, et la correction de divergence est du CODE MORT : le
    ! schema se comporte exactement comme 'vp_fbw' et tous les compteurs mdv_*
    ! sortent a zero -- ce qui se lit comme « le plafond ne mord jamais » alors
    ! que ca veut dire « rien ne s'est passe ». Sans effet de bord sur 'vp_kappa' :
    ! `kap_se` est bien calcule sous `pi_mdiv` mais `kl(j)` n'est lu que sous
    ! `pi_kappa`, donc le poids kappa y reste a 1.
    if (pi_kappa .or. pi_mdiv) then
      call compute_kappa(mesh, sol, vp, gamma_arr, boundary_2d, dt)
      call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, KAP_NV, kap_c)
      ! `mdv_c` a besoin du MEME echange que `kap_c`, et pour la meme raison :
      ! une maille FANTOME a besoin de tous ses sommets pour la divergence de
      ! son sous-element, dont certains sont hors du halo. Sans cet echange la
      ! dilatation des fantomes serait fausse et les deux rangs voisins
      ! batiraient deux Pi differents sur la meme sous-face -- donc une perte de
      ! conservation a la frontiere de rang, invisible en sequentiel.
      if (pi_mdiv .and. allocated(mdv_c)) &
        call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, KAP_NV, mdv_c)
      call kappa_scatter(mesh)
    end if

    ! --- 'vp_div' : b_pc, h_c et reste_pc, puis l'echange MPI generique. Le
    ! reste est pris sur la vitesse nodale du pas precedent, ce qui decouple le
    ! systeme nodal sans rien perdre, seul b_pc.v_p ayant besoin d'etre implicite.
    if (pi_divdirect) then
      call divdirect_prepare(mesh, vp, boundary_2d, sol, gamma_arr)
      call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 3*KAP_NV, divd_bc)
      call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, KAP_NV, divd_restc)
      call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, KAP_NV, divd_isumc)
      call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 1, divd_h)
      call divdirect_scatter(mesh)
    end if

    ! Controles en dur AVANT tout calcul, une seule fois.
    if (pi_pif .and. .not. pif_check_done) then
      pif_check_done = .true.
      if (me == 0) call pif_selfcheck(mesh)
    end if

    if (pvp_pos_sv) call sv_step_reset()
    if (pi_fbd .or. pi_fbw) call fbd_step_reset()
    if (diag_alpha .and. allocated(alpha_sal)) alpha_sal = 0.0_DOUBLE
    if (diag_force) then
      ! La 4e cible est choisie automatiquement : la maille de plus grande energie
      ! interne a l'etat INITIAL, c'est-a-dire la maille du depot. Les trois autres
      ! viennent du namelist. Choisie une seule fois, au premier pas, parce qu'apres
      ! le depot s'etale et l'argmax se deplace.
      if (iter == 1) then
        block
          integer(kind=ENTIER) :: kf, kfmax, kfv
          real(kind=DOUBLE) :: ef, efmax
          efmax = -1.0e300_DOUBLE; kfmax = -1
          do kf = 1, mesh%n_elems
            if (mesh%elem(kf)%is_ghost) cycle
            ef = sol(5, kf) - 0.5_DOUBLE*dot_product(sol(2:4, kf), sol(2:4, kf))
            if (ef > efmax) then; efmax = ef; kfmax = kf; end if
          end do
          fdiag_cells(1) = fdiag_c1; fdiag_cells(2) = fdiag_c2
          fdiag_cells(3) = fdiag_c3; fdiag_cells(4) = kfmax
          ! Si la 3e cible n'est pas donnee, on prend automatiquement une maille
          ! ADJACENTE AU PISTON : celle qui porte un sommet a vitesse imposee.
          ! Sur Saltzman c'est la seule qui compte, et l'argmax d'energie interne
          ! ne la designe pas puisque l'etat initial y est uniforme.
          if (fdiag_c3 < 0) then
            do kf = 1, mesh%n_elems
              if (mesh%elem(kf)%is_ghost) cycle
              do kfv = 1, mesh%elem(kf)%n_vert
                if (vp_is_imposed(mesh%elem(kf)%vert(kfv))) then
                  fdiag_cells(3) = kf; exit
                end if
              end do
              if (fdiag_cells(3) > 0) exit
            end do
            if (me == 0) print*, "FDIAG maille adjacente au PISTON =", fdiag_cells(3)
          end if
          if (me == 0) print*, "FDIAG maille du depot (argmax eps initial) =", kfmax, &
                               " eps", efmax
        end block
      end if
      call fdiag_reset()
    end if

    if( scheme == "vp_pif" ) then
      ! Test d'attribution : v_p vient du solveur vp, dont le rhs est JETE --
      ! seules les vitesses nodales nous interessent. Les flux restent ceux de
      ! la construction par sous-face.
      if (pif_vp_classic) then
        call compute_rhs_lagrange(mesh, sol, vp, dt, rhs, n_bc, bc_type, bc_val, &
          boundary_2d, mass, gamma_arr, vp_is_imposed, .false., grad_v, grad_p, div_v)
        rhs = 0.0_DOUBLE
      end if
      call compute_rhs_lagrange_pif(mesh, sol, vp, dt, rhs, n_bc, bc_type, &
        bc_val, boundary_2d, mass, gamma_arr, vp_is_imposed)
    else if( scheme == "vp_div" ) then
      call compute_rhs_lagrange_divdirect(mesh, sol, vp, dt, rhs, n_bc, bc_type, &
        bc_val, boundary_2d, mass, gamma_arr, vp_is_imposed)
    else if( scheme == "mean_div" ) then
      call compute_rhs_lagrange_meandiv(mesh, sol, vp, dt, rhs, n_bc, bc_type, bc_val, &
        boundary_2d, mass, gamma_arr, vp_is_imposed, grad_v, grad_p, div_v)
    else if( scheme == "classic" .or. scheme == "classic_rhoa" .or. scheme == "vp_hat" &
        .or. scheme == "classic_conv" &
        .or. scheme == "vp_hat_md" .or. scheme == "vp_facebar" &
        .or. scheme == "vp_mp" .or. scheme == "vp_kappa" &
        .or. scheme == "vp_fbd" .or. scheme == "vp_fbw" &
        .or. scheme == "vp_mdiv" ) then
      call compute_rhs_lagrange(mesh, sol, vp, &
        dt, rhs, n_bc, bc_type, bc_val, boundary_2d, mass, gamma_arr, vp_is_imposed, &
        second_order, grad_v, grad_p, div_v, alpha_p_arr)
    else if( scheme == "classic_iso" ) then
      call compute_rhs_lagrange_classic_iso(mesh, sol, vp, &
        dt, rhs, n_bc, bc_type, bc_val, boundary_2d, mass, gamma_arr, vp_is_imposed, &
        1, h_extrude)
    else if( scheme == "classic_iso_vol" ) then
      call compute_rhs_lagrange_classic_iso(mesh, sol, vp, &
        dt, rhs, n_bc, bc_type, bc_val, boundary_2d, mass, gamma_arr, vp_is_imposed, &
        iso_weight_mode, h_extrude)
    else  if( scheme == "sidil" ) then
      ! h_extrude (auto-detected from the mesh's own z-extent above), not
      ! the user-set b2d_h namelist value: relying on a manually-set
      ! parameter to match the mesh's actual extrusion thickness is a
      ! standing footgun (easy to forget/mistype per test case) now that
      ! compute_length's b2d branch genuinely depends on it being exact.
      call compute_rhs_lagrange_sidil(mesh, sol, vp, &
        dt, rhs, n_bc, bc_type, bc_val, boundary_2d, mass, method_length, h_extrude, &
        gamma_arr, vp_is_imposed, h_p_arr)
    else if( scheme == "vp_pp_pvp" .or. scheme == "vp_pp_pvp_e" ) then
      call compute_rhs_lagrange_vp_pp_pvp(mesh, sol, vp, &
        dt, rhs, n_bc, bc_type, bc_val, boundary_2d, mass, gamma_arr, vp_is_imposed)
    else if( scheme == "vp_ep" .or. scheme == "vp_ep_wip" &
        .or. scheme == "vp_pp_ep" .or. scheme == "vp_pp_ep_i" &
        .or. scheme == "vp_pp_epsp" .or. scheme == "vp_pp_eps2" &
        .or. scheme == "vp_pp_eps2m" .or. scheme == "vp_pp_epv" &
        .or. scheme == "vp_pp_vpps" &
        .or. scheme == "vp_pp_ppvn" &
        .or. scheme == "vp_pp_ep_wdiv" &
        .or. scheme == "vp_pp_ppvn_wdiv" &
        .or. scheme == "vp_pp_ep_veps" &
        .or. scheme == "vp_pp_pivp" &
        .or. scheme == "vp_pp_ppus" &
        .or. scheme == "vp_ep_i" ) then
      ! vitesse nodale + energie nodale, puis solveur 1D par sous-face.
      call compute_rhs_lagrange_vp_ep(mesh, sol, vp, &
        dt, rhs, n_bc, bc_type, bc_val, boundary_2d, mass, gamma_arr, vp_is_imposed, &
        second_order, grad_v, grad_p, div_v)
    else if( scheme == "vp_pp" .or. scheme == "nodal_vp" .or. scheme(1:min(9,len_trim(scheme))) == "vp_pp_pos" &
        .or. scheme(1:min(9,len_trim(scheme))) == "vp_pp_pur" .or. scheme == "vp_pp_psvp" ) then
      ! vitesse nodale ET pression nodale, toutes deux issues de systemes nodaux :
      ! pas de h_p, donc aucune dependance a method_length / b2d_h.
      ! "nodal_vp" est l'ancien nom du schema, conserve en alias : beaucoup de
      ! input_data.f de runs deja effectues le portent encore.
      call compute_rhs_lagrange_vp_pp(mesh, sol, vp, &
        dt, rhs, n_bc, bc_type, bc_val, boundary_2d, mass, gamma_arr, vp_is_imposed, &
        second_order, grad_v, grad_p, div_v, .false.)
    else
      print*, "No scheme !"
      error stop
    end if

    ! Phase 2 : borne de VOLUME sur le v_p qu'on vient de calculer. A l'etage 2
    ! du RK2 le maillage a deja bouge de dt/2 : reduire dt ici rendrait le demi-pas
    ! deja effectue incoherent avec le pas complet.
    if (dt_at_step_start .and. istage == 1) then
      dt_post = dt
      ! bc_p passe aussi ici : sans effet sur dt (do_acou faux) mais l'appel
      ! reinitialise dtb_acou, donc le diagnostic imprime serait sinon celui du
      ! lambda SANS la pression imposee, pas celui du schema.
      call compute_dt(mesh, sol, dt_post, cfl, vp, me, num_procs, gamma_arr, &
        only='V', bc_p=bc_p)
      if (dt_post < dt) then
        if (it_dt == 1) then
          n_dtviol = n_dtviol + 1
          dtviol_max = max(dtviol_max, dt/max(dt_post, 1.0e-300_DOUBLE))
        else
          n_dtiter2 = n_dtiter2 + 1
          dtiter2_max = max(dtiter2_max, dt/max(dt_post, 1.0e-300_DOUBLE))
        end if
        dt = dt_post
      else
        exit
      end if
    end if
    ! A l'ordre 1 le rhs ne depend pas de dt, et en mode historique il n'y a pas
    ! de phase 2 : re-predire ne changerait rien et doublerait les compteurs.
    if (.not. second_order .or. .not. dt_at_step_start) exit
    end do

    ! Source Taylor-Green : sur le maillage PRE-deplacement de l'etage courant,
    ! donc elle doit alimenter aussi le rhs du predicteur.
    if (init == 10) call add_taylor_green_source(mesh, mass, rhs)

    ! PREDICTEUR RK2 : X^{n+1/2} = X^n + dt/2 v_p^n, U^{n+1/2} = U^n + dt/2 rhs^n.
    ! L'etage 2 recalcule gradients, solveur nodal et flux sur CETTE geometrie.
    if (rk2_time .and. istage == 1) then
      call move_mesh(mesh, vp, 0.5_DOUBLE*dt)
      call mpi_memory_exchange_vert(mesh, mpi_send_recv)
      call compute_geometry_mesh(mesh, .true., boundary_2d, do_check=.false.)
      sol = sol_n + 0.5_DOUBLE*dt*rhs
      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) sol(1, i) = mesh%elem(i)%volume / mass(i)
      end do
      call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol)
    end if
    end do

    if (diag_force .and. me == 0) call fdiag_print()
    if (pi_fbw .and. do_trace .and. me == 0) print*, "FBWP pas", iter, " d2<=0 ff/cf/cc", &
      fbw_st_ff, fbw_st_cf, fbw_st_cc, " delta2min", fbw_st_d2min, &
      " (lam_w/rho a) max du pas", fbw_st_rmax
    if (pi_fbd .and. do_trace .and. me == 0) print*, "FBDP pas", iter, " |D_p|max", fbd_st_dpmax, &
      " delta2min", fbd_st_d2min, " d2<=0 ff/cf/cc", fbd_st_ff, fbd_st_cf, fbd_st_cc, &
      " lambda invalides", fbd_st_zbad, " denominateurs nuls", fbd_st_void
    if (pvp_pos_sv .and. do_trace .and. me == 0) then
      if (sv_step_n > 0) then
        print*, "POSSVP pas", iter, " intervalles VIDES", sv_st_infeas, &
                " delta2min", sv_st_d2min, " corrigees", sv_step_n, " rayon moyen", &
                sv_step_r/real(sv_step_n, DOUBLE), " de", sv_step_rmin, " a", sv_step_rmax
      else
        print*, "POSSVP pas", iter, " intervalles VIDES", sv_st_infeas, &
                " delta2min", sv_st_d2min, " corrigees", 0
      end if
    end if

    if (.not. dt_at_step_start) &
      call compute_dt(mesh, sol, dt, cfl, vp, me, num_procs, gamma_arr, bc_p=bc_p)
    ! dt = NaN rend `t + dt < t_max` FAUX : la boucle se terminerait normalement avec
    ! EXIT_CODE 0 sur un calcul mort. Trois sessions s'y sont fait prendre, d'ou l'arret
    ! immediat. dt_nan_detected etant reduit sur tous les rangs, l'arret est COLLECTIF.
    ! Vitesse nodale d ORDRE 1 NaN : aucun ordre inferieur ou se replier, donc c est une
    ! panne. Reduit pour que tous les rangs s arretent ensemble.
    inan_vp = vp_o1_nan_n
    if (num_procs > 1) call MPI_ALLREDUCE(MPI_IN_PLACE, inan_vp, 1, MPI_INTEGER, MPI_SUM, &
      MPI_COMM_WORLD, mpi_ierr)
    if (inan_vp > 0) then
      if (me == 0) print*, "PANNE : vitesse nodale d'ordre 1 NaN sur", inan_vp, &
        " noeud(s) au pas", iter, " -- Mp singuliere, 1er noeud local du rang 0 :", vp_o1_nan_node
      error stop
    end if
    if (dt /= dt .or. dt_nan_detected) then
      if (me == 0) print*, "PANNE : dt = NaN au pas", iter, " t =", t, " -- arret."
      error stop
    end if
    ! Instrumentation : les bornes rapportees par compute_dt et le dt reellement
    ! retenu se contredisaient au pas 4 de vp_fbd (bornes 3.03e-05 et 6.02e-05, dt
    ! retenu 0.99993508 = t_max - t). On imprime donc dt AVANT et APRES le clip pour
    ! savoir lequel des deux diagnostics mentait, au lieu de le deduire.
    if (diag_dtb .and. me == 0) print*, "TRACE   dt sorti de compute_dt =", dt, &
      " t =", t, " clip actif =", (t + dt > t_max)
    if( t + dt > t_max ) dt = t_max - t
    if (diag_alpha .and. me == 0) then
      block
        integer(kind=ENTIER) :: ka
        real(kind=DOUBLE), dimension(:), allocatable :: atau, apre
        allocate(atau(mesh%n_elems), apre(mesh%n_elems))
        do ka = 1, mesh%n_elems
          atau(ka) = sol(1, ka)
          apre(ka) = pressure(sol(:, ka), gamma_arr(ka))
        end do
        call alpha_report(mesh, mass, dt, iter, atau, apre, gamma_arr)
        deallocate(atau, apre)
      end block
    end if

    ! --- trace par pas INDEPENDANTE DU SCHEMA : elle sert a comparer vp et
    ! mean_div sur le meme maillage, donc elle doit etre identique pour les deux.
    ! Elle sort apres compute_dt pour que le dt imprime soit celui qui sera
    ! REELLEMENT utilise, et non celui du pas precedent.
    ! Etat initial fige au PREMIER passage, apres l'initialisation complete :
    ! c'est lui qui definit le plancher isentropique de chaque maille.
    if (diag_iso .and. .not. allocated(iso_eps0)) then
      allocate(iso_eps0(mesh%n_elems), iso_tau0(mesh%n_elems))
      do iso_k = 1, mesh%n_elems
        iso_tau0(iso_k) = sol(1, iso_k)
        iso_eps0(iso_k) = sol(5, iso_k) - 0.5_DOUBLE*dot_product(sol(2:4, iso_k), sol(2:4, iso_k))
      end do
    end if
    if (diag_iso .and. me == 0) then
      block
        integer(kind=ENTIER) :: kc, kmin, krmin
        real(kind=DOUBLE) :: ep, pl, rr, epmin, rrmin
        epmin = 1.0e300_DOUBLE; rrmin = 1.0e300_DOUBLE; kmin = -1; krmin = -1
        do kc = 1, mesh%n_elems
          if (mesh%elem(kc)%is_ghost) cycle
          if (iso_eps0(kc) <= 0.0_DOUBLE .or. iso_tau0(kc) <= 0.0_DOUBLE) cycle
          ep = sol(5, kc) - 0.5_DOUBLE*dot_product(sol(2:4, kc), sol(2:4, kc))
          if (sol(1, kc) <= 0.0_DOUBLE .or. ep /= ep) cycle
          pl = iso_eps0(kc)*(iso_tau0(kc)/sol(1, kc))**(gamma_arr(kc) - 1.0_DOUBLE)
          if (pl <= 0.0_DOUBLE) cycle
          rr = ep/pl
          if (ep < epmin) then; epmin = ep; kmin = kc; end if
          if (rr < rrmin) then; rrmin = rr; krmin = kc; end if
        end do
        if (kmin > 0) then
          ep = sol(5, kmin) - 0.5_DOUBLE*dot_product(sol(2:4, kmin), sol(2:4, kmin))
          pl = iso_eps0(kmin)*(iso_tau0(kmin)/sol(1, kmin))**(gamma_arr(kmin) - 1.0_DOUBLE)
          print*, "ISO maille eps_min", kmin, " eps", ep, " V/V0", sol(1, kmin)/iso_tau0(kmin), &
                  " eps/plancher", ep/pl
        end if
        if (krmin > 0) print*, "ISO rapport MINIMAL eps/plancher =", rrmin, " maille", krmin, &
                  " V/V0", sol(1, krmin)/iso_tau0(krmin)
        ! Maille au V/V0 MINIMAL, avec son identite : si c'est toujours la meme
        ! qui s'ecrase pendant que les autres vont bien, le mecanisme est local
        ! et on sait ou regarder ; si elle change, c'est diffus.
        block
          integer(kind=ENTIER) :: kv, kvmin, jv, nb, n50, n90, n99, nliv
          integer(kind=ENTIER) :: js, ks, idse, idsf, idf, idv
          real(kind=DOUBLE) :: vv, vvmin, cg(3), nrm(3), un, aa
          real(kind=DOUBLE) :: ssig, sabs, el, elmin, elmax, dd
          logical :: tbord
          vvmin = 1.0e300_DOUBLE; kvmin = -1
          n50 = 0; n90 = 0; n99 = 0; nliv = 0
          do kv = 1, mesh%n_elems
            if (mesh%elem(kv)%is_ghost) cycle
            if (iso_tau0(kv) <= 0.0_DOUBLE .or. sol(1, kv) <= 0.0_DOUBLE) cycle
            vv = sol(1, kv)/iso_tau0(kv)
            nliv = nliv + 1
            if (vv < 0.5_DOUBLE) n50 = n50 + 1
            if (vv < 0.9_DOUBLE) n90 = n90 + 1
            if (vv < 0.99_DOUBLE) n99 = n99 + 1
            if (vv < vvmin) then; vvmin = vv; kvmin = kv; end if
          end do
          if (kvmin > 0) then
            cg = 0.0_DOUBLE; tbord = .false.
            nb = mesh%elem(kvmin)%n_vert
            do jv = 1, nb
              cg = cg + mesh%vert(mesh%elem(kvmin)%vert(jv))%coord
              if (mesh%vert(mesh%elem(kvmin)%vert(jv))%is_bound) tbord = .true.
            end do
            cg = cg/real(nb, DOUBLE)
            ! Forme de la maille : distances entre sommets DANS LE PLAN (x,y)
            ! seulement. Le maillage est une extrusion a une maille en z, donc
            ! inclure la composante z fait ressortir la diagonale de l'extrusion
            ! (longueur ~1) et fabrique un faux rapport d'aiguille de ~80. Deux
            ! sommets superposes en (x,y) sont les deux copies d'un meme noeud 2D
            ! et sont ecartes du decompte.
            elmin = 1.0e300_DOUBLE; elmax = 0.0_DOUBLE
            do jv = 1, nb
              do js = jv + 1, nb
                el = norm2(mesh%vert(mesh%elem(kvmin)%vert(jv))%coord(1:2) &
                         - mesh%vert(mesh%elem(kvmin)%vert(js))%coord(1:2))
                if (el > 1.0e-14_DOUBLE) then
                  elmin = min(elmin, el); elmax = max(elmax, el)
                end if
              end do
            end do
            ! r = (somme signee A u_p.n) / (somme des |A u_p.n|) : la COHERENCE de la
            ! compression. La borne de volume vaut CV*V/somme|..| donc la perte de
            ! volume par pas est CV*r, independante de dt. r=1 = tous les noeuds
            ! rentrent de concert, r=0 = cisaillement pur sans changement de volume.
            ssig = 0.0_DOUBLE; sabs = 0.0_DOUBLE
            do js = 1, mesh%elem(kvmin)%n_sub_elems
              idse = mesh%elem(kvmin)%sub_elem(js)
              do ks = 1, mesh%sub_elem(idse)%n_sub_faces
                idsf = mesh%sub_elem(idse)%sub_face(ks)
                idf = mesh%sub_face(idsf)%mesh_face
                idv = mesh%sub_face(idsf)%mesh_vert
                if (mesh%face(idf)%left_neigh == kvmin) then
                  nrm = mesh%face(idf)%norm
                else
                  nrm = -mesh%face(idf)%norm
                end if
                aa = mesh%sub_face(idsf)%area
                un = dot_product(vp(:, idv), nrm)
                ssig = ssig + aa*un; sabs = sabs + aa*abs(un)
              end do
            end do
            print*, "ISOV V/V0 minimal =", vvmin, " maille", kvmin, &
                    " cdg", cg(1), cg(2), " |cdg|", norm2(cg(1:2)), &
                    " nvert", nb, " bord", tbord, " rho", 1.0_DOUBLE/sol(1, kvmin)
            print*, "ISOVD mailles vivantes", nliv, " sous 0.99 :", n99, &
                    " sous 0.90 :", n90, " sous 0.50 :", n50
            print*, "ISOVF arete min", elmin, " max", elmax, " aiguille", elmax/elmin, &
                    " r_coherence", ssig/max(sabs, 1.0e-30_DOUBLE), &
                    " vc", sol(2, kvmin), sol(3, kvmin)
            ! Les vitesses nodales de ses sommets, projetees sur la direction
            ! radiale sortante : negatif = ce noeud rentre vers le depot.
            do jv = 1, nb
              idv = mesh%elem(kvmin)%vert(jv)
              dd = norm2(mesh%vert(idv)%coord(1:2))
              print*, "ISOVN som", jv, " id", idv, " x", mesh%vert(idv)%coord(1), &
                      " y", mesh%vert(idv)%coord(2), " vp", vp(1, idv), vp(2, idv), &
                      " v_radial", (vp(1, idv)*mesh%vert(idv)%coord(1) &
                                  + vp(2, idv)*mesh%vert(idv)%coord(2)) &
                                  /max(dd, 1.0e-30_DOUBLE), " bord", mesh%vert(idv)%is_bound
            end do
          end if
        end block
      end block
    end if

    if (md_trace .and. do_trace .and. me == 0) then
      block
        integer(kind=ENTIER) :: kt, nn, nlist
        real(kind=DOUBLE) :: taumax, epsmin, ee
        real(kind=DOUBLE) :: epsmin_v, epsmin_s
        integer(kind=ENTIER) :: nnan, nvc, nvpnan
        character(len=400) :: lst
        taumax = -1.0e300_DOUBLE; epsmin = 1.0e300_DOUBLE
        nnan = 0; nlist = 0; lst = ""
        epsmin_v = 1.0e300_DOUBLE; epsmin_s = 1.0e300_DOUBLE; nvc = 0
        do kt = 1, mesh%n_elems
          if (mesh%elem(kt)%is_ghost) cycle
          if (any(sol(:, kt) /= sol(:, kt))) then
            nnan = nnan + 1
            if (nlist < 12) then
              nlist = nlist + 1
              write(lst, '(A,1X,I0)') trim(lst), kt
            end if
            cycle
          end if
          taumax = max(taumax, sol(1, kt))
          ee = sol(5, kt) - 0.5_DOUBLE*dot_product(sol(2:4, kt), sol(2:4, kt))
          epsmin = min(epsmin, ee)
          ! Les deux POPULATIONS separees : mailles portant au moins une sous-face
          ! a delta2 <= 0, et les autres. Sans cette separation, "eps_min > 0 alors
          ! que des milliers de sous-faces violent" ne dit pas si les mailles
          ! concernees sont justement celles qui vont bien -- l'hypothese de
          ! repartition est intestable sans elle.
          if (allocated(fbw_cviol)) then
            if (fbw_cviol(kt) > 0) then
              nvc = nvc + 1
              epsmin_v = min(epsmin_v, ee)
            else
              epsmin_s = min(epsmin_s, ee)
            end if
          end if
        end do
        ! max|v_p| contre max|v_c| : si v_p est d'un ordre au-dessus, c'est que
        ! le plancher de kappa est trop bas (le terme de pression de R_p n'est
        ! pas pondere, donc v_p gagne une composante en h|grad p|/kappa). C'est
        ! ce symptome qu'il faut regarder, pas l'apparition de NaN.
        ! `maxval(abs(vp))` est AVEUGLE aux NaN : gfortran rend le maximum des
        ! valeurs non-NaN, donc il a affiche 1.0000000000000000 exactement sur
        ! saltzman au pas 223 alors que 27 mailles portaient un `somme A|v_p.n|`
        ! NaN. Et la colonne `NaN` ci-dessous ne teste que `sol` : le champ NODAL
        ! n'a jamais ete sous surveillance. On le compte explicitement.
        nvpnan = count(vp /= vp)
        print*, "TRACE   max|v_p| =", maxval(abs(vp)), "  max|v_c| =", &
                maxval(abs(sol(2:4, :))), "  NaN dans v_p :", nvpnan
        if (diag_dtb) print*, "TRACE   borne acoustique", dtb_acou, " maille", dtb_cell_acou, &
                " | borne VOLUME", dtb_vol, " maille", dtb_cell_vol, &
                " V", dtb_vol_V, " somme A|u_p.n|", dtb_vol_sumav
        if (diag_dtb) print*, "TRACE   acou decomp : m", dtb_acou_m, " somme A lam", &
                dtb_acou_sal, " lam", dtb_acou_lam, " somme A", dtb_acou_area
        print*, "TRACE pas", iter, " t", t, " dt", dt, " rho_min", 1.0_DOUBLE/taumax, &
                " eps_min", epsmin, " NaN", nnan, " NaN_vp", nvpnan
        if (diag_area .and. me == 0) then
          block
          integer(kind=ENTIER) :: kf
          ! A(0) est pris au PREMIER pas, sur la geometrie initiale.
          if (.not. allocated(area_a0)) then
            allocate(area_a0(mesh%n_elems))
            call area_of_cells(mesh, area_a0)
            allocate(forme_a0(mesh%n_elems))
            do kf = 1, mesh%n_elems
              if (abs(mesh%elem(kf)%volume) > 0.0_DOUBLE .and. area_a0(kf) > 0.0_DOUBLE) then
                forme_a0(kf) = area_a0(kf)**1.5_DOUBLE/abs(mesh%elem(kf)%volume)
              else
                forme_a0(kf) = -1.0_DOUBLE
              end if
            end do
            print*, "AIRE A(0) enregistree sur", mesh%n_elems, " mailles"
          end if
          if (mod(iter, 200) == 0 .or. iter == 1) call area_report(mesh, iter)
          end block
        end if
        if (allocated(fbw_cviol)) print*, "FBWC pas", iter, " mailles portant une", &
          " sous-face violee", nvc, " leur eps_min", epsmin_v, &
          " | eps_min des AUTRES", epsmin_s
        if (nnan > 0) print*, "TRACE   premieres mailles NaN :", trim(lst)
        if (pi_meandiv) then
          print*, "TRACE   garde de vide : ", md_n_vac_step, " sous-faces sur", &
                  md_n_sf_step, " soit", &
                  100.0_DOUBLE*real(md_n_vac_step, DOUBLE) &
                  /max(1.0_DOUBLE, real(md_n_sf_step, DOUBLE)), "%"
          md_n_vac_step = 0; md_n_sf_step = 0
        end if
      end block
      if (md_watch .and. iter <= 3) then
        block
          integer(kind=ENTIER) :: kw
          integer(kind=ENTIER), parameter :: wl(12) = &
            [13780, 13785, 14121, 14122, 14123, 14124, 14125, 14127, 14128, &
             14394, 14395, 14396]
          do kw = 1, 12
            if (wl(kw) <= mesh%n_elems) &
              print*, "WATCH pas", iter, " maille", wl(kw), &
                " Pi-p_c dans [", w_dpi_min(wl(kw)), ",", w_dpi_max(wl(kw)), "]", &
                " production", w_prod(wl(kw))
          end do
        end block
      end if
      if (md_watch) then
        w_prod = 0.0_DOUBLE
        w_dpi_min = 1.0e300_DOUBLE
        w_dpi_max = -1.0e300_DOUBLE
      end if
    end if
    ! Integration de la puissance de bord avec le dt REELLEMENT utilise (celui
    ! passe a compute_rhs_* est celui du pas precedent).
    if (scheme == "vp_pp_pvp" .or. scheme == "vp_pp_pvp_e") &
      pvp_bnd_work = pvp_bnd_work - dt*pvp_bnd_power
    if (scheme == "vp_pp_pvp" .or. scheme == "vp_pp_pvp_e") &
      pvp_bnd_abs_work = pvp_bnd_abs_work + dt*pvp_bnd_abs
    ! defaut de conservation du limiteur, integre avec le dt REELLEMENT utilise
    if (pvp_pos) pos_defect_work = pos_defect_work + dt*pos_defect_step
    if (pi_hat) then
      pih_sq_work   = pih_sq_work + dt*pih_sq_step
      pih_free_work = pih_free_work + dt*pih_free_step
      pih_sq_step = 0.0_DOUBLE; pih_free_step = 0.0_DOUBLE
    end if
    ! CORRECTEUR RK2 : on repart de X^n et U^n, et on avance d'un pas ENTIER avec
    ! la vitesse nodale et le rhs de n+1/2.
    if (rk2_time) then
      do i = 1, mesh%n_vert
        mesh%vert(i)%coord = coord_n(:, i)
      end do
      sol = sol_n
    end if

    call move_mesh(mesh, vp, dt)

    call mpi_memory_exchange_vert(mesh, mpi_send_recv)

    call compute_geometry_mesh(mesh, .true., boundary_2d, do_check=.false.)
    new_sol = sol + dt * rhs
    ! Reset specific volume from geometry to enforce geometric consistency.
    ! This prevents tau from drifting away from the actual volume/mass ratio
    ! due to discretization errors, especially important for large deformations.
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) new_sol(1, i) = mesh%elem(i)%volume / mass(i)
    end do
    sol = new_sol
    call mpi_memory_exchange(mpi_send_recv, mesh%n_elems, 5, sol)
    t = t + dt
    if( mod(iter, 100) == 0 .and. me == 0 ) print*, t, dt
    if( diag_entropy .and. mod(iter, 100) == 0 ) &
      call print_entropy_diag(mesh, t, me, num_procs)

    if( t >= i_sol_vtu * t_max / real(n_sol_vtu - 1) ) then
      write(fln, *) i_sol_vtu
      write(fln, *) "output_"//trim(adjustl(fln))
      call write_sol_lag(mesh, fln, new_sol, vp, pp, gamma_arr, grad_v, grad_p, div_v, alpha_p_arr, p_limited, h_p_arr)
      call write_sol_dat_lag(mesh, fln, new_sol, gamma_arr)
      i_sol_vtu = i_sol_vtu + 1
      ! Diagnostic S PAR FENETRE, puis remise a zero. Cumule sur tout le run il
      ! n est PAS comparable entre deux schemas : un run qui meurt tot
      ! n echantillonne que les premiers instants (facteur 76 mesure entre un run
      ! de 3 pas et un run complet). C est la fenetre qui se compare.
      if (sdiag_tot > 0 .and. me == 0) then
        print*, "SDIAGW t", t, " sous-faces", sdiag_tot, &
                " ]-inf,-1e-1[", sdiag_n(1), " [-1e-1,-1e-2[", sdiag_n(2), &
                " [-1e-2,-1e-3[", sdiag_n(3), " [-1e-3,-1e-6[", sdiag_n(4), &
                " [-1e-6,0[", sdiag_n(5), " >=0", sdiag_n(6)
      end if
      sdiag_n = 0
      sdiag_tot = 0
    end if

    ! CONTROLE DE PRODUCTION de `vp_mdiv`. Les sous-faces sont visitees (mdv_n)
    ! mais la dilatation y est partout NULLE : la correction kappa*delta est
    ! alors inerte et le schema n'est qu'un `vp_fbw` deguise. C'est exactement
    ! la panne qui a eu lieu, et aucun compteur ne la disait -- `mdv_n_cap` a
    ! zero se lisait comme « le plafond ne mord jamais ». On s'arrete, plutot
    ! que de produire un balayage en theta qui multiplie un zero. Meme reflexe
    ! que `set_nml` dans `job.sbatch` : verifier qu'une injection a eu un effet
    ! au lieu de faire confiance au mecanisme.
    ! Le seuil `iter >= 3` est obligatoire : au premier pas `kap_started` est
    ! faux, donc la dilatation n'est legitimement pas encore definie.
    if (pi_mdiv .and. iter >= 3 .and. mdv_n > 0 .and. mdv_n_dltnz == 0) then
      print*, "MDIV ARRET : delta_pc identiquement NUL sur", mdv_n, &
              " sous-faces visitees, au pas", iter, &
              " -- correction inerte, vp_mdiv degenere en vp_fbw."
      print*, "MDIV   verifier que compute_kappa remplit mdv_c et que", &
              " kappa_scatter remplit mdv_se (kap_one ne doit pas les sauter)."
      error stop 7
    end if

    iter = iter + 1
  end do

  i_sol_vtu = -1
  write(fln, *) i_sol_vtu
  write(fln, *) "output_"//trim(adjustl(fln))
  call write_sol_lag(mesh, fln, new_sol, vp, pp, gamma_arr, grad_v, grad_p, div_v, alpha_p_arr, p_limited, h_p_arr)
  call write_sol_dat_lag(mesh, fln, new_sol, gamma_arr)

  if (compute_error) then
    if (init == 10) then
      ! Taylor-Green: the rate is measured on the pressure, rho being uniform here.
      call compute_error_taylor_green(mesh, new_sol, gamma_arr, h_extrude, h_err, &
        l1err, l2err, linferr)
    else
      call compute_error_vortex(mesh, new_sol, t, h_extrude, h_err, l1err, l2err, linferr)
    end if
    if (me == 0) print *, "FINAL t=", t, "h=", h_err, "L1=", l1err, "L2=", l2err, "Linf=", linferr
  end if

  if ((scheme == "vp_pp_pvp" .or. scheme == "vp_pp_pvp_e") .and. me == 0) &
    print*, "PVP_BND_WORK:", pvp_bnd_work, " PVP_BND_ABS:", pvp_bnd_abs_work
  if (diag_lamratio .and. me == 0) then
    block
      integer(kind=ENTIER) :: kk
      integer(kind=ENTIER_D) :: cum
      print*, "LAMR lambda_schema/(rho a) par sous-face, sur", lamr_n, " visites"
      print*, "LAMR   exactement 1 :", lamr_eq1, " soit", &
        100.0_DOUBLE*real(lamr_eq1, DOUBLE)/max(1.0_DOUBLE, real(lamr_n, DOUBLE)), "%"
      print*, "LAMR   > 2 :", lamr_g2, "  > 10 :", lamr_g10, "  > 1e3 :", lamr_g1e3, &
              "   max", lamr_max
      cum = 0
      do kk = 0, 16
        if (lamr_hist(kk) > 0) then
          cum = cum + lamr_hist(kk)
          print*, "LAMR   1e", kk, " a 1e", kk+1, " :", lamr_hist(kk), " cumul", &
            100.0_DOUBLE*real(cum, DOUBLE)/max(1.0_DOUBLE, real(lamr_n, DOUBLE)), "%"
        end if
      end do
    end block
  end if
  if (me == 0) print*, "DTNAN bornes NaN sur tout le run : acoustique", &
    dtb_nan_tot_acou, " volume", dtb_nan_tot_vol, &
    " -- non nul => dt n'etait plus le minimum des bornes a partir de la"
  if (me == 0 .and. dtb_nan_reel > 0) &
    print*, "DTNAN origine : mailles internes", dtb_nan_reel, &
      " | source : sol", dtb_nan_src_sol, " vitesses nodales", dtb_nan_src_vp, &
      " volume", dtb_nan_src_vol
  if (dt_at_step_start .and. me == 0) print*, "DT2PH pas ou la phase 2 a reduit dt :", &
    n_dtviol, " sur", iter-1, "  pire dt_1/dt_2", dtviol_max, &
    " | APRES re-prediction :", n_dtiter2, "  pire dt_2/dt_3", dtiter2_max, &
    " -- NB : sur ces pas le rhs tourne 2 fois, donc les compteurs de diagnostic y comptent double"
  if (pi_mdiv .and. me == 0) then
    print*, "MDIV cumul run : theta", mdiv_theta, " visites", mdv_n, &
            " plafonnees", mdv_n_cap, " soit", &
            100.0_DOUBLE*real(mdv_n_cap, DOUBLE)/max(1.0_DOUBLE, real(mdv_n, DOUBLE)), "%"
    print*, "MDIV   dilatation NON NULLE sur", mdv_n_dltnz, " visites sur", mdv_n, &
            " soit", 100.0_DOUBLE*real(mdv_n_dltnz, DOUBLE)/max(1.0_DOUBLE, real(mdv_n, DOUBLE)), &
            "% -- a ZERO, la correction est inerte et vp_mdiv vaut vp_fbw"
    print*, "MDIV   statuts : intervalle VIDE", mdv_n_empty, &
            " base DEJA HORS", mdv_n_baseout, " aucune contrainte", mdv_n_nocon
    block
      integer(kind=ENTIER) :: kk
      integer(kind=ENTIER_D) :: cum, tot
      tot = 0
      do kk = -16, 16
        tot = tot + mdv_rhist(kk)
      end do
      print*, "MDIV   DISTRIBUTION du rapport kappa_voulu/kappa_max sur", tot, &
              " visites contraintes (le max seul ne dit pas si ca mord partout) :"
      cum = 0
      do kk = -16, 16
        cum = cum + mdv_rhist(kk)
        if (mdv_rhist(kk) > 0) &
          print*, "MDIV     1e", kk, " a 1e", kk+1, " :", mdv_rhist(kk), "  cumul", &
                  100.0_DOUBLE*real(cum, DOUBLE)/max(1.0_DOUBLE, real(tot, DOUBLE)), "%"
      end do
    end block
    print*, "MDIV   kappa retenu max", mdv_kap_max, &
            "  rapport kappa_voulu/kappa_max le pire", mdv_ratio_max
    print*, "MDIV   production : ", kap_n_prod, " evaluations, ", kap_n_prod_neg, &
            " NEGATIVES (anti-dissipation), la pire ", kap_prod_min
  end if
  if (pi_fbw .and. me == 0) then
    block
      integer(kind=ENTIER) :: kk
      integer(kind=ENTIER_D) :: cum
      print*, "FBW cumul run : sous-faces", fbw_n_sf, &
              "  delta2 <= 0 : froid-froid", fbw_n_viol_ff, &
              " chaud-froid", fbw_n_viol_cf, " chaud-chaud", fbw_n_viol_cc
      print*, "FBW   delta2 MINIMAL =", fbw_d2_min, "  -- doit etre > 0"
      print*, "FBW   TAUX de violation par type (compte / visitees) :"
      print*, "FBW     froid-froid ", fbw_n_viol_ff, " /", fbw_nv_ff, " =", &
        100.0_DOUBLE*real(fbw_n_viol_ff, DOUBLE)/max(1.0_DOUBLE, real(fbw_nv_ff, DOUBLE)), "%"
      print*, "FBW     chaud-froid ", fbw_n_viol_cf, " /", fbw_nv_cf, " =", &
        100.0_DOUBLE*real(fbw_n_viol_cf, DOUBLE)/max(1.0_DOUBLE, real(fbw_nv_cf, DOUBLE)), "%"
      print*, "FBW     chaud-chaud ", fbw_n_viol_cc, " /", fbw_nv_cc, " =", &
        100.0_DOUBLE*real(fbw_n_viol_cc, DOUBLE)/max(1.0_DOUBLE, real(fbw_nv_cc, DOUBLE)), "%"
      print*, "FBW   lambda_w / rho a : DISTRIBUTION sur", fbw_n, " evaluations"
      print*, "FBW     exactement 1 (non-regression) :", fbw_n_eq1, " soit", &
              100.0_DOUBLE*real(fbw_n_eq1, DOUBLE)/max(1.0_DOUBLE, real(fbw_n, DOUBLE)), "%"
      print*, "FBW     > 10 :", fbw_n_g10, "   > 1e3 :", fbw_n_g1e3, &
              "   > 1e6 :", fbw_n_g1e6
      cum = 0
      do kk = -2, 16
        if (fbw_hist(kk) > 0) then
          cum = cum + fbw_hist(kk)
          print*, "FBW     1e", kk, " a 1e", kk+1, " :", fbw_hist(kk), &
                  "  cumul", 100.0_DOUBLE*real(cum, DOUBLE) &
                              /max(1.0_DOUBLE, real(fbw_n, DOUBLE)), "%"
        end if
      end do
    end block
  end if
  if (pi_fbd .and. me == 0) then
    print*, "FBD cumul run : sous-faces", fbd_n_sf, " noeuds avec D_p", fbd_n_node, &
            " noeuds a denominateur nul", fbd_n_node_void, " z invalide", fbd_n_zbad
    print*, "FBD   |D_p| max =", fbd_dp_absmax, &
            "  contribution CUMULEE du bord au numerateur =", fbd_num_bnd
    print*, "FBD   delta2 <= 0 (bilan ENTIER) : froid-froid", fbd_n_viol_ff, &
            " chaud-froid", fbd_n_viol_cf, " chaud-chaud", fbd_n_viol_cc, &
            " -- attendu ZERO partout"
    print*, "FBD   delta2 MINIMAL rencontre =", fbd_d2_min
    print*, "FBD   lambda invalides : NaN", fbd_n_znan, " nuls/negatifs", fbd_n_zzero
    print*, "FBD   origine : p NaN", fbd_zb_pnan, " p<=0", fbd_zb_pneg, &
            " tau NaN", fbd_zb_taunan, " tau<=0", fbd_zb_tauneg, &
            " aucun des quatre", fbd_zb_autre
  end if
  if (pvp_pos_sv .and. me == 0) then
    print*, "POSSV cumul run : sous-faces visitees", sv_n_sf, " corrigees", sv_n_corr, &
            " soit", 100.0_DOUBLE*real(sv_n_corr, DOUBLE) &
                     /max(1.0_DOUBLE, real(sv_n_sf, DOUBLE)), "%"
    print*, "POSSV   intervalle VIDE (delta2 <= 0) :", sv_n_infeas, &
            " -- toute valeur non nulle est une anomalie, rien n'a ete corrige la"
    print*, "POSSV   largeur minimale delta2 =", sv_d2_min
    print*, "POSSV   bord PHYSIQUE", sv_n_bnd, " dont corrigees", sv_n_bnd_corr, &
            " -- si 0, la question du bord ne se pose pas"
    print*, "POSSV   faces d'extrusion z", sv_n_bnd_z, " dont corrigees", sv_n_bnd_z_corr, &
            " -- plans miroir, pas des parois"
    print*, "POSSV   pp*(v_p.n) = NaN sur", sv_n_nan, " visites de sous-face"
    print*, "POSSV   dont corrigees PARCE QUE NaN :", sv_n_corr_nan, &
            " -- LAISSES PASSER, non substitues : un NaN doit rester visible"
    print*, "POSSV   ecart/largeur, sous-faces ADMISSIBLES :", sv_n_in, &
            " max", sv_ec_in_max, " moyen", &
            sv_ec_in_tot/max(1.0_DOUBLE, real(sv_n_in, DOUBLE)), " (borne a 0.5)"
    print*, "POSSV   corrections par TYPE d'interface : froid-froid", sv_n_corr_froid, &
            " chaud-chaud", sv_n_corr_chaud, " mixte", sv_n_corr_mixte, &
            " -- seuil p = 1e-6, ambiant Sedov 1e-12 et depot ~2.4e3"
    print*, "POSSV   max(p_l,p_r) sur les corrigees : du plus froid", sv_pmin_corr, &
            " au plus chaud", sv_pmax_corr
    print*, "POSSV   ecart/largeur, sous-faces CORRIGEES   :", sv_n_out, &
            " max", sv_ec_out_max, " moyen", &
            sv_ec_out_tot/max(1.0_DOUBLE, real(sv_n_out, DOUBLE))
  end if
  if ((pi_psvp .or. pi_pp_ep) .and. me == 0) then
    ! Le terme qui separe le bilan d'energie interne de ce schema de celui de
    ! classic. Un COMPTE et une somme, rapportes a la dissipation de classic :
    ! la valeur brute ne dirait pas si le terme mord.
    print*, "PSVP terme (v_c.n)(p_p - Pi_c) : sous-faces", psvp_n_sf, &
            " dont NEGATIVES", psvp_n_neg, " soit", &
            100.0_DOUBLE*real(psvp_n_neg, DOUBLE)/max(1.0_DOUBLE, real(psvp_n_sf, DOUBLE)), "%"
    print*, "PSVP somme algebrique", psvp_sum, " somme des negatives", psvp_sum_neg, &
            " pire contribution", psvp_worst
    print*, "PSVP dissipation de classic somme A lam delta^2 =", psvp_ref, &
            " ; rapport |somme|/ref =", &
            abs(psvp_sum)/max(1.0e-300_DOUBLE, psvp_ref)
    print*, "PSVP lambda de positivite : psvp_lam_iter =", psvp_lam_iter, &
            " marge", psvp_lam_marge
    print*, "PSVP eps* < 0 restants", psvp_n_epsneg, " sur", psvp_n_sf, &
            " soit", 100.0_DOUBLE*real(psvp_n_epsneg, DOUBLE) &
                     /max(1.0_DOUBLE, real(psvp_n_sf, DOUBLE)), "%"
    print*, "PSVP sous-faces HORS fenetre de faisabilite (A<=0, Ma_n>2.67)", &
            psvp_n_infeas, " soit", 100.0_DOUBLE*real(psvp_n_infeas, DOUBLE) &
                     /max(1.0_DOUBLE, real(psvp_n_sf, DOUBLE)), &
            "% -- la ou AUCUN lambda ne convenait" 
  end if
  ! Diagnostic S : production ponderee d energie interne, sous-faces INTERIEURES.
  ! Une case negative = p_p sous la moyenne QUADRATIQUE des p ponderee en 1/lam.
  ! Comptes, pas de max.
  if (pi_pp_vpps .and. me == 0) then
    print*, "VPPSTEL defaut de telescopage interieur somme A (F_l - F_r) =", &
            vpps_tel, " somme des |F| =", vpps_tel_abs, &
            " rapport =", vpps_tel/max(1.0e-300_DOUBLE, vpps_tel_abs)
  end if
  if (tang_n > 0 .and. me == 0) then
    print*, "TANG sous-faces", tang_n, " dont flux tangentiel NON NUL", tang_nact, &
            " soit", 100.0_DOUBLE*real(tang_nact, DOUBLE)/real(tang_n, DOUBLE), "%"
  end if
  if (conv_n > 0 .and. me == 0) then
    print*, "CONVSTAB sous-faces", conv_n, " dont au PLANCHER", conv_nfloor, &
            " soit", 100.0_DOUBLE*real(conv_nfloor, DOUBLE)/real(conv_n, DOUBLE), "%"
  end if
  if (pos_tot > 0 .and. me == 0) then
    print*, "POSMARG e_p - |u*|^2/2 relatif : sous-faces", pos_tot, &
            " <=0", pos_hist(1), " <1e-3", pos_hist(2), " <1e-2", pos_hist(3), &
            " <1e-1", pos_hist(4), " >=1e-1", pos_hist(5), " pire", pos_worst
  end if
  if (ppmc_tot > 0 .and. me == 0) then
    print*, "PPMC ecart |p_mc - p_moy|/p : noeuds", ppmc_tot, &
            " <1e-6", ppmc_hist(1), " <1e-4", ppmc_hist(2), &
            " <1e-2", ppmc_hist(3), " <1e-1", ppmc_hist(4), " >=1e-1", ppmc_hist(5)
  end if
  if (sdiag_tot > 0 .and. me == 0) then
    print*, "SDIAG sous-faces interieures", sdiag_tot, &
            " | S/(p.a) : <-1e-1", sdiag_n(1), " <-1e-2", sdiag_n(2), &
            " <-1e-3", sdiag_n(3), " <-1e-6", sdiag_n(4), &
            " <0", sdiag_n(5), " >=0", sdiag_n(6)
    print*, "SDIAG part de S < 0 :", &
            100.0_DOUBLE*real(sdiag_n(1)+sdiag_n(2)+sdiag_n(3) &
                              +sdiag_n(4)+sdiag_n(5), DOUBLE) &
            /max(1.0_DOUBLE, real(sdiag_tot, DOUBLE)), &
            "%  dont |S| > 1e-2 :", &
            100.0_DOUBLE*real(sdiag_n(1)+sdiag_n(2), DOUBLE) &
            /max(1.0_DOUBLE, real(sdiag_tot, DOUBLE)), "%"
  end if
  if (pvp_pos .and. me == 0) then
    print*, "POS cumul run : sous-faces ecretees", pos_n_clip, " sur", pos_n_sf, &
            " visitees, soit", 100.0_DOUBLE*real(pos_n_clip, DOUBLE)/max(1.0_DOUBLE, real(pos_n_sf, DOUBLE)), "%"
    print*, "POS cumul run : noeuds infaisables", pos_n_infeas, " repartis sur", &
            pos_n_steps_infeas, "pas ; defaut integre en temps", pos_defect_work
    print*, "POS cumul run : sous-faces de bord ecretees", pos_n_clip_bnd, &
            " ; max par NOEUD sur le run", pos_n_clip_max
    if (allocated(pos_clip_cell)) then
      block
        integer(kind=ENTIER) :: ic, kk
        real(kind=DOUBLE) :: rho_c
        print*, "POS : les 12 mailles les plus ecretees (x y z rho n_clip)"
        do kk = 1, 12
          ic = maxloc(pos_clip_cell, 1)
          if (pos_clip_cell(ic) <= 0) exit
          rho_c = 1.0_DOUBLE/sol(1, ic)
          print*, mesh%elem(ic)%coord, rho_c, pos_clip_cell(ic)
          pos_clip_cell(ic) = -1
        end do
      end block
    end if
  end if

  if (pi_hat .and. me == 0) then
    print*, "PIH theta =", theta_pi
    print*, "PIH production : part en CARRES", pih_sq_work, " part LIBRE", pih_free_work, &
            " total", pih_sq_work + pih_free_work
    print*, "PIH residu nodal max |v_p . somme A Pi n| =", pih_res_max
    if (pi_mp) then
      print*, "MP c_visc =", mp_cvisc
      print*, "MP production NODALE cumulee (doit etre >= 0) =", mp_prod_node
      print*, "MP mailles a production NEGATIVE :", mp_n_cell_neg, " sur", mp_n_cell_tot, &
              " ; la plus negative =", mp_prod_min_cell
      print*, "MP noeuds de BORD sans stabilisation :", mp_n_bnd_node, " sur", mp_n_node_tot
    end if
    if (pi_mach_damp) print*, "PIH theta_M rencontre sur le run : min =", pih_tm_min, &
      " max =", pih_tm_max
  end if

  ! Ce diagnostic doit s'imprimer pour N'IMPORTE QUEL schema, en particulier pour
  ! 'classic' ou pi_hat est faux : il etait a tort imbrique dans le bloc pi_hat.
  if (pi_meandiv .and. me == 0) then
    print*, "MD garde de vide declenchee", md_n_vac_tot, "fois sur le run"
    print*, "MD c_visc =", md_cvisc
    print*, "MD bilan nodal max |somme A Pi n| aux noeuds INTERIEURS =", md_res_max
    print*, "MD residu par direction  x =", md_res_dir(1), " y =", md_res_dir(2), &
            " z =", md_res_dir(3)
    print*, "MD somme mu (div v)^2 A h (dissipation effective) =", md_sq_div
    print*, "MD somme lambda A delta^2 (ce que vp dissiperait)  =", md_sq_delta
    print*, "MD rapport div/delta =", md_sq_div/max(1.0e-300_DOUBLE, md_sq_delta)
    print*, "MD mailles a production NEGATIVE :", md_n_cell_neg, " sur", md_n_cell_tot, &
            " ; la plus negative =", md_prod_min
  end if
  if (pi_kappa .and. me == 0) then
    print*, "KAP plancher =", kap_floor, "  kappa_one =", kap_one, &
            "  kappa_fixe =", kap_fixed
    print*, "KAP kappa rencontre : min =", kap_min, " max =", kap_max, &
            " moyenne =", kap_sum/max(1.0_DOUBLE, real(kap_n, DOUBLE))
    print*, "KAP coins au plancher :", kap_n_floor, " sur", kap_n, " soit", &
            100.0_DOUBLE*real(kap_n_floor, DOUBLE)/max(1.0_DOUBLE, real(kap_n, DOUBLE)), "%"
    print*, "KAP max |(div v)_pc| rencontre =", kap_divmax
    print*, "KAP production totale (somme A kappa lambda delta^2) =", kap_prod_tot
    print*, "KAP termes a production NEGATIVE :", kap_n_prod_neg, " sur", kap_n_prod, &
            " ; le plus negatif =", kap_prod_min
    print*, "KAP bilan nodal max |somme A Pi n| aux noeuds INTERIEURS =", kap_res_max
    if (kap_loc) then
      print*, "KAPLOC plancher local actif, mode = ", trim(kap_loc_mode), "  C =", kap_loc_c
      print*, "KAPLOC plancher par noeud : min =", kfl_min, " max =", kfl_max, &
              " moyenne =", kfl_sum/max(1.0_DOUBLE, real(kfl_n, DOUBLE))
      print*, "KAPLOC kappa EFFECTIF (= max(senseur, plancher)) : min =", keff_min, &
              " max =", keff_max, " moyenne =", keff_sum/max(1.0_DOUBLE, real(keff_n, DOUBLE))
      print*, "KAPLOC coins RELEVES par le plancher :", keff_n_raised, " sur", keff_n, &
              " soit", 100.0_DOUBLE*real(keff_n_raised, DOUBLE)/max(1.0_DOUBLE, real(keff_n, DOUBLE)), "%"
      print*, "KAPLOC noeuds a plancher SATURE (= 1) :", kfl_n_one, " sur", kfl_n, &
              " soit", 100.0_DOUBLE*real(kfl_n_one, DOUBLE)/max(1.0_DOUBLE, real(kfl_n, DOUBLE)), "%"
    end if
    ! La question qui decide : kappa etant explicite, est-il au plancher AU
    ! MOMENT ou le front passe ? La distribution globale n'y repond pas.
    print*, "KAP seuil de front |dp|/p >", kap_front_jump, " ; saut max rencontre =", &
            kap_front_jmax
    print*, "KAP kappa MINIMAL sur les coins de front =", kap_front_min
    print*, "KAP coins de front au plancher :", kap_n_front_floor, " sur", kap_n_front
  end if
  if (pi_divdirect .and. me == 0) then
    print*, "VPDIV mode = ", trim(divd_mode), "  coef =", divd_coef, "  theta =", divd_theta
    print*, "VPDIV production totale (somme sur les mailles) =", divd_prod_tot
    print*, "VPDIV mailles a production NEGATIVE :", divd_n_cell_neg, " sur", divd_n_cell, &
            " soit", 100.0_DOUBLE*real(divd_n_cell_neg, DOUBLE) &
            /max(1.0_DOUBLE, real(divd_n_cell, DOUBLE)), "%"
    print*, "VPDIV la plus negative =", divd_prod_min
    print*, "VPDIV forme non lineaire (isentropique) =", divd_nl
    print*, "VPDIV sous-faces a Pi < 0 :", divd_n_pineg, " sur", divd_n_pi, " soit", &
            100.0_DOUBLE*real(divd_n_pineg, DOUBLE)/max(1.0_DOUBLE, real(divd_n_pi, DOUBLE)), "%"
    print*, "VPDIV distribution de |M_div| : >0.1", divd_md(1), " >0.5", divd_md(2), &
            " >1", divd_md(3), " >2.857", divd_md(4), " >5", divd_md(5), &
            " >10", divd_md(6), " >20", divd_md(7), "  sur", divd_n_pi
    print*, "VPDIV sous-faces a formation de VIDE (base <= 0) :", divd_n_vide
    print*, "VPDIV noeuds a matrice SINGULIERE :", divd_n_sing, &
            "  dont bord", divd_sing_bnd, " interieur", divd_sing_int
    ! RESULTAT et non controle : la production de cette fermeture n'est pas un
    ! carre mais un carre plus un terme indefini. Un compte non nul est attendu
    ! et sa valeur mesure ce que coute la perte du carre.
  end if
  if (diag_divratio .and. me == 0) then
    print*, "DIVRATIO somme lambda A delta^2 (ce que vp dissipe)  =", dr_sq_delta
    print*, "DIVRATIO somme mu (div v)^2 V (ce qu'une divergence) =", dr_sq_div
    print*, "DIVRATIO rapport div/delta =", dr_sq_div/max(1.0e-300_DOUBLE, dr_sq_delta)
  end if

  call MPI_FINALIZE(mpi_ierr)
end program main
