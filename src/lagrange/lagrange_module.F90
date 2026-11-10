module lagrange_module
  use mpi
  use precision_module
  use mpi_module
  use mesh_module
  use mesh_reading_module
  use mesh_geometry_module
  use mesh_connectivity_module
  use linear_solver_module, only : pseudo_inverse_inplace_lapack

  implicit none

  real(kind=DOUBLE), parameter :: gamma = 7.0_DOUBLE/5.0_DOUBLE

  ! .true. = impedance acoustique pure z_c = rho_c a_c dans la matrice de coin, comme
  ! EUCCLHYD d'origine. .false. (defaut) = comportement historique de ce code.
  logical, public :: lambda_acoustic_only = .false.
  ! Mur : 'constraint' (contre-pression P_p) | 'symmetric' (defaut ; solveur nodal double par le
  ! miroir de chaque sous-face selon b = B_p/|B_p| : M_p + R M_p R, R_p + R R_p, R = I - 2 b b).
  character(len=16), public :: wall_mode = 'symmetric'
  ! > 0 : eps des poids non lineaires des gradients = grad_eps_vol * V_c (volume de la maille).
  ! 0 : ancien comportement (eps absolu, tiny() ou grad_eps_weno).
  real(kind=DOUBLE), public :: grad_eps_vol = 1.0e-8_DOUBLE
  ! Facteur du terme d'extrapolation EN TEMPS de la reconstruction (les 11 sites
  ! en 0.5*dt). 1 = MUSCL-Hancock ; 0 = reconstruction SPATIALE seule, pour RK2.
  real(kind=DOUBLE), public :: rec_dt_fac = 1.0_DOUBLE
  ! Poids non lineaire de la moyenne des gradients nodaux pour la VITESSE.
  ! 'gradp' (DEFAUT) : le meme poids que la pression, en |grad p|.
  ! 'divv'  : poids en |trace(grad v)|, c'est-a-dire en |div v|.
  character(len=8), public :: grad_weight_v = 'gradp'
  ! Terme de saut X de la 3e quantite nodale du schema 'vp_pp_ppvp' (voir
  ! compute_rhs_lagrange_vp_pp). 'vp_dv' : X = (u_p.n)(v_r - v_l), seule forme qui
  ! soit a la fois homogene a v^2 ET anti-invariante par retournement de la face.
  ! Avec ce choix q_bar = v_node * p_bar : c'est exactement la construction de pp,
  ! portee par la vitesse nodale de la face. 'dp' : le saut porte sur p au lieu de
  ! v, rendu homogene par les impedances -- conserve pour comparaison.
  character(len=16), public :: ppvp_jump_mode = 'vp_dv'
  ! .true. : ppvp est UNE valeur par noeud (moyenne en 1/lambda sur les sous-faces),
  ! comme pp. Mesure : une sous-face dont la vitesse nodale est orthogonale a la
  ! normale ne balaie aucun volume, donc son flux d'energie consistant est nul --
  ! mais la valeur nodale lui est appliquee quand meme. Sur le petit disque, 48
  ! sous-faces dans ce cas, et exactement 48 valeurs NaN : dans un ambiant a 1e-12
  ! d'energie interne, ce flux parasite la rend negative immediatement.
  ! .false. : ppvp est evalue PAR SOUS-FACE, q_bar = p_bar * v_et, qui s'annule de
  ! lui-meme quand la sous-face ne travaille pas.
  logical, public :: ppvp_nodal_average = .true.

  ! --- diagnostic de production d'entropie -----------------------------------
  ! Pour une maille c, en partant des equations semi-discretes et de Gibbs :
  !   m T ds/dt = somme_f A_f (p_c - Pi_f) (u_p - u_c).n_f      [entropy_prod]
  ! Avec la fermeture UNILATERALE Pi_f = p_c - lambda_c (u_p-u_c).n_f cela vaut
  !   somme_f A_f lambda_c [(u_p-u_c).n_f]^2 >= 0               [entropy_ref]
  ! qui est positif par construction. 'vp_pp' utilise au contraire UNE pression
  ! nodale pp commune aux deux cotes, donc (p_c - pp) n'a aucun signe impose et
  ! rien ne garantit entropy_prod >= 0. Le rapport entropy_prod/entropy_ref dit
  ! de combien on s'ecarte de la dissipation qu'on aurait avec l'unilateral.
  logical, public :: diag_entropy = .false.
  real(kind=DOUBLE), allocatable, public :: entropy_prod(:), entropy_ref(:)

  ! Energie nodale du schema 'vp_ep', une valeur par noeud (cf.
  ! compute_rhs_lagrange_vp_ep). Publique pour pouvoir etre relue/ecrite.
  real(kind=DOUBLE), allocatable, public :: ep_node(:)

  ! --- variantes 'vp_ep_wip' -------------------------------------------------
  ! vp_ep dissipe au CONTACT parce que son terme correctif lam_c (e_p - e_c) est
  ! indexe sur le saut d'ENERGIE, maximal la-bas, alors que celui de 'classic',
  ! lam_c (u_p.n - v_c), est indexe sur le saut de VITESSE, nul la-bas. Mesure sur
  ! le Sod cylindrique : de/e = 0.048 contre dv/v = 0.007 et dp/p = 0.002 de part
  ! et d'autre du contact ; largeur du contact 0.025 contre 0.010 pour classic.
  !
  ! On remplace donc lam_c par mu_c = lam_c * phi_c, avec phi_c -> 0 au contact.
  ! CRUCIAL : mu_c doit remplacer lam_c A LA FOIS dans le flux ET dans la
  ! definition de e_p, sinon la telescopie qui assure la conservation exacte est
  ! rompue. Avec la meme substitution des deux cotes,
  !   e_p = somme A[(mu_l e_l + mu_r e_r) + (p_l v_l - p_r v_r)] / somme A(mu_l+mu_r)
  ! et somme A (pv_l - pv_r) = 0 se refait a l'identique, quel que soit mu.
  !
  !   'dv'   phi = |u_p.n - v_c| / (|u_p.n - v_c| + a_c)      saut de vitesse
  !   'dp'   phi = |p_r - p_l|   / (|p_r - p_l|   + p_c)      saut de pression
  !   'none' phi = 1                                           = vp_ep d'origine
  ! Defaut 'dp' : mesure sur le Sod cylindrique, largeur du contact 0.0049 contre
  ! 0.0103 pour classic et 0.0252 pour vp_ep. 'dv' DEGRADE (0.0333) : au contact
  ! phi_dv tombe au plancher, mu -> 0, et e_p devient une indetermination 0/0
  ! dominee par le terme de flux.
  character(len=16), public :: vp_ep_wip_mode = 'dp'
  ! Plancher sur phi : a phi = 0 partout autour d'un noeud, e_p est indetermine
  ! (0/0). Le plancher garde le systeme bien pose et fixe la dissipation residuelle.
  real(kind=DOUBLE), public :: vp_ep_wip_floor = 1.0e-3_DOUBLE
  logical, public :: vp_ep_wip = .false.
  ! Plafond de positivite sur la correction d'energie. Mesure sur Noh : le terme
  ! mu_c (e_p - e_c) vide les mailles du quasi-vide (p0 = 1e-6) jusqu'a
  ! e_int = -1.07e-6 des le dump 2, AVANT tout NaN, qui n'apparaissent qu'ensuite.
  ! On borne donc |mu_c (e_p - e_c)| par clip * lam_c * e_c, homogene et local.
  ! 0 = pas de plafond, ET C'EST LE DEFAUT : teste a 1.0 et 0.2, le plafond ne
  ! debloque NI Noh NI Saltzman et CASSE sedov_quad, qui passait (189 NaN contre 0).
  ! Etant unilateral il rompt la telescopie et introduit un desequilibre nouveau.
  ! Conserve uniquement pour documenter la tentative.
  real(kind=DOUBLE), public :: vp_ep_wip_clip = 0.0_DOUBLE


  ! --- limiteur convexe du schema 'vp_pp_pvp' ---------------------------------
  ! Melange par noeud entre le flux nodal (pp, pvp) et le flux 'classic'
  ! (Pi_l, Pi_r, Pi (u_p.n)), dont la production d'entropie vaut somme A lam d^2
  ! et est donc positive par construction :
  !   flux = theta * flux_nodal + (1 - theta) * flux_classic
  ! theta = 1 redonne vp_pp_pvp, theta = 0 redonne classic.
  logical, public :: pvp_limit = .false.
  ! Construction du produit nodal (pv)_p.
  ! 'god' (defaut) : solveur acoustique par sous-face, puis systeme n (x) n.
  ! 'bal' : on annule le bilan d'energie interne du sous-element pc, projete sur
  !   les normales et somme sur le noeud,
  !       somme_c A_pc n_pc [ (pv)_p.n_pc - p_p (v_c.n_pc) ] = 0
  !   soit   N_p (pv)_p = p_p somme_c A_pc (n_pc (x) n_pc) v_c,  N_p = somme_c A_pc n_pc (x) n_pc.
  !   L'annulation sous-face par sous-face est impossible (il faudrait
  !   (pv)_p.n = p_p v_l ET = p_p v_r simultanement) ; c'est donc la projection
  !   qui est annulee. En ecoulement uniforme v_c = u et (pv)_p = p_p u, donc la
  !   consistance est preservee.
  character(len=8), public :: pvp_energy_mode = 'god'

  ! --- schema 'vp_pp_pos' -----------------------------------------------------
  ! Pas de troisieme systeme nodal. On garde v_p et p_p de vp_pp et on ECRETE le
  ! flux physique p_p (v_p.n) pour garantir eps* >= 0 des deux cotes, puis on
  ! redistribue le surplus sur le noeud pour restaurer la conservation. La ou
  ! aucune borne ne mord, c'est exactement vp_pp.
  logical, public :: pvp_pos = .false.
  ! --- schema 'vp_pp_psvp' ----------------------------------------------------
  ! Le flux est  ( -v_p.n ,  p_p n ,  Pi_c (v_p.n) ) : la quantite de mouvement
  ! de vp_pp (pression NODALE, univaluee) et l'energie de classic (pression de
  ! sous-face UNILATERALE Pi_c = p_c - lam_c (v_p - v_c).n_c).
  !
  ! Les TROIS conservations tiennent, par trois mecanismes DIFFERENTS :
  !   volume            v_p.n est univalue au noeud ;
  !   quantite de mvt   p_p est univaluee, donc somme_c A p_p n = p_p somme A n = 0,
  !                     l'identite geometrique du noeud interieur ;
  !   energie TOTALE    le flux est BIVALUE, et ce n'est pas lui qui telescope :
  !                     somme_mailles = somme_noeuds v_p . (somme_c A Pi_c n_c), et
  !                     somme_c A Pi_c n_c = 0 est le bilan nodal de classic -- que
  !                     v_p resout deja ici, le systeme en vitesse de vp_pp etant
  !                     EXACTEMENT celui de classic (meme Mp, meme Rp).
  !
  ! LE PRIX, et il n'est pas dans la conservation. En retranchant le budget
  ! cinetique, qui voit p_p alors que l'energie voit Pi_c :
  !   m dEint/dt = [-p_c dV/dt + somme A lam_c delta^2]  +  somme_f A (v_c.n)(p_p - Pi_c)
  ! Le terme en trop s'annule au REPOS -- donc Sedov a t=0 est intact -- mais n'a
  ! AUCUN signe des que la maille bouge. La production d'entropie cesse donc d'etre
  ! une somme de carres, contrairement a classic. C'est la seule chose que ce
  ! schema perd, et elle se mesure : voir les compteurs psvp_* ci-dessous.
  logical, public :: pi_psvp = .false.
  ! --- schema 'vp_pp_ep' ------------------------------------------------------
  ! Flux = ( -v_p.n , p_p n , (pv)*_s ) : la quantite de mouvement de vp_pp
  ! (pression NODALE, univaluee) et l'energie de vp_ep (flux de sous-face
  ! UNILATERAL construit sur l'energie NODALE e_p).
  !
  ! Les TROIS conservations tiennent, par trois mecanismes DIFFERENTS :
  !   volume            v_p.n univalue au noeud ;
  !   quantite de mvt   p_p univaluee, donc somme_c A p_p n = p_p somme A n = 0 ;
  !   energie TOTALE    le flux est BIVALUE, et c'est e_p qui telescope :
  !                     annuler somme_f A[(pv)*_r - (pv)*_l] donne
  !                       e_p = somme A[mu_l e_l + mu_r e_r + (p_l v_l - p_r v_r)]
  !                             / somme A (mu_l + mu_r)
  !                     et cette definition ne fait intervenir le flux de
  !                     QUANTITE DE MOUVEMENT nulle part. Remplacer Pi_c par p_p
  !                     ne touche donc pas la conservation de l'energie -- c'est
  !                     la difference structurelle avec 'vp_pp_psvp', ou c'etait
  !                     le bilan nodal de classic qui telescopait.
  !
  ! mu_s = lambda_s (phi = 1, vp_ep_wip inactif) : l'impedance elle-meme, donc
  ! aucun coefficient libre, et le denominateur somme A(lam_l+lam_r) ne peut pas
  ! s'annuler -- l'indetermination 0/0 que 'vp_ep_wip' doit border n'existe pas.
  !
  ! CE QU'IL FAUT SURVEILLER, et c'est ce qui a tue 'vp_pp_psvp' : le budget
  ! cinetique voit p_p alors que l'energie voit (pv)*_s, donc
  !   m dEint/dt = -somme A (pv)*_c + somme A p_p u_c
  ! au lieu de ... + somme A Pi_c u_c. L'ecart est somme_f A (p_p - Pi_c) u_c,
  ! exactement le terme mesure a -1.0004 fois la dissipation de classic sur psvp.
  ! Ici le flux d'energie porte sa PROPRE dissipation -mu(e_p - e_c), qui n'est
  ! pas le lam delta^2 de classic, donc rien n'impose la meme annulation -- mais
  ! on la MESURE avec les compteurs psvp_* au lieu de la supposer absente.
  logical, public :: pi_pp_ep = .false.
  ! 'vp_pp_ep_i' : la MEME construction, mais le scalaire relaxe est l'energie
  ! INTERNE eps = E - |u|^2/2 au lieu de l'energie totale E, dans le flux ET
  ! dans la definition de X_p. Le telescopage ne depend que de la forme
  ! `+-mu_s (X_p - X_s)` avec un X_p commun aux deux cotes, donc il vaut pour
  ! n'importe quel scalaire par maille : la conservation de l'energie TOTALE
  ! reste exacte au bit.
  !
  ! MOTIVATION, venue du cote NS et pas d'ici (subfv-ab, cylindre M17.6) :
  ! `lam (E_p - E_c)` relaxe l'energie TOTALE, dont l'ecart a travers une couche
  ! limite est domine par la part CINETIQUE -- environ 87 fois l'interne a
  ! M17.6. Le terme pompe donc rho a Delta(|u|^2/2) dans la maille parietale ou
  ! il devient de la chaleur, sans relaxation correspondante de la quantite de
  ! mouvement : C_H surchauffe d'un facteur 6 contre LAURA, et le cas triangulaire
  ! plante par energie interne negative.
  !
  ! PREDICTION ECRITE AVANT LA MESURE : sur Sedov et Gresho, qui sont a bas Mach
  ! la ou ce terme travaille, la substitution doit changer les resultats de
  ! QUELQUES POUR CENT seulement et le Sedov doit rester mou. Un changement
  ! massif ici invaliderait le mecanisme propose cote NS.
  logical, public :: pi_ep_interne = .false.
  ! 'vp_pp_epsp' : energie INTERNE au noeud, forme canonique. On n'impose plus
  ! e*_c = e_p (energie TOTALE commune) mais eps*_c = eps_p. La relation de saut
  ! lagrangienne -lam_c (e*_c - e_c) = (pv)*_c - p_c v_c donne alors
  !     (pv)*_c = p_c v_c -/+ lam_c ( eps_p + |v*_c|^2/2 - e_c )
  ! avec v*_c = v_c + delta_c n la vitesse etoile ACOUSTIQUE de la maille,
  ! delta_l = (p_l - p_p)/lam_l et delta_r = (p_p - p_r)/lam_r. En posant
  !     X_c = e_c - |v*_c|^2/2 = eps_c - v_c delta_c - delta_c^2/2
  ! (MEME signe des deux cotes : l'asymetrie est deja dans delta_c), la
  ! conservation nodale definit eps_p comme la moyenne ponderee des X_c, et il
  ! vient eps*_l = eps*_r = eps_p IDENTIQUEMENT.
  ! INTERET : la positivite des deux etats etoiles se reduit alors a UNE seule
  ! condition scalaire par noeud, eps_p > 0, au lieu d'une par cote de sous-face.
  ! A DISTINGUER de pi_ep_interne, qui moyenne eps_c SANS la correction de
  ! glissement : cette variante-la donne e*_l /= e*_r et perd la propriete.
  logical, public :: pi_pp_epsp = .false.
  ! 'vp_pp_eps2' / 'vp_pp_eps2m' : DEUX energies internes etoiles differentes.
  ! On n'impose plus une valeur nodale commune : chaque cote garde la sienne, a
  ! un seul scalaire nodal pres. Toute la dissipation du flux d'energie passe
  ! alors par l'energie CINETIQUE (le carre du glissement), et le terme
  ! lam_c (eps_p - eps_c) -- une conduction thermique a vitesse acoustique, qui
  ! donnait +348 % de flux parietal en NS -- DISPARAIT.
  ! En imposant eps*_c = eps_c, la relation de saut donne directement
  !     (pv)*_l = p_p v_l - D_l^2/(2 lam_l)   et   (pv)*_r = p_p v_r + D_r^2/(2 lam_r)
  ! avec D_l = p_l - p_p et D_r = p_p - p_r. Ce flux n'est PAS conservatif au
  ! noeud ; un unique scalaire kappa le rend conservatif :
  !   ADDITIF      eps*_c = eps_c + kappa, kappa = -N / somme A (lam_l + lam_r)
  !   MULTIPLICATIF eps*_c = eps_c (1 + kappa), kappa = -N / somme A (lam_l eps_l + lam_r eps_r)
  !   N = somme A [ p_p (v_r - v_l) + D_r^2/(2 lam_r) + D_l^2/(2 lam_l) ]
  ! L'ADDITIF laisse les DIFFERENCES d'energie interne inchangees, donc il ne
  ! lisse aucun gradient : c'est le candidat pour le flux de chaleur. Le
  ! MULTIPLICATIF a une positivite plus simple (kappa > -1) donc plus robuste.
  ! Dans les deux cas ep_node porte kappa et non une energie.
  logical, public :: pi_pp_eps2 = .false.
  logical, public :: pi_pp_eps2m = .false.
  ! 'vp_pp_epv' : le flux d'energie est la pression NODALE fois la vitesse
  ! ETOILE DE CHAQUE COTE, sans AUCUNE energie nodale :
  !     (pv)*_c = p_p * ( v_c.n +/- (p_c - p_p)/lam_c )
  ! Les deux energies internes etoiles sont alors DIFFERENTES et ne dependent
  ! chacune que de leur propre maille :
  !     eps*_c = eps_c + (p_p^2 - p_c^2) / (2 lam_c^2)
  ! Trois proprietes, toutes verifiees numeriquement :
  !  - CONSERVATION NODALE AUTOMATIQUE. somme A[(pv)_r - (pv)_l] = p_p somme A(v*_r - v*_l),
  !    et v*_r - v*_l = (1/lam_l + 1/lam_r)(p_p - p_bar), donc la somme s'annule
  !    EXACTEMENT parce que p_p est la moyenne ponderee en 1/lam des p_bar. C'est
  !    le solveur de pression nodale lui-meme qui conserve l'energie.
  !  - POSITIVITE AUTOMATIQUE avec le lambda CLASSIQUE : eps*_c > 0 s'ecrit
  !    lam_c^2 > (p_c^2 - p_p^2)/(2 eps_c), au pire p_c^2/(2 eps_c), qui est la
  !    condition d'Einfeldt impliquee par lam >= rho a (convexite de l'EOS).
  !    Zero violation sur 200 000 tirages. Aucune borne n'est necessaire.
  !  - PLUS DE TERME lam(e_p - e_c), donc plus de conduction thermique a vitesse
  !    acoustique : c'est elle qui donnait +346 % de flux parietal en NS.
  ! RESERVE : en DETENTE (p_p < p_c) on a eps* < eps_c, le schema refroidit. A
  ! verifier du cote de la production d'entropie.
  logical, public :: pi_pp_epv = .false.
  ! (pv)*_c = Pi_c (v_p.n) : pression de sous-face CLASSIQUE et vitesse NODALE
  ! dans le flux d ENERGIE, alors que le flux de QUANTITE DE MOUVEMENT garde la
  ! pression NODALE p_p. INCOHERENT PAR CONSTRUCTION, et volontairement :
  !  - eps* devient dependant du REFERENTIEL, au taux (Pi_c - p_p)/lam_c ;
  !  - l energie n est plus conservee, car l equilibre nodal porte sur p_p et
  !    non sur Pi_c, donc somme A Pi_c n n est pas nul.
  ! A MESURER, PAS A CROIRE : les deux diagnostics vpps_* ci-dessous sont
  ! obligatoires avec ce schema.
  logical, public :: pi_pp_vpps = .false.
  ! (pv)*_c = p_p (v_p.n) : pression NODALE et vitesse NODALE dans le flux
  ! d energie. Seule combinaison qui ait a la fois
  !   - la vitesse NODALE dans le flux d energie (structure de dissipation de
  !     classic, celle qui donne le bon flux de chaleur),
  !   - la MEME pression que le flux de quantite de mouvement, donc invariance
  !     galileenne et, dans eps*, aucun terme lineaire en la vitesse de maille,
  !   - F_l = F_r identiquement, donc telescopage et conservation triviaux.
  logical, public :: pi_pp_ppvn = .false.
  ! pi_pp_wdiv : dans la pression nodale, PONDERER les sauts de vitesse par
  ! 1/lambda comme le sont deja les pressions. Aujourd hui
  !   p_p = [somme A (p_l/lam_l + p_r/lam_r) - somme A (v_r - v_l)] / somme A invlam
  ! et le second terme, a poids UNIFORME, est exactement la divergence nodale
  ! discrete : il s annule sur un champ solenoidal, d ou un vortex de Gresho
  ! preserve EXACTEMENT -- et aucune dissipation au choc.
  ! La formule actuelle annule la MOYENNE du desaccord de vitesse etoile,
  ! somme A (u*_r - u*_l) = 0. Les MOINDRES CARRES ferment le meme systeme
  ! sur-determine en MINIMISANT somme A (u*_r - u*_l)^2 : la stationnarite donne
  ! somme A invlam (u*_r - u*_l) = 0, donc p_p = somme A invlam (S - dv) /
  ! somme A invlam^2. La ponderation en 1/lambda TOMBE de la minimisation, elle
  ! n est pas inseree, et la somme en dv n est alors plus une divergence des que
  ! les impedances varient d une sous-face a l autre. Or cette dispersion est elle-meme en O(Ma^2) dans un ecoulement
  ! bas Mach (rho et p y fluctuent en Ma^2) et d ordre UN au choc (lambda y
  ! capte les sauts). L ecart a la divergence est donc en O(Ma^2) a bas Mach et
  ! O(1) au choc, SANS aucun coefficient : mesure, le rapport ecart/Ma^2 vaut
  ! 0,21 constant de Ma = 1e-4 a 1e-1, et l ecart monte a 0,14 / 0,33 / 0,47
  ! aux Mach 1 / 3 / 10.
  ! La normalisation par <invlam> garantit la NON-REGRESSION : a impedances
  ! uniformes le poids vaut 1 et on retrouve la formule actuelle au bit pres.
  logical, public :: pi_pp_wdiv = .false.
  ! pi_pp_veps : la pression nodale fermee par une condition sur l ENERGIE
  ! INTERNE intermediaire, somme_f A (eps*_r - eps*_l) = 0, proposee par Vincent.
  ! Avec mu = lambda le facteur de relaxation vaut 1, donc e*_s = e_p des DEUX
  ! cotes et eps*_r - eps*_l = (u*_l^2 - u*_r^2)/2 : e_p s ELIMINE identiquement
  ! (verifie a 6e-14). La condition porte donc sur les VITESSES etoile,
  !   somme A (u*_r^2 - u*_l^2) = 0,
  ! et comme u*_r^2 - u*_l^2 = (u*_r - u*_l)(u*_r + u*_l), c est le desaccord de
  ! vitesse PONDERE PAR LA VITESSE ETOILE MOYENNE -- une ponderation DERIVEE d une
  ! exigence physique au noeud, et non d une minimisation.
  ! L equation est quadratique en p_p. Son terme quadratique a pour coefficient
  ! somme A (1/lam_r^2 - 1/lam_l^2), une DIFFERENCE qui s annule par symetrie, et
  ! son terme lineaire 2 somme A (u*_r/lam_r + u*_l/lam_l), une SOMME qui ne s
  ! annule pas (0,336 de son echelle). On resout donc la quadratique quand elle
  ! est franche et la lineaire sinon, et on retient la racine la plus proche de
  ! la fermeture de moyenne -- qui reste la valeur de repli si rien ne convient.
  logical, public :: pi_pp_veps = .false.
  ! pi_conv_stab : IMPEDANCE DE STABILISATION CONVECTIVE dans la pression de
  ! sous-face de `classic`. La force exercee sur une maille par ses sous-faces
  ! au noeud vaut EXACTEMENT (identite verifiee a 1e-17)
  !     somme_f A Pi_{c,f} n_f = p_c N_pc - lambda_c M_pc (v_p - v_c),
  ! avec N_pc la normale de coin et M_pc = somme_f A n(x)n. Le second terme est
  ! une TRAINEE lineaire vers la vitesse nodale, de coefficient lambda = rho a,
  ! donc a la vitesse du SON. Dans un vortex v_p - v_c = O(h grad v) ne s annule
  ! jamais : la trainee relative est O(Ma) et, rapportee a la fluctuation
  ! physique de pression O(Ma^2), elle DIVERGE en 1/Ma. C est ce qui mange le
  ! vortex de Gresho, et `vp_pp` ne l a pas du tout puisque sa force vaut
  ! p_p N_pc tout court.
  ! Remplacer l impedance de STABILISATION par le taux CONVECTIF
  !     mu_c = rho_c |v_p.n - v_c.n|
  ! rend la trainee O(Ma^2), soit l ordre de la fluctuation physique (mesure :
  ! 1,24e-5 / 1,24e-7 / 1,24e-9 aux Mach 1e-1 / 1e-2 / 1e-3, contre
  ! 3,5e-3 / 3,5e-4 / 3,5e-5 au taux acoustique). Au choc |dv| ~ a et on
  ! retrouve exactement la dissipation de `classic`. Aucun coefficient libre :
  ! c est une difference de vitesses.
  ! Les VITESSES D ONDE et le pas de temps gardent rho a : on separe les deux
  ! roles de lambda, on n affaiblit pas la borne de positivite.
  ! L equation nodale devient non lineaire en v_p ; on utilise la vitesse nodale
  ! du PAS PRECEDENT (point fixe retarde d une iteration), ce que Vincent a
  ! explicitement autorise.
  logical, public :: pi_conv_stab = .false.
  ! pi_pivp : flux d energie = Pi_c (v_p.n) + relaxation nodale, avec Pi_c la
  ! pression de sous-face CLASSIQUE et p_p dans le flux de QUANTITE DE MOUVEMENT.
  !   (pv)*_l = Pi_l (v_p.n) - mu_l (e_p - e_l)
  !   (pv)*_r = Pi_r (v_p.n) + mu_r (e_p - e_r)
  ! Le telescopage somme A[(pv)*_r - (pv)*_l] = 0 donne
  !   e_p somme A(mu_l+mu_r) = somme A[mu_l e_l + mu_r e_r] - v_p.somme A(Pi_r-Pi_l) n
  ! et le dernier terme est NUL par l equation meme de la vitesse nodale. Donc
  !   **e_p = somme A[mu_l e_l + mu_r e_r] / somme A(mu_l + mu_r)**,
  ! une COMBINAISON CONVEXE des energies de maille, bornee quel que soit mu.
  ! La singularite en mu de la version `p_c v_c` -- ou le travail de pression
  ! n avait aucun facteur mu, d ou e_p -> inf et l effondrement du maillage --
  ! DISPARAIT. Le taux convectif redevient donc utilisable.
  ! pi_pp_ppus : LE flux d energie invariant par changement de referentiel.
  !   (pv)*_c = p_p (u*_c . n) -/+ mu_c (eps_p - eps_c)
  ! avec u*_c.n = v_c + ds_c la vitesse etoile de SON PROPRE cote
  ! (ds_l = (p_l-p_p)/lam_l, ds_r = (p_p-p_r)/lam_r), celle que donne la ligne de
  ! QUANTITE DE MOUVEMENT du saut, et eps l energie INTERNE.
  ! Mesure (jouet a 4 mailles, boost de 1, 5 et 20) : l ecart sur eps*_l vaut
  ! 4e-16 pour cette forme contre 1,4e-1 / 6,8e-1 / 2,7e+0 pour pivp, ep, ppvn et
  ! ppus en energie TOTALE -- lineaire en la vitesse du repere, et eps* devient
  ! NEGATIF a grande vitesse depuis un etat au repos parfaitement sain.
  ! Les TROIS conditions sont necessaires ensemble :
  !   (a) la pression du travail = p_p, la meme que celle de la QdM ;
  !   (b) la vitesse du travail = u*_c de son cote, pas v_p.n ;
  !   (c) la quantite relaxee = eps, pas e (une energie TOTALE n est pas invariante).
  ! Conservation EXACTE et gratuite : le travail p_p (u*_c.n) telescope seul, car
  ! somme A [(v_r+ds_r) - (v_l+ds_l)] = 0 EST la definition de p_p. Donc
  !   eps_p = somme A (mu_l eps_l + mu_r eps_r) / somme A (mu_l + mu_r),
  ! combinaison convexe des energies INTERNES.
  ! vp_gcorr : REDRESSEMENT GEOMETRIQUE de la reponse en pression de la vitesse
  ! nodale. Le systeme nodal donne
  !     v_p = Mp^-1 somme A(lam_l v_l + lam_r v_r) n  +  Mp^-1 Sp,
  ! avec Sp = somme_f A (p_l - p_r) n. Pour un champ de pression LINEAIRE
  ! p = p_p + g.(x - x_p) on a Sp = -Gp g avec le tenseur de moments
  !     Gp = somme_f A n (x) (x_r - x_l),
  ! donc la reponse vaut -Mp^-1 Gp g : elle n'est correcte que si Mp^-1 Gp est
  ! ISOTROPE. Mesure sur un maillage triangulaire structure : anisotropie 2,414
  ! (= 1 + racine(2)) meme sans perturbation, et jusqu'a 11 a 30 % de
  ! perturbation ; un gradient de pression produit alors une vitesse nodale
  ! deviee de ~21 degres et de module variable d'un facteur 2.
  ! Correctif : remplacer Mp^-1 Sp par s_p Gp^-1 Sp avec s_p = tr(Mp^-1 Gp)/d,
  ! qui rend la reponse exactement parallele a g et de module independant de sa
  ! direction. La MOYENNE DE VITESSE est inchangee. Gp est PUREMENT GEOMETRIQUE
  ! (aucun lambda) : y mettre les lambda reintroduit une anisotropie (mesure :
  ! 2,52 degres au lieu de 0 sur 6 triangles reguliers a lambda variable).
  ! Redonne le schema actuel au bit pres partout ou Mp^-1 Gp est deja isotrope
  ! (1D, cartesien, polygones reguliers a lambda uniforme).
  logical, public :: vp_gcorr = .false.
  integer(kind=ENTIER), public :: gcorr_n = 0, gcorr_nfall = 0
  logical, public :: pi_pp_ppus = .false.
  logical, public :: pi_pivp = .false.
  ! pi_tang : DISSIPATION SUR LA VITESSE TANGENTIELLE. Dans toutes les
  ! fermetures precedentes le tangentiel est INCHANGE -- la relation de quantite
  ! de mouvement -lam_l(u*_l - u_l) = (p_p - p_l) n est purement normale -- donc
  ! |u*_s|^2 = (u*_n)^2 + |u_t,s|^2 avec le tangentiel intact.
  ! En imposant un tangentiel etoile COMMUN, la convexite donne
  ! |u*_t|^2 <= <|u_t,s|^2> : l energie cinetique de cisaillement DIMINUE, donc
  ! l energie interne AUGMENTE, donc la marge eps* = e_p - |u*|^2/2 s OUVRE.
  ! C est exactement ce qui manque aux ecoulements FROIDS ET EN MOUVEMENT, ou la
  ! marge vaut ~eps_0 (Noh : 11 % des sous-faces en violation).
  ! Estimation sur Noh : pour h/r >= 1e-2 le gain vaut 5 a 500 fois la marge
  ! manquante, donc ca peut fermer le trou.
  ! CONSERVATION GRATUITE : imposer nu_l(u*_t - u_t,l) + nu_r(u*_t - u_t,r) = 0
  ! rend le flux tangentiel UNIVALUE, donc il telescope exactement, et definit
  !     u*_t = (nu_l u_t,l + nu_r u_t,r)/(nu_l + nu_r).
  ! Deux taux :
  !   tang_mode 1 : nu_s = lambda_s (ACOUSTIQUE). C est le modified_three_wave de
  !     la these : positivite plus facile, mais la dissipation de shear est
  !     active partout et POURRIT LE FLUX DE CHALEUR dans la couche limite.
  !   tang_mode 2 : nu_s = rho_s |u_t,r - u_t,l| (CONVECTIF). La dissipation vaut
  !     alors rho|du_t|^3 au lieu de lambda|du_t|^2, donc elle s eteint bien plus
  !     vite dans un cisaillement lisse. Le |du_t| se simplifie dans la moyenne,
  !     qui vaut (rho_l u_t,l + rho_r u_t,r)/(rho_l+rho_r) : pas de 0/0.
  ! pi_bord_psw : traitement de BORD unilateral, porte depuis le correctif `psw`
  ! valide cote Euler. A une sous-face de BORD il n y a pas de second etat, donc
  ! le glissement (p_l - p_p)/lam_l impose par la pression nodale n est borne
  ! par rien. Mesure sur Noh : au pas 2, eps_min passe de +2,51e-6 (classic) a
  ! -2,24e-4 avec 768 NaN, et les mailles fautives sont les BORDS et les COINS
  ! (1-4 et 342-400 sur 400), pas l interieur -- une extraction valant 90 fois
  ! l energie interne disponible. rho_min y est IDENTIQUE a 16 chiffres entre
  ! les deux schemas, donc seule l energie differe.
  ! Correctif : sous-faces de bord hors de la moyenne de p_p, et flux
  ! unilateral classique (Pi n, Pi v_p.n) au bord.
  logical, public :: pi_bord_psw = .false.
  ! convex_limiting : limitation de l ordre 2 par positivite de chaque SOUS-CELLULE,
  ! dont la maille est une combinaison convexe. .false. laisse l ordre 2 pur.
  logical, public :: convex_limiting = .false.
  integer(kind=ENTIER), public :: tang_mode = 0
  integer(kind=ENTIER), public :: tang_n = 0, tang_nact = 0
  ! Plancher PUREMENT NUMERIQUE, en fraction de rho a, pour que M_p ne soit pas
  ! identiquement nulle dans un ecoulement exactement uniforme ou v_p = v_c.
  ! Ce n est pas un parametre physique : on COMPTE les sous-faces ou il mord.
  real(kind=DOUBLE), public :: conv_floor = 1.0e-10_DOUBLE
  integer(kind=ENTIER), public :: conv_n = 0, conv_nfloor = 0
  ! Diagnostic PASSIF de la MARGE DE POSITIVITE de l energie interne
  ! intermediaire. Avec mu = lambda le facteur de relaxation vaut 1, donc
  ! e*_s = e_p des DEUX cotes et
  !      eps*_s = e_p - |u*_s|^2/2,
  ! si bien que la positivite equivaut EXACTEMENT a e_p > |u*_s|^2/2 sur chaque
  ! sous-face. e_p etant fixe par la conservation, la marge
  !      r = (e_p - |u*_s|^2/2) / max(|e_p|, tiny)
  ! se mesure directement. Cases : r <= 0 (VIOLATION), < 1e-3, < 1e-2, < 1e-1, >= 1e-1.
  integer(kind=ENTIER), public :: pos_hist(5) = 0
  integer(kind=ENTIER), public :: pos_tot = 0
  real(kind=DOUBLE), public :: pos_worst = 1.0e300_DOUBLE
  ! Diagnostic PASSIF : de combien la fermeture par MOINDRES CARRES deplace-t-elle
  ! la pression nodale par rapport a la condition de MOYENNE ? Si l ecart relatif
  ! est negligeable, ce n est pas l idee qui echoue mais l amplitude du levier.
  ! Cases de |p_mc - p_moy| / max(|p_moy|, tiny) : <1e-6, <1e-4, <1e-2, <1e-1, >=1e-1.
  integer(kind=ENTIER), public :: ppmc_hist(5) = 0
  integer(kind=ENTIER), public :: ppmc_tot = 0
  ! Defaut de telescopage INTERIEUR : somme A (flux_l(5) - flux_r(5)). Pour tous
  ! les autres schemas les deux expressions sont identiques et ceci vaut 0.
  real(kind=DOUBLE), public :: vpps_tel = 0.0_DOUBLE
  real(kind=DOUBLE), public :: vpps_tel_abs = 0.0_DOUBLE
  ! Diagnostic du terme ci-dessus. On sort un COMPTE et une somme, jamais le seul
  ! extremum : un max mesure la sous-population extreme et pas la population.
  ! psvp_ref cumule la dissipation de classic somme A lam delta^2, pour que le
  ! rapport ait un sens plutot que la valeur brute.
  ! --- lambda impose par la positivite, annexe app:psvppos ---------------------
  ! eps*_c >= 0 se reorganise en un TRINOME en lambda_c :
  !     A lam^2 - B lam - C >= 0,  A = eps_c + delta u_p,  B = p_c u_p - p_p u_c,
  !     C = (p_p - p_c)^2 / 2 >= 0
  ! d'ou, quand A > 0, lam >= (B + sqrt(B^2 + 4AC))/(2A). Et A >= eps_c - u_c^2/4,
  ! donc une solution EXISTE partout ou Ma_n^2 < 4/(gamma(gamma-1)) = 7.14, soit
  ! Ma_n < 2.67 a gamma = 1.4 -- au-dela AUCUN lambda ne convient.
  !
  ! ROUTE RETENUE : monter lambda PARTOUT, systeme nodal compris. La conservation
  ! exacte de l'energie tient a somme_c A Pi_c n_c = 0, qui n'est vrai que pour le
  ! Pi_c construit avec le lambda DU SYSTEME qui a produit v_p ; ne le monter que
  ! dans le flux casserait precisement la propriete pour laquelle ce schema existe.
  ! Le prix est un point fixe -- monter lambda change v_p donc delta donc lam_req --
  ! donc psvp_lam_iter passes et AUCUNE garantie : on COMPTE les violations
  ! restantes au lieu de les proclamer absentes.
  integer(kind=ENTIER), public :: psvp_lam_iter = 0
  real(kind=DOUBLE), public :: psvp_lam_marge = 1.0_DOUBLE
  real(kind=DOUBLE), dimension(:), allocatable :: psvp_lo_l, psvp_lo_r
  logical :: psvp_ovr = .false.
  ! eps* < 0 restants, et sous-faces HORS de la fenetre de faisabilite (A <= 0),
  ! ou aucun lambda n'aurait pu convenir : c'est ce second compte qui dit si
  ! l'obstruction est la fermeture ou le cas.
  integer(kind=ENTIER), public :: psvp_n_epsneg = 0, psvp_n_infeas = 0
  integer(kind=ENTIER), public :: psvp_n_sf = 0, psvp_n_neg = 0
  real(kind=DOUBLE), public :: psvp_sum = 0.0_DOUBLE, psvp_sum_neg = 0.0_DOUBLE
  real(kind=DOUBLE), public :: psvp_ref = 0.0_DOUBLE, psvp_worst = 0.0_DOUBLE
  ! Diagnostic S : production TOTALE ponderee d energie interne de la fermeture
  ! (pv)*_c = p_p v*_c, S = 1/2[(1/lam_l+1/lam_r) p_p^2 - p_l^2/lam_l - p_r^2/lam_r].
  ! S ne depend PAS du partage de l energie entre les deux cotes : une case
  ! NEGATIVE signifie p_p sous la moyenne QUADRATIQUE des p ponderee en 1/lam,
  ! donc une fermeture qui PREND de l energie interne sur cette sous-face.
  ! Comptes par case, pas de max : le max ne mesurerait que l extreme.
  ! Cases : < -1e-1, < -1e-2, < -1e-3, < -1e-6, < 0, >= 0.
  integer(kind=ENTIER), public :: sdiag_n(6) = 0
  integer(kind=ENTIER), public :: sdiag_tot = 0
  ! 'vp_pp_pos_sv' : la variante UNIVALUEE de l'ecretage (construction de Vincent).
  ! L'ecretage existant ci-dessus prend min(t_f, cap_l) a gauche et max(t_f, flo_r)
  ! a droite, donc DEUX valeurs par sous-face, donc une perte de conservation qu'il
  ! faut ensuite redistribuer. Ici on part de l'energie interne mise a jour, qui est
  ! sans ambiguite meme en l'absence d'etat etoile unique (p_p ne satisfait aucune
  ! relation de saut avec une maille donnee, donc le volume mis a jour par v_p.n et
  ! la quantite de mouvement mise a jour par p_p impliquent des vitesses etoiles
  ! differentes) :
  !     m_c deps_c/dt = - somme_f A (pv)*_f + v_c . somme_f A p_p n
  ! Avec le pas acoustique dt <= m_c/somme(A z_c) et la part A z_c eps_c par
  ! sous-face, n sortante de c et u_k = v_k.n :
  !     Phi_haut = + z_c eps_c + p_p u_c        (cote c)
  !     Phi_bas  = - z_d eps_d + p_p u_d        (cote d)
  !     delta^2  = Phi_haut - Phi_bas = z_c eps_c + z_d eps_d + p_p (u_c - u_d)
  ! Le flux garde sa valeur naturelle p_p (v_p.n) si elle est dans l'intervalle, et
  ! vaut le MILIEU sinon -- une seule valeur dans les deux cas, calculee des deux
  ! cotes, donc la conservation est intacte. Le milieu systematique changerait le
  ! flux d'un O(Ma) relatif la ou il n'y a rien a corriger ; cette forme-ci est
  ! strictement non regressive au bas Mach, ou le flux naturel est deja admissible.
  logical, public :: pvp_pos_sv = .false.
  ! Variante de comparaison : TOUJOURS le milieu, meme quand le flux naturel est
  ! admissible. Sert a mesurer ce que coute le milieu systematique au bas Mach.
  logical, public :: pvp_pos_sv_always = .false.
  logical, public :: pvp_pos_sv_check = .false.
  ! Valeur de BASE du flux d'energie, avant ecretage. Les quatre combinaisons du
  ! plan 2x2 : chaque facteur du produit est pris nodal ou 1D de sous-face.
  !   'V0'  p_p    * (v_p.n)      'V1'  p_bar_f * (v_p.n)
  !   'V3'  p_p    * v_bar_f      'V2'  p_bar_f * v_bar_f
  ! 'V4' : le flux d'energie 1D de Godunov lui-meme, qui n'est PAS p_bar*v_bar --
  ! il porte en plus le saut -(E_r - E_l), terme dissipatif du solveur. C'est la
  ! valeur que le probleme de Riemann de sous-face donne pour le flux d'energie,
  ! celle que vp_pp_pvp projetait ensuite sur un vecteur nodal ; l'employer
  ! directement est la version sans troisieme inconnue.
  ! Les quatre sont univaluees par sous-face, donc toutes conservatives sans
  ! bilan nodal. Les BORNES ne changent pas : elles derivent de la fermeture en
  ! quantite de mouvement, qui reste p_p dans les quatre cas.
  character(len=2), public :: pos_base_mode = 'V0'
  ! Mode PUR : la cible est employee telle quelle comme flux univalue, sans
  ! ecretage ni redistribution. Les cibles etant univaluees par sous-face, elles
  ! sont conservatives PAR CONSTRUCTION ; le limiteur ne sert qu'a la positivite.
  ! Ce mode isole donc la cible du limiteur, qu'aucun de nos tests ne separait.
  logical, public :: pos_no_limit = .false.
  ! Famille W, autour de V4. Avec lambda_l = lambda_r = lambda, le flux de Godunov
  ! se reduit a (p_l v_l + p_r v_r)/2 - (lambda/2)(E_r - E_l), soit un pv centre
  ! moins un saut d'energie pondere. Les deux degres de liberte sont donc ce qu'on
  ! differencie et avec quel coefficient :
  !   'W1'  saut d'energie INTERNE eps au lieu de l'energie totale E
  !   'W2'  saut de E, coefficient alpha (pos_w_alpha)
  !   'W3'  aucun saut : pv centre seul (cas degenere, temoin de la famille)
  !   'W4'  pression centree (p_l+p_r)/2 fois la vitesse NODALE, saut de E
  real(kind=DOUBLE), public :: pos_w_alpha = 1.0_DOUBLE
  ! Compteurs. Les _step sont remis a zero a CHAQUE appel ; les cumulatifs ne le
  ! sont jamais et portent donc sur le run entier. La distinction n'est pas
  ! cosmetique : un compteur remis a zero a chaque pas ne dit rien de l'historique,
  ! et "zero au dernier appel" n'etablit pas "zero partout".
  ! --- schema 'vp_hat' : Pi_chapeau = theta*p_p + (1-theta)*p_c - lam_c*delta_pcf,
  ! LE MEME objet dans le flux de quantite de mouvement et dans celui d'energie.
  ! theta = 0 doit redonner vp au bit pres : c'est le controle de plomberie.
  ! theta = 1 est degenere -- a un noeud interieur somme(A n) = 0, donc tout
  ! scalaire constant au noeud sort du bilan nodal et le pilotage par la pression
  ! disparait du second membre ; il ne survit qu'attenue d'un facteur (1-theta).
  logical, public :: pi_hat = .false.
  ! 'vp_facebar' : Pi = p_c + beta (p_bar_f - p_c) - lambda_c delta, avec
  ! p_bar_f = (p_l + p_r)/2 une moyenne de FACE. La moyenne est 1D, la
  ! stabilisation reste multidimensionnelle par v_p, et surtout elle reste LOCALE
  ! au noeud : c'est ce qui permet au systeme nodal d'exister en d x d. Une
  ! stabilisation en (div v)_c, elle, couplerait tous les noeuds de la maille.
  ! p_bar_f est univaluee a la face, donc elle sort du bilan nodal exactement
  ! comme p_p : la matrice M_p est celle de vp, inchangee. beta = 0 redonne vp.
  logical, public :: pi_face = .false.
  ! 'vp_mp' : Pi_pcf = p_bar_f - mu_p (div v)_p, LE MEME dans les deux flux.
  ! v_p n'est PAS determine par ce flux, il vient du solveur en vitesse existant :
  ! il n'y a donc aucun systeme nodal a resoudre ici, contrairement a une
  ! stabilisation en (div v)_c qui coupleriat tous les noeuds d'une maille.
  ! (div v)_p est la divergence NODALE de acoustic_module : -(1/V_p) somme_c N_pc.v_c.
  ! Conservation de la quantite de mouvement : le terme en p_bar_f s'annule car une
  ! valeur de face est portee par deux sous-faces de normales opposees, et le terme
  ! en mu (div v)_p s'annule par somme_c N_pc = 0 -- identite que acoustic_module
  ! verifie explicitement ("Bad circ") et qui ne vaut QU'AUX NOEUDS INTERIEURS.
  ! Diagnostic INDEPENDANT du schema, actif sur vp inchange : rapport entre la
  ! dissipation qu'une stabilisation en divergence produirait, somme_c mu_c
  ! (div v)_c^2 V_c, et celle que vp produit reellement, somme_c lambda_c somme
  ! A delta^2. Si ce rapport est O(1) sur un ecoulement incompressible, les
  ! vitesses nodales ne reproduisent pas un champ a divergence nulle et
  ! l'argument bas Mach tombe -- avant d'avoir code quoi que ce soit.
  logical, public :: diag_divratio = .false.
  ! Impression unique au PREMIER appel : compare la divergence NODALE, construite
  ! sur les vitesses de MAILLE, a Vpoint_c/V_c construit sur les vitesses NODALES.
  ! Les deux ne sont pas deux variantes du meme objet : au repos la premiere est
  ! nulle par definition, la seconde est pilotee par v_p donc par le saut de
  ! pression. C'est ce qui place une forme du bon ou du mauvais cote du partage
  ! noeud-maille.
  ! --- 'mean_div' : Pi_f = p_bar^lambda_f - mu_f (div v)_f, univoque par sous-face
  logical, public :: pi_meandiv = .false.
  logical, public :: pi_fbd = .false.
  ! Les deux roles de lambda sont INDEPENDANTS et rien n'oblige a prendre le meme
  ! objet : dans `pbar_f` c'est la vitesse d'onde qui fixe la pression d'interface
  ! (il lui faut les trois termes, c'est ce qui repare le contact) ; dans le
  ! DENOMINATEUR de D_p c'est un poids de la divergence nodale, et c'est lui qui
  ! s'effondre quand 1/lambda chute. 'acou' = acoustique partout, 'sep' = trois
  ! termes dans la moyenne et acoustique au denominateur.
  character(len=4), public :: fbd_lam_mode = 'tw'
  ! 'vp_fbw' : la fermeture fbd avec une vitesse d'onde lambda_w DETERMINEE par la
  ! positivite au lieu d'etre choisie. Deux resultats la motivent : avec la seule
  ! partie centrale Pi = pbar_f l'intervalle est toujours ouvert sous la condition
  ! de Gallice lambda^2 >= p^2/(2 eps) ; et pbar_f est exactement la pression qui
  ! MAXIMISE delta^2, donc le terme nodal coute la marge quadratiquement :
  !     delta^2(pbar_f - D_p) = delta^2(pbar_f) - (D_p^2/2)(1/lam_c + 1/lam_d)
  ! Tant que lambda est le meme partout il s'elimine de la condition, parce que
  ! D_p porte somme A/lambda au denominateur donc D_p est proportionnel a lambda --
  ! c'est l'identite qui explique qu'un facteur 2 sur lambda n'ait rien change.
  ! FIGER lambda = rho a dans le SEUL denominateur de D_p casse ce lien : D_p
  ! devient une donnee et la condition redevient une vraie condition, a solution
  ! finie. Pas de point fixe : D_p ne depend pas de lambda_w, lambda_w ne depend
  ! pas de Pi, tout se calcule en une passe.
  logical, public :: pi_fbw = .false.
  ! 'vp_mdiv' : divergences MULTIPLES. La fermeture est une base corrigee par le
  ! taux de dilatation du SOUS-ELEMENT de coin :
  !     Pi_pcf = pbar_f - D_p - kappa_pc * delta_pc
  !     delta_pc = Vpoint_pc/|omega_pc|   (une frequence)
  !     kappa_pc = theta * lambda_c * h_pc ,   h_pc = |omega_pc|/A_pc
  ! Signe physique : un sous-element comprime (delta < 0) voit sa pression de
  ! sous-face MONTER.
  ! Ce que ca apporte sur tout ce qu'on a essaye : `delta_pc` est indexe par le
  ! COIN et non par le noeud, donc il varie d'une sous-face a l'autre au meme
  ! noeud et `somme A n = 0` ne l'annule PAS -- contrairement a `D_p`, qui est un
  ! scalaire par noeud et sort du bilan nodal.
  ! `delta_pc` est evalue sur les vitesses nodales du pas PRECEDENT, exactement
  ! comme `vp_kappa` le fait, et par la meme routine verifiee `kappa_subcell_div`.
  ! DIFFERENCE DE FOND AVEC `vp_kappa` : chez lui le retard etait inoffensif parce
  ! que kappa multipliait une somme de CARRES, donc un kappa imprecis ne changeait
  ! que la QUANTITE de dissipation. Ici le terme est ADDITIF, donc un retard peut
  ! changer son SIGNE et produire de l'anti-dissipation. D'ou `kap_n_prod_neg`,
  ! qui n'est pas optionnel.
  logical, public :: pi_mdiv = .false.
  real(kind=DOUBLE), public :: mdiv_theta = 1.0_DOUBLE
  ! delta_pc par sous-element, meme stockage a largeur fixe que kap_c pour
  ! traverser le meme echange MPI : une maille FANTOME a besoin de tous ses
  ! sommets pour la divergence de son sous-element, et certains sont hors du halo.
  real(kind=DOUBLE), public, dimension(:,:), allocatable :: mdv_c
  real(kind=DOUBLE), public, dimension(:), allocatable :: mdv_se
  integer(kind=ENTIER_D), public :: mdv_n = 0, mdv_n_cap = 0
  integer(kind=ENTIER_D), public :: mdv_n_empty = 0, mdv_n_baseout = 0
  integer(kind=ENTIER_D), public :: mdv_n_nocon = 0
  ! Compte les sous-faces ou la dilatation est NON NULLE. Sans ce compteur, une
  ! fermeture `vp_mdiv` entierement inerte (delta = 0 partout) se lisait comme
  ! un succes : `mdv_n` comptait ses visites et `mdv_n_cap` sortait a zero, ce
  ! qui ressemble a « le plafond ne mord jamais » alors que cela veut dire
  ! « rien ne s'est passe ». C'est l'ecart entre les deux qui est le controle.
  integer(kind=ENTIER_D), public :: mdv_n_dltnz = 0
  ! DISTRIBUTION du rapport kappa_voulu/kappa_max, par decade. `mdv_ratio_max`
  ! seul est un MAX : il decrit la sous-population extreme et ne dit pas si le
  ! plafond mord sur tout le maillage ou sur trois sous-faces. C'est l'erreur
  ! qui revient le plus souvent ici, donc on sort l'histogramme avec.
  ! Plage descendue a -16 : avec le bas a -2, repris de `fbw_hist`, la decade la
  ! plus basse est un ACCUMULATEUR DE SATURATION (le clamp y empile tout ce qui
  ! est plus petit), et elle contenait 92,9 % de la population -- donc le corps
  ! de la distribution n'etait pas resolu, seulement majore.
  integer(kind=ENTIER_D), public, dimension(-16:16) :: mdv_rhist = 0
  real(kind=DOUBLE), public :: mdv_kap_max = 0.0_DOUBLE
  real(kind=DOUBLE), public :: mdv_ratio_max = 0.0_DOUBLE
  logical, public :: diag_area = .false.
  logical, public :: diag_lamratio = .false.
  integer(kind=ENTIER_D), public :: lamr_n = 0, lamr_eq1 = 0
  integer(kind=ENTIER_D), public :: lamr_g2 = 0, lamr_g10 = 0, lamr_g1e3 = 0
  integer(kind=ENTIER_D), public, dimension(0:16) :: lamr_hist = 0
  real(kind=DOUBLE), public :: lamr_max = 0.0_DOUBLE
  logical, public :: diag_alpha = .false.
  real(kind=DOUBLE), public, dimension(:), allocatable :: alpha_sal
  real(kind=DOUBLE), public, dimension(:), allocatable :: area_a0
  real(kind=DOUBLE), public, dimension(:), allocatable :: forme_a0
  real(kind=DOUBLE), public :: fbw_marge = 1.3_DOUBLE
  ! Distribution de lambda_w / rho a, en bins log10 : un MAXIMUM ne dirait rien sur
  ! une population heterogene, et c'est ce chiffre qui decidera du pas de temps.
  integer(kind=ENTIER_D), public, dimension(-2:16) :: fbw_hist = 0
  integer(kind=ENTIER_D), public :: fbw_n = 0, fbw_n_eq1 = 0
  integer(kind=ENTIER_D), public :: fbw_n_g10 = 0, fbw_n_g1e3 = 0, fbw_n_g1e6 = 0
  integer(kind=ENTIER_D), public :: fbw_n_viol_ff = 0, fbw_n_viol_cf = 0
  integer(kind=ENTIER_D), public :: fbw_n_viol_cc = 0, fbw_n_sf = 0
  real(kind=DOUBLE), public :: fbw_d2_min = 1.0e300_DOUBLE
  integer(kind=ENTIER), public :: fbw_st_ff = 0, fbw_st_cf = 0, fbw_st_cc = 0
  real(kind=DOUBLE), public :: fbw_st_d2min = 1.0e300_DOUBLE, fbw_st_rmax = 0.0_DOUBLE
  ! Drapeau par MAILLE : porte-t-elle au moins une sous-face a delta2 <= 0 ?
  ! Sert a separer les deux populations quand on lit eps_min -- sans cette
  ! separation, "eps_min > 0 alors que 37654 sous-faces violent" est circulaire :
  ! on ne sait pas si les mailles concernees sont celles qui vont bien.
  integer(kind=ENTIER), public, dimension(:), allocatable :: fbw_cviol
  ! Denominateurs : sous-faces VISITEES par type. Un compte brut sur une population
  ! qui grandit (le souffle s'etend, donc les interfaces chaud-chaud se multiplient)
  ! ne se compare pas a lui-meme.
  integer(kind=ENTIER_D), public :: fbw_nv_ff = 0, fbw_nv_cf = 0, fbw_nv_cc = 0
  integer(kind=ENTIER_D), public :: fbd_n_sf = 0, fbd_n_node = 0
  integer(kind=ENTIER_D), public :: fbd_n_node_void = 0, fbd_n_zbad = 0
  integer(kind=ENTIER_D), public :: fbd_n_viol_ff = 0, fbd_n_viol_cf = 0
  integer(kind=ENTIER_D), public :: fbd_n_viol_cc = 0
  real(kind=DOUBLE), public :: fbd_d2_min = 1.0e300_DOUBLE
  ! Par PAS, pour savoir si |D_p| explose AVANT ou APRES le premier eps < 0 :
  ! si c'est apres, on lirait une consequence pour une cause.
  real(kind=DOUBLE), public :: fbd_st_dpmax = 0.0_DOUBLE
  real(kind=DOUBLE), public :: fbd_st_d2min = 1.0e300_DOUBLE
  integer(kind=ENTIER), public :: fbd_st_zbad = 0, fbd_st_void = 0
  integer(kind=ENTIER), public :: fbd_st_ff = 0, fbd_st_cf = 0, fbd_st_cc = 0
  ! Origine des lambda invalides. ATTENTION : `lambda <= 0` est FAUX pour un NaN,
  ! donc un lambda NaN passait mon garde-fou sans etre compte. On teste les deux
  ! separement, et on remonte l'etat (p, tau) qui l'a produit.
  integer(kind=ENTIER_D), public :: fbd_n_znan = 0, fbd_n_zzero = 0
  integer(kind=ENTIER_D), public :: fbd_zb_pnan = 0, fbd_zb_pneg = 0
  integer(kind=ENTIER_D), public :: fbd_zb_taunan = 0, fbd_zb_tauneg = 0
  integer(kind=ENTIER_D), public :: fbd_zb_autre = 0
  real(kind=DOUBLE), public :: fbd_dp_absmax = 0.0_DOUBLE
  real(kind=DOUBLE), public :: fbd_num_bnd = 0.0_DOUBLE
  real(kind=DOUBLE), public :: md_cvisc = 1.0_DOUBLE
  real(kind=DOUBLE), public :: md_sq_div = 0.0_DOUBLE, md_sq_delta = 0.0_DOUBLE
  real(kind=DOUBLE), public :: md_prod_min = 0.0_DOUBLE, md_res_max = 0.0_DOUBLE
  ! residu du bilan nodal DECOMPOSE par direction : sur un maillage 2D obtenu par
  ! extrusion les faces z sont unilaterales (idr < 0) alors que leurs noeuds ne
  ! sont PAS marques is_bound. Un residu purement en z serait donc l'artefact de
  ! l'extrusion et non une erreur d'assemblage ; un residu en x ou y le serait.
  real(kind=DOUBLE), dimension(3), public :: md_res_dir = 0.0_DOUBLE
  ! Statistiques de Pi_f au PREMIER pas : la moyenne ponderee par l'impedance
  ! p_bar^lam = (lam_d p_c + lam_c p_d)/(lam_c + lam_d) pondere p_c par l'impedance
  ! de l'AUTRE maille. Sous un saut de pression violent, lam_c >> lam_d, et elle
  ! tend vers la moyenne GEOMETRIQUE sqrt(p_c p_d), donc vers la pression BASSE.
  ! On mesure donc max p_c, max Pi_f et le nombre de Pi_f negatifs.
  logical, public :: md_first = .false., md_first_done = .false.
  real(kind=DOUBLE), public :: md_pmax = 0.0_DOUBLE, md_pimax = -1.0e300_DOUBLE
  real(kind=DOUBLE), public :: md_pimin = 1.0e300_DOUBLE
  real(kind=DOUBLE), public :: md_pimax_b = -1.0e300_DOUBLE
  ! lambda a deux ondes (TRRS de Toro) au lieu de l'acoustique rho*a : sous un
  ! saut fort lam_c/lam_d = sqrt(p_c/p_d) est enorme et p_bar^lam s'effondre sur
  ! la moyenne geometrique. Le TRRS fait voir le choc a l'impedance.
  logical, public :: md_trrs = .false.
  ! Desactive le switch de detente : mu_f actif meme quand (div v)_f > 0. VIOLE la
  ! condition 4 de Caramana -- la viscosite chaufferait alors un gaz en detente,
  ! ce qui est indefendable en general. Sert UNIQUEMENT a tester si c'est ce
  ! switch qui tue le schema au depot, ou ce dernier est en expansion violente et
  ! ou mu est donc eteint exactement la ou le premier pas est le plus raide.
  logical, public :: md_noswitch = .false.
  ! dt vaut 0 au PREMIER appel : lagrange_main l'initialise a zero et n'appelle
  ! compute_dt qu'APRES le rhs. Un diagnostic d'energie predite pris au premier
  ! appel vaut donc eps_new = eps_old et ne predit rien. On le prend au DEUXIEME.
  ! Trace PAR PAS : le schema fait une cinquantaine de pas propres puis meurt
  ! toujours au meme endroit, independamment de t_max et de c_visc. Ce qui
  ! s'accumule sur cinquante pas d'un lagrangien sur simplices, c'est la
  ! DEFORMATION du maillage. On suit donc le volume minimal, le signe des volumes,
  ! et on separe les NaN de rho de ceux de l'energie : densite d'abord = maille
  ! retournee, energie d'abord = probleme de flux.
  logical, public :: md_trace = .false.
  ! Suivi PAR MAILLE de Pi - p_c et de la production somme_f A (p_c - Pi) delta.
  ! Rempli par les DEUX schemas pour que la comparaison porte sur la meme
  ! quantite : pour vp, Pi - p_c vaut -lambda_c delta et la production est un
  ! carre ; pour mean_div elle n'a pas de signe.
  real(kind=DOUBLE), allocatable, public :: w_prod(:), w_dpi_min(:), w_dpi_max(:)
  logical, public :: md_watch = .false.
  ! Compte les sous-faces ou la garde de vide se declenche, par pas et sur le run.
  ! Sans ca, un sedov_tri qui passe ne dirait pas LEQUEL des deux schemas
  ! d'impedance a reellement tourne la ou ca se joue.
  integer(kind=ENTIER_D), public :: md_n_vac_step = 0, md_n_vac_tot = 0
  integer(kind=ENTIER_D), public :: md_n_sf_step = 0
  logical, public :: md_nan_done = .false., md_src_done = .false.
  integer(kind=ENTIER), public :: md_ncall = 0
  integer(kind=ENTIER_D), public :: md_n_eps_neg = 0
  real(kind=DOUBLE), public :: md_eps_min = 0.0_DOUBLE, md_w_at_min = 0.0_DOUBLE
  real(kind=DOUBLE), public :: md_eps_old_at_min = 0.0_DOUBLE
  integer(kind=ENTIER_D), public :: md_n_pi_neg = 0
  integer(kind=ENTIER_D), public :: md_n_cell_neg = 0, md_n_cell_tot = 0
  logical, public :: diag_firstprint = .false.
  logical, public :: dr_done = .false.
  real(kind=DOUBLE), public :: dr_max_divp = 0.0_DOUBLE, dr_max_vdot = 0.0_DOUBLE
  real(kind=DOUBLE), allocatable, public :: dr_vdot(:)
  real(kind=DOUBLE), public :: dr_sq_delta = 0.0_DOUBLE, dr_sq_div = 0.0_DOUBLE
  logical, public :: pi_mp = .false.
  real(kind=DOUBLE), public :: mp_cvisc = 1.0_DOUBLE
  ! diagnostics : la production totale par noeud est une somme de carres donc >= 0
  ! par identite, mais la positivite PAR MAILLE n'est pas garantie. On compte donc
  ! les mailles a production negative, c'est le resultat et pas un controle.
  real(kind=DOUBLE), public :: mp_prod_node = 0.0_DOUBLE
  real(kind=DOUBLE), public :: mp_prod_min_cell = 0.0_DOUBLE
  integer(kind=ENTIER_D), public :: mp_n_cell_neg = 0, mp_n_cell_tot = 0
  integer(kind=ENTIER_D), public :: mp_n_bnd_node = 0, mp_n_node_tot = 0
  real(kind=DOUBLE), public :: theta_pi = 0.0_DOUBLE
  ! theta_M = min(Ma_local, 1) applique a la part VISQUEUSE lam_c*delta : c'est
  ! la correction de type Chalons, le terme non centre amorti par le Mach local.
  ! theta_g (= theta_pi) et theta_M agissent sur deux termes DIFFERENTS : la
  ! mesure du balayage montre que theta_g atteint la cavite et pas le tourbillon,
  ! et que lam_c*delta est le terme qui tue le tourbillon en bas Mach.
  ! Le plancher dissipatif reste une somme de carres, seulement reduite :
  !   production = theta_M lam_c somme(A delta^2) + theta_g somme(A (p_c-p_p) delta)
  logical, public :: pi_mach_damp = .false.
  ! PLANCHER sur theta_M, indispensable : theta_M entre dans la MATRICE nodale
  ! M_p = somme A theta_M (lam_l+lam_r) n(x)n. Si theta_M s'annule sur toutes les
  ! sous-faces d'un noeud -- ce qui arrive des t=0 sur Sedov, ou le fluide est au
  ! repos donc Ma = 0 partout, et dans les zones au repos de Gresho -- la matrice
  ! est NULLE, et inverse_3_by_3 divise par un determinant nul sans garde.
  ! Mesure : les quatre runs a plancher nul ont rendu un champ entierement NaN.
  real(kind=DOUBLE), public :: theta_m_floor = 0.0_DOUBLE
  real(kind=DOUBLE), public :: pih_tm_min = 2.0_DOUBLE, pih_tm_max = -1.0_DOUBLE

  ! === 'vp_kappa' : vp avec un POIDS kappa par COIN sur l'impedance ==========
  !   Pi_pcf   = p_c - kappa_pc lambda_c delta_pcf      delta = (u_p - v_c).n_c
  !   kappa_pc = max(kap_floor, min(1, h_c |(div v)_pc| / a_c))
  ! ou (div v)_pc est la divergence de SOUS-MAILLE, fermee sur les vitesses
  ! NODALES et leurs interpolees -- JAMAIS sur v_c (voir kappa_subcell_div).
  !
  ! Positivite : la production vaut kappa_pc lambda_c somme(A delta^2), une somme
  ! de carres pour TOUT kappa >= 0, y compris variable par coin ET pris au pas
  ! precedent. C'est la difference structurelle avec toutes les formes en
  ! divergence essayees avant, ou la divergence etait DANS le flux : la retarder
  ! y transformait le carre en produit de deux pas de temps. Ici elle n'entre que
  ! comme multiplicateur scalaire positif, donc kappa est pris EXPLICITE.
  !
  ! PLANCHER, et il n'est pas cosmetique : M_p est proportionnelle a kappa alors
  ! que le terme de pression de R_p ne l'est pas, donc v_p gagne une composante
  ! en h |grad p| / kappa. Au bas Mach |grad p| = O(Ma^2) et kappa mesure
  ! O(h^2 Ma), ce qui porte cette composante a O(Ma/h) contre un O(Ma) physique.
  ! La condition pour que le terme de vitesse domine est kappa >~ h Ma : c'est
  ! l'echelle de kap_floor, pas une garde anti-NaN. Le symptome a surveiller est
  ! max|v_p| contre max|v_c|, pas l'apparition de NaN.
  ! Arret dur sur dt = NaN. Defaut .false. : la ligne d'alerte sort de toute
  ! facon, et passer a .true. change le comportement de TOUS les schemas, donc
  ! c'est une decision de Vincent. Declare ici et non dans le bloc de la
  ! pseudo-inverse, qui a ete retire : ce diagnostic est GENERAL.
  logical, public :: stop_on_nan_dt = .false.
  ! Pose par compute_dt quand une BORNE est NaN, et reduit sur tous les rangs. Necessaire
  ! parce qu un NaN local defait le min sans rendre dt lui-meme NaN sur les autres rangs.
  logical, public :: dt_nan_detected = .false.
  ! NaN dans la vitesse nodale d ORDRE 1, remis a zero a chaque appel de compute_rhs_lagrange.
  ! Contrairement a l ordre 2 il n y a aucun ordre inferieur ou se replier : c est une panne.
  integer(kind=ENTIER), public :: vp_o1_nan_n = 0, vp_o1_nan_node = -1
  logical, public :: pi_pif = .false.
  ! TEST D'ATTRIBUTION : prendre v_p du solveur 'vp' existant pour deplacer le
  ! maillage, et garder (Pi_f, U_f) de la construction par sous-face pour les
  ! flux. C'est licite parce que la conservation ne fait PAS intervenir v_p --
  ! mesure sur le code reel : energie a 0.0000e+00 exactement et quantite de
  ! mouvement a 6.3e-12 sur Gresho. Les deux objets sont decouples.
  !
  ! Ce n'est PAS un schema : le maillage bougerait avec le v_p de vp pendant que
  ! les flux utilisent U_f, donc la loi de conservation geometrique n'est pas
  ! satisfaite -- le dV/dt du maillage ne correspond pas a la divergence des
  ! flux. Ce desaccord existe deja avec l'ajustement actuel, qui n'est pas
  ! davantage compatible. C'est un test d'attribution et rien d'autre.
  logical, public :: pif_vp_classic = .false.
  ! IMPEDANCE DE CHOC. La single-valuedness de Pi est ce qui donne les deux
  ! conservations exactes, et elle coute un facteur ~5 sur l'entrainement du
  ! souffle : Pi est force a la moyenne acoustique (l_R p_L + l_L p_R)/(l_L+l_R),
  ! qui tend vers p_R quand lambda_L >> lambda_R -- l'interface prend la pression
  ! de la maille MOLLE. Mesure : sur un depot a 1e+05, la maille ambiante voisine
  ! recoit 20 % de la force de vp, soit un facteur 4.9, a comparer au facteur 4.3
  ! mesure sur |v_c| a temps physique egal.
  !
  ! Reparation a la Dukowicz : lambda_k = rho_k (a_k + Gamma |delta u_k|) avec
  ! Gamma = (gamma+1)/2, iteree -- c'est un point fixe local sur lambda, pas une
  ! resolution, et une seule iteration ne suffit pas.
  !
  ! Elle ne casse rien : lambda n'entre que par 1/lambda_c dans les poids, donc
  ! la structure de sous-maille, la single-valuedness et les deux conservations
  ! sont intactes. Et AU BAS MACH |delta u| -> 0 donc lambda -> rho a : la
  ! fermeture est strictement la meme dans le regime qui nous interesse et ne
  ! differe que la ou l'acoustique est fausse.
  logical, public :: pif_duk = .false.
  integer(kind=ENTIER), public :: pif_duk_it = 5
  logical, public :: pif_check_done = .false.
  logical, public :: pi_kappa = .false.
  real(kind=DOUBLE), public :: kap_floor = 1.0e-3_DOUBLE
  ! kappa = 1 partout : doit redonner 'classic' BIT A BIT, pas "a 5e-13". La
  ! formule est litteralement celle de vp, sans rearrangement arithmetique.
  logical, public :: kap_one = .false.
  ! kappa IMPOSE a une constante (> 0 pour activer), donc calcule depuis RIEN.
  ! Sert a sortir le senseur de la boucle : deux relations sont en jeu,
  !   (a) kappa = |v_p|/a       le senseur, qui lit la vibration nodale
  !   (b) |v_p| = F(kappa)      l'amortissement nodal, decroissant en kappa
  ! et la question est de savoir si kappa se verrouille sur leur intersection,
  ! c'est-a-dire si le senseur finit par mesurer sa propre vibration au lieu de
  ! l'ecoulement. En figeant kappa on sort de (a) et on mesure (b) SEULE ; il
  ! reste a verifier si l'intersection de F avec la droite |v_p| = kappa*a tombe
  ! sur le kappa_moyen observe en mode senseur.
  real(kind=DOUBLE), public :: kap_fixed = -1.0_DOUBLE
  ! kappa est porte par un tableau de MAILLE a largeur fixe, indice par la
  ! position du sommet dans elem%vert, pour pouvoir passer par l'echange MPI
  ! generique : une maille fantome a besoin de TOUS ses sommets pour sa
  ! divergence de sous-maille, et certains sont hors du halo local.
  integer(kind=ENTIER), parameter, public :: KAP_NV = 8
  real(kind=DOUBLE), allocatable, public :: kap_c(:,:)
  real(kind=DOUBLE), allocatable, public :: kap_se(:)
  logical, public :: kap_started = .false.
  real(kind=DOUBLE), public :: kap_min = 2.0_DOUBLE, kap_max = -1.0_DOUBLE
  real(kind=DOUBLE), public :: kap_sum = 0.0_DOUBLE, kap_divmax = 0.0_DOUBLE
  integer(kind=ENTIER_D), public :: kap_n = 0, kap_n_floor = 0
  ! Controle d'EXACTITUDE de l'operateur de sous-maille, une fois par run :
  ! champ nodal constant (divergence nulle) puis champ lineaire de divergence
  ! connue. C'est exactement le defaut qui avait fait sortir une erreur O(1) a
  ! l'essai precedent (v_c mis sur les faces interieures), donc il est MESURE.
  ! Production par TERME : chaque contribution vaut A kappa lambda delta^2 et est
  ! donc positive SEULE, par identite. Toute valeur negative signale un defaut
  ! d'implementation et non une propriete du schema -- c'est le controle le plus
  ! serre dont on dispose, plus serre que la positivite par maille.
  integer(kind=ENTIER_D), public :: kap_n_prod = 0, kap_n_prod_neg = 0
  real(kind=DOUBLE), public :: kap_prod_min = 0.0_DOUBLE, kap_prod_tot = 0.0_DOUBLE
  ! kappa RESTREINT AU VOISINAGE DU FRONT. La distribution globale ne repond pas
  ! a la question qui compte : kappa etant explicite, une maille devant un choc a
  ! un v_p du pas precedent ou le choc n'etait pas encore arrive, donc |div v|
  ! petit, donc kappa au plancher -- et elle recoit le choc SOUS-DISSIPEE pendant
  ! un pas. C'est exactement le mode de defaillance de la famille W, qui mourait
  ! de n'avoir rien au moment critique. On sort donc le kappa minimal parmi les
  ! mailles qui portent un vrai saut de pression, pas seulement le minimum global.
  real(kind=DOUBLE), public :: kap_front_jump = 0.1_DOUBLE
  real(kind=DOUBLE), public :: kap_front_min = 2.0_DOUBLE
  real(kind=DOUBLE), public :: kap_front_jmax = 0.0_DOUBLE
  integer(kind=ENTIER_D), public :: kap_n_front = 0, kap_n_front_floor = 0
  ! bilan nodal |somme A Pi n| aux noeuds INTERIEURS seulement
  real(kind=DOUBLE), public :: kap_res_max = 0.0_DOUBLE
  ! --- plancher ISENTROPIQUE, pour qualifier une chute d'energie interne.
  ! La production positive n'interdit PAS a eps de descendre : le terme
  ! reversible -p_c Vpoint_c est negatif en detente et vide la maille
  ! legitimement quand elle s'etend. Elle interdit seulement qu'il descende
  ! SOUS le plancher adiabatique
  !     plancher_c(t) = eps_c(0) * (V_c(0)/V_c(t))^(gamma-1)
  ! et comme la masse est constante en lagrangien, V_c est proportionnel a
  ! tau_c = sol(1,:) : le rapport ne demande que tau et eps initiaux.
  !
  ! rapport >= 1 : eps est AU plancher, la maille se detend, le schema fait ce
  !   que la preuve annonce, et une panne est alors un probleme de PAS DE TEMPS ;
  ! rapport <  1 : le plancher est VIOLE alors que chaque terme de production est
  !   positif -- la borne est semi-discrete et le pas d'Euler explicite ne la
  !   reproduit que sous condition de dt. Defaut reel, mais d'une autre nature.
  ! --- PLANCHER LOCAL sur kappa, par NOEUD ------------------------------------
  ! L'amplification parasite vient de ce que M_p est proportionnelle a kappa
  ! alors que le terme de pression de R_p ne l'est pas. Or R_p est exactement la
  ! somme des deux :
  !   R_p = somme_f A [kappa_l lam_l u_l + kappa_r lam_r u_r] n   (terme VITESSE)
  !       - somme_f A [p_r - p_l] n                               (terme PRESSION)
  ! donc la condition "le terme de vitesse domine" s'ecrit sans rien de neuf :
  !   kappa_min,p = min(1, C ||somme A (p_r-p_l) n|| / ||somme A (lam_l u_l + lam_r u_r) n||)
  ! Les deux sommes sont prises SANS kappa, ce qui evite toute circularite.
  !
  ! Comportement dans les trois regimes :
  !   zone au repos : terme vitesse -> 0, terme pression -> bruit, kappa -> 1,
  !     et c'est GRATUIT puisque delta = 0 y rend la dissipation nulle quel que
  !     soit kappa ;
  !   tourbillon    : vitesse ~ h lam |u|, pression ~ h^2 Ma^2, rapport ~ h Ma^2,
  !     donc minuscule et kappa reste celui du senseur ;
  !   choc          : les deux sont grands, rapport O(1), kappa -> 1, ce que le
  !     senseur donnait deja.
  !
  ! Calcule au noeud et applique a tous ses coins, donc SANS toucher a
  ! compute_rhs_around_node -- le controle kappa = 1 bit a bit reste intact.
  ! Coherent en MPI sans echange : le noeud d'une maille possedee a son anneau
  ! complet en local sur les deux rangs, c'est ce qui rend deja le solveur nodal
  ! identique de part et d'autre.
  logical, public :: kap_loc = .false.
  real(kind=DOUBLE), public :: kap_loc_c = 1.0_DOUBLE
  ! 'node' : plancher calcule sur les RESULTANTES vectorielles du noeud (forme
  !   d'origine). Faiblesse : ce sont des vecteurs, donc ils se COMPENSENT. Un
  !   noeud entoure de sauts de pression forts mais opposes -- exactement un
  !   noeud au coeur d'un choc, ou entre deux mailles symetriques en compression
  !   -- a une resultante faible alors que chaque terme est grand. Le plancher
  !   n'y mord pas, et c'est precisement la qu'il devrait.
  ! 'face' : meme critere PAR SOUS-FACE, en valeurs absolues scalaires, donc
  !   sans compensation possible, et pose PAR COIN (par sous-element) et non par
  !   noeud. Le coin prend le max sur ses sous-faces.
  character(len=8), public :: kap_loc_mode = 'node'
  real(kind=DOUBLE), allocatable, public :: kap_se_floor(:)
  real(kind=DOUBLE), allocatable, public :: kap_cell_floor(:)
  ! Laquelle des deux bornes de compute_dt fixe dt, et les deux valeurs. La
  ! borne acoustique vaut m_c/somme(A lambda) ; la borne en VOLUME vaut
  ! CV V_c/somme(A |u_p.n|) et depend donc des VITESSES NODALES. Si c'est elle
  ! qui mord avec un volume sain, le mode nodal entre directement dans le pas de
  ! temps -- un mecanisme distinct de l'ecrasement de maille.
  logical, public :: diag_dtb = .false.
  real(kind=DOUBLE), public :: dtb_acou = 0.0_DOUBLE, dtb_vol = 0.0_DOUBLE
  integer(kind=ENTIER), public :: dtb_cell_acou = -1, dtb_cell_vol = -1
  real(kind=DOUBLE), public :: dtb_vol_V = 0.0_DOUBLE, dtb_vol_sumav = 0.0_DOUBLE
  integer(kind=ENTIER), public :: dtb_n_nan_acou = 0, dtb_n_nan_vol = 0
  integer(kind=ENTIER), public :: dtb_i_first_nan = -1
  integer(kind=ENTIER_D), public :: dtb_nan_tot_acou = 0, dtb_nan_tot_vol = 0
  integer(kind=ENTIER_D), public :: dtb_nan_ghost = 0, dtb_nan_reel = 0
  integer(kind=ENTIER_D), public :: dtb_nan_src_sol = 0, dtb_nan_src_vp = 0
  integer(kind=ENTIER_D), public :: dtb_nan_src_vol = 0
  real(kind=DOUBLE), public :: dtb_acou_m = 0.0_DOUBLE, dtb_acou_sal = 0.0_DOUBLE
  real(kind=DOUBLE), public :: dtb_acou_lam = 0.0_DOUBLE, dtb_acou_area = 0.0_DOUBLE

  ! --- Diagnostic de FORCE NETTE sur une maille choisie.
  ! La force sur une maille est -somme A Pi n_sortante. Comme somme A n = 0 par
  ! fermeture des normales, toute partie UNIFORME de Pi ne produit aucune force :
  ! seule la VARIATION de Pi autour de la maille compte. On sort donc la force
  ! nette VECTORIELLE (pas une somme de |.|, qui masquerait justement la
  ! variation) et les moments de Pi ponderes par l'aire.
  logical, public :: diag_force = .false.
  integer(kind=ENTIER), public, dimension(4) :: fdiag_cells = -1
  real(kind=DOUBLE), public, dimension(3, 4) :: fdiag_F = 0.0_DOUBLE
  real(kind=DOUBLE), public, dimension(4) :: fdiag_sA = 0.0_DOUBLE, &
    fdiag_sAP = 0.0_DOUBLE, fdiag_sAP2 = 0.0_DOUBLE, &
    fdiag_Pmin = 1.0e300_DOUBLE, fdiag_Pmax = -1.0e300_DOUBLE
  integer(kind=ENTIER), public, dimension(4) :: fdiag_n = 0
  ! Liste brute des (Pi, maille d'en face) par sous-face, pour pouvoir verifier
  ! si une sous-face qui FAIT FACE au depot porte bien la pression du depot.
  ! Sans ca, 'Pi vaut l'ambiant partout autour de 14812' ne distingue pas un
  ! defaut de fermeture d'une maille qui ne touche simplement pas le depot.
  real(kind=DOUBLE), public, dimension(32, 4) :: fdiag_lPi = 0.0_DOUBLE
  integer(kind=ENTIER), public, dimension(32, 4) :: fdiag_lNb = -1
  real(kind=DOUBLE), allocatable, public :: kap_node_floor(:)
  real(kind=DOUBLE), public :: kfl_min = 2.0_DOUBLE, kfl_max = -1.0_DOUBLE
  real(kind=DOUBLE), public :: kfl_sum = 0.0_DOUBLE
  integer(kind=ENTIER_D), public :: kfl_n = 0, kfl_n_one = 0
  ! kappa EFFECTIF = max(senseur, plancher du noeud). Les statistiques KAP sont
  ! accumulees dans compute_kappa, donc AVANT l'application du plancher, et ne
  ! decrivent que le senseur. Sans cette seconde serie, le kappa publie n'est pas
  ! celui que le schema utilise -- et le gain mesure n'est pas interpretable.
  real(kind=DOUBLE), public :: keff_min = 2.0_DOUBLE, keff_max = -1.0_DOUBLE
  real(kind=DOUBLE), public :: keff_sum = 0.0_DOUBLE
  integer(kind=ENTIER_D), public :: keff_n = 0, keff_n_raised = 0
  ! Les termes a production negative, detailles. La production vaut
  ! A kappa lambda delta^2, positive par identite ; mais le diagnostic la calcule
  ! comme A (p_c - Pi) delta, c'est-a-dire par une DIFFERENCE DE DEUX PRESSIONS.
  ! Quand kappa*lambda*delta est negligeable devant p_c, cette difference subit
  ! une annulation catastrophique et peut sortir avec le mauvais SIGNE -- un
  ! artefact de mesure, pas un defaut du schema. On imprime donc les deux formes
  ! cote a cote pour les premiers cas et on tranche sur les chiffres.
  integer(kind=ENTIER), public :: kap_neg_print = 0
  ! === 'vp_div' : la divergence de sous-maille ENTRE DANS LA FERMETURE ========
  !   Pi_pcf = p_c - h_c (div v)_pc
  ! Base UNILATERALE p_c -- c'est elle qui decide, la mesure comparative ayant
  ! etabli que rho_min est identique a 17 chiffres entre vp et mean_div au pas 2,
  ! donc que toute la difference venait de l'energie et de la base moyennee.
  ! Base unilaterale PLUS divergence consistante est une combinaison jamais
  ! essayee.
  !
  ! COUPLAGE. (div v)_pc depend de tous les noeuds de la maille par les moyennes
  ! au centroide, donc le systeme nodal couplerait tout un anneau. On l'evite en
  ! exploitant que (div v)_pc est LINEAIRE en les vitesses nodales :
  !     (div v)_pc = ( b_pc . v_p + reste_pc ) / |omega_pc|
  ! avec b_pc purement GEOMETRIQUE. En gardant b_pc.v_p implicite et le reste au
  ! pas precedent, la matrice nodale redevient locale :
  !     M_p = somme_f A [ h_c b_pc/|om_pc| - h_d b_pd/|om_pd| ] (x) n     (d x d)
  ! non symetrique, mais locale et generiquement inversible.
  !
  ! PRODUCTION. Elle n'est PAS un carre. En ecrivant la divergence consistante
  ! comme la naive plus une correction, (div v)_pc = G_pc/|om_pc| + Delta_pc
  ! avec G_pc = (v_p - v_c).N_pc, il vient
  !   m deps/dt = -p_c Vpoint_c + (h/|om|) somme_p G_pc^2 + h somme_p Delta_pc G_pc
  ! soit un CARRE plus un terme indefini. Mesure du pair sur maillage 40x40
  ! perturbe : exactement nulle sur cisaillement pur et rotation rigide (zero
  ! coin negatif), inconditionnellement dissipative en compression uniforme
  ! (zero coin negatif), et legerement ANTI-dissipative sur un tourbillon a
  ! divergence nulle -- mais cette part converge a l'ordre 3 en h et vaut cinq
  ! ordres sous la dissipation de compression des n = 80.
  !
  ! Le compteur de coins a production negative est donc NON INFORMATIF ici : il
  ! vaut ~50 % par symetrie du bruit. C'est la SOMME PAR MAILLE et son signe qui
  ! sont le resultat, puisque c'est elle qui vide ou remplit la maille.
  logical, public :: pi_divdirect = .false.
  real(kind=DOUBLE), public :: divd_coef = 1.0_DOUBLE
  real(kind=DOUBLE), allocatable, public :: divd_b(:,:)      ! (3, n_sub_elems)
  real(kind=DOUBLE), allocatable, public :: divd_rest(:)     ! (n_sub_elems)
  real(kind=DOUBLE), allocatable, public :: divd_h(:)        ! (n_elems)
  ! COEFFICIENT, et il n'est PAS choisi. Pi = p_c - h (div v) n'est pas homogene
  ! -- div v est l'inverse d'un temps, donc h div v est une VITESSE et pas une
  ! pression : il manque une impedance. Le solveur nodal en pression la fournit
  ! deja, son coefficient etant une longueur volume-sur-surface ponderee par
  ! l'impedance. Restreint a la sous-maille, ou toutes les sous-faces de SF(pc)
  ! appartiennent a la seule maille c et portent donc la seule impedance lam_c :
  !
  !   Pi_pcf = p_c - kappa_pc (div v)_pc      kappa_pc = lam_c |om_pc| / somme A_pcf
  !
  ! Dimensionnellement correct et SANS parametre libre. Et comme
  ! (div v)_pc = (b.v_p + reste)/|om_pc|, le volume s'elimine :
  !
  !   Pi_pcf = p_c - (lam_c / somme A_pcf) (b_pc.v_p + reste_pc)
  !
  ! PIEGE deja paye deux fois : sur un quadrangle de cote h extrude de dz,
  ! |om_pc| = h^2 dz/4 et les deux demi-aretes EN PLAN donnent somme A = h dz,
  ! d'ou kappa/lam = h/4 et le dz s'elimine EXACTEMENT. En incluant les faces de
  ! normale +/-e_z, somme A = h dz + h^2/2 et le rapport redevient fonction de
  ! dz, c'est-a-dire du parametre arbitraire dont toute la construction cherche
  ! a sortir. On les exclut donc sous boundary_2d.
  real(kind=DOUBLE), allocatable, public :: divd_isum(:)     ! 1/somme A, par sous-element
  real(kind=DOUBLE), allocatable, public :: divd_isumc(:,:)
  ! 'imp' : coefficient issu de l'impedance, ci-dessus, sans parametre libre.
  ! 'h'   : la forme d'origine Pi = p_c - divd_coef h_c (div v)_pc, gardee comme
  !         point de comparaison bien qu'elle ne soit pas homogene.
  character(len=8), public :: divd_mode = 'imp'
  ! SOUS-RELAXATION de la vitesse nodale. Le reste etant retarde d'un pas,
  ! l'iteration v_p <- D^-1 (R_pression - R v_p^ancien) a un rayon spectral
  ! mesure de 5/3 EXACTEMENT sur maillage cartesien, 1.6499 a 30 % de
  ! distorsion, INDEPENDANT du nombre de mailles et de dt. Ce n'est donc pas une
  ! instabilite conditionnelle qu'un CFL rattrape : chaque pas multiplie v_p par
  ! 1.65. Les translations uniformes sont exactement preservees (1e-16), donc
  ! l'invariance galileenne est intacte et ce sont des modes parasites, pas le
  ! noyau.
  !
  ! On relaxe donc l'iteration elle-meme, v_p <- (1-theta) v_p^ancien + theta v_p^resolu,
  ! dont le rayon est celui de (1-theta) I + theta It. Le seuil est NET autour de
  ! theta ~ 0.75 (rayon 1.000 a 0.75, 0.902 a 0.5) : il ne faut pas chercher une
  ! loi d'echelle en theta mais un seuil.
  !
  ! NE PAS iterer a la place : k balayages non amortis donnent rho^k, soit 1.7,
  ! 2.8, 4.6, 12.9 -- Jacobi diverge ici, donc iterer AGGRAVE strictement.
  real(kind=DOUBLE), public :: divd_theta = 0.5_DOUBLE
  ! FORME ISENTROPIQUE. Sur quadrangle regulier kappa = lam h/4, lam = rho a et
  ! p = rho a^2/gamma, donc kappa D / p = gamma M_div / 4 avec M_div = h|div v|/a.
  ! La fermeture LINEAIRE Pi = p_c (1 - gamma M_div/4) devient donc NEGATIVE des
  ! que M_div > 4/gamma = 2.857, c'est-a-dire des que la vitesse varie de plus de
  ! 2.86 vitesses du son a travers une maille. Une pression de sous-face negative
  ! est une TENSION, et elle apparait precisement la ou la divergence est grande
  ! -- donc dans toutes les detentes et compressions fortes. Rien dans la forme
  ! lineaire ne l'en empeche.
  !
  ! Elle n'est que le developpement au premier ordre de la relation isentropique.
  ! On garde la relation :
  !     Pi_pcf = p_c * ( 1 - (gamma-1)/(2 gamma) * q )^(2 gamma/(gamma-1))
  ! avec q = kappa_pc D_pc / p_c, c'est-a-dire le decrement relatif lineaire.
  ! Au premier ordre elle redonne Pi = p_c (1 - q), donc le gain bas Mach est
  ! INTACT : l'ecart est quadratique (rapport 4 quand on divise M par deux) et
  ! n'apparait que la ou la forme lineaire est de toute facon fausse. La base ne
  ! s'annule qu'a q = 2 gamma/(gamma-1) = 7, soit M_div = 20 : sept fois plus de
  ! marge.
  ! Mise en oeuvre SANS Newton et SANS perte de conservation. La forme
  ! isentropique est la forme lineaire a coefficient reechelonne :
  !     Pi = p phi(q) = p - q p psi(q) = p_c - kappa_pc D_pc psi(q)
  !     psi(q) = (1 - phi(q))/q
  ! donc kappa_eff = kappa psi(q). Le systeme nodal reste LINEAIRE en v_p, le
  ! bilan somme A Pi n = 0 est exactement satisfait, et la conservation de
  ! quantite de mouvement est intacte. Ce n'est pas une approximation de la forme
  ! isentropique, c'est elle.
  !
  ! psi est monotone decroissante, bornee dans [0.025, 1] sur q dans [0,40], et
  ! psi(0) = 1 EXACTEMENT : le coefficient est donc inchange au bas Mach et la
  ! correction ne mord que la ou la forme lineaire etait fausse. C'est un
  ! limiteur sur le coefficient de dissipation, pas un plafond sur Pi.
  !
  ! COUT, et il ne faut pas le presenter comme gratuit : q depend de v_p, donc il
  ! est evalue a l'itere precedent -- un retard de plus. Il n'est pas de meme
  ! nature que celui qui a coute la campagne : la c'etait une contribution
  ! VECTORIELLE ADDITIVE du meme ordre que la part implicite (79 % de la norme),
  ! ici c'est un SCALAIRE MULTIPLICATIF borne et monotone sur le coefficient de
  ! diffusion. Son rayon spectral n'est pas calcule ; en cas d'oscillation sur q
  ! il faut lui appliquer le meme theta qu'au reste.
  logical, public :: divd_nl = .false.
  real(kind=DOUBLE), allocatable, public :: divd_psi(:)
  ! Compteurs. On sort une DISTRIBUTION et des COMPTES, pas un maximum : un max
  ! mesure toujours la sous-population extreme et ne dit pas si le phenomene est
  ! marginal ou massif.
  integer(kind=ENTIER_D), public :: divd_n_pi = 0, divd_n_pineg = 0
  integer(kind=ENTIER_D), public :: divd_md(7) = 0
  ! Base negative = formation de VIDE. Comptee et signalee, jamais ecrasee par un
  ! max(0,.) : on vient de mesurer ce que coute un max(0,.) silencieux sur lambda,
  ! qui transforme une pression negative en impedance NULLE au lieu d'une erreur.
  integer(kind=ENTIER_D), public :: divd_n_vide = 0
  ! Les compteurs ci-dessus sont CUMULES sur tout le run, divergence comprise :
  ! une fois la solution partie, |div v| explose et ils se remplissent seuls.
  ! Ils ne distinguent donc pas une cause d'un symptome -- exactement le piege
  ! des 119 770 matrices "singulieres", qui etaient en aval de l'effondrement.
  ! On sort donc aussi l'INCREMENT des dix premiers pas, avant que quoi que ce
  ! soit ne parte : si la tension y est deja massive, elle est causale.
  integer(kind=ENTIER), public :: divd_step = 0
  integer(kind=ENTIER_D), public :: divd_prev(9) = 0
  logical, public :: divd_check_done = .false.
  integer(kind=ENTIER), public :: divd_dbg = 0
  ! Noeuds a matrice nodale SINGULIERE. Comptes et rapportes plutot que
  ! laisses produire un NaN silencieux : rang(M_p) <= nombre de mailles au
  ! noeud, donc un maillage a noeuds pauvres les rend inevitables et il faut
  ! le voir dans le log.
  integer(kind=ENTIER_D), public :: divd_n_sing = 0
  integer(kind=ENTIER_D), public :: divd_sing_bnd = 0, divd_sing_int = 0
  integer(kind=ENTIER), public :: divd_sing_print = 0
  integer(kind=ENTIER), public :: divd_wdbg = 0, divd_wnode = -1
  ! BORNE du traceur VPDIVW. Il s'enclenche SEUL au premier noeud interieur
  ! singulier (divd_wdbg = 1, jamais remis a 0) et n'etait borne par rien : sur
  ! vp_div il s'allume des le pas 1 -- le noyau de dimension 2n en produit
  ! toujours -- puis imprime une ligne par sous-face ET par pas jusqu'a la fin,
  ! soit 75 % du log et des centaines de Mo sur un run long.
  integer(kind=ENTIER), public :: divd_wmax = 200, divd_w_printed = 0
  real(kind=DOUBLE), allocatable, public :: divd_bc(:,:), divd_restc(:,:)
  real(kind=DOUBLE), allocatable, public :: divd_prod_cell(:)
  integer(kind=ENTIER_D), public :: divd_n_cell_neg = 0, divd_n_cell = 0
  real(kind=DOUBLE), public :: divd_prod_min = 0.0_DOUBLE, divd_prod_tot = 0.0_DOUBLE
  logical, public :: diag_iso = .false.
  real(kind=DOUBLE), allocatable, public :: iso_eps0(:), iso_tau0(:)
  logical, public :: kap_check = .true.
  logical, public :: kap_check_done = .false.
  ! Decomposition de la production d'energie interne vue par la maille :
  !   p_c - Pi_chapeau = theta*(p_c - p_p) + lam_c*delta
  ! donc production = lam_c*somme(A delta^2)  +  theta*somme(A (p_c - p_p) delta)
  ! Le premier terme est une somme de carres, le second est INDEFINI. On ne
  ! dispose donc pas d'une preuve de dissipativite, mais d'un plancher plus un
  ! terme libre : les deux sont mesures separement.
  real(kind=DOUBLE), public :: pih_sq_step = 0.0_DOUBLE, pih_free_step = 0.0_DOUBLE
  real(kind=DOUBLE), public :: pih_sq_work = 0.0_DOUBLE, pih_free_work = 0.0_DOUBLE
  ! residu du bilan nodal : |v_p . somme(A Pi_chapeau n)|, qui doit etre nul.
  real(kind=DOUBLE), public :: pih_res_max = 0.0_DOUBLE
  integer(kind=ENTIER), public :: pos_n_clip_step = 0, pos_n_sf_step = 0
  integer(kind=ENTIER), public :: pos_n_infeas_step = 0
  real(kind=DOUBLE), public :: pos_defect_step = 0.0_DOUBLE
  ! cumuls sur tout le run : ENTIER_D obligatoire, le compteur de sous-faces
  ! visitees depasse 2^31 des le maillage hexa 41^3 (mesure : il est reparti
  ! negatif et le pourcentage affiche valait 8.8e9 %).
  integer(kind=ENTIER_D), public :: pos_n_clip = 0, pos_n_sf = 0
  integer(kind=ENTIER), public :: pos_n_clip_max = 0
  integer(kind=ENTIER_D), public :: pos_n_infeas = 0, pos_n_clip_bnd = 0
  integer(kind=ENTIER_D), public :: pos_n_steps_infeas = 0
  ! Compteurs de la variante univaluee. On sort des COMPTES et une distribution,
  ! pas un extremum : 'une sous-face est corrigee' et 'toutes le sont' se lisent
  ! pareil sur un max.
  integer(kind=ENTIER_D), public :: sv_n_sf = 0, sv_n_corr = 0, sv_n_bnd = 0
  integer(kind=ENTIER_D), public :: sv_n_bnd_corr = 0, sv_n_infeas = 0
  real(kind=DOUBLE), public :: sv_d2_min = 1.0e300_DOUBLE
  real(kind=DOUBLE), public :: sv_ecart_max = 0.0_DOUBLE, sv_ecart_tot = 0.0_DOUBLE
  integer(kind=ENTIER_D), public :: sv_n_bnd_z = 0, sv_n_bnd_z_corr = 0
  integer(kind=ENTIER_D), public :: sv_n_nan = 0, sv_n_ecart = 0
  real(kind=DOUBLE), public :: sv_ec_in_max = 0.0_DOUBLE, sv_ec_in_tot = 0.0_DOUBLE
  real(kind=DOUBLE), public :: sv_ec_out_max = 0.0_DOUBLE, sv_ec_out_tot = 0.0_DOUBLE
  integer(kind=ENTIER_D), public :: sv_n_in = 0, sv_n_out = 0
  integer(kind=ENTIER_D), public :: sv_n_corr_nan = 0
  integer(kind=ENTIER_D), public :: sv_n_corr_froid = 0, sv_n_corr_chaud = 0
  integer(kind=ENTIER_D), public :: sv_n_corr_mixte = 0
  real(kind=DOUBLE), public :: sv_pmax_corr = -1.0e300_DOUBLE
  real(kind=DOUBLE), public :: sv_pmin_corr = 1.0e300_DOUBLE
  ! Par PAS : le compte de corrections et le rayon moyen des sous-faces corrigees.
  ! Si la correction vit sur une couronne qui suit le front, le compte doit croitre
  ! comme le perimetre et le rayon doit se deplacer vers l'exterieur. Un compte plat
  ! refute la couronne.
  integer(kind=ENTIER), public :: sv_step_n = 0
  real(kind=DOUBLE), public :: sv_step_r = 0.0_DOUBLE
  real(kind=DOUBLE), public :: sv_step_rmin = 1.0e300_DOUBLE, sv_step_rmax = 0.0_DOUBLE
  ! Par PAS : intervalle vide et delta2 minimal. Decide si l'hypothese du theoreme
  ! (delta2 > 0) tombe sur un etat SAIN ou seulement apres degradation.
  integer(kind=ENTIER), public :: sv_st_infeas = 0
  real(kind=DOUBLE), public :: sv_st_d2min = 1.0e300_DOUBLE
  ! defaut de conservation INTEGRE en temps, avec son signe : c'est la seule forme
  ! comparable a une derive d'energie totale mesuree independamment.
  real(kind=DOUBLE), public :: pos_defect_work = 0.0_DOUBLE
  ! Ecretages subis par chaque MAILLE, cumules sur le run. Sert a dire si les
  ! mailles qui depassent la compression maximale sont celles ou le limiteur mord.
  integer(kind=ENTIER), allocatable, public :: pos_clip_cell(:)
  ! Diagnostic de conservation. Les sous-faces INTERIEURES telescopent
  ! exactement (flux_l(5) et flux_r(5) sont la meme expression), donc la seule
  ! source possible de variation d'energie totale est la sous-face de BORD, qui
  ! n'a pas de partenaire. On somme sa puissance, somme_bord A (pvp.n) ;
  ! lagrange_main l'integre en temps avec le dt REELLEMENT utilise.
  real(kind=DOUBLE), public :: pvp_bnd_power = 0.0_DOUBLE
  real(kind=DOUBLE), public :: pvp_bnd_work = 0.0_DOUBLE
  ! Somme des VALEURS ABSOLUES du meme terme. Distingue les deux facons d'obtenir
  ! zero : si abs est nul aussi, pvp.n est reellement nul au bord ; si abs est
  ! grand alors que la somme signee est nulle, il y a compensation entre
  ! sous-faces.
  real(kind=DOUBLE), public :: pvp_bnd_abs = 0.0_DOUBLE
  ! Plancher relatif sur l'energie interne predite : une maille est declaree
  ! inadmissible si e_new < eps * e_old.
  real(kind=DOUBLE), public :: pvp_eps = 1.0e-3_DOUBLE
  integer(kind=ENTIER), public :: pvp_lim_iter = 8
  real(kind=DOUBLE), public :: pvp_theta_min = 1.0_DOUBLE
  integer(kind=ENTIER), public :: pvp_n_limited = 0

  ! Sedov (init=1) : ou deposer l'energie.
  ! .false. (defaut) = comportement historique, tout dans l'UNIQUE cellule dont le
  !   centroide est le plus proche de l'origine.
  ! .true.  = si un SOMMET du maillage est a l'origine, repartir l'energie sur toutes
  !   les cellules adjacentes a ce sommet, au prorata de leur volume (l'energie
  !   specifique deposee est donc uniforme sur la zone) ; si c'est au contraire une
  !   CELLULE qui est centree sur l'origine, tout y est depose comme avant.
  ! Sans cela, sur un maillage ou l'origine est un sommet (Voronoi radial sans germe
  ! central, quart de disque a la Vilar), toute l'energie tombe dans une seule des
  ! cellules qui entourent le point, ce qui brise la symetrie du depot.
  logical, public :: sedov_nodal_deposit = .false.
  ! Energie totale deposee. <= 0 (defaut) = valeur historique, calibree pour un domaine
  ! COMPLET (front a r=1, pic 6 a t=1). Sur un quart de plan avec conditions de symetrie
  ! il faut le quart de cette valeur, sinon le choc est quatre fois trop fort.
  real(kind=DOUBLE), public :: sedov_energy = -1.0_DOUBLE

  ! Mach number of the Gresho vortex (init = 11): sets the pressure level
  ! p0 = 1/(gamma*Ma^2). Settable from the INPUT_PARAM namelist.
  real(kind=DOUBLE) :: gresho_mach = 1.0e-2_DOUBLE
  ! --- init = 13 : onde acoustique plane traversant un vortex de Gresho
  ! (Bruel, Delmas, Jung, Perrier 2019, section 4.2.2). Les deux morceaux sont
  ! activables separement : le montage se fait par etapes, le vortex seul puis
  ! l'onde seule servant de temoins a leur superposition.
  logical, public :: vac_vortex = .true.
  logical, public :: vac_wave = .true.
  ! > 0 overrides the per-init default gamma (see lagrange_main); 0 = keep default
  real(kind=DOUBLE) :: gamma_uniform = 0.0_DOUBLE

  ! Taylor-Green vortex (init = 10) -- the reference convergence test case for
  ! cell-centred Lagrangian hydrodynamics (Vilar, PhD 2012, Sec. 4.2.4). C0 has to
  ! keep the pressure positive: min[cos+cos] = -2, so C0 > rho0*U0^2/2.
  real(kind=DOUBLE), parameter :: tg_rho0 = 1.0_DOUBLE
  real(kind=DOUBLE), parameter :: tg_u0   = 1.0_DOUBLE
  real(kind=DOUBLE), parameter :: tg_c0   = 1.0_DOUBLE
  ! --- plafond de positivite du terme de divergence de SOUS-ELEMENT (md_kappa_cap)
  integer(kind=ENTIER), parameter, public :: MD_OK        = 0  ! plafond utilisable
  integer(kind=ENTIER), parameter, public :: MD_EMPTY     = 1  ! delta2(pbar_f) <= 0
  integer(kind=ENTIER), parameter, public :: MD_BASE_OUT  = 2  ! |D_p| > R
  integer(kind=ENTIER), parameter, public :: MD_NO_CONSTR = 3  ! delta_pc nul

contains
  subroutine compute_rhs_lagrange(mesh, sol, vp, dt, rhs, n_bc, bc_type, bc_val, b2d, mass, &
      gamma_arr, vp_is_imposed, second_order, grad_v, grad_p, div_v, alpha_p_arr)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), intent(in) :: dt
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mass
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    logical, dimension(mesh%n_vert), intent(in) :: vp_is_imposed
    integer(kind=ENTIER), intent(in) :: n_bc
    character(len=255), dimension(n_bc) :: bc_type
    real(kind=DOUBLE), dimension(5, n_bc) :: bc_val
    logical, intent(in) :: b2d
    logical, intent(in) :: second_order
    real(kind=DOUBLE), dimension(3, 3, mesh%n_elems), intent(in) :: grad_v
    real(kind=DOUBLE), dimension(3, mesh%n_elems), intent(in) :: grad_p
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: div_v
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(out), optional :: alpha_p_arr

    integer(kind=ENTIER) :: i, j, k, nsfn, nsen, idse, ide
    real(kind=DOUBLE), dimension(:,:), allocatable :: rhs_o1, rhs_o2
    real(kind=DOUBLE), dimension(3) :: vp_o1_i, vp_o2_i
    real(kind=DOUBLE) :: alpha_p, f_pc

    rhs = 0.0_DOUBLE
    if (present(alpha_p_arr)) alpha_p_arr = 1.0_DOUBLE
    vp_o1_nan_n = 0; vp_o1_nan_node = -1

    do i = 1, mesh%n_vert
      nsfn = mesh%vert(i)%n_sub_faces_neigh
      nsen = mesh%vert(i)%n_sub_elems_neigh
      allocate(rhs_o1(5, nsen))

      vp_o1_i = vp(:, i)
      call compute_rhs_around_node(mesh, sol, i, nsfn, nsen, vp_o1_i, vp_is_imposed(i), &
        dt, rhs_o1, n_bc, bc_type, bc_val, b2d, mass, gamma_arr, &
        .false., grad_v, grad_p, div_v)

      ! Mp peut etre singuliere (noeuds de milieu d arete) et inverse_3_by_3 divise par son
      ! determinant sans garde : le NaN d ordre 1 passait inapercu. On le COMPTE.
      if (any(rhs_o1 /= rhs_o1) .or. any(vp_o1_i /= vp_o1_i)) then
        vp_o1_nan_n = vp_o1_nan_n + 1
        if (vp_o1_nan_node < 0) vp_o1_nan_node = i
      end if

      if (second_order) then
        allocate(rhs_o2(5, nsen))
        vp_o2_i = vp(:, i)
        call compute_rhs_around_node(mesh, sol, i, nsfn, nsen, vp_o2_i, vp_is_imposed(i), &
          dt, rhs_o2, n_bc, bc_type, bc_val, b2d, mass, gamma_arr, &
          .true., grad_v, grad_p, div_v)

        ! If rhs_o2 or vp_o2_i has any NaN (negative reconstructed pressure → sqrt(neg) → NaN
        ! propagating through the nodal solver), skip blending entirely for this node.
        if (any(rhs_o2 /= rhs_o2) .or. any(vp_o2_i /= vp_o2_i)) then
          alpha_p = -0.01_DOUBLE
        else
          alpha_p = 1.0_DOUBLE
          if (convex_limiting) then
            do k = 1, nsen
              idse = mesh%vert(i)%sub_elem_neigh(k)
              ide = mesh%sub_elem(idse)%mesh_elem
              if (.not. mesh%elem(ide)%is_ghost) then
                ! rhs est divise par la masse de la MAILLE : pour tester la mise a jour de la
                ! SOUS-CELLULE, remettre a l'echelle par V_c/V_pc (densite uniforme dans la maille).
                f_pc = mesh%elem(ide)%volume &
                  /max(mesh%sub_elem(idse)%volume, 1.0e-300_DOUBLE)
                alpha_p = min(alpha_p, &
                  omega_pos(sol(:, ide) + dt*f_pc*rhs_o1(:, k), &
                            dt*f_pc*(rhs_o2(:, k) - rhs_o1(:, k))))
              end if
            end do
          end if
        end if

        if (present(alpha_p_arr)) alpha_p_arr(i) = alpha_p

        if (alpha_p > 0.0_DOUBLE) then
          vp(:, i) = vp_o1_i + alpha_p*(vp_o2_i - vp_o1_i)
          do j = 1, nsen
            idse = mesh%vert(i)%sub_elem_neigh(j)
            ide = mesh%sub_elem(idse)%mesh_elem
            rhs(:, ide) = rhs(:, ide) + rhs_o1(:,j) + alpha_p*(rhs_o2(:,j) - rhs_o1(:,j))
          end do
        else
          vp(:, i) = vp_o1_i
          do j = 1, nsen
            idse = mesh%vert(i)%sub_elem_neigh(j)
            ide = mesh%sub_elem(idse)%mesh_elem
            rhs(:, ide) = rhs(:, ide) + rhs_o1(:,j)
          end do
        end if
        deallocate(rhs_o2)
      else
        vp(:, i) = vp_o1_i
        do j = 1, nsen
          idse = mesh%vert(i)%sub_elem_neigh(j)
          ide = mesh%sub_elem(idse)%mesh_elem
          rhs(:, ide) = rhs(:, ide) + rhs_o1(:,j)
        end do
      end if

      deallocate(rhs_o1)
    end do

    ! --- finalisation du diagnostic de divergence
    if (diag_divratio) then
      block
        integer(kind=ENTIER) :: ic
        real(kind=DOUBLE) :: divc, lamc, hc, ac
        do ic = 1, mesh%n_elems
          if (mesh%elem(ic)%is_ghost) cycle
          if (mesh%elem(ic)%volume <= 0.0_DOUBLE) cycle
          divc = dr_vdot(ic)/mesh%elem(ic)%volume
          ac = sqrt(max(0.0_DOUBLE, gamma_arr(ic) &
                 *pressure(sol(:, ic), gamma_arr(ic))/sol(1, ic)))
          lamc = ac/sol(1, ic)
          ! h_c = V^(1/3). RESERVE : sur un maillage 2D obtenu par extrusion
          ! l'epaisseur en z est arbitraire et entre dans h_c, donc le rapport
          ! porte un facteur dz^(1/3) libre. Les deux sommes sont imprimees
          ! separement pour pouvoir renormaliser.
          hc = mesh%elem(ic)%volume**(1.0_DOUBLE/3.0_DOUBLE)
          dr_sq_div = dr_sq_div + lamc*hc*divc**2*mesh%elem(ic)%volume
          if (diag_firstprint .and. .not. dr_done) dr_max_vdot = max(dr_max_vdot, abs(divc))
        end do
        if (diag_firstprint .and. .not. dr_done) then
          print*, "FIRSTSTEP max_p |(div v)_p|  (sur les vitesses de MAILLE) =", dr_max_divp
          print*, "FIRSTSTEP max_c |Vpoint_c/V_c| (sur les vitesses NODALES) =", dr_max_vdot
          dr_done = .true.
        end if
        dr_vdot = 0.0_DOUBLE
      end block
    end if
  end subroutine compute_rhs_lagrange

  ! ===================================================================
  ! 'vp_pif' : l'inconnue est la PRESSION DE SOUS-FACE Pi_f, un scalaire par
  ! sous-face, exactement comme la temperature de ns_diffusion_module. Dans la
  ! sous-maille on en tire une perturbation de vitesse par la relation
  ! acoustique, et la CONTINUITE DE LA VITESSE NORMALE sur chaque sous-face
  ! ferme le systeme : une equation scalaire par sous-face, |SF(p)| inconnues.
  !
  !   U_l^(c) = v_c.n_l - (1/lam_c) somme_m (A_m/D_pc) (Pi_m - p_c) (n_m.n_l)
  !   D_pc    = moyenne des A_m sur les sous-faces de la sous-maille
  !
  ! et somme sur les DEUX sous-mailles de U_l = 0, normales sortantes de chacune.
  !
  ! PIEGE, et il n'est pas transposable depuis la diffusion : le poids y est
  ! A_l A_m/|omega| parce que le flux est K grad(theta) ; ici la perturbation de
  ! vitesse vaut Delta p / lambda, une DIFFERENCE de pression et non un gradient.
  ! D'ou le poids SANS DIMENSION A_m/D_pc. Le symptome de l'erreur est propre :
  ! la part pression sort exacte et la part vitesse exactement d'un facteur 2.
  !
  ! BORD : une sous-face sans voisin impose U_l = 0 (paroi impermeable), ce que
  ! l'assemblage a UN SEUL cote donne automatiquement -- la ligne devient
  ! U_l^(c) = 0. Rien de special a ecrire.
  !
  ! Les deux flux n'utilisent QUE (Pi_f, U_f), jamais v_p :
  !   m_c dv_c/dt = - somme A Pi_f n       m_c de_c/dt = - somme A Pi_f U_f
  ! donc les deux conservations sont exactes par univocite de Pi et continuite
  ! de U, et le defaut de v_p ne casse aucun bilan : il ne sert qu'a deplacer le
  ! maillage.
  ! ===================================================================
  subroutine compute_rhs_lagrange_pif(mesh, sol, vp, dt, rhs, n_bc, bc_type, &
      bc_val, b2d, mass, gamma_arr, vp_is_imposed)
    use linear_solver_module, only: lu_solve, inverse_3_by_3
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), intent(in) :: dt
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mass, gamma_arr
    logical, dimension(mesh%n_vert), intent(in) :: vp_is_imposed
    integer(kind=ENTIER), intent(in) :: n_bc
    character(len=255), dimension(n_bc) :: bc_type
    real(kind=DOUBLE), dimension(5, n_bc) :: bc_val
    logical, intent(in) :: b2d

    integer(kind=ENTIER) :: i, l, m, nsfn, isfl, isfp, isfl_loc, isfp_loc
    integer(kind=ENTIER) :: ifl, ifp, ise, iel, side, ide
    real(kind=DOUBLE), allocatable :: Nm(:,:), bv(:), Pif(:), Uf(:)
    real(kind=DOUBLE), allocatable :: lam_se(:), du_num(:), du_den(:)
    real(kind=DOUBLE), allocatable, save :: lam_keep(:)
    integer(kind=ENTIER) :: nsen, nit, iter_d, k, kk
    real(kind=DOUBLE) :: ac
    real(kind=DOUBLE), dimension(3) :: norml, normp, Mv
    real(kind=DOUBLE), dimension(3,3) :: Mp, Mp_inv
    real(kind=DOUBLE) :: Dpc, lam, pc, w, acc, area
    real(kind=DOUBLE), dimension(5) :: flux

    rhs = 0.0_DOUBLE
    do i = 1, mesh%n_vert
      nsfn = mesh%vert(i)%n_sub_faces_neigh
      if (nsfn <= 0) cycle
      allocate(Nm(nsfn, nsfn), bv(nsfn), Pif(nsfn), Uf(nsfn))
      Nm = 0.0_DOUBLE; bv = 0.0_DOUBLE; Pif = 0.0_DOUBLE; Uf = 0.0_DOUBLE

      ! impedance par SOUS-ELEMENT, acoustique a l'iteration 0
      nsen = mesh%vert(i)%n_sub_elems_neigh
      allocate(lam_se(nsen), du_num(nsen), du_den(nsen))
      lam_se = 0.0_DOUBLE
      do k = 1, nsen
        ise = mesh%vert(i)%sub_elem_neigh(k)
        iel = mesh%sub_elem(ise)%mesh_elem
        pc  = pressure(sol(:, iel), gamma_arr(iel))
        lam_se(mesh%sub_elem(ise)%id_loc_around_node) = &
          sqrt(max(0.0_DOUBLE, gamma_arr(iel)*pc*sol(1, iel)))/sol(1, iel)
      end do
      nit = 1
      if (pif_duk) nit = max(1, pif_duk_it)

      do iter_d = 1, nit
      Nm = 0.0_DOUBLE; bv = 0.0_DOUBLE
      do l = 1, nsfn
        isfl = mesh%vert(i)%sub_face_neigh(l)
        isfl_loc = mesh%sub_face(isfl)%id_loc_around_node
        ifl = mesh%sub_face(isfl)%mesh_face
        do side = 1, 2
          if (side == 1) then
            ise = mesh%sub_face(isfl)%left_sub_elem_neigh
            iel = mesh%face(ifl)%left_neigh
            norml = mesh%sub_face(isfl)%norm
          else
            ise = mesh%sub_face(isfl)%right_sub_elem_neigh
            iel = mesh%face(ifl)%right_neigh
            norml = -mesh%sub_face(isfl)%norm
          end if
          if (ise <= 0 .or. iel <= 0) cycle
          pc  = pressure(sol(:, iel), gamma_arr(iel))
          lam = lam_se(mesh%sub_elem(ise)%id_loc_around_node)
          if (lam <= 0.0_DOUBLE) cycle
          ! D_pc : moyenne des aires des sous-faces de la sous-maille
          Dpc = 0.0_DOUBLE
          do m = 1, mesh%sub_elem(ise)%n_sub_faces
            Dpc = Dpc + mesh%sub_face(mesh%sub_elem(ise)%sub_face(m))%area
          end do
          Dpc = Dpc/real(mesh%sub_elem(ise)%n_sub_faces, DOUBLE)
          if (Dpc <= 0.0_DOUBLE) cycle
          acc = 0.0_DOUBLE
          do m = 1, mesh%sub_elem(ise)%n_sub_faces
            isfp = mesh%sub_elem(ise)%sub_face(m)
            isfp_loc = mesh%sub_face(isfp)%id_loc_around_node
            if (mesh%sub_face(isfp)%left_sub_elem_neigh == ise) then
              normp = mesh%sub_face(isfp)%norm
            else
              normp = -mesh%sub_face(isfp)%norm
            end if
            w = (mesh%sub_face(isfp)%area/Dpc)/lam*dot_product(normp, norml)
            Nm(isfl_loc, isfp_loc) = Nm(isfl_loc, isfp_loc) + w
            acc = acc + w
          end do
          bv(isfl_loc) = bv(isfl_loc) + acc*pc &
                       + dot_product(sol(2:4, iel), norml)
        end do
        ! pression imposee : la ligne devient Pi_l = bc_val
        if (mesh%face(ifl)%right_neigh < 0) then
          if (trim(bc_type(-mesh%face(ifl)%right_neigh)) == 'pressure') then
            Nm(isfl_loc, :) = 0.0_DOUBLE
            Nm(isfl_loc, isfl_loc) = 1.0_DOUBLE
            bv(isfl_loc) = bc_val(1, -mesh%face(ifl)%right_neigh)
          end if
        end if
      end do

      call lu_solve(nsfn, Nm, Pif, bv)

      ! --- mise a jour de l'impedance de choc, sauf au dernier tour.
      ! delta u par sous-element agrege en MOYENNE QUADRATIQUE ponderee par
      ! l'aire, et non par un max : un max sur une population heterogene mesure
      ! sa sous-face la plus extreme, pas l'etat du coin.
      if (pif_duk .and. iter_d < nit) then
        du_num = 0.0_DOUBLE; du_den = 0.0_DOUBLE
        do l = 1, nsfn
          isfl = mesh%vert(i)%sub_face_neigh(l)
          isfl_loc = mesh%sub_face(isfl)%id_loc_around_node
          ifl = mesh%sub_face(isfl)%mesh_face
          do side = 1, 2
            if (side == 1) then
              ise = mesh%sub_face(isfl)%left_sub_elem_neigh
              iel = mesh%face(ifl)%left_neigh
              norml = mesh%sub_face(isfl)%norm
            else
              ise = mesh%sub_face(isfl)%right_sub_elem_neigh
              iel = mesh%face(ifl)%right_neigh
              norml = -mesh%sub_face(isfl)%norm
            end if
            if (ise <= 0 .or. iel <= 0) cycle
            k = mesh%sub_elem(ise)%id_loc_around_node
            pc  = pressure(sol(:, iel), gamma_arr(iel))
            lam = lam_se(k)
            if (lam <= 0.0_DOUBLE) cycle
            ! delta u = (Pi - p_c)/lambda, la perturbation acoustique du cote c
            du_num(k) = du_num(k) + mesh%sub_face(isfl)%area &
                      * ((Pif(isfl_loc) - pc)/lam)**2
            du_den(k) = du_den(k) + mesh%sub_face(isfl)%area
          end do
        end do
        do k = 1, nsen
          ise = mesh%vert(i)%sub_elem_neigh(k)
          iel = mesh%sub_elem(ise)%mesh_elem
          kk = mesh%sub_elem(ise)%id_loc_around_node
          if (du_den(kk) <= 0.0_DOUBLE) cycle
          pc  = pressure(sol(:, iel), gamma_arr(iel))
          ac  = sqrt(max(0.0_DOUBLE, gamma_arr(iel)*pc*sol(1, iel)))
          lam_se(kk) = (ac + 0.5_DOUBLE*(gamma_arr(iel) + 1.0_DOUBLE) &
                       *sqrt(du_num(kk)/du_den(kk)))/sol(1, iel)
        end do
      end if
      end do
      ! l'impedance convergee sert aussi au calcul de U_f, sans quoi U ne serait
      ! pas la vitesse que Pi implique et la continuite serait rompue.
      if (.not. allocated(lam_keep)) allocate(lam_keep(nsen))
      if (size(lam_keep) /= nsen) then
        deallocate(lam_keep); allocate(lam_keep(nsen))
      end if
      lam_keep = lam_se
      deallocate(lam_se, du_num, du_den)

      ! --- U_f depuis la sous-maille GAUCHE ; la continuite garantit que le
      ! cote droit donnerait la meme valeur, et c'est elle qui rend les deux
      ! bilans exacts.
      do l = 1, nsfn
        isfl = mesh%vert(i)%sub_face_neigh(l)
        isfl_loc = mesh%sub_face(isfl)%id_loc_around_node
        ifl = mesh%sub_face(isfl)%mesh_face
        ise = mesh%sub_face(isfl)%left_sub_elem_neigh
        iel = mesh%face(ifl)%left_neigh
        if (ise <= 0 .or. iel <= 0) cycle
        norml = mesh%sub_face(isfl)%norm
        pc  = pressure(sol(:, iel), gamma_arr(iel))
        lam = lam_keep(mesh%sub_elem(ise)%id_loc_around_node)
        if (lam <= 0.0_DOUBLE) cycle
        Dpc = 0.0_DOUBLE
        do m = 1, mesh%sub_elem(ise)%n_sub_faces
          Dpc = Dpc + mesh%sub_face(mesh%sub_elem(ise)%sub_face(m))%area
        end do
        Dpc = Dpc/real(mesh%sub_elem(ise)%n_sub_faces, DOUBLE)
        acc = 0.0_DOUBLE
        do m = 1, mesh%sub_elem(ise)%n_sub_faces
          isfp = mesh%sub_elem(ise)%sub_face(m)
          isfp_loc = mesh%sub_face(isfp)%id_loc_around_node
          if (mesh%sub_face(isfp)%left_sub_elem_neigh == ise) then
            normp = mesh%sub_face(isfp)%norm
          else
            normp = -mesh%sub_face(isfp)%norm
          end if
          acc = acc + (mesh%sub_face(isfp)%area/Dpc)/lam &
                *(Pif(isfp_loc) - pc)*dot_product(normp, norml)
        end do
        Uf(isfl_loc) = dot_product(sol(2:4, iel), norml) - acc
      end do

      ! --- vitesse nodale par moindres carres sur les sous-faces :
      ! M_p v_p = somme A n U_f, avec M_p = somme A n(x)n -- exactement la
      ! matrice de vp, deja dans le code, SPD et de rang plein.
      Mp = 0.0_DOUBLE; Mv = 0.0_DOUBLE
      do l = 1, nsfn
        isfl = mesh%vert(i)%sub_face_neigh(l)
        isfl_loc = mesh%sub_face(isfl)%id_loc_around_node
        norml = mesh%sub_face(isfl)%norm
        area = mesh%sub_face(isfl)%area
        Mp = Mp + area*tensor_product(norml, norml)
        Mv = Mv + area*Uf(isfl_loc)*norml
      end do
      if (.not. vp_is_imposed(i) .and. .not. pif_vp_classic &
          .and. maxval(abs(Mp)) > 0.0_DOUBLE) then
        call inverse_3_by_3(Mp, Mp_inv)
        vp(:, i) = matmul(Mp_inv, Mv)
      end if

      ! --- flux : ils n'utilisent QUE Pi et U
      do l = 1, nsfn
        isfl = mesh%vert(i)%sub_face_neigh(l)
        isfl_loc = mesh%sub_face(isfl)%id_loc_around_node
        ifl = mesh%sub_face(isfl)%mesh_face
        area = mesh%sub_face(isfl)%area
        norml = mesh%sub_face(isfl)%norm
        ise = mesh%sub_face(isfl)%left_sub_elem_neigh
        iel = mesh%face(ifl)%left_neigh
        if (ise > 0 .and. iel > 0) then
          if (mesh%sub_elem(ise)%mesh_vert == i) then
            ide = mesh%sub_elem(ise)%mesh_elem
            flux(1)   = -Uf(isfl_loc)
            flux(2:4) = Pif(isfl_loc)*norml
            flux(5)   = Pif(isfl_loc)*Uf(isfl_loc)
            rhs(:, ide) = rhs(:, ide) - area/mass(iel)*flux
            call fdiag_add(ide, area, Pif(isfl_loc), norml, mesh%face(ifl)%right_neigh)
          end if
        end if
        ise = mesh%sub_face(isfl)%right_sub_elem_neigh
        iel = mesh%face(ifl)%right_neigh
        if (ise > 0 .and. iel > 0) then
          if (mesh%sub_elem(ise)%mesh_vert == i) then
            ide = mesh%sub_elem(ise)%mesh_elem
            flux(1)   = -Uf(isfl_loc)
            flux(2:4) = Pif(isfl_loc)*norml
            flux(5)   = Pif(isfl_loc)*Uf(isfl_loc)
            rhs(:, ide) = rhs(:, ide) + area/mass(iel)*flux
            call fdiag_add(ide, area, Pif(isfl_loc), -norml, mesh%face(ifl)%left_neigh)
          end if
        end if
      end do

      deallocate(Nm, bv, Pif, Uf)
    end do
  end subroutine compute_rhs_lagrange_pif

  ! Resout le systeme de sous-face pour un etat SYNTHETIQUE (p, v, lam imposes),
  ! et rend Pi. Sert aux controles : on verifie la fermeture avant tout calcul.
  subroutine pif_solve_synth(mesh, i, pcell, vcell, lamc, Pif, nsfn)
    use linear_solver_module, only: lu_solve
    implicit none
    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: pcell, lamc
    real(kind=DOUBLE), dimension(3, mesh%n_elems), intent(in) :: vcell
    real(kind=DOUBLE), dimension(:), allocatable, intent(out) :: Pif
    integer(kind=ENTIER), intent(out) :: nsfn
    real(kind=DOUBLE), allocatable :: Nm(:,:), bv(:)
    integer(kind=ENTIER) :: l, m, isfl, isfp, isfl_loc, isfp_loc, ifl, ise, iel, side
    real(kind=DOUBLE), dimension(3) :: norml, normp
    real(kind=DOUBLE) :: Dpc, w, acc
    nsfn = mesh%vert(i)%n_sub_faces_neigh
    allocate(Nm(nsfn, nsfn), bv(nsfn), Pif(nsfn))
    Nm = 0.0_DOUBLE; bv = 0.0_DOUBLE; Pif = 0.0_DOUBLE
    do l = 1, nsfn
      isfl = mesh%vert(i)%sub_face_neigh(l)
      isfl_loc = mesh%sub_face(isfl)%id_loc_around_node
      ifl = mesh%sub_face(isfl)%mesh_face
      do side = 1, 2
        if (side == 1) then
          ise = mesh%sub_face(isfl)%left_sub_elem_neigh
          iel = mesh%face(ifl)%left_neigh
          norml = mesh%sub_face(isfl)%norm
        else
          ise = mesh%sub_face(isfl)%right_sub_elem_neigh
          iel = mesh%face(ifl)%right_neigh
          norml = -mesh%sub_face(isfl)%norm
        end if
        if (ise <= 0 .or. iel <= 0) cycle
        Dpc = 0.0_DOUBLE
        do m = 1, mesh%sub_elem(ise)%n_sub_faces
          Dpc = Dpc + mesh%sub_face(mesh%sub_elem(ise)%sub_face(m))%area
        end do
        Dpc = Dpc/real(mesh%sub_elem(ise)%n_sub_faces, DOUBLE)
        acc = 0.0_DOUBLE
        do m = 1, mesh%sub_elem(ise)%n_sub_faces
          isfp = mesh%sub_elem(ise)%sub_face(m)
          isfp_loc = mesh%sub_face(isfp)%id_loc_around_node
          if (mesh%sub_face(isfp)%left_sub_elem_neigh == ise) then
            normp = mesh%sub_face(isfp)%norm
          else
            normp = -mesh%sub_face(isfp)%norm
          end if
          w = (mesh%sub_face(isfp)%area/Dpc)/lamc(iel)*dot_product(normp, norml)
          Nm(isfl_loc, isfp_loc) = Nm(isfl_loc, isfp_loc) + w
          acc = acc + w
        end do
        bv(isfl_loc) = bv(isfl_loc) + acc*pcell(iel) &
                     + dot_product(vcell(:, iel), norml)
      end do
    end do
    call lu_solve(nsfn, Nm, Pif, bv)
    deallocate(Nm, bv)
  end subroutine pif_solve_synth

  ! CONTROLES EN DUR, avant tout calcul. Valeurs attendues EXACTES, pas des
  ! ordres de grandeur : si l'un echoue il faut s'arreter et pas passer au cas
  ! suivant -- on vient de voir ce que coute une campagne lancee sur une formule
  ! fausse.
  subroutine pif_selfcheck(mesh)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), allocatable :: pc(:), lam(:), vc(:,:), Pif(:)
    real(kind=DOUBLE), dimension(3) :: x
    real(kind=DOUBLE) :: e, emax, pi_cross, xmin, xmax
    integer(kind=ENTIER) :: i, k, nsfn, mcase, nfound

    allocate(pc(mesh%n_elems), lam(mesh%n_elems), vc(3, mesh%n_elems))
    lam = 1.0_DOUBLE

    ! --- 1. pression uniforme + trois champs de vitesse : Pi doit valoir p
    ! EXACTEMENT. Translation, rotation solide et cisaillement lineaire ont tous
    ! (v_c - v_d).n = 0 sur un maillage cartesien, donc Pi = p est solution.
    do mcase = 1, 3
      pc = 1.0_DOUBLE
      do k = 1, mesh%n_elems
        x = mesh%elem(k)%coord
        select case (mcase)
        case (1); vc(:, k) = (/ 0.3_DOUBLE, -0.2_DOUBLE, 0.0_DOUBLE /)
        case (2); vc(:, k) = (/ -x(2), x(1), 0.0_DOUBLE /)
        case (3); vc(:, k) = (/ x(2), 0.0_DOUBLE, 0.0_DOUBLE /)
        end select
      end do
      emax = 0.0_DOUBLE
      do i = 1, mesh%n_vert
        if (mesh%vert(i)%is_bound .or. mesh%vert(i)%is_ghost) cycle
        call pif_solve_synth(mesh, i, pc, vc, lam, Pif, nsfn)
        emax = max(emax, maxval(abs(Pif - 1.0_DOUBLE)))
        deallocate(Pif)
      end do
      select case (mcase)
      case (1); print*, "PIFCHK translation uniforme : max |Pi - p_c| =", emax
      case (2); print*, "PIFCHK rotation solide      : max |Pi - p_c| =", emax
      case (3); print*, "PIFCHK cisaillement lineaire: max |Pi - p_c| =", emax
      end select
    end do

    ! --- 2. Riemann 1D : p = 10 / 1, u = +1 / -1, lam = 1.
    ! Pi attendu sur les sous-faces traversantes = (p_L+p_R)/2 + lam(u_L-u_R)/2
    !                                            = 5.5 + 1 = 6.5
    xmin = 1.0e300_DOUBLE; xmax = -1.0e300_DOUBLE
    do k = 1, mesh%n_elems
      xmin = min(xmin, mesh%elem(k)%coord(1)); xmax = max(xmax, mesh%elem(k)%coord(1))
    end do
    do k = 1, mesh%n_elems
      if (mesh%elem(k)%coord(1) < 0.5_DOUBLE*(xmin + xmax)) then
        pc(k) = 10.0_DOUBLE; vc(:, k) = (/  1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE /)
      else
        pc(k) = 1.0_DOUBLE;  vc(:, k) = (/ -1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE /)
      end if
    end do
    pi_cross = 0.0_DOUBLE; nfound = 0; emax = 0.0_DOUBLE
    do i = 1, mesh%n_vert
      if (mesh%vert(i)%is_bound .or. mesh%vert(i)%is_ghost) cycle
      ! noeud sur l'interface, et seulement lui
      if (abs(mesh%vert(i)%coord(1) - 0.5_DOUBLE*(xmin + xmax)) > 1.0e-12_DOUBLE) cycle
      call pif_solve_synth(mesh, i, pc, vc, lam, Pif, nsfn)
      do k = 1, nsfn
        ! sous-faces TRAVERSANTES : normale selon x
        if (abs(abs(mesh%sub_face(mesh%vert(i)%sub_face_neigh(k))%norm(1)) &
            - 1.0_DOUBLE) > 1.0e-10_DOUBLE) cycle
        pi_cross = Pif(mesh%sub_face(mesh%vert(i)%sub_face_neigh(k))%id_loc_around_node)
        emax = max(emax, abs(pi_cross - 6.5_DOUBLE))
        nfound = nfound + 1
      end do
      deallocate(Pif)
    end do
    print*, "PIFCHK Riemann 1D : Pi traversante =", pi_cross, " ecart a 6.5 =", emax, &
            " sur", nfound, "sous-faces"
    deallocate(pc, lam, vc)
  end subroutine pif_selfcheck

  pure function kap_cross(a, b) result(c)
    implicit none
    real(kind=DOUBLE), dimension(3), intent(in) :: a, b
    real(kind=DOUBLE), dimension(3) :: c
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function kap_cross

  ! ===================================================================
  ! Divergence de SOUS-MAILLE, par Green sur le bord de omega_pc, avec v
  ! interpole depuis les vitesses NODALES seules.
  !
  ! omega_pc est EXACTEMENT le polyedre dont compute_subelem_volume somme les
  ! tetraedres : pour chaque face f de c portant le sommet p, les deux tetraedres
  ! (ce, cf, cpm, cp) et (ce, cf, cp, cpp), ou cf est le centroide de f, ce celui
  ! de c, cpm et cpp les milieux des deux aretes de f en p.
  !
  ! Son bord se reduit a quatre triangles par face :
  !   EXTERIEURS (sur la face de maillage)  (cf, cpm, cp) et (cf, cp, cpp)
  !        -- c'est exactement la decomposition de compute_subface_norm ;
  !   INTERIEURS (dans la maille)           (ce, cf, cpm) et (ce, cf, cpp).
  ! Les autres faces des tetraedres se recouvrent deux a deux A L'INTERIEUR de
  ! omega_pc et disparaissent : (ce, cp, cf) est partagee par les deux tetraedres
  ! d'une meme face, et (ce, cpm, cp) l'est avec le tetraedre de l'AUTRE face
  ! portant l'arete (p, vm). Le bord est donc clos sans rien d'autre, ce que le
  ! controle a champ constant verifie.
  !
  ! Sur un triangle plan et un champ lineaire, la valeur au centroide est la
  ! moyenne de ses trois sommets : l'integrale est EXACTE sur les champs
  ! lineaires. C'est perdu des qu'on met v_c sur les triangles interieurs --
  ! l'erreur devient O(1) -- donc v_c n'apparait nulle part ici.
  ! ===================================================================
  pure subroutine kappa_subcell_div(mesh, vn, id_se, divv, flux_tot)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: vn
    integer(kind=ENTIER), intent(in) :: id_se
    real(kind=DOUBLE), intent(out) :: divv
    real(kind=DOUBLE), dimension(3), intent(out) :: flux_tot

    integer(kind=ENTIER) :: j, k, kp, km, idv, idvp, idvm
    integer(kind=ENTIER) :: id_sf, id_f, id_el, id_vert
    real(kind=DOUBLE), dimension(3) :: ce, cf, cp, cpm, cpp
    real(kind=DOUBLE), dimension(3) :: ve, vf, vpn, vpm, vpp
    real(kind=DOUBLE), dimension(3) :: A1, A2, A3, A4
    real(kind=DOUBLE) :: S

    id_el   = mesh%sub_elem(id_se)%mesh_elem
    id_vert = mesh%sub_elem(id_se)%mesh_vert
    ce = mesh%elem(id_el)%coord

    ! v au centroide de la maille : moyenne des nodales de la maille. Pour un
    ! champ lineaire et ce = moyenne des sommets, c'est v(ce) exactement.
    ve = 0.0_DOUBLE
    do k = 1, mesh%elem(id_el)%n_vert
      ve = ve + vn(:, mesh%elem(id_el)%vert(k))
    end do
    ve = ve / real(mesh%elem(id_el)%n_vert, DOUBLE)

    S = 0.0_DOUBLE
    flux_tot = 0.0_DOUBLE
    do j = 1, mesh%sub_elem(id_se)%n_sub_faces
      id_sf = mesh%sub_elem(id_se)%sub_face(j)
      id_f  = mesh%sub_face(id_sf)%mesh_face

      kp = 1; km = 1
      do k = 1, mesh%face(id_f)%n_vert
        if (mesh%face(id_f)%vert(k) == id_vert) then
          kp = 1 + mod((k - 1) + 1, mesh%face(id_f)%n_vert)
          km = 1 + mod((k - 1) - 1 + mesh%face(id_f)%n_vert, mesh%face(id_f)%n_vert)
          exit
        end if
      end do
      idv  = id_vert
      idvp = mesh%face(id_f)%vert(kp)
      idvm = mesh%face(id_f)%vert(km)

      cf  = mesh%face(id_f)%coord
      cp  = mesh%vert(idv)%coord
      cpm = 0.5_DOUBLE*(mesh%vert(idvm)%coord + cp)
      cpp = 0.5_DOUBLE*(mesh%vert(idvp)%coord + cp)

      vpn = vn(:, idv)
      vpm = 0.5_DOUBLE*(vn(:, idvm) + vpn)
      vpp = 0.5_DOUBLE*(vn(:, idvp) + vpn)
      vf  = 0.0_DOUBLE
      do k = 1, mesh%face(id_f)%n_vert
        vf = vf + vn(:, mesh%face(id_f)%vert(k))
      end do
      vf = vf / real(mesh%face(id_f)%n_vert, DOUBLE)

      ! --- les deux triangles EXTERIEURS, orientation sortante de id_el, la
      ! meme convention que compute_subface_norm et compute_subelem_volume.
      if (mesh%face(id_f)%left_neigh == id_el) then
        A1 = 0.5_DOUBLE*kap_cross(cpm - cf, cp  - cf)
        A2 = 0.5_DOUBLE*kap_cross(cp  - cf, cpp - cf)
      else
        A1 = 0.5_DOUBLE*kap_cross(cp  - cf, cpm - cf)
        A2 = 0.5_DOUBLE*kap_cross(cpp - cf, cp  - cf)
      end if

      ! --- les deux triangles INTERIEURS. L'orientation sortante est determinee
      ! geometriquement -- cp est le quatrieme sommet du tetraedre, donc la
      ! normale sortante lui tourne le dos -- et non par la convention
      ! gauche/droite, qui ne dit rien de l'interieur de la maille.
      A3 = 0.5_DOUBLE*kap_cross(cf - ce, cpm - ce)
      if (dot_product(A3, cp - ce) > 0.0_DOUBLE) A3 = -A3
      A4 = 0.5_DOUBLE*kap_cross(cf - ce, cpp - ce)
      if (dot_product(A4, cp - ce) > 0.0_DOUBLE) A4 = -A4

      S = S + dot_product((vf + vpm + vpn)/3.0_DOUBLE, A1) &
            + dot_product((vf + vpn + vpp)/3.0_DOUBLE, A2) &
            + dot_product((ve + vf + vpm)/3.0_DOUBLE, A3) &
            + dot_product((ve + vf + vpp)/3.0_DOUBLE, A4)
      flux_tot = flux_tot + A1 + A2 + A3 + A4
    end do

    divv = 0.0_DOUBLE
    if (mesh%sub_elem(id_se)%volume > 0.0_DOUBLE) &
      divv = S / mesh%sub_elem(id_se)%volume
  end subroutine kappa_subcell_div

  ! ===================================================================
  ! kappa_pc = max(kap_floor, min(1, h_c |(div v)_pc| / a_c)), pose dans un
  ! tableau de MAILLE a largeur fixe indice par la position du sommet dans
  ! elem%vert, pour passer par l'echange MPI generique.
  !
  ! vn est la vitesse nodale du pas PRECEDENT : kappa est explicite, ce qui est
  ! licite ici et seulement ici, parce qu'il n'entre que comme multiplicateur
  ! scalaire positif devant une somme de carres.
  !
  ! h_c est une distance entre CENTROIDES, pas une racine de volume : sur un
  ! maillage 2D obtenu par extrusion, V^(1/3) porte un facteur dz^(1/3)
  ! entierement arbitraire. Les faces z de l'extrusion sont ecartees de la
  ! moyenne sous boundary_2d pour la meme raison : leur distance vaut dz.
  ! ===================================================================
  subroutine compute_kappa(mesh, sol, vn, gamma_arr, b2d, dt_prev)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: vn
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    logical, intent(in) :: b2d
    real(kind=DOUBLE), intent(in) :: dt_prev

    integer(kind=ENTIER) :: ic, k, j, id_f, id_nb, id_se, nh
    real(kind=DOUBLE) :: a_c, h_c, dv, kk, pc, pjmp, pnb
    real(kind=DOUBLE), dimension(3) :: ftot

    if (.not. allocated(kap_c))  allocate(kap_c(KAP_NV, mesh%n_elems))
    if (pi_mdiv .and. .not. allocated(mdv_c)) then
      allocate(mdv_c(KAP_NV, mesh%n_elems)); mdv_c = 0.0_DOUBLE
    end if
    if (.not. allocated(kap_se)) allocate(kap_se(mesh%n_sub_elems))
    kap_c = 1.0_DOUBLE

    if (kap_check .and. .not. kap_check_done) call kappa_selfcheck(mesh)

    ! AVANT la garde du premier appel : un kappa impose doit valoir des le pas 1,
    ! sans quoi le premier pas serait un pas de vp et la mesure de F ne porterait
    ! pas sur la constante demandee.
    if (kap_fixed > 0.0_DOUBLE) then
      kap_c  = kap_fixed
      kap_se = kap_fixed
      kap_started = .true.
      kap_min = min(kap_min, kap_fixed); kap_max = max(kap_max, kap_fixed)
      return
    end if

    ! Premier appel : aucune vitesse nodale n'existe encore (elle vaut zero au
    ! repos), donc |div v| = 0 et kappa tomberait au plancher PARTOUT, y compris
    ! dans le depot de Sedov. Le premier pas est donc un pas de vp pur, et c'est
    ! dit plutot que subi.
    if (.not. kap_started) then
      kap_started = .true.
      kap_se = 1.0_DOUBLE
      return
    end if
    ! `kap_one` neutralise kappa, mais ce `return` NE DOIT PAS sauter la boucle
    ! quand `vp_mdiv` tourne : c'est elle qui remplit `mdv_c` (plus bas). Il la
    ! sautait, et `mdv_c` -- pourtant ALLOUE -- restait identiquement nul, donc
    ! `dlt` valait 0 a chaque sous-face et la fermeture degenerait EXACTEMENT en
    ! `vp_fbw`. Aucun signal ne le disait : `mdv_n` comptait ses visites
    ! normalement, puisque seul le CONTENU de `mdv_se` manquait, pas son
    ! allocation. kappa lui-meme reste force a 1 en sortie de `kappa_scatter`,
    ! donc le flux n'est pas touche par cette traversee supplementaire.
    if (kap_one .and. .not. pi_mdiv) then
      kap_se = 1.0_DOUBLE
      return
    end if

    do ic = 1, mesh%n_elems
      if (mesh%elem(ic)%volume <= 0.0_DOUBLE) cycle
      pc  = pressure(sol(:, ic), gamma_arr(ic))
      a_c = sqrt(max(0.0_DOUBLE, gamma_arr(ic)*pc*sol(1, ic)))
      h_c = 0.0_DOUBLE; nh = 0
      do j = 1, mesh%elem(ic)%n_faces
        id_f = mesh%elem(ic)%face(j)
        if (b2d .and. abs(mesh%face(id_f)%norm(3)) > 1.0e-8_DOUBLE) cycle
        if (mesh%face(id_f)%left_neigh == ic) then
          id_nb = mesh%face(id_f)%right_neigh
        else
          id_nb = mesh%face(id_f)%left_neigh
        end if
        if (id_nb > 0) then
          h_c = h_c + norm2(mesh%elem(id_nb)%coord - mesh%elem(ic)%coord)
        else
          h_c = h_c + 2.0_DOUBLE*norm2(mesh%face(id_f)%coord - mesh%elem(ic)%coord)
        end if
        nh = nh + 1
      end do
      if (nh > 0) h_c = h_c/real(nh, DOUBLE)

      ! GARDE DURE, pas un min() : KAP_NV = 8 couvre l'hexaedre, donc tous les
      ! types de Sedov et de Gresho, mais PAS un polyedre general -- les mailles
      ! de Voronoi ont couramment plus de huit sommets, et l'etude Voronoi/Gresho
      ! est au programme. Tronquer donnerait soit un kappa au plancher sur les
      ! coins au-dela du huitieme, soit un depassement, et les deux sont
      ! silencieux. On s'arrete bruyamment plutot que de laisser relire ces
      ! conclusions dans six mois sur un maillage polyedrique.
      if (mesh%elem(ic)%n_vert > KAP_NV) then
        print*, "KAP : maille", ic, " a", mesh%elem(ic)%n_vert, &
          " sommets, au-dela de KAP_NV =", KAP_NV, &
          " -- kappa ne peut pas traverser l'echange MPI. Elargir KAP_NV. Arret."
        error stop
      end if

      ! indicateur de FRONT : le plus grand saut relatif de pression aux faces.
      pjmp = 0.0_DOUBLE
      do j = 1, mesh%elem(ic)%n_faces
        id_f = mesh%elem(ic)%face(j)
        if (mesh%face(id_f)%left_neigh == ic) then
          id_nb = mesh%face(id_f)%right_neigh
        else
          id_nb = mesh%face(id_f)%left_neigh
        end if
        if (id_nb <= 0) cycle
        pnb = pressure(sol(:, id_nb), gamma_arr(id_nb))
        pjmp = max(pjmp, abs(pnb - pc)/max(1.0e-300_DOUBLE, max(pc, pnb)))
      end do

      do k = 1, mesh%elem(ic)%n_vert
        id_se = mesh%elem(ic)%sub_elem(k)
        call kappa_subcell_div(mesh, vn, id_se, dv, ftot)
        if (a_c > 0.0_DOUBLE) then
          kk = min(1.0_DOUBLE, h_c*abs(dv)/a_c)
        else
          kk = 1.0_DOUBLE
        end if
        kap_c(k, ic) = max(kap_floor, kk)
        ! `vp_mdiv` : on garde le taux de dilatation BRUT du sous-element, sans le
        ! passer par kappa. Meme traversee, meme echange : rien de neuf a verifier
        ! cote geometrie ni cote MPI.
        if (allocated(mdv_c)) mdv_c(k, ic) = dv
        if (.not. mesh%elem(ic)%is_ghost) then
          kap_min = min(kap_min, kap_c(k, ic))
          kap_max = max(kap_max, kap_c(k, ic))
          kap_sum = kap_sum + kap_c(k, ic)
          kap_n   = kap_n + 1
          if (kap_c(k, ic) <= kap_floor) kap_n_floor = kap_n_floor + 1
          kap_divmax = max(kap_divmax, abs(dv))
          kap_front_jmax = max(kap_front_jmax, pjmp)
          if (pjmp > kap_front_jump) then
            kap_n_front = kap_n_front + 1
            kap_front_min = min(kap_front_min, kap_c(k, ic))
            if (kap_c(k, ic) <= kap_floor) kap_n_front_floor = kap_n_front_floor + 1
          end if
        end if
      end do
    end do

    if (kap_loc) then
      if (trim(kap_loc_mode) == 'face') then
        call compute_kappa_corner_floor(mesh, sol, gamma_arr)
      else if (trim(kap_loc_mode) == 'comp') then
        call compute_kappa_comp_floor(mesh, vn, dt_prev)
      else
        call compute_kappa_node_floor(mesh, sol, gamma_arr)
      end if
    end if
  end subroutine compute_kappa

  ! Plancher local par NOEUD : rapport du terme de pression au terme de vitesse
  ! du second membre nodal, tous deux pris SANS kappa.
  subroutine compute_kappa_node_floor(mesh, sol, gamma_arr)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr

    integer(kind=ENTIER) :: i, j, id_sf, id_f, idl, idr
    real(kind=DOUBLE), dimension(3) :: Rv0, Rp0, nrm
    real(kind=DOUBLE) :: area, pl, pr, vl, vr, lam_l, lam_r, kfl

    if (.not. allocated(kap_node_floor)) allocate(kap_node_floor(mesh%n_vert))
    kap_node_floor = 0.0_DOUBLE

    do i = 1, mesh%n_vert
      Rv0 = 0.0_DOUBLE
      Rp0 = 0.0_DOUBLE
      do j = 1, mesh%vert(i)%n_sub_faces_neigh
        id_sf = mesh%vert(i)%sub_face_neigh(j)
        id_f  = mesh%sub_face(id_sf)%mesh_face
        idl   = mesh%face(id_f)%left_neigh
        idr   = mesh%face(id_f)%right_neigh
        area  = mesh%sub_face(id_sf)%area
        nrm   = mesh%face(id_f)%norm
        pl    = pressure(sol(:, idl), gamma_arr(idl))
        vl    = dot_product(sol(2:4, idl), nrm)
        lam_l = sqrt(max(0.0_DOUBLE, gamma_arr(idl)*pl*sol(1, idl)))/sol(1, idl)
        if (idr > 0) then
          pr    = pressure(sol(:, idr), gamma_arr(idr))
          vr    = dot_product(sol(2:4, idr), nrm)
          lam_r = sqrt(max(0.0_DOUBLE, gamma_arr(idr)*pr*sol(1, idr)))/sol(1, idr)
        else
          ! etat miroir, le meme que celui du solveur nodal : p_r = p_l,
          ! v_r = -v_l, lam_r = lam_l. Le saut de pression y est donc nul, ce qui
          ! est le comportement voulu -- une paroi ne doit pas relever le plancher.
          pr = pl; vr = -vl; lam_r = lam_l
        end if
        Rv0 = Rv0 + area*(lam_l*vl + lam_r*vr)*nrm
        Rp0 = Rp0 + area*(pr - pl)*nrm
      end do
      kfl = min(1.0_DOUBLE, kap_loc_c*norm2(Rp0) &
            / (norm2(Rv0) + 1.0e-300_DOUBLE))
      kap_node_floor(i) = kfl
      if (.not. mesh%vert(i)%is_ghost) then
        kfl_min = min(kfl_min, kfl); kfl_max = max(kfl_max, kfl)
        kfl_sum = kfl_sum + kfl; kfl_n = kfl_n + 1
        if (kfl >= 1.0_DOUBLE) kfl_n_one = kfl_n_one + 1
      end if
    end do
  end subroutine compute_kappa_node_floor

  ! Plancher par COIN : le meme rapport, mais evalue sous-face par sous-face en
  ! valeurs absolues, donc sans compensation, et le coin prend le max.
  !
  !   plancher_pcf = min(1, C |p_r - p_l| / (lam_l |u_l| + lam_r |u_r| + eps))
  !
  ! IMPORTANT, le denominateur porte sur la NORME de la vitesse et non sur sa
  ! composante normale. Avec u.n, une sous-face dont la normale est orthogonale
  ! a l'ecoulement aurait un denominateur nul et un plancher a 1 -- or c'est
  ! exactement le cas de la plupart des sous-faces d'un tourbillon de Gresho,
  ! ou la vitesse est tangentielle. Le plancher y monterait a 1 partout et
  ! detruirait le tourbillon, c'est-a-dire l'inverse du but.
  !
  ! Entierement LOCAL a la sous-face : il ne demande que les deux mailles
  ! adjacentes, dont l'etat est echange. Il est donc calculable a l'identique
  ! pour une maille fantome, contrairement a la divergence de sous-maille qui
  ! exigeait tous les sommets de la maille. Aucun echange MPI n'est necessaire.
  subroutine compute_kappa_corner_floor(mesh, sol, gamma_arr)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr

    integer(kind=ENTIER) :: k, j, id_sf, id_f, idl, idr, id_el
    real(kind=DOUBLE), dimension(3) :: nrm
    real(kind=DOUBLE) :: pl, pr, lam_l, lam_r, ul, ur, den, kfl, kf

    if (.not. allocated(kap_se_floor)) allocate(kap_se_floor(mesh%n_sub_elems))
    kap_se_floor = 0.0_DOUBLE

    do k = 1, mesh%n_sub_elems
      id_el = mesh%sub_elem(k)%mesh_elem
      kf = 0.0_DOUBLE
      do j = 1, mesh%sub_elem(k)%n_sub_faces
        id_sf = mesh%sub_elem(k)%sub_face(j)
        id_f  = mesh%sub_face(id_sf)%mesh_face
        idl   = mesh%face(id_f)%left_neigh
        idr   = mesh%face(id_f)%right_neigh
        nrm   = mesh%face(id_f)%norm
        pl    = pressure(sol(:, idl), gamma_arr(idl))
        ul    = norm2(sol(2:4, idl))
        lam_l = sqrt(max(0.0_DOUBLE, gamma_arr(idl)*pl*sol(1, idl)))/sol(1, idl)
        if (idr > 0) then
          pr    = pressure(sol(:, idr), gamma_arr(idr))
          ur    = norm2(sol(2:4, idr))
          lam_r = sqrt(max(0.0_DOUBLE, gamma_arr(idr)*pr*sol(1, idr)))/sol(1, idr)
        else
          ! etat miroir : p_r = p_l, donc saut nul, donc une paroi ne releve
          ! pas le plancher -- le meme choix que dans la version par noeud.
          pr = pl; ur = ul; lam_r = lam_l
        end if
        den = lam_l*ul + lam_r*ur
        kfl = min(1.0_DOUBLE, kap_loc_c*abs(pr - pl)/(den + 1.0e-300_DOUBLE))
        kf = max(kf, kfl)
      end do
      kap_se_floor(k) = kf
      if (.not. mesh%elem(id_el)%is_ghost) then
        kfl_min = min(kfl_min, kf); kfl_max = max(kfl_max, kf)
        kfl_sum = kfl_sum + kf; kfl_n = kfl_n + 1
        if (kf >= 1.0_DOUBLE) kfl_n_one = kfl_n_one + 1
      end if
    end do
  end subroutine compute_kappa_corner_floor

  ! Plancher sur le TAUX DE COMPRESSION de la maille :
  !
  !   plancher_c = min(1, C |Vpoint_c| dt / V_c)
  !
  ! |Vpoint_c| dt / V_c est la fraction de volume que la maille va perdre au pas
  ! courant. Relever kappa la ou elle depasse quelques pour cent rend son
  ! amortissement au mode nodal EXACTEMENT la ou une maille est en train d'etre
  ! ecrasee, et nulle part ailleurs -- ce que ni la forme par noeud ni la forme
  ! par sous-face ne font, puisqu'elles bornent |v_p| et non sa divergence
  ! autour d'une maille. C'est une propriete du COIN et pas du noeud : des
  ! noeuds qui convergent les uns vers les autres ont chacun une vitesse modeste.
  !
  ! Sur un ecoulement incompressible, Vpoint_c = 0 et le plancher ne mord NULLE
  ! PART : le tourbillon est alors intact par construction et non par un
  ! argument d'ordre de grandeur.
  !
  ! Vpoint_c = somme_f A (v_p . n), sur les sous-faces de la maille, avec le
  ! signe de l'orientation. Pose dans kap_c AVANT l'echange MPI, donc il le
  ! traverse sans mecanisme supplementaire.
  subroutine compute_kappa_comp_floor(mesh, vn, dtc)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: vn
    real(kind=DOUBLE), intent(in) :: dtc

    integer(kind=ENTIER) :: i, j, id_sf, id_f, idl, idr, k
    real(kind=DOUBLE) :: area, v_et, kfl

    if (.not. allocated(kap_cell_floor)) allocate(kap_cell_floor(mesh%n_elems))
    kap_cell_floor = 0.0_DOUBLE

    do i = 1, mesh%n_vert
      do j = 1, mesh%vert(i)%n_sub_faces_neigh
        id_sf = mesh%vert(i)%sub_face_neigh(j)
        id_f  = mesh%sub_face(id_sf)%mesh_face
        idl   = mesh%face(id_f)%left_neigh
        idr   = mesh%face(id_f)%right_neigh
        area  = mesh%sub_face(id_sf)%area
        v_et  = dot_product(vn(:, i), mesh%face(id_f)%norm)
        kap_cell_floor(idl) = kap_cell_floor(idl) + area*v_et
        if (idr > 0) kap_cell_floor(idr) = kap_cell_floor(idr) - area*v_et
      end do
    end do

    do i = 1, mesh%n_elems
      if (mesh%elem(i)%volume > 0.0_DOUBLE) then
        kfl = min(1.0_DOUBLE, kap_loc_c*abs(kap_cell_floor(i))*dtc &
              / mesh%elem(i)%volume)
      else
        kfl = 1.0_DOUBLE
      end if
      kap_cell_floor(i) = kfl
      do k = 1, mesh%elem(i)%n_vert
        ! Le comptage se fait ICI et pas dans kappa_scatter : en mode 'comp' le
        ! plancher est pose directement dans kap_c, donc le compteur de
        ! kappa_scatter -- qui ne couvre que 'face' et 'node' -- valait
        ! structurellement zero et n'etait pas une mesure.
        if (kfl > kap_c(k, i)) then
          kap_c(k, i) = kfl
          if (.not. mesh%elem(i)%is_ghost) keff_n_raised = keff_n_raised + 1
        end if
      end do
      if (.not. mesh%elem(i)%is_ghost) then
        kfl_min = min(kfl_min, kfl); kfl_max = max(kfl_max, kfl)
        kfl_sum = kfl_sum + kfl; kfl_n = kfl_n + 1
        if (kfl >= 1.0_DOUBLE) kfl_n_one = kfl_n_one + 1
      end if
    end do
  end subroutine compute_kappa_comp_floor

  ! kap_c (par maille, indice par la position du sommet) -> kap_se (par
  ! sous-element), apres l'echange MPI. Le hot loop lit kap_se et ne cherche
  ! donc rien.
  subroutine kappa_scatter(mesh)
    implicit none
    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER) :: ic, k, id_se
    if (.not. allocated(kap_se)) allocate(kap_se(mesh%n_sub_elems))
    if (pi_mdiv .and. .not. allocated(mdv_se)) then
      allocate(mdv_se(mesh%n_sub_elems)); mdv_se = 0.0_DOUBLE
    end if
    kap_se = 1.0_DOUBLE
    if (kap_fixed > 0.0_DOUBLE) then
      kap_se = kap_fixed
      return
    end if
    ! Meme defaut qu'en amont, et c'est le second maillon : ce `return` sautait
    ! le remplissage de `mdv_se` juste en dessous. Sous `vp_mdiv` on traverse
    ! donc, et kappa est force a 1 APRES la boucle (voir la fin de la routine).
    ! `.not. kap_started` reste bloquant dans tous les cas : au premier pas il
    ! n'y a pas de vitesse nodale du pas precedent, donc pas de dilatation
    ! definie -- `mdv_se` vaut alors 0, ce qui est voulu et non subi.
    if (.not. kap_started) return
    if (kap_one .and. .not. pi_mdiv) return
    do ic = 1, mesh%n_elems
      do k = 1, mesh%elem(ic)%n_vert
        id_se = mesh%elem(ic)%sub_elem(k)
        kap_se(id_se) = kap_c(k, ic)
        if (allocated(mdv_se) .and. allocated(mdv_c)) mdv_se(id_se) = mdv_c(k, ic)
        ! Le plancher local du NOEUD releve kappa sur tous ses coins. Il
        ! s'applique apres l'echange MPI de kap_c, donc il ne le perturbe pas.
        if (kap_loc) then
          if (trim(kap_loc_mode) == 'face') then
            if (allocated(kap_se_floor)) then
              if (kap_se_floor(id_se) > kap_se(id_se)) then
                kap_se(id_se) = kap_se_floor(id_se)
                if (.not. mesh%elem(ic)%is_ghost) keff_n_raised = keff_n_raised + 1
              end if
            end if
          else if (allocated(kap_node_floor)) then
            if (kap_node_floor(mesh%sub_elem(id_se)%mesh_vert) > kap_se(id_se)) then
              kap_se(id_se) = kap_node_floor(mesh%sub_elem(id_se)%mesh_vert)
              if (.not. mesh%elem(ic)%is_ghost) keff_n_raised = keff_n_raised + 1
            end if
          end if
        end if
        if (.not. mesh%elem(ic)%is_ghost) then
          keff_min = min(keff_min, kap_se(id_se))
          keff_max = max(keff_max, kap_se(id_se))
          keff_sum = keff_sum + kap_se(id_se)
          keff_n = keff_n + 1
        end if
      end do
    end do
    ! `kap_one` : kappa doit rester STRICTEMENT neutre. On ne peut plus le faire
    ! par un `return` en amont (il emportait `mdv_se`), donc on l'impose ici, en
    ! sortie. `kap_se` n'est lu qu'en deux points du flux (kl/kr), donc cette
    ! remise a 1 suffit a garantir que la traversee supplementaire faite pour
    ! `vp_mdiv` ne change rien au schema.
    if (kap_one) kap_se = 1.0_DOUBLE
  end subroutine kappa_scatter

  ! Controle d'exactitude de l'operateur de sous-maille, une fois par run.
  ! Deux champs NODAUX synthetiques, poses depuis les coordonnees :
  !   constant   -> divergence nulle ET bord clos (somme A n = 0) ;
  !   lineaire   -> divergence connue a+e+i, dans CHAQUE sous-maille.
  ! On imprime les pires ecarts. Un echec ici invalide tout le reste, donc il
  ! est mesure avant le premier pas et non suppose.
  subroutine kappa_selfcheck(mesh)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), allocatable :: vcst(:,:), vlin(:,:)
    real(kind=DOUBLE), dimension(3) :: ftot, x
    real(kind=DOUBLE) :: dv, e_cst, e_lin, e_clos, dtrue, vmin, vmax
    integer(kind=ENTIER) :: i, k
    real(kind=DOUBLE), parameter :: AA = 0.7_DOUBLE, BB = -0.4_DOUBLE, CC = 0.25_DOUBLE
    real(kind=DOUBLE), parameter :: DD = 0.15_DOUBLE, EE = -0.3_DOUBLE, FF = 0.6_DOUBLE
    real(kind=DOUBLE), parameter :: GG = -0.2_DOUBLE, HH = 0.5_DOUBLE, II = 0.2_DOUBLE

    kap_check_done = .true.
    allocate(vcst(3, mesh%n_vert), vlin(3, mesh%n_vert))
    do i = 1, mesh%n_vert
      x = mesh%vert(i)%coord
      vcst(:, i) = (/ 1.3_DOUBLE, -0.7_DOUBLE, 0.4_DOUBLE /)
      vlin(1, i) = AA*x(1) + BB*x(2) + CC*x(3)
      vlin(2, i) = DD*x(1) + EE*x(2) + FF*x(3)
      vlin(3, i) = GG*x(1) + HH*x(2) + II*x(3)
    end do
    dtrue = AA + EE + II
    e_cst = 0.0_DOUBLE; e_lin = 0.0_DOUBLE; e_clos = 0.0_DOUBLE
    vmin = 1.0e300_DOUBLE; vmax = -1.0e300_DOUBLE
    do k = 1, mesh%n_sub_elems
      if (mesh%sub_elem(k)%volume <= 0.0_DOUBLE) cycle
      call kappa_subcell_div(mesh, vcst, k, dv, ftot)
      e_cst  = max(e_cst, abs(dv))
      e_clos = max(e_clos, maxval(abs(ftot)))
      call kappa_subcell_div(mesh, vlin, k, dv, ftot)
      e_lin = max(e_lin, abs(dv - dtrue))
      vmin = min(vmin, mesh%sub_elem(k)%volume)
      vmax = max(vmax, mesh%sub_elem(k)%volume)
    end do
    print*, "KAPCHK sous-elements :", mesh%n_sub_elems, " V dans [", vmin, ",", vmax, "]"
    print*, "KAPCHK bord clos    max|somme A n|      =", e_clos
    print*, "KAPCHK champ constant max|div|          =", e_cst
    print*, "KAPCHK champ lineaire max|div - (a+e+i)| =", e_lin, " pour div =", dtrue
    deallocate(vcst, vlin)
  end subroutine kappa_selfcheck

  ! Coefficient GEOMETRIQUE de v_p dans (div v)_pc. En reprenant la
  ! decomposition de kappa_subcell_div, les quatre triangles par face portent v_p
  ! avec les poids suivants : v_p lui-meme (1), les milieux d'arete (1/2), le
  ! centroide de face (1/n_f), celui de la maille (1/n_e). D'ou
  !   T1, T2 : (1/n_f + 1/2 + 1)/3        sur les triangles EXTERIEURS
  !   T3, T4 : (1/n_e + 1/n_f + 1/2)/3    sur les triangles INTERIEURS
  ! et b_pc = somme des A_sigma ponderes. Purement geometrique, donc calculable
  ! une fois par pas et independamment de l'ecoulement.
  pure subroutine subcell_div_b(mesh, id_se, bvec)
    implicit none
    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_se
    real(kind=DOUBLE), dimension(3), intent(out) :: bvec

    integer(kind=ENTIER) :: j, k, kp, km, idvp, idvm
    integer(kind=ENTIER) :: id_sf, id_f, id_el, id_vert, nf, ne
    real(kind=DOUBLE), dimension(3) :: ce, cf, cp, cpm, cpp, A1, A2, A3, A4
    real(kind=DOUBLE) :: c12, c34

    id_el   = mesh%sub_elem(id_se)%mesh_elem
    id_vert = mesh%sub_elem(id_se)%mesh_vert
    ce = mesh%elem(id_el)%coord
    ne = mesh%elem(id_el)%n_vert
    bvec = 0.0_DOUBLE

    do j = 1, mesh%sub_elem(id_se)%n_sub_faces
      id_sf = mesh%sub_elem(id_se)%sub_face(j)
      id_f  = mesh%sub_face(id_sf)%mesh_face
      nf    = mesh%face(id_f)%n_vert
      kp = 1; km = 1
      do k = 1, nf
        if (mesh%face(id_f)%vert(k) == id_vert) then
          kp = 1 + mod((k - 1) + 1, nf)
          km = 1 + mod((k - 1) - 1 + nf, nf)
          exit
        end if
      end do
      idvp = mesh%face(id_f)%vert(kp)
      idvm = mesh%face(id_f)%vert(km)
      cf  = mesh%face(id_f)%coord
      cp  = mesh%vert(id_vert)%coord
      cpm = 0.5_DOUBLE*(mesh%vert(idvm)%coord + cp)
      cpp = 0.5_DOUBLE*(mesh%vert(idvp)%coord + cp)

      if (mesh%face(id_f)%left_neigh == id_el) then
        A1 = 0.5_DOUBLE*kap_cross(cpm - cf, cp  - cf)
        A2 = 0.5_DOUBLE*kap_cross(cp  - cf, cpp - cf)
      else
        A1 = 0.5_DOUBLE*kap_cross(cp  - cf, cpm - cf)
        A2 = 0.5_DOUBLE*kap_cross(cpp - cf, cp  - cf)
      end if
      A3 = 0.5_DOUBLE*kap_cross(cf - ce, cpm - ce)
      if (dot_product(A3, cp - ce) > 0.0_DOUBLE) A3 = -A3
      A4 = 0.5_DOUBLE*kap_cross(cf - ce, cpp - ce)
      if (dot_product(A4, cp - ce) > 0.0_DOUBLE) A4 = -A4

      c12 = (1.0_DOUBLE/real(nf, DOUBLE) + 1.5_DOUBLE)/3.0_DOUBLE
      c34 = (1.0_DOUBLE/real(ne, DOUBLE) + 1.0_DOUBLE/real(nf, DOUBLE) &
            + 0.5_DOUBLE)/3.0_DOUBLE
      bvec = bvec + c12*(A1 + A2) + c34*(A3 + A4)
    end do
  end subroutine subcell_div_b

  ! Prepare b_pc, h_c et reste_pc. reste_pc est pris sur la vitesse nodale du
  ! PAS PRECEDENT : c'est lui qui decouple le systeme nodal, et c'est licite
  ! parce que seul b_pc.v_p a besoin d'etre implicite pour que la matrice reste
  ! locale. Tout est empaquete dans des tableaux de MAILLE a largeur fixe pour
  ! traverser l'echange MPI generique -- une maille fantome a besoin de tous ses
  ! sommets pour sa divergence, et certains sont hors du halo.
  subroutine divdirect_prepare(mesh, vn, b2d, sol, gamma_arr)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: vn
    logical, intent(in) :: b2d
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr

    integer(kind=ENTIER) :: ic, k, j, id_f, id_nb, id_se, nh, jj, id_sf2
    real(kind=DOUBLE) :: h_c, dv, sa, pcl, lam, qq
    real(kind=DOUBLE), dimension(3) :: ftot, bb

    if (.not. allocated(divd_b))     allocate(divd_b(3, mesh%n_sub_elems))
    if (.not. allocated(divd_rest))  allocate(divd_rest(mesh%n_sub_elems))
    if (.not. allocated(divd_h))     allocate(divd_h(mesh%n_elems))
    if (.not. allocated(divd_bc))    allocate(divd_bc(3*KAP_NV, mesh%n_elems))
    if (.not. allocated(divd_restc)) allocate(divd_restc(KAP_NV, mesh%n_elems))
    if (.not. allocated(divd_isum))  allocate(divd_isum(mesh%n_sub_elems))
    if (.not. allocated(divd_psi))   allocate(divd_psi(mesh%n_sub_elems))
    divd_psi = 1.0_DOUBLE
    if (.not. allocated(divd_isumc)) allocate(divd_isumc(KAP_NV, mesh%n_elems))
    if (.not. allocated(divd_prod_cell)) allocate(divd_prod_cell(mesh%n_elems))
    divd_bc = 0.0_DOUBLE; divd_restc = 0.0_DOUBLE; divd_isumc = 0.0_DOUBLE

    do ic = 1, mesh%n_elems
      if (mesh%elem(ic)%n_vert > KAP_NV) then
        print*, "VPDIV : maille", ic, " a", mesh%elem(ic)%n_vert, &
          " sommets, au-dela de KAP_NV =", KAP_NV, " -- elargir. Arret."
        error stop
      end if
      ! h_c : distance entre CENTROIDES, faces z de l'extrusion ecartees sous
      ! boundary_2d, pour la raison deja etablie -- leur distance vaut dz, qui
      ! est arbitraire.
      h_c = 0.0_DOUBLE; nh = 0
      do j = 1, mesh%elem(ic)%n_faces
        id_f = mesh%elem(ic)%face(j)
        if (b2d .and. abs(mesh%face(id_f)%norm(3)) > 1.0e-8_DOUBLE) cycle
        if (mesh%face(id_f)%left_neigh == ic) then
          id_nb = mesh%face(id_f)%right_neigh
        else
          id_nb = mesh%face(id_f)%left_neigh
        end if
        if (id_nb > 0) then
          h_c = h_c + norm2(mesh%elem(id_nb)%coord - mesh%elem(ic)%coord)
        else
          h_c = h_c + 2.0_DOUBLE*norm2(mesh%face(id_f)%coord - mesh%elem(ic)%coord)
        end if
        nh = nh + 1
      end do
      if (nh > 0) h_c = h_c/real(nh, DOUBLE)
      divd_h(ic) = h_c

      do k = 1, mesh%elem(ic)%n_vert
        id_se = mesh%elem(ic)%sub_elem(k)
        call subcell_div_b(mesh, id_se, bb)
        call kappa_subcell_div(mesh, vn, id_se, dv, ftot)
        divd_b(:, id_se) = bb
        ! S = dv * |omega_pc| ; reste = S - b.v_p
        divd_rest(id_se) = dv*mesh%sub_elem(id_se)%volume &
          - dot_product(bb, vn(:, mesh%sub_elem(id_se)%mesh_vert))
        ! somme des aires des sous-faces du COIN, faces z exclues sous b2d
        sa = 0.0_DOUBLE
        do jj = 1, mesh%sub_elem(id_se)%n_sub_faces
          id_sf2 = mesh%sub_elem(id_se)%sub_face(jj)
          if (b2d .and. abs(mesh%face(mesh%sub_face(id_sf2)%mesh_face)%norm(3)) &
              > 1.0e-8_DOUBLE) cycle
          sa = sa + mesh%sub_face(id_sf2)%area
        end do
        if (sa > 0.0_DOUBLE) then
          divd_isum(id_se) = 1.0_DOUBLE/sa
        else
          divd_isum(id_se) = 0.0_DOUBLE
        end if
        divd_bc(3*(k-1)+1:3*(k-1)+3, ic) = bb
        divd_restc(k, ic) = divd_rest(id_se)
        divd_isumc(k, ic) = divd_isum(id_se)
      end do
    end do

    ! psi(q) a l'itere precedent, UNE valeur par coin, pour que la matrice et
    ! le flux utilisent exactement le meme coefficient -- c'est cette identite
    ! qui rend le bilan nodal exact.
    if (divd_nl) then
      do ic = 1, mesh%n_elems
        pcl = pressure(sol(:, ic), gamma_arr(ic))
        if (pcl <= 0.0_DOUBLE) cycle
        lam = sqrt(max(0.0_DOUBLE, gamma_arr(ic)*pcl*sol(1, ic)))/sol(1, ic)
        do k = 1, mesh%elem(ic)%n_vert
          id_se = mesh%elem(ic)%sub_elem(k)
          qq = lam*divd_isum(id_se) &
             * (dot_product(divd_b(:, id_se), vn(:, mesh%sub_elem(id_se)%mesh_vert)) &
                + divd_rest(id_se)) / pcl
          divd_psi(id_se) = divd_psi_of(qq, gamma_arr(ic))
        end do
      end do
    end if

    if (.not. divd_check_done) then
      divd_check_done = .true.
      call divdirect_selfcheck(mesh)
    end if
  end subroutine divdirect_prepare

  ! psi(q) = (1 - phi(q))/q, avec phi(q) = (1 - (gamma-1)q/(2 gamma))^(2 gamma/(gamma-1)).
  ! Limite en zero : psi(0) = 1, posee explicitement et non laissee a un 0/0.
  ! Au-dela de q = 2 gamma/(gamma-1) la base s'annule : phi = 0, psi = 1/q, et la
  ! sous-face est en formation de VIDE -- comptee, jamais masquee.
  function divd_psi_of(q, gam) result(ps)
    implicit none
    real(kind=DOUBLE), intent(in) :: q, gam
    real(kind=DOUBLE) :: ps, base, phi
    if (abs(q) < 1.0e-12_DOUBLE) then
      ps = 1.0_DOUBLE
      return
    end if
    base = 1.0_DOUBLE - (gam - 1.0_DOUBLE)/(2.0_DOUBLE*gam)*q
    if (base <= 0.0_DOUBLE) then
      divd_n_vide = divd_n_vide + 1
      phi = 0.0_DOUBLE
    else
      phi = base**(2.0_DOUBLE*gam/(gam - 1.0_DOUBLE))
    end if
    ps = (1.0_DOUBLE - phi)/q
  end function divd_psi_of

  ! Redistribution apres l'echange MPI : les tableaux de maille reviennent vers
  ! les tableaux de sous-element que la boucle chaude lit sans chercher.
  subroutine divdirect_scatter(mesh)
    implicit none
    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER) :: ic, k, id_se
    do ic = 1, mesh%n_elems
      do k = 1, mesh%elem(ic)%n_vert
        id_se = mesh%elem(ic)%sub_elem(k)
        divd_b(:, id_se)    = divd_bc(3*(k-1)+1:3*(k-1)+3, ic)
        divd_rest(id_se)    = divd_restc(k, ic)
        divd_isum(id_se)    = divd_isumc(k, ic)
      end do
    end do
  end subroutine divdirect_scatter

  ! Controles a VALEUR ATTENDUE EXACTE, une fois par run, sur le vrai maillage.
  !   1. kappa_pc/(lam_c h_c) doit valoir 0.25 sur un quadrangle cartesien. Un
  !      ecart signale que les faces de normale +/-e_z sont entrees dans la somme
  !      de surfaces, donc que le coefficient depend de dz -- le parametre
  !      arbitraire dont toute la construction cherche a sortir.
  !   2. Sur un CISAILLEMENT PUR et sur une ROTATION RIGIDE, la production doit
  !      etre a zero machine ET sans aucun coin negatif. C'est le controle le
  !      plus discriminant : il echouerait pour n'importe quelle erreur dans la
  !      decomposition de d(omega_pc) ou dans les poids des sommets.
  !   3. En COMPRESSION UNIFORME elle doit etre strictement positive, sans coin
  !      negatif non plus.
  ! Sur un champ a divergence nulle mais non trivial (Taylor-Green), ~50 % des
  ! coins sont negatifs par simple symetrie de l'erreur : ce n'est pas un defaut
  ! et seule la SOMME compte.
  subroutine divdirect_selfcheck(mesh)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), allocatable :: vn(:,:)
    real(kind=DOUBLE), dimension(3) :: ftot, x, vc, Npc
    real(kind=DOUBLE) :: dv, g, gmin, gmax, som, pr, xx, yy
    integer(kind=ENTIER) :: k, i, ic, j, id_sf, id_f, nneg, ntot, m

    gmin = 1.0e300_DOUBLE; gmax = -1.0e300_DOUBLE
    do k = 1, mesh%n_sub_elems
      ic = mesh%sub_elem(k)%mesh_elem
      if (divd_h(ic) <= 0.0_DOUBLE) cycle
      g = mesh%sub_elem(k)%volume*divd_isum(k)/divd_h(ic)
      gmin = min(gmin, g); gmax = max(gmax, g)
    end do
    print*, "VPDIVCHK kappa/(lam h) dans [", gmin, ",", gmax, "]  attendu 0.25 (cartesien)"

    allocate(vn(3, mesh%n_vert))
    do m = 1, 3
      do i = 1, mesh%n_vert
        x = mesh%vert(i)%coord
        select case (m)
        case (1); vn(:, i) = (/  x(2), 0.0_DOUBLE, 0.0_DOUBLE /)   ! cisaillement pur
        case (2); vn(:, i) = (/ -x(2), x(1),       0.0_DOUBLE /)   ! rotation rigide
        case (3); vn(:, i) = (/ -x(1), -x(2),      0.0_DOUBLE /)   ! compression uniforme
        end select
      end do
      som = 0.0_DOUBLE; nneg = 0; ntot = 0
      do k = 1, mesh%n_sub_elems
        ic = mesh%sub_elem(k)%mesh_elem
        if (mesh%elem(ic)%is_ghost) cycle
        if (mesh%sub_elem(k)%volume <= 0.0_DOUBLE) cycle
        call kappa_subcell_div(mesh, vn, k, dv, ftot)
        ! N_pc = somme des A n sortants du coin, sur ses sous-faces
        Npc = 0.0_DOUBLE
        do j = 1, mesh%sub_elem(k)%n_sub_faces
          id_sf = mesh%sub_elem(k)%sub_face(j)
          id_f  = mesh%sub_face(id_sf)%mesh_face
          if (mesh%face(id_f)%left_neigh == ic) then
            Npc = Npc + mesh%sub_face(id_sf)%area*mesh%face(id_f)%norm
          else
            Npc = Npc - mesh%sub_face(id_sf)%area*mesh%face(id_f)%norm
          end if
        end do
        xx = mesh%elem(ic)%coord(1); yy = mesh%elem(ic)%coord(2)
        select case (m)
        case (1); vc = (/  yy, 0.0_DOUBLE, 0.0_DOUBLE /)
        case (2); vc = (/ -yy, xx,         0.0_DOUBLE /)
        case (3); vc = (/ -xx, -yy,        0.0_DOUBLE /)
        end select
        ! production du coin : kappa_pc (div v)_pc . G_pc, avec lam = 1
        pr = mesh%sub_elem(k)%volume*divd_isum(k)*dv &
           * dot_product(vn(:, mesh%sub_elem(k)%mesh_vert) - vc, Npc)
        som = som + pr; ntot = ntot + 1
        if (pr < 0.0_DOUBLE) nneg = nneg + 1
      end do
      select case (m)
      case (1); print*, "VPDIVCHK cisaillement pur   : somme", som, " coins < 0 :", nneg, "/", ntot
      case (2); print*, "VPDIVCHK rotation rigide    : somme", som, " coins < 0 :", nneg, "/", ntot
      case (3); print*, "VPDIVCHK compression unif.  : somme", som, " coins < 0 :", nneg, "/", ntot
      end select
    end do
    deallocate(vn)

    ! --- TROIS CONTROLES qui separent un defaut de FORMULE d'un defaut de CODE.
    ! Si les trois passent, le code implemente bien la formule et c'est la
    ! formule qu'il faut changer ; si l'un echoue, c'est un bug.
    block
      real(kind=DOUBLE), dimension(3) :: bb2, Npc, sN
      real(kind=DOUBLE) :: e1, e2, e3, nb, nn, sn1, rmin, rmax
      integer(kind=ENTIER) :: kk, jj, id_sf2, id_f2, ice, iv
      real(kind=DOUBLE), allocatable :: Rp_u(:,:)

      ! 1. b_pc = (3/8) N_pc ?
      e1 = 0.0_DOUBLE; sn1 = 0.0_DOUBLE
      rmin = 1.0e300_DOUBLE; rmax = -1.0e300_DOUBLE
      do kk = 1, mesh%n_sub_elems
        ice = mesh%sub_elem(kk)%mesh_elem
        call subcell_div_b(mesh, kk, bb2)
        Npc = 0.0_DOUBLE
        do jj = 1, mesh%sub_elem(kk)%n_sub_faces
          id_sf2 = mesh%sub_elem(kk)%sub_face(jj)
          id_f2  = mesh%sub_face(id_sf2)%mesh_face
          if (mesh%face(id_f2)%left_neigh == ice) then
            Npc = Npc + mesh%sub_face(id_sf2)%area*mesh%face(id_f2)%norm
          else
            Npc = Npc - mesh%sub_face(id_sf2)%area*mesh%face(id_f2)%norm
          end if
        end do
        e1 = max(e1, norm2(bb2 - 0.375_DOUBLE*Npc))
        nb = norm2(bb2); nn = norm2(Npc)
        if (nb > 0.0_DOUBLE .and. nn > 0.0_DOUBLE) then
          sn1 = max(sn1, norm2(kap_cross(bb2/nb, Npc/nn)))
          rmin = min(rmin, nb/nn); rmax = max(rmax, nb/nn)
        end if
      end do
      print*, "VPDIVV1 max |b_pc - 0.375 N_pc| =", e1, "  max sin(angle) =", sn1
      print*, "VPDIVV1 rapport |b|/|N| dans [", rmin, ",", rmax, "]"

      ! 2. somme_c N_pc = 0 a chaque noeud ?
      e2 = 0.0_DOUBLE
      do iv = 1, mesh%n_vert
        if (mesh%vert(iv)%is_bound) cycle
        sN = 0.0_DOUBLE
        do jj = 1, mesh%vert(iv)%n_sub_faces_neigh
          id_sf2 = mesh%vert(iv)%sub_face_neigh(jj)
          id_f2  = mesh%sub_face(id_sf2)%mesh_face
          if (abs(mesh%face(id_f2)%norm(3)) > 1.0e-8_DOUBLE) cycle
          sN = sN + mesh%sub_face(id_sf2)%area*mesh%face(id_f2)%norm
          if (mesh%face(id_f2)%right_neigh > 0) &
            sN = sN - mesh%sub_face(id_sf2)%area*mesh%face(id_f2)%norm
        end do
        e2 = max(e2, maxval(abs(sN)))
      end do
      print*, "VPDIVV2 max |somme_c N_pc| aux noeuds INTERIEURS =", e2

      ! 3. a pression uniforme, le second membre nodal doit etre EXACTEMENT nul.
      ! R_p = somme_f A [(p_c - p_d) - c reste_c + c reste_d] n ; avec p uniforme
      ! et vn = 0 (donc reste = 0), il se reduit a somme_f A (p-p) n = 0 terme a
      ! terme. On le verifie sur la geometrie reelle plutot que de l'admettre.
      allocate(Rp_u(3, mesh%n_vert)); Rp_u = 0.0_DOUBLE
      do iv = 1, mesh%n_vert
        do jj = 1, mesh%vert(iv)%n_sub_faces_neigh
          id_sf2 = mesh%vert(iv)%sub_face_neigh(jj)
          id_f2  = mesh%sub_face(id_sf2)%mesh_face
          if (mesh%face(id_f2)%right_neigh > 0) cycle   ! interieur : p_c - p_d = 0
          if (abs(mesh%face(id_f2)%norm(3)) > 1.0e-8_DOUBLE) cycle   ! face z : miroir
          Rp_u(:, iv) = Rp_u(:, iv) &
            + mesh%sub_face(id_sf2)%area*mesh%face(id_f2)%norm
        end do
      end do
      e3 = 0.0_DOUBLE
      do iv = 1, mesh%n_vert
        if (mesh%vert(iv)%is_bound) cycle
        e3 = max(e3, maxval(abs(Rp_u(:, iv))))
      end do
      print*, "VPDIVV3 a pression uniforme, max |R_p| aux noeuds INTERIEURS =", e3
      deallocate(Rp_u)
    end block
  end subroutine divdirect_selfcheck

  ! ===================================================================
  ! 'vp_div' : assemblage complet. UN SEUL parcours de noeuds, le decouplage
  ! etant assure par le fait que seul b_pc.v_p est implicite.
  !
  !   Pi_pcf = p_c - c_pc (b_pc.v_p + reste_pc)      c_pc = lam_c / somme A_pcf
  !
  ! Bilan nodal somme_f A (Pi_pcf - Pi_pdf) n = 0, d'ou
  !   M_p = somme_f A n (x) [c_pc b_pc - c_pd b_pd]
  !   R_p = somme_f A [ (p_c - p_d) - c_pc reste_pc + c_pd reste_pd ] n
  ! M_p est locale d x d, non symetrique, et generiquement inversible.
  !
  ! CONTROLE DE SIGNE indispensable : a l'etat uniforme au repos, p_c = p_d et
  ! les deux restes sont nuls, donc R_p doit etre EXACTEMENT nul. C'est le seul
  ! endroit ou une erreur de signe entre les deux cotes ne se verrait pas tout
  ! de suite -- elle donnerait une vitesse nodale parasite sur un fluide immobile.
  ! ===================================================================
  subroutine compute_rhs_lagrange_divdirect(mesh, sol, vp, dt, rhs, n_bc, bc_type, &
      bc_val, b2d, mass, gamma_arr, vp_is_imposed)
    use linear_solver_module, only: inverse_3_by_3
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), intent(in) :: dt
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mass, gamma_arr
    logical, dimension(mesh%n_vert), intent(in) :: vp_is_imposed
    integer(kind=ENTIER), intent(in) :: n_bc
    character(len=255), dimension(n_bc) :: bc_type
    real(kind=DOUBLE), dimension(5, n_bc) :: bc_val
    logical, intent(in) :: b2d

    integer(kind=ENTIER) :: i, j, k, id_sf, id_f, idl, idr, lse, rse, ide
    real(kind=DOUBLE), dimension(3, 3) :: Mp, Mp_inv
    real(kind=DOUBLE), dimension(3) :: Rp, Bp, nrm, vp_i, vp_old, wvec, first_wall_norm
    real(kind=DOUBLE), dimension(5) :: flux
    real(kind=DOUBLE) :: area, pl, pr, lam_l, lam_r, cl, cr, Pp
    real(kind=DOUBLE) :: Pi_l, Pi_r, v_et, vl, vr, cross_z, prod
    logical :: is_corner, have_wall

    rhs = 0.0_DOUBLE
    if (allocated(divd_prod_cell)) divd_prod_cell = 0.0_DOUBLE

    do i = 1, mesh%n_vert
      Mp = 0.0_DOUBLE; Rp = 0.0_DOUBLE; Bp = 0.0_DOUBLE
      is_corner = .false.; have_wall = .false.; first_wall_norm = 0.0_DOUBLE

      do j = 1, mesh%vert(i)%n_sub_faces_neigh
        id_sf = mesh%vert(i)%sub_face_neigh(j)
        id_f  = mesh%sub_face(id_sf)%mesh_face
        idl   = mesh%face(id_f)%left_neigh
        idr   = mesh%face(id_f)%right_neigh
        area  = mesh%sub_face(id_sf)%area
        nrm   = mesh%face(id_f)%norm
        lse   = mesh%sub_face(id_sf)%left_sub_elem_neigh
        if (lse <= 0) cycle
        pl    = pressure(sol(:, idl), gamma_arr(idl))
        lam_l = sqrt(max(0.0_DOUBLE, gamma_arr(idl)*pl*sol(1, idl)))/sol(1, idl)
        cl    = lam_l*divd_isum(lse)*divd_psi(lse)
        if (trim(divd_mode) == 'h') cl = divd_coef*divd_h(idl) &
                                       / max(mesh%sub_elem(lse)%volume, 1.0e-300_DOUBLE)

        if (idr > 0) then
          rse = mesh%sub_face(id_sf)%right_sub_elem_neigh
          if (rse <= 0) cycle
          pr    = pressure(sol(:, idr), gamma_arr(idr))
          lam_r = sqrt(max(0.0_DOUBLE, gamma_arr(idr)*pr*sol(1, idr)))/sol(1, idr)
          cr    = lam_r*divd_isum(rse)*divd_psi(rse)
          if (trim(divd_mode) == 'h') cr = divd_coef*divd_h(idr) &
                                         / max(mesh%sub_elem(rse)%volume, 1.0e-300_DOUBLE)
          wvec = cl*divd_b(:, lse) - cr*divd_b(:, rse)
          do k = 1, 3
            Mp(k, :) = Mp(k, :) + area*nrm(k)*wvec(:)
          end do
          if (divd_wdbg > 0 .and. i == divd_wnode .and. &
              divd_w_printed < divd_wmax) then
            divd_w_printed = divd_w_printed + 1
            print*, "VPDIVW  n", nrm(1), nrm(2), "  w", wvec(1), wvec(2), &
                    "  A", area, "  cl", cl, " cr", cr
          end if
          Rp = Rp + area*((pl - pr) - cl*divd_rest(lse) + cr*divd_rest(rse))*nrm
        else if (idr < 0 .and. trim(bc_type(-idr)) == 'pressure') then
          wvec = cl*divd_b(:, lse)
          do k = 1, 3
            Mp(k, :) = Mp(k, :) + area*nrm(k)*wvec(:)
          end do
          Rp = Rp + area*((pl - bc_val(1, -idr)) - cl*divd_rest(lse))*nrm
        else
          ! MIROIR -- face z de l'extrusion comme paroi. Contrairement a vp, dont
          ! la fermeture n'est pas symetrique sous le miroir et laisse subsister
          ! A(p_l + lam v_l)n, la notre l'est : p_d = p_c ET D_pd = D_pc par
          ! symetrie, donc Pi_pcf - Pi_pdf = 0 EXACTEMENT. Une face miroir ne
          ! contribue donc ni a M_p ni a R_p ; elle n'impose qu'une contrainte,
          ! v_p.n = 0, portee par B_p.
          ! Mesure qui l'a montre : R_p = (0, 0, +/-4.4648) aux noeuds interieurs
          ! de Gresho, soit exactement A p_c = 6.25e-4 x 7143, injecte par la
          ! SEULE face z du noeud -- un maillage extrude d'une couche a ses noeuds
          ! soit en z = 0 soit en z = dz, donc rien ne la compense.
          if (.true.) then
            ! Paroi : meme detection de coin que vp -- deux normales de paroi non
            ! colineaires epinglent le noeud.
            if (.not. (b2d .and. abs(nrm(3)) > 1.0e-8_DOUBLE)) then
              Bp = Bp + area*nrm
              if (.not. have_wall) then
                first_wall_norm = nrm; have_wall = .true.
              else
                cross_z = first_wall_norm(1)*nrm(2) - first_wall_norm(2)*nrm(1)
                if (abs(cross_z) > 1.0e-8_DOUBLE) is_corner = .true.
              end if
            end if
          end if
        end if
      end do

      ! Diagnostic ponctuel : la matrice nodale, son determinant et le second
      ! membre aux premiers noeuds INTERIEURS du premier appel. Sans ca on ne
      ! peut pas distinguer une matrice singuliere d'un second membre aberrant.
      if (divd_dbg < 6 .and. .not. mesh%vert(i)%is_bound .and. .not. mesh%vert(i)%is_ghost) then
        divd_dbg = divd_dbg + 1
        print*, "VPDIVDBG noeud", i, " det", &
          Mp(1,1)*(Mp(2,2)*Mp(3,3)-Mp(2,3)*Mp(3,2)) &
         -Mp(1,2)*(Mp(2,1)*Mp(3,3)-Mp(2,3)*Mp(3,1)) &
         +Mp(1,3)*(Mp(2,1)*Mp(3,2)-Mp(2,2)*Mp(3,1))
        print*, "VPDIVDBG   M ligne1", Mp(1,:)
        print*, "VPDIVDBG   M ligne2", Mp(2,:)
        print*, "VPDIVDBG   M ligne3", Mp(3,:)
        print*, "VPDIVDBG   R", Rp
      end if

      ! RANG. M_p = somme_{c,f} A n (x) w_pc, et w_pc ne depend PAS de f : la
      ! somme sur les sous-faces d'une maille se contracte en N_pc (x) w_pc, une
      ! matrice de RANG 1 par maille. Donc rang(M_p) <= nombre de mailles au
      ! noeud, et il faut au moins d mailles. C'est la difference de fond avec
      ! vp, dont la matrice somme A lam n(x)n recoit une direction propre de
      ! CHAQUE sous-face et est de rang plein des que les normales engendrent
      ! l'espace -- y compris a un noeud de coin a une seule maille.
      !
      ! Sous boundary_2d on resout donc le bloc 2x2 EN PLAN et on impose
      ! v_p,z = 0, au lieu d'inverser une 3x3 dont la 3e ligne est vide depuis
      ! que les faces z, miroir, ne contribuent plus. Le compte des mailles au
      ! noeud devient : 4 a l'interieur, 2 a une paroi (rang 2, une contrainte,
      ! un degre libre), 1 a un coin (rang 1, deux contraintes, zero degre libre
      ! -- deja traite par is_corner).
      vp_i = vp(:, i)
      vp_old = vp(:, i)
      if (.not. vp_is_imposed(i)) then
        if (is_corner) then
          vp_i = 0.0_DOUBLE
        else if (b2d) then
          block
            real(kind=DOUBLE) :: d2, a11, a12, a21, a22, r1, r2, b1, b2, s1, s2, den
            a11 = Mp(1,1); a12 = Mp(1,2); a21 = Mp(2,1); a22 = Mp(2,2)
            r1 = Rp(1); r2 = Rp(2); b1 = Bp(1); b2 = Bp(2)
            d2 = a11*a22 - a12*a21
            if (abs(d2) <= 1.0e-300_DOUBLE) then
              divd_n_sing = divd_n_sing + 1
              ! Repartition : bord ou interieur, et nombre de mailles au noeud.
              ! La borne rang(M_p) <= nb de mailles predit que seuls les noeuds
              ! pauvres soient singuliers ; si des noeuds a quatre mailles le
              ! sont, c'est autre chose et il faut le voir.
              if (mesh%vert(i)%is_bound) then
                divd_sing_bnd = divd_sing_bnd + 1
              else
                divd_sing_int = divd_sing_int + 1
                if (divd_sing_print < 3) then
                  divd_sing_print = divd_sing_print + 1
                  if (divd_wnode < 0) then
                    divd_wnode = i; divd_wdbg = 1
                  end if
                  print*, "VPDIVSING noeud INTERIEUR", i, " mailles", &
                    mesh%vert(i)%n_sub_elems_neigh, " d2", d2
                  print*, "VPDIVSING   M2", a11, a12, a21, a22
                  print*, "VPDIVSING   R2", r1, r2
                end if
              end if
              vp_i = 0.0_DOUBLE
            else
              ! M2^-1 R2
              s1 = ( a22*r1 - a12*r2)/d2
              s2 = (-a21*r1 + a11*r2)/d2
              if (abs(b1) + abs(b2) > 1.0e-8_DOUBLE) then
                ! M2^-1 B2, puis le multiplicateur de la contrainte B2.v = 0
                den = b1*(( a22*b1 - a12*b2)/d2) + b2*((-a21*b1 + a11*b2)/d2)
                if (abs(den) > 1.0e-300_DOUBLE) then
                  Pp = (b1*s1 + b2*s2)/den
                  s1 = s1 - Pp*(( a22*b1 - a12*b2)/d2)
                  s2 = s2 - Pp*((-a21*b1 + a11*b2)/d2)
                else
                  divd_n_sing = divd_n_sing + 1
                  s1 = 0.0_DOUBLE; s2 = 0.0_DOUBLE
                end if
              end if
              vp_i = (/ s1, s2, 0.0_DOUBLE /)
            end if
          end block
        else if (maxval(abs(Mp)) <= 0.0_DOUBLE) then
          divd_n_sing = divd_n_sing + 1
          vp_i = 0.0_DOUBLE
        else if (maxval(abs(Bp)) > 1e-8_DOUBLE) then
          call inverse_3_by_3(Mp, Mp_inv)
          Pp = dot_product(Rp, matmul(Mp_inv, Bp))/dot_product(Bp, matmul(Mp_inv, Bp))
          vp_i = matmul(Mp_inv, Rp - Pp*Bp)
        else
          call inverse_3_by_3(Mp, Mp_inv)
          vp_i = matmul(Mp_inv, Rp)
        end if
      end if
      ! Sous-relaxation : (1-theta) v_p^ancien + theta v_p^resolu. On l'applique
      ! APRES le traitement des coins et des parois, pour ne pas relaxer une
      ! contrainte geometrique -- un noeud epingle doit rester a zero.
      if (.not. vp_is_imposed(i) .and. .not. is_corner &
          .and. divd_theta < 1.0_DOUBLE) then
        vp_i = (1.0_DOUBLE - divd_theta)*vp_old + divd_theta*vp_i
      end if
      vp(:, i) = vp_i

      ! --- flux
      do j = 1, mesh%vert(i)%n_sub_faces_neigh
        id_sf = mesh%vert(i)%sub_face_neigh(j)
        id_f  = mesh%sub_face(id_sf)%mesh_face
        idl   = mesh%face(id_f)%left_neigh
        idr   = mesh%face(id_f)%right_neigh
        area  = mesh%sub_face(id_sf)%area
        nrm   = mesh%face(id_f)%norm
        lse   = mesh%sub_face(id_sf)%left_sub_elem_neigh
        if (lse <= 0) cycle
        v_et  = dot_product(vp_i, nrm)
        pl    = pressure(sol(:, idl), gamma_arr(idl))
        vl    = dot_product(sol(2:4, idl), nrm)
        lam_l = sqrt(max(0.0_DOUBLE, gamma_arr(idl)*pl*sol(1, idl)))/sol(1, idl)
        cl    = lam_l*divd_isum(lse)*divd_psi(lse)
        if (trim(divd_mode) == 'h') cl = divd_coef*divd_h(idl) &
                                       / max(mesh%sub_elem(lse)%volume, 1.0e-300_DOUBLE)
        Pi_l  = pl - cl*(dot_product(divd_b(:, lse), vp_i) + divd_rest(lse))
        call divd_close(pl, cl*(dot_product(divd_b(:, lse), vp_i) + divd_rest(lse)), &
                        gamma_arr(idl), Pi_l, .true.)

        flux(1)   = -v_et
        flux(2:4) = Pi_l*nrm
        flux(5)   = Pi_l*v_et
        if (mesh%sub_elem(lse)%mesh_vert == i) then
          ide = mesh%sub_elem(lse)%mesh_elem
          rhs(:, ide) = rhs(:, ide) - area/mass(idl)*flux
          if (allocated(divd_prod_cell)) divd_prod_cell(idl) = divd_prod_cell(idl) &
            + area*(pl - Pi_l)*(v_et - vl)
        end if

        if (idr > 0) then
          rse = mesh%sub_face(id_sf)%right_sub_elem_neigh
          if (rse <= 0) cycle
          pr    = pressure(sol(:, idr), gamma_arr(idr))
          vr    = dot_product(sol(2:4, idr), nrm)
          lam_r = sqrt(max(0.0_DOUBLE, gamma_arr(idr)*pr*sol(1, idr)))/sol(1, idr)
          cr    = lam_r*divd_isum(rse)*divd_psi(rse)
          if (trim(divd_mode) == 'h') cr = divd_coef*divd_h(idr) &
                                         / max(mesh%sub_elem(rse)%volume, 1.0e-300_DOUBLE)
          Pi_r  = pr - cr*(dot_product(divd_b(:, rse), vp_i) + divd_rest(rse))
          call divd_close(pr, cr*(dot_product(divd_b(:, rse), vp_i) + divd_rest(rse)), &
                          gamma_arr(idr), Pi_r, .true.)
          flux(1)   = -v_et
          flux(2:4) = Pi_r*nrm
          flux(5)   = Pi_r*v_et
          if (mesh%sub_elem(rse)%mesh_vert == i) then
            ide = mesh%sub_elem(rse)%mesh_elem
            rhs(:, ide) = rhs(:, ide) + area/mass(idr)*flux
            if (allocated(divd_prod_cell)) divd_prod_cell(idr) = divd_prod_cell(idr) &
              + area*(pr - Pi_r)*(-(v_et - vr))
          end if
        end if
      end do
    end do

    ! --- increments des dix premiers pas : cause ou symptome
    divd_step = divd_step + 1
    if (divd_step <= 10) then
      print*, "VPDIVPAS", divd_step, " Pi<0", divd_n_pineg - divd_prev(1), &
        " sur", divd_n_pi - divd_prev(2), &
        " | M_div >0.1", divd_md(1) - divd_prev(3), " >1", divd_md(3) - divd_prev(4), &
        " >2.857", divd_md(4) - divd_prev(5), " >20", divd_md(7) - divd_prev(6), &
        " | vide", divd_n_vide - divd_prev(7)
      divd_prev(1) = divd_n_pineg; divd_prev(2) = divd_n_pi
      divd_prev(3) = divd_md(1);   divd_prev(4) = divd_md(3)
      divd_prev(5) = divd_md(4);   divd_prev(6) = divd_md(7)
      divd_prev(7) = divd_n_vide
    end if

    ! --- production PAR MAILLE, avec son signe. Le compteur de COINS negatifs
    ! serait non informatif ici : la production n'est pas un carre, elle vaut un
    ! carre plus un terme indefini, et la mesure du pair donne ~50 % de coins
    ! negatifs par simple symetrie du bruit. C'est la somme par maille qui vide
    ! ou remplit la maille, donc c'est elle le RESULTAT -- pas un controle.
    if (allocated(divd_prod_cell)) then
      do i = 1, mesh%n_elems
        if (mesh%elem(i)%is_ghost) cycle
        divd_n_cell = divd_n_cell + 1
        divd_prod_tot = divd_prod_tot + divd_prod_cell(i)
        if (divd_prod_cell(i) < 0.0_DOUBLE) then
          divd_n_cell_neg = divd_n_cell_neg + 1
          divd_prod_min = min(divd_prod_min, divd_prod_cell(i))
        end if
      end do
    end if
  end subroutine compute_rhs_lagrange_divdirect

  ! Fermeture d'une sous-face : lineaire ou isentropique, plus l'instrumentation.
  ! kd = kappa_pc D_pc (le decrement absolu), pc la pression de maille.
  subroutine divd_close(pc, kd, gam, Pi, diag_only)
    implicit none
    real(kind=DOUBLE), intent(in) :: pc, kd, gam
    real(kind=DOUBLE), intent(inout) :: Pi
    logical, intent(in) :: diag_only
    real(kind=DOUBLE) :: q, base, mdiv
    integer(kind=ENTIER) :: k

    if (pc <= 0.0_DOUBLE) return
    q = kd/pc
    mdiv = 4.0_DOUBLE*q/gam          ! M_div signe
    divd_n_pi = divd_n_pi + 1
    ! distribution de |M_div|, par seuils
    if (abs(mdiv) > 0.1_DOUBLE)   divd_md(1) = divd_md(1) + 1
    if (abs(mdiv) > 0.5_DOUBLE)   divd_md(2) = divd_md(2) + 1
    if (abs(mdiv) > 1.0_DOUBLE)   divd_md(3) = divd_md(3) + 1
    if (abs(mdiv) > 2.857_DOUBLE) divd_md(4) = divd_md(4) + 1
    if (abs(mdiv) > 5.0_DOUBLE)   divd_md(5) = divd_md(5) + 1
    if (abs(mdiv) > 10.0_DOUBLE)  divd_md(6) = divd_md(6) + 1
    if (abs(mdiv) > 20.0_DOUBLE)  divd_md(7) = divd_md(7) + 1
    if (Pi < 0.0_DOUBLE) divd_n_pineg = divd_n_pineg + 1

    ! En mode psi, Pi est DEJA la forme isentropique : le coefficient a ete
    ! reechelonne en amont, donc on ne recalcule rien ici, on ne fait que mesurer.
    if (diag_only .or. .not. divd_nl) return
    base = 1.0_DOUBLE - (gam - 1.0_DOUBLE)/(2.0_DOUBLE*gam)*q
    if (base <= 0.0_DOUBLE) then
      ! VIDE. Signale et compte, pas ecrase : Pi = 0 est la pression du vide et
      ! c'est un choix physique, pas un plafond arbitraire. (negatif)^7 rendrait
      ! un NaN -- c'est exactement la faute du TRRS deja consignee.
      divd_n_vide = divd_n_vide + 1
      Pi = 0.0_DOUBLE
    else
      Pi = pc*base**(2.0_DOUBLE*gam/(gam - 1.0_DOUBLE))
    end if
  end subroutine divd_close

  ! ===================================================================
  ! 'mean_div' : Pi_f = p_bar^lambda_f - mu_f (div v)_f, UNIVOQUE par sous-face.
  !
  !   lambda_bar_f = lam_c lam_d / (lam_c + lam_d)
  !   p_bar^lam_f  = (lam_d p_c + lam_c p_d) / (lam_c + lam_d)
  !   Vpoint_c     = somme_{p dans P(c)} v_p . N_pc        (vitesses NODALES)
  !   (div v)_f    = 1/2 [ (div v)_c + (div v)_d ]
  !   mu_f         = c_visc lambda_bar_f h_f,  NUL en detente
  !
  ! Le MEME Pi_f dans les deux flux (compatibilite de Caramana). v_p vient du
  ! solveur vp existant et n'est pas re-resolu.
  !
  ! TROIS PARCOURS, UN SEUL SOLVE. Vpoint_c a besoin de TOUTES les vitesses
  ! nodales de la maille, or le code resout v_p et assemble les flux dans la
  ! meme boucle : au noeud p les v_p' des autres noeuds de c sont inconnus. Il
  ! faut donc scinder. Rien n'est mal pose et le solveur n'est pas touche, le
  ! surcout est un parcours de sous-faces, pas un solve.
  !
  ! Bords : etat miroir, p_d = p_c, lam_d = lam_c, (div v)_d = (div v)_c, d'ou
  ! p_bar^lam_f = p_c et (div v)_f = (div v)_c. La conservation vient ici de
  ! l'UNIVOCITE de Pi_f et non de somme_c N_pc = 0 -- contrairement a vp_mp, ou
  ! cette identite etait seule a la garantir et ou le bord etait donc critique.
  ! ===================================================================
  subroutine compute_rhs_lagrange_meandiv(mesh, sol, vp, dt, rhs, n_bc, bc_type, &
      bc_val, b2d, mass, gamma_arr, vp_is_imposed, grad_v, grad_p, div_v)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), intent(in) :: dt
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mass, gamma_arr
    logical, dimension(mesh%n_vert), intent(in) :: vp_is_imposed
    integer(kind=ENTIER), intent(in) :: n_bc
    character(len=255), dimension(n_bc) :: bc_type
    real(kind=DOUBLE), dimension(5, n_bc) :: bc_val
    logical, intent(in) :: b2d
    real(kind=DOUBLE), dimension(3, 3, mesh%n_elems), intent(in) :: grad_v
    real(kind=DOUBLE), dimension(3, mesh%n_elems), intent(in) :: grad_p
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: div_v

    integer(kind=ENTIER) :: i, j, k, nsfn, nsen, idse, ide
    integer(kind=ENTIER) :: id_sf, id_f, idl, idr, lse, lse_loc, rse, rse_loc
    real(kind=DOUBLE), dimension(:,:), allocatable :: rhs_o1
    real(kind=DOUBLE), dimension(3) :: vp_i, nrm, res_node
    real(kind=DOUBLE), dimension(5) :: flux
    real(kind=DOUBLE), allocatable :: vdot(:), divc(:)
    real(kind=DOUBLE) :: area, v_et, pl, pr, lam_c, lam_d, lam_bar, p_bar
    real(kind=DOUBLE) :: divf, mu_f, h_f, Pi_f, vl, vr, prod
    real(kind=DOUBLE) :: a_c, a_d, gm, zz, pst

    allocate(vdot(mesh%n_elems), source=0.0_DOUBLE)
    allocate(divc(mesh%n_elems), source=0.0_DOUBLE)

    ! ---------- PARCOURS 1 : resolution nodale, vp inchangee ----------
    do i = 1, mesh%n_vert
      nsfn = mesh%vert(i)%n_sub_faces_neigh
      nsen = mesh%vert(i)%n_sub_elems_neigh
      allocate(rhs_o1(5, nsen))
      vp_i = vp(:, i)
      call compute_rhs_around_node(mesh, sol, i, nsfn, nsen, vp_i, vp_is_imposed(i), &
        dt, rhs_o1, n_bc, bc_type, bc_val, b2d, mass, gamma_arr, &
        .false., grad_v, grad_p, div_v)
      vp(:, i) = vp_i       ! le rhs de ce parcours est JETE : seul v_p nous interesse
      deallocate(rhs_o1)
    end do

    ! ---------- PARCOURS 2 : Vpoint_c puis (div v)_c ----------
    do i = 1, mesh%n_vert
      do j = 1, mesh%vert(i)%n_sub_faces_neigh
        id_sf = mesh%vert(i)%sub_face_neigh(j)
        id_f  = mesh%sub_face(id_sf)%mesh_face
        idl   = mesh%face(id_f)%left_neigh
        idr   = mesh%face(id_f)%right_neigh
        area  = mesh%sub_face(id_sf)%area
        v_et  = dot_product(vp(:, i), mesh%face(id_f)%norm)
        vdot(idl) = vdot(idl) + area*v_et
        if (idr > 0) vdot(idr) = vdot(idr) - area*v_et
      end do
    end do
    do i = 1, mesh%n_elems
      if (mesh%elem(i)%volume > 0.0_DOUBLE) divc(i) = vdot(i)/mesh%elem(i)%volume
    end do
    ! Le NaN mesure est dans divc, pas dans la pression : p, a, lambda, p_bar et
    ! pst sont tous finis au premier Pi_f NaN, et divf vaut NaN parce que divc du
    ! voisin vaut NaN. On remonte donc a la source : v_p apres la passe 1, puis
    ! vdot apres la passe 2.
    if (.not. md_src_done) then
      block
        integer(kind=ENTIER) :: kk, nvp, nvd, firstvd
        nvp = 0; nvd = 0; firstvd = -1
        do kk = 1, mesh%n_vert
          if (any(vp(:, kk) /= vp(:, kk))) nvp = nvp + 1
        end do
        do kk = 1, mesh%n_elems
          if (vdot(kk) /= vdot(kk)) then
            nvd = nvd + 1
            if (firstvd < 0) firstvd = kk
          end if
        end do
        if (nvp > 0 .or. nvd > 0) then
          md_src_done = .true.
          print*, "MDSRC noeuds a v_p NaN :", nvp, " sur", mesh%n_vert
          print*, "MDSRC mailles a vdot NaN :", nvd, " sur", mesh%n_elems
          if (firstvd > 0) print*, "MDSRC premiere maille vdot NaN :", firstvd, &
            " ghost ?", mesh%elem(firstvd)%is_ghost, " V", mesh%elem(firstvd)%volume
        end if
      end block
    end if

    ! ---------- PARCOURS 3 : assemblage avec Pi_f ----------
    rhs = 0.0_DOUBLE
    do i = 1, mesh%n_vert
      res_node = 0.0_DOUBLE
      do j = 1, mesh%vert(i)%n_sub_faces_neigh
        id_sf = mesh%vert(i)%sub_face_neigh(j)
        id_f  = mesh%sub_face(id_sf)%mesh_face
        idl   = mesh%face(id_f)%left_neigh
        idr   = mesh%face(id_f)%right_neigh
        area  = mesh%sub_face(id_sf)%area
        nrm   = mesh%face(id_f)%norm
        v_et  = dot_product(vp(:, i), nrm)

        pl    = pressure(sol(:, idl), gamma_arr(idl))
        vl    = dot_product(sol(2:4, idl), nrm)
        a_c   = sqrt(max(0.0_DOUBLE, gamma_arr(idl)*pl*sol(1, idl)))
        lam_c = a_c/sol(1, idl)
        if (idr > 0) then
          pr    = pressure(sol(:, idr), gamma_arr(idr))
          vr    = dot_product(sol(2:4, idr), nrm)
          a_d   = sqrt(max(0.0_DOUBLE, gamma_arr(idr)*pr*sol(1, idr)))
          lam_d = a_d/sol(1, idr)
          ! --- impedance a DEUX ONDES (TRRS de Toro) : l'acoustique rho*a donne
          ! lam_c/lam_d = sqrt(p_c/p_d), qui est enorme sous saut fort, et
          ! p_bar^lam s'effondre alors sur la moyenne geometrique sqrt(p_c p_d).
          ! Le TRRS fait voir le choc a l'impedance et rend p_bar^lam proche de
          ! la pression etoile exacte.
          if (md_trrs .and. pl > 0.0_DOUBLE .and. pr > 0.0_DOUBLE) then
            gm = gamma_arr(idl)
            zz = (gm - 1.0_DOUBLE)/(2.0_DOUBLE*gm)
            ! GARDE OBLIGATOIRE : le numerateur a_c + a_d - (g-1)/2 (v_r - v_l)
            ! devient NEGATIF quand la detente est assez forte -- c'est la
            ! condition d'apparition du vide -- et (negatif)**7 rend NaN. Mesure
            ! qui l'a montre : sur sedov_tri, Pi - p_c valait deja NaN au pas 2
            ! dans les mailles du depot, alors que l'etat sol y etait encore sain.
            ! Dans ce cas il n'y a pas de choc, donc l'impedance acoustique est la
            ! bonne reponse et on ne corrige pas.
            pst = a_c + a_d - 0.5_DOUBLE*(gm - 1.0_DOUBLE)*(vr - vl)
            if (pst > 0.0_DOUBLE) then
              pst = (pst/(a_c/pl**zz + a_d/pr**zz))**(1.0_DOUBLE/zz)
            else
              pst = 0.0_DOUBLE        ! vide : pas de correction de choc
              md_n_vac_step = md_n_vac_step + 1
              md_n_vac_tot = md_n_vac_tot + 1
            end if
            md_n_sf_step = md_n_sf_step + 1
            if (pst > pl) lam_c = lam_c*sqrt(1.0_DOUBLE &
              + 0.5_DOUBLE*(gm + 1.0_DOUBLE)/gm*(pst/pl - 1.0_DOUBLE))
            if (pst > pr) lam_d = lam_d*sqrt(1.0_DOUBLE &
              + 0.5_DOUBLE*(gm + 1.0_DOUBLE)/gm*(pst/pr - 1.0_DOUBLE))
          end if
          divf  = 0.5_DOUBLE*(divc(idl) + divc(idr))
          h_f   = norm2(mesh%elem(idr)%coord - mesh%elem(idl)%coord)
        else
          pr = pl; vr = -vl; lam_d = lam_c          ! etat miroir
          divf = divc(idl)
          h_f  = 2.0_DOUBLE*abs(dot_product(mesh%vert(i)%coord &
                 - mesh%elem(idl)%coord, nrm))
        end if
        if (lam_c + lam_d <= 0.0_DOUBLE) cycle
        lam_bar = lam_c*lam_d/(lam_c + lam_d)
        p_bar   = (lam_d*pl + lam_c*pr)/(lam_c + lam_d)
        mu_f    = 0.0_DOUBLE
        if (divf < 0.0_DOUBLE .or. md_noswitch) mu_f = md_cvisc*lam_bar*h_f
        Pi_f    = p_bar - mu_f*divf
        ! Premier Pi_f NaN : on imprime TOUS ses ingredients. Deduire lequel est
        ! en cause par le raisonnement a deja echoue une fois -- la garde de vide
        ! n'a rien change au resultat, donc le NaN vient d'ailleurs.
        if (.not. md_nan_done) then
          if (Pi_f /= Pi_f) then
            md_nan_done = .true.
            print*, "MDNAN Pi_f NaN : p_l", pl, " p_r", pr
            print*, "MDNAN   a_c", a_c, " a_d", a_d, " lam_c", lam_c, " lam_d", lam_d
            print*, "MDNAN   lam_bar", lam_bar, " p_bar", p_bar, " pst", pst
            print*, "MDNAN   divf", divf, " mu_f", mu_f, " h_f", h_f
            print*, "MDNAN   divc(l)", divc(idl), " divc(r)", divc(max(idr,idl)), &
                    " V_l", mesh%elem(idl)%volume
          end if
        end if
        if (md_first .and. .not. md_first_done) then
          md_pmax  = max(md_pmax, pl)
          ! Les sous-faces de BORD ont p_d = p_c par etat miroir, donc p_bar^lam
          ! = p_c exactement : les inclure masque completement le comportement des
          ! faces interieures, qui est le seul qui nous interesse ici. Mesure qui
          ! l'a montre : max Pi_f = 254.27918625000049 = max p_c a 15 chiffres,
          ! entierement du aux faces z de l'extrusion.
          if (idr > 0) then
            md_pimax = max(md_pimax, Pi_f)
            md_pimin = min(md_pimin, Pi_f)
            if (Pi_f < 0.0_DOUBLE) md_n_pi_neg = md_n_pi_neg + 1
          else
            md_pimax_b = max(md_pimax_b, Pi_f)
          end if
        end if

        ! Le CONTROLE exclut les faces z de l'extrusion : elles sont unilaterales
        ! (idr < 0) alors que leurs noeuds ne sont pas marques is_bound, donc leur
        ! contribution ne peut pas se compenser et elle salit un bilan qui est
        ! exact dans le plan. Mesure : x = 3.0e-37, y = 6.0e-37, z = 2.5e-04.
        !
        ! RAISON EXACTE du masquage, et elle est conditionnee a boundary_2d : sous
        ! boundary_2d la composante w est contrainte a zero independamment, donc la
        ! contribution z non compensee est DYNAMIQUEMENT INERTE -- elle n'accelere
        ! rien parce qu'on l'ecarte en aval, et le flux d'energie A Pi (v_p.n_z) est
        ! nul parce que w_p = 0, pas parce que les contributions s'annulent.
        ! Ce n'est donc PAS que le flux les traite correctement : c'est qu'elles ne
        ! peuvent rien faire dans ce mode. Sur un vrai cas 3D sans boundary_2d, la
        ! meme contribution serait une erreur de quantite de mouvement REELLE, et
        ! ce masque ne doit surtout pas s'y appliquer -- d'ou le test sur b2d.
        if (.not. (b2d .and. abs(nrm(3)) > 1.0e-8_DOUBLE)) &
          res_node = res_node + area*Pi_f*nrm
        md_sq_div   = md_sq_div   + mu_f*divf**2*area*h_f
        md_sq_delta = md_sq_delta + lam_c*area*(v_et - vl)**2
        prod = area*(pl - Pi_f)*(v_et - vl)
        md_n_cell_tot = md_n_cell_tot + 1
        if (prod < 0.0_DOUBLE) then
          md_n_cell_neg = md_n_cell_neg + 1
          md_prod_min = min(md_prod_min, prod)
        end if

        if (md_watch) then
          w_dpi_min(idl) = min(w_dpi_min(idl), Pi_f - pl)
          w_dpi_max(idl) = max(w_dpi_max(idl), Pi_f - pl)
          w_prod(idl) = w_prod(idl) + area*(pl - Pi_f)*(v_et - vl)
          if (idr > 0) then
            w_dpi_min(idr) = min(w_dpi_min(idr), Pi_f - pr)
            w_dpi_max(idr) = max(w_dpi_max(idr), Pi_f - pr)
            w_prod(idr) = w_prod(idr) + area*(pr - Pi_f)*(-(v_et - vr))
          end if
        end if
        flux(1)   = -v_et
        flux(2:4) = Pi_f*nrm
        flux(5)   = Pi_f*v_et
        lse = mesh%sub_face(id_sf)%left_sub_elem_neigh
        if (mesh%sub_elem(lse)%mesh_vert == i) then
          lse_loc = mesh%sub_elem(lse)%id_loc_around_node
          rhs(:, idl) = rhs(:, idl) - area/mass(idl)*flux
        end if
        if (idr > 0) then
          if (.not. (b2d .and. abs(nrm(3)) > 1.0e-8_DOUBLE)) &
            res_node = res_node - area*Pi_f*nrm
          md_sq_delta = md_sq_delta + lam_d*area*(vr - v_et)**2
          prod = area*(pr - Pi_f)*(vr - v_et)
          md_n_cell_tot = md_n_cell_tot + 1
          if (prod < 0.0_DOUBLE) then
            md_n_cell_neg = md_n_cell_neg + 1
            md_prod_min = min(md_prod_min, prod)
          end if
          rse = mesh%sub_face(id_sf)%right_sub_elem_neigh
          if (rse > 0) then
            if (mesh%sub_elem(rse)%mesh_vert == i) &
              rhs(:, idr) = rhs(:, idr) + area/mass(idr)*flux
          end if
        end if
      end do
      ! Bilan nodal : nul par UNIVOCITE de Pi_f aux noeuds INTERIEURS. Au bord la
      ! somme vaut la force de paroi et n'a aucune raison d'etre nulle, donc on
      ! ne la controle que la ou l'identite s'applique.
      if (.not. mesh%vert(i)%is_bound) then
        md_res_max = max(md_res_max, norm2(res_node))
        md_res_dir(1) = max(md_res_dir(1), abs(res_node(1)))
        md_res_dir(2) = max(md_res_dir(2), abs(res_node(2)))
        md_res_dir(3) = max(md_res_dir(3), abs(res_node(3)))
      end if
    end do

    ! --- energie interne PREDITE au premier pas. Le mecanisme a tester est que
    ! Pi_f, etant UNIVOQUE, impose la pression du choc a une maille dont l'energie
    ! interne est quinze ordres plus petite. Mais l'ordre de grandeur est le meme
    ! sur quad, qui survit : il faut donc sortir aussi le SIGNE du travail recu et
    ! le dt, pour separer "la maille est comprimee donc eps monte", "le dt la
    ! sauve" et "la raideur seule fait basculer".
    md_ncall = md_ncall + 1
    if (md_trace) then
      block
        integer(kind=ENTIER) :: kt
        real(kind=DOUBLE) :: vmin, taumax
        integer(kind=ENTIER) :: nvneg, nnan_tau, nnan_v, nnan_e
        vmin = 1.0e300_DOUBLE; taumax = -1.0e300_DOUBLE
        nvneg = 0; nnan_tau = 0; nnan_v = 0; nnan_e = 0
        do kt = 1, mesh%n_elems
          if (mesh%elem(kt)%is_ghost) cycle
          vmin = min(vmin, mesh%elem(kt)%volume)
          if (mesh%elem(kt)%volume <= 0.0_DOUBLE) nvneg = nvneg + 1
          if (sol(1, kt) /= sol(1, kt)) then
            nnan_tau = nnan_tau + 1
          else
            taumax = max(taumax, sol(1, kt))
          end if
          if (any(sol(2:4, kt) /= sol(2:4, kt))) nnan_v = nnan_v + 1
          if (sol(5, kt) /= sol(5, kt)) nnan_e = nnan_e + 1
        end do
        print*, "MDTRACE pas", md_ncall, " dt", dt, " Vmin", vmin, " V<=0", nvneg, &
                " tau_max", taumax, " NaN tau/v/E", nnan_tau, nnan_v, nnan_e
      end block
    end if
    if (md_first .and. .not. md_first_done .and. md_ncall >= 2) then
      block
        integer(kind=ENTIER) :: kc
        real(kind=DOUBLE) :: enew, epsn, wsum
        real(kind=DOUBLE), dimension(3) :: vnew
        md_n_eps_neg = 0; md_eps_min = 1.0e300_DOUBLE
        md_w_at_min = 0.0_DOUBLE; md_eps_old_at_min = 0.0_DOUBLE
        do kc = 1, mesh%n_elems
          if (mesh%elem(kc)%is_ghost) cycle
          vnew = sol(2:4, kc) + dt*rhs(2:4, kc)
          enew = sol(5, kc) + dt*rhs(5, kc)
          epsn = enew - 0.5_DOUBLE*dot_product(vnew, vnew)
          ! somme A Pi (v_p.n) recue par la maille = -mass*rhs(5)
          wsum = -mass(kc)*rhs(5, kc)
          if (epsn < 0.0_DOUBLE) md_n_eps_neg = md_n_eps_neg + 1
          if (epsn < md_eps_min) then
            md_eps_min = epsn
            md_w_at_min = wsum
            md_eps_old_at_min = sol(5, kc) &
              - 0.5_DOUBLE*dot_product(sol(2:4, kc), sol(2:4, kc))
          end if
        end do
        print*, "MDFIRST dt =", dt
        print*, "MDFIRST mailles a eps PREDITE negative :", md_n_eps_neg, &
                " sur", mesh%n_elems
        print*, "MDFIRST eps predite la plus basse =", md_eps_min, &
                "   son eps initiale =", md_eps_old_at_min
        print*, "MDFIRST somme A Pi (v_p.n) sur cette maille =", md_w_at_min, &
                "   (positif = elle PERD de l'energie)"
      end block
      ! max_c |Vpoint_c/V_c| au PREMIER pas : construite sur les vitesses NODALES,
      ! elle n'est PAS nulle au repos puisque v_p est pilote par le saut de
      ! pression. C'est le discriminant qui separe "la base p_bar^lam tue" de
      ! "la stabilisation est absente au premier pas".
      block
        integer(kind=ENTIER) :: kc
        real(kind=DOUBLE) :: dmax
        dmax = 0.0_DOUBLE
        do kc = 1, mesh%n_elems
          if (.not. mesh%elem(kc)%is_ghost) dmax = max(dmax, abs(divc(kc)))
        end do
        print*, "MDFIRST max_c |Vpoint_c/V_c| =", dmax
      end block
      print*, "MDFIRST nombre d'appels au moment de la mesure =", md_ncall
      print*, "MDFIRST max p_c =", md_pmax, "   max Pi_f INTERIEUR =", md_pimax, &
              "   min Pi_f INTERIEUR =", md_pimin
      print*, "MDFIRST max Pi_f sur les faces MIROIR (p_d = p_c) =", md_pimax_b
      print*, "MDFIRST sous-faces a Pi_f < 0 :", md_n_pi_neg
      md_first_done = .true.
    end if
    deallocate(vdot, divc)
  end subroutine compute_rhs_lagrange_meandiv

  ! Accumule la contribution d'une sous-face a la force nette de la maille cible.
  ! 'nout' est la normale SORTANTE de la maille, donc +norm si elle est a gauche
  ! de la face et -norm si elle est a droite. Les deux chemins d'assemblage du
  ! rhs (gauche avec -A/m, droite avec +A/m) se ramenent tous deux a
  ! F = -somme A Pi n_sortante, ce qui rend l'accumulation uniforme.
  subroutine fdiag_add(icell, area, pival, nout, ineigh)
    implicit none
    integer(kind=ENTIER), intent(in) :: icell
    real(kind=DOUBLE), intent(in) :: area, pival
    real(kind=DOUBLE), dimension(3), intent(in) :: nout
    integer(kind=ENTIER), intent(in) :: ineigh
    integer(kind=ENTIER) :: k
    if (.not. diag_force) return
    do k = 1, 4
      if (fdiag_cells(k) /= icell) cycle
      if (fdiag_n(k) < 32) then
        fdiag_lPi(fdiag_n(k) + 1, k) = pival
        fdiag_lNb(fdiag_n(k) + 1, k) = ineigh
      end if
      fdiag_F(:, k) = fdiag_F(:, k) - area*pival*nout
      fdiag_sA(k)   = fdiag_sA(k)   + area
      fdiag_sAP(k)  = fdiag_sAP(k)  + area*pival
      fdiag_sAP2(k) = fdiag_sAP2(k) + area*pival*pival
      fdiag_Pmin(k) = min(fdiag_Pmin(k), pival)
      fdiag_Pmax(k) = max(fdiag_Pmax(k), pival)
      fdiag_n(k)    = fdiag_n(k) + 1
    end do
  end subroutine fdiag_add

  subroutine fdiag_reset()
    implicit none
    fdiag_F = 0.0_DOUBLE; fdiag_sA = 0.0_DOUBLE
    fdiag_sAP = 0.0_DOUBLE; fdiag_sAP2 = 0.0_DOUBLE
    fdiag_Pmin = 1.0e300_DOUBLE; fdiag_Pmax = -1.0e300_DOUBLE
    fdiag_n = 0
    fdiag_lPi = 0.0_DOUBLE; fdiag_lNb = -1
  end subroutine fdiag_reset

  subroutine fdiag_print()
    implicit none
    integer(kind=ENTIER) :: k, m
    real(kind=DOUBLE) :: pbar, var, sd
    do k = 1, 4
      if (fdiag_cells(k) < 0 .or. fdiag_n(k) == 0) cycle
      pbar = fdiag_sAP(k)/fdiag_sA(k)
      var  = fdiag_sAP2(k)/fdiag_sA(k) - pbar*pbar
      sd   = sqrt(max(0.0_DOUBLE, var))
      print*, "FDIAG maille", fdiag_cells(k), " nsf", fdiag_n(k), &
              " |F|", norm2(fdiag_F(:, k)), " F", fdiag_F(1, k), fdiag_F(2, k), &
              " Pi_moy", pbar, " ecart-type", sd, &
              " variation_rel", sd/max(abs(pbar), 1.0e-30_DOUBLE), &
              " Pi_min", fdiag_Pmin(k), " Pi_max", fdiag_Pmax(k)
      do m = 1, min(fdiag_n(k), 32)
        print*, "FDIAGSF maille", fdiag_cells(k), " sf", m, " Pi", fdiag_lPi(m, k), &
                " maille d'en face", fdiag_lNb(m, k)
      end do
    end do
  end subroutine fdiag_print

  ! --- Fermeture 'vp_fbd' : moyenne a la FACE moins divergence NODALE.
  !     Pi_pcf = pbar_f - D_p
  !     pbar_f = (p_l/z_l + p_r/z_r)/(1/z_l + 1/z_r)        LOCAL a la sous-face
  !     D_p    = somme_f A (v_nr - v_nl) / somme_f A (1/z_l + 1/z_r)   NODAL
  ! D_p est homogene a une pression -- (aire x vitesse)/(aire/impedance) -- donc
  ! AUCUNE longueur nodale a choisir, contrairement a `mu_p (div v)_p` de 'vp_mp'
  ! qui a la meme forme mais demande h_p et un coefficient libre.
  ! z est l'impedance acoustique PURE rho a. C'est un choix, et il est necessaire :
  ! le `lambda` du reste de ce code porte deux termes de saut supplementaires
  ! (pression et vitesse, ajout de type deux-chocs propre a ce code), avec lesquels
  ! la limite 1D ne redonnerait PAS le solveur acoustique. Les branches de bord
  ! existantes utilisent deja z seul, donc ce choix les rejoint au lieu de s'en
  ! ecarter.
  ! Impedance a trois termes, identique a `three_wave` (ns_euler_rs_module.F90:49)
  ! ecrite en tau : `max(a rho, sqrt(rho max(0,Dp)), -rho Du)`. L'impedance
  ! acoustique PURE etait le premier choix de cette fermeture, pour que la limite 1D
  ! redonne le solveur acoustique -- mais c'est la mauvaise cible : on a mesure que
  ! l'acoustique sous-estime la pression d'interface d'un facteur 146 a un rapport
  ! 1e5, et de 5e7 sur la force au depot de Sedov. Le terme sqrt(rho Dp) fait passer
  ! lambda du gaz FROID de 1.18e-06 a 48.99 au contact de Sedov : il cesse d'etre
  ! acoustiquement mou. A Dp -> 0 et Du -> 0 le max se reduit a rho a, donc la limite
  ! bas Mach est inchangee -- c'est le controle de non-regression de ce changement.
  pure function fbd_lambda(tau, pk, gk, p_autre, dvn) result(lam)
    implicit none
    real(kind=DOUBLE), intent(in) :: tau, pk, gk, p_autre, dvn
    real(kind=DOUBLE) :: lam
    lam = max(sqrt(max(0.0_DOUBLE, gk*pk*tau))/tau, &
              sqrt(max(0.0_DOUBLE, p_autre - pk)/tau), &
              -dvn/tau)
  end function fbd_lambda

  ! Bilan de positivite COMPLET de la sous-face. La forme par cote
  ! `|Pi - p_c| < 2.6458 p_c` suppose z_c = z_d ET eps_c = eps_d, donc deux etats
  ! IDENTIQUES : sur une interface chaud/froid elle est sans objet, parce que le
  ! terme `lambda_chaud eps_chaud` paie pour le deficit du cote froid. C'est le
  ! bilan entier qui doit etre positif, pas chaque cote separement.
  pure function fbd_delta2(lamc, epsc, pc, lamd, epsd, pd, pii) result(d2)
    implicit none
    real(kind=DOUBLE), intent(in) :: lamc, epsc, pc, lamd, epsd, pd, pii
    real(kind=DOUBLE) :: d2
    d2 = lamc*epsc + lamd*epsd &
         - (pii - pc)**2/(2.0_DOUBLE*lamc) - (pii - pd)**2/(2.0_DOUBLE*lamd)
  end function fbd_delta2

  ! Vitesse d'onde requise par la positivite pour une maille, etant donnes le D_p
  ! et le saut de vitesse Du du NOEUD. Racine positive de la condition
  ! 2 eps lam^2 + Du (p - D_p) lam - D_p^2/2 ... >= 0, d'ou la forme ci-dessous.
  ! Si D_p = 0 (etat au repos, ou champ a divergence nulle) alors lam_req = 0 et
  ! lambda_w vaut EXACTEMENT rho a : c'est la non-regression, par construction et
  ! pas par reglage.
  ! Vitesse d'onde requise par la positivite de 'vp_pp_psvp' : racine positive de
  !     A lam^2 - B lam - C >= 0,  A = eps + (u_p - u_c) u_p,  B = p u_p - p_p u_c,
  !     C = (p_p - p)^2 / 2.
  ! Rend -1 quand A <= 0 : la faisabilite est perdue et AUCUN lambda ne convient.
  ! L'appelant doit distinguer ce cas d'un lam_req simplement grand, sinon une
  ! INFAISABILITE se lit comme une exigence severe.
  pure function psvp_lambda(epsc, pc, pp_n, uc, up, lam_a, marge) result(lamw)
    implicit none
    real(kind=DOUBLE), intent(in) :: epsc, pc, pp_n, uc, up, lam_a, marge
    real(kind=DOUBLE) :: lamw, aa, bb, cc, disc
    aa = epsc + (up - uc)*up
    if (aa <= 0.0_DOUBLE) then
      lamw = -1.0_DOUBLE
      return
    end if
    bb = pc*up - pp_n*uc
    cc = 0.5_DOUBLE*(pp_n - pc)**2
    disc = bb*bb + 4.0_DOUBLE*aa*cc
    if (disc < 0.0_DOUBLE) then
      lamw = lam_a
      return
    end if
    lamw = max(lam_a, marge*(bb + sqrt(disc))/(2.0_DOUBLE*aa))
  end function psvp_lambda

  pure function fbw_lambda(tau, pk, gk, epsk, dp, du, marge) result(lamw)
    implicit none
    real(kind=DOUBLE), intent(in) :: tau, pk, gk, epsk, dp, du, marge
    real(kind=DOUBLE) :: lamw, lam_a, pcm, disc, lam_req
    lam_a = sqrt(max(0.0_DOUBLE, gk*pk*tau))/tau
    lamw = lam_a
    if (epsk <= 0.0_DOUBLE) return
    pcm = pk - dp
    disc = du*du*pcm*pcm + 8.0_DOUBLE*epsk*dp*dp
    if (disc < 0.0_DOUBLE) return
    lam_req = (-du*pcm + sqrt(disc))/(4.0_DOUBLE*epsk)
    lamw = max(lam_a, marge*lam_req)
  end function fbw_lambda

  ! --- Distribution de A(t)/A(0) par maille, ou A_c est la SOMME DES AIRES DE
  ! SOUS-FACES de la maille -- exactement la quantite qui entre dans la borne
  ! acoustique `m_c / somme(A lambda)`. On la calcule dans le code plutot que de
  ! re-deriver des aires depuis les coordonnees d'un VTU : c'est la meme somme que
  ! celle du pas de temps, au bit pres, donc la mesure et la borne parlent du meme
  ! objet. On sort une DISTRIBUTION (med, p90, p99, max et des comptes) parce qu'un
  ! etirement concentre sur 1 % du maillage et un etirement uniforme donnent la meme
  ! somme totale et ne sont pas le meme phenomene.
  subroutine area_of_cells(mesh, a)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(:), intent(out) :: a
    integer(kind=ENTIER) :: i, j, k, ise, isf
    a = 0.0_DOUBLE
    do i = 1, mesh%n_elems
      do j = 1, mesh%elem(i)%n_sub_elems
        ise = mesh%elem(i)%sub_elem(j)
        do k = 1, mesh%sub_elem(ise)%n_sub_faces
          isf = mesh%sub_elem(ise)%sub_face(k)
          a(i) = a(i) + mesh%sub_face(isf)%area
        end do
      end do
    end do
  end subroutine area_of_cells

  subroutine area_report(mesh, iter)
    implicit none
    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: iter
    real(kind=DOUBLE), dimension(:), allocatable :: at, r, fs
    integer(kind=ENTIER) :: i, n, nf
    real(kind=DOUBLE) :: st, s0
    if (.not. allocated(area_a0)) return
    allocate(at(mesh%n_elems), r(mesh%n_elems), fs(mesh%n_elems))
    call area_of_cells(mesh, at)
    n = 0; nf = 0; st = 0.0_DOUBLE; s0 = 0.0_DOUBLE
    do i = 1, mesh%n_elems
      if (mesh%elem(i)%is_ghost) cycle
      if (area_a0(i) <= 0.0_DOUBLE) cycle
      n = n + 1; r(n) = at(i)/area_a0(i)
      st = st + at(i); s0 = s0 + area_a0(i)
      ! Facteur de FORME adimensionnel A^(3/2)/V : il voit une degenerescence que
      ! ni A/A0 ni V/V0 pris separement ne montrent, parce qu'une contraction
      ! isotrope laisse A^(3/2)/V invariant alors qu'elle ecroule les deux rapports.
      if (abs(mesh%elem(i)%volume) > 0.0_DOUBLE .and. at(i) > 0.0_DOUBLE) then
        nf = nf + 1
        fs(nf) = at(i)**1.5_DOUBLE/abs(mesh%elem(i)%volume)
      end if
    end do
    if (n < 1) then; deallocate(at, r, fs); return; end if
    call tri_croissant(r, n)
    ! LES DEUX QUEUES. La premiere version de ce rapport ne sortait que la queue
    ! SUPERIEURE -- med, p90, p99, max et des comptes au-dessus de 2, 5, 12 --
    ! parce qu'elle cherchait de l'etirement. Sur noh_2d le signal est entierement
    ! dans la queue INFERIEURE : 4 mailles sur 2500 a A/A0 = 0.0024, invisibles au
    ! max comme a la somme (4/2500 = 0.16 % de deplacement). Des que la
    ! degenerescence peut aller dans les deux sens, il faut les deux queues.
    print*, "AIRE pas", iter, " n", n, " sommeA(t)/sommeA(0)", st/s0, &
            " | A/A0 BASSE : min", r(1), " p1", r(1 + n/100), " p10", r(1 + n/10), &
            " | med", r(1 + n/2), &
            " | HAUTE : p90", r(1 + (9*n)/10), " p99", r(1 + (99*n)/100), " max", r(n), &
            " | <0.5 :", count(r(1:n) < 0.5_DOUBLE), " <0.1 :", &
            count(r(1:n) < 0.1_DOUBLE), " <0.01 :", count(r(1:n) < 0.01_DOUBLE), &
            " >2 :", count(r(1:n) > 2.0_DOUBLE), " >5 :", count(r(1:n) > 5.0_DOUBLE), &
            " >12 :", count(r(1:n) > 12.0_DOUBLE)
    if (nf > 1) then
      call tri_croissant(fs, nf)
      print*, "FORME pas", iter, " A^(3/2)/V sur", nf, " mailles : min", fs(1), &
              " p1", fs(1 + nf/100), " med", fs(1 + nf/2), &
              " p99", fs(1 + (99*nf)/100), " max", fs(nf)
    end if
    ! DEGRADATION PAR MAILLE, et non rapport de deux agregats. Diviser le max
    ! global a t par le max global a l'initial n'est licite que sur un maillage
    ! UNIFORME, ou toutes les mailles partent du meme nombre. Sur le maillage
    ! triangulaire l'etendue initiale du facteur de forme est de x1.705 (min 73.59,
    ! max 125.45) et les deux max ne sont pas la meme maille : le rapport des
    ! agregats peut alors se tromper d'exactement cette etendue. On garde donc le
    ! facteur initial PAR MAILLE et on sort la distribution du rapport.
    if (allocated(forme_a0)) then
      nf = 0
      do i = 1, mesh%n_elems
        if (mesh%elem(i)%is_ghost) cycle
        if (forme_a0(i) <= 0.0_DOUBLE) cycle
        if (abs(mesh%elem(i)%volume) <= 0.0_DOUBLE .or. at(i) <= 0.0_DOUBLE) cycle
        nf = nf + 1
        fs(nf) = (at(i)**1.5_DOUBLE/abs(mesh%elem(i)%volume))/forme_a0(i)
      end do
      if (nf > 1) then
        call tri_croissant(fs, nf)
        print*, "FORMER pas", iter, " forme(t)/forme(0) PAR MAILLE sur", nf, &
                " : min", fs(1), " p1", fs(1 + nf/100), " med", fs(1 + nf/2), &
                " p75", fs(1 + (3*nf)/4), " p99", fs(1 + (99*nf)/100), " max", fs(nf), &
                " | >2 :", count(fs(1:nf) > 2.0_DOUBLE), &
                " >5 :", count(fs(1:nf) > 5.0_DOUBLE)
      end if
    end if
    deallocate(at, r, fs)
  end subroutine area_report

  ! tri par insertion sur une copie : n_elems vaut au plus quelques dizaines de
  ! milliers et le rapport n'est emis que tous les N pas, donc O(n^2) serait trop
  ! lent -- on prend un tri par tas, en place et sans allocation.
  subroutine tri_croissant(v, n)
    implicit none
    real(kind=DOUBLE), dimension(:), intent(inout) :: v
    integer(kind=ENTIER), intent(in) :: n
    integer(kind=ENTIER) :: i, j, k
    real(kind=DOUBLE) :: t
    do i = n/2, 1, -1
      call tamiser(v, i, n)
    end do
    do i = n, 2, -1
      t = v(1); v(1) = v(i); v(i) = t
      call tamiser(v, 1, i-1)
    end do
    j = 0; k = 0; t = 0.0_DOUBLE
  end subroutine tri_croissant

  subroutine tamiser(v, d, f)
    implicit none
    real(kind=DOUBLE), dimension(:), intent(inout) :: v
    integer(kind=ENTIER), intent(in) :: d, f
    integer(kind=ENTIER) :: r, c
    real(kind=DOUBLE) :: t
    r = d
    do
      c = 2*r
      if (c > f) exit
      if (c + 1 <= f) then
        if (v(c+1) > v(c)) c = c + 1
      end if
      if (v(r) >= v(c)) exit
      t = v(r); v(r) = v(c); v(c) = t
      r = c
    end do
  end subroutine tamiser

  pure function fbd_pbar(pl, pr, zl, zr) result(pb)
    implicit none
    real(kind=DOUBLE), intent(in) :: pl, pr, zl, zr
    real(kind=DOUBLE) :: pb
    pb = (pl/zl + pr/zr)/(1.0_DOUBLE/zl + 1.0_DOUBLE/zr)
  end function fbd_pbar

  pure subroutine fbd_accum(area, zl, zr, vnl, vnr, num, den)
    implicit none
    real(kind=DOUBLE), intent(in) :: area, zl, zr, vnl, vnr
    real(kind=DOUBLE), intent(inout) :: num, den
    num = num + area*(vnr - vnl)
    den = den + area*(1.0_DOUBLE/zl + 1.0_DOUBLE/zr)
  end subroutine fbd_accum

  ! Les trois controles, en dur, sur des etats et une geometrie fabriques.
  ! alpha_c = 1 - dt somme(A lambda_schema)/m_c, eq:cflexp de la these.
  ! Par PAS et jamais en cumule : un compteur cumule se remplit apres la
  ! divergence et ne dit plus a quel moment la garantie a cede.
  subroutine alpha_report(mesh, mass, dt, iter, tau, pre, gam)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(:), intent(in) :: mass, tau, pre, gam
    real(kind=DOUBLE), intent(in) :: dt
    integer(kind=ENTIER), intent(in) :: iter
    integer(kind=ENTIER) :: i, nneg, n, imin, j, k, ise, isf
    real(kind=DOUBLE) :: a, amin, sdt, sar, cc
    if (.not. allocated(alpha_sal)) return
    nneg = 0; n = 0; amin = 1.0e300_DOUBLE; imin = -1
    do i = 1, mesh%n_elems
      if (mesh%elem(i)%is_ghost) cycle
      if (mass(i) <= 0.0_DOUBLE) cycle
      a = 1.0_DOUBLE - dt*alpha_sal(i)/mass(i)
      if (a /= a) cycle
      n = n + 1
      if (a < amin) then; amin = a; imin = i; end if
      if (a < 0.0_DOUBLE) nneg = nneg + 1
    end do
    ! Et pour la maille FAUTIVE, on compare ma somme a celle que `compute_dt`
    ! construit lui-meme, par SA traversee (sous-elements -> sous-faces) et non
    ! par les voisins de face. Si les deux diffèrent, c'est mon jeu de sous-faces
    ! qui est faux et non la CFL : la borne 1 - cfl*r_c du theoreme interdit
    ! alpha_c < 1 - cfl, donc tout negatif a cfl = 0.95 est impossible.
    if (imin > 0) then
      sdt = 0.0_DOUBLE; sar = 0.0_DOUBLE
      cc = sqrt(max(0.0_DOUBLE, gam(imin)*pre(imin)*tau(imin)))
      do j = 1, mesh%elem(imin)%n_sub_elems
        ise = mesh%elem(imin)%sub_elem(j)
        do k = 1, mesh%sub_elem(ise)%n_sub_faces
          isf = mesh%sub_elem(ise)%sub_face(k)
          sdt = sdt + mesh%sub_face(isf)%area*(cc/tau(imin))
          sar = sar + mesh%sub_face(isf)%area
        end do
      end do
      print*, "ALPHAD pas", iter, " maille", imin, " MA somme(A lam)", alpha_sal(imin), &
        " somme(A rho a) facon compute_dt", sdt, " rapport", &
        alpha_sal(imin)/max(sdt, 1.0e-300_DOUBLE), " somme(A) de compute_dt", sar
    end if
    if (n > 0) print*, "ALPHA pas", iter, " mailles", n, " alpha_c < 0 :", nneg, &
      " alpha_c MINIMAL", amin
  end subroutine alpha_report

  subroutine fbd_step_reset()
    implicit none
    fbw_st_ff = 0; fbw_st_cf = 0; fbw_st_cc = 0
    if (allocated(fbw_cviol)) fbw_cviol = 0
    fbw_st_d2min = 1.0e300_DOUBLE; fbw_st_rmax = 0.0_DOUBLE
    fbd_st_dpmax = 0.0_DOUBLE; fbd_st_d2min = 1.0e300_DOUBLE
    fbd_st_zbad = 0; fbd_st_void = 0
    fbd_st_ff = 0; fbd_st_cf = 0; fbd_st_cc = 0
  end subroutine fbd_step_reset

  subroutine fbd_selfcheck()
    implicit none
    real(kind=DOUBLE) :: g, pl, pr, rl, rr, ul, ur, zl, zr
    real(kind=DOUBLE) :: num, den, dp, pi_c, pi_ref, pb
    real(kind=DOUBLE), dimension(3) :: nk, ck, ckm, vv
    real(kind=DOUBLE), dimension(2, 7) :: vx
    real(kind=DOUBLE), dimension(2, 6) :: cen
    real(kind=DOUBLE) :: th0, rad, ak, viol, pcold, phot
    integer(kind=ENTIER) :: k, kc, km, kf, kp, nviol
    character(len=22) :: nom

    g = 1.4_DOUBLE
    print*, "FBD --- controle 1 : limite 1D = solveur a DEUX CHOCS"
    ! La cible n'est PLUS l'acoustique. L'accord avec l'acoustique ne doit tenir que
    ! dans la limite Dp -> 0, Du -> 0, ou le max a trois termes se reduit a rho a.
    ! La reference a deux chocs est ecrite par le chemin des INVARIANTS de Riemann
    ! (v* puis Pi = p_l - lambda_l (v* - u_l)), c'est-a-dire la forme de
    ! ns_euler_rs_module et non une refactorisation de pbar_f - D_p.
    do k = 1, 5
      select case (k)
      case (1); nom = "Sod           "; rl=1.0_DOUBLE; pl=1.0_DOUBLE
                rr=0.125_DOUBLE; pr=0.1_DOUBLE; ul=0.0_DOUBLE; ur=0.0_DOUBLE
      case (2); nom = "choc 1e4      "; rl=1.0_DOUBLE; pl=1.0e4_DOUBLE
                rr=1.0_DOUBLE; pr=1.0_DOUBLE; ul=0.0_DOUBLE; ur=0.0_DOUBLE
      case (3); nom = "compression   "; rl=1.0_DOUBLE; pl=1.0_DOUBLE
                rr=1.0_DOUBLE; pr=1.0_DOUBLE; ul=1.0_DOUBLE; ur=-1.0_DOUBLE
      case (4); nom = "bas Mach Ma=1e-3"; rl=1.0_DOUBLE; pl=1.0_DOUBLE
                rr=1.0_DOUBLE; pr=1.0_DOUBLE
                ul=1.183216e-3_DOUBLE; ur=-1.183216e-3_DOUBLE
      case (5); nom = "contact Sedov "; rl=1.0_DOUBLE; pl=2359.484_DOUBLE
                rr=1.0_DOUBLE; pr=1.0e-12_DOUBLE; ul=0.0_DOUBLE; ur=0.0_DOUBLE
      end select
      zl = fbd_lambda(1.0_DOUBLE/rl, pl, g, pr, ur - ul)
      zr = fbd_lambda(1.0_DOUBLE/rr, pr, g, pl, ur - ul)
      num = 0.0_DOUBLE; den = 0.0_DOUBLE
      call fbd_accum(1.0_DOUBLE, zl, zr, ul, ur, num, den)
      pi_c = fbd_pbar(pl, pr, zl, zr) - num/den
      ! reference INDEPENDANTE, par les invariants
      pb = (zl*ul + zr*ur - (pr - pl))/(zl + zr)
      pi_ref = pl - zl*(pb - ul)
      ! et l'acoustique PURE, pour montrer l'ecart la ou il compte
      th0 = sqrt(g*pl*rl); rad = sqrt(g*pr*rr)
      print*, "FBD   ", nom, " Pi", pi_c, " deux chocs", pi_ref, &
              " ecart", abs(pi_c - pi_ref), " acoustique pure", &
              (rad*pl + th0*pr - th0*rad*(ur - ul))/(th0 + rad)
    end do

    print*, "FBD --- controle 2 : D_p sur un eventail de 6 triangles"
    do kc = 1, 2
      ! Eventail : noeud a l'origine, 6 sommets. Les rayons et angles sont
      ! PERTURBES au second passage, pour que l'exactitude ne tienne pas a la
      ! symetrie du maillage.
      do k = 1, 6
        th0 = 2.0_DOUBLE*3.14159265358979324_DOUBLE*real(k-1, DOUBLE)/6.0_DOUBLE
        rad = 1.0_DOUBLE
        if (kc == 2) then
          th0 = th0 + 0.17_DOUBLE*sin(3.0_DOUBLE*real(k, DOUBLE))
          rad = 1.0_DOUBLE + 0.31_DOUBLE*cos(2.0_DOUBLE*real(k, DOUBLE))
        end if
        vx(1, k) = rad*cos(th0); vx(2, k) = rad*sin(th0)
      end do
      vx(:, 7) = vx(:, 1)
      ! centre de la maille k = centroide de (noeud, sommet k, sommet k+1).
      ! PIEGE : il faut la vraie position du centroide, pas un angle suppose.
      do k = 1, 6
        cen(:, k) = (vx(:, k) + vx(:, k+1))/3.0_DOUBLE
      end do
      do kf = 1, 4
        num = 0.0_DOUBLE; den = 0.0_DOUBLE
        do k = 1, 6
          ! sous-face k : demi-arete du noeud au milieu du sommet k.
          ! aire = |v_k|/2 (epaisseur d'extrusion 1), normale perpendiculaire a v_k.
          ak = 0.5_DOUBLE*sqrt(vx(1,k)**2 + vx(2,k)**2)
          nk = 0.0_DOUBLE
          nk(1) = -vx(2,k); nk(2) = vx(1,k)
          nk = nk/norm2(nk)
          km = 1 + mod(k - 2 + 6, 6)          ! maille a gauche de la sous-face k
          kp = k                              ! maille a droite
          ckm = 0.0_DOUBLE; ckm(1:2) = cen(:, km)
          ck  = 0.0_DOUBLE; ck(1:2)  = cen(:, kp)
          call fbd_accum(ak, 1.0_DOUBLE, 1.0_DOUBLE, &
            dot_product(fbd_champ(kf, ckm), nk), &
            dot_product(fbd_champ(kf, ck), nk), num, den)
        end do
        select case (kf)
        case (1); nom = "translation           "
        case (2); nom = "rotation solide       "
        case (3); nom = "cisaillement lineaire "
        case (4); nom = "COMPRESSION (div = 2) "
        end select
        if (kc == 1) then
          print*, "FBD   regulier  ", nom, " D_p =", num/den
        else
          print*, "FBD   perturbe  ", nom, " D_p =", num/den
        end if
      end do
    end do

    print*, "FBD --- controle 3 : bilan de positivite ENTIER (delta2)"
    ! Plus de critere par cote : delta2 = lam_c eps_c + lam_d eps_d
    !                                    - (Pi-p_c)^2/(2 lam_c) - (Pi-p_d)^2/(2 lam_d)
    ! Le terme du cote CHAUD paie pour le deficit du cote froid, ce qu'un critere
    ! par cote ne peut pas voir.
    do k = 1, 4
      select case (k)
      case (1); nom = "froid-froid           "; pl=1.0e-12_DOUBLE; pr=1.0e-12_DOUBLE
      case (2); nom = "contact Sedov 2.4e3   "; pl=2359.484_DOUBLE; pr=1.0e-12_DOUBLE
      case (3); nom = "contact 1e6           "; pl=1.0e6_DOUBLE;    pr=1.0e-12_DOUBLE
      case (4); nom = "contact 1e18          "; pl=1.0e18_DOUBLE;   pr=1.0e-12_DOUBLE
      end select
      rl = 1.0_DOUBLE; rr = 1.0_DOUBLE
      zl = fbd_lambda(1.0_DOUBLE/rl, pl, g, pr, 0.0_DOUBLE)
      zr = fbd_lambda(1.0_DOUBLE/rr, pr, g, pl, 0.0_DOUBLE)
      num = 0.0_DOUBLE; den = 0.0_DOUBLE
      call fbd_accum(1.0_DOUBLE, zl, zr, 0.0_DOUBLE, 0.0_DOUBLE, num, den)
      pi_c = fbd_pbar(pl, pr, zl, zr) - num/den
      print*, "FBD   ", nom, " lambda_froid", zr, " Pi", pi_c, " delta2", &
        fbd_delta2(zl, pl/(g - 1.0_DOUBLE)/rl, pl, zr, pr/(g - 1.0_DOUBLE)/rr, pr, pi_c)
    end do

  end subroutine fbd_selfcheck

  ! Champs de vitesse du controle 2. Les trois premiers sont a divergence NULLE,
  ! le quatrieme a divergence 2 et doit donner un D_p non nul.
  pure function fbd_champ(kf, x) result(v)
    implicit none
    integer(kind=ENTIER), intent(in) :: kf
    real(kind=DOUBLE), dimension(3), intent(in) :: x
    real(kind=DOUBLE), dimension(3) :: v
    v = 0.0_DOUBLE
    select case (kf)
    case (1); v(1) = 0.7_DOUBLE; v(2) = -0.3_DOUBLE
    case (2); v(1) = -x(2);      v(2) = x(1)
    case (3); v(1) = x(2);       v(2) = 0.0_DOUBLE
    case (4); v(1) = x(1);       v(2) = x(2)
    end select
  end function fbd_champ

  ! ------------------------------------------------------------------
  ! Controle du mode HOURGLASS (sablier). Ajoute apres la mesure du
  ! 2026-10-08 sur gresho_quad a t = 0.20 : la projection du champ
  ! nodal sur le motif de coin sigma = (+1,-1,+1,-1), orthogonal aux
  ! translations, aux rotations et aux deformations lineaires, donne
  !
  !     vp_pp    max|Vp| 0.963   H med 0.030   H p75 0.138   68.8 -> 11.0 %
  !     vp_fbw   max|Vp| 71.31   H med 0.417   H p75 0.593   68.8 % > 0.3
  !
  ! soit quatorze fois plus de sablier pour vp_fbw, a vortex VIVANT dans
  ! les deux cas (pic exact 1). Le sablier n'est donc pas le prix du
  ! maintien du vortex. Et sur sedov_quad, ou classic est vivant,
  ! classic et vp_fbw ont le MEME sablier (H med 0.0073 contre 0.0049) :
  ! ce n'est pas une propriete du schema, c'est Gresho qui le revele.
  !
  ! COMPLEMENTAIRE du controle 2 de fbd_selfcheck, qui verifie D_p sur
  ! quatre champs LISSES echantillonnes aux centres de maille. Un champ
  ! lisse ne peut pas exhiber le sablier, qui est un mode NODAL a
  ! l'echelle de la maille, v_p = (-1)^(i+j) d ; et le controle 2 ne
  ! regarde que D_p, jamais la force nodale.
  !
  ! Les quatre controles ci-dessous donnent les raisons pour lesquelles
  ! aucun mecanisme de la fermeture pbar_f - D_p ne voit ce mode :
  !
  !  1  somme_f A_f n_f = 0 au noeud, sur maillage regulier ET perturbe.
  !     Donc tout Pi CONSTANT sur les sous-faces du noeud produit une
  !     force nodale nulle. Comme fbd_dp est calcule UNE SEULE FOIS par
  !     noeud dans compute_rhs_around_node et applique a toutes ses
  !     sous-faces, la part D_p de Pi = pbar_f - D_p est annulee par la
  !     fermeture QUEL QUE SOIT son contenu : D_p ne peut rien faire sur
  !     le bilan de quantite de mouvement NODAL.
  !
  !  2  mais D_p n'est PAS inerte sur la force de MAILLE, et c'est une
  !     nuance qui a failli etre ecrite a l'envers. Les sous-faces d'une
  !     maille sont visitees sous des noeuds DIFFERENTS, donc avec des
  !     D_p differents : D_p varie a travers les sous-faces d'une maille
  !     et contribue a sa force. Le controle 2 le montre en donnant un
  !     D_p distinct par noeud. Corollaire : la force de maille nulle
  !     mesuree sur le piston de Saltzman (|F| = 2.09e-25 contre
  !     1.29e-06 pour classic) ne vient PAS de la fermeture, elle vient
  !     de ce que Pi y etait reellement uniforme -- pbar_f = p_c par
  !     champ de pression uniforme et etat miroir, ET D_p = 0 a chaque
  !     noeud parce que tout est au repos et que le numerateur de D_p
  !     est bati sur les vitesses de MAILLE. Qui voudra reparer en
  !     jouant sur D_p ne doit pas croire que la maille est hors
  !     d'atteinte : elle ne l'est pas.
  !
  !  3  le mode sablier est a dV/dt nul EXACTEMENT, et le controle le
  !     donne nul sur maillage regulier ET sur maillage perturbe. Ce
  !     n'est pas un accident du maillage regulier mais une identite
  !     GEOMETRIQUE, vraie arete par arete sur n'importe quelle
  !     geometrie : la contribution d'une arete a dV/dt vaut
  !     (|e|/2)(v_a + v_b).n, et les deux extremites d'une arete de
  !     quadrangle ont des parites de damier OPPOSEES, donc v_a + v_b = 0
  !     et chaque arete contribue zero separement. J'avais predit le
  !     contraire -- un residu non nul sur maillage perturbe -- et la
  !     mesure dit mieux. PORTEE, et elle est etroite : l'argument
  !     demande que les deux bouts de chaque arete aient des parites
  !     opposees, c'est-a-dire que le graphe des noeuds soit BIPARTI.
  !     Un maillage de quadrangles l'est ; un maillage de TRIANGLES ne
  !     l'est pas (cycles de longueur 3), donc le damier n'y a pas
  !     d'analogue -- mesure sur un patch 6x6 maille des deux facons
  !     avec les MEMES noeuds : ||D v_damier|| = 0 exactement sur
  !     quadrangles contre 8.485 sur triangles. Donc ce controle ne
  !     dit rien de gresho_tri, qui est triangulaire et qui est le
  !     pire cas de la campagne.
  !     Le MECANISME survit pourtant au changement de maillage, seule
  !     la FORME du mode ne survit pas : D est de rang ligne plein sur
  !     les deux maillages, donc dim ker D = 2N - nc exactement (62 et
  !     26), et la direction la plus oscillante DANS ce noyau atteint
  !     2.0039 sur quadrangles et 2.0029 sur triangles contre 2.0057
  !     pour la plus oscillante de tout l'espace, en norme du laplacien
  !     de graphe -- soit 0.15 % pres. Une dissipation batie sur une
  !     divergence est donc aveugle aux modes maille-a-maille sur
  !     triangles AUSSI, simplement pas au damier mais a un autre
  !     element du meme noyau.
  !     Le temoin du controle sort dV/dt du champ de COMPRESSION v = x
  !     a cote, qui vaut 2 sur maillage regulier : sans lui, un zero ne
  !     prouverait pas que la mesure sait voir autre chose que zero.
  !     Consequence non verifiable ici mais mesuree : en lagrangien la
  !     pression repond au VOLUME, donc un mode a volume constant ne
  !     declenche aucune reponse thermodynamique -- sur gresho_quad,
  !     eps_min va de 17857 a 17368 pendant que max|Vp| monte a 71.3.
  !
  !  4  la force nodale sur le mode, a pression uniforme et v_c = 0 :
  !     classic donne lambda * M_p * v_p, non nulle, tandis que la part
  !     D_p de vp_fbw donne zero machine. Le terme de classic,
  !     Pi_pcf = p_c - lambda_c (v_p - v_c).n_pcf, echappe aux deux
  !     pieges a la fois : il VARIE d'une sous-face a l'autre donc il
  !     survit a la fermeture, et il n'est pas une divergence donc le
  !     sablier n'est pas dans son noyau. Il donne
  !        somme_f A_f ((v_p - v_c).n_f) n_f = M_p (v_p - v_c)/lambda,
  !     forme DEFINIE POSITIVE sur tout l'espace des vitesses nodales :
  !     aucun mode n'y est libre, sablier compris. C'est ce terme qu'on
  !     a perdu, et c'est le meme qui donne s_c = 0 dans le theoreme de
  !     positivite.
  ! ------------------------------------------------------------------
  subroutine hourglass_selfcheck()
    implicit none
    integer(kind=ENTIER), parameter :: NCD = 4            ! mailles par direction
    integer(kind=ENTIER), parameter :: NND = NCD + 1      ! noeuds par direction
    real(kind=DOUBLE), dimension(2, NND*NND) :: xn, vn
    real(kind=DOUBLE), dimension(NND*NND) :: dpn
    integer(kind=ENTIER), dimension(4, NCD*NCD) :: cel
    integer(kind=ENTIER) :: i, j, k, kc, ke, ic, ip, kg, kp
    real(kind=DOUBLE), dimension(2) :: ea, eb, nk, sn, fcl, ffb, fmp, fce
    real(kind=DOUBLE), dimension(2, 2) :: mpn
    real(kind=DOUBLE) :: ak, ll, dvm, dvc, lam, pun, dmx, pert
    character(len=12) :: nom

    lam = 2.5_DOUBLE          ! non unitaire, pour que l'egalite lam*M_p*v_p se voie
    pun = 1.0_DOUBLE          ! pression uniforme

    do kg = 1, 2              ! kg = 1 regulier, kg = 2 perturbe
      pert = 0.0_DOUBLE
      if (kg == 2) pert = 0.23_DOUBLE
      nom = "regulier    "
      if (kg == 2) nom = "perturbe    "

      ! --- geometrie : patch NCD x NCD de carres unite, noeuds interieurs
      ! deplaces d'un motif deterministe quand kg = 2. Les noeuds de BORD
      ! restent en place, sinon le patch cesse d'etre un patch.
      do j = 1, NND
        do i = 1, NND
          ip = (j - 1)*NND + i
          xn(1, ip) = real(i - 1, DOUBLE)
          xn(2, ip) = real(j - 1, DOUBLE)
          if (i > 1 .and. i < NND .and. j > 1 .and. j < NND) then
            xn(1, ip) = xn(1, ip) + pert*sin(1.7_DOUBLE*real(i, DOUBLE) &
                                           + 2.3_DOUBLE*real(j, DOUBLE))
            xn(2, ip) = xn(2, ip) + pert*cos(2.9_DOUBLE*real(i, DOUBLE) &
                                           - 1.1_DOUBLE*real(j, DOUBLE))
          end if
          ! mode SABLIER nodal : damier, (+1,-1,+1,-1) sur les coins de
          ! CHAQUE maille simultanement.
          vn(1, ip) = real((-1)**(i + j), DOUBLE)
          vn(2, ip) = 0.0_DOUBLE
          ! un D_p DIFFERENT par noeud, pour le controle 2
          dpn(ip) = 0.37_DOUBLE + 0.11_DOUBLE*real(i, DOUBLE) &
                                - 0.07_DOUBLE*real(j, DOUBLE)
        end do
      end do
      do j = 1, NCD
        do i = 1, NCD
          ic = (j - 1)*NCD + i
          cel(1, ic) = (j - 1)*NND + i
          cel(2, ic) = (j - 1)*NND + i + 1
          cel(3, ic) =  j     *NND + i + 1
          cel(4, ic) =  j     *NND + i
        end do
      end do

      ! --- controle 1 : fermeture au NOEUD. somme A n sur toutes les
      ! sous-faces qui touchent un noeud donne, cellules confondues.
      dmx = 0.0_DOUBLE
      do ip = 1, NND*NND
        i = 1 + mod(ip - 1, NND); j = 1 + (ip - 1)/NND
        if (i == 1 .or. i == NND .or. j == 1 .or. j == NND) cycle
        sn = 0.0_DOUBLE
        do ic = 1, NCD*NCD
          do kc = 1, 4
            if (cel(kc, ic) /= ip) cycle
            do ke = 1, 2
              k = kc
              if (ke == 2) k = 1 + mod(kc - 2 + 4, 4)
              ea = xn(:, cel(k, ic)); eb = xn(:, cel(1 + mod(k, 4), ic))
              ll = sqrt((eb(1) - ea(1))**2 + (eb(2) - ea(2))**2)
              nk(1) = (eb(2) - ea(2))/ll; nk(2) = -(eb(1) - ea(1))/ll
              sn = sn + 0.5_DOUBLE*ll*nk
            end do
          end do
        end do
        dmx = max(dmx, sqrt(sn(1)**2 + sn(2)**2))
      end do
      print*, "HG  ", nom, "controle 1  max|somme A n| au noeud interieur =", dmx

      ! --- controle 2 : la force de MAILLE voit-elle un D_p par noeud ?
      ! Pi_pcf = pbar_f - D_p(noeud de la sous-face), pression uniforme.
      ! Si D_p etait inerte a la maille, fce serait nulle ; elle ne l'est pas.
      ic = 1 + (NCD*NCD)/2
      fce = 0.0_DOUBLE; sn = 0.0_DOUBLE
      do kc = 1, 4
        do ke = 1, 2
          k = kc
          if (ke == 2) k = 1 + mod(kc - 2 + 4, 4)
          ea = xn(:, cel(k, ic)); eb = xn(:, cel(1 + mod(k, 4), ic))
          ll = sqrt((eb(1) - ea(1))**2 + (eb(2) - ea(2))**2)
          nk(1) = (eb(2) - ea(2))/ll; nk(2) = -(eb(1) - ea(1))/ll
          ak = 0.5_DOUBLE*ll
          sn  = sn  + ak*nk
          fce = fce - ak*(pun - dpn(cel(kc, ic)))*nk
        end do
      end do
      print*, "HG  ", nom, "controle 2  |somme A n| sur la maille =", &
        sqrt(sn(1)**2 + sn(2)**2), " |F_maille| avec un D_p par noeud =", &
        sqrt(fce(1)**2 + fce(2)**2)

      ! --- controle 3 : dV/dt de la maille sous le mode sablier,
      ! compare au champ de COMPRESSION v = x, de divergence 2.
      dvm = 0.0_DOUBLE; dvc = 0.0_DOUBLE
      do ic = 1, NCD*NCD
        do kc = 1, 4
          do ke = 1, 2
            k = kc
            if (ke == 2) k = 1 + mod(kc - 2 + 4, 4)
            ea = xn(:, cel(k, ic)); eb = xn(:, cel(1 + mod(k, 4), ic))
            ll = sqrt((eb(1) - ea(1))**2 + (eb(2) - ea(2))**2)
            nk(1) = (eb(2) - ea(2))/ll; nk(2) = -(eb(1) - ea(1))/ll
            ak = 0.5_DOUBLE*ll
            kp = cel(kc, ic)
            dvm = dvm + ak*(vn(1, kp)*nk(1) + vn(2, kp)*nk(2))
            dvc = dvc + ak*(xn(1, kp)*nk(1) + xn(2, kp)*nk(2))
          end do
        end do
        if (ic == NCD*NCD) then
          print*, "HG  ", nom, "controle 3  dV/dt SABLIER =", dvm, &
            " dV/dt compression (div=2) =", dvc
        end if
        dvm = 0.0_DOUBLE; dvc = 0.0_DOUBLE
      end do

      ! --- controle 4 : force nodale sur le mode, p uniforme et v_c = 0.
      ! classic : Pi = p_c - lam (v_p - v_c).n   -> F = lam M_p v_p
      ! vp_fbw  : Pi = pbar_f - D_p, D_p constant au noeud -> F = 0
      ip = (NND*NND + 1)/2
      fcl = 0.0_DOUBLE; ffb = 0.0_DOUBLE; mpn = 0.0_DOUBLE
      do ic = 1, NCD*NCD
        do kc = 1, 4
          if (cel(kc, ic) /= ip) cycle
          do ke = 1, 2
            k = kc
            if (ke == 2) k = 1 + mod(kc - 2 + 4, 4)
            ea = xn(:, cel(k, ic)); eb = xn(:, cel(1 + mod(k, 4), ic))
            ll = sqrt((eb(1) - ea(1))**2 + (eb(2) - ea(2))**2)
            nk(1) = (eb(2) - ea(2))/ll; nk(2) = -(eb(1) - ea(1))/ll
            ak = 0.5_DOUBLE*ll
            mpn(1, 1) = mpn(1, 1) + ak*nk(1)*nk(1)
            mpn(1, 2) = mpn(1, 2) + ak*nk(1)*nk(2)
            mpn(2, 1) = mpn(2, 1) + ak*nk(2)*nk(1)
            mpn(2, 2) = mpn(2, 2) + ak*nk(2)*nk(2)
            fcl = fcl - ak*(pun - lam*(vn(1, ip)*nk(1) + vn(2, ip)*nk(2)))*nk
            ffb = ffb - ak*(pun - dpn(ip))*nk
          end do
        end do
      end do
      fmp(1) = lam*(mpn(1, 1)*vn(1, ip) + mpn(1, 2)*vn(2, ip))
      fmp(2) = lam*(mpn(2, 1)*vn(1, ip) + mpn(2, 2)*vn(2, ip))
      print*, "HG  ", nom, "controle 4  |F classic| =", sqrt(fcl(1)**2 + fcl(2)**2), &
        " |lam M_p v_p| =", sqrt(fmp(1)**2 + fmp(2)**2), &
        " |F part D_p de vp_fbw| =", sqrt(ffb(1)**2 + ffb(2)**2)
    end do

  end subroutine hourglass_selfcheck

  ! ===========================================================================
  ! Plafond de positivite pour le terme de divergence de SOUS-ELEMENT.
  !
  ! Fermeture visee (annexe « Multiple divergences » du rapport) :
  !
  !     Pi_pcf = Pi_base_pcf - kappa_pc * delta_pc ,   kappa_pc = theta lambda_c h_pc
  !
  ! ou delta_pc est le taux de dilatation du sous-element de coin, calcule par
  ! kappa_subcell_div, qui depend de TOUS les noeuds bordant omega_pc.
  !
  ! L'admissibilite de l'etat intermediaire du probleme de Riemann de sous-face
  ! s'ecrit delta2(Pi) >= 0 avec
  !
  !     delta2(Pi) = lam_c eps_c + lam_d eps_d
  !                  - (Pi - p_c)^2/(2 lam_c) - (Pi - p_d)^2/(2 lam_d)
  !
  ! Concave en Pi, maximale a la moyenne ponderee par les impedances
  ! pbar_f = (lam_d p_c + lam_c p_d)/(lam_c + lam_d). Le carre complete donne
  !
  !     delta2(Pi) = delta2(pbar_f) - (Pi - pbar_f)^2/2 * (1/lam_c + 1/lam_d)
  !     delta2(pbar_f) = lam_c eps_c + lam_d eps_d - (p_c - p_d)^2/(2(lam_c+lam_d))
  !
  ! donc avec Pi = pbar_f - D_p - kappa delta_pc l'admissibilite est l'UNIQUE
  ! inegalite scalaire
  !
  !     ( D_p + kappa delta_pc )^2  <=  2 delta2(pbar_f) / (1/lam_c + 1/lam_d) = R^2
  !
  ! d'ou, en exigeant kappa >= 0, le plafond
  !
  !     kappa_max = ( R - D_p * sgn(delta_pc) ) / |delta_pc| .
  !
  ! CE N'EST PAS UN LIMITEUR : c'est la condition exacte d'admissibilite de l'etat
  ! intermediaire, lue comme une contrainte sur kappa. C'est la meme inegalite qui
  ! fixe lambda dans la construction de Gallice.
  !
  ! L'intervalle peut etre VIDE, et sa vacuite est un diagnostic : si
  ! delta2(pbar_f) <= 0, la fermeture de base est deja inadmissible AVANT toute
  ! correction (situation froid/froid, ou delta2 degenere en p^(3/2)). Aucun
  ! kappa ne repare ca, et le bon diagnostic est que c'est la base qu'il faut
  ! corriger. De meme si |D_p| > R : la base sort deja de l'intervalle.
  ! ===========================================================================


  ! -------------------------------------------------------------------------
  ! Plafond sur kappa_pc. kap_max n'a de sens que si status == MD_OK ;
  ! pour MD_NO_CONSTR il vaut HUGE (aucune contrainte), pour les deux autres
  ! il vaut zero et l'appelant doit traiter le cas, PAS l'ignorer.
  ! -------------------------------------------------------------------------
  pure subroutine md_kappa_cap(lam_c, eps_c, p_c, lam_d, eps_d, p_d, &
      dp_node, delta_pc, kap_max, d2_max, status)
    implicit none
    real(kind=DOUBLE), intent(in)  :: lam_c, eps_c, p_c
    real(kind=DOUBLE), intent(in)  :: lam_d, eps_d, p_d
    real(kind=DOUBLE), intent(in)  :: dp_node      ! D_p du noeud
    real(kind=DOUBLE), intent(in)  :: delta_pc     ! dilatation du sous-element
    real(kind=DOUBLE), intent(out) :: kap_max
    real(kind=DOUBLE), intent(out) :: d2_max       ! delta2(pbar_f), pour le diagnostic
    integer(kind=ENTIER), intent(out) :: status
    real(kind=DOUBLE) :: rr, sg

    kap_max = 0.0_DOUBLE
    d2_max  = lam_c*eps_c + lam_d*eps_d &
            - (p_c - p_d)**2/(2.0_DOUBLE*(lam_c + lam_d))

    if (lam_c <= 0.0_DOUBLE .or. lam_d <= 0.0_DOUBLE &
        .or. d2_max /= d2_max .or. d2_max <= 0.0_DOUBLE) then
      status = MD_EMPTY
      return
    end if

    rr = sqrt(2.0_DOUBLE*d2_max/(1.0_DOUBLE/lam_c + 1.0_DOUBLE/lam_d))

    if (abs(dp_node) > rr) then
      status = MD_BASE_OUT
      return
    end if

    ! Seuil ABSOLU sur delta_pc, et non relatif : une dilatation nulle ne
    ! contraint rien, mais un seuil relatif a rr ferait dependre le resultat
    ! de l'etat thermodynamique, ce qui n'a pas de sens ici.
    if (abs(delta_pc) <= tiny(1.0_DOUBLE)) then
      kap_max = huge(1.0_DOUBLE)
      status  = MD_NO_CONSTR
      return
    end if

    sg = sign(1.0_DOUBLE, delta_pc)
    kap_max = (rr - dp_node*sg)/abs(delta_pc)
    status  = MD_OK
  end subroutine md_kappa_cap

  ! -------------------------------------------------------------------------
  ! Controles en dur. Le controle de SERRAGE est le seul qui teste l'algebre
  ! de bout en bout : a kappa_max exactement, delta2 doit valoir zero, et de
  ! part et d'autre il doit changer de signe. Un plafond trop prudent passerait
  ! tous les autres controles sans que rien ne le signale.
  ! -------------------------------------------------------------------------
  subroutine md_posflux_selfcheck()
    implicit none
    real(kind=DOUBLE) :: lc, ec, pc, ld, ed, pd, dp, dl, km, d2m, pb, pi_v, a, b
    integer(kind=ENTIER) :: st, k
    character(len=22) :: nom

    print*, "MDPOS --- controle 1 : plafond, etat et delta2 atteint"
    do k = 1, 6
      dp = 0.0_DOUBLE; dl = 1.0_DOUBLE
      select case (k)
      case (1); nom = "Sod contact           "
        lc=1.1832_DOUBLE; ec=2.5_DOUBLE;  pc=1.0_DOUBLE
        ld=0.1323_DOUBLE; ed=2.0_DOUBLE;  pd=0.1_DOUBLE; dp=0.05_DOUBLE
      case (2); nom = "Sod, delta < 0        "
        lc=1.1832_DOUBLE; ec=2.5_DOUBLE;  pc=1.0_DOUBLE
        ld=0.1323_DOUBLE; ed=2.0_DOUBLE;  pd=0.1_DOUBLE; dp=0.05_DOUBLE; dl=-1.0_DOUBLE
      case (3); nom = "delta nul             "
        lc=1.1832_DOUBLE; ec=2.5_DOUBLE;  pc=1.0_DOUBLE
        ld=0.1323_DOUBLE; ed=2.0_DOUBLE;  pd=0.1_DOUBLE; dp=0.05_DOUBLE; dl=0.0_DOUBLE
      case (4); nom = "froid/froid           "
        lc=1.0_DOUBLE; ec=2.5e-12_DOUBLE; pc=1.0e-12_DOUBLE
        ld=1.0_DOUBLE; ed=2.5e-12_DOUBLE; pd=1.0e-12_DOUBLE
      case (5); nom = "contact Sedov         "
        lc=48.6_DOUBLE;   ec=2359.0_DOUBLE;  pc=2359.5_DOUBLE
        ld=1.18e-6_DOUBLE; ed=2.5e-12_DOUBLE; pd=1.0e-12_DOUBLE
      case (6); nom = "Gresho bas Mach       "
        lc=118.3_DOUBLE; ec=17857.0_DOUBLE; pc=7143.0_DOUBLE
        ld=118.3_DOUBLE; ed=17857.0_DOUBLE; pd=7143.0_DOUBLE; dp=0.01_DOUBLE
      end select
      call md_kappa_cap(lc, ec, pc, ld, ed, pd, dp, dl, km, d2m, st)
      pb = (ld*pc + lc*pd)/(lc + ld)
      if (st == MD_OK) then
        pi_v = pb - dp - km*dl
        print*, "MDPOS   ", nom, " kap_max", km, " delta2 atteint", &
          lc*ec + ld*ed - (pi_v-pc)**2/(2.0_DOUBLE*lc) - (pi_v-pd)**2/(2.0_DOUBLE*ld), &
          " delta2(pbar)", d2m
      else
        print*, "MDPOS   ", nom, " status", st, " delta2(pbar)", d2m
      end if
    end do

    print*, "MDPOS --- controle 2 : SYMETRIE (D_p,delta) -> (-D_p,-delta)"
    lc=1.1832_DOUBLE; ec=2.5_DOUBLE; pc=1.0_DOUBLE
    ld=0.1323_DOUBLE; ed=2.0_DOUBLE; pd=0.1_DOUBLE
    call md_kappa_cap(lc,ec,pc,ld,ed,pd,  0.05_DOUBLE,  1.0_DOUBLE, a, d2m, st)
    call md_kappa_cap(lc,ec,pc,ld,ed,pd, -0.05_DOUBLE, -1.0_DOUBLE, b, d2m, st)
    print*, "MDPOS   ", a, b, " ecart", abs(a-b)

    print*, "MDPOS --- controle 3 : SERRAGE, delta2 doit changer de signe en kap_max"
    call md_kappa_cap(lc,ec,pc,ld,ed,pd, 0.05_DOUBLE, 1.0_DOUBLE, km, d2m, st)
    pb = (ld*pc + lc*pd)/(lc + ld)
    do k = -1, 1
      pi_v = pb - 0.05_DOUBLE - km*(1.0_DOUBLE + real(k,DOUBLE)*1.0e-6_DOUBLE)
      print*, "MDPOS    eps =", real(k,DOUBLE)*1.0e-6_DOUBLE, " delta2 =", &
        lc*ec + ld*ed - (pi_v-pc)**2/(2.0_DOUBLE*lc) - (pi_v-pd)**2/(2.0_DOUBLE*ld)
    end do
  end subroutine md_posflux_selfcheck

  subroutine compute_rhs_around_node(mesh, sol, i_vert, nsfn, nsen, vp_i, vp_is_imposed_i, &
      dt, rhs, n_bc, bc_type, bc_val, b2d, mass, gamma_arr, &
      second_order, grad_v, grad_p, div_v)
    use linear_solver_module, only: lu_solve, print_mat, inverse_3_by_3
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    integer(kind=ENTIER), intent(in) :: i_vert, nsfn, nsen
    real(kind=DOUBLE), dimension(3), intent(inout) :: vp_i
    logical, intent(in) :: vp_is_imposed_i
    real(kind=DOUBLE), intent(in) :: dt
    real(kind=DOUBLE), dimension(5, nsen), intent(out) :: rhs
    integer(kind=ENTIER), intent(in) :: n_bc
    character(len=255), dimension(n_bc) :: bc_type
    real(kind=DOUBLE), dimension(5, n_bc) :: bc_val
    logical, intent(in) :: b2d
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mass
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    logical, intent(in) :: second_order
    real(kind=DOUBLE), dimension(3, 3, mesh%n_elems), intent(in) :: grad_v
    real(kind=DOUBLE), dimension(3, mesh%n_elems), intent(in) :: grad_p
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: div_v

    integer(kind=ENTIER) :: j, id_face, id_sub_face
    integer(kind=ENTIER) :: idl, idr, lse, lse_loc, rse, rse_loc
    real(kind=DOUBLE), dimension(2, nsfn) :: lambda
    real(kind=DOUBLE), dimension(5, 2, nsfn) :: sol_lr
    real(kind=DOUBLE), dimension(5) :: flux
    real(kind=DOUBLE), dimension(3, 3) :: Mp, Mp_inv
    real(kind=DOUBLE), dimension(3) :: Rp, Bp
    real(kind=DOUBLE) :: vl, vr, pl, pr, v_bar, pl_et, pr_et, v_et, Pp
    real(kind=DOUBLE) :: vn_cv, al_cv, ar_cv, mu_cv
    real(kind=DOUBLE) :: gl, gr, al, ar, taul, taur
    logical :: is_corner, have_first_wall_norm
    real(kind=DOUBLE), dimension(3) :: first_wall_norm, dx
    real(kind=DOUBLE) :: cross_z
    real(kind=DOUBLE) :: pp_h, sh_num, sh_den, invlam_h, p_bar_h, th, d_c, pbf
    real(kind=DOUBLE) :: div_p, vol_p, mu_p, h_p_loc, a_c, prod_c
    real(kind=DOUBLE) :: fbd_dp, fbd_num, fbd_den, zfl, zfr, zfd
    real(kind=DOUBLE) :: fbd_du, lam_a, fbw_rap, fbw_f5, fbw_vr
    integer(kind=ENTIER) :: kb
    real(kind=DOUBLE) :: epsc, epsd, d2sf
    real(kind=DOUBLE), dimension(3) :: Npc
    integer(kind=ENTIER) :: jj, kk, id_se, id_el, id_sf2, id_f2, n_nb
    real(kind=DOUBLE), dimension(nsfn) :: tm
    real(kind=DOUBLE), dimension(nsfn) :: kl, kr
    real(kind=DOUBLE) :: ma_l, ma_r
    real(kind=DOUBLE), dimension(3) :: res_h
    real(kind=DOUBLE), dimension(3) :: bw
    real(kind=DOUBLE), dimension(3, 3) :: Rw

    lambda   = 0.0_DOUBLE
    sol_lr   = 0.0_DOUBLE
    rhs      = 0.0_DOUBLE
    Bp = 0.0_DOUBLE
    Mp = 0.0_DOUBLE
    Rp = 0.0_DOUBLE
    is_corner = .false.
    have_first_wall_norm = .false.
    first_wall_norm = 0.0_DOUBLE
    th = 0.0_DOUBLE
    if (pi_hat) th = theta_pi
    pp_h = 0.0_DOUBLE
    sh_num = 0.0_DOUBLE
    sh_den = 0.0_DOUBLE
    res_h = 0.0_DOUBLE
    tm = 1.0_DOUBLE
    ! kappa = 1 par defaut : a ce reglage toutes les expressions ci-dessous
    ! redonnent celles de vp, le produit par 1.0 etant exact en IEEE.
    kl = 1.0_DOUBLE
    kr = 1.0_DOUBLE

    do j = 1, nsfn
      id_sub_face = mesh%vert(i_vert)%sub_face_neigh(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      idl = mesh%face(id_face)%left_neigh
      idr = mesh%face(id_face)%right_neigh
      gl = gamma_arr(idl)
      sol_lr(:, 1, j) = sol(:, idl)
      pl = pressure(sol_lr(:, 1, j), gl)
      vl = dot_product(sol_lr(2:4, 1, j), mesh%face(id_face)%norm)

      ! kappa par COIN : on le prend sur la SOUS-MAILLE de la sous-face, pas par
      ! un tableau indice au noeud. Sur une sous-face periodique le sous-element
      ! gauche porte le sommet JUMEAU, qui est bien le coin geometriquement
      ! adjacent a cette sous-face -- l'indexer par i_vert donnerait faux.
      if (pi_kappa) then
        lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
        if (lse > 0) kl(j) = kap_se(lse)
        rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
        if (rse > 0) kr(j) = kap_se(rse)
      end if

      if (second_order) then
        taul = sol_lr(1, 1, j)
        al   = sqrt(gl * pl * taul)
        dx = mesh%vert(i_vert)%coord - mesh%elem(idl)%coord
        pl = pl + dot_product(grad_p(:, idl), dx) - 0.5_DOUBLE*rec_dt_fac*dt*al**2/taul*div_v(idl)
        vl = vl + dot_product(matmul(grad_v(:,:,idl), dx), mesh%face(id_face)%norm) &
                - 0.5_DOUBLE*rec_dt_fac*dt*taul*dot_product(grad_p(:, idl), mesh%face(id_face)%norm)
      end if

      if (idr > 0) then
        gr = gamma_arr(idr)
        sol_lr(:, 2, j) = sol(:, idr)
        pr = pressure(sol_lr(:, 2, j), gr)
        vr = dot_product(sol_lr(2:4, 2, j), mesh%face(id_face)%norm)

        if (second_order) then
          taur = sol_lr(1, 2, j)
          ar   = sqrt(gr * pr * taur)
          dx = mesh%vert(i_vert)%coord - mesh%elem(idr)%coord
          pr = pr + dot_product(grad_p(:, idr), dx) - 0.5_DOUBLE*rec_dt_fac*dt*ar**2/taur*div_v(idr)
          vr = vr + dot_product(matmul(grad_v(:,:,idr), dx), mesh%face(id_face)%norm) &
                  - 0.5_DOUBLE*rec_dt_fac*dt*taur*dot_product(grad_p(:, idr), mesh%face(id_face)%norm)
        end if

        ! lambda = l'impedance de la matrice de coin M_pc. EUCCLHYD (Maire et al.,
        ! repris dans Vilar 2012 sec. 3.3.6) prend l'impedance acoustique PURE,
        ! M_pc = z_c [l-(n-@n-) + l+(n+@n+)] avec z_c = rho_c a_c, obtenue en ecrivant
        ! l'invariant de Riemann Pi_pc - P_c = -z_c (U_p - U_c).n_pc le long de chaque
        ! demi-arete. Les deux termes supplementaires ci-dessous (saut de pression et
        ! saut de vitesse) sont un ajout de type deux-chocs propre a ce code : ils
        ! n'apparaissent pas dans EUCCLHYD et ajoutent de la dissipation la ou les sauts
        ! sont forts. A noter que les branches de BORD plus bas utilisent deja z_c seul,
        ! donc l'interieur et le bord ne suivent pas la meme definition aujourd'hui.
        if (lambda_acoustic_only) then
          lambda(1, j) = sqrt(gl*pl*sol_lr(1,1,j))/sol_lr(1,1,j)
          lambda(2, j) = sqrt(gr*pr*sol_lr(1,2,j))/sol_lr(1,2,j)
        else
          lambda(1, j) = max(sqrt(gl*pl*sol_lr(1,1,j))/sol_lr(1,1,j), &
            sqrt(max(0.0_DOUBLE, pr-pl)/sol_lr(1,1,j)), -(vr-vl)/sol_lr(1,1,j))
          lambda(2, j) = max(sqrt(gr*pr*sol_lr(1,2,j))/sol_lr(1,2,j), &
            sqrt(max(0.0_DOUBLE, pl-pr)/sol_lr(1,2,j)), -(vr-vl)/sol_lr(1,2,j))
          if (pi_conv_stab) then
            ! Taux CONVECTIF EXPLICITE : le saut de MAILLE A MAILLE |v_r - v_l|,
            ! et non l ecart au noeud |v_p - v_c|. Les deux sont O(h grad v) donc
            ! donnent le meme O(Ma^2), mais celui-ci est connu AVANT de resoudre.
            ! La version au noeud fait diverger le point fixe : mu depend de v_p
            ! et v_p depend de 1/mu par l inverse de M_p, d ou petit mu -> v_p
            ! enorme -> mu enorme. Mesure : 46 % des sous-faces au plancher et
            ! max|v_p| a 1e28 sur Gresho.
            ! Ce terme n est pas neuf : c est le troisieme du max deja present
            ! dans lambda = max(rho a, sqrt(rho dp), -rho dv), pris SEUL au lieu
            ! d etre plafonne par rho a.
            vn_cv = vr - vl
            al_cv = sqrt(max(0.0_DOUBLE, gl*pl*sol_lr(1,1,j)))/sol_lr(1,1,j)
            ar_cv = sqrt(max(0.0_DOUBLE, gr*pr*sol_lr(1,2,j)))/sol_lr(1,2,j)
            mu_cv = abs(vn_cv)/sol_lr(1,1,j)
            conv_n = conv_n + 1
            if (mu_cv < conv_floor*al_cv) conv_nfloor = conv_nfloor + 1
            lambda(1, j) = max(mu_cv, conv_floor*al_cv)
            mu_cv = abs(vn_cv)/sol_lr(1,2,j)
            conv_n = conv_n + 1
            if (mu_cv < conv_floor*ar_cv) conv_nfloor = conv_nfloor + 1
            lambda(2, j) = max(mu_cv, conv_floor*ar_cv)
          end if
        end if
        ! Distribution de lambda_schema / (rho a) PAR SOUS-FACE. `compute_dt` met
        ! `lambda = rho*c` sans condition, alors que la fermeture utilise le max a
        ! trois termes, qui est >= rho a par construction. Les deux lambdas sont
        ! donc differents et c'est le PLUS PETIT qui est dans la borne de temps.
        ! On mesure l'ecart au lieu de le supposer petit ou grand : les deux queues
        ! et des comptes, parce qu'un ecart sur une poignee de sous-faces et un
        ! ecart systematique ne demandent pas la meme conclusion.
        ! Accumule somme_f (A_pcf * lambda_SCHEMA_pcf) par maille, avec le lambda
        ! du cote de la maille. Sert a evaluer, APRES que dt est choisi, le poids
        ! alpha_c = 1 - dt somme(A lambda)/m_c de eq:cflexp de la these : des qu'il
        ! est negatif, U_c^{n+1} cesse d'etre une combinaison convexe d'etats
        ! admissibles. `compute_dt` utilise rho a et non ce lambda, donc alpha_c
        ! peut etre mal evalue -- c'est le compte par MAILLE qui le dit, pas le
        ! taux par sous-face, puisque alpha_c est une somme sur la maille et qu'une
        ! seule sous-face a fort saut suffit a la dominer.
        if (diag_alpha .and. allocated(alpha_sal)) then
          alpha_sal(idl) = alpha_sal(idl) &
            + mesh%sub_face(id_sub_face)%area*lambda(1, j)
          if (idr > 0) alpha_sal(idr) = alpha_sal(idr) &
            + mesh%sub_face(id_sub_face)%area*lambda(2, j)
        end if
        if (diag_lamratio) then
          ! LES DEUX COTES. La premiere version ne sondait que `lambda(1, j)`, donc
          ! un seul cote d'une grandeur qui en a deux -- et le terme de saut de
          ! pression `sqrt(max(0, pr-pl)/tau)` ne s'allume que du cote dont le
          ! VOISIN est le plus chaud. Sur un maillage 1D oriente, `pl > pr` sur
          ! toutes les faces, donc le terme ne s'allume JAMAIS dans `lambda(1, j)`
          ! et l'instrument rendait exactement 1 sur 100 % des visites de sod_1d --
          ! un artefact d'orientation, pas une mesure. Sur sedov_quad il voyait
          ! environ la moitie des evenements, donc ses chiffres etaient des bornes
          ! INFERIEURES. Qu'un instrument rende exactement 1 sur dix millions
          ! d'echantillons sans une seule exception etait le signal : une grandeur
          ! physique ne fait pas ca.
          block
            integer(kind=ENTIER) :: kl2, kside
            real(kind=DOUBLE) :: ra, rr2
            do kside = 1, 2
              if (kside == 1) then
                ra = sqrt(max(0.0_DOUBLE, gl*pl*sol_lr(1,1,j)))/sol_lr(1,1,j)
                rr2 = lambda(1, j)
              else
                ra = sqrt(max(0.0_DOUBLE, gr*pr*sol_lr(1,2,j)))/sol_lr(1,2,j)
                rr2 = lambda(2, j)
              end if
              if (ra > 0.0_DOUBLE .and. rr2 == rr2) then
                rr2 = rr2/ra
                lamr_n = lamr_n + 1
                if (rr2 <= 1.0_DOUBLE + 1.0e-14_DOUBLE) lamr_eq1 = lamr_eq1 + 1
                if (rr2 > 2.0_DOUBLE) lamr_g2 = lamr_g2 + 1
                if (rr2 > 10.0_DOUBLE) lamr_g10 = lamr_g10 + 1
                if (rr2 > 1.0e3_DOUBLE) lamr_g1e3 = lamr_g1e3 + 1
                lamr_max = max(lamr_max, rr2)
                kl2 = min(16, max(0, int(floor(log10(max(rr2, 1.0_DOUBLE))))))
                lamr_hist(kl2) = lamr_hist(kl2) + 1
              end if
            end do
          end block
        end if

        ! p_p du schema vp_pp, ponderation en 1/lambda : le systeme nodal en
        ! quantite de mouvement ne determine PAS p_p (il sort du bilan), il faut
        ! donc le prendre d'une construction independante.
        if (pi_hat) then
          invlam_h = 1.0_DOUBLE/lambda(1,j) + 1.0_DOUBLE/lambda(2,j)
          p_bar_h = (pl/lambda(1,j) + pr/lambda(2,j) - (vr - vl))/invlam_h
          sh_num = sh_num + mesh%sub_face(id_sub_face)%area*invlam_h*p_bar_h
          sh_den = sh_den + mesh%sub_face(id_sub_face)%area*invlam_h
        end if
        if (pi_mach_damp) then
          ! Mach local par maille, sur la NORME de la vitesse, et max des deux
          ! cotes de la sous-face : theta_M -> 1 au choc, -> Ma en bas Mach.
          ma_l = sqrt(dot_product(sol(2:4, idl), sol(2:4, idl))) &
               / max(1.0e-300_DOUBLE, sqrt(max(0.0_DOUBLE, gl*pl*sol_lr(1,1,j))))
          ma_r = sqrt(dot_product(sol(2:4, idr), sol(2:4, idr))) &
               / max(1.0e-300_DOUBLE, sqrt(max(0.0_DOUBLE, gr*pr*sol_lr(1,2,j))))
          tm(j) = min(max(ma_l, ma_r, theta_m_floor), 1.0_DOUBLE)
          pih_tm_min = min(pih_tm_min, tm(j))
          pih_tm_max = max(pih_tm_max, tm(j))
        end if
        if (.not. vp_is_imposed_i) then
          ! Deux termes, deux facteurs : (1-theta_g) sur le saut de pression,
          ! theta_M sur la part visqueuse. A theta_g = 0 et theta_M = 1 on
          ! retrouve exactement v_bar de vp.
          Mp = Mp + mesh%sub_face(id_sub_face)%area &
            * (kl(j)*tm(j)*lambda(1,j) + kr(j)*tm(j)*lambda(2,j)) &
            * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
          Rp = Rp + mesh%sub_face(id_sub_face)%area &
            * ((kl(j)*tm(j)*lambda(1,j)*vl + kr(j)*tm(j)*lambda(2,j)*vr) &
               - (1.0_DOUBLE - th)*(pr - pl)) * mesh%face(id_face)%norm
        end if
      else
        if (b2d .and. abs(mesh%face(id_face)%norm(3)) > 1e-8_DOUBLE) then
          lambda(1, j) = sqrt(gl*pl*sol_lr(1,1,j))/sol_lr(1,1,j)
          if (.not. vp_is_imposed_i) then
            Mp = Mp + mesh%sub_face(id_sub_face)%area*2*kl(j)*lambda(1,j) &
              * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
          end if
        else if (idr < 0 .and. trim(bc_type(-idr)) == 'pressure') then
          ! Etat fantome rho_r = rho_l, v_r.n = v_l.n, p_r = p_bc ; meme fermeture qu'une face interieure.
          pr = bc_val(1, -idr)
          vr = vl
          taur = sol_lr(1, 1, j)
          if (lambda_acoustic_only) then
            lambda(1, j) = sqrt(gl*pl*taur)/taur
            lambda(2, j) = sqrt(max(0.0_DOUBLE, gl*pr*taur))/taur
          else
            lambda(1, j) = max(sqrt(gl*pl*taur)/taur, sqrt(max(0.0_DOUBLE, pr-pl)/taur), -(vr-vl)/taur)
            lambda(2, j) = max(sqrt(max(0.0_DOUBLE, gl*pr*taur))/taur, &
              sqrt(max(0.0_DOUBLE, pl-pr)/taur), -(vr-vl)/taur)
          end if
          if (.not. vp_is_imposed_i) then
            Mp = Mp + mesh%sub_face(id_sub_face)%area*kl(j)*(lambda(1,j) + lambda(2,j)) &
              * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
            Rp = Rp + mesh%sub_face(id_sub_face)%area &
              * (kl(j)*(lambda(1,j)*vl + lambda(2,j)*vr) - (pr - pl)) * mesh%face(id_face)%norm
          end if
        else
          lambda(1, j) = sqrt(gl*pl*sol_lr(1,1,j))/sol_lr(1,1,j)
          ! Mur symetrique : lambda de l'etat miroir (v_r.n = -v_l.n, p_r = p_l), terme -(v_r - v_l)/tau.
          if (trim(wall_mode) == 'symmetric' .and. .not. lambda_acoustic_only) &
            lambda(1, j) = max(lambda(1, j), 2.0_DOUBLE*vl/sol_lr(1,1,j))
          ! etat miroir : p_r = p_l, lam_r = lam_l, v_r = -v_l, d'ou
          ! p_bar = p_l + lam_l v_l, la meme valeur que celle du second membre.
          if (pi_mach_damp) then
            ma_l = sqrt(dot_product(sol(2:4, idl), sol(2:4, idl))) &
                 / max(1.0e-300_DOUBLE, sqrt(max(0.0_DOUBLE, gl*pl*sol_lr(1,1,j))))
            tm(j) = min(max(ma_l, theta_m_floor), 1.0_DOUBLE)
            pih_tm_min = min(pih_tm_min, tm(j))
            pih_tm_max = max(pih_tm_max, tm(j))
          end if
          if (pi_hat) then
            invlam_h = 2.0_DOUBLE/lambda(1,j)
            p_bar_h = pl + lambda(1,j)*vl
            sh_num = sh_num + mesh%sub_face(id_sub_face)%area*invlam_h*p_bar_h
            sh_den = sh_den + mesh%sub_face(id_sub_face)%area*invlam_h
          end if
          if (.not. vp_is_imposed_i) then
            if (.not. have_first_wall_norm) then
              first_wall_norm = mesh%face(id_face)%norm
              have_first_wall_norm = .true.
            else
              cross_z = first_wall_norm(1)*mesh%face(id_face)%norm(2) &
                      - first_wall_norm(2)*mesh%face(id_face)%norm(1)
              if (abs(cross_z) > 1.0e-8_DOUBLE) is_corner = .true.
            end if
            Bp = Bp + mesh%sub_face(id_sub_face)%area * mesh%face(id_face)%norm
            Mp = Mp + mesh%sub_face(id_sub_face)%area*kl(j)*tm(j)*lambda(1,j) &
              * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
            ! la part theta*p_p des sous-faces de paroi ne s'annule PAS (somme(A n)
            ! = Bp y est non nul) : elle est ajoutee en bloc apres le calcul de p_p.
            if (pi_face) then
              ! etat miroir : p_bar_f = p_l, donc beta (p_bar - p_c) s'annule et
              ! la paroi retrouve exactement le traitement de vp.
              Rp = Rp + mesh%sub_face(id_sub_face)%area &
                * (pl + kl(j)*tm(j)*lambda(1,j)*vl) * mesh%face(id_face)%norm
            else
              Rp = Rp + mesh%sub_face(id_sub_face)%area &
                * ((1.0_DOUBLE - th)*pl + kl(j)*tm(j)*lambda(1,j)*vl) &
                * mesh%face(id_face)%norm
            end if
          end if
        end if
      end if
    end do

    if (pi_hat) then
      if (sh_den > 0.0_DOUBLE) pp_h = sh_num/sh_den
      ! aux sous-faces interieures les deux cotes donnent +A theta p_p n et
      ! -A theta p_p n : la somme est nulle quel que soit le noeud. Il ne reste
      ! donc que les sous-faces de PAROI, dont somme(A n) vaut exactement Bp.
      if (.not. pi_face) Rp = Rp + th*pp_h*Bp
    end if

    if (pi_hat .and. .not. vp_is_imposed_i) then
      ! Sans cette garde l'inversion rend un NaN silencieux : inverse_3_by_3 fait
      ! mat_inv = 1/det * adj sans tester det. Mieux vaut un run qui refuse de
      ! partir qu'un champ entierement NaN sorti avec EXIT_CODE 0.
      if (maxval(abs(Mp)) <= 0.0_DOUBLE) then
        print*, "PIH : matrice nodale IDENTIQUEMENT NULLE au noeud", i_vert, &
          " -- theta_M s'annule sur toutes ses sous-faces. Relever theta_m_floor. Arret."
        error stop
      end if
    end if

    if (pi_kappa .and. .not. vp_is_imposed_i) then
      ! Meme garde que ci-dessus, meme raison : M_p = somme A kappa lambda n(x)n
      ! est identiquement nulle si kappa s'annule sur toutes les sous-faces d'un
      ! noeud, et inverse_3_by_3 ne teste pas son determinant.
      if (maxval(abs(Mp)) <= 0.0_DOUBLE) then
        print*, "KAP : matrice nodale IDENTIQUEMENT NULLE au noeud", i_vert, &
          " -- kappa s'annule sur toutes ses sous-faces. Relever kap_floor. Arret."
        error stop
      end if
    end if

    if (.not. vp_is_imposed_i) then
      if (is_corner) then
        vp_i = 0.0_DOUBLE
      else if (maxval(abs(Bp)) > 1e-8_DOUBLE .and. trim(wall_mode) == 'symmetric') then
        ! Condition de symetrie : chaque sous-face a sa jumelle miroir, la solution verifie v_p.b = 0.
        bw = Bp/norm2(Bp)
        Rw = -2.0_DOUBLE*tensor_product(bw, bw)
        Rw(1,1) = Rw(1,1) + 1.0_DOUBLE; Rw(2,2) = Rw(2,2) + 1.0_DOUBLE; Rw(3,3) = Rw(3,3) + 1.0_DOUBLE
        Mp = Mp + matmul(Rw, matmul(Mp, Rw))
        Rp = Rp + matmul(Rw, Rp)
        call inverse_3_by_3(Mp, Mp_inv)
        vp_i = matmul(Mp_inv, Rp)
      else if (maxval(abs(Bp)) > 1e-8_DOUBLE) then
        call inverse_3_by_3(Mp, Mp_inv)
        Pp = dot_product(Rp, matmul(Mp_inv, Bp)) / dot_product(Bp, matmul(Mp_inv, Bp))
        vp_i = matmul(Mp_inv, Rp - Pp*Bp)
      else
        call inverse_3_by_3(Mp, Mp_inv)
        vp_i = matmul(Mp_inv, Rp)
      end if
    end if

    ! --- 'vp_mp' : divergence NODALE et viscosite de volume au noeud
    div_p = 0.0_DOUBLE; mu_p = 0.0_DOUBLE; vol_p = 0.0_DOUBLE
    if (pi_mp) then
      mp_n_node_tot = mp_n_node_tot + 1
      ! acoustic_module ne definit la divergence nodale QUE si .not. is_bound, et
      ! pour une raison de fond : somme_c N_pc = 0 n'est vraie qu'a l'interieur.
      ! L'etendre au bord casserait la conservation de quantite de mouvement, que
      ! cette identite est seule a garantir. On laisse donc zero au bord, et on
      ! COMPTE les noeuds concernes au lieu de le passer sous silence.
      if (mesh%vert(i_vert)%is_bound) then
        mp_n_bnd_node = mp_n_bnd_node + 1
      else
        a_c = 0.0_DOUBLE; h_p_loc = 0.0_DOUBLE; n_nb = 0
        do jj = 1, mesh%vert(i_vert)%n_sub_elems_neigh
          id_se = mesh%vert(i_vert)%sub_elem_neigh(jj)
          id_el = mesh%sub_elem(id_se)%mesh_elem
          Npc = 0.0_DOUBLE
          do kk = 1, mesh%sub_elem(id_se)%n_sub_faces
            id_sf2 = mesh%sub_elem(id_se)%sub_face(kk)
            id_f2 = mesh%sub_face(id_sf2)%mesh_face
            if (mesh%sub_face(id_sf2)%left_elem_neigh == id_el) then
              Npc = Npc + mesh%sub_face(id_sf2)%area*mesh%face(id_f2)%norm
            else
              Npc = Npc - mesh%sub_face(id_sf2)%area*mesh%face(id_f2)%norm
            end if
          end do
          div_p = div_p - dot_product(Npc, sol(2:4, id_el))
          vol_p = vol_p + mesh%sub_elem(id_se)%volume
          a_c = a_c + sqrt(max(0.0_DOUBLE, gamma_arr(id_el) &
                 *pressure(sol(:, id_el), gamma_arr(id_el))/sol(1, id_el)))/sol(1, id_el)
          h_p_loc = h_p_loc + norm2(mesh%vert(i_vert)%coord - mesh%elem(id_el)%coord)
          n_nb = n_nb + 1
        end do
        if (vol_p > 0.0_DOUBLE .and. n_nb > 0) then
          div_p = div_p/vol_p
          ! longueur nodale : distance moyenne aux centres des mailles voisines.
          ! Definie independamment de l'extrusion, contrairement a V_p^(1/3).
          h_p_loc = h_p_loc/real(n_nb, DOUBLE)
          ! mu_p homogene a rho x vitesse x longueur, et NUL en detente
          if (div_p < 0.0_DOUBLE) mu_p = mp_cvisc*(a_c/real(n_nb, DOUBLE))*h_p_loc
          mp_prod_node = mp_prod_node + mu_p*vol_p*div_p**2
          if (diag_firstprint .and. .not. dr_done) &
            dr_max_divp = max(dr_max_divp, abs(div_p))
        else
          div_p = 0.0_DOUBLE
        end if
      end if
    end if

    ! --- D_p nodal de la fermeture 'vp_fbd', assemble AVANT la boucle de flux
    ! parce qu'il est global au noeud. Contrairement a `div_p` de 'vp_mp', il n'est
    ! PAS mis a zero au bord : l'exclusion y etait imposee par `somme_c N_pc = 0`,
    ! identite qui garantit seule la conservation de quantite de mouvement pour
    ! cette forme-la. Ici Pi est unique par sous-face, donc la conservation est
    ! automatique et ne depend pas de D_p. La contribution du bord est instrumentee
    ! au lieu d'etre supposee inoffensive.
    ! RESERVE : D_p est bati sur les moyennes de maille, donc l'extrapolation GRP
    ! demi-pas de l'ordre 2 ne lui est PAS appliquee.
    fbd_dp = 0.0_DOUBLE; fbd_du = 0.0_DOUBLE
    if (pi_fbd .or. pi_fbw) then
      block
        integer(kind=ENTIER) :: jf, jdl, jdr, jsf, jff
        real(kind=DOUBLE) :: za, zb, van, vbn, aa
        real(kind=DOUBLE), dimension(3) :: nrm
        fbd_num = 0.0_DOUBLE; fbd_den = 0.0_DOUBLE
        do jf = 1, nsfn
          jsf = mesh%vert(i_vert)%sub_face_neigh(jf)
          jff = mesh%sub_face(jsf)%mesh_face
          jdl = mesh%face(jff)%left_neigh
          jdr = mesh%face(jff)%right_neigh
          if (jdl <= 0) cycle
          aa = mesh%sub_face(jsf)%area
          nrm = mesh%face(jff)%norm
          van = dot_product(sol(2:4, jdl), nrm)
          if (jdr > 0) then
            vbn = dot_product(sol(2:4, jdr), nrm)
            ! Du = le plus grand saut de vitesse normale du noeud. Il entre dans
            ! lambda_w, qui est par maille mais alimente par des quantites NODALES.
            fbd_du = max(fbd_du, abs(vbn - van))
            if (fbd_lam_mode == 'tw' .and. .not. pi_fbw) then
              za = fbd_lambda(sol(1, jdl), pressure(sol(:, jdl), gamma_arr(jdl)), &
                     gamma_arr(jdl), pressure(sol(:, jdr), gamma_arr(jdr)), vbn - van)
              zb = fbd_lambda(sol(1, jdr), pressure(sol(:, jdr), gamma_arr(jdr)), &
                     gamma_arr(jdr), pressure(sol(:, jdl), gamma_arr(jdl)), vbn - van)
            else
              ! pi_fbw : lambda FIGE a rho a dans le denominateur de D_p. C'est ce
              ! gel qui casse la proportionnalite D_p ~ lambda et rend la condition
              ! de positivite soluble.
              ! 'acou' et 'sep' : acoustique PURE au denominateur de D_p
              za = fbd_lambda(sol(1, jdl), pressure(sol(:, jdl), gamma_arr(jdl)), &
                     gamma_arr(jdl), -1.0_DOUBLE, 0.0_DOUBLE)
              zb = fbd_lambda(sol(1, jdr), pressure(sol(:, jdr), gamma_arr(jdr)), &
                     gamma_arr(jdr), -1.0_DOUBLE, 0.0_DOUBLE)
            end if
            if (za <= 0.0_DOUBLE .or. zb <= 0.0_DOUBLE) cycle
          else
            fbd_du = max(fbd_du, 2.0_DOUBLE*abs(van))
            if (fbd_lam_mode == 'tw' .and. .not. pi_fbw) then
              za = fbd_lambda(sol(1, jdl), pressure(sol(:, jdl), gamma_arr(jdl)), &
                     gamma_arr(jdl), pressure(sol(:, jdl), gamma_arr(jdl)), &
                     -2.0_DOUBLE*van)
            else
              za = fbd_lambda(sol(1, jdl), pressure(sol(:, jdl), gamma_arr(jdl)), &
                     gamma_arr(jdl), -1.0_DOUBLE, 0.0_DOUBLE)
            end if
            if (za <= 0.0_DOUBLE) cycle
            ! etat miroir : p_r = p_l donc z_r = z_l, et v_nr = -v_nl.
            zb = za; vbn = -van
            fbd_num_bnd = fbd_num_bnd + aa*(vbn - van)
          end if
          call fbd_accum(aa, za, zb, van, vbn, fbd_num, fbd_den)
        end do
        if (fbd_den > 0.0_DOUBLE) then
          fbd_dp = fbd_num/fbd_den
          fbd_n_node = fbd_n_node + 1
          if (abs(fbd_dp) == abs(fbd_dp)) then
            fbd_dp_absmax = max(fbd_dp_absmax, abs(fbd_dp))
            fbd_st_dpmax = max(fbd_st_dpmax, abs(fbd_dp))
          end if
        else
          fbd_n_node_void = fbd_n_node_void + 1
          fbd_st_void = fbd_st_void + 1
        end if
      end block
    end if

    do j = 1, nsfn
      id_sub_face = mesh%vert(i_vert)%sub_face_neigh(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      idl = mesh%face(id_face)%left_neigh
      idr = mesh%face(id_face)%right_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      gl = gamma_arr(idl)

      v_et = dot_product(vp_i, mesh%face(id_face)%norm)
      vl = dot_product(sol_lr(2:4, 1, j), mesh%face(id_face)%norm)
      pl = pressure(sol_lr(:, 1, j), gl)
      if (second_order) then
        taul = sol_lr(1, 1, j)
        al   = sqrt(gl * pl * taul)
        dx = mesh%vert(i_vert)%coord - mesh%elem(idl)%coord
        pl = pl + dot_product(grad_p(:, idl), dx) - 0.5_DOUBLE*rec_dt_fac*dt*al**2/taul*div_v(idl)
        vl = vl + dot_product(matmul(grad_v(:,:,idl), dx), mesh%face(id_face)%norm) &
                - 0.5_DOUBLE*rec_dt_fac*dt*taul*dot_product(grad_p(:, idl), mesh%face(id_face)%norm)
      end if
      if (pi_fbw .or. pi_mdiv) then
        ! lambda_w par maille, alimente par le D_p et le Du du NOEUD.
        epsc = sol(5, idl) - 0.5_DOUBLE*dot_product(sol(2:4, idl), sol(2:4, idl))
        zfl = fbw_lambda(sol(1, idl), pl, gl, epsc, fbd_dp, fbd_du, fbw_marge)
        if (idr > 0) then
          zfr = pressure(sol_lr(:, 2, j), gamma_arr(idr))
          epsd = sol(5, idr) - 0.5_DOUBLE*dot_product(sol(2:4, idr), sol(2:4, idr))
          zfd = fbw_lambda(sol(1, idr), zfr, gamma_arr(idr), epsd, fbd_dp, fbd_du, &
                           fbw_marge)
          ! `vr` de la routine n'est affecte que PLUS BAS dans la boucle, apres le
          ! flux : le lire ici donnerait la valeur de la sous-face PRECEDENTE. On
          ! calcule donc la vitesse normale du cote droit dans sa propre variable.
          fbw_vr = dot_product(sol_lr(2:4, 2, j), mesh%face(id_face)%norm)
        else
          ! etat miroir d'une sous-face de bord
          zfr = pl; epsd = epsc; zfd = zfl; fbw_vr = -vl
        end if
        if (zfl > 0.0_DOUBLE .and. zfd > 0.0_DOUBLE) then
          pbf = fbd_pbar(pl, zfr, zfl, zfd)
        else
          pbf = pl
          fbd_n_zbad = fbd_n_zbad + 1
        end if
        pl_et = pbf - fbd_dp
        ! --- 'vp_mdiv' : correction par la dilatation du SOUS-ELEMENT, plafonnee
        ! par la condition exacte d'admissibilite de l'etat intermediaire. Ce
        ! n'est pas un limiteur : `md_kappa_cap` rend la meme inegalite que celle
        ! qui fixe lambda chez Gallice.
        if (pi_mdiv) then
          block
            integer(kind=ENTIER) :: ks, ksf, kst, rse_mdv, kbm
            real(kind=DOUBLE) :: dlt, apc, hpc, kapw, kapm, d2m, prd, rap
            if (lse > 0 .and. allocated(mdv_se)) then
              ! SYMETRISATION : moyenne arithmetique des dilatations des DEUX
              ! sous-elements de la sous-face. Sans elle, Pi etait bati sur le
              ! seul cote gauche, donc le schema n'etait pas invariant par
              ! retournement de face : echanger les etiquettes gauche/droite --
              ! une convention du fichier de maillage, sans contenu physique --
              ! changeait Pi de kappa(delta_pc - delta_pd), c'est-a-dire de
              ! l'ordre de la correction elle-meme, les deux dilatations
              ! differant en O(1) au choc. Meme defaut que la parite sous
              ! retournement de face de 'vp_pp_ppvp'.
              ! Pi reste UNIVALUE (pr_et = pl_et), donc la conservation a la
              ! face est intacte. Au BORD il n'y a pas de sous-element droit et
              ! on garde le gauche seul, ce qui est la seule valeur definie.
              rse_mdv = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
              if (rse_mdv > 0) then
                dlt = 0.5_DOUBLE*(mdv_se(lse) + mdv_se(rse_mdv))
              else
                dlt = mdv_se(lse)
              end if
              apc = 0.0_DOUBLE
              do ks = 1, mesh%sub_elem(lse)%n_sub_faces
                ksf = mesh%sub_elem(lse)%sub_face(ks)
                apc = apc + mesh%sub_face(ksf)%area
              end do
              if (apc > 0.0_DOUBLE .and. abs(mesh%sub_elem(lse)%volume) > 0.0_DOUBLE) then
                hpc = abs(mesh%sub_elem(lse)%volume)/apc
                kapw = mdiv_theta*zfl*hpc
                ! plafond de positivite. On teste le STATUT et non la valeur :
                ! elle vaut huge() pour MD_NO_CONSTR, ce qui est correct dans un
                ! min, mais ZERO pour MD_EMPTY et MD_BASE_OUT, ou le zero est un
                ! signal que la base est deja inadmissible et pas une valeur.
                call md_kappa_cap(zfl, epsc, pl, zfd, epsd, zfr, fbd_dp, dlt, &
                                  kapm, d2m, kst)
                mdv_n = mdv_n + 1
                if (dlt /= 0.0_DOUBLE) mdv_n_dltnz = mdv_n_dltnz + 1
                if (kst == MD_EMPTY) then
                  mdv_n_empty = mdv_n_empty + 1
                else if (kst == MD_BASE_OUT) then
                  mdv_n_baseout = mdv_n_baseout + 1
                else if (kst == MD_NO_CONSTR) then
                  mdv_n_nocon = mdv_n_nocon + 1
                  pl_et = pl_et - kapw*dlt
                else
                  ! Histogramme sur TOUTES les visites contraintes, et pas
                  ! seulement sur celles qui sont plafonnees : c'est la seule
                  ! facon de distinguer « le plafond mord partout » de « le
                  ! plafond mord sur trois sous-faces et le max vient de la ».
                  rap = kapw/max(kapm, 1.0e-300_DOUBLE)
                  kbm = min(16, max(-16, int(floor(log10(max(rap, 1.0e-16_DOUBLE))))))
                  mdv_rhist(kbm) = mdv_rhist(kbm) + 1
                  if (kapw > kapm) then
                    mdv_n_cap = mdv_n_cap + 1
                    mdv_ratio_max = max(mdv_ratio_max, rap)
                    kapw = kapm
                  end if
                  pl_et = pl_et - kapw*dlt
                end if
                mdv_kap_max = max(mdv_kap_max, kapw)
                ! Signe de la production. Le retard d'un pas sur delta_pc peut en
                ! changer le signe et produire de l'ANTI-dissipation : c'est le
                ! seul risque propre a cette construction, et il se lit ici.
                prd = mesh%sub_face(id_sub_face)%area*(pl - pl_et)*(v_et - vl)
                kap_n_prod = kap_n_prod + 1
                kap_prod_tot = kap_prod_tot + prd
                if (prd < 0.0_DOUBLE) then
                  kap_n_prod_neg = kap_n_prod_neg + 1
                  kap_prod_min = min(kap_prod_min, prd)
                end if
              end if
            end if
          end block
        end if
        fbw_n_sf = fbw_n_sf + 1
        ! distribution de lambda_w / rho a, en bins log10
        if (zfl > 0.0_DOUBLE) then
          lam_a = sqrt(max(0.0_DOUBLE, gl*pl*sol(1, idl)))/sol(1, idl)
          if (lam_a > 0.0_DOUBLE) then
            fbw_rap = zfl/lam_a
            fbw_n = fbw_n + 1
            if (fbw_rap <= 1.0_DOUBLE + 1.0e-14_DOUBLE) fbw_n_eq1 = fbw_n_eq1 + 1
            if (fbw_rap > 10.0_DOUBLE) fbw_n_g10 = fbw_n_g10 + 1
            if (fbw_rap > 1.0e3_DOUBLE) fbw_n_g1e3 = fbw_n_g1e3 + 1
            if (fbw_rap > 1.0e6_DOUBLE) fbw_n_g1e6 = fbw_n_g1e6 + 1
            kb = min(16, max(-2, int(floor(log10(max(fbw_rap, 1.0e-2_DOUBLE))))))
            fbw_hist(kb) = fbw_hist(kb) + 1
            fbw_st_rmax = max(fbw_st_rmax, fbw_rap)
          end if
        end if
        ! delta^2 ENTIER, avec lambda_w et le Pi retenu
        if (zfl > 0.0_DOUBLE .and. zfd > 0.0_DOUBLE) then
          d2sf = fbd_delta2(zfl, epsc, pl, zfd, epsd, zfr, pl_et)
          if (d2sf == d2sf) then
            fbw_d2_min = min(fbw_d2_min, d2sf)
            fbw_st_d2min = min(fbw_st_d2min, d2sf)
          end if
          ! visitees par type : les denominateurs des taux
          if (max(pl, zfr) < 1.0e-6_DOUBLE) then
            fbw_nv_ff = fbw_nv_ff + 1
          else if (min(pl, zfr) < 1.0e-6_DOUBLE) then
            fbw_nv_cf = fbw_nv_cf + 1
          else
            fbw_nv_cc = fbw_nv_cc + 1
          end if
          if (d2sf <= 0.0_DOUBLE) then
            if (max(pl, zfr) < 1.0e-6_DOUBLE) then
              fbw_n_viol_ff = fbw_n_viol_ff + 1; fbw_st_ff = fbw_st_ff + 1
            else if (min(pl, zfr) < 1.0e-6_DOUBLE) then
              fbw_n_viol_cf = fbw_n_viol_cf + 1; fbw_st_cf = fbw_st_cf + 1
            else
              fbw_n_viol_cc = fbw_n_viol_cc + 1; fbw_st_cc = fbw_st_cc + 1
            end if
            if (allocated(fbw_cviol)) then
              fbw_cviol(idl) = 1
              if (idr > 0) fbw_cviol(idr) = 1
            end if
          end if
        end if
      else if (pi_fbd) then
        ! Pi = pbar_f - D_p, pbar_f PONDERE PAR L'IMPEDANCE (et non la moyenne
        ! arithmetique de 'vp_mp'/'vp_facebar'), avec l'impedance a trois termes.
        if (idr > 0) then
          zfr = pressure(sol_lr(:, 2, j), gamma_arr(idr))
          if (fbd_lam_mode == 'acou') then
            zfl = fbd_lambda(sol_lr(1, 1, j), pl, gl, -1.0_DOUBLE, 0.0_DOUBLE)
            zfd = fbd_lambda(sol(1, idr), zfr, gamma_arr(idr), -1.0_DOUBLE, 0.0_DOUBLE)
          else
            zfl = fbd_lambda(sol_lr(1, 1, j), pl, gl, zfr, vr - vl)
            zfd = fbd_lambda(sol(1, idr), zfr, gamma_arr(idr), pl, vr - vl)
          end if
          if (zfl > 0.0_DOUBLE .and. zfd > 0.0_DOUBLE) then
            pbf = fbd_pbar(pl, zfr, zfl, zfd)
          else
            pbf = pl
            fbd_n_zbad = fbd_n_zbad + 1
            fbd_st_zbad = fbd_st_zbad + 1
            if (zfl /= zfl .or. zfd /= zfd) then
              fbd_n_znan = fbd_n_znan + 1
            else
              fbd_n_zzero = fbd_n_zzero + 1
            end if
            ! d'ou vient-il : un etat deja casse en amont, ou le max lui-meme ?
            if (pl /= pl .or. zfr /= zfr) then
              fbd_zb_pnan = fbd_zb_pnan + 1
            else if (pl <= 0.0_DOUBLE .or. zfr <= 0.0_DOUBLE) then
              fbd_zb_pneg = fbd_zb_pneg + 1
            else if (sol_lr(1, 1, j) /= sol_lr(1, 1, j) .or. sol(1, idr) /= sol(1, idr)) then
              fbd_zb_taunan = fbd_zb_taunan + 1
            else if (sol_lr(1, 1, j) <= 0.0_DOUBLE .or. sol(1, idr) <= 0.0_DOUBLE) then
              fbd_zb_tauneg = fbd_zb_tauneg + 1
            else
              fbd_zb_autre = fbd_zb_autre + 1
            end if
          end if
        else
          ! etat miroir : p_r = p_l et lambda_r = lambda_l, donc pbar_f = p_l.
          pbf = pl
          zfl = -1.0_DOUBLE; zfd = -1.0_DOUBLE
        end if
        pl_et = pbf - fbd_dp
        fbd_n_sf = fbd_n_sf + 1
        ! Positivite : le BILAN ENTIER de la sous-face, pas un critere par cote.
        ! Classe par type d'interface, parce qu'un compte global melange deux
        ! populations dont une seule est en question.
        if (idr > 0 .and. zfl > 0.0_DOUBLE .and. zfd > 0.0_DOUBLE) then
          epsc = sol(5, idl) - 0.5_DOUBLE*dot_product(sol(2:4, idl), sol(2:4, idl))
          epsd = sol(5, idr) - 0.5_DOUBLE*dot_product(sol(2:4, idr), sol(2:4, idr))
          d2sf = fbd_delta2(zfl, epsc, pl, zfd, epsd, zfr, pl_et)
          fbd_d2_min = min(fbd_d2_min, d2sf)
          if (d2sf == d2sf) fbd_st_d2min = min(fbd_st_d2min, d2sf)
          if (d2sf <= 0.0_DOUBLE) then
            if (max(pl, zfr) < 1.0e-6_DOUBLE) then
              fbd_n_viol_ff = fbd_n_viol_ff + 1; fbd_st_ff = fbd_st_ff + 1
            else if (min(pl, zfr) < 1.0e-6_DOUBLE) then
              fbd_n_viol_cf = fbd_n_viol_cf + 1; fbd_st_cf = fbd_st_cf + 1
            else
              fbd_n_viol_cc = fbd_n_viol_cc + 1; fbd_st_cc = fbd_st_cc + 1
            end if
          end if
        end if
      else if (pi_mp) then
        ! Pi = p_bar_f - mu_p (div v)_p, UNE SEULE valeur par sous-face, donc la
        ! meme des deux cotes : c'est elle qui rend la conservation gratuite.
        if (idr > 0) then
          pbf = 0.5_DOUBLE*(pressure(sol_lr(:, 1, j), gl) &
                          + pressure(sol_lr(:, 2, j), gamma_arr(idr)))
        else
          pbf = pressure(sol_lr(:, 1, j), gl)
        end if
        pl_et = pbf - mu_p*div_p
      else if (pi_face) then
        ! p_bar_f = (p_l + p_r)/2 ; au bord l'etat miroir donne p_r = p_l.
        if (idr > 0) then
          pbf = 0.5_DOUBLE*(pressure(sol_lr(:, 1, j), gl) + pressure(sol_lr(:, 2, j), gamma_arr(idr)))
        else
          pbf = pressure(sol_lr(:, 1, j), gl)
        end if
        pl_et = pl + th*(pbf - pl) - tm(j)*lambda(1, j)*(v_et - vl)
      else
        pl_et = th*pp_h + (1.0_DOUBLE - th)*pl &
              - kl(j)*tm(j)*lambda(1, j)*(v_et - vl)
      end if
      if (pi_hat) then
        d_c = v_et - vl
        pih_sq_step   = pih_sq_step &
          + tm(j)*mesh%sub_face(id_sub_face)%area*lambda(1, j)*d_c**2
        pih_free_step = pih_free_step + th*mesh%sub_face(id_sub_face)%area*(pl - pp_h)*d_c
        res_h = res_h + mesh%sub_face(id_sub_face)%area*pl_et*mesh%face(id_face)%norm
      end if

      if (md_watch) then
        w_dpi_min(idl) = min(w_dpi_min(idl), pl_et - pl)
        w_dpi_max(idl) = max(w_dpi_max(idl), pl_et - pl)
        w_prod(idl) = w_prod(idl) &
          + mesh%sub_face(id_sub_face)%area*(pl - pl_et)*(v_et - vl)
      end if
      if (pi_kappa) then
        prod_c = mesh%sub_face(id_sub_face)%area*(pl - pl_et)*(v_et - vl)
        kap_n_prod = kap_n_prod + 1
        kap_prod_tot = kap_prod_tot + prod_c
        if (prod_c < 0.0_DOUBLE) then
          kap_n_prod_neg = kap_n_prod_neg + 1
          kap_prod_min = min(kap_prod_min, prod_c)
          if (kap_neg_print < 14) then
            kap_neg_print = kap_neg_print + 1
            print*, "KAPNEG maille", idl, " ghost ?", mesh%elem(idl)%is_ghost
            print*, "KAPNEG   kappa", kl(j), " lambda", lambda(1,j), " delta", v_et - vl, &
                    " A", mesh%sub_face(id_sub_face)%area
            print*, "KAPNEG   p_c", pl, " Pi", pl_et, " p_c - Pi", pl - pl_et
            print*, "KAPNEG   prod par DIFFERENCE", prod_c, " prod par CARRE", &
                    mesh%sub_face(id_sub_face)%area*kl(j)*tm(j)*lambda(1,j)*(v_et - vl)**2
          end if
        end if
        ! Le bilan nodal exclut les faces z de l'extrusion sous boundary_2d,
        ! pour la raison ecrite au long dans compute_rhs_lagrange_meandiv : leur
        ! contribution y est dynamiquement inerte, pas compensee.
        if (.not. (b2d .and. abs(mesh%face(id_face)%norm(3)) > 1.0e-8_DOUBLE)) &
          res_h = res_h + mesh%sub_face(id_sub_face)%area*pl_et &
                * mesh%face(id_face)%norm
      end if
      if (diag_divratio) then
        dr_vdot(idl) = dr_vdot(idl) + mesh%sub_face(id_sub_face)%area*v_et
        dr_sq_delta = dr_sq_delta &
          + lambda(1, j)*mesh%sub_face(id_sub_face)%area*(v_et - vl)**2
        if (idr > 0) then
          dr_vdot(idr) = dr_vdot(idr) - mesh%sub_face(id_sub_face)%area*v_et
          dr_sq_delta = dr_sq_delta &
            + lambda(2, j)*mesh%sub_face(id_sub_face)%area &
              *(dot_product(sol_lr(2:4, 2, j), mesh%face(id_face)%norm) - v_et)**2
        end if
      end if
      if (pi_mp) then
        prod_c = mesh%sub_face(id_sub_face)%area*(pl - pl_et)*(v_et - vl)
        mp_n_cell_tot = mp_n_cell_tot + 1
        if (prod_c < 0.0_DOUBLE) then
          mp_n_cell_neg = mp_n_cell_neg + 1
          mp_prod_min_cell = min(mp_prod_min_cell, prod_c)
        end if
      end if
      flux(1)   = -v_et
      flux(2:4) = pl_et * mesh%face(id_face)%norm
      flux(5)   = pl_et * v_et
      if (pi_fbw) then
        ! Flux d'energie = MILIEU PUR de l'intervalle de positivite, bornes
        ! evaluees avec lambda_w et le Pi retenu, n sortante de la maille GAUCHE :
        !   Phi_haut = +lam_c eps_c + Pi u_c      Phi_bas = -lam_d eps_d + Pi u_d
        ! Aucune branche, aucun test d'appartenance : le flux naturel pl_et*v_et
        ! n'est jamais lu. Une seule valeur par sous-face, donc la conservation
        ! reste automatique -- c'est elle qui rachete la positivite.
        fbw_f5 = 0.5_DOUBLE*((zfl*epsc + pl_et*vl) + (-zfd*epsd + pl_et*fbw_vr))
        flux(5) = fbw_f5
      end if
      if (mesh%sub_elem(lse)%mesh_vert == i_vert) then
        rhs(:, lse_loc) = rhs(:, lse_loc) - mesh%sub_face(id_sub_face)%area/mass(idl) * flux
        call fdiag_add(idl, mesh%sub_face(id_sub_face)%area, pl_et, &
                       mesh%face(id_face)%norm, mesh%face(id_face)%right_neigh)
      end if

      if (idr > 0) then
        gr = gamma_arr(idr)
        vr = dot_product(sol_lr(2:4, 2, j), mesh%face(id_face)%norm)
        pr = pressure(sol_lr(:, 2, j), gr)
        if (second_order) then
          taur = sol_lr(1, 2, j)
          ar   = sqrt(gr * pr * taur)
          dx = mesh%vert(i_vert)%coord - mesh%elem(idr)%coord
          pr = pr + dot_product(grad_p(:, idr), dx) - 0.5_DOUBLE*rec_dt_fac*dt*ar**2/taur*div_v(idr)
          vr = vr + dot_product(matmul(grad_v(:,:,idr), dx), mesh%face(id_face)%norm) &
                  - 0.5_DOUBLE*rec_dt_fac*dt*taur*dot_product(grad_p(:, idr), mesh%face(id_face)%norm)
        end if
        if (pi_fbd .or. pi_fbw .or. pi_mp) then
          pr_et = pl_et            ! meme Pi des deux cotes
        else if (pi_face) then
          pr_et = pr + th*(pbf - pr) + tm(j)*lambda(2, j)*(v_et - vr)
        else
          pr_et = th*pp_h + (1.0_DOUBLE - th)*pr &
                + kr(j)*tm(j)*lambda(2, j)*(v_et - vr)
        end if
        if (pi_hat) then
          ! delta_r = (u_p - v_r).n_r avec n_r = -n, donc -(v_et - vr).
          d_c = -(v_et - vr)
          pih_sq_step   = pih_sq_step &
            + tm(j)*mesh%sub_face(id_sub_face)%area*lambda(2, j)*d_c**2
          pih_free_step = pih_free_step + th*mesh%sub_face(id_sub_face)%area*(pr - pp_h)*d_c
          res_h = res_h - mesh%sub_face(id_sub_face)%area*pr_et*mesh%face(id_face)%norm
        end if
        if (md_watch) then
          ! delta_r = (u_p - v_r).n_r avec n_r = -n
          w_dpi_min(idr) = min(w_dpi_min(idr), pr_et - pr)
          w_dpi_max(idr) = max(w_dpi_max(idr), pr_et - pr)
          w_prod(idr) = w_prod(idr) &
            + mesh%sub_face(id_sub_face)%area*(pr - pr_et)*(-(v_et - vr))
        end if
        if (pi_kappa) then
          ! delta_r = (u_p - v_r).n_r avec n_r = -n
          prod_c = mesh%sub_face(id_sub_face)%area*(pr - pr_et)*(-(v_et - vr))
          kap_n_prod = kap_n_prod + 1
          kap_prod_tot = kap_prod_tot + prod_c
          if (prod_c < 0.0_DOUBLE) then
            kap_n_prod_neg = kap_n_prod_neg + 1
            kap_prod_min = min(kap_prod_min, prod_c)
          end if
          if (.not. (b2d .and. abs(mesh%face(id_face)%norm(3)) > 1.0e-8_DOUBLE)) &
            res_h = res_h - mesh%sub_face(id_sub_face)%area*pr_et &
                  * mesh%face(id_face)%norm
        end if
        flux(1)   = -v_et
        flux(2:4) = pr_et * mesh%face(id_face)%norm
        flux(5)   = pr_et * v_et
        ! meme valeur des deux cotes : le flux d'energie est UNIQUE par sous-face
        if (pi_fbw) flux(5) = fbw_f5
        rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
        if (rse > 0) then
        if (mesh%sub_elem(rse)%mesh_vert == i_vert) then
          rse_loc = mesh%sub_elem(rse)%id_loc_around_node
          rhs(:, rse_loc) = rhs(:, rse_loc) + mesh%sub_face(id_sub_face)%area/mass(idr) * flux
          call fdiag_add(idr, mesh%sub_face(id_sub_face)%area, pr_et, &
                         -mesh%face(id_face)%norm, mesh%face(id_face)%left_neigh)
        end if
        end if
      end if
    end do

    if (pi_hat) then
      ! bilan nodal : v_p . somme(A Pi_chapeau n) doit etre nul, c'est ce qui rend
      ! la construction conservative en energie. On garde le pire ecart du run.
      pih_res_max = max(pih_res_max, abs(dot_product(vp_i, res_h)))
    end if
    if (pi_kappa .and. .not. mesh%vert(i_vert)%is_bound) then
      ! Aux noeuds INTERIEURS seulement : somme_c N_pc = 0 n'est vraie que la,
      ! donc le bilan n'a pas a etre nul au bord et l'y exiger serait un faux
      ! controle.
      kap_res_max = max(kap_res_max, maxval(abs(res_h)))
    end if
  end subroutine compute_rhs_around_node

  subroutine compute_gradients(mesh, sol, gamma_arr, dt, grad_v, grad_p, div_v, p_limited, weno, &
      gg_corrected, eps_weno_in, bound_skip, bound_interior, b2d)
    use linear_solver_module, only: tensor_product_3, inverse_3_by_3
    ! Nodal WENO Green-Gauss gradients weighted by ||grad_p||.
    ! div_v = trace(grad_v).
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    real(kind=DOUBLE), intent(in) :: dt
    real(kind=DOUBLE), dimension(3, 3, mesh%n_elems), intent(out) :: grad_v
    real(kind=DOUBLE), dimension(3, mesh%n_elems), intent(out) :: grad_p
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(out) :: div_v
    integer(kind=ENTIER), dimension(mesh%n_elems), intent(out), optional :: p_limited
    ! .false. -> plain volume-weighted average of the nodal gradients (linear
    ! scheme). Default .true. = the WENO-style weight below, which is nonlinear
    ! even in smooth flow and is what the scheme has always used.
    logical, intent(in), optional :: weno
    ! .true. -> CORRECTED Green-Gauss, as aho_gg does: grad = mat^-1 . (sum A n dphi)
    ! with mat = sum A n (x) dx, instead of dividing the raw flux sum by the nodal
    ! volume Omega_p. The uncorrected form is only exact for a linear field when
    ! mat = Omega_p * I, which holds on a Cartesian hex mesh by symmetry but not in
    ! general. Default .false. = the historical behaviour.
    logical, intent(in), optional :: gg_corrected
    ! eps in the nonlinear nodal-gradient weight w_p = Omega_p/(eps + ||S_p||^4).
    ! Default tiny() makes the weighting fully nonlinear even in smooth flow; a
    ! larger eps (e.g. 1e-2) pushes it back towards the plain volume-weighted
    ! (linear) combination wherever ||S_p||^4 is small compared to eps.
    real(kind=DOUBLE), intent(in), optional :: eps_weno_in
    ! Diagnostic switch, DEFAULT .false. -> behaviour unchanged.
    ! .true. drops boundary vertices from the nodal gradient fit entirely, the way
    ! ns_euler_recon_module::compute_nodal_grad does (it guards the fit with
    ! `if (re > 0)` and then sets nodal_grad = 0 at is_bound vertices). The current
    ! path instead MIRRORS the missing neighbour (pr = pl, vr = vl - 2(vl.n)n). The
    ! Taylor-Green wall excess is order-2-only -- 15x the deep-interior error at o2,
    ! none at all at o1, stable over three refinements -- and o1 and o2 share the same
    ! nodal solve, so the o2 gradient reconstruction at the wall is the suspect and
    ! this switch is what isolates it. Skipping the whole boundary VERTEX (rather than
    ! individual boundary sub-faces) subsumes the `re > 0` guard and sidesteps the
    ! boundary_2d z-cap trap: with b2d every vertex has z-cap sub-faces.
    logical, intent(in), optional :: bound_skip
    ! Noeud de bord : GG corrige sur les sous-faces interieures seules (exact pour un champ lineaire),
    ! faces z de b2d en miroir, noeud a matrice singuliere exclu de la moyenne.
    logical, intent(in), optional :: bound_interior, b2d

    integer(kind=ENTIER) :: i, j, id_sub_face, id_sub_elem, id_face, id_elem, idl, idr
    real(kind=DOUBLE) :: omega_p, pl, pr, area, w_p, norm_Sp, p0
    real(kind=DOUBLE) :: w_v, tr_Sgv
    real(kind=DOUBLE) :: eps_weno, eps_loc
    real(kind=DOUBLE), dimension(3) :: S_p, vl_vec, vr_vec, dv_vec, n, dx
    real(kind=DOUBLE), dimension(3, 3) :: S_p_gv, mat, mat_inv
    real(kind=DOUBLE), dimension(:), allocatable :: sum_omega_p, sum_omega_gv
    logical :: use_weno, use_corr, mat_ok, skip_bound
    real(kind=DOUBLE) :: det_mat
    real(kind=DOUBLE), dimension(3, 3) :: mat_inv_t
    logical :: bnd_int, is_b2d, has_wall
    real(kind=DOUBLE), dimension(3) :: bwv, dxg, nw1
    logical :: wcorner
    real(kind=DOUBLE), dimension(3, 3) :: Rg, mat_w, sgv_w

    use_weno = .true.
    if (present(weno)) use_weno = weno
    eps_weno = tiny(1.0_DOUBLE)
    if (present(eps_weno_in)) then
      if (eps_weno_in > 0.0_DOUBLE) eps_weno = eps_weno_in
    end if
    use_corr = .false.
    if (present(gg_corrected)) use_corr = gg_corrected

    grad_v = 0.0_DOUBLE
    grad_p = 0.0_DOUBLE
    div_v  = 0.0_DOUBLE
    if (present(p_limited)) p_limited = 0

    allocate(sum_omega_p(mesh%n_elems), sum_omega_gv(mesh%n_elems))
    sum_omega_p  = 0.0_DOUBLE
    sum_omega_gv = 0.0_DOUBLE

    skip_bound = .false.
    if (present(bound_skip)) skip_bound = bound_skip
    bnd_int = .false.
    if (present(bound_interior)) bnd_int = bound_interior
    is_b2d = .false.
    if (present(b2d)) is_b2d = b2d

    do i = 1, mesh%n_vert
      if (skip_bound) then
        if (mesh%vert(i)%is_bound) cycle
      end if
      omega_p = mesh%vert(i)%volume
      S_p    = 0.0_DOUBLE
      S_p_gv = 0.0_DOUBLE
      mat    = 0.0_DOUBLE
      has_wall = .false.
      bwv = 0.0_DOUBLE; mat_w = 0.0_DOUBLE; sgv_w = 0.0_DOUBLE
      nw1 = 0.0_DOUBLE; wcorner = .false.

      do j = 1, mesh%vert(i)%n_sub_faces_neigh
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        idl    = mesh%sub_face(id_sub_face)%left_elem_neigh
        idr    = mesh%sub_face(id_sub_face)%right_elem_neigh
        area   = mesh%sub_face(id_sub_face)%area
        n      = mesh%sub_face(id_sub_face)%norm
        if (bnd_int .and. idr <= 0 .and. trim(wall_mode) /= 'symmetric') then
          if (.not. (is_b2d .and. abs(n(3)) > 1.0e-8_DOUBLE)) then
            has_wall = .true.
            cycle
          end if
        end if
        pl     = pressure(sol(:, idl), gamma_arr(idl))
        vl_vec = sol(2:4, idl)

        if (idr > 0) then
          pr     = pressure(sol(:, idr), gamma_arr(idr))
          vr_vec = sol(2:4, idr)
          mat = mat + tensor_product_3(area * n, &
            mesh%elem(idr)%coord - mesh%elem(idl)%coord)
        else
          pr     = pressure(sol(:, idl), gamma_arr(idl))
          vr_vec = sol(2:4, idl) - 2.0_DOUBLE*dot_product(sol(2:4, idl), n)*n
          if (trim(wall_mode) == 'symmetric') then
            ! Maille fantome = image miroir par rapport au mur : decalage 2 ((x_f - x_c).n) n.
            dxg = 2.0_DOUBLE*dot_product(mesh%face(id_face)%coord - mesh%elem(idl)%coord, n)*n
            if (.not. (is_b2d .and. abs(n(3)) > 1.0e-8_DOUBLE)) then
              ! Face de mur = face partagee entre la maille et son image : comptee une fois.
              if (norm2(nw1) == 0.0_DOUBLE) then
                nw1 = n
              else if (norm2(n - dot_product(n, nw1)*nw1) > 1.0e-8_DOUBLE) then
                wcorner = .true.
              end if
              bwv = bwv + area*n
              mat_w = mat_w + tensor_product_3(area * n, dxg)
              sgv_w = sgv_w + area * tensor_product_3(vr_vec - vl_vec, n)
              cycle
            end if
            mat = mat + tensor_product_3(area * n, dxg)
          else
            mat = mat + tensor_product_3(area * n, &
              2.0_DOUBLE*(mesh%face(id_face)%coord - mesh%elem(idl)%coord))
          end if
        end if

        dv_vec = vr_vec - vl_vec
        S_p    = S_p    + area * n * (pr - pl)
        S_p_gv = S_p_gv + area * tensor_product_3(dv_vec, n)
      end do

      ! Coin (deux murs non paralleles) : un seul miroir ne suffit pas, le noeud est exclu de la moyenne.
      if (trim(wall_mode) == 'symmetric' .and. wcorner) cycle
      ! Mur symetrique : gabarit double par le miroir R = I - 2 b b de toutes les autres faces.
      if (trim(wall_mode) == 'symmetric' .and. norm2(bwv) > 0.0_DOUBLE) then
        bwv = bwv/norm2(bwv)
        Rg = -2.0_DOUBLE*tensor_product_3(bwv, bwv)
        Rg(1,1) = Rg(1,1) + 1.0_DOUBLE; Rg(2,2) = Rg(2,2) + 1.0_DOUBLE; Rg(3,3) = Rg(3,3) + 1.0_DOUBLE
        S_p    = S_p + matmul(Rg, S_p)
        S_p_gv = S_p_gv + matmul(Rg, matmul(S_p_gv, Rg)) + sgv_w
        mat    = mat + matmul(Rg, matmul(mat, Rg)) + mat_w
        has_wall = .true.
      end if

      if (use_corr .or. has_wall) then
        ! S_p(i)      = mat(i,j) dp/dx_j            -> grad_p = mat^-1 S_p
        ! S_p_gv(i,j) = J(i,k) mat(j,k)             -> J      = S_p_gv . (mat^-1)^T
        ! Singular/ill-conditioned mat (degenerate vertex stencils) falls back to
        ! the uncorrected form rather than producing garbage.
        det_mat = mat(1,1)*(mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)) &
                - mat(1,2)*(mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1)) &
                + mat(1,3)*(mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1))
        mat_ok = abs(det_mat) > 1.0e-12_DOUBLE*max(omega_p, tiny(1.0_DOUBLE))**3
        if (mat_ok) then
          call inverse_3_by_3(mat, mat_inv)
          mat_inv_t = transpose(mat_inv)
          S_p    = matmul(mat_inv, S_p)
          S_p_gv = matmul(S_p_gv, mat_inv_t)
        else if (has_wall) then
          cycle
        else
          S_p    = S_p    / omega_p
          S_p_gv = S_p_gv / omega_p
        end if
      else
        S_p    = S_p    / omega_p
        S_p_gv = S_p_gv / omega_p
      end if
      norm_Sp = norm2(S_p)
      tr_Sgv = S_p_gv(1,1) + S_p_gv(2,2) + S_p_gv(3,3)

      do j = 1, mesh%vert(i)%n_sub_elems_neigh
        id_sub_elem = mesh%vert(i)%sub_elem_neigh(j)
        id_elem     = mesh%sub_elem(id_sub_elem)%mesh_elem
        ! eps = grad_eps_vol*V_c, commun aux noeuds de la maille : poids proportionnel a Omega_p
        ! dans tous les regimes (avec eps ~ Omega_p, un champ 1D devenait 2D au mur symmetric).
        eps_loc = eps_weno
        if (grad_eps_vol > 0.0_DOUBLE) eps_loc = grad_eps_vol*mesh%elem(id_elem)%volume
        if (use_weno) then
          w_p = omega_p / (eps_loc + norm_Sp**4)
        else
          w_p = omega_p
        end if
        ! Poids de la VITESSE : celui de la pression, ou indexe sur |trace(grad v)| en mode 'divv'.
        w_v = w_p
        if (use_weno .and. trim(grad_weight_v) == 'divv') w_v = omega_p / (eps_loc + abs(tr_Sgv)**4)
        grad_p(:, id_elem)    = grad_p(:, id_elem)    + w_p * S_p
        grad_v(:, :, id_elem) = grad_v(:, :, id_elem) + w_v * S_p_gv
        sum_omega_p(id_elem)  = sum_omega_p(id_elem)  + w_p
        sum_omega_gv(id_elem) = sum_omega_gv(id_elem) + w_v
      end do
    end do

    do i = 1, mesh%n_elems
      if (sum_omega_p(i)  > 0.0_DOUBLE) grad_p(:, i)    = grad_p(:, i)    / sum_omega_p(i)
      if (sum_omega_gv(i) > 0.0_DOUBLE) grad_v(:, :, i) = grad_v(:, :, i) / sum_omega_gv(i)
      div_v(i) = grad_v(1, 1, i) + grad_v(2, 2, i) + grad_v(3, 3, i)
    end do

    deallocate(sum_omega_p, sum_omega_gv)

    ! Positivity limiter: if the full reconstruction of p at any vertex gives p ≤ 0,
    ! zero out the gradients so that reconstruction falls back to cell-centred p > 0.
    do i = 1, mesh%n_vert
      do j = 1, mesh%vert(i)%n_sub_elems_neigh
        id_sub_elem = mesh%vert(i)%sub_elem_neigh(j)
        id_elem     = mesh%sub_elem(id_sub_elem)%mesh_elem
        if (mesh%elem(id_elem)%is_ghost) cycle
        dx = mesh%vert(i)%coord - mesh%elem(id_elem)%coord
        p0 = pressure(sol(:, id_elem), gamma_arr(id_elem))
        ! Full reconstruction: p0 + grad_p·dx - 0.5*dt*(γ*p0)*div_v
        if (p0 + dot_product(grad_p(:, id_elem), dx) &
            - 0.5_DOUBLE*rec_dt_fac*dt*gamma_arr(id_elem)*p0*div_v(id_elem) <= 0.0_DOUBLE) then
          grad_p(:, id_elem)    = 0.0_DOUBLE
          grad_v(:, :, id_elem) = 0.0_DOUBLE
          div_v(id_elem)        = 0.0_DOUBLE
          if (present(p_limited)) p_limited(id_elem) = 1
        end if
      end do
    end do
  end subroutine compute_gradients

  subroutine compute_rhs_lagrange_classic_iso(mesh, sol, vp, dt, rhs, n_bc, bc_type, &
      bc_val, b2d, mass, gamma_arr, vp_is_imposed, iso_weight_mode, h_extrude)
    ! "classic_iso" family: same node-based nodal-velocity solve as the classic
    ! scheme, but every INTERNAL in-plane sub-face flux is blended with a 1D
    ! acoustic face flux:
    !
    !     F = w_pcf * F_classic + (1 - w_pcf) * F_face
    !
    ! - w_pcf(j) = scale_a / area_j  (clipped to [0,1]) over the sub-faces around
    !   the node, in the spirit of the Euler scheme multi_point_iso.
    !     iso_weight_mode = 1  ("classic_iso")     : scale_a = min_k area_k
    !     iso_weight_mode = 2  ("classic_iso_vol") : scale_a = (min_k V_sub,k)^(2/3)
    !       with V_sub,k the sub-element volumes around the vertex.
    ! - F_face uses classic 1D Lagrangian star values (p_f, v_f) from an acoustic
    !   Riemann solver, as done for the Lagrange part of the LS1D schemes.
    ! Boundary sub-faces and, in boundary_2d, the top/bottom z-faces keep the pure
    ! classic flux (w_pcf = 1, unchanged BC handling). The nodal system uses the
    ! matching weight A_j * w_pcf(j) so it stays the discrete adjoint of the flux.
    use linear_solver_module, only: lu_solve, print_mat, inverse_3_by_3
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), intent(in) :: dt
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mass
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    logical, dimension(mesh%n_vert), intent(in) :: vp_is_imposed
    integer(kind=ENTIER), intent(in) :: n_bc
    character(len=255), dimension(n_bc) :: bc_type
    real(kind=DOUBLE), dimension(5, n_bc) :: bc_val
    logical, intent(in) :: b2d
    integer(kind=ENTIER), intent(in) :: iso_weight_mode
    real(kind=DOUBLE), intent(in) :: h_extrude   ! extrusion height (b2d meshes)

    integer(kind=ENTIER) :: i, j, id_face, id_sub_face
    integer(kind=ENTIER) :: nsfn, idl, idr, var
    real(kind=DOUBLE), dimension(:, :), allocatable :: lambda
    real(kind=DOUBLE), dimension(:, :, :), allocatable :: sol_lr
    real(kind=DOUBLE), dimension(:), allocatable :: w_pcf
    real(kind=DOUBLE), dimension(5) :: flux, flux_c, flux_f
    real(kind=DOUBLE), dimension(3, 3) :: Mp, Mp_inv
    real(kind=DOUBLE), dimension(3) :: Rp, Bp
    real(kind=DOUBLE) :: vl, vr, pl, pr, v_bar, pl_et, pr_et, v_et, Pp
    real(kind=DOUBLE) :: gl, gr
    real(kind=DOUBLE) :: tau_l, tau_r, a_l, a_r, rho_m, a_m, p_f, v_f, wj, scale_a, aw
    real(kind=DOUBLE) :: p_w
    integer(kind=ENTIER) :: n_ip
    logical :: is_corner, have_first_wall_norm
    real(kind=DOUBLE), dimension(3) :: first_wall_norm
    real(kind=DOUBLE) :: cross_z

    rhs = 0.0_DOUBLE
    do i=1, mesh%n_vert
      if (.not. vp_is_imposed(i)) vp(:, i) = 0.0_DOUBLE
    end do

    do i=1, mesh%n_vert
      nsfn = mesh%vert(i)%n_sub_faces_neigh
      allocate(lambda(2, nsfn))
      lambda = 0.0_DOUBLE
      allocate(sol_lr(5, 2, nsfn))
      sol_lr = 0.0_DOUBLE
      allocate(w_pcf(nsfn))

      ! ---- Sub-face weight w_pcf(j) = min(1, (scale_a / A_j)^p_w) --------------
      ! Recipes for the reference area scale scale_a (and exponent p_w), chosen
      ! by iso_weight_mode.  In boundary_2d the top/bottom z-faces are excluded
      ! from every statistic and given w_pcf = 1 (pure classic).
      !   1  min_k A_k                                (p=1)  "classic_iso"
      !   2  b2d: sqrt(A_2D_min)*h ; 3D: V_min^(2/3)   (p=1)  "classic_iso_vol"
      !   3  arithmetic mean of A_k                   (p=1)
      !   4  min_k A_k                                (p=1/2)
      !   5  (min_k A_k / max_k A_k)                  (p=1)   node aspect ratio
      !   6  harmonic mean of A_k                     (p=1)
      !   7  b2d: sqrt(A_2D_min)*h ; 3D: V_min^(2/3)   (p=1/2)
      p_w = 1.0_DOUBLE
      if (iso_weight_mode == 2 .or. iso_weight_mode == 7) then
        scale_a = huge(1.0_DOUBLE)
        do j = 1, mesh%vert(i)%n_sub_elems_neigh
          scale_a = min(scale_a, &
            abs(mesh%sub_elem(mesh%vert(i)%sub_elem_neigh(j))%volume))
        end do
        if (b2d) then
          scale_a = sqrt(scale_a / h_extrude) * h_extrude
        else
          scale_a = scale_a ** (2.0_DOUBLE/3.0_DOUBLE)
        end if
        if (iso_weight_mode == 7) p_w = 0.5_DOUBLE
      else
        w_pcf = huge(1.0_DOUBLE)
        n_ip = 0
        do j=1, nsfn
          id_sub_face = mesh%vert(i)%sub_face_neigh(j)
          id_face = mesh%sub_face(id_sub_face)%mesh_face
          if( b2d .and. abs(mesh%face(id_face)%norm(3)) > 1e-8_DOUBLE ) cycle
          w_pcf(j) = mesh%sub_face(id_sub_face)%area
          n_ip = n_ip + 1
        end do
        select case (iso_weight_mode)
        case (3)
          scale_a = sum(w_pcf, mask = (w_pcf < huge(1.0_DOUBLE))) / max(1, n_ip)
        case (5)
          scale_a = minval(w_pcf) &
            / maxval(w_pcf, mask = (w_pcf < huge(1.0_DOUBLE)))
        case (6)
          scale_a = real(max(1, n_ip), DOUBLE) &
            / sum(1.0_DOUBLE / w_pcf, mask = (w_pcf < huge(1.0_DOUBLE)))
        case (4)
          scale_a = minval(w_pcf)
          p_w = 0.5_DOUBLE
        case default   ! 1
          scale_a = minval(w_pcf)
        end select
      end if
      do j=1, nsfn
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        idr = mesh%face(id_face)%right_neigh
        if( idr > 0 .and. .not. (b2d .and. &
            abs(mesh%face(id_face)%norm(3)) > 1e-8_DOUBLE) ) then
          if (iso_weight_mode == 5) then
            w_pcf(j) = min(1.0_DOUBLE, scale_a)
          else
            w_pcf(j) = min(1.0_DOUBLE, &
              (scale_a / mesh%sub_face(id_sub_face)%area) ** p_w)
          end if
        else
          ! boundary faces and z-faces: pure classic flux, keep full area weight.
          w_pcf(j) = 1.0_DOUBLE
        end if
      end do

      Bp = 0.0_DOUBLE
      Mp = 0.0_DOUBLE
      Rp = 0.0_DOUBLE
      is_corner = .false.
      have_first_wall_norm = .false.
      first_wall_norm = 0.0_DOUBLE
      do j=1, nsfn
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        idl = mesh%face(id_face)%left_neigh
        idr = mesh%face(id_face)%right_neigh
        gl = gamma_arr(idl)
        sol_lr(:, 1, j) = sol(:, idl)
        pl = pressure(sol_lr(:, 1, j), gl)
        vl = dot_product(sol_lr(2:4, 1, j), mesh%face(id_face)%norm)
        ! The nodal system is the discrete adjoint of the ACTUAL sub-face flux.
        ! The classic (multipoint) part of that flux carries the factor
        ! A_j * w_pcf(j); the (1-w_pcf) face-flux part cancels between left and
        ! right and drops out. So every A_j below is replaced by A_j * w_pcf(j)
        ! -> on internal in-plane faces the area cancels (constant min_area),
        ! on boundary / z-faces w_pcf(j) = 1 and the area weight is unchanged.
        aw = mesh%sub_face(id_sub_face)%area * w_pcf(j)

        if( idr > 0 ) then
          gr = gamma_arr(idr)
          sol_lr(:, 2, j) = sol(:, idr)
          pr = pressure(sol_lr(:, 2, j), gr)
          vr = dot_product(sol_lr(2:4, 2, j), mesh%face(id_face)%norm)

          lambda(1, j) = max(sqrt(gl*pl*sol_lr(1, 1, j))/sol_lr(1, 1, j), &
            sqrt(max(0.0_DOUBLE, pr-pl)/sol_lr(1, 1, j)), &
            -(vr-vl)/sol_lr(1, 1, j))
          lambda(2, j) = max(sqrt(gr*pr*sol_lr(1, 2, j))/sol_lr(1, 2, j), &
            sqrt(max(0.0_DOUBLE, pl-pr)/sol_lr(1, 2, j)), &
            -(vr-vl)/sol_lr(1, 2, j))

          if (.not. vp_is_imposed(i)) then
            v_bar = (lambda(1, j)*vl + lambda(2, j)*vr - (pr - pl)) &
              /(lambda(2, j) + lambda(1, j))
            Mp = Mp + aw*(lambda(1,j)+lambda(2,j)) &
              * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
            Rp = Rp + aw*(lambda(1,j)+lambda(2,j)) &
              *mesh%face(id_face)%norm*v_bar
          end if

        else
          if( b2d .and. abs(mesh%face(id_face)%norm(3)) > 1e-8_DOUBLE ) then
            lambda(1, j) = sqrt(gl*pl*sol_lr(1, 1, j))/sol_lr(1, 1, j)
            if (.not. vp_is_imposed(i)) then
              Mp = Mp + aw*2*lambda(1,j) &
                * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
            end if
          else if (idr < 0 .and. trim(bc_type(-idr)) == 'pressure') then
            lambda(1, j) = sqrt(gl*pl*sol_lr(1, 1, j))/sol_lr(1, 1, j)
            if (.not. vp_is_imposed(i)) then
              Mp = Mp + aw*lambda(1,j) &
                * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
              ! Etat fantome p_r = p_bc, v_r = v_l : avec Mp = A lambda, la solution a une face
              ! doit etre v.n = v_l + (p_l - p_bc)/lambda. p_l manquait, p_bc avait le signe inverse.
              Rp = Rp + aw &
                *(pl - bc_val(1, -idr) + lambda(1,j)*vl) &
                *mesh%face(id_face)%norm
            end if
          else
            lambda(1, j) = sqrt(gl*pl*sol_lr(1, 1, j))/sol_lr(1, 1, j)
            if (.not. vp_is_imposed(i)) then
              if (.not. have_first_wall_norm) then
                first_wall_norm = mesh%face(id_face)%norm
                have_first_wall_norm = .true.
              else
                cross_z = first_wall_norm(1)*mesh%face(id_face)%norm(2) &
                        - first_wall_norm(2)*mesh%face(id_face)%norm(1)
                if (abs(cross_z) > 1.0e-8_DOUBLE) is_corner = .true.
              end if
              Bp = Bp + aw*mesh%face(id_face)%norm
              Mp = Mp + aw*lambda(1,j) &
                * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
              Rp = Rp + aw*(pl + lambda(1,j)*vl) &
                *mesh%face(id_face)%norm
            end if
          end if
        end if
      end do

      if (.not. vp_is_imposed(i)) then
        if (is_corner) then
          vp(:, i) = 0.0_DOUBLE
        else if( maxval(abs(Bp)) > 1e-8_DOUBLE ) then
          call inverse_3_by_3(Mp, Mp_inv)
          Pp = dot_product(Rp, matmul(Mp_inv,Bp))/dot_product(Bp, matmul(Mp_inv, Bp))
          vp(:, i) = matmul(Mp_inv, Rp-Pp*BP)
        else
          call inverse_3_by_3(Mp, Mp_inv)
          vp(:, i) = matmul(Mp_inv, Rp)
        end if
      end if

      do j=1, nsfn
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        idl = mesh%face(id_face)%left_neigh
        idr = mesh%face(id_face)%right_neigh
        gl = gamma_arr(idl)

        v_et = dot_product(vp(:, i), mesh%face(id_face)%norm)
        vl = dot_product(sol_lr(2:4, 1, j), mesh%face(id_face)%norm)
        pl = pressure(sol_lr(:, 1, j), gl)
        pl_et = pl - lambda(1, j)*(v_et - vl)

        ! Blend the 1D acoustic face flux only on INTERNAL IN-PLANE sub-faces.
        ! Boundary faces and, in boundary_2d, the top/bottom z-faces stay pure
        ! classic (wj = 1), so their handling is identical to compute_rhs_lagrange.
        wj = 1.0_DOUBLE
        flux_f = 0.0_DOUBLE
        if( idr > 0 .and. .not. (b2d .and. &
            abs(mesh%face(id_face)%norm(3)) > 1e-8_DOUBLE) ) then
          wj = w_pcf(j)
          gr = gamma_arr(idr)
          vr = dot_product(sol_lr(2:4, 2, j), mesh%face(id_face)%norm)
          pr = pressure(sol_lr(:, 2, j), gr)
          tau_l = sol_lr(1, 1, j)
          tau_r = sol_lr(1, 2, j)
          a_l = sqrt(max(0.0_DOUBLE, gl*pl*tau_l))
          a_r = sqrt(max(0.0_DOUBLE, gr*pr*tau_r))
          rho_m = 0.5_DOUBLE*(1.0_DOUBLE/tau_l + 1.0_DOUBLE/tau_r)
          a_m   = max(a_l, a_r, 1.0e-10_DOUBLE)
          v_f = 0.5_DOUBLE*(vl + vr) - (pr - pl)/(2.0_DOUBLE*rho_m*a_m)
          p_f = 0.5_DOUBLE*(pl + pr) - 0.5_DOUBLE*rho_m*a_m*(vr - vl)
          flux_f(1)   = -v_f
          flux_f(2:4) = p_f*mesh%face(id_face)%norm
          flux_f(5)   = p_f*v_f
        end if

        ! left cell
        flux_c(1)   = -v_et
        flux_c(2:4) = pl_et*mesh%face(id_face)%norm
        flux_c(5)   = pl_et*v_et
        flux = wj*flux_c + (1.0_DOUBLE - wj)*flux_f
        do var=1, 5
          rhs(var, idl) = rhs(var, idl) &
            - mesh%sub_face(id_sub_face)%area/mass(idl) * flux(var)
        end do

        if( idr > 0 ) then
          gr = gamma_arr(idr)
          vr = dot_product(sol_lr(2:4, 2, j), mesh%face(id_face)%norm)
          pr = pressure(sol_lr(:, 2, j), gr)
          pr_et = pr + lambda(2, j)*(v_et - vr)
          flux_c(1)   = -v_et
          flux_c(2:4) = pr_et*mesh%face(id_face)%norm
          flux_c(5)   = pr_et*v_et
          flux = wj*flux_c + (1.0_DOUBLE - wj)*flux_f
          do var=1, 5
            rhs(var, idr) = rhs(var, idr) &
              + mesh%sub_face(id_sub_face)%area/mass(idr) * flux(var)
          end do
        end if
      end do

      deallocate(lambda, sol_lr, w_pcf)
    end do
  end subroutine compute_rhs_lagrange_classic_iso

  subroutine compute_rhs_lagrange_sidil(mesh, sol, vp, dt, rhs, &
      n_bc, bc_type, bc_val, b2d, mass, method, b2d_h, gamma_arr, vp_is_imposed, h_p_arr)
    use linear_solver_module, only: lu_solve, print_mat, inverse_3_by_3
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), intent(in) :: dt
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mass
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    logical, dimension(mesh%n_vert), intent(in) :: vp_is_imposed
    integer(kind=ENTIER), intent(in) :: n_bc, method
    character(len=255), dimension(n_bc) :: bc_type
    real(kind=DOUBLE), dimension(5, n_bc) :: bc_val
    logical, intent(in) :: b2d
    real(kind=DOUBLE), intent(in) :: b2d_h
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(out), optional :: h_p_arr

    integer(kind=ENTIER) :: i, j, id_elem, id_face, id_sub_face
    integer(kind=ENTIER) :: id_sub_elem
    integer(kind=ENTIER) :: nsfn, idl, idr, var
    real(kind=DOUBLE), dimension(:, :), allocatable :: lambda
    real(kind=DOUBLE), dimension(:, :, :), allocatable :: sol_lr
    real(kind=DOUBLE), dimension(5) :: flux, sol_w
    real(kind=DOUBLE), dimension(3, 3) :: Mp, Mp_inv
    real(kind=DOUBLE), dimension(3) :: Rp, Bp
    real(kind=DOUBLE) :: vl, vr, pl, pr, v_bar, pl_et, pr_et, v_et, pp, a_p
    logical :: corner

    rhs = 0.0_DOUBLE
    do i=1, mesh%n_vert
      if (.not. vp_is_imposed(i)) vp(:, i) = 0.0_DOUBLE
    end do

    do i=1, mesh%n_vert

      if (.not. vp_is_imposed(i)) then
        call compute_nodal_velocity_sidil(mesh, i, sol, vp(:, i), method, b2d, b2d_h, gamma_arr, &
          n_bc, bc_type, bc_val)
      end if
      call compute_nodal_pressure_sidil(mesh, i, sol, pp, method, b2d, b2d_h, gamma_arr, &
        n_bc, bc_type, bc_val)
      if (present(h_p_arr)) h_p_arr(i) = compute_length(mesh, i, method, b2d, b2d_h)

      a_p = 0.0_DOUBLE
      do j=1, mesh%vert(i)%n_sub_elems_neigh
        id_sub_elem = mesh%vert(i)%sub_elem_neigh(j)
        id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
        sol_w = lag_to_primit(sol(:, id_elem), gamma_arr(id_elem))
        a_p = a_p + mesh%sub_elem(id_sub_elem)%volume &
          * sqrt(gamma_arr(id_elem)*sol_w(5)/sol_w(1))
      end do
      a_p = a_p / mesh%vert(i)%volume

      nsfn = mesh%vert(i)%n_sub_faces_neigh
      do j=1, nsfn
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        idl = mesh%face(id_face)%left_neigh
        idr = mesh%face(id_face)%right_neigh

        v_et = dot_product(vp(:, i), mesh%face(id_face)%norm)

        flux(1)   = -v_et
        flux(2:4) = pp*mesh%face(id_face)%norm
        flux(5)   = pp*v_et
        do var=1, 5
          rhs(var, idl) = rhs(var, idl) &
            - mesh%sub_face(id_sub_face)%area/mass(idl) * flux(var)
        end do

        if( idr > 0 ) then
          flux(1)   = -v_et
          flux(2:4) = pp*mesh%face(id_face)%norm
          flux(5)   = pp*v_et
          do var=1, 5
            rhs(var, idr) = rhs(var, idr) &
              + mesh%sub_face(id_sub_face)%area/mass(idr) * flux(var)
          end do
        end if
      end do

    end do
  end subroutine compute_rhs_lagrange_sidil

  ! Extrapolation GRP d'un demi-etat de la moyenne de maille vers le NOEUD, avec le
  ! demi-pas de temps. Meme forme que dans 'classic' (compute_rhs_around_node), ou
  ! elle est ecrite deux fois en ligne ; ici elle sert quatre fois (systeme de
  ! vitesse et systeme de pression, cote gauche et cote droit), d'ou la routine.
  ! a est evalue sur la valeur MOYENNE de p, avant reconstruction, comme dans
  ! 'classic'. p et v sont modifies sur place.
  subroutine grp_half_step(mesh, ic, i_vert, nrm, dt, tau, g, grad_v, grad_p, div_v, p, v)
    implicit none
    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: ic, i_vert
    real(kind=DOUBLE), dimension(3), intent(in) :: nrm
    real(kind=DOUBLE), intent(in) :: dt, tau, g
    real(kind=DOUBLE), dimension(3, 3, mesh%n_elems), intent(in) :: grad_v
    real(kind=DOUBLE), dimension(3, mesh%n_elems), intent(in) :: grad_p
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: div_v
    real(kind=DOUBLE), intent(inout) :: p, v
    real(kind=DOUBLE), dimension(3) :: dx
    real(kind=DOUBLE) :: a

    a  = sqrt(g*p*tau)
    dx = mesh%vert(i_vert)%coord - mesh%elem(ic)%coord
    v = v + dot_product(matmul(grad_v(:,:,ic), dx), nrm) &
          - 0.5_DOUBLE*rec_dt_fac*dt*tau*dot_product(grad_p(:, ic), nrm)
    p = p + dot_product(grad_p(:, ic), dx) - 0.5_DOUBLE*rec_dt_fac*dt*a**2/tau*div_v(ic)
  end subroutine grp_half_step

  ! ---------------------------------------------------------------------------
  ! 'vp_pp' (anciennement 'nodal_vp'): nodal velocity AND nodal pressure, both obtained from NODAL
  ! SYSTEMS, so that no nodal length h_p is needed anywhere.
  !
  ! Motivation. 'sidil' gets its two nodal quantities by correcting volume
  ! averages with length-scaled terms,
  !     v_p = <v> - (h_p/(2 rho_p a_p)) grad_p,   p_p = <p> - (h_p/2) rho_p a_p div_v,
  ! which makes the scheme depend on an ad-hoc h_p (compute_length). Here both
  ! come out of the same Riemann-solver data as the fluxes:
  !
  !   velocity   [sum_f A_f (lam_l+lam_r) n(x)n] u_p = sum_f A_f (lam_l+lam_r) vbar_f n
  !              vbar_f = (lam_l v_l + lam_r v_r - (p_r-p_l))/(lam_l+lam_r)
  !              i.e. the EUCCLHYD nodal solver, as in 'classic' (sub-face
  !              based, one impedance per side of each sub-face -- not GLACE,
  !              which builds corner vectors instead).
  !
  !   pressure   p_p = sum_f A_f L_f pbar_f / sum_f A_f L_f,   L_f = 1/lam_l + 1/lam_r
  !              pbar_f = (p_l/lam_l + p_r/lam_r - (v_r-v_l)) / L_f
  !              i.e. the impedance-weighted nodal pressure of the multi-point
  !              pressure solver (Del Grosso, Barsukow, Loubere & Maire 2026),
  !              whose velocity-jump term is the multi-d divergence.
  !
  ! Unlike 'classic', which closes with the one-sided pressures
  ! p_l - lam_l (u_p.n - v_l), a SINGLE p_p is used on both sides of a sub-face,
  ! so the flux is symmetric and the scheme conservative by construction.
  !
  ! Boundary faces carry no right state, so they enter the pressure average with
  ! weight 1/lam_l and the one-sided value p_l - lam_l (u_p.n - v_l), evaluated
  ! after the velocity system has been solved. b2d z-normal faces are excluded
  ! from the pressure average entirely (they are not physical faces).
  ! Bornes de positivite du flux d'energie, forme UNIVALUEE. Isolee dans sa propre
  ! routine pour que le controle en dur teste EXACTEMENT le code que le solveur
  ! execute : un controle qui reimplemente la formule ne verifie rien.
  ! n est la normale SORTANTE de la maille c, et u_c comme u_d sont projetes sur
  ! CETTE normale, pas chacun sur la sienne.
  pure subroutine pos_sv_bounds(lam_c, eps_c, u_c, lam_d, eps_d, u_d, pp, &
                                phi_hi, phi_lo, d2)
    implicit none
    real(kind=DOUBLE), intent(in) :: lam_c, eps_c, u_c, lam_d, eps_d, u_d, pp
    real(kind=DOUBLE), intent(out) :: phi_hi, phi_lo, d2
    phi_hi =  lam_c*eps_c + pp*u_c
    phi_lo = -lam_d*eps_d + pp*u_d
    d2     = phi_hi - phi_lo
  end subroutine pos_sv_bounds

  ! Le flux retenu : naturel s'il est admissible, milieu sinon. UNE valeur, donc
  ! conservatif. 'corr' dit si la correction a mordu, 'infeas' si l'intervalle est
  ! vide -- auquel cas on ne corrige RIEN et on compte, parce qu'un intervalle vide
  ! est une anomalie a signaler et pas a absorber.
  pure subroutine pos_sv_flux(t_nat, phi_hi, phi_lo, d2, always_mid, fl, corr, &
                              infeas, isnan)
    implicit none
    real(kind=DOUBLE), intent(in) :: t_nat, phi_hi, phi_lo, d2
    logical, intent(in) :: always_mid
    real(kind=DOUBLE), intent(out) :: fl
    logical, intent(out) :: corr, infeas, isnan
    corr = .false.; infeas = .false.; isnan = .false.
    if (t_nat /= t_nat) isnan = .true.
    if (d2 <= 0.0_DOUBLE) infeas = .true.
    if (always_mid) then
      ! MILIEU PUR : aucune branche, aucun test d'appartenance, le flux naturel
      ! n'est PAS lu. Les deux drapeaux ci-dessus restent leves pour le comptage,
      ! mais ils ne gouvernent rien -- et c'est la difference de fond avec la
      ! version a limiteur : prendre le milieu sans consulter un flux naturel NaN
      ! n'est pas un blanchiment, c'est ne pas utiliser une entree cassee. Le NaN
      ! reste visible par le compteur, ce qui etait tout l'enjeu.
      fl = 0.5_DOUBLE*(phi_hi + phi_lo); corr = .true.; return
    end if
    ! Version a LIMITEUR. Ici le flux naturel EST lu, donc un NaN doit passer au
    ! lieu d'etre remplace par un milieu fini : 219 des 295 "corrections" mesurees
    ! sur sedov_tri etaient exactement ce blanchiment.
    if (isnan .or. infeas) then
      fl = t_nat; return
    end if
    if (t_nat > phi_lo .and. t_nat < phi_hi) then
      fl = t_nat
    else
      fl = 0.5_DOUBLE*(phi_hi + phi_lo); corr = .true.
    end if
  end subroutine pos_sv_flux

  ! Les quatre controles de la construction, en dur, sur des etats fabriques.
  ! RESERVE a lire avec les resultats : p_p est ici remplace par la pression
  ! d'interface acoustique 1D, alors que dans vp_pp c'est la sortie d'un systeme
  ! NODAL. Les controles 1 et 2 sont algebriques et ne dependent pas de ce choix ;
  ! les controles 3 et 4 en dependent, et leur version qui compte est le compte de
  ! sous-faces corrigees mesure sur les vrais calculs.
  subroutine sv_step_reset()
    implicit none
    sv_step_n = 0; sv_step_r = 0.0_DOUBLE
    sv_step_rmin = 1.0e300_DOUBLE; sv_step_rmax = 0.0_DOUBLE
    sv_st_infeas = 0; sv_st_d2min = 1.0e300_DOUBLE
  end subroutine sv_step_reset

  subroutine pos_sv_selfcheck()
    implicit none
    real(kind=DOUBLE) :: g, rc, pc, rd, pd, zc, zd, ec, ed, uc, ud, pp, vn
    real(kind=DOUBLE) :: hi, lo, d2, hi2, lo2, d22, fl, w, ma, ac, rap
    real(kind=DOUBLE) :: d2min, ecart
    integer(kind=ENTIER) :: k, nout
    logical :: corr, infeas, isnan
    g = 1.4_DOUBLE

    print*, "POSSV --- controle 1 : largeur de l'intervalle sur un saut de pression"
    d2min = 1.0e300_DOUBLE
    do k = 1, 6
      rap = 10.0_DOUBLE**real(k, DOUBLE)
      rc = 1.0_DOUBLE; pc = rap; rd = 1.0_DOUBLE; pd = 1.0_DOUBLE
      zc = sqrt(g*pc*rc); zd = sqrt(g*pd*rd)
      ec = pc/((g - 1.0_DOUBLE)*rc); ed = pd/((g - 1.0_DOUBLE)*rd)
      uc = 0.0_DOUBLE; ud = 0.0_DOUBLE
      pp = (zd*pc + zc*pd)/(zc + zd)
      call pos_sv_bounds(zc, ec, uc, zd, ed, ud, pp, hi, lo, d2)
      d2min = min(d2min, d2)
      print*, "POSSV   rapport", rap, " Phi_haut", hi, " Phi_bas", lo, " delta2", d2
    end do
    print*, "POSSV   largeur MINIMALE =", d2min, &
            "  -- doit etre > 0, sinon intervalle vide"

    print*, "POSSV --- controle 2 : invariance galileenne de delta2"
    rc = 1.0_DOUBLE; pc = 100.0_DOUBLE; rd = 1.0_DOUBLE; pd = 1.0_DOUBLE
    zc = sqrt(g*pc*rc); zd = sqrt(g*pd*rd)
    ec = pc/((g - 1.0_DOUBLE)*rc); ed = pd/((g - 1.0_DOUBLE)*rd)
    pp = (zd*pc + zc*pd)/(zc + zd)
    call pos_sv_bounds(zc, ec, 0.3_DOUBLE, zd, ed, -0.2_DOUBLE, pp, hi, lo, d2)
    do k = 1, 4
      select case (k)
      case (1); w = 0.0_DOUBLE
      case (2); w = 1.0_DOUBLE
      case (3); w = 5.0_DOUBLE
      case (4); w = 20.0_DOUBLE
      end select
      call pos_sv_bounds(zc, ec, 0.3_DOUBLE + w, zd, ed, -0.2_DOUBLE + w, pp, &
                         hi2, lo2, d22)
      print*, "POSSV   w =", w, " delta2 =", d22, " ecart a w=0 :", abs(d22 - d2)
    end do

    print*, "POSSV --- controle 3 : bas Mach, le flux naturel doit etre admissible"
    ! RESERVE SUR LA PREMIERE VERSION DE CE CONTROLE, gardee ici en garde-fou :
    ! avec p_c = p_d, rho_c = rho_d et u_d = -u_c, on a milieu = p_p (u_c+u_d)/2 = 0
    ! et flux naturel = p_p v_n = 0, donc les deux coincident IDENTIQUEMENT et le
    ! controle sort 0 corrigee quelle que soit la construction. Il n'avait aucun
    ! pouvoir de detection. La configuration ci-dessous est asymetrique des deux
    ! cotes -- pressions, densites et vitesses differentes -- pour que l'ecart au
    ! milieu soit non nul et porte la vraie dependance en Mach.
    nout = 0
    do k = 1, 5
      select case (k)
      case (1); ma = 0.5_DOUBLE
      case (2); ma = 0.2_DOUBLE
      case (3); ma = 0.05_DOUBLE
      case (4); ma = 0.01_DOUBLE
      case (5); ma = 0.002_DOUBLE
      end select
      rc = 1.0_DOUBLE; pc = 1.0_DOUBLE
      rd = 1.3_DOUBLE; pd = 1.0_DOUBLE + 0.7_DOUBLE*g*ma**2
      ac = sqrt(g*pc/rc)
      zc = sqrt(g*pc*rc); zd = sqrt(g*pd*rd)
      ec = pc/((g - 1.0_DOUBLE)*rc); ed = pd/((g - 1.0_DOUBLE)*rd)
      uc = 0.4_DOUBLE*ma*ac; ud = -1.1_DOUBLE*ma*ac
      pp = (zd*pc + zc*pd)/(zc + zd) - zc*zd/(zc + zd)*(ud - uc)
      vn = (zc*uc + zd*ud - (pd - pc))/(zc + zd)
      call pos_sv_bounds(zc, ec, uc, zd, ed, ud, pp, hi, lo, d2)
      call pos_sv_flux(pp*vn, hi, lo, d2, .false., fl, corr, infeas, isnan)
      ecart = abs(pp*vn - 0.5_DOUBLE*(hi + lo))/max(d2, 1.0e-30_DOUBLE)
      if (corr) nout = nout + 1
      print*, "POSSV   Ma", ma, " corrigee", corr, " ecart au milieu / largeur", ecart
    end do
    print*, "POSSV   sous-faces corrigees au bas Mach :", nout, " sur 5 -- doit valoir 0"

    print*, "POSSV --- controle 4 : choc, la correction doit mordre"
    do k = 1, 4
      select case (k)
      case (1); rap = 1.0e2_DOUBLE
      case (2); rap = 1.0e4_DOUBLE
      case (3); rap = 1.0e6_DOUBLE
      case (4); rap = 1.0e8_DOUBLE
      end select
      rc = 1.0_DOUBLE; pc = rap; rd = 1.0_DOUBLE; pd = 1.0_DOUBLE
      zc = sqrt(g*pc*rc); zd = sqrt(g*pd*rd)
      ec = pc/((g - 1.0_DOUBLE)*rc); ed = pd/((g - 1.0_DOUBLE)*rd)
      uc = 0.0_DOUBLE; ud = 0.0_DOUBLE
      pp = (zd*pc + zc*pd)/(zc + zd)
      vn = (zc*uc + zd*ud - (pd - pc))/(zc + zd)
      call pos_sv_bounds(zc, ec, uc, zd, ed, ud, pp, hi, lo, d2)
      call pos_sv_flux(pp*vn, hi, lo, d2, .false., fl, corr, infeas, isnan)
      ecart = abs(pp*vn - 0.5_DOUBLE*(hi + lo))/max(d2, 1.0e-30_DOUBLE)
      print*, "POSSV   rapport", rap, " flux naturel", pp*vn, " milieu", &
              0.5_DOUBLE*(hi + lo), " corrigee", corr, &
              " ecart au milieu / largeur", ecart
    end do
  end subroutine pos_sv_selfcheck

  subroutine compute_rhs_lagrange_vp_pp(mesh, sol, vp, dt, rhs, &
      n_bc, bc_type, bc_val, b2d, mass, gamma_arr, vp_is_imposed, &
      second_order, grad_v, grad_p, div_v, use_ppvp)
    use linear_solver_module, only: inverse_3_by_3
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), intent(in) :: dt
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mass
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    logical, dimension(mesh%n_vert), intent(in) :: vp_is_imposed
    integer(kind=ENTIER), intent(in) :: n_bc
    character(len=255), dimension(n_bc) :: bc_type
    real(kind=DOUBLE), dimension(5, n_bc) :: bc_val
    logical, intent(in) :: b2d
    logical, intent(in) :: second_order
    real(kind=DOUBLE), dimension(3, 3, mesh%n_elems), intent(in) :: grad_v
    real(kind=DOUBLE), dimension(3, mesh%n_elems), intent(in) :: grad_p
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: div_v
    ! .true. = schema 'vp_pp_ppvp' : une TROISIEME quantite nodale ppvp porte le
    ! flux d'energie a la place du produit pp*(u_p.n) de deux quantites resolues
    ! separement. Construite exactement comme pp : meme ponderation en 1/lambda,
    ! meme exclusion des faces z, meme etat miroir aux bords.
    logical, intent(in) :: use_ppvp

    integer(kind=ENTIER) :: i, j, id_elem, id_face, id_sub_face, nsfn
    integer(kind=ENTIER) :: idl, idr, var
    real(kind=DOUBLE), dimension(3, 3) :: Mp, Mp_inv
    real(kind=DOUBLE), dimension(3) :: Rp, Bp, first_wall_norm
    real(kind=DOUBLE), dimension(5) :: flux
    real(kind=DOUBLE) :: vl, vr, pl, pr, taul, taur, gl, gr
    real(kind=DOUBLE) :: lam_l, lam_r, v_bar, p_bar, invlam
    real(kind=DOUBLE) :: pp, s_num, s_den, v_et, Pw, cross_z, area, v_target
    real(kind=DOUBLE) :: ppvp, q_num, q_den, q_bar, ql, qr, v_node, xjump
    real(kind=DOUBLE) :: psvp_pil, psvp_pir, psvp_t
    real(kind=DOUBLE), dimension(3,3) :: Gp, Gp_inv, MinvG
    real(kind=DOUBLE), dimension(3) :: Sp, dxc
    real(kind=DOUBLE) :: s_gc, detG
    integer(kind=ENTIER) :: nd_gc
    integer(kind=ENTIER) :: psvp_it, psvp_nit
    real(kind=DOUBLE) :: psvp_epss, psvp_lreq
    real(kind=DOUBLE), dimension(:), allocatable :: ppvp_face
    real(kind=DOUBLE), dimension(:), allocatable :: pv_l_f, pv_r_f
    real(kind=DOUBLE) :: sv_hi, sv_lo, sv_d2, sv_fl, sv_ecart, sv_ray
    logical :: sv_corr, sv_infeas, sv_isnan
    real(kind=DOUBLE) :: d_l, d_r, q_visc, tr_d, h_f, mu_v
    real(kind=DOUBLE), dimension(3, 3) :: G_t, D_t
    real(kind=DOUBLE), dimension(3) :: tau_n
    real(kind=DOUBLE) :: eps_l, eps_r, B_l, B_r, cap_l, flo_r, t_f, p_bar_f, v_bar_f, E_l_s, E_r_s
    real(kind=DOUBLE) :: av_l, av_r, room, S_sur, R_room, reduc, flux5_l, flux5_r
    integer(kind=ENTIER) :: n_clip_step
    logical :: is_corner, have_first_wall_norm, is_zface

    rhs = 0.0_DOUBLE
    if (diag_entropy) then
      if (.not. allocated(entropy_prod)) then
        allocate(entropy_prod(mesh%n_elems), entropy_ref(mesh%n_elems))
      end if
      entropy_prod = 0.0_DOUBLE
      entropy_ref  = 0.0_DOUBLE
    end if
    if (pvp_pos) then
      if (.not. allocated(pos_clip_cell)) then
        allocate(pos_clip_cell(mesh%n_elems)); pos_clip_cell = 0
      end if
      allocate(pv_l_f(maxval(mesh%vert(:)%n_sub_faces_neigh)))
      allocate(pv_r_f(maxval(mesh%vert(:)%n_sub_faces_neigh)))
      pos_n_clip_step = 0; pos_n_infeas_step = 0
      pos_defect_step = 0.0_DOUBLE; pos_n_sf_step = 0
    end if
    if (use_ppvp .and. .not. ppvp_nodal_average) then
      allocate(ppvp_face(maxval(mesh%vert(:)%n_sub_faces_neigh)))
    end if
    if (pi_psvp) then
      if (.not. allocated(psvp_lo_l)) then
        allocate(psvp_lo_l(maxval(mesh%vert(:)%n_sub_faces_neigh)))
        allocate(psvp_lo_r(maxval(mesh%vert(:)%n_sub_faces_neigh)))
      end if
    end if
    do i = 1, mesh%n_vert
      if (.not. vp_is_imposed(i)) vp(:, i) = 0.0_DOUBLE
    end do

    do i = 1, mesh%n_vert
      nsfn = mesh%vert(i)%n_sub_faces_neigh

      ! --- 'vp_pp_psvp' : boucle de point fixe sur lambda. A psvp_lam_iter = 0,
      ! psvp_ovr reste faux et les passes 1 et 2 sont INCHANGEES au bit pres --
      ! c'est le controle de plomberie de cette boucle.
      psvp_nit = 0
      if (pi_psvp) psvp_nit = psvp_lam_iter
      psvp_ovr = .false.
      do psvp_it = 0, psvp_nit

      ! ---- pass 1: nodal velocity system ----------------------------------
      Mp = 0.0_DOUBLE
      Rp = 0.0_DOUBLE
      Gp = 0.0_DOUBLE
      Sp = 0.0_DOUBLE
      Bp = 0.0_DOUBLE
      is_corner = .false.
      have_first_wall_norm = .false.
      first_wall_norm = 0.0_DOUBLE

      do j = 1, nsfn
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        idl = mesh%face(id_face)%left_neigh
        idr = mesh%face(id_face)%right_neigh
        area = mesh%sub_face(id_sub_face)%area

        gl = gamma_arr(idl)
        taul = sol(1, idl)
        pl = pressure(sol(:, idl), gl)
        vl = dot_product(sol(2:4, idl), mesh%face(id_face)%norm)
        if (second_order) call grp_half_step(mesh, idl, i, mesh%face(id_face)%norm, &
          dt, taul, gl, grad_v, grad_p, div_v, pl, vl)

        if (idr > 0) then
          gr = gamma_arr(idr)
          taur = sol(1, idr)
          pr = pressure(sol(:, idr), gr)
          vr = dot_product(sol(2:4, idr), mesh%face(id_face)%norm)
          if (second_order) call grp_half_step(mesh, idr, i, mesh%face(id_face)%norm, &
            dt, taur, gr, grad_v, grad_p, div_v, pr, vr)

          lam_l = max(sqrt(gl*pl*taul)/taul, &
            sqrt(max(0.0_DOUBLE, pr-pl)/taul), -(vr-vl)/taul)
          lam_r = max(sqrt(gr*pr*taur)/taur, &
            sqrt(max(0.0_DOUBLE, pl-pr)/taur), -(vr-vl)/taur)
          if (psvp_ovr) then
            lam_l = psvp_lo_l(j); lam_r = psvp_lo_r(j)
          end if

          if (.not. vp_is_imposed(i)) then
            v_bar = (lam_l*vl + lam_r*vr - (pr - pl))/(lam_l + lam_r)
            Mp = Mp + area*(lam_l + lam_r) &
              * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
            Rp = Rp + area*(lam_l + lam_r)*mesh%face(id_face)%norm*v_bar
            if (vp_gcorr) then
              ! Sp = somme A (p_l - p_r) n  : la part de pression de Rp.
              ! Gp = somme A n (x) (x_r - x_l)  : purement geometrique.
              Sp = Sp + area*(pl - pr)*mesh%face(id_face)%norm
              dxc = mesh%elem(idr)%coord - mesh%elem(idl)%coord
              Gp = Gp + area*tensor_product(mesh%face(id_face)%norm, dxc)
            end if
          end if
        else
          lam_l = sqrt(gl*pl*taul)/taul
          is_zface = (b2d .and. abs(mesh%face(id_face)%norm(3)) > 1e-8_DOUBLE)
          if (.not. vp_is_imposed(i)) then
            if (is_zface) then
              Mp = Mp + area*2.0_DOUBLE*lam_l &
                * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
            else if (idr < 0 .and. trim(bc_type(-idr)) == 'pressure') then
              Mp = Mp + area*lam_l &
                * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
              ! Etat fantome p_r = p_bc, v_r = v_l : avec Mp = A lambda, v.n = v_l + (p_l - p_bc)/lambda.
              ! p_l manquait, p_bc avait le signe inverse.
              Rp = Rp + area*(pl - bc_val(1,-idr) + lam_l*vl)*mesh%face(id_face)%norm
            else
              if (.not. have_first_wall_norm) then
                first_wall_norm = mesh%face(id_face)%norm
                have_first_wall_norm = .true.
              else
                cross_z = first_wall_norm(1)*mesh%face(id_face)%norm(2) &
                        - first_wall_norm(2)*mesh%face(id_face)%norm(1)
                if (abs(cross_z) > 1.0e-8_DOUBLE) is_corner = .true.
              end if
              Bp = Bp + area*mesh%face(id_face)%norm
              Mp = Mp + area*lam_l &
                * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
              Rp = Rp + area*(pl + lam_l*vl)*mesh%face(id_face)%norm
            end if
          end if
        end if
      end do

      if (.not. vp_is_imposed(i)) then
        if (is_corner) then
          vp(:, i) = 0.0_DOUBLE
        else if (maxval(abs(Bp)) > 1e-8_DOUBLE) then
          call inverse_3_by_3(Mp, Mp_inv)
          Pw = dot_product(Rp, matmul(Mp_inv, Bp))/dot_product(Bp, matmul(Mp_inv, Bp))
          vp(:, i) = matmul(Mp_inv, Rp - Pw*Bp)
        else
          call inverse_3_by_3(Mp, Mp_inv)
          vp(:, i) = matmul(Mp_inv, Rp)
          if (vp_gcorr) then
            ! en boundary_2d le probleme est plan : on force la 3e direction a
            ! l'identite pour que Gp soit inversible et que v_p,z reste nul.
            nd_gc = 3
            if (b2d) then
              Gp(3,:) = 0.0_DOUBLE; Gp(:,3) = 0.0_DOUBLE; Gp(3,3) = 1.0_DOUBLE
              nd_gc = 2
            end if
            detG = Gp(1,1)*(Gp(2,2)*Gp(3,3)-Gp(2,3)*Gp(3,2)) &
                 - Gp(1,2)*(Gp(2,1)*Gp(3,3)-Gp(2,3)*Gp(3,1)) &
                 + Gp(1,3)*(Gp(2,1)*Gp(3,2)-Gp(2,2)*Gp(3,1))
            gcorr_n = gcorr_n + 1
            if (abs(detG) > 1.0e-30_DOUBLE) then
              call inverse_3_by_3(Gp, Gp_inv)
              MinvG = matmul(Mp_inv, Gp)
              s_gc = 0.0_DOUBLE
              do j = 1, nd_gc
                s_gc = s_gc + MinvG(j,j)
              end do
              s_gc = s_gc / real(nd_gc, DOUBLE)
              ! on retire la reponse en pression actuelle et on remet la redressee
              vp(:, i) = matmul(Mp_inv, Rp - Sp) + s_gc*matmul(Gp_inv, Sp)
              if (b2d) vp(3, i) = 0.0_DOUBLE
            else
              ! Gp singulier (noeud de bord, mailles d un seul cote) : repli
              gcorr_nfall = gcorr_nfall + 1
            end if
          end if
        end if
      end if

      ! ---- pass 2: nodal pressure system ----------------------------------
      s_num = 0.0_DOUBLE
      s_den = 0.0_DOUBLE
      do j = 1, nsfn
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        idl = mesh%face(id_face)%left_neigh
        idr = mesh%face(id_face)%right_neigh
        area = mesh%sub_face(id_sub_face)%area

        ! z-cap exclusion, strict form (|n_z| ~ 1), as in compute_nodal_pressure_sidil
        ! after the z-cap leak fix -- not the loose |n_z|>1e-8 used elsewhere.
        if (b2d .and. abs(mesh%face(id_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE) cycle

        gl = gamma_arr(idl)
        taul = sol(1, idl)
        pl = pressure(sol(:, idl), gl)
        vl = dot_product(sol(2:4, idl), mesh%face(id_face)%norm)
        if (second_order) call grp_half_step(mesh, idl, i, mesh%face(id_face)%norm, &
          dt, taul, gl, grad_v, grad_p, div_v, pl, vl)

        if (idr > 0) then
          gr = gamma_arr(idr)
          taur = sol(1, idr)
          pr = pressure(sol(:, idr), gr)
          vr = dot_product(sol(2:4, idr), mesh%face(id_face)%norm)
          if (second_order) call grp_half_step(mesh, idr, i, mesh%face(id_face)%norm, &
            dt, taur, gr, grad_v, grad_p, div_v, pr, vr)
        else
          ! Boundary: build the MIRROR ghost state and then use the very same
          ! two-state formula, which is how both precedents in this code do it --
          ! compute_nodal_pressure_sidil here, and ns_euler_recon_module's
          ! BC_EULER_WALL on the Euler side. No special branch in the solver, and
          ! no dependence on the nodal velocity. A piston wall simply mirrors
          ! about its own velocity instead of about zero.
          v_target = 0.0_DOUBLE
          if (idr < 0) then
            if (-idr <= n_bc) then
              if (trim(bc_type(-idr)) == 'piston') then
                v_target = dot_product(bc_val(2:4, -idr), mesh%face(id_face)%norm)
              end if
            end if
          end if
          gr = gl
          taur = taul
          pr = pl
          vr = vl + 2.0_DOUBLE*(v_target - vl)
        end if

        lam_l = max(sqrt(gl*pl*taul)/taul, &
          sqrt(max(0.0_DOUBLE, pr-pl)/taul), -(vr-vl)/taul)
        lam_r = max(sqrt(gr*pr*taur)/taur, &
          sqrt(max(0.0_DOUBLE, pl-pr)/taur), -(vr-vl)/taur)
        if (psvp_ovr) then
          lam_l = psvp_lo_l(j); lam_r = psvp_lo_r(j)
        end if
        invlam = 1.0_DOUBLE/lam_l + 1.0_DOUBLE/lam_r
        p_bar = (pl/lam_l + pr/lam_r - (vr - vl))/invlam

        s_num = s_num + area*invlam*p_bar
        s_den = s_den + area*invlam
      end do

      if (s_den > 0.0_DOUBLE) then
        pp = s_num/s_den
      else
        pp = 0.0_DOUBLE
      end if


      ! --- relevement de lambda : on a maintenant v_p ET p_p de cette passe, donc
      ! A, B et C sont calculables. On monte lambda DES DEUX COTES et on refait la
      ! passe ; la conservation est intacte parce que le systeme nodal voit le meme
      ! lambda que Pi_c.
      if (psvp_it < psvp_nit) then
        do j = 1, nsfn
          id_sub_face = mesh%vert(i)%sub_face_neigh(j)
          id_face = mesh%sub_face(id_sub_face)%mesh_face
          idl = mesh%face(id_face)%left_neigh
          idr = mesh%face(id_face)%right_neigh
          v_et = dot_product(vp(:, i), mesh%face(id_face)%norm)
          gl = gamma_arr(idl); taul = sol(1, idl)
          pl = pressure(sol(:, idl), gl)
          vl = dot_product(sol(2:4, idl), mesh%face(id_face)%norm)
          eps_l = sol(5, idl) - 0.5_DOUBLE*dot_product(sol(2:4, idl), sol(2:4, idl))
          if (idr > 0) then
            gr = gamma_arr(idr); taur = sol(1, idr)
            pr = pressure(sol(:, idr), gr)
            vr = dot_product(sol(2:4, idr), mesh%face(id_face)%norm)
            eps_r = sol(5, idr) - 0.5_DOUBLE*dot_product(sol(2:4, idr), sol(2:4, idr))
          else
            gr = gl; taur = taul; pr = pl; vr = -vl; eps_r = eps_l
          end if
          lam_l = max(sqrt(max(0.0_DOUBLE, gl*pl*taul))/taul, &
            sqrt(max(0.0_DOUBLE, pr-pl)/taul), -(vr-vl)/taul)
          lam_r = max(sqrt(max(0.0_DOUBLE, gr*pr*taur))/taur, &
            sqrt(max(0.0_DOUBLE, pl-pr)/taur), -(vr-vl)/taur)
          ! n stockee sortante de idl : le cote droit voit -n, donc u_c et u_p
          ! changent de signe.
          psvp_t = psvp_lambda(eps_l, pl, pp, vl, v_et, lam_l, psvp_lam_marge)
          psvp_lo_l(j) = lam_l
          if (psvp_t > 0.0_DOUBLE) psvp_lo_l(j) = psvp_t
          psvp_t = psvp_lambda(eps_r, pr, pp, -vr, -v_et, lam_r, psvp_lam_marge)
          psvp_lo_r(j) = lam_r
          if (psvp_t > 0.0_DOUBLE) psvp_lo_r(j) = psvp_t
        end do
        psvp_ovr = .true.
      end if
      end do

      ! ---- pass 3 (schema vp_pp_ppvp): quantite nodale ppvp ----------------
      ! q = p*(u.n) est la densite de flux d'energie. On la moyenne sur le MEME
      ! modele que pp (ponderation en 1/lambda) au lieu de former le produit
      ! pp*(u_p.n) de deux quantites resolues separement.
      !
      ! Le terme de saut X doit satisfaire DEUX contraintes :
      !   - homogeneite : q/lambda est en v^2, donc X aussi ;
      !   - parite : sous retournement de la face (echange gauche/droite et n -> -n),
      !     q change de signe donc q_bar doit etre ANTI-invariant. Or delta_q est
      !     INVARIANT (les deux changements de signe se compensent), tout comme
      !     p_bar*delta_v et v_bar*delta_p. Il faut un nombre IMPAIR de facteurs
      !     anti-invariants : (u_p.n)*delta_v convient, delta_q non.
      ppvp = 0.0_DOUBLE
      if (use_ppvp) then
        if (allocated(ppvp_face)) ppvp_face = 0.0_DOUBLE
        q_num = 0.0_DOUBLE
        q_den = 0.0_DOUBLE
        do j = 1, nsfn
          id_sub_face = mesh%vert(i)%sub_face_neigh(j)
          id_face = mesh%sub_face(id_sub_face)%mesh_face
          idl = mesh%face(id_face)%left_neigh
          idr = mesh%face(id_face)%right_neigh
          area = mesh%sub_face(id_sub_face)%area

          if (b2d .and. abs(mesh%face(id_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE) cycle

          gl = gamma_arr(idl)
          taul = sol(1, idl)
          pl = pressure(sol(:, idl), gl)
          vl = dot_product(sol(2:4, idl), mesh%face(id_face)%norm)
          if (second_order) call grp_half_step(mesh, idl, i, mesh%face(id_face)%norm, &
            dt, taul, gl, grad_v, grad_p, div_v, pl, vl)

          if (idr > 0) then
            gr = gamma_arr(idr)
            taur = sol(1, idr)
            pr = pressure(sol(:, idr), gr)
            vr = dot_product(sol(2:4, idr), mesh%face(id_face)%norm)
            if (second_order) call grp_half_step(mesh, idr, i, mesh%face(id_face)%norm, &
              dt, taur, gr, grad_v, grad_p, div_v, pr, vr)
          else
            ! meme etat miroir que la passe de pression
            v_target = 0.0_DOUBLE
            if (idr < 0) then
              if (-idr <= n_bc) then
                if (trim(bc_type(-idr)) == 'piston') then
                  v_target = dot_product(bc_val(2:4, -idr), mesh%face(id_face)%norm)
                end if
              end if
            end if
            gr = gl
            taur = taul
            pr = pl
            vr = vl + 2.0_DOUBLE*(v_target - vl)
          end if

          lam_l = max(sqrt(gl*pl*taul)/taul, &
            sqrt(max(0.0_DOUBLE, pr-pl)/taul), -(vr-vl)/taul)
          lam_r = max(sqrt(gr*pr*taur)/taur, &
            sqrt(max(0.0_DOUBLE, pl-pr)/taur), -(vr-vl)/taur)
          invlam = 1.0_DOUBLE/lam_l + 1.0_DOUBLE/lam_r

          ! q = p * (u_p.n) : le produit porte la vitesse NODALE, pas celle des
          ! mailles. Avec les vitesses de maille, q_l = q_r = 0 a t=0 (fluide au
          ! repos) et ppvp serait identiquement nul quelle que soit la ponderation
          ! -- le souffle ne fournirait aucun travail et le schema ne demarrerait
          ! pas. C'est aussi ce que dit le nom : ppvp = pp * vp.
          v_node = dot_product(vp(:, i), mesh%face(id_face)%norm)
          ql = pl*v_node
          qr = pr*v_node
          if (trim(ppvp_jump_mode) == 'dp') then
            ! variante : le saut porte sur p (donc sur q a v_node fixe), rendu
            ! homogene par les impedances.
            xjump = v_node*(pr - pl)*invlam*0.5_DOUBLE
          else
            xjump = v_node*(vr - vl)
          end if
          q_bar = (ql/lam_l + qr/lam_r - xjump)/invlam

          if (ppvp_nodal_average) then
            q_num = q_num + area*invlam*q_bar
            q_den = q_den + area*invlam
          else
            ! variante par sous-face : on stocke q_bar la ou le flux le relira
            ppvp_face(j) = q_bar
          end if
        end do
        if (ppvp_nodal_average .and. q_den > 0.0_DOUBLE) ppvp = q_num/q_den
      end if

      ! ---- passe 'vp_pp_pos' : ecretage du flux d'energie puis redistribution
      ! Relations de saut -> energie interne des etats intermediaires :
      !   eps*_s = eps_s - (sgn_s/lam_s)[(pv)*_s - p_bar v_n_s] - (p_bar-p_s)^2/(2 lam_s^2)
      ! eps*_s >= 0 donne, avec B_s = lam_s eps_s - (p_bar-p_s)^2/(2 lam_s) :
      !   (pv)*_l <= p_bar v_n_l + B_l        (plafond)
      !   (pv)*_r >= p_bar v_n_r - B_r        (plancher)
      ! L'ecretage rend le flux BIVALUE, donc non conservatif ; on redistribue le
      ! surplus sur le noeud, pondere par le mou disponible de chaque sous-face.
      ! Ponderer par le mou n'est pas un choix : reduire la ou il n'y en a pas
      ! ferait repasser sous la borne.
      if (pvp_pos) then
        S_sur = 0.0_DOUBLE
        R_room = 0.0_DOUBLE
        n_clip_step = 0
        do j = 1, nsfn
          id_sub_face = mesh%vert(i)%sub_face_neigh(j)
          id_face = mesh%sub_face(id_sub_face)%mesh_face
          idl = mesh%face(id_face)%left_neigh
          idr = mesh%face(id_face)%right_neigh
          area = mesh%sub_face(id_sub_face)%area
          v_et = dot_product(vp(:, i), mesh%face(id_face)%norm)

          gl = gamma_arr(idl)
          taul = sol(1, idl)
          pl = pressure(sol(:, idl), gl)
          vl = dot_product(sol(2:4, idl), mesh%face(id_face)%norm)
          eps_l = sol(5, idl) - 0.5_DOUBLE*dot_product(sol(2:4, idl), sol(2:4, idl))
          if (idr > 0) then
            gr = gamma_arr(idr)
            taur = sol(1, idr)
            pr = pressure(sol(:, idr), gr)
            vr = dot_product(sol(2:4, idr), mesh%face(id_face)%norm)
            eps_r = sol(5, idr) - 0.5_DOUBLE*dot_product(sol(2:4, idr), sol(2:4, idr))
            E_r_s = sol(5, idr)
          else
            gr = gl; taur = taul; pr = pl; vr = -vl; eps_r = eps_l
            E_r_s = sol(5, idl)
          end if
          E_l_s = sol(5, idl)
          lam_l = max(sqrt(max(0.0_DOUBLE, gl*pl*taul))/taul, &
            sqrt(max(0.0_DOUBLE, pr-pl)/taul), -(vr-vl)/taul)
          lam_r = max(sqrt(max(0.0_DOUBLE, gr*pr*taur))/taur, &
            sqrt(max(0.0_DOUBLE, pl-pr)/taur), -(vr-vl)/taur)
          if (lam_l <= 0.0_DOUBLE .or. lam_r <= 0.0_DOUBLE) then
            pv_l_f(j) = pp*v_et; pv_r_f(j) = pp*v_et; cycle
          end if

          if (pvp_pos_sv) then
            ! vl et vr sont tous deux projetes sur mesh%face(id_face)%norm, qui est
            ! la normale sortante de idl : c'est exactement la convention de la
            ! construction (un seul n, celui de c). Sur une sous-face de BORD, l'etat
            ! miroir donne eps_r = eps_l et vr = -vl, donc Phi_bas = -Phi_haut :
            ! l'intervalle est symetrique et le flux naturel d'une paroi, ou
            ! v_p.n = 0, en est le milieu exact. C'est ce que sv_n_bnd_corr verifie
            ! au lieu de le supposer.
            call pos_sv_bounds(lam_l, eps_l, vl, lam_r, eps_r, vr, pp, &
                               sv_hi, sv_lo, sv_d2)
            call pos_sv_flux(pp*v_et, sv_hi, sv_lo, sv_d2, pvp_pos_sv_always, &
                             sv_fl, sv_corr, sv_infeas, sv_isnan)
            pv_l_f(j) = sv_fl; pv_r_f(j) = sv_fl
            sv_n_sf = sv_n_sf + 1
            if (sv_infeas) then
              sv_n_infeas = sv_n_infeas + 1
              sv_st_infeas = sv_st_infeas + 1
            end if
            if (sv_d2 == sv_d2) sv_st_d2min = min(sv_st_d2min, sv_d2)
            if (sv_corr) sv_n_corr = sv_n_corr + 1
            ! Un flux naturel NaN echoue le test d'appartenance a l'intervalle, donc
            ! il est 'corrige' et remplace par le milieu, qui est fini. Le limiteur
            ! BLANCHIT alors le NaN au lieu de le laisser voir. C'est a compter
            ! separement : une correction de positivite et un NaN maquille ne sont
            ! pas le meme evenement, et seul le second est une panne masquee.
            if (sv_isnan) sv_n_corr_nan = sv_n_corr_nan + 1
            ! Test du mecanisme annonce : la correction ne mordrait que sur les
            ! sous-faces separant deux mailles FROIDES, ou la borne vaut
            ! z_froid eps_froid et ou les quantites nodales, fixees par la maille
            ! chaude, injectent un flux d'ordre 1. On classe donc chaque correction
            ! par la pression de ses DEUX cotes au lieu de le supposer.
            if (sv_corr) then
              sv_step_n = sv_step_n + 1
              sv_ray = norm2(mesh%vert(mesh%sub_face(id_sub_face)%mesh_vert)%coord(1:2))
              sv_step_r = sv_step_r + sv_ray
              sv_step_rmin = min(sv_step_rmin, sv_ray)
              sv_step_rmax = max(sv_step_rmax, sv_ray)
              sv_pmax_corr = max(sv_pmax_corr, max(pl, pr))
              sv_pmin_corr = min(sv_pmin_corr, max(pl, pr))
              if (max(pl, pr) < 1.0e-6_DOUBLE) then
                sv_n_corr_froid = sv_n_corr_froid + 1
              else if (min(pl, pr) > 1.0e-6_DOUBLE) then
                sv_n_corr_chaud = sv_n_corr_chaud + 1
              else
                sv_n_corr_mixte = sv_n_corr_mixte + 1
              end if
            end if
            if (idr <= 0) then
              ! Sur un maillage 2D par extrusion, les deux faces z sont des faces
              ! de BORD au sens de idr <= 0 alors que ce sont des plans miroir et
              ! pas des parois. Les compter ensemble donnerait ici 659721 bords sur
              ! 1310817 sous-faces, soit la moitie du maillage, et rendrait le compte
              ! de bords corriges illisible. On les separe.
              if (b2d .and. abs(mesh%face(id_face)%norm(3)) > 1.0e-8_DOUBLE) then
                sv_n_bnd_z = sv_n_bnd_z + 1
                if (sv_corr) sv_n_bnd_z_corr = sv_n_bnd_z_corr + 1
              else
                sv_n_bnd = sv_n_bnd + 1
                if (sv_corr) sv_n_bnd_corr = sv_n_bnd_corr + 1
              end if
            end if
            ! pp*v_et peut etre NaN : vp_pp en produit sur sedov_tri. On les COMPTE
            ! au lieu de les sommer, sinon un seul NaN rend la moyenne inutilisable
            ! -- et le max, lui, les ignore et n'en dit rien.
            if (pp*v_et /= pp*v_et) then
              sv_n_nan = sv_n_nan + 1
            else if (sv_d2 > 0.0_DOUBLE) then
              sv_d2_min = min(sv_d2_min, sv_d2)
              sv_ecart = abs(pp*v_et - 0.5_DOUBLE*(sv_hi + sv_lo))/sv_d2
              ! POPULATION BIMODALE : une sous-face admissible a par definition un
              ! ecart <= 0.5, une sous-face corrigee peut etre a 1e14 largeurs. Une
              ! moyenne sur les deux ne mesure que la seconde diluee par la premiere
              ! -- la premiere version de ce diagnostic sortait 1.23e9 de moyenne
              ! alors que 99.98 % de la population est sous 0.5. On separe.
              if (sv_corr) then
                sv_ec_out_max = max(sv_ec_out_max, sv_ecart)
                sv_ec_out_tot = sv_ec_out_tot + sv_ecart
                sv_n_out = sv_n_out + 1
              else
                sv_ec_in_max = max(sv_ec_in_max, sv_ecart)
                sv_ec_in_tot = sv_ec_in_tot + sv_ecart
                sv_n_in = sv_n_in + 1
              end if
            end if
            cycle
          end if

          ! valeurs 1D de la sous-face, eq. (68) : meme solveur acoustique que
          ! celui qui alimente les systemes nodaux.
          p_bar_f = (lam_r*pl + lam_l*pr - lam_l*lam_r*(vr-vl))/(lam_l+lam_r)
          v_bar_f = (lam_l*vl + lam_r*vr - (pr-pl))/(lam_l+lam_r)
          select case (pos_base_mode)
          case ('V1'); t_f = p_bar_f*v_et
          case ('V2'); t_f = p_bar_f*v_bar_f
          case ('V3'); t_f = pp*v_bar_f
          case ('V4')
            ! E miroir = E de la maille (mirror_state conserve sol(5)), donc le
            ! saut est nul sur une sous-face de bord.
            t_f = (pl*vl/lam_l + pr*vr/lam_r - (E_r_s - E_l_s)) &
                 /(1.0_DOUBLE/lam_l + 1.0_DOUBLE/lam_r)
          case ('V5')
            ! Le flux de Godunov 1D de V4, mais differenciant l'energie INTERNE
            ! au lieu de l'energie TOTALE. Rien d'autre ne change : memes poids
            ! 1/lam_l et 1/lam_r que le probleme de Riemann de sous-face, aucun
            ! coefficient libre, aucun ecretage. C'est strictement moins
            ! dissipatif que V4 -- la part CINETIQUE du saut n'est plus dissipee --
            ! et ca se derive au lieu de se regler : W1 fait la meme substitution
            ! mais avec une moyenne ARITHMETIQUE (lam_l+lam_r)/2, qui elle n'est
            ! pas celle du solveur et est donc le choix arbitraire de cette forme.
            t_f = (pl*vl/lam_l + pr*vr/lam_r - (eps_r - eps_l)) &
                 /(1.0_DOUBLE/lam_l + 1.0_DOUBLE/lam_r)
          case ('W1')
            t_f = 0.5_DOUBLE*(pl*vl + pr*vr) &
                - 0.5_DOUBLE*0.5_DOUBLE*(lam_l+lam_r)*(eps_r - eps_l)
          case ('W2')
            ! coefficient reglable sur la forme QUI GAGNE, c'est-a-dire le saut
            ! d'energie interne de W1, et non le saut d'energie totale de V4.
            t_f = 0.5_DOUBLE*(pl*vl + pr*vr) &
                - pos_w_alpha*0.5_DOUBLE*0.5_DOUBLE*(lam_l+lam_r)*(eps_r - eps_l)
          case ('W3')
            t_f = 0.5_DOUBLE*(pl*vl + pr*vr)
          case ('W6')
            ! controle complementaire de W1 : E = eps + |v|^2/2 donne la
            ! decomposition EXACTE  (E_r-E_l) = (eps_r-eps_l) + (ek_r-ek_l).
            ! W1 ne dissipe que la premiere part, W6 que la seconde, sur la meme
            ! base centree -- les deux termes de jump s'additionnent donc en
            ! celui de la forme totale. ek est pris comme E - eps et non recalcule
            ! depuis sol(2:4), pour que l'identite soit vraie au bit pres, et il
            ! s'annule sur une sous-face de bord comme (E_r - E_l).
            t_f = 0.5_DOUBLE*(pl*vl + pr*vr) &
                - 0.5_DOUBLE*0.5_DOUBLE*(lam_l+lam_r) &
                  *((E_r_s - eps_r) - (E_l_s - eps_l))
          case ('W5')
            ! combinaison de W1 et W4 : les deux modifications portent sur des
            ! termes DISTINCTS du flux -- W1 change ce qu'on differencie (eps au
            ! lieu de E, donc la part cinetique du saut n'est plus dissipee), W4
            ! change le facteur advectif du terme centre (vitesse nodale, donc le
            ! couplage multidimensionnel est conserve). Rien ne les oppose
            ! structurellement.
            t_f = 0.5_DOUBLE*(pl + pr)*v_et &
                - pos_w_alpha*0.5_DOUBLE*0.5_DOUBLE*(lam_l+lam_r)*(eps_r - eps_l)
          case ('W4')
            t_f = 0.5_DOUBLE*(pl + pr)*v_et &
                - 0.5_DOUBLE*0.5_DOUBLE*(lam_l+lam_r)*(E_r_s - E_l_s)
          case ('W7')
            ! forme COMPATIBLE au sens de Caramana (1998) -- meme pression p_p que
            ! dans le flux de quantite de mouvement, donc V0 -- AUGMENTEE du saut
            ! d'energie interne, celui que W1/W6 ont designe comme le seul terme
            ! dissipatif utile. On n'y remplace pas la base, on lui ajoute.
            t_f = pp*v_et &
                - 0.5_DOUBLE*0.5_DOUBLE*(lam_l+lam_r)*(eps_r - eps_l)
          case ('W8')
            ! idem W7 mais avec le saut TOTAL : mesure l'ecart interne/total sur
            ! la base compatible, comme W1 contre V4 le mesure sur la base centree.
            t_f = pp*v_et &
                - 0.5_DOUBLE*0.5_DOUBLE*(lam_l+lam_r)*(E_r_s - E_l_s)
          case ('W9')
            ! dissipation en SAUT DE VITESSE et non en saut d'energie : elle est
            ! NULLE a t=0 partout, donc elle ne touche pas au depot de Sedov, au
            ! contraire du saut d'energie interne qui y est maximal.
            ! ATTENTION, defaut connu et mesure : ce terme est PAIR sous le
            ! reetiquetage (gauche<->droite, n -> -n), alors que la base pp*v_et
            ! est IMPAIRE. t_f n'est donc pas independant de l'orientation
            ! stockee de la face. Le schema reste conservatif -- une seule valeur
            ! t_f alimente les deux mailles -- mais il depend de l'orientation du
            ! maillage et non de la physique. W9B corrige exactement ce point.
            t_f = pp*v_et &
                - 0.5_DOUBLE*0.5_DOUBLE*(lam_l+lam_r)*(vr - vl)**2
          case ('WB')
            ! W9B : meme idee, parite rendue IMPAIRE en multipliant le saut de
            ! vitesse par la vitesse nodale normale au lieu de l'elever au carre.
            ! Sous le reetiquetage, v_et change de signe et (vr-vl) est invariant,
            ! donc le produit change de signe comme la base. Toujours nul a t=0.
            t_f = pp*v_et &
                - 0.5_DOUBLE*0.5_DOUBLE*(lam_l+lam_r)*v_et*(vr - vl)
          case ('WA')
            ! W10 : forme UNILATERALE a la maniere de vp, un carre par cote.
            ! Elle est traitee apres le select car elle donne deux valeurs
            ! distinctes ; t_f ne sert ici que de base.
            t_f = pp*v_et
          case ('WD')
            ! W11 : forme de MANUEL -- la dissipation entre comme une PRESSION
            ! VISQUEUSE q, scalaire positif, et non comme un terme ajoute au flux.
            ! q est PAIR sous le reetiquetage et v_et IMPAIR, donc le produit a la
            ! bonne parite, et q >= 0 par construction. C'est la seule facon de
            ! concilier "positif" et "impair", que W9 (paire) et W9B (de signe
            ! quelconque) ratent chacune par un bout.
            ! RESERVE a lire avec le resultat : q >= 0 garantit la condition 1 de
            ! Caramana dans la construction COMPATIBLE, ou q entre aussi dans le
            ! flux de quantite de mouvement. Ici il n'entre que dans l'energie, le
            ! travail et la force ne derivent donc plus du meme objet, et la
            ! positivite de la production n'est plus garantie par q >= 0 seul.
            ! Elle doit etre MESUREE.
            q_visc = 0.5_DOUBLE*0.5_DOUBLE*(lam_l+lam_r)*max(0.0_DOUBLE, -(vr - vl))
            t_f = (pp + q_visc)*v_et
          case ('WE')
            ! W12 : viscosite TENSORIELLE de type Navier-Stokes, premiere des deux
            ! solutions de Caramana au probleme de la compression uniforme.
            ! Trois des quatre conditions sont satisfaites PAR CONSTRUCTION :
            !   compression uniforme : grad_v = a*I -> D = a*I -> Dp = 0
            !   rotation rigide      : grad_v antisymetrique -> D = 0
            !   compression nulle    : v uniforme -> grad_v = 0 -> tau = 0
            ! La condition 1 (dissipativite) ne l'est PAS : tau n'entre que dans
            ! le flux d'energie et pas dans celui de quantite de mouvement, donc
            ! la production n'a pas de contrepartie et son signe n'est pas acquis.
            G_t = 0.5_DOUBLE*(grad_v(:, :, idl) + grad_v(:, :, max(idr, idl)))
            D_t = 0.5_DOUBLE*(G_t + transpose(G_t))
            if (b2d) then
              ! maillage 2D par extrusion : la trace deviatorique porte sur x,y
              ! seulement, et la direction z est entierement retiree du tenseur.
              D_t(3, :) = 0.0_DOUBLE; D_t(:, 3) = 0.0_DOUBLE
              tr_d = 0.5_DOUBLE*(D_t(1, 1) + D_t(2, 2))
              D_t(1, 1) = D_t(1, 1) - tr_d
              D_t(2, 2) = D_t(2, 2) - tr_d
            else
              tr_d = (D_t(1, 1) + D_t(2, 2) + D_t(3, 3))/3.0_DOUBLE
              D_t(1, 1) = D_t(1, 1) - tr_d
              D_t(2, 2) = D_t(2, 2) - tr_d
              D_t(3, 3) = D_t(3, 3) - tr_d
            end if
            ! longueur locale : distance entre centres de mailles, definie
            ! independamment de l'extrusion, contrairement a sqrt(aire).
            if (idr > 0) then
              h_f = norm2(mesh%elem(idr)%coord - mesh%elem(idl)%coord)
            else
              h_f = 2.0_DOUBLE*abs(dot_product(mesh%vert(i)%coord &
                    - mesh%elem(idl)%coord, mesh%face(id_face)%norm))
            end if
            mu_v = pos_w_alpha*0.5_DOUBLE*(lam_l + lam_r)*h_f
            tau_n = 2.0_DOUBLE*mu_v*matmul(D_t, mesh%face(id_face)%norm)
            t_f = pp*v_et + dot_product(tau_n, vp(:, i))
          case default; t_f = pp*v_et
          end select

          if (pos_base_mode == 'WA' .or. pos_base_mode == 'WC') then
            ! NON CONSERVATIF PAR CONSTRUCTION : pv_r_f - pv_l_f est une somme de
            ! deux carres, donc de signe constant. Le defaut est exactement celui
            ! que S_sur/pos_defect_work mesure deja pour l'ecretage, on le laisse
            ! donc s'accumuler et se lire dans le diagnostic.
            ! delta_c = (u_p - v_c).n_c avec n_c la normale SORTANTE de la maille c,
            ! donc delta_l = v_et - vl pour la gauche et delta_r = vr - v_et pour la
            ! droite : les deux sont negatifs en compression. WC (= W10b) n'allume
            ! le terme que la, condition 4 de Caramana -- nul en detente.
            d_l = v_et - vl
            d_r = vr - v_et
            if (pos_base_mode == 'WC') then
              if (d_l >= 0.0_DOUBLE) d_l = 0.0_DOUBLE
              if (d_r >= 0.0_DOUBLE) d_r = 0.0_DOUBLE
            end if
            pv_l_f(j) = t_f - 0.5_DOUBLE*lam_l*d_l**2
            pv_r_f(j) = t_f + 0.5_DOUBLE*lam_r*d_r**2
            pos_defect_step = pos_defect_step + area*(pv_r_f(j) - pv_l_f(j))
            cycle
          end if

          if (pos_no_limit) then
            pv_l_f(j) = t_f; pv_r_f(j) = t_f
            cycle
          end if

          B_l = lam_l*eps_l - (pp-pl)**2/(2.0_DOUBLE*lam_l)
          B_r = lam_r*eps_r - (pp-pr)**2/(2.0_DOUBLE*lam_r)
          cap_l = pp*vl + B_l
          flo_r = pp*vr - B_r

          pv_l_f(j) = min(t_f, cap_l)
          if (idr > 0) then
            pv_r_f(j) = max(t_f, flo_r)
          else
            ! sous-face de bord : un seul cote, une seule borne, et aucun
            ! partenaire pour la conservation -> hors redistribution, comptee.
            pv_r_f(j) = pv_l_f(j)
            if (pv_l_f(j) < t_f) pos_n_clip_bnd = pos_n_clip_bnd + 1
            cycle
          end if
          if (pv_l_f(j) < t_f .or. pv_r_f(j) > t_f) then
            n_clip_step = n_clip_step + 1
            if (allocated(pos_clip_cell)) then
              pos_clip_cell(idl) = pos_clip_cell(idl) + 1
              pos_clip_cell(idr) = pos_clip_cell(idr) + 1
            end if
          end if

          S_sur = S_sur + area*(pv_r_f(j) - pv_l_f(j))
          av_l = cap_l - pv_l_f(j)
          av_r = pv_r_f(j) - flo_r
          R_room = R_room + area*(av_l + av_r)
        end do

        if (S_sur > 0.0_DOUBLE) then
          if (S_sur > R_room) then
            ! INFAISABLE : aucun choix ne satisfait bornes ET conservation. On fait
            ! au mieux (tout le mou) et on COMPTE le defaut, plutot que de choisir
            ! en silence entre monter lambda et lacher la conservation.
            pos_n_infeas_step = pos_n_infeas_step + 1
            pos_defect_step = pos_defect_step + (S_sur - R_room)
          end if
          do j = 1, nsfn
            id_sub_face = mesh%vert(i)%sub_face_neigh(j)
            id_face = mesh%sub_face(id_sub_face)%mesh_face
            if (mesh%face(id_face)%right_neigh <= 0) cycle
            idl = mesh%face(id_face)%left_neigh
            idr = mesh%face(id_face)%right_neigh
            v_et = dot_product(vp(:, i), mesh%face(id_face)%norm)
            gl = gamma_arr(idl); taul = sol(1, idl)
            pl = pressure(sol(:, idl), gl)
            vl = dot_product(sol(2:4, idl), mesh%face(id_face)%norm)
            eps_l = sol(5, idl) - 0.5_DOUBLE*dot_product(sol(2:4, idl), sol(2:4, idl))
            gr = gamma_arr(idr); taur = sol(1, idr)
            pr = pressure(sol(:, idr), gr)
            vr = dot_product(sol(2:4, idr), mesh%face(id_face)%norm)
            eps_r = sol(5, idr) - 0.5_DOUBLE*dot_product(sol(2:4, idr), sol(2:4, idr))
            lam_l = max(sqrt(max(0.0_DOUBLE, gl*pl*taul))/taul, &
              sqrt(max(0.0_DOUBLE, pr-pl)/taul), -(vr-vl)/taul)
            lam_r = max(sqrt(max(0.0_DOUBLE, gr*pr*taur))/taur, &
              sqrt(max(0.0_DOUBLE, pl-pr)/taur), -(vr-vl)/taur)
            if (lam_l <= 0.0_DOUBLE .or. lam_r <= 0.0_DOUBLE) cycle
            if (pos_no_limit) then
            pv_l_f(j) = t_f; pv_r_f(j) = t_f
            cycle
          end if

          B_l = lam_l*eps_l - (pp-pl)**2/(2.0_DOUBLE*lam_l)
            B_r = lam_r*eps_r - (pp-pr)**2/(2.0_DOUBLE*lam_r)
            av_l = (pp*vl + B_l) - pv_l_f(j)
            av_r = pv_r_f(j) - (pp*vr - B_r)
            room = av_l + av_r
            if (room <= 0.0_DOUBLE) cycle
            if (S_sur > R_room) then
              reduc = room
            else
              reduc = S_sur*room/R_room
            end if
            pv_l_f(j) = pv_l_f(j) + reduc*av_l/room
            pv_r_f(j) = pv_r_f(j) - reduc*av_r/room
          end do
        end if
        pos_n_clip_step = pos_n_clip_step + n_clip_step
        pos_n_sf_step = pos_n_sf_step + nsfn
        if (n_clip_step > pos_n_clip_max) pos_n_clip_max = n_clip_step
      end if

      ! ---- flux assembly: one nodal pressure, one nodal velocity ----------
      do j = 1, nsfn
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        idl = mesh%face(id_face)%left_neigh
        idr = mesh%face(id_face)%right_neigh

        v_et = dot_product(vp(:, i), mesh%face(id_face)%norm)

        if (diag_entropy) then
          ! contribution de cette sous-face aux deux mailles. La normale sortante
          ! est +n pour idl et -n pour idr, d'ou le changement de signe.
          call add_entropy_diag(mesh, idl, id_face, &
            mesh%sub_face(id_sub_face)%area, 1.0_DOUBLE, vp(:, i), sol, gamma_arr, pp)
          if (idr > 0) call add_entropy_diag(mesh, idr, id_face, &
            mesh%sub_face(id_sub_face)%area, -1.0_DOUBLE, vp(:, i), sol, gamma_arr, pp)
        end if

        flux(1)   = -v_et
        flux(2:4) = pp*mesh%face(id_face)%norm
        if (use_ppvp) then
          if (ppvp_nodal_average) then
            flux(5) = ppvp
          else
            flux(5) = ppvp_face(j)
          end if
        else
          flux(5) = pp*v_et
        end if
        ! 'vp_pp_pos' : le flux d'energie est BIVALUE apres ecretage.
        flux5_l = flux(5)
        flux5_r = flux(5)
        if (pvp_pos) then
          flux5_l = pv_l_f(j)
          flux5_r = pv_r_f(j)
        end if
        ! 'vp_pp_psvp' : l'energie prend la pression de sous-face UNILATERALE de
        ! classic. Le flux est donc bivalue, et la conservation vient du bilan
        ! nodal somme_c A Pi_c n_c = 0 et non de l'univocite -- voir pi_psvp.
        if (pi_psvp) then
          gl = gamma_arr(idl)
          taul = sol(1, idl)
          pl = pressure(sol(:, idl), gl)
          vl = dot_product(sol(2:4, idl), mesh%face(id_face)%norm)
          if (second_order) call grp_half_step(mesh, idl, i, mesh%face(id_face)%norm, &
            dt, taul, gl, grad_v, grad_p, div_v, pl, vl)
          if (idr > 0) then
            gr = gamma_arr(idr)
            taur = sol(1, idr)
            pr = pressure(sol(:, idr), gr)
            vr = dot_product(sol(2:4, idr), mesh%face(id_face)%norm)
            if (second_order) call grp_half_step(mesh, idr, i, mesh%face(id_face)%norm, &
              dt, taur, gr, grad_v, grad_p, div_v, pr, vr)
          else
            ! meme etat miroir que la passe nodale : un seul cote, pas de voisin.
            gr = gl; taur = taul; pr = pl; vr = -vl
          end if
          eps_l = sol(5, idl) - 0.5_DOUBLE*dot_product(sol(2:4, idl), sol(2:4, idl))
          if (idr > 0) then
            eps_r = sol(5, idr) - 0.5_DOUBLE*dot_product(sol(2:4, idr), sol(2:4, idr))
          else
            eps_r = eps_l
          end if
          lam_l = max(sqrt(max(0.0_DOUBLE, gl*pl*taul))/taul, &
            sqrt(max(0.0_DOUBLE, pr-pl)/taul), -(vr-vl)/taul)
          lam_r = max(sqrt(max(0.0_DOUBLE, gr*pr*taur))/taur, &
            sqrt(max(0.0_DOUBLE, pl-pr)/taur), -(vr-vl)/taur)
          ! Le flux DOIT voir le meme lambda que le systeme nodal, sinon
          ! somme_c A Pi_c n_c cesse d'etre nulle et la conservation de l'energie
          ! -- la raison d'etre de ce schema -- est perdue.
          if (psvp_ovr) then
            lam_l = psvp_lo_l(j); lam_r = psvp_lo_r(j)
          end if
          ! La normale stockee est SORTANTE de idl, donc n_r = -n : le signe du
          ! terme en lambda s'inverse pour le cote droit.
          psvp_pil = pl - lam_l*(v_et - vl)
          psvp_pir = pr + lam_r*(v_et - vr)
          ! --- positivite effectivement atteinte : on COMPTE, on ne proclame pas.
          ! eps*_c = eps_c - (1/lam)[(pv)* - p_p u_c] - (p_p - p_c)^2/(2 lam^2)
          psvp_epss = eps_l - (psvp_pil*v_et - pp*vl)/lam_l &
                    - (pp - pl)**2/(2.0_DOUBLE*lam_l**2)
          if (psvp_epss < 0.0_DOUBLE) psvp_n_epsneg = psvp_n_epsneg + 1
          if (eps_l + (v_et - vl)*v_et <= 0.0_DOUBLE) psvp_n_infeas = psvp_n_infeas + 1
          if (idr > 0) then
            psvp_epss = eps_r - (psvp_pir*(-v_et) - pp*(-vr))/lam_r &
                      - (pp - pr)**2/(2.0_DOUBLE*lam_r**2)
            if (psvp_epss < 0.0_DOUBLE) psvp_n_epsneg = psvp_n_epsneg + 1
            if (eps_r + (-v_et + vr)*(-v_et) <= 0.0_DOUBLE) psvp_n_infeas = psvp_n_infeas + 1
          end if
          flux5_l = psvp_pil*v_et
          flux5_r = psvp_pir*v_et
          ! --- diagnostic : le terme qui separe ce bilan de celui de classic.
          area = mesh%sub_face(id_sub_face)%area
          psvp_t = area*vl*(pp - psvp_pil)
          psvp_n_sf = psvp_n_sf + 1
          psvp_sum = psvp_sum + psvp_t
          if (psvp_t < 0.0_DOUBLE) then
            psvp_n_neg = psvp_n_neg + 1
            psvp_sum_neg = psvp_sum_neg + psvp_t
            if (psvp_t < psvp_worst) psvp_worst = psvp_t
          end if
          psvp_ref = psvp_ref + area*lam_l*(v_et - vl)**2
          if (idr > 0) then
            psvp_t = area*(-vr)*(pp - psvp_pir)
            psvp_n_sf = psvp_n_sf + 1
            psvp_sum = psvp_sum + psvp_t
            if (psvp_t < 0.0_DOUBLE) then
              psvp_n_neg = psvp_n_neg + 1
              psvp_sum_neg = psvp_sum_neg + psvp_t
              if (psvp_t < psvp_worst) psvp_worst = psvp_t
            end if
            psvp_ref = psvp_ref + area*lam_r*(v_et - vr)**2
          end if
        end if

        do var = 1, 4
          rhs(var, idl) = rhs(var, idl) &
            - mesh%sub_face(id_sub_face)%area/mass(idl)*flux(var)
        end do
        rhs(5, idl) = rhs(5, idl) &
          - mesh%sub_face(id_sub_face)%area/mass(idl)*flux5_l
        if (idr > 0) then
          do var = 1, 4
            rhs(var, idr) = rhs(var, idr) &
              + mesh%sub_face(id_sub_face)%area/mass(idr)*flux(var)
          end do
          rhs(5, idr) = rhs(5, idr) &
            + mesh%sub_face(id_sub_face)%area/mass(idr)*flux5_r
        end if
      end do
    end do

    if (pvp_pos) then
      pos_n_clip = pos_n_clip + pos_n_clip_step
      pos_n_sf = pos_n_sf + pos_n_sf_step
      pos_n_infeas = pos_n_infeas + pos_n_infeas_step
      if (pos_n_infeas_step > 0) pos_n_steps_infeas = pos_n_steps_infeas + 1
    end if
  end subroutine compute_rhs_lagrange_vp_pp

  ! Contribution d'une sous-face a la production d'entropie de la maille ic.
  ! sgn = +1 si la normale stockee est sortante pour ic, -1 sinon.
  subroutine add_entropy_diag(mesh, ic, id_face, area, sgn, u_p, sol, gamma_arr, pp)
    implicit none
    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: ic, id_face
    real(kind=DOUBLE), intent(in) :: area, sgn, pp
    real(kind=DOUBLE), dimension(3), intent(in) :: u_p
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    real(kind=DOUBLE), dimension(3) :: nout
    real(kind=DOUBLE) :: pc, tauc, lamc, du

    nout = sgn*mesh%face(id_face)%norm
    pc   = pressure(sol(:, ic), gamma_arr(ic))
    tauc = sol(1, ic)
    lamc = sqrt(max(0.0_DOUBLE, gamma_arr(ic)*pc*tauc))/tauc
    du   = dot_product(u_p - sol(2:4, ic), nout)
    entropy_prod(ic) = entropy_prod(ic) + area*(pc - pp)*du
    entropy_ref(ic)  = entropy_ref(ic)  + area*lamc*du*du
  end subroutine add_entropy_diag

  ! Resume du diagnostic, appele depuis la boucle en temps.
  subroutine print_entropy_diag(mesh, t, me, num_procs)
    implicit none
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), intent(in) :: t
    integer(kind=ENTIER), intent(in) :: me, num_procs
    integer(kind=ENTIER) :: i, nneg, ntot, ierr
    real(kind=DOUBLE) :: sprod, sref, worst, r_worst, rr

    if (.not. allocated(entropy_prod)) return
    nneg = 0; ntot = 0; sprod = 0.0_DOUBLE; sref = 0.0_DOUBLE
    worst = 0.0_DOUBLE; r_worst = -1.0_DOUBLE
    do i = 1, mesh%n_elems
      if (mesh%elem(i)%is_ghost) cycle
      ntot = ntot + 1
      sprod = sprod + entropy_prod(i)
      sref  = sref  + entropy_ref(i)
      if (entropy_prod(i) < 0.0_DOUBLE) then
        nneg = nneg + 1
        if (entropy_prod(i) < worst) then
          worst = entropy_prod(i)
          rr = norm2(mesh%elem(i)%coord(:2))
          r_worst = rr
        end if
      end if
    end do
    if (num_procs > 1) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, nneg, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, ntot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, sprod, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, sref, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, worst, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    end if
    if (me == 0) then
      print '(a,es12.5,a,i7,a,i7,a,f7.3,a,es12.5,a,es12.5,a,f7.4)', &
        "[entropie] t=", t, "  mailles<0: ", nneg, " /", ntot, " = ", &
        100.0_DOUBLE*real(nneg, DOUBLE)/max(real(ntot, DOUBLE), 1.0_DOUBLE), &
        " %  somme_prod=", sprod, "  somme_ref=", sref, "  ratio=", &
        sprod/max(abs(sref), 1.0e-300_DOUBLE)
    end if
  end subroutine print_entropy_diag

  ! ---------------------------------------------------------------------------
  ! 'vp_ep' : vitesse nodale vp ET energie nodale ep, calculees d'abord, puis
  ! injectees dans un solveur de Riemann 1D applique sous-face par sous-face.
  !
  ! Trois etages nettement separes, pour pouvoir travailler sur l'un sans toucher
  ! aux autres :
  !
  !   PASSE 1  vitesse nodale u_p. Identique a celle de 'vp_pp' : c'est le
  !            systeme nodal d'EUCCLHYD
  !              [somme_f A_f (lam_l+lam_r) n@n] u_p
  !                = somme_f A_f {p_l - p_r + lam_l v_l + lam_r v_r} n
  !            avec le traitement des parois, des faces z (boundary_2d) et des
  !            coins inchange.
  !
  !   PASSE 2  energie nodale e_p, rangee dans ep_node(i). << A DEFINIR >> :
  !            la routine nodal_energy ci-dessous est un PLACEHOLDER.
  !
  !   PASSE 3  assemblage des flux. Chaque sous-face appelle rs_vp_ep, qui recoit
  !            sol_w_l, sol_w_r, la normale, PLUS u_p.n et e_p, et rend les deux
  !            flux unilateraux. << A DEFINIR >> : rs_vp_ep est aussi un
  !            PLACEHOLDER, cf. sa propre en-tete.
  !
  ! Convention de signe des flux, identique au reste du module :
  !     rhs(:, idl) -= A_f/mass(idl) * flux_l
  !     rhs(:, idr) += A_f/mass(idr) * flux_r
  ! et sol = (tau, u, v, w, E) avec E l'energie totale MASSIQUE.
  subroutine compute_rhs_lagrange_vp_ep(mesh, sol, vp, dt, rhs, &
      n_bc, bc_type, bc_val, b2d, mass, gamma_arr, vp_is_imposed, &
      second_order, grad_v, grad_p, div_v)
    use linear_solver_module, only: inverse_3_by_3
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), intent(in) :: dt
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mass
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    logical, dimension(mesh%n_vert), intent(in) :: vp_is_imposed
    integer(kind=ENTIER), intent(in) :: n_bc
    character(len=255), dimension(n_bc) :: bc_type
    real(kind=DOUBLE), dimension(5, n_bc) :: bc_val
    logical, intent(in) :: b2d
    logical, intent(in) :: second_order
    real(kind=DOUBLE), dimension(3, 3, mesh%n_elems), intent(in) :: grad_v
    real(kind=DOUBLE), dimension(3, mesh%n_elems), intent(in) :: grad_p
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: div_v

    integer(kind=ENTIER) :: i, j, id_face, id_sub_face, nsfn, idl, idr, var
    real(kind=DOUBLE), dimension(3, 3) :: Mp, Mp_inv
    real(kind=DOUBLE), dimension(3) :: Rp, Bp, first_wall_norm
    real(kind=DOUBLE), dimension(5) :: flux_l, flux_r
    real(kind=DOUBLE) :: vl, vr, pl, pr, taul, taur, gl, gr
    real(kind=DOUBLE) :: lam_l, lam_r, v_bar, Pw, cross_z, area, v_target, v_et
    logical :: is_corner, have_first_wall_norm, is_zface
    real(kind=DOUBLE) :: pp_ep, s_num, s_den, p_bar, invlam
    real(kind=DOUBLE) :: s_presw, s_divw, s_den2
    real(kind=DOUBLE) :: q_a, q_b, q_c, disc, r1, r2, pp_lin
    real(kind=DOUBLE) :: mbar

    rhs = 0.0_DOUBLE
    if (.not. allocated(ep_node)) allocate(ep_node(mesh%n_vert))
    ep_node = 0.0_DOUBLE
    do i = 1, mesh%n_vert
      if (.not. vp_is_imposed(i)) vp(:, i) = 0.0_DOUBLE
    end do

    do i = 1, mesh%n_vert
      nsfn = mesh%vert(i)%n_sub_faces_neigh

      ! ---- PASSE 1 : vitesse nodale (identique a vp_pp) ----------------------
      Mp = 0.0_DOUBLE
      Rp = 0.0_DOUBLE
      Bp = 0.0_DOUBLE
      is_corner = .false.
      have_first_wall_norm = .false.
      first_wall_norm = 0.0_DOUBLE

      do j = 1, nsfn
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        idl = mesh%face(id_face)%left_neigh
        idr = mesh%face(id_face)%right_neigh
        area = mesh%sub_face(id_sub_face)%area

        gl = gamma_arr(idl)
        taul = sol(1, idl)
        pl = pressure(sol(:, idl), gl)
        vl = dot_product(sol(2:4, idl), mesh%face(id_face)%norm)
        if (second_order) call grp_half_step(mesh, idl, i, mesh%face(id_face)%norm, &
          dt, taul, gl, grad_v, grad_p, div_v, pl, vl)

        if (idr > 0) then
          gr = gamma_arr(idr)
          taur = sol(1, idr)
          pr = pressure(sol(:, idr), gr)
          vr = dot_product(sol(2:4, idr), mesh%face(id_face)%norm)
          if (second_order) call grp_half_step(mesh, idr, i, mesh%face(id_face)%norm, &
            dt, taur, gr, grad_v, grad_p, div_v, pr, vr)

          lam_l = max(sqrt(gl*pl*taul)/taul, &
            sqrt(max(0.0_DOUBLE, pr-pl)/taul), -(vr-vl)/taul)
          lam_r = max(sqrt(gr*pr*taur)/taur, &
            sqrt(max(0.0_DOUBLE, pl-pr)/taur), -(vr-vl)/taur)

          if (.not. vp_is_imposed(i)) then
            v_bar = (lam_l*vl + lam_r*vr - (pr - pl))/(lam_l + lam_r)
            Mp = Mp + area*(lam_l + lam_r) &
              * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
            Rp = Rp + area*(lam_l + lam_r)*mesh%face(id_face)%norm*v_bar
          end if
        else
          lam_l = sqrt(gl*pl*taul)/taul
          is_zface = (b2d .and. abs(mesh%face(id_face)%norm(3)) > 1e-8_DOUBLE)
          if (.not. vp_is_imposed(i)) then
            if (is_zface) then
              Mp = Mp + area*2.0_DOUBLE*lam_l &
                * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
            else if (idr < 0 .and. trim(bc_type(-idr)) == 'pressure') then
              Mp = Mp + area*lam_l &
                * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
              ! Etat fantome p_r = p_bc, v_r = v_l : avec Mp = A lambda, v.n = v_l + (p_l - p_bc)/lambda.
              ! p_l manquait, p_bc avait le signe inverse.
              Rp = Rp + area*(pl - bc_val(1,-idr) + lam_l*vl)*mesh%face(id_face)%norm
            else
              if (.not. have_first_wall_norm) then
                first_wall_norm = mesh%face(id_face)%norm
                have_first_wall_norm = .true.
              else
                cross_z = first_wall_norm(1)*mesh%face(id_face)%norm(2) &
                        - first_wall_norm(2)*mesh%face(id_face)%norm(1)
                if (abs(cross_z) > 1.0e-8_DOUBLE) is_corner = .true.
              end if
              Bp = Bp + area*mesh%face(id_face)%norm
              Mp = Mp + area*lam_l &
                * tensor_product(mesh%face(id_face)%norm, mesh%face(id_face)%norm)
              Rp = Rp + area*(pl + lam_l*vl)*mesh%face(id_face)%norm
            end if
          end if
        end if
      end do

      if (.not. vp_is_imposed(i)) then
        if (is_corner) then
          vp(:, i) = 0.0_DOUBLE
        else if (maxval(abs(Bp)) > 1e-8_DOUBLE) then
          call inverse_3_by_3(Mp, Mp_inv)
          Pw = dot_product(Rp, matmul(Mp_inv, Bp))/dot_product(Bp, matmul(Mp_inv, Bp))
          vp(:, i) = matmul(Mp_inv, Rp - Pw*Bp)
        else
          call inverse_3_by_3(Mp, Mp_inv)
          vp(:, i) = matmul(Mp_inv, Rp)
        end if
      end if


      ! ---- PASSE 1bis ('vp_pp_ep') : pression nodale p_p ---------------------
      ! Strictement la passe 2 de compute_rhs_lagrange_vp_pp : moyenne des
      ! pressions d'interface acoustiques, ponderee en 1/lambda, avec la MEME
      ! exclusion des faces z. Toute divergence entre les deux rendrait les
      ! colonnes vp_pp et vp_pp_ep incomparables.
      pp_ep = 0.0_DOUBLE
      if (pi_pp_ep) then
        s_num = 0.0_DOUBLE; s_den = 0.0_DOUBLE
        s_presw = 0.0_DOUBLE; s_divw = 0.0_DOUBLE; s_den2 = 0.0_DOUBLE
        q_a = 0.0_DOUBLE; q_b = 0.0_DOUBLE; q_c = 0.0_DOUBLE
        do j = 1, nsfn
          id_sub_face = mesh%vert(i)%sub_face_neigh(j)
          id_face = mesh%sub_face(id_sub_face)%mesh_face
          idl = mesh%face(id_face)%left_neigh
          idr = mesh%face(id_face)%right_neigh
          area = mesh%sub_face(id_sub_face)%area
          if (b2d .and. abs(mesh%face(id_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE) cycle
          ! correctif psw : les sous-faces de BORD sortent de la moyenne de p_p,
          ! comme les faces z, et pour la meme raison -- pas de second etat.
          if (pi_bord_psw .and. mesh%face(id_face)%right_neigh <= 0) cycle
          gl = gamma_arr(idl); taul = sol(1, idl)
          pl = pressure(sol(:, idl), gl)
          vl = dot_product(sol(2:4, idl), mesh%face(id_face)%norm)
          if (second_order) call grp_half_step(mesh, idl, i, mesh%face(id_face)%norm, &
            dt, taul, gl, grad_v, grad_p, div_v, pl, vl)
          if (idr > 0) then
            gr = gamma_arr(idr); taur = sol(1, idr)
            pr = pressure(sol(:, idr), gr)
            vr = dot_product(sol(2:4, idr), mesh%face(id_face)%norm)
            if (second_order) call grp_half_step(mesh, idr, i, mesh%face(id_face)%norm, &
              dt, taur, gr, grad_v, grad_p, div_v, pr, vr)
          else
            v_target = 0.0_DOUBLE
            if (idr < 0) then
              if (-idr <= n_bc) then
                if (trim(bc_type(-idr)) == 'piston') &
                  v_target = dot_product(bc_val(2:4, -idr), mesh%face(id_face)%norm)
              end if
            end if
            gr = gl; taur = taul; pr = pl
            vr = vl + 2.0_DOUBLE*(v_target - vl)
          end if
          lam_l = max(sqrt(max(0.0_DOUBLE, gl*pl*taul))/taul, &
            sqrt(max(0.0_DOUBLE, pr-pl)/taul), -(vr-vl)/taul)
          lam_r = max(sqrt(max(0.0_DOUBLE, gr*pr*taur))/taur, &
            sqrt(max(0.0_DOUBLE, pl-pr)/taur), -(vr-vl)/taur)
          invlam = 1.0_DOUBLE/lam_l + 1.0_DOUBLE/lam_r
          p_bar = (pl/lam_l + pr/lam_r - (vr - vl))/invlam
          s_num = s_num + area*invlam*p_bar
          s_den = s_den + area*invlam
          ! parts separees, pour la variante ponderee (voir pi_pp_wdiv)
          ! coefficients de somme A (u*_r^2 - u*_l^2) = 0, developpee en p_p
          ! avec u*_l = vl + (pl - p_p)/lam_l et u*_r = vr + (p_p - pr)/lam_r.
          q_a = q_a + area*(1.0_DOUBLE/(lam_r*lam_r) - 1.0_DOUBLE/(lam_l*lam_l))
          q_b = q_b + area*2.0_DOUBLE*( vr/lam_r - pr/(lam_r*lam_r) &
                                      + vl/lam_l + pl/(lam_l*lam_l) )
          q_c = q_c + area*( vr*vr - 2.0_DOUBLE*vr*pr/lam_r + pr*pr/(lam_r*lam_r) &
                           - vl*vl - 2.0_DOUBLE*vl*pl/lam_l - pl*pl/(lam_l*lam_l) )
          ! sommes des MOINDRES CARRES (voir pi_pp_wdiv)
          s_presw = s_presw + area*invlam*(pl/lam_l + pr/lam_r)
          s_divw  = s_divw  + area*invlam*(vr - vl)
          s_den2  = s_den2  + area*invlam*invlam
        end do
        if (s_den > 0.0_DOUBLE) pp_ep = s_num/s_den
        ! diagnostic passif : ecart entre les deux fermetures, TOUJOURS calcule
        if (s_den > 0.0_DOUBLE .and. s_den2 > 0.0_DOUBLE) then
          mbar = abs((s_presw - s_divw)/s_den2 - s_num/s_den) &
                 / max(abs(s_num/s_den), 1.0e-300_DOUBLE)
          ppmc_tot = ppmc_tot + 1
          if (mbar < 1.0e-6_DOUBLE) then
            ppmc_hist(1) = ppmc_hist(1) + 1
          else if (mbar < 1.0e-4_DOUBLE) then
            ppmc_hist(2) = ppmc_hist(2) + 1
          else if (mbar < 1.0e-2_DOUBLE) then
            ppmc_hist(3) = ppmc_hist(3) + 1
          else if (mbar < 1.0e-1_DOUBLE) then
            ppmc_hist(4) = ppmc_hist(4) + 1
          else
            ppmc_hist(5) = ppmc_hist(5) + 1
          end if
        end if
        if (pi_pp_veps .and. abs(q_b) > 0.0_DOUBLE) then
          ! quadratique franche ? on compare le terme en p_p^2 au terme lineaire
          ! a l echelle de la pression de repli.
          if (abs(q_a)*abs(pp_ep) > 1.0e-8_DOUBLE*abs(q_b)) then
            disc = q_b*q_b - 4.0_DOUBLE*q_a*q_c
            if (disc >= 0.0_DOUBLE) then
              r1 = (-q_b + sqrt(disc))/(2.0_DOUBLE*q_a)
              r2 = (-q_b - sqrt(disc))/(2.0_DOUBLE*q_a)
              ! racine la plus proche de la fermeture de moyenne, et POSITIVE
              if (r1 > 0.0_DOUBLE .and. &
                  (r2 <= 0.0_DOUBLE .or. abs(r1-pp_ep) <= abs(r2-pp_ep))) then
                pp_ep = r1
              else if (r2 > 0.0_DOUBLE) then
                pp_ep = r2
              end if
            end if
          else
            pp_lin = -q_c/q_b
            if (pp_lin > 0.0_DOUBLE) pp_ep = pp_lin
          end if
        end if
        if (pi_pp_wdiv .and. s_den2 > 0.0_DOUBLE) then
          ! MOINDRES CARRES : p_p minimise somme A (u*_r - u*_l)^2, alors que la
          ! formule ci-dessus annule la MOYENNE somme A (u*_r - u*_l). Les deux
          ! ferment le meme systeme sur-determine, aucune n'a de parametre, et
          ! elles COINCIDENT EXACTEMENT a impedances uniformes.
          ! La stationnarite donne somme A invlam (u*_r - u*_l) = 0, d ou
          !   p_p = somme A invlam (S - dv) / somme A invlam^2,  S = p_l/lam_l + p_r/lam_r.
          ! C est la ponderation par invlam -- et non une insertion a la main --
          ! qui fait que le terme en dv n est plus une divergence exacte.
          pp_ep = (s_presw - s_divw)/s_den2
        end if
      end if

      ! ---- PASSE 2 : energie nodale ------------------------------------------
      ep_node(i) = nodal_energy(mesh, i, nsfn, sol, gamma_arr, vp(:, i), b2d, &
        n_bc, bc_type, bc_val, dt, second_order, grad_v, grad_p, div_v, pp_ep)

      ! ---- PASSE 3 : assemblage des flux -------------------------------------
      do j = 1, nsfn
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        idl = mesh%face(id_face)%left_neigh
        idr = mesh%face(id_face)%right_neigh
        v_et = dot_product(vp(:, i), mesh%face(id_face)%norm)

        ! Les faces z (boundary_2d) sont exclues, exactement comme dans
        ! nodal_energy. C'est indispensable et pas cosmetique : le flux d'energie
        !     pv_c = p_c v_c -/+ lam_c (e_p - e_c)
        ! contient un terme SANS facteur vitesse, qui ne s'annule donc pas quand
        ! u_p.n = 0 -- ce qui est le cas sur une face z. Le flux de 'classic',
        ! Pi*(u_p.n), s'annule lui de lui-meme, d'ou sa conservation. Si ces faces
        ! entrent dans le flux mais pas dans la moyenne definissant e_p, la
        ! telescopie qui rend le schema conservatif est rompue : mesure, -4.5 %
        ! d'energie totale sur le disque, -3.0 % sur le quad. Avec l'exclusion,
        ! la conservation est EXACTE.
        if (b2d .and. abs(mesh%face(id_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE) cycle

        if (idr > 0) then
          if (second_order) then
            pl = pressure(sol(:, idl), gamma_arr(idl))
            vl = dot_product(sol(2:4, idl), mesh%face(id_face)%norm)
            pr = pressure(sol(:, idr), gamma_arr(idr))
            vr = dot_product(sol(2:4, idr), mesh%face(id_face)%norm)
            call grp_half_step(mesh, idl, i, mesh%face(id_face)%norm, dt, &
              sol(1, idl), gamma_arr(idl), grad_v, grad_p, div_v, pl, vl)
            call grp_half_step(mesh, idr, i, mesh%face(id_face)%norm, dt, &
              sol(1, idr), gamma_arr(idr), grad_v, grad_p, div_v, pr, vr)
            if (pi_pp_ep) then
              call rs_vp_ep(sol(:, idl), sol(:, idr), gamma_arr(idl), gamma_arr(idr), &
                mesh%face(id_face)%norm, v_et, ep_node(i), flux_l, flux_r, &
                pl, vl, pr, vr, pp_nodal=pp_ep)
            else
              call rs_vp_ep(sol(:, idl), sol(:, idr), gamma_arr(idl), gamma_arr(idr), &
                mesh%face(id_face)%norm, v_et, ep_node(i), flux_l, flux_r, &
                pl, vl, pr, vr)
            end if
          else
            if (pi_pp_ep) then
              call rs_vp_ep(sol(:, idl), sol(:, idr), gamma_arr(idl), gamma_arr(idr), &
                mesh%face(id_face)%norm, v_et, ep_node(i), flux_l, flux_r, &
                pp_nodal=pp_ep)
            else
              call rs_vp_ep(sol(:, idl), sol(:, idr), gamma_arr(idl), gamma_arr(idr), &
                mesh%face(id_face)%norm, v_et, ep_node(i), flux_l, flux_r)
            end if
          end if
        else
          ! bord : etat miroir, comme dans la passe de pression de vp_pp
          v_target = 0.0_DOUBLE
          if (idr < 0) then
            if (-idr <= n_bc) then
              if (trim(bc_type(-idr)) == 'piston') then
                v_target = dot_product(bc_val(2:4, -idr), mesh%face(id_face)%norm)
              end if
            end if
          end if
          if (pi_pp_ep) then
            call rs_vp_ep(sol(:, idl), mirror_state(sol(:, idl), &
              mesh%face(id_face)%norm, v_target), gamma_arr(idl), gamma_arr(idl), &
              mesh%face(id_face)%norm, v_et, ep_node(i), flux_l, flux_r, &
              pp_nodal=pp_ep, bord=.true.)
          else
            call rs_vp_ep(sol(:, idl), mirror_state(sol(:, idl), &
              mesh%face(id_face)%norm, v_target), gamma_arr(idl), gamma_arr(idl), &
              mesh%face(id_face)%norm, v_et, ep_node(i), flux_l, flux_r, &
              bord=.true.)
          end if
        end if

        do var = 1, 5
          rhs(var, idl) = rhs(var, idl) &
            - mesh%sub_face(id_sub_face)%area/mass(idl)*flux_l(var)
        end do
        if (idr > 0) then
          do var = 1, 5
            rhs(var, idr) = rhs(var, idr) &
              + mesh%sub_face(id_sub_face)%area/mass(idr)*flux_r(var)
          end do
        end if
      end do
    end do
  end subroutine compute_rhs_lagrange_vp_ep

  ! ===========================================================================
  ! << A DEFINIR >>  Energie nodale e_p du schema 'vp_ep'.
  !
  ! Appelee une fois par noeud, APRES que la vitesse nodale u_p soit connue.
  ! Elle a acces a toutes les sous-faces du noeud, donc aux etats des mailles
  ! voisines, a leurs impedances et a u_p.
  !
  ! Formule retenue (cf. corps) : moyenne ponderee en LAMBDA de l'energie totale
  ! massique, corrigee par le saut du FLUX d'energie p*v_n. C'est la duale exacte
  ! de la vitesse nodale. L'energie utilisee est sol(5), l'energie TOTALE
  ! massique ; si l'energie INTERNE etait voulue, c'est la seule ligne a changer.
  !
  ! Rappels utiles :
  !   sol(:,c) = (tau_c, u_c, v_c, w_c, E_c)   E = energie totale MASSIQUE
  !   e_int    = E - |u|^2/2
  !   p        = (g-1)/tau * (E - |u|^2/2)      [fonction pressure]
  !   lambda_c = sqrt(g p tau)/tau              [impedance rho a]
  ! Les faces z sont exclues en boundary_2d comme dans la passe de pression de
  ! 'vp_pp' (test strict |n_z| > 1 - 1e-8).
  ! ===========================================================================
  function nodal_energy(mesh, i_vert, nsfn, sol, gamma_arr, u_p, b2d, &
      n_bc, bc_type, bc_val, dt, second_order, grad_v, grad_p, div_v, pp_in) result(e_p)
    implicit none
    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: i_vert, nsfn, n_bc
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    real(kind=DOUBLE), dimension(3), intent(in) :: u_p
    logical, intent(in) :: b2d
    character(len=255), dimension(n_bc) :: bc_type
    real(kind=DOUBLE), dimension(5, n_bc) :: bc_val
    real(kind=DOUBLE), intent(in) :: dt
    logical, intent(in) :: second_order
    real(kind=DOUBLE), dimension(3, 3, mesh%n_elems), intent(in) :: grad_v
    real(kind=DOUBLE), dimension(3, mesh%n_elems), intent(in) :: grad_p
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: div_v
    ! pression nodale, requise par 'vp_pp_epsp' pour le glissement acoustique
    real(kind=DOUBLE), intent(in), optional :: pp_in
    real(kind=DOUBLE) :: e_p

    integer(kind=ENTIER) :: j, id_sub_face, id_face, idl, idr
    real(kind=DOUBLE) :: area, gl, gr, taul, taur, pl, pr, lam_l, lam_r
    real(kind=DOUBLE) :: e_bar, num, den, v_target, vl, vr, mu_l, mu_r
    real(kind=DOUBLE) :: dsl, dsr, xl, xr
    real(kind=DOUBLE), dimension(5) :: solr

    num = 0.0_DOUBLE
    den = 0.0_DOUBLE
    do j = 1, nsfn
      id_sub_face = mesh%vert(i_vert)%sub_face_neigh(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      idl = mesh%face(id_face)%left_neigh
      idr = mesh%face(id_face)%right_neigh
      area = mesh%sub_face(id_sub_face)%area
      if (b2d .and. abs(mesh%face(id_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE) cycle
      ! Les sous-faces de BORD sortent aussi, pour la meme raison que les faces z
      ! et avec la meme consequence si on l'oublie. Le flux d'energie porte le
      ! terme `-mu (e_p - e_c)`, SANS facteur vitesse, qui ne s'annule donc pas
      ! quand u_p.n = 0 -- ce qui est le cas a une paroi impermeable. Une
      ! sous-face de bord n'ayant PAS de partenaire pour telescoper, ce terme est
      ! une fuite nette d'energie hors du domaine, et a une paroi isotherme c'est
      ! de surcroit un faux flux de chaleur non visqueux.
      ! MESURE, sedov_quad_n21 en domaine ferme : derive de l'energie totale
      ! 1.7e-06 pour vp_ep et 6.3e-07 pour vp_pp_ep, contre 6e-16 pour classic.
      ! Le flux de bord est donc ramene a (pression de quantite de mouvement)
      ! x (u_p.n), nul a une paroi, et ces sous-faces sortent d'ici pour que le
      ! telescopage INTERIEUR reste exact -- inclure dans l'un sans l'autre est
      ! precisement ce qui le romprait.
      if (idr <= 0) cycle

      gl = gamma_arr(idl)
      taul = sol(1, idl)
      pl = pressure(sol(:, idl), gl)
      if (idr > 0) then
        solr = sol(:, idr)
        gr = gamma_arr(idr)
      else
        v_target = 0.0_DOUBLE
        if (idr < 0) then
          if (-idr <= n_bc) then
            if (trim(bc_type(-idr)) == 'piston') &
              v_target = dot_product(bc_val(2:4, -idr), mesh%face(id_face)%norm)
          end if
        end if
        solr = mirror_state(sol(:, idl), mesh%face(id_face)%norm, v_target)
        gr = gl
      end if
      taur = solr(1)
      pr = pressure(solr, gr)
      vl = dot_product(sol(2:4, idl), mesh%face(id_face)%norm)
      vr = dot_product(solr(2:4), mesh%face(id_face)%norm)
      ! ordre 2 : meme extrapolation GRP que dans la passe de vitesse. Seuls p et
      ! v.n sont reconstruits (ce sont les seuls gradients disponibles) ; e_l et
      ! e_r restent les moyennes de maille, comme dans 'classic' dont la
      ! fermeture ne fait jamais intervenir l'energie de maille.
      if (second_order) then
        call grp_half_step(mesh, idl, i_vert, mesh%face(id_face)%norm, &
          dt, taul, gl, grad_v, grad_p, div_v, pl, vl)
        if (idr > 0) call grp_half_step(mesh, idr, i_vert, mesh%face(id_face)%norm, &
          dt, taur, gr, grad_v, grad_p, div_v, pr, vr)
      end if

      if (lambda_acoustic_only) then
        lam_l = sqrt(max(0.0_DOUBLE, gl*pl*taul))/taul
        lam_r = sqrt(max(0.0_DOUBLE, gr*pr*taur))/taur
      else
        lam_l = max(sqrt(max(0.0_DOUBLE, gl*pl*taul))/taul, &
          sqrt(max(0.0_DOUBLE, pr-pl)/taul), -(vr-vl)/taul)
        lam_r = max(sqrt(max(0.0_DOUBLE, gr*pr*taur))/taur, &
          sqrt(max(0.0_DOUBLE, pl-pr)/taur), -(vr-vl)/taur)
      end if
      if (lam_l + lam_r <= 0.0_DOUBLE) cycle

      ! Energie nodale, duale de la vitesse nodale : ponderation en LAMBDA (et
      ! non en 1/lambda comme la pression), et terme de saut porte par le FLUX
      ! d'energie p*v_n.
      !     v_bar = (lam_l v_l + lam_r v_r - (p_r - p_l))       / (lam_l+lam_r)
      !     e_bar = (lam_l e_l + lam_r e_r - (p_r v_r - p_l v_l))/ (lam_l+lam_r)
      ! Les deux contraintes qui avaient elimine les formes naives sont ici
      ! satisfaites : lam*e est en rho v^3 comme Delta(p v) (homogeneite), et e
      ! comme Delta(p v) sont INVARIANTS par retournement de la face, comme doit
      ! l'etre e_bar (parite).
      ! mu = lam * phi (phi = 1 hors mode wip). La MEME substitution est faite
      ! dans rs_vp_ep, ce qui preserve la telescopie donc la conservation exacte.
      mu_l = lam_l*wip_phi(dot_product(u_p, mesh%face(id_face)%norm), vl, taul, pl, gl, pr - pl)
      mu_r = lam_r*wip_phi(dot_product(u_p, mesh%face(id_face)%norm), vr, taur, pr, gr, pr - pl)
      if (mu_l + mu_r <= 0.0_DOUBLE) cycle
      if ((pi_pp_eps2 .or. pi_pp_eps2m) .and. present(pp_in)) then
        ! ep_node porte kappa : num accumule -N, den la ponderation choisie.
        dsl = (pl - pp_in)/lam_l
        dsr = (pp_in - pr)/lam_r
        num = num - area*(pp_in*(vr - vl) &
              + 0.5_DOUBLE*dsr*dsr*lam_r + 0.5_DOUBLE*dsl*dsl*lam_l)
        if (pi_pp_eps2m) then
          xl = sol(5, idl) - 0.5_DOUBLE*dot_product(sol(2:4, idl), sol(2:4, idl))
          xr = solr(5) - 0.5_DOUBLE*dot_product(solr(2:4), solr(2:4))
          den = den + area*(lam_l*xl + lam_r*xr)
        else
          den = den + area*(lam_l + lam_r)
        end if
        cycle
      end if
      if (pi_pp_ppus) then
        ! eps_p : combinaison convexe des energies INTERNES. Le travail de
        ! pression p_p (u*_c.n) telescope seul, il ne doit donc PAS figurer ici.
        num = num + area*(mu_l*(sol(5, idl) &
                   - 0.5_DOUBLE*dot_product(sol(2:4, idl), sol(2:4, idl))) &
                   + mu_r*(solr(5) &
                   - 0.5_DOUBLE*dot_product(solr(2:4), solr(2:4))))
        den = den + area*(mu_l + mu_r)
        cycle
      end if
      if (pi_pivp) then
        ! e_p = somme A[mu_l e_l + mu_r e_r] / somme A(mu_l+mu_r) : le travail de
        ! pression a disparu du numerateur (il telescope seul), donc e_p est une
        ! combinaison convexe et ne peut pas diverger quand mu devient petit.
        num = num + area*(mu_l*sol(5,idl) + mu_r*solr(5))
        den = den + area*(mu_l + mu_r)
        cycle
      end if
      if (pi_pp_epsp .and. present(pp_in)) then
        ! X_c = eps_c - v_c delta_c - delta_c^2/2, le meme signe des deux cotes.
        ! Le glissement porte lam (relation de quantite de mouvement) et non mu.
        dsl = (pl - pp_in)/lam_l
        dsr = (pp_in - pr)/lam_r
        xl = sol(5, idl) - 0.5_DOUBLE*dot_product(sol(2:4, idl), sol(2:4, idl)) &
             - vl*dsl - 0.5_DOUBLE*dsl*dsl
        xr = solr(5) - 0.5_DOUBLE*dot_product(solr(2:4), solr(2:4)) &
             - vr*dsr - 0.5_DOUBLE*dsr*dsr
        e_bar = (mu_l*xl + mu_r*xr - (pr*vr - pl*vl)) / (mu_l + mu_r)
      else if (pi_ep_interne) then
        ! eps = E - |u|^2/2, des DEUX cotes. La meme substitution est faite dans
        ! rs_vp_ep, sans quoi le telescopage tombe.
        e_bar = (mu_l*(sol(5, idl) - 0.5_DOUBLE*dot_product(sol(2:4, idl), sol(2:4, idl))) &
               + mu_r*(solr(5) - 0.5_DOUBLE*dot_product(solr(2:4), solr(2:4))) &
               - (pr*vr - pl*vl)) / (mu_l + mu_r)
      else
        e_bar = (mu_l*sol(5, idl) + mu_r*solr(5) - (pr*vr - pl*vl)) / (mu_l + mu_r)
      end if

      num = num + area*(mu_l + mu_r)*e_bar
      den = den + area*(mu_l + mu_r)
    end do

    if (den > 0.0_DOUBLE) then
      e_p = num/den
    else
      e_p = 0.0_DOUBLE
    end if
  end function nodal_energy

  ! ===========================================================================
  ! << A DEFINIR >>  Solveur de Riemann 1D du schema 'vp_ep'.
  !
  ! Entrees : les deux etats de maille sol_w_l, sol_w_r et leurs gamma, la
  ! normale de la face, PLUS les deux quantites nodales deja calculees --
  ! vn_nodal = u_p.n  et  e_nodal = e_p.
  !
  ! Sorties : les deux flux unilateraux, avec la convention du module
  !     rhs(:, idl) -= A/mass(idl) * flux_l
  !     rhs(:, idr) += A/mass(idr) * flux_r
  ! donc flux = (-u_p.n, Pi*n, Pi*u_p.n) pour une fermeture en pression Pi.
  !
  ! Le PLACEHOLDER ci-dessous est la fermeture UNILATERALE de 'classic',
  !     Pi_l = p_l - lam_l (u_p.n - v_l),   Pi_r = p_r + lam_r (u_p.n - v_r)
  ! qui rend le schema conservatif (l'equilibre nodal est deja satisfait par la
  ! vitesse de la passe 1) et a production d'entropie positive. e_nodal N'EST
  ! PAS ENCORE UTILISE : c'est la ligne a ecrire. Tel quel, 'vp_ep' reproduit
  ! donc 'classic' -- c'est voulu, ca donne une reference qui tourne.
  ! ===========================================================================
  subroutine rs_vp_ep(sol_w_l, sol_w_r, g_l, g_r, n, vn_nodal, e_nodal, &
      flux_l, flux_r, pl_rec, vl_rec, pr_rec, vr_rec, pp_nodal, bord)
    implicit none
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), intent(in) :: g_l, g_r, vn_nodal, e_nodal
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), dimension(5), intent(out) :: flux_l, flux_r
    ! ordre 2 : p et v.n deja extrapoles au noeud par l'appelant. Absents = ordre 1.
    real(kind=DOUBLE), intent(in), optional :: pl_rec, vl_rec, pr_rec, vr_rec
    ! 'vp_pp_ep' : present => la quantite de mouvement prend la pression NODALE
    ! au lieu des deux Pi unilateraux. L'energie, elle, ne change pas.
    real(kind=DOUBLE), intent(in), optional :: pp_nodal
    ! Sous-face de BORD : le terme en (e_p - e) saute et le flux d'energie
    ! redevient (pression de quantite de mouvement) x (u_p.n), nul a une paroi
    ! impermeable et egal a p v_piston sur un piston. Voir nodal_energy.
    logical, intent(in), optional :: bord

    real(kind=DOUBLE) :: pl, pr, taul, taur, vl, vr, lam_l, lam_r, pi_l, pi_r
    real(kind=DOUBLE) :: pvl, pvr, el, er, mu_l, mu_r, cl, cr, cap
    real(kind=DOUBLE) :: dsl, dsr
    real(kind=DOUBLE) :: s_prod, s_ref
    real(kind=DOUBLE) :: pm_l, pm_r, pm_rel
    real(kind=DOUBLE), dimension(3) :: utl, utr, ust, ft
    real(kind=DOUBLE) :: dut, nu_l, nu_r
    logical :: bord_psw_actif
    logical :: bord_loc

    bord_psw_actif = .false.
    taul = sol_w_l(1)
    taur = sol_w_r(1)
    pl = pressure(sol_w_l, g_l)
    pr = pressure(sol_w_r, g_r)
    vl = dot_product(sol_w_l(2:4), n)
    vr = dot_product(sol_w_r(2:4), n)
    if (present(pl_rec)) pl = pl_rec
    if (present(vl_rec)) vl = vl_rec
    if (present(pr_rec)) pr = pr_rec
    if (present(vr_rec)) vr = vr_rec
    ! meme definition de lambda que 'classic' (et meme bascule
    ! lambda_acoustic_only), pour que le placeholder soit une vraie reference
    if (lambda_acoustic_only) then
      lam_l = sqrt(max(0.0_DOUBLE, g_l*pl*taul))/taul
      lam_r = sqrt(max(0.0_DOUBLE, g_r*pr*taur))/taur
    else
      lam_l = max(sqrt(max(0.0_DOUBLE, g_l*pl*taul))/taul, &
        sqrt(max(0.0_DOUBLE, pr-pl)/taul), -(vr-vl)/taul)
      lam_r = max(sqrt(max(0.0_DOUBLE, g_r*pr*taur))/taur, &
        sqrt(max(0.0_DOUBLE, pl-pr)/taur), -(vr-vl)/taur)
    end if

    pi_l = pl - lam_l*(vn_nodal - vl)
    pi_r = pr + lam_r*(vn_nodal - vr)

    el = sol_w_l(5)
    er = sol_w_r(5)
    if (pi_pp_epsp .and. present(pp_nodal)) then
      ! LA MEME substitution que dans nodal_energy, sans quoi le telescopage
      ! tombe et la conservation avec lui.
      dsl = (pl - pp_nodal)/lam_l
      dsr = (pp_nodal - pr)/lam_r
      el = el - 0.5_DOUBLE*dot_product(sol_w_l(2:4), sol_w_l(2:4)) &
           - vl*dsl - 0.5_DOUBLE*dsl*dsl
      er = er - 0.5_DOUBLE*dot_product(sol_w_r(2:4), sol_w_r(2:4)) &
           - vr*dsr - 0.5_DOUBLE*dsr*dsr
    else if (pi_ep_interne) then
      el = el - 0.5_DOUBLE*dot_product(sol_w_l(2:4), sol_w_l(2:4))
      er = er - 0.5_DOUBLE*dot_product(sol_w_r(2:4), sol_w_r(2:4))
    end if

    mu_l = lam_l*wip_phi(vn_nodal, vl, taul, pl, g_l, pr - pl)
    mu_r = lam_r*wip_phi(vn_nodal, vr, taur, pr, g_r, pr - pl)
    if (pi_pp_ppus .and. present(pp_nodal)) then
      ! pvl = p_p (v_l + ds_l) - mu_l (eps_p - eps_l). el/er portent deja
      ! l energie INTERNE (branche pi_ep_interne plus haut), e_nodal aussi.
      dsl = (pl - pp_nodal)/lam_l
      dsr = (pp_nodal - pr)/lam_r
      cl = (pl - pp_nodal)*vl - pp_nodal*dsl + mu_l*(e_nodal - el)
      cr = (pp_nodal - pr)*vr + pp_nodal*dsr + mu_r*(e_nodal - er)
    else if (pi_pivp .and. present(pp_nodal)) then
      ! pvl = Pi_l (v_p.n) - mu_l (e_p - e_l), avec pvl = pl*vl - cl
      cl = pl*vl - pi_l*vn_nodal + mu_l*(e_nodal - el)
      cr = pi_r*vn_nodal - pr*vr + mu_r*(e_nodal - er)
    else if (pi_pp_ppvn .and. present(pp_nodal)) then
      ! pvl = pvr = p_p (v_p.n)
      cl = pl*vl - pp_nodal*vn_nodal
      cr = pp_nodal*vn_nodal - pr*vr
    else if (pi_pp_vpps .and. present(pp_nodal)) then
      ! pvl = Pi_l (v_p.n) et pvr = Pi_r (v_p.n), avec pvl = pl*vl - cl et
      ! pvr = pr*vr + cr.
      cl = pl*vl - pi_l*vn_nodal
      cr = pi_r*vn_nodal - pr*vr
    else if (pi_pp_epv .and. present(pp_nodal)) then
      ! (pv)*_c = p_p v*_c, avec v*_c la vitesse etoile ACOUSTIQUE de la maille.
      ! Aucune energie nodale n'intervient : e_nodal n'est pas utilise.
      dsl = (pl - pp_nodal)/lam_l
      dsr = (pp_nodal - pr)/lam_r
      cl = (pl - pp_nodal)*vl - pp_nodal*dsl
      cr = (pp_nodal - pr)*vr + pp_nodal*dsr
    else if ((pi_pp_eps2 .or. pi_pp_eps2m) .and. present(pp_nodal)) then
      ! e_nodal porte KAPPA. (pv)*_l = p_p v_l - D_l^2/(2 lam_l) - lam_l X_l kappa
      ! est obtenu en posant cl tel que pvl = p_l v_l - cl.
      dsl = (pl - pp_nodal)/lam_l
      dsr = (pp_nodal - pr)/lam_r
      if (pi_pp_eps2m) then
        cl = (pl - pp_nodal)*vl + 0.5_DOUBLE*dsl*dsl*lam_l &
             + lam_l*e_nodal*(sol_w_l(5) - 0.5_DOUBLE*dot_product(sol_w_l(2:4), sol_w_l(2:4)))
        cr = (pp_nodal - pr)*vr + 0.5_DOUBLE*dsr*dsr*lam_r &
             + lam_r*e_nodal*(sol_w_r(5) - 0.5_DOUBLE*dot_product(sol_w_r(2:4), sol_w_r(2:4)))
      else
        cl = (pl - pp_nodal)*vl + 0.5_DOUBLE*dsl*dsl*lam_l + lam_l*e_nodal
        cr = (pp_nodal - pr)*vr + 0.5_DOUBLE*dsr*dsr*lam_r + lam_r*e_nodal
      end if
    else
      cl = mu_l*(e_nodal - el)
      cr = mu_r*(e_nodal - er)
    end if
    if (vp_ep_wip .and. vp_ep_wip_clip > 0.0_DOUBLE) then
      cap = vp_ep_wip_clip*lam_l*max(el - 0.5_DOUBLE*dot_product(sol_w_l(2:4), sol_w_l(2:4)), 0.0_DOUBLE)
      cl = sign(min(abs(cl), cap), cl)
      cap = vp_ep_wip_clip*lam_r*max(er - 0.5_DOUBLE*dot_product(sol_w_r(2:4), sol_w_r(2:4)), 0.0_DOUBLE)
      cr = sign(min(abs(cr), cap), cr)
    end if
    pvl = pl*vl - cl
    pvr = pr*vr + cr

    ! --- dissipation tangentielle (voir tang_mode) -----------------------
    if (tang_mode > 0) then
      utl = sol_w_l(2:4) - dot_product(sol_w_l(2:4), n)*n
      utr = sol_w_r(2:4) - dot_product(sol_w_r(2:4), n)*n
      dut = sqrt(dot_product(utr - utl, utr - utl))
      if (tang_mode == 1) then
        nu_l = lam_l ; nu_r = lam_r
      else
        nu_l = dut/taul ; nu_r = dut/taur
      end if
      tang_n = tang_n + 1
      if (nu_l + nu_r > 0.0_DOUBLE) then
        ust = (nu_l*utl + nu_r*utr)/(nu_l + nu_r)
        ft  = nu_l*(ust - utl)          ! = -nu_r*(ust - utr), donc UNIVALUE
        if (dot_product(ft, ft) > 0.0_DOUBLE) tang_nact = tang_nact + 1
      else
        ust = 0.5_DOUBLE*(utl + utr) ; ft = 0.0_DOUBLE
      end if
    else
      ft = 0.0_DOUBLE ; ust = 0.0_DOUBLE
    end if

    flux_l(1)   = -vn_nodal
    flux_l(2:4) = pi_l*n - ft
    flux_l(5)   = pvl - dot_product(ft, ust)

    flux_r(1)   = -vn_nodal
    flux_r(2:4) = pi_r*n - ft
    flux_r(5)   = pvr - dot_product(ft, ust)

    if (present(bord)) then
      if (bord) then
        if (pi_bord_psw) then
          flux_l(5) = pi_l*vn_nodal
          flux_r(5) = pi_r*vn_nodal
          bord_psw_actif = .true.
        else if (present(pp_nodal)) then
          flux_l(5) = pp_nodal*vn_nodal
          flux_r(5) = pp_nodal*vn_nodal
        else
          flux_l(5) = pi_l*vn_nodal
          flux_r(5) = pi_r*vn_nodal
        end if
      end if
    end if

    ! --- marge de positivite (passif, voir pos_hist) ---------------------
    if (present(pp_nodal)) then
      dsl = (pl - pp_nodal)/lam_l
      dsr = (pp_nodal - pr)/lam_r
      pm_l = e_nodal - 0.5_DOUBLE*(vl + dsl)**2
      pm_r = e_nodal - 0.5_DOUBLE*(vr + dsr)**2
      pm_rel = min(pm_l, pm_r)/max(abs(e_nodal), 1.0e-300_DOUBLE)
      pos_tot = pos_tot + 1
      pos_worst = min(pos_worst, pm_rel)
      if (pm_rel <= 0.0_DOUBLE) then
        pos_hist(1) = pos_hist(1) + 1
      else if (pm_rel < 1.0e-3_DOUBLE) then
        pos_hist(2) = pos_hist(2) + 1
      else if (pm_rel < 1.0e-2_DOUBLE) then
        pos_hist(3) = pos_hist(3) + 1
      else if (pm_rel < 1.0e-1_DOUBLE) then
        pos_hist(4) = pos_hist(4) + 1
      else
        pos_hist(5) = pos_hist(5) + 1
      end if
    end if
    if (present(pp_nodal)) then
      if (.not. bord_psw_actif) then
        flux_l(2:4) = pp_nodal*n - ft
        flux_r(2:4) = pp_nodal*n - ft
      end if
      ! Le terme qui separe ce bilan d'energie interne de celui de vp_ep :
      ! somme_f A (p_p - Pi_c) u_c, nul au repos et sans signe sinon.
      psvp_n_sf = psvp_n_sf + 2
      psvp_sum = psvp_sum + (pp_nodal - pi_l)*vl - (pp_nodal - pi_r)*vr
      if ((pp_nodal - pi_l)*vl < 0.0_DOUBLE) psvp_n_neg = psvp_n_neg + 1
      if (-(pp_nodal - pi_r)*vr < 0.0_DOUBLE) psvp_n_neg = psvp_n_neg + 1
      psvp_ref = psvp_ref + lam_l*(vn_nodal - vl)**2 + lam_r*(vn_nodal - vr)**2
      ! Defaut de telescopage interieur. Nul pour toute fermeture ou les deux
      ! flux d energie sont la meme expression ; non nul pour pi_pp_vpps.
      vpps_tel = vpps_tel + (flux_l(5) - flux_r(5))
      vpps_tel_abs = vpps_tel_abs + abs(flux_l(5)) + abs(flux_r(5))
      ! --- diagnostic S, sous-faces INTERIEURES seulement : sur une paroi ou une
      ! face z, v* = 0 et S serait negatif mecaniquement sans rien mesurer.
      bord_loc = .false.
      if (present(bord)) bord_loc = bord
      if (.not. bord_loc) then
        sdiag_tot = sdiag_tot + 1
        s_prod = 0.5_DOUBLE*((1.0_DOUBLE/lam_l + 1.0_DOUBLE/lam_r)*pp_nodal**2 &
                 - pl*pl/lam_l - pr*pr/lam_r)
        s_ref = max(pl, pr)*0.5_DOUBLE*(lam_l*taul + lam_r*taur)
        if (s_ref > 0.0_DOUBLE) then
          s_prod = s_prod/s_ref
        else
          s_prod = 0.0_DOUBLE
        end if
        if (s_prod < -1.0e-1_DOUBLE) then
          sdiag_n(1) = sdiag_n(1) + 1
        else if (s_prod < -1.0e-2_DOUBLE) then
          sdiag_n(2) = sdiag_n(2) + 1
        else if (s_prod < -1.0e-3_DOUBLE) then
          sdiag_n(3) = sdiag_n(3) + 1
        else if (s_prod < -1.0e-6_DOUBLE) then
          sdiag_n(4) = sdiag_n(4) + 1
        else if (s_prod < 0.0_DOUBLE) then
          sdiag_n(5) = sdiag_n(5) + 1
        else
          sdiag_n(6) = sdiag_n(6) + 1
        end if
      end if
    end if
  end subroutine rs_vp_ep

  ! Etat miroir d'une paroi : meme tau, meme energie interne, composante normale
  ! de la vitesse reflechie autour de v_target (0 pour une paroi fixe).
  function mirror_state(sol_c, n, v_target) result(solm)
    implicit none
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_c
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), intent(in) :: v_target
    real(kind=DOUBLE), dimension(5) :: solm
    real(kind=DOUBLE) :: vn

    vn = dot_product(sol_c(2:4), n)
    solm = sol_c
    solm(2:4) = sol_c(2:4) + 2.0_DOUBLE*(v_target - vn)*n
    ! E est inchange en norme de vitesse (reflexion), donc sol(5) reste valable
  end function mirror_state

  ! Facteur phi_c du schema 'vp_ep_wip' : il doit s'annuler au CONTACT, ou la
  ! vitesse et la pression sont continues et seule l'energie saute. Utilise a
  ! l'identique dans nodal_energy et dans rs_vp_ep, faute de quoi la conservation
  ! exacte est perdue.
  function wip_phi(vn_nodal, v_c, tau_c, p_c, g_c, dp) result(phi)
    implicit none
    real(kind=DOUBLE), intent(in) :: vn_nodal, v_c, tau_c, p_c, g_c, dp
    real(kind=DOUBLE) :: phi, a_c, num

    if (.not. vp_ep_wip) then
      phi = 1.0_DOUBLE
      return
    end if
    a_c = sqrt(max(0.0_DOUBLE, g_c*p_c*tau_c))
    select case (trim(vp_ep_wip_mode))
    case ('dv')
      num = abs(vn_nodal - v_c)
      phi = num/(num + max(a_c, 1.0e-30_DOUBLE))
    case ('dp')
      num = abs(dp)
      phi = num/(num + max(abs(p_c), 1.0e-30_DOUBLE))
    case default
      phi = 1.0_DOUBLE
    end select
    phi = max(phi, vp_ep_wip_floor)
  end function wip_phi

  function tensor_product(a,b) result(c)
    implicit none

    real(kind=DOUBLE), dimension(:) :: a, b
    real(kind=DOUBLE), dimension(size(a), size(b)) :: c

    integer(kind=ENTIER) :: i, j

    do i=1, size(a)
      do j=1, size(b)
        c(i, j) = a(i)*b(j)
      end do
    end do
  end function tensor_product

  ! Impedance du schema : acoustique, saut de pression, saut de vitesse.
  ! Partagee par le flux et compute_dt pour qu'ils ne puissent pas diverger.
  ! `only` : 'A' = borne acoustique seule (sans v_p, calculable en tete de pas),
  ! 'V' = borne de volume seule, absent = les deux.
  subroutine compute_dt(mesh, sol, dt, cfl, vp, me, num_procs, gamma_arr, only, bc_p)
    use, intrinsic :: ieee_arithmetic
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: vp
    real(kind=DOUBLE), intent(in) :: cfl
    real(kind=DOUBLE), intent(inout) :: dt
    integer(kind=ENTIER) :: me, num_procs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    character(len=1), intent(in), optional :: only
    ! Pression imposee par condition de bord (< 0 si ce n'est pas une condition 'pressure').
    real(kind=DOUBLE), dimension(:), intent(in), optional :: bc_p
    logical :: do_acou, do_vol

    ! Coefficient de la borne en VITESSE dt <= CV*omega_c/somme_f(a_pcf |u_p.n|),
    ! celle qui limite le changement relatif de volume d'une maille par pas.
    ! Passe de 0.8 a 0.2 : la reecriture de cette borne (volume complet et SOMME
    ! sur les sous-faces, au lieu du volume moyen de sous-element et de l'aire
    ! d'UNE face) l'a rendue ~1.9x plus permissive a CV egal. Sur Saltzman cette
    ! marge suffisait a faire TRAVERSER les parois y=0 et y=0.1 au maillage
    ! (y de -0.031 a 0.108 au lieu de rester dans [0, 0.1]) -- sans aucun NaN ni
    ! code d'erreur, donc invisible au depouillement habituel. Des 0.4 la geometrie
    ! est retablie et les autres cas sont inchanges a 0.1 % pres (Sedov, Noh, Sod
    ! bi-materiau, Sod cylindrique) ; 0.2 prend une marge supplementaire, le pic de
    ! Saltzman n'etant pas encore convergie en temps a 0.4 (21.3 a CV=0.4, 18.9
    ! en divisant encore le pas).
    real(kind=DOUBLE), parameter :: CV = 0.2_DOUBLE
    integer(kind=ENTIER) :: i, j, k, mpi_ierr
    real(kind=DOUBLE) :: bnd, dt_nan_val
    logical :: dt_nan_seen
    integer :: inan
    integer(kind=ENTIER) :: id_sub_face, id_face
    integer(kind=ENTIER) :: id_vert, id_sub_elem
    real(kind=DOUBLE) :: dx, c, rho, lambda, mass_c, sum_al, sum_av
    ! lambda par SOUS-FACE, celui du schema (sauts compris) et non rho*a seul.
    integer(kind=ENTIER) :: id_nb
    real(kind=DOUBLE) :: p_i, p_nb, vn_i, vn_nb, lam_sf
    real(kind=DOUBLE), dimension(5) :: w
    real(kind=DOUBLE), dimension(3) :: norm

    do_acou = .true.; do_vol = .true.
    if (present(only)) then
      do_acou = (only /= 'V'); do_vol = (only /= 'A')
    end if
    dt = 1e6
    dt_nan_seen = .false.; dt_nan_val = 0.0_DOUBLE
    if (diag_dtb) then
      dtb_acou = 1e300_DOUBLE; dtb_vol = 1e300_DOUBLE
      dtb_n_nan_acou = 0; dtb_n_nan_vol = 0; dtb_i_first_nan = -1
    end if
    do i=1, mesh%n_elems
      c = sqrt(max(0.0_DOUBLE, gamma_arr(i)*pressure(sol(:, i), gamma_arr(i))*sol(1, i)))
      c = max(c, 1.0e-10_DOUBLE)
      ! Lagrangian wave speed (acoustic impedance), homogeneous to rho*a.
      rho = 1.0_DOUBLE/max(sol(1, i), 1.0e-30_DOUBLE)
      lambda = rho*c
      p_i    = pressure(sol(:, i), gamma_arr(i))
      mass_c = rho*abs(mesh%elem(i)%volume)
      ! Accumulate over ALL sub-faces of the cell. The previous version took a min
      ! of avg_sub_vol/(A_pcf*c) over sub-faces, i.e. omega_c/MAX_f(A_pcf*c) with the
      ! sub-element volume replaced by its cell average. Both are non-conservative on
      ! an irregular mesh (they coincide with the sum only when every sub-face has the
      ! same area), and let a node be ejected on Voronoi cells while the run still
      ! ends with EXIT_CODE 0.
      sum_al = 0.0_DOUBLE
      sum_av = 0.0_DOUBLE
      do j=1, mesh%elem(i)%n_sub_elems
        id_sub_elem = mesh%elem(i)%sub_elem(j)
        do k=1, mesh%sub_elem(id_sub_elem)%n_sub_faces
          id_sub_face = mesh%sub_elem(id_sub_elem)%sub_face(k)
          id_face = mesh%sub_face(id_sub_face)%mesh_face
          id_vert = mesh%sub_face(id_sub_face)%mesh_vert
          if( mesh%face(id_face)%left_neigh == i ) then
            norm = mesh%face(id_face)%norm
            id_nb = mesh%face(id_face)%right_neigh
          else
            norm = -mesh%face(id_face)%norm
            id_nb = mesh%face(id_face)%left_neigh
          end if
          ! Au BORD on garde rho*a seul, comme la branche de bord du flux.
          if (lambda_acoustic_only .or. id_nb <= 0) then
            lam_sf = lambda
            ! Face de pression : meme lambda que l'etat fantome du flux (p_r = p_bc, v_r = v_l).
            if (.not. lambda_acoustic_only .and. id_nb < 0 .and. present(bc_p)) then
              if (bc_p(-id_nb) >= 0.0_DOUBLE) &
                lam_sf = max(lambda, sqrt(max(0.0_DOUBLE, bc_p(-id_nb) - p_i)/sol(1, i)))
            end if
            ! Mur symetrique : meme lambda que le flux (etat miroir) ; faces z du faux 2D exclues.
            if (.not. lambda_acoustic_only .and. id_nb < 0 .and. trim(wall_mode) == 'symmetric' &
                .and. abs(norm(3)) < 0.5_DOUBLE .and. present(bc_p)) then
              if (bc_p(-id_nb) < 0.0_DOUBLE) &
                lam_sf = max(lambda, 2.0_DOUBLE*dot_product(sol(2:4, i), norm)/sol(1, i))
            end if
          else
            p_nb   = pressure(sol(:, id_nb), gamma_arr(id_nb))
            vn_i   = dot_product(sol(2:4, i), norm)
            vn_nb  = dot_product(sol(2:4, id_nb), norm)
            lam_sf = max(lambda, sqrt(max(0.0_DOUBLE, p_nb-p_i)/sol(1, i)), &
              -(vn_nb-vn_i)/sol(1, i))
          end if
          sum_al = sum_al + mesh%sub_face(id_sub_face)%area*lam_sf
          sum_av = sum_av + mesh%sub_face(id_sub_face)%area &
            *abs(dot_product(vp(:, id_vert), norm))
        end do
      end do
      ! dt <= omega_c / sum_f (a_pcf * lambda_pcf), written with the mass so that the
      ! impedance form is explicit: m_c/sum(A*rho*a) == omega_c/sum(A*a).
      ! Compte les bornes NaN AVANT de les passer a min(), pour etablir la cause
      ! de la perte du minimum au lieu de la supposer.
      ! Compte CUMULATIF sur le run, hors diag_dtb : des qu'une borne est NaN,
      ! `dt = min(dt, borne)` cesse de rendre le minimum (min avec NaN rend le second
      ! operande), donc le run tourne avec un dt faux a partir de la. C'est le seul
      ! moyen de savoir par cas lesquels sont affectes, sans rebalayer les logs.
      if (mass_c/max(sum_al, 1.0e-30_DOUBLE) /= mass_c/max(sum_al, 1.0e-30_DOUBLE)) &
        dtb_nan_tot_acou = dtb_nan_tot_acou + 1
      if (CV*abs(mesh%elem(i)%volume)/(1.0e-8_DOUBLE + sum_av) &
          /= CV*abs(mesh%elem(i)%volume)/(1.0e-8_DOUBLE + sum_av)) &
        dtb_nan_tot_vol = dtb_nan_tot_vol + 1
      ! ORIGINE de la borne NaN. Cette boucle parcourt TOUTES les mailles, fantomes
      ! incluses, alors que le compteur de NaN de la trace les SAUTE
      ! (`if (is_ghost) cycle`). Une maille fantome a l'etat casse produit donc une
      ! borne NaN, defait le `min`, et reste invisible au depouillement. Mesure sur
      ! saltzman : au pas 223 l'interieur est sain -- eps_min 1.479e-06, NaN 0,
      ! max|v_p| exactement 1.0 -- et `dt` passe quand meme de 8.8e-14 a 0.728.
      ! On classe donc chaque borne NaN par sa source au lieu de la supposer.
      if (mass_c /= mass_c .or. sum_al /= sum_al .or. sum_av /= sum_av) then
        if (mesh%elem(i)%is_ghost) then
          dtb_nan_ghost = dtb_nan_ghost + 1
        else
          dtb_nan_reel = dtb_nan_reel + 1
        end if
        if (any(sol(:, i) /= sol(:, i))) dtb_nan_src_sol = dtb_nan_src_sol + 1
        if (sum_av /= sum_av) dtb_nan_src_vp = dtb_nan_src_vp + 1
        if (mesh%elem(i)%volume /= mesh%elem(i)%volume) &
          dtb_nan_src_vol = dtb_nan_src_vol + 1
      end if
      if (diag_dtb) then
        if (mass_c/max(sum_al, 1.0e-30_DOUBLE) /= mass_c/max(sum_al, 1.0e-30_DOUBLE)) &
          dtb_n_nan_acou = dtb_n_nan_acou + 1
        if (CV*abs(mesh%elem(i)%volume)/(1.0e-8_DOUBLE + sum_av) &
            /= CV*abs(mesh%elem(i)%volume)/(1.0e-8_DOUBLE + sum_av)) &
          dtb_n_nan_vol = dtb_n_nan_vol + 1
      end if
      ! min(dt, NaN) rend dt : une borne NaN etait donc SILENCIEUSEMENT ignoree et le
      ! run continuait avec un dt faux. On la retient au lieu de la perdre.
      if (do_acou) then
        bnd = mass_c/max(sum_al, 1.0e-30_DOUBLE)
        if (bnd /= bnd) then
          dt_nan_seen = .true.; dt_nan_val = bnd
        else
          dt = min(dt, bnd)
        end if
      end if
      if (do_vol) then
        bnd = CV*abs(mesh%elem(i)%volume)/(1.0e-8_DOUBLE + sum_av)
        if (bnd /= bnd) then
          dt_nan_seen = .true.; dt_nan_val = bnd
        else
          dt = min(dt, bnd)
        end if
      end if
      if (diag_dtb .and. dt /= dt .and. dtb_i_first_nan < 0) dtb_i_first_nan = i
      if (diag_dtb) then
        if (mass_c/max(sum_al, 1.0e-30_DOUBLE) < dtb_acou) then
          dtb_acou = mass_c/max(sum_al, 1.0e-30_DOUBLE); dtb_cell_acou = i
          ! Decomposition de la borne acoustique m_c / somme(A lambda). La masse est
          ! constante en lagrangien, donc si la borne s'effondre c'est que
          ! somme(A lambda) croit -- et il faut savoir si c'est lambda (donc la
          ! thermodynamique, donc la fermeture) ou somme(A) (donc la GEOMETRIE, des
          ! mailles qui s'etirent a volume constant). La dichotomie
          ! "borne volume = geometrie / borne acoustique = fermeture" ne suffit pas :
          ! la geometrie peut entrer par somme(A) a V constant.
          dtb_acou_m = mass_c
          dtb_acou_sal = sum_al
          dtb_acou_lam = lambda
          dtb_acou_area = sum_al/max(lambda, 1.0e-300_DOUBLE)
        end if
        if (CV*abs(mesh%elem(i)%volume)/(1.0e-8_DOUBLE + sum_av) < dtb_vol) then
          dtb_vol = CV*abs(mesh%elem(i)%volume)/(1.0e-8_DOUBLE + sum_av)
          dtb_cell_vol = i
          dtb_vol_V = abs(mesh%elem(i)%volume); dtb_vol_sumav = sum_av
        end if
      end if
    end do
    ! Instrumentation : `dt` sorti de la boucle ne correspondait pas au minimum
    ! rapporte par dtb_acou/dtb_vol (dt 1317.83 contre des bornes a 3.03e-05).
    ! Hypothese a tester : `min(dt, x)` est defait par un NaN -- les compteurs
    ! dtb utilisent un `<` explicite, FAUX pour NaN, donc eux gardent le vrai
    ! minimum, alors que `min` avec NaN peut rendre le second operande et donc
    ! REMETTRE dt a la valeur de la maille courante. On compte les bornes NaN et on
    ! recalcule le minimum avec un `<` explicite pour comparer.
    if (diag_dtb) then
      block
        integer(kind=ENTIER) :: ii, nnan
        real(kind=DOUBLE) :: ba, bv, dmin
        nnan = 0; dmin = 1e300_DOUBLE
        do ii = 1, mesh%n_elems
          ba = 0.0_DOUBLE; bv = 0.0_DOUBLE
          if (dtb_acou == dtb_acou) ba = dtb_acou
          if (dtb_vol == dtb_vol) bv = dtb_vol
        end do
        dmin = min(dtb_acou, dtb_vol)
        print*, "DTBUG dt de la boucle =", dt, " minimum des deux bornes =", dmin, &
                " rapport =", dt/max(dmin, 1.0e-300_DOUBLE), " cfl =", cfl, &
                " bornes NaN : acoustique", dtb_n_nan_acou, " volume", dtb_n_nan_vol, &
                " 1re maille ou dt devient NaN", dtb_i_first_nan, &
                " sur", mesh%n_elems
        if (dt /= dt) print*, "DTBUG dt est NaN en sortie de boucle"
        if (dt > 2.0_DOUBLE*dmin) print*, &
          "DTBUG ANOMALIE : dt retenu DEPASSE le minimum des bornes -- le min a ete perdu"
      end block
    end if
    dt = cfl * dt

    if( num_procs > 1 ) call MPI_ALLREDUCE(MPI_IN_PLACE, &
      dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, mpi_ierr)

    ! MPI_MIN ne transporte pas un NaN de facon fiable : on reduit le DRAPEAU, pour que
    ! tous les rangs le sachent et s arretent ENSEMBLE plutot qu un seul.
    inan = 0
    if (dt_nan_seen) inan = 1
    if( num_procs > 1 ) call MPI_ALLREDUCE(MPI_IN_PLACE, inan, 1, MPI_INTEGER, &
      MPI_MAX, MPI_COMM_WORLD, mpi_ierr)
    dt_nan_detected = (inan > 0)
    if (dt_nan_detected) then
      if (dt_nan_seen) then
        dt = dt_nan_val
      else
        dt = ieee_value(dt, ieee_quiet_nan)
      end if
    end if
  end subroutine compute_dt

  subroutine init_sol(mesh, sol, sol_uniform, init, me, num_procs, gamma_arr, boundary_2d)
    use mpi
    use mpi_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: sol
    real(kind=DOUBLE) :: r_gre, p_gre, p0_gre
    real(kind=DOUBLE) :: p_inf_ac, a_inf_ac, rho_ac, u_ac, rho_tot_ac, x_ac
    real(kind=DOUBLE), parameter :: pi_const = 4.0_DOUBLE*atan(1.0_DOUBLE)
    real(kind=DOUBLE), dimension(2) :: v_gre
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_uniform
    integer(kind=ENTIER), intent(in) :: init
    integer(kind=ENTIER), intent(in) :: me, num_procs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(inout) :: gamma_arr
    logical, intent(in) :: boundary_2d

    integer(kind=ENTIER) :: i, imin, mpi_ierr
    real(kind=DOUBLE) :: tot_vol, r, rmin, rmin_share, e_dep, vol_dep, rv, rvt, rv_share
    integer(kind=ENTIER) :: iv0, j, k
    logical :: node_at_origin
    !real(kind=DOUBLE), parameter :: r_sedov = 2.4_DOUBLE/20.0_DOUBLE
    real(kind=DOUBLE), parameter :: r_sedov = 0.011_DOUBLE
    real(kind=DOUBLE), parameter :: r_sedov_3d = 0.12_DOUBLE

    if ( init == 0 ) then
      do i=1, mesh%n_elems
        if( mesh%elem(i)%coord(1) < 0.5_DOUBLE ) then
          sol(1, i) = 1.0_DOUBLE
          sol(2:4, i) = 0.0_DOUBLE
          sol(5, i) = 1.0_DOUBLE
        else
          sol(1, i) = 0.125_DOUBLE
          sol(2:4, i) = 0.0_DOUBLE
          sol(5, i) = 0.1_DOUBLE
        end if
        !sol(:, i) = sol_uniform
        sol(1, i) = 1.0_DOUBLE/sol(1, i)
        sol(5, i) = sol(5, i)/((gamma_arr(i) - 1.0_DOUBLE)/sol(1, i)) &
          + 0.5_DOUBLE*norm2(sol(2:4, i))**2
      end do
    else if( init == 1) then
      ! Find the cell nearest to the origin (same as Euler init_sedov).
      rmin = huge(1.0_DOUBLE)
      imin = -1
      do i = 1, mesh%n_elems
        if (mesh%elem(i)%is_ghost) cycle
        if (boundary_2d) then
          if (norm2(mesh%elem(i)%coord(:2)) < rmin) then
            rmin = norm2(mesh%elem(i)%coord(:2))
            imin = i
          end if
        else
          if (norm2(mesh%elem(i)%coord) < rmin) then
            rmin = norm2(mesh%elem(i)%coord)
            imin = i
          end if
        end if
      end do


      ! Rayon minimal GLOBAL. Une seule reduction, et en MPI_DOUBLE_PRECISION : la
      ! version precedente passait MPI_DOUBLE, qui est la constante du binding C et
      ! n'est pas le type Fortran attendu ici.
      rmin_share = rmin
      if (num_procs > 1) then
        call MPI_ALLREDUCE(MPI_IN_PLACE, rmin_share, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
          MPI_COMM_WORLD, mpi_ierr)
        if (abs(rmin - rmin_share) > 1e-14_DOUBLE) imin = -1
      end if

      do i = 1, mesh%n_elems
        sol(1, i) = 1.0_DOUBLE
        sol(2:4, i) = 0.0_DOUBLE
        sol(5, i) = 1e-12_DOUBLE/(gamma_arr(i) - 1.0_DOUBLE)
      end do

      ! energie totale a deposer
      e_dep = sedov_energy
      if (e_dep <= 0.0_DOUBLE) then
        e_dep = merge(0.984042_DOUBLE, 0.851072_DOUBLE, boundary_2d)
      end if

      ! Depot nodal : si un sommet du maillage tombe sur l'origine, l'energie va dans
      ! TOUTES les cellules qui le touchent, au prorata du volume.
      iv0 = -1
      node_at_origin = .false.
      if (sedov_nodal_deposit) then
        ! variable dediee : ne PAS reutiliser rmin_share ici, il porte le rayon
        ! minimal des CENTROIDES et sert au test d'egalite plus bas.
        rv = huge(1.0_DOUBLE)
        do i = 1, mesh%n_vert
          if (boundary_2d) then
            rvt = norm2(mesh%vert(i)%coord(:2))
          else
            rvt = norm2(mesh%vert(i)%coord)
          end if
          if (rvt < rv) then
            rv = rvt
            iv0 = i
          end if
        end do
        ! Le test "un sommet est-il sur l'origine ?" doit etre GLOBAL. En local il
        ! est faux sur tous les rangs qui ne possedent pas l'origine : ceux-la
        ! tombaient alors dans la branche centroide et redeposaient de l'energie
        ! sur leur propre minimum local, tout en contribuant a un ALLREDUCE que
        ! le rang de l'origine faisait sur une autre quantite.
        rv_share = rv
        if (num_procs > 1) then
          call MPI_ALLREDUCE(MPI_IN_PLACE, rv_share, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
            MPI_COMM_WORLD, mpi_ierr)
        end if
        node_at_origin = (rv_share <= 1.0e-10_DOUBLE)
        ! iv0 ne reste positif que sur les rangs qui touchent reellement l'origine
        ! (un sommet d'interface est duplique : chacun sommera ses mailles non-ghost).
        if (.not. node_at_origin .or. rv > 1.0e-10_DOUBLE) iv0 = -1
      end if

      if (node_at_origin) then
        vol_dep = 0.0_DOUBLE
        if (iv0 > 0) then
          do j = 1, mesh%vert(iv0)%n_elems_neigh
            k = mesh%vert(iv0)%elem_neigh(j)
            if (k > 0) then
              if (.not. mesh%elem(k)%is_ghost) vol_dep = vol_dep + mesh%elem(k)%volume
            end if
          end do
        end if
        ! collective prise par TOUS les rangs : node_at_origin est global.
        if (num_procs > 1) then
          call MPI_ALLREDUCE(MPI_IN_PLACE, vol_dep, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
            MPI_COMM_WORLD, mpi_ierr)
        end if
        if (iv0 > 0 .and. vol_dep > 0.0_DOUBLE) then
          do j = 1, mesh%vert(iv0)%n_elems_neigh
            k = mesh%vert(iv0)%elem_neigh(j)
            if (k > 0) then
              if (.not. mesh%elem(k)%is_ghost) sol(5, k) = e_dep/vol_dep
            end if
          end do
        end if
      else
        ! Pas de sommet a l'origine : deposer dans TOUTES les cellules a egale distance
        ! minimale, pas dans une seule. Le cas se produit des que l'origine tombe sur une
        ! arete ou une face -- sur le nid d'abeilles par exemple, deux centroides sont a
        ! 0.02143 exactement, et l'ancien code tranchait l'egalite par l'ordre de boucle,
        ! ce qui deposait toute l'energie d'un seul cote et brisait la symetrie des le
        ! premier pas. La tolerance est relative au rayon minimal lui-meme.
        ! Ce partage ne depend PAS de sedov_nodal_deposit, qui ne concerne que le cas
        ! d un SOMMET sur l origine : le verrou etait mal place. Sans lui, vol_dep restait
        ! nul au defaut et chaque rang a egalite deposait e_dep en entier, soit k e_dep.
        vol_dep = 0.0_DOUBLE
        do i = 1, mesh%n_elems
          if (mesh%elem(i)%is_ghost) cycle
          if (boundary_2d) then
            r = norm2(mesh%elem(i)%coord(:2))
          else
            r = norm2(mesh%elem(i)%coord)
          end if
          if (r <= rmin_share*(1.0_DOUBLE + 1.0e-8_DOUBLE) + 1.0e-14_DOUBLE) &
            vol_dep = vol_dep + mesh%elem(i)%volume
        end do
        if (num_procs > 1) then
          call MPI_ALLREDUCE(MPI_IN_PLACE, vol_dep, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
            MPI_COMM_WORLD, mpi_ierr)
        end if
        if (vol_dep > 0.0_DOUBLE) then
          do i = 1, mesh%n_elems
            if (mesh%elem(i)%is_ghost) cycle
            if (boundary_2d) then
              r = norm2(mesh%elem(i)%coord(:2))
            else
              r = norm2(mesh%elem(i)%coord)
            end if
            if (r <= rmin_share*(1.0_DOUBLE + 1.0e-8_DOUBLE) + 1.0e-14_DOUBLE) &
              sol(5, i) = e_dep/vol_dep
          end do
        else if (imin > 0) then
          sol(5, imin) = e_dep/mesh%elem(imin)%volume
        end if
      end if
    else if( init == 13) then

      tot_vol = 0.0_DOUBLE
      do i = 1, mesh%n_elems
        if( .not. mesh%elem(i)%is_ghost ) then
          if (norm2(mesh%elem(i)%coord) < r_sedov_3d) then
            tot_vol = tot_vol + mesh%elem(i)%volume
          end if
        end if
      end do

      if( num_procs > 1 ) call MPI_ALLREDUCE(MPI_IN_PLACE, &
        tot_vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)

      if (tot_vol < 1e-10_DOUBLE) then
        print *, "[-] r_sedov_3d too small, tot_vol_init:", tot_vol
        error stop
      end if

      do i = 1, mesh%n_elems
        if (norm2(mesh%elem(i)%coord) < r_sedov_3d) then
          sol(1, i) = 1.0_DOUBLE
          sol(2, i) = 0.0_DOUBLE
          sol(3, i) = 0.0_DOUBLE
          sol(4, i) = 0.0_DOUBLE
          !sol(5,i) = 0.311357_DOUBLE / tot_vol
          !sol(5, i) = 0.244816_DOUBLE/tot_vol
          !sol(5, i) = 0.851072_DOUBLE/tot_vol
          sol(5, i) = 0.983909_DOUBLE/tot_vol
        else
          sol(1, i) = 1.0_DOUBLE
          sol(2, i) = 0.0_DOUBLE
          sol(3, i) = 0.0_DOUBLE
          sol(4, i) = 0.0_DOUBLE
          sol(5, i) = 1e-8_DOUBLE/(gamma_arr(i) - 1.0_DOUBLE)
        end if
      end do
    else if( init == 2) then
      do i=1, mesh%n_elems
        call sol_isentropic_vortex(mesh%elem(i)%coord, sol(:, i), 0.0_DOUBLE)
        sol(5, i) = sol(5, i)/(sol(1, i)*(gamma_arr(i)-1.0_DOUBLE)) &
          + 0.5_DOUBLE*norm2(sol(2:4, i))**2
        sol(1, i) = 1.0_DOUBLE/sol(1, i)
      end do
    else if( init == 10) then
      ! Taylor-Green. Initialised with the exact CELL AVERAGE of the conserved
      ! variables (quadrature over the cell), not the centroid point value as
      ! init=2 does: a centroid initialisation injects an O(h^2) error of its own,
      ! which is the same order as what a second-order scheme is being measured
      ! against. rho is uniform here, so the mass average is the volume average.
      do i=1, mesh%n_elems
        call taylor_green_cell_average(mesh, i, sol(:, i))
      end do
    else if( init == 11 ) then
      ! Gresho vortex, same profile as the Euler solver's sol_gresho_mach:
      ! rho = 1, purely azimuthal velocity v_theta = 5r (r<0.2), 2-5r
      ! (0.2<r<0.4), 0 beyond, and the pressure that balances it. It is a
      ! STEADY solution, so any departure from the initial state is scheme
      ! error -- which is what makes it the low-Mach test.
      ! gresho_mach sets the pressure level p0 = 1/(gamma*Ma^2); the smaller
      ! it is, the harder the low-Mach regime.
      do i=1, mesh%n_elems
        r_gre = norm2(mesh%elem(i)%coord(:2))
        p0_gre = 1.0_DOUBLE/(gamma_arr(i)*gresho_mach**2)
        if (r_gre < 0.2_DOUBLE) then
          v_gre(1) = -5.0_DOUBLE*mesh%elem(i)%coord(2)
          v_gre(2) =  5.0_DOUBLE*mesh%elem(i)%coord(1)
          p_gre = p0_gre + 12.5_DOUBLE*r_gre**2
        else if (r_gre < 0.4_DOUBLE) then
          v_gre(1) = (5.0_DOUBLE - 2.0_DOUBLE/r_gre)*mesh%elem(i)%coord(2)
          v_gre(2) = (-5.0_DOUBLE + 2.0_DOUBLE/r_gre)*mesh%elem(i)%coord(1)
          p_gre = p0_gre + 12.5_DOUBLE*r_gre**2 &
            + 4.0_DOUBLE - 20.0_DOUBLE*r_gre + 4.0_DOUBLE*log(5.0_DOUBLE*r_gre)
        else
          v_gre(1) = 0.0_DOUBLE
          v_gre(2) = 0.0_DOUBLE
          p_gre = p0_gre - 2.0_DOUBLE + 4.0_DOUBLE*log(2.0_DOUBLE)
        end if
        sol(1, i) = 1.0_DOUBLE            ! specific volume, rho = 1
        sol(2, i) = v_gre(1)
        sol(3, i) = v_gre(2)
        sol(4, i) = 0.0_DOUBLE
        sol(5, i) = p_gre/(gamma_arr(i) - 1.0_DOUBLE) &
          + 0.5_DOUBLE*(v_gre(1)**2 + v_gre(2)**2)
      end do
    else if( init == 14 ) then
      ! Onde acoustique plane traversant un vortex de Gresho.
      ! NOTE : numerote 14 et non 13 -- init = 13 est DEJA pris plus haut par un
      ! depot de type Sedov, et une seconde branche 13 etait inatteignable.
      ! Le vortex est centre a l'origine et de rayon 0.4 ; le paquet acoustique
      ! occupe -0.8 < x < -0.6, donc les deux supports sont DISJOINTS a t = 0 et
      ! la superposition est sans ambiguite : la ou l'onde est, l'ecoulement de
      ! base est l'etat uniforme du champ lointain.
      do i = 1, mesh%n_elems
        p0_gre = 1.0_DOUBLE/(gamma_arr(i)*gresho_mach**2)
        p_inf_ac = p0_gre - 2.0_DOUBLE + 4.0_DOUBLE*log(2.0_DOUBLE)
        a_inf_ac = sqrt(gamma_arr(i)*p_inf_ac)          ! rho_inf = 1

        ! --- vortex de Gresho, identique a init = 11
        r_gre = norm2(mesh%elem(i)%coord(:2))
        if (vac_vortex .and. r_gre < 0.2_DOUBLE) then
          v_gre(1) = -5.0_DOUBLE*mesh%elem(i)%coord(2)
          v_gre(2) =  5.0_DOUBLE*mesh%elem(i)%coord(1)
          p_gre = p0_gre + 12.5_DOUBLE*r_gre**2
        else if (vac_vortex .and. r_gre < 0.4_DOUBLE) then
          v_gre(1) = (5.0_DOUBLE - 2.0_DOUBLE/r_gre)*mesh%elem(i)%coord(2)
          v_gre(2) = (-5.0_DOUBLE + 2.0_DOUBLE/r_gre)*mesh%elem(i)%coord(1)
          p_gre = p0_gre + 12.5_DOUBLE*r_gre**2 &
            + 4.0_DOUBLE - 20.0_DOUBLE*r_gre + 4.0_DOUBLE*log(5.0_DOUBLE*r_gre)
        else
          v_gre(1) = 0.0_DOUBLE
          v_gre(2) = 0.0_DOUBLE
          p_gre = p_inf_ac
        end if

        ! --- paquet acoustique, demi-cosinus nul a ses deux bords
        x_ac = mesh%elem(i)%coord(1)
        rho_ac = 0.0_DOUBLE
        if (vac_wave .and. x_ac > -0.8_DOUBLE .and. x_ac < -0.6_DOUBLE) then
          rho_ac = gresho_mach*cos(5.0_DOUBLE*pi_const*(x_ac + 0.7_DOUBLE))
        end if
        rho_tot_ac = 1.0_DOUBLE + rho_ac
        ! invariant de Riemann d'une onde simple : a(rho) = a_inf rho^((g-1)/2)
        u_ac = 2.0_DOUBLE/(gamma_arr(i) - 1.0_DOUBLE) &
             * a_inf_ac*(rho_tot_ac**(0.5_DOUBLE*(gamma_arr(i) - 1.0_DOUBLE)) - 1.0_DOUBLE)
        if (rho_ac /= 0.0_DOUBLE) p_gre = p_inf_ac*rho_tot_ac**gamma_arr(i)

        sol(1, i) = 1.0_DOUBLE/rho_tot_ac          ! volume specifique
        sol(2, i) = v_gre(1) + u_ac
        sol(3, i) = v_gre(2)
        sol(4, i) = 0.0_DOUBLE
        sol(5, i) = p_gre*sol(1, i)/(gamma_arr(i) - 1.0_DOUBLE) &
          + 0.5_DOUBLE*(sol(2, i)**2 + sol(3, i)**2)
      end do
    else if( init == 42 ) then
      do i=1, mesh%n_elems
        sol(1, i) = 1.0_DOUBLE
        !sol(2:4, i) = 0.0_DOUBLE
        sol(2:4, i) = 0.0_DOUBLE
        sol(2, i) = mesh%elem(i)%coord(1)
        !sol(5, i) = (10+mesh%elem(i)%coord(1))/(gamma-1.0_DOUBLE)
        sol(5, i) = 10/(gamma_arr(i)-1.0_DOUBLE)
      end do
    else if( init == 4 ) then
      do i=1, mesh%n_elems
        if( norm2(mesh%elem(i)%coord(:2)) < 0.5_DOUBLE ) then
          sol(1, i) = 1.0_DOUBLE
          sol(2:4, i) = 0.0_DOUBLE
          sol(5, i) = 1.0_DOUBLE
        else
          sol(1, i) = 0.125_DOUBLE
          sol(2:4, i) = 0.0_DOUBLE
          sol(5, i) = 0.1_DOUBLE
        end if
        sol(1, i) = 1.0_DOUBLE/sol(1, i)
        sol(5, i) = sol(5, i)/((gamma_arr(i) - 1.0_DOUBLE)/sol(1, i)) &
          + 0.5_DOUBLE*norm2(sol(2:4, i))**2
      end do
    else if( init == 12 ) then
      ! Etat uniforme au repos, pris tel quel dans sol_uniform = (tau, u, v, w, E).
      ! Sert au test d'odd-even decoupling : le gaz est au repos et c'est le piston
      ! qui cree le choc, donc aucune structure n'est imposee a l'initialisation.
      do i=1, mesh%n_elems
        sol(:, i) = sol_uniform
      end do
    else if( init == 3 ) then
      ! Noh 2D cylindrical: uniform density rho=1, radial inward unit velocity, p~0
      do i=1, mesh%n_elems
        sol(1, i) = 1.0_DOUBLE
        r = sqrt(mesh%elem(i)%coord(1)**2 + mesh%elem(i)%coord(2)**2)
        if (r < 1.0e-10_DOUBLE) then
          sol(2, i) = 0.0_DOUBLE
          sol(3, i) = 0.0_DOUBLE
        else
          sol(2, i) = -mesh%elem(i)%coord(1) / r
          sol(3, i) = -mesh%elem(i)%coord(2) / r
        end if
        sol(4, i) = 0.0_DOUBLE
        ! E = p/((gamma-1)*rho) + 0.5*|v|^2, with p=1e-6 and |v|=1 away from center
        sol(5, i) = 1.0e-6_DOUBLE / (gamma_arr(i) - 1.0_DOUBLE) &
          + 0.5_DOUBLE * (sol(2, i)**2 + sol(3, i)**2)
      end do
    else if( init == 5 ) then
      ! Saltzman piston: uniform gas at rest (rho=1, p~0); piston BC at left wall sets v=(1,0,0)
      do i=1, mesh%n_elems
        sol(1, i) = 1.0_DOUBLE
        sol(2:4, i) = 0.0_DOUBLE
        sol(5, i) = 1.0e-6_DOUBLE / (gamma_arr(i) - 1.0_DOUBLE)
      end do
    else if( init == 6 ) then
      ! Sod bi-materiaux: left (x<0.5) gamma=1.4, right (x>=0.5) gamma=1.5
      do i=1, mesh%n_elems
        if( mesh%elem(i)%coord(1) < 0.5_DOUBLE ) then
          gamma_arr(i) = 1.4_DOUBLE
          sol(1, i) = 1.0_DOUBLE
          sol(2:4, i) = 0.0_DOUBLE
          sol(5, i) = 1.0_DOUBLE
        else
          gamma_arr(i) = 1.5_DOUBLE
          sol(1, i) = 0.125_DOUBLE
          sol(2:4, i) = 0.0_DOUBLE
          sol(5, i) = 0.1_DOUBLE
        end if
        sol(1, i) = 1.0_DOUBLE/sol(1, i)
        sol(5, i) = sol(5, i)/((gamma_arr(i) - 1.0_DOUBLE)/sol(1, i)) &
          + 0.5_DOUBLE*norm2(sol(2:4, i))**2
      end do
    else if( init == 7 ) then
      ! Triple point 3D (LLNL version): three states, two materials
      ! Left  (x<1):                     rho=1,   p=1,   gamma=1.5
      ! Right, outside cyl (r_yz>1.5):  rho=0.1, p=0.1, gamma=1.5
      ! Right, inside  cyl (r_yz<=1.5): rho=1,   p=0.1, gamma=1.4
      ! where r_yz = sqrt(y^2 + z^2) is the radius from the x-axis
      do i=1, mesh%n_elems
        r = sqrt(mesh%elem(i)%coord(2)**2 + mesh%elem(i)%coord(3)**2)
        if (mesh%elem(i)%coord(1) < 1.0_DOUBLE) then
          gamma_arr(i) = 1.5_DOUBLE
          sol(1, i)    = 1.0_DOUBLE
          sol(5, i)    = 1.0_DOUBLE
        else if (r > 1.5_DOUBLE) then
          gamma_arr(i) = 1.5_DOUBLE
          sol(1, i)    = 0.1_DOUBLE
          sol(5, i)    = 0.1_DOUBLE
        else
          gamma_arr(i) = 1.4_DOUBLE
          sol(1, i)    = 1.0_DOUBLE
          sol(5, i)    = 0.1_DOUBLE
        end if
        sol(2:4, i) = 0.0_DOUBLE
        sol(1, i) = 1.0_DOUBLE / sol(1, i)
        sol(5, i) = sol(5, i) / ((gamma_arr(i) - 1.0_DOUBLE) / sol(1, i)) &
          + 0.5_DOUBLE * norm2(sol(2:4, i))**2
      end do
    else if( init == 8 ) then
      ! Gradient test: linear pressure p = 1 + x + y, zero velocity, rho = 1.
      ! Exact grad_p = [1, 1, 0]; exact grad_v = 0; exact div_v = 0.
      do i=1, mesh%n_elems
        sol(1, i) = 1.0_DOUBLE
        sol(2:4, i) = 0.0_DOUBLE
        sol(5, i) = (1.0_DOUBLE + mesh%elem(i)%coord(1) + mesh%elem(i)%coord(2)) &
                    / (gamma_arr(i) - 1.0_DOUBLE)
      end do
    else if( init == 18 ) then
      ! Test de gradient : p = 1 + x + y et v = (0,3x + 0,2y, -0,1x + 0,4y, 0), rho = 1.
      ! grad_p = (1,1,0), grad_vx = (0,3;0,2;0), grad_vy = (-0,1;0,4;0), div_v = 0,7.
      do i=1, mesh%n_elems
        sol(1, i) = 1.0_DOUBLE
        sol(2, i) = 0.3_DOUBLE*mesh%elem(i)%coord(1) + 0.2_DOUBLE*mesh%elem(i)%coord(2)
        sol(3, i) = -0.1_DOUBLE*mesh%elem(i)%coord(1) + 0.4_DOUBLE*mesh%elem(i)%coord(2)
        sol(4, i) = 0.0_DOUBLE
        sol(5, i) = (1.0_DOUBLE + mesh%elem(i)%coord(1) + mesh%elem(i)%coord(2)) &
                    / (gamma_arr(i) - 1.0_DOUBLE) + 0.5_DOUBLE*dot_product(sol(2:4, i), sol(2:4, i))
      end do
    else if( init == 9 ) then
      ! Shock-bubble interaction (Haas & Sturtevant 1987, He+28% air bubble),
      ! Ms=1.22 planar shock in air. Geometry/gas data from Razmi et al.
      ! (2019, JAFM 12(2), Fig. 2 / Table 1) and Quirk & Karni (1996, JFM
      ! 318, Fig. 6). Full domain (no symmetry reduction), full circular
      ! bubble centered at mid-height (xc_bub, yc_bub), radius r_bub.
      ! The whole domain starts quiescent (rho=p=1, gamma=1.4) except the
      ! bubble, which is in mechanical/thermal equilibrium with the ambient
      ! air (same p, T), giving density R_air/R_bubble and gamma=1.648
      ! (Table 1 of the reference: He contaminated 28% by mass with air).
      ! The Ms=1.22 shock is generated by a 'piston' BC at the right wall
      ! (x=Lx) imposing the post-shock velocity u2 there (see input_data.f);
      ! this reproduces the reference's inflow BC on side CD ("exact flow
      ! conditions behind the incident shock wave") and matches the normal-
      ! shock relations for gamma=1.4: rho2/rho1=1.376417, p2/p1=1.569800,
      ! u2=-(2/(gamma+1))*(Ms-1/Ms)*sqrt(gamma*p1/rho1)=-0.394664.
      block
        real(kind=DOUBLE), parameter :: xc_bub = 0.350_DOUBLE
        real(kind=DOUBLE), parameter :: yc_bub = 0.0445_DOUBLE
        real(kind=DOUBLE), parameter :: r_bub = 0.025_DOUBLE
        real(kind=DOUBLE), parameter :: gamma_air = 1.4_DOUBLE
        real(kind=DOUBLE), parameter :: gamma_bub = 1.648_DOUBLE
        real(kind=DOUBLE), parameter :: rho_bub = 0.287_DOUBLE / 1.578_DOUBLE
        real(kind=DOUBLE) :: rb
        do i = 1, mesh%n_elems
          rb = sqrt((mesh%elem(i)%coord(1) - xc_bub)**2 + (mesh%elem(i)%coord(2) - yc_bub)**2)
          if (rb <= r_bub) then
            gamma_arr(i) = gamma_bub
            sol(1, i) = rho_bub
            sol(2, i) = 0.0_DOUBLE
            sol(5, i) = 1.0_DOUBLE
          else
            gamma_arr(i) = gamma_air
            sol(1, i) = 1.0_DOUBLE
            sol(2, i) = 0.0_DOUBLE
            sol(5, i) = 1.0_DOUBLE
          end if
          sol(3:4, i) = 0.0_DOUBLE
          sol(1, i) = 1.0_DOUBLE / sol(1, i)
          sol(5, i) = sol(5, i) / ((gamma_arr(i) - 1.0_DOUBLE) / sol(1, i)) &
            + 0.5_DOUBLE * norm2(sol(2:4, i))**2
        end do
      end block
    end if

  end subroutine init_sol

  function pressure(u, g)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: u
    real(kind=DOUBLE), intent(in) :: g
    real(kind=DOUBLE) :: pressure

    pressure = (g-1.0_DOUBLE)/u(1)*(u(5) - 0.5_DOUBLE*norm2(u(2:4))**2)
  end function pressure

  subroutine mpi_memory_exchange_vert(mesh, mpi_send_recv)
    use mpi
    implicit none

    type(mesh_type), intent(inout) :: mesh
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    integer :: mpi_ierr
    integer(kind=ENTIER) :: i, k, j
    integer(kind=ENTIER) :: id_elem, n_vert_tot, id_vert

    do i = 1, mpi_send_recv%n_mpi_send_neigh

      n_vert_tot = 0
      do k = 1, mpi_send_recv%mpi_send_neigh(i)%n_elems
        id_elem = mpi_send_recv%mpi_send_neigh(i)%elem_id(k)
        n_vert_tot = n_vert_tot + mesh%elem(id_elem)%n_vert
      end do

      allocate(mpi_send_recv%mpi_send_neigh(i)%sol(3, n_vert_tot))

      n_vert_tot = 1
      do k = 1, mpi_send_recv%mpi_send_neigh(i)%n_elems
        id_elem = mpi_send_recv%mpi_send_neigh(i)%elem_id(k)
        do j=1, mesh%elem(id_elem)%n_vert
          id_vert = mesh%elem(id_elem)%vert(j)
          mpi_send_recv%mpi_send_neigh(i)%sol(:, n_vert_tot) = mesh%vert(id_vert)%coord
          n_vert_tot = n_vert_tot + 1
        end do
      end do

      call mpi_isend(mpi_send_recv%mpi_send_neigh(i)%sol(1, 1), &
        3*(n_vert_tot-1), MPI_DOUBLE, &
        mpi_send_recv%mpi_send_neigh(i)%partition_id, &
        0, MPI_COMM_WORLD, mpi_send_recv%mpi_reqsend(i), mpi_ierr)
    end do

    do i = 1, mpi_send_recv%n_mpi_recv_neigh
      n_vert_tot = 0
      do k = 1, mpi_send_recv%mpi_recv_neigh(i)%n_elems
        id_elem = mpi_send_recv%mpi_recv_neigh(i)%elem_id(k)
        n_vert_tot = n_vert_tot + mesh%elem(id_elem)%n_vert
      end do
      allocate(mpi_send_recv%mpi_recv_neigh(i)%sol(3, n_vert_tot))
      call mpi_irecv(mpi_send_recv%mpi_recv_neigh(i)%sol(1, 1), &
        3*n_vert_tot, MPI_DOUBLE, &
        mpi_send_recv%mpi_recv_neigh(i)%partition_id, &
        MPI_ANY_TAG, MPI_COMM_WORLD, mpi_send_recv%mpi_reqrecv(i), mpi_ierr)
    end do

    call mpi_waitall(mpi_send_recv%n_mpi_send_neigh, &
      mpi_send_recv%mpi_reqsend, mpi_send_recv%mpi_sendstat, mpi_ierr)
    call mpi_waitall(mpi_send_recv%n_mpi_recv_neigh, &
      mpi_send_recv%mpi_reqrecv, mpi_send_recv%mpi_recvstat, mpi_ierr)

    do i = 1, mpi_send_recv%n_mpi_recv_neigh
      n_vert_tot = 1
      do k = 1, mpi_send_recv%mpi_recv_neigh(i)%n_elems
        id_elem = mpi_send_recv%mpi_recv_neigh(i)%elem_id(k)
        do j = 1, mesh%elem(id_elem)%n_vert
          id_vert = mesh%elem(id_elem)%vert(j)
          mesh%vert(id_vert)%coord = mpi_send_recv%mpi_recv_neigh(i)%sol(:, n_vert_tot)
          n_vert_tot = n_vert_tot + 1
        end do
      end do
    end do

    do i = 1, mpi_send_recv%n_mpi_send_neigh
      deallocate(mpi_send_recv%mpi_send_neigh(i)%sol)
    end do
    do i = 1, mpi_send_recv%n_mpi_recv_neigh
      deallocate(mpi_send_recv%mpi_recv_neigh(i)%sol)
    end do
  end subroutine mpi_memory_exchange_vert

  pure subroutine sol_isentropic_vortex(coord, w, t)
    use lagrange_global_data_module
    implicit none

    real(kind=DOUBLE), intent(in) :: t
    real(kind=DOUBLE), dimension(3), intent(in) :: coord
    real(kind=DOUBLE), dimension(5), intent(inout) :: w

    real(kind=DOUBLE), dimension(3) :: center_coord, vel
    real(kind=DOUBLE) :: beta, r

    beta = 5.0_DOUBLE
    vel(:) = (/0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE/)
    center_coord(:) = (/0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE/) + vel*t
    r = (coord(1) - center_coord(1))**2 + (coord(2) - center_coord(2))**2
    w(1) = 1.0_DOUBLE*(1.0_DOUBLE - ((gamma - 1)*beta**2)/(8.0_DOUBLE*gamma*pi**2)* &
      exp(1.0_DOUBLE - r))**(1.0_DOUBLE/(gamma - 1.0_DOUBLE))
    w(2) = vel(1) &
      - (coord(2) - center_coord(2))*beta/(2.0_DOUBLE*pi)*exp(0.5_DOUBLE*(1.0_DOUBLE - r))
    w(3) = vel(2) &
      + (coord(1) - center_coord(1))*beta/(2.0_DOUBLE*pi)*exp(0.5_DOUBLE*(1.0_DOUBLE - r))
    w(4) = 0.0_DOUBLE
    w(5) = w(1)**gamma
  end subroutine sol_isentropic_vortex

  subroutine compute_nodal_velocity_sidil(mesh, id_vert, sol, vp, method, b2d, b2d_h, gamma_arr, &
      n_bc, bc_type, bc_val)
    use lagrange_global_data_module, only : boundary_2d
    use linear_solver_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert, method
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3), intent(inout) :: vp
    real(kind=DOUBLE), intent(in) :: b2d_h
    logical, intent(in) :: b2d
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    integer(kind=ENTIER), intent(in), optional :: n_bc
    character(len=255), dimension(:), intent(in), optional :: bc_type
    real(kind=DOUBLE), dimension(:, :), intent(in), optional :: bc_val

    integer(kind=ENTIER) :: j, k, le, re
    integer(kind=ENTIER) :: id_sub_elem, id_elem
    integer(kind=ENTIER) :: id_sub_face, id_face
    real(kind=DOUBLE) :: rho_p, a_p, h_p
    real(kind=DOUBLE), dimension(3) :: grad_p, Bp, v_target
    real(kind=DOUBLE), dimension(5) :: sol_w, sol_l, sol_r
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(3,3) :: mat

    vp = 0.0_DOUBLE
    rho_p = 0.0_DOUBLE
    a_p = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_w = lag_to_primit(sol(:, id_elem), gamma_arr(id_elem))
      if( sol_w(5) < 0.0_DOUBLE ) then
        print*,"Error pos", id_elem, mesh%elem(id_elem)%coord
      end if
      vp = vp &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(2:4)
      rho_p = rho_p &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(1)
      a_p = a_p &
        + mesh%sub_elem(id_sub_elem)%volume &
        * sqrt(gamma_arr(id_elem)*sol_w(5)/sol_w(1))
    end do
    rho_p = rho_p / mesh%vert(id_vert)%volume
    a_p = a_p / mesh%vert(id_vert)%volume
    vp = vp / mesh%vert(id_vert)%volume

    grad_p = 0.0_DOUBLE
    mat = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      ! In boundary_2d, z-normal "cap" faces (top/bottom of the single
      ! prism/hex layer) carry no physical wall, exactly like in
      ! boundary_normal/corner_normal above -- but unlike those, this loop
      ! never excluded them. If the fluid ever picks up even a tiny
      ! out-of-plane z-velocity (numerical noise, since nothing in this
      ! scheme forces it to stay exactly 0 the way vp_z is forced to for
      ! the grid velocity), the cap-face mirror flips its sign, and the
      ! resulting jump gets weighted by the cap's *full 2D footprint area*
      ! (large, ~1e-4) instead of a lateral face's thin sliver (~1e-7) --
      ! turning a negligible z-velocity into a wildly amplified spurious
      ! grad_p/div_v contribution. Confirmed empirically: an interior
      ! vertex with genuinely uniform, at-rest neighboring cells still blew
      ! up to div_v~500 purely from this effect.
      if (b2d .and. abs(mesh%face(id_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE) cycle
      le = mesh%face(id_face)%left_neigh
      re = mesh%face(id_face)%right_neigh
      sol_l = sol(:, le)
      if( re > 0 ) then
        sol_r = sol(:, re)
      else
        ! Generic boundary mirror targets zero normal velocity (a static
        ! wall). A piston BC moves with an imposed velocity instead: its
        ! mirror must target that velocity's normal component, not zero,
        ! or the reconstructed boundary state (and hence grad_p, which
        ! depends on it through the kinetic-energy term in pressure())
        ! silently assumes the piston isn't moving at all.
        v_target = 0.0_DOUBLE
        if (present(n_bc) .and. re < 0) then
          if (-re <= n_bc) then
            if (trim(bc_type(-re)) == 'piston') v_target = bc_val(2:4, -re)
          end if
        end if
        sol_r = sol(:, le)
        sol_r(2:4) = sol_l(2:4) + 2.0_DOUBLE*(dot_product(v_target, mesh%face(id_face)%norm) &
          - dot_product(sol_l(2:4), mesh%face(id_face)%norm))*mesh%face(id_face)%norm
      end if
      sol_w_l = lag_to_primit(sol_l, gamma_arr(le))
      sol_w_r = lag_to_primit(sol_r, gamma_arr(le))
      grad_p = grad_p + (sol_w_r(5) - sol_w_l(5)) &
        * mesh%sub_face(id_sub_face)%area&
        * mesh%face(id_face)%norm
    end do
    grad_p = grad_p / mesh%vert(id_vert)%volume

    if( mesh%vert(id_vert)%is_bound ) then
      Bp = boundary_normal(mesh, id_vert, b2d)
      Bp = Bp/norm2(Bp)
      grad_p = grad_p - dot_product(grad_p, Bp)*Bp
    end if

    h_p = compute_length(mesh, id_vert, method, b2d, b2d_h)
    vp = vp - 0.5_DOUBLE*h_p/(rho_p*a_p)*grad_p

    ! The correction above only ever adds a tangential contribution (grad_p
    ! was already projected against Bp), but vp's base value -- the
    ! volume-weighted average of the neighboring cells' velocities, before
    ! any correction -- was never projected: a general flow field has no
    ! reason for that average to be tangential to the wall already, so the
    ! final vp used to move the boundary node could still carry a spurious
    ! wall-normal (penetrating) component.
    if( mesh%vert(id_vert)%is_bound ) then
      vp = vp - dot_product(vp, Bp)*Bp
    end if
  end subroutine compute_nodal_velocity_sidil

  subroutine compute_nodal_pressure_sidil(mesh, id_vert, sol, pp, method, b2d, b2d_h, gamma_arr, &
      n_bc, bc_type, bc_val)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert, method
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), intent(inout) :: pp
    real(kind=DOUBLE), intent(in) :: b2d_h
    logical :: b2d
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    integer(kind=ENTIER), intent(in), optional :: n_bc
    character(len=255), dimension(:), intent(in), optional :: bc_type
    real(kind=DOUBLE), dimension(:, :), intent(in), optional :: bc_val

    integer(kind=ENTIER) :: j, id_elem, id_sub_elem, le, re
    integer(kind=ENTIER) :: id_face, id_sub_face
    real(kind=DOUBLE), dimension(5) :: sol_w
    real(kind=DOUBLE) :: rho_p, a_p, h_p, div_v
    real(kind=DOUBLE), dimension(3) :: v_target
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r, sol_l, sol_r

    pp = 0.0_DOUBLE
    rho_p = 0.0_DOUBLE
    a_p = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      sol_w = lag_to_primit(sol(:, id_elem), gamma_arr(id_elem))
      pp = pp &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(5)
      rho_p = rho_p &
        + mesh%sub_elem(id_sub_elem)%volume * sol_w(1)
      a_p = a_p &
        + mesh%sub_elem(id_sub_elem)%volume &
        * sqrt(gamma_arr(id_elem)*sol_w(5)/sol_w(1))
    end do
    pp = pp / mesh%vert(id_vert)%volume
    rho_p = rho_p / mesh%vert(id_vert)%volume
    a_p = a_p / mesh%vert(id_vert)%volume

    div_v = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      ! See the identical exclusion + comment in compute_nodal_velocity_sidil.
      if (b2d .and. abs(mesh%face(id_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE) cycle
      le = mesh%face(id_face)%left_neigh
      re = mesh%face(id_face)%right_neigh
      sol_l = sol(:, le)
      if( re > 0 ) then
        sol_r = sol(:, re)
      else
        ! Same piston-vs-static-wall mirror fix as compute_nodal_velocity_sidil.
        v_target = 0.0_DOUBLE
        if (present(n_bc) .and. re < 0) then
          if (-re <= n_bc) then
            if (trim(bc_type(-re)) == 'piston') v_target = bc_val(2:4, -re)
          end if
        end if
        sol_r = sol(:, le)
        sol_r(2:4) = sol_l(2:4) + 2.0_DOUBLE*(dot_product(v_target, mesh%face(id_face)%norm) &
          - dot_product(sol_l(2:4), mesh%face(id_face)%norm))*mesh%face(id_face)%norm
      end if
      sol_w_l = lag_to_primit(sol_l, gamma_arr(le))
      sol_w_r = lag_to_primit(sol_r, gamma_arr(le))
      div_v = div_v &
        + dot_product(sol_w_r(2:4) - sol_w_l(2:4), &
        mesh%sub_face(id_sub_face)%area&
        * mesh%face(id_face)%norm)
    end do
    div_v = div_v / mesh%vert(id_vert)%volume
    h_p = compute_length(mesh, id_vert, method, b2d, b2d_h)
    pp = pp - 0.5_DOUBLE*h_p*rho_p*a_p*div_v
    if( pp < 0.0_DOUBLE ) then
      print*, "Neg prex in sidil nodal"
      error stop
    end if
  end subroutine compute_nodal_pressure_sidil

  function compute_length(mesh, id_vert, method, b2d, b2d_h) result(h_p)
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert, method
    logical, intent(in) :: b2d
    real(kind=DOUBLE), intent(in) :: b2d_h
    real(kind=DOUBLE) :: h_p

    integer(kind=ENTIER) :: j, k
    integer(kind=ENTIER) :: id_sub_elem, id_sub_face
    integer(kind=ENTIER) :: id_elem, id_face
    real(kind=DOUBLE) :: area_sum, vol_2d, length_scale, coeff

    area_sum = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
      if( b2d ) then
        ! The lateral faces' area scales with the (purely numerical, not
        ! physical) extrusion thickness -- only the z-normal "cap" faces
        ! (top/bottom of the single layer) carry a genuine 2D footprint
        ! area, independent of that thickness. Sum it as a plain scalar
        ! area (not corner_normal's vector-summed normal, which mixes
        ! several differently-oriented faces): each sub_elem contributes
        ! its footprint area once (top and bottom caps have equal area, so
        ! either one alone already represents it -- summing both would
        ! just double-count the same physical quantity).
        do k=1, mesh%sub_elem(id_sub_elem)%n_sub_faces
          id_sub_face = mesh%sub_elem(id_sub_elem)%sub_face(k)
          id_face = mesh%sub_face(id_sub_face)%mesh_face
          if( abs(mesh%face(id_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE ) then
            area_sum = area_sum + mesh%sub_face(id_sub_face)%area
            exit
          end if
        end do
      else
        area_sum = area_sum + norm2(corner_normal(mesh, id_sub_elem, b2d))
      end if
    end do

    if( b2d ) then
      ! mesh%vert(...)%volume is a genuine 3D volume (physical 2D area
      ! times the artificial extrusion thickness): dividing by b2d_h (the
      ! same thickness, set by the user to match the mesh's actual
      ! extrusion) recovers the physical 2D area. area_sum is now itself a
      ! physical area (see above), so sqrt(area_sum) is the length scale to
      ! divide by, keeping h_p homogeneous to a length either way.
      vol_2d = mesh%vert(id_vert)%volume / b2d_h
      length_scale = sqrt(area_sum)
    else
      vol_2d = mesh%vert(id_vert)%volume
      length_scale = area_sum
    end if

    if( method == 0 ) then
      coeff = 0.0_DOUBLE
    else if( method == 1) then
      coeff = 1.0_DOUBLE
    else if( method == 2) then
      coeff = 2.0_DOUBLE
    else if( method == 3) then
      coeff = 4.0_DOUBLE
    else if( method == 4) then
      coeff = 8.0_DOUBLE
    end if
    h_p = coeff*vol_2d/length_scale
  end function compute_length

  function boundary_normal(mesh, id_vert, boundary_2d) result(Bp)
    ! boundary_2d is passed in explicitly (matching every other routine in
    ! this file, e.g. compute_nodal_velocity_sidil's b2d argument) rather
    ! than pulled from lagrange_global_data_module: lagrange_main.F90
    ! declares its OWN local boundary_2d (read from the namelist, default
    ! .true.), which shadows the module's own (unrelated, always-.FALSE.
    ! default, never assigned) variable of the same name -- using the
    ! module one here silently always saw .FALSE., regardless of the actual
    ! run's boundary_2d setting.
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    logical, intent(in) :: boundary_2d

    real(kind=DOUBLE), dimension(3) :: Bp
    integer(kind=ENTIER) :: j
    integer(kind=ENTIER) :: id_sub_face, id_face

    Bp = 0.0_DOUBLE
    do j=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      if( mesh%face(id_face)%right_neigh <= 0 ) then
        ! In boundary_2d, the mesh is a single prism/hex layer: the z-normal
        ! "cap" faces (top/bottom of that one layer) are also boundary
        ! faces (right_neigh<=0) but carry no physical wall -- including
        ! them here mixes a spurious z-component into Bp, corrupting the
        ! wall-tangential projection for every genuinely 2D lateral wall.
        if (.not. (boundary_2d .and. abs(mesh%face(id_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE)) then
          Bp = Bp &
            + mesh%sub_face(id_sub_face)%area*mesh%face(id_face)%norm
        end if
      end if
    end do
  end function boundary_normal

  function corner_normal(mesh, id_sub_elem, boundary_2d) result(norm)
    ! boundary_2d passed in explicitly -- same reasoning as boundary_normal
    ! above: lagrange_global_data_module's own boundary_2d is a dead,
    ! always-.FALSE. variable, shadowed everywhere else by lagrange_main's
    ! local one (the one actually read from the namelist and threaded as
    ! the b2d argument). Pulling it from the module here silently disabled
    ! the entire z-cap-face exclusion below (.not. .FALSE. .or. ... is
    ! always true), so area_sum in compute_length mixed in the z-normal cap
    ! faces' area alongside the genuine lateral faces -- those scale with
    ! the artificial extrusion thickness (h_extrude, e.g. 1e-4) completely
    ! differently than the physical 2D geometry, corrupting the nodal
    ! characteristic length h_p for every method that goes through
    ! compute_length (the "classic" and "sidil" schemes).
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_sub_elem
    logical, intent(in) :: boundary_2d
    real(kind=DOUBLE), dimension(3) :: norm

    integer(kind=ENTIER) :: j, id_sub_face, id_face

    norm = 0.0_DOUBLE
    do j=1, mesh%sub_elem(id_sub_elem)%n_sub_faces
      id_sub_face = mesh%sub_elem(id_sub_elem)%sub_face(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      if( .not. boundary_2d .or. &
        ( boundary_2d .and. abs(mesh%face(id_face)%norm(3)) < 1e-12_DOUBLE) ) then
        if( mesh%face(id_face)%left_neigh  &
          == mesh%sub_elem(id_sub_elem)%mesh_elem ) then
          norm = norm + mesh%sub_face(id_sub_face)%area*mesh%face(id_face)%norm
        else
          norm = norm - mesh%sub_face(id_sub_face)%area*mesh%face(id_face)%norm
        end if
      end if
    end do
  end function corner_normal

  function lag_to_primit(u, g) result(w)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: u
    real(kind=DOUBLE), intent(in) :: g
    real(kind=DOUBLE), dimension(5) :: w

    w(1) = 1.0_DOUBLE/u(1)
    w(2:4) = u(2:4)
    w(5) = pressure(u, g)
  end function lag_to_primit

  subroutine move_mesh(mesh, vp, dt)
    implicit none

    type(mesh_type), intent(inout) :: mesh
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: vp
    real(kind=DOUBLE), intent(in) :: dt

    integer(kind=ENTIER) :: i

    do i=1, mesh%n_vert
      mesh%vert(i)%coord = mesh%vert(i)%coord + dt*vp(:, i)
    end do
  end subroutine move_mesh  
  
  function omega_pos(V, deltaV)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: V, deltaV
    real(kind=DOUBLE) :: omega_pos

    real(kind=DOUBLE) :: newEps
    real(kind=DOUBLE), dimension(5) :: newV

    if (any(deltaV /= deltaV)) then
      omega_pos = 0.0_DOUBLE
      return
    end if

    omega_pos = 1.0_DOUBLE
    newV = V + omega_pos * deltaV
    newEps = newV(5) - 0.5_DOUBLE*dot_product(newV(2:4),newV(2:4))
    do while (.not. (newEps >= 0.0_DOUBLE) .and. omega_pos > 1e-14_DOUBLE)
      omega_pos = 0.9_DOUBLE*omega_pos
      newV = V + omega_pos * deltaV
      newEps = newV(5) - 0.5_DOUBLE*dot_product(newV(2:4),newV(2:4))
    end do
    if (.not. (newEps >= 0.0_DOUBLE)) omega_pos = 0.0_DOUBLE
  end function omega_pos

  ! L1/L2/Linf density error against the steady isentropic vortex (init=2).
  ! Lagrangian specific: the mesh has MOVED, so the exact cell average is
  ! re-quadratured on the deformed cell at the current vertex positions.
  ! The vortex has zero background velocity, so it is an exact steady state
  ! of the Euler equations and the exact field is the t=0 one at the moved
  ! cell's own location (t is kept in the signature for symmetry).
  subroutine compute_error_vortex(mesh, sol, t, h_extrude, h, l1err, l2err, linferr)
    use quadrature_module, only: volume_quad_pts
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), intent(in) :: t, h_extrude
    real(kind=DOUBLE), intent(out) :: h, l1err, l2err, linferr

    integer(kind=ENTIER) :: i, k, q, n_v, n_inner, mpi_ierr
    real(kind=DOUBLE) :: err1, err2, errinf, area_tot, area_i, rho_num, rho_ex, vol_q
    real(kind=DOUBLE), dimension(5) :: wq
    real(kind=DOUBLE), dimension(:, :), allocatable :: vcoords, pts
    real(kind=DOUBLE), dimension(:), allocatable :: wts
    real(kind=DOUBLE), dimension(4) :: buf

    err1     = 0.0_DOUBLE
    err2     = 0.0_DOUBLE
    errinf   = 0.0_DOUBLE
    area_tot = 0.0_DOUBLE
    n_inner  = 0

    do i = 1, mesh%n_elems
      if (mesh%elem(i)%is_ghost) cycle

      n_v = mesh%elem(i)%n_vert
      rho_ex = 0.0_DOUBLE
      if (n_v == 4 .or. n_v == 5 .or. n_v == 6 .or. n_v == 8) then
        allocate(vcoords(3, n_v))
        do k = 1, n_v
          vcoords(:, k) = mesh%vert(mesh%elem(i)%vert(k))%coord
        end do
        call volume_quad_pts(n_v, vcoords, 3_ENTIER, pts, wts)
        deallocate(vcoords)
        vol_q = 0.0_DOUBLE
        do q = 1, size(wts)
          call sol_isentropic_vortex(pts(:, q), wq, t)
          rho_ex = rho_ex + wts(q)*wq(1)
          vol_q  = vol_q  + wts(q)
        end do
        deallocate(pts, wts)
        if (vol_q > 0.0_DOUBLE) then
          rho_ex = rho_ex/vol_q
        else
          call sol_isentropic_vortex(mesh%elem(i)%coord, wq, t)
          rho_ex = wq(1)
        end if
      else
        call sol_isentropic_vortex(mesh%elem(i)%coord, wq, t)
        rho_ex = wq(1)
      end if

      rho_num = 1.0_DOUBLE/sol(1, i)
      area_i  = mesh%elem(i)%volume/h_extrude

      err1     = err1 + area_i*abs(rho_num - rho_ex)
      err2     = err2 + area_i*(rho_num - rho_ex)**2
      errinf   = max(errinf, abs(rho_num - rho_ex))
      area_tot = area_tot + area_i
      n_inner  = n_inner + 1
    end do

    buf = (/err1, err2, area_tot, real(n_inner, DOUBLE)/)
    call MPI_ALLREDUCE(MPI_IN_PLACE, buf, 3, MPI_DOUBLE_PRECISION, MPI_SUM, &
      MPI_COMM_WORLD, mpi_ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, buf(4), 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
      MPI_COMM_WORLD, mpi_ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, errinf, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
      MPI_COMM_WORLD, mpi_ierr)
    err1     = buf(1)
    err2     = buf(2)
    area_tot = buf(3)

    if (buf(4) > 0.0_DOUBLE .and. area_tot > 0.0_DOUBLE) then
      h     = sqrt(area_tot/buf(4))
      l1err = err1/area_tot
      l2err = sqrt(err2/area_tot)
    else
      h = 0.0_DOUBLE; l1err = 0.0_DOUBLE; l2err = 0.0_DOUBLE
    end if
    linferr = errinf
  end subroutine compute_error_vortex

  ! 2D footprint area of a cell on an extruded "2D" mesh. Takes the area of the
  ! z-normal face directly, the way euler_ho's compute_error_vortex does, instead
  ! of volume/h_extrude: on a MOVING Lagrangian mesh the latter silently assumes
  ! the z-extent is still the one auto-detected at t=0, and it hides a degenerate
  ! cell (a NaN or zero volume shows up as a plausible area). Falls back to
  ! volume/h_extrude when no z-face is found (a genuinely 3D mesh).
  function cell_area_2d(mesh, i, h_extrude) result(area)
    use quadrature_module, only: face_quad_pts
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE), intent(in) :: h_extrude
    real(kind=DOUBLE) :: area

    integer(kind=ENTIER) :: kf, iface_loc, n_fv, k
    real(kind=DOUBLE), dimension(:, :), allocatable :: fc, qpts_v
    real(kind=DOUBLE), dimension(:), allocatable :: qwts_v

    area = -1.0_DOUBLE
    do kf = 1, mesh%elem(i)%n_faces
      iface_loc = mesh%elem(i)%face(kf)
      if (abs(abs(mesh%face(iface_loc)%norm(3)) - 1.0_DOUBLE) < 1.0e-6_DOUBLE) then
        n_fv = mesh%face(iface_loc)%n_vert
        if (n_fv >= 3 .and. allocated(mesh%face(iface_loc)%vert)) then
          allocate(fc(3, n_fv))
          do k = 1, n_fv
            fc(:, k) = mesh%vert(mesh%face(iface_loc)%vert(k))%coord
          end do
          call face_quad_pts(n_fv, fc, 3_ENTIER, qpts_v, qwts_v)
          deallocate(fc)
          area = sum(qwts_v)
          deallocate(qpts_v, qwts_v)
          exit
        end if
      end if
    end do
    if (area <= 0.0_DOUBLE) area = mesh%elem(i)%volume/h_extrude
  end function cell_area_2d

  ! ---------------------------------------------------------------------------
  ! Taylor-Green vortex, following Vilar (PhD 2012, Sec. 4.2.4), the standard
  ! convergence test case for cell-centred Lagrangian hydrodynamics.
  !
  !   domain [0,1]^2,  rho = rho0  (uniform),
  !   U = U0 ( sin(pi x) cos(pi y), -cos(pi x) sin(pi y) )   -- divergence free,
  !   P = rho0 U0^2 / 4 * [cos(2 pi x) + cos(2 pi y)] + C0.
  !
  ! The velocity field being divergence free, the volume equation holds identically,
  ! and that P is exactly what balances the inertia in the momentum equation. The
  ! total energy equation does NOT close on its own: it needs the source term in
  ! taylor_green_source (Vilar eq. 4.102-4.103). WITH that source the Eulerian
  ! fields are steady, so the exact solution at any time is the initial one, which
  ! is what compute_error_taylor_green relies on. WITHOUT it the computed solution
  ! drifts and the error measurement is meaningless.
  !
  ! Note Vilar measures the rate on the PRESSURE, never the density -- here rho is
  ! uniform and carries no signal at all.
  ! ---------------------------------------------------------------------------
  pure subroutine sol_taylor_green(coord, w)
    use lagrange_global_data_module
    implicit none

    real(kind=DOUBLE), dimension(3), intent(in) :: coord
    real(kind=DOUBLE), dimension(5), intent(inout) :: w

    real(kind=DOUBLE) :: x, y

    x = coord(1)
    y = coord(2)
    w(1) = tg_rho0
    w(2) =  tg_u0*sin(PI*x)*cos(PI*y)
    w(3) = -tg_u0*cos(PI*x)*sin(PI*y)
    w(4) = 0.0_DOUBLE
    w(5) = 0.25_DOUBLE*tg_rho0*tg_u0**2 &
      *(cos(2.0_DOUBLE*PI*x) + cos(2.0_DOUBLE*PI*y)) + tg_c0
  end subroutine sol_taylor_green

  ! Energy source per unit VOLUME that closes the total energy equation, Vilar
  ! eq. (4.103). Re-derived independently and it matches the printed formula.
  pure function taylor_green_source(coord) result(src)
    use lagrange_global_data_module
    implicit none

    real(kind=DOUBLE), dimension(3), intent(in) :: coord
    real(kind=DOUBLE) :: src, x, y

    x = coord(1)
    y = coord(2)
    src = PI*tg_rho0*tg_u0**3/(4.0_DOUBLE*(gamma - 1.0_DOUBLE)) &
      *(cos(3.0_DOUBLE*PI*x)*cos(PI*y) - cos(3.0_DOUBLE*PI*y)*cos(PI*x))
  end function taylor_green_source

  ! Exact cell average of the Lagrangian conserved state (tau, u, E) over cell i.
  subroutine taylor_green_cell_average(mesh, i, u_avg)
    use quadrature_module, only: volume_quad_pts
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE), dimension(5), intent(out) :: u_avg

    integer(kind=ENTIER) :: k, q, n_v
    real(kind=DOUBLE) :: vol_q, e_tot
    real(kind=DOUBLE), dimension(5) :: w
    real(kind=DOUBLE), dimension(:, :), allocatable :: vcoords, pts
    real(kind=DOUBLE), dimension(:), allocatable :: wts

    u_avg = 0.0_DOUBLE
    n_v = mesh%elem(i)%n_vert
    if (n_v == 4 .or. n_v == 5 .or. n_v == 6 .or. n_v == 8) then
      allocate(vcoords(3, n_v))
      do k = 1, n_v
        vcoords(:, k) = mesh%vert(mesh%elem(i)%vert(k))%coord
      end do
      call volume_quad_pts(n_v, vcoords, 3_ENTIER, pts, wts)
      deallocate(vcoords)
      vol_q = 0.0_DOUBLE
      do q = 1, size(wts)
        call sol_taylor_green(pts(:, q), w)
        e_tot = w(5)/(w(1)*(gamma - 1.0_DOUBLE)) + 0.5_DOUBLE*norm2(w(2:4))**2
        u_avg(2:4) = u_avg(2:4) + wts(q)*w(2:4)
        u_avg(5)   = u_avg(5)   + wts(q)*e_tot
        vol_q = vol_q + wts(q)
      end do
      deallocate(pts, wts)
      if (vol_q > 0.0_DOUBLE) then
        u_avg(2:5) = u_avg(2:5)/vol_q
      end if
    else
      call sol_taylor_green(mesh%elem(i)%coord, w)
      u_avg(2:4) = w(2:4)
      u_avg(5)   = w(5)/(w(1)*(gamma - 1.0_DOUBLE)) + 0.5_DOUBLE*norm2(w(2:4))**2
    end if
    u_avg(1) = 1.0_DOUBLE/tg_rho0     ! rho is uniform, so tau is exact
  end subroutine taylor_green_cell_average

  ! Adds the Taylor-Green energy source to the RHS. rhs(5,:) is dE/dt PER UNIT
  ! MASS (new_sol = sol + dt*rhs with sol(5) the specific total energy), while the
  ! source is per unit volume, hence the integral over the cell divided by its mass.
  subroutine add_taylor_green_source(mesh, mass, rhs)
    use quadrature_module, only: volume_quad_pts
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mass
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs

    integer(kind=ENTIER) :: i, k, q, n_v
    real(kind=DOUBLE) :: src_int
    real(kind=DOUBLE), dimension(:, :), allocatable :: vcoords, pts
    real(kind=DOUBLE), dimension(:), allocatable :: wts

    do i = 1, mesh%n_elems
      if (mesh%elem(i)%is_ghost) cycle
      if (mass(i) <= 0.0_DOUBLE) cycle
      n_v = mesh%elem(i)%n_vert
      src_int = 0.0_DOUBLE
      if (n_v == 4 .or. n_v == 5 .or. n_v == 6 .or. n_v == 8) then
        allocate(vcoords(3, n_v))
        do k = 1, n_v
          vcoords(:, k) = mesh%vert(mesh%elem(i)%vert(k))%coord
        end do
        call volume_quad_pts(n_v, vcoords, 3_ENTIER, pts, wts)
        deallocate(vcoords)
        do q = 1, size(wts)
          src_int = src_int + wts(q)*taylor_green_source(pts(:, q))
        end do
        deallocate(pts, wts)
      else
        src_int = mesh%elem(i)%volume*taylor_green_source(mesh%elem(i)%coord)
      end if
      rhs(5, i) = rhs(5, i) + src_int/mass(i)
    end do
  end subroutine add_taylor_green_source

  ! Error on the PRESSURE, as Vilar does. The exact reference is the pressure
  ! re-quadratured on the DEFORMED cell at the final time: the Eulerian field is
  ! steady, so the exact solution is the initial one evaluated at the cell's
  ! current position.
  subroutine compute_error_taylor_green(mesh, sol, gamma_arr, h_extrude, h, l1err, l2err, linferr)
    use quadrature_module, only: volume_quad_pts
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    real(kind=DOUBLE), intent(in) :: h_extrude
    real(kind=DOUBLE), intent(out) :: h, l1err, l2err, linferr

    integer(kind=ENTIER) :: i, k, q, n_v, n_inner, mpi_ierr
    real(kind=DOUBLE) :: err1, err2, errinf, area_tot, area_i, p_num, p_ex, vol_q
    real(kind=DOUBLE), dimension(5) :: wq
    real(kind=DOUBLE), dimension(:, :), allocatable :: vcoords, pts
    real(kind=DOUBLE), dimension(:), allocatable :: wts
    real(kind=DOUBLE), dimension(4) :: buf

    err1 = 0.0_DOUBLE; err2 = 0.0_DOUBLE; errinf = 0.0_DOUBLE
    area_tot = 0.0_DOUBLE; n_inner = 0

    do i = 1, mesh%n_elems
      if (mesh%elem(i)%is_ghost) cycle

      n_v = mesh%elem(i)%n_vert
      p_ex = 0.0_DOUBLE
      if (n_v == 4 .or. n_v == 5 .or. n_v == 6 .or. n_v == 8) then
        allocate(vcoords(3, n_v))
        do k = 1, n_v
          vcoords(:, k) = mesh%vert(mesh%elem(i)%vert(k))%coord
        end do
        call volume_quad_pts(n_v, vcoords, 3_ENTIER, pts, wts)
        deallocate(vcoords)
        vol_q = 0.0_DOUBLE
        do q = 1, size(wts)
          call sol_taylor_green(pts(:, q), wq)
          p_ex  = p_ex + wts(q)*wq(5)
          vol_q = vol_q + wts(q)
        end do
        deallocate(pts, wts)
        if (vol_q > 0.0_DOUBLE) then
          p_ex = p_ex/vol_q
        else
          call sol_taylor_green(mesh%elem(i)%coord, wq)
          p_ex = wq(5)
        end if
      else
        call sol_taylor_green(mesh%elem(i)%coord, wq)
        p_ex = wq(5)
      end if

      p_num  = pressure(sol(:, i), gamma_arr(i))
      area_i = cell_area_2d(mesh, i, h_extrude)

      err1     = err1 + area_i*abs(p_num - p_ex)
      err2     = err2 + area_i*(p_num - p_ex)**2
      errinf   = max(errinf, abs(p_num - p_ex))
      area_tot = area_tot + area_i
      n_inner  = n_inner + 1
    end do

    buf = (/err1, err2, area_tot, real(n_inner, DOUBLE)/)
    call MPI_ALLREDUCE(MPI_IN_PLACE, buf, 4, MPI_DOUBLE_PRECISION, MPI_SUM, &
      MPI_COMM_WORLD, mpi_ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, errinf, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
      MPI_COMM_WORLD, mpi_ierr)
    err1 = buf(1); err2 = buf(2); area_tot = buf(3)

    ! A diverged run makes area_tot NaN, and "NaN > 0" is false, so the plain
    ! else branch below would report h=0, L1=0, L2=0 -- a blow-up disguised as a
    ! zero error. Flag it explicitly instead.
    if (.not. (area_tot == area_tot)) then
      h = -1.0_DOUBLE; l1err = -1.0_DOUBLE; l2err = -1.0_DOUBLE; linferr = -1.0_DOUBLE
      return
    end if
    if (buf(4) > 0.0_DOUBLE .and. area_tot > 0.0_DOUBLE) then
      h     = sqrt(area_tot/buf(4))
      l1err = err1/area_tot
      l2err = sqrt(err2/area_tot)
    else
      h = 0.0_DOUBLE; l1err = 0.0_DOUBLE; l2err = 0.0_DOUBLE
    end if
    linferr = errinf
  end subroutine compute_error_taylor_green

  subroutine solve_vp_pp_pvp_systems(nsfn, sol_lr, gammas, norms, areas, vp, pp, pvp)
    implicit none

    integer(kind=ENTIER), intent(in) :: nsfn
    real(kind=DOUBLE), dimension(5, 2, nsfn), intent(in) :: sol_lr
    real(kind=DOUBLE), dimension(3, nsfn), intent(in) :: norms
    real(kind=DOUBLE), dimension(nsfn), intent(in) :: areas
    real(kind=DOUBLE), dimension(2, nsfn), intent(in) :: gammas
    real(kind=DOUBLE), dimension(3), intent(inout) :: vp, pvp
    real(kind=DOUBLE), intent(inout) :: pp

    integer(kind=ENTIER) :: i
    real(kind=DOUBLE) :: vnl, vnr, lambda_l, lambda_r, taul, taur, pl, pr, el, er
    real(kind=DOUBLE) :: pgod, vgod, pvgod, ppdenom

    real(kind=DOUBLE), dimension(3, 3) :: mat_vp, mat_pvp
    real(kind=DOUBLE), dimension(3) :: rhs_vp, rhs_pvp, acc_bal

    vp = 0.0_DOUBLE
    mat_vp = 0.0_DOUBLE
    rhs_vp = 0.0_DOUBLE
    pp = 0.0_DOUBLE
    ppdenom = 0.0_DOUBLE
    pvp = 0.0_DOUBLE
    mat_pvp = 0.0_DOUBLE
    rhs_pvp = 0.0_DOUBLE
    acc_bal = 0.0_DOUBLE
    do i=1, nsfn
      taul = sol_lr(1, 1, i)
      vnl = dot_product(sol_lr(2:4, 1, i), norms(:, i))
      el = sol_lr(5, 1, i)
      pl = pressure(sol_lr(:, 1, i), gammas(1, i))

      taur = sol_lr(1, 2, i)
      vnr = dot_product(sol_lr(2:4, 2, i), norms(:, i))
      er = sol_lr(5, 2, i)
      pr = pressure(sol_lr(:, 2, i), gammas(2, i))

      lambda_l = max(sqrt(max(0.0_DOUBLE, gammas(1, i)*pl*taul))/taul, &
        sqrt(max(0.0_DOUBLE, pr-pl)/taul), -(vnr-vnl)/taul)
      lambda_r = max(sqrt(max(0.0_DOUBLE, gammas(2, i)*pr*taur))/taur, &
        sqrt(max(0.0_DOUBLE, pl-pr)/taur), -(vnr-vnl)/taur)

      pgod = (pl/lambda_l + pr/lambda_r - (vnr-vnl))&
        /(1.0_DOUBLE/lambda_l+1.0_DOUBLE/lambda_r)
      
      vgod = (lambda_l*vnl + lambda_r*vnr - (pr - pl))&
        /(lambda_l + lambda_r)

      pvgod = (pl*vnl/lambda_l + pr*vnr/lambda_r - (er - el))&
        /(1.0_DOUBLE/lambda_l+1.0_DOUBLE/lambda_r)

      mat_vp = mat_vp &
        + areas(i)*tensor_product(norms(:, i), norms(:, i))&
        *(lambda_l + lambda_r)
      rhs_vp = rhs_vp &
        + areas(i)*(lambda_l+lambda_r)*vgod*norms(:, i)

      pp = pp &
        + areas(i)*(1.0_DOUBLE/lambda_l+1.0_DOUBLE/lambda_r)*pgod
      ppdenom = ppdenom &
        + areas(i)*(1.0_DOUBLE/lambda_l+1.0_DOUBLE/lambda_r)

      mat_pvp = mat_pvp &
        + areas(i)*(1.0_DOUBLE/lambda_l+1.0_DOUBLE/lambda_r)&
        *tensor_product(norms(:, i), norms(:, i))
      rhs_pvp = rhs_pvp &
        + areas(i)*(1.0_DOUBLE/lambda_l+1.0_DOUBLE/lambda_r)*pvgod*norms(:, i)

      ! mode 'bal' : la sous-face porte DEUX sous-elements, c = l et c = r, de
      ! normales opposees -- n (x) n est identique, seule v_c change.
      acc_bal = acc_bal + areas(i)*(vnl + vnr)*norms(:, i)
    end do

    pp = pp / ppdenom
    call pseudo_inverse_inplace_lapack(3, mat_vp)
    vp = matmul(mat_vp, rhs_vp)
    if (trim(pvp_energy_mode) == 'bal') then
      mat_pvp = 0.0_DOUBLE
      do i = 1, nsfn
        mat_pvp = mat_pvp + 2.0_DOUBLE*areas(i)*tensor_product(norms(:, i), norms(:, i))
      end do
      rhs_pvp = pp*acc_bal
    end if
    call pseudo_inverse_inplace_lapack(3, mat_pvp)
    pvp = matmul(mat_pvp, rhs_pvp)
  end subroutine solve_vp_pp_pvp_systems

  subroutine assemble_vp_pp_pvp_rhs(mesh, sol, vp, rhs, n_bc, bc_type, bc_val, &
      b2d, mass, gamma_arr, pp_node, pvp_node, theta)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: vp
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(out) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mass, gamma_arr
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: pp_node, theta
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: pvp_node
    integer(kind=ENTIER), intent(in) :: n_bc
    character(len=255), dimension(n_bc) :: bc_type
    real(kind=DOUBLE), dimension(5, n_bc) :: bc_val
    logical, intent(in) :: b2d

    integer(kind=ENTIER) :: i, j, id_sub_face, id_face, idl, idr, var
    real(kind=DOUBLE), dimension(5) :: flux_l, flux_r, solr
    real(kind=DOUBLE), dimension(3) :: nrm
    real(kind=DOUBLE) :: taul, taur, pl, pr, vl, vr, gl, gr, lam_l, lam_r
    real(kind=DOUBLE) :: vn, pi_l, pi_r, th, v_target, area

    rhs = 0.0_DOUBLE
    pvp_bnd_power = 0.0_DOUBLE
    pvp_bnd_abs = 0.0_DOUBLE
    do i = 1, mesh%n_vert
      th = theta(i)
      do j = 1, mesh%vert(i)%n_sub_faces_neigh
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        nrm = mesh%face(id_face)%norm
        if (b2d .and. abs(nrm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE) cycle
        idl = mesh%face(id_face)%left_neigh
        idr = mesh%face(id_face)%right_neigh
        area = mesh%sub_face(id_sub_face)%area
        vn = dot_product(vp(:, i), nrm)

        gl = gamma_arr(idl)
        taul = sol(1, idl)
        pl = pressure(sol(:, idl), gl)
        vl = dot_product(sol(2:4, idl), nrm)
        if (idr > 0) then
          solr = sol(:, idr)
          gr = gamma_arr(idr)
        else
          v_target = 0.0_DOUBLE
          if (idr < 0) then
            if (-idr <= n_bc) then
              if (trim(bc_type(-idr)) == 'piston') &
                v_target = dot_product(bc_val(2:4, -idr), nrm)
            end if
          end if
          solr = mirror_state(sol(:, idl), nrm, v_target)
          gr = gl
        end if
        taur = solr(1)
        pr = pressure(solr, gr)
        vr = dot_product(solr(2:4), nrm)

        lam_l = max(sqrt(max(0.0_DOUBLE, gl*pl*taul))/taul, &
          sqrt(max(0.0_DOUBLE, pr-pl)/taul), -(vr-vl)/taul)
        lam_r = max(sqrt(max(0.0_DOUBLE, gr*pr*taur))/taur, &
          sqrt(max(0.0_DOUBLE, pl-pr)/taur), -(vr-vl)/taur)
        pi_l = pl - lam_l*(vn - vl)
        pi_r = pr + lam_r*(vn - vr)

        flux_l(1)   = -vn
        flux_l(2:4) = (th*pp_node(i) + (1.0_DOUBLE-th)*pi_l)*nrm
        flux_l(5)   = th*dot_product(pvp_node(:, i), nrm) + (1.0_DOUBLE-th)*pi_l*vn
        flux_r(1)   = -vn
        flux_r(2:4) = (th*pp_node(i) + (1.0_DOUBLE-th)*pi_r)*nrm
        flux_r(5)   = th*dot_product(pvp_node(:, i), nrm) + (1.0_DOUBLE-th)*pi_r*vn

        do var = 1, 5
          rhs(var, idl) = rhs(var, idl) - area/mass(idl)*flux_l(var)
        end do
        if (idr > 0) then
          do var = 1, 5
            rhs(var, idr) = rhs(var, idr) + area/mass(idr)*flux_r(var)
          end do
        else
          pvp_bnd_power = pvp_bnd_power + area*flux_l(5)
          pvp_bnd_abs = pvp_bnd_abs + abs(area*flux_l(5))
        end if
      end do
    end do
  end subroutine assemble_vp_pp_pvp_rhs

  subroutine compute_rhs_lagrange_vp_pp_pvp(mesh, sol, vp, dt, rhs, &
      n_bc, bc_type, bc_val, b2d, mass, gamma_arr, vp_is_imposed)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), intent(in) :: dt
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: vp
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: mass
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: gamma_arr
    logical, dimension(mesh%n_vert), intent(in) :: vp_is_imposed
    integer(kind=ENTIER), intent(in) :: n_bc
    character(len=255), dimension(n_bc) :: bc_type
    real(kind=DOUBLE), dimension(5, n_bc) :: bc_val
    logical, intent(in) :: b2d

    integer(kind=ENTIER) :: i, j, k, it, nsfn, id_sub_face, id_face, idl, idr
    integer(kind=ENTIER) :: n_bad
    real(kind=DOUBLE), dimension(:, :, :), allocatable :: sol_lr
    real(kind=DOUBLE), dimension(:, :), allocatable :: norms, gammas
    real(kind=DOUBLE), dimension(:), allocatable :: areas
    real(kind=DOUBLE), dimension(:), allocatable :: pp_node, theta
    real(kind=DOUBLE), dimension(:, :), allocatable :: pvp_node
    logical, dimension(:), allocatable :: bad
    real(kind=DOUBLE), dimension(3) :: vp_loc, pvp_loc, u_new
    real(kind=DOUBLE) :: pp_loc, v_target, e_old, e_new

    allocate(pp_node(mesh%n_vert), pvp_node(3, mesh%n_vert), theta(mesh%n_vert))
    allocate(bad(mesh%n_elems))
    pp_node = 0.0_DOUBLE
    pvp_node = 0.0_DOUBLE
    theta = 1.0_DOUBLE

    do i = 1, mesh%n_vert
      nsfn = mesh%vert(i)%n_sub_faces_neigh
      if (nsfn <= 0) cycle
      allocate(sol_lr(5, 2, nsfn), norms(3, nsfn), areas(nsfn), gammas(2, nsfn))
      k = 0
      do j = 1, nsfn
        id_sub_face = mesh%vert(i)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        if (b2d .and. abs(mesh%face(id_face)%norm(3)) > 1.0_DOUBLE - 1e-8_DOUBLE) cycle
        idl = mesh%face(id_face)%left_neigh
        idr = mesh%face(id_face)%right_neigh
        k = k + 1
        norms(:, k) = mesh%face(id_face)%norm
        areas(k) = mesh%sub_face(id_sub_face)%area
        sol_lr(:, 1, k) = sol(:, idl)
        gammas(1, k) = gamma_arr(idl)
        if (idr > 0) then
          sol_lr(:, 2, k) = sol(:, idr)
          gammas(2, k) = gamma_arr(idr)
        else
          v_target = 0.0_DOUBLE
          if (idr < 0) then
            if (-idr <= n_bc) then
              if (trim(bc_type(-idr)) == 'piston') &
                v_target = dot_product(bc_val(2:4, -idr), mesh%face(id_face)%norm)
            end if
          end if
          sol_lr(:, 2, k) = mirror_state(sol(:, idl), mesh%face(id_face)%norm, v_target)
          gammas(2, k) = gamma_arr(idl)
        end if
      end do
      if (k > 0) then
        vp_loc = vp(:, i)
        call solve_vp_pp_pvp_systems(k, sol_lr(:, :, 1:k), gammas(:, 1:k), &
          norms(:, 1:k), areas(1:k), vp_loc, pp_loc, pvp_loc)
        if (.not. vp_is_imposed(i)) vp(:, i) = vp_loc
        pp_node(i) = pp_loc
        pvp_node(:, i) = pvp_loc
      end if
      deallocate(sol_lr, norms, areas, gammas)
    end do

    call assemble_vp_pp_pvp_rhs(mesh, sol, vp, rhs, n_bc, bc_type, bc_val, &
      b2d, mass, gamma_arr, pp_node, pvp_node, theta)

    if (pvp_limit) then
      ! theta est attache au NOEUD, jamais a la sous-face : la conservation de la
      ! part 'classic' du melange repose sur l'equilibre nodal
      !   somme_f A (Pi_r - Pi_l) n = 0,
      ! qui est l'equation meme definissant vp. Avec un theta par sous-face la
      ! somme devient somme_f A (1-theta_f)(Pi_r - Pi_l) n et ne s'annule plus.
      ! Avec un theta par noeud il sort de la somme et la conservation est exacte.
      do it = 1, pvp_lim_iter
        bad = .false.
        n_bad = 0
        do i = 1, mesh%n_elems
          if (mesh%elem(i)%is_ghost) cycle
          e_old = sol(5, i) - 0.5_DOUBLE*dot_product(sol(2:4, i), sol(2:4, i))
          u_new = sol(2:4, i) + dt*rhs(2:4, i)
          e_new = (sol(5, i) + dt*rhs(5, i)) - 0.5_DOUBLE*dot_product(u_new, u_new)
          if (e_new /= e_new .or. e_new < pvp_eps*max(e_old, 0.0_DOUBLE)) then
            bad(i) = .true.
            n_bad = n_bad + 1
          end if
        end do
        ! n_bad reste LOCAL : theta n'est pas echange entre rangs, donc un noeud
        ! partage peut recevoir deux theta differents et la conservation tombe a
        ! l'interface MPI. Le limiteur n'est valide qu'en sequentiel pour l'instant.
        if (n_bad == 0) exit
        do i = 1, mesh%n_vert
          do j = 1, mesh%vert(i)%n_sub_faces_neigh
            id_sub_face = mesh%vert(i)%sub_face_neigh(j)
            id_face = mesh%sub_face(id_sub_face)%mesh_face
            idl = mesh%face(id_face)%left_neigh
            idr = mesh%face(id_face)%right_neigh
            if (idl > 0) then
              if (bad(idl)) theta(i) = 0.0_DOUBLE
            end if
            if (idr > 0) then
              if (bad(idr)) theta(i) = 0.0_DOUBLE
            end if
          end do
        end do
        if (it < pvp_lim_iter) then
          do i = 1, mesh%n_vert
            if (theta(i) == 0.0_DOUBLE) cycle
            theta(i) = 0.5_DOUBLE*theta(i)
          end do
        end if
        call assemble_vp_pp_pvp_rhs(mesh, sol, vp, rhs, n_bc, bc_type, bc_val, &
          b2d, mass, gamma_arr, pp_node, pvp_node, theta)
      end do
      pvp_theta_min = minval(theta)
      pvp_n_limited = count(theta < 1.0_DOUBLE)
    end if

    deallocate(pp_node, pvp_node, theta, bad)
  end subroutine compute_rhs_lagrange_vp_pp_pvp
end module lagrange_module
