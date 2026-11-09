# Maillages de Voronoï pour subfv

## Pourquoi CGNS et pas gmsh

Une cellule de Voronoï est un polygone quelconque (ici 4 à 8 côtés). Le format
`.msh` de gmsh n'a pas d'élément polygonal général : il n'y a donc pas de moyen
d'écrire un maillage de Voronoï en `.msh`. On passe par CGNS avec les sections
`NGON_n` (faces) / `NFACE_n` (cellules), qui décrivent des polyèdres arbitraires.

Côté code, aucune modification n'est nécessaire : `read_mesh_msh`
(`src/mesh/mesh_reading_module.F90:157`) est en fait un dispatcher, qui route
toute extension `.cgns` vers `read_mesh_cgns_dist`
(`src/partition/dist_mesh_module.F90:452`). Chaque rang lit le fichier par blocs
et le partitionnement est fait dans le code par ParMETIS — il n'y a donc **pas**
de pré-découpage `subfv-gmsh` à faire (et d'ailleurs `subfv-gmsh` ne sait pas
lire de `.geo` : il attend un `.msh` déjà maillé).

Tous les exécutables qui passent par ce dispatcher en profitent : `subfvns`
(`src/ns/ns_main.F90:42`) comme `subfvlagrange`
(`src/lagrange/lagrange_main.F90:84`).

Il faut en revanche que le binaire soit construit avec `-DSUBFV_ENABLE_CGNS=ON`,
sinon le lecteur s'arrête avec un message explicite. C'est le cas en local et,
depuis le rebuild du 2026-09-29, sur Curta (`cmake.log` doit afficher
« subfv: CGNS reader with ParMETIS partitioning »).

## Générateurs

- `voronoi_box.py` — Voronoï d'un rectangle `[x0,x1]x[y0,y1]`, extrudé d'une
  seule couche en z (donc utilisable avec `boundary_2d = .true.`).
- `voronoi_disc.py` — Voronoï radial d'un disque, topologie du maillage de Sedov
  polygonal de Vilar (thèse 2012, fig. 4.13a).
- `voronoi_halfcyl.py` — Voronoï du demi-cylindre (cas flux de chaleur).
- `voronoi_halfcyl_bl.py` — même domaine, mais couche limite structurée en quads
  sur la paroi et Voronoï au loin, conforme à l'interface.

### Dépendances

`numpy`, `scipy` (pour `scipy.spatial.Voronoi`) et `pycgns`. **Attention :
pycgns et scipy ne sont installés ni dans le python système local, ni sur
Curta.** Un venv suffit en local :

```sh
python3 -m venv --system-site-packages vor-venv
./vor-venv/bin/pip install scipy pycgns
```

Comme le cluster n'a pas ces modules, **les `.cgns` sont générés en local puis
envoyés par rsync** avec le reste du cas test (un maillage 1600 cellules pèse
830 Ko, c'est négligeable).

### Exemple : le maillage du Gresho lagrangien

Le domaine doit reproduire exactement celui de `gresho_quad.geo` : boîte
`[-0.5,0.5]^2`, une couche en z d'épaisseur `lc = 1/N`, et toutes les surfaces
étiquetées `wall` (y compris les chapeaux en z, exclus ensuite par
`boundary_2d`) :

```sh
./vor-venv/bin/python tools/voronoi/voronoi_box.py gresho_voro.cgns \
    --n 40 --box -0.5 0.5 -0.5 0.5 --z 0 0.025 \
    --jitter 0.3 --lloyd 3 --bc wall --caps-bc wall
```

Ce qui donne 1600 cellules — **le même nombre que le maillage quad 40x40**, pour
que la comparaison soit à résolution égale — de 4 à 8 côtés (5,70 en moyenne), et
une somme des aires égale à l'aire de la boîte à 1e-12 près (contrôle de
non-recouvrement).

`--lloyd` lisse les cellules (0 = Voronoï jitteré brut, plus dur) ; `--periodic`
construit des germes symétriques par rapport aux deux axes médians, de sorte que
les côtés opposés se correspondent exactement, pour `periodic_mesh = .true.`.

## Coût : le pas de temps

À nombre de cellules égal, une cellule de Voronoï a un rayon inscrit plus petit
qu'un quadrangle, donc un `dt` limité par CFL plus petit. Mesuré sur le Gresho
lagrangien `classic` à Ma = 1e-2 : `dt = 1.30e-5` sur Voronoï contre
`2.48e-5` sur quad, soit **un facteur 1,9 de temps de calcul en plus** pour le
même `t_max`. Les `job.sbatch` des cas Voronoï prévoient donc le double du
walltime de leur équivalent quad/tri.

## Vérifier qu'un maillage est bien lu

Le lecteur CGNS est bavard, et ces lignes suffisent à valider l'ingestion :

```
[cgns] BC 'wall' (family '') matches no bc_name entry -> tag 0
[cgns] boundary faces: 3360 in BCs, 3360 attached to a cell
[cgns] 1600 cells, 6088 nodes, partitioner parmetis, 1 ghost layer(s)
```

- les deux nombres de « boundary faces » doivent être **égaux** (toute face de
  bord déclarée est bien rattachée à une cellule) ;
- `matches no bc_name entry -> tag 0` est normal quand `n_bc = 0` (c'est le cas
  du Gresho lagrangien, dont le `.geo` étiquette aussi tout en `wall`) ; avec
  `n_bc > 0`, le nom passé à `--bc` doit figurer dans la liste `bc_name`, sinon
  la condition sera silencieusement remplacée par le tag 0.

## Le piège du mode périodique (à lire avant d'étendre aux schémas Euler)

Les deux familles de cas Gresho n'ont **pas** les mêmes conditions de bord, et
c'est ce qui décide des options du générateur :

| cas | `periodic_mesh` | nom de BC | options |
|---|---|---|---|
| `test/lagrange_gresho` (lagrangien, `n_bc = 0`) | absent, donc `.false.` | `wall` | `--bc wall --caps-bc wall` |
| `test/gresho_quad` (Euler / `subfvns`) | **`.true.`** | `all_surf` | `--periodic --bc all_surf --caps-bc all_surf`, `--n` pair |

`stitch_periodic_faces` / `build_periodic_stencil` codent en dur la boîte
`[-0.5,0.5]^3`, exigent que les faces et sommets opposés se correspondent à
1e-12, et apparient aussi les chapeaux en z. Un maillage de Voronoï construit
par miroir n'est **pas** périodique et plante là avec un indice -1 dans
`build_periodic_stencil` : d'où `--periodic`, qui engendre les germes dans un
seul quadrant puis les réfléchit à travers les deux axes médians, de sorte que
les sommets de bord opposés sont bit-à-bit identiques. **Conséquence à
mentionner dans tout rapport : le maillage périodique est alors
miroir-symétrique par rapport à x = 0 et y = 0**, ce n'est pas un Voronoï
générique. (Une vraie tessellation périodique découpée sur la boîte produisait
des cellules en éclats d'aire quasi nulle, piste abandonnée.)

Enfin `--caps-bc` : le `all_surf` de `test/gresho_quad` englobe aussi les deux
chapeaux en z (surfaces 5 et 27 du `.geo`). Sans `--caps-bc`, ils restent au
tag 0 et le lecteur le signale par une ligne d'information.

## Régularité des cellules : mettre du Lloyd pour le lagrangien

`--jitter 0.3` seul donne un rapport d'aires min/max d'environ 0,4.
Mesuré sur la boîte du Gresho à `--n 40` : `--lloyd 3` donne 0,658 et
`--lloyd 20` donne 0,707 (min edge 0,052 h dans les deux cas).

Pour un solveur **lagrangien**, les cellules se déforment au cours du calcul :
partir d'un maillage déjà régulier évite de confondre une dégénérescence
géométrique avec un défaut de schéma. On lance donc les deux (`voro` =
`--lloyd 3`, `voro_l20` = `--lloyd 20`) pour pouvoir attribuer un éventuel échec
à l'un ou à l'autre.

À `num_procs == 1` le lecteur court-circuite complètement ParMETIS : un build
CGNS suffit, ParMETIS n'est pas requis. Les maillages périodiques sont de toute
façon séquentiels uniquement.

## Un piège indépendant du Voronoï

Sur le Gresho, `scheme = multi_point_iso` avec `second_order = .true.` et
`method = 2` en périodique meurt sur « BAD OMEGA » (NaN) après ~35 itérations,
**identiquement sur le `.msh` quad d'origine (iter 36), sur le même maillage
quad réécrit en NFACE (34) et sur le Voronoï périodique (38)**. Cette
configuration diverge donc sur le Gresho quel que soit le maillage : ce n'est pas
un symptôme du Voronoï. En parois (`periodic_mesh = .false.`), la même
configuration passe 6000 itérations sans problème.

## Si on passe le Gresho Euler en parois : corriger bc_type

`test/gresho_quad` et `test_low_mach_compar_resume/gresho_quad` déclarent
`n_bc = 1`, `bc_name = 'all_surf', ...` et `bc_type = 'outflowsupersonic', ...`.
Tant que `periodic_mesh = .true.`, ce `bc_type` est **inoffensif** : les faces
sont recousues, la condition n'est jamais appliquée. C'est le cas des cas
`gresho_voro` / `convergence_gresho_voro`, qui restent périodiques.

En revanche, si on met `periodic_mesh = .false.`, il faut **impérativement**
passer `bc_type(1)` à `'wall'`. Laissé à `'outflowsupersonic'`, le résultat est
faux sans planter : `u_t` max passe de 1,0 à 8,1 dès t = 1e-3 (mesuré sur
Voronoï, et identique en np1 et np4, donc ce n'est pas un problème MPI).

Avec `bc_type = 'wall'`, quad 50x50 et Voronoï se comportent de la même façon
sur 3000 itérations (`multi_point_iso`, ordre 2, method 2) : `u_t` max
0,993 -> 0,844 sur quad contre 0,998 -> 0,802 sur Voronoï, et `rho` reste à
1,000000 sur les deux.

Le template **lagrangien** (parois, `n_bc = 0`) n'est pas concerné.

## Demi-cylindre hybride : couche limite en quads + lointain en Voronoï

`voronoi_halfcyl_bl.py` produit le maillage du demi-cylindre NS avec une **couche
limite structurée en quadrangles** sur la paroi et un **lointain en Voronoï**,
conformes à l'interface, extrudés d'une couche en z — le pendant direct des
maillages `half_cylinder_quad_ns` et `half_cylinder_tri_ns` existants.

```sh
./vor-venv/bin/python tools/voronoi/voronoi_halfcyl_bl.py half_cylinder_voro_ns.cgns
```

### Les valeurs par defaut sont calees sur le maillage TRI, pas sur le quad

C'est le point a comprendre avant de toucher aux parametres.
`half_cylinder_tri_ns.geo` **est deja un hybride** : quadrangles structures sur
$r\in[1;1{,}04]$ (`Recombine`, `nlr=16`, `Beta_HWall 1e-4`, 128 azimutales) puis
triangles jusqu'a $r=1{,}7$ avec une taille $\ell \propto r$
($\ell_2=\pi r_{bl}/128$ a l'interface, $\ell_3=\pi r_{out}/128$ au loin).
C'est donc lui le modele direct, et les defauts le transposent :
`--r-bl 1.04 --nt 128 --nr 16 --h-wall 1e-4`.

Se caler sur le maillage **quad** conduit au contraire a un lointain trop
grossier : ce maillage-la n'a que **13 mailles radiales au-dela de $r=1{,}2$**.

### Longueur caracteristique ou densite de cellules ?

Les deux ne coincident pas. A $h$ egal, une cellule de Voronoi couvre environ
$0{,}866\,h^2$ contre $0{,}433\,h^2$ pour un triangle : **deux fois plus**. Donc

- reprendre $\ell_3$ tel quel egalise la **resolution** mais donne moitie moins
  de cellules (4608 au total, contre 8326 pour le tri) ;
- diviser $\ell_3$ par $\sqrt2$ egalise la **densite**, donc le nombre de
  cellules --- c'est la convention suivie depuis le debut de l'etude (1600
  cellules face au quad $40\times40$, etc.), et c'est le defaut.

Comme le premier anneau de Voronoi doit avoir exactement $n_t$ germes pour la
conformite, sa taille est bloquee a $\pi r_{bl}/n_t$. Pour descendre en dessous
sans casser la conformite, `--iface-refine M` **subdivise la face externe de
chaque maille quad en $M$ faces** : la maille devient un polygone a $3+M$ cotes,
ce qui est parfaitement licite en NGON. Le defaut `M=2` amene la taille a
l'interface a $1{,}28\cdot10^{-2}$, soit tres pres de la derniere maille
radiale quad ($1{,}08\cdot10^{-2}$) : la transition ne presente donc de saut ni
dans un sens ni dans l'autre.

### Ce que donne la configuration par defaut

```
quad BL: 2048 cells (16 x 128), growth ratio 1.3663, first 1.000e-04, last 1.080e-02
Voronoi: 7055 cells, sides min/mean/max 4/5.89/8
interface: 256 shared faces (expected 256, 2 per quad face), max snap 1.05e-15
total 9103 cells, 32094 nodes, 43355 faces; area 2.969070743445 (exact 2.968805057642)
BCs: surf_cyl: 128, surf_in: 182, surf_out: 98
```

Soit **9103 cellules** face aux **8326** du tri (+9 %) et aux 4096 du quad. L'ecart
d'aire tombe a $+0{,}009\,\%$.

### Comment l'interface est rendue conforme (exactement)

C'est le seul point délicat, et il a une solution exacte plutôt qu'un découpage.

L'anneau externe du bloc quad a ses nœuds en $(R_I, \theta_k)$, avec
$\theta_k = \pi/2 + k\,\delta\theta$ et $\delta\theta = \pi/n_t$ ; ses faces
externes sont donc les cordes $[\theta_k, \theta_{k+1}]$.

On place le premier anneau de germes de Voronoï aux angles **milieux** des faces,
$\varphi_k = \theta_k + \delta\theta/2$, et on réfléchit chaque germe radialement
à travers le cercle de rayon

$$R_m = R_I \cos(\delta\theta/2).$$

La bissectrice d'un germe et de son miroir est alors la **tangente** à ce cercle
en $\varphi_k$. Deux tangentes voisines se coupent, par symétrie, sur le rayon
$\theta_k$ à la distance $R_m/\cos(\delta\theta/2) = R_I$ — c'est-à-dire
exactement sur le nœud quad. Le même calcul place les extrémités sur $x=0$ à la
distance exactement $R_I$. Chaque cellule de Voronoï de l'interface possède donc
**une arête droite coïncidant avec une face quad**, une pour une.

Le script **mesure** l'écart résiduel et **refuse d'écrire** s'il dépasse
`--snap-tol` (1e-9 par défaut), plutôt que de recaler un défaut sous le tapis.
Il vérifie aussi que le nombre de faces partagées à l'interface vaut bien $n_t$.
Mesuré sur la configuration par défaut : écart maximal **1,2e-15**, 64 faces
partagées pour 64 attendues.

### Ce que donne la configuration par défaut

```
quad BL: 3072 cells (48 x 64), growth ratio 1.1672, first 2.000e-05, last 2.867e-02
Voronoi: 436 cells, sides min/mean/max 4/5.63/8
interface: 64 shared faces (expected 64), max snap 1.18e-15
total 3508 cells, 7962 nodes, 14504 faces
BCs: surf_cyl: 64, surf_in: 61, surf_out: 110
```

3508 cellules, à comparer aux 4096 du quad d'origine ($64\times64$) : résolutions
comparables.

### Deux points a savoir

**L'aire n'est pas exacte, et c'est normal.** Les arcs sont approches par des
segments : le cercle interieur est *inscrit* (noeuds sur le cercle, cordes
droites) et l'arc exterieur *circonscrit* (les tangentes du miroir), deux effets
qui augmentent l'aire. Le `.msh` de gmsh fait la meme approximation. L'ecart
decroit en $\delta\theta^2$ : il vaut $+0{,}056\,\%$ a $n_t=64$ et
$+0{,}009\,\%$ a $n_t=128$.

**Ou placer l'interface.** `--r-bl` vaut 1,04 par defaut, la valeur du maillage
tri. Elle place l'interface **tres pres de la paroi**, donc le choc de culot
(decollement de l'ordre de $0{,}15\,R$ a Mach 27, soit $r\approx1{,}15$) tombe
dans la region **de Voronoi**, exactement comme il tombe dans la region
triangulaire du maillage tri. C'est voulu : le bloc structure sert a la couche
limite et au flux de chaleur parietal, pas a capturer le choc. Si l'on veut au
contraire le choc dans les quads, il faut monter `--r-bl` vers 1,2 et augmenter
`--nr` en consequence --- mais on ne compare alors plus la meme chose que le tri.
