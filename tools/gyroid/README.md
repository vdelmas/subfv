# Maillage tétraédrique conforme d'une TPMS (gyroïde) dans un cube

`gyroid_mesh.py` construit, dans le cube `[0,L]^3`, les **trois** régions
séparées par les deux surfaces décalées `f = +delta` et `f = -delta` :

| région  | définition           | `Physical Volume` |
|---------|----------------------|-------------------|
| `solid` | `-delta < f <= +delta` (la nappe, épaisseur `~ 2*delta/|grad f|`) | 1 |
| `airA`  | `f > +delta`  (premier labyrinthe d'air)  | 2 |
| `airB`  | `f <= -delta` (second labyrinthe d'air)   | 3 |

avec, pour la gyroïde,
`f = sin(kx)cos(ky) + sin(ky)cos(kz) + sin(kz)cos(kx)`, `k = 2*pi/period`.

Le maillage est **conforme** : les trois régions partagent les mêmes nœuds et
les mêmes faces aux interfaces (vérifié par `--check`).

## Utilisation

    ./gyroid_mesh.py --nper 1 --n 48 --delta 0.35 -o gyroide --check
    ./gyroid_mesh.py --nper 2 --n 96 --delta 0.5 --size 0.05 -o gyroide2
    ./gyroid_mesh.py --surface schwarz_p --n 48 --delta 0.4 -o schwarz
    ./gyroid_mesh.py --n 32 --no-mesh          # complexe frontière seulement

Sorties : `<out>_surface.msh` (complexe frontière), `<out>.geo` (lisible et
éditable), `<out>.msh` (maillage volumique, msh 4.1 par défaut).

Options principales :

| option | rôle |
|--------|------|
| `--n` | cellules de la grille de fond par direction = **résolution des nappes** |
| `--delta` | demi-décalage des deux lignes de niveau |
| `--nper` | nombre de périodes par direction (`L = nper*period`) |
| `--snap` | recalage des nœuds sur les niveaux (défaut 0.10), voir plus bas |
| `--size` | `Mesh.MeshSizeMax` pour le remplissage (0 = gmsh prolonge les tailles du bord) |
| `--algo3d` | `Mesh.Algorithm3D` (1 = Delaunay, 10 = HXT) |
| `--surface` | `gyroid`, `schwarz_p`, `diamond` |
| `--check` | relit le maillage, vérifie volumes et conformité |

Groupes physiques de surface : `xmin xmax ymin ymax zmin zmax` (chacun couvrant
les trois régions sur cette face), `sheet_plus`, `sheet_minus`.

## Pourquoi pas un simple `.geo`

gmsh ne maille pas de surface implicite et ses noyaux CAO ne construisent pas de
TPMS. Le script construit donc la **frontière** des trois régions et gmsh ne
fait que le remplissage 3D.

1. grille structurée, chaque voxel coupé en 6 tétraèdres (subdivision de
   Kuhn/Freudenthal). La triangulation induite sur une face de voxel ne dépend
   que de cette face : elle est identique vue des deux côtés, donc étanche, et
   **aucun cas ambigu** n'est à lever (contrairement au marching cubes) ;
2. marching tetrahedra à `+delta` et `-delta` -> les deux nappes ; chaque point
   d'intersection est créé sur une arête de grille et mémorisé dans une table
   globale indexée par `(arête, niveau)` ;
3. les 6 faces du cube sont découpées par les mêmes lignes de niveau **en
   réutilisant ces mêmes nœuds d'arête** : la trace des nappes sur les faces
   coïncide exactement avec la triangulation des faces, sans tolérance ni
   recollement ;
4. tout est écrit dans un seul `.msh` avec **une seule liste de nœuds** : le
   complexe est conforme par construction ;
5. les composantes connexes de chaque région sont détectées et chacune devient
   un `Volume` gmsh distinct, regroupés dans un même `Physical Volume`.
   **Nécessaire** : si plusieurs coques déconnectées sont mises dans le même
   `Surface Loop`, gmsh en maille silencieusement une seule et abandonne les
   autres (avec `--nper 1` la nappe est tranchée par le cube et donne 3
   composantes, dont deux éclats de coin de volume `3e-5`) ;
6. le `.geo` définit les `Surface Loop` / `Volume` / `Physical`, `gmsh -3`
   remplit.

## `--snap` : qualité des tétraèdres contre erreur géométrique

Les nœuds dont la valeur est à moins de `snap` fois la variation maximale de `f`
sur une arête sont recalés exactement sur le niveau. Les coupes dégénérées
s'effondrent alors sur des nœuds existants au lieu de produire des triangles
dégénérés, qui sont ce qui détruit la qualité des tétraèdres.

Mesuré sur `--n 32 --delta 0.35` (gyroïde, 1 période) :

| `snap` | tri. par nappe | tets `q<0.1` | distance à la surface exacte (rms / max) | volume solide |
|--------|---------------|--------------|------------------------------------------|---------------|
| 0      | 27648 | 3791 | 0.013 h / 0.057 h | 22.733 % |
| 0.10   | 21240 |  150 | 0.019 h / 0.099 h | 22.777 % |
| 0.15   | 18828 |   44 | —                 | 22.779 % |
| 0.25   | 13200 |   23 | 0.069 h / 0.233 h | 22.704 % |

Le script affiche cette distance à la surface exacte à chaque exécution : le
compromis reste donc chiffré. Défaut `0.10`.

## Vérifications faites par le script

- chaque coque de région est fermée (toute arête partagée par exactement deux
  triangles) et d'orientation cohérente ;
- volume enclos par chaque coque > 0, et somme des trois régions = `L^3`
  (erreur de fermeture `~1e-16`) ;
- après `gmsh -3` (`--check`) : volume maillé par groupe physique, total = `L^3`
  à la précision machine, faces partagées aux interfaces, absence de nœuds
  dupliqués et de face intérieure sans région voisine.

Si une de ces vérifications échoue, le script le dit et **n'appelle pas gmsh**
(ou retourne un code d'erreur) au lieu de produire un maillage douteux.

## Convergence

Volume du solide (gyroïde, 1 période, `delta = 0.35`, `snap = 0`) :
`n = 16 -> 23.181 %`, `n = 32 -> 22.733 %`, `n = 64 -> 22.622 %`.
Les écarts successifs sont dans un rapport `4.4`, c'est-à-dire la convergence
`O(h^2)` attendue pour une isosurface linéaire par morceaux ; extrapolation de
Richardson `~22.58 %`.

## Coût indicatif

`--nper 1 --n 48 --delta 0.35` : 128 736 triangles de frontière, 2 s pour le
complexe, 50 s pour `gmsh -3`, 943 000 éléments, 135 000 nœuds.

## Configurations validées

| cas | éléments | `gmsh -3` | somme des volumes | conformité |
|-----|----------|-----------|-------------------|------------|
| gyroïde, `--n 48 --delta 0.35` (1 période) | 943 000 | 50 s | `L^3` à 2.7e-15 | 48 576 faces partagées de chaque côté, 0 nœud dupliqué |
| gyroïde, `--nper 2 --n 64 --delta 0.4` | 353 000 nœuds | 108 s | `8` à 4.8e-14 | 165 504 faces partagées de chaque côté, 0 nœud dupliqué |
| Schwarz-P, `--n 32 --delta 0.4` | 37 829 nœuds | 4.6 s | `L^3` à 1.4e-14 | 14 904 faces partagées de chaque côté, 0 nœud dupliqué |

Dans tous les cas `airA / airB` partagent **0** face : les deux labyrinthes sont
bien séparés par la nappe, ce qui est la topologie attendue d'une TPMS.

Sur Schwarz-P les deux régions d'air ont le même volume (38.603 %) mais pas le
même nombre de triangles de frontière (24 420 contre 16 368) : elles sont
congruentes par la translation `(L/2, L/2, L/2)`, qui est une symétrie du réseau
mais pas du cube, donc les faces du cube les découpent différemment.
