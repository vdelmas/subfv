# all_mach_number_scheme

Suite de cas tests autonomes (indépendants de ctest/CMake) pour valider le
schéma `WIP` (voir `../tex/wip.tex` pour la présentation du schéma et
`../tex_all_mach_number_scheme/` pour le rapport de cette étude), en local ou
sur le cluster Curta.

Ces dossiers sont extraits de `../test/<case>/` : mêmes maillages, mêmes
paramètres physiques (`input.template`), mais `scheme='WIP'` est figé en dur
(plus de boucle multi-schémas) et ordre 1 uniquement (`WIP_o2_*` non repris
ici, cf. `../tex/wip.tex` qui ne référence que `outputs/WIP`).

## Cas tests

| Dossier | Binaire | NPROC par défaut | Remarque |
|---|---|---|---|
| `test_shear` | `subfvshear` | 1 | 4 maillages (`mesh_shear_1/2/4/8`) lancés en séquence |
| `gresho_quad` | `subfvns` | 1 | `periodic_mesh=.true.` force NPROC=1 (maillage quad périodique non partitionnable) |
| `gresho_tri` | `subfvns` | 1 | idem — `test/gresho_tri/run.sh` n'avait pas ce garde-fou (contrairement à `gresho_quad`), ce qui plantait bien sur le cluster à 32 rangs ("Error periodic boundaries !", `mesh_connectivity_module.F90` : le rattachement périodique exige des comptes de faces gauche/droite égaux **sur chaque rang**, non garanti par un partitionnement MPI quelconque). Garde-fou ajouté ici (2026-09-15). |
| `convergence_gresho_quad` | `subfvgreshoconv` | 1 (imposé) | balaie plusieurs Mach en un seul run |
| `convergence_gresho_tri` | `subfvgreshoconv` | 1 (imposé) | idem |
| `sedov_tri` | `subfvns` | 4 | |
| `sedov_hex` | `subfvns` | 4 | maillage ~1.07M éléments / 129 Mo — job cluster avec `--mem=32G` |
| `half_cylinder_quad_euler` | `subfvns` | 4 | |
| `half_cylinder_quad_ns` | `subfvns` | 4 | produit aussi `plot_cp.pdf`/`plot_st.pdf` (comparaison LAURA/fun3D) |
| `half_cylinder_tri_euler` | `subfvns` | 4 | |
| `half_cylinder_tri_ns` | `subfvns` | 4 | idem quad_ns |
| `gresho_summary` | — | — | combine les 4 cas gresho/convergence en un `summary.pdf` (pdflatex, **local uniquement**, Curta n'a pas pdflatex) |

Chaque dossier de cas contient `run.sh` (calcul), `visu.sh`
(post-traitement pvbatch/gnuplot) et `job.sbatch` (gabarit Slurm Curta qui
enchaîne les deux). Relancer un cas :

```bash
cd all_mach_number_scheme/<case>
./run.sh && ./visu.sh          # en local
# ou, sur le cluster :
sbatch job.sbatch
```

## Build dédié

Ces scripts pointent par défaut vers `../../build_all_mach_number_scheme`
(variable `SUBFV_BUILD`, surchargeable), **jamais** vers `../../build` /
`/scratch/vdelmas/subfv/build` — ce dossier partagé est utilisé par d'autres
études en cours (ex. `aho_study`) avec des jobs actifs ; ne jamais y faire de
`rm -rf`/reconfiguration. Construire le build dédié (une fois, ou après
modification du code) :

```bash
cd /scratch/vdelmas/subfv   # sur le cluster
mkdir -p build_all_mach_number_scheme && cd build_all_mach_number_scheme
module load gcc/11.2.0 lapack mpi blas
cmake ..
make -j subfvns subfvshear subfvgreshoconv subfv-gmsh
```

## Déploiement cluster

1. Sync (depuis le poste local, dans `subfv/subfv/`) :
   ```bash
   rsync -avz --exclude='test/' --exclude='large_test/' \
     ./ vdelmas@curta.mcia.fr:/scratch/vdelmas/subfv/
   ```
   (`all_mach_number_scheme/` est inclus automatiquement — seuls `test/` et
   `large_test/` sont exclus.)
2. Construire `build_all_mach_number_scheme` (voir ci-dessus).
3. Depuis `/scratch/vdelmas/subfv/all_mach_number_scheme/` sur le cluster :
   `./submit_all.sh` (soumet les 11 cas ; régénère systématiquement
   maillage/partition à chaque soumission).
4. Vérifier que chaque job se termine réellement (`run.timestamp` présent +
   `sacct -j <id>` = `COMPLETED`, pas seulement "il y a un fichier de sortie" —
   voir la mémoire cluster sur les faux positifs de jobs `TIMEOUT`).
5. Rapatrier les résultats vers le poste local :
   ```bash
   rsync -avz --exclude='*.msh' --exclude='*_[0-9]*.msh' \
     vdelmas@curta.mcia.fr:/scratch/vdelmas/subfv/all_mach_number_scheme/*/outputs/ \
     ./<case>/outputs/   # par cas, ou boucle équivalente
   ```
6. Lancer `gresho_summary/combine.sh` en local une fois les 4 cas gresho
   rapatriés.
