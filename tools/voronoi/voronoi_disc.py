#!/usr/bin/env python3
"""Voronoi radial d'un disque, extrude d'une seule couche en z, ecrit en CGNS NGON/NFACE.

Reproduit la topologie du maillage de Sedov polygonal de Vilar (these 2012, fig. 4.13a) :
germes places sur des anneaux concentriques, le nombre de germes par anneau croissant avec
le rayon pour garder des cellules de taille a peu pres constante. Vilar travaille sur un
QUART de disque (775 cellules, rayon 1.2) avec des conditions de symetrie sur les axes ;
ici on construit le quart puis on le MIROITE dans les quatre quadrants, ce qui donne un
disque complet exactement symetrique d'ordre 4.

Deux ecarts assumes par rapport a la figure de Vilar, tous deux documentes :
  - un germe est ajoute a l'ORIGINE. Sans lui, les cellules du premier anneau se rejoignent
    au centre et le depot d'energie de Sedov (qui va dans la cellule la plus proche de
    l'origine) tomberait arbitrairement dans l'une d'elles, brisant la symetrie.
  - le bord exterieur est le polygone engendre par l'inversion des germes dans le cercle,
    donc une approximation polygonale de l'arc -- exactement comme sur sa figure.

Usage :
  ./vor-venv/bin/python tools/voronoi/voronoi_disc.py sedov_disc.cgns --R 1.2 --nr 32
"""
import argparse
import numpy as np
from scipy.spatial import Voronoi
import CGNS.PAT.cgnslib as CL
import CGNS.PAT.cgnskeywords as CK
import CGNS.MAP as CM

ap = argparse.ArgumentParser()
ap.add_argument("out")
ap.add_argument("--R", type=float, default=1.2, help="rayon du disque")
ap.add_argument("--nr", type=int, default=32, help="nombre d'anneaux sur le rayon")
ap.add_argument("--z", type=float, nargs=2, default=[0.0, 1.0], help="z0 z1 de l'unique couche")
ap.add_argument("--quarter", action="store_true", help="ne garder que le quart (comme Vilar)")
ap.add_argument("--no-center-seed", action="store_true", help="pas de germe a l'origine")
ap.add_argument("--bc", default="all_surf", help="nom de la CL du bord lateral")
ap.add_argument("--caps-bc", default="", help="nom de CL des chapeaux z (vide = non tague)")
a = ap.parse_args()

R, nr = a.R, a.nr
dr = R / nr
h = dr                                   # taille caracteristique, utilisee par la queue

# --- germes : anneaux concentriques, quart de disque puis miroir ---------------
quarter = []
if a.quarter and not a.no_center_seed:
    # Sur un QUART, un germe a l'origine est une erreur geometrique : sa cellule de
    # Voronoi est un losange centre sur l'origine qui deborde dans les trois autres
    # quadrants (verifie : sommets a +/-0.053 en x ET en y pour nr=8). Les cellules
    # miroir voisines etant jetees, ses aretes deviennent des faces de BORD de
    # normales a +/-45 deg, donc non colineaires : le solveur lagrangien les classe
    # toutes en 'coin' et epingle a v=0 les quatre sommets de la maille qui porte
    # toute l'energie de Sedov -- le calcul reste fige (|v| ~ 1e-22 partout).
    # Vilar n'a pas de germe au coin non plus. On desactive donc d'office.
    a.no_center_seed = True
    print("quart : germe central desactive (il deborderait du quadrant)")
if not a.no_center_seed:
    quarter.append(np.zeros((1, 2)))      # cellule centrale (voir docstring)
for k in range(nr):
    rk = (k + 0.5) * dr
    nq = max(int(round(2.0 * np.pi * rk / dr / 4.0)), 1)   # germes sur le quart d'anneau
    th = (np.arange(nq) + 0.5) * (0.5 * np.pi / nq)        # decale : aucun germe sur un axe
    quarter.append(np.c_[rk * np.cos(th), rk * np.sin(th)])
q = np.vstack(quarter)

if a.quarter:
    seeds = q
else:
    # miroir dans les quatre quadrants ; le germe central n'est pas duplique
    c = q[np.all(np.abs(q) < 1e-14, axis=1)]
    p = q[~np.all(np.abs(q) < 1e-14, axis=1)]
    seeds = np.vstack([c, p, p * [-1, 1], p * [1, -1], -p])
ns = len(seeds)

# --- Voronoi avec miroir dans le CERCLE ---------------------------------------
def voronoi(s):
    """Inversion des germes par rapport au cercle de rayon R : un germe au rayon r a
    son image au rayon 2R-r sur le meme rayon, donc la bissectrice des deux tombe
    exactement a r = R. Les aretes de Voronoi epousent ainsi le bord du disque.
    Pour le quart, on ajoute les miroirs par rapport aux deux axes."""
    rad = np.linalg.norm(s, axis=1)
    keep = rad > 1e-14
    u = np.zeros_like(s)
    u[keep] = s[keep] / rad[keep, None]
    m = [u[keep] * (2.0 * R - rad[keep])[:, None]]
    if a.quarter:
        # Le germe a l'ORIGINE est son propre miroir sous les trois symetries
        # d'axes : l'inclure envoie quatre copies de (0,0) a scipy, ce qui produit
        # un maillage geometriquement plausible mais degenere (la cellule centrale
        # se retrouve isolee et le solveur lagrangien se fige, toutes les vitesses
        # nodales a ~1e-22). On ne miroite donc que les germes hors origine.
        sq = s[keep]
        m += [sq * [-1, 1], sq * [1, -1], -sq]
        mm = m[0]
        m += [mm * [-1, 1], mm * [1, -1], -mm]
    return Voronoi(np.vstack([s] + m))


def polygon(vor, i):
    reg = vor.regions[vor.point_region[i]]
    if -1 in reg:
        raise SystemExit("cellule interieure non bornee")
    return reg


vor = voronoi(seeds)
V = vor.vertices
used = sorted({v for i in range(ns) for v in polygon(vor, i)})
parent = {v: v for v in used}


def find(v):
    while parent[v] != v:
        parent[v] = parent[parent[v]]
        v = parent[v]
    return v


for (p0, p1), (v0, v1) in zip(vor.ridge_points, vor.ridge_vertices):
    if (p0 < ns or p1 < ns) and v0 >= 0 and v1 >= 0 and np.linalg.norm(V[v0] - V[v1]) < 0.05 * h:
        parent[find(v0)] = find(v1)
roots = sorted({find(v) for v in used})
vid = {r: k for k, r in enumerate(roots)}
V2 = np.array([V[r] for r in roots])

# Recalage radial : les miroirs placent les sommets de bord sur le cercle a l'arrondi
# pres, et la fusion de sommets ci-dessus peut en tirer un legerement vers l'interieur.
rad = np.linalg.norm(V2, axis=1)
onb = np.abs(rad - R) < 0.02 * h
V2[onb] *= (R / rad[onb])[:, None]
if a.quarter:                              # et sur les deux axes de symetrie
    for c in (0, 1):
        V2[np.abs(V2[:, c]) < 0.02 * h, c] = 0.0

polys = []
for i in range(ns):
    reg = [vid[find(v)] for v in polygon(vor, i)]
    poly = [v for k, v in enumerate(reg) if v != reg[k - 1]]
    x, y = V2[poly, 0], V2[poly, 1]
    if 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y) < 0:
        poly = poly[::-1]
    polys.append(poly)

edge = {}
for c, poly in enumerate(polys):
    for k in range(len(poly)):
        p, q = poly[k], poly[(k + 1) % len(poly)]
        edge.setdefault((min(p, q), max(p, q)), []).append((c, p, q))

nv = len(V2)
faces, cell_faces, side_ids, cap_ids = [], [[] for _ in polys], [], []
for c, poly in enumerate(polys):
    faces.append([v + 1 for v in poly[::-1]])  # z0 cap, outward -z
    cell_faces[c].append(len(faces)); cap_ids.append(len(faces))
    faces.append([v + 1 + nv for v in poly])   # z1 cap, outward +z
    cell_faces[c].append(len(faces)); cap_ids.append(len(faces))
for owners in edge.values():
    c, p, q = owners[0]                          # p->q is c's CCW order: quad is outward for c
    faces.append([p + 1, q + 1, q + 1 + nv, p + 1 + nv])
    fid = len(faces)
    cell_faces[c].append(fid)
    if len(owners) == 2:
        cell_faces[owners[1][0]].append(-fid)
    elif len(owners) == 1:
        side_ids.append(fid)
    else:
        raise SystemExit("non-manifold edge")

nodes = np.vstack([np.c_[V2, np.full(nv, a.z[0])], np.c_[V2, np.full(nv, a.z[1])]])


def csr(lists):
    off = np.cumsum([0] + [len(l) for l in lists]).astype(np.int64)
    return np.concatenate([np.asarray(l, dtype=np.int64) for l in lists]), off


nf, nc = len(faces), len(cell_faces)
T = CL.newCGNSTree()
zone = CL.newZone(CL.newBase(T, "Base", 3, 3), "Zone", np.array([[len(nodes), nc, 0]], dtype=np.int64),
                  CK.Unstructured_s)
gc = CL.newGridCoordinates(zone, CK.GridCoordinates_s)
for k, name in enumerate([CK.CoordinateX_s, CK.CoordinateY_s, CK.CoordinateZ_s]):
    CL.newDataArray(gc, name, np.ascontiguousarray(nodes[:, k]))
conn, off = csr(faces)
CL.newElements(zone, "Faces", CK.NGON_n_s, np.array([1, nf], dtype=np.int64), conn, off)
conn, off = csr(cell_faces)
CL.newElements(zone, "Cells", CK.NFACE_n_s, np.array([nf + 1, nf + nc], dtype=np.int64), conn, off)
zbc = CL.newZoneBC(zone)
groups = {a.bc: list(side_ids)}
if a.caps_bc:
    groups.setdefault(a.caps_bc, []).extend(cap_ids)
for name, ids in groups.items():
    bc = CL.newBoundary(zbc, name, np.array([ids], dtype=np.int64), "BCWall", pttype=CK.PointList_s)
    CL.newGridLocation(bc, CK.FaceCenter_s)
CM.save(a.out, T)

areas = np.array([0.5 * abs(np.sum(V2[p, 0] * np.roll(V2[p, 1], -1) - np.roll(V2[p, 0], -1) * V2[p, 1])) for p in polys])

sides = np.array([len(p) for p in polys])
areas = np.array([abs(0.5 * np.sum(V2[p, 0] * np.roll(V2[p, 1], -1)
                                   - np.roll(V2[p, 0], -1) * V2[p, 1])) for p in polys])
elen = [np.linalg.norm(V2[p[k]] - V2[p[k - 1]]) for p in polys for k in range(len(p))]
exact = 0.25 * np.pi * R * R if a.quarter else np.pi * R * R
print(f"longueur d'arete min {min(elen) / h:.3f} h")
print(f"{nc} cellules, cotes min/moy/max {sides.min()}/{sides.mean():.2f}/{sides.max()}, "
      f"aire totale {areas.sum():.9f} (disque exact {exact:.9f}, ecart "
      f"{100 * abs(areas.sum() / exact - 1):.4f} %), "
      f"rapport d'aires min/max {areas.min() / areas.max():.3f}, {len(side_ids)} faces laterales")
