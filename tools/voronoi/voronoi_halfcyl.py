"""Voronoi mesh of the half_cylinder_tri_euler domain (r in [1, 2.5], x <= 0),
extruded by one layer in z into polygonal prisms, written as CGNS (NGON_n/NFACE_n,
named FaceCenter BCs). Needs numpy, scipy, pycgns. Usage: voronoi_halfcyl.py out.cgns [h_in h_out seed]"""
import sys
import numpy as np
from scipy.spatial import Voronoi

R_IN, R_OUT, HZ = 1.0, 2.5, 1e-3
out = sys.argv[1]
h_in = float(sys.argv[2]) if len(sys.argv) > 2 else 0.04
h_out = float(sys.argv[3]) if len(sys.argv) > 3 else 0.08
rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 1)


def hloc(r):
    return h_in + (h_out - h_in) * (r - R_IN) / (R_OUT - R_IN)


# --- jittered seeds on rings, left half annulus -----------------------------
seeds = []
r = R_IN + 0.5 * hloc(R_IN)
while r < R_OUT - 0.25 * hloc(r):
    h = hloc(r)
    nt = max(int(round(np.pi * r / h)), 2)
    th = np.pi / 2 + (np.arange(nt) + 0.5) * np.pi / nt
    rr = r + rng.uniform(-0.3, 0.3, nt) * h
    th = th + rng.uniform(-0.3, 0.3, nt) * h / r
    rr = np.clip(rr, R_IN + 0.15 * h, R_OUT - 0.15 * h)
    th = np.clip(th, np.pi / 2 + 0.15 * h / r, 1.5 * np.pi - 0.15 * h / r)
    seeds += list(zip(rr * np.cos(th), rr * np.sin(th)))
    r += h
seeds = np.array(seeds)
ns = len(seeds)

# --- mirror seeds: across x=0 and radially across both arcs (+ corners) -----
def refl_r(p, R):
    rr = np.hypot(p[:, 0], p[:, 1])
    return p * ((2 * R - rr) / rr)[:, None]

rad = np.hypot(seeds[:, 0], seeds[:, 1])
hh = hloc(rad)
near_x = -seeds[:, 0] < 3 * hh
near_in = rad - R_IN < 3 * hh
near_out = R_OUT - rad < 3 * hh
mx = seeds[near_x] * [-1, 1]
mirrors = [mx, refl_r(seeds[near_in], R_IN), refl_r(seeds[near_out], R_OUT),
           refl_r(seeds[near_x & near_in] * [-1, 1], R_IN),
           refl_r(seeds[near_x & near_out] * [-1, 1], R_OUT)]
pts = np.vstack([seeds] + mirrors)
vor = Voronoi(pts)

# --- interior cells, merged near-duplicate vertices -------------------------
used = sorted({v for i in range(ns) for v in vor.regions[vor.point_region[i]]})
if -1 in used:
    sys.exit("unbounded interior Voronoi cell")
V = vor.vertices
parent = {v: v for v in used}

def find(v):
    while parent[v] != v:
        parent[v] = parent[parent[v]]
        v = parent[v]
    return v

# merge vertices closer than 5% of the local size (slivers from near-cocircular seeds)
for (a, b), rv in zip(vor.ridge_points, vor.ridge_vertices):
    if a >= ns and b >= ns:
        continue
    v0, v1 = rv
    if v0 < 0 or v1 < 0:
        continue
    if np.linalg.norm(V[v0] - V[v1]) < 0.05 * hloc(np.hypot(*V[v0])):
        parent[find(v0)] = find(v1)

roots = sorted({find(v) for v in used})
vid = {rt: k for k, rt in enumerate(roots)}
V2 = np.array([V[rt] for rt in roots])

polys = []
for i in range(ns):
    reg = [vid[find(v)] for v in vor.regions[vor.point_region[i]]]
    poly = [v for k, v in enumerate(reg) if v != reg[k - 1]]  # drop collapsed edges
    x, y = V2[poly, 0], V2[poly, 1]
    if 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y) < 0:
        poly = poly[::-1]
    polys.append(poly)

# --- 2D edges -> owner/neighbour ---------------------------------------------
edge = {}
for c, poly in enumerate(polys):
    for k in range(len(poly)):
        a, b = poly[k], poly[(k + 1) % len(poly)]
        edge.setdefault((min(a, b), max(a, b)), []).append((c, a, b))

def bc_of(a, b):
    m = 0.5 * (V2[a] + V2[b])
    d = {"cyl_surf": abs(np.hypot(*m) - R_IN), "in_surf": abs(np.hypot(*m) - R_OUT),
         "out_surf": abs(m[0])}
    return min(d, key=d.get)

# --- extrusion ---------------------------------------------------------------
nv = len(V2)
faces, cell_faces = [], [[] for _ in polys]
bcs = {"cyl_surf": [], "in_surf": [], "out_surf": []}
for c, poly in enumerate(polys):
    faces.append([v + 1 for v in poly[::-1]])          # bottom, outward -z
    cell_faces[c].append(len(faces))
    faces.append([v + 1 + nv for v in poly])           # top, outward +z
    cell_faces[c].append(len(faces))
for key, owners in edge.items():
    c, a, b = owners[0]                                 # a->b is c's CCW order
    faces.append([a + 1, b + 1, b + 1 + nv, a + 1 + nv])
    fid = len(faces)
    cell_faces[c].append(fid)
    if len(owners) == 2:
        cell_faces[owners[1][0]].append(-fid)
    elif len(owners) == 1:
        bcs[bc_of(a, b)].append(fid)
    else:
        sys.exit("non-manifold edge")

nodes = np.vstack([np.c_[V2, np.zeros(nv)], np.c_[V2, np.full(nv, HZ)]])

# --- CGNS output (pyCGNS): NGON_n faces, NFACE_n cells, FaceCenter BCs ------
import CGNS.PAT.cgnslib as CL
import CGNS.PAT.cgnskeywords as CK
import CGNS.MAP as CM


def csr(lists):
    off = np.cumsum([0] + [len(l) for l in lists]).astype(np.int64)
    return np.concatenate([np.asarray(l, dtype=np.int64) for l in lists]), off


nf, nc = len(faces), len(cell_faces)
T = CL.newCGNSTree()
base = CL.newBase(T, "Base", 3, 3)
zone = CL.newZone(base, "Zone", np.array([[len(nodes), nc, 0]], dtype=np.int64), CK.Unstructured_s)
gc = CL.newGridCoordinates(zone, CK.GridCoordinates_s)
for k, name in enumerate([CK.CoordinateX_s, CK.CoordinateY_s, CK.CoordinateZ_s]):
    CL.newDataArray(gc, name, np.ascontiguousarray(nodes[:, k]))
conn, off = csr(faces)
CL.newElements(zone, "Faces", CK.NGON_n_s, np.array([1, nf], dtype=np.int64), conn, off)
conn, off = csr(cell_faces)
CL.newElements(zone, "Cells", CK.NFACE_n_s, np.array([nf + 1, nf + nc], dtype=np.int64), conn, off)
zbc = CL.newZoneBC(zone)
for name, ids in bcs.items():
    bc = CL.newBoundary(zbc, name, np.array([ids], dtype=np.int64), "BCWall",
                        pttype=CK.PointList_s)
    CL.newGridLocation(bc, CK.FaceCenter_s)
CM.save(out, T)

nsides = np.array([len(p) for p in polys])
print(f"{ns} cells, {nv} 2D vertices, {len(faces)} faces; polygon sides "
      f"min/mean/max {nsides.min()}/{nsides.mean():.2f}/{nsides.max()}; "
      + ", ".join(f"{k}: {len(v)}" for k, v in bcs.items()))
