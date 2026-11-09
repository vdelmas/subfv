"""Hybrid half-cylinder mesh: structured QUAD boundary layer on the wall + VORONOI
farfield, conformal at the interface, extruded by one layer in z, written as CGNS
(NGON_n/NFACE_n, named FaceCenter BCs). Needs numpy, scipy, pycgns. See --help.

Domain (matching half_cylinder_*_ns.geo): half annulus r in [R_IN, R_OUT], x <= 0,
i.e. theta in [pi/2, 3pi/2]. BCs surf_cyl (r=R_IN), surf_in (r=R_OUT), surf_out (x=0).
The z caps are left untagged, exactly as the .geo does (boundary_2d excludes them).

HOW THE INTERFACE IS MADE EXACTLY CONFORMAL
-------------------------------------------
The quad block's outer ring has nodes at (R_I, theta_k), theta_k = pi/2 + k*dth,
dth = pi/n_t, so its outer faces are the chords [theta_k, theta_k+1].

Put the first Voronoi seed ring at the FACE MIDPOINT angles phi_k = theta_k + dth/2,
and mirror each seed radially across the circle of radius

    R_m = R_I * cos(dth/2).

The bisector of a seed and its mirror is then the tangent line to that circle at
angle phi_k. Two neighbouring tangents (phi_k-1, phi_k) meet, by symmetry, on the
ray theta_k at radius R_m/cos(dth/2) = R_I -- exactly the quad node. The same
algebra puts the end vertices on x=0 at radius exactly R_I. So each interface
Voronoi cell has ONE straight edge coinciding with ONE quad face: one-to-one
conformity, no clipping and no hanging node.

This is exact in real arithmetic; the script measures the residual (snap distance)
and REFUSES to write if it is not at round-off level, rather than snapping a
mismatch out of sight.
"""
import argparse
import sys
import numpy as np
from scipy.spatial import Voronoi
import CGNS.PAT.cgnslib as CL
import CGNS.PAT.cgnskeywords as CK
import CGNS.MAP as CM

ap = argparse.ArgumentParser()
ap.add_argument("out")
ap.add_argument("--r-in", type=float, default=1.0, help="cylinder radius")
ap.add_argument("--r-bl", type=float, default=1.04, help="interface radius (quad/Voronoi)")
ap.add_argument("--r-out", type=float, default=1.7, help="farfield radius")
ap.add_argument("--nt", type=int, default=128, help="azimuthal cells over the half annulus")
ap.add_argument("--nr", type=int, default=16, help="radial cells in the quad BL block")
ap.add_argument("--h-wall", type=float, default=1e-4, help="first cell height at the wall")
ap.add_argument("--iface-refine", type=int, default=2, metavar="M",
                help="M Voronoi cells per quad face at the interface. The outer edge of "
                     "each outermost quad cell is subdivided into M faces (the cell "
                     "becomes a 3+M-gon, legitimate in NGON), which decouples the "
                     "farfield size from nt. M=1 = one-to-one.")
ap.add_argument("--h-far", type=float, default=0.0,
                help="Voronoi target size at r_out (0 = pi*r_out/(nt*sqrt2): the tri "
                     "mesh ls3 divided by sqrt2, which matches CELL DENSITY rather than "
                     "characteristic length, a Voronoi cell covering ~2x a triangle)")
ap.add_argument("--z", type=float, nargs=2, default=[0.0, 1e-3], help="z0 z1 of the single layer")
ap.add_argument("--jitter", type=float, default=0.25, help="seed perturbation, fraction of spacing")
ap.add_argument("--lloyd", type=int, default=15, help="Lloyd iterations on the free seeds")
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--snap-tol", type=float, default=1e-9,
                help="max tolerated interface mismatch before refusing to write")
a = ap.parse_args()

R_IN, R_I, R_OUT = a.r_in, a.r_bl, a.r_out
Z0, Z1 = a.z
n_t, n_r = a.nt, a.nr
if R_IN >= R_I or R_I >= R_OUT:
    sys.exit("need r_in < r_bl < r_out")
m_ref = max(a.iface_refine, 1)
n_f = n_t * m_ref                                     # Voronoi faces on the interface
dth = np.pi / n_t                                     # quad azimuthal step
dphi = np.pi / n_f                                    # interface face step
th_k = np.pi / 2 + np.arange(n_t + 1) * dth           # quad node angles (inner rings)
th_f = np.pi / 2 + np.arange(n_f + 1) * dphi          # interface node angles
R_m = R_I * np.cos(dphi / 2)                          # mirror circle (see docstring)
rng = np.random.default_rng(a.seed)

# --- radial grading of the quad block: geometric, first cell = h_wall --------
L = R_I - R_IN
if a.h_wall * n_r >= L:
    sys.exit("h_wall * nr >= BL thickness: no room to grade")


def total(q):
    return a.h_wall * n_r if abs(q - 1) < 1e-14 else a.h_wall * (q ** n_r - 1) / (q - 1)


lo, hi = 1.0 + 1e-12, 2.0
for _ in range(200):                                   # bisection on the growth ratio
    mid = 0.5 * (lo + hi)
    if total(mid) < L:
        lo = mid
    else:
        hi = mid
q = 0.5 * (lo + hi)
dr = a.h_wall * q ** np.arange(n_r)
dr *= L / dr.sum()                                     # exact fit to the BL thickness
r_i = np.concatenate([[R_IN], R_IN + np.cumsum(dr)])
r_i[-1] = R_I                                          # exact, no round-off drift

h_interface = R_I * dphi                               # Voronoi size at the interface
# default: the tri mesh's law (ls proportional to r, ls3 = pi*r_out/nt) divided by
# sqrt(2), so that the CELL COUNT matches rather than the characteristic length --
# at equal h a Voronoi cell covers ~0.866 h^2 against ~0.433 h^2 for a triangle.
h_far = a.h_far if a.h_far > 0 else np.pi * R_OUT / (n_t * np.sqrt(2.0))


def hloc(r):
    return h_interface + (h_far - h_interface) * (r - R_I) / (R_OUT - R_I)


# --- Voronoi seeds ----------------------------------------------------------
# Ring 1 is PINNED at the face midpoint angles and NOT jittered: it is what makes
# the interface conformal. Its radius sits half a local size outside R_I.
phi = np.pi / 2 + (np.arange(n_f) + 0.5) * dphi
r1 = R_I + 0.5 * h_interface
ring1 = np.c_[r1 * np.cos(phi), r1 * np.sin(phi)]

free = []
r = r1 + h_interface
while r < R_OUT - 0.25 * hloc(r):
    h = hloc(r)
    nt = max(int(round(np.pi * r / h)), 2)
    t = np.pi / 2 + (np.arange(nt) + 0.5) * np.pi / nt
    rr = r + rng.uniform(-a.jitter, a.jitter, nt) * h
    t = t + rng.uniform(-a.jitter, a.jitter, nt) * h / r
    rr = np.clip(rr, R_I + 0.6 * h_interface, R_OUT - 0.15 * h)
    t = np.clip(t, np.pi / 2 + 0.15 * h / r, 1.5 * np.pi - 0.15 * h / r)
    free += list(zip(rr * np.cos(t), rr * np.sin(t)))
    r += h
free = np.array(free) if free else np.zeros((0, 2))


def refl_r(p, R):
    """reflect radially across the circle of radius R"""
    rr = np.hypot(p[:, 0], p[:, 1])
    return p * ((2 * R - rr) / rr)[:, None]


def build(seeds):
    """Voronoi of `seeds` plus the ghost seeds that force edges onto the boundaries."""
    rad = np.hypot(seeds[:, 0], seeds[:, 1])
    hh = hloc(np.clip(rad, R_I, R_OUT))
    near_x = -seeds[:, 0] < 3 * hh
    near_out = R_OUT - rad < 3 * hh
    near_i = rad - R_I < 3 * hh
    mir = [seeds[near_x] * [-1, 1],                                  # flat ends, x=0
           refl_r(seeds[near_out], R_OUT),                           # farfield arc
           refl_r(seeds[near_i], R_m),                               # interface
           refl_r(seeds[near_x & near_out] * [-1, 1], R_OUT),
           refl_r(seeds[near_x & near_i] * [-1, 1], R_m)]
    return Voronoi(np.vstack([seeds] + [m for m in mir if len(m)]))


seeds = np.vstack([ring1, free])
n_ring1 = len(ring1)

# Lloyd relaxation on the FREE seeds only; ring 1 must stay pinned.
for _ in range(a.lloyd):
    v = build(seeds)
    new = seeds.copy()
    for i in range(n_ring1, len(seeds)):
        reg = v.regions[v.point_region[i]]
        if -1 in reg or not reg:
            continue
        p = v.vertices[reg]
        x, y = p[:, 0], p[:, 1]
        cr = x * np.roll(y, -1) - np.roll(x, -1) * y
        A = 0.5 * cr.sum()
        if abs(A) < 1e-14:
            continue
        new[i] = [((x + np.roll(x, -1)) * cr).sum() / (6 * A),
                  ((y + np.roll(y, -1)) * cr).sum() / (6 * A)]
    rr = np.hypot(new[:, 0], new[:, 1])
    keep = (rr > R_I + 0.35 * h_interface) & (rr < R_OUT - 0.1 * hloc(np.clip(rr, R_I, R_OUT)))
    keep[:n_ring1] = False
    seeds[keep] = new[keep]

vor = build(seeds)
ns = len(seeds)

# --- Voronoi polygons, near-duplicate vertices merged -----------------------
used = sorted({v for i in range(ns) for v in vor.regions[vor.point_region[i]]})
if -1 in used:
    sys.exit("unbounded interior Voronoi cell: seeds too close to a boundary")
V = vor.vertices
parent = {v: v for v in used}


def find(v):
    while parent[v] != v:
        parent[v] = parent[parent[v]]
        v = parent[v]
    return v


for (i, j), rv in zip(vor.ridge_points, vor.ridge_vertices):
    if i >= ns and j >= ns:
        continue
    v0, v1 = rv
    if v0 < 0 or v1 < 0:
        continue
    rloc = np.clip(np.hypot(*V[v0]), R_I, R_OUT)
    if np.linalg.norm(V[v0] - V[v1]) < 0.05 * hloc(rloc):
        parent[find(v0)] = find(v1)

roots = sorted({find(v) for v in used})
Vv = {rt: V[rt] for rt in roots}

# --- shared node table ------------------------------------------------------
nodes2d, index = [], {}


def node(p):
    key = (round(p[0], 10), round(p[1], 10))
    if key not in index:
        index[key] = len(nodes2d)
        nodes2d.append(np.array([p[0], p[1]]))
    return index[key]


# quad block first, so the interface ring exists before any Voronoi vertex snaps
qid = [[node((r * np.cos(t), r * np.sin(t))) for t in th_k] for r in r_i[:-1]]
qid.append([node((R_I * np.cos(t), R_I * np.sin(t))) for t in th_f])   # subdivided ring
iface_pts = np.array([[R_I * np.cos(t), R_I * np.sin(t)] for t in th_f])

# Voronoi vertices: snap onto an interface node when within tolerance
snap_max = 0.0
vmap = {}
for rt in roots:
    p = Vv[rt]
    d = np.hypot(*(iface_pts - p).T)
    k = int(np.argmin(d))
    if d[k] < 0.25 * h_interface:
        snap_max = max(snap_max, d[k])
        vmap[rt] = qid[-1][k]
    else:
        vmap[rt] = node(p)

if snap_max > a.snap_tol:
    sys.exit(f"interface mismatch {snap_max:.3e} > tol {a.snap_tol:.1e}: the conformal "
             f"construction did not hold, refusing to write a non-conformal mesh")

polys = []
for c in range(n_r):                                            # quad cells (CCW)
    for k in range(n_t):
        if c == n_r - 1:
            # outer edge subdivided into m_ref faces -> a (3+m_ref)-gon, valid in NGON
            outer = [qid[c + 1][k * m_ref + j] for j in range(m_ref + 1)]
            polys.append([qid[c][k]] + outer + [qid[c][k + 1]])
        else:
            polys.append([qid[c][k], qid[c + 1][k], qid[c + 1][k + 1], qid[c][k + 1]])
n_quad = len(polys)

for i in range(ns):                                             # Voronoi cells
    reg = [vmap[find(v)] for v in vor.regions[vor.point_region[i]]]
    poly = [v for k, v in enumerate(reg) if v != reg[k - 1]]
    if len(poly) < 3:
        sys.exit("degenerate Voronoi cell after vertex merging")
    P = np.array([nodes2d[v] for v in poly])
    x, y = P[:, 0], P[:, 1]
    if 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y) < 0:
        poly = poly[::-1]
    polys.append(poly)

V2 = np.array(nodes2d)

# --- 2D edges -> owner/neighbour --------------------------------------------
edge = {}
for c, poly in enumerate(polys):
    for k in range(len(poly)):
        u, v = poly[k], poly[(k + 1) % len(poly)]
        edge.setdefault((min(u, v), max(u, v)), []).append((c, u, v))


def bc_of(u, v):
    m = 0.5 * (V2[u] + V2[v])
    d = {"surf_cyl": abs(np.hypot(*m) - R_IN), "surf_in": abs(np.hypot(*m) - R_OUT),
         "surf_out": abs(m[0])}
    return min(d, key=d.get)


# --- extrusion ---------------------------------------------------------------
nv = len(V2)
faces, cell_faces = [], [[] for _ in polys]
bcs = {"surf_cyl": [], "surf_in": [], "surf_out": []}
for c, poly in enumerate(polys):
    faces.append([v + 1 for v in poly[::-1]])
    cell_faces[c].append(len(faces))
    faces.append([v + 1 + nv for v in poly])
    cell_faces[c].append(len(faces))
n_iface = 0
for key, own in edge.items():
    c, u, v = own[0]
    faces.append([u + 1, v + 1, v + 1 + nv, u + 1 + nv])
    fid = len(faces)
    cell_faces[c].append(fid)
    if len(own) == 2:
        cell_faces[own[1][0]].append(-fid)
        if (own[0][0] < n_quad) != (own[1][0] < n_quad):
            n_iface += 1
    elif len(own) == 1:
        bcs[bc_of(u, v)].append(fid)
    else:
        sys.exit("non-manifold edge")

if n_iface != n_f:
    sys.exit(f"interface has {n_iface} shared faces, expected {n_f}: not conformal")

nodes = np.vstack([np.c_[V2, np.full(nv, Z0)], np.c_[V2, np.full(nv, Z1)]])

# --- CGNS output -------------------------------------------------------------


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
CM.save(a.out, T)

# --- report ------------------------------------------------------------------
area = 0.0
for poly in polys:
    P = V2[poly]
    x, y = P[:, 0], P[:, 1]
    area += 0.5 * abs(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))
exact = 0.5 * np.pi * (R_OUT ** 2 - R_IN ** 2)
nsides = np.array([len(p) for p in polys[n_quad:]])
print(f"quad BL: {n_quad} cells ({n_r} x {n_t}), growth ratio {q:.4f}, "
      f"first {dr[0]:.3e}, last {dr[-1]:.3e}, azimuthal at interface {h_interface:.3e}")
print(f"Voronoi: {ns} cells, sides min/mean/max {nsides.min()}/{nsides.mean():.2f}/{nsides.max()}")
print(f"interface: {n_iface} shared faces (expected {n_f}, {m_ref} per quad face), "
      f"max snap {snap_max:.2e}, Voronoi size there {h_interface:.4e}")
print(f"total {nc} cells, {len(nodes)} nodes, {nf} faces; area {area:.12f} (exact {exact:.12f})")
print("BCs: " + ", ".join(f"{k}: {len(v)}" for k, v in bcs.items()))
