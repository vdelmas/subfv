"""Voronoi mesh of a rectangle [x0,x1]x[y0,y1], extruded by one layer in z, written as CGNS
(NGON_n/NFACE_n, named FaceCenter BCs). Needs numpy, scipy, pycgns. See --help."""
import argparse
import numpy as np
from scipy.spatial import Voronoi
import CGNS.PAT.cgnslib as CL
import CGNS.PAT.cgnskeywords as CK
import CGNS.MAP as CM

ap = argparse.ArgumentParser()
ap.add_argument("out")
ap.add_argument("--n", type=int, default=50, help="seeds per unit length along x")
ap.add_argument("--box", type=float, nargs=4, default=[-0.5, 0.5, -0.5, 0.5], help="x0 x1 y0 y1")
ap.add_argument("--z", type=float, nargs=2, default=[-0.5, 0.5], help="z0 z1 of the single layer")
ap.add_argument("--jitter", type=float, default=0.3, help="seed perturbation, fraction of spacing")
ap.add_argument("--lloyd", type=int, default=0, help="Lloyd iterations (0 = raw jittered Voronoi)")
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--bc", default="all_surf", help="one BC name for the 4 sides")
ap.add_argument("--caps-bc", default="", help="BC name for the z caps ('' = untagged)")
ap.add_argument("--hex", action="store_true",
                help="regular honeycomb: seeds on a TRIANGULAR lattice, unperturbed. "
                     "The Voronoi of a triangular lattice is a hexagonal tiling, so the "
                     "interior cells are regular hexagons; only the boundary rows are cut "
                     "into pentagons/quads. Ignores --jitter and --lloyd.")
ap.add_argument("--rotate", type=float, default=0.0, metavar="DEG",
                help="rotate the --hex lattice by DEG degrees about the box centre. "
                     "Separates the two candidate sources of anisotropy in a square box: "
                     "if an artefact turns with the lattice it comes from the tiling, if it "
                     "stays aligned with the walls it does not. 30 deg maps a flat-top "
                     "honeycomb onto a pointy-top one.")
ap.add_argument("--periodic", action="store_true",
                help="centre-line-symmetric seeds: opposite sides match exactly (for periodic_mesh)")
a = ap.parse_args()

x0, x1, y0, y1 = a.box
h = (x1 - x0) / a.n
nx, ny = a.n, max(int(round((y1 - y0) / h)), 1)
rng = np.random.default_rng(a.seed)
if a.hex:
    # Triangular lattice: rows spaced dy = a*sqrt(3)/2, odd rows offset by a/2.
    # ny is rounded to fit the box, so the hexagons are stretched by the printed
    # aspect ratio (a fraction of a percent for a well-chosen --n).
    dx = (x1 - x0) / nx
    ny = max(int(round((y1 - y0) / (dx * np.sqrt(3.0) / 2.0))), 1)
    dy = (y1 - y0) / ny
    if a.rotate == 0.0:
        rows = []
        for j in range(ny):
            yy = y0 + (j + 0.5) * dy
            if j % 2 == 0:                   # nx seeds, centred in each column
                xx = x0 + (np.arange(nx) + 0.5) * dx
            else:                            # nx-1 seeds, offset by half a column
                xx = x0 + (np.arange(nx - 1) + 1.0) * dx
            rows.append(np.c_[xx, np.full(len(xx), yy)])
        seeds = np.vstack(rows)
    else:
        # Lay the lattice over a disc large enough to cover the rotated box, turn it,
        # then keep the seeds well inside the walls so the mirrors stay clean. The
        # boundary row is no longer commensurate with the walls, so the edge cells are
        # less regular than at rotate=0 -- that is the price of the test.
        th = np.radians(a.rotate)
        R = 0.75 * np.hypot(x1 - x0, y1 - y0) + 2.0 * dx
        xc, yc = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
        m = int(np.ceil(R / dx)) + 2
        n = int(np.ceil(R / dy)) + 2
        pts = []
        for j in range(-n, n + 1):
            off = 0.0 if j % 2 == 0 else 0.5 * dx
            pts.append(np.c_[xc + off + np.arange(-m, m + 1) * dx,
                             np.full(2 * m + 1, yc + j * dy)])
        q = np.vstack(pts) - [xc, yc]
        c, si = np.cos(th), np.sin(th)
        q = np.c_[q[:, 0] * c - q[:, 1] * si, q[:, 0] * si + q[:, 1] * c] + [xc, yc]
        keep = ((q[:, 0] > x0 + 0.30 * dx) & (q[:, 0] < x1 - 0.30 * dx) &
                (q[:, 1] > y0 + 0.30 * dy) & (q[:, 1] < y1 - 0.30 * dy))
        seeds = q[keep]
    print(f"honeycomb: {nx} columns x {ny} rows, spacing {dx:.6f}, row pitch {dy:.6f}, "
          f"hexagon aspect {dy / (dx * np.sqrt(3.0) / 2.0):.5f} (1 = perfectly regular)")
else:
    gx, gy = np.meshgrid(x0 + (np.arange(nx) + 0.5) * (x1 - x0) / nx, y0 + (np.arange(ny) + 0.5) * (y1 - y0) / ny)
    seeds = np.c_[gx.ravel(), gy.ravel()] + rng.uniform(-a.jitter, a.jitter, (nx * ny, 2)) * h
Lx, Ly = x1 - x0, y1 - y0
xm, ym = 0.5 * (x0 + x1), 0.5 * (y0 + y1)


def symmetrize(q):
    # quadrant seeds (x<xm, y<ym) reflected through both centre lines: the box-mirror
    # construction is then also periodic (left side = right side, bottom = top)
    q = np.c_[np.minimum(q[:, 0], xm - 0.05 * h), np.minimum(q[:, 1], ym - 0.05 * h)]
    return np.vstack([q, q * [-1, 1] + [2 * xm, 0], q * [1, -1] + [0, 2 * ym], -q + [2 * xm, 2 * ym]])


if a.hex:
    pass                                      # lattice already fits the box, no clipping
elif a.periodic:
    if nx % 2 or ny % 2:
        raise SystemExit("--periodic needs an even number of seeds per side")
    gxq, gyq = np.meshgrid(x0 + (np.arange(nx // 2) + 0.5) * Lx / nx, y0 + (np.arange(ny // 2) + 0.5) * Ly / ny)
    q = np.c_[gxq.ravel(), gyq.ravel()] + rng.uniform(-a.jitter, a.jitter, (gxq.size, 2)) * h
    q = np.c_[np.maximum(q[:, 0], x0 + 0.05 * h), np.maximum(q[:, 1], y0 + 0.05 * h)]
    seeds = symmetrize(q)
else:
    seeds[:, 0] = np.clip(seeds[:, 0], x0 + 0.05 * h, x1 - 0.05 * h)
    seeds[:, 1] = np.clip(seeds[:, 1], y0 + 0.05 * h, y1 - 0.05 * h)
ns = len(seeds)


def voronoi(s):
    # mirrors across the 4 sides and the 4 corners: cell edges land exactly on the box
    m = [s * [-1, 1] + [2 * x0, 0], s * [-1, 1] + [2 * x1, 0], s * [1, -1] + [0, 2 * y0], s * [1, -1] + [0, 2 * y1]]
    for cx in (x0, x1):
        for cy in (y0, y1):
            m.append(-s + [2 * cx, 2 * cy])
    return Voronoi(np.vstack([s] + m))


def polygon(vor, i):
    reg = vor.regions[vor.point_region[i]]
    if -1 in reg:
        raise SystemExit("unbounded interior cell")
    return reg


for _ in range(a.lloyd):
    vor = voronoi(seeds)
    for i in range(ns):
        p = vor.vertices[polygon(vor, i)]
        x, y = p[:, 0], p[:, 1]
        cr = x * np.roll(y, -1) - np.roll(x, -1) * y
        A = 0.5 * cr.sum()
        seeds[i] = [np.sum((x + np.roll(x, -1)) * cr) / (6 * A), np.sum((y + np.roll(y, -1)) * cr) / (6 * A)]
    if a.periodic:
        seeds = symmetrize(seeds[: ns // 4])

vor = voronoi(seeds)
V = vor.vertices
used = sorted({v for i in range(ns) for v in polygon(vor, i)})
# merge vertices closer than 5% of the spacing (slivers from near-cocircular seeds)
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

# Snap boundary vertices exactly onto the box. The vertex merge above can pull a
# vertex that belongs on a wall slightly inward -- on a regular lattice the wall
# carries exact 4-fold degeneracies (seed + mirror + two neighbours equidistant),
# so scipy emits near-duplicate vertices there and the merge shifts them. Left
# alone this loses volume (measured 8e-5 of the box on a 56x65 honeycomb) even
# though the mesh stays watertight and its bounding box stays exact.
_snap = 0.02 * h
for _k, (_lo, _hi) in enumerate(((x0, x1), (y0, y1))):
    V2[np.abs(V2[:, _k] - _lo) < _snap, _k] = _lo
    V2[np.abs(V2[:, _k] - _hi) < _snap, _k] = _hi
# snap boundary vertices exactly onto the box (mirrors put them there up to round-off)
for c, lo, hi in ((0, x0, x1), (1, y0, y1)):
    V2[np.abs(V2[:, c] - lo) < 1e-9 * h, c] = lo
    V2[np.abs(V2[:, c] - hi) < 1e-9 * h, c] = hi

polys = []
for i in range(ns):
    reg = [vid[find(v)] for v in polygon(vor, i)]
    poly = [v for k, v in enumerate(reg) if v != reg[k - 1]]
    x, y = V2[poly, 0], V2[poly, 1]
    if 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y) < 0:
        poly = poly[::-1]
    polys.append(poly)

if a.periodic:
    # opposite sides are mirror images through the centre lines: make them bit-identical
    for c, lo, hi in ((0, x0, x1), (1, y0, y1)):
        o = 1 - c
        il = np.where(V2[:, c] == lo)[0]; ir = np.where(V2[:, c] == hi)[0]
        il = il[np.argsort(V2[il, o])]; ir = ir[np.argsort(V2[ir, o])]
        if len(il) != len(ir) or np.abs(V2[il, o] - V2[ir, o]).max() > 1e-8 * h:
            raise SystemExit("periodic sides do not match")
        V2[ir, o] = V2[il, o]

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
elen = [np.linalg.norm(V2[p[k]] - V2[p[k - 1]]) for p in polys for k in range(len(p))]
print(f"min edge length {min(elen) / h:.3f} h")
print(f"{nc} cells, sides min/mean/max {sides.min()}/{sides.mean():.2f}/{sides.max()}, "
      f"area sum {areas.sum():.12f} (box {(x1 - x0) * (y1 - y0):.12f}), "
      f"cell area min/max ratio {areas.min() / areas.max():.3f}, {len(side_ids)} side faces")
