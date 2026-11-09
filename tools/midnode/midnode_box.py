"""Cartesian quad (2D, extruded one layer) or hex (3D) mesh with an extra node at the middle of every
edge (and, in 3D, at every face centre): degenerate octagons / 26-node hexes, written as CGNS
NGON_n/NFACE_n. Needs numpy, pycgns. See --help."""
import argparse
import numpy as np
import CGNS.PAT.cgnslib as CL
import CGNS.PAT.cgnskeywords as CK
import CGNS.MAP as CM

ap = argparse.ArgumentParser()
ap.add_argument("out")
ap.add_argument("--dim", type=int, choices=(2, 3), default=2, help="2 = quads extruded one layer, 3 = hexes")
ap.add_argument("--n", type=int, nargs="+", required=True, help="cells nx ny [nz] (nz only for --dim 3)")
ap.add_argument("--box", type=float, nargs="+", required=True, help="x0 x1 y0 y1 z0 z1 (2D: z0 z1 = extrusion)")
ap.add_argument("--plain", action="store_true", help="no extra nodes (reference quad/hex mesh, same CGNS path)")
ap.add_argument("--bc", default="wall", help="BC name for every boundary side")
ap.add_argument("--side", action="append", default=[], help="override one side, e.g. xmin=left_surf "
                "(sides: xmin xmax ymin ymax zmin zmax; in 2D zmin/zmax are the caps)")
ap.add_argument("--no-caps-bc", action="store_true", help="2D: leave the z caps untagged")
a = ap.parse_args()

nx, ny = a.n[0], a.n[1]
nz = a.n[2] if a.dim == 3 else 1
x0, x1, y0, y1, z0, z1 = a.box
step = 1 if a.plain else 2          # lattice points per cell edge
side_bc = {s: a.bc for s in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")}
for s in a.side:
    k, v = s.split("=")
    side_bc[k] = v

# --- nodes on the (half-)lattice: keep points with at most 2 odd indices (3D) / 1 odd index (2D plane)
LX, LY = step * nx + 1, step * ny + 1
LZ = step * nz + 1 if a.dim == 3 else 2   # 2D: two z planes only (no mid-height node)
nid = -np.ones((LX, LY, LZ), dtype=np.int64)
coords = []
for k in range(LZ):
    for j in range(LY):
        for i in range(LX):
            odd = (i % 2) + (j % 2) + ((k % 2) if a.dim == 3 else 0)
            if step == 2 and odd > (2 if a.dim == 3 else 1):
                continue
            nid[i, j, k] = len(coords) + 1
            z = z0 + (z1 - z0) * k / (LZ - 1)
            coords.append((x0 + (x1 - x0) * i / (LX - 1), y0 + (y1 - y0) * j / (LY - 1), z))
coords = np.array(coords)

faces, face_of, cell_faces, bc_faces = [], {}, [], {}


def add_face(cell, pts, centre, side):
    """pts: lattice index triples of a planar polygon; oriented outward w.r.t. the cell centre."""
    ids = [int(nid[p]) for p in pts]
    X = coords[np.array(ids) - 1]
    n = np.zeros(3)
    for m in range(len(X)):
        n += np.cross(X[m], X[(m + 1) % len(X)])
    if np.dot(n, X.mean(0) - centre) < 0:
        ids = ids[::-1]
    key = tuple(sorted(ids))
    if key in face_of:
        fid, owner = face_of[key]
        assert owner != cell
        cell_faces[cell].append(-fid)
        return
    faces.append(ids)
    fid = len(faces)
    face_of[key] = (fid, cell)
    cell_faces[cell].append(fid)
    if side:
        bc_faces.setdefault(side_bc[side], []).append(fid)


def ring(i0, j0, k):
    """CCW lattice points around the 2D cell with lower corner (i0, j0) in plane k."""
    s = step
    pts = [(i0 + t, j0, k) for t in range(s)] + [(i0 + s, j0 + t, k) for t in range(s)] + \
          [(i0 + s - t, j0 + s, k) for t in range(s)] + [(i0, j0 + s - t, k) for t in range(s)]
    return pts


def patches(fixed_axis, fixed, a0, b0):
    """Sub-quads of one hex face: lattice patch of size step x step in the two free axes."""
    out = []
    for da in range(step):
        for db in range(step):
            q = []
            for (u, v) in ((0, 0), (1, 0), (1, 1), (0, 1)):
                p = [0, 0, 0]
                free = [ax for ax in range(3) if ax != fixed_axis]
                p[fixed_axis] = fixed
                p[free[0]] = a0 + da + u
                p[free[1]] = b0 + db + v
                q.append(tuple(p))
            out.append(q)
    return out


for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            c = len(cell_faces)
            cell_faces.append([])
            I, J = step * i, step * j
            if a.dim == 2:
                centre = np.array([x0 + (i + .5) * (x1 - x0) / nx, y0 + (j + .5) * (y1 - y0) / ny, .5 * (z0 + z1)])
                caps = None if a.no_caps_bc else "zmin"
                add_face(c, ring(I, J, 0), centre, caps)
                add_face(c, ring(I, J, 1), centre, None if a.no_caps_bc else "zmax")
                r = ring(I, J, 0)
                for m in range(len(r)):
                    p, q = r[m], r[(m + 1) % len(r)]
                    side = None
                    if p[0] == q[0] == 0: side = "xmin"
                    elif p[0] == q[0] == LX - 1: side = "xmax"
                    elif p[1] == q[1] == 0: side = "ymin"
                    elif p[1] == q[1] == LY - 1: side = "ymax"
                    add_face(c, [p, q, (q[0], q[1], 1), (p[0], p[1], 1)], centre, side)
            else:
                K = step * k
                centre = np.array([x0 + (i + .5) * (x1 - x0) / nx, y0 + (j + .5) * (y1 - y0) / ny,
                                   z0 + (k + .5) * (z1 - z0) / nz])
                for axis, lo, hi, idx, names in ((0, I, I + step, (J, K), ("xmin", "xmax")),
                                                 (1, J, J + step, (I, K), ("ymin", "ymax")),
                                                 (2, K, K + step, (I, J), ("zmin", "zmax"))):
                    last = (LX, LY, LZ)[axis] - 1
                    for fixed, name in ((lo, names[0]), (hi, names[1])):
                        side = name if fixed in (0, last) else None
                        for q in patches(axis, fixed, *idx):
                            add_face(c, q, centre, side)

# --- checks: closed cells (sum of area vectors), positive volumes, interior faces shared twice
vols = np.zeros(len(cell_faces))
for c, cf in enumerate(cell_faces):
    s = np.zeros(3)
    for f in cf:
        ids = faces[abs(f) - 1] if f > 0 else faces[abs(f) - 1][::-1]
        X = coords[np.array(ids) - 1]
        for m in range(1, len(X) - 1):
            cr = np.cross(X[m] - X[0], X[m + 1] - X[0])
            s += 0.5 * cr
            vols[c] += np.dot(cr, X[0]) / 6.0
    assert np.linalg.norm(s) < 1e-12 * max(1.0, abs(vols[c]) ** (2 / 3)), "open cell"
assert vols.min() > 0, "negative volume"
nb = sum(len(v) for v in bc_faces.values())
n_owned = sum(1 for f in face_of.values())


def csr(lists):
    off = np.cumsum([0] + [len(l) for l in lists]).astype(np.int64)
    return np.concatenate([np.asarray(l, dtype=np.int64) for l in lists]), off


nf, nc = len(faces), len(cell_faces)
T = CL.newCGNSTree()
zone = CL.newZone(CL.newBase(T, "Base", 3, 3), "Zone", np.array([[len(coords), nc, 0]], dtype=np.int64),
                  CK.Unstructured_s)
gc = CL.newGridCoordinates(zone, CK.GridCoordinates_s)
for k, name in enumerate([CK.CoordinateX_s, CK.CoordinateY_s, CK.CoordinateZ_s]):
    CL.newDataArray(gc, name, np.ascontiguousarray(coords[:, k]))
conn, off = csr(faces)
CL.newElements(zone, "Faces", CK.NGON_n_s, np.array([1, nf], dtype=np.int64), conn, off)
conn, off = csr(cell_faces)
CL.newElements(zone, "Cells", CK.NFACE_n_s, np.array([nf + 1, nf + nc], dtype=np.int64), conn, off)
zbc = CL.newZoneBC(zone)
for name, ids in bc_faces.items():
    bc = CL.newBoundary(zbc, name, np.array([ids], dtype=np.int64), "BCWall", pttype=CK.PointList_s)
    CL.newGridLocation(bc, CK.FaceCenter_s)
CM.save(a.out, T)

box = (x1 - x0) * (y1 - y0) * (z1 - z0)
npc = sorted({len({v for f in cf for v in faces[abs(f) - 1]}) for cf in cell_faces})
print(f"{nc} cells, {len(coords)} nodes, {nf} faces ({nb} boundary), nodes per cell {npc}, "
      f"faces per cell {sorted({len(cf) for cf in cell_faces})}; volume {vols.sum():.12g} (box {box:.12g}); "
      + ", ".join(f"{k}: {len(v)}" for k, v in bc_faces.items()))
