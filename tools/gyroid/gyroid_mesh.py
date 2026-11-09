#!/usr/bin/env python3
"""
Conformal 3-region tetrahedral mesh of a TPMS (gyroid) sheet-solid in a cube.

The implicit function
    f(x,y,z) = sin(kx)cos(ky) + sin(ky)cos(kz) + sin(kz)cos(kx),   k = 2*pi/period
splits the cube [0,L]^3 into three regions

    airA  :  f >  +delta          (first air labyrinth)
    solid : -delta < f <= +delta  (the sheet, thickness ~ 2*delta/|grad f|)
    airB  :  f <= -delta          (second air labyrinth)

separated by the two offset surfaces f = +delta and f = -delta.

Why not a pure .geo file
------------------------
gmsh has no implicit-surface mesher and its CAD kernels cannot build a TPMS.
So this script builds the *boundary* of the three regions and gmsh only does
the 3D filling (which is what it is good at):

 1. structured node grid over the cube; every voxel is split into 6 tetrahedra
    (Kuhn / Freudenthal subdivision).  The triangulation it induces on a voxel
    face depends only on that face, so it is identical seen from either side:
    the result is watertight and, unlike marching cubes, there is no ambiguous
    case to disambiguate.
 2. marching tetrahedra at level +delta and -delta -> the two nappes.
 3. the 6 cube faces are cut by the same two level sets, reusing the *same*
    grid-edge intersection nodes, so the trace of the nappes on the cube faces
    matches the face triangulation exactly.
 4. everything is written as one .msh with a single node list, hence conformal
    by construction (no node welding, no tolerance).  Connected components of
    each region are detected and each one becomes its own gmsh Volume -- gmsh
    silently drops all but one component if several are put in the same
    Surface Loop, which would lose part of the solid.
 5. the generated .geo defines the Surface Loops, the Volumes, and the
    Physical groups; gmsh -3 fills them with tetrahedra.

--snap moves node values that lie within a fraction of the per-edge field
variation onto the level set.  Degenerate cuts then collapse onto existing grid
nodes instead of producing sliver triangles, which is what kills tet quality.

Examples
--------
  ./gyroid_mesh.py --nper 1 --n 48 --delta 0.35 -o gyroid --check
  ./gyroid_mesh.py --nper 2 --n 96 --delta 0.5 --size 0.05 -o gyroid2
  ./gyroid_mesh.py --n 32 --delta 0.3 --no-mesh          # surfaces only
"""

import argparse
import collections
import os
import subprocess
import sys
import time

import numpy as np


# --------------------------------------------------------------------------
# implicit functions
# --------------------------------------------------------------------------
def gyroid(x, y, z, period):
    k = 2.0 * np.pi / period
    return (np.sin(k * x) * np.cos(k * y)
            + np.sin(k * y) * np.cos(k * z)
            + np.sin(k * z) * np.cos(k * x))


def schwarz_p(x, y, z, period):
    k = 2.0 * np.pi / period
    return np.cos(k * x) + np.cos(k * y) + np.cos(k * z)


def diamond(x, y, z, period):
    k = 2.0 * np.pi / period
    sx, sy, sz = np.sin(k * x), np.sin(k * y), np.sin(k * z)
    cx, cy, cz = np.cos(k * x), np.cos(k * y), np.cos(k * z)
    return sx * sy * sz + sx * cy * cz + cx * sy * cz + cx * cy * sz


SURFACES = {"gyroid": gyroid, "schwarz_p": schwarz_p, "diamond": diamond}

REGIONS = ("solid", "airA", "airB")
OUT = "OUT"


# --------------------------------------------------------------------------
# global node pool: grid nodes + level-set cut nodes sitting on grid edges
# --------------------------------------------------------------------------
class NodePool:
    """Grid nodes [0, nbase) followed by cut nodes created on grid edges.

    A cut node is keyed by (grid edge, level), so tetrahedra and cube-face
    triangles always share it -> the complex is conformal.  When the
    intersection falls exactly on a grid node (t = 0 or 1, which --snap
    deliberately provokes) that grid node is returned instead of a duplicate;
    the resulting degenerate triangles are dropped afterwards.
    """

    def __init__(self, P, F):
        self.Pb = P                 # (nbase,3) grid coordinates
        self.Fb = F                 # (nbase,)  grid field values
        self.nbase = P.shape[0]
        self.count = self.nbase
        self.eP = []                # extra node coordinates
        self.eF = []                # extra node field value (= its level)
        self.eH = []                # extra node host grid edge (a,b)
        self.cache = {}             # level -> {edge key: node id}

    # -- accessors ---------------------------------------------------------
    def val(self, i):
        return self.Fb[i] if i < self.nbase else self.eF[i - self.nbase]

    def host(self, i):
        """(node,) for a grid node, (a,b) for a node sitting on grid edge a-b"""
        return (i,) if i < self.nbase else self.eH[i - self.nbase]

    def coords(self):
        if self.eP:
            return np.vstack([self.Pb, np.array(self.eP, dtype=float)])
        return self.Pb.copy()

    # -- creation ----------------------------------------------------------
    def _make(self, lo, hi, lev, cache):
        Fa, Fb = self.Fb[lo], self.Fb[hi]
        t = (lev - Fa) / (Fb - Fa)
        if t <= 0.0:
            nid = int(lo)
        elif t >= 1.0:
            nid = int(hi)
        else:
            nid = self.count
            self.count += 1
            self.eP.append(tuple(self.Pb[lo] + t * (self.Pb[hi] - self.Pb[lo])))
            self.eF.append(lev)
            self.eH.append((int(lo), int(hi)))
        cache[lo * self.nbase + hi] = nid
        return nid

    def get_one(self, a, b, lev):
        lo, hi = (a, b) if a < b else (b, a)
        cache = self.cache.setdefault(lev, {})
        nid = cache.get(lo * self.nbase + hi)
        return nid if nid is not None else self._make(lo, hi, lev, cache)

    def get_many(self, a, b, lev):
        a = np.asarray(a, dtype=np.int64)
        b = np.asarray(b, dtype=np.int64)
        key = (np.minimum(a, b) * self.nbase + np.maximum(a, b))
        uk, inv = np.unique(key.ravel(), return_inverse=True)
        cache = self.cache.setdefault(lev, {})
        ids = np.empty(uk.size, dtype=np.int64)
        for i, kk in enumerate(uk):
            kk = int(kk)
            nid = cache.get(kk)
            if nid is None:
                nid = self._make(kk // self.nbase, kk % self.nbase, lev, cache)
            ids[i] = nid
        return ids[inv].reshape(key.shape)

    # -- convex clipping of a polygon by a level set -----------------------
    def seg_edge(self, u, v):
        """grid edge hosting segment u-v, or None if the segment is a chord"""
        hu, hv = self.host(u), self.host(v)
        if len(hu) == 1 and len(hv) == 1:
            return (hu[0], hv[0])
        if len(hu) == 2 and len(hv) == 2:
            return hu if (hu[0] == hv[0] and hu[1] == hv[1]) else None
        if len(hu) == 1:
            return hv if hu[0] in hv else None
        return hu if hv[0] in hu else None

    def clip(self, poly, lev, keep_above):
        """Sutherland-Hodgman.  Tie break: val == lev counts as *below*."""
        if len(poly) < 3:
            return []
        out = []
        m = len(poly)
        for i in range(m):
            u = poly[i]
            v = poly[(i + 1) % m]
            iu = (self.val(u) > lev) if keep_above else (self.val(u) <= lev)
            iv = (self.val(v) > lev) if keep_above else (self.val(v) <= lev)
            if iu:
                out.append(u)
            if iu != iv:
                e = self.seg_edge(u, v)
                if e is None:
                    # the chord is an isoline of the other level, it cannot be
                    # crossed by this one; reaching this point means a bug
                    raise RuntimeError("clip: straddling chord with no host edge")
                out.append(self.get_one(e[0], e[1], lev))
        clean = []
        for x in out:
            if not clean or clean[-1] != x:
                clean.append(x)
        while len(clean) > 1 and clean[0] == clean[-1]:
            clean.pop()
        return clean if len(clean) >= 3 else []


# --------------------------------------------------------------------------
# geometry helpers
# --------------------------------------------------------------------------
def drop_degenerate(tris, P, scale):
    if len(tris) == 0:
        return np.zeros((0, 3), dtype=np.int64)
    t = np.asarray(tris, dtype=np.int64).reshape(-1, 3)
    t = t[(t[:, 0] != t[:, 1]) & (t[:, 1] != t[:, 2]) & (t[:, 0] != t[:, 2])]
    if len(t) == 0:
        return t
    n = np.cross(P[t[:, 1]] - P[t[:, 0]], P[t[:, 2]] - P[t[:, 0]])
    return t[0.5 * np.linalg.norm(n, axis=1) > 1e-10 * scale * scale]


def orient_towards(tris, P, ref):
    n = np.cross(P[tris[:, 1]] - P[tris[:, 0]], P[tris[:, 2]] - P[tris[:, 0]])
    bad = np.einsum('ij,ij->i', n, ref - P[tris[:, 0]]) < 0.0
    tris[bad] = tris[bad][:, [0, 2, 1]]
    return tris


def shell_volume(tris, P):
    """volume enclosed by an oriented closed triangulation (outward normals)"""
    a, b, c = P[tris[:, 0]], P[tris[:, 1]], P[tris[:, 2]]
    return float(np.einsum('ij,ij->i', a, np.cross(b, c)).sum() / 6.0)


def edge_pairs(tris):
    """(key-sorted triangle pairs sharing an edge, number of unpaired edges)"""
    nt = len(tris)
    M = int(tris.max()) + 1
    e = np.vstack([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
    tid = np.tile(np.arange(nt, dtype=np.int64), 3)
    key = (np.minimum(e[:, 0], e[:, 1]).astype(np.int64) * M
           + np.maximum(e[:, 0], e[:, 1]))
    o = np.argsort(key, kind='stable')
    key, tid = key[o], tid[o]
    _, cnt = np.unique(key, return_counts=True)
    same = key[1:] == key[:-1]
    pairs = np.stack([tid[:-1][same], tid[1:][same]], axis=1)
    return pairs, int((cnt != 2).sum())


def shell_components(tris):
    """(number of components, per-triangle component label)"""
    nt = len(tris)
    pairs, nbad = edge_pairs(tris)
    parent = np.arange(nt, dtype=np.int64)

    def find(x):
        r = x
        while parent[r] != r:
            r = parent[r]
        while parent[x] != r:
            parent[x], x = r, parent[x]
        return r

    for a, b in pairs:
        ra, rb = find(int(a)), find(int(b))
        if ra != rb:
            parent[ra] = rb
    roots = np.array([find(x) for x in range(nt)])
    uniq, lab = np.unique(roots, return_inverse=True)
    return len(uniq), lab.reshape(-1), nbad


# --------------------------------------------------------------------------
# construction of the boundary complex
# --------------------------------------------------------------------------
FACES = [(0, 0, 'xmin'), (0, 1, 'xmax'),
         (1, 0, 'ymin'), (1, 1, 'ymax'),
         (2, 0, 'zmin'), (2, 1, 'zmax')]


def build(args):
    period = args.period
    L = args.nper * period
    n = args.n
    N = n + 1
    dp, dm = args.delta, -args.delta
    fun = SURFACES[args.surface]
    h = L / n

    t0 = time.time()
    xs = np.linspace(0.0, L, N)
    Xg, Yg, Zg = np.meshgrid(xs, xs, xs, indexing='ij')
    F = np.ascontiguousarray(fun(Xg, Yg, Zg, period).ravel())   # (i*N+j)*N+k
    P = np.stack([Xg.ravel(), Yg.ravel(), Zg.ravel()], axis=1)
    del Xg, Yg, Zg
    print("grid %d^3 nodes, h = %.4g, cube L = %.4g (%d period(s)), delta = %.4g"
          % (N, h, L, args.nper, args.delta))

    # ---- snapping: remove the degenerate cuts at their source ------------
    var = max(np.abs(np.diff(F.reshape(N, N, N), axis=a)).max()
              for a in (0, 1, 2))
    thr = args.snap * var
    if thr > 0.0:
        d1, d2 = np.abs(F - dp), np.abs(F - dm)
        m1 = (d1 < thr) & (d1 <= d2)
        m2 = (d2 < thr) & (d2 < d1)
        F[m1], F[m2] = dp, dm
        print("snap: |f - level| < %.3e (%.0f%% of the max per-edge variation "
              "%.3e) -> %d + %d nodes moved onto the level sets"
              % (thr, 100 * args.snap, var, m1.sum(), m2.sum()))

    pool = NodePool(P, F)

    # ---- Kuhn subdivision: 6 tetrahedra per voxel ------------------------
    def off(a, b, c):
        return (a * N + b) * N + c

    ii, jj, kk = np.meshgrid(np.arange(n), np.arange(n), np.arange(n),
                             indexing='ij')
    base = ((ii * N + jj) * N + kk).ravel().astype(np.int64)
    del ii, jj, kk

    tets = []
    for pm in [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]:
        cur, nds = [0, 0, 0], [off(0, 0, 0)]
        for ax in pm:
            cur[ax] += 1
            nds.append(off(*cur))
        tets.append(np.array(nds, dtype=np.int64))

    # ---- marching tetrahedra --------------------------------------------
    OTH = np.array([[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]])
    raw = {dp: [], dm: []}
    for lev in (dp, dm):
        for shift in tets:
            tv = base[:, None] + shift[None, :]
            S = F[tv] > lev
            cnt = S.sum(axis=1)

            for target in (1, 3):
                sel = np.flatnonzero(cnt == target)
                if sel.size == 0:
                    continue
                tt, ss = tv[sel], S[sel]
                mino = ss if target == 1 else ~ss
                mi = np.argmax(mino, axis=1)
                ar = np.arange(len(tt))
                m = tt[ar, mi]
                o = np.take_along_axis(tt, OTH[mi], axis=1)
                c = pool.get_many(np.repeat(m, 3), o.ravel(), lev).reshape(-1, 3)
                w = ss.astype(float)
                above = (P[tt] * w[:, :, None]).sum(1) / w.sum(1)[:, None]
                raw[lev].append((c, above))

            sel = np.flatnonzero(cnt == 2)
            if sel.size:
                tt, ss = tv[sel], S[sel]
                order = np.argsort(~ss, axis=1, kind='stable')
                ar = np.arange(len(tt))
                a1, a2 = tt[ar, order[:, 0]], tt[ar, order[:, 1]]
                b1, b2 = tt[ar, order[:, 2]], tt[ar, order[:, 3]]
                q0 = pool.get_many(a1, b1, lev)
                q1 = pool.get_many(a1, b2, lev)
                q2 = pool.get_many(a2, b2, lev)
                q3 = pool.get_many(a2, b1, lev)
                above = 0.5 * (P[a1] + P[a2])
                raw[lev].append((np.stack([q0, q1, q2], 1), above))
                raw[lev].append((np.stack([q0, q2, q3], 1), above))
            del tv, S, cnt

    Pall = pool.coords()
    nappe = {}
    for lev in (dp, dm):
        tr = np.vstack([c for c, _ in raw[lev]])
        rf = np.vstack([a for _, a in raw[lev]])
        keep = ((tr[:, 0] != tr[:, 1]) & (tr[:, 1] != tr[:, 2])
                & (tr[:, 0] != tr[:, 2]))
        tr, rf = tr[keep], rf[keep]
        nrm = np.cross(Pall[tr[:, 1]] - Pall[tr[:, 0]],
                       Pall[tr[:, 2]] - Pall[tr[:, 0]])
        keep = 0.5 * np.linalg.norm(nrm, axis=1) > 1e-10 * h * h
        nappe[lev] = orient_towards(tr[keep], Pall, rf[keep])  # normal -> f up
        print("nappe f = %+.4g : %7d triangles" % (lev, keep.sum()))
    del raw

    # ---- how far are the nappe vertices from the true level set? --------
    #  distance ~ (f_true - level) / |grad f_true|, in units of h.  This is
    #  what --snap trades against tetrahedron quality, so measure it.
    eps = 1e-6 * period
    for lev in (dp, dm):
        v = np.unique(nappe[lev].ravel())
        q = Pall[v]
        ft = fun(q[:, 0], q[:, 1], q[:, 2], period)
        g = np.stack([(fun(q[:, 0] + (i == 0) * eps, q[:, 1] + (i == 1) * eps,
                           q[:, 2] + (i == 2) * eps, period)
                       - fun(q[:, 0] - (i == 0) * eps, q[:, 1] - (i == 1) * eps,
                             q[:, 2] - (i == 2) * eps, period)) / (2 * eps)
                      for i in range(3)], axis=1)
        d = np.abs(ft - lev) / np.maximum(np.linalg.norm(g, axis=1), 1e-30)
        print("    distance to the exact level set: rms %.3g h, max %.3g h"
              % (np.sqrt((d ** 2).mean()) / h, d.max() / h))

    # ---- the 6 cube faces, cut by the same two level sets ----------------
    # The Kuhn subdivision induces on every voxel face the diagonal joining
    # the (0,0) and (1,1) corners of the two in-plane axes taken in increasing
    # order -- hence the templates below.
    face_tris = {}
    for d, side, fname in FACES:
        u_ax, v_ax = [a for a in (0, 1, 2) if a != d]
        fixed = 0 if side == 0 else n
        sgn = -1.0 if side == 0 else 1.0
        idx = [0, 0, 0]

        def gid(iu, iv):
            idx[d], idx[u_ax], idx[v_ax] = fixed, iu, iv
            return off(*idx)

        oriented = []
        for tri in (((0, 0), (1, 0), (1, 1)), ((0, 0), (0, 1), (1, 1))):
            a, b, c = [np.array(t, dtype=float) for t in tri]
            e1, e2 = np.zeros(3), np.zeros(3)
            e1[u_ax], e1[v_ax] = (b - a)
            e2[u_ax], e2[v_ax] = (c - a)
            oriented.append(tri if np.cross(e1, e2)[d] * sgn > 0
                            else (tri[0], tri[2], tri[1]))

        for iu in range(n):
            for iv in range(n):
                for tri in oriented:
                    poly = [gid(iu + o[0], iv + o[1]) for o in tri]
                    for rname, pg in (
                            ('solid', pool.clip(pool.clip(poly, dp, False),
                                                dm, True)),
                            ('airA', pool.clip(poly, dp, True)),
                            ('airB', pool.clip(poly, dm, False))):
                        if len(pg) < 3:
                            continue
                        lst = face_tris.setdefault((fname, rname), [])
                        for q in range(1, len(pg) - 1):
                            lst.append((pg[0], pg[q], pg[q + 1]))

    Pall = pool.coords()
    for key in list(face_tris):
        t = drop_degenerate(face_tris[key], Pall, h)
        if len(t):
            face_tris[key] = t
        else:
            del face_tris[key]
    print("cube faces : %7d triangles in %d patches"
          % (sum(len(t) for t in face_tris.values()), len(face_tris)))

    # ---- items: oriented patches, with the region on each side -----------
    items = [dict(label='sheet_plus', face=None, tris=nappe[dp],
                  neg='solid', pos='airA'),
             dict(label='sheet_minus', face=None, tris=nappe[dm],
                  neg='airB', pos='solid')]
    for _, _, fname in FACES:
        for rname in REGIONS:
            if (fname, rname) in face_tris:
                items.append(dict(label='%s_%s' % (fname, rname), face=fname,
                                  tris=face_tris[(fname, rname)],
                                  neg=rname, pos=OUT))
    print("surface complex built in %.1f s (%d nodes)"
          % (time.time() - t0, Pall.shape[0]))
    return Pall, items, L, h


# --------------------------------------------------------------------------
# split each region into its connected components
# --------------------------------------------------------------------------
def partition(Pall, items, L):
    def outward(it, R):
        return it['tris'] if it['neg'] == R else it['tris'][:, [0, 2, 1]]

    labels = {}          # region -> {item index: per-triangle component label}
    ncomp = {}
    volumes = {}         # (region, comp) -> enclosed volume
    ok = True
    print("\nregion shells:")
    for R in REGIONS:
        parts, owner = [], []
        for i, it in enumerate(items):
            if R in (it['neg'], it['pos']):
                parts.append(outward(it, R))
                owner.append((i, len(it['tris'])))
        if not parts:
            ncomp[R], labels[R] = 0, {}
            print("  %-6s empty" % R)
            continue
        tr = np.vstack(parts)
        nc, lab, nbad = shell_components(tr)
        if nbad:
            print("  !! %s: %d edges not shared by exactly two triangles "
                  "(shell not closed)" % (R, nbad))
            ok = False
        d, o = {}, 0
        for i, m in owner:
            d[i] = lab[o:o + m]
            o += m
        labels[R], ncomp[R] = d, nc
        tot = 0.0
        for c in range(nc):
            v = shell_volume(tr[lab == c], Pall)
            volumes[(R, c)] = v
            tot += v
            if v <= 0.0:
                print("  !! %s component %d has non-positive enclosed volume "
                      "%.3e: it is a cavity inside another region, which this "
                      "script does not handle -- inspect the geometry" % (R, c, v))
                ok = False
        print("  %-6s %8d tri  %2d component(s)  volume %.6g (%.3f%% of L^3)"
              % (R, len(tr), nc, tot, 100 * tot / L ** 3))
        if nc > 1:
            sizes = " ".join("%.3g" % volumes[(R, c)] for c in range(nc))
            print("           component volumes: %s" % sizes)

    tot = sum(volumes.values())
    print("  total volume %.9g  (L^3 = %.9g, closure error %.2e)"
          % (tot, L ** 3, abs(tot - L ** 3) / L ** 3))
    if abs(tot - L ** 3) > 1e-8 * L ** 3:
        print("  !! the three regions do not tile the cube -- do not mesh this")
        ok = False

    # ---- surface entities: one per (patch, component pair) ---------------
    ents = []            # dict(tag, tris, item, label, face, neg, pos, cneg, cpos)
    for i, it in enumerate(items):
        m = len(it['tris'])
        z = np.zeros(m, dtype=np.int64)
        ln = labels[it['neg']][i] if it['neg'] in REGIONS else z
        lp = labels[it['pos']][i] if it['pos'] in REGIONS else z
        key = ln * (max(ncomp.values()) + 1) + lp
        for k in np.unique(key):
            sel = key == k
            ents.append(dict(tag=len(ents) + 1, tris=it['tris'][sel],
                             label=it['label'], face=it['face'],
                             neg=it['neg'], pos=it['pos'],
                             cneg=int(ln[sel][0]), cpos=int(lp[sel][0])))
    for e in ents:
        suffix = ''
        if ncomp.get(e['neg'], 1) > 1 or ncomp.get(e['pos'], 1) > 1:
            suffix = '_c%d_%d' % (e['cneg'], e['cpos'])
        e['name'] = e['label'] + suffix
    print("  %d surface entities" % len(ents))

    # ---- one Volume per region component ---------------------------------
    vols = []            # dict(tag, region, comp, loop=[signed surface tags])
    for R in REGIONS:
        for c in range(ncomp[R]):
            loop = []
            for e in ents:
                if e['neg'] == R and e['cneg'] == c:
                    loop.append(e['tag'])
                elif e['pos'] == R and e['cpos'] == c:
                    loop.append(-e['tag'])
            vols.append(dict(tag=len(vols) + 1, region=R, comp=c, loop=loop))
    return ents, vols, ok


# --------------------------------------------------------------------------
# output
# --------------------------------------------------------------------------
def write_msh22_surface(path, P, ents):
    used = np.zeros(P.shape[0], dtype=bool)
    for e in ents:
        used[e['tris'].ravel()] = True
    idx = np.flatnonzero(used)
    ren = np.zeros(P.shape[0], dtype=np.int64)
    ren[idx] = np.arange(1, idx.size + 1)
    ne = sum(len(e['tris']) for e in ents)
    with open(path, 'w') as fh:
        fh.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        fh.write("$Nodes\n%d\n" % idx.size)
        fh.writelines("%d %.17g %.17g %.17g\n" % (m, P[i, 0], P[i, 1], P[i, 2])
                      for m, i in enumerate(idx, start=1))
        fh.write("$EndNodes\n$Elements\n%d\n" % ne)
        k = 0
        for e in ents:
            t = e['tag']
            for tr in ren[e['tris']]:
                k += 1
                fh.write("%d 2 2 %d %d %d %d %d\n"
                         % (k, t, t, tr[0], tr[1], tr[2]))
        fh.write("$EndElements\n")
    print("wrote %s  (%d nodes, %d triangles)" % (path, idx.size, ne))


def write_geo(path, surfmsh, ents, vols, size, algo3d):
    with open(path, 'w') as fh:
        fh.write('// gyroid sheet-solid + the two air labyrinths in a cube\n')
        fh.write('// generated by gyroid_mesh.py -- do not edit by hand\n')
        fh.write('// the boundary triangulation is imported as is; gmsh only\n')
        fh.write('// fills the volumes with tetrahedra\n\n')
        fh.write('Merge "%s";\n\n' % os.path.basename(surfmsh))
        fh.write('Geometry.Tolerance = 1e-12;\n')
        fh.write('Mesh.Algorithm3D = %d;\n' % algo3d)
        fh.write('Mesh.Optimize = 1;\n')
        fh.write('Mesh.MeshSizeFromPoints = 0;\n')
        fh.write('Mesh.MeshSizeFromCurvature = 0;\n')
        fh.write('Mesh.MeshSizeExtendFromBoundary = 1;\n')
        if size:
            fh.write('Mesh.MeshSizeMax = %.17g;\n' % size)
        fh.write('\n')
        for v in vols:
            fh.write('Surface Loop(%d) = {%s};   // %s component %d\n'
                     % (v['tag'], ', '.join(str(t) for t in v['loop']),
                        v['region'], v['comp']))
            fh.write('Volume(%d) = {%d};\n' % (v['tag'], v['tag']))
        fh.write('\n')
        for i, R in enumerate(REGIONS, start=1):
            tags = [str(v['tag']) for v in vols if v['region'] == R]
            if tags:
                fh.write('Physical Volume("%s", %d) = {%s};\n'
                         % (R, i, ', '.join(tags)))
        fh.write('\n')
        groups = {}
        for e in ents:
            key = e['face'] if e['face'] else e['label']
            groups.setdefault(key, []).append(str(e['tag']))
        for i, (key, tags) in enumerate(sorted(groups.items()), start=1):
            fh.write('Physical Surface("%s", %d) = {%s};\n'
                     % (key, i, ', '.join(tags)))
    print("wrote %s" % path)


# --------------------------------------------------------------------------
# a posteriori check of the produced volume mesh
# --------------------------------------------------------------------------
def read_msh22(path):
    """-> {node id: xyz}, {physical tag: [tet node quadruples]}"""
    nodes, tets = {}, {}
    with open(path) as fh:
        line = fh.readline()
        while line:
            if line.startswith('$Nodes'):
                for _ in range(int(fh.readline())):
                    p = fh.readline().split()
                    nodes[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
            elif line.startswith('$Elements'):
                for _ in range(int(fh.readline())):
                    p = fh.readline().split()
                    if p[1] != '4':
                        continue
                    nt = int(p[2])
                    tets.setdefault(int(p[3]), []).append(
                        tuple(int(q) for q in p[3 + nt:7 + nt]))
            line = fh.readline()
    return nodes, tets


def volume_per_group(nodes, tets):
    out = {}
    for ph, tt in tets.items():
        s = 0.0
        for a, b, c, d in tt:
            pa, pb, pc, pd = (np.array(nodes[a]), np.array(nodes[b]),
                              np.array(nodes[c]), np.array(nodes[d]))
            s += abs(np.dot(pb - pa, np.cross(pc - pa, pd - pa))) / 6.0
        out[ph] = s
    return out


def check_conformity(nodes, tets, L):
    """the interfaces must be face-to-face: same triangles seen from both sides"""
    bnd = {}
    for ph, tt in tets.items():
        cnt = collections.Counter()
        for t in tt:
            for f in ((t[0], t[1], t[2]), (t[0], t[1], t[3]),
                      (t[0], t[2], t[3]), (t[1], t[2], t[3])):
                cnt[tuple(sorted(f))] += 1
        bnd[ph] = {f for f, k in cnt.items() if k == 1}

    def on_box(f):
        for d in range(3):
            if (all(abs(nodes[v][d]) < 1e-9 for v in f)
                    or all(abs(nodes[v][d] - L) < 1e-9 for v in f)):
                return True
        return False

    print("\nconformity of the volume mesh:")
    keys = sorted(bnd)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            print("  physical %d / %d : %d shared faces"
                  % (keys[i], keys[j], len(bnd[keys[i]] & bnd[keys[j]])))
    orphan = 0
    for ph in keys:
        for f in bnd[ph]:
            if on_box(f):
                continue
            if not any(f in bnd[o] for o in keys if o != ph):
                orphan += 1
    pos = collections.Counter(tuple(round(x, 12) for x in v)
                              for v in nodes.values())
    dup = sum(k - 1 for k in pos.values() if k > 1)
    print("  %d nodes, %d duplicated positions, %d interior faces without a "
          "neighbour region" % (len(nodes), dup, orphan))
    if dup or orphan:
        print("  !! the mesh is NOT conformal across the interfaces")
        return False
    return True


def voxel_fractions(args):
    L, n = args.nper * args.period, 4 * args.n
    xs = (np.arange(n) + 0.5) * L / n
    cnt = dict.fromkeys(REGIONS, 0)
    tot = 0
    for x0 in np.array_split(xs, 16):
        Xg, Yg, Zg = np.meshgrid(x0, xs, xs, indexing='ij')
        f = SURFACES[args.surface](Xg, Yg, Zg, args.period)
        cnt['airA'] += int((f > args.delta).sum())
        cnt['airB'] += int((f <= -args.delta).sum())
        cnt['solid'] += int(((f > -args.delta) & (f <= args.delta)).sum())
        tot += f.size
    return {k: v / tot for k, v in cnt.items()}


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--surface', default='gyroid', choices=sorted(SURFACES))
    ap.add_argument('--period', type=float, default=1.0,
                    help='period of the TPMS (default 1)')
    ap.add_argument('--nper', type=int, default=1,
                    help='periods per direction (cube side = nper*period)')
    ap.add_argument('--n', type=int, default=48,
                    help='background grid cells per direction; this sets the '
                         'resolution of the two nappes')
    ap.add_argument('--delta', type=float, default=0.35,
                    help='half offset of the two level sets f = +/- delta')
    ap.add_argument('--snap', type=float, default=0.10,
                    help='node snapping, as a fraction of the max per-edge '
                         'field variation (0 disables it)')
    ap.add_argument('--size', type=float, default=0.0,
                    help='Mesh.MeshSizeMax for the volume filling '
                         '(0 = let gmsh extend the boundary sizes)')
    ap.add_argument('--algo3d', type=int, default=1,
                    help='gmsh Mesh.Algorithm3D (1 = Delaunay, 10 = HXT)')
    ap.add_argument('-o', '--out', default='gyroid', help='output prefix')
    ap.add_argument('--format', default='msh4', choices=['msh4', 'msh22'],
                    help='final mesh format (subfv reads msh 4.1)')
    ap.add_argument('--no-mesh', action='store_true',
                    help='write the boundary complex and the .geo, stop there')
    ap.add_argument('--check', action='store_true',
                    help='re-read the volume mesh and compare the volumes with '
                         'a fine voxel count')
    ap.add_argument('--gmsh', default='gmsh')
    args = ap.parse_args()

    Pall, items, L, h = build(args)
    ents, vols, ok = partition(Pall, items, L)

    surfmsh = args.out + '_surface.msh'
    geo = args.out + '.geo'
    print()
    write_msh22_surface(surfmsh, Pall, ents)
    write_geo(geo, surfmsh, ents, vols, args.size, args.algo3d)

    if not ok:
        print("\n!! the boundary complex failed its checks, see the messages "
              "above; not calling gmsh")
        return 2
    if args.no_mesh:
        return 0

    out = args.out + '.msh'
    fmt = 'msh4' if args.format == 'msh4' else 'msh22'
    cmd = [args.gmsh, geo, '-3', '-format', fmt, '-o', out]
    print("\n$ " + " ".join(cmd))
    t0 = time.time()
    r = subprocess.run(cmd)
    print("gmsh returned %d in %.1f s" % (r.returncode, time.time() - t0))
    if r.returncode != 0:
        return r.returncode

    if args.check:
        chk = out
        if fmt != 'msh22':
            # convert, do not mesh again
            chk = args.out + '_check22.msh'
            subprocess.run([args.gmsh, out, '-save', '-format', 'msh22',
                            '-o', chk], stdout=subprocess.DEVNULL)
        nodes, tets = read_msh22(chk)
        mv = volume_per_group(nodes, tets)
        ok_conf = check_conformity(nodes, tets, L)
        ref = voxel_fractions(args)
        print("\nmeshed volume per physical group vs fine voxel count:")
        tot = 0.0
        for i, R in enumerate(REGIONS, start=1):
            v = mv.get(i, 0.0)
            tot += v
            print("  %-6s %.6g  (%.3f%% of L^3, voxel reference %.3f%%)"
                  % (R, v, 100 * v / L ** 3, 100 * ref[R]))
        print("  total  %.9g  (L^3 = %.9g, relative error %.2e)"
              % (tot, L ** 3, abs(tot - L ** 3) / L ** 3))
        if abs(tot - L ** 3) > 1e-6 * L ** 3:
            print("  !! gmsh did not fill the whole cube: some volume is "
                  "missing, do not use this mesh as is")
            return 3
        if not ok_conf:
            return 4
    return 0


if __name__ == '__main__':
    sys.exit(main())
