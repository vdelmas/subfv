#!/usr/bin/env pvpython
"""Mesh images (surface with edges) and the geometric properties that matter for these cases.

For a hypersonic half cylinder the numbers that decide what the solver can resolve are the
wall-normal size of the FIRST cell (the boundary layer is only captured if that is small
enough), how fast the spacing grows away from the wall, and the tangential resolution along
the shock. So rather than a plain cell count this reports, per mesh:

  cells / vertices / cell type
  r_wall, r_outer          -- the O-grid's radial extent
  first cell height        -- 2 * (min cell-centroid radius - r_wall), i.e. the wall-normal
                              size of the wall-adjacent cell
  last cell height         -- same at the outer boundary, giving the growth over the layer
  tangential size at wall  -- r_wall * (angular span / number of tangential cells)
  aspect ratio at wall     -- tangential / wall-normal, the stretching the solver must cope with

usage: pvpython mesh_report.py <out_dir> <label>:<file.pvtu> [<label>:<file.pvtu> ...]
"""
import sys, os
import numpy as np
from paraview.simple import *
from paraview import servermanager as sm
from paraview.numpy_support import vtk_to_numpy

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)

H_PX = 1000   # the camera is fitted per mesh: the Euler O-grid runs out to r=2.5 and the
              # viscous one only to r=1.7, so a single hard-coded framing would crop one of them


def blocks_of(reader):
    out = []

    def walk(o):
        if o.IsA('vtkCompositeDataSet'):
            it = o.NewIterator(); it.InitTraversal()
            while not it.IsDoneWithTraversal():
                walk(it.GetCurrentDataObject()); it.GoToNextItem()
        else:
            out.append(o)
    walk(sm.Fetch(reader))
    return out


rows = []
for spec in sys.argv[2:]:
    label, path = spec.split(':', 1)
    r = XMLPartitionedUnstructuredGridReader(FileName=[path]); r.UpdatePipeline()
    bl = blocks_of(r)
    cen = np.concatenate([vtk_to_numpy(b.GetCellData().GetArray('Centroid')) for b in bl])
    ncell = cen.shape[0]
    nvert = sum(b.GetNumberOfPoints() for b in bl)

    # Everything below is measured from the actual cell geometry rather than from centroid
    # "shells": subfvns writes its cells as VTK polyhedra and the tri meshes are unstructured,
    # so any layer-detection heuristic on centroid radii gives nonsense on 3 of the 4 meshes.
    cell_r = []      # (r_min, r_max, point thetas, point radii) per cell
    n_by_type = {}
    allpts = np.concatenate([vtk_to_numpy(b.GetPoints().GetData()) for b in bl])
    xspan = allpts[:, 0].max() - allpts[:, 0].min()
    yspan = allpts[:, 1].max() - allpts[:, 1].min()
    for b in bl:
        pts = vtk_to_numpy(b.GetPoints().GetData())
        pr = np.hypot(pts[:, 0], pts[:, 1])
        pth = np.arctan2(pts[:, 1], pts[:, 0])
        for c in range(b.GetNumberOfCells()):
            cell = b.GetCell(c)
            n_by_type[cell.GetCellType()] = n_by_type.get(cell.GetCellType(), 0) + 1
            ids = [cell.GetPointId(k) for k in range(cell.GetNumberOfPoints())]
            rr = pr[ids]
            cell_r.append((rr.min(), rr.max(), pth[ids], rr))
    r_wall = min(c[0] for c in cell_r)
    r_out = max(c[1] for c in cell_r)
    tol = 1e-9 + 1e-6 * (r_out - r_wall)
    wall_cells = [c for c in cell_r if c[0] <= r_wall + tol]

    # Wall-normal size of the wall-adjacent cell: its radial extent. This is the number the
    # boundary layer lives or dies on, and it is what Beta_HWall in the .geo targets.
    h_first = float(np.median([c[1] - c[0] for c in wall_cells]))

    # Tangential size at the wall: arc spanned by the cell's points that actually sit ON the
    # wall. Cells touching the wall at a single vertex (common on the triangular meshes)
    # span nothing and would drag a median to zero, so they are excluded rather than counted.
    def arc(c):
        th = c[2][c[3] <= r_wall + tol]
        if th.size < 2:
            return np.nan
        th = np.unwrap(np.sort(th))
        return r_wall * (th[-1] - th[0])
    arcs = np.array([arc(c) for c in wall_cells], dtype=float)
    arcs = arcs[np.isfinite(arcs) & (arcs > 0)]
    h_tan = float(np.median(arcs)) if arcs.size else float('nan')

    # No outer-layer or growth-ratio figure is reported: the viscous meshes' outer boundary is
    # NOT a circle (measured x in [-1.7, 0] but y in [-3.575, 3.575]), so "the last radial
    # cell" is not a well-defined single layer and any growth ratio built on it is meaningless.
    npc = int(np.median([len(c[3]) for c in cell_r]))
    ctype = {8: 'hexahedron', 6: 'prism', 4: 'tet'}.get(npc, '%d-point cell' % npc)
    n_wall = int(arcs.size)
    rows.append((label, ncell, nvert, ctype, r_wall, xspan, yspan, n_wall,
                 h_first, h_tan, h_tan / h_first if h_first else float('nan')))

    # --- mesh image, camera fitted to this mesh's own bounds ---
    xmin, xmax = cen[:, 0].min(), cen[:, 0].max()
    ymin, ymax = cen[:, 1].min(), cen[:, 1].max()
    pad = 0.04 * max(xmax - xmin, ymax - ymin)
    xmin -= pad; xmax += pad; ymin -= pad; ymax += pad
    y_half = 0.5 * (ymax - ymin)
    w_px = max(1, int(round((xmax - xmin) / (ymax - ymin) * H_PX)))
    view = CreateView("RenderView")
    view.ViewSize = [w_px, H_PX]
    view.CameraParallelProjection = 1
    view.OrientationAxesVisibility = 0
    view.UseColorPaletteForBackground = 0
    view.Background = [1, 1, 1]; view.Background2 = [1, 1, 1]
    d = Show(r, view)
    d.Representation = "Surface With Edges"
    d.DisableLighting = 1
    ColorBy(d, ('POINTS', ''))   # solid colour: the mesh lines are the subject, not a field
    d.AmbientColor = [1, 1, 1]; d.DiffuseColor = [1, 1, 1]
    d.EdgeColor = [0.1, 0.1, 0.1]
    d.LineWidth = 0.4
    d.SetScalarBarVisibility(view, False)
    Render()
    cx, cy = 0.5 * (xmin + xmax), 0.5 * (ymin + ymax)
    view.CameraFocalPoint = [cx, cy, 0.0]
    view.CameraPosition = [cx, cy, 10.0]
    view.CameraViewUp = [0, 1, 0]
    view.CameraParallelScale = y_half
    png = os.path.join(OUT, "mesh_%s.png" % label)
    SaveScreenshot(png, view, ImageResolution=[w_px, H_PX],
                   OverrideColorPalette="WhiteBackground")
    Delete(view); Delete(r)
    print("saved " + png, flush=True)

print()
hdr = ("mesh", "cells", "verts", "cell type", "r_wall", "Lx", "Ly", "n_wall",
       "h_first", "h_tan", "AR_wall")
print("%-12s %8s %8s %-12s %7s %7s %7s %7s %11s %11s %8s" % hdr)
for r_ in rows:
    print("%-12s %8d %8d %-12s %7.3f %7.3f %7.3f %7d %11.4e %11.4e %8.2f" % r_)
