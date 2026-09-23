#!/usr/bin/env pvpython
"""Print the shared colour ranges to render every panel of the comparison with.

Panels are only comparable if they share a range, so the range is taken over ALL the runs
passed in: the total-enthalpy min/max across every cell of every run, and the schlieren
99.9th percentile (not the max -- one carbuncle cell would otherwise set the scale and wash
out the whole shock).

usage: pvpython ranges.py <file.pvtu> [more.pvtu ...]   ->  "<schl_vmax> <H_min> <H_max>"
"""
import sys
import numpy as np
from paraview.simple import *
from paraview import servermanager as sm
from paraview.numpy_support import vtk_to_numpy

H, S = [], []
for f in sys.argv[1:]:
    r = XMLPartitionedUnstructuredGridReader(FileName=[f]); r.UpdatePipeline()
    blocks = []

    def walk(o):
        if o.IsA('vtkCompositeDataSet'):
            it = o.NewIterator(); it.InitTraversal()
            while not it.IsDoneWithTraversal():
                walk(it.GetCurrentDataObject()); it.GoToNextItem()
        else:
            blocks.append(o)
    walk(sm.Fetch(r))
    H.append(np.concatenate([vtk_to_numpy(b.GetCellData().GetArray('H')) for b in blocks]))
    # Same two naming conventions as render.py -- see the comment there.
    if blocks[0].GetCellData().GetArray('grad_rho_cell') is not None:
        g = np.concatenate([vtk_to_numpy(b.GetCellData().GetArray('grad_rho_cell'))
                            for b in blocks])
    else:
        g = np.concatenate([vtk_to_numpy(b.GetPointData().GetArray('Nodal_Grad_Density'))
                            for b in blocks])
    S.append(np.log(np.linalg.norm(g, axis=1) + 1))
    Delete(r)

H = np.concatenate(H); S = np.concatenate(S)
print("%.4f %.4f %.4f" % (np.quantile(S, 0.999), H.min(), H.max()))
