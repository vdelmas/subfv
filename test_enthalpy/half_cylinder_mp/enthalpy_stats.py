#!/usr/bin/env pvpython
"""Total-enthalpy uniformity statistics for the half-cylinder runs.

For a steady inviscid solution fed by a uniform inflow, h = e + p/rho is a Riemann invariant
along streamlines, so the exact solution has h == h_inf in the whole domain, bow shock
included. This script measures how far each run is from that, which is the entire point of
flux_scheme='three_wave_enthalpy' versus the baseline 'three_wave'.

usage: pvpython enthalpy_stats.py <h_inf> <file.pvtu> [more.pvtu ...]

Reported per run: the distribution of |H - h_inf|/h_inf over the cells. The MEDIAN and the
upper quantiles matter more than the max -- on this case the 1D solvers develop a carbuncle
kink at one point of the shock, and that single unsteady spot dominates the max for both
schemes while saying nothing about the enthalpy property itself.
"""
import sys, os
import numpy as np
from paraview.simple import *
from paraview import servermanager as sm
from paraview.numpy_support import vtk_to_numpy


def cell_arrays(pvtu, names):
    r = XMLPartitionedUnstructuredGridReader(FileName=[pvtu])
    r.UpdatePipeline()
    blocks = []

    def walk(o):
        if o.IsA('vtkCompositeDataSet'):
            it = o.NewIterator(); it.InitTraversal()
            while not it.IsDoneWithTraversal():
                walk(it.GetCurrentDataObject()); it.GoToNextItem()
        else:
            blocks.append(o)
    walk(sm.Fetch(r))
    out = {n: np.concatenate([vtk_to_numpy(b.GetCellData().GetArray(n)) for b in blocks])
           for n in names}
    Delete(r)
    return out


def main():
    hinf = float(sys.argv[1])
    print("h_inf = %.10g" % hinf)
    print("%-34s %8s %10s %10s %10s %10s %10s" %
          ("run", "ncells", "median", "q90", "q99", "rms", "max"))
    for pvtu in sys.argv[2:]:
        A = cell_arrays(pvtu, ['H', 'Centroid'])
        H = A['H']
        dev = np.abs(H - hinf) / hinf
        tag = os.path.basename(os.path.dirname(pvtu))
        print("%-34s %8d %10.3e %10.3e %10.3e %10.3e %10.3e" %
              (tag, H.size, np.median(dev), np.quantile(dev, 0.9),
               np.quantile(dev, 0.99), np.sqrt((dev ** 2).mean()), dev.max()))
        c = A['Centroid']
        i = int(np.argmax(dev))
        print("%-34s   worst cell at r=%.3f theta=%.1f deg, H=%.6g" %
              ("", np.hypot(c[i, 0], c[i, 1]), np.degrees(np.arctan2(c[i, 1], c[i, 0])), H[i]))


main()
