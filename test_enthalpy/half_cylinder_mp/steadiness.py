#!/usr/bin/env pvpython
"""How close to a steady state is each run?

The enthalpy-preservation property is a STEADY-state property: h == h_inf only holds for a
steady solution. This measures ||rho(t2) - rho(t1)|| / ||rho(t1)|| between the last two
solution dumps, i.e. whether the run has actually settled.

usage: pvpython steadiness.py <outputs_dir> [more dirs ...]
"""
import sys, os, re
import numpy as np
from paraview.simple import *
from paraview import servermanager as sm
from paraview.numpy_support import vtk_to_numpy


def rho(pvtu):
    r = XMLPartitionedUnstructuredGridReader(FileName=[pvtu]); r.UpdatePipeline()
    blocks = []

    def walk(o):
        if o.IsA('vtkCompositeDataSet'):
            it = o.NewIterator(); it.InitTraversal()
            while not it.IsDoneWithTraversal():
                walk(it.GetCurrentDataObject()); it.GoToNextItem()
        else:
            blocks.append(o)
    walk(sm.Fetch(r))
    # euler_ho calls it 'rho', ns calls it 'Density'.
    dens = 'rho' if blocks[0].GetCellData().GetArray('rho') is not None else 'Density'
    out = {'rho': np.concatenate([vtk_to_numpy(b.GetCellData().GetArray(dens)) for b in blocks]),
           'H': np.concatenate([vtk_to_numpy(b.GetCellData().GetArray('H')) for b in blocks])}
    Delete(r)
    return out


def last_two(d):
    """The final dump (output_-1) and the numbered dump just before it."""
    n = sorted(int(m.group(1))
               for m in (re.match(r'output_(\d+)\.pvtu$', f) for f in os.listdir(d)) if m)
    if not n:
        raise SystemExit("no numbered dump in " + d)
    return os.path.join(d, 'output_%d.pvtu' % n[-1]), os.path.join(d, 'output_-1.pvtu')


print("%-34s %12s %12s" % ("run", "d(rho)_rel", "d(H)_rel"))
for d in sys.argv[1:]:
    f1, f2 = last_two(d)
    a, b = rho(f1), rho(f2)
    dr = np.linalg.norm(b['rho'] - a['rho']) / np.linalg.norm(a['rho'])
    dh = np.linalg.norm(b['H'] - a['H']) / np.linalg.norm(a['H'])
    print("%-34s %12.3e %12.3e" % (os.path.basename(d.rstrip('/')), dr, dh))
