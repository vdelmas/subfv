"""L2(rho-1) for each Mach of a Lagrangian Gresho convergence sweep.

Same definition as ns_gresho_convergence_main.F90: sqrt(sum_i V_i (rho_i-1)^2)
over cells at t_final, so the numbers can be put next to the Euler ones.
"""
import os, sys, glob
from paraview.simple import *
from paraview import servermanager as sm

root = sys.argv[1]
for case in sorted(os.listdir(root)):
    d = os.path.join(root, case)
    if not os.path.isdir(d): continue
    print(f"=== {case} ===")
    for ma in ("1e-1", "1e-2", "1e-3", "1e-4"):
        f = os.path.join(d, f"Ma{ma}", "output_-1.pvtu")
        if not os.path.exists(f):
            print(f"  Ma={ma}: (absent)"); continue
        r = OpenDataFile(f); r.UpdatePipeline()
        c = Calculator(Input=r); c.AttributeType = 'Cell Data'
        c.ResultArrayName = 'dr2'; c.Function = '(Density-1)^2'
        iv = IntegrateVariables(Input=c); iv.DivideCellDataByVolume = 0
        val = sm.Fetch(iv).GetCellData().GetArray('dr2').GetValue(0)
        print(f"  Ma={ma}: L2_rho = {val**0.5:.5e}")
        Delete(iv); Delete(c); Delete(r)
