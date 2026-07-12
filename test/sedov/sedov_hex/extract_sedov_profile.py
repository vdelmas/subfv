"""
pvbatch extract_sedov_profile.py
Run from outputs/SCHEME/ — reads output_-1.pvtu (final state).
Writes density_profile.dat.
Columns: r  Density
"""

import os, sys

try:
    import paraview.simple as pv
except ImportError:
    print("ERROR: must be run with pvbatch")
    sys.exit(1)

from paraview import servermanager

pv._DisableFirstRenderCameraReset()

pvtu_file = "output_-1.pvtu"
dat_file  = "density_profile.dat"

if not os.path.isfile(pvtu_file):
    print(f"SKIP (missing): {pvtu_file}")
    sys.exit(0)

reader = pv.XMLPartitionedUnstructuredGridReader(FileName=[pvtu_file])
pv.UpdatePipeline()
data = servermanager.Fetch(reader)

r_arr       = data.GetCellData().GetArray("r")
density_arr = data.GetCellData().GetArray("Density")

if r_arr is None or density_arr is None:
    print(f"ERROR: 'r' or 'Density' not found in {pvtu_file}")
    sys.exit(1)

n = data.GetNumberOfCells()
print(f"  {pvtu_file}: {n} cells -> {dat_file}")
with open(dat_file, "w") as f:
    f.write("# r  Density\n")
    for i in range(n):
        r   = r_arr.GetValue(i)
        rho = density_arr.GetValue(i)
        f.write(f"{r:.10e}  {rho:.10e}\n")

pv.Delete(reader)
print("Done.")
