"""
pvbatch extract_velocity_profile.py
Run from outputs/SCHEME/ — reads output_0.pvtu (initial) and output_-1.pvtu (final).
Writes velocity_profile_initial.dat and velocity_profile_final.dat.
Columns: Rxy  |V|   (radius from vortex center (0.5, 0.5), velocity norm)
"""

import os, sys

try:
    import paraview.simple as pv
except ImportError:
    print("ERROR: must be run with pvbatch")
    sys.exit(1)

from paraview import servermanager

pv._DisableFirstRenderCameraReset()

CX, CY = 0.5, 0.5

def extract(pvtu_file, dat_file):
    if not os.path.isfile(pvtu_file):
        print(f"SKIP (missing): {pvtu_file}")
        return

    reader = pv.XMLPartitionedUnstructuredGridReader(FileName=[pvtu_file])
    pv.UpdatePipeline()
    data = servermanager.Fetch(reader)

    centroid_arr = data.GetCellData().GetArray("Centroid")
    velocity_arr = data.GetCellData().GetArray("Velocity")

    if centroid_arr is None or velocity_arr is None:
        print(f"ERROR: Centroid or Velocity not found in {pvtu_file}")
        return

    n = data.GetNumberOfCells()
    print(f"  {pvtu_file}: {n} cells -> {dat_file}")
    with open(dat_file, "w") as f:
        f.write("# Rxy  |V|\n")
        for i in range(n):
            cx, cy, _ = centroid_arr.GetTuple(i)
            vx, vy, _ = velocity_arr.GetTuple(i)
            Rxy = ((cx - CX)**2 + (cy - CY)**2)**0.5
            V   = (vx**2 + vy**2)**0.5
            f.write(f"{Rxy:.10e}  {V:.10e}\n")

    pv.Delete(reader)

extract("output_0.pvtu",  "velocity_profile_initial.dat")
extract("output_-1.pvtu", "velocity_profile_final.dat")
print("Done.")
