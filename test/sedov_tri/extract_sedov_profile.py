"""
pvbatch extract_sedov_profile.py
Run from outputs/SCHEME/ — reads output_-1.pvtu (final state).
Writes density_profile.dat and velocity_profile.dat.
"""

import os, sys, math

try:
    import paraview.simple as pv
except ImportError:
    print("ERROR: must be run with pvbatch")
    sys.exit(1)

from paraview import servermanager

pv._DisableFirstRenderCameraReset()

pvtu_file = "output_-1.pvtu"

if not os.path.isfile(pvtu_file):
    print(f"SKIP (missing): {pvtu_file}")
    sys.exit(0)

reader = pv.XMLPartitionedUnstructuredGridReader(FileName=[pvtu_file])
pv.UpdatePipeline()
data = servermanager.Fetch(reader)

rxy_arr     = data.GetCellData().GetArray("rxy")
density_arr = data.GetCellData().GetArray("Density")
vel_arr     = data.GetCellData().GetArray("Velocity")

if rxy_arr is None or density_arr is None:
    print(f"ERROR: 'rxy' or 'Density' not found in {pvtu_file}")
    sys.exit(1)

n = data.GetNumberOfCells()
print(f"  {pvtu_file}: {n} cells")

with open("density_profile.dat", "w") as f:
    f.write("# rxy  Density\n")
    for i in range(n):
        f.write(f"{rxy_arr.GetValue(i):.10e}  {density_arr.GetValue(i):.10e}\n")

if vel_arr is not None:
    with open("velocity_profile.dat", "w") as f:
        f.write("# rxy  v_radial\n")
        for i in range(n):
            rxy = rxy_arr.GetValue(i)
            vx  = vel_arr.GetValue(3*i)
            vy  = vel_arr.GetValue(3*i + 1)
            # project velocity onto radial direction
            vr = (vx*vx + vy*vy)**0.5 if rxy < 1e-14 else (vx * data.GetCell(i).GetBounds()[0] + vy * data.GetCell(i).GetBounds()[2]) / rxy
            # simpler: radial speed = sqrt(vx^2+vy^2) signed by dot with r_hat
            # use cell centre from bounds midpoint
            b = data.GetCell(i).GetBounds()
            cx = 0.5*(b[0]+b[1]); cy = 0.5*(b[2]+b[3])
            vr = (vx*cx + vy*cy) / rxy if rxy > 1e-14 else 0.0
            f.write(f"{rxy:.10e}  {vr:.10e}\n")
    print("  -> density_profile.dat, velocity_profile.dat")
else:
    print("  -> density_profile.dat (no Velocity field found)")

pv.Delete(reader)
print("Done.")
