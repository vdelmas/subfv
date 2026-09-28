#!/usr/bin/env pvpython
"""Export the wall-coefficient surface written by subfvns to CSV.

Pure data extraction -- no Render(), no SaveScreenshot -- so it runs on an ordinary compute
node without a GL context, unlike the field-rendering part of the project's own
visualize_fields.py that this is carved out of.

Produces coeffs_export.csv with the columns the existing half_cylinder_*_ns plot scripts
expect: 7 = Theta3 (angle along the wall), 8 = Cp, 9 = wall heat flux q.

usage: pvpython export_coeffs.py [coeffs_-1.pvtu] [coeffs_export.csv]
"""
import sys
from paraview.simple import *

src = sys.argv[1] if len(sys.argv) > 1 else "coeffs_-1.pvtu"
out = sys.argv[2] if len(sys.argv) > 2 else "coeffs_export.csv"

reader = XMLPartitionedUnstructuredGridReader(FileName=[src])
SaveData(out, proxy=reader, Precision=10)
print("wrote " + out)
