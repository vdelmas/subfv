"""
pvbatch visu_corners.py
Génère une image PNG haute résolution du maillage pour chacun des 4 cas aux coins.
Fond blanc, sans axes d'orientation.
Usage: pvbatch visu_corners.py
"""

import os, sys

try:
    from paraview.simple import *
except ImportError:
    print("ERROR: must be run with pvbatch")
    sys.exit(1)

paraview.simple._DisableFirstRenderCameraReset()

ROOT = os.path.dirname(os.path.abspath(__file__))
CORNERS_DIR = os.path.join(ROOT, "corners_visu")

RESOLUTION = [7200, 5400]  # haute résolution 4:3, ~600dpi sur 30x22cm

CORNERS = {
    "top_left":     "FF grossier / CL fine",
    "top_right":    "FF grossier / CL grossiere",
    "bottom_left":  "FF fin / CL fine",
    "bottom_right": "FF fin / CL grossiere",
}

for label, description in CORNERS.items():
    pvtu_file = os.path.join(CORNERS_DIR, label, "output_0.pvtu")
    png_file  = os.path.join(CORNERS_DIR, label + "_mesh.png")

    if not os.path.isfile(pvtu_file):
        print(f"SKIP (missing): {pvtu_file}")
        continue

    print(f"Processing: {label} ({description})")

    # Create view
    view = CreateView("RenderView")
    view.ViewSize = RESOLUTION
    view.OrientationAxesVisibility = 0

    # Load solution
    reader = XMLPartitionedUnstructuredGridReader(FileName=[pvtu_file])
    UpdatePipeline()

    # Display: surface with edges (mesh visualization)
    display = Show(reader, view)
    display.Representation = "Surface With Edges"
    display.EdgeColor = [0.0, 0.0, 0.0]
    display.AmbientColor = [0.9, 0.9, 0.9]
    display.DiffuseColor = [0.9, 0.9, 0.9]
    display.SetScalarColoring(None, 0)

    ResetCamera(view)
    Render(view)

    SaveScreenshot(png_file, view, ImageResolution=RESOLUTION,
                   OverrideColorPalette='WhiteBackground')
    print(f"  -> {png_file}  ({RESOLUTION[0]}x{RESOLUTION[1]})")

    Delete(display)
    Delete(reader)
    Delete(view)

print("Done.")
