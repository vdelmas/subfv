from paraview.simple import *
import os

os.makedirs("images", exist_ok=True)

reader = XMLPartitionedUnstructuredGridReader(FileName=["output_-1.pvtu"])

# Crinkle clip at z=0 (keeps full hex cells intersecting the plane)
sliceFilter = Clip(Input=reader)
sliceFilter.ClipType = "Plane"
sliceFilter.ClipType.Origin = [0.0, 0.0, 0.0]
sliceFilter.ClipType.Normal = [0.0, 0.0, 1.0]
sliceFilter.Crinkleclip = 1

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [2560, 2560]
view.CameraParallelProjection = 1
view.Background = [1, 1, 1]
view.OrientationAxesVisibility = 0

display = Show(sliceFilter, view)
display.Representation = "Surface With Edges"
display.DisableLighting = 1

ColorBy(display, ("CELLS", "Density"))

lut = GetColorTransferFunction("Density")
lut.ApplyPreset("Fast", True)
display.RescaleTransferFunctionToDataRange()
display.SetScalarBarVisibility(view, False)

# Top-down view (XY plane)
view.CameraPosition = [0, 0, 1]
view.CameraFocalPoint = [0, 0, 0]
view.CameraViewUp = [0, 1, 0]
ResetCamera()

Render()
view.CameraParallelScale = 1.2

SaveScreenshot(
    "images/sedov_hex_density.png",
    view,
    ImageResolution=[2560, 2560],
    OverrideColorPalette="WhiteBackground"
)

print("Saved: images/sedov_hex_density.png")
