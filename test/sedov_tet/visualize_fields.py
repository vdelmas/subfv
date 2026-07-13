from paraview.simple import *
import os

os.makedirs("images", exist_ok=True)

reader = XMLPartitionedUnstructuredGridReader(FileName=["output_-1.pvtu"])

# Slice at z=0 (XY plane)
sliceFilter = Slice(Input=reader)
sliceFilter.SliceType = "Plane"
sliceFilter.SliceType.Origin = [0.0, 0.0, 0.0]
sliceFilter.SliceType.Normal = [0.0, 0.0, 1.0]

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
    "images/sedov_tet_density.png",
    view,
    ImageResolution=[2560, 2560],
    OverrideColorPalette="WhiteBackground"
)

print("Saved: images/sedov_tet_density.png")
