from paraview.simple import *
import os

os.makedirs("images", exist_ok=True)

reader = XMLPartitionedUnstructuredGridReader(FileName=["output_-1.pvtu"])
UpdatePipeline()

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [1280, 1280]
view.CameraParallelProjection = 1
view.Background = [1, 1, 1]
view.OrientationAxesVisibility = 0

display = Show(reader, view)
display.Representation = "Surface With Edges"
display.DisableLighting = 1

ColorBy(display, ("CELLS", "Velocity"))
lut = GetColorTransferFunction("Velocity")
lut.VectorMode = "Magnitude"
display.RescaleTransferFunctionToDataRange(True)
lut.ApplyPreset("Fast", True)
display.SetScalarBarVisibility(view, False)

# Top-down view, domain [0,1]x[0,1]
view.CameraPosition = [0.5, 0.5, 1]
view.CameraFocalPoint = [0.5, 0.5, 0]
view.CameraViewUp = [0, 1, 0]
ResetCamera()

Render()
view.CameraParallelScale = 0.5

SaveScreenshot(
    "images/gresho_tri_velocity.png",
    view,
    ImageResolution=[1280, 1280],
    OverrideColorPalette="WhiteBackground"
)
print("Saved: images/gresho_tri_velocity.png")
