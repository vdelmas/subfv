from paraview.simple import *
import os

os.makedirs("images", exist_ok=True)

reader = XMLPartitionedUnstructuredGridReader(FileName=["output_-1.pvtu"])
UpdatePipeline()

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [2560, 1280]
view.CameraParallelProjection = 1
view.Background = [1, 1, 1]
view.OrientationAxesVisibility = 0

display = Show(reader, view)
display.Representation = "Surface"
display.DisableLighting = 1

ColorBy(display, ("CELLS", "Density"))
lut = GetColorTransferFunction("Density")
lut.ApplyPreset("Cool to Warm (Extended)", True)
lut.NumberOfTableValues = 1024
display.RescaleTransferFunctionToDataRange(True)
display.SetScalarBarVisibility(view, False)

# Top-down view, domain [0,2]x[-0.5,0.5]
view.CameraPosition = [1.0, 0, 1]
view.CameraFocalPoint = [1.0, 0, 0]
view.CameraViewUp = [0, 1, 0]
ResetCamera()

Render()
view.CameraParallelScale = 0.5

SaveScreenshot(
    "images/kelvin_helmholtz_density.png",
    view,
    ImageResolution=[2560, 1280],
    OverrideColorPalette="WhiteBackground"
)
print("Saved: images/kelvin_helmholtz_density.png")
