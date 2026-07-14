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

view.CameraPosition = [0.5, 0.5, 1]
view.CameraFocalPoint = [0.5, 0.5, 0]
view.CameraViewUp = [0, 1, 0]
ResetCamera()
view.CameraParallelScale = 0.5

def save_field(field, is_vector, filename, preset="Fast"):
    ColorBy(display, ("CELLS", field))
    lut = GetColorTransferFunction(field)
    if is_vector:
        lut.VectorMode = "Magnitude"
    display.RescaleTransferFunctionToDataRange(True)
    lut.ApplyPreset(preset, True)
    display.SetScalarBarVisibility(view, False)
    Render()
    SaveScreenshot(filename, view, ImageResolution=[1280, 1280],
                   OverrideColorPalette="WhiteBackground")
    print(f"Saved: {filename}")

save_field("Velocity", True,  "images/gresho_tri_velocity.png")
save_field("Density",  False, "images/gresho_tri_density.png")
save_field("Pressure", False, "images/gresho_tri_pressure.png")
