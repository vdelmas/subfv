from paraview.simple import *
import os

os.makedirs("images", exist_ok=True)

# Domain [0,7] x [0,3]. VIEW_SIZE matches that 7:3 aspect and PARALLEL_SCALE
# is the domain half-height, so the domain fills the whole frame ("plein ecran").
VIEW_SIZE = [1680, 720]
PARALLEL_SCALE = 1.5
GAMMA = 1.4

reader = XMLPartitionedUnstructuredGridReader(FileName=["output_-1.pvtu"])
UpdatePipeline()

d3 = D3(Input=reader)
UpdatePipeline()

# Internal energy per unit mass: e = p / ((gamma - 1) * rho)
calc = Calculator(Input=d3)
calc.AttributeType = "Cell Data"
calc.ResultArrayName = "InternalEnergy"
calc.Function = "Pressure/(%g*Density)" % (GAMMA - 1.0)
UpdatePipeline()

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = VIEW_SIZE
view.CameraParallelProjection = 1
view.Background = [1, 1, 1]
view.OrientationAxesVisibility = 0

display = Show(calc, view)
display.Representation = "Surface"
display.DisableLighting = 1


def save(field, fname):
    ColorBy(display, ("CELLS", field))
    lut = GetColorTransferFunction(field)
    lut.ApplyPreset("Viridis", True)
    lut.NumberOfTableValues = 1024
    display.RescaleTransferFunctionToDataRange(True)
    display.SetScalarBarVisibility(view, False)
    # Frame the domain: explicit camera (don't rely on ResetCamera), and set
    # CameraParallelScale AFTER Render() because Render() re-fits the camera.
    view.CameraFocalPoint = [3.5, 1.5, 0.005]
    view.CameraPosition = [3.5, 1.5, 10]
    view.CameraViewUp = [0, 1, 0]
    Render()
    view.CameraParallelScale = PARALLEL_SCALE
    SaveScreenshot("images/" + fname, view,
                   ImageResolution=VIEW_SIZE,
                   OverrideColorPalette="WhiteBackground")
    print("Saved: images/" + fname)


save("Density", "triple_point_density.png")
save("Pressure", "triple_point_pressure.png")
save("InternalEnergy", "triple_point_internal_energy.png")
