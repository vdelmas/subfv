from paraview.simple import *
import sys, time

# usage: pvbatch render_cyl.py <rundir> <out_dir> <frame_start> <frame_end> <scale>
rundir, out_dir, frame_start, frame_end, scale = sys.argv[1:6]
frame_start, frame_end = int(frame_start), int(frame_end)
scale = float(scale)

Lx, Ly = 4.0, 2.0
W = int(round(1600 * scale))
H = int(round(800 * scale))

import os
os.makedirs(f"{out_dir}/rho", exist_ok=True)
os.makedirs(f"{out_dir}/schlieren", exist_ok=True)

reader = XMLPartitionedUnstructuredGridReader(FileName=[f"{rundir}/output_{frame_start}.pvtu"])
reader.UpdatePipeline()

calc = Calculator(Input=reader)
calc.AttributeType = 'Cell Data'
calc.ResultArrayName = 'schlieren'
calc.Function = "log(mag(grad_rho_cell)*0.3+1)"

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [W, H]
view.CameraParallelProjection = 1
view.CameraPosition = [Lx / 2.0, Ly / 2.0, 1]
view.CameraFocalPoint = [Lx / 2.0, Ly / 2.0, 0]
view.CameraViewUp = [0, 1, 0]
view.OrientationAxesVisibility = 0
view.Background = [1, 1, 1]

display = Show(calc, view)
display.Representation = "Surface"
display.DisableLighting = 1

t0 = time.time()
for i in range(frame_start, frame_end + 1):
    reader.FileName = [f"{rundir}/output_{i}.pvtu"]
    reader.UpdatePipeline()

    ColorBy(display, ("CELLS", "rho"))
    lut = GetColorTransferFunction("rho")
    lut.ApplyPreset("Viridis", True)
    display.RescaleTransferFunctionToDataRange(False)
    Render()
    view.CameraPosition = [Lx / 2.0, Ly / 2.0, 1]
    view.CameraFocalPoint = [Lx / 2.0, Ly / 2.0, 0]
    view.CameraViewUp = [0, 1, 0]
    view.CameraParallelScale = Ly / 2.0 * 0.995
    SaveScreenshot(f"{out_dir}/rho/frame_{i:03d}.png", view, ImageResolution=[W, H])

    ColorBy(display, ("CELLS", "schlieren"))
    lut2 = GetColorTransferFunction("schlieren")
    lut2.ApplyPreset("Grayscale", True)
    lut2.InvertTransferFunction()
    display.RescaleTransferFunctionToDataRange(False)
    Render()
    view.CameraPosition = [Lx / 2.0, Ly / 2.0, 1]
    view.CameraFocalPoint = [Lx / 2.0, Ly / 2.0, 0]
    view.CameraViewUp = [0, 1, 0]
    view.CameraParallelScale = Ly / 2.0 * 0.995
    SaveScreenshot(f"{out_dir}/schlieren/frame_{i:03d}.png", view, ImageResolution=[W, H])

    elapsed = time.time() - t0
    print(f"frame {i} done, elapsed={elapsed:.1f}s", flush=True)

print("ALL_FRAMES_DONE", flush=True)
