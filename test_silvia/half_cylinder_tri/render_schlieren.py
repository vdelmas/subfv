# Numerical schlieren of the final state of each order: log(1+|grad rho|/rho_inf),
# inverted grayscale, one colour range shared by all orders (p99.5 of the
# pooled values) so the images are directly comparable.
# Usage: pvbatch render_schlieren.py 1 2 3   (from test_silvia/half_cylinder_tri)
#        order:idx picks outputs/o<order>/output_<idx>.pvtu (default idx -1 = final)
import sys
import numpy as np
from paraview.simple import *
from paraview import servermanager as sm
from vtk.util.numpy_support import vtk_to_numpy

orders = []
for arg in sys.argv[1:] or ["1", "2", "3"]:
    o, _, idx = arg.partition(":")
    orders.append((o, idx or "-1"))
W, H = 1300, 670


def schlieren_source(order, idx):
    reader = XMLPartitionedUnstructuredGridReader(
        FileName=[f"outputs/o{order}/output_{idx}.pvtu"])
    c2p = CellDatatoPointData(Input=reader)
    grad = Gradient(Input=c2p)
    grad.ScalarArray = ["POINTS", "rho"]
    grad.ResultArrayName = "grad_rho"
    calc = Calculator(Input=grad)
    calc.ResultArrayName = "schlieren"
    calc.Function = "ln(1+mag(grad_rho)/1e-3)"  # grad normalised by rho_inf
    # Domain outline: boundary edges of a mid-thickness slice
    cut = Slice(Input=reader)
    cut.SliceType = "Plane"
    cut.SliceType.Origin = [0.0, 0.0, 5.0e-4]
    cut.SliceType.Normal = [0.0, 0.0, 1.0]
    edges = FeatureEdges(Input=cut)
    edges.BoundaryEdges = 1
    edges.FeatureEdges = 0
    edges.NonManifoldEdges = 0
    edges.ManifoldEdges = 0
    return calc, edges


sources = {o: schlieren_source(o, idx) for o, idx in orders}

pooled = []
for o, (src, _) in sources.items():
    data = sm.Fetch(src)
    pooled.append(vtk_to_numpy(data.GetPointData().GetArray("schlieren")))
pooled = np.concatenate(pooled)
vmin, vmax = float(pooled.min()), float(np.percentile(pooled, 99.5))
print(f"schlieren range: [{vmin:.4g}, {vmax:.4g}]")

view = CreateView("RenderView")
view.ViewSize = [W, H]
view.OrientationAxesVisibility = 0
view.UseColorPaletteForBackground = 0
view.Background = [1, 1, 1]
view.Background2 = [1, 1, 1]
view.CameraParallelProjection = 1

lut = GetColorTransferFunction("schlieren")
lut.ApplyPreset("Grayscale", True)
lut.InvertTransferFunction()
lut.RescaleTransferFunction(vmin, vmax)
lut.AutomaticRescaleRangeMode = "Never"

for o, (src, edges) in sources.items():
    disp = Show(src, view)
    ColorBy(disp, ("POINTS", "schlieren"))
    disp.LookupTable = lut
    disp.SetScalarBarVisibility(view, False)
    disp.DisableLighting = 1
    edisp = Show(edges, view)
    edisp.ColorArrayName = [None, ""]
    edisp.AmbientColor = edisp.DiffuseColor = [0.5, 0.5, 0.5]
    edisp.LineWidth = 2.0
    Render(view)
    # Camera AFTER the first Render (Render resets it). Half annulus x<=0,
    # |y|<=2.5: put the free stream (coming from -x) at the top.
    view.CameraFocalPoint = [-1.25, 0.0, 0.0]
    view.CameraPosition = [-1.25, 0.0, 10.0]
    view.CameraViewUp = [-1.0, 0.0, 0.0]
    view.CameraParallelScale = 1.30
    SaveScreenshot(f"outputs/o{o}/schlieren.png", view, ImageResolution=[W, H])
    print("saved", f"outputs/o{o}/schlieren.png")
    Hide(src, view)
    Hide(edges, view)
