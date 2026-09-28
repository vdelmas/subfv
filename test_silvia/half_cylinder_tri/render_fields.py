# Cell fields of the final state in viridis, colour bar in a white margin on
# the right. Usage: pvbatch render_fields.py <outdir> [fields...]
#   e.g. pvbatch render_fields.py outputs/o1_fine rho p T
import sys
from paraview.simple import *

outdir = sys.argv[1]
fields = sys.argv[2:] or ["rho", "p", "T"]
H, field_w, bar_w = 670, 1300, 260
W = field_w + bar_w
scale = 1.30                      # half height of the view in world units
ly = 2.0 * scale
shift = (bar_w / H) * ly / 2.0    # move the domain left to free the margin

reader = XMLPartitionedUnstructuredGridReader(FileName=[f"{outdir}/output_-1.pvtu"])

cut = Slice(Input=reader)
cut.SliceType = "Plane"
cut.SliceType.Origin = [0.0, 0.0, 5.0e-4]
cut.SliceType.Normal = [0.0, 0.0, 1.0]
edges = FeatureEdges(Input=cut)
edges.BoundaryEdges = 1
edges.FeatureEdges = 0
edges.NonManifoldEdges = 0
edges.ManifoldEdges = 0

view = CreateView("RenderView")
view.ViewSize = [W, H]
view.OrientationAxesVisibility = 0
view.UseColorPaletteForBackground = 0
view.Background = [1, 1, 1]
view.Background2 = [1, 1, 1]
view.CameraParallelProjection = 1

disp = Show(reader, view)
disp.DisableLighting = 1
edisp = Show(edges, view)
edisp.ColorArrayName = [None, ""]
edisp.AmbientColor = edisp.DiffuseColor = [0.5, 0.5, 0.5]
edisp.LineWidth = 2.0

for field in fields:
    ColorBy(disp, ("CELLS", field))
    lut = GetColorTransferFunction(field)
    lut.ApplyPreset("Viridis", True)
    disp.RescaleTransferFunctionToDataRange(False, True)
    bar = GetScalarBar(lut, view)
    bar.Orientation = "Vertical"
    bar.WindowLocation = "Any Location"
    bar.Position = [field_w / W + 0.02, 0.08]
    bar.ScalarBarLength = 0.85
    bar.Title = field
    bar.ComponentTitle = ""
    bar.TitleColor = [0, 0, 0]
    bar.LabelColor = [0, 0, 0]
    bar.TitleFontSize = 22
    bar.LabelFontSize = 18
    disp.SetScalarBarVisibility(view, True)
    Render(view)
    # Camera AFTER Render (Render resets it); free stream (from -x) on top.
    # Screen right = +y here: shifting the view centre by +y moves the
    # domain to the left of the image.
    view.CameraFocalPoint = [-1.25, shift, 0.0]
    view.CameraPosition = [-1.25, shift, 10.0]
    view.CameraViewUp = [-1.0, 0.0, 0.0]
    view.CameraParallelScale = scale
    out = f"{outdir}/viridis_{field}.png"
    SaveScreenshot(out, view, ImageResolution=[W, H])
    print("saved", out)
    disp.SetScalarBarVisibility(view, False)
