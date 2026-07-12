from paraview.simple import *
import os
import numpy as np

os.makedirs("images", exist_ok=True)

reader = XMLPartitionedUnstructuredGridReader(FileName=["output_-1.pvtu"])
UpdatePipeline()

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [1920, 640]
view.CameraParallelProjection = 1
view.Background = [1, 1, 1]
view.OrientationAxesVisibility = 0

d3 = D3(Input=reader)
UpdatePipeline()

cell2point = CellDatatoPointData(Input=d3)
UpdatePipeline()

# Density range via Fetch (serial-aggregated)
data = servermanager.Fetch(cell2point)
rho_min, rho_max = data.GetPointData().GetArray("Density").GetRange()

# Mach range via DataInformation (parallel-safe, no Fetch needed)
pt_info = cell2point.GetDataInformation().GetPointDataInformation()
mach_min, mach_max = pt_info.GetArrayInformation("Mach").GetComponentRange(0)

bounds = cell2point.GetDataInformation().GetBounds()
z_mid = 0.5 * (bounds[4] + bounds[5])


def make_iso_display(source, field, n_iso, z_offset):
    """Build a contour → slice → translate pipeline and return its display."""
    values = list(np.linspace(
        *source.GetDataInformation().GetPointDataInformation()
              .GetArrayInformation(field).GetComponentRange(0),
        n_iso + 2
    )[1:-1])
    contour = Contour(Input=source)
    contour.ContourBy = ["POINTS", field]
    contour.Isosurfaces = values
    UpdatePipeline()

    slicer = Slice(Input=contour)
    slicer.SliceType = "Plane"
    slicer.SliceType.Origin = [1.5, 0.5, z_mid]
    slicer.SliceType.Normal = [0, 0, 1]
    UpdatePipeline()

    xfm = Transform(Input=slicer)
    xfm.Transform.Translate = [0, 0, z_offset]
    UpdatePipeline()

    disp = Show(xfm, view)
    disp.Representation = "Surface"
    disp.AmbientColor = [0, 0, 0]
    disp.DiffuseColor = [0, 0, 0]
    disp.LineWidth = 1.0
    disp.ColorArrayName = ['POINTS', '']
    return xfm, disp


xfm_rho, iso_rho = make_iso_display(cell2point, "Density", 50, 0.1)
xfm_mach, iso_mach = make_iso_display(cell2point, "Mach", 30, 0.1)


def apply_diverging_centered(lut, vcenter, preset="Cool to Warm (Extended)"):
    """Resample preset so that vcenter maps to the neutral midpoint colour.

    Must be called AFTER RescaleTransferFunctionToDataRange so the LUT already
    spans [vmin, vmax] in data space.  RGBPoints after ApplyPreset are in data
    space, not [0,1], so we normalise before applying the two-slope mapping.
    """
    lut.ApplyPreset(preset, True)
    pts = list(lut.RGBPoints)
    pt_min = pts[0]
    pt_max = pts[-4]
    pt_range = pt_max - pt_min if pt_max != pt_min else 1.0
    new_pts = []
    for i in range(0, len(pts), 4):
        t, r, g, b = pts[i], pts[i+1], pts[i+2], pts[i+3]
        t_norm = (t - pt_min) / pt_range  # normalise to [0, 1]
        if t_norm <= 0.5:
            val = pt_min + t_norm / 0.5 * (vcenter - pt_min)
        else:
            val = vcenter + (t_norm - 0.5) / 0.5 * (pt_max - vcenter)
        new_pts += [val, r, g, b]
    lut.RGBPoints = new_pts
    lut.NumberOfTableValues = 1024


display = Show(reader, view)
display.Representation = "Surface"
display.DisableLighting = 1

# Top-down view, domain [0,3]x[0,1]
view.CameraPosition = [1.5, 0.5, 1]
view.CameraFocalPoint = [1.5, 0.5, 0]
view.CameraViewUp = [0, 1, 0]
ResetCamera()
view.CameraParallelScale = 0.5

# ================================================================
# IMAGE 1 — Density (Fast colormap) + density isolines
# ================================================================
Hide(xfm_mach, view)

ColorBy(display, ("CELLS", "Density"))
lut_density = GetColorTransferFunction("Density")
lut_density.ApplyPreset("Fast", True)
lut_density.NumberOfTableValues = 1024
display.RescaleTransferFunctionToDataRange(True)
display.SetScalarBarVisibility(view, False)

Render()

SaveScreenshot(
    "images/forward_facing_density.png",
    view,
    ImageResolution=[1920, 640],
    OverrideColorPalette="WhiteBackground"
)
print("Saved: images/forward_facing_density.png")


# ================================================================
# IMAGE 2 — Mach (Cool to Warm Extended, centred on Mach=1) + Mach isolines
# ================================================================
Hide(xfm_rho, view)
Show(xfm_mach, view)

ColorBy(display, ("CELLS", "Mach"))
lut_mach = GetColorTransferFunction("Mach")
# Let ParaView aggregate the true range across all MPI ranks, then
# apply the two-slope diverging map with Mach=1 at the neutral centre.
display.RescaleTransferFunctionToDataRange()
apply_diverging_centered(lut_mach, vcenter=1.0)

display.SetScalarBarVisibility(view, False)
Render()

SaveScreenshot(
    "images/forward_facing_mach.png",
    view,
    ImageResolution=[1920, 640],
    OverrideColorPalette="WhiteBackground"
)
print("Saved: images/forward_facing_mach.png")
