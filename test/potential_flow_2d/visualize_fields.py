from paraview.simple import *
import os

os.makedirs("images", exist_ok=True)

W, H = 1920, 1920  # square: annular domain is symmetric

reader = XMLPartitionedUnstructuredGridReader(FileName=["output_-1.pvtu"])
UpdatePipeline()

bounds   = reader.GetDataInformation().GetBounds()
x_center = 0.5 * (bounds[0] + bounds[1])
y_center = 0.5 * (bounds[2] + bounds[3])
z_mid    = 0.5 * (bounds[4] + bounds[5])
height   = bounds[3] - bounds[2]
width    = bounds[1] - bounds[0]
aspect   = W / H
cam_scale = max(height / 2.0, width / (2.0 * aspect))

# ================================================================
# Pipeline
# ================================================================
d3  = D3(Input=reader)
c2p = CellDatatoPointData(Input=d3)

# Velocity magnitude via Calculator
calc = Calculator(Input=c2p)
calc.ResultArrayName = "Vmag"
calc.Function = "mag(Velocity)"
UpdatePipeline()

# Fetch range of Vmag across the domain
data  = servermanager.Fetch(calc)
arr   = data.GetPointData().GetArray("Vmag")
vmin, vmax = arr.GetRange()
print(f"Vmag range: [{vmin:.6g}, {vmax:.6g}]")

n_iso = 40
iso_vals = [vmin + i * (vmax - vmin) / (n_iso - 1) for i in range(n_iso)]

# Contour lines of velocity magnitude
contour = Contour(Input=calc)
contour.ContourBy   = ["POINTS", "Vmag"]
contour.Isosurfaces = iso_vals

# Slice at z_mid to project isolines onto 2D plane
sl_iso = Slice(Input=contour)
sl_iso.SliceType        = "Plane"
sl_iso.SliceType.Origin = [x_center, y_center, z_mid]
sl_iso.SliceType.Normal = [0, 0, 1]

# Lift isolines slightly above z=0 so they appear on top of the mesh
xf_iso = Transform(Input=sl_iso)
xf_iso.Transform.Translate = [0, 0, 0.1]
UpdatePipeline()

# Mesh slice for wireframe display
sl_mesh = Slice(Input=d3)
sl_mesh.SliceType        = "Plane"
sl_mesh.SliceType.Origin = [x_center, y_center, z_mid]
sl_mesh.SliceType.Normal = [0, 0, 1]
UpdatePipeline()

# ================================================================
# View
# ================================================================
view = GetActiveViewOrCreate("RenderView")
view.ViewSize                  = [W, H]
view.CameraParallelProjection  = 1
view.Background                = [1, 1, 1]
view.OrientationAxesVisibility = 0
view.CameraPosition            = [x_center, y_center, 1]
view.CameraFocalPoint          = [x_center, y_center, 0]
view.CameraViewUp              = [0, 1, 0]
view.CameraParallelScale       = cam_scale

# Wireframe mesh (grey edges on white)
disp_mesh = Show(sl_mesh, view)
disp_mesh.Representation = "Wireframe"
disp_mesh.AmbientColor   = [0.7, 0.7, 0.7]
disp_mesh.DiffuseColor   = [0.7, 0.7, 0.7]
disp_mesh.LineWidth      = 0.5
disp_mesh.ColorArrayName = ["POINTS", ""]

# Velocity iso-contours (black lines on top)
disp_iso = Show(xf_iso, view)
disp_iso.Representation = "Surface"
disp_iso.AmbientColor   = [0, 0, 0]
disp_iso.DiffuseColor   = [0, 0, 0]
disp_iso.LineWidth      = 1.5
disp_iso.ColorArrayName = ["POINTS", ""]

Render()

view.CameraPosition      = [x_center, y_center, 1]
view.CameraFocalPoint    = [x_center, y_center, 0]
view.CameraViewUp        = [0, 1, 0]
view.CameraParallelScale = cam_scale

SaveScreenshot(
    os.path.abspath("images/velocity_iso.png"),
    view,
    ImageResolution=[W, H],
    OverrideColorPalette="WhiteBackground"
)
print("Saved: images/velocity_iso.png")
