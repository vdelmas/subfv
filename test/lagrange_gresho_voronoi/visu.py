"""Wireframe of the FINAL Lagrangian mesh for each Gresho run.

Run with:  pvbatch visu.py
Produces:  outputs/<scheme>_<mesh>/wireframe_final.png   (and wireframe_init.png)

The point of the figure: the Gresho vortex is a steady solution, so in a
Lagrangian scheme the cells are dragged around by a velocity field that should
stay put. The shear between r<0.2 (solid-body rotation) and the outer ring is
what tangles the mesh, so the final wireframe shows directly how much spurious
deformation each nodal solver produces.
"""
import os
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

BASE_PX = 1600
ROOT = os.path.dirname(os.path.abspath(__file__))

CASES = [(s, m) for m in ("gresho_quad", "gresho_tri")
                for s in ("classic", "sidil")]


def get_img_size(reader, base=BASE_PX):
    reader.UpdatePipeline()
    b = reader.GetDataInformation().GetBounds()
    dx, dy = b[1] - b[0], b[3] - b[2]
    if dx >= dy:
        return base, max(1, round(base * dy / dx))
    return max(1, round(base * dx / dy)), base


def top_view(view, reader, W, H):
    reader.UpdatePipeline()
    b = reader.GetDataInformation().GetBounds()
    xc, yc = (b[0] + b[1]) / 2.0, (b[2] + b[3]) / 2.0
    dx, dy = b[1] - b[0], b[3] - b[2]
    ps = max(dy / 2.0, dx / (2.0 * (W / float(H)))) * 1.03
    cam = view.GetActiveCamera()
    cam.SetPosition(xc, yc, 1.0)
    cam.SetFocalPoint(xc, yc, 0.0)
    cam.SetViewUp(0.0, 1.0, 0.0)
    cam.SetParallelProjection(1)
    cam.SetParallelScale(ps)


view = GetActiveViewOrCreate('RenderView')
view.UseColorPaletteForBackground = 0
view.Background = [1.0, 1.0, 1.0]
view.Background2 = [1.0, 1.0, 1.0]
view.BackgroundColorMode = 0
view.OrientationAxesVisibility = 0

for scheme, mesh in CASES:
    d = os.path.join(ROOT, "outputs", f"{scheme}_{mesh}")
    for tag, pvtu in (("final", "output_-1.pvtu"), ("init", "output_0.pvtu")):
        src = os.path.join(d, pvtu)
        if not os.path.exists(src):
            print(f"-- absent: {scheme}/{mesh} {pvtu}")
            continue
        r = OpenDataFile(src)
        rep = Show(r, view)
        rep.Representation = 'Wireframe'
        rep.AmbientColor = [0, 0, 0]
        rep.LineWidth = 1.0
        W, H = get_img_size(r)
        top_view(view, r, W, H)
        RenderAllViews()
        out = os.path.join(d, f"wireframe_{tag}.png")
        SaveScreenshot(out, view, ImageResolution=[W, H])
        print(f"-> {out}  [{W}x{H}]")
        Delete(r)
        del r
