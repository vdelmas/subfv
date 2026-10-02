# Render one half-cylinder .pvtu to a schlieren PNG and a total-enthalpy PNG.
#
# Both color ranges are supplied by the caller so that every panel of the
# three_wave / three_wave_enthalpy x quad / tri comparison is directly comparable.
#
# usage: pvbatch render.py <file.pvtu> <out_prefix> <schlieren_vmax> <H_min> <H_max>
#                          [<rho_min> <rho_max> <p_min> <p_max>]
#
# Density and pressure are optional: pass their ranges to get the extra panels.
#
# Layout follows the project's standing figure style: field full-bleed on the left with no
# axis box, native vertical ParaView scalar bar in a reserved white margin on the right,
# black text, rendered fully offscreen (DISPLAY must be unset).
from paraview.simple import *
import sys, os

src, out_prefix = sys.argv[1], sys.argv[2]
SCHL_VMAX = float(sys.argv[3])
H_MIN, H_MAX = float(sys.argv[4]), float(sys.argv[5])
EXTRA = [float(v) for v in sys.argv[6:10]] if len(sys.argv) >= 10 else None
# "panels" mode: render each field on its own, with NO reserved colour-bar margin, so a LaTeX
# table can place one original image per cell without anything being cropped afterwards. The
# colour bars are then rendered separately (see BARS below), also as original renders: the
# scalar bar is a 2D overlay, so pushing the camera off the domain leaves the bar alone on
# white without touching the pixels of the bar itself.
PANELS = "panels" in sys.argv
BARS = "bars" in sys.argv

# Domain: left half-annulus, r in [1, 2.5] -> x in [-2.5, 0], y in [-2.5, 2.5].
X_LO, X_HI, Y_HALF = -2.5, 0.0, 2.5
H_PX = 1000
FIELD_PX = int(round((X_HI - X_LO) / (2 * Y_HALF) * H_PX))   # 500
BAR_PX = 0 if PANELS else 460
W_PX = FIELD_PX + BAR_PX if not BARS else 420
# Push the camera centre right by half the reserved margin's world width, so the domain
# renders flush in the left FIELD_PX pixels and the margin stays blank for the colour bar.
CAM_X = 0.5 * (X_LO + X_HI) + 0.5 * (BAR_PX / H_PX) * (2 * Y_HALF)

os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)

reader = XMLPartitionedUnstructuredGridReader(FileName=[src])
reader.UpdatePipeline()

# The two solvers name their density gradient differently, so pick whichever this file has:
#   euler_ho (subfveulerho): 'grad_rho_cell', cell data. The vertex twin grad_rho_vert is
#     zeroed on every boundary vertex by its write_vtu, which would draw a spurious white
#     ring along the cylinder wall -- exactly where the interesting physics is -- so the
#     cell array is the one to use.
#   ns (subfvns): 'Nodal_Grad_Density', point data.
_cd = reader.CellData.keys()
if 'grad_rho_cell' in _cd:
    GRAD, SCHL_ASSOC, CALC_ATTR = 'grad_rho_cell', 'CELLS', 'Cell Data'
    RHO, PRES = 'rho', 'p'
else:
    GRAD, SCHL_ASSOC, CALC_ATTR = 'Nodal_Grad_Density', 'POINTS', 'Point Data'
    RHO, PRES = 'Density', 'Pressure'

calc = Calculator(Input=reader)
calc.AttributeType = CALC_ATTR
calc.ResultArrayName = "schlieren"
calc.Function = "log(mag(%s)+1)" % GRAD
calc.UpdatePipeline()


def render(source, assoc, field, preset, vmin, vmax, bar_title, out_png, invert=False,
           labels=None):
    view = CreateView("RenderView")
    view.ViewSize = [W_PX, H_PX]
    view.CameraParallelProjection = 1
    view.OrientationAxesVisibility = 0
    view.UseColorPaletteForBackground = 0
    view.Background = [1, 1, 1]
    view.Background2 = [1, 1, 1]

    d = Show(source, view)
    d.Representation = "Surface"
    d.DisableLighting = 1
    ColorBy(d, (assoc, field))
    lut = GetColorTransferFunction(field)
    lut.ApplyPreset(preset, True)
    if invert:
        lut.InvertTransferFunction()
    lut.RescaleTransferFunction(vmin, vmax)
    lut.NumberOfTableValues = 1024

    d.SetScalarBarVisibility(view, not PANELS or BARS)
    bar = GetScalarBar(lut, view)
    bar.Title = bar_title
    bar.ComponentTitle = ""
    bar.Orientation = "Vertical"
    bar.WindowLocation = "Any Location"
    bar.Position = [0.10 if BARS else FIELD_PX / W_PX + 0.02, 0.08]
    bar.ScalarBarLength = 0.85
    bar.TitleColor = [0, 0, 0]
    bar.LabelColor = [0, 0, 0]
    bar.TitleFontSize = 26
    bar.LabelFontSize = 20
    bar.AutomaticLabelFormat = 1
    if labels is not None:
        # ParaView 6.1 picks the tick count from the bar length alone, which gives ~19 labels
        # on a bar this tall; pin them instead so the two panels carry the same ticks.
        bar.UseCustomLabels = 1
        bar.CustomLabels = labels

    # Render() re-fits the camera, so the camera must be set AFTER it and immediately
    # before SaveScreenshot (which renders again without re-triggering the reset).
    Render()
    # In BARS mode the camera is pushed far off the domain so nothing but the 2D scalar-bar
    # overlay remains in frame.
    cam_x = CAM_X + 1000.0 if BARS else CAM_X
    view.CameraFocalPoint = [cam_x, 0.0, 0.0]
    view.CameraPosition = [cam_x, 0.0, 10.0]
    view.CameraViewUp = [0, 1, 0]
    view.CameraParallelScale = Y_HALF
    SaveScreenshot(out_png, view, ImageResolution=[W_PX, H_PX],
                   OverrideColorPalette="WhiteBackground")
    Delete(view)
    print("saved " + out_png, flush=True)


render(calc, SCHL_ASSOC, "schlieren", "Grayscale", 0.0, SCHL_VMAX,
       "schlieren", out_prefix + "_schlieren.png", invert=True,
       labels=[round(SCHL_VMAX * f, 2) for f in (0.0, 0.25, 0.5, 0.75, 1.0)])
render(reader, "CELLS", "H", "Viridis", H_MIN, H_MAX,
       "H", out_prefix + "_enthalpy.png",
       labels=[round(H_MIN + (H_MAX - H_MIN) * f, 2) for f in (0.0, 0.25, 0.5, 0.75, 1.0)])

if EXTRA is not None:
    r_lo, r_hi, p_lo, p_hi = EXTRA
    render(reader, "CELLS", RHO, "Viridis", r_lo, r_hi,
           "rho", out_prefix + "_density.png",
           labels=[round(r_lo + (r_hi - r_lo) * f, 2) for f in (0.0, 0.25, 0.5, 0.75, 1.0)])
    render(reader, "CELLS", PRES, "Viridis", p_lo, p_hi,
           "p", out_prefix + "_pressure.png",
           labels=[round(p_lo + (p_hi - p_lo) * f, 2) for f in (0.0, 0.25, 0.5, 0.75, 1.0)])
