from paraview.simple import *
import os

input_files = ["output_-1.pvtu"]

fields = [
    "Pressure",
    "Temperature"
]

output_dir = "images"
os.makedirs(output_dir, exist_ok=True)


def setup_view():
    view = GetActiveViewOrCreate("RenderView")

    view.ViewSize = [1920, 1080]
    view.CameraParallelProjection = 1

    view.CameraPosition = [0, 0, 1]
    view.CameraFocalPoint = [0, 0, 0]
    view.CameraViewUp = [0, 1, 0]
    view.CameraParallelScale = 1.

    view.Background = [1, 1, 1]
    view.OrientationAxesVisibility = 0

    cam = view.GetActiveCamera()
    cam.Azimuth(0)
    cam.Elevation(0)
    cam.Roll(-90)

    return view


def setup_field(display, field, view):
    ColorBy(display, ("CELLS", field))

    display.Representation = "Surface"
    display.DisableLighting = 1

    lut = GetColorTransferFunction(field)
    pwf = GetOpacityTransferFunction(field)

    lut.ApplyPreset("Fast", True)

    display.RescaleTransferFunctionToDataRange()
    lut.RescaleTransferFunctionToDataRange()

    scalarBar = GetScalarBar(lut, view)
    scalarBar.Orientation = "Horizontal"
    scalarBar.WindowLocation = "Lower Center"
    scalarBar.ScalarBarLength = 0.6
    scalarBar.ScalarBarThickness = 25

    scalarBar.Title = field
    scalarBar.ComponentTitle = ""

    scalarBar.DrawTickMarks = 0
    scalarBar.DrawTickLabels = 0

    return lut


for f in input_files:

    print(f"Processing: {f}")

    reader = XMLPartitionedUnstructuredGridReader(FileName=[f])

    view = setup_view()

    display = Show(reader, view)

    Render()

    view.CameraParallelScale = 1.6
    view.CameraPosition[0] = view.CameraPosition[0] - 0.6
    view.CameraFocalPoint[0] = view.CameraFocalPoint[0] - 0.6

    # =========================================================
    # PLOT OVER LINE
    # =========================================================

    d3 = D3(Input=reader)

    c2p = CellDatatoPointData(Input=d3)

    plotLine = PlotOverLine(Input=c2p)

    plotLine.Point1 = [-3.0, 0.0, 0.0]
    plotLine.Point2 = [-1.0, 0.0, 0.0]

    plotLine.SamplingPattern = 'Sample At Cell Boundaries'

    SaveData(
        "line_export.csv",
        proxy=plotLine,
        Precision=10
    )

    display.Representation = "Surface With Edges"

    Render()

    base = os.path.splitext(os.path.basename(f))[0]

    mesh_out = os.path.join(output_dir, f"{base}_mesh.png")

    SaveScreenshot(
        mesh_out,
        view,
        ImageResolution=[1920, 1080],
        OverrideColorPalette="WhiteBackground"
    )

    print("Saved mesh:", mesh_out)

    for field in fields:

        print(f"  -> Field: {field}")

        setup_field(display, field, view)

        Render()

        # =========================================================
        # TEST AJOUTS : ISOLINES
        # =========================================================

        try:

            # conversion cell -> point
            d3_iso = D3(Input=reader)

            c2p_iso = CellDatatoPointData(Input=d3_iso)

            # contour
            contour = Contour(Input=c2p_iso)

            contour.ContourBy = ['POINTS', field]

            # =========================================================
            # MIN / MAX ROBUSTE (PVBATCH SAFE)
            # =========================================================
            
            data = servermanager.Fetch(c2p_iso)
            
            array = data.GetPointData().GetArray(field)
            
            vmin, vmax = array.GetRange()


            # 30 iso valeurs
            n_iso = 30

            contour.Isosurfaces = [
                vmin + i * (vmax - vmin) / (n_iso - 1)
                for i in range(n_iso)
            ]

            # slice selon +Z
            slice1 = Slice(Input=contour)

            slice1.SliceType = 'Plane'

            slice1.SliceType.Origin = [0.0, 0.0, 0.0]

            slice1.SliceType.Normal = [0.0, 0.0, 1.0]

            # translation +Z
            transform = Transform(Input=slice1)

            transform.Transform = 'Transform'

            transform.Transform.Translate = [0.0, 0.0, 1.0]

            # affichage noir
            isoDisplay = Show(transform, view)

            isoDisplay.Representation = 'Surface'

            isoDisplay.ColorArrayName = [None, '']

            isoDisplay.DiffuseColor = [0.0, 0.0, 0.0]

            isoDisplay.LineWidth = 2.0

            Render()

        except Exception as e:

            print("Isoline error:", e)

        # =========================================================
        # FIN TEST AJOUTS
        # =========================================================

        base = os.path.splitext(os.path.basename(f))[0]

        out_file = os.path.join(
            output_dir,
            f"{base}_{field}.png"
        )

        SaveScreenshot(
            out_file,
            view,
            ImageResolution=[1920, 1080],
            OverrideColorPalette="WhiteBackground"
        )

        print("     saved:", out_file)

        display.SetScalarBarVisibility(view, False)

        # =========================================================
        # TEST AJOUTS : CLEANUP ISOLINES
        # =========================================================

        try:

            Delete(isoDisplay)
            Delete(transform)
            Delete(slice1)
            Delete(contour)
            Delete(c2p_iso)
            Delete(d3_iso)

        except:
            pass

        # =========================================================
        # FIN TEST AJOUTS
        # =========================================================

    Delete(reader)

    Delete(display)

    del reader
    del display


coeff_reader = XMLPartitionedUnstructuredGridReader(
    FileName=["coeffs_-1.pvtu"]
)

SaveData(
    "coeffs_export.csv",
    proxy=coeff_reader,
    Precision=10
)
