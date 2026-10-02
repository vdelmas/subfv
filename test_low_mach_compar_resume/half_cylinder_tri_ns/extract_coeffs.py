"""
pvbatch extract_coeffs.py
Run from outputs/WIP/ — reads output_-1.pvtu and coeffs_-1.pvtu.
Writes line_export.csv and coeffs_export.csv. No Show()/Render() calls, so
this works even on compute nodes with no OpenGL/OSMesa (unlike
visualize_fields.py's screenshot step).
"""
from paraview.simple import *

reader = XMLPartitionedUnstructuredGridReader(FileName=["output_-1.pvtu"])
UpdatePipeline()

d3 = D3(Input=reader)
c2p = CellDatatoPointData(Input=d3)
plotLine = PlotOverLine(Input=c2p)
plotLine.Point1 = [-3.0, 0.0, 0.0]
plotLine.Point2 = [-1.0, 0.0, 0.0]
plotLine.SamplingPattern = 'Sample At Cell Boundaries'

SaveData("line_export.csv", proxy=plotLine, Precision=10)
print("Saved: line_export.csv")

coeff_reader = XMLPartitionedUnstructuredGridReader(FileName=["coeffs_-1.pvtu"])
SaveData("coeffs_export.csv", proxy=coeff_reader, Precision=10)
print("Saved: coeffs_export.csv")
