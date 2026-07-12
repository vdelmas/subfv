from paraview.simple import *

coeff_reader = XMLPartitionedUnstructuredGridReader(FileName=["coeffs_-1.pvtu"])
SaveData("coeffs_export.csv", proxy=coeff_reader, Precision=10)
