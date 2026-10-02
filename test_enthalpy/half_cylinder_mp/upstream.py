import sys, numpy as np
from paraview.simple import *
from paraview import servermanager as sm
from paraview.numpy_support import vtk_to_numpy
HINF=283.5
for f in sys.argv[1:]:
    r=XMLPartitionedUnstructuredGridReader(FileName=[f]); r.UpdatePipeline()
    bl=[]
    def w(o):
        if o.IsA('vtkCompositeDataSet'):
            it=o.NewIterator(); it.InitTraversal()
            while not it.IsDoneWithTraversal(): w(it.GetCurrentDataObject()); it.GoToNextItem()
        else: bl.append(o)
    w(sm.Fetch(r))
    get=lambda n: np.concatenate([vtk_to_numpy(b.GetCellData().GetArray(n)) for b in bl])
    H=get('H'); c=get('Centroid'); rho=get('Density')
    rad=np.hypot(c[:,0],c[:,1])
    # "undisturbed upstream": cells where the density is still the freestream value to 0.1%
    up = np.abs(rho-1.0) < 1e-3
    dev=np.abs(H-HINF)/HINF
    print("%-34s cells upstream(rho==1 to 0.1%%)=%4d  median|dH| there=%.3e  max=%.3e"%(
        f.split('/')[1], up.sum(), np.median(dev[up]) if up.any() else np.nan,
        dev[up].max() if up.any() else np.nan))
    print("%-34s  rho range in the 'freestream' shell r>2.2 : %.10f .. %.10f"%(
        "", rho[rad>2.2].min(), rho[rad>2.2].max()))
    Delete(r)
