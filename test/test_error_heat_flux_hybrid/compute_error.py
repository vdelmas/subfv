import numpy as np
import csv
import os

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
COEFFS_FILE = "coeffs_export.csv"

RHO_INF = 1e-3
V_INF   = 5000.0

ref_st = np.loadtxt(os.path.join(SCRIPT_DIR, "fun3D_st.csv"), delimiter=',')
ref_cp = np.loadtxt(os.path.join(SCRIPT_DIR, "fun3D_cp.csv"), delimiter=',')

theta_ref_st, st_ref = ref_st[:, 0], ref_st[:, 1]
theta_ref_cp, cp_ref = ref_cp[:, 0], ref_cp[:, 1]

theta_sim = []
q_sim     = []
cp_sim    = []
with open(COEFFS_FILE) as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            theta_sim.append(float(row['Theta3']))
            q_sim.append(float(row['q']))
            cp_sim.append(float(row['Cp']))
        except (ValueError, KeyError):
            pass

theta_sim = np.abs(np.array(theta_sim))
q_sim     = np.array(q_sim)
cp_sim    = np.array(cp_sim)
st_sim    = 2.0 * q_sim / (RHO_INF * V_INF**3)

mask      = (theta_sim > 0) & (theta_sim <= np.pi / 2)
theta_sim = theta_sim[mask]
st_sim    = st_sim[mask]
cp_sim    = cp_sim[mask]

idx       = np.argsort(theta_sim)
theta_sim = theta_sim[idx]
st_sim    = st_sim[idx]
cp_sim    = cp_sim[idx]

st_ref_interp = np.interp(theta_sim, theta_ref_st, st_ref)
cp_ref_interp = np.interp(theta_sim, theta_ref_cp, cp_ref)

err_st = np.sqrt(np.mean((st_sim - st_ref_interp)**2)) / np.mean(np.abs(st_ref_interp)) * 100
err_cp = np.sqrt(np.mean((cp_sim - cp_ref_interp)**2)) / np.mean(np.abs(cp_ref_interp)) * 100

print(f"{err_st:.6e} {err_cp:.6e}")
