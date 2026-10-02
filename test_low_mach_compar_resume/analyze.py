#!/usr/bin/env python3
"""Score every ZB_*_* variant (+ WIP baseline) against the 4 WIP2 acceptance
criteria, from the already-rsynced outputs/<scheme>/ directories. Prints a
plain-text table; not a paper figure (see gnuplot plots for those)."""
import numpy as np
import os
import re

ROOT = "/run/media/delmasv/scratch/Codes/subfv/subfv/test_low_mach_compar_resume"
SCHEMES = ["multi_point", "three_wave", "two_wave", "modified_three_wave",
           "multi_point_pressure", "WIP", "WIP2_NOLM",
           "WIP2_EJUMP", "WIP2_EDIV", "WIP2_EMAX", "WIP2_EPOS", "WIP2_EPOS2",
           "ZB_ARMD_LPP", "ZB_ARMDU_LPP", "ZB_ARMDMAT_LPP", "ZB_ARMDUMAT_LPP",
           "WIP2_ED1", "WIP2_ED2", "WIP2_MD", "WIP2_MD5", "WIP2_MDH"]


def read_two_col(path):
    if not os.path.exists(path):
        return None
    data = np.loadtxt(path, comments="#")
    return data


def gresho_preservation(case):
    out = {}
    for s in SCHEMES:
        d = os.path.join(ROOT, case, "outputs", s)
        fi = read_two_col(os.path.join(d, "velocity_profile_initial.dat"))
        ff = read_two_col(os.path.join(d, "velocity_profile_final.dat"))
        if fi is None or ff is None:
            out[s] = None
            continue
        # match by nearest radius (files share the same cell sampling/order)
        vi = fi[:, 1]
        vf = ff[:, 1]
        n = min(len(vi), len(vf))
        rmse = float(np.sqrt(np.mean((vf[:n] - vi[:n]) ** 2)))
        vmax = float(np.max(vi)) if len(vi) else 1.0
        out[s] = rmse / vmax if vmax > 0 else rmse
    return out


def convergence_order(case):
    out = {}
    for s in SCHEMES:
        log = os.path.join(ROOT, case, "outputs", s, "log.txt")
        if not os.path.exists(log):
            out[s] = None
            continue
        mas, l2s = [], []
        with open(log) as f:
            for line in f:
                m = re.search(r"Ma=\s*([\d.eE+-]+)\s+L2_rho=\s*([\d.eE+-]+)", line)
                if m:
                    mas.append(float(m.group(1)))
                    l2s.append(float(m.group(2)))
        if len(mas) < 2:
            out[s] = None
            continue
        mas = np.array(mas)
        l2s = np.array(l2s)
        order = np.polyfit(np.log10(mas), np.log10(l2s), 1)[0]
        out[s] = order
    return out


def sedov_status(case):
    out = {}
    analytic = os.path.join(ROOT, case,
                             "analytic_sedov_3D.dat" if "hex" in case else "analytic_sedov_2D.dat")
    an = None
    if os.path.exists(analytic):
        an = np.loadtxt(analytic)
    for s in SCHEMES:
        d = os.path.join(ROOT, case, "outputs", s)
        prof = os.path.join(d, "density_profile.dat")
        if not os.path.exists(prof):
            # distinguish "never ran" from "ran and died": a crashed run still
            # leaves log.txt behind, a job that never started leaves nothing
            if not os.path.exists(os.path.join(d, "log.txt")):
                out[s] = ("NOT-RUN", None)
            else:
                out[s] = ("CRASHED", None)
            continue
        data = read_two_col(prof)
        if data is None or len(data) == 0:
            out[s] = ("CRASHED", None)
            continue
        if an is None:
            out[s] = ("OK", None)
            continue
        r_an, rho_an = an[:, 1], an[:, 2]
        r_sim, rho_sim = data[:, 0], data[:, 1]
        mask = (r_sim >= r_an.min()) & (r_sim <= r_an.max())
        if mask.sum() == 0:
            out[s] = ("OK", None)
            continue
        rho_an_interp = np.interp(r_sim[mask], r_an, rho_an)
        rmse = float(np.sqrt(np.mean((rho_sim[mask] - rho_an_interp) ** 2)))
        out[s] = ("OK", rmse)
    return out


def half_cylinder_heatflux(case):
    ref_path = os.path.join(ROOT, case, "fun3D_st.csv")
    ref = np.loadtxt(ref_path, delimiter=",")
    x_ref = np.concatenate([ref[:, 0], -ref[:, 0]])
    y_ref = np.concatenate([ref[:, 1], ref[:, 1]])
    order = np.argsort(x_ref)
    x_ref, y_ref = x_ref[order], y_ref[order]

    def st(x):
        return 2.0 * x / (1e-3 * 5000.0 ** 3)

    out = {}
    for s in SCHEMES:
        f = os.path.join(ROOT, case, "outputs", s, "coeffs_export.csv")
        if not os.path.exists(f):
            out[s] = None
            continue
        data = np.genfromtxt(f, delimiter=",", names=True)
        theta = data["Theta3"]
        q = data["q"]
        st_sim = st(q)
        order2 = np.argsort(theta)
        theta_s, st_s = theta[order2], st_sim[order2]
        mask = (theta_s >= x_ref.min()) & (theta_s <= x_ref.max())
        if mask.sum() == 0:
            out[s] = None
            continue
        ref_interp = np.interp(theta_s[mask], x_ref, y_ref)
        rmse = float(np.sqrt(np.mean((st_s[mask] - ref_interp) ** 2)))
        out[s] = rmse
    return out


if __name__ == "__main__":
    print("=== 1) Gresho velocity preservation (RMSE(|V|final-|V|initial)/Vmax, lower=better) ===")
    for case in ["gresho_quad", "gresho_tri"]:
        res = gresho_preservation(case)
        print(f"-- {case} --")
        for s, v in res.items():
            print(f"  {s:16s} {v if v is None else f'{v:.4e}'}")

    print("\n=== 2) Convergence order (slope of log10(L2_rho) vs log10(Ma), target: -1 min, -2 ideal) ===")
    for case in ["convergence_gresho_quad", "convergence_gresho_tri", "convergence_gresho_voro"]:
        res = convergence_order(case)
        print(f"-- {case} --")
        for s, v in res.items():
            print(f"  {s:16s} {v if v is None else f'{v:.3f}'}")

    print("\n=== 3) Sedov robustness (status, density-profile RMSE vs analytic) ===")
    for case in ["sedov_tri", "sedov_hex"]:
        res = sedov_status(case)
        print(f"-- {case} --")
        for s, (status, rmse) in res.items():
            rmse_s = "" if rmse is None else f"{rmse:.4e}"
            print(f"  {s:16s} {status:10s} {rmse_s}")

    print("\n=== 4) Half-cylinder Mach-27 heat flux RMSE vs fun3D/LAURA reference (lower=better) ===")
    for case in ["half_cylinder_quad_ns", "half_cylinder_tri_ns", "half_cylinder_voro_ns"]:
        res = half_cylinder_heatflux(case)
        print(f"-- {case} --")
        for s, v in res.items():
            print(f"  {s:16s} {v if v is None else f'{v:.4e}'}")
