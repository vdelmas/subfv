#!/usr/bin/env python3
"""Mesh-convergence study of subfvsilvia on the stationary isentropic vortex.

For each mesh size ls: mesh vortex.geo with gmsh, write input_data.f from
input.template, run subfvsilvia, read the final error printed by the code,
then tabulate the L1 errors and the observed orders.

The code prints err(rho,u,v,w,p) = sum over cells with |x|,|y| < 3 of
volume*|w_num - w_exact|, and the total volume of those cells; the mean L1
error is err/volume. The mesh size is h = sqrt(L^2 / n_cells), n_cells being
the number of 3D cells in the .msh (one layer).

Examples (from test_silvia/isentropic_vortex_convergence):
  python3 run_convergence.py --order 1 --scheme multi_point
  python3 run_convergence.py --order 2 --scheme two_wave --mesh quad --np 4
  python3 run_convergence.py --order 3 --quad gauss --ls 0.5 0.25 0.125 0.0625
"""
import argparse
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "..", "..", "build", "subfvsilvia")


def count_cells(msh):
    """Number of 3D elements in a MSH4 file (sum of the dim-3 blocks)."""
    with open(msh) as f:
        lines = f.read().splitlines()
    i = lines.index("$Elements") + 1
    n_blocks = int(lines[i].split()[0])
    i += 1
    n3 = 0
    for _ in range(n_blocks):
        dim, _tag, _etype, n = map(int, lines[i].split())
        if dim == 3:
            n3 += n
        i += 1 + n
    return n3


def read_error(log):
    """Last 'volume=' and 'err(rho,u,v,w,p)=' lines printed by the code."""
    with open(log) as f:
        text = f.read()
    vols = re.findall(r"volume=\s*(\S+)", text)
    errs = re.findall(r"err\(rho,u,v,w,p\)=\s*(.+)", text)
    if not vols or not errs:
        raise RuntimeError(f"no error found in {log}")
    volume = float(vols[-1])
    err = [float(x) for x in errs[-1].split()]
    return volume, err


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--order", type=int, default=1)
    ap.add_argument("--scheme", default="multi_point",
                    help="multi_point, two_wave, three_wave, ...")
    ap.add_argument("--quad", default="simpson", help="subface quadrature: simpson or gauss")
    def parse_bool(value):
        if isinstance(value, bool):
            return value
        normalized = value.lower()
        if normalized in {"true", "t", "1", "yes", "y"}:
            return True
        if normalized in {"false", "f", "0", "no", "n"}:
            return False
        raise argparse.ArgumentTypeError("expected true or false")

    ap.add_argument("--subface_quad", "--subface-quad", type=parse_bool, nargs="?", const=True,
                    default=True, help="use subfaces (true) or full faces (false)")
    ap.add_argument("--mesh", default="tri", choices=["tri", "quad"])
    ap.add_argument("--ls", type=float, nargs="+", default=[0.4, 0.2, 0.1],
                    help="target mesh sizes (finest last)")
    ap.add_argument("--L", type=float, default=10.0, help="domain side")
    ap.add_argument("--tmax", type=float, default=1.0)
    ap.add_argument("--cfl", type=float, default=0.5)
    ap.add_argument("--np", type=int, default=1, help="MPI ranks")
    ap.add_argument("--tag", default="", help="suffix of the output directory")
    args = ap.parse_args()

    if args.subface_quad==True:
        name = f"{args.scheme}_sf_o{args.order}_{args.quad}_{args.mesh}{args.tag}"
    else:
        name = f"{args.scheme}_f_o{args.order}_{args.quad}_{args.mesh}{args.tag}"
    outroot = os.path.join(HERE, "outputs", name)
    os.makedirs(outroot, exist_ok=True)
    with open(os.path.join(HERE, "input.template")) as f:
        template = f.read()

    rows = []
    for ls in args.ls:
        run = os.path.join(outroot, f"ls_{ls:g}")
        os.makedirs(run, exist_ok=True)
        inp = (template.replace("ORDER_PLACEHOLDER", str(args.order))
               .replace("SCHEME_PLACEHOLDER", args.scheme)
               .replace("SUBFACE_QUAD_PLACEHOLDER", ".true." if args.subface_quad else ".false.")
               .replace("QUAD_PLACEHOLDER", args.quad)
               .replace("CFL_PLACEHOLDER", repr(args.cfl))
               .replace("TMAX_PLACEHOLDER", repr(args.tmax)))
        with open(os.path.join(run, "input_data.f"), "w") as f:
            f.write(inp)

        subprocess.run(["gmsh", "-3", os.path.join(HERE, "vortex.geo"),
                        "-setnumber", "ls", repr(ls), "-setnumber", "L", repr(args.L),
                        "-setnumber", "quad", "1" if args.mesh == "quad" else "0",
                        "-o", "vortex.msh"],
                       cwd=run, check=True, stdout=subprocess.DEVNULL)
        n_cells = count_cells(os.path.join(run, "vortex.msh"))
        if args.np > 1:
            subprocess.run(["subfv-gmsh", "-3", "vortex.msh", "-part", str(args.np),
                            "-part_split", "-part_ghosts"],
                           cwd=run, check=True, stdout=subprocess.DEVNULL)

        print(f"[{name}] ls={ls:g}: {n_cells} cells ...", flush=True)
        with open(os.path.join(run, "log.txt"), "w") as log:
            res = subprocess.run(["mpirun", "-np", str(args.np), BIN, "input_data.f"],
                                 cwd=run, stdout=log, stderr=subprocess.STDOUT)
        if res.returncode != 0:
            sys.exit(f"subfvsilvia failed for ls={ls:g} (exit {res.returncode}), see {run}/log.txt")

        volume, err = read_error(os.path.join(run, "log.txt"))
        h = math.sqrt(args.L**2 / n_cells)
        rows.append((ls, h, n_cells, [e / volume for e in err]))

    # table: mean L1 error of rho, u, v, p and observed order between meshes
    names = ["rho", "u", "v", "p"]
    comp = [0, 1, 2, 4]
    dat = os.path.join(outroot, "convergence.dat")
    with open(dat, "w") as f:
        f.write(f"# {name}  tmax={args.tmax:g} cfl={args.cfl:g}  mean L1 error over |x|,|y|<3\n")
        f.write("# ls h n_cells " + " ".join(f"err_{n}" for n in names)
                + " " + " ".join(f"order_{n}" for n in names) + "\n")
        print(f"\n{name}  (tmax={args.tmax:g}, cfl={args.cfl:g})")
        print(f"{'h':>10} {'cells':>8} " + " ".join(f"{'L1 ' + n:>11}" for n in names)
              + "  " + " ".join(f"{'ord ' + n:>7}" for n in names))
        for k, (ls, h, n, e) in enumerate(rows):
            errs = [e[c] for c in comp]
            if k == 0:
                orders = [float("nan")] * len(comp)
            else:
                h0, e0 = rows[k - 1][1], [rows[k - 1][3][c] for c in comp]
                orders = [math.log(a / b) / math.log(h0 / h) if a > 0 and b > 0 else float("nan")
                          for a, b in zip(e0, errs)]
            f.write(f"{ls:g} {h:.6e} {n} " + " ".join(f"{x:.6e}" for x in errs)
                    + " " + " ".join(f"{x:.3f}" for x in orders) + "\n")
            print(f"{h:10.4e} {n:8d} " + " ".join(f"{x:11.4e}" for x in errs)
                  + "  " + " ".join(f"{x:7.2f}" for x in orders))
    print(f"\nwritten {dat}")


if __name__ == "__main__":
    main()
