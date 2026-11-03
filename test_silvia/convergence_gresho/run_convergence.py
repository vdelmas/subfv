#!/usr/bin/env python3
"""Mesh-convergence study for the stationary Gresho vortex."""

import argparse
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "..", "..", "build", "subfvsilvia")
GEO = os.path.join(HERE, "gresho.geo")


def count_cells(msh):
    """Count 3D elements in a Gmsh MSH4 file."""
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
    with open(log) as f:
        text = f.read()
    vols = re.findall(r"volume=\s*(\S+)", text)
    errs = re.findall(r"err\(rho,u,v,w,p\)=\s*(.+)", text)
    if not vols or not errs:
        raise RuntimeError(f"no error found in {log}")
    return float(vols[-1]), [float(x) for x in errs[-1].split()]


def parse_bool(value):
    normalized = value.lower()
    if normalized in {"true", "t", "1", "yes", "y"}:
        return True
    if normalized in {"false", "f", "0", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError("expected true or false")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--order", type=int, default=2)
    ap.add_argument("--scheme", default="multi_point")
    ap.add_argument("--quad", default="gauss", choices=["simpson", "gauss"])
    ap.add_argument("--subface_quad", "--subface-quad", type=parse_bool, nargs="?", const=True,
                    default=True, help="use subfaces (true) or full faces (false)")
    ap.add_argument("--ls", type=float, nargs="+", default=[0.04, 0.02, 0.01])
    ap.add_argument("--tmax", type=float, default=1.0)
    ap.add_argument("--cfl", type=float, default=0.5)
    ap.add_argument("--np", type=int, default=1)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()

    mode = "sf" if args.subface_quad else "f"
    name = f"{args.scheme}_{mode}_o{args.order}_{args.quad}{args.tag}"
    outroot = os.path.join(HERE, "outputs", name)
    os.makedirs(outroot, exist_ok=True)
    with open(os.path.join(HERE, "input.template")) as f:
        template = f.read()

    rows = []
    previous_cell_count = None
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

        subprocess.run(["gmsh", "-3", GEO, "-setnumber", "ls", repr(ls), "-o", "gresho.msh"],
                       cwd=run, check=True, stdout=subprocess.DEVNULL)
        n_cells = count_cells(os.path.join(run, "gresho.msh"))
        if n_cells == previous_cell_count:
            sys.exit(f"ls={ls:g} produced the same {n_cells} cells as the previous mesh; "
                     "check that Gmsh changed the mesh resolution")
        previous_cell_count = n_cells
        if args.np > 1:
            subprocess.run(["subfv-gmsh", "-3", "gresho.msh", "-part", str(args.np),
                            "-part_split", "-part_ghosts"],
                           cwd=run, check=True, stdout=subprocess.DEVNULL)

        print(f"[{name}] ls={ls:g}: {n_cells} cells ...", flush=True)
        with open(os.path.join(run, "log.txt"), "w") as log:
            result = subprocess.run(["mpirun", "-np", str(args.np), BIN, "input_data.f"],
                                    cwd=run, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            sys.exit(f"subfvsilvia failed for ls={ls:g}; see {run}/log.txt")

        volume, error = read_error(os.path.join(run, "log.txt"))
        h = math.sqrt(1.0 / n_cells)
        rows.append((ls, h, n_cells, [value / volume for value in error]))

    components = [0, 1, 2, 4]
    labels = ["rho", "u", "v", "p"]
    datafile = os.path.join(outroot, "convergence.dat")
    with open(datafile, "w") as f:
        f.write(f"# {name}  tmax={args.tmax:g} cfl={args.cfl:g}\n")
        f.write("# ls h n_cells " + " ".join(f"err_{v}" for v in labels)
                + " " + " ".join(f"order_{v}" for v in labels) + "\n")
        print(f"\n{name}  (tmax={args.tmax:g}, cfl={args.cfl:g})")
        print(f"{'h':>10} {'cells':>8} " + " ".join(f"{'L1 ' + v:>11}" for v in labels)
              + "  " + " ".join(f"{'ord ' + v:>7}" for v in labels))
        for index, (ls, h, cells, errors) in enumerate(rows):
            values = [errors[c] for c in components]
            if index == 0:
                orders = [float("nan")] * len(components)
            else:
                h0 = rows[index - 1][1]
                prev = [rows[index - 1][3][c] for c in components]
                orders = [math.log(a / b) / math.log(h0 / h) if a > 0 and b > 0 else float("nan")
                          for a, b in zip(prev, values)]
            f.write(f"{ls:g} {h:.6e} {cells} " + " ".join(f"{v:.6e}" for v in values)
                    + " " + " ".join(f"{v:.3f}" for v in orders) + "\n")
            print(f"{h:10.4e} {cells:8d} " + " ".join(f"{v:11.4e}" for v in values)
                  + "  " + " ".join(f"{v:7.2f}" for v in orders))
    print(f"\nwritten {datafile}")


if __name__ == "__main__":
    main()
