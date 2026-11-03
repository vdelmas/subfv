#!/usr/bin/env python3
"""Run the 2D Cartesian Sedov blast test with subfvsilvia."""
import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "..", "..", "build", "subfvsilvia")
GEO = os.path.join(HERE, "sedov.geo")


def run_quiet(command, cwd):
    result = subprocess.run(command, cwd=cwd, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        sys.stderr.write(result.stdout)
        result.check_returncode()


def parse_bool(value):
    normalized = value.lower()
    if normalized in {"true", "t", "1", "yes", "y"}:
        return True
    if normalized in {"false", "f", "0", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError("expected true or false")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--order", type=int, default=1)
    ap.add_argument("--scheme", default="multi_point")
    ap.add_argument("--quad", default="simpson", choices=["simpson", "gauss"])
    ap.add_argument("--subface-quad", "--subface_quad", type=parse_bool, nargs="?",
                    const=True, default=True,
                    help="use subfaces (true) or full faces (false)")
    ap.add_argument("--no-subface-quad", dest="subface_quad", action="store_false")
    ap.add_argument("--tmax", type=float, default=1.0)
    ap.add_argument("--cfl", type=float, default=0.5)
    ap.add_argument("--np", type=int, default=8)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()

    mode = "sf" if args.subface_quad else "f"
    name = f"{args.scheme}_{mode}_o{args.order}_{args.quad}{args.tag}"
    run = os.path.join(HERE, "outputs", name)
    os.makedirs(run, exist_ok=True)
    with open(os.path.join(HERE, "input.template")) as f:
        template = f.read()
    replacements = {
        "ORDER_PLACEHOLDER": str(args.order),
        "SCHEME_PLACEHOLDER": args.scheme,
        "QUAD_PLACEHOLDER": args.quad,
        "SUBFACE_QUAD_PLACEHOLDER": ".true." if args.subface_quad else ".false.",
        "CFL_PLACEHOLDER": repr(args.cfl),
        "TMAX_PLACEHOLDER": repr(args.tmax),
    }
    for key in sorted(replacements, key=len, reverse=True):
        template = template.replace(key, replacements[key])
    with open(os.path.join(run, "input_data.f"), "w") as f:
        f.write(template)

    run_quiet(["gmsh", "-3", GEO, "-o", "sedov.msh"], cwd=run)
    if args.np > 1:
        run_quiet(["subfv-gmsh", "-3", "sedov.msh", "-part", str(args.np),
                   "-part_split", "-part_ghosts"], cwd=run)

    command = ([BIN, "input_data.f"] if args.np == 1 else
               ["mpirun", "-np", str(args.np), BIN, "input_data.f"])
    with open(os.path.join(run, "log.txt"), "w") as log:
        result = subprocess.run(command, cwd=run, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        sys.exit(f"subfvsilvia failed; see {run}/log.txt")
    print(f"completed Sedov: {run}")


if __name__ == "__main__":
    main()
