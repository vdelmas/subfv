#!/usr/bin/env python3
"""Run the Mach-dependent odd-even shock test on its baseline mesh."""
import argparse
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "..", "..", "build", "subfvsilvia")
GEO = os.path.join(HERE, "oddeven.geo")
GAMMA = 7.0 / 5.0
DX = 0.5
PERTURBATION = 1.0e-6


def shock_state(mach):
    """Post-shock primitive state for a right-moving shock into (1,0,0,0,1)."""
    m2 = mach * mach
    rho = (GAMMA + 1.0) * m2 / ((GAMMA - 1.0) * m2 + 2.0)
    pressure = 1.0 + 2.0 * GAMMA / (GAMMA + 1.0) * (m2 - 1.0)
    speed = mach * math.sqrt(GAMMA)
    velocity = speed * (1.0 - 1.0 / rho)
    return rho, velocity, pressure


def perturb_centre_row(msh):
    """Apply the alternating x displacement to MSH4 nodes at y=0."""
    with open(msh) as f:
        lines = f.readlines()
    start = lines.index("$Nodes\n")
    end = lines.index("$EndNodes\n", start)
    header = lines[start + 1].split()
    n_blocks = int(header[0])
    index = start + 2
    changed = 0
    for _ in range(n_blocks):
        dim, _entity, parametric, n_nodes = map(int, lines[index].split())
        index += 1
        tags_read = 0
        while tags_read < n_nodes:
            tags_read += len(lines[index].split())
            index += 1
        for _node in range(n_nodes):
            fields = lines[index].split()
            x, y = float(fields[0]), float(fields[1])
            if abs(y) < 1.0e-10:
                column = round(x / DX)
                if 0 < column < round(400.0 / DX):
                    x += PERTURBATION if column % 2 == 0 else -PERTURBATION
                    fields[0] = f"{x:.16g}"
                    lines[index] = " ".join(fields) + "\n"
                    changed += 1
            index += 1
            if parametric:
                # MSH4 parametric nodes append `dim` parameters on this line.
                pass
    if changed == 0:
        raise RuntimeError(f"no y=0 nodes were perturbed in {msh}")
    with open(msh, "w") as f:
        f.writelines(lines)
    print(f"perturbed {changed} centre-row mesh nodes by +/-{PERTURBATION:g}")


def run_quiet(command, cwd):
    """Keep successful mesh-tool output quiet, but show it if a tool fails."""
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
    ap.add_argument("--mach", type=float, default=6.0)
    ap.add_argument("--order", type=int, default=1)
    ap.add_argument("--scheme", default="multi_point")
    ap.add_argument("--quad", default="simpson", choices=["simpson", "gauss"])
    ap.add_argument("--subface-quad", "--subface_quad", type=parse_bool, nargs="?",
                    const=True, default=True,
                    help="use subfaces (true) or full faces (false)")
    ap.add_argument("--no-subface-quad", dest="subface_quad", action="store_false")
    ap.add_argument("--tmax", type=float, default=50.0)
    ap.add_argument("--cfl", type=float, default=0.5)
    ap.add_argument("--np", type=int, default=8)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    if args.mach <= 1.0:
        ap.error("--mach must be greater than 1")

    mode = "sf" if args.subface_quad else "f"
    name = f"{args.scheme}_{mode}_o{args.order}_{args.quad}_M{args.mach:g}{args.tag}"
    run = os.path.join(HERE, "outputs", name)
    os.makedirs(run, exist_ok=True)
    rho, velocity, pressure = shock_state(args.mach)
    with open(os.path.join(HERE, "input.template")) as f:
        template = f.read()
    replacements = {
        "ORDER_PLACEHOLDER": str(args.order), "SCHEME_PLACEHOLDER": args.scheme,
        "QUAD_PLACEHOLDER": args.quad,
        "SUBFACE_QUAD_PLACEHOLDER": ".true." if args.subface_quad else ".false.",
        "CFL_PLACEHOLDER": repr(args.cfl), "TMAX_PLACEHOLDER": repr(args.tmax),
        "RHO2_PLACEHOLDER": f"{rho:.16g}", "U2_PLACEHOLDER": f"{velocity:.16g}",
        "P2_PLACEHOLDER": f"{pressure:.16g}",
    }
    for key in sorted(replacements, key=len, reverse=True):
        value = replacements[key]
        template = template.replace(key, value)
    with open(os.path.join(run, "input_data.f"), "w") as f:
        f.write(template)
    run_quiet(["gmsh", "-3", GEO, "-o", "oddeven.msh"], cwd=run)
    perturb_centre_row(os.path.join(run, "oddeven.msh"))
    if args.np > 1:
        run_quiet(["subfv-gmsh", "-3", "oddeven.msh", "-part", str(args.np),
                   "-part_split", "-part_ghosts"], cwd=run)
    print(f"[{name}] RH inflow: rho={rho:.8g}, u={velocity:.8g}, p={pressure:.8g}", flush=True)
    command = [BIN, "input_data.f"] if args.np == 1 else ["mpirun", "-np", str(args.np), BIN, "input_data.f"]
    with open(os.path.join(run, "log.txt"), "w") as log:
        result = subprocess.run(command, cwd=run, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        sys.exit(f"subfvsilvia failed; see {run}/log.txt")
    print(f"completed: {run}")


if __name__ == "__main__":
    main()
