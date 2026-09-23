#!/usr/bin/env python3
"""Assemble the eight rendered panels into one comparison figure.

Grid: 2 rows (schlieren / total enthalpy) x 4 columns (quad three_wave, quad
three_wave_enthalpy, tri three_wave, tri three_wave_enthalpy). Every panel was rendered by
render.py on the SAME colour range, so the colour bar is drawn once per row, on the right.
render.py lays each PNG out as [field | reserved white margin holding the bar], so the first
three columns are cropped to the field and only the last keeps its margin.
"""
import os
from PIL import Image, ImageDraw, ImageFont

FIELD_W, BAR_W, H = 500, 460, 1000
# All four solvers, both meshes, in one strip. Every run in outputs/ comes from the SAME
# executable (subfvns) with the same residual-driven stopping, so the columns really are
# comparable -- that is why three_wave / three_wave_enthalpy are re-run here rather than
# reused from ../half_cylinder, which marches in time through euler_ho instead.
SCHEMES = ["three_wave", "three_wave_enthalpy", "multi_point", "multi_point_enthalpy"]
COLS = [("%s_%s" % (m, s), "%s - %s" % (m, s)) for m in ("quad", "tri") for s in SCHEMES]
ROWS = [("schlieren", "schlieren  log(|grad rho|+1)"),
        ("enthalpy", "total enthalpy H   (h_inf = 283.5)")]
TITLE_H, ROWLAB_W, PAD = 54, 46, 8


def font(sz):
    for p in ("/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf",
              "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
              "/usr/share/fonts/dejavu-sans-fonts/DejaVuSans.ttf"):
        if os.path.exists(p):
            return ImageFont.truetype(p, sz)
    return ImageFont.load_default()


W = ROWLAB_W + len(COLS) * FIELD_W + BAR_W + 2 * PAD
HT = TITLE_H + 2 * (H + PAD)
out = Image.new("RGB", (W, HT), "white")
d = ImageDraw.Draw(out)
f_col, f_row = font(26), font(26)

for c, (tag, label) in enumerate(COLS):
    x = ROWLAB_W + c * FIELD_W
    d.text((x + FIELD_W // 2, TITLE_H // 2), label, fill="black", font=f_col, anchor="mm")

for r, (field, rowlabel) in enumerate(ROWS):
    y = TITLE_H + r * (H + PAD)
    for c, (tag, _) in enumerate(COLS):
        im = Image.open("figures/%s_%s.png" % (tag, field))
        box = (0, 0, FIELD_W, H) if c < len(COLS) - 1 else (0, 0, FIELD_W + BAR_W, H)
        out.paste(im.crop(box), (ROWLAB_W + c * FIELD_W, y))
    lab = Image.new("RGB", (H, ROWLAB_W), "white")
    ImageDraw.Draw(lab).text((H // 2, ROWLAB_W // 2), rowlabel, fill="black",
                             font=f_row, anchor="mm")
    out.paste(lab.rotate(90, expand=True), (0, y))

out.save("figures/comparison.png")
print("saved figures/comparison.png  (%dx%d)" % out.size)
