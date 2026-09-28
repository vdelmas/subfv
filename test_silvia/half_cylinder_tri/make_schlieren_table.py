# Assembles outputs/o<k>/schlieren.png into one labelled table image.
# Usage: python3 make_schlieren_table.py 1 2 3   (after render_schlieren.py)
import sys
from PIL import Image, ImageDraw, ImageFont

orders = sys.argv[1:] or ["1", "2", "3"]
imgs = [Image.open(f"outputs/o{o}/schlieren.png").convert("RGB") for o in orders]
w, h = imgs[0].size
header, gap = 70, 20

try:
    font = ImageFont.truetype("DejaVuSans.ttf", 40)
except OSError:
    font = ImageFont.load_default()

table = Image.new("RGB", (len(imgs) * w + (len(imgs) - 1) * gap, h + header), "white")
draw = ImageDraw.Draw(table)
for k, (o, im) in enumerate(zip(orders, imgs)):
    x0 = k * (w + gap)
    table.paste(im, (x0, header))
    label = f"ordre {o}"
    tw = draw.textlength(label, font=font)
    draw.text((x0 + (w - tw) / 2, 12), label, fill="black", font=font)

table.save("schlieren_table.png")
print("saved schlieren_table.png", table.size)
