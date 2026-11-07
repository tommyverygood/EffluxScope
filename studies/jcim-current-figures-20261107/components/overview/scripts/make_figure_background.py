"""Embed the approved illustration with editable, centered panel headings.

The PNG is a conceptual illustration, not experimental data. This wrapper
keeps its native pixels and aspect ratio at 6.5 inches of manuscript width.
Only the PDF heading layer is reflowed; the raster source remains unchanged.
Requires reportlab and Pillow; no model execution or data fitting is involved.
"""
from pathlib import Path
from PIL import Image
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'assets/Figure_background_approved.png'
OUTPUT = ROOT / 'outputs/Figure_background.pdf'


def main():
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with Image.open(SOURCE) as image:
        pixel_width, pixel_height = image.size
    width = 6.5 * 72
    height = width * pixel_height / pixel_width
    c = canvas.Canvas(str(OUTPUT), pagesize=(width, height), pageCompression=1)
    c.setTitle('Reliability challenges in ABCB1 efflux prediction')
    c.drawImage(ImageReader(str(SOURCE)), 0, 0, width=width, height=height,
                preserveAspectRatio=True)
    # The original illustration is retained intact underneath this PDF text
    # layer. Cover only old headings, restore border segments, then place
    # letter + title as one centered bold group.
    font_path = Path('/System/Library/Fonts/Supplemental/Arial Bold.ttf')
    font = 'Helvetica-Bold'
    if font_path.exists():
        pdfmetrics.registerFont(TTFont('Arial-Bold', str(font_path)))
        font = 'Arial-Bold'
    scale = width / pixel_width
    def rect_px(x, y, w, h):
        c.rect(x * scale, height - (y + h) * scale,
               w * scale, h * scale, fill=1, stroke=0)
    def line_px(x0, y0, x1, y1):
        c.line(x0 * scale, height - y0 * scale,
               x1 * scale, height - y1 * scale)
    panels = [
        (14, 620, 84, 168, 'A', ['ABCB1-mediated efflux']),
        (676, 1261, 84, 195, 'B', ['Different evidence,', 'different scores']),
        (1327, 2050, 84, 168, 'C', ['Two evaluation questions']),
    ]
    title_size = 8.5
    for x0, x1, top, bottom, letter, lines in panels:
        c.setFillColorRGB(1, 1, 1)
        rect_px(x0 - 8, top, x1 - x0 + 16, bottom - top)
        c.setStrokeColorRGB(0.45, 0.61, 0.94)
        c.setLineWidth(0.65)
        line_px(x0, 147, x1, 147)
        line_px(x0, 147, x0, bottom + 1)
        line_px(x1, 147, x1, bottom + 1)
        line1 = f'{letter}  {lines[0]}'
        text_w = max([pdfmetrics.stringWidth(line1, font, title_size)] +
                     [pdfmetrics.stringWidth(line, font, title_size) for line in lines[1:]])
        center = (x0 + x1) * scale / 2
        # A white badge interrupts the top rule, matching the original panels.
        c.setFillColorRGB(1, 1, 1)
        badge_y0 = 94
        badge_y1 = 184 if len(lines) > 1 else 166
        c.rect(center - text_w / 2 - 3, height - badge_y1 * scale,
               text_w + 6, (badge_y1 - badge_y0) * scale, fill=1, stroke=0)
        c.setFillColorRGB(0, 0, 0)
        c.setFont(font, title_size)
        baseline = 134 if len(lines) > 1 else 150
        c.drawCentredString(center, height - baseline * scale, line1)
        for line_index, line in enumerate(lines[1:], start=1):
            c.drawCentredString(center, height - (baseline + 38 * line_index) * scale, line)
    c.showPage()
    c.save()
    print(f'{OUTPUT.relative_to(ROOT)}: {width:.2f} x {height:.2f} pt; {pixel_width / 6.5:.1f} ppi')


if __name__ == '__main__':
    main()
