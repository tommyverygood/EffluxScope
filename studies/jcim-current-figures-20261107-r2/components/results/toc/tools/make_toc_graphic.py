#!/usr/bin/env python3
"""Table-of-contents graphic for the EffluxScope method-evaluation manuscript.

ACS specification (toc_abstract_graphics_guidelines.pdf, 28 Feb 2024): at most 3.25 x 1.75 in,
sans-serif type no smaller than 6 pt, TIF at 300 dpi (colour) or EPS (RGB, fonts embedded).
The graphic is a schematic: bar lengths are read from the frozen equal-count metrics
(data_for_figures/metrics.csv) but no numbers are printed, and it does not reuse any panel of
Figures 1-5. Colours follow the project palette (palette.SEMANTIC): weak-label purple,
endpoint-aligned green, accepted/policy blue, abstained grey.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon, Circle
from matplotlib.path import Path as MPath

TOOLS = Path(__file__).resolve().parent
ROOT = TOOLS.parent
sys.path.insert(0, str(TOOLS))
import palette
from palette import C, tint

OUT = ROOT / "toc_graphic"
W, H = 3.25, 1.75                                  # inches, ACS maximum
INK, MID, LIGHT = "#262626", "#666666", "#9A9A9A"  # neutral greys (layout colours)
CARD, CARD_EDGE = "#F4F4F4", "#D6D6D6"

FS_HEAD, FS_BODY, FS_SUB = 6.8, 6.4, 6.1           # all >= 6 pt (ACS floor)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.unicode_minus": False,
    "savefig.bbox": "standard", "savefig.pad_inches": 0.0,   # keep the exact 3.25 x 1.75 in page
})

# ---- values behind the bar lengths (frozen equal-count comparison, 453 retained) -------------
M = pd.read_csv(ROOT / "data_for_figures" / "metrics.csv", encoding="utf-8-sig")
sel = M[(M.analysis_role == "primary_matched_count") & (M.version == "identity_audited_v2")]
bacc = dict(zip(sel.policy, sel.balanced_accuracy))
V01, V03, FULL = bacc["v01_margin"], bacc["v03_margin"], bacc["full_v04_pool_primary_margin"]
assert V01 < FULL < V03, (V01, FULL, V03)          # ordering the graphic depicts


def rbox(ax, x, y, w, h, fc, ec="none", lw=0.0, r=0.045, z=1):
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                       fc=fc, ec=ec, lw=lw, zorder=z)
    ax.add_patch(p); return p


def poisson_ellipse(n, cx, cy, rx, ry, dmin, seed):
    rng = np.random.default_rng(seed); pts = []; tries = 0
    while len(pts) < n and tries < 40000:
        tries += 1
        x, y = rng.uniform(-1, 1, 2)
        if x * x + y * y > 1:
            continue
        p = (cx + x * rx, cy + y * ry)
        if all((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 >= dmin ** 2 for q in pts):
            pts.append(p)
    return pts


def build(adc=False):
    blue, purple, green, grey = C["blue"], C["purple"], C["green"], C["grey"]
    palette.assert_palette([blue, purple, green, grey, tint("blue", 0.30)])
    fig = plt.figure(figsize=(W, H), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")

    # ---- cards ------------------------------------------------------------------------------
    cy0, ch = 0.335, 1.365
    cards = {"A": (0.06, 0.94), "B": (1.09, 0.90), "C": (2.08, 1.11)}
    for k, (x0, w) in cards.items():
        rbox(ax, x0, cy0, w, ch, CARD, CARD_EDGE, 0.5, r=0.05, z=0)
    ytop, ysub = 1.585, 1.465
    heads = {"A": ("2,624 molecules", "measured efflux"),
             "B": ("Abstention", "composite policy"),
             "C": ("Equal retained count", "balanced accuracy")}
    for k, (x0, w) in cards.items():
        ax.text(x0 + w / 2, ytop, heads[k][0], ha="center", va="center", fontsize=FS_HEAD,
                fontweight="bold", color=INK)
        ax.text(x0 + w / 2, ysub, heads[k][1], ha="center", va="center", fontsize=FS_SUB, color=MID)

    # ---- card A: molecules -----------------------------------------------------------------
    ax.add_patch(mpl.patches.Ellipse((0.53, 0.925), 0.80, 0.80, fc="#EAEAEA", ec="none", zorder=1))
    cloud = poisson_ellipse(54, 0.53, 0.925, 0.355, 0.365, 0.064, seed=4)
    for (x, y) in cloud:
        ax.add_patch(Circle((x, y), 0.0165, fc="#5B5B5B", ec="none", zorder=2))
    teal = C["teal"]
    if adc:                                          # optional: five ADC payload structures of the case study
        palette.assert_palette([teal])
        rng = np.random.default_rng(11)
        pick = [cloud[i] for i in rng.choice([i for i, (x, y) in enumerate(cloud) if abs(x - 0.53) < 0.25 and abs(y - 0.925) < 0.26], 5, replace=False)]
        for (x, y) in pick:
            ax.add_patch(Circle((x, y), 0.032, fc="none", ec=teal, lw=1.0, zorder=3))
        ax.add_patch(Circle((0.255, 0.395), 0.032, fc="none", ec=teal, lw=1.0, zorder=3))
        ax.add_patch(Circle((0.255, 0.395), 0.0165, fc="#5B5B5B", ec="none", zorder=3))
        ax.text(0.31, 0.395, "ADC payloads", ha="left", va="center", fontsize=FS_SUB, color=MID)

    # ---- card B: policy funnel ---------------------------------------------------------------
    X0, XN, XE, YC = 1.15, 1.62, 1.91, 0.93
    HM, HN = 0.38, 0.07                              # half-heights at the mouth and at the neck
    fx = [(X0, YC + HM), (XN, YC + HN), (XE, YC + HN), (XE, YC - HN), (XN, YC - HN), (X0, YC - HM)]
    ax.add_patch(Polygon(fx, closed=True, fc=tint("blue", 0.30), ec=blue, lw=0.9, joinstyle="round", zorder=2))
    def half(x):                                     # half-height of the cone at x
        return HM + (HN - HM) * (x - X0) / (XN - X0)
    # (x, t): t in [-1, 1] is the relative offset from the centre line; |t| small = passes the neck
    passing = [(1.23, 0.10), (1.30, -0.12), (1.37, 0.14), (1.44, -0.10), (1.51, 0.08)] + [(1.67, 0), (1.76, 0), (1.85, 0)]
    blocked = [(1.24, 0.72), (1.31, 0.58), (1.40, 0.66), (1.36, 0.34),
               (1.24, -0.72), (1.32, -0.60), (1.40, -0.68), (1.30, -0.92)]
    path = MPath(fx)
    for (x, t) in blocked:
        y = YC + t * half(x) if x < XN else YC + t * HN
        assert path.contains_point((x, y)), (x, y)
        ax.add_patch(Circle((x, y), 0.0165, fc=LIGHT, ec="none", zorder=3))
    for (x, t) in passing:
        y = YC + t * half(x) if x < XN else YC + t * HN
        assert path.contains_point((x, y)), (x, y)
        ax.add_patch(Circle((x, y), 0.0165, fc=blue, ec="none", zorder=3))
    if adc:
        for (x, t) in (passing[1], passing[6]) + (blocked[1], blocked[4], blocked[6]):
            y = YC + t * half(x) if x < XN else YC + t * HN
            ax.add_patch(Circle((x, y), 0.032, fc="none", ec=teal, lw=1.0, zorder=4))
    ax.text(1.775, 1.075, "retain", ha="center", va="center", fontsize=FS_SUB, color=blue, fontweight="bold")
    ax.text(1.40, 0.465, "abstain", ha="center", va="center", fontsize=FS_SUB, color=MID, fontweight="bold")

    # ---- arrows between cards ----------------------------------------------------------------
    for xa, xb in ((1.008, 1.082), (1.998, 2.072)):
        ax.annotate("", xy=(xb, 0.93), xytext=(xa, 0.93),
                    arrowprops=dict(arrowstyle="-|>", color=MID, lw=1.0, mutation_scale=6.5,
                                    shrinkA=0, shrinkB=0), zorder=4)

    # ---- card C: ranking at equal retained count ---------------------------------------------
    x0c, L = 2.17, 1.30                             # baseline (chance) and length of 0.5 -> 1.0
    rows = [("Assay-aware score", V03, green, 1.325), ("Weak-label + policy", FULL, blue, 1.025),
            ("Weak-label score", V01, purple, 0.725)]
    xpol = x0c + (FULL - 0.5) / 0.5 * L
    ax.plot([x0c, x0c], [0.50, 1.27], color=LIGHT, lw=0.7, zorder=1.5)
    for lab, v, col, yl in rows:
        yb = yl - 0.115
        ax.text(x0c + 0.03, yl, lab, ha="left", va="center", fontsize=FS_BODY, color=INK)
        rbox(ax, x0c, yb - 0.045, (v - 0.5) / 0.5 * L, 0.09, col, r=0.04, z=2)
        # guide at the level reached by the policy, drawn only beside the bars (never through text)
        ax.plot([xpol, xpol], [yb - 0.075, yb + 0.075], color=blue, lw=0.9, ls=(0, (1.1, 1.3)), zorder=3)
    ax.text(x0c + 0.03, 0.435, "chance", ha="left", va="center", fontsize=FS_SUB, color=MID)

    # ---- take-home banner ---------------------------------------------------------------------
    rbox(ax, 0.06, 0.055, 3.13, 0.21, blue, r=0.05, z=0)
    ax.text(W / 2, 0.16, "Abstention did not replace endpoint alignment", ha="center", va="center",
            fontsize=6.9, fontweight="bold", color="white")
    return fig


def text_overlaps(fig):
    fig.canvas.draw(); r = fig.canvas.get_renderer(); items = []
    for t in fig.findobj(mpl.text.Text):
        if t.get_text().strip() and t.get_visible():
            items.append((t.get_text(), t.get_window_extent(r)))
    bad = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            if items[i][1].overlaps(items[j][1]):
                bad.append((items[i][0], items[j][0]))
    small = [t.get_text() for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_fontsize() < 6.0]
    return bad, small


def save_all(fig, stem):
    import io
    from PIL import Image
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.eps", format="eps")
    fig.savefig(OUT / f"{stem}_preview.png", dpi=600)
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=300); buf.seek(0)
    im = Image.open(buf).convert("RGB")                 # RGB, no alpha channel
    im.save(OUT / f"{stem}.tif", dpi=(300, 300), compression="tiff_lzw")
    return im.size


if __name__ == "__main__":
    for adc, stem in ((False, "TOC_graphic"), (True, "TOC_graphic_ADC")):
        f = build(adc=adc); bad, small = text_overlaps(f)
        print(stem, "| text overlaps:", bad, "| text below 6 pt:", small, "| TIF px:", save_all(f, stem))
    print("BAcc behind the bar lengths (not printed on the graphic):", round(V01, 3), round(FULL, 3), round(V03, 3))
