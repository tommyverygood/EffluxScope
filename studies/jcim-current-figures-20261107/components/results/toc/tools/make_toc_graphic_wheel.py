#!/usr/bin/env python3
"""TOC graphic, wheel style (layout chosen by the authors from a reference they like):
a central three-sector wheel with pale-blue/pink fills and a white hub, four corner panels (two tinted) with black titles,
hexagon badges, blue-bordered circles and thick curved arrows. Only the layout idea is borrowed;
all artwork is drawn here, without logos (ACS does not allow them in TOC graphics).

Content: the three ways compared at an equal retained count (weak-label score; weak-label score + abstention
policy; assay-aware score ranked by its own margin) and the resulting balanced-accuracy rank (low / mid / high).
The rank order is asserted against data_for_figures/metrics.csv at build time; no numeric result is printed.

Specification (ACS, 28 Feb 2024): at most 3.25 x 1.75 in, type >= 6 pt, TIF 300 dpi or EPS.
Colours: the Figure 2 palette and pale tints (weak-label red, assay-aware light pink, policy blue);
PANEL and application circles use decorative pale-blue fills.
Model icon: Morgan-bit vector -> class-balanced logistic score p -> own margin |p - 0.5| (Methods, 'EffluxScope model').
Layout checks (run in __main__): text overlaps, type size floor, wheel labels inside their wedge and clear of the hub,
all artists inside the canvas margin; checked at 100, 300 and 600 dpi (hinted and unhinted text metrics).
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle, Polygon, Rectangle, FancyArrowPatch, FancyBboxPatch

TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS))
import palette
from palette import C
import make_toc_graphic as v1                      # shared constants, metrics and exporters

W, H = v1.W, v1.H
INK, MID, LIGHT = "#000000", "#000000", palette.C["grey"]
FS_TITLE, FS_HEAD, FS_BODY, FS_SUB = 7.0, 6.8, 6.3, 6.1
V01, V03, FULL = v1.V01, v1.V03, v1.FULL
STYLE = palette.FIGURE_STYLE
PANEL = STYLE["panel"]                                   # decorative tint of two corner panels


def hexagon(ax, x, y, r, fc="white", ec=INK, lw=0.8, z=6, pointy_top=False):
    ang = np.deg2rad(60 * np.arange(6) + (30 if pointy_top else 0))
    ax.add_patch(Polygon(np.c_[x + r * np.cos(ang), y + r * np.sin(ang)], closed=True, fc=fc, ec=ec, lw=lw,
                         joinstyle="round", zorder=z))


def pol(cx, cy, r, deg):
    return cx + r * np.cos(np.deg2rad(deg)), cy + r * np.sin(np.deg2rad(deg))


def molecule(ax, smiles, box, lw=0.75):
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    m = Chem.MolFromSmiles(smiles); rdDepictor.SetPreferCoordGen(True); rdDepictor.Compute2DCoords(m)
    Chem.Kekulize(m, clearAromaticFlags=True)
    xy = np.array([list(m.GetConformer().GetAtomPosition(i))[:2] for i in range(m.GetNumAtoms())])
    x0, y0, x1, y1 = box
    sx = (x1 - x0) / (np.ptp(xy[:, 0]) or 1); sy = (y1 - y0) / (np.ptp(xy[:, 1]) or 1); s = min(sx, sy)
    xy = (xy - xy.mean(axis=0)) * s + np.array([(x0 + x1) / 2, (y0 + y1) / 2])
    ring_c = {}
    ri = m.GetRingInfo()
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx(); p, q = xy[i], xy[j]
        ax.plot([p[0], q[0]], [p[1], q[1]], color=INK, lw=lw, solid_capstyle="round", zorder=3)
        if b.GetBondType() == Chem.BondType.DOUBLE:
            d = q - p; n = np.array([-d[1], d[0]]) / (np.linalg.norm(d) or 1)
            cen = None
            for ring in ri.AtomRings():
                if i in ring and j in ring: cen = xy[list(ring)].mean(axis=0)
            off = 0.16 * np.linalg.norm(d)
            if cen is not None and np.dot(n, cen - (p + q) / 2) < 0: n = -n
            if cen is None: n = n * 1.0
            a2 = p + 0.18 * d + n * off; b2 = q - 0.18 * d + n * off
            ax.plot([a2[0], b2[0]], [a2[1], b2[1]], color=INK, lw=lw, solid_capstyle="round", zorder=3)
    col = {"N": C["blue"], "O": C["red"], "S": C["amber"]}
    for a in m.GetAtoms():
        c = col.get(a.GetSymbol())
        if c: ax.add_patch(Circle(xy[a.GetIdx()], 0.020, fc=c, ec="white", lw=0.4, zorder=4))
    return xy


def build():
    assert V01 < FULL < V03, f"rank order behind the hexagons changed: {V01:.3f}, {FULL:.3f}, {V03:.3f}"
    blue, purple, green, amber, red, grey = C["blue"], C["purple"], C["green"], C["amber"], C["red"], C["grey"]
    assay_c = palette.tint("red", 0.45)                  # assay-aware score: lighter member of the red (weak-label) family
    palette.assert_palette([blue, purple, green, amber, red, grey, assay_c])
    fig = plt.figure(figsize=(W, H), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
    cx, cy, R, rh = W / 2, 0.965, 0.675, 0.17
    gy0 = 0.225

    # ---- four corner panels -------------------------------------------------------------------
    for (x0, x1, y0, y1, fc) in ((0, W / 2, cy, H, PANEL), (W / 2, W, cy, H, "white"),
                                 (0, W / 2, gy0, cy, "white"), (W / 2, W, gy0, cy, PANEL)):
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc=fc, ec="none", zorder=0))
    ttl = dict(fontsize=FS_TITLE, fontweight="bold", color="#000000", va="center")
    ax.text(0.07, 1.68, "Input / Sources", ha="left", **ttl)
    ax.text(W - 0.07, 1.68, "Models / Framework", ha="right", **ttl)
    ax.text(0.07, 0.31, "Applications", ha="left", **ttl)
    ax.text(W - 0.07, 0.31, "Outputs", ha="right", **ttl)

    # ---- wheel: three ways of choosing which predictions to keep, compared at an equal retained count --
    for name, t1, t2, col in (("policy", 90, 210, STYLE["blue_fill"]), ("assay", 210, 330, STYLE["aa_fill"]), ("weak", 330, 450, STYLE["red_fill"])):
        ax.add_patch(Wedge((cx, cy), R, t1, t2, fc=col, ec="white", lw=1.5, zorder=2))
    ax.add_patch(Circle((cx, cy), R, fc="none", ec="white", lw=2.2, zorder=3))
    ax.add_patch(Circle((cx, cy), rh, fc="white", ec=STYLE["border"], lw=1.0, zorder=5))
    ax.text(cx, cy, "Equal\ncount", ha="center", va="center", fontsize=FS_BODY, fontweight="bold", color=INK,
            linespacing=1.0, zorder=6)
    lab = dict(ha="center", va="center", fontsize=FS_BODY, fontweight="bold", linespacing=1.08, zorder=7)
    ax.text(cx - 0.335, cy + 0.235, "Weak-label\n+ abstain", color=INK, **lab)
    ax.text(cx + 0.335, cy + 0.235, "Weak-label\nscore", color=INK, **lab)
    ax.text(cx, cy - 0.40, "Assay-aware\nscore", color=INK, **lab)

    # ---- Input / Sources ---------------------------------------------------------------------------
    molecule(ax, "COc1ccc2nccc(C(O)C3CC4CCN3CC4C=C)c2c1", (0.13, 1.14, 0.87, 1.56))
    ax.text(0.50, 1.045, "2,624 molecules", ha="center", va="center", fontsize=FS_HEAD, fontweight="bold", color=INK)

    # ---- Models / Framework: fingerprint bits -> class-balanced logistic score p -> own margin |p - 0.5| ----------
    nx = 2.88
    for i, bit in enumerate((1, 0, 1, 1, 0, 0, 1, 0, 1)):                      # a fingerprint bit vector
        ax.add_patch(Rectangle((nx - 0.227 + i * 0.052, 1.525), 0.038, 0.038, fc=INK if bit else "white", ec=INK, lw=0.5,
                               zorder=3))
    ax.add_patch(FancyArrowPatch((nx, 1.506), (nx, 1.466), arrowstyle="-|>", color=MID, lw=0.8, mutation_scale=5,
                                 shrinkA=0, shrinkB=0, zorder=3))
    px0, px1, py0, py1, pmid = 2.70, 3.06, 1.20, 1.45, 1.325                    # logistic curve with the 0.5 line
    ax.plot([px0, px0, px1], [py1, py0, py0], color=LIGHT, lw=0.7, zorder=2, solid_joinstyle="miter")
    ax.plot([px0, px1], [pmid, pmid], color=LIGHT, lw=0.5, ls=(0, (2, 1.5)), zorder=2)
    sig = lambda x: py0 + (py1 - py0) / (1 + np.exp(-(x - nx) * 18))
    xx = np.linspace(px0 + 0.01, px1, 160); ax.plot(xx, sig(xx), color=INK, lw=1.1, zorder=3)
    xm = 2.985; ym = float(sig(xm))                                             # margin = distance of p from 0.5
    ax.plot([xm, xm], [pmid, ym], color=INK, lw=0.9, zorder=4)
    for yy in (pmid, ym): ax.plot([xm - 0.02, xm + 0.02], [yy, yy], color=INK, lw=0.9, zorder=4)
    ax.add_patch(Circle((xm, ym), 0.016, fc=INK, ec="white", lw=0.4, zorder=5))
    ax.text(nx, 1.125, "|s \u2212 0.5|", ha="center", va="center", fontsize=FS_HEAD, fontweight="bold", color=INK)
    ax.text(nx, 1.035, "own margin", ha="center", va="center", fontsize=FS_SUB, color=MID)

    # ---- Applications ------------------------------------------------------------------------------
    for (x, y, t) in ((0.245, 0.68, "Efflux\nrisk"), (0.705, 0.68, "ADC\npayloads")):
        ax.add_patch(Circle((x, y), 0.21, fc=STYLE["panel"], ec=STYLE["border"], lw=1.5, zorder=3))
        ax.text(x, y, t, ha="center", va="center", fontsize=FS_SUB, fontstyle="italic", fontweight="bold",
                color=INK, linespacing=1.05, zorder=4)

    # ---- Outputs: balanced-accuracy rank at the equal retained count (qualitative, no numbers) ----
    ax.text(2.785, 0.745, "Balanced accuracy", ha="center", va="center", fontsize=FS_SUB, fontweight="bold", color=INK)
    for x, col, word in ((2.53, red, "low"), (2.785, blue, "mid"), (3.04, assay_c, "high")):
        hexagon(ax, x, 0.52, 0.125, fc=STYLE["panel"] if col == blue else STYLE["panel_pink"], ec=col, lw=1.5, z=4, pointy_top=True)
        ax.text(x, 0.52, word, ha="center", va="center", fontsize=FS_SUB, fontweight="bold", color=INK, zorder=5)

    # ---- thick curved arrows ---------------------------------------------------------------------
    ar = dict(arrowstyle="-|>", color="black", lw=1.7, mutation_scale=9, shrinkA=0, shrinkB=0, zorder=8)
    ax.add_patch(FancyArrowPatch((1.10, 1.665), (2.18, 1.665), connectionstyle="arc3,rad=-0.08", **ar))
    ax.add_patch(FancyArrowPatch((2.43, 1.20), (2.43, 0.84), connectionstyle="arc3,rad=-0.55", **ar))
    ax.add_patch(FancyArrowPatch((2.15, 0.292), (1.10, 0.292), connectionstyle="arc3,rad=-0.07", **ar))

    # Replace the redundant take-home footer with a centered black subject title.
    # Retain the original physical artwork scale: crop the now-empty lower band,
    # translate the existing drawing downward and reserve a white title strip.
    ax.set_position([0, -0.15 / H, 1, 1])
    fig.text(0.5, 0.971, "ABCB1 efflux prediction at matched coverage",
             ha="center", va="top", fontsize=7.0, fontweight="bold", color="#000000")
    return fig


def layout_report(fig, dpi=600):
    """Check (at the export resolution, unhinted) that wheel labels sit inside their wedge and clear of the hub."""
    fig.set_dpi(dpi); fig.canvas.draw(); r = fig.canvas.get_renderer(); ax = fig.axes[0]; inv = ax.transData.inverted()
    wedges = [p for p in ax.patches if isinstance(p, Wedge)]
    hub = [p for p in ax.patches if isinstance(p, Circle) and abs(p.get_radius() - 0.17) < 1e-6][0]
    out = []
    for t in ax.texts:
        s_ = t.get_text()
        if s_ in ("Weak-label\n+ abstain", "Weak-label\nscore", "Assay-aware\nscore"):
            bb = t.get_window_extent(r); (x0, y0), (x1, y1) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
            cs = [(x0, y0), (x0, y1), (x1, y0), (x1, y1)]
            ok = any(all(w.get_path().transformed(w.get_patch_transform()).contains_point(c, radius=-0.012) for c in cs)
                     for w in wedges)
            d_hub = min(np.hypot(c[0] - hub.center[0], c[1] - hub.center[1]) for c in cs) - hub.get_radius()
            d_rim = max(np.hypot(c[0] - hub.center[0], c[1] - hub.center[1]) for c in cs)
            out.append((s_.replace("\n", " "), bool(ok), round(float(d_hub), 3), round(float(d_rim), 3)))
    return out


def bounds_report(fig, margin=0.02, dpi=600):
    """Artists whose extent leaves the canvas margin (inches)."""
    fig.set_dpi(dpi); fig.canvas.draw(); r = fig.canvas.get_renderer(); ax = fig.axes[0]; inv = ax.transData.inverted(); out = []
    items = [(t.get_text().replace("\n", " "), t.get_window_extent(r)) for t in ax.texts]
    items += [(type(p).__name__, p.get_window_extent(r)) for p in ax.patches if not isinstance(p, Rectangle)]
    for name, bb in items:
        (x0, y0), (x1, y1) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
        if x0 < margin or y0 < margin or x1 > W - margin or y1 > H - margin:
            out.append((name, round(float(x0), 3), round(float(y0), 3), round(float(x1), 3), round(float(y1), 3)))
    return out


if __name__ == "__main__":
    f = build()
    for dpi in (100, 300, 600):
        f.set_dpi(dpi); bad, small = v1.text_overlaps(f)
        print(f"{dpi:>3} dpi | text overlaps: {bad} | text below 6 pt: {small} | wedge/hub: "
              f"{[(a, b) for a, b, c, d in layout_report(f, dpi)]} | outside margin: {bounds_report(f, dpi=dpi)}")
    f.set_dpi(100)
    print("TIF px:", v1.save_all(f, "TOC_graphic_wheel"))
    print("BAcc ranks behind the hexagons (not printed on the graphic):", round(V01, 3), round(FULL, 3), round(V03, 3))
