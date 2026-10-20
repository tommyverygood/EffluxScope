#!/usr/bin/env python3
"""Render the original Figure 1 with stronger lines for the JCIM manuscript.

Run with an existing Python environment containing matplotlib and numpy.
The original generator and released input data are used without modification.
Only this manuscript's figure assets are written.
"""
from pathlib import Path
import runpy
import subprocess

from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.collections import Collection
from matplotlib.text import Annotation
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / "figures_src" / "make_figure1_v2.py"


def main():
    source = runpy.run_path(str(SOURCE))
    source["verify_palette"]()
    grid, _ = source["load"]()
    frozen = SOURCE.resolve().parents[2] / "数据/processed/abcb1_public_predictor_metrics_v1.csv"
    differences = source["verify_against_frozen"](grid, frozen)
    if differences:
        raise AssertionError(differences)
    source["palette"].apply_base_style({"savefig.bbox": None})
    figure = source["draw"]()
    # The unchanged canvas is included at approximately 0.956 scale. A 0.65-pt
    # source stroke leaves margin above the 0.5-pt final-size recommendation.
    for artist in figure.findobj():
        if isinstance(artist, Annotation) and artist.arrow_patch is not None:
            artist.arrow_patch.set_linewidth(max(0.65, artist.arrow_patch.get_linewidth()))
        if isinstance(artist, (Line2D, Patch)):
            width = artist.get_linewidth()
            if 0 < width < 0.65:
                artist.set_linewidth(0.65)
        elif isinstance(artist, Collection):
            widths = artist.get_linewidths()
            artist.set_linewidths([max(w, 0.65) if w > 0 else w for w in widths])
    stem = ROOT / "figures/Figure1_endpoint_shift"
    figure.savefig(stem.with_suffix(".pdf"), bbox_inches=None)
    figure.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches=None)
    plt.close(figure)
    subprocess.run([
        "gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=eps2write",
        f"-sOutputFile={stem.with_suffix('.eps')}", str(stem.with_suffix('.pdf')),
    ], check=True)
    print("Figure 1 rendered; frozen numerical metrics verified unchanged.")


if __name__ == "__main__":
    main()
