"""Readability adjustments for JCIM Figure 3 / 4 only.

The archived data and original plotting generators remain unchanged.  These
adjustments increase minimum font and stroke sizes before PDF/EPS export, use
full-size digits for logarithm bases, and reflow long explanatory text.
"""
import re
import textwrap
import shutil
import subprocess

from matplotlib.text import Text, Annotation
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.collections import Collection

MIN_FONT_PT = 7.2
MIN_STROKE_PT = 0.85


def apply_jcim_readability(fig):
    """Retain scientific values, labels and panel order; change typography."""
    for text in fig.findobj(match=Text):
        value = text.get_text()
        # Same logarithm bases, without a small separately sized math script.
        value = re.sub(r"\\log_(?:\{10\}|10)", r"\\mathrm{log10}", value)
        value = re.sub(r"\\log_(?:\{2\}|2)", r"\\mathrm{log2}", value)
        value = value.replace(r"\mathrm{log10} 2", r"\mathrm{log10}(2)")
        value = value.replace("log$_{10}$", "log10")
        if value.startswith("Results are reproduced") or value.startswith("Only the two endpoint-aligned"):
            value = "\n".join(textwrap.wrap(value.replace("\n", " "), width=112,
                break_long_words=False, break_on_hyphens=False))
        if value.startswith("criterion: point estimate"):
            value = value.replace(" and bootstrap", "\nand bootstrap")
            text.set_bbox(dict(facecolor="white", edgecolor="none", pad=0.3))
        if value.startswith("peak $"):
            value = value.replace("  aligned", "\naligned")
        text.set_text(value)
        text.set_fontsize(max(MIN_FONT_PT, text.get_fontsize()))
        text.set_fontfamily("DejaVu Sans")
    # Increase strokes that are actually drawn; retain deliberate fill-only
    # elements and their zero-width/non-existent outlines.
    for artist in fig.findobj():
        if isinstance(artist, Annotation) and artist.arrow_patch is not None:
            artist.arrow_patch.set_linewidth(MIN_STROKE_PT)
        if isinstance(artist, Line2D):
            if artist.get_linewidth() > 0:
                artist.set_linewidth(max(MIN_STROKE_PT, artist.get_linewidth()))
            if artist.get_markeredgewidth() > 0:
                artist.set_markeredgewidth(max(MIN_STROKE_PT, artist.get_markeredgewidth()))
        elif isinstance(artist, Patch):
            if artist.get_linewidth() > 0 and artist.get_edgecolor()[-1] > 0:
                artist.set_linewidth(max(MIN_STROKE_PT, artist.get_linewidth()))
        elif isinstance(artist, Collection):
            artist.set_linewidths([max(MIN_STROKE_PT, width) if width > 0 else width
                                  for width in artist.get_linewidths()])
    for ax in fig.axes:
        ax.tick_params(which="both", width=MIN_STROKE_PT)
    fig.canvas.draw()
    return fig


def export_eps_from_pdf(pdf_path):
    """Preserve transparency appearance when converting vector PDF to EPS."""
    ghostscript = shutil.which("gs")
    if not ghostscript:
        raise RuntimeError("Ghostscript is required for the matching EPS export")
    subprocess.run([ghostscript, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE",
        "-sDEVICE=eps2write", "-r600", f"-sOutputFile={pdf_path.with_suffix('.eps')}",
        str(pdf_path)], check=True)
