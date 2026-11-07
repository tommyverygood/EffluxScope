"""Redraw current manuscript Figures 4, 5 and 6 from bundled frozen inputs.

No model is refitted. Figure 6 checks aggregate confusion counts; statistical
interval recomputation is a separate supported entrypoint in the full release.
Run with Python 3.12, NumPy, pandas, Matplotlib; Ghostscript supplies EPS export.
"""
from pathlib import Path
import subprocess
import sys
from matplotlib.font_manager import findfont
try:
    findfont('Arial', fallback_to_default=False)
except ValueError as exc:
    raise SystemExit('Exact approved-figure rendering requires a separately installed Arial font; font files are not redistributed.') from exc
ROOT = Path(__file__).resolve().parent
TOOLS = ROOT / '论文主稿/latex_jcim/tools'
(ROOT / '论文主稿/latex_jcim/figures').mkdir(parents=True, exist_ok=True)
for name in ['make_figure3_jcim.py', 'make_figure2_jcim.py', 'make_figure4_jcim.py']:
    subprocess.run([sys.executable, str(TOOLS / name)], cwd=ROOT, check=True)

subprocess.run([sys.executable, str(ROOT / "toc/tools/make_toc_graphic_wheel.py")], cwd=ROOT, check=True)
