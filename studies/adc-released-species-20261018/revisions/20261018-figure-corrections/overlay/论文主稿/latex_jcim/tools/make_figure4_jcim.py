"""Generate the current Figure 4 using its single maintained source."""
from pathlib import Path
import runpy

for ancestor in Path(__file__).resolve().parents:
    generator = ancestor / "figures_src/make_figure4.py"
    if generator.is_file():
        runpy.run_path(str(generator), run_name="__main__")
        break
else:
    raise FileNotFoundError("Cannot find figures_src/make_figure4.py")
