from pathlib import Path
import subprocess,sys
root=Path(__file__).resolve().parent
(root/'论文主稿/latex_jcim/figures').mkdir(parents=True,exist_ok=True)
subprocess.run([sys.executable,str(root/'论文主稿/latex_jcim/tools/make_figure1_jcim.py')],check=True,cwd=root)
