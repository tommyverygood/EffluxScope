"""Verify inputs, redraw seven graphics, and recompute the stated Figure 6 intervals."""
from pathlib import Path
import argparse,hashlib,importlib.metadata,json,os,platform,shutil,subprocess,sys
from PIL import Image,ImageChops
from pypdf import PdfReader
from matplotlib.font_manager import findfont
ROOT=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output-dir',type=Path,default=Path('verification'));args=ap.parse_args()
    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    for binary in ['gs','pdftoppm']:
        if shutil.which(binary) is None:raise RuntimeError(f'Required executable missing: {binary}')
    findfont('Arial',fallback_to_default=False)
    lines=(ROOT/'SHA256SUMS.txt').read_text().splitlines()
    for line in lines:
        h,n=line.split('  ',1);path=ROOT/n
        if not path.is_file() or sha(path)!=h:raise AssertionError(f'Frozen supplied-file checksum mismatch: {n}')
    env=dict(os.environ);env.pop('PYTHONPATH',None);env['PYTHONNOUSERSITE']='1';env['MPLBACKEND']='Agg';env['PYTHONDONTWRITEBYTECODE']='1'
    commands=[['components/overview/scripts/make_figure_background.py'],['components/overview/scripts/make_figure_overview.py'],['components/overview/scripts/verify_overview.py'],['components/endpoint/run.py'],['components/results/run_results_figures.py'],['components/uncertainty/recompute_intervals.py','--output',str(out/'figure6_intervals.json')]]
    for i,cmd in enumerate(commands,1):
        with (out/f'command_{i}.log').open('w') as log:subprocess.run([sys.executable,str(ROOT/cmd[0]),*cmd[1:]],cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    ci=json.loads((out/'figure6_intervals.json').read_text());plot=json.loads((ROOT/'components/results/handoff/fig4_data.json').read_text())
    n=0
    for key,data in ci['result']['b'].items():
        for field,value in data.items():assert value==plot['b'][key][field];n+=1
    for field,value in ci['result']['band'].items():assert value==plot['band'][field];n+=1
    assert n==100 and ci['status']=='PASS'
    checks=[]
    for i,m in enumerate(json.loads((ROOT/'figure_map.json').read_text()),1):
        reference=ROOT/m['reference'];actual=ROOT/m['generated'];d=out/f'figure_{i}';d.mkdir(exist_ok=True)
        for name,path in [('expected',reference),('actual',actual)]:subprocess.run(['pdftoppm','-r','150','-png','-singlefile',str(path),str(d/name)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        a=Image.open(d/'expected.png').convert('RGB');b=Image.open(d/'actual.png').convert('RGB');ra=PdfReader(reference);rb=PdfReader(actual)
        c={'reference':m['reference'],'generated':m['generated'],'reference_sha256':sha(reference),'generated_sha256':sha(actual),'page_geometry_equal':len(ra.pages)==len(rb.pages) and all(list(x.mediabox)==list(y.mediabox) for x,y in zip(ra.pages,rb.pages)),'pdf_text_equal':len(ra.pages)==len(rb.pages) and all(x.extract_text()==y.extract_text() for x,y in zip(ra.pages,rb.pages)),'pixels_equal_150dpi':a.size==b.size and ImageChops.difference(a,b).getbbox() is None}
        checks.append(c)
        assert all(c[k] for k in ['page_geometry_equal','pdf_text_equal','pixels_equal_150dpi']),c
    report={'status':'PASS','release':'jcim-figures-20261107','scope':'Frozen-input redraw of six main figures and TOC; Figure 6 interval recomputation from fixed scaffold sufficient statistics, not regeneration of service/model predictions.','python':platform.python_version(),'platform':platform.system()+' '+platform.machine(),'packages':{n:importlib.metadata.version(n) for n in ['numpy','pandas','matplotlib','Pillow','rdkit','scipy','joblib','reportlab','pypdf']},'external_tools':{'ghostscript':subprocess.run(['gs','--version'],text=True,capture_output=True,check=True).stdout.strip(),'poppler':subprocess.run(['pdftoppm','-v'],text=True,capture_output=True).stderr.splitlines()[0]},'verified_supplied_files':len(lines),'figures':checks,'figure6_intervals':{k:ci[k] for k in ['status','seed','replicates','clusters','numeric_fields_checked','max_absolute_error']},'interval_fields_match_plotted_data':n,'model_refitting':False,'new_biological_validation':False}
    (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'status':report['status'],'figures':len(checks),'exact_pixels_and_text':True,'recomputed_interval_fields':n,'max_interval_error':ci['max_absolute_error']},indent=2))
if __name__=='__main__':main()
