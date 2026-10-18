from pathlib import Path
import urllib.request,json,re,concurrent.futures
from pypdf import PdfReader
BASE=Path(__file__).resolve().parent
endpoints={
'OBI992_PMC_metadata':'https://pmc.ncbi.nlm.nih.gov/articles/PMC11791482/',
'Weng2023_PubMed_metadata':'https://pubmed.ncbi.nlm.nih.gov/36693125/?format=pubmed',
'OBI992_Crossref_metadata':'https://api.crossref.org/works/10.1158/1535-7163.MCT-24-0588',
'Weng2023_Crossref_metadata':'https://api.crossref.org/works/10.1158/2159-8290.CD-22-1368',
'Maji2025_Crossref_metadata':'https://api.crossref.org/works/10.1038/s41598-025-86919-7'
}
def f(x):
 k,u=x
 try:
  r=urllib.request.urlopen(urllib.request.Request(u,headers={'User-Agent':'Mozilla/5.0 citation verification'}),timeout=25);d=r.read();p=BASE/'source_extracts'/(k+('.json' if 'Crossref' in k else '.html'));p.write_bytes(d)
  return {'id':k,'url':u,'status':r.status,'path':str(p)}
 except Exception as e:return {'id':k,'url':u,'error':repr(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:res=list(pool.map(f,endpoints.items()))
(BASE/'bib_metadata_fetch_receipts.json').write_text(json.dumps(res,ensure_ascii=False,indent=2))
for r in res:print(r['id'],r.get('status'),r.get('error',''))
for file in ['PMDA_Padcev_current.pdf','FDA_Blenrep2020_review.pdf','ISSX2025_abstracts.pdf']:
 p=PdfReader(BASE/'source_extracts'/file);print('\nCOVER',file,p.pages[0].extract_text()[:2500])
