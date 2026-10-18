from pathlib import Path
import hashlib,json,concurrent.futures,urllib.request,html,re
from pypdf import PdfReader
BASE=Path(__file__).resolve().parent/'source_extracts'
sources={
'ISSX2025_abstracts':'https://cdn.ymaws.com/www.issx.org/resource/resmgr/meetings_-_general/2025/ISSX_2025_Complete_Abstract_.pdf',
'EMA_Elahere_assessment':'https://www.ema.europa.eu/en/documents/assessment-report/elahere-epar-public-assessment-report_en.pdf',
'EMA_Elahere_RMP':'https://www.ema.europa.eu/en/documents/rmp/elahere-epar-risk-management-plan_en.pdf',
'Maji2025_PMC':'https://pmc.ncbi.nlm.nih.gov/articles/PMC11814085/',
'AU2012244675A1':'https://patents.google.com/patent/AU2012244675A1/en',
'Blenrep2025_emc':'https://www.medicines.org.uk/emc/product/100782/smpc',
'Rock2015_publisher':'https://www.sciencedirect.com/science/article/pii/S0090955624112925',
'SLC46A3_2026_publisher':'https://pubs.acs.org/doi/10.1021/acs.molpharmaceut.6c00535',
'Tomabechi2022_PMC':'https://pmc.ncbi.nlm.nih.gov/articles/PMC9896951/',
'GSK117159_protocol':'https://cdn.clinicaltrials.gov/large-docs/87/NCT02064387/Prot_000.pdf',
}
def fetch(item):
 key,url=item;out={'id':key,'url':url}
 try:
  r=urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 research source verification'}),timeout=45);content=r.read();out.update(status=r.status,final_url=r.url,content_type=r.headers.get('Content-Type'),bytes=len(content),sha256=hashlib.sha256(content).hexdigest())
  ispdf=content.startswith(b'%PDF'); path=BASE/(key+('.pdf' if ispdf else '.html'));path.write_bytes(content);out['saved_path']=str(path)
  if ispdf:
   p=PdfReader(path);tx='\n'.join(f'\n===== PDF page {i+1} =====\n'+(page.extract_text() or '') for i,page in enumerate(p.pages));(BASE/(key+'.txt')).write_text(tx);out['pdf_pages']=len(p.pages)
  else:
   tx=html.unescape(re.sub(r'<[^>]+>', ' ', content.decode('utf8',errors='replace')));(BASE/(key+'.txt')).write_text(tx);out['looks_blocked']=any(x in tx[:2000].lower() for x in ['checking your browser','access denied','captcha','forbidden'])
 except Exception as e:out['error']=repr(e)
 return out
if __name__ == '__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:res=list(pool.map(fetch,sources.items()))
 (BASE/'fetch_receipts.json').write_text(json.dumps(res,ensure_ascii=False,indent=2))
 for r in res:print(r['id'],r.get('status'),r.get('pdf_pages'),r.get('looks_blocked'),r.get('error',''))
