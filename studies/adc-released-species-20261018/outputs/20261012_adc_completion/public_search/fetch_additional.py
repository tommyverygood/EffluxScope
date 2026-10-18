from fetch_sources import fetch,BASE
import json,concurrent.futures,urllib.request,urllib.parse
more={
'FDA_Blenrep2020_review':'https://www.accessdata.fda.gov/drugsatfda_docs/nda/2020/761158Orig1s000MultidisciplineR.pdf',
'PMDA_Padcev_current':'https://www.pmda.go.jp/drugs/2021/P20211011002/800126000_30300AMX00454000_I100_1.pdf',
'PMDA_Blenrep_current':'https://www.pmda.go.jp/files/000280862.pdf',
'Maji2025_methods_Mahon2017':'https://doi.org/10.1021/acs.jmedchem.7b00565',
}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:res=list(pool.map(fetch,more.items()))
(BASE/'fetch_additional_receipts.json').write_text(json.dumps(res,ensure_ascii=False,indent=2))
for r in res:print(r['id'],r.get('status'),r.get('pdf_pages'),r.get('error',''))
terms=['("MMAF" AND ("cys-mcMMAF" OR "cysteine")) AND ("P-glycoprotein" OR "MDR1" OR "ABCB1")','("DM1" AND "Lys-MCC-DM1") AND ("transport" OR "permeability")','("DM4" AND "S-methyl-DM4") AND ("MDR1" OR "P-glycoprotein")','"TransPorter evaluation of payloads"']
receipts=[]
for i,t in enumerate(terms):
 url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'+urllib.parse.urlencode({'db':'pubmed','term':t,'retmax':100,'retmode':'json'})
 try:
  r=urllib.request.urlopen(url,timeout=30);d=json.loads(r.read());(BASE/f'PubMed_query_{i+1}.json').write_text(json.dumps(d,indent=2));receipts.append({'query':t,'url':url,'count':d['esearchresult']['count'],'ids':d['esearchresult']['idlist']})
 except Exception as e:receipts.append({'query':t,'error':repr(e)})
(BASE/'PubMed_search_receipts.json').write_text(json.dumps(receipts,ensure_ascii=False,indent=2))
print(receipts)
