from pathlib import Path
import json,re,hashlib
BASE=Path(__file__).resolve().parent
entries=[]
def esc(s):return s.replace('–','--').replace('—','---').replace('&',r'\&')
def entry(kind,key,fields):
 entries.append('@'+kind+'{'+key+',\n'+',\n'.join('  '+k+' = {'+esc(str(v))+'}' for k,v in fields.items())+'\n}\n')
entry('techreport','ema2017besponsa',{
'author':'{European Medicines Agency}',
 'title':'{BESPONSA} (inotuzumab ozogamicin): Assessment report',
 'institution':'European Medicines Agency, Committee for Medicinal Products for Human Use',
 'year':'2017','number':'EMA/289046/2017',
 'note':'21 April 2017; Procedure No. EMEA/H/C/004119/0000',
 'url':'https://www.ema.europa.eu/en/documents/assessment-report/besponsa-epar-public-assessment-report_en.pdf'})
for key,file,year,note in [
('obi2025payload','OBI992_Crossref_metadata.json','2025','First published online 30 December 2024; 2025 issue. Supplementary Methods and Tables, Table S5.'),
('weng2023exatecan','Weng2023_Crossref_metadata.json','2023','Supporting Information: Figures S1--S19, DOI: 10.1158/2159-8290.22553455.v1; Figure S1E.')]:
 d=json.loads((BASE/'source_extracts'/file).read_text())['message']
 fields={'author':' and '.join(a['family']+', '+a['given'] for a in d['author']),'title':d['title'][0],'journal':d['container-title'][0],'year':year,'volume':d['volume'],'number':d['issue'],'pages':d['page'].replace('-','--'),'doi':d['DOI'],'url':'https://doi.org/'+d['DOI'],'note':note}
 fields['title']=fields['title'].replace('OBI-992','{OBI-992}').replace('TROP2','{TROP2}')
 entry('article',key,fields)
entry('techreport','pmda2021padcev',{
'author':'{Astellas Pharma Inc.}',
 'title':'Enfortumab vedotin: Pharmacokinetics written summary and tabulated summaries',
 'institution':'Pharmaceuticals and Medical Devices Agency',
 'year':'2021',
 'note':'Japanese public CTD, Sections 2.6.4--2.6.5; Study 8004, Section 2.6.5.15.1, printed p. 39 (PDF p. 76)',
 'url':'https://www.pmda.go.jp/drugs/2021/P20211011002/800126000_30300AMX00454000_I100_1.pdf'})
entry('techreport','fda2020blenrep',{
'author':'{U.S. Food and Drug Administration}',
 'title':'{BLA} Multi-disciplinary Review and Evaluation: Belantamab mafodotin',
 'institution':'Center for Drug Evaluation and Research, U.S. Food and Drug Administration',
 'year':'2020','number':'761158Orig1s000',
 'note':'Reference ID: 4651784; Section 6.3.2.4, printed p. 98 (PDF p. 99)',
 'url':'https://www.accessdata.fda.gov/drugsatfda_docs/nda/2020/761158Orig1s000MultidisciplineR.pdf'})
entry('inproceedings','bao2025transporters',{
'author':'Bao, Yifan and Park, Seong Hee and Shi, Yifan and Yang, Bin and Evers, Raymond',
 'title':'Transporter Evaluation of Payloads in Antibody-Drug Conjugates for Cancer Treatment',
 'booktitle':'Abstracts from the 14th International ISSX Meeting',
 'year':'2025','address':'Chicago, Illinois, USA','pages':'140',
 'note':'Abstract P178; 21--24 September 2025. Conference abstract; no quantitative efflux ratios reported.',
 'url':'https://cdn.ymaws.com/www.issx.org/resource/resmgr/meetings_-_general/2025/ISSX_2025_Complete_Abstract_.pdf'})
entry('article','maji2025auristatin',{
'author':'Maji, Dolonchampa and Wolke, Malerie and Khaja, Shamim and Savaryn, John P. and Kalvass, John C. and Jenkins, Gary J.',
 'title':'Application of drug-induced growth rate inhibition and intracellular drug exposures for comprehensive evaluation of cellular drug sensitivity',
 'journal':'Scientific Reports','year':'2025','volume':'15','number':'1','pages':'5064',
 'doi':'10.1038/s41598-025-86919-7','url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC11814085/',
 'note':'Table 1; article number 5064'})
(BASE/'new_references.bib').write_text('\n'.join(entries))
keys=re.findall(r'@\w+\{([^,]+),', '\n'.join(entries));assert len(keys)==7 and len(set(keys))==7
sources={
'ema2017besponsa':{'basis':'Official EMA original PDF cover','metadata_file':'../family_expansion/sources/EMA_Besponsa_assessment.txt','checked_fields':['title','date','number','procedure','agency']},
'obi2025payload':{'basis':'Publisher-deposited Crossref metadata plus official publisher indexed title/PMC supplement page','metadata_file':'source_extracts/OBI992_Crossref_metadata.json','note':'Formal issue year2025 differs from online2024. All13 authors from metadata.'},
'weng2023exatecan':{'basis':'Publisher-deposited Crossref metadata plus AACR Figshare SI','metadata_file':'source_extracts/Weng2023_Crossref_metadata.json','SI_metadata':'../family_expansion/sources/Weng2023_figshare_metadata.json','note':'All23 authors from metadata; issue13(4)950-973; SI figureS1-S19 deposit DOI10.1158/2159-8290.22553455.v1.'},
'pmda2021padcev':{'basis':'Original Japanese public sponsor CTD PDF cover and p76; English title descriptive translation','metadata_file':'source_extracts/PMDA_Padcev_current.pdf','note':'Astellas sponsor author printed on CTD; PMDA public host; year2021 from official release URL.'},
'fda2020blenrep':{'basis':'FDA original review cover and printed98','metadata_file':'source_extracts/FDA_Blenrep2020_review.pdf','note':'Institutional author; CDER/application/referenceID verified.'},
'bao2025transporters':{'basis':'Original meeting abstract PDF cover and p140','metadata_file':'source_extracts/ISSX2025_abstracts.pdf','note':'All5 authors; conference place/date verified cover; abstractP178.'},
'maji2025auristatin':{'basis':'Original PMC HTML header and PubMed PMID39934176','metadata_file':'source_extracts/Maji2025_PMC.html','note':'Crossref request429; did not infer fields from failed request.'}
}
(BASE/'new_references_verification.json').write_text(json.dumps({'date_client':'2026-10-12','keys':keys,'sources':sources,'bib_sha256':hashlib.sha256((BASE/'new_references.bib').read_bytes()).hexdigest(),'manuscript_modified':False},ensure_ascii=False,indent=2))
print('PASS:',', '.join(keys))
