"""Extract published Yusoh 2025 summaries; do not infer raw dose responses.

S2 is a raster image, initially OCRed by installed tesseract. Fixed row order
and three explicitly documented corrections were verified against the image.
No Bliss/ZIP synergy or independent statistical tests are calculated.
"""
from pathlib import Path
import csv, json, re, unicodedata, hashlib
import pandas as pd
from docx import Document

BASE = Path(__file__).parent
LIBRARY = BASE.parent.parent / 'upload/L1021-DiscoveryProbe-FDA-approved-Drug-Library (1).xlsx'
DOCX = BASE / 'MOL2-19-3387-s001.docx'
SOURCE = 'https://pmc-oa-opendata.s3.amazonaws.com/PMC12591312.1/MOL2-19-3387-s001.docx'
ARTICLE = 'https://febs.onlinelibrary.wiley.com/doi/10.1002/1878-0261.70093'

# The original OCR intermediates were not retained in the delivery archive.
# Do not silently substitute derived CI tables for the source extraction.
required = [LIBRARY, DOCX, BASE / 'PMC12591312.1.json']
required.extend(BASE / f'S2_block{block}.txt' for block in range(1, 5))
missing = [str(path) for path in required if not path.is_file()]
if missing:
    raise FileNotFoundError(
        'Missing published-supplement extraction inputs:\n  '
        + '\n  '.join(missing)
        + '\nThe delivery archive retains the extracted result tables but '
        'omits the supplement, PMC metadata, and OCR text intermediates. '
        'See docs/DRUG_REPRODUCIBILITY.md. Recreating the OCR inputs '
        'requires checking their row order and values against the source image.'
    )

def name_key(value):
    """Normalize typography/parenthetical synonyms only, retaining salts."""
    s = unicodedata.normalize('NFKC', str(value)).casefold()
    s = re.sub(r'\([^)]*\)', '', s)
    s = s.replace('‐', '-').replace('–', '-')
    return re.sub(r'\s+', ' ', s).strip()

def value_parts(raw):
    s = re.sub(r'\s+', ' ', raw.replace('µ', 'u')).strip()
    # S5 has missing plus-minus separators in its original cells; leave raw.
    if s == 'ND': return {'raw': raw, 'qualifier':'not_determined', 'mean':None, 'sd':None}
    if s.startswith('>'):
        return {'raw':raw,'qualifier':'greater_than','mean':float(s[1:].strip()),'sd':None}
    if '±' in s:
        a,b=s.split('±');return {'raw':raw,'qualifier':'equal','mean':float(a.strip()),'sd':float(b.strip())}
    if re.fullmatch(r'\d+(?:\.\d+)?',s):
        return {'raw':raw,'qualifier':'equal','mean':float(s),'sd':None}
    return {'raw':raw,'qualifier':'ambiguous_original_separator','mean':None,'sd':None}

def write_csv(name, rows):
    if not rows: return
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with (BASE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

d=Document(DOCX)
s1=[]
for row in d.tables[0].rows[1:]:
    plate,well,nsc,cas,name,mw=[c.text.strip() for c in row.cells]
    s1.append(dict(screen_id=f'{plate}/{well[0]}{int(well[1:])}',plate=plate,well=well,nsc=nsc,cas=cas,drug=name,mw_raw=mw))
assert len(s1)==166 and len({r['screen_id'] for r in s1})==166

# OCR parsed numeric columns only. The source TIFF and cropped OCR outputs
# remain alongside the outputs for audit. Labels are constructed by the
# source's explicit ordered NCI plate grid, then checked against S1.
manual_corrections={
    (1,8): [0.40,0.86,2.57],  # OCR's 2.0f; visually 2.57, 4893/A3.
    (2,18): [3.18,0.65,0.17], # OCR missed decimal points, 4893/C9.
    (4,38): [4.61,3.96,3.22], # OCR's 461; visually 4.61, 4894/G11.
}
ci=[]
for block in range(1,6):
    if block==5:
        numbers=[[6.43,3.55,1.46],[4.29,2.84,0.26],[1.11,0.95,3.64],
                 [2.55,0.71,1.05],[7.54,7.68,1.91],[1.48,1.48,0.21]]
        # Six rows visually transcribed; plain tesseract returned no text.
        plate='4895';start=2
    else:
        lines=[s for s in (BASE/f'S2_block{block}.txt').read_text().splitlines() if s.strip()]
        assert len(lines)==40,(block,len(lines))
        numbers=[]
        for index,line in enumerate(lines):
            if (block,index) in manual_corrections:
                numbers.append(manual_corrections[block,index]);continue
            tokens=re.findall(r'(?<![\w])\d+[.,]\d+(?![\w])',line)
            assert len(tokens)>=3,(block,index,line,tokens)
            numbers.append([float(x.replace(',','.')) for x in tokens[-3:]])
        plate='4893' if block<=2 else '4894';start=2 if block in (1,3) else 7
    for index,values in enumerate(numbers):
        key=f'{plate}/{chr(65+index%8)}{start+index//8}'
        ci.append(dict(screen_id=key,ci_0_1uM=values[0],ci_1uM=values[1],ci_10uM=values[2],
                       n_ci_lt_0_9=sum(v<0.9 for v in values),n_ci_le_0_9=sum(v<=0.9 for v in values),
                       ci_min=min(values),ci_extraction='OCR_plus_visual_check' if block<5 else 'visual_transcription'))
assert len(ci)==166 and {r['screen_id'] for r in ci}=={r['screen_id'] for r in s1}
s1_by={r['screen_id']:r for r in s1}
for r in ci:r.update(drug=s1_by[r['screen_id']]['drug'],cas=s1_by[r['screen_id']]['cas'])

s3=[]
for row in d.tables[1].rows[2:]:
    ident,drug,mono,combo,fold=[c.text.strip() for c in row.cells]
    # S3 uses plate serials with an added '/49'; S1 and S2 do not.
    key=ident.replace('/49_','/')
    entry=dict(screen_id=key,drug=drug,cell_line='MDA-MB-231',treatment_h=24,
               olaparib_uM=10,independent_experiments=2,
               mono_ic50_raw=mono,combo_ic50_raw=combo,fold_change_raw=fold)
    for prefix,value in [('mono_ic50',mono),('combo_ic50',combo),('fold_change',fold)]:
        for field,result in value_parts(value).items():
            if field!='raw':entry[f'{prefix}_{field}']=result
    assert key in s1_by,(key,drug)
    entry.update(cas=s1_by[key]['cas'])
    entry.update({k:v for k,v in next(r for r in ci if r['screen_id']==key).items() if k.startswith('ci_') or k.startswith('n_ci')})
    s3.append(entry)

s5=[]
for row in d.tables[3].rows[2:]:
    model,mono,combo,fold=[c.text.strip() for c in row.cells]
    r=dict(cell_line=model,drug='Exemestane',olaparib_uM=10,treatment_h=24,
           mono_ic50_raw=mono,combo_ic50_raw=combo,fold_change_raw=fold,independent_experiments=2)
    for prefix,value in [('mono_ic50',mono),('combo_ic50',combo),('fold_change',fold)]:
        for field,result in value_parts(value).items():
            if field!='raw':r[f'{prefix}_{field}']=result
    s5.append(r)

lib=pd.read_excel(LIBRARY,sheet_name='Chemical Data').fillna('')
lib_names={}
for index,row in lib.iterrows():lib_names.setdefault(name_key(row['Item Name']),[]).append(index)
# Manual equivalent identifiers verified by FDA GSRS / NCI, preserving CAS.
synonyms={'acalabrutinib':'acp-196','omacetaxine mepesuccinate':'homoharringtonine'}
synonym_sources={
    'acalabrutinib':'https://www.cancer.gov/publications/dictionaries/cancer-drug/def/acalabrutinib',
    'omacetaxine mepesuccinate':'https://precision.fda.gov/uniisearch/srs/unii/6FG8041S5B',
}
matches=[]
for r in s1:
    key=name_key(r['drug'])
    candidates={i:'exact_CAS' for i in lib.index[lib['CAS Number'].astype(str).str.strip()==r['cas']] if r['cas']}
    for i in lib_names.get(key,[]):candidates.setdefault(i,'normalized_name_retaining_salt')
    for i in lib_names.get(synonyms.get(key,'__none__'),[]):candidates.setdefault(i,'verified_synonym')
    for i,method in candidates.items():
        q=lib.loc[i]
        match={**r,'catalog_number':q['CatalogNumber'],'library_name':q['Item Name'],
               'library_cas':q['CAS Number'],'rack':q['Rack Number'],'well_library':q['Plate Location'],
               'match_method':method,'synonym_source':synonym_sources.get(key,'') if method=='verified_synonym' else '',
               'vendor_target_annotation':q['Target']}
        match['identity_status']='exact_CAS' if method=='exact_CAS' else 'provisional_name_match_CAS_conflict'
        match.update({k:v for k,v in next(z for z in ci if z['screen_id']==r['screen_id']).items() if k.startswith('ci_') or k.startswith('n_ci')})
        hit=next((z for z in s3 if z['screen_id']==r['screen_id']),None)
        match['author_selected_hit']=hit is not None
        if hit:match.update({k:v for k,v in hit.items() if k.startswith('mono_') or k.startswith('combo_') or k.startswith('fold_')})
        matches.append(match)

write_csv('S1_drug_identities.csv',s1)
write_csv('S2_CI_numeric.csv',ci)
write_csv('S3_IC50_numeric.csv',s3)
write_csv('S5_exemestane_models.csv',s5)
write_csv('L1021_screen_overlap.csv',matches)
write_csv('L1021_author_hit_overlap.csv',[r for r in matches if r['author_selected_hit']])
matched_ids={r['screen_id'] for r in matches}
exact_cas_ids={r['screen_id'] for r in matches if r['match_method']=='exact_CAS'}
verified_ci=[r for r in ci if r['screen_id'] in exact_cas_ids]
docx_md5=hashlib.md5(DOCX.read_bytes()).hexdigest()
pmc_metadata=json.loads((BASE/'PMC12591312.1.json').read_text())
expected_md5=next(u.split('md5=')[1] for u in pmc_metadata['media_urls'] if 's001.docx' in u)
assert docx_md5==expected_md5
summary=dict(source=SOURCE,article=ARTICLE,screen_drugs=166,ci_rows=166,ci_values=498,
             supplement_md5=docx_md5,supplement_md5_matches_PMC_metadata=True,
             library_catalog_rows=len(lib),matched_screen_compounds=len(matched_ids),
             exact_CAS_matched_screen_compounds=len(exact_cas_ids),
             additional_provisional_name_matched_compounds=len(matched_ids-exact_cas_ids),
             matched_library_catalog_rows=len({r['catalog_number'] for r in matches}),
             author_hits=10,author_hits_matched=len({r['screen_id'] for r in matches if r['author_selected_hit']}),
             author_hit_percent_of_166=100*10/166,ci_lt_0_9_all_3=sum(r['n_ci_lt_0_9']==3 for r in ci),
             ci_lt_0_9_at_least_2=sum(r['n_ci_lt_0_9']>=2 for r in ci),
             ci_le_0_9_at_least_2=sum(r['n_ci_le_0_9']>=2 for r in ci),
             confirmed_matches_ci_lt_0_9_at_least_2=sum(r['n_ci_lt_0_9']>=2 for r in verified_ci),
             confirmed_matches_ci_lt_0_9_all_3=sum(r['n_ci_lt_0_9']==3 for r in verified_ci),
             secondary_OCR_verification='Independent numeric-only OCR recovered 166 rows; differences were checked against original image and first-pass values retained. Three first-pass numeric corrections and six final-block rows manually visually verified.',
             model='parental MDA-MB-231; 24h; 10uM olaparib; not MDA-MB-468',
             missing_selected_hits=[r['drug'] for r in s3 if r['screen_id'] not in matched_ids],
             limitations=[
                 'Published Chou-Talalay CI summary; no raw viability matrix; no Bliss/ZIP recomputation.',
                 'Printed CI rounded to two decimal places; borderline 0.90 values cannot resolve underlying threshold.',
                 'Initial drug concentrations 0.1/1/10uM; selected-hit IC50 from additional dose ranges.',
                 'Author selected 10 hits; S2 alone cannot reconstruct the >50% inhibition condition.',
                 '166 library compounds versus reported 10=16.9% denominator discrepancy.',
                 'S3 IC50 >100 values are right censored; fold changes with > are lower bounds, not point estimates.',
                 'S5 original cells lack plus-minus separators for several mean/SD entries; retained raw, not parsed.',
                 'Vendor target annotations are not accepted as validated molecular mechanisms.',
                 'One additional screen compound (vemurafenib) matches name but not CAS and remains provisional.',
                 'The S2 CI table includes olaparib added to olaparib; interpretation of those scores needs caution.',
                 'CAS/name matches identify chemicals; uploaded list is not proof of compounds physically available.',
             ])
(BASE/'extraction_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
print('Author-hit overlap:')
for r in matches:
    if r['author_selected_hit']:
        print(r['drug'],r['catalog_number'],r['library_name'],r.get('mono_ic50_raw'),r.get('combo_ic50_raw'),r.get('fold_change_raw'),[r['ci_0_1uM'],r['ci_1uM'],r['ci_10uM']])
