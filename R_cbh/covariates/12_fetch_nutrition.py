#!/usr/bin/env python3
"""Cache official regional WHO-standard nutrition indicators; no microdata read/upload."""
import csv, hashlib, json, urllib.parse, urllib.request
from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
out = Path('data/derived_cbh/regional_adjustment/planned17_audit')
out.mkdir(parents=True, exist_ok=True)
with open('data/derived_dhs/survey_registry.csv') as f:
    surveys = sorted({r['SurveyId'] for r in csv.DictReader(f) if r['SurveyId'] and r['SurveyId'] != 'NA'})
ids = 'CN_NUTS_C_HA2,CN_NUTS_C_WH2'
def fetch(endpoint, params, name):
    page, rows, manifest = 1, [], []
    while True:
        url = 'https://api.dhsprogram.com/rest/dhs/' + endpoint + '?' + urllib.parse.urlencode(dict(f='json',perpage=1000,page=page,**params))
        path = out / f'{name}_page{page}.json'
        if not path.exists():
            response = urllib.request.urlopen(url, timeout=60)
            assert urllib.parse.urlparse(response.geturl()).hostname == 'api.dhsprogram.com'
            payload = response.read()
            json.loads(payload)
            path.write_bytes(payload)
        x = json.loads(path.read_bytes())
        assert x['RecordsReturned'] == len(x['Data'])
        rows.extend(x['Data'])
        manifest.append(dict(url=url,file=str(path),md5=hashlib.md5(path.read_bytes()).hexdigest(),retrieved_utc=datetime.fromtimestamp(path.stat().st_mtime,timezone.utc).isoformat()))
        if page >= x['TotalPages']: break
        page += 1
    assert len(rows) == x['RecordCount']
    return rows, manifest
metadata, provenance = fetch('indicators',dict(indicatorIds=ids),'nutrition_metadata')
assert len(metadata)==2 and all('WHO standard' in x['Definition'] for x in metadata)
batches = [surveys[i:i+20] for i in range(0,len(surveys),20)]
def batch(item):
    i, ss=item
    return fetch('data',dict(indicatorIds=ids,surveyIds=','.join(ss),breakdown='subnational'),f'nutrition_{i}')
rows=[]
with ThreadPoolExecutor(max_workers=3) as pool:
    for records, manifest in pool.map(batch,enumerate(batches,1)):
        rows.extend(records); provenance.extend(manifest)
assert rows and all(x['IndicatorId'] in ids.split(',') and x['SurveyId'] in surveys for x in rows)
for filename, values in [('nutrition_published.csv',rows),('nutrition_source_manifest.csv',provenance)]:
    fields = list(dict.fromkeys(k for row in values for k in row))
    with (out/filename).open('w',newline='') as f:
        writer=csv.DictWriter(f,fields,lineterminator='\n');writer.writeheader();writer.writerows(values)
print(f'Cached {len(rows)} official regional nutrition observations for {len({x["SurveyId"] for x in rows})} surveys.')
