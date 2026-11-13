import csv
import hashlib
import io

def parse_samples(text):
    if len(text.encode()) > 1500000: raise ValueError('CSV 不得超过 1.5MB')
    reader=csv.DictReader(io.StringIO(text))
    fieldnames=reader.fieldnames or []
    if not {'plant_id','greenhouse','symptoms'}.issubset(fieldnames):
        raise ValueError('需要 plant_id,greenhouse,symptoms 列，可选 unexamined 列记录未查病征，症状以 | 分隔')
    track_missing='unexamined' in fieldnames
    rows=[]; ids=set(); universe=set()
    for raw in reader:
        row={k:(raw.get(k) or '').strip() for k in ('plant_id','greenhouse','symptoms')}
        if not row['plant_id'] or not row['greenhouse']: raise ValueError('植株和温室不能为空')
        if row['plant_id'] in ids: raise ValueError('植株编号不能重复')
        items=sorted(set(x.strip() for x in row['symptoms'].split('|') if x.strip()))
        if any(len(x)>30 for x in items): raise ValueError('症状名称最长30字')
        record={'plant_id':row['plant_id'],'greenhouse':row['greenhouse'],'symptoms':items}
        universe.update(items)
        if track_missing:
            missed=sorted(set(x.strip() for x in (raw.get('unexamined') or '').split('|') if x.strip()))
            if any(len(x)>30 for x in missed): raise ValueError('症状名称最长30字')
            both=sorted(set(items)&set(missed))
            if both: raise ValueError(f"病征不能既是阳性又是未查：{'、'.join(both)}")
            universe.update(missed)
            record['unexamined']=missed
        ids.add(row['plant_id'])
        rows.append(record)
        if len(rows)>5000: raise ValueError('最多5000株')
    if len(rows)<20: raise ValueError('至少20株样本')
    if not 2<=len(universe)<=12: raise ValueError('症状种类应为2至12种')
    return rows

def checksum(rows):
    import json
    return hashlib.sha256(json.dumps(rows,sort_keys=True,ensure_ascii=False).encode()).hexdigest()
