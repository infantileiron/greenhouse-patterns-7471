"""Apriori + blocked permutation, with BH over the prespecified full rule family.

每株每个病征为阳性/阴性/未查三态。每条规则只在相关病征全部查过的有效植株上
计算支持度、置信度、提升度；置换同样限于这些植株并保留温室分层。有效样本不足
MIN_EFFECTIVE 的规则不进入检验族，连同排除原因一并返回供核对。
"""
from collections import defaultdict
from itertools import combinations
import random
from src.utils.csv_data import checksum

MIN_EFFECTIVE=20  # 单条规则的最少有效植株数，与调查批次下限一致

def bh_adjust(pvalues):
    m=len(pvalues); result=[1.0]*m; running=1.0
    for rank,index in reversed(list(enumerate(sorted(range(m),key=lambda i:pvalues[i]),1))):
        running=min(running,pvalues[index]*m/rank)
        result[index]=min(1.0,running)
    return result

def apriori(items, pos_masks, unx_masks, full, min_support, max_length=3):
    """频繁项集按逐项集完整观测分母判定：只数相关病征都查过的植株。"""
    if not 0<min_support<=1: raise ValueError('支持度须在(0,1]')
    if full.bit_count()<20: raise ValueError('至少20株样本')
    def eligible(itemset):
        mask=0
        for x in itemset: mask|=unx_masks[x]
        return full&~mask
    def count(itemset, within):
        mask=within
        for x in itemset: mask&=pos_masks[x]
        return mask.bit_count()
    frequent={}
    previous={frozenset([x]) for x in items if (lambda e:e.bit_count()>0 and count((x,),e)/e.bit_count()>=min_support)(eligible((x,)))}
    frequent.update({s:count(s,eligible(s)) for s in previous})
    for k in range(2,max_length+1):
        candidates={a|b for a in previous for b in previous if len(a|b)==k}
        candidates={s for s in candidates if all(frozenset(z) in previous for z in combinations(s,k-1))}
        previous=set()
        for s in candidates:
            e=eligible(s); n=e.bit_count()
            if n and count(s,e)/n>=min_support: previous.add(s)
        frequent.update({s:count(s,eligible(s)) for s in previous})
        if not previous: break
    return frequent

def analyze(rows, params):
    support=float(params['min_support']); confidence=float(params['min_confidence'])
    permutations=int(params['permutations']); seed=int(params['seed'])
    if not .02<=support<=.8 or not .1<=confidence<=1: raise ValueError('支持度0.02–0.8，置信度0.1–1')
    if not 99<=permutations<=1999: raise ValueError('置换次数99–1999')
    n=len(rows)
    if not 20<=n<=5000: raise ValueError('样本量20–5000')
    # 旧数据没有 unexamined 键，按完整观测处理。
    positives=[frozenset(r['symptoms']) for r in rows]
    unexamined=[frozenset(r.get('unexamined') or ()) for r in rows]
    items=sorted(set().union(*positives,*unexamined))
    if not 2<=len(items)<=12: raise ValueError('症状种类应为2–12')
    full=(1<<n)-1
    pos_masks={x:sum(1<<i for i,t in enumerate(positives) if x in t) for x in items}
    unx_masks={x:sum(1<<i for i,u in enumerate(unexamined) if x in u) for x in items}
    blocks=defaultdict(list)
    for i,row in enumerate(rows): blocks[row['greenhouse']].append(i)
    # Determine the complete family before applying support/confidence display filters.
    # Rules whose relevant symptoms were never examined together on enough plants are
    # excluded from the test family and reported with their reason.
    family=[]; excluded=[]
    for size in (1,2):
        for antecedent in combinations(items,size):
            a_mask=full
            for x in antecedent: a_mask&=pos_masks[x]
            for consequent in items:
                if consequent in antecedent: continue
                ineligible=0
                for x in (*antecedent,consequent): ineligible|=unx_masks[x]
                eligible=full&~ineligible; ne=eligible.bit_count()
                if ne<MIN_EFFECTIVE:
                    excluded.append({'antecedent':list(antecedent),'consequent':consequent,'effective_n':ne,
                        'reason':'无共同有效样本：相关病征没有同时查过的植株' if ne==0 else f'有效样本仅{ne}株，不足{MIN_EFFECTIVE}株下限'})
                    continue
                ac=(a_mask&eligible).bit_count(); bc=(pos_masks[consequent]&eligible).bit_count()
                joint=(a_mask&pos_masks[consequent]&eligible).bit_count()
                conf=joint/ac if ac else 0
                family.append({'antecedent':list(antecedent),'consequent':consequent,'count':joint,
                    'antecedent_count':ac,'consequent_count':bc,'effective_n':ne,'excluded_n':n-ne,
                    'excluded_by':[x for x in (*antecedent,consequent) if unx_masks[x]],
                    'support':joint/ne,'confidence':conf,'lift':conf/(bc/ne) if ac and bc else 0,
                    '_mask':a_mask&eligible,'_eligible':eligible,'_exceed':0})
    # Permute the consequent within each greenhouse, restricted to each rule's eligible
    # plants. Rules sharing consequent and eligible set share one permutation stream.
    groups=defaultdict(list)
    for rule in family: groups[(rule['consequent'],rule['_eligible'])].append(rule)
    rng=random.Random(seed)
    for consequent in items:
        own=sorted(((key[1],rs) for key,rs in groups.items() if key[0]==consequent),key=lambda t:t[0])
        for eligible,rules in own:
            strata=[]
            for indices in blocks.values():
                ids=[i for i in indices if eligible>>i&1]
                strata.append((ids,sum(consequent in positives[i] for i in ids)))
            for _ in range(permutations):
                mask=0
                for ids,count in strata:
                    for idx in rng.sample(ids,count): mask|=1<<idx
                for rule in rules:
                    if (rule['_mask']&mask).bit_count()>=rule['count']: rule['_exceed']+=1
    for rule in family:
        rule['p_value']=(1+rule['_exceed'])/(permutations+1)
        for key in ('_mask','_eligible','_exceed'): rule.pop(key)
    adjusted=bh_adjust([r['p_value'] for r in family])
    for row,q in zip(family,adjusted): row['q_value']=q
    frequent=apriori(items,pos_masks,unx_masks,full,support,3)
    rules=[r for r in family if r['support']>=support and r['confidence']>=confidence]
    rules.sort(key=lambda r:(r['q_value'],-r['lift'],r['antecedent'],r['consequent']))
    for rank,rule in enumerate(rules,1): rule['rule_id']=f'R{rank:03}'
    frequencies=[{'item':x,'count':pos_masks[x].bit_count(),'examined':n-unx_masks[x].bit_count(),
        'unexamined':unx_masks[x].bit_count(),
        'rate':pos_masks[x].bit_count()/(n-unx_masks[x].bit_count()) if unx_masks[x].bit_count()<n else 0} for x in items]
    heatmap=[]
    for a in items:
        for b in items:
            both=full&~(unx_masks[a]|unx_masks[b])
            heatmap.append({'a':a,'b':b,'count':(pos_masks[a]&pos_masks[b]&both).bit_count(),'effective_n':both.bit_count()})
    observation={'total_cells':n*len(items),'unexamined_cells':sum(len(u) for u in unexamined),
        'complete_plants':sum(1 for u in unexamined if not u),'partial_plants':sum(1 for u in unexamined if u),
        'by_item':[{'item':x,'examined':n-unx_masks[x].bit_count(),'unexamined':unx_masks[x].bit_count()} for x in items],
        'by_greenhouse':[{'name':g,'plants':len(ids),'complete_plants':sum(1 for i in ids if not unexamined[i]),
            'unexamined_cells':sum(len(unexamined[i]) for i in ids)} for g,ids in sorted(blocks.items())]}
    return {'sample_count':n,'item_count':len(items),'family_size':len(family),'frequent_itemsets':len(frequent),
        'rules':rules,'significant_count':sum(r['q_value']<=.05 for r in rules),'frequencies':frequencies,
        'greenhouses':[{'name':g,'count':len(ids)} for g,ids in sorted(blocks.items())], 'heatmap':heatmap,
        'observation':observation,'excluded_rules':excluded,'min_effective':MIN_EFFECTIVE,
        'checksum':checksum(rows),'params':params,
        'method':'逐规则完整观测分母；温室内置换后件且限于相关病征均查过的植株；完整单后件规则族先做BH，再按支持度/置信度展示筛选',
        'caveat':'探索性关联，不代表因果。未查病征不计入阴性，各规则有效样本量不同；BH依赖独立或正相关假设，相关症状规则存在依赖，q值仅供复查排序，需独立样本验证。'}
