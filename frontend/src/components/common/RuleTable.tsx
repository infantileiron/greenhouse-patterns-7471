import type {Rule} from '../../types/domain';
import {pct} from '../../utils/format';
export function RuleTable({rules,selected=[],onToggle,sampleN}:{rules:Rule[];selected?:string[];onToggle?:(id:string)=>void;sampleN?:number}){
 const effective=(r:Rule)=>r.effective_n??sampleN;
 return <div className="table-wrap"><table><thead><tr>{onToggle&&<th>选取</th>}<th>病征组合 → 后件</th><th>支持数</th><th>有效样本</th><th>置信度</th><th>提升度</th><th>q 值</th></tr></thead><tbody>{rules.map(r=><tr key={r.rule_id}>{onToggle&&<td><input aria-label={`选择${r.rule_id}`} type="checkbox" checked={selected.includes(r.rule_id)} onChange={()=>onToggle(r.rule_id)}/></td>}<td><small>{r.rule_id}</small><b>{r.antecedent.join(' + ')} <em>→</em> {r.consequent}</b></td><td>{r.count}</td><td>{effective(r)==null?'—':effective(r)}{(r.excluded_n??0)>0&&<small>排除 {r.excluded_n}{r.excluded_by?.length?`（${r.excluded_by.join('、')}未查）`:''}</small>}</td><td>{pct(r.confidence)}</td><td>{r.lift.toFixed(2)}×</td><td className={r.q_value<=.05?'good':''}>{r.q_value.toFixed(4)}</td></tr>)}</tbody></table></div>
}
