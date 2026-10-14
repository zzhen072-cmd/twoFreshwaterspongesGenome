"""Validate Merqury's installed QV output and distinct reliable-k-mer counts."""
import json, math, re, sys
from pathlib import Path

root=Path(sys.argv[1])
def stats(name):
    values={}
    for line in (root/name).read_text().splitlines():
        m=re.match(r'\s*(\d+)\s+(unique|distinct|total)\b',line,re.I)
        if m:values[m[2].lower()]=int(m[1])
        m=re.match(r'\s*(unique|distinct|present|total)\s+(\d+)\b',line,re.I)
        if m:values['total' if m[1].lower()=='present' else m[1].lower()]=int(m[2])
    if not {'distinct','total'}<=values.keys():
        raise ValueError('Unrecognized meryl statistics format: '+name)
    return values

solid=stats('reads_solid.stats')['distinct']
represented=stats('represented_v2.stats')['distinct']
lost=stats('lost_solid.stats')['distinct']
asm=stats('new_assembly.stats')['total']
errors=stats('new_assembly_only.stats')['total']
assert solid==185575451, 'Reliable read database differs from the original analysis'
assert represented+lost==144606286, 'New representation does not reconcile to frozen original'
assert 0<errors<asm
error=-math.expm1(math.log1p(-errors/asm)/19)
qv=-10*math.log10(error)
q=(root/'S04_v2.qv').read_text().strip().split('\t')
assert len(q)==5 and int(q[1])==errors and int(q[2])==asm
assert abs(float(q[3])-qv)<0.001 and abs(float(q[4])-error)<1e-7
result={'assembly':'S04_decontaminated_v2','k':19,'read_solid_threshold':'multiplicity > 3',
        'QV':qv,'estimated_error_probability':error,'estimated_error_percent':100*error,
        'assembly_kmer_total':asm,'assembly_only_kmer_total':errors,
        'reliable_distinct_read_kmers':solid,'represented_reliable_distinct_kmers':represented,
        'lost_reliable_distinct_kmers':lost,'original_represented_reliable_distinct_kmers':144606286,
        'completeness_percent':100*represented/solid,'reconciliation':'PASS'}
(root/'S04_v2_merqury.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
