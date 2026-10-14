"""Ensure post-processing changed only the requested contig membership."""
import collections, json, re
from pathlib import Path
ROOT=Path('/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010/manuscript_recheck_20261011')
OUT=ROOT/'repeats'
remove=set((ROOT/'scripts/removed_short_ids.txt').read_text().split())
def annotations(path,drop=False):
    rows=collections.Counter();intervals=collections.defaultdict(list)
    for line in open(path):
        q=line.split()
        if len(q)<15 or not q[0].isdigit():continue
        if drop and q[4] in remove:continue
        assert q[4] not in remove
        # IDs are renumbered by ProcessRepeats; compare every other field.
        rows[tuple(q[:14]+q[15:])]+=1
        intervals[q[4]].append((int(q[5]),int(q[6])))
    union=0
    for ints in intervals.values():
        a=b=None
        for st,en in sorted(ints):
            if a is None:a,b=st,en
            elif st<=b+1:b=max(b,en)
            else:union+=b-a+1;a,b=st,en
        if a is not None:union+=b-a+1
    return rows,union
old,oldunion=annotations('/work1/sponge/genome_annotation/S04/repeatmasker/S04_shortids.fasta.out',True)
new,newunion=annotations(OUT/'S04_v2_shortids.fasta.out')
assert old==new, 'Postprocessed repeat coordinates/classes changed beyond contig filtering'
tbl=(OUT/'S04_v2_shortids.fasta.tbl').read_text()
masked=int(re.search(r'bases masked:\s+(\d+)',tbl)[1])
bp=int(re.search(r'total length:\s+(\d+)',tbl)[1])
assert bp==267354623 and masked==newunion==oldunion
result={'retained_repeat_rows':sum(new.values()),'all_coordinates_and_classes_unchanged':True,
        'assembly_bp':bp,'masked_union_bp':masked,'masked_percent_exact':100*masked/bp,
        'original_masked_bp':136372399,'masked_bp_removed':136372399-masked}
(OUT/'repeat_output_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
