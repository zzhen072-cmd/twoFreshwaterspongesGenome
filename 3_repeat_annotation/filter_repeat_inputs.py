"""Filter complete contigs only; retain original repeat search/alignment data."""
from pathlib import Path
import collections, csv, gzip, hashlib, json, re, sys

ROOT=Path('/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010/manuscript_recheck_20261011')
BASE=Path('/work1/sponge/genome_annotation')
OUT=ROOT/'repeats';OUT.mkdir(exist_ok=True)
remove=set((ROOT/'scripts/removed_short_ids.txt').read_text().split())
assert len(remove)==21
mapping={r['short_id']:r for r in csv.DictReader(open(ROOT/'scripts/S04_contig_id_mapping.tsv'),delimiter='\t')}
assert len(mapping)==2118
def fa(p):
    h=None;seq=[]
    with open(p) as f:
        for l in f:
            if l.startswith('>'):
                if h:yield h,''.join(seq).upper()
                h=l[1:].split()[0];seq=[]
            else:seq.append(l.strip())
        if h:yield h,''.join(seq).upper()
count=length=0;seen=set()
with open(OUT/'S04_v2_shortids.fasta','x') as f:
    for name,seq in fa(BASE/'S04/S04_shortids.fasta'):
        assert name in mapping and name not in seen
        seen.add(name)
        assert hashlib.sha256(seq.encode()).hexdigest()==mapping[name]['sequence_sha256']
        if name in remove:continue
        count+=1;length+=len(seq)
        f.write('>'+name+'\n')
        for j in range(0,len(seq),80):f.write(seq[j:j+80]+'\n')
assert len(seen)==2118 and count==2097 and length==267354623

header=re.compile(r'^\s*\d+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+(\S+)\s+\d+\s+\d+\s+\(')
metadata=re.compile(r'^(##\s*Total\s+(Sequences|Length|NonMask[^:]*|NonSub[^:]*):\s*)(\d+)')
def filter_align(src,dest,is_cat=False):
    opener=gzip.open if str(src).endswith('.gz') else open
    counts=collections.Counter();keep=True;started=False
    with opener(src,'rt') as inp, (gzip.open(dest,'wt',compresslevel=3) if str(dest).endswith('.gz') else open(dest,'x')) as out:
        for line in inp:
            if line.startswith('##'):
                m=metadata.match(line)
                if m:
                    val=count if m[2]=='Sequences' else length
                    line=line[:m.start(3)]+str(val)+line[m.end(3):]
                m=re.match(r'^##\s+(S04_contig_\d+)\s+',line)
                if m and m[1] in remove:continue
                out.write(line);continue
            m=header.match(line)
            if m:
                name=m[1];started=True
                if name in mapping:keep=name not in remove
                elif name=='refinement' and is_cat:keep=True
                else:raise ValueError('Unexpected alignment query: '+name)
                counts['retained_blocks' if keep else 'removed_blocks']+=1
                counts['refinement_blocks']+=name=='refinement'
            if keep:out.write(line)
    assert started and counts['removed_blocks']>0
    return dict(counts)

summary={'genome_contigs':count,'genome_bp':length,'removed_contigs':len(remove)}
summary['cat']=filter_align(BASE/'S04/repeatmasker/S04_shortids.fasta.cat.gz',OUT/'S04_v2_shortids.fasta.cat.gz',True)
summary['align']=filter_align(BASE/'repeat_landscapes_20260922/S04/S04_shortids.fasta.align',OUT/'S04_v2.filtered_original.align')
(OUT/'input_filter_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
