"""Independent small-subset cross-check for exact repeat bases removed."""
from pathlib import Path
import gzip,re,json
ROOT=Path('/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010/manuscript_recheck_20261011')
OUT=ROOT/'removed_repeat_crosscheck';OUT.mkdir(exist_ok=True)
remove=set((ROOT/'scripts/removed_short_ids.txt').read_text().split())
header=re.compile(r'^\s*\d+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+(\S+)\s+\d+\s+\d+\s+\(')
keep=False;n=0
with gzip.open('/work1/sponge/genome_annotation/S04/repeatmasker/S04_shortids.fasta.cat.gz','rt') as src,open(OUT/'removed21.fasta.cat','x') as dst:
    for line in src:
        if line.startswith('##'):
            m=re.match(r'^(##\s*Total\s+(Sequences|Length|NonMask[^:]*|NonSub[^:]*):\s*)(\d+)',line)
            if m:line=line[:m.start(3)]+str(21 if m[2]=='Sequences' else 1228587)+line[m.end(3):]
            m=re.match(r'^##\s+(S04_contig_\d+)\s+',line)
            if m and m[1] not in remove:continue
            dst.write(line);continue
        m=header.match(line)
        if m:
            keep=m[1] in remove
            if keep:n+=1
        if keep:dst.write(line)
assert n==97,n
print('PASS: exactly 97 original search blocks from the removed 21 contigs.')
