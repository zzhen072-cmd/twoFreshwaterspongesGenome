#!/usr/bin/env python3
"""Identify contig renaming by exact uppercase sequence SHA256; export CDS union."""
import csv,gzip,hashlib,json,collections
from pathlib import Path
ROOT=Path('/work1/sponge/genome_annotation')
OUT=ROOT/'manuscript_validation_20260923/annotation'
OUT.mkdir(parents=True,exist_ok=True)
def fa_hashes(p):
    records=[];k=None;h=None;n=0
    with open(p) as f:
        for l in f:
            if l.startswith('>'):
                if k:records.append((k,n,h.hexdigest()))
                k=l[1:].split()[0];h=hashlib.sha256();n=0
            else:
                s=l.strip().upper();h.update(s.encode());n+=len(s)
    if k:records.append((k,n,h.hexdigest()))
    return records
results={}
for sp,filename in [('gzx','gzx_cleaned.fasta'),('S04','S04_clean.fasta')]:
    short=fa_hashes(ROOT/sp/(sp+'_shortids.fasta'))
    original=fa_hashes(Path('/public/home/sponge/sci_data')/sp/'00_input'/filename)
    index=collections.defaultdict(list)
    for k,n,h in original:index[(n,h)].append(k)
    mapping={}
    for k,n,h in short:
        assert len(index[(n,h)])==1,('nonunique_or_missing_contig',sp,k)
        mapping[k]=index[(n,h)][0]
    assert len(set(mapping.values()))==len(original)==len(short)
    with (OUT/(sp+'_contig_id_mapping.tsv')).open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['short_id','original_id','length_bp','sequence_sha256'])
        for k,n,h in short:w.writerow([k,mapping[k],n,h])
    qc='qc_busco_protein' if sp=='gzx' else 'qc_busco_longest'
    rp=ROOT/'braker4_ep'/sp/qc/(sp+'_braker.longest_isoform.aa')
    reps={l[1:].split()[0] for l in rp.open() if l.startswith('>')}
    segments=collections.defaultdict(list)
    with gzip.open(ROOT/'braker4_ep'/sp/'output'/sp/'results/braker.gff3.gz','rt') as f:
        for l in f:
            if l.startswith('#'):continue
            x=l.rstrip().split('\t')
            if x[2]!='CDS':continue
            a=dict(s.split('=',1) for s in x[8].split(';') if '=' in s)
            if any(t in reps for t in a.get('Parent','').split(',')):segments[mapping[x[0]]].append((int(x[3])-1,int(x[4])))
    merged=[]
    for contig in sorted(segments):
        for st,en in sorted(segments[contig]):
            if merged and merged[-1][0]==contig and st<=merged[-1][2]:merged[-1][2]=max(en,merged[-1][2])
            else:merged.append([contig,st,en])
    with (OUT/(sp+'_representative_CDS_union.original_ids.bed')).open('w') as f:
        for row in merged:f.write('\t'.join(map(str,row))+'\n')
    results[sp]={'contigs':len(short),'all_original_and_short_contigs_bijectively_sequence_identical':True,'representative_models':len(reps),'merged_cds_intervals':len(merged),'merged_cds_bp':sum(en-st for _,st,en in merged),'mapping':str(OUT/(sp+'_contig_id_mapping.tsv')),'cds_union_bed':str(OUT/(sp+'_representative_CDS_union.original_ids.bed'))}
(OUT/'contig_mapping_validation.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
