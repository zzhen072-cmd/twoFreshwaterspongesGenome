#!/usr/bin/env python3
"""Read-only sequence audit; optional genome reconstruction. Python >=3.9 stdlib."""
import argparse, collections, csv, gzip, hashlib, itertools, json, statistics
from pathlib import Path
from urllib.parse import unquote

def op(p):
    return gzip.open(p, 'rt') if str(p).endswith('.gz') else open(p)

def fasta(p):
    d={}; k=None
    with op(p) as f:
        for l in f:
            if l.startswith('>'):
                k=l[1:].split()[0]
                if k in d: raise ValueError('Duplicate FASTA identifier '+k)
                d[k]=[]
            elif k: d[k].append(l.strip())
    return {k:''.join(v).upper() for k,v in d.items()}

BASES='TCAG'
CODONS=[a+b+c for a in BASES for b in BASES for c in BASES]
TABLE=dict(zip(CODONS,'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
RC=str.maketrans('ACGTNRYKMSWBDHV','TGCANYRMKSWVHDB')
IUPAC={'A':'A','C':'C','G':'G','T':'T','N':'ACGT','R':'AG','Y':'CT','K':'GT','M':'AC','S':'CG','W':'AT','B':'CGT','D':'AGT','H':'ACT','V':'ACG'}
def codon_aa(c):
    if c not in TABLE:
        aas={TABLE[''.join(x)] for x in itertools.product(*(IUPAC.get(b,'ACGT') for b in c))}
        TABLE[c]=next(iter(aas)) if len(aas)==1 else 'X'
    return TABLE[c]
def trans(s):return ''.join(codon_aa(s[i:i+3]) for i in range(0,len(s)-2,3))
def digest(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def audit(args):
    proteins=fasta(args.proteins);cds=fasta(args.cds); reps=fasta(args.representative)
    genes={}; tx={}; seg=collections.defaultdict(lambda:collections.defaultdict(list)); ids=set(); bad=[]; featurecounts=collections.Counter()
    with op(args.gff) as f:
        for n,l in enumerate(f,1):
            if l.startswith('#') or not l.strip():continue
            q=l.rstrip('\n').split('\t'); attr=dict((unquote(k),unquote(v)) for k,v in (a.split('=',1) for a in q[8].split(';') if '=' in a))
            st,en=int(q[3]),int(q[4]);featurecounts[q[2]]+=1
            if 'ID' in attr:
                if attr['ID'] in ids:bad.append({'line':n,'kind':'duplicate_feature_id','id':attr['ID']})
                ids.add(attr['ID'])
            if st>en:bad.append({'line':n,'kind':'reversed_interval','id':attr.get('ID'),'parent':attr.get('Parent')})
            r=(q[0],st,en,q[6],q[7])
            if q[2]=='gene':genes[attr['ID']]=r
            elif q[2] in ('mRNA','transcript'):tx[attr['ID']]=(attr['Parent'],r)
            else:
                for p in attr.get('Parent','').split(','):
                    if p:seg[p][q[2]].append(r)
    genome=fasta(args.genome) if args.genome else None
    stats=collections.Counter();rows=[]
    for t,(g,tr) in tx.items():
        p=proteins.get(t); c=cds.get(t); ss=sorted(seg[t]['CDS'],key=lambda x:x[1],reverse=tr[3]=='-')
        ex=sorted(set((x[1],x[2]) for x in seg[t]['exon']))
        intr=sorted(set((x[1],x[2]) for x in seg[t]['intron']))
        gaps=[(a[1]+1,b[0]-1) for a,b in zip(ex,ex[1:]) if a[1]+1<=b[0]-1]
        rawlen=sum(x[2]-x[1]+1 for x in ss)
        ph=int(ss[0][4]) if ss and ss[0][4] in '012' else None
        phasebad=0; cumulative=0
        for i,s in enumerate(ss):
            expected=ph if i==0 else (-cumulative)%3
            phasebad+=str(expected)!=s[4]
            cumulative+=s[2]-s[1]+1-(ph if i==0 and ph is not None else 0)
        r={'sample':args.sample,'transcript_id':t,'gene_id':g,'representative':int(t in reps),'gff_cds_bp':rawlen,'first_cds_phase':ph,'phase_inconsistency_count':phasebad,'explicit_introns_match_exon_gaps':int(intr==gaps),'protein_present':int(p is not None),'cds_present':int(c is not None)}
        if c is not None and p is not None:
            aa=p.rstrip('*'); translated=trans(c)
            r.update(cds_bp=len(c),protein_aa=len(aa),cds_modulo3=len(c)%3,cds_starts_ATG=int(c.startswith('ATG')),cds_ends_stop=int(c[-3:] in ('TAA','TAG','TGA')),cds_internal_stop_count=translated.rstrip('*').count('*'),protein_internal_stop_count=aa.count('*'),unknown_aa_count=aa.count('X'),translation_matches_protein=int(translated.rstrip('*')==aa),gff_minus_cds_bp=rawlen-len(c),representative_sequence_matches=int(t not in reps or reps[t].rstrip('*')==aa))
        if genome:
            chunks=[];bounds=True
            for seq,st,en,strand,phase in ss:
                if seq not in genome or st<1 or en>len(genome[seq]) or st>en:bounds=False;break
                s=genome[seq][st-1:en];chunks.append(s if strand=='+' else s.translate(RC)[::-1])
            recon=''.join(chunks) if bounds else None
            r['genome_bounds_valid']=int(bounds)
            if recon is not None and c is not None:
                r['raw_gff_extracted_cds_matches']=int(recon==c)
                r['phase_trimmed_gff_cds_matches']=int(recon[ph or 0:]==c)
                r['phase_trimmed_remove_terminal_stop_matches']=int(recon[ph or 0:][:-3]==c and recon[-3:] in ('TAA','TAG','TGA'))
                r['phase_trimmed_translation_matches_protein']=int(trans(recon[ph or 0:]).rstrip('*')==p.rstrip('*')) if p is not None else 0
                # AUGUSTUS/BRAKER FASTA export preserves partial codons by adding
                # terminal Ns; terminal stop codons are omitted from codingseq.
                leading_padding=(3-(ph or 0))%3
                exported='N'*leading_padding+recon
                trailing_padding=(-len(exported))%3
                exported+='N'*trailing_padding
                terminal_stop_removed=exported[-3:] in ('TAA','TAG','TGA')
                if terminal_stop_removed:exported=exported[:-3]
                r.update(leading_N_padding=leading_padding,trailing_N_padding=trailing_padding,terminal_stop_removed=int(terminal_stop_removed),braker_export_convention_matches_cds=int(exported==c),braker_export_convention_matches_protein=int(trans(exported)==p.rstrip('*')))
        rows.append(r)
    metrics=['phase_inconsistency_count','explicit_introns_match_exon_gaps','cds_modulo3','cds_starts_ATG','cds_ends_stop','cds_internal_stop_count','protein_internal_stop_count','unknown_aa_count','translation_matches_protein','gff_minus_cds_bp','representative_sequence_matches','raw_gff_extracted_cds_matches','phase_trimmed_gff_cds_matches','phase_trimmed_translation_matches_protein','genome_bounds_valid','leading_N_padding','trailing_N_padding','terminal_stop_removed','braker_export_convention_matches_cds','braker_export_convention_matches_protein']
    summary={'sample':args.sample,'gff_genes':len(genes),'gff_transcripts':len(tx),'cds_records':len(cds),'protein_records':len(proteins),'representative_records':len(reps),'gff_transcript_ids_equal_cds_ids':set(tx)==set(cds),'gff_transcript_ids_equal_protein_ids':set(tx)==set(proteins),'representative_ids_subset_all_proteins':set(reps)<=set(proteins),'representative_gene_count':len({tx[t][0] for t in reps}),'gff_feature_counts':dict(featurecounts),'gff_feature_issues':bad,'all_transcripts':{},'representatives':{}}
    for population in ['all_transcripts','representatives']:
        group=rows if population=='all_transcripts' else [r for r in rows if r['representative']]
        for metric in metrics:
            if any(metric in r for r in group):summary[population][metric]=dict(collections.Counter(str(r.get(metric,'NA')) for r in group))
    summary['source_files']={k:{'path':str(getattr(args,k)),'sha256':digest(getattr(args,k))} for k in ['gff','cds','proteins','representative']+(['genome'] if args.genome else [])}
    return {'summary':summary,'transcripts':rows}

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    for name in ['sample','gff','cds','proteins','representative']:ap.add_argument('--'+name,required=True)
    ap.add_argument('--genome')
    a=ap.parse_args();print(json.dumps(audit(a),separators=(',',':')))
