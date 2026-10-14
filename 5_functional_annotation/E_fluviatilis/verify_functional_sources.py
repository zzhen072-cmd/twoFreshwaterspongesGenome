"""Count final annotation evidence directly in the filtered raw tool outputs."""
import csv,gzip,json
from pathlib import Path
ROOT=Path('/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010')
DATA=ROOT/'release_v2/functional_annotation'
OUT=ROOT/'manuscript_recheck_20261011'
rep={l[1:].split()[0] for l in open(ROOT/'release_v2/annotation/S04_decontaminated_v2.representative_proteins.faa') if l.startswith('>')}
assert len(rep)==29180
with gzip.open(DATA/'S04_v2_complete_annotation.tsv.gz','rt') as f:records=list(csv.DictReader(f,delimiter='\t'))
byid={r['protein_id']:r for r in records};assert set(byid)==rep
flags={k:{r['protein_id'] for r in records if r[k]=='true'} for k in ['eggNOG','InterPro','InterProScan_any_signature','SwissProt','any_upset_source']}
egg=set();go=set();ko=set();cog=set();pathway=set();header=None
with open(DATA/'S04_v2.longest.emapper.annotations') as f:
    for line in f:
        if line.startswith('#query'):
            header=line.rstrip('\n').lstrip('#').split('\t');continue
        if line.startswith('#') or not line.strip():continue
        assert header is not None
        row=dict(zip(header,line.rstrip('\n').split('\t')))
        q=row['query'];assert q in rep
        good=lambda name:row[name] not in ('','-','NA')
        if good('seed_ortholog'):egg.add(q)
        if good('GOs'):go.add(q)
        if good('KEGG_ko'):ko.add(q)
        if good('COG_category'):cog.add(q)
        if good('KEGG_Pathway'):pathway.add(q)
ipr=set();signatures=set()
with open(DATA/'S04_v2.longest.interproscan.tsv') as f:
    for line in f:
        if not line.strip():continue
        q=line.rstrip('\n').split('\t');assert len(q)>=13 and q[0] in rep
        signatures.add(q[0])
        assert q[11]=='-' or q[11].startswith('IPR')
        if q[11]!='-':ipr.add(q[0])
swiss=[]
with open(DATA/'S04_v2_vs_swissprot.best.tsv') as f:
    for line in f:
        if line.strip() and not line.startswith('#'):swiss.append(line.split('\t',1)[0])
assert len(swiss)==len(set(swiss)) and set(swiss)<=rep
assert egg==flags['eggNOG']
assert ipr==flags['InterPro'] and signatures==flags['InterProScan_any_signature']
assert set(swiss)==flags['SwissProt']
assert egg|ipr|set(swiss)==flags['any_upset_source']
for key,observed in [('GO_eggNOG',go),('KEGG_KO_union',ko),('COG_categories',cog)]:
    expected={r['protein_id'] for r in records if r[key] not in ('','-','NA','[]')}
    assert observed==expected, key
result={'status':'PASS','representative_proteins':len(rep),'eggNOG_seed_ortholog':len(egg),
 'InterPro_IPR_accession':len(ipr),'InterProScan_any_signature':len(signatures),
 'SwissProt_first_hit':len(swiss),'any_of_three_sources':len(egg|ipr|set(swiss)),
 'eggNOG_GO':len(go),'eggNOG_KEGG_KO':len(ko),'eggNOG_COG':len(cog),
 'eggNOG_KEGG_Pathway':len(pathway),
 'all_protein_id_sets_match_frozen_annotation_flags':True}
(OUT/'raw_functional_source_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
