#!/usr/bin/env python3
import csv
import gzip
import os
from collections import Counter
from pathlib import Path

ROOT = Path(os.environ.get("S04_AUDIT_ROOT", "/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010"))
REPO_ROOT = Path(__file__).resolve().parents[2]
REMOVAL_LIST = Path(os.environ.get(
    "E_FLUVIATILIS_REMOVAL_LIST",
    REPO_ROOT / "2_contamination_curation/E_fluviatilis/E_fluviatilis_removed_21_contigs.txt",
))
RELEASE = ROOT / "release_v2"
QC = RELEASE / "qc"
TABLES = ROOT / "tables"
FASTA = RELEASE / "S04_decontaminated_v2.fasta"
GFF = RELEASE / "annotation/S04_decontaminated_v2.annotation.gff3.gz"
REP = RELEASE / "annotation/S04_decontaminated_v2.representative_proteins.faa"
REMOVED_IDS = set(x.strip() for x in open(REMOVAL_LIST) if x.strip())
OLD_ALL = Path("/work1/sponge/genome_annotation/deep_validation_20260923/platform_taxonomy/tables/S04_all_contigs_integrated_evidence.tsv")


def fasta_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


ids = fasta_ids(FASTA)
id_set = set(ids)
assert len(ids) == len(id_set) == 2097
assert not (id_set & REMOVED_IDS)

gff_seqids = set()
gff_genes = 0
with gzip.open(GFF, "rt") as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t", 3)
        gff_seqids.add(f[0])
        if len(f) > 2 and f[2] == "gene":
            gff_genes += 1
assert gff_seqids <= id_set
assert not (gff_seqids & REMOVED_IDS)
assert gff_genes == 29180

rep_ids = fasta_ids(REP)
assert len(rep_ids) == len(set(rep_ids)) == 29180

tiara_counts = Counter()
recommendation_counts = Counter()
retained_review_rows = []
with open(OLD_ALL) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    fields = reader.fieldnames
    for row in reader:
        if row["contig_id"] not in id_set:
            continue
        tiara_counts[row["classification"]] += 1
        recommendation_counts[row["recommendation"]] += 1
        if row["recommendation"] != "not_in_140_candidate_set":
            retained_review_rows.append(row)

with open(TABLES / "S04_v2_remaining_review_flags.tsv", "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(retained_review_rows)

kraken_status = Counter()
kraken_names = Counter()
kraken_path = QC / "S04_decontaminated_v2.kraken2.out"
if kraken_path.exists():
    with open(kraken_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                kraken_status[f[0]] += 1
                if f[0] == "C":
                    kraken_names[f[2]] += 1

with open(QC / "S04_v2_validation_summary.tsv", "w") as out:
    out.write("section\tmetric\tvalue\n")
    base = {
        "assembly_contigs": len(ids),
        "removed_ids_present": len(id_set & REMOVED_IDS),
        "gff_distinct_seqids": len(gff_seqids),
        "gff_genes": gff_genes,
        "representative_proteins": len(rep_ids),
        "remaining_review_flag_contigs": len(retained_review_rows),
    }
    for key, value in base.items():
        out.write(f"validation\t{key}\t{value}\n")
    for key, value in sorted(tiara_counts.items()):
        out.write(f"existing_tiara\t{key}\t{value}\n")
    for key, value in sorted(recommendation_counts.items()):
        out.write(f"existing_review_recommendation\t{key}\t{value}\n")
    for key, value in sorted(kraken_status.items()):
        out.write(f"new_kraken_status\t{key}\t{value}\n")
    for key, value in kraken_names.most_common(30):
        out.write(f"new_kraken_top_assignment\t{key}\t{value}\n")

print(f"PASS v2_validation contigs={len(ids)} genes={gff_genes} proteins={len(rep_ids)} review_flags={len(retained_review_rows)}")
