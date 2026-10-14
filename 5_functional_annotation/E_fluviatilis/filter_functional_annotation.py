#!/usr/bin/env python3
import csv
import gzip
import hashlib
import os
import re
from pathlib import Path

ROOT = Path(os.environ.get("S04_AUDIT_ROOT", "/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010"))
OUT = ROOT / "release_v2/functional_annotation"
OUT.mkdir(parents=True, exist_ok=True)
REP = ROOT / "release_v2/annotation/S04_decontaminated_v2.representative_proteins.faa"

SOURCES = {
    "eggnog_annotations": Path("/work1/sponge/genome_annotation/eggnog_mapper/S04/results/S04.longest.emapper.annotations"),
    "eggnog_hits": Path("/work1/sponge/genome_annotation/eggnog_mapper/S04/results/S04.longest.emapper.hits"),
    "eggnog_seed_orthologs": Path("/work1/sponge/genome_annotation/eggnog_mapper/S04/results/S04.longest.emapper.seed_orthologs"),
    "interproscan": Path("/work1/sponge/genome_annotation/eggnog_mapper/S04/interpro/S04.longest.interproscan.tsv"),
    "swissprot_best": Path("/work1/sponge/genome_annotation/eggnog_mapper/S04/swissprot/S04_vs_swissprot.best.tsv"),
    "swissprot_all": Path("/work1/sponge/genome_annotation/eggnog_mapper/S04/swissprot/S04_vs_swissprot.tsv"),
    "complete_annotation": Path("/work1/sponge/genome_annotation/deep_validation_20260923/platform_taxonomy/source/frozen/S04_complete_annotation.tsv.gz"),
}

OUTPUT_NAMES = {
    "eggnog_annotations": "S04_v2.longest.emapper.annotations",
    "eggnog_hits": "S04_v2.longest.emapper.hits",
    "eggnog_seed_orthologs": "S04_v2.longest.emapper.seed_orthologs",
    "interproscan": "S04_v2.longest.interproscan.tsv",
    "swissprot_best": "S04_v2_vs_swissprot.best.tsv",
    "swissprot_all": "S04_v2_vs_swissprot.tsv",
    "complete_annotation": "S04_v2_complete_annotation.tsv.gz",
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


protein_ids = set()
gene_ids = set()
with open(REP) as fh:
    for line in fh:
        if not line.startswith(">"):
            continue
        header = line[1:].strip()
        protein_ids.add(header.split()[0])
        match = re.search(r"(?:^|\s)gene=([^\s]+)", header)
        if match:
            gene_ids.add(match.group(1))
assert len(protein_ids) == len(gene_ids) == 29180


def filter_first_column(source, destination):
    source_rows = kept_rows = 0
    source_queries = set()
    kept_queries = set()
    with open(source) as src, open(destination, "w") as dst:
        for line in src:
            if line.startswith("#") or not line.strip():
                dst.write(line)
                continue
            source_rows += 1
            query = line.split("\t", 1)[0]
            source_queries.add(query)
            if query in protein_ids:
                dst.write(line)
                kept_rows += 1
                kept_queries.add(query)
    return source_rows, kept_rows, len(source_queries), len(kept_queries)


stats = []
for role in ("eggnog_annotations", "eggnog_hits", "eggnog_seed_orthologs", "interproscan", "swissprot_best", "swissprot_all"):
    values = filter_first_column(SOURCES[role], OUT / OUTPUT_NAMES[role])
    stats.append((role,) + values)

src_path = SOURCES["complete_annotation"]
dst_path = OUT / OUTPUT_NAMES["complete_annotation"]
source_rows = kept_rows = 0
boolean_counts_before = {}
boolean_counts_after = {}
with gzip.open(src_path, "rt") as src, gzip.open(dst_path, "wt", compresslevel=9) as dst:
    reader = csv.DictReader(src, delimiter="\t")
    writer = csv.DictWriter(dst, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    bool_fields = ["eggNOG", "InterProScan_any_signature", "InterPro", "Pfam", "SwissProt", "any_upset_source"]
    boolean_counts_before = {field: 0 for field in bool_fields}
    boolean_counts_after = {field: 0 for field in bool_fields}
    for row in reader:
        source_rows += 1
        for field in bool_fields:
            boolean_counts_before[field] += int(row.get(field, "").lower() == "true")
        if row["protein_id"] in protein_ids and row["gene_id"] in gene_ids:
            writer.writerow(row)
            kept_rows += 1
            for field in bool_fields:
                boolean_counts_after[field] += int(row.get(field, "").lower() == "true")
assert source_rows == 29824 and kept_rows == 29180, (source_rows, kept_rows)
stats.append(("complete_annotation", source_rows, kept_rows, source_rows, kept_rows))

with open(OUT / "S04_v2_functional_filter_summary.tsv", "w") as out:
    out.write("file_role\tsource_rows\tretained_rows\tsource_unique_queries\tretained_unique_queries\n")
    for row in stats:
        out.write("\t".join(map(str, row)) + "\n")
    out.write("\nmetric\tsource_count\tretained_count\n")
    out.write(f"representative_models\t29824\t{len(protein_ids)}\n")
    for field in boolean_counts_before:
        out.write(f"complete_annotation_{field}_true\t{boolean_counts_before[field]}\t{boolean_counts_after[field]}\n")

with open(OUT / "SOURCE_AND_OUTPUT_SHA256.tsv", "w") as out:
    out.write("role\tsource_path\tsource_sha256\tfiltered_path\tfiltered_sha256\n")
    for role, source in SOURCES.items():
        destination = OUT / OUTPUT_NAMES[role]
        out.write(f"{role}\t{source}\t{sha256(source)}\t{destination}\t{sha256(destination)}\n")

print(f"PASS functional_filter retained_models={len(protein_ids)}")
