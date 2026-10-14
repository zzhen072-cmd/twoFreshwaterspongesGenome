#!/usr/bin/env python3
import csv
import gzip
import hashlib
import os
import re
from pathlib import Path

ROOT = Path(os.environ.get("S04_AUDIT_ROOT", "/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010"))
SCRIPT_DIR = Path(__file__).resolve().parent
REMOVAL_LIST = Path(os.environ.get(
    "E_FLUVIATILIS_REMOVAL_LIST",
    SCRIPT_DIR / "E_fluviatilis_removed_21_contigs.txt",
))
INPUT = ROOT / "inputs"
TABLES = ROOT / "tables"
RELEASE = ROOT / "release_v2"
ANNOT = RELEASE / "annotation"
QC = RELEASE / "qc"
for directory in (RELEASE, ANNOT, QC):
    directory.mkdir(parents=True, exist_ok=True)

SOURCE_FASTA = Path("/public/home/sponge/sci_data/S04/00_input/S04_clean.fasta")
SOURCE_GFF = Path("/work1/sponge/genome_annotation/deep_validation_20260923/platform_taxonomy/source/frozen/S04.annotation.original_contig_ids.gff3.gz")
SOURCE_REP_PROTEIN = Path("/work1/sponge/genome_annotation/braker4_ep/S04/qc_busco_longest/S04_braker.longest_isoform.aa")
EXPECTED_FASTA_SHA256 = "f665f1c20f2f4a322a4f9afe99798b737269629fbabfb9d8020875278d355678"

HOST_FASTA = RELEASE / "S04_decontaminated_v2.fasta"
REMOVED_FASTA = RELEASE / "S04_removed_21_nonhost_candidates.fasta"
HOST_GFF = ANNOT / "S04_decontaminated_v2.annotation.gff3"
HOST_GFF_GZ = Path(str(HOST_GFF) + ".gz")
HOST_REP_PROTEIN = ANNOT / "S04_decontaminated_v2.representative_proteins.faa"
REMOVED_REP_PROTEIN = ANNOT / "S04_removed_21.representative_proteins.faa"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def fasta_records(path):
    name = None
    header = None
    parts = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, header, "".join(parts)
                header = line[1:].rstrip("\n")
                name = header.split()[0]
                parts = []
            else:
                parts.append(line.strip())
        if name is not None:
            yield name, header, "".join(parts)


def write_record(fh, header, seq):
    fh.write(f">{header}\n")
    for i in range(0, len(seq), 80):
        fh.write(seq[i:i + 80] + "\n")


remove_ids = [x.strip() for x in open(REMOVAL_LIST) if x.strip()]
remove_set = set(remove_ids)
assert len(remove_ids) == len(remove_set) == 21
assert sha256(SOURCE_FASTA) == EXPECTED_FASTA_SHA256

host_n = host_bp = removed_n = removed_bp = 0
seen = set()
with open(HOST_FASTA, "w") as host_out, open(REMOVED_FASTA, "w") as removed_out:
    for name, header, seq in fasta_records(SOURCE_FASTA):
        if name in remove_set:
            write_record(removed_out, header, seq)
            removed_n += 1
            removed_bp += len(seq)
            seen.add(name)
        else:
            write_record(host_out, header, seq)
            host_n += 1
            host_bp += len(seq)

assert seen == remove_set, sorted(remove_set - seen)
assert removed_n == 21 and removed_bp == 1228587, (removed_n, removed_bp)
assert host_n == 2097 and host_bp == 267354623, (host_n, host_bp)

removed_genes = set()
retained_genes = set()
feature_counts = {"retained": {}, "removed": {}}

def bump(group, feature):
    feature_counts[group][feature] = feature_counts[group].get(feature, 0) + 1

with open_text(SOURCE_GFF) as src, open(HOST_GFF, "w") as dst:
    for line in src:
        if line.startswith("##sequence-region"):
            fields = line.rstrip("\n").split()
            if len(fields) >= 2 and fields[1] in remove_set:
                continue
            dst.write(line)
            continue
        if line.startswith("#") or not line.strip():
            dst.write(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        group = "removed" if fields[0] in remove_set else "retained"
        bump(group, fields[2])
        if fields[2] == "gene":
            match = re.search(r"(?:^|;)ID=([^;]+)", fields[8])
            if match:
                (removed_genes if group == "removed" else retained_genes).add(match.group(1))
        if group == "retained":
            dst.write(line)

assert removed_genes.isdisjoint(retained_genes)
assert len(removed_genes) == 644, len(removed_genes)

with open(HOST_GFF, "rb") as src, gzip.open(HOST_GFF_GZ, "wb", compresslevel=9) as dst:
    for block in iter(lambda: src.read(1024 * 1024), b""):
        dst.write(block)

rep_total = rep_host = rep_removed = 0
with open(HOST_REP_PROTEIN, "w") as host_out, open(REMOVED_REP_PROTEIN, "w") as removed_out:
    for name, header, seq in fasta_records(SOURCE_REP_PROTEIN):
        rep_total += 1
        match = re.search(r"(?:^|\s)gene=([^\s]+)", header)
        gene = match.group(1) if match else name.split(".", 1)[0]
        if gene in removed_genes:
            write_record(removed_out, header, seq)
            rep_removed += 1
        else:
            write_record(host_out, header, seq)
            rep_host += 1

assert rep_total == 29824, rep_total
assert rep_removed == 644 and rep_host == 29180, (rep_removed, rep_host)

with open(TABLES / "S04_21_full_integrated.tsv") as fh, open(TABLES / "S04_21_final_decisions.tsv", "w", newline="") as out:
    reader = csv.DictReader(fh, delimiter="\t")
    columns = reader.fieldnames + ["final_decision", "decision_confidence", "decision_basis"]
    writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    count = 0
    for row in reader:
        count += 1
        blast = float(row["blastn_union_percent"] or 0)
        coding = float(row["prodigal_coding_percent"] or 0)
        strong = int(row["strong_nonhost_genes"] or 0)
        viral = int(row["viral_keyword_genes"] or 0)
        confidence = "very_high" if (blast >= 10 or coding >= 80 or strong + viral >= 4) else "high"
        basis = (
            f"Tiara={row['tiara_class']}; FCS={row['fcs_division'] or 'none'}; "
            f"BLASTN_union={blast:.4f}% ({row['blastn_phyla'] or 'no_nt_hit'}); "
            f"strong_nonhost_genes={strong}; strong_metazoa_genes={row['strong_metazoa_genes']}; "
            f"viral_keyword_genes={viral}; Prodigal_coding={coding:.4f}%; "
            f"Illumina_Q20_median={row['illumina_q20_window_median']}x, "
            f"Illumina_Q20_max={row['illumina_q20_window_max']}x, "
            f"ONT_Q20_median={row['ont_q20_window_median']}x"
        )
        row.update(final_decision="remove_from_host_assembly_and_archive", decision_confidence=confidence, decision_basis=basis)
        writer.writerow(row)
    assert count == 21

with open(TABLES / "S04_v2_build_summary.tsv", "w") as out:
    out.write("metric\tvalue\n")
    values = {
        "source_contigs": host_n + removed_n,
        "source_bp": host_bp + removed_bp,
        "removed_contigs": removed_n,
        "removed_bp": removed_bp,
        "retained_contigs": host_n,
        "retained_bp": host_bp,
        "source_representative_proteins": rep_total,
        "removed_representative_proteins": rep_removed,
        "retained_representative_proteins": rep_host,
        "source_gff_genes": len(retained_genes) + len(removed_genes),
        "removed_gff_genes": len(removed_genes),
        "retained_gff_genes": len(retained_genes),
    }
    for key, value in values.items():
        out.write(f"{key}\t{value}\n")
    for group in ("retained", "removed"):
        for feature, value in sorted(feature_counts[group].items()):
            out.write(f"{group}_gff_{feature}\t{value}\n")

with open(TABLES / "S04_removed_gene_ids.txt", "w") as out:
    out.write("\n".join(sorted(removed_genes)) + "\n")
with open(TABLES / "S04_retained_gene_ids.txt", "w") as out:
    out.write("\n".join(sorted(retained_genes)) + "\n")

with open(RELEASE / "README.txt", "w") as out:
    out.write("S04 post hoc decontaminated release v2\n")
    out.write("The frozen source assembly is unchanged. Exactly 21 whole contigs were removed from the host assembly and archived separately.\n")
    out.write(f"Source assembly: {SOURCE_FASTA}\n")
    out.write(f"Source SHA256: {EXPECTED_FASTA_SHA256}\n")
    out.write("Decision table: ../tables/S04_21_final_decisions.tsv\n")

manifest_paths = [HOST_FASTA, REMOVED_FASTA, HOST_GFF, HOST_GFF_GZ, HOST_REP_PROTEIN, REMOVED_REP_PROTEIN,
                  TABLES / "S04_21_final_decisions.tsv", TABLES / "S04_v2_build_summary.tsv"]
with open(RELEASE / "SHA256SUMS.build.txt", "w") as out:
    for path in manifest_paths:
        out.write(f"{sha256(path)}  {path}\n")

print(f"PASS build_v2 retained={host_n} contigs/{host_bp} bp removed={removed_n} contigs/{removed_bp} bp rep_proteins={rep_host}")
