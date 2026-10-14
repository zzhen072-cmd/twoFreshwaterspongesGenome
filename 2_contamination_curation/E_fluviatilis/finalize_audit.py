#!/usr/bin/env python3
import csv
import gzip
import os
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(os.environ.get("S04_AUDIT_ROOT", "/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010"))
TABLES = ROOT / "tables"


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


base_rows = read_tsv(TABLES / "S04_21_precoverage_integrated.tsv")
by_id = {r["contig_id"]: r for r in base_rows}


def parse_mosdepth(label):
    means = defaultdict(list)
    breadth = defaultdict(lambda: defaultdict(list))
    with gzip.open(TABLES / f"{label}.regions.bed.gz", "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            contig, start, end, mean = f[0], int(f[1]), int(f[2]), float(f[4])
            means[contig].append((start, end, mean))
    with gzip.open(TABLES / f"{label}.thresholds.bed.gz", "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            contig, start, end = f[0], int(f[1]), int(f[2])
            length = end - start
            vals = [int(x) for x in f[4:]]
            for threshold, value in zip(header[4:], vals):
                breadth[contig][threshold].append(value / length if length else 0.0)
    return means, breadth


platforms = {}
for label in ("illumina_q20_1kb", "ont_q20_1kb"):
    platforms[label] = parse_mosdepth(label)

prodigal = defaultdict(lambda: {"genes": 0, "intervals": []})
gff = TABLES / "S04_21_candidates.prodigal.gff"
if gff.exists():
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5 and f[2] == "CDS":
                prodigal[f[0]]["genes"] += 1
                prodigal[f[0]]["intervals"].append((int(f[3]), int(f[4])))


def merged_bp(intervals):
    if not intervals:
        return 0
    iv = sorted(intervals)
    a, b = iv[0]
    total = 0
    for x, y in iv[1:]:
        if x <= b + 1:
            b = max(b, y)
        else:
            total += b - a + 1
            a, b = x, y
    return total + b - a + 1


kraken = {}
kpath = TABLES / "S04_21_candidates.kraken2.out"
if kpath.exists():
    with open(kpath) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5:
                kraken[f[1]] = {"status": f[0], "taxid": f[2], "kraken_length": f[3]}

extra_cols = []
for prefix in ("illumina_q20", "ont_q20"):
    extra_cols.extend([
        f"{prefix}_weighted_mean", f"{prefix}_window_median", f"{prefix}_window_min", f"{prefix}_window_max",
        f"{prefix}_windows_breadth_ge90pct", f"{prefix}_windows_breadth_le10pct",
    ])
extra_cols += ["prodigal_genes", "prodigal_coding_bp", "prodigal_coding_percent", "kraken_status_new", "kraken_taxid_new"]
columns = list(base_rows[0]) + extra_cols

with open(TABLES / "S04_21_full_integrated.tsv", "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in base_rows:
        contig = row["contig_id"]
        length = int(row["length_bp"])
        for label, prefix in (("illumina_q20_1kb", "illumina_q20"), ("ont_q20_1kb", "ont_q20")):
            means, breadth = platforms[label]
            windows = means.get(contig, [])
            values = [x[2] for x in windows]
            weighted = sum((end - start) * mean for start, end, mean in windows) / length if windows else 0.0
            b1 = breadth[contig].get("1X", [])
            row[f"{prefix}_weighted_mean"] = f"{weighted:.6f}"
            row[f"{prefix}_window_median"] = f"{statistics.median(values):.6f}" if values else ""
            row[f"{prefix}_window_min"] = f"{min(values):.6f}" if values else ""
            row[f"{prefix}_window_max"] = f"{max(values):.6f}" if values else ""
            row[f"{prefix}_windows_breadth_ge90pct"] = sum(x >= 0.9 for x in b1)
            row[f"{prefix}_windows_breadth_le10pct"] = sum(x <= 0.1 for x in b1)
        coding_bp = merged_bp(prodigal[contig]["intervals"])
        row["prodigal_genes"] = prodigal[contig]["genes"]
        row["prodigal_coding_bp"] = coding_bp
        row["prodigal_coding_percent"] = f"{100 * coding_bp / length:.4f}"
        row["kraken_status_new"] = kraken.get(contig, {}).get("status", "")
        row["kraken_taxid_new"] = kraken.get(contig, {}).get("taxid", "")
        writer.writerow(row)

print(f"PASS\tfinalized {len(base_rows)} candidates")
