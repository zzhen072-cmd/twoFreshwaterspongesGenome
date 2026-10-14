#!/usr/bin/env python3
import csv
import gzip
import hashlib
import math
import os
import statistics
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(os.environ.get("S04_AUDIT_ROOT", "/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010"))
SCRIPT_DIR = Path(__file__).resolve().parent
REMOVAL_LIST = Path(os.environ.get(
    "E_FLUVIATILIS_REMOVAL_LIST",
    SCRIPT_DIR / "E_fluviatilis_removed_21_contigs.txt",
))
INPUT = ROOT / "inputs"
TABLES = ROOT / "tables"
SEQDIR = ROOT / "sequences"
TABLES.mkdir(parents=True, exist_ok=True)
SEQDIR.mkdir(parents=True, exist_ok=True)

FASTA = Path("/public/home/sponge/sci_data/S04/00_input/S04_clean.fasta")
FAI = Path(str(FASTA) + ".fai")
GFF = Path("/work1/sponge/genome_annotation/deep_validation_20260923/platform_taxonomy/source/frozen/S04.annotation.original_contig_ids.gff3.gz")
CONTIG_TABLE = Path("/work1/sponge/genome_annotation/deep_validation_20260923/platform_taxonomy/tables/S04_140_candidate_decisions.tsv")
GENE_TABLE = Path("/work1/sponge/genome_annotation/deep_validation_20260923/platform_taxonomy/tables/S04_candidate_gene_evidence.tsv")
EXPECTED_FASTA_SHA256 = "f665f1c20f2f4a322a4f9afe99798b737269629fbabfb9d8020875278d355678"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def merge_bp(intervals):
    if not intervals:
        return 0
    intervals = sorted((min(a, b), max(a, b)) for a, b in intervals)
    total = 0
    left, right = intervals[0]
    for a, b in intervals[1:]:
        if a <= right + 1:
            right = max(right, b)
        else:
            total += right - left + 1
            left, right = a, b
    return total + right - left + 1


candidates = [x.strip() for x in open(REMOVAL_LIST) if x.strip()]
candidate_set = set(candidates)
assert len(candidates) == 21 and len(candidate_set) == 21
actual_sha = sha256(FASTA)
assert actual_sha == EXPECTED_FASTA_SHA256, (actual_sha, EXPECTED_FASTA_SHA256)

lengths = {}
with open(FAI) as fh:
    for line in fh:
        fields = line.rstrip("\n").split("\t")
        lengths[fields[0]] = int(fields[1])
missing = candidate_set - lengths.keys()
assert not missing, sorted(missing)

with open(TABLES / "candidate_windows_1kb.bed", "w") as out:
    for contig in candidates:
        for start in range(0, lengths[contig], 1000):
            end = min(start + 1000, lengths[contig])
            out.write(f"{contig}\t{start}\t{end}\t{contig}:{start + 1}-{end}\n")

seqs = {}
name = None
parts = []
with open(FASTA) as fh:
    for line in fh:
        if line.startswith(">"):
            if name in candidate_set:
                seqs[name] = "".join(parts).upper()
            name = line[1:].split()[0]
            parts = []
        elif name in candidate_set:
            parts.append(line.strip())
    if name in candidate_set:
        seqs[name] = "".join(parts).upper()
assert set(seqs) == candidate_set

with open(SEQDIR / "S04_21_candidates.fasta", "w") as out:
    for contig in candidates:
        out.write(f">{contig}\n")
        seq = seqs[contig]
        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")

seqstats = {}
for contig, seq in seqs.items():
    counts = Counter(seq)
    acgt = sum(counts[b] for b in "ACGT")
    gc = (counts["G"] + counts["C"]) / acgt if acgt else 0.0
    probs = [counts[b] / acgt for b in "ACGT" if counts[b]]
    entropy = -sum(p * math.log2(p) for p in probs)
    seqstats[contig] = {
        "sequence_length": len(seq),
        "sequence_gc_percent": 100 * gc,
        "sequence_n_bases": counts["N"],
        "sequence_acgt_entropy": entropy,
    }

blast_intervals = defaultdict(list)
blast_by_phylum = defaultdict(lambda: defaultdict(list))
blast_rows = defaultdict(int)
blast_max_identity = defaultdict(float)
blast_max_bitscore = defaultdict(float)
blast_accessions = defaultdict(set)
with open(INPUT / "S04_21_blastn_hsps_annotated.tsv") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        contig = row["contig_id"]
        if contig not in candidate_set:
            continue
        a, b = int(row["qstart"]), int(row["qend"])
        blast_intervals[contig].append((a, b))
        blast_by_phylum[contig][row["phylum"]].append((a, b))
        blast_rows[contig] += 1
        blast_max_identity[contig] = max(blast_max_identity[contig], float(row["pident"]))
        blast_max_bitscore[contig] = max(blast_max_bitscore[contig], float(row["bitscore"]))
        blast_accessions[contig].add(row["subject_accession"])

gene_coords = {}
with open_text(GFF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] not in candidate_set or f[2] != "gene":
            continue
        attrs = dict(item.split("=", 1) for item in f[8].split(";") if "=" in item)
        gene_coords[attrs["ID"]] = (f[0], int(f[3]), int(f[4]))

gene_summary = defaultdict(Counter)
evidence_positions = defaultdict(list)
swiss_groups = defaultdict(Counter)
with open(GENE_TABLE) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        contig = row["contig_id"]
        if contig not in candidate_set:
            continue
        s = gene_summary[contig]
        s["genes"] += 1
        s["single_exon"] += int(row["exons"] == "1")
        s["protein_hint"] += int(row["protein_hint_support"] == "1")
        s["eggnog_metazoa"] += int(row["eggNOG_metazoa_assignment"] == "true")
        s["interpro"] += int(row["InterPro_IPR"] == "true")
        group = row["Swiss_broad_group"] or "Unresolved"
        swiss_groups[contig][group] += int(bool(row["Swiss_taxon"]))
        strong = row["Swiss_strong_operational_hit"] == "1"
        nonhost = strong and group not in {"Metazoa", "Unresolved"}
        metazoa = strong and group == "Metazoa"
        viral = row["phage_or_viral_keyword"] == "1"
        s["strong_nonhost"] += int(nonhost)
        s["strong_metazoa"] += int(metazoa)
        s["viral_keyword"] += int(viral)
        if nonhost or viral:
            coord = gene_coords.get(row["gene_id"])
            if coord:
                evidence_positions[contig].append((coord[1], coord[2]))

contig_rows = {}
with open(CONTIG_TABLE) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        if row["contig_id"] in candidate_set:
            contig_rows[row["contig_id"]] = row
assert set(contig_rows) == candidate_set

columns = [
    "contig_id", "length_bp", "gc_percent", "tiara_class", "tiara_bac", "tiara_arc", "tiara_euk",
    "fcs_result", "fcs_division", "fcs_best_taxon", "fcs_prok_cov_bp", "fcs_animal_cov_bp", "fcs_virus_cov_bp",
    "illumina_mean_depth", "illumina_breadth_1x_percent", "ont_mean_depth", "ont_breadth_1x_percent",
    "blastn_union_bp", "blastn_union_percent", "blastn_phyla", "blastn_max_identity", "blastn_max_bitscore", "blastn_accessions",
    "genes", "single_exon_genes", "protein_hint_genes", "eggnog_metazoa_genes", "interpro_genes",
    "strong_nonhost_genes", "strong_metazoa_genes", "viral_keyword_genes", "swiss_groups",
    "evidence_span_start", "evidence_span_end", "evidence_span_percent", "sequence_n_bases", "sequence_acgt_entropy",
]
with open(TABLES / "S04_21_precoverage_integrated.tsv", "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for contig in candidates:
        r = contig_rows[contig]
        gs = gene_summary[contig]
        pos = evidence_positions[contig]
        ev_start = min((a for a, _ in pos), default=0)
        ev_end = max((b for _, b in pos), default=0)
        blast_bp = merge_bp(blast_intervals[contig])
        row = {
            "contig_id": contig,
            "length_bp": lengths[contig],
            "gc_percent": f"{seqstats[contig]['sequence_gc_percent']:.4f}",
            "tiara_class": r["tiara_parent_first_stage"],
            "tiara_bac": r["tiara_probability_bac"],
            "tiara_arc": r["tiara_probability_arc"],
            "tiara_euk": r["tiara_probability_euk"],
            "fcs_result": r["fcs_retained_interval_results"],
            "fcs_division": r["fcs_retained_interval_divisions"],
            "fcs_best_taxon": r["fcs_best_taxon_names"],
            "fcs_prok_cov_bp": r["fcs_prokaryotic_division_coverage_bp"],
            "fcs_animal_cov_bp": r["fcs_animal_division_coverage_bp"],
            "fcs_virus_cov_bp": r["fcs_prokaryotic_virus_division_coverage_bp"],
            "illumina_mean_depth": r["Illumina_mean_depth"],
            "illumina_breadth_1x_percent": r["Illumina_percent_ge1x"],
            "ont_mean_depth": r["ONT_mean_depth"],
            "ont_breadth_1x_percent": r["ONT_percent_ge1x"],
            "blastn_union_bp": blast_bp,
            "blastn_union_percent": f"{100 * blast_bp / lengths[contig]:.4f}",
            "blastn_phyla": ";".join(f"{p}:{merge_bp(v)}" for p, v in sorted(blast_by_phylum[contig].items())),
            "blastn_max_identity": f"{blast_max_identity[contig]:.3f}" if blast_rows[contig] else "",
            "blastn_max_bitscore": f"{blast_max_bitscore[contig]:.1f}" if blast_rows[contig] else "",
            "blastn_accessions": ";".join(sorted(blast_accessions[contig])),
            "genes": gs["genes"],
            "single_exon_genes": gs["single_exon"],
            "protein_hint_genes": gs["protein_hint"],
            "eggnog_metazoa_genes": gs["eggnog_metazoa"],
            "interpro_genes": gs["interpro"],
            "strong_nonhost_genes": gs["strong_nonhost"],
            "strong_metazoa_genes": gs["strong_metazoa"],
            "viral_keyword_genes": gs["viral_keyword"],
            "swiss_groups": ";".join(f"{k}:{v}" for k, v in sorted(swiss_groups[contig].items()) if v),
            "evidence_span_start": ev_start or "",
            "evidence_span_end": ev_end or "",
            "evidence_span_percent": f"{100 * (ev_end - ev_start + 1) / lengths[contig]:.4f}" if pos else "0.0000",
            "sequence_n_bases": seqstats[contig]["sequence_n_bases"],
            "sequence_acgt_entropy": f"{seqstats[contig]['sequence_acgt_entropy']:.6f}",
        }
        writer.writerow(row)

with open(ROOT / "source_manifest.tsv", "w") as out:
    out.write("role\tpath\tsha256\n")
    for role, path in [
        ("frozen_assembly", FASTA), ("candidate_ids", REMOVAL_LIST),
        ("existing_blastn_hsps", INPUT / "S04_21_blastn_hsps_annotated.tsv"),
        ("frozen_gff", GFF), ("integrated_contig_evidence", CONTIG_TABLE), ("gene_evidence", GENE_TABLE),
    ]:
        out.write(f"{role}\t{path}\t{sha256(path)}\n")

print(f"PASS\tprepared {len(candidates)} candidates\t{sum(lengths[x] for x in candidates)} bp\tassembly_sha256={actual_sha}")
