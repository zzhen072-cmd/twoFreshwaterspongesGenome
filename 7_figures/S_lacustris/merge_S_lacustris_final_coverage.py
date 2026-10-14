#!/usr/bin/env python3
"""Merge final-reference mosdepth coverage into an audited S. lacustris plotting table.

This script is specific to the verified 1,411-contig GZXS202303 final assembly. It
validates its file MD5, count and span, and requires an exact ID/length match
between the FASTA, existing plotting table and mosdepth summary. The source
table and existing outputs are never overwritten. No alignments are performed.

Example:
  python3 merge_S_lacustris_final_coverage.py \
    --table output/gzx_blobplot/gzx_final_plot_data.tsv \
    --summary gzx.final.ONT.mosdepth.summary.txt \
    --fasta gzx_cleaned.fasta \
    --output-prefix output/S_lacustris_blobplot_final_coverage/S_lacustris_final_plot_data

Outputs: <prefix>.tsv and <prefix>.audit.json. Python 3.9+, standard library only.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import json
import math
import sys
from pathlib import Path


EXPECTED_MD5 = "097a27663506e2db23fb61c67e9825a6"
EXPECTED_COUNT = 1411
EXPECTED_SPAN = 264956519
ADDED_FIELDS = (
    "previous_ont_coverage", "previous_coverage_reference",
    "previous_coverage_status", "coverage_reference_md5", "coverage_method",
    "mosdepth_bases", "mosdepth_reported_mean",
)


def fail(message):
    raise ValueError(message)


def nonnegative_integer(value, context):
    try:
        result = int(value)
    except (ValueError, TypeError):
        fail(f"{context}: expected an integer, got {value!r}")
    if result < 0:
        fail(f"{context}: negative value")
    return result


def nonnegative_number(value, context):
    try:
        result = float(value)
    except (ValueError, TypeError):
        fail(f"{context}: expected a number, got {value!r}")
    if not math.isfinite(result) or result < 0:
        fail(f"{context}: must be finite and nonnegative")
    return result


def read_final_fasta(path):
    """Read lengths and hash the original FASTA bytes in one pass."""
    lengths = {}
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    current = None
    with path.open("rb") as handle:
        for line_number, raw in enumerate(handle, start=1):
            md5.update(raw)
            sha256.update(raw)
            line = raw.strip()
            if not line:
                continue
            if line.startswith(b">"):
                fields = line[1:].split()
                if not fields:
                    fail(f"FASTA line {line_number}: empty identifier")
                current = fields[0].decode("ascii")
                if current in lengths:
                    fail(f"FASTA line {line_number}: duplicate ID {current}")
                lengths[current] = 0
            else:
                if current is None:
                    fail(f"FASTA line {line_number}: sequence before first header")
                if len(line.split()) != 1:
                    fail(f"FASTA line {line_number}: whitespace within sequence")
                lengths[current] += len(line)
    observed_md5 = md5.hexdigest()
    if observed_md5 != EXPECTED_MD5:
        fail(f"Final FASTA MD5 mismatch: expected {EXPECTED_MD5}, got {observed_md5}")
    if len(lengths) != EXPECTED_COUNT or sum(lengths.values()) != EXPECTED_SPAN:
        fail(f"Final FASTA count/span mismatch: {len(lengths)} contigs, "
             f"{sum(lengths.values())} bp")
    if any(length <= 0 for length in lengths.values()):
        fail("Final FASTA contains an empty sequence")
    return lengths, {"md5": observed_md5, "sha256": sha256.hexdigest()}


def check_ids(observed, expected, source):
    missing = sorted(set(expected) - set(observed))
    extra = sorted(set(observed) - set(expected))
    if missing or extra:
        fail(f"{source}: ID mismatch; missing={len(missing)} {missing[:5]}, "
             f"extra={len(extra)} {extra[:5]}")


def read_plot_table(path, lengths):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        required = {"contig_id", "length_bp", "gc_fraction", "ont_coverage", "phylum"}
        if len(fields) != len(set(fields)) or not required.issubset(fields):
            fail("Plot table has duplicate headers or is missing required plotting columns")
        if "previous_ont_coverage" in fields:
            fail("Plot table already contains previous_ont_coverage; use the original source table")
        rows = []
        seen = set()
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                fail(f"Plot table line {line_number}: inconsistent column count")
            ident = row["contig_id"]
            if not ident or ident in seen:
                fail(f"Plot table line {line_number}: empty or duplicate ID {ident!r}")
            seen.add(ident)
            length = nonnegative_integer(row["length_bp"], f"Plot table {ident} length")
            if ident in lengths and length != lengths[ident]:
                fail(f"Plot table {ident}: length {length} differs from FASTA {lengths[ident]}")
            nonnegative_number(row["ont_coverage"], f"Plot table {ident} old coverage")
            gc = nonnegative_number(row["gc_fraction"], f"Plot table {ident} GC")
            if gc > 1:
                fail(f"Plot table {ident}: GC exceeds 1")
            if not row["phylum"].strip():
                fail(f"Plot table {ident}: blank phylum")
            rows.append(row)
    check_ids(seen, lengths, "Plot table")
    return fields, rows


def read_mosdepth(path, lengths):
    """Read standard whole-contig summary; never infer zero from a missing row."""
    records = {}
    total = None
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or not {"chrom", "length", "bases", "mean"}.issubset(fields):
            fail("mosdepth summary requires unique chrom, length, bases and mean headers")
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                fail(f"mosdepth line {line_number}: inconsistent column count")
            ident = row["chrom"]
            length = nonnegative_integer(row["length"], f"mosdepth {ident} length")
            bases = nonnegative_integer(row["bases"], f"mosdepth {ident} bases")
            mean = nonnegative_number(row["mean"], f"mosdepth {ident} mean")
            record = {"length": length, "bases": bases, "reported_mean": row["mean"],
                      "reported_mean_numeric": mean}
            if ident == "total":
                if total is not None:
                    fail("mosdepth summary has duplicate total rows")
                total = record
                continue
            if not ident or ident in records:
                fail(f"mosdepth line {line_number}: empty or duplicate ID {ident!r}")
            if ident in lengths and length != lengths[ident]:
                fail(f"mosdepth {ident}: length {length} differs from FASTA {lengths[ident]}")
            records[ident] = record
    check_ids(records, lengths, "mosdepth summary")
    observed_span = sum(record["length"] for record in records.values())
    observed_bases = sum(record["bases"] for record in records.values())
    if total is not None and (total["length"] != observed_span or total["bases"] != observed_bases):
        fail("mosdepth total length/bases does not equal the sum of contig rows")
    return records, total


def fingerprint(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {"path": str(path.resolve()), "bytes": path.stat().st_size,
            "sha256": digest.hexdigest()}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args(argv)
    output_table = Path(str(args.output_prefix) + ".tsv")
    output_audit = Path(str(args.output_prefix) + ".audit.json")
    try:
        for path in (output_table, output_audit):
            if path.exists() or path.is_symlink():
                fail(f"Refusing to overwrite existing output: {path}")
        lengths, fasta_hashes = read_final_fasta(args.fasta)
        fields, rows = read_plot_table(args.table, lengths)
        coverage, total = read_mosdepth(args.summary, lengths)
        for field in (*ADDED_FIELDS, "coverage_reference", "coverage_status", "assembly_scope"):
            if field not in fields:
                fields.append(field)
        zero_ids = []
        for row in rows:
            ident = row["contig_id"]
            record = coverage[ident]
            row["previous_ont_coverage"] = row["ont_coverage"]
            row["previous_coverage_reference"] = row.get("coverage_reference", "")
            row["previous_coverage_status"] = row.get("coverage_status", "")
            row["ont_coverage"] = repr(record["bases"] / record["length"])
            row["coverage_reference"] = args.fasta.name
            row["coverage_reference_md5"] = fasta_hashes["md5"]
            row["coverage_status"] = "computed_on_final_reference"
            row["assembly_scope"] = "final_clean_assembly_with_final_reference_ONT_coverage"
            row["coverage_method"] = "mosdepth whole-contig depth: bases / length"
            row["mosdepth_bases"] = record["bases"]
            row["mosdepth_reported_mean"] = record["reported_mean"]
            if record["bases"] == 0:
                zero_ids.append(ident)
        audit = {
            "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "script": Path(__file__).name,
            "assembly": {"path": str(args.fasta.resolve()), **fasta_hashes,
                         "contig_count": len(lengths), "span_bp": sum(lengths.values())},
            "inputs": {"plot_table": fingerprint(args.table), "mosdepth_summary": fingerprint(args.summary)},
            "validation": {"expected_fasta_md5": EXPECTED_MD5,
                           "expected_count": EXPECTED_COUNT, "expected_span_bp": EXPECTED_SPAN,
                           "all_three_id_sets_and_lengths_match": True,
                           "mosdepth_total_row_present": total is not None},
            "coverage": {"method": "bases / length, without rounding the reported summary mean",
                         "denominator": "full contig length including N positions",
                         "total_depth_bases": sum(r["bases"] for r in coverage.values()),
                         "zero_coverage_count": len(zero_ids),
                         "zero_coverage_span_bp": sum(lengths[ident] for ident in zero_ids),
                         "zero_coverage_ids": zero_ids,
                         "min": min(float(r["ont_coverage"]) for r in rows),
                         "max": max(float(r["ont_coverage"]) for r in rows),
                         "previous_depth_column": "previous_ont_coverage"},
            "provenance_note": "The cluster mapping run must establish that the BAM was aligned to this exact FASTA. A mosdepth summary contains IDs and lengths, but not reference sequence checksums or mapping parameters.",
            "plot_note": "Zero depth is retained as 0. The existing log-scale plotting script excludes zero-depth points and reports their count/span; missing rows are rejected rather than converted to zero.",
            "outputs": {"table": str(output_table.resolve()), "audit": str(output_audit.resolve())},
        }
        output_table.parent.mkdir(parents=True, exist_ok=True)
        # Exclusive creation is a second safeguard against output-name collisions.
        with output_table.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        with output_audit.open("x", encoding="utf-8") as handle:
            json.dump(audit, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        print(f"Merged {len(rows):,} contigs; span {sum(lengths.values()):,} bp.")
        print(f"Zero coverage: {len(zero_ids):,} contigs; {sum(lengths[q] for q in zero_ids):,} bp.")
        print(f"Table: {output_table}\nAudit: {output_audit}")
    except (OSError, ValueError, UnicodeError, csv.Error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
