#!/usr/bin/env python3
"""Independently validate E. fluviatilis blob-plot tables and rendered deliverables."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import subprocess
from collections import defaultdict
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--pdf-text-audit", type=Path, required=True)
    args = parser.parse_args()
    root = args.root

    source = read_tsv(root / "E_fluviatilis_final_plot_data.tsv")
    points = read_tsv(root / "E_fluviatilis_final_blobplot.plotted_points.tsv")
    summary = read_tsv(root / "E_fluviatilis_final_blobplot.category_summary.tsv")
    bins = read_tsv(root / "E_fluviatilis_final_blobplot.histogram_bins.tsv")
    source_by_id = {row["contig_id"]: row for row in source}
    points_by_id = {row["contig_id"]: row for row in points}
    checks = []

    def check(name: str, passed: bool, **details) -> None:
        checks.append({"name": name, "passed": bool(passed), **details})

    check(
        "unique_and_identical_contig_ids",
        len(source_by_id) == len(source)
        and len(points_by_id) == len(points)
        and set(source_by_id) == set(points_by_id),
        source_count=len(source),
        point_count=len(points),
    )
    fields_match = all(
        source_by_id[ident][field] == points_by_id[ident][field]
        for ident in source_by_id
        for field in ("length_bp", "gc_fraction", "ont_coverage", "phylum")
    )
    check("plotted_fields_match_source_table", fields_match)

    positive = [row for row in source if float(row["ont_coverage"]) > 0]
    zero = [row for row in source if float(row["ont_coverage"]) == 0]
    check(
        "log_axis_exclusion_is_explicit",
        sum(row["plotted"] == "True" for row in points) == len(positive)
        and sum(row["plotted"] == "False" for row in points) == len(zero),
        plotted_contigs=len(positive),
        zero_coverage_contigs=len(zero),
        zero_coverage_length_bp=sum(int(row["length_bp"]) for row in zero),
    )

    positive_length = sum(int(row["length_bp"]) for row in positive)
    histogram_totals = defaultdict(int)
    for row in bins:
        histogram_totals[row["axis"]] += int(row["sum_length_bp"])
    check(
        "marginal_histograms_preserve_positive_coverage_length",
        histogram_totals["gc_fraction"] == positive_length
        and histogram_totals["ont_coverage"] == positive_length,
        expected_length_bp=positive_length,
        gc_histogram_length_bp=histogram_totals["gc_fraction"],
        coverage_histogram_length_bp=histogram_totals["ont_coverage"],
    )

    total_row = next(row for row in summary if row["plot_group"] == "Total")
    total_length = sum(int(row["length_bp"]) for row in source)
    check(
        "category_summary_total",
        int(total_row["contig_count"]) == len(source)
        and int(total_row["sum_length_bp"]) == total_length
        and int(total_row["plotted_contig_count"]) == len(positive),
        contig_count=len(source),
        sum_length_bp=total_length,
    )

    collision = json.loads(
        (root / "E_fluviatilis_final_blobplot.collision-audit.json").read_text(
            encoding="utf-8"
        )
    )
    check(
        "rendered_collision_audit",
        collision.get("verdict") == "PASS"
        and collision.get("summary", {}).get("fail") == 0
        and collision.get("summary", {}).get("warn") == 0,
        verdict=collision.get("verdict"),
        summary=collision.get("summary"),
    )
    alignment = json.loads(
        (root / "E_fluviatilis_final_blobplot.alignment.json").read_text(encoding="utf-8")
    )
    check(
        "single_composite_panel_alignment_gate",
        alignment.get("verdict") == "NOT APPLICABLE"
        and len(alignment.get("layout", {}).get("panels", [])) == 1,
        verdict=alignment.get("verdict"),
        rationale=(
            "The scientific display is one blob plot; top/right axes are linked marginals "
            "and the fourth axis is a legend table, not comparable panels."
        ),
    )

    text_result = subprocess.run(
        [
            "python",
            str(args.pdf_text_audit),
            str(root / "E_fluviatilis_final_blobplot.pdf"),
            "--min-pt",
            "5",
            "--json",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    pdf_text = json.loads(text_result.stdout)
    check(
        "pdf_text_glyph_floor",
        pdf_text.get("auditable")
        and pdf_text.get("below_minimum_count") == 0
        and pdf_text.get("minimum_found_pt", 0) >= 5,
        minimum_found_pt=pdf_text.get("minimum_found_pt"),
        text_run_count=pdf_text.get("text_run_count"),
    )

    output_files = {}
    for name in (
        "E_fluviatilis_final_blobplot.png",
        "E_fluviatilis_final_blobplot.svg",
        "E_fluviatilis_final_blobplot.pdf",
        "E_fluviatilis_final_blobplot.category_summary.tsv",
        "E_fluviatilis_final_blobplot.histogram_bins.tsv",
        "E_fluviatilis_final_blobplot.plotted_points.tsv",
    ):
        path = root / name
        output_files[name] = {
            "exists": path.exists(),
            "bytes": path.stat().st_size if path.exists() else 0,
            "sha256": sha256(path) if path.exists() else "",
        }
    check("all_delivery_outputs_exist", all(item["exists"] for item in output_files.values()))

    report = {
        "status": "passed" if all(item["passed"] for item in checks) else "failed",
        "checked_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "scope": "Independent numeric and rendered-output QA; source data were not modified.",
        "figure_contract": {
            "core_conclusion": (
                "The final Ephydatia fluviatilis (NJ151201) assembly is dominated by "
                "Porifera-assigned sequence, while "
                "minor taxonomic groups and no-hit contigs occupy distinct GC/depth positions."
            ),
            "archetype": "asymmetric mixed-modality figure",
            "backend": "Python/Matplotlib",
            "display_role": "assembly contamination and taxonomic-QC visualization",
        },
        "counts": {
            "all_contigs": len(source),
            "plotted_positive_coverage_contigs": len(positive),
            "zero_coverage_contigs": len(zero),
            "assembly_span_bp": total_length,
        },
        "output_files": output_files,
        "checks": checks,
    }
    (root / "plot_qa.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
