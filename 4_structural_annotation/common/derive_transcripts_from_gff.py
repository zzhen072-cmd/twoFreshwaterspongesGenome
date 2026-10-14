#!/usr/bin/env python3
"""Derive transcript FASTA from an audited genome+GFF3 and validate its CDS.

The output is a deterministic gzip stream (mtime=0). Transcript sequence is the
strand-aware concatenation of GFF3 exon features. The same extraction is applied
to CDS features and compared record-by-record with an existing all-isoform CDS
FASTA, which guards against contig, strand, coordinate and annotation-version
mismatches.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import unquote


RC_TABLE = str.maketrans(
    "ACGTRYMKBDHVNacgtrymkbdhvn",
    "TGCAYRKMVHDBNtgcayrkmvhdbn",
)


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_attributes(raw: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for field in raw.strip().split(";"):
        if not field:
            continue
        if "=" not in field:
            raise ValueError(f"Malformed GFF3 attribute: {field!r}")
        key, value = field.split("=", 1)
        parsed[unquote(key)] = unquote(value)
    return parsed


def build_fasta_index(path: Path) -> dict[str, dict[str, int]]:
    index: dict[str, dict[str, int]] = {}
    current = None
    current_id = None
    saw_short_line = False
    with path.open("rb") as handle:
        while True:
            line_start = handle.tell()
            line = handle.readline()
            if not line:
                break
            if line.startswith(b">"):
                header_id = line[1:].strip().split(None, 1)[0].decode("utf-8")
                if not header_id:
                    raise ValueError(f"Empty FASTA identifier at byte {line_start}")
                if header_id in index:
                    raise ValueError(f"Duplicate genome FASTA identifier: {header_id}")
                current_id = header_id
                current = {
                    "length": 0,
                    "offset": handle.tell(),
                    "line_bases": 0,
                    "line_bytes": 0,
                }
                index[header_id] = current
                saw_short_line = False
                continue
            if current is None:
                raise ValueError("Sequence data found before first FASTA header")
            stripped = line.rstrip(b"\r\n")
            if not stripped:
                continue
            bases = len(stripped)
            if current["line_bases"] == 0:
                current["line_bases"] = bases
                current["line_bytes"] = len(line)
            else:
                if saw_short_line:
                    raise ValueError(f"Nonterminal short FASTA line for {current_id}")
                if bases > current["line_bases"]:
                    raise ValueError(f"Inconsistent FASTA wrapping for {current_id}")
                if bases < current["line_bases"]:
                    saw_short_line = True
            current["length"] += bases
    if not index:
        raise ValueError("Genome FASTA is empty")
    return index


class IndexedFasta:
    def __init__(self, path: Path):
        self.path = path
        self.index = build_fasta_index(path)
        self.handle = path.open("rb")

    def close(self):
        self.handle.close()

    def fetch(self, seqid: str, start: int, end: int) -> str:
        if seqid not in self.index:
            raise KeyError(f"GFF3 seqid absent from genome: {seqid}")
        entry = self.index[seqid]
        if start < 1 or end < start or end > entry["length"]:
            raise ValueError(
                f"Out-of-range interval {seqid}:{start}-{end}; length={entry['length']}"
            )
        position = start - 1
        remaining = end - start + 1
        chunks: list[bytes] = []
        while remaining:
            line_index, within = divmod(position, entry["line_bases"])
            take = min(remaining, entry["line_bases"] - within)
            byte_offset = entry["offset"] + line_index * entry["line_bytes"] + within
            self.handle.seek(byte_offset)
            chunk = self.handle.read(take)
            if len(chunk) != take or b"\n" in chunk or b"\r" in chunk:
                raise ValueError(f"FASTA random-access failure at {seqid}:{position + 1}")
            chunks.append(chunk)
            position += take
            remaining -= take
        return b"".join(chunks).decode("ascii").upper()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(RC_TABLE)[::-1]


def parse_gff(path: Path):
    transcripts: dict[str, dict] = {}
    order: list[str] = []
    feature_counts: dict[str, int] = {}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GFF3 line {line_number} has {len(fields)} columns")
            seqid, _source, feature, start_raw, end_raw, _score, strand, phase_raw, attrs_raw = fields
            feature_counts[feature] = feature_counts.get(feature, 0) + 1
            attrs = parse_attributes(attrs_raw)
            start, end = int(start_raw), int(end_raw)
            if feature in {"mRNA", "transcript"}:
                transcript_id = attrs.get("ID")
                if not transcript_id:
                    raise ValueError(f"GFF3 transcript without ID at line {line_number}")
                if transcript_id in transcripts and transcripts[transcript_id].get("defined"):
                    raise ValueError(f"Duplicate transcript ID: {transcript_id}")
                record = transcripts.setdefault(
                    transcript_id, {"exons": [], "cds": [], "defined": False}
                )
                record.update(
                    {
                        "defined": True,
                        "seqid": seqid,
                        "strand": strand,
                        "start": start,
                        "end": end,
                        "gene_id": attrs.get("Parent", ""),
                    }
                )
                order.append(transcript_id)
            elif feature in {"exon", "CDS"}:
                parents = attrs.get("Parent", "")
                if not parents:
                    raise ValueError(f"GFF3 {feature} without Parent at line {line_number}")
                for transcript_id in parents.split(","):
                    record = transcripts.setdefault(
                        transcript_id, {"exons": [], "cds": [], "defined": False}
                    )
                    record["exons" if feature == "exon" else "cds"].append(
                        (start, end, seqid, strand, phase_raw)
                    )

    for transcript_id, record in transcripts.items():
        if not record.get("defined"):
            raise ValueError(f"Feature Parent lacks transcript record: {transcript_id}")
        if record["strand"] not in {"+", "-"}:
            raise ValueError(f"Unsupported strand for {transcript_id}: {record['strand']}")
        if not record["exons"]:
            raise ValueError(f"Transcript lacks exon features: {transcript_id}")
        if not record["cds"]:
            raise ValueError(f"Transcript lacks CDS features: {transcript_id}")
        for feature_name in ("exons", "cds"):
            intervals = record[feature_name]
            for start, end, seqid, strand, _phase in intervals:
                if seqid != record["seqid"] or strand != record["strand"]:
                    raise ValueError(
                        f"{feature_name} seqid/strand mismatch for {transcript_id}"
                    )
                if start < record["start"] or end > record["end"]:
                    raise ValueError(
                        f"{feature_name} outside transcript bounds for {transcript_id}"
                    )
            ordered = sorted(intervals)
            for left, right in zip(ordered, ordered[1:]):
                if right[0] <= left[1]:
                    raise ValueError(f"Overlapping {feature_name} for {transcript_id}")
    if len(order) != len(transcripts):
        raise ValueError("Transcript order/count mismatch")
    return transcripts, order, feature_counts


def parse_fasta_digests(path: Path):
    digests: dict[str, str] = {}
    current_id = None
    digest = None
    count = 0
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if current_id is not None:
                    digests[current_id] = digest.hexdigest()
                current_id = line[1:].strip().split(None, 1)[0]
                if not current_id or current_id in digests:
                    raise ValueError(f"Empty or duplicate FASTA ID in {path}: {current_id}")
                digest = hashlib.sha256()
                count += 1
            else:
                if digest is None:
                    raise ValueError(f"Sequence before header in {path}")
                digest.update(line.strip().upper().encode("ascii"))
        if current_id is not None:
            digests[current_id] = digest.hexdigest()
    return digests, count


def write_wrapped(handle, sequence: str, width: int = 60):
    for offset in range(0, len(sequence), width):
        handle.write(sequence[offset : offset + width] + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--cds", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--expected-transcripts", type=int, required=True)
    args = parser.parse_args()

    transcripts, order, feature_counts = parse_gff(args.gff)
    if len(order) != args.expected_transcripts:
        raise ValueError(
            f"Expected {args.expected_transcripts} transcripts, found {len(order)}"
        )
    expected_cds, expected_cds_count = parse_fasta_digests(args.cds)
    if expected_cds_count != args.expected_transcripts:
        raise ValueError(
            f"Expected {args.expected_transcripts} CDS records, found {expected_cds_count}"
        )
    if set(expected_cds) != set(transcripts):
        missing = sorted(set(transcripts) - set(expected_cds))[:10]
        extra = sorted(set(expected_cds) - set(transcripts))[:10]
        raise ValueError(f"GFF/CDS ID-set mismatch; missing={missing}, extra={extra}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    genome = IndexedFasta(args.genome)
    generated_cds: dict[str, str] = {}
    terminal_stop_codons_omitted_by_audited_cds = 0
    transcripts_with_leading_n_padding = 0
    transcripts_with_trailing_n_padding = 0
    transcript_lengths: list[int] = []
    uncompressed_hash = hashlib.sha256()
    try:
        with args.output.open("wb") as raw_output:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw_output, mtime=0) as gzip_output:
                with io.TextIOWrapper(gzip_output, encoding="ascii", newline="\n") as output:
                    class HashingWriter:
                        def write(self, text: str):
                            uncompressed_hash.update(text.encode("ascii"))
                            output.write(text)

                    writer = HashingWriter()
                    for transcript_id in order:
                        record = transcripts[transcript_id]
                        exon_sequence = "".join(
                            genome.fetch(record["seqid"], start, end)
                            for start, end, _seqid, _strand, _phase in sorted(record["exons"])
                        )
                        cds_sequence = "".join(
                            genome.fetch(record["seqid"], start, end)
                            for start, end, _seqid, _strand, _phase in sorted(record["cds"])
                        )
                        if record["strand"] == "-":
                            exon_sequence = reverse_complement(exon_sequence)
                            cds_sequence = reverse_complement(cds_sequence)
                        cds_in_transcription_order = sorted(
                            record["cds"],
                            key=lambda interval: interval[0],
                            reverse=record["strand"] == "-",
                        )
                        first_phase_raw = cds_in_transcription_order[0][4]
                        if first_phase_raw not in {"0", "1", "2"}:
                            raise ValueError(
                                f"Invalid first CDS phase for {transcript_id}: {first_phase_raw}"
                            )
                        first_phase = int(first_phase_raw)
                        leading_padding = (3 - first_phase) % 3
                        exported_cds = "N" * leading_padding + cds_sequence
                        trailing_padding = (-len(exported_cds)) % 3
                        exported_cds += "N" * trailing_padding
                        if leading_padding:
                            transcripts_with_leading_n_padding += 1
                        if trailing_padding:
                            transcripts_with_trailing_n_padding += 1
                        if exported_cds[-3:] in {"TAA", "TAG", "TGA"}:
                            exported_cds = exported_cds[:-3]
                            terminal_stop_codons_omitted_by_audited_cds += 1
                        cds_hash = hashlib.sha256(exported_cds.encode("ascii")).hexdigest()
                        if cds_hash != expected_cds[transcript_id]:
                            raise ValueError(
                                "Extracted CDS differs from audited CDS FASTA under the "
                                f"documented BRAKER export convention: {transcript_id}"
                            )
                        generated_cds[transcript_id] = cds_hash
                        transcript_lengths.append(len(exon_sequence))
                        header = (
                            f">{transcript_id} gene={record['gene_id']} "
                            f"seqid={record['seqid']} strand={record['strand']}\n"
                        )
                        writer.write(header)
                        for offset in range(0, len(exon_sequence), 60):
                            writer.write(exon_sequence[offset : offset + 60] + "\n")
    finally:
        genome.close()

    if generated_cds != expected_cds:
        raise ValueError("Generated and audited CDS digest maps are not identical")
    report = {
        "schema_version": "1.0",
        "method": "strand-aware concatenation of final GFF3 exon features",
        "genome_path": str(args.genome),
        "genome_sha256": sha256_file(args.genome),
        "gff3_path": str(args.gff),
        "gff3_sha256": sha256_file(args.gff),
        "audited_cds_path": str(args.cds),
        "audited_cds_sha256": sha256_file(args.cds),
        "output_path": str(args.output),
        "output_gzip_sha256": sha256_file(args.output),
        "output_uncompressed_sha256": uncompressed_hash.hexdigest(),
        "transcript_records": len(order),
        "unique_transcript_ids": len(set(order)),
        "gff3_mrna_records": feature_counts.get("mRNA", 0)
        + feature_counts.get("transcript", 0),
        "gff3_exon_records": feature_counts.get("exon", 0),
        "gff3_cds_records": feature_counts.get("CDS", 0),
        "cds_record_id_set_matches_gff3": True,
        "all_extracted_cds_sequences_match_audited_cds_fasta": True,
        "transcripts_with_leading_n_padding_in_audited_cds_fasta": transcripts_with_leading_n_padding,
        "transcripts_with_trailing_n_padding_in_audited_cds_fasta": transcripts_with_trailing_n_padding,
        "terminal_stop_codons_omitted_by_audited_cds_fasta": terminal_stop_codons_omitted_by_audited_cds,
        "transcript_length_min": min(transcript_lengths),
        "transcript_length_max": max(transcript_lengths),
        "transcript_length_sum": sum(transcript_lengths),
        "transcript_length_mean": sum(transcript_lengths) / len(transcript_lengths),
        "status": "PASS",
    }
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
