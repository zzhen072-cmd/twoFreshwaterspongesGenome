#!/bin/sh
set -eu

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 INPUT_GX_DIR OUTPUT_DIR" >&2
    exit 2
fi

gx_dir=$1
out_dir=$2
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

cleaned="$gx_dir/Cleaned.fasta"
actions="$gx_dir/Galaxy20-[NCBI FCS GX on dataset 18_ Action report].tabular"
excluded="$gx_dir/Galaxy21-[NCBI FCS GX on dataset 18 and 20_ Fasta for EXCLUDE entries].fasta"

for input in "$cleaned" "$actions" "$excluded"; do
    test -s "$input" || { echo "Missing input: $input" >&2; exit 1; }
done

mkdir -p "$out_dir"

clean_out="$out_dir/S04.manual_clean.candidate.fasta"
removed_out="$out_dir/S04.removed_review.fasta"
review_out="$out_dir/S04.manual_review.fasta"
decisions_out="$out_dir/S04.cleaning_decisions.tsv"

: > "$clean_out"
: > "$removed_out"
: > "$review_out"
: > "$decisions_out"

awk \
    -v clean_out="$clean_out" \
    -v removed_out="$removed_out" \
    -v review_out="$review_out" \
    -v decisions_out="$decisions_out" \
    -v action_file="$actions" \
    -v cleaned_file="$cleaned" \
    -v excluded_file="$excluded" \
    -f "$script_dir/E_fluviatilis_manual_cleanup.awk" \
    "$actions" "$cleaned" "$excluded"

for output in "$clean_out" "$removed_out" "$review_out" "$decisions_out"; do
    test -s "$output" || { echo "Empty output: $output" >&2; exit 1; }
done

awk '/^>/{n++;next}{bp+=length($0)}END{print "candidate_seqs\t" n "\ncandidate_bp\t" bp}' "$clean_out" > "$out_dir/S04.manual_clean.stats.tsv"
awk '/^>/{n++;next}{bp+=length($0)}END{print "removed_seqs\t" n "\nremoved_bp\t" bp}' "$removed_out" >> "$out_dir/S04.manual_clean.stats.tsv"
awk '/^>/{n++;next}{bp+=length($0)}END{print "review_seqs\t" n "\nreview_bp\t" bp}' "$review_out" >> "$out_dir/S04.manual_clean.stats.tsv"

echo "Created outputs in $out_dir"
