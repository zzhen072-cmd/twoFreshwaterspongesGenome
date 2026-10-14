#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 4 ]; then
  echo "Usage: $0 CLEANED_FASTA FCS_GX_ACTION_REPORT TIARA_TSV OUTPUT_DIR" >&2
  exit 2
fi

cleaned=$1
action=$2
tiara=$3
out_dir=$4

for input in "$cleaned" "$action" "$tiara"; do
  test -s "$input" || { echo "Missing input: $input" >&2; exit 1; }
done

mkdir -p "$out_dir"

awk -F '\t' 'BEGIN{OFS="\t"}
  NR==FNR {
    if ($0 !~ /^#/ && $5=="REVIEW") {
      id=$1; sub(/^lcl\|/, "", id)
      len[id]=$4; div[id]=$6; cov[id]=$7; tax[id]=$8
    }
    next
  }
  NR>1 {
    split($1,a," "); id=a[1]
    if ((id in len) && $2=="bacteria")
      print "lcl|" id,len[id],div[id],cov[id],tax[id],$2,$4,$5,$6,$7
  }' "$action" "$tiara" | sort -t $'\t' -k2,2nr > "$out_dir/removed_review_bacteria.tsv.tmp"

{
  printf 'seq_id\tlength\tfcs_div\tfcs_agg_cont_cov\tfcs_top_tax_name\ttiara_class\ttiara_org\ttiara_bac\ttiara_arc\ttiara_euk\n'
  awk '{print}' "$out_dir/removed_review_bacteria.tsv.tmp"
} > "$out_dir/removed_review_bacteria.tsv"
rm "$out_dir/removed_review_bacteria.tsv.tmp"

awk -F '\t' 'NR==FNR{if (FNR>1) drop[$1]=1; next}
  /^>/ {
    id=substr($0,2); sub(/[[:space:]].*$/, "", id)
    keep=!(id in drop)
  }
  keep {print}' "$out_dir/removed_review_bacteria.tsv" "$cleaned" \
  > "$out_dir/gzx.Cleaned.review_filtered.fasta"

awk -F '\t' 'NR==FNR{if (FNR>1) drop[$1]=1; next}
  /^>/ {
    id=substr($0,2); sub(/[[:space:]].*$/, "", id)
    emit=(id in drop)
  }
  emit {print}' "$out_dir/removed_review_bacteria.tsv" "$cleaned" \
  > "$out_dir/gzx.removed_review_bacteria.fasta"

for fasta in "$cleaned" "$out_dir/gzx.Cleaned.review_filtered.fasta" "$out_dir/gzx.removed_review_bacteria.fasta"; do
  awk -v file="$fasta" '
    /^>/ {if (len) {n++; bp+=len}; len=0; next}
    {gsub(/[[:space:]]/,""); len+=length($0)}
    END {if (len) {n++; bp+=len}; printf "%s\t%d\t%d\n",file,n,bp}' "$fasta"
done > "$out_dir/fasta_stats.tsv"

printf 'Created %s\n' "$out_dir/gzx.Cleaned.review_filtered.fasta"
