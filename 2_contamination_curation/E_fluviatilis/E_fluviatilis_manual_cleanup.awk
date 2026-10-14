BEGIN {
    FS = "\t"
    OFS = "\t"
    print "seq_id", "fcs_action", "manual_action", "coordinates", "reason" > decisions_out
}

# Pass 1: read the FCS action report.
FILENAME == action_file {
    if ($0 ~ /^#/) next
    if ($5 == "REVIEW") review[$1] = 1
    next
}

# Pass 2: stream Cleaned.fasta.
FILENAME == cleaned_file {
    if ($0 ~ /^>/) {
        flush_cleaned()
        header = $0
        split(substr($0, 2), fields, /[[:space:]]+/)
        id = fields[1]
        seq = ""
        next
    }
    seq = seq $0
    next
}

# Pass 3: stream the FCS EXCLUDE fasta and recover the eukaryotic part of
# tig00018876 only.
FILENAME == excluded_file {
    if (!cleaned_flushed) {
        flush_cleaned()
        cleaned_flushed = 1
    }
    if ($0 ~ /^>/) {
        flush_excluded()
        excluded_header = $0
        split(substr($0, 2), fields, /[[:space:]]+/)
        excluded_id = fields[1]
        excluded_seq = ""
        next
    }
    excluded_seq = excluded_seq $0
    next
}

END {
    flush_excluded()
}

function write_fasta(out, h, s,    i) {
    print h > out
    for (i = 1; i <= length(s); i += 80)
        print substr(s, i, 80) > out
}

function flush_cleaned(    keepseq) {
    if (header == "") return

    if (id == "lcl|tig00003180_pilon_1_np1212") {
        write_fasta(review_out, header, seq)
        keepseq = substr(seq, 228692)
        write_fasta(clean_out,
            ">lcl|tig00003180_pilon_1_np1212_trimmed_228692_280325 length=51634 source=S04_FCS_REVIEW",
            keepseq)
        print id, "REVIEW", "SPLIT", "remove:1-228691;keep:228692-280325", "FCS bacterial first block; Tiara eukarya" > decisions_out
    } else if (id == "lcl|tig00016872_pilon_1_np1212") {
        write_fasta(review_out, header, seq)
        print id, "REVIEW", "HOLD_OUT", "entire_contig", "Tiara plastid; not placed in nuclear assembly" > decisions_out
    } else if (id in review) {
        write_fasta(removed_out, header, seq)
        print id, "REVIEW", "REMOVE", "entire_contig", "FCS prokaryote/virus review candidate" > decisions_out
    } else {
        write_fasta(clean_out, header, seq)
    }

    header = ""
    seq = ""
}

function flush_excluded(    keepseq) {
    if (excluded_header == "") return
    if (excluded_id == "lcl|tig00018876_pilon_1_np1212") {
        write_fasta(review_out, excluded_header, excluded_seq)
        keepseq = substr(excluded_seq, 122176)
        write_fasta(clean_out,
            ">lcl|tig00018876_pilon_1_np1212_recovered_122176_420595 length=298420 source=S04_FCS_EXCLUDE",
            keepseq)
        print excluded_id, "EXCLUDE", "RECOVER_SPLIT", "remove:1-122175;keep:122176-420595", "FCS/Taxonomy report supports bacterial-eukaryotic chimera" > decisions_out
    }
    excluded_header = ""
    excluded_seq = ""
}
