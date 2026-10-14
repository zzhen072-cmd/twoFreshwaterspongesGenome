#!/usr/bin/env python3
"""Reproduce existing Swiss-Prot best-hit TSV: first DIAMOND output row/query.
Usage: python3 first_hit.py input.tsv > output.best.tsv
"""
import sys
seen=set()
with open(sys.argv[1]) as source:
    for line in source:
        query=line.split('\t',1)[0]
        if query not in seen:
            seen.add(query)
            sys.stdout.write(line)
