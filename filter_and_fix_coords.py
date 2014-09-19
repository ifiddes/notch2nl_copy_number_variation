#!/usr/bin/env/python

"""supposed to be used with full pipeline to filter pysam records
because pysam sucks"""

import pysam, sys, os

header = {"HD": {"VN": "1.3"},
            "SQ": [{"LN": 250522664, "SN": "chr1"}]}

outfile = pysam.Samfile(sys.argv[2], "wb", header=header)

samfile = pysam.Samfile(sys.argv[1], "rb")

for record in samfile:
    if not record.is_unmapped:
        chrom, span = samfile.getrname(record.tid).split(":")
        start, end = map(int, span.split("-"))
        record.pos = record.pos + start - 1
        outfile.write(record)

outfile.close()