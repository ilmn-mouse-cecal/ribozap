#!/usr/bin/env python3

from Bio import SeqIO
import csv
import sys

fasta_path = sys.argv[1]
output_csv = sys.argv[2]

with open(output_csv, "w", newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Probe ID", "Sequence", "Length"])

    for record in SeqIO.parse(fasta_path, "fasta"):
        writer.writerow([record.id, str(record.seq), len(record.seq)])
