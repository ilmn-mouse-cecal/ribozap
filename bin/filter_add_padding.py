#!/usr/bin/env python3
import argparse
from Bio import SeqIO

# Argument parsing
parser = argparse.ArgumentParser(description="Extract antisense hits from blast result table and add padding.")
parser.add_argument("-i", "--input", required=True, help="Input alignment TSV file")
parser.add_argument("-o", "--output", required=True, help="Output TSV file")
parser.add_argument("-p", "--padding", type=int, default=50, help="Padding around aligned region (default: 50)")
parser.add_argument("-f", "--fasta", required=True, help="Reference fasta file")

args = parser.parse_args()

# Load reference fasta
fasta_dict = SeqIO.index(args.fasta, "fasta")
print(f"Loaded {len(fasta_dict)} sequences from {args.fasta}")

with open(args.input, "r") as fin, open(args.output, "w") as fout:
    for line in fin:
        fields = line.strip().split('\t')

        # Check antisense alignment
        if int(fields[8]) > int(fields[9]):
            start = max(1, int(fields[9]) - args.padding)
            end = min(int(fields[8]) + args.padding, len(fasta_dict[fields[1]].seq))

            identity = float(fields[2])
            aln_len = float(fields[3])
            if (identity * aln_len) / 50 >= 80:
                fout.write(f"{fields[1]}\t{start}\t{end}\n")
