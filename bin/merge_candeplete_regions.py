#!/usr/bin/env python3
import argparse

# Argument parsing
parser = argparse.ArgumentParser(description="Merge overlapping regions from a sorted deplete BED file.")
parser.add_argument("-i", "--input", required=True, help="Input sorted BED file")
parser.add_argument("-o", "--output", required=True, help="Output file with merged regions")
args = parser.parse_args()

with open(args.input, "r") as fin, open(args.output, "w") as fout:
    prev_genome = None
    start = end = None

    for line in fin:
        genome, s, e = line.strip().split("\t")
        s, e = int(s), int(e)

        if prev_genome is None:
            # First line setup
            prev_genome = genome
            start = s
            end = e
        elif genome != prev_genome:
            # New genome, write previous and reset
            fout.write(f"{prev_genome}\t{start}\t{end}\n")
            prev_genome = genome
            start = s
            end = e
        elif s <= end:  # Overlapping
            end = max(end, e)
        else:
            # Non-overlapping, write current block and start new one
            fout.write(f"{prev_genome}\t{start}\t{end}\n")
            start = s
            end = e

    # Writing last block
    if prev_genome is not None:
        fout.write(f"{prev_genome}\t{start}\t{end}\n")
