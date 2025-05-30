#!/usr/bin/env python3
import argparse
from pathlib import Path

"""
This script reads in a list of high-coverage blocks (from idenify_blocks.py) and merges blocks that are within 3 bases of each other on the same genome.

Purpose:
To consolidate adjacent high-coverage regions into unified blocks . . .
"""

# Argument parsing
parser = argparse.ArgumentParser(description="Merge nearby high-coverage blocks allowing up to 3bp gap.")
parser.add_argument("-s", "--sample", required=True, help="Sample name (used for input/output filenames)")
parser.add_argument("-c", "--cov-blocks-tsv", required=True, help="High coverage blocks tsv file")

args = parser.parse_args()

# Paths
input_path = Path(args.cov_blocks_tsv)
output_path = Path(f"{args.sample}_high_coverage_blocks_gap_merged.tsv")

# Initialize state
prev_genome = ""
block_start = block_end = block_max = -1

with input_path.open("r") as fin, output_path.open("w") as fout:
    fout.write("Genome Name\tBlock Start\tBlock End\tBlock Max Coverage\n")
    next(fin)  # skip header

    for line in fin:
        genome, start, end, cov = line.strip().split('\t')
        start, end, cov = int(start), int(end), int(cov)

        if genome == prev_genome:
            # Check if current block is within 3bp of previous block
            if block_end + 3 >= start:
                block_end = end
                block_max = max(block_max, cov)
            else:
                fout.write(f"{prev_genome}\t{block_start}\t{block_end}\t{block_max}\n")
                block_start, block_end, block_max = start, end, cov
        else:
            if prev_genome != "":
                fout.write(f"{prev_genome}\t{block_start}\t{block_end}\t{block_max}\n")
            block_start, block_end, block_max = start, end, cov

        prev_genome = genome

    # Writing final block
    fout.write(f"{prev_genome}\t{block_start - 1}\t{block_end}\t{block_max}\n")
