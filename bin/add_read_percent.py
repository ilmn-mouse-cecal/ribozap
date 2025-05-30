#!/usr/bin/env python3
import argparse
from pathlib import Path

# Argument parsing
parser = argparse.ArgumentParser(description="Add block length and percent of total coverage to sorted coverage BED file.")
parser.add_argument("-s", "--sample", required=True, help="Sample name (used in file naming)")
parser.add_argument("-c", "--cov-sorted-bed", required=True, help="High coverage blocks and merged gaps BED file")
parser.add_argument("-n", "--num_reads", type=float, required=True, help="Number of reads in the original FASTQ file")

args = parser.parse_args()

# Input and output paths
input_path = Path(args.cov_sorted_bed)
output_path = Path(f"{args.sample}_high_coverage_blocks_gap_merged_cov_sorted_percentage_added.bed")

with input_path.open("r") as fin, output_path.open("w") as fout:
    for line in fin:
        fields = line.strip().split('\t')
        chrom, start, end, cov = fields[0], int(fields[1]), int(fields[2]), float(fields[3])
        length = end - start + 1
        percentage = (cov / args.num_reads) * 100
        fout.write(f"{chrom}\t{start}\t{end}\t{length}\t{cov}\t{percentage:.2f}\n")
