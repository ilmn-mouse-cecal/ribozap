#!/usr/bin/env python3

import argparse
from pathlib import Path

# Argument parsing
parser = argparse.ArgumentParser(description="Generate 50bp probe windows with 25bp gaps from FASTA index (.fai).")
parser.add_argument("fasta_fai", help="Path to the FASTA index (.fai) file")
parser.add_argument("output_bed", help="Path to output BED file")
args = parser.parse_args()

fai_path = Path(args.fasta_fai)
output_path = Path(args.output_bed)

with output_path.open("w") as fout, fai_path.open("r") as fin:
    for line in fin:
        chrom, length = line.strip().split('\t')[:2]
        length = int(length)
        probe_num = 1
        i = length
        while i > 50:
            start = i - 50
            end = i
            fout.write(f"{chrom}\t{start}\t{end}\t{chrom}_Probe{probe_num}\t1\t-\n")
            i -= 75
            probe_num += 1
