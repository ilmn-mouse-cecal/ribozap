#!/usr/bin/env python3
import argparse
from pathlib import Path

# Argument parsing
parser = argparse.ArgumentParser(description="Merge overlapping or nested top coverage blocks from multiple samples.")
parser.add_argument("-t", "--top", type=int,  default=50, help="Top N string (e.g. '50')")
parser.add_argument("-i", "--input", required=True, help="Input BED file path")
parser.add_argument("-o", "--output", required=True, help="Output BED file path")
args = parser.parse_args()

# Paths
top = args.top
input_path = Path(args.input)
output_path = Path(args.output)

# Store merged regions
merged = []

with input_path.open("r") as fin:
    for line in fin:
        fields = line.strip().split('\t')
        genome, start, end, cov = fields[0], int(fields[1]), int(fields[2]), int(float(fields[4]))
        updated = False

        for item in merged:
            if item['Genome'] != genome:
                continue

            # Check overlap
            if (start <= item['End'] and end >= item['Start']):
                item['Start'] = min(item['Start'], start)
                item['End'] = max(item['End'], end)
                item['Cov'] = max(item['Cov'], cov)
                updated = True
                break

        if not updated:
            merged.append({'Genome': genome, 'Start': start, 'End': end, 'Cov': cov})

# ouput
with output_path.open("w") as fout:
    for region in merged:
        length = region['End'] - region['Start'] + 1
        if length >= top:
            fout.write(f"{region['Genome']}\t{region['Start']}\t{region['End']}\t{length}\t{region['Cov']}\n")
