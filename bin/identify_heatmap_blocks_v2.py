#!/usr/bin/env python3

import argparse
import json
from pathlib import Path
from Bio import SeqIO

def filter_hits(hit_dict):
    return [hit for hit, identity in hit_dict.items() if identity >= 80]

def get_hits_for_region(region_name, hits_json):
    for entry in hits_json:
        if entry['Name'] == region_name:
            return filter_hits(entry['Hits'])
    return []

def get_sequence(region_name, fasta_path):
    for record in SeqIO.parse(fasta_path, "fasta"):
        if record.id == region_name:
            return record.seq
    return None

# Parse arguments
parser = argparse.ArgumentParser(description="Select non-redundant regions with ≥80% similarity for probe design.")
parser.add_argument("-t", "--top", type=str, default="50", help="Top N string (e.g. '50')")
parser.add_argument("-b", "--blast_json", required=True, help="Input BLAST result JSON file")
parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file with region sequences")
parser.add_argument("-o", "--output", required=True, help="Output FASTA file for selected regions")
args = parser.parse_args()

# Read BLAST JSON
with open(args.blast_json, "r") as f:
    blast_data = json.load(f)

# Build region hit info
hits_info = [{"Name": entry["Name"], "HitsNum": len(entry["Hits"])} for entry in blast_data]
hits_info.sort(key=lambda x: x["HitsNum"], reverse=True)

# Select non-overlapping regions
ignore_list = set()
keep_list = []

for info in hits_info:
    region = info["Name"]
    if region in ignore_list:
        continue
    keep_list.append(region)
    ignore_list.update(get_hits_for_region(region, blast_data))

# Write selected sequences to output
with open(args.output, "w") as fout:
    for region in keep_list:
        seq = get_sequence(region, args.fasta)
        if seq:
            fout.write(f">{region}\n{seq}\n")
