#!/usr/bin/env python3

import argparse
import json
from pathlib import Path
from Bio import SeqIO

def filter_hits(hit_dict):
    # Return regions with >=80% identity
    return [hit for hit, identity in hit_dict.items() if identity >= 80]

def get_hits_for_region(region_name, hits_json):
    # Get all regions hit by region_name with >=80% identity
    for entry in hits_json:
        if entry['Name'] == region_name:
            return filter_hits(entry['Hits'])
    return []

def get_sequence(region_name, fasta_dict):
    record = fasta_dict.get(region_name)
    return record.seq if record else None

def main():
    parser = argparse.ArgumentParser(description="Select non-redundant regions with ≥80% similarity for probe design.")
    parser.add_argument("-t", "--top", type=int, default=50, help="Top N regions (not used currently)")
    parser.add_argument("-b", "--blast_json", required=True, help="Input BLAST result JSON file")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file with region sequences")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file for selected regions")
    args = parser.parse_args()

    # Read BLAST JSON
    with open(args.blast_json, "r") as f:
        blast_data = json.load(f)

    # Build region hit info and sort
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

    # TODO: Filter --top N regions
    # keep_list = keep_list[:args.top]

    fasta_dict = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

    # Write selected sequences to output
    with open(args.output, "w") as fout:
        for region in keep_list:
            seq = get_sequence(region, fasta_dict)
            if seq:
                fout.write(f">{region}\n{seq}\n")

if __name__ == "__main__":
    main()
