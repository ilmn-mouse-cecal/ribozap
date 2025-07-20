#!/usr/bin/env python3
import argparse
from pathlib import Path

def extract_high_cov_blocks(input_path, output_path, high_cov):
    prev_genome = ""
    block_flag = False
    block_start = block_end = max_cov = -1

    with output_path.open("w") as fout:
        fout.write("Genome Name\tBlock Start\tBlock End\tBlock Max Coverage\n")

        with input_path.open("r") as fin:
            for line in fin:
                try:
                    genome, start, end, cov = line.strip().split('\t')
                    start, end, cov = int(start), int(end), int(cov)
                except Exception as e:
                    print(f"Skipping malformed line: {line.strip()}")
                    continue

                if genome == prev_genome:
                    if block_flag:
                        if cov >= high_cov:
                            block_end = end
                            max_cov = max(max_cov, cov)
                        else:
                            fout.write(f"{prev_genome}\t{block_start}\t{block_end}\t{max_cov}\n")
                            block_flag = False
                    else:
                        if cov >= high_cov:
                            block_start, block_end = start, end
                            max_cov = cov
                            block_flag = True
                else:
                    if block_flag:
                        fout.write(f"{prev_genome}\t{block_start}\t{block_end}\t{max_cov}\n")
                    block_flag = False
                    if cov >= high_cov:
                        block_start, block_end = start, end
                        max_cov = cov
                        block_flag = True

                prev_genome = genome

        if block_flag:
            fout.write(f"{prev_genome}\t{block_start}\t{block_end}\t{max_cov}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract high-coverage blocks from a genomeCoverage.bed file.")
    parser.add_argument("-s", "--sample", required=True, help="Sample name (used for file naming)")
    parser.add_argument("-c", "--cov", required=True, help="Genome coverage file")
    parser.add_argument("--high", type=int, required=True, help="Minimum coverage for high coverage block")
    parser.add_argument("-o", "--output", help="Output file name (optional)")

    args = parser.parse_args()

    input_path = Path(args.cov)
    output_path = Path(args.output) if args.output else Path(f"{args.sample}_high_coverage_blocks.tsv")

    print(f"Identifying Blocks for {args.sample}.\nHigh Cov: {args.high}+")
    extract_high_cov_blocks(input_path, output_path, args.high)
