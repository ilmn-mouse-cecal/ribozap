import argparse
from pathlib import Path
from .rewrite import validate_and_rewrite
from .docker_launcher import run_docker
import logging
import re
import sys

MIN_CPUS = 2
MIN_MEM = 10

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)

def main():
    parser = argparse.ArgumentParser(
        prog="ribozap",
        description="Run sample prep and containerized workflow for probe design.",
    )

    parser.add_argument(
        "--analysis-name",
        type=str,
        required=True,
        help="Name of the analysis"
    )

    parser.add_argument(
        "--sample-sheet",
        type=Path,
        required=True,
        help="Path to the input sample sheet (CSV)",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory to write outputs and mount in Docker",
    )

    parser.add_argument(
        "--cpus",
        default=4,
        type=int,
        help="Limit number of CPUs for Docker container (e.g. 4)",
    )

    parser.add_argument(
        "--memory", default=16, type=int, help="Limit memory for Docker container (e.g. 16)"
    )

    parser.add_argument(
        "--resume",
        action="store_true",
        help="Enable Nextflow -resume mode to continue from previous work directory",
    )

    parser.add_argument(
        "--num-cov-regions",
        default=50,
        type=int,
        help="Enter the number of top high-coverage regions from the BED file to keep (default: 50)"
    )

    parser.add_argument(
        "--image",
        default="ribozap",
        type=str,
        help="Enter the image name you would like to use. Default is 'ribozap'"
    )

    parser.add_argument(
        "--image-tag",
        default='latest',
        type=str,
        help="Enter the docker image tag. Default is 'latest'"
    )

    args = parser.parse_args()
    
    if args.cpus < MIN_CPUS:
        logging.error(f"Minimum required CPUs is {MIN_CPUS}. You provided {args.cpus}.")
        sys.exit(1)
    
    if args.memory < MIN_MEM:
        logging.error(f"Minimum required memory is {MIN_MEM} GB. You provided {args.memory} GB.")
        sys.exit(1)

    invalid = re.search(r"[^A-Za-z0-9_.-]", args.analysis_name)
    if invalid:
        logging.error("Invalid analysis name: only letters, numbers, dashes, underscores, and dots are allowed . . .")
        sys.exit(1)

    high_cpus = args.cpus
    high_memory = args.memory
    out_dir = args.output_dir
    analysis_name = args.analysis_name
    out_dir.mkdir(parents=True, exist_ok=True)
    num_coverage_regions = args.num_cov_regions

    rewritten_path = out_dir / "rewritten_sample_sheet.csv"
    mount_path = out_dir / "docker_mounts.txt"

    # Step 1: Rewrite sample sheet
    validate_and_rewrite(args.sample_sheet, rewritten_path, mount_path)

    resume_flag = "-resume" if args.resume else ""

    # Step 2: Run Docker
    run_docker(
        image=f"{args.image}:{args.image_tag}",
        mount_file=mount_path,
        out_dir=out_dir,
        analysis_name=analysis_name,
        container_cmd=f"nextflow run main.nf -work-dir /app/{out_dir.name}/{analysis_name}/work/ --sample_sheet /app/{out_dir.name}/{rewritten_path.name} --outdir /app/{out_dir.name}/{analysis_name} --trace_dir /app/{out_dir.name}/{analysis_name}/trace_dir --top_coverage_regions {num_coverage_regions} --cpus {high_cpus} --memory '{high_memory} GB' {resume_flag}",
        cpus=args.cpus,
        memory=args.memory
    )


if __name__ == "__main__":
    main()
