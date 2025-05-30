import csv
import sys
import argparse
import logging
from pathlib import Path
from collections import OrderedDict


logger = logging.getLogger(__name__)

def parse_sample_sheet(sheet_path):
    sheet_path = Path(sheet_path).resolve()
    with sheet_path.open() as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        if not header or "sample_id" not in header or "read1" not in header:
            raise ValueError(
                "Sample sheet must have at least 'sample_id' and 'read1' columns"
            )
        is_paired = "read2" in header
        rows = []
        for row in reader:
            cleaned_row = {k: v.strip() for k, v in row.items()}
            rows.append(cleaned_row)
    return rows, is_paired


def validate_sample_row(row, is_paired):
    sample_id = row["sample_id"]
    read1 = Path(row["read1"]).expanduser().resolve()
    if not read1.is_file():
        raise FileNotFoundError(f"read1 file not found: {read1}")

    read2 = None
    if is_paired:
        read2 = Path(row["read2"]).expanduser().resolve()
        if not read2.is_file():
            raise FileNotFoundError(f"read2 file not found: {read2}")
        if read1.parent != read2.parent:
            raise ValueError(
                f"read1 and read2 for sample '{sample_id}' must be in the same directory."
            )

    return sample_id, read1, read2


def get_mount_path(sample_id):
    return f"/mnt/data/{sample_id}"


def rewrite_row(sample_id, read1, read2, container_dir, is_paired):
    row = {"sample_id": sample_id, "read1": f"{container_dir}/{read1.name}"}
    if is_paired:
        row["read2"] = f"{container_dir}/{read2.name}"
    return row


def write_sample_sheet(output_path, rows, is_paired):
    with output_path.open("w", newline="") as fout:
        fieldnames = ["sample_id", "read1"] + (["read2"] if is_paired else [])
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def validate_and_rewrite(sheet_path, out_sheet, mount_file):
    logger.info(f"Validating and rewriting sample sheet: {sheet_path}")
    rows, is_paired = parse_sample_sheet(sheet_path)

    rewritten_rows = []
    mounts = OrderedDict()
    sample_dirs = {}

    for row in rows:
        sample_id, read1, read2 = validate_sample_row(row, is_paired)
        logger.info(f"Sample '{sample_id}' passed validation")
        host_dir = read1.parent

        if sample_id in sample_dirs and sample_dirs[sample_id] != host_dir:
            raise ValueError(
                f"Sample '{sample_id}' has reads from multiple directories:\n"
                f"- {sample_dirs[sample_id]}\n"
                f"- {host_dir}\n"
                f"Please consolidate files or use unique sample IDs."
            )

        sample_dirs[sample_id] = host_dir
        container_dir = get_mount_path(sample_id)
        mounts[str(host_dir)] = container_dir
        rewritten_rows.append(
            rewrite_row(sample_id, read1, read2, container_dir, is_paired)
        )

    write_sample_sheet(out_sheet, rewritten_rows, is_paired)
    logger.info(f"Rewritten sample sheet to: {out_sheet}")

    with mount_file.open("w") as mf:
        for host, container in mounts.items():
            mf.write(f'-v "{host}":"{container}"\n')
    
    logger.info(f"Docker mount file written to: {mount_file}")
    return out_sheet, mount_file


def main():
    parser = argparse.ArgumentParser(
        description="Rewrite sample sheet and generate Docker mount paths."
    )
    parser.add_argument(
        "--sample-sheet",
        required=True,
        type=Path,
        help="Path to input sample sheet (CSV)"
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("output"),
        help="Directory to store outputs [default: ./output]",
    )
    parser.add_argument(
        "--out-sample-sheet",
        default="rewritten_sample_sheet.csv",
        help="Name of rewritten sample sheet",
    )
    parser.add_argument(
        "--mount-file", default="docker_mounts.txt", help="Name of Docker mounts file"
    )

    args = parser.parse_args()

    try:
        args.out_dir.mkdir(parents=True, exist_ok=True)
        out_sheet = args.out_dir / args.out_sample_sheet
        mount_file = args.out_dir / args.mount_file

        rewritten_sheet, mount_file = validate_and_rewrite(
            args.sample_sheet, out_sheet, mount_file
        )

        logger.info(f"Rewritten sample sheet saved to: {rewritten_sheet}")
        logger.info(f"Docker mount flags saved to: {mount_file}")
        logger.info("Docker -v arguments:")
        with mount_file.open() as f:
            for line in f:
                logger.info(line.strip())

    except Exception as e:
        logger.error(f"❌ Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
