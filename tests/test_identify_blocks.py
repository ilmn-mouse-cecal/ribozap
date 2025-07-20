import subprocess
from pathlib import Path

def test_script_runs(tmp_path):
    # Create a sample BED file
    bed = tmp_path / "test.bed"
    bed.write_text("chr1\t0\t100\t10\nchr1\t100\t200\t20\nchr1\t200\t300\t5\n")
    output = Path("mysample_high_coverage_blocks.tsv")
    cmd = [
        "python", "bin/identify_blocks.py",
        "-s", "mysample",
        "-c", str(bed),
        "--high", "15"
    ]
    subprocess.run(cmd, check=True)
    lines = output.read_text().splitlines()
    assert lines[0] == "Genome Name\tBlock Start\tBlock End\tBlock Max Coverage"
    assert lines[1] == "chr1\t100\t200\t20"
