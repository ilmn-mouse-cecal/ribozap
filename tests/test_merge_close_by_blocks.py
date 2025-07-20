import subprocess
from pathlib import Path

def test_merge_high_coverage_blocks(tmp_path):
    tsv = tmp_path / "blocks.tsv"
    tsv.write_text(
        "Genome Name\tBlock Start\tBlock End\tBlock Max Coverage\n"
        "chr1\t0\t50\t10\n"
        "chr1\t53\t100\t15\n"
        "chr1\t105\t150\t9\n"
        "chr2\t0\t20\t8\n"
        "chr2\t23\t50\t12\n"
    )
    sample_name = "mysample"
    output = Path(f"{sample_name}_high_coverage_blocks_gap_merged.tsv")
    cmd = [
        "python", "bin/merge_close_by_blocks.py",
        "-s", sample_name,
        "-c", str(tsv)
    ]
    subprocess.run(cmd, check=True)
    out = output.read_text().splitlines()
    assert out[0] == "Genome Name\tBlock Start\tBlock End\tBlock Max Coverage"
    # chr1: 0-100 (merged), 105-150 (separate)
    assert out[1] == "chr1\t0\t100\t15"
    assert out[2] == "chr1\t105\t150\t9"
    # chr2: 0-50 (merged)
    assert out[3] == "chr2\t0\t50\t12"
    assert len(out) == 4
