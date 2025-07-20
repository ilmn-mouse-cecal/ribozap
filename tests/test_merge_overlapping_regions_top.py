import subprocess
from pathlib import Path

def test_merge_overlapping_top_blocks(tmp_path):
    bed = tmp_path / "test.bed"
    bed.write_text(
        "chr1\t0\t40\tprobe1\t10\n"
        "chr1\t30\t90\tprobe2\t20\n"
        "chr1\t95\t120\tprobe3\t30\n"
        "chr2\t0\t30\tprobeA\t15\n"
        "chr2\t28\t80\tprobeB\t12\n"
    )
    out_bed = tmp_path / "out.bed"
    script_path = "bin/merge_overlapping_regions_top.py"

    cmd = [
        "python", script_path,
        "-i", str(bed),
        "-o", str(out_bed),
        "-t", "50"
    ]
    subprocess.run(cmd, check=True)
    out = out_bed.read_text().splitlines()

    # Should merge chr1: 0-90 (max cov 20), chr1:95-120 (30), chr2:0-80 (15)
    assert any(line.startswith("chr1\t0\t90\t91\t20") for line in out)
    assert not any(line.startswith("chr1\t95\t120") for line in out)
    assert any(line.startswith("chr2\t0\t80\t81\t15") for line in out)
    # All regions should be at least 50bp (top)
    for line in out:
        fields = line.split("\t")
        assert int(fields[3]) >= 50

    # Output should not include any merged region shorter than 50bp
    assert all(int(line.split("\t")[3]) >= 50 for line in out)