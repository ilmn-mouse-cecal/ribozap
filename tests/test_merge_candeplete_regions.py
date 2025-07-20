import subprocess

def test_bed_merge(tmp_path):
    input_bed = tmp_path / "in.bed"
    input_bed.write_text(
        "chr1\t10\t50\n"
        "chr1\t40\t100\n"
        "chr1\t120\t150\n"
        "chr2\t5\t25\n"
        "chr2\t20\t30\n"
        "chr2\t40\t60\n"
    )
    output_bed = tmp_path / "out.bed"
    cmd = [
        "python", "bin/merge_candeplete_regions.py",
        "-i", str(input_bed),
        "-o", str(output_bed)
    ]
    subprocess.run(cmd, check=True)
    out = output_bed.read_text().splitlines()
    # chr1: 10-100 (overlap), 120-150
    # chr2: 5-30 (overlap), 40-60
    assert out[0] == "chr1\t10\t100"
    assert out[1] == "chr1\t120\t150"
    assert out[2] == "chr2\t5\t30"
    assert out[3] == "chr2\t40\t60"
    assert len(out) == 4
