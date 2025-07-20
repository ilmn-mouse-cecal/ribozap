import subprocess
from pathlib import Path

def test_cli_probe_generation(tmp_path):
    # Prepare a fake .fai input
    fai = tmp_path / "test.fai"
    fai.write_text("chr1\t200\nchr2\t100\n")
    out_bed = tmp_path / "probes.bed"

    # cmd: python bin/make_bed_25bp_gap_v2.py <fai> <output_bed>
    cmd = [
        "python", "bin/make_bed_25bp_gap_v2.py",
        str(fai),
        str(out_bed)
    ]

    # Run the script
    result = subprocess.run(cmd, check=True)

    # Now check the output BED
    lines = out_bed.read_text().splitlines()
    assert any(line.startswith("chr1\t150\t200") for line in lines)
    assert any(line.startswith("chr2\t50\t100") for line in lines)
    assert len(lines) == 3
    assert result.returncode == 0
