import subprocess
from pathlib import Path
import csv

def test_fasta_to_csv(tmp_path):
    fasta = tmp_path / "test.fa"
    fasta.write_text(
        ">probe1\nACGTACGT\n"
        ">probe2\nTTTTTT\n"
    )
    out_csv = tmp_path / "probes.csv"
    script_path = "bin/summarize_probes.py"

    # Run script
    cmd = ["python", script_path, str(fasta), str(out_csv)]
    subprocess.run(cmd, check=True)

    # Check CSV output
    rows = list(csv.reader(open(out_csv)))
    assert rows[0] == ["Probe ID", "Sequence", "Length"]
    assert ["probe1", "ACGTACGT", "8"] in rows
    assert ["probe2", "TTTTTT", "6"] in rows
    assert len(rows) == 3
    