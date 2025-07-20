import subprocess
from pathlib import Path

def test_filter_add_padding(tmp_path):
    # Create a reference fasta
    fasta_content = ">ref1\n" + "A" * 200 + "\n"
    fasta_file = tmp_path / "ref.fa"
    fasta_file.write_text(fasta_content)

    # Columns: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    tsv_content = "query1\tref1\t90.0\t100\t0\t0\t1\t100\t150\t50\t1e-10\t200\n"
    tsv_file = tmp_path / "blast.tsv"
    tsv_file.write_text(tsv_content)

    # Output file
    out_file = tmp_path / "output.tsv"

    script_path = Path("bin/filter_add_padding.py")
    subprocess.run([
        "python3", str(script_path),
        "-i", str(tsv_file),
        "-o", str(out_file),
        "-p", "10",
        "-f", str(fasta_file)
    ], check=True)

    # Read and check output
    lines = out_file.read_text().strip().split("\n")
    assert lines[0] == "ref1\t40\t160"

