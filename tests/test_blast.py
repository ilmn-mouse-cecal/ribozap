import subprocess
import json
from pathlib import Path


def test_blast_script(tmp_path):
    # Write test FASTA
    fasta_content = ">seq1\nACTGACTGACTG\n"
    fasta_file = tmp_path / "test.fa"
    fasta_file.write_text(fasta_content)

    # Make BLAST DB
    blast_db = tmp_path / "blastdb"
    db_fasta = tmp_path / "db.fa"
    db_fasta.write_text(fasta_content)
    subprocess.run([
        "makeblastdb",
        "-in", str(db_fasta),
        "-dbtype", "nucl",
        "-out", str(blast_db),
    ], check=True)

    # Create output dir
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    script_path = Path("bin/blast.py")
    subprocess.run([
        "python3", str(script_path),
        "-f", str(fasta_file),
        "-d", str(blast_db),
        "-o", str(output_dir)
    ], check=True)

    # Check outputs
    seq_list = output_dir / "seqList.txt"
    assert seq_list.exists()
    lines = seq_list.read_text().splitlines()
    assert "seq1" in lines

    blast_json = output_dir / "blastHitResult.json"
    assert blast_json.exists()
    results = json.loads(blast_json.read_text())
    assert isinstance(results, list)
    assert results[0]["Name"] == "seq1"
    assert isinstance(results[0]["Hits"], dict)
    assert len(results[0]["Hits"]) == 0 
