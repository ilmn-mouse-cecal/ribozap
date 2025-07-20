import json
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
from identify_heatmap_blocks_v2 import filter_hits, get_hits_for_region, get_sequence

def test_filter_hits():
    hit_dict = {"regA": 80, "regB": 79, "regC": 95}
    assert sorted(filter_hits(hit_dict)) == ["regA", "regC"]

def test_get_hits_for_region():
    hits_json = [
        {"Name": "R1", "Hits": {"R2": 90, "R3": 75}},
        {"Name": "R2", "Hits": {"R1": 91}},
    ]
    assert get_hits_for_region("R1", hits_json) == ["R2"]
    assert get_hits_for_region("R2", hits_json) == ["R1"]
    assert get_hits_for_region("missing", hits_json) == []

def test_get_sequence(tmp_path):
    # Create a sample FASTA
    fasta_file = tmp_path / "seqs.fa"
    records = [
        SeqRecord(Seq("ACGT"), id="foo", description=""),
        SeqRecord(Seq("GGGG"), id="bar", description=""),
    ]
    SeqIO.write(records, fasta_file, "fasta")
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    assert str(get_sequence("foo", fasta_dict)) == "ACGT"
    assert str(get_sequence("bar", fasta_dict)) == "GGGG"
    assert get_sequence("nope", fasta_dict) is None

def test_full_pipeline(tmp_path):
    # Create small BLAST JSON
    blast_json = [
        {"Name": "A", "Hits": {"B": 85, "C": 70}},
        {"Name": "B", "Hits": {"A": 82}},
        {"Name": "C", "Hits": {}},
    ]
    blast_json_path = tmp_path / "blast.json"
    blast_json_path.write_text(json.dumps(blast_json))

    # Create sample FASTA
    records = [
        SeqRecord(Seq("AAAA"), id="A", description=""),
        SeqRecord(Seq("BBBB"), id="B", description=""),
        SeqRecord(Seq("CCCC"), id="C", description=""),
    ]
    fasta_path = tmp_path / "seqs.fa"
    SeqIO.write(records, fasta_path, "fasta")

    # Build region hit info and sort
    hits_info = [{"Name": entry["Name"], "HitsNum": len(entry["Hits"])} for entry in blast_json]
    hits_info.sort(key=lambda x: x["HitsNum"], reverse=True)

    # Select non-overlapping regions
    ignore_list = set()
    keep_list = []

    for info in hits_info:
        region = info["Name"]
        if region in ignore_list:
            continue
        keep_list.append(region)
        ignore_list.update(get_hits_for_region(region, blast_json))

    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

    # Write to output
    output_fasta = tmp_path / "output.fa"
    with open(output_fasta, "w") as fout:
        for region in keep_list:
            seq = get_sequence(region, fasta_dict)
            if seq:
                fout.write(f">{region}\n{seq}\n")

    # Check output
    written = list(SeqIO.parse(output_fasta, "fasta"))
    ids = [rec.id for rec in written]
    assert "A" in ids
    assert "C" in ids
    assert "B" not in ids
