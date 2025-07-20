import subprocess

def test_cov_percentage(tmp_path):
    bed_content = "chr1\t100\t200\t500\nchr1\t300\t350\t100\n"
    num_reads = 1000

    # Write test input file
    in_file = tmp_path / "input.bed"
    in_file.write_text(bed_content)

    expected = [
        "chr1\t100\t200\t101\t500.0\t50.00",
        "chr1\t300\t350\t51\t100.0\t10.00",
    ]
    
    out_file = tmp_path / "sample_high_coverage_blocks_gap_merged_cov_sorted_percentage_added.bed"
    subprocess.run([
        "python3", "bin/add_read_percent.py", 
        "-s", "sample", 
        "-c", str(in_file), 
        "-n", str(num_reads)
    ], check=True)

    # Read and check output
    lines = [line.strip() for line in out_file.read_text().splitlines()]
    assert lines == expected
