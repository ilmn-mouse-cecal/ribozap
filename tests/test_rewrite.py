import csv
import sys
from pathlib import Path
import tempfile
import pytest
from ribozap.rewrite import validate_and_rewrite, write_mounts


def _write_sample_sheet(path, rows, header=("sample_id", "read1", "read2")):
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def test_valid_paired_sample_sheet():
    test_dir = Path("tests/data/Sample1")
    read1 = test_dir / "sample1_R1.fastq"
    read2 = test_dir / "sample1_R2.fastq"

    assert read1.exists(), f"Missing test file: {read1}"
    assert read2.exists(), f"Missing test file: {read2}"

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        sample_sheet = tmp_path / "samples.csv"
        out_sheet = tmp_path / "rewritten_sample_sheet.csv"
        mount_file = tmp_path / "docker_mounts.txt"

        _write_sample_sheet(
            sample_sheet, [["Sample1", str(read1.resolve()), str(read2.resolve())]]
        )

        rewritten_sheet, mounts_txt = validate_and_rewrite(
            sample_sheet, out_sheet, mount_file
        )

        assert rewritten_sheet.exists()
        assert mount_file.exists()


def test_missing_read1_file():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        sample_sheet = tmp_path / "bad_samples.csv"
        out_sheet = tmp_path / "rewritten.csv"
        mount_file = tmp_path / "mounts.txt"

        _write_sample_sheet(
            sample_sheet,
            [
                [
                    "Sample1",
                    str(tmp_path / "missing_R1.fastq"),
                    str(tmp_path / "read2.fastq"),
                ]
            ],
        )

        with pytest.raises(FileNotFoundError, match="read1 file not found"):
            validate_and_rewrite(sample_sheet, out_sheet, mount_file)


def test_missing_read2_file():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        read1 = tmp_path / "read1.fastq"
        read1.write_text("@SEQ_ID\nGATTACA\n+\n!!!!!!\n")

        sample_sheet = tmp_path / "bad_samples.csv"
        out_sheet = tmp_path / "rewritten.csv"
        mount_file = tmp_path / "mounts.txt"

        _write_sample_sheet(
            sample_sheet,
            [["Sample1", str(read1.resolve()), str(tmp_path / "missing_R2.fastq")]],
        )

        with pytest.raises(FileNotFoundError, match="read2 file not found"):
            validate_and_rewrite(sample_sheet, out_sheet, mount_file)


def test_read1_and_read2_from_different_dirs():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        dir1 = tmp_path / "A"
        dir2 = tmp_path / "B"
        dir1.mkdir()
        dir2.mkdir()

        read1 = dir1 / "read1.fastq"
        read2 = dir2 / "read2.fastq"
        read1.write_text("@SEQ_ID\nGATTACA\n+\n!!!!!!\n")
        read2.write_text("@SEQ_ID\nGATTACA\n+\n!!!!!!\n")

        sample_sheet = tmp_path / "samples.csv"
        out_sheet = tmp_path / "rewritten.csv"
        mount_file = tmp_path / "mounts.txt"

        _write_sample_sheet(sample_sheet, [["SampleX", str(read1), str(read2)]])

        with pytest.raises(ValueError, match="must be in the same directory"):
            validate_and_rewrite(sample_sheet, out_sheet, mount_file)


def test_same_sample_multiple_dirs_error():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        dir1 = tmp_path / "dir1"
        dir2 = tmp_path / "dir2"
        dir1.mkdir()
        dir2.mkdir()

        r1a = dir1 / "read1.fastq"
        r2a = dir1 / "read2.fastq"
        r1b = dir2 / "read1.fastq"
        r2b = dir2 / "read2.fastq"

        for f in [r1a, r2a, r1b, r2b]:
            f.write_text("@SEQ_ID\nGATTACA\n+\n!!!!!!\n")

        sample_sheet = tmp_path / "samples.csv"
        out_sheet = tmp_path / "rewritten.csv"
        mount_file = tmp_path / "mounts.txt"

        _write_sample_sheet(
            sample_sheet, [["dup", str(r1a), str(r2a)], ["dup", str(r1b), str(r2b)]]
        )

        with pytest.raises(ValueError, match="multiple directories"):
            validate_and_rewrite(sample_sheet, out_sheet, mount_file)

def test_write_mounts():
    mounts = {
        "/host/path1": "/container/path1",
        "/host/path2": "/container/path2"
    }
    with tempfile.TemporaryDirectory() as tmpdir:
        mount_file = Path(tmpdir) / "mounts.txt"
        write_mounts(mount_file, mounts)
        with mount_file.open() as f:
            lines = f.readlines()
        assert lines == [
            '-v "/host/path1":"/container/path1"\n',
            '-v "/host/path2":"/container/path2"\n'
        ]