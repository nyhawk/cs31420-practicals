from pathlib import Path

import pytest

import practical1a


DATASETS = Path(__file__).parent / "test_datasets"


@pytest.mark.parametrize(
    ("filename", "expected"),
    [
        (
            "multiple_sequences.fasta",
            [("seq1", "ATGC"), ("seq2", "GGTA")],
        ),
        (
            "special_characters.fasta",
            [
                ("seq|one: description (sample)[1]", "ACGTNRYKMSWBDHV-."),
                ("seq_two", "acgtnrykmswbdhv"),
            ],
        ),
        ("empty_header.fasta", [("", "ATGC"), ("", "TTAA")]),
        (
            "blank_lines.fasta",
            [
                ("seq1_blank_lines", "ATGC"),
                ("seq2_blank_lines", "TTAA"),
            ],
        ),
        ("mixed_case.fasta", [("mixed_case", "acgTnRy")]),
        (
            "empty_sequence.fasta",
            [("empty_sequence", ""), ("nonempty_sequence", "ATGC")],
        ),
    ],
)
def test_reads_valid_dataset_records(filename, expected):
    assert practical1a.read_fasta(DATASETS / filename) == expected


@pytest.mark.parametrize(
    ("filename", "error"),
    [
        ("no_header.fasta", "sequence data before first FASTA header"),
        ("sequence_before_header.fasta", "sequence data before first FASTA header"),
        ("unexpected_header_marker.fasta", "unexpected '>'"),
        ("invalid_sequence_symbol.fasta", "invalid sequence character"),
    ],
)
def test_rejects_malformed_dataset_records(filename, error):
    with pytest.raises(ValueError, match=error):
        practical1a.read_fasta(DATASETS / filename)


def test_empty_file_returns_no_records(tmp_path):
    fasta_file = tmp_path / "empty.fasta"
    fasta_file.write_bytes(b"")

    assert practical1a.read_fasta(fasta_file) == []


@pytest.mark.parametrize("line_ending", [b"\n", b"\r\n", b"\r"])
def test_reads_different_line_endings_and_missing_final_newline(
    tmp_path, line_ending
):
    fasta_file = tmp_path / "line_endings.fasta"
    fasta_file.write_bytes(
        line_ending.join([b">first", b"AC", b"GT", b">second", b"TTAA"])
    )

    assert practical1a.read_fasta(fasta_file) == [
        ("first", "ACGT"),
        ("second", "TTAA"),
    ]


def test_reads_utf8_bom(tmp_path):
    fasta_file = tmp_path / "bom.fasta"
    fasta_file.write_bytes(b"\xef\xbb\xbf>bom_header\nACGT\n")

    assert practical1a.read_fasta(fasta_file) == [("bom_header", "ACGT")]
