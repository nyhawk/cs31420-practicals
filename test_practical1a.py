from pathlib import Path

import pytest

import practical1a


ROOT = Path(__file__).parent
DATASETS = ROOT / "test_datasets"


@pytest.mark.parametrize(
    ("filename", "sequence_type", "expected"),
    [
        (
            "multiple_sequences.fasta",
            "gene",
            [("seq1", "ATGC"), ("seq2", "GGTA")],
        ),
        (
            "special_characters.fasta",
            "gene",
            [
                ("seq|one: description (sample)[1]", "ACGTNRYKMSWBDHV-."),
                ("seq_two", "acgtnrykmswbdhv"),
            ],
        ),
        ("empty_header.fasta", "gene", [("", "ATGC"), ("", "TTAA")]),
        (
            "blank_lines.fasta",
            "gene",
            [("seq1_blank_lines", "ATGC"), ("seq2_blank_lines", "TTAA")],
        ),
        ("mixed_case.fasta", "gene", [("mixed_case", "acgTnRy")]),
        (
            "empty_sequence.fasta",
            "gene",
            [("empty_sequence", ""), ("nonempty_sequence", "ATGC")],
        ),
    ],
)
def test_reads_fasta_datasets(filename, sequence_type, expected):
    assert practical1a.read_fasta(DATASETS / filename, sequence_type) == expected


@pytest.mark.parametrize(
    ("sequence_type", "sequence"),
    [
        ("gene", "ACGTRYSWKMBDHVN-."),
        ("rna", "ACGUTRYSWKMBDHVN-."),
        ("protein", "ACDEFGHIKLMNPQRSTVWYBJZXOU*-."),
    ],
)
def test_accepts_sequence_alphabet(tmp_path, sequence_type, sequence):
    fasta_file = tmp_path / "alphabet.fasta"
    fasta_file.write_text(f">sequence\n{sequence}\n")

    assert practical1a.read_fasta(fasta_file, sequence_type) == [
        ("sequence", sequence)
    ]


@pytest.mark.parametrize(
    ("filename", "sequence_type", "message"),
    [
        ("no_header.fasta", "gene", "sequence data before first FASTA header"),
        (
            "sequence_before_header.fasta",
            "gene",
            "sequence data before first FASTA header",
        ),
        (
            "unexpected_header_marker.fasta",
            "gene",
            "unexpected '>'",
        ),
        (
            "invalid_sequence_symbol.fasta",
            "gene",
            "invalid gene sequence symbol",
        ),
    ],
)
def test_rejects_invalid_fasta_datasets(filename, sequence_type, message):
    with pytest.raises(ValueError, match=message):
        practical1a.read_fasta(DATASETS / filename, sequence_type)


@pytest.mark.parametrize("sequence_type", ["gene", "rna", "protein"])
def test_reads_project_biological_fastas(sequence_type):
    for gene in ("OPN1LW_datasets", "OPN1MW_datasets", "OPN1SW_datasets"):
        suffix = "protein.faa" if sequence_type == "protein" else (
            "rna.fna" if sequence_type == "rna" else "gene.fna"
        )
        fasta_file = (
            ROOT / gene / "ncbi_dataset" / "data" / suffix
        )
        assert practical1a.read_fasta(fasta_file, sequence_type)


def test_empty_file_returns_no_records(tmp_path):
    fasta_file = tmp_path / "empty.fasta"
    fasta_file.write_bytes(b"")

    assert practical1a.read_fasta(fasta_file, "gene") == []


@pytest.mark.parametrize("line_ending", [b"\n", b"\r\n", b"\r"])
def test_reads_line_endings_and_missing_final_newline(tmp_path, line_ending):
    fasta_file = tmp_path / "line_endings.fasta"
    fasta_file.write_bytes(
        line_ending.join([b">first", b"AC", b"GT", b">second", b"TTAA"])
    )

    assert practical1a.read_fasta(fasta_file, "gene") == [
        ("first", "ACGT"),
        ("second", "TTAA"),
    ]


def test_reads_utf8_bom(tmp_path):
    fasta_file = tmp_path / "bom.fasta"
    fasta_file.write_bytes(b"\xef\xbb\xbf>bom_header\nACGT\n")

    assert practical1a.read_fasta(fasta_file, "gene") == [
        ("bom_header", "ACGT")
    ]


def test_rejects_unknown_sequence_type(tmp_path):
    fasta_file = tmp_path / "sequence.fasta"
    fasta_file.write_text(">sequence\nACGT\n")

    with pytest.raises(ValueError, match="sequence_type must be one of"):
        practical1a.read_fasta(fasta_file, "unknown")
