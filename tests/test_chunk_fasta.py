from sat.scripts.chunk_fasta import chunk_seq


def test_chunk_seq_not_overlapping():
    seq = "ABCDEFGHIJKLMNO"
    expected = ["ABCDEF", "GHIJKL", "MNO"]
    assert chunk_seq(seq, max_length=6, output_overlapping=False) == expected


def test_chunk_seq_overlapping():
    seq = "ABCDEFGHIJKLMNO"
    expected = ["ABCDEF", "DEFGHI", "GHIJKL", "JKLMNO"]
    assert chunk_seq(seq, max_length=6, output_overlapping=True) == expected


def test_chunk_seq_overlapping_ODD():
    # round(7/2) gives 4
    seq = "ABCDEFGHIJKLMNO"
    expected = ["ABCDEFG", "DEFGHIJ", "GHIJKLM", "JKLMNO"]
    assert chunk_seq(seq, max_length=7, output_overlapping=True) == expected
