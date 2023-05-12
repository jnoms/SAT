import pytest


from sat.scripts.seq_multimerize import multimerize_seqs


# ------------------------------------------------------------------------------------ #
# Individual function tests
# ------------------------------------------------------------------------------------ #
def test_multimerize_seqs_normal_case():
    sequences = ["ATG", "TAG", "GTA"]
    cardinality = [2, 3, 1]
    expected = "ATG:ATG:TAG:TAG:TAG:GTA"
    assert multimerize_seqs(sequences, cardinality) == expected


def test_multimerize_seqs_dimer():
    sequences = ["ATG"]
    cardinality = [2]
    expected = "ATG:ATG"
    assert multimerize_seqs(sequences, cardinality) == expected


def test_multimerize_seqs_empty_case():
    sequences = []
    cardinality = []
    expected = ""
    assert multimerize_seqs(sequences, cardinality) == expected


def test_multimerize_seqs_mismatched_lengths():
    sequences = ["ATG", "TAG", "GTA"]
    cardinality = [2, 3]
    with pytest.raises(ValueError, match="sequences is a dif size than cardinality!"):
        multimerize_seqs(sequences, cardinality)


def test_multimerize_seqs_different_sequence_lengths():
    sequences = ["ATG", "TAGC", "G"]
    cardinality = [2, 2, 2]
    expected = "ATG:ATG:TAGC:TAGC:G:G"
    assert multimerize_seqs(sequences, cardinality) == expected
