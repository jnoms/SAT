from sat.scripts.remove_redundant_domains import is_overlapping


def test_is_overlapping():

    # Seq2 within seq1
    seq1 = "ABCDEFGHIJKLMNO"
    seq2 = "DEFGH"
    assert is_overlapping(seq1, seq2, 5) == True

    # No overlap
    seq1 = "ABCDEFGHIJKLxxMNO"
    seq2 = "GHIJKLMNOPQRST"
    assert is_overlapping(seq1, seq2, 5) == False

    # same string
    seq1 = "THIS_IS_A_SEQ"
    seq2 = "THIS_IS_A_SEQ"
    assert is_overlapping(seq1, seq2, 5) == True

    # Seq1 is an extended part of seq2... need to trim off 3 to see it. 2 doesn't work.
    seq1 = "ABCDEFGHIJKLMNO"
    seq2 = "DEFGHIJKLMNOPQRST"
    assert is_overlapping(seq1, seq2, 3) == True
    assert is_overlapping(seq1, seq2, 2) == False
