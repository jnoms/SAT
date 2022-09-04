from sat.scripts.chunk_fasta import chunk_seq, chunk_fasta_main
from sat.scripts.utils.misc import read_fasta_to_memory


# ------------------------------------------------------------------------------------ #
# Individual function tests
# ------------------------------------------------------------------------------------ #
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


# ------------------------------------------------------------------------------------ #
# Whole-script tests
# ------------------------------------------------------------------------------------ #
def test_chunk_fasta_nonoverlapping(tmp_path):

    # Define inputs
    class args:
        pass

    args.in_fasta = "tests/test_data/t.fasta"
    args.out_fasta = f"{tmp_path}/observed.fasta"
    args.max_seq_length = 50
    args.minimum_sequence_output_size = 20
    args.overlapping_chunks = False
    args.individual = False

    # Run script
    chunk_fasta_main(args)

    # Check that the output is expected
    expected = {
        "PART0_100AA_seq": "MSSATGEGSQGARATYRAALNNEKRHDHVALTVPCCGTEAKVTALSPWFM",
        "PART1_100AA_seq": "DGMLAYETVKEMLLKGEQLLFAPSNLSGYIKFLPGPRVYLVERLTGGTYS",
    }
    assert read_fasta_to_memory(args.out_fasta) == expected


def test_chunk_fasta_nonoverlapping_two(tmp_path):

    # Define inputs
    class args:
        pass

    args.in_fasta = "tests/test_data/t.fasta"
    args.out_fasta = f"{tmp_path}/observed.fasta"
    args.max_seq_length = 50
    args.minimum_sequence_output_size = 10
    args.overlapping_chunks = False
    args.individual = False

    # Run script
    chunk_fasta_main(args)

    # Check that the output is expected
    expected = {
        "17AA_seq": "DPKGKFSNKLYKKLCGG",
        "PART0_100AA_seq": "MSSATGEGSQGARATYRAALNNEKRHDHVALTVPCCGTEAKVTALSPWFM",
        "PART1_100AA_seq": "DGMLAYETVKEMLLKGEQLLFAPSNLSGYIKFLPGPRVYLVERLTGGTYS",
    }
    assert read_fasta_to_memory(args.out_fasta) == expected


def test_chunk_fasta_overlapping(tmp_path):

    # Define inputs
    class args:
        pass

    args.in_fasta = "tests/test_data/t.fasta"
    args.out_fasta = f"{tmp_path}/observed.fasta"
    args.max_seq_length = 50
    args.minimum_sequence_output_size = 20
    args.overlapping_chunks = True
    args.individual = False

    # Run script
    chunk_fasta_main(args)

    # Check that the output is expected
    expected = {
        "PART0_100AA_seq": "MSSATGEGSQGARATYRAALNNEKRHDHVALTVPCCGTEAKVTALSPWFM",
        "PART1_100AA_seq": "HDHVALTVPCCGTEAKVTALSPWFMDGMLAYETVKEMLLKGEQLLFAPSN",
        "PART2_100AA_seq": "DGMLAYETVKEMLLKGEQLLFAPSNLSGYIKFLPGPRVYLVERLTGGTYS",
    }
    assert read_fasta_to_memory(args.out_fasta) == expected


def test_chunk_fasta_overlapping_individual(tmp_path):

    # Define inputs
    class args:
        pass

    args.in_fasta = "tests/test_data/t.fasta"
    args.out_fasta = f"{tmp_path}/observed"
    args.max_seq_length = 50
    args.minimum_sequence_output_size = 20
    args.overlapping_chunks = True
    args.individual = True

    # Run script
    chunk_fasta_main(args)

    # Check that the output is expected
    one = read_fasta_to_memory(f"{tmp_path}/observed/PART0_100AA_seq.fasta")
    two = read_fasta_to_memory(f"{tmp_path}/observed/PART1_100AA_seq.fasta")
    three = read_fasta_to_memory(f"{tmp_path}/observed/PART2_100AA_seq.fasta")

    assert one == {
        "PART0_100AA_seq": "MSSATGEGSQGARATYRAALNNEKRHDHVALTVPCCGTEAKVTALSPWFM"
    }
    assert two == {
        "PART1_100AA_seq": "HDHVALTVPCCGTEAKVTALSPWFMDGMLAYETVKEMLLKGEQLLFAPSN"
    }
    assert three == {
        "PART2_100AA_seq": "DGMLAYETVKEMLLKGEQLLFAPSNLSGYIKFLPGPRVYLVERLTGGTYS"
    }
