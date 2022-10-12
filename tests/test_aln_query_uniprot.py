from sat.scripts.aln_query_uniprot import aln_query_uniprot_main
from sat.scripts.utils.uniprot import uniprot_object


def test_aln_query_uniprot_main_no_cache(tmp_path):

    # Define args
    class args:
        pass

    args.infile = "tests/test_data/foldseek_related/af2_foldseek_result_small.m8"
    args.uniprot_lookup_output = f"{tmp_path}/result.tsv"
    args.infile_col = 1
    args.uniprot_cache = ""

    # Run script
    aln_query_uniprot_main(args)

    # Parse the expected and observed into a dictionary of uniprot objects
    expected_uniprot_object_dict = dict()
    with open(
        "tests/test_data/foldseek_related/af2_foldseek_result_small_uniprot_lookup.tsv"
    ) as infile:
        for line in infile:
            line = line.rstrip("\n").rstrip("\t").split("\t")
            acc, geneName, fullName = line
            uo = uniprot_object(acc)
            uo.geneName = geneName
            uo.fullName = fullName
            expected_uniprot_object_dict[acc] = uo

    observed_uniprot_object_dict = dict()
    with open(args.uniprot_lookup_output) as infile:
        for line in infile:
            line = line.rstrip("\n").rstrip("\t").split("\t")
            acc, geneName, fullName = line
            uo = uniprot_object(acc)
            uo.geneName = geneName
            uo.fullName = fullName
            observed_uniprot_object_dict[acc] = uo

    # Make sure the values are teh same
    for acc, expected_uo in expected_uniprot_object_dict.items():
        observed_uo = observed_uniprot_object_dict[acc]
        assert expected_uo.__eq__(observed_uo)


def test_aln_query_uniprot_main_write_cache(tmp_path):

    # Define args
    class args:
        pass

    args.infile = "tests/test_data/foldseek_related/af2_foldseek_result_small.m8"
    args.uniprot_lookup_output = f"{tmp_path}/result.tsv"
    args.infile_col = 1
    args.uniprot_cache = f"{tmp_path}/new_cache.pkl"

    # Run script
    aln_query_uniprot_main(args)

    # Parse the expected and observed into a dictionary of uniprot objects
    expected_uniprot_object_dict = dict()
    with open(
        "tests/test_data/foldseek_related/af2_foldseek_result_small_uniprot_lookup.tsv"
    ) as infile:
        for line in infile:
            line = line.rstrip("\n").rstrip("\t").split("\t")
            acc, geneName, fullName = line
            uo = uniprot_object(acc)
            uo.geneName = geneName
            uo.fullName = fullName
            expected_uniprot_object_dict[acc] = uo

    observed_uniprot_object_dict = dict()
    with open(args.uniprot_lookup_output) as infile:
        for line in infile:
            line = line.rstrip("\n").rstrip("\t").split("\t")
            acc, geneName, fullName = line
            uo = uniprot_object(acc)
            uo.geneName = geneName
            uo.fullName = fullName
            observed_uniprot_object_dict[acc] = uo

    # Make sure the values are teh same
    for acc, expected_uo in expected_uniprot_object_dict.items():
        observed_uo = observed_uniprot_object_dict[acc]
        assert expected_uo.__eq__(observed_uo)


def test_aln_query_uniprot_main_read_cache(tmp_path):

    # Define args
    class args:
        pass

    args.infile = "tests/test_data/foldseek_related/af2_foldseek_result_small.m8"
    args.uniprot_lookup_output = f"{tmp_path}/result.tsv"
    args.infile_col = 1
    args.uniprot_cache = (
        "tests/test_data/foldseek_related/af2_foldseek_result_small_uniprot_cache.pkl"
    )

    # Run script
    aln_query_uniprot_main(args)

    # Parse the expected and observed into a dictionary of uniprot objects
    expected_uniprot_object_dict = dict()
    with open(
        "tests/test_data/foldseek_related/af2_foldseek_result_small_uniprot_lookup.tsv"
    ) as infile:
        for line in infile:
            line = line.rstrip("\n").rstrip("\t").split("\t")
            acc, geneName, fullName = line
            uo = uniprot_object(acc)
            uo.geneName = geneName
            uo.fullName = fullName
            expected_uniprot_object_dict[acc] = uo

    observed_uniprot_object_dict = dict()
    with open(args.uniprot_lookup_output) as infile:
        for line in infile:
            line = line.rstrip("\n").rstrip("\t").split("\t")
            acc, geneName, fullName = line
            uo = uniprot_object(acc)
            uo.geneName = geneName
            uo.fullName = fullName
            observed_uniprot_object_dict[acc] = uo

    # Make sure the values are teh same
    for acc, expected_uo in expected_uniprot_object_dict.items():
        observed_uo = observed_uniprot_object_dict[acc]
        assert expected_uo.__eq__(observed_uo)
