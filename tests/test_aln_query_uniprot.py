from sat.scripts.aln_query_uniprot import aln_query_uniprot_main
from sat.scripts.utils.uniprot import uniprot_object


def compare_uniprot_objects(expected_uniprot_object_dict, observed_uniprot_object_dict):

    # Compare expected and observed. Need to handle the fact that uniprot can update
    # over time, so my downloaded test data will eventually get old. As long as at least
    # 1/2 of the entries are exactly the same lets pass this test.

    total = len(expected_uniprot_object_dict)
    test_data_probably_old_count = 0

    for acc, expected_uo in expected_uniprot_object_dict.items():
        observed_uo = observed_uniprot_object_dict[acc]

        try:
            assert expected_uo.__eq__(observed_uo)

        # Handle cases where uniprot has updated and the lookup information varies from my
        # old downloaded test data
        except AssertionError:
            assert expected_uo.accession == observed_uo.accession
            assert len(expected_uo.__dict__) == len(observed_uo.__dict__)
            test_data_probably_old_count += 1

    # Make sure at least half of the test data haven't updated.
    assert test_data_probably_old_count / total < 0.5


# def test_aln_query_uniprot_main_no_cache(tmp_path):

#     # Define args
#     class args:
#         pass

#     args.infile = "tests/test_data/foldseek_related/af2_foldseek_result_small.m8"
#     args.uniprot_lookup_output = f"{tmp_path}/result.tsv"
#     args.infile_col = 1
#     args.uniprot_cache = ""

#     # Run script
#     aln_query_uniprot_main(args)

#     # Parse the expected and observed into a dictionary of uniprot objects
#     expected_uniprot_object_dict = dict()
#     with open(
#         "tests/test_data/foldseek_related/af2_foldseek_result_small_uniprot_lookup.tsv"
#     ) as infile:
#         for line in infile:
#             line = line.rstrip("\n").rstrip("\t").split("\t")
#             acc, geneName, fullName = line
#             uo = uniprot_object(acc)
#             uo.geneName = geneName
#             uo.fullName = fullName
#             expected_uniprot_object_dict[acc] = uo

#     observed_uniprot_object_dict = dict()
#     with open(args.uniprot_lookup_output) as infile:
#         for line in infile:
#             line = line.rstrip("\n").rstrip("\t").split("\t")
#             acc, geneName, fullName = line
#             uo = uniprot_object(acc)
#             uo.geneName = geneName
#             uo.fullName = fullName
#             observed_uniprot_object_dict[acc] = uo

#     compare_uniprot_objects(expected_uniprot_object_dict, observed_uniprot_object_dict)


# def test_aln_query_uniprot_main_write_cache(tmp_path):

#     # Define args
#     class args:
#         pass

#     args.infile = "tests/test_data/foldseek_related/af2_foldseek_result_small.m8"
#     args.uniprot_lookup_output = f"{tmp_path}/result.tsv"
#     args.infile_col = 1
#     args.uniprot_cache = f"{tmp_path}/new_cache.pkl"

#     # Run script
#     aln_query_uniprot_main(args)

#     # Parse the expected and observed into a dictionary of uniprot objects
#     expected_uniprot_object_dict = dict()
#     with open(
#         "tests/test_data/foldseek_related/af2_foldseek_result_small_uniprot_lookup.tsv"
#     ) as infile:
#         for line in infile:
#             line = line.rstrip("\n").rstrip("\t").split("\t")
#             acc, geneName, fullName = line
#             uo = uniprot_object(acc)
#             uo.geneName = geneName
#             uo.fullName = fullName
#             expected_uniprot_object_dict[acc] = uo

#     observed_uniprot_object_dict = dict()
#     with open(args.uniprot_lookup_output) as infile:
#         for line in infile:
#             line = line.rstrip("\n").rstrip("\t").split("\t")
#             acc, geneName, fullName = line
#             uo = uniprot_object(acc)
#             uo.geneName = geneName
#             uo.fullName = fullName
#             observed_uniprot_object_dict[acc] = uo

#     compare_uniprot_objects(expected_uniprot_object_dict, observed_uniprot_object_dict)


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

    compare_uniprot_objects(expected_uniprot_object_dict, observed_uniprot_object_dict)
