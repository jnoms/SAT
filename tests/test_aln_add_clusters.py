from sat.scripts.utils.Foldseek_Dataset import (
    Foldseek_Dataset,
    compare_foldseek_datasets,
)
from sat.scripts.aln_add_clusters import (
    aln_add_clusters_main,
    get_percentage_of_members_with_alignments,
    Super_cluster,
)


def test_get_percentage_of_members_with_alignments__all_align():

    members_1 = set(["a", "b"])
    members_2 = set(["e", "f", "g"])
    alignment_dict = {
        "a": set("e"),
        "b": set(["f"]),
        "c": set(["g", "f"]),
        "d": set(["f"]),
        "e": set("b"),
        "f": set("c"),
    }

    observed = get_percentage_of_members_with_alignments(
        members_1, members_2, alignment_dict
    )
    expected = 1.0

    assert observed == expected


def test_get_percentage_of_members_with_alignments__none_align():

    members_1 = set(["a", "b"])
    members_2 = set(["e", "f", "g"])
    alignment_dict = {
        "a": set("z"),
        "b": set(["y"]),
        "c": set(["x", "q"]),
        "d": set(["f"]),
        "e": set("b"),
        "f": set("c"),
    }

    observed = get_percentage_of_members_with_alignments(
        members_1, members_2, alignment_dict
    )
    expected = 0

    assert observed == expected


def test_get_percentage_of_members_with_alignments__overlapping_members():

    members_1 = set(["a", "b"])
    members_2 = set(["a", "e", "f", "g"])
    alignment_dict = {
        "a": set("e"),
        "b": set(["f"]),
        "c": set(["g", "f"]),
        "d": set(["f"]),
        "e": set("b"),
        "f": set("c"),
    }

    observed = get_percentage_of_members_with_alignments(
        members_1, members_2, alignment_dict
    )
    expected = 1.0

    assert observed == expected


def test_get_percentage_of_members_with_alignments__half_align():

    members_1 = set(["a", "b"])
    members_2 = set(["a", "e", "f", "g"])
    alignment_dict = {
        "a": set("e"),
        "b": set(["z"]),
        "c": set(["g", "f"]),
        "d": set(["f"]),
        "e": set("b"),
        "f": set("c"),
    }

    observed = get_percentage_of_members_with_alignments(
        members_1, members_2, alignment_dict
    )
    expected = 0.5

    assert observed == expected


def test_all_are_linked__all_linked():
    one = Super_cluster("a")
    one.add_member("b")
    one.add_member("c")

    two = Super_cluster("x")
    two.add_member("y")
    two.add_member("z")

    cluster_linkages = {
        "a": set(["x", "y", "z", "q"]),
        "b": set(["x", "y", "z", "w"]),
        "c": set(["x", "y", "z", "e"]),
        "x": set(["a", "b", "c", "r"]),
        "y": set(["a", "b", "c", "t"]),
        "z": set(["a", "b", "c", "y"]),
    }

    assert one.all_are_linked(two, cluster_linkages)


def test_all_are_linked__only_not_reciprocal():

    one = Super_cluster("a")
    one.add_member("b")
    one.add_member("c")

    two = Super_cluster("x")
    two.add_member("y")
    two.add_member("z")

    # HERE - x from sc "two" is not linked with b and c from sc "one"
    cluster_linkages = {
        "a": set(["x", "y", "z", "q"]),
        "b": set(["x", "y", "z", "w"]),
        "c": set(["x", "y", "z", "e"]),
        "x": set(["a", "r"]),
        "y": set(["a", "b", "c", "t"]),
        "z": set(["a", "b", "c", "y"]),
    }

    assert not one.all_are_linked(two, cluster_linkages)


def test_all_are_linked__none_linked():

    one = Super_cluster("a")
    one.add_member("b")
    one.add_member("c")

    two = Super_cluster("x")
    two.add_member("y")
    two.add_member("z")

    cluster_linkages = {
        "a": set(["q"]),
        "b": set(["w"]),
        "c": set(["e"]),
        "x": set(["r"]),
        "y": set(["t"]),
        "z": set(["y"]),
    }
    assert not one.all_are_linked(two, cluster_linkages)


def test_all_are_linked__same_members():

    one = Super_cluster("a")
    one.add_member("b")
    one.add_member("c")

    cluster_linkages = {}

    assert one.all_are_linked(one, cluster_linkages)


# ------------------------------------------------------------------------------------ #
# Whole-script tests
# ------------------------------------------------------------------------------------ #
# def test_aln_add_clusters_rep_out(tmp_path):

#     # Define inputs
#     class args:
#         pass

#     args.alignment_file = "tests/test_data/foldseek_related/clusters_alignment.m8"
#     args.cluster_file = "tests/test_data/foldseek_related/clusters.tsv"
#     args.rep_out = f"{tmp_path}/rep_out.m8"
#     args.alignment_fields = (
#         "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,"
#         "evalue,bits,alntmscore"
#     ).split(",")
#     args.all_out = ""

#     # Run the script
#     aln_add_clusters_main(args)

#     # Validate the outputs
#     expected = Foldseek_Dataset()
#     expected.parse_alignment("tests/test_data/foldseek_related/rep_out.m8")

#     observed = Foldseek_Dataset()
#     observed.parse_alignment(args.rep_out)

#     compare_foldseek_datasets(expected, observed)


# def test_aln_add_clusters_all(tmp_path):

#     # Define inputs
#     class args:
#         pass

#     args.alignment_file = "tests/test_data/foldseek_related/clusters_alignment.m8"
#     args.cluster_file = "tests/test_data/foldseek_related/clusters.tsv"
#     args.all_out = f"{tmp_path}/all_out.m8"
#     args.rep_out = ""
#     args.alignment_fields = (
#         "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,"
#         "evalue,bits,alntmscore"
#     ).split(",")

#     # Run the script
#     aln_add_clusters_main(args)

#     # Read in expected and observed
#     expected = Foldseek_Dataset()
#     expected.parse_alignment("tests/test_data/foldseek_related/all_out.m8")

#     observed = Foldseek_Dataset()
#     observed.parse_alignment(args.all_out)

#     compare_foldseek_datasets(expected, observed)

#     # Also compare the input and observed - should have pretty much the same alignments
#     # except for removal of self-self alignments
#     infile_aln_set = set()
#     with open(args.alignment_file) as infile:
#         for line in infile:
#             line = line.rstrip("\n").split("\t")
#             if line[0] == line[1]:
#                 continue
#             if line[0] == "query":
#                 continue
#             aln = (line[0], line[1])
#             infile_aln_set.add(aln)

#     observed_aln_set = set()
#     with open(args.all_out) as infile:
#         for line in infile:
#             line = line.rstrip("\n").split("\t")
#             if line[0] == line[1]:
#                 continue
#             if line[0] == "query":
#                 continue
#             aln = (line[0], line[1])
#             observed_aln_set.add(aln)

#     assert infile_aln_set == observed_aln_set


# def test_aln_add_clusters_rep_out_some_queries_removed(tmp_path):

#     # Define inputs
#     class args:
#         pass

#     args.alignment_file = (
#         "tests/test_data/foldseek_related/clusters_alignment_some_queries_removed.m8"
#     )
#     args.cluster_file = "tests/test_data/foldseek_related/clusters.tsv"
#     args.rep_out = f"{tmp_path}/rep_some_removed.m8"
#     args.alignment_fields = (
#         "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,"
#         "evalue,bits,alntmscore"
#     ).split(",")
#     args.all_out = ""

#     # Run the script
#     aln_add_clusters_main(args)
#     fields = args.alignment_fields + [
#         "cluster_ID",
#         "cluster_count",
#         "cluster_rep",
#     ]

#     # Validate the outputs
#     expected = Foldseek_Dataset()
#     expected.parse_alignment(
#         "tests/test_data/foldseek_related/rep_some_removed.m8",
#         fields,
#     )

#     observed = Foldseek_Dataset()
#     observed.parse_alignment(args.rep_out, fields)

#     compare_foldseek_datasets(expected, observed)
