from sat.scripts.utils.Foldseek_Dataset import (
    Foldseek_Dataset,
    compare_foldseek_datasets,
)
from sat.scripts.aln_add_clusters import aln_add_clusters_main


# ------------------------------------------------------------------------------------ #
# Whole-script tests
# ------------------------------------------------------------------------------------ #
def test_aln_add_clusters_rep_out(tmp_path):

    # Define inputs
    class args:
        pass

    args.alignment_file = "tests/test_data/foldseek_related/clusters_alignment.m8"
    args.cluster_file = "tests/test_data/foldseek_related/clusters.tsv"
    args.rep_out = f"{tmp_path}/rep_out.m8"
    args.alignment_fields = (
        "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,"
        "evalue,bits,alntmscore"
    ).split(",")
    args.all_out = ""

    # Run the script
    aln_add_clusters_main(args)

    # Validate the outputs
    expected = Foldseek_Dataset()
    expected.parse_alignment("tests/test_data/foldseek_related/rep_out.m8")

    observed = Foldseek_Dataset()
    observed.parse_alignment(args.rep_out)

    compare_foldseek_datasets(expected, observed)


def test_aln_add_clusters_all(tmp_path):

    # Define inputs
    class args:
        pass

    args.alignment_file = "tests/test_data/foldseek_related/clusters_alignment.m8"
    args.cluster_file = "tests/test_data/foldseek_related/clusters.tsv"
    args.all_out = f"{tmp_path}/all_out.m8"
    args.rep_out = ""
    args.alignment_fields = (
        "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,"
        "evalue,bits,alntmscore"
    ).split(",")

    # Run the script
    aln_add_clusters_main(args)

    # Read in expected and observed
    expected = Foldseek_Dataset()
    expected.parse_alignment("tests/test_data/foldseek_related/all_out.m8")

    observed = Foldseek_Dataset()
    observed.parse_alignment(args.all_out)

    compare_foldseek_datasets(expected, observed)

    # Also compare the input and observed - should have pretty much the same alignments
    # except for removal of self-self alignments
    infile_aln_set = set()
    with open(args.alignment_file) as infile:
        for line in infile:
            line = line.rstrip("\n").split("\t")
            if line[0] == line[1]:
                continue
            if line[0] == "query":
                continue
            aln = (line[0], line[1])
            infile_aln_set.add(aln)

    observed_aln_set = set()
    with open(args.all_out) as infile:
        for line in infile:
            line = line.rstrip("\n").split("\t")
            if line[0] == line[1]:
                continue
            if line[0] == "query":
                continue
            aln = (line[0], line[1])
            observed_aln_set.add(aln)

    assert infile_aln_set == observed_aln_set


def test_aln_add_clusters_rep_out_some_queries_removed(tmp_path):

    # Define inputs
    class args:
        pass

    args.alignment_file = (
        "tests/test_data/foldseek_related/clusters_alignment_some_queries_removed.m8"
    )
    args.cluster_file = "tests/test_data/foldseek_related/clusters.tsv"
    args.rep_out = f"{tmp_path}/rep_some_removed.m8"
    args.alignment_fields = (
        "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,"
        "evalue,bits,alntmscore"
    ).split(",")
    args.all_out = ""

    # Run the script
    aln_add_clusters_main(args)
    fields = args.alignment_fields + [
        "cluster_ID",
        "cluster_count",
        "cluster_rep",
    ]

    # Validate the outputs
    expected = Foldseek_Dataset()
    expected.parse_alignment(
        "tests/test_data/foldseek_related/rep_some_removed.m8",
        fields,
    )

    observed = Foldseek_Dataset()
    observed.parse_alignment(args.rep_out, fields)

    compare_foldseek_datasets(expected, observed)
