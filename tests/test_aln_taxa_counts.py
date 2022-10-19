from sat.scripts.aln_taxa_counts import aln_taxa_counts_main


def test_aln_taxa_counts_main(tmp_path):

    # Define inputs
    class args:
        pass

    args.alignment_file = (
        "tests/test_data/foldseek_related/top_query_per_cluster_tax.m8"
    )
    args.output_file = f"{tmp_path}/top_query_per_cluster_tax_counts.tsv"
    args.taxonomy_levels = "superkingdom,phylum,class,order,family,genus,species"

    # Run script
    aln_taxa_counts_main(args)

    # Read observed and expected into a very nested dictionary
    expected = dict()
    with open(
        "tests/test_data/foldseek_related/top_query_per_cluster_tax_counts.tsv"
    ) as infile:
        for line in infile:
            line = line.rstrip("\n").split("\t")
            cluster_ID, cluster_rep, level, superkingdom, taxon, count = line

            if cluster_ID not in expected:
                expected[cluster_ID] = dict()

            if level not in expected[cluster_ID]:
                expected[cluster_ID][level] = dict()

            if taxon not in expected[cluster_ID][level]:
                expected[cluster_ID][level][taxon] = count

    observed = dict()
    with open(args.output_file) as infile:
        for line in infile:
            line = line.rstrip("\n").split("\t")
            cluster_ID, cluster_rep, level, superkingdom, taxon, count = line

            if cluster_ID not in observed:
                observed[cluster_ID] = dict()

            if level not in observed[cluster_ID]:
                observed[cluster_ID][level] = dict()

            if taxon not in observed[cluster_ID][level]:
                observed[cluster_ID][level][taxon] = count

    # Compare them
    assert expected == observed
