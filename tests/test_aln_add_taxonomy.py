from sat.scripts.aln_add_taxonomy import aln_add_taxonomy_main
from sat.scripts.utils.Foldseek_Dataset import Foldseek_Dataset


def test_aln_add_taxonomy_main(tmp_path):

    # Define inputs
    class args:
        pass

    args.alignment_file = "tests/test_data/foldseek_related/top_query_per_cluster.m8"
    args.output_file = f"{tmp_path}/top_query_per_cluster_tax.m8"
    args.alignment_fields = (
        "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,"
        "evalue,bits,alntmscore,cluster_ID,cluster_count,top_query"
    ).split(",")
    args.taxonID_finder_delimiter = "__"
    args.taxonID_finder_pos = -1
    args.taxonomy_levels = "superkingdom,phylum,class,order,family,genus,species"

    # Run program
    aln_add_taxonomy_main(args)

    # Compare result with expected
    taxonomy_cols = [f"query_{level}" for level in args.taxonomy_levels]
    taxonomy_cols += [f"target_{level}" for level in args.taxonomy_levels]
    taxonomy_cols += ["query_taxonID", "target_taxonID"]
    result_fields = args.alignment_fields + taxonomy_cols

    expected = Foldseek_Dataset()
    expected.parse_alignment(
        "tests/test_data/foldseek_related/top_query_per_cluster_tax.m8",
        result_fields,
    )

    observed = Foldseek_Dataset()
    observed.parse_alignment(args.output_file, result_fields)

    # Check content
    for query, alignment_group in expected.alignment_groups.items():
        assert query in observed.alignment_groups
        for alignment in alignment_group.alignments:
            assert alignment in observed.alignment_groups[query].alignments

            observed_alignment = [
                a for a in expected.alignment_groups[query].alignments if a == alignment
            ][0]
            assert alignment.__eq__(observed_alignment)
