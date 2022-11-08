from sat.scripts.aln_add_taxonomy import aln_add_taxonomy_main
from sat.scripts.utils.Foldseek_Dataset import (
    Foldseek_Dataset,
    compare_foldseek_datasets,
)


def test_aln_add_taxonomy_main(tmp_path):

    # Define inputs
    class args:
        pass

    args.alignment_file = "tests/test_data/foldseek_related/rep_out.m8"
    args.output_file = f"{tmp_path}/rep_out.tax.m8"
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
    expected = Foldseek_Dataset()
    expected.parse_alignment("tests/test_data/foldseek_related/rep_out.tax.m8")

    observed = Foldseek_Dataset()
    observed.parse_alignment(args.output_file)

    # Check content
    compare_foldseek_datasets(expected, observed)
