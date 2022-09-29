from sat.scripts.add_uniprot_information_to_alignments import (
    add_uniprot_information_to_alignments_main,
)
from sat.scripts.utils.Foldseek_Dataset import Foldseek_Dataset


def test_add_uniprot_information_to_alignments(tmp_path):

    # Inputs
    class args:
        pass

    args.alignment_file = (
        "tests/test_data/foldseek_related/af2_foldseek_result_small.m8"
    )
    args.uniprot_information = (
        "tests/test_data/foldseek_related/af2_foldseek_result_small_uniprot_lookup.tsv"
    )
    args.outfile = f"{tmp_path}/af2_foldseek_result_small.uniprot.m8"
    args.alignment_fields = (
        "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,"
        "bits,alntmscore,taxid"
    )

    # Running script
    add_uniprot_information_to_alignments_main(args)

    # Parse observed and expected
    expected = Foldseek_Dataset()
    expected.parse_alignment(
        "tests/test_data/foldseek_related/af2_foldseek_result_small.uniprot.m8"
    )
    observed = Foldseek_Dataset()
    observed.parse_alignment(args.outfile)

    # Compare them
    for query, alignment_group in expected.alignment_groups.items():
        assert query in observed.alignment_groups
        for alignment in alignment_group.alignments:
            assert alignment in observed.alignment_groups[query].alignments

            observed_alignment = [
                a for a in expected.alignment_groups[query].alignments if a == alignment
            ][0]
            assert alignment.__eq__(observed_alignment)
