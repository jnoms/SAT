from sat.scripts.utils.alignments import Alignment_group, Alignment_object
from sat.scripts.filter_alignments import filter_alignments_main
from sat.scripts.utils.Foldseek_Dataset import Foldseek_Dataset


def test_filter_alignments_func():
    aln1 = Alignment_object(["q1", "4.2E-01"], ["query", "alntmscore"])
    aln2 = Alignment_object(["q2", "1.2E-01"], ["query", "alntmscore"])
    aln3 = Alignment_object(["q3", "6.0E-01"], ["query", "alntmscore"])
    aln4 = Alignment_object(["q4", "0.1E-01"], ["query", "alntmscore"])
    aln5 = Alignment_object(["q5", "3.0E-01"], ["query", "alntmscore"])

    aln_group = Alignment_group("query")
    aln_group.alignments = [aln1, aln2, aln3, aln4, aln5]

    aln_group.filter_alignments("alntmscore", 1, 0.3)

    assert len(aln_group.alignments) == 3
    for q in ["q4", "q2"]:
        assert q not in [aln.query for aln in aln_group.alignments]


def test_keep_top_N_alignments():
    aln1 = Alignment_object(["q1", "4.2E-01"], ["query", "alntmscore"])
    aln2 = Alignment_object(["q2", "1.2E-01"], ["query", "alntmscore"])
    aln3 = Alignment_object(["q3", "6.0E-01"], ["query", "alntmscore"])
    aln4 = Alignment_object(["q4", "0.1E-01"], ["query", "alntmscore"])
    aln5 = Alignment_object(["q5", "3.0E-01"], ["query", "alntmscore"])

    aln_group = Alignment_group("query")
    aln_group.alignments = [aln1, aln2, aln3, aln4, aln5]

    aln_group.keep_top_N_alignments("alntmscore", 3)

    assert len(aln_group.alignments) == 3
    for q in ["q1", "q3", "q5"]:
        assert q in [aln.query for aln in aln_group.alignments]


def test_filter_alignments_script(tmp_path):

    # Define inputs
    class args:
        pass

    args.alignment_file = (
        "tests/test_data/foldseek_related/top_query_per_cluster_tax.m8"
    )
    args.output_file = f"{tmp_path}/test_filtered.m8"

    args.alignment_fields = ""
    args.N = 10
    args.filter_field = "alntmscore"
    args.min_val_filter_field = 0.4
    args.max_val_filter_field = 1

    # Run script
    filter_alignments_main(args)

    # Parse expected and observed
    expected = Foldseek_Dataset()
    expected.parse_alignment(
        "tests/test_data/foldseek_related/top_query_per_cluster_tax.filtered.m8"
    )
    observed = Foldseek_Dataset()
    observed.parse_alignment(args.output_file)

    # Compare with expected file
    for query, alignment_group in observed.alignment_groups.items():
        assert query in expected.alignment_groups
        for alignment in alignment_group.alignments:
            assert alignment in expected.alignment_groups[query].alignments

            expected_alignment = [
                a for a in expected.alignment_groups[query].alignments if a == alignment
            ][0]
            assert alignment.__eq__(expected_alignment)

    # Also make sure all of the alntmscores are appropriate and N's are correct
    for _, alignment_group in observed.alignment_groups.items():
        assert len(alignment_group.alignments) <= args.N
        for aln in alignment_group.alignments:
            assert float(aln.alntmscore) >= args.min_val_filter_field
