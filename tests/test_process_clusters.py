from sat.scripts.utils.alignments import Alignment_group, Alignment_object
from sat.scripts.utils.clusters import Cluster

aln1 = Alignment_object(["seq1", "seq2", "0.9"], ["query", "target", "alntmscore"])
aln2 = Alignment_object(["seq1", "seq3", "0.8"], ["query", "target", "alntmscore"])
aln3 = Alignment_object(["seq1", "seq4", "0.9"], ["query", "target", "alntmscore"])
aln4 = Alignment_object(["seq1", "seq1", "1"], ["query", "target", "alntmscore"])

aln_group1 = Alignment_group("seq1")
aln_group1.add_alignment(aln1)
aln_group1.add_alignment(aln2)
aln_group1.add_alignment(aln3)
aln_group1.add_alignment(aln4)

aln5 = Alignment_object(["seq5", "seq6", "0.9"], ["query", "target", "alntmscore"])
aln6 = Alignment_object(["seq5", "seq7", "0.1"], ["query", "target", "alntmscore"])
aln7 = Alignment_object(["seq5", "seq8", "0.1"], ["query", "target", "alntmscore"])
aln8 = Alignment_object(["seq5", "seq9", "0.1"], ["query", "target", "alntmscore"])

aln_group2 = Alignment_group("seq5")
aln_group2.add_alignment(aln5)
aln_group2.add_alignment(aln6)
aln_group2.add_alignment(aln7)
aln_group2.add_alignment(aln8)

aln_group3 = Alignment_group("seq5")
aln_group3.add_alignment(aln5)
aln_group3.add_alignment(aln6)
aln_group3.add_alignment(aln7)

# ------------------------------------------------------------------#
# CLUSTER1 - seq5 wins because it has more alignments
# ------------------------------------------------------------------#
cluster1 = Cluster("seq1")
cluster1.alignment_groups.append(aln_group1)
cluster1.alignment_groups.append(aln_group2)
for i in range(9):
    i = i + 1
    cluster1.cluster_members.add(f"seq{i}")


def test_Cluster_add_top_query_more_alignments():
    cluster1.add_top_query()
    assert cluster1.top_query == "seq5"


# ------------------------------------------------------------------#
# CLUSTER2 - seq1 wins because higher average TM score
# ------------------------------------------------------------------#
cluster2 = Cluster("seq1")
cluster2.alignment_groups.append(aln_group1)
cluster2.alignment_groups.append(aln_group3)
for i in range(6):
    i = i + 1
    cluster2.cluster_members.add(f"seq{i}")


def test_Cluster_add_top_query_higher_avg_TM():
    cluster2.add_top_query()
    assert cluster2.top_query == "seq1"


# ------------------------------------------------------------------#
# CLUSTER3 - no alignment_group had alignments. So,top_query should
# be equal to just the cluster id
# ------------------------------------------------------------------#
aln_group4 = Alignment_group("seq10")
cluster3 = Cluster("seq10")
cluster3.alignment_groups.append(aln_group4)


def test_Cluster_add_top_query_no_alignments():
    cluster3.add_top_query()
    assert cluster3.top_query == "seq10"


# ------------------------------------------------------------------#
# Confirm that all of the top query's alignments have been written out
# appropriately.
# ------------------------------------------------------------------#
def test_Cluster_write_top_query_alignments():
    cluster1.add_top_query()
    cluster2.add_top_query()
    cluster3.add_top_query()

    cluster1_expected = [
        "seq5\tseq6\t0.9",
        "seq5\tseq7\t0.1",
        "seq5\tseq8\t0.1",
        "seq5\tseq9\t0.1",
    ]
    cluster1_observed = (
        cluster1.write_top_query_alignments(["query", "target", "alntmscore"])
        .rstrip("\n")
        .split("\n")
    )
    assert cluster1_observed == cluster1_expected


def test_Cluster_write_all_nonredundant_alignments():
    cluster1_expected = [
        "seq1\tseq2\t0.9",
        "seq1\tseq3\t0.8",
        "seq1\tseq4\t0.9",
        "seq5\tseq6\t0.9",
        "seq5\tseq7\t0.1",
        "seq5\tseq8\t0.1",
        "seq5\tseq9\t0.1",
    ]
    cluster1_observed = (
        cluster1.write_all_nonredundant_alignments(["query", "target", "alntmscore"])
        .rstrip("\n")
        .split("\n")
    )

    assert cluster1_observed == cluster1_expected
