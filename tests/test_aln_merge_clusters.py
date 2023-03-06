# from sat.scripts.aln_merge_clusters import aln_merge_clusters_main
# from sat.scripts.utils.Foldseek_Dataset import (
#     Foldseek_Dataset,
#     compare_foldseek_datasets,
# )

from sat.scripts.aln_merge_clusters import get_superclusters_no_reciprocal


def test_get_superclusters_no_reciprocal__all_linked():
    linkage_dict = dict()
    linkage_dict["1"] = set(["2", "3"])
    linkage_dict["2"] = set(["2", "3"])
    linkage_dict["3"] = set(["2", "3"])

    observed = get_superclusters_no_reciprocal(linkage_dict)
    expected = [set(["1", "2", "3"])]
    assert expected == observed


def test_get_superclusters_no_reciprocal__two_linked():
    linkage_dict = dict()
    linkage_dict["1"] = set(["1", "2"])
    linkage_dict["2"] = set(["1", "2"])
    linkage_dict["3"] = set(["6"])

    observed = get_superclusters_no_reciprocal(linkage_dict)
    print(observed)
    expected = [set(["1", "2"]), set(["3", "6"])]
    assert expected == observed
