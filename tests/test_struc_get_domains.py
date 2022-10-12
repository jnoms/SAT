from sat.scripts.struc_get_domains import (
    filter_clusters,
    smooth_array,
    plddt_trim_clusters,
)
import numpy as np


def test_filter_clusters_two_good():

    # Two clusters are all good
    cluster_coords = [(1, 2, 3, 4), (6, 7, 8, 9, 10)]  # 0 indexed
    plddt_array = np.array([70, 70, 70, 70, 70, 70, 60, 70, 50, 90])
    desired = [(1, 2, 3, 4), (6, 7, 8, 9, 10)]
    observed = filter_clusters(cluster_coords, plddt_array, 3, 60)
    assert desired == observed


def test_filter_clusters_one_below_plddt():
    # One cluster below plddt
    cluster_coords = [(1, 2, 3, 4), (6, 7, 8, 9, 10)]  # 0 indexed
    plddt_array = np.array([10, 10, 10, 70, 70, 70, 60, 70, 50, 90])
    desired = [(6, 7, 8, 9, 10)]
    observed = filter_clusters(cluster_coords, plddt_array, 3, 60)
    assert desired == observed


def test_filter_clusters_both_below_plddt():
    # Both clusters below plddt
    cluster_coords = [(1, 2, 3, 4), (6, 7, 8, 9, 10)]  # 0 indexed
    plddt_array = np.array([10, 10, 10, 70, 70, 70, 60, 10, 10, 100])
    desired = []
    observed = filter_clusters(cluster_coords, plddt_array, 3, 60)
    assert desired == observed


def test_filter_clusters_three_clust_one_filtered():
    # three clusters - one filtered
    cluster_coords = [(1, 2, 3, 4), (6, 7, 8, 9, 10), (12, 13, 14, 15)]  # 0 indexed
    plddt_array = np.array(
        [10, 70, 150, 70, 70, 70, 60, 10, 10, 100, 70, 70, 70, 60, 60]
    )
    desired = [(1, 2, 3, 4), (12, 13, 14, 15)]
    observed = filter_clusters(cluster_coords, plddt_array, 3, 60)
    assert desired == observed


def test_smooth_array_simple_case():
    in_arr = [10, 10, 3, 3, 3, 10, 3, 3, 10]
    in_arr = np.array(in_arr)
    expected = [10, 10, 1, 1, 1, 1, 1, 1, 10]
    assert list(smooth_array(in_arr)) == expected


def test_smooth_array_no_smoothing():
    in_arr = [10, 10, 6, 6, 6, 10, 6, 6, 10]
    in_arr = np.array(in_arr)
    expected = [10, 10, 6, 6, 6, 10, 6, 6, 10]
    assert list(smooth_array(in_arr)) == expected


def test_smooth_array_smooth_start():
    in_arr = [2, 2, 4, 4, 4, 10, 2, 4, 10]
    in_arr = np.array(in_arr)
    expected = [1, 1, 1, 1, 1, 1, 1, 1, 10]
    assert list(smooth_array(in_arr)) == expected


def test_smooth_array_smooth_end():
    in_arr = [10, 10, 4, 4, 4, 10, 2, 2, 2]
    in_arr = np.array(in_arr)
    expected = [10, 10, 1, 1, 1, 1, 1, 1, 1]
    assert list(smooth_array(in_arr)) == expected


def test_smooth_array_n_long_enough():
    in_arr = [10, 10, 3, 3, 3, 10, 10, 10, 3, 3, 10]
    in_arr = np.array(in_arr)
    expected = [10, 10, 1, 1, 1, 10, 10, 10, 1, 1, 10]
    assert list(smooth_array(in_arr, n=2)) == expected


def test_plddt_trim_clusters_no_trimming():
    clusters = [frozenset([1, 2, 3, 4, 5]), frozenset([6, 7, 8, 9, 10])]
    plddt_array = [70 for i in range(10)]
    plddt_array = np.array(plddt_array)
    expected = [frozenset([1, 2, 3, 4, 5]), frozenset([6, 7, 8, 9, 10])]
    assert plddt_trim_clusters(clusters, plddt_array, min_avg_plddt=60) == expected


def test_plddt_trim_clusters_everything_trimmed():
    clusters = [frozenset([1, 2, 3, 4, 5]), frozenset([6, 7, 8, 9, 10])]
    plddt_array = [10 for i in range(10)]
    plddt_array = np.array(plddt_array)
    expected = [frozenset(), frozenset()]
    assert plddt_trim_clusters(clusters, plddt_array, min_avg_plddt=60) == expected


def test_plddt_trim_clusters_N_trimmed():
    clusters = [frozenset([1, 2, 3, 4, 5]), frozenset([6, 7, 8, 9, 10])]
    plddt_array = [10 for i in range(4)] + [70 for i in range(6)]
    plddt_array = np.array(plddt_array)
    expected = [frozenset([5]), frozenset([6, 7, 8, 9, 10])]
    assert plddt_trim_clusters(clusters, plddt_array, min_avg_plddt=60) == expected


def test_plddt_trim_clusters_C_trimmed():
    clusters = [frozenset([1, 2, 3, 4, 5]), frozenset([6, 7, 8, 9, 10])]
    plddt_array = [60 for i in range(4)] + [10 for i in range(6)]
    plddt_array = np.array(plddt_array)
    expected = [frozenset([1, 2, 3, 4]), frozenset()]
    assert plddt_trim_clusters(clusters, plddt_array, min_avg_plddt=60) == expected
