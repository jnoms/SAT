from sat.scripts.get_domains import filter_clusters
import numpy as np


def test_filter_clusters():

    # Two clusters are all good
    cluster_coords = [(1, 2, 3, 4), (6, 7, 8, 9, 10)]  # 0 indexed
    plddt_array = np.array([70, 70, 70, 70, 70, 70, 60, 70, 50, 90])
    desired = [(1, 2, 3, 4), (6, 7, 8, 9, 10)]
    observed = filter_clusters(cluster_coords, plddt_array, 3, 60)
    assert desired == observed

    # One cluster below plddt
    cluster_coords = [(1, 2, 3, 4), (6, 7, 8, 9, 10)]  # 0 indexed
    plddt_array = np.array([10, 10, 10, 70, 70, 70, 60, 70, 50, 90])
    desired = [(6, 7, 8, 9, 10)]
    observed = filter_clusters(cluster_coords, plddt_array, 3, 60)
    assert desired == observed

    # Both clusters below plddt
    cluster_coords = [(1, 2, 3, 4), (6, 7, 8, 9, 10)]  # 0 indexed
    plddt_array = np.array([10, 10, 10, 70, 70, 70, 60, 10, 10, 100])
    desired = []
    observed = filter_clusters(cluster_coords, plddt_array, 3, 60)
    assert desired == observed

    # three clusters - one filtered
    cluster_coords = [(1, 2, 3, 4), (6, 7, 8, 9, 10), (12, 13, 14, 15)]  # 0 indexed
    plddt_array = np.array(
        [10, 70, 150, 70, 70, 70, 60, 10, 10, 100, 70, 70, 70, 60, 60]
    )
    desired = [(1, 2, 3, 4), (12, 13, 14, 15)]
    observed = filter_clusters(cluster_coords, plddt_array, 3, 60)
    assert desired == observed
