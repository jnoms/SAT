from sat.scripts.struc_get_domains import (
    filter_clusters,
    smooth_array,
    plddt_trim_clusters,
    struc_get_domains_main,
)
from sat.scripts.utils.structure import pdb_to_structure_object, compare_structures
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


def test_struc_get_domains_three_domain_structure(tmp_path):
    # Inputs
    class args:
        pass

    args.structure_file_path = (
        "tests/test_data/structure_related/get_domains/inputs/three_domains.pdb"
    )
    args.pae_path = (
        "tests/test_data/structure_related/get_domains/inputs/three_domains.scores.json"
    )
    args.output_dir = f"{tmp_path}/domains"
    args.pae_power = 1
    args.pae_cutoff = 5
    args.graph_resolution = 1
    args.min_domain_length = 50
    args.min_domain_plddt = 60
    args.smooth_n = 0
    args.plddt_report = f"{tmp_path}/plddt_report.tsv"
    args.pae_report = f"{tmp_path}/pae_report.tsv"

    # Run script
    struc_get_domains_main(args)

    # Compare outputs - domains
    expected_domain_1 = pdb_to_structure_object(
        "tests/test_data/structure_related/get_domains/outputs/three_domains_domain-1.pdb"
    )
    expected_domain_2 = pdb_to_structure_object(
        "tests/test_data/structure_related/get_domains/outputs/three_domains_domain-2.pdb"
    )
    expected_domain_3 = pdb_to_structure_object(
        "tests/test_data/structure_related/get_domains/outputs/three_domains_domain-3.pdb"
    )

    observed_domain_1 = pdb_to_structure_object(
        f"{tmp_path}/domains/three_domains_domain-1.pdb"
    )
    observed_domain_2 = pdb_to_structure_object(
        f"{tmp_path}/domains/three_domains_domain-2.pdb"
    )
    observed_domain_3 = pdb_to_structure_object(
        f"{tmp_path}/domains/three_domains_domain-3.pdb"
    )

    compare_structures(expected_domain_1, observed_domain_1)
    compare_structures(expected_domain_2, observed_domain_2)
    compare_structures(expected_domain_3, observed_domain_3)

    # Check the pae and plddt reports
    expected_plddt_report = open(
        "tests/test_data/structure_related/get_domains/outputs/three_domains_plddt_report.tsv"
    )
    observed_plddt_report = open(args.plddt_report)
    assert expected_plddt_report.read() == observed_plddt_report.read()
    expected_plddt_report.close()
    observed_plddt_report.close()

    # Check the pae and plddt reports
    expected_pae_report = open(
        "tests/test_data/structure_related/get_domains/outputs/three_domains_pae_report.tsv"
    )
    observed_pae_report = open(args.pae_report)
    assert expected_pae_report.read() == observed_pae_report.read()
    expected_pae_report.close()
    observed_pae_report.close()


def test_struc_get_domains_whole_structure_used_as_domain(tmp_path):
    # Inputs
    class args:
        pass

    args.structure_file_path = "tests/test_data/structure_related/get_domains/inputs/full_struc_used_as_domain.pdb"
    args.pae_path = "tests/test_data/structure_related/get_domains/inputs/full_struc_used_as_domain.scores.json"
    args.output_dir = f"{tmp_path}/domains"
    args.pae_power = 1
    args.pae_cutoff = 5
    args.graph_resolution = 1
    args.min_domain_length = 50
    args.min_domain_plddt = 60
    args.smooth_n = 0
    args.plddt_report = f"{tmp_path}/plddt_report.tsv"
    args.pae_report = f"{tmp_path}/pae_report.tsv"

    # Run script
    struc_get_domains_main(args)

    # Compare outputs - domains
    expected_domain_1 = pdb_to_structure_object(
        "tests/test_data/structure_related/get_domains/outputs/full_struc_used_as_domain_domain-1.pdb"
    )
    observed_domain_1 = pdb_to_structure_object(
        f"{tmp_path}/domains/full_struc_used_as_domain_domain-1.pdb"
    )
    compare_structures(expected_domain_1, observed_domain_1)

    # Check the pae and plddt reports
    expected_plddt_report = open(
        "tests/test_data/structure_related/get_domains/outputs/full_struc_used_as_domain_plddt_report.tsv"
    )
    observed_plddt_report = open(args.plddt_report)
    assert expected_plddt_report.read() == observed_plddt_report.read()
    expected_plddt_report.close()
    observed_plddt_report.close()

    # Check the pae and plddt reports
    expected_pae_report = open(
        "tests/test_data/structure_related/get_domains/outputs/full_struc_used_as_domain_pae_report.tsv"
    )
    observed_pae_report = open(args.pae_report)
    assert expected_pae_report.read() == observed_pae_report.read()
    expected_pae_report.close()
    observed_pae_report.close()
