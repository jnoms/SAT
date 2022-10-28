# ------------------------------------------------------------------------------------ #
# Import dependencies
# ------------------------------------------------------------------------------------ #

import json
import numpy as np
import networkx as nx
from networkx.algorithms import community
from itertools import groupby
import matplotlib.pyplot as plt
import os

from .utils.misc import make_output_dir, talk_to_me
from .utils.structure import pdb_to_structure_object, write_structure_subset

# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #


def parse_json_file(pae_json_file):
    """
    This function was adapted from
    https://github.com/tristanic/pae_to_domains/blob/main/pae_to_domains.py
    which is under the MIT license.
    """
    # Returns json information as a dictionary
    with open(pae_json_file, "rt") as f:
        data = json.load(f)

    # Colabfold format only
    if "pae" not in data.keys() or "plddt" not in data.keys():
        msg = "This script is only configured for the colabfold PAE format."
        raise ValueError(msg)

    pae_matrix = np.array(data["pae"], dtype=np.float64)
    plddt_array = data["plddt"]

    return pae_matrix, plddt_array


def domains_from_pae_matrix_networkx(
    pae_matrix, pae_power=1, pae_cutoff=5, graph_resolution=1
):
    """
    This function was adapted from
    https://github.com/tristanic/pae_to_domains/blob/main/pae_to_domains.py
    under the MIT license. Because this function is used as-is, this function
    is still liscenced under the original MIT license of that source. See
    https://github.com/tristanic/pae_to_domains/blob/main/LICENSE.

    Takes a predicted aligned error (PAE) matrix representing the predicted error in
    distances between each pair of residues in a model, and uses a graph-based community
    clustering algorithm to partition the model into approximately rigid groups.
    Arguments:
        * pae_matrix: a (n_residues x n_residues) numpy array. Diagonal elements should
        be set to some non-zero value to avoid divide-by-zero warnings
        * pae_power (optional, default=1): each edge in the graph will be weighted
        proportional to (1/pae**pae_power)
        * pae_cutoff (optional, default=5): graph edges will only be created for
        residue pairs with pae<pae_cutoff
        * graph_resolution (optional, default=1): regulates how aggressively the
        clustering algorithm is. Smaller values lead to larger clusters. Value should
        be larger than zero, and values larger than 5 are unlikely to be useful.
    Returns: a series of lists, where each list contains the indices of residues
    belonging to one cluster.
    """
    weights = 1 / pae_matrix**pae_power

    g = nx.Graph()
    size = weights.shape[0]
    g.add_nodes_from(range(size))
    edges = np.argwhere(pae_matrix < pae_cutoff)
    sel_weights = weights[edges.T[0], edges.T[1]]
    wedges = [(i, j, w) for (i, j), w in zip(edges, sel_weights)]
    g.add_weighted_edges_from(wedges)
    clusters = community.greedy_modularity_communities(
        g, weight="weight", resolution=graph_resolution
    )
    return clusters


def smooth_array(in_arr, threshold=5, n=20, block_replace=1):
    """
    This function takes in a 1 dimensional array of values, where lower values are
    'good'. It finds groups of values that are lower than the threshold value. If
    a 'good' group ends but the subequent 'bad' group (e.g. has values > threshold) is
    not long (is less long than n), it now considers that formerly-bad group to be
    good. Then, it overrides all good values as block_replace.

    In terms of taking in an array from a PAE matrix, this function:
    - finds regions of the vector that have a PAE below threshold
    - if two regions of PAE below threshold are separated by a region of PAE
      above threshold that is less long than n, the regions are merged.
    - all pae regions below threshold (including the converted, formerly-above-threshold
      regions) are converted to block_replace.
    """
    status = in_arr < threshold

    # (bool, count). Note that grouped_status is essentially 1-indexed.
    grouped_status = [(k, sum(1 for i in g)) for k, g in groupby(status)]

    pos_to_smooth = set()  # set with indices of positions to smooth
    i = 0
    for (status, length) in grouped_status:

        i += length

        # If the stretch is at the start, won't correct a False but will correct a
        # true.
        if i - length == 0 and not status:
            continue

        # If the stretch is at the end, won't correct a False but will correct a
        # true.
        if i == len(in_arr) and not status:
            continue

        # If true or if false but under length, will smooth it
        if status:
            pos_to_smooth.update(range(i - length, i))
        elif length < n and not status:
            pos_to_smooth.update(range(i - length, i))

    out_arr = []
    for i, val in enumerate(in_arr):
        if i in pos_to_smooth:
            out_arr.append(block_replace)
        else:
            out_arr.append(val)
    out_arr = np.array(out_arr)
    return out_arr


def visualize_PAE(pae_matrix):
    """
    Given a PAE matrix, returns a matplotlib heatmap of the residues.
    """
    plt.imshow(np.asarray(pae_matrix))
    plt.xlabel("Residue 1")
    plt.ylabel("Residue 2")
    return plt


def filter_clusters(clusters, plddt_array, min_length, min_avg_plddt):
    """
    Filters clusters to keep those that are longer than a minimum length and
    have an average pLDDT higher than a minimum value.

    clusters is a list, where each item is a cluster. Each cluster is represented as
    a FrozenSet of positions present in the domain. Notably, the coordinates are
    1-indexed!

    plddt_array is a numpy array wich the plddt at each residue.
    """

    filtered_clusters = []
    for cluster in clusters:

        # Length filter
        cluster_length = len(cluster)
        if cluster_length < min_length:
            continue

        # pLDDT filter
        plddt_sum = 0
        num_residues = 0
        for pos in cluster:

            # Because the positions in clusters are 1-indexed but I am subsetting a
            # vector, need to -1 to convert it to 0-indexed.
            plddt = plddt_array[pos - 1]
            plddt_sum += plddt
            num_residues += 1

        pLDDT_avg = plddt_sum / num_residues
        if pLDDT_avg < min_avg_plddt:
            continue

        # Sanity check
        if num_residues != cluster_length:
            msg = "Something went wrong - number of residues parsed for plddt"
            msg += " is not the same length as the cluster. Cluser in question is"
            msg += f"{cluster}"
            raise ValueError(msg)

        filtered_clusters.append(cluster)

    return filtered_clusters


def plddt_trim_clusters(clusters, plddt_array, min_avg_plddt):
    """
    clusters is a list of Frozensets, where each Frozenset contains the 1-indexed
    positions of residues belonging to a domain.

    Unlike filter_clusters (which just removes clusters with an average pLDDT below
    the threshold), this function trims residues on the N and C sides of the cluster
    that are below the threshold.

    Remember, the positions in clusters are 1-indexed! So need to -1 to get 0 index
    when looking up in the plddt_array.
    """
    out_clusters = []
    for cluster in clusters:
        cluster = list(cluster)
        to_trim = []

        # Working from the start, add positions to the to_trim until there is a residue
        # that has sufficiently high plddt.
        for pos in cluster:
            i = pos - 1
            plddt = plddt_array[i]
            if plddt < min_avg_plddt:
                to_trim.append(pos)
            else:
                break

        # Working from the end, add positions to the to_trim until there is a residue
        # that has sufficiently high plddt.
        for pos in cluster[::-1]:
            i = pos - 1
            plddt = plddt_array[i]
            if plddt < min_avg_plddt:
                to_trim.append(pos)
            else:
                break

        cluster = set(cluster) - set(to_trim)
        cluster = frozenset(cluster)
        out_clusters.append(cluster)

    return out_clusters


def get_avg_plddt(cluster, plddt_array):
    """
    This function returns the average plddt for all of the residues in the cluster.
    cluster is a frozenset of 1-indexed positions of the cluster. plddt_array is a numpy
    index of the plddts at all positions in the input structure.
    """

    running_sum = 0
    for i in cluster:

        # i in 1-indexed, so need to adjust to 0 indexing
        plddt = plddt_array[i - 1]
        running_sum += plddt

    return running_sum / len(cluster)


def get_avg_pae(cluster, pae_matrix):
    """
    Returns the average pae between every pair of residues in the cluster.
    cluster is a frozenset of 1-indexed positions of the cluster. plddt_array is a numpy
    index of the plddts at all positions in the input structure.
    """
    running_sum = 0
    n = 0
    for r1 in cluster:
        for r2 in cluster:
            pae = pae_matrix[r1, r2]
            running_sum += pae
            n += 1

    return running_sum / n


def struc_get_domains_main(args):

    talk_to_me("Reading PAE json file.")
    pae_matrix, plddt_array = parse_json_file(args.pae_path)

    talk_to_me("Parsing structure.")
    structure = pdb_to_structure_object(args.structure_file_path)

    if args.smooth_n != 0:
        talk_to_me("Smoothing out PAE matrix.")
        pae_matrix = np.apply_along_axis(smooth_array, 0, pae_matrix, n=args.smooth_n)
        pae_matrix = np.apply_along_axis(smooth_array, 1, pae_matrix, n=args.smooth_n)

    talk_to_me("Finding clusters from PAE matrix.")
    clusters = domains_from_pae_matrix_networkx(
        pae_matrix,
        pae_power=args.pae_power,
        pae_cutoff=args.pae_cutoff,
        graph_resolution=args.graph_resolution,
    )

    talk_to_me("Filtering cluster coordinates by length and pLDDT.")
    clusters = filter_clusters(
        clusters,
        plddt_array,
        min_length=args.min_domain_length,
        min_avg_plddt=args.min_domain_plddt,
    )

    talk_to_me("Trimming cluster coordinates to remove low-pLDDT-ends.")
    clusters = plddt_trim_clusters(clusters, plddt_array, args.min_domain_plddt)

    if clusters == []:
        talk_to_me(
            "No domain made it past initial filtering, so using whole structure."
        )
        structure_length = len([i for i in structure.get_residues()])
        clusters = [frozenset([i + 1 for i in range(structure_length)])]

        talk_to_me("Trimming cluster coordinates to remove low-pLDDT-ends.")
        clusters = plddt_trim_clusters(clusters, plddt_array, args.min_domain_plddt)

        talk_to_me("Filtering cluster coordinates by length and pLDDT.")
        clusters = filter_clusters(
            clusters,
            plddt_array,
            min_length=args.min_domain_length,
            min_avg_plddt=args.min_domain_plddt,
        )

    if clusters == []:
        talk_to_me(
            "No domains were found, and the full structure didn't pass filtration."
        )
        exit(0)

    talk_to_me("Writing domains to output pdb files.")
    make_output_dir(args.output_dir, is_dir=True)
    if args.plddt_report != "":
        make_output_dir(args.plddt_report, is_dir=False)
    if args.pae_report != "":
        make_output_dir(args.pae_report, is_dir=False)

    basename = os.path.basename(args.structure_file_path).rstrip(".pdb")
    i = 0
    for cluster in clusters:
        i += 1
        outfile_name = f"{basename}_domain-{i}.pdb"
        write_structure_subset(structure, cluster, f"{args.output_dir}/{outfile_name}")

        if args.plddt_report != "":
            avg_plddt = get_avg_plddt(cluster, plddt_array)
            with open(args.plddt_report, "a") as outfile:
                out = f"{outfile_name}\t{avg_plddt}\n"
                outfile.write(out)

        if args.pae_report != "":
            avg_pae = get_avg_pae(cluster, pae_matrix)
            with open(args.pae_report, "a") as outfile:
                out = f"{outfile_name}\t{avg_pae}\n"
                outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
