# ------------------------------------------------------------------------------------ #
# Import dependencies
# ------------------------------------------------------------------------------------ #

import json
import numpy as np
import networkx as nx
from networkx.algorithms import community
from Bio.PDB import *

from .utils.misc import make_output_dir, talk_to_me
from .utils.structure import pdb_to_structure_object

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
            msg = f"Something went wrong - number of residues parsed for plddt"
            msg += f" is not the same length as the cluster. Cluser in question is"
            msg += f" c_start: {c_start}, c_end: {c_end}"
            raise ValueErorr(msg)

        filtered_clusters.append(cluster)

    return filtered_clusters


def write_structure_subset(structure, residues_to_keep, outfile):
    """
    Writes a pdb outfile that is just a subset of the input structure from the
    residue start to the residue end. Note that the pdb file is 1-indexed, so the
    coordinates are expected to be 1-indexed as well.

    - structure: biopython structure object
    - residues_to_keep: some iterable/list of numbers, where each number is the position
        of one of the residues to keep.
    """

    class ResSelect(Select):
        def accept_residue(self, res):
            if res.id[1] in residues_to_keep:
                return True
            else:
                return False

    io = PDBIO()
    io.set_structure(structure)
    io.save(outfile, ResSelect())


def get_domains_main(args):

    # Find clusters from json file.
    talk_to_me("Reading PAE json file.")
    pae_matrix, plddt_array = parse_json_file(args.pae_path)

    talk_to_me("Finding clusters from PAE matrix.")
    clusters = domains_from_pae_matrix_networkx(
        pae_matrix,
        pae_power=args.pae_power,
        pae_cutoff=args.pae_cutoff,
        graph_resolution=args.graph_resolution,
    )

    # Filter clusters based on length and average plddt
    talk_to_me("Filtering cluster coordinates by length and pLDDT.")
    clusters = filter_clusters(
        clusters,
        plddt_array,
        min_length=args.min_domain_length,
        min_avg_plddt=args.min_domain_plddt,
    )

    # Parse structure
    talk_to_me("Parsing structure.")
    structure = pdb_to_structure_object(args.structure_file_path)

    # Generate output directory if necessary
    make_output_dir(args.output_prefix, is_dir=False)

    # Write to output directory
    talk_to_me("Writing domains to output pdb files.")
    i = 0
    for cluster in clusters:
        i += 1
        write_structure_subset(
            structure, cluster, f"{args.output_prefix}_domain-{i}.pdb"
        )


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
