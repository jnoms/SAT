import re
from collections import Counter

from .utils.Foldseek_Dataset import Foldseek_Dataset
from .utils.misc import talk_to_me, make_output_dir
from .utils.clusters import Cluster_information


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def parse_ecod_descriptions(ecod_path):
    """
    Parses the ECOD description file into a nested dictionary of format:
    ecod_dict[ecod_accession][level] = value
    Where level is A, X, H, T, or F and the value is the classification (such as, e.g.
    aspartice protease)
    """

    ecod_dict = dict()
    with open(ecod_path) as infile:
        for line in infile:

            line = line.rstrip("\n").strip(">").split(" | ")

            ecod_accession = line[0].split(".")[0]
            ecod_classifications = line[3]

            # This regex is necessary because some field values have commas >:(
            # the [3:] removes the initial preceeding "A: "
            ecod_classifications = re.split(", [AXHTF]: ", ecod_classifications[3:])
            A, X, H, T, F = [i.replace(" ", "_") for i in ecod_classifications]

            ecod_dict[ecod_accession] = dict()
            ecod_dict[ecod_accession]["A"] = A
            ecod_dict[ecod_accession]["X"] = X
            ecod_dict[ecod_accession]["H"] = H
            ecod_dict[ecod_accession]["T"] = T
            ecod_dict[ecod_accession]["F"] = F

    return ecod_dict


def parse_domain_list(domain_list_path):
    """
    This parses the domain list file into a simple dictionary of structure
    domain: cluster_ID
    """

    domain_list_dict = dict()
    with open(domain_list_path) as infile:
        for line in infile:
            line = line.rstrip("\n").split("\t")
            query = line[0].rstrip(".pdb")
            cluster_ID = line[1]
            domain_list_dict[query] = cluster_ID

    return domain_list_dict


def alignment_to_ecod_count(data, cluster_info, ecod_dict):
    """
    This function iterates through each of the queries int he alignments and finds
    the target (which is, of course, an ECOD entry). It uses the ecod_dict to look up
    the ECOD entry, providing a dictionary of structure
    ECOD_info[level] = value.

    The ultimate output of this function is a nested dictionary/counter of structure
    cluster_ID --> level ---> COUNTER of level: count

    Inputs:
    - data: Foldseek_dataset object
    - cluster_info: Cluster_informatoin object
    - ecod_dict: dictionary of structure ecod_dict[ecod_accession][level] = value. This
        basically lets you connect the ECOD accession to the values at each level
        (level being A, X, H, T, F - the ECOD level classifications).

    Output:
    Tuple of (cluster_info_dict, aligned_domains)
    - cluster_info dict is the main output, described above
    - aligned domains is a set of the domains that were output - this is just for
        error checking later. Also, lets you add back information for queries that
        don't have any alignments.
    """

    cluster_info_dict = dict()
    aligned_domains = set()

    for _, alignment_group in data.alignment_groups.items():
        for alignment in alignment_group.alignments:

            # Because ECOD queries/targets dont end in .pdb, need to add them
            query = f"{alignment.query}.pdb"

            try:
                cluster_ID = cluster_info.cluster_member_to_ID[query]
            except KeyError:
                # This happens when a domain isn't classified into a cluster_ID because
                # it had no alignments during all-by-all alignment - not even to itself.
                # This can happen for very short domains.
                print(alignment.query)
                continue

            # keep a record of which domains were successfully aligned
            aligned_domains.add(query)

            ECOD_info = ecod_dict[alignment.target]

            if cluster_ID not in cluster_info_dict:
                cluster_info_dict[cluster_ID] = dict()

            for level in ["A", "X", "H", "T", "F"]:
                if level not in cluster_info_dict[cluster_ID]:
                    cluster_info_dict[cluster_ID][level] = Counter()

                cluster_info_dict[cluster_ID][level][ECOD_info[level]] += 1

    return cluster_info_dict, aligned_domains


def add_unaligned_entires_to_cluster_info(
    cluster_info, cluster_info_dict, aligned_domains
):
    """
    This function takes in the cluster_info object (which has all members),
    aligned_domains set (which just has the queries that aligned to ECOD, and are
    therefore already present in cluster_info_dict), and the cluster_info_dict, and
    adds an UNALLIGNED entry to each level for each unaligned cluster member.
    """
    cluster_info_dict = cluster_info_dict.copy()
    for domain, cluster_ID in cluster_info.cluster_member_to_ID.items():
        if domain not in aligned_domains:

            if cluster_ID not in cluster_info_dict:
                cluster_info_dict[cluster_ID] = dict()

            for level in ["A", "X", "H", "T", "F"]:
                if level not in cluster_info_dict[cluster_ID]:
                    cluster_info_dict[cluster_ID][level] = Counter()
                cluster_info_dict[cluster_ID][level]["UNALLIGNED"] += 1

    return cluster_info_dict


def cluster_info_dict_to_output(cluster_info_dict, cluster_ID_to_rep):
    out = ["cluster_ID", "cluster_rep", "level", "value", "count"]
    out = "\t".join(out) + "\n"

    for cluster_ID, ECOD_dict in cluster_info_dict.items():
        cluster_rep = cluster_ID_to_rep[cluster_ID]
        for level, level_counter in ECOD_dict.items():
            for value, count in level_counter.items():
                out_line = [str(cluster_ID), cluster_rep, level, value, str(count)]
                out_line = "\t".join(out_line) + "\n"
                out += out_line

    return out


# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def aln_ecod_purity_main(args):

    talk_to_me("Reading in ECOD description and clusters.")
    ecod_dict = parse_ecod_descriptions(args.ECOD_key)

    cluster_info = Cluster_information()
    cluster_info.parse_cluster_file(args.cluster_file, args.cluster_file_fields)
    cluster_info.add_cluster_ID()
    cluster_info.generate_member_and_rep_to_cluster_ID()
    cluster_info.cluster_ID_to_rep = {
        v: k for k, v in cluster_info.cluster_rep_to_ID.items()
    }

    talk_to_me("Reading in alignments.")
    data = Foldseek_Dataset()
    data.parse_alignment(args.alignment_file, args.alignment_fields)

    talk_to_me("Formating cluster information and ECOD counts")
    cluster_info_dict, aligned_domains = alignment_to_ecod_count(
        data, cluster_info, ecod_dict
    )
    cluster_info_dict = add_unaligned_entires_to_cluster_info(
        cluster_info, cluster_info_dict, aligned_domains
    )

    talk_to_me("Generating and writing output")
    out = cluster_info_dict_to_output(cluster_info_dict, cluster_info.cluster_ID_to_rep)
    make_output_dir(args.output_file)
    with open(args.output_file, "w") as outfile:
        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
