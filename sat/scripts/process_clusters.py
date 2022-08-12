from .utils.alignments import parse_alignment
from .utils.misc import make_output_dir, talk_to_me

from collections import OrderedDict, Counter


# ------------------------------------------------------------------------------------ #
# Classes
# ------------------------------------------------------------------------------------ #
class Alignment_cluster:
    def __init__(self, alignment_rep):
        self.alignment_rep = alignment_rep
        self.alignments = []

    def remove_redundant_alignments(self):
        """
        We only need one alignment for every query-target pair. However, because
        clusters were generated with an all-by-all search, each query-target pair is
        present as query-target as well as target-query.

        This function keeps only one of the alignments.
        """
        alignments = self.alignments

        seen = []
        filtered_alignments = []
        for alignment in alignments:

            pair1 = [alignment.query, alignment.target]
            pair2 = [alignment.target, alignment.query]

            if pair1 not in seen and pair2 not in seen:
                seen.append(pair1)
                seen.append(pair2)
                filtered_alignments.append(alignment)
            else:
                continue

        self.alignments = filtered_alignments


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #


def parse_cluster_files(cluster_file_path: str):
    """
    Given a path to the cluster file, reads in clusters to the format of
    cluster_rep:set(cluster_members). Note that the cluster_rep will be in the
    cluster_members.
    """
    clusters = dict()
    with open(cluster_file_path) as infile:
        for line in infile:
            line = line.rstrip("\n").split("\t")
            cluster_rep, cluster_member = line

            if cluster_rep not in clusters:
                clusters[cluster_rep] = [cluster_member]
            else:
                clusters[cluster_rep].append(cluster_member)

    return clusters


def generate_cluster_lookup_dict(clusters: dict):
    """
    clusters -  a dict of structure cluster_rep:[list of cluster members]

    This function outputs a dict of the following format:
    member:(cluster_ID, cluster_rep, #_members_in_cluster). cluster_ID is the 1-indexed
    order of clusters, where cluster 1 has the most members, cluster 2 has the next
    most, etc.
    """

    # First, count and make sure that every item is only in one cluster
    members = Counter()
    for cluster_rep, cluster_members in clusters.items():
        for member in cluster_members:
            members[member] += 1
    if members.most_common()[0][1] > 1:
        msg = "There seems to be an entry that is a member of multiple clusters! This "
        msg += "isn't allowed. Something is wrong."
        raise ValueError(msg)

    # Convert clusters to an OrderedDict, ordered by the number of members of each
    # cluster
    clusters = OrderedDict(
        sorted(clusters.items(), key=lambda x: len(x[1]), reverse=True)
    )

    cluster_lookup = dict()
    for i, (cluster_rep, cluster_members) in enumerate(clusters.items()):
        cluster_ID = i + 1
        number_in_cluster = len(cluster_members)

        for cluster_member in cluster_members:
            cluster_lookup[cluster_member] = (
                cluster_ID,
                cluster_rep,
                number_in_cluster,
            )

    return cluster_lookup


def load_alignment_clusters(alignment_groups, cluster_lookup):

    alignment_clusters = dict()

    for query, alignment_group in alignment_groups.items():

        for alignment in alignment_group.alignments:

            # Determine if query or target is not in a cluster - in this case, skip
            if (
                alignment.query not in cluster_lookup
                or alignment.target not in cluster_lookup
            ):
                continue

            # Determine if query and target are in the same cluster. If not, skip
            query_cluster_info = cluster_lookup[alignment.query]
            target_cluster_info = cluster_lookup[alignment.target]
            if query_cluster_info != target_cluster_info:
                continue

            # Add the cluster information into the alignment object
            cluster_ID, cluster_rep, cluster_count = query_cluster_info
            alignment.cluster_ID = cluster_ID
            alignment.cluster_rep = cluster_rep
            alignment.cluster_count = cluster_count

            # Add to alignment_cluster object
            if alignment.cluster_rep not in alignment_clusters:
                alignment_clusters[alignment.cluster_rep] = Alignment_cluster(
                    alignment.cluster_rep
                )
                alignment_clusters[
                    alignment.cluster_rep
                ].cluster_ID = alignment.cluster_ID
                alignment_clusters[
                    alignment.cluster_rep
                ].cluster_rep = alignment.cluster_rep
                alignment_clusters[
                    alignment.cluster_rep
                ].cluster_count = alignment.cluster_count

            alignment_clusters[alignment.cluster_rep].alignments.append(alignment)

    return alignment_clusters


# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def process_clusters_main(args):

    # Parse alignment and cluster files
    talk_to_me("Parsing input files.")
    alignment_groups = parse_alignment(args.alignment_file, args.alignment_fields)
    clusters = parse_cluster_files(args.cluster_file)
    cluster_lookup = generate_cluster_lookup_dict(clusters)

    # Separate alignments into alignment cluster objects
    talk_to_me("Separating alignments into clusters.")
    alignment_clusters = load_alignment_clusters(alignment_groups, cluster_lookup)

    # Final processing and save to outfile
    talk_to_me("Final processing and writing output file.")

    progress = 0
    total = len(alignment_clusters)

    out = ""
    for rep, alignment_cluster in alignment_clusters.items():
        alignment_cluster.remove_redundant_alignments()
        for alignment in alignment_cluster.alignments:
            output_fields = alignment.alignment_fields + [
                "cluster_ID",
                "cluster_rep",
                "cluster_count",
            ]
            out += alignment.write_output(output_fields)
        progress += 1
        if progress % 100 == 0:
            talk_to_me(f"Progress: {progress}/{total}")

    make_output_dir(args.output_file)
    with open(args.output_file, "w") as outfile:
        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
