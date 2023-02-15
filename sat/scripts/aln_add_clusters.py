from .utils.misc import talk_to_me, make_output_dir
from .utils.Foldseek_Dataset import Foldseek_Dataset


# ------------------------------------------------------------------------------------ #
# Functions and objects
# ------------------------------------------------------------------------------------ #
class Super_cluster:
    def __init__(self, starting_member):
        self.members = set()
        self.members.add(starting_member)

    def add_member(self, member):
        self.members.add(member)

    def all_are_linked(self, other, cluster_linkages):
        """
        Compares the members of this supercluster with the other supercluster. It asks
        if every member of this supercluster is linked to every member of the other
        supercluster.

        cluster_linkages is a dictionary of structure
        cluster: {set of other clusters it's linked to} where the clusters are
        'named' by the cluster representative
        """

        # Handle situation where the superclusters have the same members - aka
        # they are the same
        if self.members == other.members:
            return True

        for m1 in self.members:
            m1_linkages = cluster_linkages[m1]
            for m2 in other.members:
                if m2 not in m1_linkages:
                    return False

        for m2 in other.members:
            m2_linkages = cluster_linkages[m2]
            for m1 in self.members:
                if m1 not in m2_linkages:
                    return False

        return True

    def merge(self, other):
        merged_members = self.members.union(other.members)
        self.members = merged_members
        other.members = merged_members

    def get_all_members(self, cluster_dict):
        """
        Currently, the members slot is full of the cluster representatives.

        This function creates several new slots:
        - all_members: a set of all members of the clusters that are in this
          super cluster
        - largest_subcluster_rep_name: the (foldseek-derived) representative name
          of the cluster in this supercluster with the most members.
        - total_supercluster_members: an interger value indicating the length of
          the all_members slot.


        I will make
        an all_members slot that contains all members, representatives and non-
        representatives. all_members slot is a set.

        Keeps track of largest subcluster and saves the subcluster rep name in the
        largest_subcluster_rep_name slot.

        Finally, keeps track of total number of members in the supercluster in the
        total_supercluster_members slot.
        """
        all_members = set()
        largest_subcluster_rep_name = ""
        largest_subcluster_n = 0

        total_supercluster_members = 0

        for rep in self.members:
            subcluster_members = cluster_dict[rep]
            subcluster_size = len(subcluster_members)
            total_supercluster_members += subcluster_size

            if subcluster_size > largest_subcluster_n:
                largest_subcluster_rep_name = rep
                largest_subcluster_n = subcluster_size

            all_members = all_members.union(subcluster_members)

        self.all_members = all_members
        self.largest_subcluster_rep_name = largest_subcluster_rep_name
        self.total_supercluster_members = total_supercluster_members

        if len(self.all_members) != total_supercluster_members:
            msg = "Somehow the all_members slot does not equal the length of the "
            msg += "total_supercluster_members slot. That's weird."
            raise ValueError(msg)


def parse_cluster_file(cluster_file):
    """
    Returns two lookup dictionaries.
    1) cluster_rep: {set of members}
    2) member: cluster_rep
    """

    cluster_rep_to_members = dict()
    cluster_member_to_rep = dict()

    with open(cluster_file) as infile:
        for line in infile:
            line = line.rstrip("\n").split("\t")
            rep, member = line
            if rep not in cluster_member_to_rep:
                cluster_rep_to_members[rep] = set()
            cluster_rep_to_members[rep].add(member)
            cluster_member_to_rep[member] = rep

    return cluster_rep_to_members, cluster_member_to_rep


def get_alignment_dict(alignments):
    """
    Takes in a foldseek dataset and returns a dictionary of structure
    query:set of targets
    """
    alignment_dict = dict()
    for _, alignment_group in alignments.alignment_groups.items():
        alignment_dict[alignment_group.query] = set(
            [alignment.target for alignment in alignment_group.alignments]
        )

    return alignment_dict


def get_percentage_of_members_with_alignments(
    cluster_1_members, cluster_2_members, alignment_dict
):
    """
    The purpose of this function is to determine the percentage of members of cluster_1
    that have an alignment to any member of cluster 2.

    cluster_1_members and cluster_2_members are sets of the members of each of the
    two clusters.

    The alignment_dict is a dict of structure member:{set of other members it maps to}.
    """
    c1_count = 0
    for m in cluster_1_members:
        try:
            alignment_targets = alignment_dict[m]
        except KeyError:
            # This happens where the member didn't have any alignments except to
            # itself.
            continue
        if len(alignment_targets.intersection(cluster_2_members)) > 0:
            c1_count += 1

    c1_p = c1_count / len(cluster_1_members)
    return c1_p


def get_cluster_linkages(cluster_dict, alignment_dict, linkage_threshold):
    """
    This function identifies, for every cluster, the clusters it is linked to.
    Two clusters are considered linked if at least linkage_threshold fraction of members
    of at least one of the clusters have an alignment against targets in the other
    cluster.

    Input:
    - cluster_dict - dictionary of structure cluster_rep:{set of cluster members}
    - alignment_dict - dictionary of structure query:{set of targets} for every query
      in the alignment

    Output:
    - cluster_linkages - dictionary of structure cluster_rep:{set of cluster_reps it is
    linked to}.
    """

    cluster_linkages = dict()
    compared_clusters = set()

    # For a counter report
    total_comparisons = len(cluster_dict) * len(cluster_dict)
    current_count = 0
    previous_progress_report = 0

    for rep1, cluster_members1 in cluster_dict.items():
        for rep2, cluster_members2 in cluster_dict.items():

            current_count += 1

            # Don't repeat effort
            comparison = frozenset([rep1, rep2])
            if comparison in compared_clusters:
                continue
            else:
                compared_clusters.add(comparison)

            # Here, are determining the percentage of members of cluster 1 that have at
            # least one alignment against members of cluster2. Am using set intersection
            # which should be faster. Then do the reverse.
            c1_p = get_percentage_of_members_with_alignments(
                cluster_members1, cluster_members2, alignment_dict
            )
            c2_p = get_percentage_of_members_with_alignments(
                cluster_members2, cluster_members1, alignment_dict
            )

            # If cluster 1 or cluster 2 has the sufficient fraction of members with
            # cross-cluster alignments, they're considered linked.
            if c1_p >= linkage_threshold or c2_p >= linkage_threshold:
                if rep1 not in cluster_linkages:
                    cluster_linkages[rep1] = set()
                if rep2 not in cluster_linkages:
                    cluster_linkages[rep2] = set()

                cluster_linkages[rep1].add(rep2)
                cluster_linkages[rep2].add(rep1)

            comparison_progress = (current_count / total_comparisons) * 100
            comparison_progress = round(comparison_progress, 0)
            if (
                comparison_progress % 1 == 0
                and comparison_progress > previous_progress_report
            ):
                msg = "get_cluster_linkages - have completed: "
                msg += f"{comparison_progress}% of comparisons..."
                print(msg)
                previous_progress_report = comparison_progress

    return cluster_linkages


def get_superclusters(cluster_linkages):
    """
    Given the cluster linkages, this function determines which groups of linked clusters
    are superclusters (sc's). sc's are defined as groups of clusters whom are all linked
    to one another. Incorporation of a cluster into a sc requires that the cluster is
    linked to all other clusters currently in the sc.

    cluster_linkages - Dictionary of format cluster:{set of clusters that are linked
    to this cluster.}

    The output is a dictionary of structure cluster: supercluster object, where the
    supercluster object has a members slot that lists all clusters in the supercluster.

    Note that the output cluster_to_supercluster is pretty redundant - identical
    Super_cluster objects are present for all members in the same supercluster.
    """

    cluster_to_supercluster = dict()

    for cluster, linked_clusters in cluster_linkages.items():

        if cluster not in cluster_to_supercluster:
            cluster_to_supercluster[cluster] = Super_cluster(cluster)
        csc = cluster_to_supercluster[cluster]  # csc - cluster supercluster

        for linked_cluster in linked_clusters:
            if linked_cluster not in cluster_to_supercluster:
                cluster_to_supercluster[linked_cluster] = Super_cluster(linked_cluster)

            linked_csc = cluster_to_supercluster[linked_cluster]

            # Determine if all members of both sc's are linked. If so, merge
            if csc.all_are_linked(linked_csc, cluster_linkages):
                csc.merge(linked_csc)

    return cluster_to_supercluster


def get_unique_super_clusters(cluster_to_supercluster):
    """
    cluster_to_supercluster is a dictionary of cluster_rep to the supercluster object.
    Notably, the same supercluster object is present multiple times (once per member).
    This function returns a list of just the unique individual supercluster objects.
    """
    seen_scs = set()
    unique_scs = []

    for _, sc in cluster_to_supercluster.items():
        members = frozenset(sc.members)
        if members in seen_scs:
            continue
        else:
            seen_scs.add(members)
            unique_scs.append(sc)

    return unique_scs


def generate_supercluster_outputs(
    fs_dataset_object, cluster_member_to_rep, cluster_to_supercluster
):
    """
    This function iterates through all alignments and adds cluster information. It
    returns two very large strings, each of them containing tab-delimited output
    information ready to be written to output files. The output will have the same
    columns as the input alignment columns in addition to the columns cluster_ID,
    cluster_rep, and cluster_count (where the 'cluster' is the super_cluster determined
    here).

    Note that cluster information is provided based on each QUERY. Some alignments
    will have queries and targets in different clusters.

    Inputs:
    - fs_dsataset_object: Foldseek_dataset object containing the alignments
    - cluster_member_to_rep: dictionary to convert an item to it's (foldseek-derived)
      cluster representative.
    - cluster_to_supercluster: dictionary to convert an item's (foldseek-derived)
      cluster representative to the appropriate Super_cluster object.
    """

    nr_out = fs_dataset_object.input_alignment_fields + [
        "cluster_ID",
        "cluster_rep",
        "cluster_count",
    ]
    nr_out = "\t".join(nr_out) + "\n"
    top_out = nr_out

    # track progress
    count = 0
    total = len(fs_dataset_object.alignment_groups)

    for _, alignment_group in fs_dataset_object.alignment_groups.items():
        for alignment in alignment_group.alignments:
            subcluster_rep = cluster_member_to_rep[alignment.query]
            supercluster = cluster_to_supercluster[subcluster_rep]
            supercluster_rep = supercluster.largest_subcluster_rep_name

            # Generate output
            output_line = alignment.write_output(
                fs_dataset_object.input_alignment_fields
            )
            output_line = output_line.rstrip("\n")

            # Add the cluster information
            cluster_information = [
                supercluster.id,
                supercluster.largest_subcluster_rep_name,
                supercluster.total_supercluster_members,
            ]
            cluster_information = [str(x) for x in cluster_information]
            cluster_information = "\t".join(cluster_information)
            output_line = output_line + "\t" + cluster_information + "\n"

            # Write to nr and/or top as appropriate
            nr_out += output_line
            if alignment.query == supercluster_rep:
                top_out += output_line

        # Handle progress
        count += 1
        if count % 1000 == 0:
            msg = "generate_supercluster_outputs - Current Progress: "
            msg += f"{count}/{total}."
            print(msg)

    return nr_out, top_out


# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def aln_add_clusters_main(args):

    talk_to_me("Parsing cluster file and alignments")
    cluster_rep_to_members, cluster_member_to_rep = parse_cluster_file(
        args.cluster_file
    )
    data = Foldseek_Dataset()
    data.parse_alignment(args.alignment_file, args.alignment_fields)
    alignment_dict = get_alignment_dict(data)

    talk_to_me("Determining linked clusters and grouping to superclusters.")
    cluster_linkages = get_cluster_linkages(
        cluster_rep_to_members, alignment_dict, args.linkage_threshold
    )
    cluster_to_supercluster = get_superclusters(cluster_linkages)

    unique_scs = get_unique_super_clusters(cluster_to_supercluster)
    for sc in unique_scs:
        sc.get_all_members(cluster_rep_to_members)

    # Assign super_clusters an ID/ranking and remake cluster_to_supercluster
    talk_to_me("Assigning supercluster IDs.")
    ordered_scs = sorted(unique_scs, key=lambda x: x.total_supercluster_members)[::-1]
    i = 0
    cluster_to_supercluster = dict()
    for sc in ordered_scs:
        i += 1
        sc.id = i

        # Have added the ID above - now remake the cluster_to_supercluster dictionary
        # with the updated Super_cluster objects
        for member in sc.members:
            cluster_to_supercluster[member] = sc

    # Write output
    talk_to_me("Generating and writing output.")
    nr_out, top_out = generate_supercluster_outputs(
        data, cluster_member_to_rep, cluster_to_supercluster
    )

    make_output_dir(args.nr_out)
    make_output_dir(args.rep_out)
    with open(args.nr_out, "w") as outfile:
        outfile.write(nr_out)
    with open(args.rep_out, "w") as outfile:
        outfile.write(top_out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
