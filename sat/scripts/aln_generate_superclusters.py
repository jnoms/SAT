from .utils.misc import talk_to_me, make_output_dir
from .utils.Foldseek_Dataset import Foldseek_Dataset
from .utils.clusters import Cluster_information


# ------------------------------------------------------------------------------------ #
# Functions and objects
# ------------------------------------------------------------------------------------ #
class Super_cluster:
    def __init__(self, starting_members):

        if isinstance(starting_members, str):
            self.members = set([starting_members])
        elif isinstance(starting_members, set):
            self.members = starting_members
        else:
            msg = "starting_members must be a string or a set."
            raise ValueError(msg)

    def add_member(self, member):
        self.members.add(member)

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

            # Initialize items in cluster_linkages if necessary
            if rep1 not in cluster_linkages:
                cluster_linkages[rep1] = set()
            if rep2 not in cluster_linkages:
                cluster_linkages[rep2] = set()

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


def all_are_linked(s1, s2, cluster_linkages):
    """
    Given two superclusters (s1 and s2), it checks if all members of
    s1 are linked to mebers of s2, and if all members of s2 are linked to members
    of s1. If all are linked, this returns True. Otherwise, it returns False.
    """
    # Handle situation where the superclusters have the same members - aka
    # they are the same
    if s1 == s2:
        return True

    for m1 in s1:
        m1_linkages = cluster_linkages[m1]
        for m2 in s2:
            if m2 not in m1_linkages:
                return False

    for m2 in s2:
        m2_linkages = cluster_linkages[m2]
        for m1 in s1:
            if m1 not in m2_linkages:
                return False

    return True


def get_superclusters(cluster_linkages):
    """
    Takes in linked_clusters, which is a dictionary of format
    cluster: {set of clusters it is linked to}.

    Returns superclusters, which is a list of sets, where each set is the members
    of the supercluster.
    """

    scs = []
    progress = 0
    total = len(cluster_linkages)
    for c1, linked_clusters in cluster_linkages.items():

        progress += 1
        if progress % 100 == 0:
            print(f"get_superclusters - progress: {progress}/{total}")

        if linked_clusters == set():
            c1_sc = set([c1])
            for sc in scs:
                if c1 in sc:
                    c1_sc = sc
            if c1_sc == set([c1]):
                scs.append(c1_sc)
                continue

        for c2 in linked_clusters:
            c1_sc = set([c1])
            c2_sc = set([c2])
            for sc in scs:
                if c1 in sc:
                    c1_sc = sc
                if c2 in sc:
                    c2_sc = sc
                if c1_sc != set([c1]) and c2_sc != set([c2]):
                    break

            if c1_sc == set([c1]) and c1_sc not in scs:
                scs.append(c1_sc)
            if c2_sc == set([c2]) and c2_sc not in scs:
                scs.append(c2_sc)

            if all_are_linked(c1_sc, c2_sc, cluster_linkages):
                merged_sc = c1_sc.union(c2_sc)
                if c1_sc == c2_sc:
                    continue
                scs.remove(c1_sc)
                scs.remove(c2_sc)
                scs.append(merged_sc)
    return scs


def generate_output(ordered_scs, cluster_rep_to_members, all_members):
    """
    Takes in ordered_scs, which is a list of supercluster objects. Then iterates through
    each subcluster rep (e.g. sc.member), converts them to their members via
    cluster_rep_to_members, and returns an output string with the columns
    cluster_rep, cluster_member, cluster_ID, subcluster_rep.

    all_members is a set of all members that should be observed - this is a sanity check
    """

    out = [
        "cluster_rep",
        "cluster_member",
        "cluster_ID",
        "cluster_count",
        "subcluster_rep",
    ]
    out = "\t".join(out) + "\n"
    total_members_observed = 0
    for sc in ordered_scs:
        cluster_rep = sc.largest_subcluster_rep_name
        cluster_id = sc.id
        cluster_count = sc.total_supercluster_members
        for subcluster_rep in sc.members:
            subcluster_members = cluster_rep_to_members[subcluster_rep]
            total_members_observed += len(subcluster_members)
            for member in subcluster_members:
                out += f"{cluster_rep}\t{member}\t{cluster_id}\t"
                out += f"{cluster_count}\t{subcluster_rep}\n"

    if len(all_members) != total_members_observed:
        msg = "Total number of cluster members/original alignment items in the cluster "
        msg += "file doesn't match the number of members observed in superclusters! "
        msg += "Total clustered members in input cluster file: "
        msg += f"{len(all_members)}, total observed members: "
        msg += f"{total_members_observed}"
        raise ValueError(msg)
    return out


# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def aln_generate_superclusters_main(args):

    talk_to_me("Parsing cluster file")
    cluster_info = Cluster_information()
    cluster_info.parse_cluster_file(args.cluster_file, args.cluster_file_fields)

    talk_to_me("Parsing alignment")
    data = Foldseek_Dataset()
    data.parse_alignment(args.alignment_file, args.alignment_fields)
    alignment_dict = get_alignment_dict(data)

    talk_to_me("Determining linked clusters and grouping to superclusters.")
    cluster_linkages = get_cluster_linkages(
        cluster_info.cluster_rep_to_members, alignment_dict, args.linkage_threshold
    )
    superclusters = get_superclusters(cluster_linkages)

    # Write to super_cluster objects and order them by total # members
    scs = []
    total_clusters_obseved = 0
    for sc in superclusters:
        total_clusters_obseved += len(sc)
        sc_object = Super_cluster(sc)
        sc_object.get_all_members(cluster_info.cluster_rep_to_members)
        scs.append(sc_object)
    ordered_scs = sorted(scs, key=lambda x: x.total_supercluster_members)[::-1]

    # Label each supercluster with its ID
    i = 0
    for sc in ordered_scs:
        i += 1
        sc.id = i

    # Sanity check
    if total_clusters_obseved != len(cluster_info.cluster_rep_to_members):
        msg = "Total number of clusters doesn't match the number of clusters that are "
        msg += "present in a supercluster!"
        raise ValueError(msg)

    out = generate_output(
        ordered_scs, cluster_info.cluster_rep_to_members, cluster_info.all_members
    )

    make_output_dir(args.outfile)
    with open(args.outfile, "w") as outfile:
        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
