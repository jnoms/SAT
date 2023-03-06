from .utils.Foldseek_Dataset import Foldseek_Dataset
from .utils.misc import talk_to_me, make_output_dir
from .utils.clusters import Cluster_information
from .aln_generate_superclusters import Super_cluster


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def generate_basic_linkage_dict(alignments):
    """
    alignments is a foldseek_dataset object. This function iterates through the
    alignments and generates an output dictionary of linkages between members. Here,
    two members are considered linked if they have any one alignment.

    Output is:
    member: set(members it aligns to)

    Note that this will generate linkages for every member, query or target, rep or
    not rep, present in the alignment_file. It isn't even looking at the cluster file.
    """

    linkage_dict = dict()
    for _, alignment_group in alignments.alignment_groups.items():
        for alignment in alignment_group.alignments:
            if alignment.query not in linkage_dict:
                linkage_dict[alignment.query] = set()
            if alignment.target not in linkage_dict:
                linkage_dict[alignment.target] = set()

            linkage_dict[alignment.query].add(alignment.target)
            linkage_dict[alignment.target].add(alignment.query)
    return linkage_dict


def get_superclusters_no_reciprocal(cluster_linkages):
    """
    Takes in linked_clusters, which is a dictionary of format
    cluster: {set of clusters it is linked to}, which each cluster is represented by
    the cluster representatives.

    Returns superclusters, which is a list of sets, where each set is the members
    of the supercluster.

    Note that this function puts all reps that are linked into the same supercluster -
    it doesn't check for linkage/alignments between cluster members.
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

            merged_sc = c1_sc.union(c2_sc)
            if c1_sc == c2_sc:
                continue
            scs.remove(c1_sc)
            scs.remove(c2_sc)
            scs.append(merged_sc)
    return scs


def generate_output(ordered_scs, cluster_information):
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
        "old_rep",
    ]
    out = "\t".join(out) + "\n"
    total_members_observed = 0
    for sc in ordered_scs:
        cluster_rep = sc.largest_subcluster_rep_name
        cluster_id = sc.id
        cluster_count = sc.total_supercluster_members
        for old_rep in sc.members:
            subcluster_members = cluster_information.cluster_rep_to_members[old_rep]
            total_members_observed += len(subcluster_members)
            for member in subcluster_members:
                out += f"{cluster_rep}\t{member}\t{cluster_id}\t"
                out += f"{cluster_count}\t{old_rep}\n"

    if len(cluster_information.all_members) != total_members_observed:
        msg = "Total number of cluster members/original alignment items in the cluster "
        msg += "file doesn't match the number of members observed in superclusters! "
        msg += "Total clustered members in input cluster file: "
        msg += f"{len(cluster_information.all_members)}, total observed members: "
        msg += f"{total_members_observed}"
        raise ValueError(msg)
    return out


# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def aln_merge_clusters_main(args):

    talk_to_me("Parsing inputs.")
    data = Foldseek_Dataset()
    data.parse_alignment(args.alignment_file, args.alignment_fields)
    clusters = Cluster_information()
    clusters.parse_cluster_file(args.cluster_file, args.cluster_file_fields)

    talk_to_me("Generating linkage dict and finding superclusters.")
    linkage_dict = generate_basic_linkage_dict(data)
    superclusters = get_superclusters_no_reciprocal(linkage_dict)

    # Write to super_cluster objects
    talk_to_me("Generating supercluster objects.")
    scs = []
    seen_reps = set()
    total_clusters_obseved = 0
    for sc in superclusters:
        sc_object = Super_cluster(sc)
        sc_object.get_all_members(clusters.cluster_rep_to_members)
        scs.append(sc_object)
        for cluster in sc:
            total_clusters_obseved += 1
            seen_reps.add(cluster)

    # Generate supercluster objects for those cluster entries that didn't have
    # alignments
    for cluster in set(clusters.clusters.keys()) - seen_reps:
        total_clusters_obseved += 1
        sc_object = Super_cluster(cluster)
        sc_object.get_all_members(clusters.cluster_rep_to_members)
        scs.append(sc_object)

    # Sanity check
    if total_clusters_obseved != len(clusters.clusters):
        msg = "Have not encountered all clusters that are present in the input cluster"
        msg += f" file! Observed: {total_clusters_obseved}, in cluster file: "
        msg += f"{len(clusters.cluster_rep_to_members)}"
        raise ValueError(msg)

    # Order by total # of members
    ordered_scs = sorted(scs, key=lambda x: x.total_supercluster_members)[::-1]

    # Label each supercluster with its ID
    i = 0
    for sc in ordered_scs:
        i += 1
        sc.id = i

    talk_to_me("Writing output.")
    out = generate_output(ordered_scs, clusters)
    make_output_dir(args.output_file)
    with open(args.output_file, "w") as outfile:
        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
