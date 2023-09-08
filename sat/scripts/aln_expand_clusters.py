import os


from .utils.clusters import Cluster_information
from .utils.misc import talk_to_me, make_output_dir


def merge_clusters_and_subclusters_to_string(clusters, subclusters):
    """
    clusters and subclusters are each Cluster_info() objects. It is assumed that all
    members of clusters are the reps of the subclusters. This script generates a tab-
    delimited string of structure cluster_rep, subcluster_rep, cluster_member.

    Because each member of each cluster in the clusters object is a subcluster rep,
    this is basically adding in all of the subcluster members.

    This string should be output and later parsed into a new Cluster_info() object...
    this is useful because during parsing all of the information will be distribubed
    correctly into Cluster() and Cluster_member() objects.
    """

    intermediate_out = ["cluster_rep", "subcluster_rep", "cluster_member"]
    intermediate_out = "\t".join(intermediate_out) + "\n"

    for cluster_rep, cluster_members in clusters.cluster_rep_to_members.items():
        cluster_rep = cluster_rep.rstrip(".pdb")
        for cluster_member in cluster_members:
            cluster_member = cluster_member.rstrip(".pdb")

            if cluster_member not in subclusters.cluster_rep_to_members:
                msg = (
                    "For this script, it is expected that each member of the clusters "
                    "should be a REPRESENTATIVE of the subclusters. Currently, cannot "
                    f"find the cluster member {cluster_member} in the subclusters "
                    "cluster_rep_to_members dictionary."
                )
                raise ValueError(msg)
            subcluster_members = subclusters.cluster_rep_to_members[cluster_member]
            for subcluster_member in subcluster_members:
                out_line = [cluster_rep, cluster_member, subcluster_member]
                out_line = "\t".join(out_line) + "\n"
                intermediate_out += out_line

    return intermediate_out


def aln_expand_clusters_main(args):
    talk_to_me("Parsing cluster and subcluster files")
    clusters = Cluster_information()
    clusters.parse_cluster_file(args.cluster_file, args.cluster_file_fields)

    subclusters = Cluster_information()
    subclusters.parse_cluster_file(args.subcluster_file, args.subcluster_file_fields)

    talk_to_me("Merging clusters and subclusters, and writing an intermediate file")
    intermediate_out = merge_clusters_and_subclusters_to_string(clusters, subclusters)
    intermediate_file_path = args.out_file + "_TEMP"
    make_output_dir(args.out_file)
    with open(intermediate_file_path, "w") as outfile:
        outfile.write(intermediate_out)

    talk_to_me("Parsing the intermediate file, adding cluster_IDs")
    final_clusters = Cluster_information()
    final_clusters.parse_cluster_file(intermediate_file_path, cluster_file_fields="")
    final_clusters.add_cluster_ID()

    talk_to_me("Writing output file")
    out = [
        "cluster_ID",
        "cluster_rep",
        "subcluster_rep",
        "cluster_member",
        "cluster_count",
    ]
    out = "\t".join(out) + "\n"
    for cluster_rep, cluster_obj in final_clusters.clusters.items():
        cluster_count = len(cluster_obj.cluster_members)
        for cluster_member in cluster_obj.cluster_members:
            cluster_member_name = cluster_member.cluster_member
            subcluster_rep = cluster_member.subcluster_rep
            cluster_ID = cluster_member.cluster_ID

            out_line = [
                cluster_ID,
                cluster_rep,
                subcluster_rep,
                cluster_member_name,
                cluster_count,
            ]
            out_line = [str(x) for x in out_line]
            out_line = "\t".join(out_line) + "\n"
            out += out_line
    with open(args.out_file, "w") as outfile:
        outfile.write(out)

    # clean up intermediate file
    os.remove(intermediate_file_path)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
