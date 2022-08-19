from .utils.alignments import parse_alignment
from .utils.misc import talk_to_me, make_output_dir
from .utils.clusters import load_cluster_objects, add_cluster_ID


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def format_args(args):
    args.alignment_fields = args.alignment_fields.split(",")
    if args.all_nonredundant_out == "" and args.top_query_per_cluster_out == "":
        msg = "At least one of -1 --top_query_per_cluster_out or"
        msg += " -2 --all_nonredundant_out must be specified."
        raise ValueError(msg)
    return args


# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def process_clusters_main(args):
    args = format_args(args)

    talk_to_me("Reading in alignment file.")
    alignment_groups = parse_alignment(args.alignment_file, args.alignment_fields)

    # Load alignment_groups into Cluster objects. Will define the cluster object
    # by the foldseek_cluster_rep
    talk_to_me("Loading cluster objects.")
    cluster_objects = load_cluster_objects(args.cluster_file, alignment_groups)

    # Load cluster_ID (ranked clusters by # of items/queries), cluster count, and
    # top query.
    add_cluster_ID(cluster_objects)
    for cluster in cluster_objects:
        cluster.add_top_query()
        cluster.add_cluster_count()

    # For output fields, will be outputting the alignments using args.alignment_fields.
    # But also want to include some items from the cluster objects.
    cluster_fields = ["cluster_ID", "cluster_count", "top_query"]

    # Determine the top query (e.g. the query that has the highest number of alignments
    # or, if there is a tie, the highest average TMscore). Then output the alignments
    # for that query.
    if args.top_query_per_cluster_out != "":
        talk_to_me("Writing alignments for the top query per cluster")
        out = ""
        for cluster in cluster_objects:
            out += cluster.write_top_query_alignments(
                args.alignment_fields, cluster_fields
            )

        make_output_dir(args.top_query_per_cluster_out)
        with open(args.top_query_per_cluster_out, "w") as outfile:
            outfile.write(out)

    if args.all_nonredundant_out != "":
        talk_to_me("Writing out all non-redundant alignments.")
        out = ""
        for cluster in cluster_objects:
            out += cluster.write_all_nonredundant_alignments(
                args.alignment_fields, cluster_fields
            )

        make_output_dir(args.all_nonredundant_out)
        with open(args.all_nonredundant_out, "w") as outfile:
            outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
