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

    # For output fields, will be outputting the alignments using args.alignment_fields.
    # But also want to include some items from the cluster objects.
    cluster_fields = ["cluster_ID", "cluster_count", "top_query"]

    # Load cluster_ID (ranked clusters by # of items/queries)
    add_cluster_ID(cluster_objects)

    # Make any output directories that are necessary, and open output files
    if args.top_query_per_cluster_out != "":
        make_output_dir(args.top_query_per_cluster_out)
        top_query_per_cluster_out_outfile = open(args.top_query_per_cluster_out, "w")

    if args.all_nonredundant_out != "":
        make_output_dir(args.all_nonredundant_out)
        all_nonredundant_out_outfile = open(args.all_nonredundant_out, "w")

    # Keeping track...
    total = len(cluster_objects)
    progress = 0

    talk_to_me("Processing clusters and writing to output file(s)")
    for cluster in cluster_objects:

        # Adding the top_query and the cluser count
        cluster.add_top_query()
        cluster.add_cluster_count()

        # Writing out the alignments of the top query for each cluster - e.g. the query
        # with the most alignments or, if a tie, the query with the highest average
        # TMscore.
        if args.top_query_per_cluster_out != "":
            out = cluster.write_top_query_alignments(
                args.alignment_fields, cluster_fields
            )
            top_query_per_cluster_out_outfile.write(out)

        # Writing out all non-redundant alignments that make it into a cluster object.
        # Note that alignments with a query and target not in the same cluster are
        # excluded.
        if args.all_nonredundant_out != "":
            out = cluster.write_all_nonredundant_alignments(
                args.alignment_fields, cluster_fields
            )
            all_nonredundant_out_outfile.write(out)

        # Update progress
        progress += 1
        if progress % 100 == 0:
            talk_to_me(f"Progress: {progress}/{total}")

    # Close the output files
    if args.top_query_per_cluster_out != "":
        top_query_per_cluster_out_outfile.close()
    if args.all_nonredundant_out != "":
        all_nonredundant_out_outfile.close()


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
