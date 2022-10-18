from .utils.misc import talk_to_me, make_output_dir
from .utils.Foldseek_Dataset import Foldseek_Dataset


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def format_args(args):
    if args.all_nonredundant_out == "" and args.top_query_per_cluster_out == "":
        msg = "At least one of -1 --top_query_per_cluster_out or"
        msg += " -2 --all_nonredundant_out must be specified."
        raise ValueError(msg)
    return args


# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def aln_add_clusters_main(args):
    args = format_args(args)

    talk_to_me("Reading in alignment file.")
    data = Foldseek_Dataset()
    data.parse_alignment(args.alignment_file, args.alignment_fields)

    # Load alignment_groups into Cluster objects. Will define the cluster object
    # by the foldseek_cluster_rep
    talk_to_me("Loading cluster objects.")
    data.load_cluster_objects(args.cluster_file)

    # Load cluster_ID (ranked clusters by # of items/queries)
    data.add_cluster_ID()

    output_fields = data.input_alignment_fields
    cluster_fields = ["cluster_ID", "cluster_count", "top_query"]

    # Generate outputs
    talk_to_me("Generating and writing outputs.")
    out_header = "\t".join(data.input_alignment_fields + cluster_fields) + "\n"
    out_dict = data.write_out_cluster_alignments(
        output_fields, cluster_fields=cluster_fields, top_or_nonredundant="both"
    )

    # Write outputs
    if args.top_query_per_cluster_out != "":
        make_output_dir(args.top_query_per_cluster_out)
        with open(args.top_query_per_cluster_out, "w") as outfile:
            out = out_header + out_dict["top"]
            outfile.write(out)

    if args.all_nonredundant_out != "":
        make_output_dir(args.all_nonredundant_out)
        with open(args.all_nonredundant_out, "w") as outfile:
            out = out_header + out_dict["nr"]
            outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
