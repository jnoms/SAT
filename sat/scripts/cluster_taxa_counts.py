from .utils.Foldseek_Dataset import Foldseek_Dataset
from .utils.misc import make_output_dir, talk_to_me


def format_args(args):
    args.taxonomy_levels = args.taxonomy_levels.split(",")
    return args


def cluster_taxa_counts_main(args):
    args = format_args(args)

    talk_to_me("Reading in alignments")
    data = Foldseek_Dataset()
    data.parse_alignment(args.alignment_file)

    talk_to_me("Counting taxa")
    out = data.write_out_cluster_taxa_count(args.taxonomy_levels)

    talk_to_me("Writing output")
    make_output_dir(args.outfile)
    with open(args.outfile, "w") as outfile:
        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
