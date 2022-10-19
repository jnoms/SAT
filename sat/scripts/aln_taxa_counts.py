from .utils.Foldseek_Dataset import Foldseek_Dataset
from .utils.misc import make_output_dir, talk_to_me


def format_args(args):
    args.taxonomy_levels = args.taxonomy_levels.split(",")
    return args


def aln_taxa_counts_main(args):

    args = format_args(args)

    talk_to_me("Reading in alignments")
    data = Foldseek_Dataset()
    data.parse_alignment(args.alignment_file)

    talk_to_me("Adding taxon objects to every alignment")
    data.add_taxon_to_alignments(args.taxonomy_levels)

    talk_to_me("Grouping alignments into clusters")
    data.load_cluster_objects_from_alignments()

    talk_to_me("Writing output.")
    make_output_dir(args.output_file)
    out = (
        "\t".join(
            ["cluster_ID", "top_query", "level", "superkingdom", "taxon", "count"]
        )
        + "\n"
    )
    for cluster in data.clusters:
        cluster.add_taxa_counts(args.taxonomy_levels)
        out += cluster.write_taxa_counts_output()
    with open(args.output_file, "w") as outfile:
        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
