from .utils.misc import make_output_dir, talk_to_me
from .utils.Foldseek_Dataset import Foldseek_Dataset


def validate_and_format_args(args):
    # Format args
    args.alignment_fields = args.alignment_fields.split(",")
    args.taxonomy_levels = args.taxonomy_levels.split(",")

    return args


def add_taxonomy_to_alignments_main(args):
    args = validate_and_format_args(args)

    # Read in alignment file
    talk_to_me("Reading in alignment file.")
    data = Foldseek_Dataset()
    data.parse_alignment(args.alignment_file, args.alignment_fields)

    # Prepare the output file header
    taxonomy_cols = []
    taxonomy_cols += [f"query_{level}" for level in args.taxonomy_levels]
    taxonomy_cols += [f"target_{level}" for level in args.taxonomy_levels]
    taxonomy_cols += ["query_taxonID", "target_taxonID"]
    out = "\t".join(args.alignment_fields + taxonomy_cols) + "\n"

    # Label with taxonomy and generate an output string
    data.add_taxon_to_alignments(
        args.taxonomy_levels, args.taxonID_finder_delimiter, args.taxonID_finder_pos
    )

    # Note that query_taxon and target_taxon will yield a Taxon() object in the
    # alignment object, which will then be parsed for the canonical lineage.
    # The query_taxonID and target_taxonID refer to the flat taxonIDs.
    out += data.write_out_alignments(
        args.alignment_fields
        + ["query_taxon", "target_taxon", "query_taxonID", "target_taxonID"]
    )

    # Write output file
    talk_to_me("Writing output file.")
    make_output_dir(args.output_file)
    with open(args.output_file, "w") as outfile:
        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
