from .utils.alignments import parse_alignment
from .utils.misc import make_output_dir, talk_to_me


def validate_and_format_args(args):
    taxID_location_options = set([0, 1, 2])
    if args.query_taxid_location not in taxID_location_options:
        msg = "args.query_taxID_location must be set to 0, 1, or 2. You entered "
        msg += f"{args.query_taxID_location}"
        raise ValueError(msg)
    if args.target_taxid_location not in taxID_location_options:
        msg = "args.target_taxID_location must be set to 0, 1, or 2. You entered "
        msg += f"{args.target_taxID_location}"
        raise ValueError(msg)

    if args.query_taxid_location == 2:
        msg = (
            "You indicated that the query_taxID is located in the taxID field. This is "
            "weird, because that is probably target_taxID :/"
        )
        raise ValueError(msg)
    if args.target_taxid_location == 2 and "taxid" not in args.alignment_fields:
        msg = (
            "You indicated that the target_taxid_location is in the taxid field of the "
            "alignment. However, taxid was not specified as an alignment field."
        )
        raise ValueError(msg)

    # Format args
    args.alignment_fields = args.alignment_fields.split(",")
    args.taxonomy_levels = args.taxonomy_levels.split(",")

    return args


def add_taxonomy_to_alignments_main(args):
    args = validate_and_format_args(args)

    # Read in alignment file
    talk_to_me("Reading in alignment file.")
    alignment_groups = parse_alignment(args.alignment_file, args.alignment_fields)

    # Prepare the output file header
    taxonomy_cols = []
    if args.query_taxid_location != 0:
        taxonomy_cols += [f"query_{level}" for level in args.taxonomy_levels]
    if args.target_taxid_location != 0:
        taxonomy_cols += [f"target_{level}" for level in args.taxonomy_levels]
    out = "\t".join(args.alignment_fields + taxonomy_cols) + "\n"

    # Label with taxonomy and generate an output string
    # Have a lineage lookup dict to vastly reduce the search speed by holding relevant
    # taxonIDs in memory.
    talk_to_me("Labeling taxonomies.")
    lineage_lookup = dict()
    progress = 0
    total = len(alignment_groups)
    for query, alignment_group in alignment_groups.items():
        for alignment in alignment_group.alignments:
            output_fields = args.alignment_fields

            # Deal with query lineage if specified
            if args.query_taxid_location != 0:
                alignment.add_taxid("query", args.query_taxid_location)
                if alignment.query_taxid in lineage_lookup:
                    alignment.query_lineage = lineage_lookup[alignment.query_taxid]
                else:
                    alignment.add_query_lineage(args.taxonomy_levels)
                    lineage_lookup[alignment.query_taxid] = alignment.query_lineage
                output_fields = output_fields + ["query_lineage"]

            # Deal with target lineage if specified
            if args.target_taxid_location != 0:
                alignment.add_taxid("target", args.target_taxid_location)
                if alignment.target_taxid in lineage_lookup:
                    alignment.target_lineage = lineage_lookup[alignment.target_taxid]
                else:
                    alignment.add_target_lineage(args.taxonomy_levels)
                    lineage_lookup[alignment.target_taxid] = alignment.target_lineage
                output_fields = output_fields + ["target_lineage"]

            # Write output
            out += alignment.write_output(output_fields)

        progress += 1
        if progress % 100 == 0:
            talk_to_me(f"Progress: {progress}/{total}")

    # Write output file
    talk_to_me("Writing output file.")
    make_output_dir(args.output_file)
    with open(args.output_file, "w") as outfile:
        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
