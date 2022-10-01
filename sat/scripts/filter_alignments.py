from .utils.Foldseek_Dataset import Foldseek_Dataset
from .utils.misc import talk_to_me, make_output_dir


def validate_and_format_args(args):
    # Format args
    if args.alignment_fields != "":
        args.alignment_fields = args.alignment_fields.split(",")
    return args


def filter_alignments_main(args):
    args = validate_and_format_args(args)

    talk_to_me("Reading in alignments.")
    data = Foldseek_Dataset()
    data.parse_alignment(args.alignment_file, args.alignment_fields)

    if args.alignment_fields == "":
        output_fields = data.input_alignment_fields
        out = "\t".join(data.input_alignment_fields) + "\n"
    else:
        output_fields = args.alignment_fields
        out = "\t".join(args.alignment_fields) + "\n"

    talk_to_me("Filtering and generating output string.")
    for _, aln_group in data.alignment_groups.items():
        aln_group.filter_alignments(
            args.filter_field, args.max_val_filter_field, args.min_val_filter_field
        )
        aln_group.keep_top_N_alignments(args.filter_field, args.N)
        for aln in aln_group.alignments:
            out += aln.write_output(output_fields)

    talk_to_me("Writing output file.")
    make_output_dir(args.output_file)
    with open(args.output_file, "w") as outfile:
        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
