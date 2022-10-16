from .utils.Foldseek_Dataset import Foldseek_Dataset
from .utils.misc import talk_to_me, make_output_dir


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def get_merged_output_fields(aln1_fields, aln2_fields):
    """
    Takes in two lists of alignment fields. Keeps the basic order of aln1_fields, but
    adds all of the fields that are only present in one of the lists to the end.
    """
    # Copy to avoid side effects
    aln1_fields = aln1_fields.copy()
    aln2_fields = aln2_fields.copy()

    # Find which fields are present in one but not the other
    unique_fields = []
    unique_fields.extend(list(set(aln1_fields) - set(aln2_fields)))
    unique_fields.extend(list(set(aln2_fields) - set(aln1_fields)))

    # Remove those from the aln1_fields - will add them back at the end of the list
    for f in unique_fields:
        aln1_fields.remove(f)

    return aln1_fields + unique_fields


# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def aln_merge_main(args):
    talk_to_me(f"Reading in datasets.")
    data1 = Foldseek_Dataset()
    data1.parse_alignment(args.aln1, args.aln1_fields)
    talk_to_me(f"Dataset1 has {data1.count_alignments()} alignments.")

    data2 = Foldseek_Dataset()
    data2.parse_alignment(args.aln2, args.aln2_fields)
    talk_to_me(f"Dataset2 has {data2.count_alignments()} alignments.")

    talk_to_me(f"Merging.")
    data1.merge(data2)
    talk_to_me(f"Merged dataset has {data1.count_alignments()} alignments.")

    talk_to_me(f"Writing output.")
    output_fields = get_merged_output_fields(
        data1.input_alignment_fields, data2.input_alignment_fields
    )
    make_output_dir(args.output)
    out = "\t".join(output_fields) + "\n"
    out += data1.write_out_alignments(output_fields)
    with open(args.output, "w") as outfile:
        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
