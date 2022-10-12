# ------------------------------------------------------------------------------------ #
# Import dependencies
# ------------------------------------------------------------------------------------ #
from .utils.Foldseek_Dataset import Foldseek_Dataset
from .utils.uniprot import uniprot_object
from .utils.misc import talk_to_me, make_output_dir


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def process_args(args):
    if args.alignment_fields != "":
        args.alignment_fields = args.alignment_fields.split(",")

    return args


def parse_uniprot_information(uniprot_file_path):
    """
    Takes in the path to the uniprot file. The file should be tab-delimited and have
    the columns accession, geneName, and fullName.

    This will return a dictionary of structure
    accession: uniprot_object,
    where the uniprot_objects have the slots accession, geneName, and fullName
    """

    uniprot_objects = dict()
    with open(uniprot_file_path) as infile:

        for line in infile:
            line = line.rstrip("\n").strip("\t").split("\t")
            if len(line) == 3:
                accession, geneName, fullName = line
            elif len(line) == 1:
                accession = line[0]
            else:
                msg = (
                    "The uniprot information for line f{line} is shorter than"
                    "expected! Something is werid."
                )
                raise ValueError(msg)
            uo = uniprot_object(accession)
            uo.geneName = geneName
            uo.fullName = fullName
            uniprot_objects[accession] = uo

    return uniprot_objects


def get_uniprot_accession_from_alignment(alignment, field="target"):
    """
    Given an alignment object, will search in the indicated field for the uniprotID
    and will return it.
    """
    ID = alignment.__dict__[field]
    if "-" in ID:
        ID = ID.split("-")[1]

    return ID


# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def aln_add_uniprot_main(args):
    args = process_args(args)

    talk_to_me("Reading in alignments.")
    data = Foldseek_Dataset()
    data.parse_alignment(args.alignment_file, args.alignment_fields)

    talk_to_me("Reading in uniprot information.")
    uniprot_dict = parse_uniprot_information(args.uniprot_information)

    # For each alignment, add the uniprot information and generate output
    talk_to_me("Adding uniprot information to each alignment.")
    if args.alignment_fields != "":
        output_fields = args.alignment_fields + [
            "uniprot_accession",
            "uniprot_geneName",
            "uniprot_fullName",
        ]
    else:
        output_fields = data.input_alignment_fields + [
            "uniprot_accession",
            "uniprot_geneName",
            "uniprot_fullName",
        ]
    out = "\t".join(output_fields) + "\n"
    for query, alignment_group in data.alignment_groups.items():

        for alignment in alignment_group.alignments:

            accession = get_uniprot_accession_from_alignment(alignment)
            uo = uniprot_dict[accession]
            alignment.uniprot_accession = uo.accession
            alignment.uniprot_geneName = uo.geneName
            alignment.uniprot_fullName = uo.fullName
            out += alignment.write_output(output_fields)

    talk_to_me("Writing output.")
    make_output_dir(args.outfile)
    with open(args.outfile, "w") as outfile:
        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
