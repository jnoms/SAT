from .utils.structure import pdb_to_structure_object, structure_to_seq
from .utils.misc import make_output_dir


def structure_to_seq_main(args):
    structure = pdb_to_structure_object(args.structure_file, args.header)
    seq = structure_to_seq(structure)

    # Save or print
    if args.out_file == "":
        print(seq)
    else:
        if args.header == "":
            msg = "If --out_file is specified, you must fill in --header."
            raise ValueError(msg)
        make_output_dir(args.out_file)
        with open(args.out_file, "a") as outfile:
            out = f">{args.header}\n{seq}\n"
            outfile.write(out)
