from .utils.structure import (
    pdb_to_structure_object,
    write_structure_to_pdb,
    extract_chains,
)
from .utils.misc import make_output_dir, talk_to_me


def format_args(args):
    args.chains = set(args.chains.split(","))
    return args


def struc_extract_chains_main(args):
    args = format_args(args)

    talk_to_me("Parsing input file.")
    structure = pdb_to_structure_object(args.structure_file)

    talk_to_me("Extracting chain(s).")
    output_structure = extract_chains(structure, args.chains)

    talk_to_me("Writing output file.")
    make_output_dir(args.output_file)
    write_structure_to_pdb(output_structure, args.output_file)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
