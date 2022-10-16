from .utils.structure import (
    pdb_to_structure_object,
    write_structure_to_pdb,
    struc_rebase,
)
from .utils.misc import make_output_dir, talk_to_me


def struc_rebase_main(args):
    talk_to_me("Reading in structure")
    structure = pdb_to_structure_object(args.structure_file)

    talk_to_me("Rebasing.")
    struc_rebase(structure)

    make_output_dir(args.out_file)
    write_structure_to_pdb(structure, args.out_file)
    talk_to_me("Script completed.")


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
