from .utils.structure import (
    pdb_to_structure_object,
    write_structure_to_pdb,
    rebase_structure,
)
from .utils.misc import make_output_dir


def rebase_structure_main(args):
    structure = pdb_to_structure_object(args.structure_file)
    rebase_structure(structure)
    make_output_dir(args.out_file)
    write_structure_to_pdb(structure, args.out_file)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
