import os

from .utils.structure import pdb_to_structure_object, structure_to_pLDDT


def struc_qc_main(args):
    structure = pdb_to_structure_object(args.structure_file_path)
    plddts = structure_to_pLDDT(structure, "l")

    total = 0
    passing = 0
    for r in plddts:
        total += 1
        if r >= args.plddt_cutoff:
            passing += 1

    structure_base = os.path.basename(args.structure_file_path)

    out = [structure_base, total, passing, round(passing / total, 2)]
    out = [str(i) for i in out]
    out = "\t".join(out)
    print(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
