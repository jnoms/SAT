import os

from .utils.structure import pdb_to_structure_object, structure_to_pLDDT
from .utils.misc import make_output_dir


def struc_to_plddt_main(args):
    structure = pdb_to_structure_object(args.structure_file)
    plddts = structure_to_pLDDT(structure)

    count = 0
    plddt_sum = 0
    for _, plddt in plddts.items():
        count += 1
        plddt_sum += plddt
    plddt = round(plddt_sum / count, ndigits=2)

    # Save or print
    if args.out_file == "":
        print(plddt)
    else:
        make_output_dir(args.out_file)
        base = os.path.basename(args.structure_file)
        out = f"{base}\t{plddt}\n"
        with open(args.out_file, "a") as outfile:
            outfile.write(out)
