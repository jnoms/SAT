import os

from .utils.structure import pdb_to_structure_object, structure_to_pLDDT
from .utils.misc import make_output_dir, talk_to_me


def find_disorder(plddts, cutoff, n_sequential):
    """
    plddts is a list of residue-wise pLDDT values. This function returns a list of
    0-indexed positions of residues that are considered disordered.

    A residue is considered disordered if it is in a stretch of at least n_sequential
    residues that have a pLDDT of <= cutoff.
    """

    stretches = []
    current_stretch = []

    for i, value in enumerate(plddts):
        if value <= cutoff:
            current_stretch.append(i)
        else:
            if len(current_stretch) >= n_sequential:
                stretches.append(current_stretch)
            current_stretch = []

    # Check for a valid stretch at the end of the list
    if len(current_stretch) >= n_sequential:
        stretches.append(current_stretch)

    residues = [item for sublist in stretches for item in sublist]

    return residues


def find_order(plddts, cutoff, n_sequential):
    """
    plddts is a list of residue-wise pLDDT values. This function returns a list of
    0-indexed positions of residues that are considered ordered.

    A residue is considered ordered if it is in a stretch of at least n_sequential
    residues that have a pLDDT of >= cutoff.
    """
    stretches = []
    current_stretch = []

    for i, value in enumerate(plddts):
        if value >= cutoff:
            current_stretch.append(i)
        else:
            if len(current_stretch) >= n_sequential:
                stretches.append(current_stretch)
            current_stretch = []

    # Check for a valid stretch at the end of the list
    if len(current_stretch) >= n_sequential:
        stretches.append(current_stretch)

    residues = [item for sublist in stretches for item in sublist]

    return residues


def struc_disorder_main(args):
    # Parse the structure
    talk_to_me("Parsing structure file.")
    struc = pdb_to_structure_object(args.structure_file)
    plddts = structure_to_pLDDT(struc, format="l")

    talk_to_me("Counting order and disorder.")
    ordered_positions = find_order(
        plddts, cutoff=args.order_cutoff, n_sequential=args.n_sequential
    )
    order = len(ordered_positions)
    disordered_positions = find_disorder(
        plddts, cutoff=args.order_cutoff, n_sequential=args.n_sequential
    )
    disorder = len(disordered_positions)
    intermediate = len(plddts) - order - disorder
    total = len(plddts)

    file_basename = os.path.basename(args.structure_file)

    out = [file_basename, order, disorder, intermediate, total]
    out = [str(x) for x in out]
    out = "\t".join(out) + "\n"

    make_output_dir(args.out_file)
    with open(args.out_file, "w") as outfile:
        outfile.write(out)
