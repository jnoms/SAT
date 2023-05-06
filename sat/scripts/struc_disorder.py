import os

from .utils.structure import pdb_to_structure_object, structure_to_pLDDT
from .utils.misc import make_output_dir, talk_to_me


def find_disorder(plddts, cutoff, n_sequential):
    """
    plddts is a list of residue-wise pLDDT values. This function returns a list lists,
    where each sublist contains the 0-indexed positions of a stretch of residues that
    is disordered.

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

    return stretches


def find_order(plddts, cutoff, n_sequential):
    """
    plddts is a list of residue-wise pLDDT values. This function returns a list lists,
    where each sublist contains the 0-indexed positions of a stretch of residues that
    is ordered.

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

    return stretches


def struc_disorder_main(args):
    # Parse the structure
    talk_to_me("Parsing structure file.")
    struc = pdb_to_structure_object(args.structure_file)
    plddts = structure_to_pLDDT(struc, format="l")

    talk_to_me("Counting order and disorder.")
    ordered_stretches = find_order(
        plddts, cutoff=args.order_cutoff, n_sequential=args.n_sequential
    )
    ordered_residues = [item for sublist in ordered_stretches for item in sublist]
    total_ordered = len(ordered_residues)

    disordered_stretches = find_disorder(
        plddts, cutoff=args.order_cutoff, n_sequential=args.n_sequential
    )
    disordered_residues = [item for sublist in disordered_stretches for item in sublist]
    total_disordered = len(disordered_residues)
    intermediate = len(plddts) - total_ordered - total_disordered
    total = len(plddts)

    # Check if there is at least one ordered stretch that is at least
    # args.check_for_domain_len in size
    there_is_a_domain = "no"
    for stretch in ordered_stretches:
        if len(stretch) >= args.check_for_domain_len:
            there_is_a_domain = "yes"
            break

    file_basename = os.path.basename(args.structure_file)

    out = [
        file_basename,
        total_ordered,
        total_disordered,
        intermediate,
        total,
        there_is_a_domain,
    ]
    out = [str(x) for x in out]
    out = "\t".join(out) + "\n"

    make_output_dir(args.out_file)
    with open(args.out_file, "w") as outfile:
        outfile.write(out)
