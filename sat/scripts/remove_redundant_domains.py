# ------------------------------------------------------------------------------------ #
# Import dependencies
# ------------------------------------------------------------------------------------ #
from .utils.misc import make_output_dir, talk_to_me
from .utils.structure import (
    pdb_to_structure_object,
    structure_to_seq,
    structure_to_pLDDT,
    write_structure_to_pdb,
)
from glob import glob
from statistics import fmean
import os


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def is_overlapping(seq1, seq2, end_overlap=15):
    """
    Determines if two strings, seq1 and seq2, either:
    1) are the same,
    2) contain one or the other, or
    3) if, after trimming the first and last 15 nucleotides of one of the sequences, it
       fits in the other sequence.

    If so, returns True.
    """

    # Check if they are equal
    if seq1 == seq2:
        return True

    # Check if one contains the other
    if seq1 in seq2 or seq2 in seq1:
        return True

    # Determine if either a sequence, when trimmed by 15 nucleotides on either side,
    # then fits into the other sequence
    if seq1[end_overlap : len(seq1) - end_overlap] in seq2:
        return True
    if seq2[end_overlap : len(seq2) - end_overlap] in seq1:
        return True

    return False


# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def remove_redundant_domains_main(args):

    # Read in the domains
    talk_to_me("Reading in structures.")
    domains = []
    for structure_path in glob(args.input_structure_glob):
        file_basename = os.path.basename(structure_path)
        domain = pdb_to_structure_object(structure_path, file_basename)
        domain.seq = structure_to_seq(domain)
        domain.pLDDT = fmean(structure_to_pLDDT(domain, "l"))
        domains.append(domain)

    # Determine if any domain overlaps and loses to another domain
    talk_to_me("Detecting overlap and filtering..")
    loosers = set()
    for i in range(len(domains)):
        domain1 = domains[i]
        for domain2 in domains[i + 1 :]:
            if is_overlapping(domain1.seq, domain2.seq):
                if len(domain1.seq) > len(domain2.seq):
                    loosers.add(domain2)
                    continue
                elif len(domain2.seq) > len(domain1.seq):
                    loosers.add(domain1)
                    continue
                if domain1.pLDDT > domain2.pLDDT:
                    loosers.add(domain1)
                    continue
                elif domain2.pLDDT > domain1.pLDDT:
                    loosers.add(domain2)
                    continue

                # Final case - if length and pLDDT are the same, need to make sure that
                # just one of them is added to the loosers.
                if (
                    len(domain1.seq) == len(domain2.seq)
                    and domain2.pLDDT == domain1.pLDDT
                ):
                    if domain1 in loosers:
                        continue
                    elif domain2 in loosers:
                        continue
                    else:
                        loosers.add(domain1)

    filtered_domains = [domain for domain in domains if domain not in loosers]

    # Write to output to the same basename
    talk_to_me("Writing output structures.")
    make_output_dir(args.output_dir, is_dir=True)
    for domain in filtered_domains:
        write_structure_to_pdb(domain, f"{args.output_dir}/{domain.id}")


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
