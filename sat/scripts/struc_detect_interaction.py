import os

from .struc_get_domains import parse_json_file, domains_from_pae_matrix_networkx
from .utils.structure import pdb_to_structure_object, find_near_residues
from .utils.misc import talk_to_me, make_output_dir


def classify_members(members, boundary):
    """
    Members is an iterable containing 0-indexed cluster members, usually from clustering
    the PAE matrix. Boundary is a 1-indexed interger value indicating the length of
    chain A. This indicates the number of residues in members that are in chain A and
    chain B as an output dictionary, where keys are A or B and values are integer
    counts.
    """
    result = {"A": 0, "B": 0}

    for member in members:

        # Adjust to 1-indexing
        member = member + 1

        if member <= boundary:
            result["A"] += 1
        else:
            result["B"] += 1
    return result


def is_cross_boundary(classification):
    """
    Takes in the output from classify_members()
    """
    if classification["A"] > 0 and classification["B"] > 0:
        return True
    else:
        return False


def get_components_from_structure_path(structure_path, delimiter):
    """
    Gets the name of the individual molecules being compared - it's assumed that these
    are stored in the structure path with a delimiter. E.g. structure path could be
    molecule1__molecule2.pdb, with __ being the delimiter.
    """
    base = os.path.basename(structure_path).rstrip(".pdb")
    return base.split(delimiter)


def struc_detect_interaction_main(args):
    talk_to_me("Parsing input files")
    struc = pdb_to_structure_object(args.structure)
    if len(list(struc[0].get_chains())) != 2:
        msg = (
            "The input structure must have two chains! Not less, not more. The input "
            f"structure has {len(list(struc[0].get_chains()))} chains."
        )
        raise ValueError(msg)
    chain_A = struc[0]["A"]
    chain_A_residue_count = len(list(chain_A.get_residues()))
    pae = parse_json_file(args.pae)[0]

    talk_to_me("Clustering PAE matrix")
    domains = domains_from_pae_matrix_networkx(pae, pae_cutoff=5)

    # Looking for an interaction. But also, store the cross-chian domains incase needed
    talk_to_me("Detecting interaction...")
    interaction = False
    cross_chain_domains = []
    for domain in domains:

        # Classification gives a dictionary with keys "A"/"B", with values being
        # the number of residues within Chain A or Chain B within the PAE cluster
        classification = classify_members(domain, chain_A_residue_count)
        if is_cross_boundary(classification):
            interaction = True
            cross_chain_domains.append(domain)
    talk_to_me(f"Interaction status: {interaction}")

    # Calculate the number of residues within distance cutoff
    talk_to_me(f"Calculating residues within {args.distance_cutoff} angstroms")
    near_resdiues_A = find_near_residues(struc, "A", "B", args.distance_cutoff)
    near_resdiues_B = find_near_residues(struc, "B", "A", args.distance_cutoff)
    talk_to_me(f"Chain A: {len(near_resdiues_A)} residues")
    talk_to_me(f"Chain B: {len(near_resdiues_B)} residues")

    # Write output
    talk_to_me("Writing output")
    make_output_dir(args.out_file)
    with open(args.out_file, "w") as outfile:
        member1, member2 = get_components_from_structure_path(
            args.structure, args.delimiter
        )
        out = (
            f"{member1}\t{member2}\t{interaction}\t{len(near_resdiues_A)}"
            f"\t{len(near_resdiues_B)}\n"
        )
        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
