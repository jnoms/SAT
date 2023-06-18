import os
import numpy as np

from .struc_get_domains import parse_json_file, domains_from_pae_matrix_networkx
from .utils.structure import pdb_to_structure_object
from .utils.misc import talk_to_me, make_output_dir


def get_chain_ids(structure):
    """
    Input: Biopython structure object
    Output: List of strings indicating the chain IDs present in the structure
    """
    chain_ids = []
    for model in structure:
        for chain in model:
            chain_ids.append(chain.id)
    return chain_ids


def classify_members(members, boundary, chains):
    """
    Members is an iterable containing 0-indexed cluster members, usually from clustering
    the PAE matrix. Boundary is a 1-indexed interger value indicating the length of
    chain 1. This indicates the number of residues in members that are in chain 1 and
    chain 2 as an output dictionary, where keys are A or B and values are integer
    counts.

    chains is a list of chain IDs
    """
    result = {chains[0]: 0, chains[1]: 0}

    for member in members:

        # Adjust to 1-indexing
        member = member + 1

        if member <= boundary:
            result[chains[0]] += 1
        else:
            result[chains[1]] += 1
    return result


def is_cross_boundary(classification, chains):
    """
    Takes in the output from classify_members().

    Chains is a list of chain IDs
    """
    if classification[chains[0]] > 0 and classification[chains[1]] > 0:
        return True
    else:
        return False


def domain_chian_counts(domain, chain_1_residue_count):
    """
    Returns a tuple with the number of residues in chain1 and chain2.

    Notes on indexing:
    - domain is a frozen set with intergers, where the integers are 0-indexed positions.
    - SO, to do this calculation relative to the residue count need to either add
      one to each member of domain or subtract 1 from the chain_1_residue_count.

    Because the values in domain are 0-indexed, subtract 1 from chain_1_residue_count
    """
    chain1_count = len([num for num in domain if num <= chain_1_residue_count - 1])
    chain2_count = len([num for num in domain if num > chain_1_residue_count - 1])
    return chain1_count, chain2_count


def atoms_interact(
    atom1, atom2, vdw_radii={"H": 1.2, "C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8}
):
    """
    Determines if two atoms interact, where an interaction is defined as them being
    closer than the sum of their van der Waals radii plus 0.5 angstroms
    """
    d = atom1 - atom2
    atom1_vdw = vdw_radii.get(atom1.element, 0)
    atom2_vdw = vdw_radii.get(atom2.element, 0)
    if d < atom1_vdw + atom2_vdw + 0.5:
        return True
    return False


def find_interacting_residues(chain_1, chain_2):
    """
    chain_1 and chain_2 are biopython chain objects. This function runs the
    atoms_interact() function on every pair of atoms in each residue of each chain. If
    one residue of one chain has an atom that interacts with one atom of one residue in
    the other chain, those two residues are considered interacting.

    The output is a set of tuples, where each tuple is two biopython residue objects.
    The first is the reisdue from chain 1 and the second is the residue from chain 2
    that is interacting with that residue.

    Because this iterates over every pair of atoms, this can take some time (10s of
    seconds for a small protein pair, probably 1min + for larger proteins).
    """
    # a set to hold the interacting residue pairs
    interacting_residues = set()

    # for each pair of atoms, one from chain A and the other from chain B
    for atom_a in chain_1.get_atoms():
        for atom_b in chain_2.get_atoms():
            if atoms_interact(atom_a, atom_b):
                interacting_residues.add((atom_a.get_parent(), atom_b.get_parent()))

    return interacting_residues


def interacting_residues_cross_chain_pae(r1, r2, chain_1_residue_count, pae):
    """
    r1 and r2 are biopython residue objects, where r.id[1] corresponds to the
    1-indexed residue position of a residue. r1 is a residue from chain 1, and r2
    is a residue from chain 2.

    This finds the pae between those two residues.
    """
    return pae[r1.id[1] - 1, r2.id[1] - 1 + chain_1_residue_count]


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
    chains = get_chain_ids(struc)
    if len(chains) != 2:
        msg = (
            "The input structure must have two chains! Not less, not more. The input "
            f"structure has {len(chains)} chains."
        )
        raise ValueError(msg)
    chain_1 = struc[0][chains[0]]
    chain_1_residue_count = len(list(chain_1.get_residues()))
    chain_2 = struc[0][chains[1]]
    chain_2_residue_count = len(list(chain_2.get_residues()))
    pae = parse_json_file(args.pae)[0]

    # Determine if there is a cross-chain cluster
    talk_to_me("Detecting cross-chain clusters...")
    domains = domains_from_pae_matrix_networkx(pae, pae_cutoff=5)
    cross_chain_cluster = False
    for domain in domains:
        # Classification gives a dictionary with keys "A"/"B", with values being
        # the number of residues within Chain A or Chain B within the PAE cluster
        classification = classify_members(domain, chain_1_residue_count, chains)
        if is_cross_boundary(classification, chains):
            cross_chain_cluster = False
    talk_to_me(f"Cross-chain clusters status: {cross_chain_cluster}")

    talk_to_me(f"Identifying interacting residues")
    interacting_residues = find_interacting_residues(chain_1, chain_2)
    paes = []
    for r1, r2 in interacting_residues:
        paes.append(
            interacting_residues_cross_chain_pae(r1, r2, chain_1_residue_count, pae)
        )
    avg_interaction_pae = np.average(paes)
    number_of_interactions = len(interacting_residues)
    interface_residue_count = len(
        {item for sublist in interacting_residues for item in sublist}
    )

    # Write output
    talk_to_me("Writing output")
    make_output_dir(args.out_file)
    with open(args.out_file, "w") as outfile:
        try:
            member1, member2 = get_components_from_structure_path(
                args.structure, args.delimiter
            )
        except ValueError:
            msg = (
                "There may be a delimiter problem... couldn't split "
                f"{args.structure} by {args.delimiter}. Or, this input file is e.g. a"
                " dimer and you don't want to parse the component names anyway. "
                "Will simply put the structure basename as member1 and member2 in the "
                "output file."
            )
            member1 = os.path.basename(args.structure).rstrip(".pdb")
            member2 = member1
        out = [
            member1,
            member2,
            avg_interaction_pae,
            number_of_interactions,
            interface_residue_count,
            cross_chain_cluster,
            chain_1_residue_count,
            chain_2_residue_count,
        ]
        out = [str(x) for x in out]
        out = "\t".join(out) + "\n"
        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
