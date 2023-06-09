import os

from .struc_get_domains import parse_json_file, domains_from_pae_matrix_networkx
from .utils.structure import pdb_to_structure_object, find_near_residues
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

    talk_to_me("Clustering PAE matrix")
    domains = domains_from_pae_matrix_networkx(pae, pae_cutoff=5)

    # Looking for an interaction. But also, store the cross-chian domains incase needed
    talk_to_me("Detecting interaction...")
    interaction = False
    cross_chain_domains = []

    # Keep track of total number of residues in each chain that are present in domains
    # that are cross-chain.
    chain1_cross_chain_count = 0
    chain2_cross_chain_count = 0

    for domain in domains:

        # Classification gives a dictionary with keys "A"/"B", with values being
        # the number of residues within Chain A or Chain B within the PAE cluster
        classification = classify_members(domain, chain_1_residue_count, chains)
        if is_cross_boundary(classification, chains):
            interaction = True
            cross_chain_domains.append(domain)
            domain_chain1_count, domain_chain2_count = domain_chian_counts(
                domain, chain_1_residue_count
            )
            chain1_cross_chain_count += domain_chain1_count
            chain2_cross_chain_count += domain_chain2_count
    talk_to_me(f"Interaction status: {interaction}")
    talk_to_me(f"Chain 1: {chain1_cross_chain_count} residues in cross-chain clusters.")
    talk_to_me(f"Chain 2: {chain2_cross_chain_count} residues in cross-chain clusters.")

    # Calculate the number of residues within distance cutoff
    talk_to_me(f"Calculating residues within {args.distance_cutoff} angstroms")
    near_residues_1 = find_near_residues(
        struc, chains[0], chains[1], args.distance_cutoff
    )
    near_residues_2 = find_near_residues(
        struc, chains[1], chains[0], args.distance_cutoff
    )
    talk_to_me(
        f"Chain 1: {len(near_residues_1)} residues within {args.distance_cutoff} "
        "angstroms of chain 2"
    )
    talk_to_me(
        f"Chain 2: {len(near_residues_2)} residues within {args.distance_cutoff} "
        "angstroms of chain 1"
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
            interaction,
            chain_1_residue_count,
            chain_2_residue_count,
            chain1_cross_chain_count,
            chain2_cross_chain_count,
            chain1_cross_chain_count / chain_1_residue_count,
            chain2_cross_chain_count / chain_2_residue_count,
            len(near_residues_1),
            len(near_residues_2),
        ]
        out = [str(x) for x in out]
        out = "\t".join(out) + "\n"
        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
