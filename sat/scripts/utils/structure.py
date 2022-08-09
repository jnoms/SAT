from Bio.PDB import PDBParser, PDBIO
from Bio.SeqUtils import seq1


def pdb_to_structure_object(pdb_file_path, structure_name="structure"):
    """
    Given the path to a pdb file and the name of the structure, returns a biopython
    structure object.
    """
    parser = PDBParser()
    structure = parser.get_structure(structure_name, pdb_file_path)
    return structure


def structure_to_seq(structure):
    """
    Takes in a biopython structure object and returns the amino acid sequence as a
    string.
    """

    chains = {
        chain.id: seq1("".join(residue.resname for residue in chain))
        for chain in structure.get_chains()
    }

    if len(chains) > 1:
        msg = (
            "This function is designed for AF2 or Colabfold-generated structures with"
            " a single chain. The input has multiple chains!"
        )
        raise ValueError(msg)

    for chain, seq in chains.items():
        return seq


def structure_to_pLDDT(structure, format="d"):
    """
    Takes in a biopython structure object and returns pLDDT.

    If format is specified as 'd', will return a dictionary of structure pos:pLDDT,
    where pos is the 1-indexed residue position.

    If format is specified as 'l', will return a list with numeric pLDDTs.
    """

    if format not in set(["l", "d"]):
        msg = f"The 'format' parameter must be 'd' or 'l'. Received {format}."
        raise ValueError(msg)

    def get_bfactor_from_residue(residue):
        for atom in residue:
            return atom.get_bfactor()

    pLDDTs = dict()
    for residue in structure.get_residues():

        # 1-indexed position
        pos = residue.get_full_id()[3][1]

        # Get pLDDT, which is stored under bfactor of each atom
        pLDDT = get_bfactor_from_residue(residue)

        pLDDTs[pos] = pLDDT

    if format == "d":
        return pLDDTs
    elif format == "l":
        return list(pLDDTs.values())


def write_structure_to_pdb(structure, path):
    io = PDBIO()
    io.set_structure(structure)
    io.save(path)


if __name__ == "__main__":
    msg = "This script has utilities and functions. Don't call it directly!"
    raise ValueError(msg)
