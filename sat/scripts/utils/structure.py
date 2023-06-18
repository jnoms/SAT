from Bio.PDB import PDBParser, PDBIO, Select, Selection, NeighborSearch
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.SeqUtils import seq1


def pdb_to_structure_object(pdb_file_path, structure_name="structure"):
    """
    Given the path to a pdb file and the name of the structure, returns a biopython
    structure object.
    """
    parser = PDBParser()
    structure = parser.get_structure(structure_name, pdb_file_path)
    return structure


def struc_to_seq(structure):
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


def write_structure_subset(structure, residues_to_keep, outfile):
    """
    Writes a pdb outfile that is just a subset of the input structure from the
    residue start to the residue end. Note that the pdb file is 1-indexed, so the
    coordinates are expected to be 1-indexed as well.

    - structure: biopython structure object
    - residues_to_keep: some iterable/list of numbers, where each number is the position
        of one of the residues to keep.
    """

    class ResSelect(Select):
        def accept_residue(self, res):
            if res.id[1] in residues_to_keep:
                return True
            else:
                return False

    io = PDBIO()
    io.set_structure(structure)
    io.save(outfile, ResSelect())


def struc_rebase(structure):
    """
    Input is a biopython structure object. This function renumbers all residues such
    that the first residue starts at 1 and all residues are sequential - so, it takes
    out number gaps.

    Note that this function acts in-place. It doesn't return anything - it just edits
    the structure that is input.
    """
    for i, residue in enumerate(structure.get_residues()):
        i += 1
        res_id = list(residue.id)
        res_id[1] = i
        residue.id = tuple(res_id)


def compare_structures(structure1, structure2):
    """
    This makes sure that all residues are the same in the two input structures. The
    inputs should be biopython structure objects. This is used for testing.
    """
    s1_residues = [r for r in structure1.get_residues()]
    s2_residues = [r for r in structure2.get_residues()]

    assert len(s1_residues) == len(s2_residues)
    for i in range(len(s1_residues)):
        r1 = s1_residues[i]
        r2 = s2_residues[i]
        assert r1 == r2


def extract_chains(structure, chains):
    """
    This function makes a new structure object containing the chain or chains that
    are specified.

    structure - a structure object
    chains - an iterable of chain IDs you want to keep in the new structure object
    """

    # Create a new structure and model
    new_structure = Structure("new_structure")
    new_model = Model(0)
    new_structure.add(new_model)

    # Get the chain from the original structure
    for model in structure:
        for chain in model:
            if chain.get_id() in chains:
                # Create a new chain and copy residues from the original chain
                new_chain = Chain(chain.get_id())
                for residue in chain:
                    new_chain.add(residue.copy())
                # Add the new chain to the new model
                new_model.add(new_chain)
                break

    # Check that all desired chains are present in the final structure
    final_chains = [chain.get_id() for chain in new_structure.get_chains()]
    for chain in chains:
        if chain not in final_chains:
            msg = (
                f"Cannot find the chain, {chain}. The input 'chains' iterable was "
                f"{chains}. Make sure this chain is in the input structure."
            )
            raise ValueError(msg)

    return new_structure


if __name__ == "__main__":
    msg = "This script has utilities and functions. Don't call it directly!"
    raise ValueError(msg)
