from sat.scripts.rebase_structure import rebase_structure_main
from sat.scripts.utils.structure import (
    pdb_to_structure_object,
    structure_to_pLDDT,
    structure_to_seq,
)


# ------------------------------------------------------------------------------------ #
# Whole-script tests
# ------------------------------------------------------------------------------------ #
def test_rebase_structure(tmp_path):

    # Define inputs
    class args:
        pass

    args.structure_file = (
        "tests/test_data/structure_related/discontinuous_structure.pdb"
    )
    args.out_file = f"{tmp_path}/rebased.pdb"

    # Run script
    rebase_structure_main(args)

    # Load expected and observed structures
    expected = pdb_to_structure_object("tests/test_data/structure_related/rebased.pdb")
    observed = pdb_to_structure_object(args.out_file)

    # Make sure the primary sequence and pLDDT is still the same as expected
    assert structure_to_seq(expected) == structure_to_seq(observed)
    assert structure_to_pLDDT(expected) == structure_to_pLDDT(observed)

    # Make the result starts at one and doens't skip any gaps
    current = 1
    for residue in observed.get_residues():
        pos = residue.get_full_id()[3][1]
        assert pos == current
        current += 1
