from glob import glob

from sat.scripts.struc_download import struc_download_main


# ------------------------------------------------------------------------------------ #
# Whole-script tests
# ------------------------------------------------------------------------------------ #
def teststruc_download_main(tmp_path):

    # Define inputs
    class args:
        pass

    args.infile = "tests/test_data/structure_related/struc_download_input.txt"
    args.output_dir = f"{tmp_path}/output"

    args.infile_columns = "uniprotID,taxonID"
    args.additional_field_delimiter = "__"

    # Run the script
    struc_download_main(args)

    # Validate the outputs - make sure there are 3 pdb and 3 pae (json) files
    assert len(glob(f"{args.output_dir}/*pdb")) == 3
    assert len(glob(f"{args.output_dir}/*json")) == 3
