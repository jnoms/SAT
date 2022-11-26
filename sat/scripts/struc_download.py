from urllib.request import urlopen, urlretrieve
import json
import os

from .utils.misc import make_output_dir, talk_to_me


def format_and_validate_args(args):
    args.infile_columns = args.infile_columns.split(",")
    return args


class AF2_downloader:
    """
    Class for downloading from alphafold2
    """

    def __init__(self, input, input_fields):
        """
        This takes in an input list of content. It also takes in
        input_fields, which is a list that describes this content. One of the content
        fields should be uniprotID which will be used for downloading from the AF2
        database. If multiple fields are input, those additional fields will be stored
        in the object and can later be appended to outfile names.
        """

        # Check that input_fields are valid
        if len(input) > 1 and input_fields == "uniprotID":
            msg = (
                "The input file has multiple columns - you need to use the"
                "--infile_columns input to indicate what they are."
            )
            raise ValueError(msg)
        if len(input) != len(input_fields):
            msg = (
                "The line of the infile doesn't have the same number of columns as "
                f"indicated in the input fields. Current line: {input}\ninput_fields: "
                f"{input_fields}"
            )
            raise ValueError(msg)
        if "uniprotID" not in input_fields:
            msg = (
                "One column must be titled uniprotID. Current infile_columns are set as"
                f" {input_fields}."
            )
            raise ValueError(msg)

        # Parse the input
        for i in range(len(input)):
            self.__dict__[input_fields[i]] = input[i]

    def __str__(self):
        return self.uniprotID

    def __repr__(self):
        return self.uniprotID

    def query_af2_api(self):
        """
        Reads in json format from the AF2 api. This will include paths to the pdb and
        PAE files.
        """
        base = "https://alphafold.ebi.ac.uk/api/prediction"
        with urlopen(f"{base}/{self.uniprotID}") as url:
            data = json.load(url)[0]
        self.data = data

    def format_outfile_path(
        self, url, output_dir, additional_field_delimiter, infile_columns
    ):
        """
        Takes the basename of the url file (with will be a .pbd or pae .json file),
        adds the output_dir path, and also appends the additional input fields.
        """
        # Remove file extension but save it for later
        file_extension = url.split(".")[-1]
        url = url.rstrip(f".{file_extension}")
        basename = os.path.basename(url)

        # Add aditional input files with the delimiter
        if len(infile_columns) > 1:
            for field in infile_columns:
                if field == "uniprotID":
                    continue
                basename += additional_field_delimiter + str(self.__dict__[field])

        # get final output path
        return f"{output_dir}/{basename}.{file_extension}"

    def download_pdb(self, output_dir, additional_field_delimiter, infile_columns):
        url = self.data["pdbUrl"]
        destination = self.format_outfile_path(
            url, output_dir, additional_field_delimiter, infile_columns
        )
        urlretrieve(url, destination)

    def download_pae(self, output_dir, additional_field_delimiter, infile_columns):
        url = self.data["paeDocUrl"]
        destination = self.format_outfile_path(
            url, output_dir, additional_field_delimiter, infile_columns
        )
        urlretrieve(url, destination)


def struc_download_main(args):
    args = format_and_validate_args(args)
    make_output_dir(args.output_dir, is_dir=True)
    talk_to_me("Starting script.")

    progress = 0
    with open(args.infile) as infile:
        for line in infile:
            line = line.rstrip("\n").split("\t")
            AF2 = AF2_downloader(line, args.infile_columns)
            AF2.query_af2_api()
            AF2.download_pdb(
                args.output_dir, args.additional_field_delimiter, args.infile_columns
            )
            AF2.download_pae(
                args.output_dir, args.additional_field_delimiter, args.infile_columns
            )
            progress += 1
            if progress % 10 == 0:
                print(progress)

    talk_to_me("Download complete.")


if __name__ == "__main__":
    msg = "This script has utilities and functions. Don't call it directly!"
    raise ValueError(msg)
