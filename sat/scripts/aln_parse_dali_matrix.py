import re

from .utils.misc import talk_to_me, make_output_dir


def validate_args(args):
    if args.tree == "" and args.matrix == "":
        msg = "Must specify --tree and/or --matrix."
        raise ValueError(msg)
    if args.tree != "" and args.tree_out == "":
        msg = "You specified --tree, so you must specify --tree_out!"
        raise ValueError(msg)
    if args.matrix != "" and args.matrix_out == "":
        msg = "You specified --matrix, so you must specify --matrix_out!"
        raise ValueError(msg)


def parse_key(key_path):
    key = dict()
    with open(key_path) as infile:
        for line in infile:
            line = line.rstrip("\n").split(",,")
            key[line[1]] = line[0].rstrip(".pdb")
    return key


def format_nwk(nwk, key):

    # Regex to find potential keys (assuming keys end with a letter followed by a colon and a number)
    pattern = re.compile(r"(\w+):\d+(\.\d+)?")
    keys_found = set(pattern.findall(nwk))
    keys_found = {k[0] for k in keys_found}

    # Replace keys in the nwk string if the base key is in the dictionary
    for full_key in keys_found:
        base_key = full_key[:-1]
        nwk = nwk.replace(full_key, key[base_key])

    return nwk


def parse_matrix(matrix_path):
    """
    Parses the matrix file as a list of lists
    """
    matrix = []
    first_line = True
    with open(matrix_path) as infile:
        for line in infile:
            # Ignore the first line, which is an integer
            if first_line:
                first_line = False
                continue
            row = line.split("\t")
            matrix.append(row)

    return matrix


def format_matrix(matrix, key):
    """
    matrix is a list of lists, where each list is a row. The first item in each row
    should be the ID. This script will convert the ID using the key.
    """

    # Replace the IDs
    for row in matrix:
        row[0] = key[row[0][:-1]]
    return matrix


def aln_parse_dali_matrix_main(args):
    validate_args(args)

    talk_to_me("Parsing key")
    key = parse_key(args.key)

    if args.tree != "":

        talk_to_me("Processing tree")

        with open(args.tree) as infile:
            nwk = infile.read()
            formatted_nwk = format_nwk(nwk, key)

        make_output_dir(args.tree_out)
        with open(args.tree_out, "w") as outfile:
            outfile.write(formatted_nwk)

    if args.matrix != "":

        talk_to_me("Processing matrix")
        matrix = parse_matrix(args.matrix)
        matrix = format_matrix(matrix, key)

        make_output_dir(args.matrix_out)
        with open(args.matrix_out, "w") as outfile:
            out = ""
            for row in matrix:
                row = "\t".join(row)
                out += row
            outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
