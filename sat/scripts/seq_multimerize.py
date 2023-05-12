# ------------------------------------------------------------------------------------ #
# Import dependencies
# ------------------------------------------------------------------------------------ #
from .utils.misc import read_fasta_to_memory, make_output_dir, talk_to_me


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
# Parse infiles to seqeunces
def parse_infiles(infiles):
    """
    Infiles is a list of paths to fastas
    """
    seqs = []
    for infile in infiles:
        fasta = read_fasta_to_memory(infile)
        if len(fasta) != 1:
            msg = (
                "There should be only one sequence in the fasta! This script doesn't "
                f"know how to handle more than one. Detected {len(fasta)} sequences."
            )
            raise ValueError(msg)
        header = list(fasta.keys())[0]
        seq = fasta[header]
        seqs.append(seq)

    return seqs


def multimerize_seqs(sequences, cardinality):
    """
    sequences is a list of sequences, cardinality is a list of integers indicating the
    number of times a sequence should occur.
    """
    if len(sequences) != len(cardinality):
        msg = "sequences is a dif size than cardinality!"
        raise ValueError(msg)

    out_seq = ""
    for i, seq in enumerate(sequences):
        for n in range(cardinality[i]):
            out_seq += seq + ":"

    return out_seq.rstrip(":")


# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def seq_multimerize_main(args):

    # Parse the input fasta(s) and cardinality
    talk_to_me("Parsing input files")
    infiles = args.infiles.split(",,")
    sequences = parse_infiles(infiles)

    cardinality = args.cardinality.split(",,")
    try:
        cardinality = [int(n) for n in cardinality]
    except ValueError:
        msg = (
            "Did you use a single-comma delimiter for cardinality? That's a cause for"
            " this error! You need to use double-comma delimiter"
        )
        raise ValueError(msg)

    if len(infiles) != len(cardinality):
        msg = (
            f"There are {len(infiles)} infiles, but only cardinality information for "
            f"{len(cardinality)}. The input for infiles is {args.infiles} and the input"
            f" for cardinality is {cardinality}"
        )
        raise ValueError(msg)

    talk_to_me("Generating output sequences and file")
    out_seq = multimerize_seqs(sequences, cardinality)

    out = f">{args.header}\n{out_seq}"
    make_output_dir(args.out_file)
    with open(args.out_file, "w") as outfile:
        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
