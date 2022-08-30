import pathlib
import os
import sys
import gzip
from Bio import SeqIO
import argparse


# ------------------------------------------------------------------------------------ #
# Misc
# ------------------------------------------------------------------------------------ #
def talk_to_me(msg):
    """
    Prints message to the screen, but prefixes it with the parent pyscript name.
    """
    base = os.path.basename(sys.argv[0])
    print(f"{base}: {msg}")


def arg_str2bool(v):
    """
    For use as an argparse argument type. Makes it easy to use boolean flags.
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


# ------------------------------------------------------------------------------------ #
# Read/write
# ------------------------------------------------------------------------------------ #
def make_output_dir(path, is_dir=False):
    """
    Finds the dirname of the path, and makes the directory if it doesn't exist.
    Notably, if the path is the a directory it skips the dirname finding and just makes
    that directory
    """
    if not is_dir:
        out_dir = os.path.dirname(path)
        pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
    else:
        pathlib.Path(path).mkdir(parents=True, exist_ok=True)


def read_fasta_to_memory(input_fasta):
    """
    Reads fasta into a memory as a dictionary with header:sequence.
    This function can handle .gzip files, but input_fasta needs to
    end with .gz
    """
    fasta_dict = dict()

    if not input_fasta.endswith(".gz"):
        for seq_record in SeqIO.parse(input_fasta, "fasta"):
            fasta_dict[seq_record.id] = str(seq_record.seq)

    elif input_fasta.endswith(".gz"):
        with gzip.open(input_fasta, "rt") as handle:
            for seq_record in SeqIO.parse(handle, "fasta"):
                fasta_dict[seq_record.id] = str(seq_record.seq)

    return fasta_dict


if __name__ == "__main__":
    msg = "This script has utilities and functions. Don't call it directly!"
    raise ValueError(msg)
