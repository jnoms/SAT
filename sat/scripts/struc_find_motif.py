import re
import itertools

from .utils.structure import pdb_to_structure_object, struc_to_seq
from .utils.misc import read_fasta_to_memory


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def generate_motif_options(motif, AAs="ARNDBCEQZGHILKMFPSTWYV"):
    in_options = False
    pos = 0
    pos_options = dict()
    for c in motif:
        if c == "[":
            in_options = True
            pos += 1
            continue

        if c == "]":
            in_options = False
            continue

        if in_options:
            c = c.upper()
            if pos not in pos_options:
                pos_options[pos] = set()
            pos_options[pos].add(c)

        if not in_options:
            pos += 1
            if pos not in pos_options:
                pos_options[pos] = set()
            if c.lower() != "x":
                c = c.upper()
                pos_options[pos].add(c)
            else:
                for AA in AAs:
                    pos_options[pos].add(AA)

    return pos_options


def generate_all_possible_motifs(pos_options):
    if len(pos_options) > 10:
        msg = "Max of 10 spaces in a motif! Sorry"
        raise ValueError(msg)

    if len(pos_options) == 1:
        options = set(itertools.product(pos_options[1]))
    elif len(pos_options) == 2:
        options = set(itertools.product(pos_options[1], pos_options[2]))
    elif len(pos_options) == 3:
        options = set(itertools.product(pos_options[1], pos_options[2], pos_options[3]))
    elif len(pos_options) == 4:
        options = set(
            itertools.product(
                pos_options[1], pos_options[2], pos_options[3], pos_options[4]
            )
        )
    elif len(pos_options) == 5:
        options = set(
            itertools.product(
                pos_options[1],
                pos_options[2],
                pos_options[3],
                pos_options[4],
                pos_options[5],
            )
        )
    elif len(pos_options) == 6:
        options = set(
            itertools.product(
                pos_options[1],
                pos_options[2],
                pos_options[3],
                pos_options[4],
                pos_options[5],
                pos_options[6],
            )
        )
    elif len(pos_options) == 7:
        options = set(
            itertools.product(
                pos_options[1],
                pos_options[2],
                pos_options[3],
                pos_options[4],
                pos_options[5],
                pos_options[6],
                pos_options[7],
            )
        )
    elif len(pos_options) == 8:
        options = set(
            itertools.product(
                pos_options[1],
                pos_options[2],
                pos_options[3],
                pos_options[4],
                pos_options[5],
                pos_options[6],
                pos_options[7],
                pos_options[8],
            )
        )
    elif len(pos_options) == 9:
        options = set(
            itertools.product(
                pos_options[1],
                pos_options[2],
                pos_options[3],
                pos_options[4],
                pos_options[5],
                pos_options[6],
                pos_options[7],
                pos_options[8],
                pos_options[9],
            )
        )
    elif len(pos_options) == 10:
        options = set(
            itertools.product(
                pos_options[1],
                pos_options[2],
                pos_options[3],
                pos_options[4],
                pos_options[5],
                pos_options[6],
                pos_options[7],
                pos_options[8],
                pos_options[9],
                pos_options[10],
            )
        )

    options = set("".join(i) for i in options)
    return options


def get_match(option, seq):
    out = ""
    if re.search(option, seq) is not None:
        for result in re.finditer(option, seq):
            out_result = [result[0], result.start() + 1, result.end()]
            out_result = [str(i) for i in out_result]
            out_result = "\t".join(out_result) + "\n"
            out += out_result
        return out
    else:
        return ""


# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def struc_find_motif_main(args):
    if args.structure == "" and args.seq == "" and args.fasta == "":
        msg = "You must enter either a structure, a fasta, or a seq..."
        msg += " Call the -h flag for help and instructions."
        raise ValueError(msg)

    # Parse seq
    if args.fasta != "":
        fasta = read_fasta_to_memory(args.fasta)
        seq = ""
        for _, fasta_seq in fasta.items():
            if seq != "":
                msg = "There is more than one sequence in the fasta file - stop that."
                raise ValueError(msg)
            else:
                seq = fasta_seq

    if args.structure != "":
        structure = pdb_to_structure_object(args.structure)
        seq = struc_to_seq(structure)

    if args.seq != "":
        seq = args.seq

    seq = seq.upper()

    pos_options = generate_motif_options(args.motif)
    options = generate_all_possible_motifs(pos_options)

    out = ["match", "start", "end"]
    out = "\t".join(out) + "\n"
    for option in options:
        if option in seq:
            out += get_match(option, seq)
    print(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
