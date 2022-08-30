# ------------------------------------------------------------------------------------ #
# Import dependencies
# ------------------------------------------------------------------------------------ #
from .utils.misc import read_fasta_to_memory, make_output_dir, talk_to_me


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def chunk_seq(seq: str, max_length: int, output_overlapping: bool):
    """
    Given an input sequence, splits it into chunks of size max_length. If
    output_overlapping is set to true, the chunks will overlap by round(max_length/2).
    """
    if output_overlapping:
        out = [
            seq[i : i + max_length]
            for i in range(0, len(seq), max_length - round(max_length / 2))
        ]
    else:
        out = [seq[i : i + max_length] for i in range(0, len(seq), max_length)]

    # If the last chunk is just a part of the preceeding chunk, remove it.
    if out[-1] in out[-2]:
        del out[-1]

    return out


# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def chunk_fasta_main(args):
    fasta_dict = read_fasta_to_memory(args.in_fasta)
    out_fasta_dict = dict()

    for header, seq in fasta_dict.items():

        # If the sequence is below the max sequence length but longer than the
        # minimum seq length, write it out without chunking it
        if len(seq) <= args.max_seq_length:
            if len(seq) < args.minimum_sequence_output_size:
                continue
            out_fasta_dict[header] = seq
            continue

        # Split the sequence into chunks
        chunks = chunk_seq(seq, args.max_seq_length, args.overlapping_chunks)

        for i, chunk in enumerate(chunks):
            # Ditch chunks that fall below a minimum cutoff size - colabfold
            # fails with very small inputs.
            if len(chunk) < args.minimum_sequence_output_size:
                continue

            new_header = "PART{}_{}".format(str(i), header)
            out_fasta_dict[new_header] = chunk

    # Write output
    if not args.individual:
        talk_to_me("Writing output to a single fasta file.")
        make_output_dir(args.out_fasta)
        with open(args.out_fasta, "w") as outfile:
            for header, seq in out_fasta_dict.items():
                out = ">{}\n{}\n".format(header, seq)
                outfile.write(out)
    else:
        talk_to_me("Writing output to separate fasta files.")
        if not args.out_fasta.endswith("/"):
            args.out_fasta += "/"
        make_output_dir(args.out_fasta, is_dir=True)
        for header, seq in out_fasta_dict.items():
            outfile = args.out_fasta + header + ".fasta"
            with open(outfile, "w") as outfile:
                out = ">{}\n{}\n".format(header, seq)
                outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
