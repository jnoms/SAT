from Bio import SeqIO

from .utils.misc import talk_to_me, make_output_dir


def pad_to_five_digits(number):
    """
    Takes an integer and returns it as a string of at least 5 characters,
    padded with leading zeros.

    e.g.
    in: 5, out: 00005
    in: 1432, out: 01432
    in: 24132, out: 24132
    """
    if number > 99999:
        msg = "You broke everything! Nice nice. Never thought this would happen."
        raise ValueError(msg)

    # Convert the number to a string and pad it with zeros to make it at least 5 characters long
    padded_string = str(number).zfill(5)
    return padded_string


def remove_bad_characters(
    input_string,
    bad_characters=[
        ":",
        ";",
        "_",
        "\\",
        "/",
        "?",
        "{",
        "}",
        "[",
        "]",
        "*",
        "!",
        "#",
        "+",
        " ",
        "!",
    ],
):
    for char in bad_characters:
        input_string = input_string.replace(char, "-")
    return input_string


def parse_gb_record(record, delimiter, accs=set()):
    """
    record is a biopython gb record object.mro

    acc_list is an iterable containing nucleotide accessions you want to keep. If
    specified, only those will be written out.

    Returns out_fasta, out_csv:
    - out_fasta: fasta containing all protein information
    - out_csv: csv information about each protein
    """

    genome_acc = record.name
    organism_name = remove_bad_characters(record.name)

    protein_order_tracker = 1

    out_fasta = ""
    out_csv = ""

    if accs != set():
        if genome_acc not in accs:
            return out_fasta, out_csv

    # Loop over each feature in the record
    for feature in record.features:

        # Only keep CDS'
        if feature.type != "CDS":
            continue

        # Get feature details

        # A dictionary with information such as gene name, protein ID, etc.
        qualifiers = feature.qualifiers

        location = feature.location
        protein_id = qualifiers.get("protein_id", ["X"])[0]
        locus_tag = qualifiers.get("locus_tag", ["X"])[0]
        protein_name = qualifiers.get("product", ["X"])[0]
        protein_order = pad_to_five_digits(protein_order_tracker)
        start = int(location.start)
        end = int(location.end)
        strand = location.strand
        protein_sequence = qualifiers.get("translation", ["None"])[0]

        # Remove bad characters
        locus_tag = remove_bad_characters(locus_tag)
        protein_name = remove_bad_characters(protein_name)

        # Output name!
        output_name = f"{genome_acc}{delimiter}{protein_id}{delimiter}{locus_tag}{delimiter}{protein_order}"

        # Write output fasta
        entry = f">{output_name}\n{protein_sequence}\n"
        out_fasta += entry

        # Write output table
        entry = f"{output_name},{genome_acc},{locus_tag},{protein_id},{start},{end},{strand},{protein_order},{organism_name},{protein_name}\n"
        out_csv += entry

        protein_order_tracker += 1

    return out_fasta, out_csv


def seq_parse_genbank_main(args):
    out_csv = "name,scaffold,locus_tag,protein_id,start,end,strand,gene_order,organism_name,protein_name\n"
    out_fasta = ""

    # Parse a file including only nucleotide accessions that are desired (optional)
    accs = set()
    if args.nuc_accs_to_keep != "":
        with open(args.nuc_accs_to_keep) as infile:
            for line in infile:
                accs.add(line.rstrip("\n"))

    talk_to_me("Parsing genbank file...")
    # This try/except stuff is needed to parse malformed genbank entries.
    try:
        with open(args.gb, "r") as file:
            while True:
                try:
                    # Parse the next record
                    record = next(SeqIO.parse(file, "genbank"))
                    out_fasta_entry, out_csv_entry = parse_gb_record(
                        record, args.delimiter, accs
                    )
                    out_csv += out_csv_entry
                    out_fasta += out_fasta_entry
                except StopIteration:
                    # If no more records, stop the loop
                    break
                except Exception as e:
                    msg = (
                        "An entry has an error. This record is just before the "
                        f"one that has an error: {record.id}. Thus, the errored entry "
                        "and the entry directly following it won't be output."
                    )
                    print(msg)
                    # Skip to the next record by continuing the while loop
                    continue
    except Exception as e:
        print(f"Failed to read file: {str(e)}")

    talk_to_me("Writing output files")
    make_output_dir(args.out_csv)
    make_output_dir(args.out_fasta)
    with open(args.out_csv, "w") as outfile:
        outfile.write(out_csv)
    with open(args.out_fasta, "w") as outfile:
        outfile.write(out_fasta)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
