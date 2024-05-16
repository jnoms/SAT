from .utils.misc import talk_to_me, make_output_dir


def parse_structure_key(structure_key_file, delim=",,", existing_dict=dict()):
    """
    Takes in a path to a file with the first column the structure, and second column
    the 4-digit identifier, and returns a dictionary of format
    identifier:structure.

    Add an existing_dict if you just want to update that dictionary - this lets you
    have multiple structure_key files.
    """

    key_to_structure = existing_dict.copy()
    with open(structure_key_file) as infile:
        for line in infile:
            line = line.rstrip("\n")
            structure, key = line.split(delim)

            if len(key) != 4:
                msg = "The key is expected to be a four-digit identifier! This key, "
                msg += f" {key}, is a length of {len(key)}!"
                raise ValueError(msg)

            if key in key_to_structure:
                msg = f"Have obseved a key, {key}, that is already present in"
                msg += " key_to_structure! This means it may be present multiple times!"
                raise ValueError(msg)

            key_to_structure[key] = structure

    return key_to_structure


def aln_parse_dali_main(args):

    # Parse structure key. Multiple keys are allowed if they are a comma-delimited
    # string.
    talk_to_me("Parsing structure key.")
    if "," not in args.structure_key:
        structure_key = parse_structure_key(
            args.structure_key, delim=args.structure_key_delim
        )
    else:
        structure_key = dict()
        for structure_key_file in args.structure_key.split(","):
            structure_key = parse_structure_key(
                structure_key_file,
                delim=args.structure_key_delim,
                existing_dict=structure_key,
            )

    talk_to_me("Parsing the alignment file.")
    query_id = ""
    query = ""
    qlen = 0
    out = [
        "query",
        "target",
        "query_id",
        "target_id",
        "alnlen",
        "qlen",
        "tlen",
        "cov",
        "pident",
        "rmsd",
        "z",
    ]
    out = "\t".join(out) + "\n"

    with open(args.alignment_file) as infile:

        for line in infile:
            line = line.rstrip("\n")

            if line == "":
                continue

            # Need to stop looping after the summary segment
            if line.startswith("# Structural equivalences"):
                break

            # Find query
            if line.startswith("# Query: "):
                query_id = line.replace("# Query: ", "")[:-1]

                # If the query name was a commandline paramater, go with that.
                # Otherwise, look in structure key
                if args.query_name != "":
                    talk_to_me(f"query_name has been specified as {args.query_name}")
                    query = args.query_name
                else:
                    try:
                        query = structure_key[query_id]
                    except KeyError:
                        msg = f"Cannot find the query_id, {query_id}, in structure_key!"
                        msg += " Continuing anyway."
                        talk_to_me(msg)
                continue

            # Find the colnames and the other header row
            if line.startswith("#"):
                continue

            # Otherwise, we're at alignments
            line = line.split()
            try:
                _, chain, z, rmsd, alnlen, tlen, pident, *_ = line
            except ValueError:
                msg = (
                    f"Failed to parse LINE: {line}. This probably happened because the"
                    "RMSD is above 100 so its column merged with z score column."
                )
                print(msg)

            # Removing the chain identifier (-A for chain A, etc)
            target_id = chain[:-2]
            try:
                target = structure_key[target_id]
            except KeyError:
                msg = f"Cannot find the target_id, {target_id}, in structure_key!"
                msg += " Continuing anyway."
                target = ""
                talk_to_me(msg)

            z = float(z)
            rmsd = float(rmsd)
            alnlen = float(alnlen)
            tlen = float(tlen)
            pident = float(pident)

            # Update qlen if we get a self alignment
            if query_id == target_id:
                qlen = tlen

            # Generate a "cov" column that is alnlen/(max(qlen, tlen)). If there is no
            # labeled qlen, it will just end up alnlen/tlen.
            cov = alnlen / max(qlen, tlen)

            # Do filtration
            if (
                cov < args.min_cov
                or z < args.min_z
                or alnlen < args.min_alnlen
                or rmsd < args.min_rsmd
            ):
                continue

            # Write output
            out_line = [
                query,
                target,
                query_id,
                target_id,
                alnlen,
                qlen,
                tlen,
                cov,
                pident,
                rmsd,
                z,
            ]
            out_line = [str(i) for i in out_line]
            out_line = "\t".join(out_line) + "\n"
            out += out_line

    make_output_dir(args.output_file)
    talk_to_me("Writing output file.")
    with open(args.output_file, "w") as outfile:
        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
