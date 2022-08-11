from ete3 import NCBITaxa

ncbi = NCBITaxa()

canonical_prefixes = {
    "superkingdom": "sk__",
    "kingdom": "k__",
    "phylum": "p__",
    "class": "c__",
    "order": "o__",
    "family": "f__",
    "genus": "g__",
    "species": "s__",
}


def name_to_taxID(name: str):
    """
    Given name, returns taxonID. If the name taxonID can't be found,
    will try to strip any following numbers from the name. Otherwise,
    will return X.
    """
    try:
        return ncbi.get_name_translator([name])[name][0]
    except:
        # Try to repair by stripping any numbers from the end of the name
        name = name.rstrip("0123456789").rstrip(" ")
        try:
            return ncbi.get_name_translator([name])[name][0]
        except:
            return "X"


def get_level(taxonID: int):
    """
    Given a single taxonID, returns the taxonomic level.
    """
    level = list(ncbi.get_rank([taxonID]).values())

    # Unknown taxonID would yield [], which can't be indexed by [0] to get the string
    if level == []:
        level = "UNKNOWN"
    else:
        level = level[0]
    return level


def get_name(taxonID: int):
    """
    Given a single taxonID, returns the name of the taxon.
    """
    name = list(ncbi.get_taxid_translator([taxonID]).values())

    # Unknown taxonID would yield [], which can't be indexed by [0] to get the string
    if name == []:
        name = "UNKNOWN"
    else:
        name = name[0]

    name = name.replace(" ", "_")
    return name


def get_lineage(taxonID: int):
    """
    Given a taxonID, returns a list of all taxonIDs in its lineage.
    """
    try:
        lineage = ncbi.get_lineage(taxonID)
    except ValueError:
        print("Cannot find taxonID " + str(taxonID))
        lineage = [taxonID]
    if lineage is None:
        lineage = [0]

    return lineage


def get_cannonical_lineage(
    taxonID: int,
    desired_levels: int = list(canonical_prefixes.keys()),
    prefix_dictionary: dict = canonical_prefixes,
):
    """
    taxonID: interger taxonID
    desired_levels: list of desired levels - e.g. ['superkingdom', 'species']
    prefix_dictionary: A dict with the desired levels as the keys and an abbreviation as
        the value. This abbreviation will be appended to the start of the taxonomy.
        E.g. 'species':'s__' will append s__ to the species name.

    This function outputs the lineage as the NAMES of each taxon.
    """

    for level in desired_levels:
        if level not in prefix_dictionary:
            prefix_dictionary[level] = ""

    # Get a list of lineage and a list of their ranks
    lineage = get_lineage(taxonID)
    levels = [get_level(taxon) for taxon in lineage]

    cannonical_lineage = []
    # Iterate over each of the levels in prefix dictionary
    for level, prefix in prefix_dictionary.items():

        if level not in desired_levels:
            continue

        # If the level isn't here, it is unknown
        if level not in set(levels):
            cannonical_lineage.append("")
            continue

        # Get the taxon name
        index = levels.index(level)
        taxon = lineage[index]
        name = get_name(taxon)

        # Report it out with the lineage
        cannonical_lineage.append(prefix + name)

    return cannonical_lineage


if __name__ == "__main__":
    msg = "This script has utilities and functions. Don't call it directly!"
    raise ValueError(msg)
