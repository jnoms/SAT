import pickle


from .utils.misc import talk_to_me, make_output_dir
from .utils.uniprot import (
    uniprot_object,
    submit_id_mapping,
    check_id_mapping_results_ready,
    get_id_mapping_results_link,
    get_id_mapping_results_search,
)

# ------------------------------------------------------------------------------------ #
# Constants
# ------------------------------------------------------------------------------------ #
output_delimiter = "\t"
output_fields = ["accession", "geneName", "fullName"]


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def read_in_uniprot_IDs(infile, infile_col=0):
    """
    Takes in the path to a file and extracts the uniprot IDs. If the IDs are not in the
    first column, specify infile_col.

    infile_col is 0-indexed!

    This script will remove unneccesary stuff from the front and end of each entry by
    keeping the second item delimited by a dash (-)... e.g. AF-K0EZQ3-F1-model_v2.pdb.gz
    will turn into K0EZQ3
    """

    IDs = set()
    with open(infile) as infile:
        for line in infile:
            if line.startswith("query"):
                continue
            line = line.split("\t")
            ID = line[infile_col]
            if "-" in ID:
                ID = ID.split("-")[1]
            IDs.add(ID)

    return IDs


def reformat_uniprot_lookup_results(results):
    """
    Currently, the results is a terrible dictionary of structure
    results["results"] --> [{from:{...}, to:{...}}, ....]

    Where results["results"][i]["from"] is the uniprotID.

    I will reformat the dict such that is is structured as

    uniprotID:{the dict that was initially in "to"}
    """

    reformatted = dict()

    for i in range(len(results["results"])):

        res = results["results"][i]

        ID = res["from"]
        to = res["to"]
        reformatted[ID] = to

    return reformatted


def parse_uniprot_results(results):
    """
    Takes in the results highly-nested dictionary from uniprot parsing, and returns
    a dictionary of accession:uniprot_object.

    The uniprot_object will be loaded with the slots
    -accession
    -fullName
    -geneName

    If a slot can't be parsed, it will be entered as a blank string.
    """
    res = dict()
    for acc, results_dict in results.items():

        # Try various lookups to get the gene name
        fullName = ""
        if fullName == "":
            try:
                fullName = results_dict["proteinDescription"]["recommendedName"][
                    "fullName"
                ]["value"]
            except KeyError:
                pass
        if fullName == "":
            try:
                fullName = results_dict["proteinDescription"]["submissionNames"][0][
                    "fullName"
                ]["value"]
            except KeyError:
                pass

        # Look up the geneID
        geneName = ""
        if geneName == "":
            try:
                geneName = results_dict["genes"][0]["geneName"]["value"]
            except KeyError:
                pass
        if geneName == "":
            try:
                geneName = results_dict["genes"][0]["orfNames"][0]["value"]
            except KeyError:
                pass

        # Generate a uniprot object
        u = uniprot_object(acc)
        u.fullName = fullName
        u.geneName = geneName

        res[acc] = u

    return res


def query_uniprot_main(args):
    talk_to_me("Reading in file to get uniprot IDs")
    IDs = read_in_uniprot_IDs(args.infile, args.infile_col)

    if args.uniprot_cache != "":

        # Read in cache
        try:
            talk_to_me("Reading in cache file")
            cache = pickle.load(open(args.uniprot_cache, "rb"))
            talk_to_me(f"Found {len(cache.keys())} entires in the cache.")
        except FileNotFoundError:
            msg = (
                f"There is currently no cache file {args.uniprot_cache}."
                " So, a cache will be made at this path after this script runs."
            )
            talk_to_me(msg)
            cache = dict()

        # Remove ID's that are already in the cache
        IDs = IDs - cache.keys()

    # Look up IDs!
    if len(IDs) > 0:
        talk_to_me(f"Looking up {len(IDs)} IDs from uniprot")
        job_id = submit_id_mapping(
            from_db="UniProtKB_AC-ID", to_db="UniProtKB", ids=IDs
        )

        if check_id_mapping_results_ready(job_id):
            link = get_id_mapping_results_link(job_id)
            results = get_id_mapping_results_search(link)

        talk_to_me("Reformatting the uniprot results")
        results = reformat_uniprot_lookup_results(results)
    else:
        if args.uniprot_cache == "":
            msg = "For some reason, there are no IDs present despite there being no"
            msg += " cache file. This isn't right!"
            raise ValueError(msg)

        talk_to_me("All IDs are in the cache, so not querying uniprot.")
        results = dict()

    if args.uniprot_cache != "":

        # Add the cache to the new results dict
        talk_to_me("Combining the results with the cache")
        results.update(cache)

        # Write the entire thing as a new cache if IDs were added
        if len(IDs) > 0:
            talk_to_me("Writing a new cache")
            with open(args.uniprot_cache, "wb") as cache_handle:
                pickle.dump(results, cache_handle)

    talk_to_me("Generating uniprot objects from the results")
    uniprot_object_dict = parse_uniprot_results(results)

    talk_to_me("Writing flat output file")
    make_output_dir(args.uniprot_lookup_output)
    with open(args.uniprot_lookup_output, "w") as outfile:
        out = ""
        for _, uniprot in uniprot_object_dict.items():
            out += uniprot.write_out(output_delimiter, output_fields)

        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
