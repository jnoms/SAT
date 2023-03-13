from .utils.Foldseek_Dataset import Foldseek_Dataset
from .utils.misc import make_output_dir, talk_to_me
from .utils.clusters import Cluster_information
from .utils.ete3_taxonomy import taxon_list_to_lineage_counts, Taxon


def format_args(args):
    args.taxonomy_levels = args.taxonomy_levels.split(",")
    return args


def get_taxonID(member, delimiter, pos):
    """
    member is a string, which is probably a domain or structure name. This pulls out the
    taxonID, which presumably is stored somewhere in the file name. The delimiter and
    pos indicate what to split the member string by and where to find the taxonID. This
    function returns the taxonID
    """
    taxonID = member.split(delimiter)[pos].rstrip(".pdb")
    while "_" in taxonID:
        taxonID = taxonID[:-1]
    return taxonID


def get_cluster_rep_to_taxonIDs(clusters, taxonID_finder_delimiter, taxonID_finder_pos):
    """
    This takes clusters, the Cluster_information object, and iterates through each
    member of the cluster. For each member, it uses the delimiter and position to
    pull out the taxonID. This function then creates a dictionary of structure
    cluster_ID:set of taxonIDs present in that cluster
    """
    cluster_rep_to_taxonIDs = dict()
    for rep, cluster in clusters.clusters.items():
        cluster_taxids = set()
        cluster_rep = cluster.cluster_rep
        cluster_members = cluster.cluster_members
        for member in cluster_members:
            taxonID = get_taxonID(
                member.cluster_member, taxonID_finder_delimiter, taxonID_finder_pos
            )
            cluster_taxids.add(taxonID)

        if cluster_rep in cluster_rep_to_taxonIDs:
            msg = f"I alreadu have cluster_rep {cluster_rep} in cluster_rep_to_taxonIDs!! "
            msg += "weird! Something is wrong."
            raise ValueError(msg)

        cluster_rep_to_taxonIDs[cluster_rep] = cluster_taxids
    return cluster_rep_to_taxonIDs


def aln_taxa_counts_main(args):

    args = format_args(args)

    talk_to_me("Reading in cluster file")
    clusters = Cluster_information()
    clusters.parse_cluster_file(args.cluster_file, args.cluster_file_fields)
    clusters.generate_member_and_rep_to_cluster_ID()

    talk_to_me("Extracting taxonIDs")
    cluster_rep_to_taxonIDs = get_cluster_rep_to_taxonIDs(
        clusters, args.taxonID_finder_delimiter, args.taxonID_finder_pos
    )

    # If there is an alignment file, will all the target taxonIDs to the cluster_ID of
    # each query cluster_ID.
    if args.alignment_file != "":
        talk_to_me("Parsing alignment file")
        data = Foldseek_Dataset()
        data.parse_alignment(args.alignment_file, args.alignment_fields)
        data.add_taxon_to_alignments(args.taxonomy_levels)

        # Find each query's cluster_ID, and then add each target taxonID to the
        # cluster_ID_to_taxonIDs
        for query, alignment_group in data.alignment_groups.items():
            query_cluster_rep = clusters.cluster_member_to_rep[query]
            for alignment in alignment_group.alignments:
                cluster_rep_to_taxonIDs[query_cluster_rep].add(alignment.target_taxonID)

    # Convert the structure from cluster_ID:set of taxonIDs TO
    # cluster_ID:set of Taxon objects.
    talk_to_me("Generating Taxon objects")
    taxonID_to_taxon = dict()
    cluster_rep_to_taxons = dict()
    for cluster_rep, taxonID_list in cluster_rep_to_taxonIDs.items():
        for taxonID in taxonID_list:
            if taxonID in taxonID_to_taxon:
                taxon_object = taxonID_to_taxon[taxonID]
            else:
                taxon_object = Taxon(taxonID, args.taxonomy_levels)
                taxonID_to_taxon[taxonID] = taxon_object

            if cluster_rep not in cluster_rep_to_taxons:
                cluster_rep_to_taxons[cluster_rep] = set()

            cluster_rep_to_taxons[cluster_rep].add(taxon_object)

    # Generate the output
    talk_to_me("Generating output")
    out = [
        "cluster_ID",
        "cluster_rep",
        "cluster_count",
        "superkingdom",
        "level",
        "taxon",
        "count",
    ]
    out = "\t".join(out) + "\n"
    for cluster_rep, taxon_list in cluster_rep_to_taxons.items():
        cluster_ID = clusters.cluster_rep_to_ID[cluster_rep]
        cluster_count = len(clusters.clusters[cluster_rep].cluster_members)
        count_dict, taxa_to_superkingdom = taxon_list_to_lineage_counts(
            taxon_list, args.taxonomy_levels
        )
        for level, taxa_dict in count_dict.items():
            for taxa, count in taxa_dict.items():
                superkingdom = taxa_to_superkingdom[taxa]
                out_line = [
                    cluster_ID,
                    cluster_rep,
                    cluster_count,
                    superkingdom,
                    level,
                    taxa,
                    count,
                ]
                out_line = [str(i) for i in out_line]
                out_line = "\t".join(out_line) + "\n"
                out += out_line

    make_output_dir(args.output_file)
    with open(args.output_file, "w") as outfile:
        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
