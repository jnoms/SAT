from .alignments import Alignment_group, Alignment_object
from .misc import talk_to_me
from .clusters import Cluster
from .ete3_taxonomy import Taxon, taxonID_list_to_lineage_counts

# ------------------------------------------------------------------------------------ #
# Classes
# ------------------------------------------------------------------------------------ #
class Foldseek_Dataset:
    """
    This class is designed to store foldseek alignment information from a foldseek
    tabular output file. Alignments are organized in the following way:
    Foldseek_Dataset --> Alignment_group --> Alignment_object

    Furthermore, cluster objects can also be stored in this object. The cluster
    objects will themselves contain references to the alignment_groups (and therefore
    specific alignment_objects) that contain queries in the cluster.
    Foldseek_Dataset --> Cluster --> Alignment_group --> Alignment_object

    Working upwards:
    - Alignment_object: contains all information from the foldseek alignment, with
      each column as a discrete slot. Contains functions for writing out the alignment.
      ---> Can also contain a Taxon object which has helpful methods for handling
            taxonomy information.
    - Alignment_group: Groups of alignment objects that have the same query. Thus, note
      that in all-by-all alignments a specific item will be present as both a query
      and a target (and thus be present in multiple alignment_groups).
    - Cluster: Contains alignment groups for each query that is present in the same
      cluster (determined by Foldseek cluster command).
    - Foldseek_Dataset: Groups all of the information together, with helpful methods
      for reading in data and such.
    """

    def parse_alignment(self, alignment_file_path, alignment_fields=""):
        """
        Given an input path, parses the alignment file into a dictionary of
        structure query:Alingment_group, where each Alignment_group contains
        a list of Alignment objects.

        Stores this dictionary in the alignment_groups slot of the Foldseek_Dataset
        object.
        """
        alignments_dict = dict()
        with open(alignment_file_path) as infile:

            # Check if the first line is a header - if so, and alignment_fields is
            # an empty string, use those as the alignment_fields
            if alignment_fields == "":
                for line in infile:
                    if line.startswith("query"):
                        alignment_fields = line.rstrip("\n").split("\t")
                        break
                    else:
                        msg = "alignment_fields has not been passed to parse_alignment,"
                        msg += " which is only allowed when the first line has headers!"
                        msg += " (e.g. first line should start with 'query')"
                        raise ValueError(msg)

            for line in infile:
                line = line.rstrip("\n").split("\t")

                if len(line) != len(alignment_fields):
                    msg = "The line and alignment_fields don't have the same "
                    msg += "number of entries!"
                    raise ValueError(msg)

                alignment = Alignment_object(line, alignment_fields)

                # Save to output dictionary
                if alignment.query not in alignments_dict:
                    alignments_dict[alignment.query] = Alignment_group(alignment.query)
                    alignments_dict[alignment.query].add_alignment(alignment)
                else:
                    alignments_dict[alignment.query].add_alignment(alignment)

        self.alignment_groups = alignments_dict

    def load_cluster_objects(self, cluster_file_path: str):
        """
        Given the path to the cluster file and a dictionary containing alignment_groups,
        will return a list of cluster_objects. Here, each cluster object contains the
        alignment_groups that have a cluster member as a query.

        Inputs:
        - cluster_file_path: path to the foldseek cluster file.
        - alignment_groups: dictionary of structure query:alignment_group

        Output:
        - List of cluster objects into the clusters slot
        """

        # structure - foldseek_cluster_rep:set of member names
        clusters = dict()

        with open(cluster_file_path) as infile:
            for line in infile:
                # Skip header if there is one
                if line.startswith("query"):
                    continue

                line = line.rstrip("\n").split("\t")
                foldseek_cluster_rep, member = line
                if foldseek_cluster_rep not in clusters:
                    clusters[foldseek_cluster_rep] = set()
                clusters[foldseek_cluster_rep].add(member)

        cluster_objects = []
        for foldseek_cluster_rep, cluster_members in clusters.items():

            # Skip clusters that only have one member (aka the cluster_rep).
            if len(cluster_members) == 1:
                continue

            cluster = Cluster(foldseek_cluster_rep)

            for member in cluster_members:

                # Look up the alignment_group for the cluster member - e.g.
                # all of the alignments with that member as the query
                member_alignment_group = self.alignment_groups[member]

                # Remove any alignments whose targets are not in the same cluster
                for i, alignment in enumerate(member_alignment_group.alignments):
                    if alignment.target not in cluster_members:
                        msg = (
                            f"Found a cluster with a query, {alignment.query}, who has"
                            f" an alignment against a target, {alignment.target}, that "
                            " is not in the same cluster. The foldseek_cluster_rep is "
                            f"{foldseek_cluster_rep}"
                        )
                        print(msg)
                        del member_alignment_group.alignments[i]

                cluster.alignment_groups.append(member_alignment_group)
                cluster.cluster_members.add(member)

            cluster_objects.append(cluster)

        self.clusters = cluster_objects

    def add_cluster_ID(self):
        """
        This function orders the clusters by the number of items they have (aka the
        number of alignment_groups, aka the number of queries) and adds a cluster_ID
        slot. Lower numbers indicate more members of the cluster. cluster_ID starts at
        1.

        This function ACTS IN PLACE, and modifies each cluster object in the input.
        """

        # Make a list of sets of structure [(cluster_object, count_of_alignments), ...]
        unordered = [
            (cluster, len(cluster.alignment_groups)) for cluster in self.clusters
        ]

        # Order it and then add the cluster_ID
        ordered = sorted(unordered, key=lambda x: x[1])[::-1]
        for i, (cluster, count) in enumerate(ordered):
            i += 1
            cluster.cluster_ID = i

    def write_out_cluster_alignments(
        self,
        alignment_fields,
        top_or_nonredundant="both",
        cluster_fields=["cluster_ID", "cluster_count", "top_query"],
    ):
        """
        alignment_fields are the input alignments fields present in each
        alignment_object. cluster_fields are the fields in the Cluster object

        top_or_nonredundant: options are 'top', 'nr' (non-redundant), or 'both'.
        Dictates if will be writing out only the alignments for the top query in each
        cluster, all non-redundant alignments, or both.

        The output will be a dict with the keys top and/or nr depending on what is
        selected.
        """
        if top_or_nonredundant not in ["top", "nr", "both"]:
            msg = "top_or_nonredundant must be top, nr, or both. You "
            msg += f"entered {top_or_nonredundant}."
            raise ValueError(msg)

        # Keeping track...
        total = len(self.clusters)
        progress = 0

        out = dict()
        out["top"] = ""
        out["nr"] = ""

        for cluster in self.clusters:

            # Adding the top_query and the cluser count
            cluster.add_top_query()
            cluster.add_cluster_count()

            # Writing out the alignments of the top query for each cluster - e.g. the
            # query with the most alignments or, if a tie, the query with the highest
            # average TMscore.
            if top_or_nonredundant in ["top", "both"]:
                out["top"] += cluster.write_top_query_alignments(
                    alignment_fields, cluster_fields
                )

            # Writing out all non-redundant alignments that make it into a cluster
            # object. Note that alignments with a query and target not in the same
            # cluster are excluded.
            if top_or_nonredundant in ["nr", "both"]:
                out["nr"] += cluster.write_all_nonredundant_alignments(
                    alignment_fields, cluster_fields
                )

            # Update progress
            progress += 1
            if progress % 100 == 0:
                talk_to_me(f"Progress: {progress}/{total}")

        return out

    def add_taxon_to_alignments(self, taxonomy_levels, delimiter="__", pos=-1):
        """
        Iterates through all alignment_groups and all alignment_objects and adds the
        query and target taxonIDs as a Taxon object in the query_taxonID and
        target_taxonID slots of the alignment_objects.
        """

        # Keep track of Taxa I've already made
        seen_taxa = dict()

        for query, alignment_group in self.alignment_groups.items():
            for alignment in alignment_group.alignments:
                alignment.add_query_taxonID(delimiter=delimiter, pos=pos)
                alignment.add_target_taxonID(delimiter=delimiter, pos=pos)

                if alignment.query_taxonID in seen_taxa:
                    alignment.query_taxon = seen_taxa[alignment.query_taxonID]
                else:
                    alignment.query_taxon = Taxon(
                        alignment.query_taxonID, taxonomy_levels
                    )
                    seen_taxa[alignment.query_taxonID] = alignment.query_taxon

                if alignment.target_taxonID in seen_taxa:
                    alignment.target_taxon = seen_taxa[alignment.target_taxonID]
                else:
                    alignment.target_taxon = Taxon(
                        alignment.target_taxonID, taxonomy_levels
                    )
                    seen_taxa[alignment.target_taxonID] = alignment.target_taxon

    def write_out_alignments(self, alignment_fields):
        """
        This differs from write_out_cluster_alignments because that one will iterate
        over cluster objects, while this one iterates directly over the alignments.

        This returns a string with the output from all of the alignment_objects, where
        alignment_fields are the slots (in the correct order!) that you want in the
        output. Output will be delimited by \\t.
        """
        out = ""
        for query, alignment_group in self.alignment_groups.items():
            for alignment in alignment_group.alignments:
                out += alignment.write_output(alignment_fields)
        return out

    def write_out_cluster_taxa_count(self, taxonomy_levels):
        """
        This function iterates through each alignment and counts the taxons at every
        level in each cluster. Result is a tab-delimited string of columns cluster_ID,
        top_query, level, taxon, count.

        Notably, this function requires that each alignment has the slots
        cluster_ID, top_query, query_taxonID, and target_taxonID. Typically this
        will come by parsing an alignment file that has been processed through
        process_clusters and then add_taxonomy_to_alignments.
        """

        cluster_taxa_dict = dict()
        for query, alignment_group in self.alignment_groups.items():
            for alignment in alignment_group.alignments:
                cluster = alignment.cluster_ID
                top_query = alignment.top_query
                key = tuple([cluster, top_query])

                query_taxonID = alignment.query_taxonID
                target_taxonID = alignment.target_taxonID

                if key not in cluster_taxa_dict:
                    cluster_taxa_dict[key] = set()

                cluster_taxa_dict[key].add(query_taxonID)
                cluster_taxa_dict[key].add(target_taxonID)

        # Outputing observed_taxa each time to avoid repetitive taxonID lookups, which
        # can be very slow.
        observed_taxa = dict()
        out = "\t".join(["cluster_ID", "top_query", "level", "taxon", "count"]) + "\n"
        for (cluster_ID, top_query), taxonIDs in cluster_taxa_dict.items():
            counts_dict, observed_taxa = taxonID_list_to_lineage_counts(
                taxonIDs, taxonomy_levels, observed_taxa
            )
            for level, c in counts_dict.items():
                for taxon, count in c.items():
                    line = (
                        "\t".join(
                            [cluster_ID, top_query, level, str(taxon), str(count)]
                        )
                        + "\n"
                    )
                    out += line

        return out


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #

if __name__ == "__main__":
    msg = "This script has utilities and functions. Don't call it directly!"
    raise ValueError(msg)
