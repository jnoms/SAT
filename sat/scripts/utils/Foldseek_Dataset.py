from .alignments import Alignment_group, Alignment_object
from .misc import talk_to_me
from .clusters import Cluster
from .ete3_taxonomy import Taxon


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

        Note that this function checks for a header in the input file. If a header is
        present and alignment_fields is not provided, it will use the header as the
        alignment fields and will store the alignment_fields list in the
        Foldseek_Dataset object's input_alignment_fields attribute.
        """
        # Format alignment fields from comma-delimited string to list if neceesary
        if alignment_fields != "" and not isinstance(alignment_fields, list):
            alignment_fields = alignment_fields.split(",")

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

            # Record alignment fields for later use
            self.input_alignment_fields = alignment_fields

            for line in infile:
                # It's okay if there is a header, but need to remove it
                if line.startswith("query"):
                    continue

                line = line.rstrip("\n").split("\t")

                if len(line) != len(alignment_fields):
                    msg = "The line and alignment_fields don't have the same "
                    msg += "number of entries!"
                    msg += f"Current line is: {line}.\n"
                    msg += f"Alignment fields are: {alignment_fields}"
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

            # # Skip clusters that only have one member (aka the cluster_rep).
            # if len(cluster_members) == 1:
            #     continue

            cluster = Cluster(foldseek_cluster_rep)

            for member in cluster_members:

                # Look up the alignment_group for the cluster member - e.g.
                # all of the alignments with that member as the query
                try:
                    member_alignment_group = self.alignment_groups[member]
                    cluster.alignment_groups.append(member_alignment_group)
                    cluster.cluster_members.add(member)
                except KeyError:
                    # Sometimes alignments have been filtered out before we get here -
                    # this can happen if all the alignments with a query of this cluster
                    # member were filtered out.
                    continue

            # only want to keep cluster objects that have alignment groups. If we're
            # using a cluster file from one foldseek run to add cluster information to
            # another foldseek run that was against a different database, it's possible
            # that some clusters won't have an alignment against any cluster member.
            if cluster.alignment_groups != []:
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
        rep_or_all="both",
        cluster_fields=["cluster_ID", "cluster_count", "top_query"],
    ):
        """
        alignment_fields are the input alignments fields present in each
        alignment_object. cluster_fields are the fields in the Cluster object

        rep_or_nonredundant: options are 'top', 'all', or 'both'.
        Dictates if will be writing out only the alignments for the top query in each
        cluster, all non-redundant alignments, or both.

        The output will be a dict with the keys top and/or nr depending on what is
        selected.
        """
        if rep_or_all not in ["top", "all", "both"]:
            msg = "top_or_nonredundant must be top, all, or both. You "
            msg += f"entered {rep_or_all}."
            raise ValueError(msg)

        # Keeping track...
        total = len(self.clusters)
        progress = 0

        out = dict()
        out["top"] = ""
        out["all"] = ""

        for cluster in self.clusters:

            # Adding the top_query and the cluser count
            cluster.add_cluster_count()

            # Writing out the alignments of the top query for each cluster - e.g. the
            # query with the most alignments or, if a tie, the query with the highest
            # average TMscore.
            if rep_or_all in ["top", "both"]:
                out["top"] += cluster.write_rep_alignments(
                    alignment_fields, cluster_fields
                )

            # Writing out all non-redundant alignments that make it into a cluster
            # object. Note that alignments with a query and target not in the same
            # cluster are excluded.
            if rep_or_all in ["all", "both"]:
                out["all"] += cluster.write_all_alignments(
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

    def load_cluster_objects_from_alignments(self):
        """
        This function iterates through all alignments and generates cluster objects that
        contain all members of a cluster. This requires that alignments already have the
        attributes (e.g. column in the alignment file) cluster_ID, cluster_count, and
        top_query.

        Cluster objects will be added to the Foldseek_dataset object as a list of
        cluster objects in the .clusters attribute.
        """
        required_attributes = ["cluster_ID", "cluster_count", "top_query"]
        cluster_lookup_dict = dict()

        for query, alignment_group in self.alignment_groups.items():
            for alignment in alignment_group.alignments:

                for a in required_attributes:
                    if not hasattr(alignment, a):
                        msg = (
                            f"Cannot find the attribute {a} in the alignments. "
                            f"Required attributes are {required_attributes} and "
                            "these must be present as columns in the input file."
                        )
                        raise ValueError(msg)

                cluster_ID = alignment.cluster_ID
                cluster_count = alignment.cluster_count
                top_query = alignment.top_query

                # Initialize a cluster object
                if cluster_ID not in cluster_lookup_dict:
                    cluster_lookup_dict[cluster_ID] = Cluster(cluster_ID)
                    cluster_lookup_dict[cluster_ID].cluster_count = cluster_count
                    cluster_lookup_dict[cluster_ID].top_query = top_query
                    cluster_lookup_dict[cluster_ID].taxon_set = set()
                    cluster_lookup_dict[cluster_ID].alignments = []

                # Store the alignment and have a set of taxons that are in the cluster
                cluster_lookup_dict[cluster_ID].alignments.append(alignment)
                cluster_lookup_dict[cluster_ID].taxon_set.add(alignment.query_taxon)
                cluster_lookup_dict[cluster_ID].taxon_set.add(alignment.target_taxon)

        clusters = [cluster_object for _, cluster_object in cluster_lookup_dict.items()]
        self.clusters = clusters

    def merge(self, other):
        """
        other should be another Foldseek_dataset.

        This is a simple function that reads in all alignment_groups/alignments from
        other and adds them to this foldseek dataset. If a query is already present
        in this foldseek dataset it's alignments will be added to the existing alignment
        group for that query.
        """

        for query, other_alignment_group in other.alignment_groups.items():

            # If the query isn't already present, add it and continue
            if query not in self.alignment_groups:
                self.alignment_groups[query] = other_alignment_group
                continue

            # If the query is already present, add all of the alignments to the
            # pertinent alingment_group of that query
            if query in self.alignment_groups:
                self.alignment_groups[query].alignments.extend(
                    other_alignment_group.alignments
                )

    def count_alignments(self):
        """
        Returns an interger indicating the number of alignments in this dataset
        """
        aln_count = 0
        for _, alignment_group in self.alignment_groups.items():
            aln_count += len(alignment_group.alignments)
        return aln_count

    def add_field_to_alignments(self, field: str, value: any):
        """
        Adds the field to each alignment object with the indicated value. Typical use-
        case is labeling alignments that come from different sources when merging
        multiple Foldseek_datasets.
        """
        for _, alignment_group in self.alignment_groups.items():
            for alignment in alignment_group.alignments:
                alignment.__dict__[field] = value


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def compare_foldseek_datasets(dataset1, dataset2):
    """
    For testing. This compares the contents and alignments for two datasets.
    """
    try:
        assert len(dataset1.alignment_groups) == len(dataset1.alignment_groups)
    except AssertionError:
        msg = (
            "The two datasets have different numbers of alignment groups."
            f"Dataset1: {len(dataset1.alignment_groups)} alignment groups, "
            f"Dataset2: {len(dataset1.alignment_groups)} alignment groups"
        )
        raise AssertionError(msg)

    for query, alignment_group in dataset2.alignment_groups.items():
        try:
            assert query in dataset1.alignment_groups
        except AssertionError:
            msg = (
                f"Cannot find the dataset2 query, {query}, as a query in dataset1 "
                "alignment groups"
            )
            raise AssertionError(msg)
        for alignment in alignment_group.alignments:
            try:
                assert alignment in dataset1.alignment_groups[query].alignments
            except AssertionError:
                msg = (
                    f"Dataset2 has an alignment, query - {alignment.query}, "
                    f"target - {alignment.target} in the alignment group with query "
                    f"{query}. However, this alignment cannot be found in dataset1's "
                    f"alignment group with query {query}."
                )
                raise AssertionError(msg)

            dataset1_alignment = [
                a for a in dataset1.alignment_groups[query].alignments if a == alignment
            ][0]
            assert alignment.__eq__(dataset1_alignment)


if __name__ == "__main__":
    msg = "This script has utilities and functions. Don't call it directly!"
    raise ValueError(msg)
