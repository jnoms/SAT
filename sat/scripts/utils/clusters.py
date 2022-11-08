from collections import Counter

from .ete3_taxonomy import taxon_list_to_lineage_counts


# ------------------------------------------------------------------------------------ #
# Classes
# ------------------------------------------------------------------------------------ #
class Cluster:
    def __init__(self, cluster_rep):
        self.cluster_rep = cluster_rep
        self.alignment_groups = []
        self.cluster_members = set()

    def __str__(self):
        return self.cluster_rep

    def __repr__(self):
        return self.cluster_rep

    def add_cluster_count(self):
        """
        Counts the number of queries in the cluster
        """
        self.cluster_count = len(self.alignment_groups)

    def write_output(self, fields: list, sep: str = "\t"):
        """
        Returns the desired fields as a tab-delimited string followed by a \\n.
        """
        out = []
        for field in fields:
            val = self.__dict__[field]
            if type(val) == list:
                val = [str(v) for v in val]
                val = sep.join(val)
            elif type(val) in {float, int}:
                val = str(val)
            out.append(val)
        return sep.join(out) + "\n"

    def write_rep_alignments(self, alignment_fields, cluster_fields=[]):
        """
        Writes the alignments of the foldseek cluster representative as a string.

        The alignment_fields come from the alignment objects present in the
        alignment groups in the cluster. The cluster_fields come directly from the
        cluster object.
        """
        rep = self.cluster_rep
        out = ""
        for alignment_group in self.alignment_groups:
            if alignment_group.query == rep:
                for alignment in alignment_group.alignments:
                    if cluster_fields != []:
                        alignment_out = alignment.write_output(alignment_fields)
                        cluster_out = self.write_output(cluster_fields)
                        out += alignment_out.rstrip("\n") + "\t" + cluster_out
                    else:
                        alignment_out = alignment.write_output(alignment_fields)
                        out += alignment_out
        return out

    def write_all_alignments(self, alignment_fields, cluster_fields=[]):
        """
        Writes all alignments as a string. Note that upon loading into alignment_objects
        self alignments were removed.

        The alignment_fields come from the alignment objects present in the
        alignment groups in the cluster. The cluster_fields come directly from the
        cluster object.
        """
        out = ""
        for alignment_group in self.alignment_groups:
            for alignment in alignment_group.alignments:
                if cluster_fields != []:
                    alignment_out = alignment.write_output(alignment_fields)
                    cluster_out = self.write_output(cluster_fields)
                    out += alignment_out.rstrip("\n") + "\t" + cluster_out
                else:
                    alignment_out = alignment.write_output(alignment_fields)
                    out += alignment_out
        return out

    def add_taxa_counts(self, taxonomy_levels):
        """
        Looks for the cluster object's taxon_set attribute and yields counts and
        superkingdoms. This function loads the .taxa_count_dict and
        .taxa_superkingdom_dict attributes.

        taxon_set is simply a set of taxon objects that are present in any alignment in
        the cluster.
        """
        if not hasattr(self, "taxon_set"):
            msg = (
                "This cluster does not have the taxon_set attribute. This is required."
            )
            raise ValueError(msg)

        taxon_count_dict, taxon_superkingdom_dict = taxon_list_to_lineage_counts(
            self.taxon_set, taxonomy_levels
        )
        self.taxon_count_dict = taxon_count_dict
        self.taxon_superkingdom_dict = taxon_superkingdom_dict

    def write_taxa_counts_output(self):
        """
        This outputs a tab-delimited string with the entires cluster_ID, cluster_rep,
        level, superkingdom, taxon, and count.
        """
        out = ""
        for level, c in self.taxon_count_dict.items():
            for taxon, count in c.items():
                cluster_ID = self.cluster_ID
                superkingdom = self.taxon_superkingdom_dict[taxon]
                cluster_rep = self.cluster_rep
                line = (
                    "\t".join(
                        [
                            cluster_ID,
                            cluster_rep,
                            level,
                            superkingdom,
                            str(taxon),
                            str(count),
                        ]
                    )
                    + "\n"
                )
                out += line

        return out


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
