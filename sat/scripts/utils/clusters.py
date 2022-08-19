from collections import Counter


# ------------------------------------------------------------------------------------ #
# Classes
# ------------------------------------------------------------------------------------ #
class Cluster:
    def __init__(self, id):
        self.id = id
        self.alignment_groups = []
        self.cluster_members = set()

    def add_top_query(self):
        """
        This function goes through each query's alignment group. It finds the query with
        the highest number of alignments. If there are multiple queries with the same
        number of alignments, the one with the highest average TMscore is chosen. The
        name of the query is added to the self.top_query slot of the cluster object.
        """
        query_counts = Counter()
        query_avg_TMscore = dict()
        for alignment_group in self.alignment_groups:
            query = alignment_group.query

            # If a query doesn't have any alignments, this means that is was actually a
            # self alignment that was removed. But, in any case, want to continue
            if alignment_group.alignments == []:
                continue

            # Count the # of alignmnets the query has
            if query in query_counts:
                msg = (
                    "There seems to be multiple of the same query in the"
                    "alignment_group!"
                )
                raise ValueError(msg)

            query_counts[query] = len(alignment_group.alignments)

            # Get the average TMscore of the query's alignments
            query_avg_TMscore[query] = 0
            for alignment in alignment_group.alignments:
                query_avg_TMscore[query] += float(alignment.alntmscore)
            query_avg_TMscore[query] = query_avg_TMscore[query] / query_counts[query]

        # If no alignment_groups had alignments, I assume there is only one
        # alignment_group and that the group doesn't have any alignments. Thus, the
        # representative is the query
        if query_counts == Counter():
            if len(self.alignment_groups) != 1:
                msg = "There are no alignment_groups with alignments, yet "
                msg += "there are multiple alignments. Weird!"
                raise ValueError(msg)
            if self.alignment_groups[0].query != self.id:
                msg = "There are no alignment_groups with alignments, yet "
                msg += "the only alignment_group has a query that is not the cluster ID"
                msg += ". Weird!"
                raise ValueError(msg)
            self.top_query = self.id
            return

        # Now, want to pick the query with the most alignments. If multiple queries have
        # the same # alignments, pick the one with the highest average TMscore.
        max_count = query_counts.most_common(1)[0][1]
        representative = ""
        representative_avg_TM = 0
        for query, query_count in query_counts.most_common():

            # Can stop the loop if the count is not max_counts
            if query_count < max_count:
                break

            # At this point we are a query that has the max_counts (there may be
            # multiple). Need to compare with average TMscore.
            query_avg_TM = query_avg_TMscore[query]
            if query_avg_TM > representative_avg_TM:
                representative = query
                representative_avg_TM = query_avg_TM

        self.top_query = representative

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

    def write_top_query_alignments(self, alignment_fields, cluster_fields=[]):
        """
        Writes the alignments of the top query as a string.

        The alignment_fields come from the alignment objects present in the
        alignment groups in the cluster. The cluster_fields come directly from the
        cluster object.
        """
        top_query = self.top_query
        out = ""
        for alignment_group in self.alignment_groups:
            if alignment_group.query == top_query:
                for alignment in alignment_group.alignments:
                    if cluster_fields != []:
                        alignment_out = alignment.write_output(alignment_fields)
                        cluster_out = self.write_output(cluster_fields)
                        out += alignment_out.rstrip("\n") + "\t" + cluster_out
                    else:
                        alignment_out = alignment.write_output(alignment_fields)
                        out += alignment_out
        return out

    def write_all_nonredundant_alignments(self, alignment_fields, cluster_fields=[]):
        """
        Writes all alignments as a string, but only one alignment per
        query-target pair (in all-by-all- searches, each item will be listed as a
        query and will match its partner as a target... so each alignment will be
        present twice!). Also note that upon loading into alignment_objects self
        alignments were removed.

        The alignment_fields come from the alignment objects present in the
        alignment groups in the cluster. The cluster_fields come directly from the
        cluster object.
        """
        out = ""
        seen = set()
        for alignment_group in self.alignment_groups:
            for alignment in alignment_group.alignments:
                lookup = frozenset([alignment.query, alignment.target])
                if lookup in seen:
                    continue
                else:
                    seen.add(lookup)
                    if cluster_fields != []:
                        alignment_out = alignment.write_output(alignment_fields)
                        cluster_out = self.write_output(cluster_fields)
                        out += alignment_out.rstrip("\n") + "\t" + cluster_out
                    else:
                        alignment_out = alignment.write_output(alignment_fields)
                        out += alignment_out
        return out


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def load_cluster_objects(cluster_file_path: str, alignment_groups: dict):
    """
    Given the path to the cluster file and a dictionary containing alignment_groups,
    will return a list of cluster_objects. Here, each cluster object contains the
    alignment_groups that have a cluster member as a query.

    Inputs:
    - cluster_file_path: path to the foldseek cluster file.
    - alignment_groups: dictionary of structure query:alignment_group

    Output:
    - list of cluster objects
    """

    # structure - foldseek_cluster_rep:set of member names
    clusters = dict()

    with open(cluster_file_path) as infile:
        for line in infile:
            line = line.rstrip("\n").split("\t")
            foldseek_cluster_rep, member = line
            if foldseek_cluster_rep not in clusters:
                clusters[foldseek_cluster_rep] = set()
            clusters[foldseek_cluster_rep].add(member)

    cluster_objects = []
    for foldseek_cluster_rep, cluster_members in clusters.items():
        cluster = Cluster(foldseek_cluster_rep)

        for member in cluster_members:
            member_alignment_group = alignment_groups[member]
            cluster.alignment_groups.append(member_alignment_group)
            cluster.cluster_members.add(member)

        cluster_objects.append(cluster)

    return cluster_objects


def add_cluster_ID(cluster_objects: list):
    """
    Input: cluster_objects, a list of cluster objects.

    This function orders them by the number of items they have (aka the number of
    alignment_groups, aka the number of queries) and adds a cluster_ID slot. Lower
    numbers indicate more members of the cluster. cluster_ID starts at 1.

    This function ACTS IN PLACE, and modifies each cluster object in the input.
    """

    # Make a list of sets of structure [(cluster_object, count_of_alignments), ...]
    unordered = [
        (cluster, len(cluster.alignment_groups)) for cluster in cluster_objects
    ]

    # Order it and then add the cluster_ID
    ordered = sorted(unordered, key=lambda x: x[1])[::-1]
    for i, (cluster, count) in enumerate(ordered):
        i += 1
        cluster.cluster_ID = i


if __name__ == "__main__":
    msg = "This script has utilities and functions. Don't call it directly!"
    raise ValueError(msg)
