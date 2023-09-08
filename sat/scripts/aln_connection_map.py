from collections import Counter

from .utils.clusters import Cluster_information
from .utils.misc import talk_to_me, make_output_dir, print_args


def aln_connection_map_main(args):

    print_args(args)

    talk_to_me("Parsing cluster file")
    cluster_info = Cluster_information()
    cluster_info.parse_cluster_file(args.cluster_file, args.cluster_file_fields)

    if "family" not in cluster_info.input_cluster_file_fields:
        msg = "Can't find family in the cluster file fields. Does "
        msg += "the cluster file have taxonomy information?"
        raise ValueError(msg)

    if "species" not in cluster_info.input_cluster_file_fields:
        msg = "Can't find species in the cluster file fields. Does "
        msg += "the cluster file have taxonomy information?"
        raise ValueError(msg)

    talk_to_me("Recording connections")
    fam_connections = Counter()
    for cluster_rep, cluster in cluster_info.clusters.items():

        # Only want to count one connection per cluster
        seen_connections = set()

        for m1 in cluster.cluster_members:
            for m2 in cluster.cluster_members:

                m1_fam = m1.family
                m2_fam = m2.family

                if m1_fam == m2_fam:
                    continue

                if m1_fam == "" or m2_fam == "":
                    continue

                connection = [m1_fam, m2_fam]
                connection.sort()
                connection = frozenset(connection)

                if connection in seen_connections:
                    continue

                seen_connections.add(connection)
                fam_connections[connection] += 1

    talk_to_me("Counting how many clusterIDs are in each family")
    fam_cluster_totals = Counter()
    for cluster_rep, cluster in cluster_info.clusters.items():

        # Record how many species each fam has in this cluster
        fam_species_counts = dict()

        for m in cluster.cluster_members:

            fam = m.family
            species = m.species

            if fam == "":
                continue

            if fam not in fam_species_counts:
                fam_species_counts[fam] = set()
            fam_species_counts[fam].add(species)

        for fam, species in fam_species_counts.items():
            count = len(species)
            fam_cluster_totals[fam] += 1

    talk_to_me("Generating output connections file")
    out = ["f1", "f2", "count", "f1_total", "f2_total", "jaccard"]
    out = "\t".join(out) + "\n"

    for entry in fam_connections.most_common():

        connection, count = entry
        f1, f2 = connection
        f1_total = fam_cluster_totals[f1]
        f2_total = fam_cluster_totals[f2]
        jarccard = count / (f1_total + f2_total - count)

        out_line = [f1, f2, count, f1_total, f2_total, jarccard]
        out_line = [str(x) for x in out_line]
        out_line = "\t".join(out_line) + "\n"
        out += out_line

    make_output_dir(args.outfile)
    with open(args.outfile, "w") as outfile:
        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
