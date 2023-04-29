from .utils.Foldseek_Dataset import Foldseek_Dataset
from .utils.misc import talk_to_me, make_output_dir


def aln_cluster_main(args):

    talk_to_me("Parsing input alignment file")
    data = Foldseek_Dataset()
    data.parse_alignment(args.alignment_file, args.alignment_fields)

    talk_to_me("Generating clusters")
    clusters = []
    seen_members = set()
    for query, alignment_group in data.alignment_groups.items():
        for alignment in alignment_group.alignments:
            query = alignment.query
            target = alignment.target

            # Find current cluster for query and target. Will also set the defaults here
            query_cluster = {query}
            target_cluster = {target}
            for cluster in clusters:
                if query in cluster:
                    query_cluster = cluster
                if target in cluster:
                    target_cluster = cluster
                if query_cluster != {query} and target_cluster != {target}:
                    break

            # If they're alread in the same cluster, we're done here
            if query_cluster == target_cluster:
                continue

            # Need to merge their clusters
            merged_cluster = query_cluster.union(target_cluster)

            # If the query cluster or the target cluster is in the clusters, need to
            # remove it
            if query_cluster in clusters:
                clusters.remove(query_cluster)
            if target_cluster in clusters:
                clusters.remove(target_cluster)

            # Then replace it with the merged cluster
            clusters.append(merged_cluster)

            # Keep track of all seen members
            seen_members.add(query)
            seen_members.add(target)

    # If there is an all_inputs file: For those inputs that don't have an alignment,
    # add them as a single-member cluster.
    if args.all_inputs != "":
        talk_to_me("Making sure all alignment inputs are present in a cluster...")
        all_inputs = set()
        with open(args.all_inputs) as infile:
            for line in infile:
                line = line.rstrip("\n")
                all_inputs.add(line)

        not_seen = all_inputs - seen_members

        for member in not_seen:
            clusters.append(set([member]))

    talk_to_me("Generating output")
    out = ["cluster_rep", "cluster_member"]
    out = "\t".join(out) + "\n"
    for cluster in clusters:

        # Pick the cluster rep... it will be just the first member alphabetically. Will
        # make a sorted list so it's reproducible (compared to picking from a set)
        cluster_list = list(cluster)
        cluster_list.sort()
        cluster_rep = cluster_list[0]

        for cluster_member in cluster:
            out_line = [cluster_rep, cluster_member]
            out_line = "\t".join(out_line)
            out += out_line + "\n"

    make_output_dir(args.outfile)
    with open(args.outfile, "w") as outfile:
        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
