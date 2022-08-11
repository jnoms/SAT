#!/usr/bin/env python

import argparse


def main():

    # Top-level parser
    parser = argparse.ArgumentParser(
        description=(
            "SAT - Structural Analysis Toolkit. A python package for manipulating "
            "predicted structures and structural alignments."
        ),
        usage="""python SAT.py <subcommand> [options]""",
    )
    subparsers = parser.add_subparsers(
        title="Subcommands",
        required=True,
    )

    # -------------------------------------------------------------------------------- #
    # Parser for get_domains subcommand
    # -------------------------------------------------------------------------------- #
    parser_get_domains = subparsers.add_parser(
        "get_domains",
        help="""
        Extract domains from a structure using PAE information. Notably, this script is
        designed for colabfold-generated pae json files. json files generated by other
        structural prediction software likely will not work.
        """,
    )
    parser_get_domains.add_argument(
        "-s",
        "--structure_file_path",
        type=str,
        required=True,
        default="",
        help="""
        Path to the input structure in .pdb format.
        """,
    )
    parser_get_domains.add_argument(
        "-p",
        "--pae_path",
        type=str,
        required=True,
        default="",
        help="""
        Path to the pae .json file.
        """,
    )
    parser_get_domains.add_argument(
        "-o",
        "--output_prefix",
        type=str,
        required=True,
        default="",
        help="""
        Prefix of resultant output files. Files will be labled {prefix}_domain-{i}.pdb.
        Note that the domain number will be 1-indexed.
        """,
    )
    parser_get_domains.add_argument(
        "-1",
        "--pae_power",
        type=int,
        required=False,
        default=1,
        help="""
        Default: 1
        Each edge in the graph will be weighted proportional to (1/pae**pae_power)
        """,
    )
    parser_get_domains.add_argument(
        "-2",
        "--pae_cutoff",
        type=int,
        required=False,
        default=5,
        help="""
        Default: 5
        Graph edges will only be created for residue pairs with pae<pae_cutoff. Lowering
        this will make domain identification more stringent by reducing the amount of
        error allowed.
        """,
    )
    parser_get_domains.add_argument(
        "-3",
        "--graph_resolution",
        type=int,
        required=False,
        default=1,
        help="""
        Default: 1
        Regulates how aggressively the clustering algorithm is. Smaller values lead to
        larger clusters. Value should be larger than zero, and values larger than 5 are
        unlikely to be useful.
        """,
    )
    parser_get_domains.set_defaults(func=call_get_domains)

    # -------------------------------------------------------------------------------- #
    # Parser for remove_redundant_domains subcommand
    # -------------------------------------------------------------------------------- #
    parser_remove_redundant_domains = subparsers.add_parser(
        "remove_redundant_domains",
        help=(
            """Remove PDB files that overlap. If a structure has a primary amino acid
            sequence that overlaps another structure in the input, the structure with
            the longer length will be output. If structures are the same length, the
            structure with the highest average pLDDT will be output. If a structure has
            no overlap with any other structures, it will be output."""
        ),
    )
    parser_remove_redundant_domains.add_argument(
        "-i",
        "--input_structure_glob",
        type=str,
        required=True,
        default="",
        help="""
        Glob specifying the structures to be compared. Remember to wrap the glob in
        quotes!
        """,
    )
    parser_remove_redundant_domains.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=True,
        default="",
        help="""
        Path to the output directory in which the filtered files will be saved.
        """,
    )
    parser_remove_redundant_domains.set_defaults(func=call_remove_redundant_domains)

    # -------------------------------------------------------------------------------- #
    # Parser for process_clusters subcommand
    # -------------------------------------------------------------------------------- #
    parser_process_clusters = subparsers.add_parser(
        "process_clusters",
        help=(
            """
            Labels a foldseek all-by-all alignment output with clustering information
            from foldseek cluster. Only alignments that are present in a cluster are
            kept, and only one alignment for every query-target pair. The output file
            is an alignment file with three added columns:
            - cluster_ID: ID of the cluster, starting at 1. A lower number indicates a
                larger cluster
            - cluster_rep: The name of the structure that is the cluster representative
                chosen by foldseek cluster
            - cluster_count: The number of structures in the cluster.
            """
        ),
    )
    parser_process_clusters.add_argument(
        "-a",
        "--alignment_file",
        type=str,
        required=True,
        default="",
        help="""
        Path to the foldseek alignment file.
        """,
    )
    parser_process_clusters.add_argument(
        "-c",
        "--cluster_file",
        type=str,
        required=True,
        default="",
        help="""
        Path to the foldseek cluster tsv file.
        """,
    )
    parser_process_clusters.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        default="",
        help="""
        Path to the output file.
        """,
    )
    parser_process_clusters.add_argument(
        "-f",
        "--alignment_fields",
        type=str,
        required=False,
        default=(
            "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,"
            "evalue,bits,alntmscore"
        ),
        help="""
        A comma-delimited string of the fields in the input foldseek alignment file.
        Make sure to wrap in quotes!

        Default: 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,
        evalue,bits,alntmscore'
        """,
    )
    parser_process_clusters.set_defaults(func=call_process_clusters)

    # -------------------------------------------------------------------------------- #
    # Parser for add_taxonomy_to_alignments subcommand
    # -------------------------------------------------------------------------------- #
    parser_add_taxonomy_to_alignments = subparsers.add_parser(
        "add_taxonomy_to_alignments",
        help=(
            """
            Takes in a foldseek alignment file and adds taxonomy information for the
            query and target. Notably, the query taxonID  is assumed to be
            query_name.stripl('.pdb').split('__')[-1] - e.g. the last value in the
            double-underscore-delimited list. In contrast, the target taxonID is assumed
            to be at that location OR can be present in the taxonid field outputted
            by foldseek.
            """
        ),
    )
    parser_add_taxonomy_to_alignments.add_argument(
        "-a",
        "--alignment_file",
        type=str,
        required=True,
        default="",
        help="""
        Path to the foldseek alignment file.
        """,
    )
    parser_add_taxonomy_to_alignments.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        default="",
        help="""
        Path to the output file.
        """,
    )
    parser_add_taxonomy_to_alignments.add_argument(
        "-f",
        "--alignment_fields",
        type=str,
        required=False,
        default=(
            "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,"
            "evalue,bits,alntmscore,cluster_ID,cluster_rep,cluster_count"
        ),
        help="""
        Fields present in the alignment file.
        Default: query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,
        evalue,bits,alntmscore,cluster_ID,cluster_rep,cluster_count
        """,
    )
    parser_add_taxonomy_to_alignments.add_argument(
        "-q",
        "--query_taxid_location",
        type=int,
        required=False,
        default=1,
        choices=[0, 1],
        help="""
        Default: 1
        Where is the query_taxid located? Options:
        - 0: Indicates that the query_taxid is not present or not desired in the output.
        - 1: Indicates that the query_taxid is present in the query name of the
            alignment as query_name.stripl('.pdb').split('__')[-1]
        """,
    )
    parser_add_taxonomy_to_alignments.add_argument(
        "-t",
        "--target_taxid_location",
        type=int,
        required=False,
        default=1,
        choices=[0, 1, 2],
        help="""
        Default: 1
        Where is the target_taxid located? Options:
        - 0: Indicates that the target_taxid is not present or not desired in the output
        - 1: Indicates that the target_taxid is present in the query name of the
            alignment as target_name.stripl('.pdb').split('__')[-1]
        -2:  Indicates that the target_taxid is present in the taxid field of the
            alignment file.
        """,
    )
    parser_add_taxonomy_to_alignments.add_argument(
        "-T",
        "--taxonomy_levels",
        type=str,
        required=False,
        default="superkingdom,phylum,class,order,family,genus,species",
        help="""
        Default: 1
        Where is the target_taxid located? Options:
        - 0: Indicates that the target_taxid is not present or not desired in the output
        - 1: Indicates that the target_taxid is present in the query name of the
            alignment as target_name.stripl('.pdb').split('__')[-1]
        -2:  Indicates that the target_taxid is present in the taxid field of the
            alignment file.
        """,
    )
    parser_add_taxonomy_to_alignments.set_defaults(func=call_add_taxonomy_to_alignments)

    # Parse the args and call the function associated with the subcommand
    args = parser.parse_args()
    args.func(args)


def call_get_domains(args):
    from scripts.get_domains import get_domains_main

    get_domains_main(args)


def call_remove_redundant_domains(args):
    from scripts.remove_redundant_domains import remove_redundant_domains_main

    remove_redundant_domains_main(args)


def call_process_clusters(args):
    from scripts.process_clusters import process_clusters_main

    args.alignment_fields = args.alignment_fields.split(",")
    process_clusters_main(args)


def call_add_taxonomy_to_alignments(args):
    from scripts.add_taxonomy_to_alignments import add_taxonomy_to_alignments_main

    add_taxonomy_to_alignments_main(args)


if __name__ == "__main__":
    main()
