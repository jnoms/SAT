#!/usr/bin/env python

import argparse


def main():

    # Top-level parser
    parser = argparse.ArgumentParser(
        description=(
            "SAT - Structural Analysis Toolkit. A python package for manipulating "
            "predicted structures and structural alignments."
        ),
        usage="""sat.py <subcommand>""",
    )
    subparsers = parser.add_subparsers(
        title="Subcommands",
        required=True,
    )

    # -------------------------------------------------------------------------------- #
    # Parser for struc_get_domains subcommand
    # -------------------------------------------------------------------------------- #
    parser_struc_get_domains = subparsers.add_parser(
        "struc_get_domains",
        help="""
        Extract domains from a structure using PAE information. Notably, this script is
        designed for colabfold-generated pae json files. json files generated by other
        structural prediction software likely will not work.
        """,
    )
    parser_struc_get_domains.add_argument(
        "-s",
        "--structure_file_path",
        type=str,
        required=True,
        default="",
        help="""
        Path to the input structure in .pdb format.
        """,
    )
    parser_struc_get_domains.add_argument(
        "-p",
        "--pae_path",
        type=str,
        required=True,
        default="",
        help="""
        Path to the pae .json file.
        """,
    )
    parser_struc_get_domains.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=True,
        default="",
        help="""
        Directory of the output domains. Files will be labled
        {output_dir}/{basename of structure}_domain-{i}.pdb.
        Note that the domain number will be 1-indexed.
        """,
    )
    parser_struc_get_domains.add_argument(
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
    parser_struc_get_domains.add_argument(
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
    parser_struc_get_domains.add_argument(
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
    parser_struc_get_domains.add_argument(
        "-l",
        "--min_domain_length",
        type=int,
        required=False,
        default=50,
        help="""
        Default: 50
        A domain must be at least this long to be called.
        """,
    )
    parser_struc_get_domains.add_argument(
        "-P",
        "--min_domain_plddt",
        type=int,
        required=False,
        default=60,
        help="""
        Default: 60
        A domain must have an average pLDDT of at least this to be called.
        """,
    )
    parser_struc_get_domains.add_argument(
        "-n",
        "--smooth_n",
        type=int,
        required=False,
        default=0,
        help="""
        Default: 0
        If set to non-zero value, will smooth out the PAE matrix. If a region of high
        PAE (>5) is less than n residues long and is surrounded by a region of low PAE,
        it will be overridden with a low PAE. Furthermore, all of the low-PAE regions
        will be overridden with a uniformely low PAE of 1.
        """,
    )
    parser_struc_get_domains.add_argument(
        "-r",
        "--plddt_report",
        type=str,
        required=False,
        default="",
        help="""
        If specified, for each domain will write the average plddt to this file. 
        This file will be APPENDED to. The columns are domain_name, average_plddt.
        """,
    )
    parser_struc_get_domains.add_argument(
        "-R",
        "--pae_report",
        type=str,
        required=False,
        default="",
        help="""
        If specified, for each domain will write the average PAE to this file. 
        This file will be APPENDED to. The columns are domain_name, average_PAE.
        """,
    )
    parser_struc_get_domains.set_defaults(func=call_struc_get_domains)

    # -------------------------------------------------------------------------------- #
    # Parser for struc_remove_redundant subcommand
    # -------------------------------------------------------------------------------- #
    parser_struc_remove_redundant = subparsers.add_parser(
        "struc_remove_redundant",
        help=(
            """Remove PDB files that overlap. If a structure has a primary amino acid
            sequence that overlaps another structure in the input, the structure with
            the longer length will be output. If structures are the same length, the
            structure with the highest average pLDDT will be output. If a structure has
            no overlap with any other structures, it will be output."""
        ),
    )
    parser_struc_remove_redundant.add_argument(
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
    parser_struc_remove_redundant.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=True,
        default="",
        help="""
        Path to the output directory in which the filtered files will be saved.
        """,
    )
    parser_struc_remove_redundant.set_defaults(func=call_struc_remove_redundant)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_add_clusters subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_add_clusters = subparsers.add_parser(
        "aln_add_clusters",
        help=(
            """
            This subcommand incorporates clustering information from foldseek cluster
            into the foldseek alignment tabular output file. Notably, there are two
            possible outputs from this script:
            1) An output alignment file where, for each cluster, only the 'top' query
                (aka the one with the highest number of alignments or, if there is a
                tie, the one with the highest average TMscore).
            2) An output alignment file containing all non-redundant alingments. Here,
                all self-self alignments are removed. Furthermore, because this file
                arose from all-by-all alignments, there will be two alignments for each
                pair of items. Only one will be present in the output file. Finally,
                Plese NOTE!! that only alignments with both the query and target in the
                SAME CLUSTER are kept - cross-cluster alignments are removed.

            This script adds the fields "cluster_ID", "cluster_count", and
            "top_query" to the output files.
            """
        ),
    )
    parser_aln_add_clusters.add_argument(
        "-a",
        "--alignment_file",
        type=str,
        required=True,
        default="",
        help="""
        Path to the foldseek alignment file.
        """,
    )
    parser_aln_add_clusters.add_argument(
        "-c",
        "--cluster_file",
        type=str,
        required=True,
        default="",
        help="""
        Path to the foldseek cluster tsv file.
        """,
    )
    parser_aln_add_clusters.add_argument(
        "-1",
        "--top_query_per_cluster_out",
        type=str,
        required=False,
        default="",
        help="""
        Path to the output file.
        """,
    )
    parser_aln_add_clusters.add_argument(
        "-2",
        "--all_nonredundant_out",
        type=str,
        required=False,
        default="",
        help="""
        Path to the output file.
        """,
    )
    parser_aln_add_clusters.add_argument(
        "-f",
        "--alignment_fields",
        type=str,
        required=False,
        default="",
        help="""
        A comma-delimited string of the fields in the input foldseek alignment file.
        Make sure to wrap in quotes!

        Default: ''
        """,
    )
    parser_aln_add_clusters.set_defaults(func=call_aln_add_clusters)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_add_taxonomy subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_add_taxonomy = subparsers.add_parser(
        "aln_add_taxonomy",
        help=(
            """
            Takes in a foldseek alignment file and adds taxonomy information for the
            query and target. Notably, the query taxonID  is assumed to be
            query_name.rstrip('.pdb').split('__')[-1] - e.g. the last value in the
            double-underscore-delimited list. In contrast, the target taxonID is assumed
            to be at that location OR can be present in the taxonid field outputted
            by foldseek.
            """
        ),
    )
    parser_aln_add_taxonomy.add_argument(
        "-a",
        "--alignment_file",
        type=str,
        required=True,
        default="",
        help="""
        Path to the foldseek alignment file.
        """,
    )
    parser_aln_add_taxonomy.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        default="",
        help="""
        Path to the output file.
        """,
    )
    parser_aln_add_taxonomy.add_argument(
        "-f",
        "--alignment_fields",
        type=str,
        required=False,
        default="",
        help="""
        Fields present in the alignment file.
        Default: ''
        """,
    )
    parser_aln_add_taxonomy.add_argument(
        "-d",
        "--taxonID_finder_delimiter",
        type=str,
        required=False,
        default="__",
        help="""
        This script will try to find the taxonID in the query or target strings of each
        alignment (although if there is a taxid column in the alignment_object, it will
        use that taxonID for the target taxonID). To parse for the taxonID, this is
        the delimiter.
        """,
    )
    parser_aln_add_taxonomy.add_argument(
        "-p",
        "--taxonID_finder_pos",
        type=int,
        required=False,
        default=-1,
        help="""
        This is where in the delimited string the taxonID is located.
        """,
    )
    parser_aln_add_taxonomy.add_argument(
        "-T",
        "--taxonomy_levels",
        type=str,
        required=False,
        default="superkingdom,phylum,class,order,family,genus,species",
        help="""
        These are the taxonomic levels to include in the output file.
        """,
    )
    parser_aln_add_taxonomy.set_defaults(func=call_aln_add_taxonomy)

    # -------------------------------------------------------------------------------- #
    # Parser for struc_to_seq subcommand
    # -------------------------------------------------------------------------------- #
    parser_struc_to_seq = subparsers.add_parser(
        "struc_to_seq",
        help=(
            """
            Simple subcommand that returns the amino-acid sequence of a specified
            structure (in pdb format). The AA sequence will be APPENDED to the outfile
            if specified, or printed to the screen if -o --out_file is not specified.
            """
        ),
    )
    parser_struc_to_seq.add_argument(
        "-s",
        "--structure_file",
        type=str,
        required=True,
        help="""
        Path to the structure file in pdb format.
        """,
    )
    parser_struc_to_seq.add_argument(
        "-o",
        "--out_file",
        type=str,
        required=False,
        default="",
        help="""
        Path to a file the sequence will be appended to. If left blank, sequence will be
        printed to the screen.
        """,
    )
    parser_struc_to_seq.add_argument(
        "-H",
        "--header",
        type=str,
        required=False,
        default="",
        help="""
        Header of the entry if writing to a fasta. Only required if -o --out_file is
        specified.
        """,
    )
    parser_struc_to_seq.set_defaults(func=call_struc_to_seq)

    # -------------------------------------------------------------------------------- #
    # Parser for struc_rebase subcommand
    # -------------------------------------------------------------------------------- #
    parser_struc_rebase = subparsers.add_parser(
        "struc_rebase",
        help=(
            """
            Simple subcommand that renumbers all residues in a structure such that
            the first residue is #1 and all residues are sequential (e.g. it takes out
            numeric gaps in residue numbers).
            """
        ),
    )
    parser_struc_rebase.add_argument(
        "-s",
        "--structure_file",
        type=str,
        required=True,
        help="""
        Path to the structure file in pdb format.
        """,
    )
    parser_struc_rebase.add_argument(
        "-o",
        "--out_file",
        type=str,
        required=True,
        help="""
        Path to the output structure file.
        """,
    )
    parser_struc_rebase.set_defaults(func=call_struc_rebase_main)

    # -------------------------------------------------------------------------------- #
    # Parser for seq_chunk subcommand
    # -------------------------------------------------------------------------------- #
    parser_seq_chunk = subparsers.add_parser(
        "seq_chunk",
        help=(
            """
            Tool to split sequences in an input fasta into overlapping or not
            overlapping chunks. Can write all resultant sequences to a single output
            file, or to separate output files (one per header). The header of
            chunks are >PART{N}_{header}.

            For example, if a sequence is 2300AA long but you desire sequences
            of 1000 max, this script can either generate the following:
            1) 1-1000, 1001-2000, 2001-2300. (if -v is NOT specified)
            2) 1-1000, 501-1500, 1001-2000, 1501-2300, 2001-2300 (if -v is specified)
            """
        ),
    )
    parser_seq_chunk.add_argument(
        "-i",
        "--in_fasta",
        type=str,
        required=True,
        help="""
        Path to the input fasta.
        """,
    )
    parser_seq_chunk.add_argument(
        "-o",
        "--out_fasta",
        type=str,
        required=True,
        help="""
        Path to the output fasta. If -n is specified, will put every fasta entry
        into a separate file and the entry in the -o switch will specify the
        base path - should be a directory in this case.
        """,
    )
    parser_seq_chunk.add_argument(
        "-m",
        "--max_seq_length",
        type=int,
        required=True,
        help="""
        The maximum sequence length.
        """,
    )
    parser_seq_chunk.add_argument(
        "-s",
        "--minimum_sequence_output_size",
        type=int,
        required=False,
        default=50,
        help="""
        Will not output any sequences below this size. This is because colabfold
        fails when sequences are very small. [50]
        """,
    )
    parser_seq_chunk.add_argument(
        "-v",
        "--overlapping_chunks",
        type=arg_str2bool,
        required=False,
        default=False,
        nargs="?",
        const=True,
        help="""
        If specified, each chunk will overlap by max_seq_length/2. If not
        specified, chunks will be just
        1:max_seq_length, max_seq_length+1:max_seq_length*2, etc
        """,
    )
    parser_seq_chunk.add_argument(
        "-n",
        "--individual",
        type=arg_str2bool,
        required=False,
        default=False,
        nargs="?",
        const=True,
        help="""
        If specified, each fasta entry will be passed to a separate file
        named by its header. Furthermore, it will assume that args.out_fasta is
        the base path.
        """,
    )
    parser_seq_chunk.set_defaults(func=call_seq_chunk_main)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_taxa_counts subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_taxa_counts = subparsers.add_parser(
        "aln_taxa_counts",
        help=(
            """
            This takes in a processed alignment file (typically generated from a
            foldseek alignment that was then processed through aln_add_clusters and
            add_taxonomy_to_alignment) and returns, for each cluster, the number of
            taxa at each taxonomic level and their names. The output file has the
            following columns:
            cluster_ID, cluster_rep, level, superkingdom, taxon, count.
            """
        ),
    )
    parser_aln_taxa_counts.add_argument(
        "-a",
        "--alignment_file",
        type=str,
        required=True,
        help="""
        Path to the alignment file.
        """,
    )
    parser_aln_taxa_counts.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        help="""
        Path to the output count file. This file will have the following columns:
        cluster_ID, cluster_rep, level, superkingdom, taxon, count.
        """,
    )
    parser_aln_taxa_counts.add_argument(
        "-t",
        "--taxonomy_levels",
        type=str,
        required=False,
        default="superkingdom,phylum,class,order,family,genus,species",
        help="""
        Taxonomy levels to count and output
        """,
    )
    parser_aln_taxa_counts.set_defaults(func=call_parser_aln_taxa_counts_main)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_query_uniprot subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_query_uniprot = subparsers.add_parser(
        "aln_query_uniprot",
        help=(
            """
            This subcommand takes in a file that contains uniprotIDs (or alphafold IDs,
            even if they are formatted like e.g. AF-K0EZQ3-F1-model_v2.pdb.gz). You just
            need to specify the 0-indexed column position of the ID.

            This script uses the uniprot REST API to download information for each
            uniprotID. It offers the funcionality to save a cache containing the raw
            download information to prevent unncessary downloading.

            Besides making or updating a cache file, the main purpose of this script
            is to output a flat, tab-delimited file with the columns uniprotID,
            geneName, and fullName (fullName is a descriptive protein name).
            """
        ),
    )
    parser_aln_query_uniprot.add_argument(
        "-i",
        "--infile",
        type=str,
        required=True,
        help="""
        Path to the input file that contains the uniprot IDs. If there are multiple
        columns, this must be tab-delimited.
        """,
    )
    parser_aln_query_uniprot.add_argument(
        "-o",
        "--uniprot_lookup_output",
        type=str,
        required=True,
        help="""
        This is the main output file. Will be a tab-delimited file with the columns
        uniprotID, geneName, fullName
        """,
    )
    parser_aln_query_uniprot.add_argument(
        "-c",
        "--infile_col",
        type=int,
        required=False,
        default=1,
        help="""
        Default: 1
        This is the 0-indexed column holding the uniprot of alphafold IDs.
        """,
    )
    parser_aln_query_uniprot.add_argument(
        "-u",
        "--uniprot_cache",
        type=str,
        required=False,
        default="",
        help="""
        A pkl file containing a cache from former uniprot downloads. This is optional.
        If specified, this file will be read in and will be updated with this script's
        downloads.
        """,
    )
    parser_aln_query_uniprot.set_defaults(func=call_aln_query_uniprot_main)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_add_uniprot subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_add_uniprot = subparsers.add_parser(
        "aln_add_uniprot",
        help=(
            """
            This subcommand takes in uniprot information generated by aln_query_uniprot
            and adds it to a foldseek alignment.
            """
        ),
    )
    parser_aln_add_uniprot.add_argument(
        "-a",
        "--alignment_file",
        type=str,
        required=True,
        help="""
        Path to the alignment file. It is OK if the first row is the header, as long as
        the first column is 'query'.
        """,
    )
    parser_aln_add_uniprot.add_argument(
        "-u",
        "--uniprot_information",
        type=str,
        required=True,
        help="""
        Path to the uniprot information file generated by aln_query_uniprot.
        """,
    )
    parser_aln_add_uniprot.add_argument(
        "-o",
        "--outfile",
        type=str,
        required=True,
        help="""
        Path to the output file.
        """,
    )
    parser_aln_add_uniprot.add_argument(
        "-f",
        "--alignment_fields",
        type=str,
        required=False,
        default="",
        help="""
        Alignment fields present in the input file. Note that if there are labeled
        column headers in the input alignment file, they can be automatically parsed
        and you can leave the alignment_fields command empty.
        """,
    )
    parser_aln_add_uniprot.set_defaults(func=call_aln_add_uniprot_main)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_filter subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_filter = subparsers.add_parser(
        "aln_filter",
        help=(
            """
            This subcommand filters a foldseek alignment file to keep only those
            alignments with a value below/above the specified value in a field
            (alntmscore is a common one). It also only outputs a maximum of N alignments
            for each query.
            """
        ),
    )
    parser_aln_filter.add_argument(
        "-a",
        "--alignment_file",
        type=str,
        required=True,
        help="""
        Path to the alignment file. It is OK if the first row is the header, as long as
        the first column is 'query'.
        """,
    )
    parser_aln_filter.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        help="""
        Path to the output alignment file.
        """,
    )
    parser_aln_filter.add_argument(
        "-f",
        "--alignment_fields",
        type=str,
        required=False,
        default="",
        help="""
        Comma-delimited string of alignment fields. Leave blank if there is a header in
        the input file.
        """,
    )
    parser_aln_filter.add_argument(
        "-x",
        "--filter_field",
        type=str,
        required=False,
        default="alntmscore",
        help="""
        Default: 'alntmscore'
        String indicating the alignment field which contains the value to be filtered/
        sorted with.
        """,
    )
    parser_aln_filter.add_argument(
        "-N",
        "--N",
        type=int,
        required=False,
        default=10,
        help="""
        Default: 10
        This is the maximum number of alignments to output for each query. Set to 0 if
        you want to return all alignments.
        """,
    )
    parser_aln_filter.add_argument(
        "-m",
        "--min_val_filter_field",
        type=float,
        required=False,
        default=0.4,
        help="""
        Default: 0.4
        An alignment must have at least this minimum value in the filter_field for the
        alignment to be output.
        """,
    )
    parser_aln_filter.add_argument(
        "-M",
        "--max_val_filter_field",
        type=float,
        required=False,
        default=1,
        help="""
        Default: 1
        An alignment must have this value or less in the filter_field for the
        alignment to be output.
        """,
    )
    parser_aln_filter.set_defaults(func=call_aln_filter_main)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_merge subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_merge = subparsers.add_parser(
        "aln_merge",
        help=(
            """
            This subcommand merges to foldseek alignments. Columns that are present
            in one but not present in the other will be carried over.
            """
        ),
    )
    parser_aln_merge.add_argument(
        "-1",
        "--aln1",
        type=str,
        required=True,
        help="""
        Path the first alignment file.
        """,
    )
    parser_aln_merge.add_argument(
        "-2",
        "--aln2",
        type=str,
        required=True,
        help="""
        Path the second alignment file.
        """,
    )
    parser_aln_merge.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="""
        Path the output file.
        """,
    )
    parser_aln_merge.add_argument(
        "-f",
        "--aln1_fields",
        type=str,
        required=False,
        default="",
        help="""
        Alignment fields present in the first alignment file. Can leave blank if the
        alignment file has headers.
        """,
    )
    parser_aln_merge.add_argument(
        "-F",
        "--aln2_fields",
        type=str,
        required=False,
        default="",
        help="""
        Alignment fields present in the second alignment file. Can leave blank if the
        alignment file has headers.
        """,
    )
    parser_aln_merge.add_argument(
        "-s",
        "--aln1_source",
        type=str,
        required=False,
        default="",
        help="""
        String indicating the source of alignment 1. This value will be added to each 
        alignment that originated from alignment 1. You must use either neither 
        aln1_source or aln2_source or both. Leave both blank if not desired. 
        """,
    )
    parser_aln_merge.add_argument(
        "-S",
        "--aln2_source",
        type=str,
        required=False,
        default="",
        help="""
        String indicating the source of alignment 2. This value will be added to each 
        alignment that originated from alignment 2. You must use either neither 
        aln1_source or aln2_source or both. Leave both blank if not desired. 
        """,
    )

    parser_aln_merge.set_defaults(func=call_aln_merge_main)

    # Parse the args and call the function associated with the subcommand
    args = parser.parse_args()
    args.func(args)


def arg_str2bool(v):
    """
    For use as an argparse argument type. Makes it easy to use boolean flags.
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


def call_struc_get_domains(args):
    from scripts.struc_get_domains import struc_get_domains_main

    struc_get_domains_main(args)


def call_struc_remove_redundant(args):
    from scripts.struc_remove_redundant import struc_remove_redundant_main

    struc_remove_redundant_main(args)


def call_aln_add_clusters(args):
    from scripts.aln_add_clusters import aln_add_clusters_main

    aln_add_clusters_main(args)


def call_aln_add_taxonomy(args):
    from scripts.aln_add_taxonomy import aln_add_taxonomy_main

    aln_add_taxonomy_main(args)


def call_struc_to_seq(args):
    from scripts.struc_to_seq import struc_to_seq_main

    struc_to_seq_main(args)


def call_struc_rebase_main(args):
    from scripts.struc_rebase import struc_rebase_main

    struc_rebase_main(args)


def call_seq_chunk_main(args):
    from scripts.seq_chunk import seq_chunk_main

    seq_chunk_main(args)


def call_parser_aln_taxa_counts_main(args):
    from scripts.aln_taxa_counts import aln_taxa_counts_main

    aln_taxa_counts_main(args)


def call_aln_query_uniprot_main(args):
    from scripts.aln_query_uniprot import aln_query_uniprot_main

    aln_query_uniprot_main(args)


def call_aln_add_uniprot_main(args):
    from scripts.aln_add_uniprot import (
        aln_add_uniprot_main,
    )

    aln_add_uniprot_main(args)


def call_aln_filter_main(args):
    from scripts.aln_filter import aln_filter_main

    aln_filter_main(args)


def call_aln_merge_main(args):
    from scripts.aln_merge import aln_merge_main

    aln_merge_main(args)


#
#
#
if __name__ == "__main__":
    main()
