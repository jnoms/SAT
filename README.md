# SAT
Structural Analysis Toolkit - A python package for manipulation of structural data and structural alignments.  

![Tests](https://github.com/jnoms/SAT/actions/workflows/main.yml/badge.svg)

# Installation
There are two methods to install this package - in a poetry environment, or within a conda environment.

## Installation in a poetry environment
1. Make sure the [Poetry package manager](https://python-poetry.org/) is installed.
1. Clone this repository
2. Enter the SAT directory, and install the package via poetry. This will download all dependencies into a poetry virtual environment.
```
poetry install
```
3. Now, you can run SAT in the following way:
```
poetry run sat.py <subcommand>
```

## Installation in a conda environment
1. Make sure the [conda package manager](https://docs.conda.io/en/latest/miniconda.html) is installed.
1. Clone this repository
2. Create a conda environment that contains the poetry package manager
```
conda create --name SAT -c conda-forge poetry
```
3. Activate environment. Enter the SAT directory and download the dependencies using poetry. Dependencies will be downloaded specifically into that conda environment.
```
conda activate SAT
poetry install
```
4. The SAT conda environment will now contain all dependencies. You can run SAT in the following way:
```
conda activate SAT #if conda enviornment is not active  

sat.py <subcommand>
```

## Test installation
Navigate to the sat directory and enter `pytest` (if using a conda environment) or `poetry run pytest` if using a poetry environment. This will make sure that all tests pass and sat is properly installed.  

## Note on ETE3 
When you run the tests or the first time you run any taxonomy-related script, ete3 will download a taxonomy database to ~/.etetoolkit/. **This database is 560MB** as of October 2022. Capability to specify the database download location is on the to-do list, but in the interim you can make a symlink from ~/.etetoolkit/ to wherever you want the database to reside (see below). This is particularly important if taxonomy related queries are going slowly, as that probably means your home directory has slow IO so you should symlink to somewhere with faster IO.    
```ln -s /desired/ete/database/location ~/.etetoolkit```  


# List of subcommands

## Structure-focused
`sat.py struc_get_domains` - Uses PAE information to extract well-folded domains from an input structure.  
`sat.py struc_remove_redundant` - Removes domains that have strongly overlapping primary amino-acid sequences.   
`sat.py struc_find_motif` - Checks if there is a motif in a structure or sequence input.  
`sat.py struc_to_seq` - Prints the primary amino acid sequence of a structure to the screen or appends to a specified file in fasta format.  
`sat.py struc_to_plddt` - Prints the average pLDDT of a structure to the screen or appends to a specified file.  
`sat.py struc_rebase` - Rebases an input structure such that the first residue is residue #1, and all subsequent residues are sequential (e.g. removes numeric gaps present in discontinuous domains).  
`sat.py struc_qc` - Get information on the fraction of residues that are at least a specified pLDDT - this can be good for filtration.  
`sat.py struc_disorder` - Get information on the number of residues in an input structure that are considered disordered and ordered.  

## Alignment-focused
`sat.py aln_add_clusters` - Adds foldseek clustering information to the foldseek tabular alignment file.  
`sat.py aln_add_taxonomy` - Adds specified taxonomic levels for the query and/or target of foldseek alignments.  
`sat.py aln_taxa_counts` - Returns counts at desired taxonomic levels within each foldseek cluster.  
`sat.py aln_query_uniprot` - Lets you look up alphafold or uniprot IDs using the Uniprot REST API, and get the geneName and fullName (an informative protein name) for each.  
`sat.py aln_add_uniprot` - After retreiving the uniprot unformation using aln_query_uniprot, adds the information as columns to the alignment file.  
`sat.py aln_filter` - This filters for alignments below/above a specified value in a specified column, and can also filter to keep a maximum number of queries per alignment.  
`sat.py aln_merge` - This merges two alignment files.  
`sat.py aln_cluster` - This lets you do connected-component clustering, similar to foldseek/mmseqs cluster mode 1, but lets you have more control on filtering the alignments prior.  
`sat.py aln_expand_clusters` - When you use the cluster representatives from one clustering (often sequence-based) to cluster using another clustering (often structure based), will merge them.  
`sat.py aln_merge_clusters` - This takes in a cluster file and an alignment file of alignments between cluster representatives, and merges clusters whose representatives align.  
`sat.py aln_parse_dali` - This parses a Dalilite alignment output into a tab-delimited format. It can also filter based on various alignment statistics.  
`sat.py aln_generate_superclusters` - This takes in an all-by-all foldseek alignment and foldseek-reported clusters and identifies clusters that are linked (e.g. - some specified number of members of each cluster align to members of the other cluster). These clusters are then merged into a supercluster.  
`sat.py aln_ecod_purity` - This takes in a cluster file and the alignments between the cluster members and the ECOD HMM database and counts, for each ECOD level, the number of members per cluster.

## Sequence-focused  
`sat.py seq_chunk` - Splits a fasta file into overlapping or non-overlapping chunks.

# SAT struc_get_domains
Extract separate domain structures from a predicted structure.  
This uses the PAE information to cluster residues that likely fall into linear domains. Notably, the script is currently only configured to process colabfold-generated PAE files. 
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_get_domains -h`](.github/img/struc_get_domains.png)  


# SAT struc_remove_redundant
Given a glob specifying multiple structure (often domains), will remove structures that have an overlapping pirmary amino acid sequence.   
Priority is given to the longer structure or, if the sequences are the same length, the structure with the highest pLDDT
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_remove_redundant -h`](.github/img/struc_remove_redundant.png)  

# SAT struc_find_motif
Given a motif of structure [OPTIONS]xxx[OPTIONS]xx, where x indicates any amino acid and [] indicate any of the amino acids present within the brackets, this returns the match and position start/end of the motif present in the input sequence.  

The input can be a structure file, a fasta, or just a sequence. The output is tab-delimited and printed to the screen, with the columns  
- match  
- start   
- end   
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_find_motif -h`](.github/img/struc_find_motif.png) 

# SAT struc_qc
Given a structure, determines the percentage of residues that have at least the specified pLDDT. The output is returned to STDOUT!! It is tab-delimited and has the following columns:   
- structure_file (the basename of the file, including suffix)  
- number of residues  
- number of residues that pass the pLDDT threshold  
- proportion of residues that pass the pLDDT threshold (this will be a decimal between 0 and 1)
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_qc -h`](.github/img/struc_qc.png) 

# SAT struc_disorder
This takes an input structure and calculates the number of residues that are considered ordered, disordered, or intermediate. A residue is considered ordered if it is in a stretch of at least n_sequential residues that have a pLDDT of >= order_cutoff. A residue is considered disordered if it is in a stretech of at least n_sequential residues <= disorder_cutoff.  

This returns an output file with the following columns:  
- basename of the input structure  
- number of ordered residues  
- number of disordered residues   
- number of intermediate residues (neither ordered or disordered)  
- total number of residues  
- there_is_a_domain: yes or no. This checks that there is at least one stretech of continuous residues that have ordered pLDDTs. The required stretch size is args.check_for_domain_len  
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_disorder -h`](.github/img/struc_disorder.png) 

# SAT aln_add_taxonomy
Adds taxonomy information to a foldseek alignment  
- Taxonomy for each query can be built in to the query_name as the final element of the double-underscore-delimited list.   
- Taxonomy for each target can either be built in to the target name in a similar manner, or present as the taxid field in the foldseek output.  
- This script adds the taxonomic names of the query and/or target at the specified levels.   
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_add_taxonomy -h`](.github/img/aln_add_taxonomy.png)  

# SAT struc_to_seq
Simple subcommand to produce the amino-acid sequence from a structure file.  

Can append the sequence to an outfile if provided, or will print to screen.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_to_seq -h`](.github/img/struc_to_seq.png)  


# SAT struc_to_plddt
Simple subcommand that returns the average plddt of the input structure file. If --out_file is not specified, the average plddt is simply printed to the screen. If --out_file is specified, the output file will be APPENDED to with the following: [basename input structure_file]\\t[plddt]\\n
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_to_plddt -h`](.github/img/struc_to_plddt.png)  


# SAT struc_rebase
Simple subcommand that renumbers all residues in a structure such that the first residue is #1 and all residues are sequential (e.g. it takes out numeric gaps in residue numbers).
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_rebase -h`](.github/img/struc_rebase.png)  


# SAT seq_chunk
Splits entries into a fasta into overlapping or non-overlapping chunks. This is helpful when you want to split up sequences that are too long to effectively use for structure prediction. This subcommand is able to generate overlapping sequences. 
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py seq_chunk -h`](.github/img/seq_chunk.png)  

# SAT struc_download
This subcommand takes in a file of uniprot IDs and downloads the AF2 database pdb and pae files to the indicated directory. Furthermore, if any additional information is present in the tabular infile it will be appended to the output files - this is a good way to lable the files with information like taxonomyID, etc.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_download -h`](.github/img/struc_download.png)  


# SAT aln_taxa_counts
This takes in a cluster file (required columns are cluster_ID, cluster_rep, cluster_member, and cluster_count) and tallies up the taxons for each cluster. It makes a tidy file for each cluster where, for every taxon at every level, it specifies the count. The cluster file is assumed to be generated from an all-by-all alignment, perhaps with some additional merging steps. If you are also interested in adding taxonomy count information for the targets of a search of the cluster members against a separate database, you can enter an alignment file to this script. In the event an alignment file is provided, taxonIDs from the TARGET will be added to the cluster_ID of the QUERY.  
            
The output file has the following columns:  
cluster_ID, cluster_rep, cluster_count, superkingdom, level, taxon, count.  

<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_taxa_counts -h`](.github/img/aln_taxa_counts.png)  

# SAT aln_query_uniprot
This script takes alphafold IDs (or raw uniprot IDs) and uses the Uniprot REST API to get information on the geneName and fullName (an informative name of the protein) for each ID. You can specify in which column of the infile the IDs live.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_query_uniprot -h`](.github/img/aln_query_uniprot.png)  


# SAT aln_add_uniprot
This script adds the uniprot information garther from aln_query_uniprot to a foldseek alignment file.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_add_uniprot -h`](.github/img/aln_add_uniprot.png)  


# SAT aln_filter
This subcommand filters a foldseek alignment file to keep only those alignments with a value below/above the specified value in a field (alntmscore is a common one). It also only outputs a maximum of N alignments for each query.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_filter -h`](.github/img/aln_filter.png)  

# SAT aln_parse_dali
This subcommand reads in a DALI alignment output file and formats it as a tab-delimited file. This script will written to the specified output file. There is also functionality to filter the alignments by zscore, alnlen, coverage, or rmsd.  
There are two main inputs:  
1) alignment_file: This is the DALI alignment file. Notably, the first output field MUST BE the 'summary' and the second output field MUST BE 'equivalences'
2) structure_key: DALI only processes files that have a 4-digit identifier. The structure key must be of format structure[delimiter]identifier, and lets you convert the identifiers back to the actual structure name. Note that the structure_key identifiers should not have the DALI segment (e.g. A, B, C...) at the end - this will be taken care of.  

The qlen field of the output is dependent on their being a self alignment in the alignment file, as then the qlen=tlen. If not present, qlen will be listed as 0.  
            
Note also the coverage is determined by alnlen/max(qlen, tlen)  

The output file is a .m8 file (e.g. tab delimited) and has the following columns: query, target, query_id, target_id, alnlen, qlen, tlen, cov, pident, rmsd, z
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_parse_dali -h`](.github/img/aln_parse_dali.png)  


# SAT aln_merge
This subcommand is used to merge two foldseek alignment files.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_merge -h`](.github/img/aln_merge.png)  

# SAT aln_cluster
This subcommand generates clusters from an input alignment file, where every query-target pair will be put into the same cluster.  

This subcommand basically does what mmseqs/foldseek cluster mode 1 does (e.g. connected-compontent clustering). Here, any two members that are aligned will end up in the same cluster. Because of this strong clustering, the alignment file should be strinctly filtered to only keep those alignments with high coverage and high confidence (e.g. high TMscore from foldseek or high z score from DALI).  

The output file is essentially a foldseek/mmseqs cluster file with two columns: cluster_rep, cluster_member.  

The optional --all_inputs switch can be used to provide information for all members that initially were input to the alignment. If provided, the output cluster file will include those members that aren't present in the alignment file as a cluster with only one member (themselves). This is very useful because the alignment file should be strictly filtered prior to using this script, so many of the items inputted to foldseek or mmseqs won't be present in the alignment file.  
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_cluster -h`](.github/img/aln_cluster.png)  

# SAT aln_expand_clusters
This subcommand is used to merge cluster files from a teired pair of clustering, often first a sequence clustering and then a structure
clustering. It it assumed that the 'subcluster file' representatives were the members subsequently clustered in the 'cluster file'. This script simply adds the members of each subcluster to the parent cluster in the cluster
file.  

Ideally, there will not be a cluster_ID column in the input cluster files. If there is, those cluster_IDs will be used. Otherwise, cluster_IDs will be generated which rank clusters by total members.  

The output file has the following columns:  
- cluster_ID  
- cluster_rep  
- subcluster_rep  
- cluster_member  
- cluster_count  

<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_expand_clusters -h`](.github/img/aln_expand_clusters.png)  

# SAT aln_merge_clusters
This subcommand takes in a cluster file and alignments between the REPRESENTATIVES of the clusters, and merges clusters whose representatives align together.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_merge_clusters -h`](.github/img/aln_merge_clusters.png)  

# SAT aln_ecod_purity
This subcommand takes in a cluster file and an alignment file of those same members aligned (using an HMM approach) to the ECOD HMM database.  It takes in an ECOD information file that connects each ECOD accession to it's classification at various annotation levels. This script returns a tidy-format output file with, for each cluster, the counts of members with alignments against each ECOD entry.  
The output columns are as follows:  
- cluster_ID  
- cluster_rep  
- level  
- value  
- count  

Note that this assumes that each member only has ONE alignment - e.g. the best ECOD alignment.  
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_ecod_purity -h`](.github/img/aln_ecod_purity.png)  


# Planned improvements
struc_get_domains
- Add functionality to parse PAE json files from additional sources

ete3  
- Add ability to specify where the ete3 taxonomy database is downloaded.
