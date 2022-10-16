# SAT
Structural Analysis Toolkit - A python package for manipulation of structural data and identification of shared structural features.  

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
`sat struc_get_domains` - Uses PAE information to extract well-folded domains from an input structure.  
`sat struc_remove_redundant` - Removes domains that have strongly overlapping primary amino-acid sequences.   
`sat struc_to_seq` - Prints the primary amino acid sequence of a structure to the screen or appends to a specified file in fasta format.  
`sat struc_rebase` - Rebases an input structure such that the first residue is residue #1, and all subsequent residues are sequential (e.g. removes numeric gaps present in discontinuous domains).  

## Foldseek-focused
`sat aln_clusters` - Adds foldseek clustering information to the foldseek tabular alignment file.  
`sat aln_add_taxonomy` - Adds specified taxonomic levels for the query and/or target of foldseek alignments.  
`sat aln_taxa_counts` - Returns counts at desired taxonomic levels within each foldseek cluster.  
`sat aln_query_uniprot` - Lets you look up alphafold or uniprot IDs using the Uniprot REST API, and get the geneName and fullName (an informative protein name) for each.  
`sat aln_add_uniprot` - After retreiving the uniprot unformation using aln_query_uniprot, adds the information as columns to the alignment file.  
`sat aln_filter` - This filters for alignments below/above a specified value in a specified column, and can also filter to keep a maximum number of queries per alignment.  

## Sequence-focused  
`sat seq_chunk` - Splits a fasta file into overlapping or non-overlapping chunks.

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


# SAT aln_clusters
Adds cluster information from foldseek cluster into the foldseek alignment information.  
Currently, foldseek cluster can generate clusters, and the alignment can output alignments, but it will be helpful to annotate each alginment with their cluster. Furthermore, all-by-all searchs  result in a redundant query-target and target-query alignment for each entry. These redundant alignments are removed.  
This script adds the following fields to the input alignment file:  
1) cluster_ID: ID of the cluster, starting at 1. A lower number indicates a larger cluster  
2) cluster_rep: The name of the structure that is the cluster representative chosen by foldseek  
3) cluster_count: The number of structures in the cluster.  
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_clusters -h`](.github/img/aln_clusters.png)  

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

# SAT struc_rebase
Simple subcommand that renumbers all residues in a structure such that the first residue is #1 and all residues are sequential (e.g. it takes out numeric gaps in residue numbers).
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_rebase -h`](.github/img/struc_rebase.png)  


# SAT seq_chunk
Splits entries into a fasta into overlapping or non-overlapping chunks. This is helpful when you want to split up sequences that are too long to effectively use for structure prediction. This subcommand is able to generate overlapping sequences. 
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py seq_chunk -h`](.github/img/seq_chunk.png)  


# SAT aln_taxa_counts
This takes in a processed alignment file (typically generated from a foldseek alignment that was then processed through aln_clusters and  add_taxonomy_to_alignment) and returns, for each cluster, the number of  taxa at each taxonomic level and their names. The output file has the  following columns: cluster_ID, cluster_rep, level, taxon, count.
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


# SAT aln_merge
This subcommand is used to merge two foldseek alignment files.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_merge -h`](.github/img/aln_merge.png)  


# Planned improvements
struc_get_domains
- Add functionality to parse PAE json files from additional sources

ete3  
- Add ability to specify where the ete3 taxonomy database is downloaded.
