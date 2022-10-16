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
## Extract separate domain structures from a predicted structure.
### This uses the PAE information to cluster residues that likely fall into linear domains. Notably, the script is currently only configured to process colabfold-generated PAE files. 

Inputs:
1) structure (PDB format)
2) PAE information (json format)

Output:
1) Separate pdb structure files, one per domain.

## Usage
```
sat.py struc_get_domains \
-s sturucture.pdb \
-p pae.json \
-o output/structure_prefix
```

## Required Parameters
`-s --structure_file_path`: Path to the input structure. 

`-p --pae_path`: Path to the colabfold-generated PAE json file. Notably, this script is designed to parse the PAE json specifically from colabfold. 

`-o --output_prefix`: Path specifying the output prefix for the output structure. Files will be labled {prefix}_domain-{i}.pdb. Note that the domain number will be 1-indexed.

## Optional Parameters
`-1 --pae_power [1]`: Used when clustering residues using PAE. EAch edge in the graph will be weighted proportional to (1/pae**pae_power).

`-2 --pae_cutoff [5]`: Used when clustering residues using PAE. Lowering this will make domain identification more stringent by reducing the amount of error allowed.

`-3 --graph_resolution [1]`: Regulates how aggressively the clustering algorithm is. Smaller values lead to larger clusters. Value should be larger than zero, and values larger than 5 are unlikely to be useful.

# SAT struc_remove_redundant
## Given a glob specifying multiple structure (often domains), will remove structures that have an overlapping pirmary amino acid sequence. 
### Priority is given to the longer structure or, if the sequences are the same length, the structure with the highest pLDDT

Inputs: 
1) A glob specifying the path to multiple structure files (that may be redundant). 

Output: 
1) Structure files, with redundant structure files removed.

## Usage
```
sat.py struc_remove_redundant \
-i "path/to/structures/*pdb" \
-o output_directory
```
## Required Parameters
`-i --input_structure_glob`: Glob indicating the structures you want to compare. Make sure to wrap in quotes!

`-o --output_dir`: Path to the directory that will contain the output, filtered structure files. This script will make the directory if it doesn't exist.

# SAT aln_clusters
## Adds cluster information from foldseek cluster into the foldseek alignment information.
### Currently, foldseek cluster can generate clusters, and the alignment can output alignments, but it will be helpful to annotate each alginment with their cluster. Furthermore, all-by-all searchs  result in a redundant query-target and target-query alignment for each entry. These redundant alignments are removed.
### This script adds the following fields to the input alignment file:
### 1) cluster_ID: ID of the cluster, starting at 1. A lower number indicates a larger cluster
### 2) cluster_rep: The name of the structure that is the cluster representative chosen by foldseek
### 3) cluster_count: The number of structures in the cluster.

Inputs:
1) A foldseek alignment file
2) A foldseek cluster file

Output:
1) The foldseek alignment file with removed redundant alignments and labled with cluster information.

## Usage
```
sat.py aln_clusters \
-a alignment.m8 \
-c clusters.tsv \
-o result.m8
```

## Required Parameters
`-a --alignment_file`: Path to the foldseek alignment file.

`-c --cluster_file`: Path to the foldseek cluster file.

`-o --output_file`: Path to the output file.

## Optional Parameters
`-f --alignment_fields ["query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore"]`: A comma-delimited string of the fields in the input foldseek alignment file. Make sure to wrap in quotes!

# SAT aln_add_taxonomy
## Adds taxonomy information to a foldseek alignment
### Taxonomy for each query can be built in to the query_name as the final element of the double-underscore-delimited list. 
### Taxonomy for each target can either be built in to the target name in a similar manner, or present as the taxid field in the foldseek output.
### This script adds the taxonomic names of the query and/or target at the specified levels. 

Inputs:
1) A foldseek alignment file

Output:
1) The foldseek alignment file with taxonomy information

## Usage
```
sat.py aln_add_taxonomy \
-a alignment.m8 \
-o result.m8
```

## Required Parameters
`-a --alignment_file`: Path to the foldseek alignment file.

`-o --output_file`: Path to the output file.

## Optional Parameters
`-f --alignment_fields ["query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore,cluster_ID,cluster_rep,cluster_count"]`: A comma-delimited string of the fields in the input foldseek alignment file. Make sure to wrap in quotes!

`-q --query_taxid_location [1]`: Options are {0, 1}.

    0: Indicates that the query_taxid is not present or not desired in the output.

    1: Indicates that the query_taxid is present in the query name of the alignment as query_name.rstrip('.pdb').split('__')[-1]

`-t --target_taxid_location [1]`: Options are {0, 1, 2}.

    0: Indicates that the target_taxid is not present or not desired in the output.

    1: Indicates that the target_taxid is present in the query name of the alignment as query_name.rstrip('.pdb').split('__')[-1]

    2: Indicates that the target_taxid is present in the taxid field of the alignment file.

`-T --taxonomy_levels ["superkingdom,phylum,class,order,family,genus,species"]`: Comma-delimited string of the levels you want as taxonomy information in the output. Make sure to wrap in quotes.

# SAT struc_to_seq
## Simple subcommand to produce the amino-acid sequence from a structure file.
### Can append the sequence to an outfile if provided, or will print to screen.

Inputs:
1) A structure file in pdb format.

Output:
1) The amino acid sequence printed to the screen or appended to an output file.

## Usage
```
sat.py struc_to_seq \
-s structure.pdb \
-o sequence.fasta \
-H header_of_the_sequence
```

## Required Parameters
`-s --structure_file`: Path to the structure file.


## Optional Parameters
`-o --output_file`: Path to the output fasta that will be APPENDED to. If not specified, the amino acid sequence will simply be written to the screen.  
`-H --header`: String that is used as the header in the fasta (no need to include the >). This is only required if -o is specified.

# SAT struc_rebase
## Simple subcommand that renumbers all residues in a structure such that the first residue is #1 and all residues are sequential (e.g. it takes out numeric gaps in residue numbers).

Inputs:
1) A structure file in pdb format.

Output:
1) Path to the rebased structure output file.

## Usage
```
sat.py struc_rebase \
-s structure.pdb \
-o rebased_structure.pdb
```

## Required Parameters
`-s --structure_file`: Path to the structure file.  
`-o --output_file`: Path to the output, rebased structure.


# SAT seq_chunk
## Splits entries into a fasta into overlapping or non-overlapping chunks. This is helpful when you want to split up sequences that are too long to effectively use for structure prediction. This subcommand is able to generate overlapping sequences. 

Inputs:
1) A fasta containing one or more sequences.

Output:
1) Either a single fasta or multiple fastas (named by header) containing overlapping or not overlapping sequences with a maximum and minimum length.

## Usage
```
sat.py seq_chunk \
-i input.fasta \ 
-o output.fasta # or output_dir/ if -n is specified
```

## Required Parameters
`-i --in_fasta`: Path to the input fasta file.  
`-o --out_fasta`: Path to a singular fasta output file. If -n is specified and multiple fastas are desired, the entry here will be used as the BASE directory and output files will be {BASE}/{header}.fasta
`-m --max_seq_length`: Maximum length of an output sequence. If an input sequence is longer than this, it will be chunked.  

## Optional Parameters
`-s --minimum_sequence_output_size [50]`: Minimum length of an output sequence. If a sequence or chunk is below this size it will be discarded.  
`-v --overlapping_chunks [False]`: This is a boolean flag. If you specify -v in the command line, each chunk will overlap by round(max_seq_length/2). If not specified, chunks will not be overlapping.  
`-n --individual [False]`: This is a boolean flag. If you specify -n in the command line, each resultant fasta will be passed to a separate file named by its header. It will assume that --out_fasta specifies the path to the output directory... individual filenames will be named by their header.

# SAT aln_taxa_counts
## This takes in a processed alignment file (typically generated from a foldseek alignment that was then processed through aln_clusters and  add_taxonomy_to_alignment) and returns, for each cluster, the number of  taxa at each taxonomic level and their names. The output file has the  following columns: cluster_ID, cluster_rep, level, taxon, count.

Inputs:
1) A foldseek alignment file (typically processed through aln_clusters and add_taxonomy_to_alignment).

Output:
1) A tidy output file with counts, for each cluster, for the number of members of each cluster in every taxon at each desired level.

## Usage
```
sat.py aln_taxa_counts \
-a alignment.m8 \ 
-o result.m8
```

## Required Parameters
`-a --alignment_file`: Path to the alignment file.
`-o --outfile`: Path to the output count file.

## Optional Parameters
`-t --taxonomy_levels ["superkingdom,phylum,class,order,family,genus,species"]`: Taxonomy levels to count and output

# SAT aln_query_uniprot
## This script takes alphafold IDs (or raw uniprot IDs) and uses the Uniprot REST API to get information on the geneName and fullName (an informative name of the protein) for each ID. You can specify in which column of the infile the IDs live.

Inputs:
1) A single-column or tab-delimited file, with one of the columns containing alphafold IDs or uniprot IDs. Note that the ID will be found as the second item in a dash-delimited list if dashes are present (as they are in alphafold IDs).
2) Optional, but recommended: A cache file or file path. This will save the downloaded uniprot contents to a pkl binary file. Subsequent future runs will only query uniprot if the information is not in the cache.

Output:
1) A tab-delimited file with the columns uniprotID, geneName, and fullName

## Usage
```
sat.py aln_query_uniprot \
-i some_file.m8 \ 
-o result.tsv \
-c 1 \
-u cache.pkl
```

## Required Parameters
`-i --infile`: Path to the input file that contains the uniprot IDs. If there are multiple columns, this must be tab-delimited.
`-o --uniprot_lookup_output`: This is the main output file. Will be a tab-delimited file with the columns uniprotID, geneName, fullName.

## Optional Parameters
`-c --infile_col [1]`: This is the 0-indexed column holding the uniprot of alphafold IDs.
`-u --uniprot_cache`: A pkl file containing a cache from former uniprot downloads. This is optional. If specified, this file will be read in and will be updated with this script's downloads.

# SAT aln_add_uniprot
## This script adds the uniprot information garther from aln_query_uniprot to a foldseek alignment file.

Inputs:
1) A foldseek alignment file.
2) Uniprot information produced by aln_query_uniprot.

Output:
1) The alignment file with uniprot information added.

## Usage
```
sat.py aln_add_uniprot \
-a alignment.m8 \
-u uniprot_info.tsv \
-o alignment_w_uniprot.m8
```

## Required Parameters
`-a --alignment_file`: Path to the alignment file. It is OK if the first row is the header, as long as the first column is 'query'.  
`-u --uniprot_information`: Path to the uniprot information file generated by aln_query_uniprot.  
`-o --output_file`: Path to the output file.  

## Optional Parameters
`-f --alignment_fields [""]`: Comma-delimited list of fields in the input alignment files (aka colnames). If the colnames are the first line in the input file and the first column is 'query', you can leave this blank.  

# SAT aln_filter
## This subcommand filters a foldseek alignment file to keep only those alignments with a value below/above the specified value in a field (alntmscore is a common one). It also only outputs a maximum of N alignments for each query.

Inputs:
1) A foldseek alignment file.

Output:
1) The filtered alignment file

## Usage
```
sat.py aln_filter \
-a alignment.m8 \
-o alignment.filtered.m8 \
-x alntmscore \
-N 10 \
-m 0.4
```
The above filters for alignments that have an alntmscore of at least 0.4, and keeps a max of 10 alignments per query.

## Required Parameters
`-a --alignment_file`: Path to the alignment file. It is OK if the first row is the header, as long as the first column is 'query'.   
`-o --output_file`: Path to the output file.  

## Optional Parameters
`-f --alignment_fields [""]`: Comma-delimited list of fields in the input alignment files (aka colnames). If the colnames are the first line in the input file and the first column is 'query', you can leave this blank.  
`-x --filter_field ["alntmscore"]`: String indicating the alignment field which contains the value to be filtered/sorted with.  
`-N --N [10]`: This is the maximum number of alignments to output for each query. Set to 0 if you want to return all alignments.  
`-m --min_val_filter_field [0.4]`: An alignment must have at least this minimum value in the filter_field for the alignment to be output.  
`-M --max_val_filter_field [1]`: An alignment must have this value or less in the filter_field for the alignment to be output.  


# SAT aln_merge
## This subcommand is used to merge two foldseek alignment files.

## Usage
![`poetry run sat.py aln_merge -h`](.github/img/aln_merge.png)


# Planned improvements
struc_get_domains
- Add functionality to parse PAE json files from additional sources

ete3  
- Add ability to specify where the ete3 taxonomy database is downloaded.