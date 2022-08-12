# SAT
Structural Analysis Toolkit - A python package for manipulation of structural data and identification of shared structural features.

# Installation
1. Clone this repository
2. Create a conda environment that contains the poetry package manager
```
conda create --name SAT -c conda-forge poetry
```
3. Activate environment. Enter the SAT directory that contains the pyproject.toml file, and download the dependencies using poetry!
```
conda activate SAT
poetry install
```
4. Finished! The SAT conda environment will now contain all dependencies.

# SAT get_domains
## Extract separate domain structures from a predicted structure.
### This uses the PAE information to cluster residues that likely fall into linear domains. Notably, the script is currently only configured to process colabfold-generated PAE files. 

Inputs:
1) structure (PDB format)
2) PAE information (json format)

Output:
1) Separate pdb structure files, one per domain.

## Usage
```
python sat.py get_domains \
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

# SAT remove_redundant_domains
## Given a glob specifying multiple structure (often domains), will remove structures that have an overlapping pirmary amino acid sequence. 
### Priority is given to the longer structure or, if the sequences are the same length, the structure with the highest pLDDT

Inputs: 
1) A glob specifying the path to multiple structure files (that may be redundant). 

Output: 
1) Structure files, with redundant structure files removed.

## Usage
```
python sat.py remove_redundant_domains \
-i "path/to/structures/*pdb" \
-o output_directory
```
## Required Parameters
`-i --input_structure_glob`: Glob indicating the structures you want to compare. Make sure to wrap in quotes!

`-o --output_dir`: Path to the directory that will contain the output, filtered structure files. This script will make the directory if it doesn't exist.

# SAT process_clusters
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
python sat.py process_clusters \
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

# SAT add_taxonomy_to_alignments
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
python sat.py add_taxonomy_to_alignments \
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

# SAT structure_to_seq
## Simple subcommand to produce the amino-acid sequence from a structure file.
### Can append the sequence to an outfile if provided, or will print to screen.

Inputs:
1) A structure file in pdb format.

Output:
1) The amino acid sequence printed to the screen or appended to an output file.

## Usage
```
python sat.py structure_to_seq \
-s structure.pdb \
-o sequence.fasta \
-H header_of_the_sequence
```

## Required Parameters
`-s --structure_file`: Path to the structure file.


## Optional Parameters
`-o --output_file`: Path to the output fasta that will be APPENDED to. If not specified, the amino acid sequence will simply be written to the screen.  
`-H --header`: String that is used as the header in the fasta (no need to include the >). This is only required if -o is specified.

# Planned improvements
get_domains
- Add functionality to parse PAE json files from additional sources