# SAT
Structural Analysis Toolkit - A python package for manipulation of structural data and identification of shared structural features.

# Installation
to do

# Quick start

## SAT get_domains
### Extract separate domain structures from a predicted structure.

Inputs:
1) structure (PDB format)
2) PAE information (json format)

Output:
1) Separate pdb structure files, one per domain.

### Usage
```
python sat.py get_domains \
-s sturucture.pdb \
-p pae.json \
-o output/structure_prefix
```

### Required Parameters
`-s --structure_file_path`: Path to the input structure. 

`-p --pae_path`: Path to the colabfold-generated PAE json file. Notably, this script is designed to parse the PAE json specifically from colabfold. 

`-o --output_prefix`: Path specifying the output prefix for the output structure. Files will be labled {prefix}_domain-{i}.pdb. Note that the domain number will be 1-indexed.

### Optional Parameters
`-1 --pae_power [1]`: Used when clustering residues using PAE. EAch edge in the graph will be weighted proportional to (1/pae**pae_power).

`-2 --pae_cutoff [5]`: Used when clustering residues using PAE. Lowering this will make domain identification more stringent by reducing the amount of error allowed.

`-3 --graph_resolution [1]`: Regulates how aggressively the clustering algorithm is. Smaller values lead to larger clusters. Value should be larger than zero, and values larger than 5 are unlikely to be useful.

# Planned improvements
get_domains
- Add functionality to parse PAE json files from additional sources