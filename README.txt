This repository contains scripts for preparing and splitting the benchmark data we were tasked with examining in our project.
The file csvPrep_dataAug.py contains functions for splitting datasets and augmenting data with noncanonical SMILES.
The file SmallTanimoto.py contains functions for computing tanimoto similarity and cosine similarity from lists of SMILES/embeddings
The file embeddings.ipynb contains code to take MoLFormer embeddings and use them as input in SVMs/RFs
The files finetune_pubchem_light_commented.ipynb and rotary_commented.ipynb are more heavily commented versions of scripts that appear in the original MoLFormer GitHub
The Datasets folder contains datasets used in this paper: https://pubs.acs.org/doi/10.1021/acs.jcim.8b00542
The benchmark_data contains activity cliffs datasets used in this paper: https://pubs.acs.org/doi/10.1021/acs.jcim.2c01073
original molformer github can be found here: https://github.com/IBM/molformer
