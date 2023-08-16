This repository contains scripts for preparing and splitting the benchmark data we were tasked with examining in our project.
The file csvPrep_dataAug.py contains functions for splitting datasets and augmenting data with noncanonical SMILES.
The file SmallTanimoto.py contains functions for computing tanimoto similarity and cosine similarity from lists of SMILES/embeddings
The file embeddings.ipynb contains code to take MoLFormer embeddings and use them as input in SVMs/RFs
The files finetune_pubchem_light_commented.ipynb and rotary_commented.ipynb are more heavily commented versions of scripts that appear in the original MoLFormer GitHub
original molformer github can be found here: https://github.com/IBM/molformer
