import pandas as pd
import os
import random
import numpy as np
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import rdFingerprintGenerator

"""
Collection of functions to determine Tanimoto similarity of train/test sets of fingerprints
    -find_most_similar_training_mol takes a set of training smiles and a test smiles and returns the training smile most
    similar the test smiles, as well as the pair's similarity score
    -score_test_smiles essentially does that for every test smile in a set, and returns dictionaries that map
    test smiles to their closest training smiles/associated scores
"""

def find_most_similar_training_mol(training_smiles, test_smile):
    """
    Searches through a list of smiles in the training data to find the most tanimoto-similar one to the input smile
    :param training_smiles: list of smiles in the training data
    :param test_smile: smile for which you want to find the closest training smile
    :return: the closest training smile to the input smile, as well as their tanimoto similarity
    """

    t_mols = [Chem.MolFromSmiles(smi) for smi in training_smiles]
    test_mol = Chem.MolFromSmiles(test_smile)
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    t_fingies = [mfpgen.GetFingerprint(mol) for mol in t_mols]
    test_fingie = mfpgen.GetFingerprint(test_mol)
    scores = Chem.DataStructs.BulkTanimotoSimilarity(test_fingie, t_fingies)
    score_i = np.argmax(scores)
    return training_smiles[score_i], scores[score_i]


def score_test_smiles(training_smiles, test_smiles):
    """
    Finds the most closely related training smile for every test smile, and gives their Tanimoto similarity
    :param training_smiles:
    :param test_smiles:
    :return: two maps; one of test smiles to their closest training smile, and one to their Tanimoto similarity
    """
    score_map = {}
    trainer_map = {}
    for i in range(len(test_smiles)):
        smile = test_smiles[i]
        print(str(i) + '/' + str(len(test_smiles)))
        trainer_map[smile], score_map[smile] = find_most_similar_training_mol(training_smiles, smile)
    return trainer_map, score_map

def pairwise_tanimoto_similarity(smiles):
    """
    Takes a list of smiles of length i and computes pairwise tanimoto similarity for all of them
    :param smiles: a list of smiles
    :return: a dictionary where key (i,j) maps to the tanimoto similarity between the smile at index i and
    the smile at index j
    """
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    mols = [Chem.MolFromSmiles(smile) for smile in smiles]
    fingies = [mfpgen.GetFingerprint(mol) for mol in mols]
    results = {}
    for i in range(len(smiles)):
        for j in range(i, len(smiles)):
            results[(i, j)] = Chem.DataStructs.TanimotoSimilarity(fingies[i], fingies[j])
            if i != j:
                results[(j, i)] = results[(i, j)]
    return results

def pairwise_euclidean_distance(embeddings):
    """
    Takes a 2D numpy array and computes the pairwise euclidean distance between the row vectors
    :param embeddings: a 2D numpy array
    :return: a dictionary where the key (i,j) maps to the Euclidean distance between the ith and jth
    row vectors of the embeddings array
    """
    results = {}
    for i in range(len(embeddings)):
        for j in range(i, len(embeddings)):
            results[(i, j)] = np.linalg.norm(embeddings[i] - embeddings[j])
            if i != j:
                results[(j, i)] = results[(i, j)]
    return results

def pairwise_cosine_similarity(embeddings):
    """
    Takes a 2D numpy array and computes the pairwise cosine similarity between the row vectors
    :param embeddings: a 2D numpy array
    :return: a dictionary where the key (i,j) maps to the cosine similarity between the ith and jth
    row vectors of the embeddings array
    """
    results = {}
    for i in range(len(embeddings)):
        for j in range(i, len(embeddings)):
            results[(i, j)] = np.dot(embeddings[i], embeddings[j])/(np.linalg.norm(embeddings[i]) * np.linalg.norm(embeddings[j]))
            if i != j:
                results[(j, i)] = results[(i, j)]
    return results

rand_smiles = ['CC(c1ccccc1)N1CC[C@H]1[C@@H](N)c1cccc(Cl)c1', 'N[C@@H](c1ccc(Cl)cc1)[C@@H]1CCN1C(c1ccccc1)c1ccccc1',
               'CCNC(=O)c1ccc(C(=C2CC3CCC(C2)N3c2cccn2C)c2ccccc2)cc1',
               'CC(C)C[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@@H]1CCCN1C(=O)CNC(=O)[C@@H](N)Cc1ccc(O)cc1)C(N)=O']
#
# rand_tester = 'CCNC(=O)c1ccc(C(=C2CC3CCC(C2)N3CCO)c2ccccc2)cc1'
#
# print(find_most_similar_training_mol(rand_smiles, rand_tester))
#print(pairwise_cosine_similarity(rand_smiles).items())
# train_smis = list(pd.read_csv('219_prep/train.csv')['smiles'])
# test_smis = list(pd.read_csv('219_prep/test.csv')['smiles'])
# results = score_test_smiles(train_smis, test_smis)
# for smi in test_smis:
#     print(smi, results[0][smi], results[1][smi])