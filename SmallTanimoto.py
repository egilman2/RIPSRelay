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


# rand_smiles = ['CC(c1ccccc1)N1CC[C@H]1[C@@H](N)c1cccc(Cl)c1', 'N[C@@H](c1ccc(Cl)cc1)[C@@H]1CCN1C(c1ccccc1)c1ccccc1',
#                'CCNC(=O)c1ccc(C(=C2CC3CCC(C2)N3c2cccn2C)c2ccccc2)cc1',
#                'CC(C)C[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@@H]1CCCN1C(=O)CNC(=O)[C@@H](N)Cc1ccc(O)cc1)C(N)=O']
#
# rand_tester = 'CCNC(=O)c1ccc(C(=C2CC3CCC(C2)N3CCO)c2ccccc2)cc1'
#
# print(find_most_similar_training_mol(rand_smiles, rand_tester))

# train_smis = list(pd.read_csv('219_prep/train.csv')['smiles'])
# test_smis = list(pd.read_csv('219_prep/test.csv')['smiles'])
# results = score_test_smiles(train_smis, test_smis)
# for smi in test_smis:
#     print(smi, results[0][smi], results[1][smi])