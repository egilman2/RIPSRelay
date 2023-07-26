import pandas as pd
import os
import random

from rdkit import Chem
from sklearn.model_selection import train_test_split
from rdkit.Chem import PandasTools

"""
Collection of utility functions for data prep
"""

def prep_csvs(path, dir_name):
    """
    Reads in datasets from MoleculeACE and splits them
    :param path: path to dataset csv
    :param dir_name: name of folder to spit out splits
    :return: splits the .csv and puts the splits in your folder
    """
    bigCSV = pd.read_csv(path)
    #make dirs
    os.mkdir(dir_name)
    train_1 = bigCSV[bigCSV['split']=='train']
    test = bigCSV[bigCSV['split']=='test']
    train, valid = train_test_split(train_1, test_size=.1)
    train.to_csv(dir_name+'/train.csv')
    test.to_csv(dir_name + '/test.csv')
    valid.to_csv(dir_name + '/valid.csv')


def prep_sdf(path):
    """
    Preps and splits data from the other paper
    :param path: path to csv
    :return: splits them
    """
    imp_sdf = PandasTools.LoadSDF(path)
    if 'ROMol' not in imp_sdf.columns:
        print('no romol')
        return None
    imp_sdf['smiles'] = [Chem.MolToSmiles(mol) for mol in imp_sdf['ROMol']]
    if 'ChEMBL_ID' in imp_sdf.columns:
        imp_sdf = imp_sdf.drop(columns=['ChEMBL_ID'])
    if 'ID' in imp_sdf.columns:
        imp_sdf = imp_sdf.drop(columns=['ID'])
    train_1, test = train_test_split(imp_sdf, test_size=.1)
    train, valid = train_test_split(train_1, test_size=.1)
    dir_name = path.split('/')[-1][:-4] + '_split'
    os.makedirs(dir_name)
    train.to_csv(dir_name + '/train.csv')
    test.to_csv(dir_name + '/test.csv')
    valid.to_csv(dir_name + '/valid.csv')


# stolen from moleculeACE


def random_smiles(mol):
    """ Generate a random non-canonical SMILES string from a molecule"""
    # https://github.com/michael1788/virtual_libraries/blob/master/experiments/do_data_processing.py
    mol.SetProp("_canonicalRankingNumbers", "True")
    idxs = list(range(0, mol.GetNumAtoms()))
    random.shuffle(idxs)
    for i, v in enumerate(idxs):
        mol.GetAtomWithIdx(i).SetProp("_canonicalRankingNumber", str(v))
    return Chem.MolToSmiles(mol)


def smile_augmentation(smile: str, augmentation: int, max_len: int = 200, max_tries: int = 1000):
    """Generate n random non-canonical SMILES strings from a SMILES string with length constraints"""
    # https://github.com/michael1788/virtual_libraries/blob/master/experiments/do_data_processing.py
    mol = Chem.MolFromSmiles(smile)
    s = set()
    for i in range(max_tries):
        if len(s) == augmentation:
            break

        smiles = random_smiles(mol)
        if len(smiles) <= max_len:
            # tokens = smi_tokenizer(smiles)
            # if all([tok in smiles_encoding['token_indices'] for tok in tokens]):
            s.add(smiles)

    return list(s)


def data_augment(df: pd.DataFrame, num):
    """
    Augments the data in an activity cliff dataset by throwing in non-canonical copies of SMILES
    :param df: base data
    :param num: how many new smiles per o.g. smile
    :return: a dataframe that has been augmented with new smiles
    """
    new_rows = []
    for index, row in df.iterrows():
        smiles = row['smiles']
        y_val = row['y']
        new_smiles = smile_augmentation(smiles, num)
        for smile in new_smiles:
            new_rows.append({'smiles': smile, 'y': y_val})
    return pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)


