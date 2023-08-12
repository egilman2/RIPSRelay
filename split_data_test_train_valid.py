import pandas as pd
import os
import random

from rdkit import Chem
from sklearn.model_selection import train_test_split
from rdkit.Chem import PandasTools

def prep_csvs(path, dir_name):
    bigCSV = pd.read_csv(path)
    #make dirs
    os.mkdir(dir_name)
    train_1 = bigCSV[bigCSV['split']=='train']
    test = bigCSV[bigCSV['split']=='test']
    train, valid = train_test_split(train_1, test_size=.1)
    train.to_csv(dir_name+'/train.csv')
    test.to_csv(dir_name + '/test.csv')
    valid.to_csv(dir_name + '/valid.csv')

def main():
    fs = os.listdir('/home/joiaz/molformer/data/data_7') # returns list of filenames in directory

    for f in fs:
        path = '/home/joiaz/molformer/data/data_7/' + f
        dir_name = '/home/joiaz/molformer/data/data_7_split/' + f.replace('.csv', '')
        prep_csvs(path, dir_name)
    
if __name__ == '__main__':
    main()






    