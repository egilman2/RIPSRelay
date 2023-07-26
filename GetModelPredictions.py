import pandas as pd
import os
import random
import SmallTanimoto
import numpy as np
import matplotlib

pred_df = pd.read_csv('predictions465.csv')
analysis_df = pred_df.copy()
#print(analysis_df['preds'])
analysis_df['abs_error'] = analysis_df.apply(lambda x: np.abs(x['preds'] - x['actual']), axis=1)
test_df = pd.read_csv('219_prep/test.csv')
train_df = pd.read_csv('219_prep/train.csv')
train_smiles = list(train_df['smiles'])
test_smiles = list(test_df['smiles'])
score_map = SmallTanimoto.score_test_smiles(train_smiles, test_smiles)[1]
scores = []
for smile in test_smiles:
    scores.append(score_map[smile])
analysis_df['smiles'] = test_smiles
analysis_df['scores'] = scores
analysis_df.to_csv('219_analysis2.csv')

an_df = pd.read_csv('219_analysis2.csv')
an_df['cliff_mol'] = test_df['cliff_mol']
an_df.to_csv('219_analysis_wcliff2.csv')





