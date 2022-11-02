#rdkit
from rdkit import Chem
from rdkit.Chem import PandasTools
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report
from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import BaggingClassifier
from sklearn.model_selection import cross_val_score

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

norm_datasets = pd.read_csv("NR-Aromatase_209_with_activity.csv") #NR-Aromatase_412_with_activity.csv NR-Aromatase_309_with_activity.csv
columns = norm_datasets.columns

PandasTools.AddMoleculeColumnToFrame(frame=norm_datasets, smilesCol = 'smiles')

original = pd.read_csv( '../Quantum_Chemistry/tox21.csv')

# fingerprint
rdkit_fp = []
for mol in norm_datasets.ROMol:
    mol = Chem.AddHs(mol) 
    fp = [x for x in Chem.RDKFingerprint(mol)]
    rdkit_fp.append(fp)
    
calc_datasets = norm_datasets.iloc[:,3:-1] #3D
rdkit_fp_df = pd.DataFrame(rdkit_fp) #2D
calc_sets = pd.concat([calc_datasets,rdkit_fp_df], axis=1) #2D+3D

norm_datasets['activity'].value_counts()

#Only Rdkit Finger Print
#Calculation Results + Rdkit Finger Print
forest1 = GradientBoostingClassifier(n_estimators=500,max_depth=7)
bagging_forest1 = BaggingClassifier(forest1,n_estimators=60) 
scores1 = cross_val_score(bagging_forest1, rdkit_fp_df, norm_datasets.activity, cv=10, scoring="average_precision") 
print('=======Only Rdkit Finger Print========')
print('scores1=',scores1)

#Calculation Results + Rdkit Finger Print
forest2 = GradientBoostingClassifier(n_estimators=500,max_depth=7)
bagging_forest2 = BaggingClassifier(forest2,n_estimators=60) 
scores2 = cross_val_score(bagging_forest2, calc_sets, norm_datasets.activity, cv=10, scoring="average_precision")
print('=======Calculation Results + Rdkit Finger Print========')
print('scores2=',scores2)

#Only Calculation Results
forest3 = GradientBoostingClassifier(n_estimators=500,max_depth=7)
bagging_forest3 = BaggingClassifier(forest3,n_estimators=60) 
scores3 = cross_val_score(bagging_forest3, calc_datasets, norm_datasets.activity, cv=10, scoring="average_precision") 
print('=======Only Calculation Results========')
print('scores3=',scores3)