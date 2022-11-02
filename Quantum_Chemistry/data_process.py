from rdkit import rdBase, Chem
from rdkit.Chem import PandasTools
import pandas as pd
import os 

df = pd.read_csv( 'tox21_process.csv', header= 0, sep = ',')
PandasTools.AddMoleculeColumnToFrame(frame=df, smilesCol = 'smiles') # mol objects
df[ 'MOL'] = df.ROMol.map(lambda x: False if x == None else True) 
del_index = df[df.MOL == False].index
df1 = df.drop(del_index)
df1.drop('ROMol',axis=1,inplace=True) 
df1.drop('MOL',axis=1,inplace=True) 

not_caculate = pd.read_csv('not_caculate.csv',sep=',',header=None)
not_caculate_list = not_caculate.iloc[0,0:].tolist()
print(len(not_caculate_list))
for i in not_caculate_list:
    df1 = df1.drop(i-1)
    
df1.to_csv("tox21_process.csv")
