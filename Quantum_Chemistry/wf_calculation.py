from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem.rdDistGeom import ETKDGv3, EmbedMolecule
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMolecule
import psi4

import pandas as pd
import time
from timeout_decorator import timeout, TimeoutError

import os
os.environ["CUDA_VISIBLE_DEVICES"] = "0"

original_dataset = pd.read_csv( "tox21.csv")

data_num = len(original_dataset.index)
print('start with {} samples'.format(data_num))

psi4.set_num_threads(nthread=5)
psi4.set_memory("10GB")
timeout_sec = 60*60 

@timeout(timeout_sec)
def opt(level, molecule):
    energy, wave_function = psi4.optimize(level, molecule=molecule, return_wfn=True)
    return energy, wave_function

rm_num = 0 
calc_num = 0 
task_num = 0 

while len(calc_datasets.index) > 0: 

    loop_datasets = calc_datasets.copy()
    
    for i, (ID, smile) in enumerate(zip(loop_datasets['ID'], loop_datasets['smiles'])):
        try:
            #logfile
            print('ID: {}, smiles: {}'.format(ID, smile))
            
            try:
                psi4.set_output_file("../log_files/{}.log".format(ID))
            except SystemError:
                calc_datasets.drop(ID-1, inplace = True) 
                rm_num += 1
                task_num = data_num - rm_num - calc_num
                print('Remaining molecules:', task_num)
                psi4.core.opt_clean()
                
                with open('not_caculate.csv','a') as f:
                    f.write("{}".format(ID)) 
                    f.write(",")
                print("ID=",ID)
                
                print('skiped this molecule: ', smile)
                break 
            
            mol = Chem.MolFromSmiles(smile)   
            mol = Chem.AddHs(mol)   # plus H
            params = ETKDGv3()
            params.randomseed = 1
            EmbedMolecule(mol, params) 
            
            MMFFOptimizeMolecule(mol)
            conf = mol.GetConformer()
            
            plus = smile.count('+')
            minus = smile.count('-')
            charge = plus - minus + 0
            
            mol_input = str(charge) + " 1" # 0 1
            print(str(i+1) + ": " + mol_input + ": "+ smile)
            
            print(mol_input)
            
            for atom in mol.GetAtoms():
                mol_input += "\n" + atom.GetSymbol() + " " + str(conf.GetAtomPosition(atom.GetIdx()).x)\
                + " " + str(conf.GetAtomPosition(atom.GetIdx()).y)\
                + " " + str(conf.GetAtomPosition(atom.GetIdx()).z)
            molecule = psi4.geometry(mol_input) 
            level = "HF/sto-3g"
            psi4.set_options({'geom_maxiter': 1000})
            try:
                print("wavefunction caculation begin!")
                energy, wave_function = opt(level, molecule) 
                print("wavefunction caculation finish!")
            except psi4.OptimizationConvergenceError:
                calc_datasets.drop(ID-1, inplace = True)
                rm_num += 1
                task_num = data_num - rm_num - calc_num
                print('Remaining molecules:', task_num)
                psi4.core.opt_clean()
                with open('not_caculate.csv','a') as f:
                    f.write("{}".format(ID)) # 或者f.write(str(ID))
                    f.write(",")
                print("ID=",ID)
                print('skiped this molecule: ', smile)
                break
            except TimeoutError:
                calc_datasets.drop(ID-1, inplace = True)
                rm_num += 1
                task_num = data_num - rm_num - calc_num
                print('Remaining molecules:', task_num)
                psi4.core.opt_clean()
                with open('not_caculate.csv','a') as f:
                    f.write("{}".format(ID)) 
                    f.write(",")
                print("ID=",ID)
                print('skiped this molecule: ', smile)
                break
        except Exception as e:
            print('Unexpected Error')
            print(e)
            calc_datasets.drop(ID-1, inplace = True)
            rm_num += 1
            task_num = data_num - rm_num - calc_num
            print('Remaining molecules:', task_num)
            with open('not_caculate.csv','a') as f:
                f.write("{}".format(ID))
                f.write(",")
            print("ID=",ID)
            print('skiped this molecule: ', smile)
            psi4.core.opt_clean()
            break
        except Exception:
            print('Rise Other Error')
            calc_datasets.drop(ID-1, inplace = True)
            rm_num += 1
            task_num = data_num - rm_num - calc_num
            print('Remaining molecules:', task_num)
            with open('not_caculate.csv','a') as f:
                f.write("{}".format(ID))
                f.write(",")
            print("ID=",ID)
            print('skiped this molecule: ', smile)
            psi4.core.opt_clean() 
            break
            
        wave_function.to_file("./WaveFunctions/ID_{}_wavefunction".format(ID)) 
        calc_datasets.drop(ID-1, inplace = True)
        calc_num += 1
        task_num = data_num - rm_num - calc_num
        print("completed caluculations: ", calc_num)
        print('Remaining molecules:', task_num)
        
print("Calculation was finished!!")