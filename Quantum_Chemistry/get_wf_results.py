# get values from wavefunctions

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, PandasTools
from rdkit.Chem.rdDistGeom import ETKDGv3, EmbedMolecule
import psi4
import pandas as pd
import numpy as np
import glob
import re

flies = []
for file in glob.glob("./WaveFunctions/ID_*_wavefunction"):
    files.append(file)
data_num = len(files)

# Settings of PSI4
psi4.set_num_threads(nthread=5)
psi4.set_memory("10GB")

ID_list = []
wf_list = []
for file in files:
    pre_ID = re.findall('./WaveFunctions/ID_(.*)_wavefunction.npy', file)
    ID = pre_ID[0]
    ID_list.append(int(ID))
    wave_function = psi4.core.Wavefunction.from_file(file)
    wf_list.append(wave_function) 

data_sets = pd.DataFrame({'ID': ID_list, 'wavefunction': wf_list})
data_sets = data_sets.sort_values('ID')
sort_data_sets = data_sets.reset_index(drop=True)
original_df = pd.read_csv( 'tox21.csv')

smiles_list = []
NR_AR_list = []
NR_AR_LBD_list = []
NR_AhR_list = []
NR_Aromatase_list = []
NR_ER_list = []
NR_ER_LBD_list = []
NR_PPAR_gamma_list = []
SR_ARE_list = []
SR_ATAD5_list = []
SR_HSE_list = []
SR_MMP_list = []
SR_p53_list = []

for ID in sort_data_sets.ID :
    same_df = original_df.query('ID == @ID')
    smile = same_df['smiles'].values[0]
    NR_AR = same_df['NR-AR'].values[0]
    NR_AR_LBD = same_df['NR-AR-LBD'].values[0]
    NR_AhR = same_df['NR-AhR'].values[0]
    NR_Aromatase = same_df['NR-Aromatase'].values[0]
    NR_ER = same_df['NR-ER'].values[0]
    NR_ER_LBD = same_df['NR-ER-LBD'].values[0]
    NR_PPAR_gamma = same_df['NR-PPAR-gamma'].values[0]
    SR_ARE = same_df['SR-ARE'].values[0]
    SR_ATAD5 = same_df['SR-ATAD5'].values[0]
    SR_HSE = same_df['SR-HSE'].values[0]
    SR_MMP = same_df['SR-MMP'].values[0]
    SR_p53 = same_df['SR-p53'].values[0]
    
    smiles_list.append(smile)
    NR_AR_list.append(NR_AR)
    NR_AR_LBD_list.append(NR_AR_LBD)
    NR_AhR_list.append(NR_AhR)
    NR_Aromatase_list.append(NR_Aromatase)
    NR_ER_list.append(NR_ER)
    NR_ER_LBD_list.append(NR_ER_LBD)
    NR_PPAR_gamma_list.append(NR_PPAR_gamma)
    SR_ARE_list.append(SR_ARE)
    SR_ATAD5_list.append(SR_ATAD5)
    SR_HSE_list.append(SR_HSE)
    SR_MMP_list.append(SR_MMP)
    SR_p53_list.append(SR_p53)
    
sort_data_sets['smiles'] = smiles_list
sort_data_sets['NR-AR'] = NR_AR_list
sort_data_sets['NR-AR-LBD'] = NR_AR_LBD_list
sort_data_sets['NR-AhR'] = NR_AhR_list
sort_data_sets['NR-Aromatase'] = NR_Aromatase_lis
sort_data_sets['NR-ER'] = NR_ER_list
sort_data_sets['NR-ER-LBD'] = NR_ER_LBD_list
sort_data_sets['NR-PPAR-gamma'] = NR_PPAR_gamma_list
sort_data_sets['SR-ARE'] = SR_ARE_list
sort_data_sets['SR-ATAD5'] = SR_ATAD5_list
sort_data_sets['SR-HSE'] = SR_HSE_list
sort_data_sets['SR-MMP'] = SR_MMP_list
sort_data_sets['SR-p53'] = SR_p53_list


# Make data_list
energy_list = []
homo_list = []
lumo_list = []
HLgap_list = []
charge_ave_list = []
charge_var_list = []
dipole_list = []
index_list = []
mass_list = []
vol_list = []
dens_list = []

# Energy，homo，lumo，HLgap，The mean and variance of the charges，dipole，mass，volume，density are extracted from wavefuntions
for wave_function,ID,smile in zip(sort_data_sets.wavefunction,sort_data_sets.ID,sort_data_sets.smiles):
    rdmol = Chem.AddHs(Chem.MolFromSmiles(smile))  
    AllChem.EmbedMolecule(rdmol) 
    vol = AllChem.ComputeMolVolume(rdmol)  
    vol_list.append(vol)
    energy = wave_function.energy()
    energy_list.append(energy)
    LUMO_idx = wave_function.nalpha()
    HOMO_idx = LUMO_idx - 1
    homo =wave_function.epsilon_a_subset("AO", "ALL").np[HOMO_idx]
    lumo = wave_function.epsilon_a_subset("AO", "ALL").np[LUMO_idx]
    homo_list.append(homo)
    lumo_list.append(lumo)
    HLgap = lumo - homo
    HLgap_list.append(HLgap)
    psi4.oeprop(wave_function, "MULLIKEN_CHARGES")
    charges = np.array(wave_function.atomic_point_charges())
    charge_ave_list.append(np.average(charges)) 
    charge_var_list.append(np.var(charges)) 
    dipole_x, dipole_y, dipole_z = wave_function.variable('SCF DIPOLE X'), wave_function.variable('SCF DIPOLE Y'), wave_function.variable('SCF DIPOLE Z')
    dipole_moment = np.sqrt(dipole_x ** 2 + dipole_y ** 2 + dipole_z ** 2)
    dipole_list.append(dipole_moment)
    mol = wave_function.molecule() 
    atom_num = mol.natom() 
    atom_mass_list = []
    for n in range(0,atom_num):
        mass = mol.mass(n)
        atom_mass_list.append(mass)
    mol_mass = sum(atom_mass_list)
    mass_list.append(mol_mass)
    atom_mass_list.clear()
    dens = mol_mass/vol
    dens_list.append(dens)
sort_data_sets['energy'] = energy_list
sort_data_sets['homo'] = homo_list
sort_data_sets['lumo'] = lumo_list
sort_data_sets['HLgap'] = HLgap_list
sort_data_sets['charge_ave'] = charge_ave_list
sort_data_sets['charge_var'] = charge_var_list
sort_data_sets['dipole'] = dipole_list
sort_data_sets['mass'] = mass_list
sort_data_sets['vol'] = vol_list
sort_data_sets['dens'] = dens_list

print_data_sets = sort_data_sets.drop(columns='wavefunction')
print_data_sets.to_csv("Tox21_wf_results_activity_{}samples.csv".format(data_num),index=False)