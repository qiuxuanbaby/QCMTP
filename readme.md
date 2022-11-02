# QuantumMTP: Utilizing quantum chemistry with Ensemble Learning for Molecular Toxicity Prediction

## Requirements
RDKit
PSi4
sklearn
numpy
matplotlib
pandas
seaborn

To install RDKit, please follow the instructions: http://www.rdkit.org/docs/Install.html

## Quantum_Chemistry 

a. wf_calculation.py contains code for calculating wavefunctions.
b. The WaveFunctions folder contains some examples of wavefunctions
c. data_process.py contains code for processing molecules that cannot calculate wavefunctions.
d. get_wf_results.py contains code for getting values from wavefunctions.

## Ensemble_Learning
Ablation+Different_ML contains code for ablation studies and comparisons with different machine learning models.
imbalance.py contains code for experiments with imbalanced data sets.
