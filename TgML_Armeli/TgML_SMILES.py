from __future__ import annotations

''' this script includes code snippets for loading pickle files and returning Tg prediction
authored by Gianluca Armeli'''

''' AZ notes: This script variant was derived from TgML_minimal.py, but adapted (simplified) by Andi Zuend  
    for faster loading of the specific pickle file and a list of SMILES from file for a 
    SMILES-only-input model (model_name = ['sm_no_tm']) prediction. 
    It generates an output file that lists the Tg and corresponding SMILES, one per line'''

''' input formats for feat '''
# sm_no_tm: feat = [smiles1, smiles2, ...]

print("... importing modules for TgML_SMILES \n")

import pickle
import numpy as np
#import logging
import sys
import os.path
from rdkit import Chem
from rdkit.Chem import Descriptors, MolFromSmiles #, AllChem
from rdkit import DataStructs
from deepchem.utils.typing import RDKitMol
from deepchem.feat.base_classes import MolecularFeaturizer

print("... end of importing modules TgML_SMILES \n")
    

# -----------------------------------------------------
# classes and functions for SMILES processing
# -----------------------------------------------------
class RDKitDescriptors(MolecularFeaturizer):
    def __init__(self, use_fragment=True, ipc_avg=True):
        self.use_fragment = use_fragment
        self.ipc_avg = ipc_avg
        self.descriptors = []
        self.descList = []
        
    def _featurize(self, mol: RDKitMol) -> np.ndarray:
        # initialize
        if len(self.descList) == 0:
            try:
                for descriptor, function in Descriptors.descList:
                    if self.use_fragment is False and descriptor.startswith('fr_'):
                        continue
                    self.descriptors.append(descriptor)
                    self.descList.append((descriptor, function))
            except ModuleNotFoundError:
                raise ImportError("This class requires RDKit to be installed.")
            
        # check initialization
        assert len(self.descriptors) == len(self.descList)
        features = []
        for desc_name, function in self.descList:
            if desc_name == 'Ipc' and self.ipc_avg:
                feature = function(mol, avg=True)
            else:
                feature = function(mol)
            features.append(feature)
        return np.asarray(features)


class InvalidSmilesError(Exception):
    """Exception raised when a SMILES string is chemically invalid."""
    pass

def rd_descriptor_list(list_of_smiles, is_valid_smiles):
    fingerprints = []
    featurizer = RDKitDescriptors()
    #for smiles in list_of_smiles:
    for ind, smiles in enumerate(list_of_smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            fp = featurizer.featurize(mol)
            fp = fp.reshape((208,))
            fingerprints.append(fp)
        else:
            is_valid_smiles[ind] = False
            # As workaround, use 'C' simply so that there will be a valid fingerprint; 
            # is_valid_smiles will be used to flag this entry as false
            mol = Chem.MolFromSmiles('C')   
            fp = featurizer.featurize(mol)
            fp = fp.reshape((208,))
            fingerprints.append(fp)
            
    return fingerprints, is_valid_smiles 


def load(model_name):
    pickle_in = open(picklepath.format(model_name),'rb')
    model = pickle.load(pickle_in)
    return model
#------------------------------------------------------


# -----------------------------------------------------
# Main program
# -----------------------------------------------------

# receive arguments about files from command line input:
print(' ')
cmdline = sys.argv[0:]

if len(cmdline) > 1:
    smiles_file = sys.argv[1]
    outputfile  = sys.argv[2]
else:   # for tests and debugging
    smiles_file = "InputFiles/input_1000_SMILES.txt" 
    outputfile  = "OutputFiles/output_1000_Tg.txt"

# determine the relative path to output folder using the location of this *.py file as the starting point.
# use of the os.path functions is necessary to ensure proper path strings when calling this .py file 
# from a non-local directory, e.g. the AIOMFAC program;
inpfile_path = os.path.abspath(smiles_file)
pyfilepath = os.path.abspath(cmdline[0])
locpath = os.path.dirname(pyfilepath)
outpath = os.path.abspath(locpath +'/OutputFiles')
picklepath = os.path.abspath(locpath + '/pickle/{}')    # absolute path to open the pickle

# load the specific pickle file needed:
model = load('sm_no_tm')

list_of_smiles = []

# open and read input file:
with open(inpfile_path, 'r', newline='') as file1:
    list_of_smiles = [line.strip() for line in file1]

# make sure each array entry is a nonzero string or otherwise remove:
list_of_smiles = [s for s in list_of_smiles if len(s.strip()) > 0]

is_valid_smiles = [True]*len(list_of_smiles) 

# run the prediction model for the whole list of SMILES:
fp, is_valid_smiles = rd_descriptor_list(list_of_smiles, is_valid_smiles)
try: 
    Tg = model.predict(fp)
except:     #shouldn't ever end up here
    print(f"ERROR: at least one SMILES submitted is invalid.")
    Tg = [-99]*len(list_of_smiles)

# output Tg in [K] and SMILES string to file2 (line by line):
with open(outputfile,'w') as file2:
    for index, smiles in enumerate(list_of_smiles):
        if is_valid_smiles[index]:
            file2.write(f'{Tg[index]:.2f}' + '   ' + smiles + '\n')
        else:   # to indicate an issue to the calling program
            file2.write(f'{-99.00:.2f}' + '   ' + smiles + '\n')

print(f'done with processing {len(list_of_smiles)} SMILES... \n')

if (all(is_valid_smiles)):
    print(f'Note: all SMILES were confirmed to be valid.\n')
else:
    print(f'Note: at least one SMILES submitted is invalid; Tg values of -99.00 indicate such invalid SMILES in the output file.\n')
# -----------------------------------------------------