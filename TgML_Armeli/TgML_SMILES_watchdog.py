

''' this script includes code snippets for loading pickle files and returning Tg prediction
authored by Gianluca Armeli'''

''' AZ notes: This script variant was derived from TgML_minimal.py, but adapted (simplified) by Andi Zuend  
    for faster loading of the specific pickle file and a list of SMILES from file for a 
    SMILES-only-input model (model_name = ['sm_no_tm']) prediction. 
    It generates an output file that lists the Tg and corresponding SMILES, one per line.
    Further, this program includes a watchdog event handler so it can be run on a server 
    to monitor and act when a new SMILES input file is added for Tg calculations.'''

print("... importing modules\n")

# modules and libraries for TgML script:
import pickle
import numpy as np
import sys
import os.path
import os
import time
from pathlib import Path
import platform
import re
import logging
from logging.handlers import RotatingFileHandler

from rdkit import Chem
from rdkit.Chem import Descriptors, MolFromSmiles #, AllChem
from rdkit import DataStructs
from deepchem.utils.typing import RDKitMol
from deepchem.feat.base_classes import MolecularFeaturizer

# in addition, for watchdog monitoring script: 
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

print("... end of importing modules\n")


# -----------------------------------------------------
# classes and functions for watchdog event handler
# -----------------------------------------------------
class MyHandler(FileSystemEventHandler):
    # def on_created(self, event):
    #     super().on_created(event)
    #     what = 'directory' if event.is_directory else 'file'
    #     print(f"Created {what}: {event.src_path}")

    def on_created(self, event):
        super().on_created(event)
        # if a new SMILES input file was added, process it to predict Tg:
        if not event.is_directory:
            file_path = event.src_path
            if file_path.lower().endswith('smiles.txt'):
                #print(f"Detected new {'.txt'} file: {os.path.basename(file_path)}")
                n_smiles_proc = TgML_SMILES_prediction(file_path)
                logging.info(f"processed {n_smiles_proc} SMILES from file {os.path.basename(file_path)}")
 
# -----------------------------------------------------


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
            logging.info(f"ERROR: The SMILES '{smiles}' failed RDKit validation.")
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


# -----------------------------------------------------
# Function for TgML calculation on event
# -----------------------------------------------------
def TgML_SMILES_prediction(file_path):

    smiles_file = os.path.basename(file_path)
    outputfile = re.sub(r"input", "output", smiles_file, flags=re.IGNORECASE)
    outfile_path  = './OutputFiles/' + outputfile

    # determine the relative path to output folder using the location of this *.py file as the starting point.
    # use of the os.path functions is necessary to ensure proper path strings when calling this .py file 
    # from a non-local directory, e.g. the AIOMFAC program;
    inpfile_path = './InputFiles/'+ smiles_file

    list_of_smiles = []
    # open and read input file:
    try:
        with open(inpfile_path, 'r', newline='') as file1:
            list_of_smiles = [line.strip() for line in file1]
    except PermissionError:
        logging.info(f"ERROR: no permission to open smiles input file at {inpfile_path}")
        print(f"ERROR: no permission to open smiles input file at {inpfile_path}")
    except FileNotFoundError:
        logging.info(f"ERROR: smiles input file at {inpfile_path} not found")
        print(f"ERROR: smiles input file at {inpfile_path} not found")
    except OSError as e:
        logging.info(f"Cannot open file: {e}")
        print(f"Cannot open file: {e}")

    # make sure each array entry is a nonzero string or otherwise remove:
    list_of_smiles = [s for s in list_of_smiles if len(s.strip()) > 0]

    # remove lines starting with '[##]' since that indicates a remnant line of a 
    # previously appended text from this program:
    list_of_smiles = [item for item in list_of_smiles if not item.startswith('[##]')] 
    is_valid_smiles = [True]*len(list_of_smiles) 

    # run the prediction model for the whole list of SMILES:
    fp, is_valid_smiles = rd_descriptor_list(list_of_smiles, is_valid_smiles)
    try: 
        Tg = model.predict(fp)
    except:
        logging.info(f"ERROR: at least one SMILES submitted is invalid.")
        Tg = [-99]*len(list_of_smiles)


    # output Tg in [K] and SMILES string to file2 (line by line):
    with open(outfile_path,'w') as file2:
        for index, smiles in enumerate(list_of_smiles):
            if is_valid_smiles[index]:
                file2.write(f'{Tg[index]:.2f}' + '   ' + smiles + '\n')
            else:   # to indicate an issue to the calling program
                file2.write(f'{-99.00:.2f}' + '   ' + smiles + '\n')
                
    if platform.system() == "Windows":
        is_linux = False
        # do nothing for now...
    else:
        is_linux = True
        # logging.info(f"OS is Linux.")
        # Linux RHEL: Standard octal for read, write, execute for Owner, Group, and Others
        os.chmod(outfile_path, 0o777)
    
    # for watchdog mode, append a line of text to the input file to increase its size as hint to the 
    # calling (Fortran) program to get the go-ahead for opening the output file.
    with open(inpfile_path, 'a', newline='') as file1:
        file1.write(' \n')
        file1.write(f'[##] NOTE: the input file has been processed and output written to {outfile_path}' + '\n')
              

    return len(list_of_smiles)
# -----------------------------------------------------



# -----------------------------------------------------
# Main watchdog monitoring program
# -----------------------------------------------------
if __name__ == "__main__":
    
    # get and save the process ID in a file as well as a heartbeat file for later lookup 
    # from Fortran program to verify that watchdog program is running:
    pid_file = "../Auxiliary/watchdog.pid"
    with open(pid_file, "w") as f:
        f.write(str(os.getpid()))
    
    # get the directory to monitor from command line arguments
    # if no argument is provided, use the current directory
    monitor_folder_path = sys.argv[1] if len(sys.argv) > 1 else "./InputFiles"

    if not os.path.isdir(monitor_folder_path):
        print(f"Error: The specified path '{monitor_folder_path}' is not a valid directory.")
        sys.exit(1)

    locpath = os.path.dirname(monitor_folder_path)
    picklepath = os.path.abspath(locpath + '/pickle/{}')    # absolute path to open the TgML pickle

    # load the specific pickle file needed:
    model = load('sm_no_tm')

    # configure the logging system:
    log_file = 'logTgMLprog.log'

    handler = RotatingFileHandler(
        log_file, maxBytes=10485760, backupCount=5, encoding="utf-8"
    )

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[handler],
    )

    # configure the watchdog system:
    event_handler = MyHandler()
    observer = Observer()

    # schedule the event handler to watch the specified path recursively
    # recursive=True means it will also watch subdirectories
    observer.schedule(event_handler, monitor_folder_path, recursive=False)

    print(f"[*] Starting to monitor directory: {os.path.abspath(monitor_folder_path)}\n")
    print(f"See the {log_file} file for information on events.\n")
    print("[*] Press Ctrl+C to stop monitoring.\n")
    logging.info(f"Monitoring of {monitor_folder_path} started.\n")
    
    hb_tmp = "../Auxiliary/watchdog.heartbeat.tmp"
    hb  = "../Auxiliary/watchdog.heartbeat"
    
    # start the observer thread
    observer.start()    
    last_heartbeat = 0.0
    try:
        # keep the main thread alive by sleeping
        while True:
            now = time.time()
            # update date/time modified every ~ 5 seconds
            if now - last_heartbeat > 5.0:  
                # update epoch time in .heartbeat file
                with open(hb_tmp, "w") as f:
                    f.write(f"{time.time():.6f}\n")
                    f.flush()
                    os.fsync(f.fileno())
                # try to replace file, do so in a loop in case the Fortran program 
                # happens to access the file at the exact same time
                for _ in range(50):
                    try:
                        os.replace(hb_tmp, hb)
                        break
                    except PermissionError:
                        time.sleep(0.1)
                else:
                    raise RuntimeError("Could not replace file")
                
                last_heartbeat = now
                
            time.sleep(1) 
            
    except KeyboardInterrupt:
        # stop the observer gracefully when Ctrl+C is pressed
        print("\n[*] Stopping monitor...")
        observer.stop()

    # wait for the observer thread to finish
    observer.join()
    print("[*] Monitoring stopped.\n")
    
# ---- end of watchdog script ----